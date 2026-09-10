"""B3-A — doyma təzyiqindən aşağı Bo düzəlişi.

TAPILMA TARİXÇƏSİ (bax `ISH_HESABATI.md` → Seans 5 və 7): sahibkar
proqramı açıb PVT modelini işə saldı; simulyasiya t = 0 gündə yığılmadı
(RUN-002). Testlərin heç biri bunu tutmurdu, çünki panelin DEFOLT
kombinasiyası (Pb = 240 bar + tam implicit mühərrik) heç bir testdə
yox idi.

KÖK SƏBƏB. Korrelyasiya cədvəlində Bo doyma təzyiqinə qədər ARTIR
(qaz həll olur, neft şişir), ondan yuxarı AZALIR (sıxılma). Neft
tənliyinin təzyiq üzrə diaqonalı `−So·B'o/Bo²`-yə mütənasib olduğuna
görə işarə dəyişən yerdə SIFIRDAN keçir → Jakobian təkləşir → Nyuton
addımı absurd böyüklüyə (ölçüldü: 656 bar) qalxır, qoruyucular onu
kəsir və iterasiya DONUR.

Bu, ədədi qüsur DEYİL: doymuş budaq ayrılan qazı təsvir edir, iki
fazalı model isə o qazı modelləşdirmir. Ona görə düzəliş yalnız QAZ
SÜTUNLARI OLMAYAN cədvələ tətbiq olunur.
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter


def _table(bubble_point=240.0, include_gas=False):
    return build_pvt_table(pressure_min=1.0, pressure_max=400.0, n_points=40,
                           bubble_point_bar=bubble_point,
                           include_gas=include_gas)


def _model(bubble_point=240.0, nx=11, ny=11):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=20.0, dy=20.0, dz=10.0,
        porosity=0.22, permx_base=150.0)
    return ReservoirModelBuilder().build(
        geological_model=geology, wells=five_spot(geology.grid),
        scal=default_scal(), pvt_table=_table(bubble_point),
        name="B3-A sınağı")


def _service():
    return ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(),
        engine_factory=FullyImplicitEngine)


# ═══════════════════════════ cədvəlin özü ═══════════════════════════════
def test_raw_correlation_table_really_is_non_monotone():
    """Problemin ƏSASI — bu olmasaydı düzəlişə ehtiyac qalmazdı."""
    table = _table()
    bo = np.asarray(table.oil_fvf)
    slopes = np.diff(bo) / np.diff(np.asarray(table.pressure))
    assert slopes.min() < 0.0 and slopes.max() > 0.0, \
        "korrelyasiya cədvəlində Bo meyli işarə dəyişməlidir"


def test_correction_makes_bo_non_increasing_in_pressure():
    """Düzəlişdən sonra `dBo/dp ≤ 0` — yəni neft sıxılması MÜSBƏTDİR."""
    table = _table()
    provider = BlackOilPVTProvider(table, dead_oil_below_bubble_point=True)
    assert provider.dead_oil_corrected

    bo = provider.oil_fvf(table.pressure)
    slopes = np.diff(bo) / np.diff(np.asarray(table.pressure))
    assert slopes.max() <= 1e-12, "düzəlişdən sonra Bo artmamalıdır"


def test_provider_stays_a_pure_interpolator_by_default():
    """Defolt davranış DƏYİŞMİR — provider cədvəli olduğu kimi verir."""
    table = _table()
    provider = BlackOilPVTProvider(table)
    assert not provider.dead_oil_corrected
    assert np.allclose(provider.oil_fvf(table.pressure), table.oil_fvf)


def test_gas_table_is_never_touched():
    """Üç fazalı modeldə doymuş budaq DÜZGÜNDÜR — qaz tənliyi onu
    kompensasiya edir, ona görə düzəliş tətbiq OLUNMAMALIDIR."""
    table = _table(include_gas=True)
    provider = BlackOilPVTProvider(table, dead_oil_below_bubble_point=True)
    assert not provider.dead_oil_corrected
    assert np.allclose(provider.oil_fvf(table.pressure), table.oil_fvf)


def test_table_without_a_saturated_branch_is_unchanged():
    """Doyma təzyiqi cədvəlin altındadırsa artan budaq yoxdur."""
    table = build_pvt_table(pressure_min=250.0, pressure_max=400.0,
                            n_points=30, bubble_point_bar=100.0)
    provider = BlackOilPVTProvider(table, dead_oil_below_bubble_point=True)
    assert not provider.dead_oil_corrected


# ═══════════════════════════ mühərrik ═══════════════════════════════════
@pytest.mark.parametrize("bubble_point", [220.0, 240.0, 300.0])
def test_simulation_converges_above_the_old_breaking_point(bubble_point):
    """ƏSAS TEST — bu üç dəyərin HAMISI əvvəl yığılmırdı."""
    result = _service().run(_model(bubble_point),
                            SimulationConfig(end_time=300.0))
    assert result.converged, result.message
    assert result.steps > 0
    assert result.series.recovery_factor[-1] > 1.0


def test_users_actual_failing_case_now_runs():
    """Sahibkarın RUN-002-si: 41×41, panelin defolt Pb = 240 bar.

    Ölçülmüşdü: t = 0.0 gündə "zaman addımı minimal həddə də yığılmadı".
    """
    result = _service().run(_model(240.0, nx=41, ny=41),
                            SimulationConfig(end_time=1500.0))
    assert result.converged, result.message
    assert result.series.recovery_factor[-1] > 5.0


def test_low_bubble_point_result_is_unchanged_bit_for_bit():
    """Heç bir hüceyrə Pb-dən aşağı düşmürsə nəticə DƏYİŞMƏMƏLİDİR.

    Bu, düzəlişin köhnə modelləri sındırmadığının zəmanətidir.
    """
    model = _model(60.0)
    with_fix = _service().run(model, SimulationConfig(end_time=300.0))

    service = _service()
    service.pvt_provider = BlackOilPVTProvider(_table(60.0))   # düzəlişsiz
    engine = FullyImplicitEngine(
        model=_model(60.0), config=SimulationConfig(end_time=300.0),
        relperm=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(),
        pvt=BlackOilPVTProvider(_table(60.0)))
    without_fix = engine.run()

    assert with_fix.converged and without_fix.converged
    assert with_fix.steps == without_fix.steps
    assert np.isclose(with_fix.series.recovery_factor[-1],
                      without_fix.series.recovery_factor[-1], rtol=1e-9)


def test_recovery_grows_with_bubble_point():
    """Fiziki sağlamlıq: doyma təzyiqi yüksəldikcə neft daha çox şişir,
    hasilat da artmalıdır — düzəliş bu meyli pozmamalıdır."""
    recoveries = []
    for bubble_point in (150.0, 240.0, 300.0):
        result = _service().run(_model(bubble_point),
                                SimulationConfig(end_time=300.0))
        assert result.converged, result.message
        recoveries.append(result.series.recovery_factor[-1])
    assert recoveries == sorted(recoveries), recoveries
