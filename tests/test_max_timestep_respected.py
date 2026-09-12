"""Maksimal zaman addımı (`max_dt`) istifadəçinin dediyi kimi işləməlidir.

TAPILAN SƏHV (sahibkarın bildirişi, `ISH_HESABATI.md` → Seans 17).
İnterfeysdə "Maks. Δt = 20 gün" qoyulurdu, mühərrik isə 30 günlük
addımlar atırdı. Səbəb `three_phase_engine._time_config`-də idi:

    max_dt=max(stepping.max_dt, 30.0)      # ← süni döşəmə

Ölçüldü: 5 və 20 gün TAM EYNİ nəticə verirdi (31 addım, maks Δt 30.00),
çünki hər ikisi eyni həddə qaldırılırdı.

⚠️ BU SƏHV BİR DƏFƏ ARTIQ GERİ QAYIDIB. O, əvvəlcə iki fazalı
mühərrikdə tapılıb düzəldilmişdi (bax `engine.py::_time_config`
sənədi), lakin üç fazalı yol A7 bərpasında (B2) git tarixçəsindən
qaytarıldığı üçün köhnə kodu özü ilə geri gətirdi.

BU FAYLIN İŞİ: hər İKİ mühərrik üçün kilid qoymaq ki, üçüncü dəfə
qayıtmasın.
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.config import (SimulationConfig, TimeSteppingConfig)
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.scal import GasCoreyParameters
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter


def _model(with_gas: bool, bubble_point=240.0):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=8, ny=8, dx=20.0, dy=20.0, dz=10.0,
        porosity=0.22, permx_base=150.0)
    return ReservoirModelBuilder().build(
        geological_model=geology, wells=five_spot(geology.grid),
        scal=default_scal(),
        gas_scal=GasCoreyParameters() if with_gas else None,
        pvt_table=build_pvt_table(pressure_min=1.0, pressure_max=400.0,
                                  n_points=40, bubble_point_bar=bubble_point,
                                  include_gas=with_gas),
        name="Δt sınağı")


def _service():
    return ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)


def _config(max_dt: float, end_time=400.0):
    return SimulationConfig(end_time=end_time,
                            time_stepping=TimeSteppingConfig(max_dt=max_dt))


def _steps_taken(result) -> np.ndarray:
    """Faktiki atılan Δt-lər — `series.time`-ın fərqi."""
    time = np.asarray(result.series.time, float)
    return np.diff(np.concatenate(([0.0], time)))


# ═══════════════════ konfiqurasiya səviyyəsi (sürətli) ════════════════

@pytest.mark.parametrize("with_gas", [False, True])
@pytest.mark.parametrize("max_dt", [0.5, 5.0, 20.0, 50.0])
def test_engine_config_passes_max_dt_through_untouched(with_gas, max_dt):
    """`_time_config` dəyəri OLDUĞU KİMİ ötürməlidir — döşəmə YOXDUR."""
    service = _service()
    config = _config(max_dt)
    engine = service.create_engine(_model(with_gas), config)
    assert engine._time_config(config).max_dt == pytest.approx(max_dt)


def test_both_engines_agree_on_the_limit():
    """İki və üç fazalı mühərrik EYNİ həddi qurmalıdır."""
    service = _service()
    config = _config(7.0)
    two = service.create_engine(_model(False), config)._time_config(config)
    three = service.create_engine(_model(True), config)._time_config(config)
    assert two.max_dt == three.max_dt == pytest.approx(7.0)


# ═══════════════════ real qaçış — həddi AŞMIR ════════════════════════

@pytest.mark.parametrize("with_gas", [False, True])
@pytest.mark.parametrize("max_dt", [5.0, 20.0])
def test_no_step_exceeds_the_requested_limit(with_gas, max_dt):
    """ƏSAS TƏLƏB: heç bir addım istifadəçinin həddindən böyük olmamalıdır."""
    result = _service().run(_model(with_gas), _config(max_dt))
    assert result.converged, result.message

    steps = _steps_taken(result)
    assert steps.max() <= max_dt + 1e-9, \
        f"ən böyük addım {steps.max():.4f} > hədd {max_dt}"


def test_smaller_limit_really_produces_smaller_steps():
    """SƏHVİN SİMPTOMU: 5 və 20 gün eyni nəticə verirdi.

    Döşəmə hər ikisini 30-a qaldırdığı üçün fərq YOX İDİ. İndi
    fərqlənməlidirlər.
    """
    fine = _service().run(_model(True), _config(5.0))
    coarse = _service().run(_model(True), _config(20.0))
    assert fine.converged and coarse.converged

    assert _steps_taken(fine).max() < _steps_taken(coarse).max()
    assert fine.steps > coarse.steps, (fine.steps, coarse.steps)


def test_three_phase_no_longer_inflates_small_limits():
    """Üç fazalı mühərrik kiçik həddi 30-a QALDIRMAMALIDIR."""
    result = _service().run(_model(True), _config(5.0))
    assert result.converged, result.message
    assert _steps_taken(result).max() <= 5.0 + 1e-9


# ═══════════════════ fizika: iki və üç fazalı UYĞUNLUĞU ══════════════

def test_three_phase_matches_two_phase_when_no_gas_is_liberated():
    """Qaz ayrılmayan rejimdə iki mühərrik EYNİ nəticə verməlidir.

    Pb = 100 bar, ən aşağı lay təzyiqi ~150 bar → heç bir hüceyrə
    doyma təzyiqinə çatmır, sərbəst qaz yaranmır. Belədə üç fazalı
    mühərrik faktiki olaraq iki fazalı məsələni həll edir.

    ÖLÇÜLMÜŞ: iki fazalı RF 62.861505, üç fazalı RF 62.861429 —
    fərq 7.6·10⁻⁵.

    ⚠️ Bu uyğunluq `max_dt` düzəlişindən ƏVVƏL YOX İDİ: üç fazalı
    30 günlük addım atdığı üçün RF 62.72 çıxırdı. Yəni köhnə etalon
    səhvin öz dəyərini kilidləyirdi.
    """
    config = _config(20.0)
    two = _service().run(_model(False, bubble_point=100.0), config)
    three = _service().run(_model(True, bubble_point=100.0), config)
    assert two.converged and three.converged

    peak_gas = max(float(np.max(snapshot.gas_saturation))
                   for snapshot in three.snapshots
                   if snapshot.gas_saturation is not None)
    assert peak_gas < 1e-9, "bu sınaq qazsız rejim tələb edir"

    assert three.series.recovery_factor[-1] == pytest.approx(
        two.series.recovery_factor[-1], abs=1e-3)
    assert three.steps == two.steps
