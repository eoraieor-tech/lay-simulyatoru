"""B3-B — üç fazalı mühərrik yüksək doyma təzyiqində (Pb ≥ 240 bar).

PROBLEM (ölçülüb, `ISH_HESABATI.md` → Seans 10). Üç fazalı mühərrik
Pb = 240 və 300 bar-da yığılmırdı ("zaman addımı minimal həddə də
yığılmadı"), halbuki Pb ≤ 220 işləyirdi. B3-A-nın iki fazalı düzəlişi
buna KÖMƏK ETMİRDİ, çünki o, qaz sütunlu cədvələ qəsdən tətbiq
olunmur.

KÖK SƏBƏB. `build_fluid()` Bo-nu HƏR hüceyrədə `oil_fvf(p)` ilə, yəni
cədvəlin DOYMUŞ qolundan oxuyurdu. Doymamış hüceyrədə bu yanlışdır:
orada neftin tərkibi sabitdir (`Rs` sərbəst dəyişəndir, `Rs < Rs_sat(p)`),
ona görə Bo həmin Rs-in doyma təzyiqindən başlayan SIXILMA qoluna
aiddir. İki nəticəsi vardı:

  1. `dBo/dp` doymuş qolda Pb-də sıfırdan keçir → neft tənliyinin
     təzyiq diaqonalı itir, Jakobian kilidlənir;
  2. `∂N_o/∂Rs = 0` — neft tənliyi 3-cü dəyişəndən qopmuşdu, ona görə
     qaz tənliyi 1-ci problemi kompensasiya EDƏ BİLMİRDİ.

HƏLL. Doymamış qol sənaye standartına (Eclipse `PVTO`) uyğun bərpa
olundu:

    Bo(p, Rs) = Bo_sat(Pb(Rs)) · exp(c_o · (Pb(Rs) − p))

`Rs`-in özü TOXUNULMUR — qaz kütlə balansı pozulmur. Bax
`BlackOilPVTProvider.oil_fvf_undersaturated`.
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
from imex2d.domain.scal import GasCoreyParameters
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.implicit.three_phase_state import ThreePhaseState
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter


def _table(bubble_point=240.0, include_gas=True):
    return build_pvt_table(pressure_min=1.0, pressure_max=400.0, n_points=40,
                           bubble_point_bar=bubble_point,
                           include_gas=include_gas)


def _model(bubble_point=240.0, nx=8, ny=8):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=20.0, dy=20.0, dz=10.0,
        porosity=0.22, permx_base=150.0)
    return ReservoirModelBuilder().build(
        geological_model=geology, wells=five_spot(geology.grid),
        scal=default_scal(), gas_scal=GasCoreyParameters(),
        pvt_table=_table(bubble_point), name=f"B3-B Pb={bubble_point}")


def _service():
    return ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(),
        engine_factory=FullyImplicitEngine)


def _engine(model, end_time=400.0):
    return _service().create_engine(model, SimulationConfig(end_time=end_time))


# ═══════════════════════ PVT — doymamış qol ═══════════════════════════

@pytest.mark.parametrize("bubble_point", [200.0, 240.0, 300.0])
def test_branch_is_continuous_at_the_saturation_boundary(bubble_point):
    """Doymuş hüceyrədə Rs = Rs_sat(p) → Pb(Rs) = p → Bo = Bo_sat(p).

    Dəyişən keçidin (variable switching) İŞLƏMƏSİ buna bağlıdır:
    hüceyrə doymuşdan doymamışa keçəndə Bo sıçramamalıdır.
    """
    provider = BlackOilPVTProvider(_table(bubble_point))
    pressure = np.linspace(20.0, 380.0, 40)
    rs_saturated = provider.solution_gor(pressure)

    branch = provider.oil_fvf_undersaturated(pressure, rs_saturated)

    # ⚠️ ŞƏRT DƏYİŞDİ (Seans 36, lövbər düzəlişi). Əvvəl `atol=1e-5` tələb
    # olunurdu və KEÇİRDİ, çünki lövbər cədvəlin ÖZ interpolyasiyası ilə
    # hesablanırdı — yəni test fiziki dəqiqliyi yox, ÖZ-ÖZÜNƏ UYĞUNLUĞU
    # yoxlayırdı. İndi lövbər doymamış qolun fit-indən gəlir (μ üçün 2.27 %
    # meyl aradan qalxdı), ona görə xam sütunla eynilik ARTIQ DÜZGÜN ŞƏRT
    # DEYİL: sütun Pb ətrafında sınığı XƏTTİ kəsir, qol isə hamardır.
    #
    # ÖLÇÜLDÜ: cədvəl DÜYÜNLƏRİNDƏ fərq TAM SIFIRDIR; sapma yalnız son
    # düyünlə Pb arasındakı intervalda yaranır və maksimumu 0.0195-dir.
    # DƏQİQLİK YALNIZ DOYMUŞ ZONADA gözlənilir: Pb-dən yuxarı lövbər
    # qəsdən doymuş meylin davamıdır (törəmə uyğunluğu üçün — bax
    # `_build_saturated_grid`), ona görə orada xam sütunla eynilik şərt deyil.
    nodes = provider.table.pressure
    nodes = nodes[(nodes >= 20.0) & (nodes < provider.table.bubble_point)]
    assert np.allclose(provider.oil_fvf_undersaturated(
        nodes, provider.solution_gor(nodes)), provider.oil_fvf(nodes), rtol=1e-12)
    assert np.allclose(branch, provider.oil_fvf(pressure), atol=2.5e-2)


@pytest.mark.parametrize("bubble_point", [200.0, 240.0, 300.0])
def test_above_the_bubble_point_nothing_changes(bubble_point):
    """GERİYƏ UYĞUNLUQ: Pb-dən yuxarı cədvəlin öz qolu qalır.

    Orada `Rs = Rs_sat(Pb_cədvəl)`, yəni `Pb(Rs) = Pb_cədvəl` və düstur
    cədvəlin qurulduğu düsturun eynisinə çevrilir.
    """
    provider = BlackOilPVTProvider(_table(bubble_point))
    # Pb + 15: cədvəlin Pb-dən sonrakı İLK düyününü keçirik — orada
    # cədvəlin öz xətti interpolyasiyası sınıq nöqtəni kəsir və
    # düsturla deyil, ŞƏBƏKƏ SIXLIĞI ilə bağlı fərq verir.
    pressure = np.linspace(bubble_point + 15.0, 395.0, 20)
    rs_plateau = np.full_like(pressure, float(np.max(
        provider.table.solution_gor)))

    # ⚠️ HƏDD DƏYİŞDİ (Seans 36): lövbər artıq fit-dən gəlir və Pb-dən
    # yuxarı qol HAMAR analitik əyridir, cədvəl sütunu isə düyünlər arasında
    # DÜZ XƏTTDİR. ÖLÇÜLDÜ: maksimal fərq 0.0195 (Pb = 300 halında).
    # Köhnə 1e-5 həddi yalnız ona görə keçirdi ki, lövbər də həmin xətti
    # interpolyasiyadan alınırdı — yəni iki tərəf eyni səhvi bölüşürdü.
    assert np.allclose(provider.oil_fvf_undersaturated(pressure, rs_plateau),
                       provider.oil_fvf(pressure), atol=2.5e-2)


@pytest.mark.parametrize("bubble_point", [240.0, 300.0])
def test_pressure_derivative_never_changes_sign(bubble_point):
    """ƏSL DÜZƏLİŞ: `dBo/dp` doymamış qolda HƏMİŞƏ mənfidir.

    Doymuş qolda o, Pb-də sıfırdan keçirdi — Jakobianın kilidlənməsi
    məhz oradan gəlirdi.
    """
    provider = BlackOilPVTProvider(_table(bubble_point))
    pressure = np.linspace(10.0, 390.0, 60)
    rs = np.full_like(pressure, float(provider.solution_gor(
        np.array([bubble_point * 0.7]))[0]))

    d_dp, _ = provider.oil_fvf_undersaturated_derivatives(pressure, rs)
    assert np.all(d_dp < 0.0), "sıxılma müsbət olmalıdır (dBo/dp < 0)"


def test_saturated_branch_really_does_change_sign():
    """NƏZARƏT: problemin mövcudluğunu təsdiqləyir.

    Bu test uğursuz olsa, yuxarıdakı düzəlişin səbəbi yox olub demək
    olar — yəni cədvəl artıq degenerasiya etmir.
    """
    table = _table(240.0)
    slope = np.diff(table.oil_fvf) / np.diff(table.pressure)
    assert np.any(np.sign(slope[:-1]) != np.sign(slope[1:])), \
        "doymuş qolda dBo/dp işarə dəyişməlidir — problemin mənbəyi budur"


@pytest.mark.parametrize("bubble_point", [240.0, 300.0])
def test_rs_derivative_is_not_zero(bubble_point):
    """`∂Bo/∂Rs ≠ 0` — neft tənliyini 3-cü dəyişənə bağlayan hədd."""
    provider = BlackOilPVTProvider(_table(bubble_point))
    # Rs-in doymuş platosundan AŞAĞI aralıq — doymamış rejim
    rs_max = float(np.max(provider.table.solution_gor))
    rs = np.linspace(0.15 * rs_max, 0.85 * rs_max, 20)
    pressure = np.full_like(rs, bubble_point * 0.9)

    _, d_drs = provider.oil_fvf_undersaturated_derivatives(pressure, rs)
    assert np.all(d_drs > 0.0), "daha çox həll olmuş qaz → daha şişkin neft"


def test_derivatives_match_finite_differences():
    """Analitik törəmələr ↔ mərkəzi fərq."""
    provider = BlackOilPVTProvider(_table(240.0))
    pressure = np.linspace(40.0, 360.0, 25)
    rs = np.full_like(pressure, float(provider.solution_gor(
        np.array([150.0]))[0]))

    step = 1e-4
    numeric = (provider.oil_fvf_undersaturated(pressure + step, rs)
               - provider.oil_fvf_undersaturated(pressure - step, rs)) / (2 * step)
    analytic, _ = provider.oil_fvf_undersaturated_derivatives(pressure, rs)
    assert np.allclose(analytic, numeric, rtol=1e-6)


def test_saturation_pressure_inverts_the_solution_gor():
    """Pb(Rs_sat(p)) = p — tərs DƏQİQ olmalıdır."""
    provider = BlackOilPVTProvider(_table(240.0))
    pressure = np.linspace(20.0, 230.0, 30)
    recovered = provider.saturation_pressure(provider.solution_gor(pressure))
    assert np.allclose(recovered, pressure, atol=1e-6)


# ═══════════════════════ Jakobian ════════════════════════════════════

@pytest.mark.parametrize("bubble_point", [240.0, 300.0])
def test_oil_equation_is_coupled_to_the_third_variable(bubble_point):
    """KİLİDLƏNMƏNİN QARŞISINI ALAN HƏDD: `∂R_neft/∂Rs ≠ 0`.

    Əvvəl bu blok TAM SIFIR idi (`blocks[:, 1, 2] = 0`) — neft tənliyi
    doymamış hüceyrədə 3-cü dəyişəndən ümumiyyətlə asılı deyildi.
    """
    engine = _engine(_model(bubble_point, nx=4, ny=4), end_time=50.0)
    state = engine.state
    assert not state.is_saturated.all(), "sınaq doymamış hüceyrə tələb edir"

    newton = engine.newton
    fluid = newton.build_fluid(state)
    assert fluid.bo_rs is not None
    matrix = newton.jacobian.assemble(state, fluid, 1.0, state.pressure)

    dense = matrix.toarray()
    oil_rows = np.arange(1, state.ncell * 3, 3)
    rs_columns = np.arange(2, state.ncell * 3, 3)
    coupling = np.diag(dense[np.ix_(oil_rows, rs_columns)])
    assert np.all(np.abs(coupling) > 1e-9)


def test_jacobian_matches_finite_differences_when_undersaturated():
    """Analitik Jakobian ↔ tam qalığın sonlu fərqi — AKKUMULYASIYA.

    Axın və quyu blokları ayrıca (`test_three_phase_residual.py`)
    yoxlanılır; burada məqsəd yeni `∂Bo/∂Rs` hədlərinin akkumulyasiyada
    düzgün çıxarıldığını təsdiqləməkdir.
    """
    engine = _engine(_model(240.0, nx=4, ny=4), end_time=50.0)
    newton, state = engine.newton, engine.state
    accumulator = newton.accumulator

    def totals(vector):
        probe = ThreePhaseState.from_vector(vector, state.is_saturated)
        return np.array(accumulator.accumulation(probe, newton.build_fluid(probe)))

    base = state.to_vector()
    analytic = newton.jacobian.accumulation_jacobian.blocks(state, newton.build_fluid(state))

    for cell in range(state.ncell):
        for variable in range(3):
            index = cell * 3 + variable
            step = 1e-6 * max(1.0, abs(base[index]))
            forward, backward = base.copy(), base.copy()
            forward[index] += step
            backward[index] -= step
            numeric = (totals(forward)[:, cell] - totals(backward)[:, cell]) / (2 * step)
            assert np.allclose(analytic[cell, :, variable], numeric,
                               rtol=2e-3, atol=1e-6), \
                f"hüceyrə {cell}, dəyişən {variable}"


# ═══════════════════════ mühərrik — ƏSAS TƏLƏB ════════════════════════

@pytest.mark.parametrize("bubble_point", [240.0, 300.0])
def test_three_phase_engine_converges_at_high_bubble_point(bubble_point):
    """B3-B-nin ƏSL TƏLƏBİ — bu iki dəyər əvvəl yığılmırdı."""
    result = _service().run(_model(bubble_point), SimulationConfig(end_time=400.0))
    assert result.converged, result.message
    assert result.steps > 1
    assert result.series.recovery_factor[-1] > 10.0


@pytest.mark.parametrize("bubble_point", [240.0, 300.0])
def test_gas_is_liberated_at_high_bubble_point(bubble_point):
    """M4 yüksək Pb-də də işləməlidir: sərbəst qaz yaranır, GOR qalxır."""
    result = _service().run(_model(bubble_point), SimulationConfig(end_time=400.0))
    assert result.converged, result.message

    peak_gas = max(float(np.max(snapshot.gas_saturation))
                   for snapshot in result.snapshots
                   if snapshot.gas_saturation is not None)
    assert peak_gas > 1e-3, "doyma təzyiqindən aşağı qaz ayrılmalıdır"

    initial_rs = float(np.mean(_engine(_model(bubble_point)).state.third_variable))
    assert result.series.gas_oil_ratio[-1] > initial_rs


def test_more_dissolved_gas_means_more_free_gas():
    """Fiziki monotonluq: Pb artdıqca ayrılan qaz da artır."""
    peaks = []
    for bubble_point in (200.0, 240.0, 300.0):
        result = _service().run(_model(bubble_point),
                                SimulationConfig(end_time=400.0))
        assert result.converged, result.message
        peaks.append(max(float(np.max(snapshot.gas_saturation))
                         for snapshot in result.snapshots
                         if snapshot.gas_saturation is not None))

    assert peaks[0] < peaks[1] < peaks[2], f"monoton deyil: {peaks}"


def test_low_bubble_point_results_are_unchanged():
    """REQRESSİYA QORUMASI — düzəliş köhnə rejimə TOXUNMAMALIDIR.

    Pb = 100 bar-da heç bir hüceyrə doyma təzyiqinə çatmır, bütün
    hüceyrələr `Rs = Rs_sat(Pb)` platosunda qalır, yəni Bo cədvəlin
    öz doymamış qoludur.

    ⚠️ ETALON DƏYİŞDİ: 62.72 → 62.86 (Seans 17, `max_dt` düzəlişi).
    Köhnə 62.72 SƏHVİN ÖZ DƏYƏRİNİ kilidləyirdi — üç fazalı mühərrik
    istifadəçinin 20 günlük həddini 30-a qaldırırdı, ona görə nəticə
    daha kobud idi. İndi iki fazalı mühərriklə üst-üstə düşür
    (62.861505 vs 62.861429) — bax
    `tests/test_max_timestep_respected.py`.
    """
    result = _service().run(_model(100.0), SimulationConfig(end_time=400.0))
    assert result.converged, result.message
    # QEYD (Seans 36): G2b-nin ilk variantında bu dəyər 62.73-ə düşmüşdü və
    # mən onu "düzgün fizika" kimi izah etmişdim — SƏHV İDİ. Ölçmə göstərdi
    # ki, fərq doyma nöqtəsindəki LÖVBƏR qüsurundan gəlirdi (μ 2.27 % şişik).
    # Lövbər düzəldiləndən sonra dəyər 62.8612-yə qayıtdı, yəni iki və üç
    # fazalı mühərriklərin qazsız rejimdəki uyğunluğu da bərpa olundu.
    assert result.series.recovery_factor[-1] == pytest.approx(62.86, abs=0.05)

    peak_gas = max(float(np.max(snapshot.gas_saturation))
                   for snapshot in result.snapshots
                   if snapshot.gas_saturation is not None)
    assert peak_gas < 1e-6, "Pb = 100 bar-da qaz ayrılmamalıdır"


def test_two_phase_path_is_untouched():
    """Qaz sütunsuz cədvəl hələ də B3-A-nın ölü-neft düzəlişini alır."""
    provider = BlackOilPVTProvider(_table(240.0, include_gas=False),
                                   dead_oil_below_bubble_point=True)
    assert provider.dead_oil_corrected

    gas_provider = BlackOilPVTProvider(_table(240.0, include_gas=True),
                                       dead_oil_below_bubble_point=True)
    assert not gas_provider.dead_oil_corrected, \
        "üç fazalı cədvəl ölü-neftə çevrilməməlidir — qaz balansı pozulardı"
