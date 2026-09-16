"""G2 (2-ci hissə) — doymamış özlülüyün MÜHƏRRİYƏ qoşulması.

G2a provider-i genişləndirdi, bu addım isə onu flüid vəziyyətinə və
Jakobiana bağlayır. Neft mobilliyi `kro/(μo·Bo)`-dur, ona görə doymamış
hüceyrədə 3-cü dəyişənə (Rs) görə törəmə İKİ hədddən ibarətdir:

    ∂mob_o/∂Rs = −mob_o · ( Bo'_Rs/Bo  +  μo'_Rs/μo )

İkinci hədd bu commit-də əlavə olundu. Testlər həm DÜZGÜNLÜYÜNÜ (sonlu
fərq), həm də ZƏRURİLİYİNİ (onsuz uyğunluq pozulur) yoxlayır.
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder, five_spot
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.scal import CoreyParameters, GasCoreyParameters
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.implicit.three_phase_state import ThreePhaseState
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter

WATER_OIL = CoreyParameters()
GAS = GasCoreyParameters()


def _engine(end_time: float = 60.0, nx: int = 4, ny: int = 4):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    model = ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=WATER_OIL, gas_scal=GAS,
        pvt_table=build_pvt_table(pressure_min=1.0, pressure_max=400.0,
                                  n_points=40, bubble_point_bar=240.0,
                                  include_gas=True),
        name="G2b sınağı")
    service = ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)
    return service.create_engine(model, SimulationConfig(end_time=end_time))


def _mixed_state(engine, pressure_value: float = 300.0):
    """Yarısı doymuş (3-cü dəyişən Sg), yarısı doymamış (Rs)."""
    n = engine.model.ncell
    pressure = np.full(n, pressure_value)
    saturated = np.zeros(n, bool)
    saturated[: n // 2] = True
    rs_sat = np.asarray(engine.newton.pvt.solution_gor(pressure))
    third = np.where(saturated, 0.08, 0.6 * rs_sat)
    return ThreePhaseState(pressure, np.full(n, 0.3), third, saturated)


def _saturated_state(engine, pressure_value: float = 200.0):
    n = engine.model.ncell
    return ThreePhaseState(np.full(n, pressure_value), np.full(n, 0.3),
                           np.full(n, 0.08), np.ones(n, bool))


def _column_kind_errors(engine, state, zero_viscosity_rs=False):
    """Sütun NÖVÜ üzrə (p / Sw / 3-cü) maksimal nisbi xəta.

    NİYƏ NÖV ÜZRƏ: ölçüldü ki, TƏZYİQ sütununda bu modeldə ~0.19 səviyyəsində
    MÖVCUD (G2b-dən əvvəlki) qeyri-dəqiqlik var — ümumi maksimum onunla
    üstələnir və yeni həddin təsirini gizlədir. Bax modul sənədləşməsi.
    """
    errors = {}
    for kind in range(3):
        columns = range(kind, state.ncell * 3, 3)
        errors[("p", "Sw", "third")[kind]] = _jacobian_error(
            engine, state, columns=columns, zero_viscosity_rs=zero_viscosity_rs)
    return errors


def _jacobian_error(engine, state, columns=None, zero_viscosity_rs=False):
    """Yığılmış Jakobianın tam qalığın SONLU FƏRQİ ilə uyğunluğu.

    `zero_viscosity_rs` — yeni həddi SÜNİ ŞƏKİLDƏ söndürür; onda uyğunluq
    pozulmalıdır (həddin zəruriliyinin sübutu).
    """
    newton, dt = engine.newton, 1.0
    previous = state.copy()
    previous_fluid = newton.build_fluid(previous)

    def residual(s):
        return newton.compute_residual(s, previous, previous_fluid, dt)[0]

    fluid = newton.build_fluid(state)
    if zero_viscosity_rs:
        fluid.mu_o_rs = np.zeros_like(fluid.mu_o_rs)
    jacobian = newton.jacobian.assemble(state, fluid, dt)

    vector = state.to_vector()
    worst = 0.0
    for column in (columns if columns is not None else range(len(vector))):
        step = 1e-6 * max(1.0, abs(vector[column]))
        forward, backward = vector.copy(), vector.copy()
        forward[column] += step
        backward[column] -= step
        numeric = (residual(ThreePhaseState.from_vector(forward, state.is_saturated))
                   - residual(ThreePhaseState.from_vector(backward, state.is_saturated))
                   ) / (2 * step)
        analytic = np.asarray(jacobian[:, column].todense()).ravel()
        worst = max(worst, np.max(np.abs(numeric - analytic))
                    / max(1.0, np.max(np.abs(numeric))))
    return worst


# ═══════════════════════ flüid vəziyyəti ═════════════════════════════

def test_fluid_state_carries_the_viscosity_derivatives():
    engine = _engine()
    state = _mixed_state(engine)
    fluid = engine.newton.build_fluid(state)
    assert fluid.mu_o_p is not None and fluid.mu_o_rs is not None


def test_saturated_cells_keep_the_saturated_branch():
    engine = _engine()
    state = _mixed_state(engine)
    fluid = engine.newton.build_fluid(state)
    saturated = state.is_saturated
    mu_sat = np.asarray(engine.newton.pvt.oil_viscosity(state.pressure))
    assert np.allclose(fluid.mu_o[saturated], mu_sat[saturated])
    assert np.all(fluid.mu_o_rs[saturated] == 0.0), "doymuşda 3-cü dəyişən Sg-dir"


def test_undersaturated_cells_are_thicker_than_the_saturated_branch():
    """ÖLÇÜLDÜ: doymuş qoldan oxumaq özlülüyü ~1.6 dəfə AZ göstərirdi."""
    engine = _engine()
    state = _mixed_state(engine)
    fluid = engine.newton.build_fluid(state)
    under = ~state.is_saturated
    mu_sat = np.asarray(engine.newton.pvt.oil_viscosity(state.pressure))
    ratio = fluid.mu_o[under] / mu_sat[under]
    assert np.all(ratio > 1.3), ratio[:3]
    assert np.all(np.abs(fluid.mu_o_rs[under]) > 0.0)


# ═══════════════════════ Jakobian — sonlu fərq ═══════════════════════

def test_saturation_columns_match_finite_difference_in_a_mixed_state():
    """ƏSAS YOXLAMA: Sw və 3-cü dəyişən sütunları qalıqla UYĞUNDUR.

    ÖLÇÜLDÜ: hər ikisi ~10⁻¹⁰ səviyyəsindədir, yəni analitik törəmələr
    praktik olaraq dəqiqdir.
    """
    engine = _engine()
    errors = _column_kind_errors(engine, _mixed_state(engine))
    assert errors["Sw"] < 1e-8, errors
    assert errors["third"] < 1e-8, errors


def test_pressure_column_inaccuracy_is_pre_existing_not_from_this_change():
    """TƏZYİQ sütunundakı qeyri-dəqiqlik G2b-DƏN GƏLMİR.

    ÖLÇÜLDÜ (eyni vəziyyət, eyni model):

        G2b (yeni)                        0.1936
        yalnız ∂μo/∂p köhnəyə qaytarılmış 0.1936   ← eynidir
        tam köhnə davranış (G2b-dən əvvəl) 0.5672   ← DAHA PİS

    YERLƏŞDİRİLDİ: ən pis element 8-ci hüceyrənin QAZ tənliyindədir və
    həmin hüceyrədə QUYU YOXDUR — yəni mənbə quyu həddi deyil, AXIN
    Jakobianıdır (ehtimal: cazibə həddində sıxlığın təzyiqdən asılılığının
    və upstream seçiminin diferensiallaşdırılmaması — hər ikisi sənədlərdə
    açıq göstərilmiş mövcud sadələşdirmələrdir).

    Yəni bu qeyri-dəqiqlik G2b-dən ƏVVƏLKİ texniki borcdur və G2b onu
    YAXŞILAŞDIRIB. Test bunu kilidləyir ki, gələcəkdə kimsə onu G2b-nin
    üstünə yazmasın.
    """
    engine = _engine()
    state = _mixed_state(engine)
    errors = _column_kind_errors(engine, state)
    same_without_term = _column_kind_errors(engine, state, zero_viscosity_rs=True)
    assert errors["p"] == pytest.approx(same_without_term["p"], rel=1e-9)
    assert errors["p"] < 1.0


def test_the_new_viscosity_term_is_necessary():
    """ZƏRURİLİK SÜBUTU: yeni həddi söndürəndə 3-cü sütun DAĞILIR.

    ÖLÇÜLDÜ: hədd ilə 3.6×10⁻¹¹, hədd olmadan 4.2×10⁻¹ — yəni on milyard
    dəfə fərq. Müqayisə YALNIZ 3-cü sütun üzrədir, çünki ümumi maksimum
    təzyiq sütunundakı mövcud qeyri-dəqiqliklə üstələnir.
    """
    engine = _engine()
    state = _mixed_state(engine)
    with_term = _column_kind_errors(engine, state)["third"]
    without_term = _column_kind_errors(engine, state,
                                       zero_viscosity_rs=True)["third"]
    assert without_term > with_term * 1000.0, (with_term, without_term)


def test_saturated_only_state_is_unaffected():
    """GERİYƏ UYĞUNLUQ: hamısı doymuş olanda yeni hədd SIFIRDIR."""
    engine = _engine()
    state = _saturated_state(engine)
    fluid = engine.newton.build_fluid(state)
    assert np.all(fluid.mu_o_rs == 0.0)
    errors = _column_kind_errors(engine, state)
    assert errors["Sw"] < 1e-6 and errors["third"] < 1e-6, errors


# ═══════════════════════ uc-uca ═══════════════════════════════════════

def test_run_converges_with_the_undersaturated_branch():
    engine = _engine(end_time=120.0)
    result = engine.run()
    assert result.converged, result.message


def test_engine_now_reads_the_undersaturated_viscosity():
    """G2a-dakı staging testinin TƏRSİ: qoşulma artıq baş verib."""
    import inspect
    from imex2d.simulation.implicit import three_phase_newton
    source = inspect.getsource(three_phase_newton.ThreePhaseNewtonSolver)
    assert "oil_viscosity_undersaturated" in source
    assert "mu_o_p=mu_o_p, mu_o_rs=mu_o_rs" in source
