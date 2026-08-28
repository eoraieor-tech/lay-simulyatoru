"""Quyu tənlikləri — OPM tipli quyu modeli, MƏRHƏLƏ 2-4-6.

v69 addım 4b: testlər İKİ FAZALI (neft-su) modelə uyğunlaşdırıldı.
Qaza aid testlər (sərbəst+həll olmuş qaz, üçüncü dəyişən x) silindi;
doymuş/doymamış rejim cütlüyü isə İKİ FƏRQLİ DOYUMLULUQ (aşağı/yuxarı
Sw) ilə əvəz olundu — metodoloji qayda qalır: törəmələr YALNIZ bir
rejimdə deyil, mühərrikin işlədiyi hər rejimdə yoxlanılmalıdır
(bax `A7_PLAN.md`, v60 dərsi).

Debitlər BHP naməlumundan hesablanır; ən vacib xüsusiyyət —
KƏSMƏ YOXDUR, ona görə funksiya sıfır-debit nöqtəsində də hamardır.
"""

import numpy as np

from helpers import default_scal
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.domain.scal import CoreyParameters
from imex2d.domain.wells import ControlMode
from imex2d.simulation.implicit.residual import FluidState
from imex2d.simulation.implicit.standard_well import StandardWellModel
from imex2d.simulation.implicit.state import ReservoirState
from imex2d.simulation.implicit.well_state import WellUnknowns
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.well_model import PeacemanWellModel


def _providers():
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0))
    relperm = CoreyRelativePermeabilityAdapter(CoreyParameters())
    return pvt, relperm


def _fluid_state(pvt, relperm, state: ReservoirState) -> FluidState:
    p, sw = state.pressure, state.water_saturation
    mu_w = np.asarray(pvt.water_viscosity(p), float)
    mu_o = np.asarray(pvt.oil_viscosity(p), float)
    return FluidState(
        mu_w=mu_w, mu_o=mu_o,
        bw=np.asarray(pvt.water_fvf(p), float),
        bo=np.asarray(pvt.oil_fvf(p), float),
        lam_w=np.asarray(relperm.krw(sw), float) / mu_w,
        lam_o=np.asarray(relperm.kro(sw), float) / mu_o)


def _setup(nx=5, ny=5, pressure=213.5, sw=0.35):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    model = ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                          scal=default_scal())
    connections = PeacemanWellModel().build_connections(model)
    n = model.ncell

    pvt, relperm = _providers()
    p = np.full(n, pressure)
    state = ReservoirState(p, np.full(n, sw))
    fluid = _fluid_state(pvt, relperm, state)
    model_wells = StandardWellModel(connections, n)
    wells = WellUnknowns.from_connections(connections, p)
    return model, connections, state, fluid, model_wells, wells


def _with_bhp(wells, name, value):
    bhp = wells.bhp.copy()
    bhp[wells.index_of(name)] = value
    return WellUnknowns(list(wells.names), bhp)


# ── HAMARLIQ — bu mərhələnin əsas məqsədi ───────────────────────────
def test_rate_passes_through_zero_smoothly():
    """Sıfır-debit nöqtəsində SINIQ OLMAMALIDIR.

    Köhnə modeldə `min(q, 0)` kəsməsi burada törəməni 0-dan 1-ə
    sıçradırdı və Nyuton osilyasiya edirdi (bax `A7_PLAN.md`).
    Yeni modeldə debit sadəcə xəttidir.
    """
    model, connections, state, fluid, well_model, wells = _setup()
    cell_pressure = state.pressure[0]

    rates = []
    for offset in (-2.0, -1.0, 0.0, 1.0, 2.0):
        modified = _with_bhp(wells, "PROD-1", cell_pressure + offset)
        result = well_model.perforation_rates(state, fluid, modified)
        rates.append(result.per_well_oil["PROD-1"])

    # bərabər addımlarda debit də BƏRABƏR artmalıdır (xəttilik)
    differences = np.diff(rates)
    assert np.allclose(differences, differences[0], rtol=1e-6)


def test_rate_is_exactly_zero_when_bhp_equals_cell_pressure():
    model, connections, state, fluid, well_model, wells = _setup()
    modified = _with_bhp(wells, "PROD-1", state.pressure[0])
    rates = well_model.perforation_rates(state, fluid, modified)
    assert abs(rates.per_well_oil["PROD-1"]) < 1e-12


def test_no_clamping_allows_sign_reversal():
    """Kəsmə olmadığı üçün istismarçı BHP-si hüceyrədən yuxarı olanda
    debit MÜSBƏTƏ keçir (çarpaz axın) — OPM də belə modelləşdirir."""
    model, connections, state, fluid, well_model, wells = _setup()
    below = _with_bhp(wells, "PROD-1", state.pressure[0] - 20.0)
    above = _with_bhp(wells, "PROD-1", state.pressure[0] + 20.0)
    assert well_model.perforation_rates(state, fluid, below
                                        ).per_well_oil["PROD-1"] < 0
    assert well_model.perforation_rates(state, fluid, above
                                        ).per_well_oil["PROD-1"] > 0


# ── işarə konvensiyası ──────────────────────────────────────────────
def test_producer_rates_are_negative_at_its_target():
    model, connections, state, fluid, well_model, wells = _setup()
    rates = well_model.perforation_rates(state, fluid, wells)
    assert rates.per_well_oil["PROD-1"] < 0
    assert rates.per_well_water["PROD-1"] < 0


def test_injector_rate_is_positive_at_its_target():
    model, connections, state, fluid, well_model, wells = _setup()
    rates = well_model.perforation_rates(state, fluid, wells)
    assert rates.per_well_water["INJ-1"] > 0


def test_water_injector_produces_no_oil():
    """Su vurucusu yalnız su verir — vurulan fazanın mobilliyi."""
    model, connections, state, fluid, well_model, wells = _setup()
    rates = well_model.perforation_rates(state, fluid, wells)
    assert abs(rates.per_well_oil["INJ-1"]) < 1e-12


# ── mobilliyin upstream seçimi ──────────────────────────────────────
def test_injector_switches_to_reservoir_mobility_on_cross_flow():
    """Vurucunun BHP-si hüceyrədən AŞAĞI olsa, axın laydan çıxır —
    onda hüceyrənin öz mobillikləri işlədilməlidir (yalnız su yox)."""
    model, connections, state, fluid, well_model, wells = _setup()
    reversed_ = _with_bhp(wells, "INJ-1", state.pressure[0] - 30.0)
    rates = well_model.perforation_rates(state, fluid, reversed_)
    assert rates.per_well_water["INJ-1"] < 0
    assert rates.per_well_oil["INJ-1"] < 0        # neft də çıxır


# ── idarəetmə tənlikləri ────────────────────────────────────────────
def test_bhp_control_residual_is_zero_at_the_target():
    model, connections, state, fluid, well_model, wells = _setup()
    rates = well_model.perforation_rates(state, fluid, wells)
    residuals = well_model.control_residuals(rates, wells)
    assert np.allclose(residuals, 0.0)


def test_bhp_control_residual_measures_the_offset():
    model, connections, state, fluid, well_model, wells = _setup()
    offset = _with_bhp(wells, "PROD-1", 150.0 + 7.0)
    rates = well_model.perforation_rates(state, fluid, offset)
    residuals = well_model.control_residuals(rates, offset)
    assert abs(residuals[offset.index_of("PROD-1")] - 7.0) < 1e-9


def test_rate_control_residual_uses_liquid_rate():
    """RATE hədəfi maye (su+neft) debitidir — A6 konvensiyası."""
    from imex2d.simulation.well_model import WellConnection

    model, connections, state, fluid, _, _ = _setup()
    producer = next(c for c in connections if not c.is_injector)
    rate_connection = WellConnection(
        well_name=producer.well_name, cell=producer.cell,
        well_index=producer.well_index, is_injector=False,
        mode=ControlMode.RATE, target=50.0)
    well_model = StandardWellModel([rate_connection], model.ncell)
    wells = WellUnknowns.from_connections([rate_connection], state.pressure)

    rates = well_model.perforation_rates(state, fluid, wells)
    residuals = well_model.control_residuals(rates, wells)
    liquid = (rates.per_well_water[producer.well_name]
              + rates.per_well_oil[producer.well_name])
    expected = (liquid - (-50.0)) / producer.well_index
    assert abs(residuals[0] - expected) < 1e-9


def test_rate_control_residual_is_scaled_to_pressure_units():
    """RATE qalığı m³/gün-dədir, BHP qalığı bar — eyni vektorda çox
    fərqli böyüklüklər pis şərtlənmə yaradır, ona görə miqyaslanır."""
    from imex2d.simulation.well_model import WellConnection

    model, connections, state, fluid, _, _ = _setup()
    producer = next(c for c in connections if not c.is_injector)
    rate_connection = WellConnection(
        well_name=producer.well_name, cell=producer.cell,
        well_index=producer.well_index, is_injector=False,
        mode=ControlMode.RATE, target=50.0)
    well_model = StandardWellModel([rate_connection], model.ncell)
    wells = WellUnknowns.from_connections([rate_connection], state.pressure)
    rates = well_model.perforation_rates(state, fluid, wells)
    residuals = well_model.control_residuals(rates, wells)
    # miqyaslanmamış qalıq yüzlərlə olardı; miqyaslanmış onlarla
    assert abs(residuals[0]) < 500.0


def test_control_residual_count_matches_well_count():
    model, connections, state, fluid, well_model, wells = _setup()
    rates = well_model.perforation_rates(state, fluid, wells)
    assert well_model.control_residuals(rates, wells).size == wells.count


def test_is_bhp_controlled_reports_the_mode():
    model, connections, state, fluid, well_model, wells = _setup()
    assert well_model.is_bhp_controlled("PROD-1") is True


# ── cəm uyğunluğu ───────────────────────────────────────────────────
def test_per_well_totals_match_the_cell_arrays():
    model, connections, state, fluid, well_model, wells = _setup()
    rates = well_model.perforation_rates(state, fluid, wells)
    for name in wells.names:
        cells = [c.cell for c in connections if c.well_name == name]
        assert abs(rates.per_well_oil[name]
                   - sum(rates.oil[c] for c in cells)) < 1e-9


# ══════════════════════════════════════════════════════════════════
# MƏRHƏLƏ 3 — Jakobian
# ══════════════════════════════════════════════════════════════════

def _jacobian_setup(sw: float):
    """Törəmələr İKİ FƏRQLİ doyumluluqda yoxlanılır.

    DƏRS (A7_PLAN.md, v60): əvvəllər Jakobian yalnız BİR rejimdə
    doğrulanırdı və digər rejimdəki səhv uzun müddət gizli qaldı.
    İki fazalı modeldə "rejim" doyumluluq səviyyəsidir: aşağı Sw-də
    kro böyük/krw kiçik, yuxarı Sw-də əksinə — törəmələrin miqyası
    tamamilə fərqlidir.
    """
    from imex2d.simulation.implicit.standard_well import StandardWellJacobian

    model, connections, base, _, well_model, wells = _setup()
    n = model.ncell
    pvt, relperm = _providers()
    state = ReservoirState(base.pressure, np.full(n, sw))

    def build_fluid(s):
        return _fluid_state(pvt, relperm, s)

    jacobian = StandardWellJacobian(well_model, pvt, relperm)
    return (model, connections, well_model, wells, state, build_fluid,
            jacobian)


def _check_reservoir_column(sw: float, variable: int, tolerance=1e-6):
    (model, connections, well_model, wells, state, build_fluid,
     jacobian) = _jacobian_setup(sw)
    cell = next(c for c in connections if not c.is_injector).cell
    blocks = jacobian.blocks(state, build_fluid(state), wells)

    step = 1e-5
    vector = state.to_vector()
    forward = vector.copy(); forward[cell * 2 + variable] += step
    backward = vector.copy(); backward[cell * 2 + variable] -= step
    state_f = ReservoirState.from_vector(forward)
    state_b = ReservoirState.from_vector(backward)

    rates_f = well_model.perforation_rates(state_f, build_fluid(state_f), wells)
    rates_b = well_model.perforation_rates(state_b, build_fluid(state_b), wells)
    numeric = np.array([
        (rates_f.water[cell] - rates_b.water[cell]) / (2 * step),
        (rates_f.oil[cell] - rates_b.oil[cell]) / (2 * step)])
    analytic = blocks.rate_wrt_reservoir[cell, :, variable]
    error = np.max(np.abs(numeric - analytic)) / max(1.0, np.abs(numeric).max())
    assert error < tolerance, f"xəta {error}"


def test_rate_pressure_derivative_at_low_saturation():
    _check_reservoir_column(sw=0.30, variable=0)


def test_rate_pressure_derivative_at_high_saturation():
    _check_reservoir_column(sw=0.60, variable=0)


def test_rate_water_saturation_derivative_at_low_saturation():
    _check_reservoir_column(sw=0.30, variable=1)


def test_rate_water_saturation_derivative_at_high_saturation():
    _check_reservoir_column(sw=0.60, variable=1)


def _check_bhp_column(sw: float):
    (model, connections, well_model, wells, state, build_fluid,
     jacobian) = _jacobian_setup(sw)
    producer = next(c for c in connections if not c.is_injector)
    cell = producer.cell
    fluid = build_fluid(state)
    blocks = jacobian.blocks(state, fluid, wells)

    step = 1e-5
    position = wells.index_of(producer.well_name)
    up = wells.bhp.copy(); up[position] += step
    down = wells.bhp.copy(); down[position] -= step
    rates_f = well_model.perforation_rates(
        state, fluid, WellUnknowns(list(wells.names), up))
    rates_b = well_model.perforation_rates(
        state, fluid, WellUnknowns(list(wells.names), down))
    numeric = np.array([
        (rates_f.water[cell] - rates_b.water[cell]) / (2 * step),
        (rates_f.oil[cell] - rates_b.oil[cell]) / (2 * step)])
    error = (np.max(np.abs(numeric - blocks.rate_wrt_bhp[cell]))
             / max(1.0, np.abs(numeric).max()))
    assert error < 1e-6, f"xəta {error}"


def test_rate_bhp_derivative_at_low_saturation():
    """`∂q/∂p_bhp` — köhnə modeldə YOX idi (BHP orada sabit idi).
    Məhz bu bağlantı quyunu sistemin bir hissəsinə çevirir."""
    _check_bhp_column(sw=0.30)


def test_rate_bhp_derivative_at_high_saturation():
    _check_bhp_column(sw=0.60)


def test_bhp_control_derivative_is_one():
    """BHP idarəsində `R = p_bhp − hədəf` → ∂R/∂p_bhp = 1."""
    (model, connections, well_model, wells, state, build_fluid,
     jacobian) = _jacobian_setup(sw=0.35)
    blocks = jacobian.blocks(state, build_fluid(state), wells)
    assert np.allclose(blocks.control_wrt_bhp, 1.0)


def test_bhp_control_has_no_reservoir_coupling():
    """BHP idarəsində idarəetmə tənliyi rezervuardan ASILI DEYİL."""
    (model, connections, well_model, wells, state, build_fluid,
     jacobian) = _jacobian_setup(sw=0.35)
    blocks = jacobian.blocks(state, build_fluid(state), wells)
    for coupling in blocks.control_wrt_reservoir.values():
        assert coupling == {}


def test_rate_control_derivative_matches_finite_difference():
    """RATE idarəsində idarəetmə tənliyi HƏM BHP-dən, HƏM rezervuardan
    asılıdır — hər iki bağlantı yoxlanılır."""
    from imex2d.simulation.implicit.standard_well import StandardWellJacobian
    from imex2d.simulation.well_model import WellConnection

    model, connections, base, _, _, _ = _setup()
    n = model.ncell
    producer = next(c for c in connections if not c.is_injector)
    rate_connection = WellConnection(
        well_name=producer.well_name, cell=producer.cell,
        well_index=producer.well_index, is_injector=False,
        mode=ControlMode.RATE, target=50.0)
    well_model = StandardWellModel([rate_connection], n)
    wells = WellUnknowns.from_connections([rate_connection], base.pressure)

    pvt, relperm = _providers()
    state = ReservoirState(base.pressure, base.water_saturation)
    fluid = _fluid_state(pvt, relperm, state)

    jacobian = StandardWellJacobian(well_model, pvt, relperm)
    blocks = jacobian.blocks(state, fluid, wells)

    step = 1e-5
    up = wells.bhp.copy(); up[0] += step
    down = wells.bhp.copy(); down[0] -= step
    residual_f = well_model.control_residuals(
        well_model.perforation_rates(state, fluid,
                                     WellUnknowns(list(wells.names), up)),
        WellUnknowns(list(wells.names), up))
    residual_b = well_model.control_residuals(
        well_model.perforation_rates(state, fluid,
                                     WellUnknowns(list(wells.names), down)),
        WellUnknowns(list(wells.names), down))
    numeric = (residual_f[0] - residual_b[0]) / (2 * step)
    assert abs(numeric - blocks.control_wrt_bhp[0]) < 1e-6 * max(
        1.0, abs(numeric))


# ══════════════════════════════════════════════════════════════════
# MƏRHƏLƏ 4 — birləşmiş Nyuton
# ══════════════════════════════════════════════════════════════════

def _coupled_setup(nx=4, ny=4):
    from imex2d.simulation.discretization import TwoPointFluxDiscretization
    from imex2d.simulation.implicit.coupled_newton import CoupledNewtonSolver
    from imex2d.simulation.implicit.linear import NewtonLinearSolver
    from imex2d.simulation.implicit.well_state import CoupledState

    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    model = ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                          scal=default_scal())
    grid = TwoPointFluxDiscretization().build(model)
    connections = PeacemanWellModel().build_connections(model)
    pvt, relperm = _providers()
    solver = CoupledNewtonSolver(model, relperm, pvt, NewtonLinearSolver(),
                                 grid, connections)
    n = model.ncell
    reservoir = ReservoirState(np.linspace(200.0, 225.0, n), np.full(n, 0.35))
    state = CoupledState(reservoir,
                         WellUnknowns.from_connections(connections,
                                                       reservoir.pressure))
    return solver, state


def test_coupled_system_size_includes_wells():
    solver, state = _coupled_setup()
    residual, _, _ = solver.compute_residual(
        state, state, solver.build_fluid(state.reservoir), 1.0)
    assert residual.size == state.size
    assert state.size == state.reservoir.ncell * 2 + state.wells.count


def test_coupled_jacobian_matches_finite_difference():
    """BİRLƏŞMİŞ Jakobian — rezervuar, quyu VƏ onların qarşılıqlı
    bloklarını birlikdə yoxlayır. Ayrı-ayrı hissələrin doğru olması
    onların düzgün BİRLƏŞDİYİNİ təmin etmir."""
    from imex2d.simulation.implicit.well_state import CoupledState

    solver, state = _coupled_setup()
    previous = state.copy()
    previous_fluid = solver.build_fluid(previous.reservoir)
    dt = 1.0

    def residual_of(candidate):
        value, _, _ = solver.compute_residual(candidate, previous,
                                              previous_fluid, dt)
        return value

    _, fluid, _ = solver.compute_residual(state, previous, previous_fluid, dt)
    jacobian = solver.assemble_jacobian(state, fluid, dt)

    step = 1e-6
    vector = state.to_vector()
    names = state.wells.names
    worst = 0.0
    for column in range(vector.size):
        forward = vector.copy(); forward[column] += step
        backward = vector.copy(); backward[column] -= step
        numeric = (residual_of(CoupledState.from_vector(forward, names))
                   - residual_of(CoupledState.from_vector(backward, names))
                   ) / (2 * step)
        analytic = np.asarray(jacobian[:, column].todense()).ravel()
        worst = max(worst, np.max(np.abs(numeric - analytic))
                    / max(1.0, np.abs(numeric).max()))
    assert worst < 1e-5, f"birləşmiş Jakobian xətası: {worst}"


def test_coupled_solver_converges_on_a_single_step():
    solver, state = _coupled_setup()
    result = solver.solve(state, dt=1.0)
    assert result.converged, result.status


def test_bhp_controlled_wells_stay_at_their_target():
    """BHP idarəsində quyu tənliyi `p_bhp = hədəf` deməkdir — həll
    sonunda BHP dəqiq hədəfdə olmalıdır."""
    solver, state = _coupled_setup()
    original = state.wells.bhp.copy()
    result = solver.solve(state, dt=1.0)
    assert result.converged
    assert np.allclose(result.state.wells.bhp, original, atol=1e-6)


def test_coupled_solver_never_raises():
    """Mövcud həllediciylə eyni zəmanət (bax `NewtonSolver`)."""
    solver, state = _coupled_setup()

    class Broken:
        def __getattr__(self, name):
            raise RuntimeError("qəsdən sındırılmış")

    solver.reservoir = Broken()
    result = solver.solve(state, dt=1.0)
    assert result is not None
    assert not result.converged


# ══════════════════════════════════════════════════════════════════
# MƏRHƏLƏ 6 — zaman addımının quyu debitinə görə hədd
# ══════════════════════════════════════════════════════════════════

def _stepper_setup():
    from imex2d.simulation.discretization import TwoPointFluxDiscretization
    from imex2d.simulation.implicit.coupled_newton import CoupledNewtonSolver
    from imex2d.simulation.implicit.linear import NewtonLinearSolver

    model, connections, reservoir, fluid, well_model, wells = _setup()
    grid = TwoPointFluxDiscretization().build(model)
    pvt, relperm = _providers()
    solver = CoupledNewtonSolver(model, relperm, pvt, NewtonLinearSolver(),
                                 grid, connections)
    return solver, model, connections, reservoir, wells


def test_max_stable_dt_is_bounded_by_the_fastest_draining_cell():
    """1250 m³ məsamə həcmi, 654 m³/gün debit → dt=1.0 hüceyrənin
    yarısını çıxarır BİR ADDIMDA (kök səbəb, bax A7_PLAN.md)."""
    from imex2d.simulation.implicit.well_state import CoupledState

    solver, model, connections, reservoir, wells = _stepper_setup()
    state = CoupledState(reservoir, wells)
    limit = solver.max_stable_dt(state, alpha=0.2)
    assert 0.0 < limit < 1.0, f"gözlənilməz hədd: {limit}"


def test_max_stable_dt_scales_with_alpha():
    from imex2d.simulation.implicit.well_state import CoupledState

    solver, model, connections, reservoir, wells = _stepper_setup()
    state = CoupledState(reservoir, wells)
    narrow = solver.max_stable_dt(state, alpha=0.1)
    wide = solver.max_stable_dt(state, alpha=0.4)
    assert wide > narrow


def test_max_stable_dt_is_infinite_without_active_wells():
    from imex2d.simulation.implicit.well_state import (CoupledState,
                                                       WellUnknowns)

    solver, model, connections, reservoir, wells = _stepper_setup()
    # hər iki quyunun BHP-si hüceyrə təzyiqinə bərabər -> debit sıfır
    stationary = WellUnknowns(list(wells.names), reservoir.pressure[
        [connections[0].cell, connections[-1].cell]])
    state = CoupledState(reservoir, stationary)
    assert solver.max_stable_dt(state) == float("inf")


def test_adaptive_stepper_respects_the_well_rate_limit():
    """Ümumi (generic) `AdaptiveTimeStepper` `max_stable_dt`-i
    DUCK-TYPING ilə tanıyır — A6-nın iki fazalı `NewtonSolver`-inə
    HEÇ BİR TƏSİR ETMİR (o, bu metodu tanımır)."""
    from imex2d.simulation.implicit.time_stepping import (
        AdaptiveTimeStepConfig, AdaptiveTimeStepper)
    from imex2d.simulation.implicit.well_state import CoupledState

    solver, model, connections, reservoir, wells = _stepper_setup()
    state = CoupledState(reservoir, wells)
    well_limit = solver.max_stable_dt(state)

    stepper = AdaptiveTimeStepper(
        solver, AdaptiveTimeStepConfig(initial_dt=10.0, max_dt=30.0))
    _, dt, _ = stepper.advance(state, 0.0, 30.0)
    assert dt <= well_limit + 1e-9


def test_old_two_phase_engine_is_unaffected_by_the_well_limit():
    """DUCK-TYPING TƏHLÜKƏSİZLİYİ: A6-nın `NewtonSolver`-ində
    `max_stable_dt` yoxdur — `getattr(..., None)` bunu düzgün aşkarlayır
    və köhnə mühərrik heç bir dəyişiklik hiss etmir."""
    from imex2d.simulation.implicit.newton import NewtonSolver

    assert not hasattr(NewtonSolver, "max_stable_dt")
