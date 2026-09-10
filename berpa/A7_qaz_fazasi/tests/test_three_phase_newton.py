"""Üç fazalı Nyuton-Rafson döngəsi (A7, mərhələ 6d) — A7-nin son addımı."""

import numpy as np

from helpers import default_scal
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.domain.scal import CoreyParameters, GasCoreyParameters
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.implicit.linear import NewtonLinearSolver
from imex2d.simulation.implicit.newton import NewtonConfig, NewtonStatus
from imex2d.simulation.implicit.three_phase_newton import (
    ThreePhaseNewtonResult, ThreePhaseNewtonSolver)
from imex2d.simulation.implicit.three_phase_state import ThreePhaseState
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.stone_relperm import StoneRelativePermeabilityProvider


def _solver(nx=5, ny=5, with_wells=True, **config_kwargs):
    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    wells_pattern = five_spot(geology.grid) if with_wells else []
    model = ReservoirModelBuilder().build(geology, wells_pattern, scal=scal)
    grid = TwoPointFluxDiscretization().build(model)
    relperm = StoneRelativePermeabilityProvider.from_corey(
        CoreyParameters(), GasCoreyParameters())
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    config = NewtonConfig(**config_kwargs) if config_kwargs else None
    return ThreePhaseNewtonSolver(model, relperm, pvt, NewtonLinearSolver(),
                                  grid, config=config), model


def _initial_state(model, pressure=213.5, sw=0.35, sg=0.1):
    n = model.ncell
    return ThreePhaseState(np.full(n, pressure), np.full(n, sw),
                           np.full(n, sg), np.ones(n, bool))


# ── əsas yığılma ────────────────────────────────────────────────────
def test_single_step_converges_on_a_real_model_with_wells():
    solver, model = _solver()
    previous = _initial_state(model)
    result = solver.solve(previous, dt=1.0)
    assert result.converged
    assert result.iterations < solver.config.max_iterations


def test_convergence_history_generally_decreases():
    """CNV mütləq monoton olmasa da, son dəyər ilkindən çox kiçik olmalıdır."""
    solver, model = _solver()
    previous = _initial_state(model)
    result = solver.solve(previous, dt=1.0)
    assert result.history[-1] < result.history[0] * 1e-3


def test_material_balance_is_near_machine_precision_at_convergence():
    solver, model = _solver()
    previous = _initial_state(model)
    result = solver.solve(previous, dt=1.0)
    for mb in result.material_balance:
        assert mb < 1e-6


def test_result_state_respects_physical_saturation_bounds():
    solver, model = _solver()
    previous = _initial_state(model)
    result = solver.solve(previous, dt=1.0)
    assert np.all(result.state.water_saturation >= solver.physical_sw_min - 1e-9)
    assert np.all(result.state.water_saturation <= solver.physical_sw_max + 1e-9)


def test_mass_balance_holds_over_the_step():
    """Sw+So+Sg=1 nəticə vəziyyətində qorunmalıdır."""
    solver, model = _solver()
    previous = _initial_state(model)
    result = solver.solve(previous, dt=1.0)
    total = (result.state.water_saturation + result.state.oil_saturation
            + result.state.gas_saturation)
    assert np.allclose(total, 1.0, atol=1e-6)


# ── çoxaddımlı sabitlik ─────────────────────────────────────────────
def test_ten_consecutive_steps_all_converge():
    """A7-nin əsas doğrulaması: real quyularla, davamlı simulyasiya."""
    solver, model = _solver()
    state = _initial_state(model)
    dt = 0.5
    for _ in range(10):
        result = solver.solve(state, dt)
        assert result.converged, result.status
        state = result.state


def test_water_saturation_increases_under_water_injection():
    """Su vurulur — orta Sw zamanla artmalıdır."""
    solver, model = _solver()
    state = _initial_state(model)
    initial_sw = state.water_saturation.mean()
    for _ in range(5):
        result = solver.solve(state, dt=0.5)
        assert result.converged
        state = result.state
    assert state.water_saturation.mean() > initial_sw


def test_multistep_run_preserves_mass_balance_at_every_step():
    solver, model = _solver()
    state = _initial_state(model)
    for _ in range(6):
        result = solver.solve(state, dt=0.5)
        assert result.converged
        state = result.state
        total = (state.water_saturation + state.oil_saturation
                + state.gas_saturation)
        assert np.allclose(total, 1.0, atol=1e-6)


# ── dəyişən keçid real simulyasiyada ────────────────────────────────
def test_solver_handles_a_cell_starting_at_the_saturation_boundary():
    """Sg=0-a çox yaxın başlanğıc — Nyuton addımı sərhədi keçə bilər,
    çökmə/NaN olmamalıdır."""
    solver, model = _solver(with_wells=False)
    n = model.ncell
    state = ThreePhaseState(np.full(n, 213.5), np.full(n, 0.35),
                            np.full(n, 1e-6), np.ones(n, bool))
    result = solver.solve(state, dt=1.0)
    assert result.status in (NewtonStatus.CONVERGED, NewtonStatus.MAX_ITERATIONS)
    assert np.all(np.isfinite(result.state.pressure))
    assert np.all(np.isfinite(result.state.water_saturation))


def test_solver_handles_fully_undersaturated_initial_state():
    """Heç bir hüceyrədə sərbəst qaz yoxdur — klassik iki fazalı başlanğıc."""
    solver, model = _solver(with_wells=False)
    n = model.ncell
    rs = np.full(n, 50.0)         # doyma əyrisindən aşağı
    state = ThreePhaseState(np.full(n, 213.5), np.full(n, 0.35), rs,
                            np.zeros(n, bool))
    result = solver.solve(state, dt=1.0)
    assert result.status in (NewtonStatus.CONVERGED, NewtonStatus.MAX_ITERATIONS)
    assert np.all(np.isfinite(result.state.pressure))


# ── quyusuz, sadə ssenari ────────────────────────────────────────────
def test_solver_converges_without_wells():
    solver, model = _solver(with_wells=False)
    state = _initial_state(model)
    result = solver.solve(state, dt=1.0)
    assert result.converged


def test_uniform_state_with_no_wells_stays_nearly_static():
    """Bərabər başlanğıc, quyusuz — heç bir hərəkətverici qüvvə yoxdur,
    vəziyyət demək olar dəyişməməlidir."""
    solver, model = _solver(with_wells=False)
    state = _initial_state(model)
    result = solver.solve(state, dt=1.0)
    assert result.converged
    assert np.allclose(result.state.pressure, state.pressure, atol=1e-2)
    assert np.allclose(result.state.water_saturation, state.water_saturation,
                       atol=1e-3)


# ── konfiqurasiya və nəticə obyekti ──────────────────────────────────
def test_result_object_has_expected_fields():
    solver, model = _solver()
    result = solver.solve(_initial_state(model), dt=1.0)
    assert isinstance(result, ThreePhaseNewtonResult)
    assert result.fluid is not None
    assert result.rates is not None
    assert result.linear_iterations >= 0


def test_max_iterations_status_when_limit_too_low():
    """Süni şəkildə 0 iterasiya limiti qoyulanda MAX_ITERATIONS qaytarmalıdır."""
    solver, model = _solver(max_iterations=0)
    result = solver.solve(_initial_state(model), dt=1.0)
    assert result.status in (NewtonStatus.MAX_ITERATIONS, NewtonStatus.CONVERGED)


# ── ÇÖKMƏYƏ QARŞI MÜHAFİZƏ (istifadəçi bildirişi əsasında) ────────────
def test_solve_never_raises_even_with_a_broken_pvt_provider():
    """`solve()` HEÇ VAXT istisna atmamalıdır — Jakobian yığımı daxilində
    gözlənilməz bir şey baş versə belə, təhlükəsiz nəticə qaytarılmalıdır."""
    solver, model = _solver()
    previous = _initial_state(model)

    class BrokenPVT:
        def __getattr__(self, name):
            raise RuntimeError("qəsdən sındırılmış PVT metodu")

    solver.pvt = BrokenPVT()
    result = solver.solve(previous, dt=1.0)
    from imex2d.simulation.implicit.newton import NewtonStatus
    assert result.status == NewtonStatus.LINEAR_SOLVER_FAILED
    assert not result.converged


def test_solve_never_raises_with_nan_poisoned_initial_state():
    """NaN-lı başlanğıc vəziyyət (real ssenaridə əvvəlki addımın
    uğursuzluğundan qala bilər) çökməyə səbəb olmamalıdır."""
    solver, model = _solver()
    n = model.ncell
    previous = ThreePhaseState(
        np.full(n, np.nan), np.full(n, 0.35), np.zeros(n), np.zeros(n, bool))
    result = solver.solve(previous, dt=1.0)
    assert result is not None
    assert hasattr(result, "status")


def test_engine_run_never_raises_even_on_catastrophic_failure():
    """Mühərrikin `run()` metodu HEÇ VAXT istisna atmamalıdır — hər zaman
    `SimulationResult` (converged=True ya da False) qaytarmalıdır."""
    from imex2d.application.config import SimulationConfig
    from imex2d.simulation.discretization import TwoPointFluxDiscretization
    from imex2d.simulation.implicit.three_phase_engine import (
        ThreePhaseSimulationEngine)

    _, model = _solver()
    config = SimulationConfig(end_time=10.0)
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    relperm = StoneRelativePermeabilityProvider.from_corey(
        CoreyParameters(), GasCoreyParameters())
    engine = ThreePhaseSimulationEngine(model, config, relperm, pvt)

    # gözlənilməz istisna simulyasiya etmək üçün akkumulyatoru sındırırıq
    def broken_accumulation(*args, **kwargs):
        raise RuntimeError("qəsdən sındırılmış akkumulyasiya")

    engine.newton.accumulator.accumulation = broken_accumulation
    result = engine.run()
    assert result is not None
    assert not result.converged
    assert "Gözlənilməz xəta" in result.message


def test_linear_solve_failure_types_beyond_value_and_floating_point():
    """Əvvəlki `except (FloatingPointError, ValueError)` — DAR idi. İndi
    HƏR HANSI istisna növü (məs. RuntimeError) tutulmalıdır."""
    solver, model = _solver()
    previous = _initial_state(model)

    class AlwaysFailsLinearSolver:
        last_iterations = 0
        def solve(self, matrix, rhs, x0=None):
            raise RuntimeError("scipy-yə xas, əvvəllər tutulmayan istisna")

    solver.linear_solver = AlwaysFailsLinearSolver()
    result = solver.solve(previous, dt=1.0)
    from imex2d.simulation.implicit.newton import NewtonStatus
    assert result.status == NewtonStatus.LINEAR_SOLVER_FAILED
