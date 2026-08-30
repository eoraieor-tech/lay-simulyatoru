"""Adaptiv zaman addımı və fully implicit mühərrik (A6, mərhələ 4)."""

import numpy as np

from helpers import default_scal, five_spot_model, make_service
from imex2d.application.config import (OutputConfig, SimulationConfig,
                                       TimeSteppingConfig)
from imex2d.domain.wells import ControlMode, WellControl
from imex2d.interfaces.services import ISimulationEngine
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.implicit.newton import (NewtonConfig, NewtonResult,
                                               NewtonStatus)
from imex2d.simulation.implicit.state import ReservoirState
from imex2d.simulation.implicit.time_stepping import (AdaptiveTimeStepConfig,
                                                      AdaptiveTimeStepper)
from imex2d.simulation.initialization.equilibrium import EquilibriumInitializationProvider
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from test_implicit_newton import _initial, _rate_controlled, _solver


def _engine(model, scal, config=None, **kwargs):
    config = config or SimulationConfig(end_time=400.0,
                                        output=OutputConfig(snapshot_count=4))
    return FullyImplicitEngine(model, config,
                               CoreyRelativePermeabilityAdapter(scal), **kwargs)


# ── adaptiv addım ─────────────────────────────────────────────────────
def test_time_step_grows_when_newton_converges_quickly():
    scal = default_scal()
    model = _rate_controlled(nx=13, fraction=0.02, scal=scal)
    solver, _ = _solver(model, scal)
    stepper = AdaptiveTimeStepper(solver, AdaptiveTimeStepConfig(initial_dt=1.0))

    state = _initial(model)
    steps = []
    time = 0.0
    for _ in range(6):
        state, dt, _ = stepper.advance(state, time, 1e6)
        assert dt > 0
        steps.append(dt)
        time += dt
    assert steps[-1] > steps[0] * 2, steps


def test_time_step_is_cut_and_repeated_when_newton_fails():
    """Yığılmayan addım kəsilib TƏKRARLANMALIDIR, atılmamalıdır."""
    scal = default_scal()
    model = _rate_controlled(nx=13, fraction=0.2, scal=scal)
    solver, _ = _solver(model, scal, config=NewtonConfig(max_iterations=3))
    stepper = AdaptiveTimeStepper(
        solver, AdaptiveTimeStepConfig(initial_dt=300.0, min_dt=0.05))

    state, dt, result = stepper.advance(_initial(model), 0.0, 1e6)
    assert dt > 0, "Addım heç cür yığılmadı"
    assert dt < 300.0, "Addım kəsilmədi"
    assert stepper.history[-1].repeats > 0


def test_saturation_change_limit_is_respected():
    scal = default_scal()
    model = _rate_controlled(nx=13, fraction=0.05, scal=scal)
    solver, _ = _solver(model, scal)
    stepper = AdaptiveTimeStepper(
        solver, AdaptiveTimeStepConfig(initial_dt=200.0,
                                       max_saturation_change=0.05))
    state = _initial(model)
    for _ in range(4):
        new_state, dt, _ = stepper.advance(state, 0.0, 1e6)
        assert dt > 0
        change = np.abs(new_state.water_saturation
                        - state.water_saturation).max()
        assert change <= 0.05 + 1e-9, change
        state = new_state


def test_time_step_never_exceeds_configured_bounds():
    scal = default_scal()
    model = _rate_controlled(nx=13, fraction=0.02, scal=scal)
    solver, _ = _solver(model, scal)
    config = AdaptiveTimeStepConfig(initial_dt=1.0, max_dt=25.0)
    stepper = AdaptiveTimeStepper(solver, config)
    state = _initial(model)
    time = 0.0
    for _ in range(10):
        state, dt, _ = stepper.advance(state, time, 1e6)
        assert config.min_dt <= dt <= config.max_dt + 1e-9
        time += dt


def test_remaining_time_limits_the_step():
    scal = default_scal()
    model = _rate_controlled(nx=11, fraction=0.02, scal=scal)
    solver, _ = _solver(model, scal)
    stepper = AdaptiveTimeStepper(solver,
                                  AdaptiveTimeStepConfig(initial_dt=100.0))
    _, dt, _ = stepper.advance(_initial(model), 0.0, remaining=7.0)
    assert abs(dt - 7.0) < 1e-9


def test_summary_reports_statistics():
    scal = default_scal()
    model = _rate_controlled(nx=11, fraction=0.02, scal=scal)
    solver, _ = _solver(model, scal)
    stepper = AdaptiveTimeStepper(solver)
    state = _initial(model)
    time = 0.0
    for _ in range(5):
        state, dt, _ = stepper.advance(state, time, 1e6)
        time += dt
    summary = stepper.summary()
    assert summary["addım"] == 5
    assert summary["orta Δt"] > 0
    assert summary["orta iterasiya"] > 0


class _AlmostConvergedSolver:
    """Saxta Nyuton — CNV/MB HƏMİŞƏ sərt meyarın bir az üstündə, lakin
    yumşaq hədlərin altındadır. `AdaptiveTimeStepper`-in yumşaq-qəbul
    məntiqini fiziki modeldən asılı olmadan, dəqiq idarə edilən
    dəyərlərlə sınamaq üçün."""

    def __init__(self, cnv: float, mb: float):
        self._cnv = cnv
        self._mb = mb

    def solve(self, previous, dt):
        return NewtonResult(NewtonStatus.MAX_ITERATIONS, previous.copy(), 12,
                            [self._cnv], (self._mb, self._mb))


def test_soft_failure_accepts_a_step_within_relaxed_tolerances():
    """TAPILAN SƏHV: `soft_failure_cnv_tolerance` heç vaxt işləmirdi —
    kod `getattr(result, "history", None)` axtarırdı, halbuki
    `NewtonResult`-un sahəsi `cnv_history`-dir. Səssizcə `None`
    qaytarılırdı, şərt HƏMİŞƏ False olurdu, funksiya heç vaxt işə
    düşmürdü (istifadəçi tərəfindən aktivləşdirilsə belə)."""
    solver = _AlmostConvergedSolver(cnv=1e-4, mb=1e-6)
    config = AdaptiveTimeStepConfig(initial_dt=1.0, min_dt=1.0,
                                    soft_failure_cnv_tolerance=1e-2,
                                    soft_failure_mb_tolerance=1e-4)
    stepper = AdaptiveTimeStepper(solver, config)
    state = ReservoirState(np.array([200.0]), np.array([0.3]))

    _, dt, result = stepper.advance(state, 0.0, 10.0)
    assert dt > 0, "Yumşaq qəbul işləmədi (reqressiya)"
    assert stepper.history[-1].soft_failure


def test_stall_detection_stops_after_repeated_soft_failures():
    """Yumşaq qəbul EYNI nöqtədə ardıcıl-ardıcıl işə düşəndə (davamlı,
    keçici olmayan çətinlik) simulyasiya sonsuz kiçik addımlarla
    SÜRÜNMƏK əvəzinə təmiz dayanmalıdır — "sonsuz sürünmə" nəticəsiz
    CPU yeyən, istifadəçini aldadan (t irəliləmir) ən pis haldır."""
    solver = _AlmostConvergedSolver(cnv=1e-4, mb=1e-6)
    config = AdaptiveTimeStepConfig(initial_dt=1.0, min_dt=1.0,
                                    soft_failure_cnv_tolerance=1e-2,
                                    soft_failure_mb_tolerance=1e-4,
                                    max_consecutive_soft_failures=5)
    stepper = AdaptiveTimeStepper(solver, config)
    state = ReservoirState(np.array([200.0]), np.array([0.3]))

    for _ in range(5):
        state, dt, result = stepper.advance(state, 0.0, 10.0)
        assert dt > 0

    _, dt, result = stepper.advance(state, 0.0, 10.0)
    assert dt <= 0.0, "Hədddən sonra da qəbul edilməyə davam edir"
    assert not result.converged


# ── mühərrik ──────────────────────────────────────────────────────────
def test_engine_implements_the_simulation_interface():
    scal = default_scal()
    engine = _engine(_rate_controlled(nx=11, scal=scal), scal)
    assert isinstance(engine, ISimulationEngine)
    assert hasattr(engine, "pressure") and hasattr(engine, "sw")


def test_engine_produces_a_complete_result():
    scal = default_scal()
    result = _engine(_rate_controlled(nx=13, scal=scal), scal).run()
    assert result.converged
    assert result.steps > 0
    assert result.ooip > 0

    series = result.series
    length = len(series.time)
    for values in (series.oil_rate, series.water_rate, series.water_cut,
                   series.cumulative_oil, series.average_pressure,
                   series.recovery_factor):
        assert len(values) == length
    assert result.snapshots
    assert result.snapshots[-1].water_saturation.shape == \
           result.snapshots[0].pressure.shape


def test_engine_uses_far_fewer_steps_than_impes():
    scal = default_scal()
    config = SimulationConfig(end_time=900.0,
                              output=OutputConfig(snapshot_count=3))
    impes = make_service(scal).run(_rate_controlled(nx=15, scal=scal), config)
    implicit = _engine(_rate_controlled(nx=15, scal=scal), scal, config).run()

    assert implicit.converged
    assert implicit.steps < impes.steps / 10, \
        f"IMPES {impes.steps} vs implicit {implicit.steps}"


def test_engine_agrees_with_impes_on_recovery():
    scal = default_scal()
    config = SimulationConfig(end_time=600.0,
                              output=OutputConfig(snapshot_count=3))
    impes = make_service(scal).run(_rate_controlled(nx=15, scal=scal), config)
    implicit = _engine(_rate_controlled(nx=15, scal=scal), scal, config).run()
    assert abs(implicit.final_recovery_factor
               - impes.final_recovery_factor) < 1.0


def test_result_is_insensitive_to_the_maximum_time_step():
    """Adaptiv nəzarət nəticəni Δt seçimindən asılı etməməlidir."""
    scal = default_scal()
    recoveries = []
    for max_dt in (10.0, 120.0):
        config = SimulationConfig(
            end_time=600.0, output=OutputConfig(snapshot_count=3),
            time_stepping=TimeSteppingConfig(max_dt=max_dt))
        result = _engine(_rate_controlled(nx=13, scal=scal), scal, config).run()
        assert result.converged
        recoveries.append(result.final_recovery_factor)
    assert abs(recoveries[0] - recoveries[1]) < 0.3, recoveries


def test_saturation_stays_within_limits_throughout():
    scal = default_scal()
    result = _engine(_rate_controlled(nx=13, fraction=0.1, scal=scal),
                     scal).run()
    for snapshot in result.snapshots:
        assert snapshot.water_saturation.min() >= scal.swc - 1e-9
        assert snapshot.water_saturation.max() <= 1.0 - scal.sor + 1e-9


def test_cumulative_production_is_monotonic():
    scal = default_scal()
    series = _engine(_rate_controlled(nx=13, scal=scal), scal).run().series
    assert np.all(np.diff(series.cumulative_oil) >= -1e-9)
    assert np.all(np.diff(series.cumulative_water) >= -1e-9)


def test_engine_works_in_three_dimensions():
    from imex2d.application.model_builder import ReservoirModelBuilder
    from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                              five_spot)

    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=7, ny=7, dx=30.0, dy=30.0, dz=5.0, porosity=0.2,
        permx_base=150.0, nz=3, top_depth=2000.0)
    model = ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                          scal=scal)
    rate = model.pore_volume()[0] * 0.05
    for well in model.wells:
        well.control = (WellControl(ControlMode.RATE, rate) if well.is_injector
                        else WellControl(ControlMode.BHP, 150.0))
    result = _engine(model, scal).run()
    assert result.converged
    assert result.snapshots[-1].water_saturation.shape == model.grid.shape


def test_engine_works_with_equilibration_provider():
    import dataclasses

    scal = default_scal()
    model = _rate_controlled(nx=11, scal=scal)
    surface = np.add.outer(np.zeros(11), np.arange(11) * 4.0) + 2000.0
    model.geometry = dataclasses.replace(model.geometry,
                                         top_depth_map=surface)
    model.initial_conditions.use_equilibration = True
    model.initial_conditions.datum_depth = 2000.0
    model.initial_conditions.oil_water_contact = 2035.0

    engine = _engine(model, scal,
                     initialization=EquilibriumInitializationProvider())
    assert np.ptp(engine.pressure) > 0
    result = engine.run()
    assert result.converged


def test_engine_works_with_pvt_above_bubble_point():
    scal = default_scal()
    result = _engine(_rate_controlled(nx=11, scal=scal), scal,
                     pvt=BlackOilPVTProvider(
                         build_pvt_table(bubble_point_bar=150.0))).run()
    assert result.converged


def test_failed_run_is_reported_not_silently_truncated():
    """Yığılmayan model səssizcə dayanmamalı, mesajla xəbər verməlidir."""
    scal = default_scal()
    model = _rate_controlled(nx=11, fraction=0.5, scal=scal)
    config = SimulationConfig(end_time=400.0,
                              output=OutputConfig(snapshot_count=2))
    engine = FullyImplicitEngine(
        model, config, CoreyRelativePermeabilityAdapter(scal),
        newton_config=NewtonConfig(max_iterations=2),
        time_step_config=AdaptiveTimeStepConfig(initial_dt=50.0, min_dt=20.0,
                                                max_repeats=2))
    result = engine.run()
    if not result.converged:
        assert "yığılmadı" in result.message
