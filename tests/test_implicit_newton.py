"""Nyuton döngəsi (A6, mərhələ 3)."""

import numpy as np

from helpers import default_scal, five_spot_model, make_service
from imex2d.application.config import OutputConfig, SimulationConfig
from imex2d.domain.wells import ControlMode, WellControl
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.implicit.linear import (NewtonLinearSolver,
                                               NewtonLinearSolverConfig)
from imex2d.simulation.implicit.newton import (NewtonConfig, NewtonSolver,
                                               NewtonStatus)
from imex2d.simulation.implicit.residual import ResidualAssembler
from imex2d.simulation.implicit.state import ReservoirState
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.well_model import PeacemanWellModel


def _solver(model, scal, wells=True, pvt=None, config=None):
    grid = TwoPointFluxDiscretization().build(model)
    connections = PeacemanWellModel().build_connections(model) if wells else []
    residual = ResidualAssembler(model, grid, connections,
                                 CoreyRelativePermeabilityAdapter(scal), pvt=pvt)
    return NewtonSolver(residual, config=config), residual


def _initial(model):
    ic = model.initial_conditions
    return ReservoirState(np.full(model.ncell, ic.datum_pressure),
                          np.full(model.ncell, ic.water_saturation))


def _rate_controlled(nx=15, fraction=0.05, scal=None):
    """Debit quyu hüceyrəsinin məsamə həcminə nisbətdə verilir."""
    scal = scal or default_scal()
    model = five_spot_model(nx=nx, ny=nx, scal=scal)
    rate = model.pore_volume()[0] * fraction
    for well in model.wells:
        well.control = (WellControl(ControlMode.RATE, rate) if well.is_injector
                        else WellControl(ControlMode.BHP, 150.0))
    return model


# ── xətti həlledici ───────────────────────────────────────────────────
def test_linear_solver_handles_unsymmetric_matrix():
    """Jakobian simmetrik deyil — KQ işləməz, BiCGStab/LU lazımdır."""
    import scipy.sparse as sp

    matrix = sp.csr_matrix(np.array([[4.0, 1.0, 0.0],
                                     [0.0, 3.0, 1.0],
                                     [2.0, 0.0, 5.0]]))
    rhs = np.array([1.0, 2.0, 3.0])
    solution = NewtonLinearSolver().solve(matrix, rhs)
    assert np.allclose(matrix @ solution, rhs, atol=1e-8)


def test_linear_solver_uses_iterative_path_for_large_systems():
    import scipy.sparse as sp

    size = 400
    matrix = sp.diags([np.full(size - 1, -1.0), np.full(size, 4.0),
                       np.full(size - 1, -1.2)], [-1, 0, 1], format="csr")
    rhs = np.ones(size)
    solver = NewtonLinearSolver(
        NewtonLinearSolverConfig(direct_threshold=10))
    solution = solver.solve(matrix, rhs)
    assert solver.last_method.startswith("BiCGStab")
    assert np.allclose(matrix @ solution, rhs, atol=1e-6)


# ── konvergensiya ─────────────────────────────────────────────────────
def test_converges_in_one_iteration_at_equilibrium():
    """Quyusuz, tarazlıqdakı model artıq həlldir."""
    scal = default_scal()
    model = five_spot_model(nx=7, ny=7, scal=scal)
    solver, _ = _solver(model, scal, wells=False)
    result = solver.solve(_initial(model), dt=10.0)
    assert result.converged
    assert result.iterations == 0


def test_converges_for_a_wide_range_of_time_steps_without_wells():
    """Quyusuz axın məsələsi istənilən Δt-də bir addımda həll olunur."""
    scal = default_scal()
    model = five_spot_model(nx=7, ny=7, scal=scal)
    solver, _ = _solver(model, scal, wells=False)
    state = ReservoirState(
        np.linspace(230.0, 270.0, model.ncell),
        np.full(model.ncell, 0.45))
    for dt in (1.0, 100.0, 1000.0):
        result = solver.solve(state, dt)
        assert result.converged, f"dt={dt}: {result.status}"
        assert result.iterations <= 2


def test_converges_with_wells_at_moderate_time_steps():
    scal = default_scal()
    model = _rate_controlled(scal=scal)
    solver, _ = _solver(model, scal)
    state = _initial(model)
    for dt in (1.0, 10.0, 30.0):
        result = solver.solve(state, dt)
        assert result.converged, f"dt={dt}: {result.status.value}"


def test_iteration_count_grows_with_time_step():
    scal = default_scal()
    model = _rate_controlled(scal=scal)
    solver, _ = _solver(model, scal)
    state = _initial(model)
    small = solver.solve(state, 1.0).iterations
    large = solver.solve(state, 30.0).iterations
    assert large > small


def test_residual_decreases_monotonically_in_easy_cases():
    scal = default_scal()
    model = _rate_controlled(fraction=0.02, scal=scal)
    solver, _ = _solver(model, scal)
    result = solver.solve(_initial(model), dt=10.0)
    assert result.converged
    history = result.cnv_history
    assert all(later < earlier * 1.05
               for earlier, later in zip(history, history[1:]))


def test_reports_max_iterations_when_step_is_too_large():
    """Çox böyük Δt-də Nyuton yığılmır — bu, adaptiv Δt üçün siqnaldır."""
    scal = default_scal()
    model = _rate_controlled(fraction=0.2, scal=scal)
    solver, _ = _solver(model, scal, config=NewtonConfig(max_iterations=5))
    result = solver.solve(_initial(model), dt=365.0)
    assert not result.converged
    assert result.status is NewtonStatus.MAX_ITERATIONS


def test_convergence_measures_are_dimensionless():
    """CNV və MB grid ölçüsündən asılı olmamalıdır."""
    scal = default_scal()
    measures = []
    for size in (7, 13):
        model = five_spot_model(nx=size, ny=size, scal=scal)
        solver, residual = _solver(model, scal, wells=False)
        state = _initial(model)
        perturbed = ReservoirState(state.pressure.copy(),
                                   state.water_saturation + 0.02)
        vector, fluid, _ = residual.residual(perturbed, state, 5.0)
        measures.append(solver.convergence_measures(vector, perturbed,
                                                    fluid, 5.0))
    assert abs(measures[0][0] - measures[1][0]) / measures[0][0] < 1e-9


def test_material_balance_is_tighter_than_local_residual():
    """MB toleransı CNV-dən qat-qat sıx olmalıdır."""
    config = NewtonConfig()
    assert config.material_balance_tolerance < config.cnv_tolerance / 100.0


# ── fiziki doğruluq ───────────────────────────────────────────────────
def test_solution_is_independent_of_the_time_step():
    """Eyni məsələ müxtəlif Δt-lərdə eyni cavaba yaxınlaşmalıdır."""
    scal = default_scal()
    end_time = 300.0
    recoveries = []
    for dt in (10.0, 60.0):
        model = _rate_controlled(scal=scal)
        solver, residual = _solver(model, scal)
        state = _initial(model)
        initial_oil = residual.accumulation(
            state, residual.fluid_state(state))[1].sum()
        time = 0.0
        while time < end_time - 1e-9:
            step = min(dt, end_time - time)
            result = solver.solve(state, step)
            assert result.converged, f"dt={dt} t={time}"
            state = result.state
            time += step
        final_oil = residual.accumulation(
            state, residual.fluid_state(state))[1].sum()
        recoveries.append((initial_oil - final_oil) / initial_oil * 100.0)
    assert abs(recoveries[0] - recoveries[1]) < 0.1, recoveries


def test_agrees_with_impes_within_discretisation_error():
    """İki sxem eyni fizikanı həll edir — fərq yalnız diskretizasiyadadır."""
    scal = default_scal()
    end_time = 400.0

    impes = make_service(scal).run(
        _rate_controlled(scal=scal),
        SimulationConfig(end_time=end_time,
                         output=OutputConfig(snapshot_count=2)))

    model = _rate_controlled(scal=scal)
    solver, residual = _solver(model, scal)
    state = _initial(model)
    initial_oil = residual.accumulation(state,
                                        residual.fluid_state(state))[1].sum()
    time = 0.0
    while time < end_time - 1e-9:
        step = min(30.0, end_time - time)
        result = solver.solve(state, step)
        assert result.converged
        state = result.state
        time += step
    final_oil = residual.accumulation(state,
                                      residual.fluid_state(state))[1].sum()
    newton_rf = (initial_oil - final_oil) / initial_oil * 100.0

    assert abs(newton_rf - impes.final_recovery_factor) < 1.0, \
        f"IMPES {impes.final_recovery_factor:.3f} vs Newton {newton_rf:.3f}"


def test_saturation_stays_within_physical_limits():
    scal = default_scal()
    model = _rate_controlled(fraction=0.1, scal=scal)
    solver, _ = _solver(model, scal)
    state = _initial(model)
    for _ in range(10):
        result = solver.solve(state, 20.0)
        assert result.converged
        state = result.state
        assert state.water_saturation.min() >= scal.swc - 1e-12
        assert state.water_saturation.max() <= 1.0 - scal.sor + 1e-12


def test_works_with_pvt_above_the_bubble_point():
    """PVT ilə Nyuton model doyma təzyiqinin ÜSTÜNDƏ qalanda yığılır."""
    scal = default_scal()
    model = _rate_controlled(scal=scal)
    solver, _ = _solver(model, scal,
                        pvt=BlackOilPVTProvider(
                            build_pvt_table(bubble_point_bar=150.0)))
    state = _initial(model)
    for dt in (5.0, 20.0, 60.0):
        result = solver.solve(state, dt)
        assert result.converged, f"dt={dt}: {result.status.value}"


def test_crossing_the_bubble_point_now_converges():
    """ƏVVƏLKİ MƏLUM MƏHDUDİYYƏT DÜZƏLDİLDİ: Pb-də Bo/μo-nun DOYMUŞ
    qoldan DOYMAMIŞ qola keçməsi ∂Bo/∂p-nin işarəsini dəyişdirirdi və
    Nyutonu ossilyasiyaya salırdı. Qaz fazası (mühərrikdə HEÇ YERDƏ
    istifadə olunmayan Rs xaric) onsuz da modelləşdirilmədiyi üçün
    (bax `ReservoirModel.diagnose()`-un xəbərdarlığı: "nəticələr nikbin
    ola bilər"), Bo/μo indi Pb-nin HƏR İKİ tərəfində EYNİ (doymamış
    maye) düsturu ilə HAMAR davam edir — bax `test_pvt.py`-də
    `test_oil_fvf_is_smooth_across_the_bubble_point`. Tam variable
    switching (real qaz fazası) hələ də A7-nin işi olaraq qalır, lakin
    bu artıq YIĞILMA üçün lazım deyil.
    """
    scal = default_scal()
    model = _rate_controlled(scal=scal)
    solver, _ = _solver(model, scal,
                        pvt=BlackOilPVTProvider(
                            build_pvt_table(bubble_point_bar=240.0)))
    result = solver.solve(_initial(model), dt=20.0)
    assert result.converged, result.status.value


def test_line_search_prevents_infinite_oscillation_near_a_well():
    """TAPILAN SƏHV: geri-izləmə olmadan Nyuton quyu hüceyrəsi ətrafında
    bir neçə vəziyyət arasında SONSUZ DÖVRƏYƏ düşürdü (upstream
    mobillik seçimi iterasiyadan-iterasiyaya dəyişəndə) — CNV heç
    vaxt ~1.4×10⁻³-dən aşağı düşmədən dəqiq təkrarlanırdı (ölçüldü:
    [...0.00151, 0.00861, 0.00292, 0.00151, 0.00254...] sonsuza qədər).
    Kəsmə (Appleyard) YALNIZ addımın uzunluğunu məhdudlaşdırır,
    istiqamətin faydalı olub-olmadığını yoxlamır;
    `CoupledNewtonSolver`-də artıq doğrulanmış geri-izləmə (bax onun
    sənədləşməsi) bunu `NewtonSolver`-ə də əlavə edir — nəticədə CNV
    MONOTON azalır (dövr etmir) və tolerantlığa çox yaxınlaşır.

    Kiçik Δt-nin BÖYÜKDƏN daha çətin olması gözlənilənin TƏRSİNƏdir
    (implicit sxem şərti-sabit olmalıdır) — məhz bu, problemi üzə
    çıxarır: Δt=2.0 həmin nöqtəni bir sıçrayışda keçir, Δt=0.25 isə
    məhz oscillasiya zolağında ilişib qalır.
    """
    scal = default_scal()
    model = five_spot_model(nx=15, scal=scal)
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=200.0))
    solver, _ = _solver(model, scal, pvt=pvt,
                        config=NewtonConfig(max_iterations=20))
    state = _initial(model)
    dt = 0.25

    for _ in range(40):
        result = solver.solve(state, dt)
        if not result.converged:
            break
        state = result.state
    else:
        return  # bütün addımlar yığılıb — oscillasiya artıq baş vermir

    # Dövr edən (unfixed) hal 0.00137-dən aşağı heç vaxt düşmürdü;
    # geri-izləmə ilə monoton azalma tolerantlığın (1e-3) lap yaxınına
    # çatır (~1.0035e-3). Hər ikisini ayıran hədd.
    assert min(result.cnv_history) < 0.0012, (
        "CNV gözlənilən qədər aşağı düşmür — geri-izləmə reqressiyası?")


def test_compressibility_allows_production_beyond_injected_volume():
    """Sıxılma olmadan hasilat yalnız vurulan həcmlə məhdudlaşardı.

    Bu, IMPES ilə uyğunsuzluğun səbəbi idi: IMPES `ct` ilə sıxılmanı
    modelləşdirirdi, implicit sxem isə B-ni sabit saxlayırdı.
    """
    scal = default_scal()
    model = five_spot_model(nx=9, ny=9, scal=scal)
    solver, residual = _solver(model, scal)
    state = _initial(model)

    initial_oil = residual.accumulation(state,
                                        residual.fluid_state(state))[1].sum()
    result = solver.solve(state, dt=5.0)
    assert result.converged

    produced = initial_oil - residual.accumulation(
        result.state, residual.fluid_state(result.state))[1].sum()
    injected = result.rates.water.sum() * 5.0
    assert produced > 0
    # sıxılma hesabına çıxarılan neft vurulan sudan çox ola bilər
    assert produced > injected * 0.5
