"""CPR ön-şərtçisi (A6, mərhələ 5)."""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from helpers import default_scal
from imex2d.simulation.implicit.cpr import CprConfig, CprPreconditioner
from imex2d.simulation.implicit.linear import (NewtonLinearSolver,
                                               NewtonLinearSolverConfig)
from imex2d.simulation.implicit.state import VARIABLES_PER_CELL
from test_implicit_newton import _initial, _rate_controlled, _solver


def _jacobian(nx=15, dt=20.0):
    """Real Nyuton addımından alınmış Jakobian və sağ tərəf."""
    scal = default_scal()
    model = _rate_controlled(nx=nx, scal=scal)
    solver, residual = _solver(model, scal)
    state = _initial(model)
    result = solver.solve(state, dt)
    vector, fluid, _ = residual.residual(result.state, state, dt)
    return solver.J.assemble(result.state, fluid, dt).copy(), -vector, model


def _bicgstab(matrix, rhs, preconditioner, max_iterations=500):
    counter = {"n": 0}
    solution, info = spla.bicgstab(
        matrix, rhs, rtol=1e-8, atol=0.0, maxiter=max_iterations,
        M=preconditioner, callback=lambda _: counter.__setitem__("n", counter["n"] + 1))
    error = np.linalg.norm(matrix @ solution - rhs) / np.linalg.norm(rhs)
    return solution, counter["n"], info, error


# ── ötürmə operatorları ───────────────────────────────────────────────
def test_transfer_operator_shapes():
    matrix, _, model = _jacobian(nx=9)
    cpr = CprPreconditioner(matrix)
    assert cpr.ncell == model.ncell
    assert cpr.restriction.shape == (model.ncell, model.ncell * VARIABLES_PER_CELL)
    assert cpr.prolongation.shape == (model.ncell * VARIABLES_PER_CELL, model.ncell)


def test_prolongation_maps_to_pressure_variables_only():
    matrix, _, model = _jacobian(nx=7)
    cpr = CprPreconditioner(matrix)
    vector = cpr.prolongation @ np.ones(model.ncell)
    assert np.allclose(vector[0::VARIABLES_PER_CELL], 1.0)
    assert np.allclose(vector[1::VARIABLES_PER_CELL], 0.0)


def test_decoupling_removes_the_saturation_derivative():
    """CPR-in mahiyyəti: çəkili cəmdə ∂/∂Sw yox olmalıdır.

    Bu pozulsa, təzyiq alt-sistemi elliptik olmaz və ön-şərtçi
    mənasını itirər.
    """
    matrix, _, model = _jacobian(nx=9)
    cpr = CprPreconditioner(matrix)

    cells = np.arange(model.ncell)
    rows = cells * VARIABLES_PER_CELL
    saturation_columns = rows + 1

    weights = cpr.restriction.toarray()
    d_water = np.asarray(matrix[rows, saturation_columns]).ravel()
    d_oil = np.asarray(matrix[rows + 1, saturation_columns]).ravel()

    combined = (weights[cells, rows] * d_water
                + weights[cells, rows + 1] * d_oil)
    scale = np.maximum(np.abs(d_water) + np.abs(d_oil), 1e-30)
    assert np.max(np.abs(combined) / scale) < 1e-10


def test_pressure_system_is_smaller_and_nearly_symmetric():
    """Dekuplinq düzgündürsə, təzyiq matrisi elliptik operatordur."""
    matrix, _, model = _jacobian(nx=13)
    cpr = CprPreconditioner(matrix)
    assert cpr.pressure_matrix.shape == (model.ncell, model.ncell)

    asymmetry = abs(cpr.pressure_matrix - cpr.pressure_matrix.T).max()
    assert asymmetry / abs(cpr.pressure_matrix).max() < 0.1


# ── hamarlayıcı ───────────────────────────────────────────────────────
def test_block_jacobi_smoother_inverts_the_diagonal_blocks():
    matrix, _, model = _jacobian(nx=9)
    cpr = CprPreconditioner(matrix, CprConfig(smoother="block_jacobi"))
    assert cpr.smoother_kind == "blok-Jakobi"

    cell = model.ncell // 2
    rows = np.array([cell * VARIABLES_PER_CELL,
                     cell * VARIABLES_PER_CELL + 1])
    block = matrix[np.ix_(rows, rows)].toarray()

    unit = np.zeros(matrix.shape[0])
    unit[rows[0]] = 1.0
    applied = cpr._apply_block_jacobi(unit)[rows]
    assert np.allclose(block @ applied, np.array([1.0, 0.0]), atol=1e-9)


def test_block_jacobi_is_always_available():
    """ILU böyük sistemlərdə qurula bilmir — blok-Jakobi həmişə işləyir."""
    matrix, _, _ = _jacobian(nx=15)
    cpr = CprPreconditioner(matrix, CprConfig(smoother="block_jacobi"))
    assert cpr._smoother is not None
    assert cpr.smoother_kind == "blok-Jakobi"


def test_singular_blocks_do_not_break_the_smoother():
    size = 6
    matrix = sp.eye(size, format="csr").tolil()
    matrix[2, 2] = 0.0          # tək blok yaradır
    matrix[2, 3] = 0.0
    matrix[3, 2] = 0.0
    matrix[3, 3] = 0.0
    cpr = CprPreconditioner(matrix.tocsr(), CprConfig(smoother="block_jacobi"))
    result = cpr._apply_block_jacobi(np.ones(size))
    assert np.all(np.isfinite(result))


# ── həll keyfiyyəti ───────────────────────────────────────────────────
def test_cpr_solves_the_newton_system():
    matrix, rhs, _ = _jacobian(nx=15)
    _, iterations, info, error = _bicgstab(matrix, rhs,
                                           CprPreconditioner(matrix))
    assert info == 0
    assert error < 1e-7
    assert iterations < 20


def test_cpr_with_block_jacobi_still_converges():
    matrix, rhs, _ = _jacobian(nx=15)
    cpr = CprPreconditioner(matrix, CprConfig(smoother="block_jacobi"))
    _, iterations, info, error = _bicgstab(matrix, rhs, cpr)
    assert info == 0
    assert error < 1e-7
    assert iterations < 50


def test_cpr_beats_no_preconditioner_by_a_wide_margin():
    matrix, rhs, _ = _jacobian(nx=15)
    _, plain_iterations, _, _ = _bicgstab(matrix, rhs, None, 3000)
    _, cpr_iterations, _, _ = _bicgstab(matrix, rhs, CprPreconditioner(matrix))
    assert cpr_iterations * 20 < plain_iterations, \
        f"ön-şərtçisiz {plain_iterations}, CPR {cpr_iterations}"


def _preconditioner_memory(matrix, model):
    """(ILU nnz, CPR nnz) — ön-şərtçilərin saxladığı qeyri-sıfır sayı."""
    ilu = spla.spilu(matrix.tocsc(), drop_tol=1e-4, fill_factor=10)
    cpr = CprPreconditioner(matrix, CprConfig(smoother="block_jacobi"))
    pressure_ilu = spla.spilu(cpr.pressure_matrix.tocsc(),
                              drop_tol=1e-4, fill_factor=8)
    return ilu.nnz, pressure_ilu.nnz + 4 * model.ncell


def test_cpr_uses_less_memory_than_strong_ilu():
    """CPR-in ƏSAS üstünlüyü yaddaşdadır, sürətdə deyil.

    Blok-Jakobi hamarlayıcı ilə CPR yalnız təzyiq sistemi üçün
    faktorlaşdırma saxlayır (N ölçü), tam sistem üçün deyil (2N).
    """
    matrix, _, model = _jacobian(nx=25)
    ilu_memory, cpr_memory = _preconditioner_memory(matrix, model)
    assert cpr_memory < ilu_memory * 0.7, f"ILU {ilu_memory}, CPR {cpr_memory}"


def test_memory_advantage_grows_with_model_size():
    """Ölçü artdıqca CPR-in yaddaş üstünlüyü genişlənir.

    Ölçülmüş nisbətlər: 225 hüceyrədə 0.63, 3025 hüceyrədə 0.44.
    Səbəb: ILU-nun doldurulması sistem ölçüsü ilə superxətti artır,
    CPR isə yarıya bölünmüş sistemdə faktorlaşdırma aparır.
    """
    ratios = []
    for size in (15, 40):
        matrix, _, model = _jacobian(nx=size)
        ilu_memory, cpr_memory = _preconditioner_memory(matrix, model)
        ratios.append(cpr_memory / ilu_memory)
    assert ratios[1] < ratios[0], ratios


# ── həlledici ilə inteqrasiya ─────────────────────────────────────────
def test_solver_selects_ilu_for_small_systems():
    matrix, rhs, _ = _jacobian(nx=11)
    solver = NewtonLinearSolver(
        NewtonLinearSolverConfig(direct_threshold=0, preconditioner="auto",
                                 cpr_threshold=10**9))
    solver.solve(matrix, rhs)
    assert solver.last_preconditioner == "ILU"


def test_solver_selects_cpr_above_the_threshold():
    matrix, rhs, _ = _jacobian(nx=11)
    solver = NewtonLinearSolver(
        NewtonLinearSolverConfig(direct_threshold=0, preconditioner="auto",
                                 cpr_threshold=1))
    solver.solve(matrix, rhs)
    assert solver.last_preconditioner.startswith("CPR")


def test_solver_result_is_the_same_with_either_preconditioner():
    """Ön-şərtçi seçimi HƏLLİ dəyişməməlidir — yalnız yola təsir edir."""
    matrix, rhs, _ = _jacobian(nx=13)
    solutions = []
    for preconditioner in ("ilu", "cpr"):
        solver = NewtonLinearSolver(
            NewtonLinearSolverConfig(direct_threshold=0,
                                     preconditioner=preconditioner))
        solutions.append(solver.solve(matrix, rhs))
    difference = np.linalg.norm(solutions[0] - solutions[1])
    assert difference / np.linalg.norm(solutions[0]) < 1e-5


def test_engine_runs_with_cpr_and_matches_ilu_result():
    from imex2d.application.config import OutputConfig, SimulationConfig
    from imex2d.simulation.implicit.engine import FullyImplicitEngine
    from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter

    scal = default_scal()
    config = SimulationConfig(end_time=400.0,
                              output=OutputConfig(snapshot_count=2))
    recoveries = []
    for preconditioner in ("ilu", "cpr"):
        model = _rate_controlled(nx=13, scal=scal)
        engine = FullyImplicitEngine(
            model, config, CoreyRelativePermeabilityAdapter(scal))
        engine.newton.linear_solver = NewtonLinearSolver(
            NewtonLinearSolverConfig(direct_threshold=0,
                                     preconditioner=preconditioner))
        result = engine.run()
        assert result.converged, preconditioner
        recoveries.append(result.final_recovery_factor)
    assert abs(recoveries[0] - recoveries[1]) < 0.05, recoveries
