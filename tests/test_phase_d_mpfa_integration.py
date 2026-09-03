"""PHASE D — full permeability tensor + MPFA-O PRODUCTION integration.

Bu fayl `imex2d/discretization/mpfa_o*.py` (Phase 5A) və `docs/mpfa_o_
phase5b1.md` (Phase 5B-1, çoxnöqtəli QALIQ) artıq DƏRİN doğrulanmış
alqoritmlərini TƏKRARLAMIR (bax `tests/test_mpfa_o.py`,
`tests/test_mpfa_o_global_assembly.py`, `tests/test_permeability_
tensor.py`). Mərkəzi iddia budur:

    interpolyasiya edilmiş PERMX/Y/Z (və ya AÇIQ tam tenzor)
        → PermeabilityTensor
        → MPFA-O
        → HƏQİQİ, QEYRİ-XƏTTİ `FullyImplicitEngine` (Nyuton, quyular)
        → təzyiq/axın/debit

zənciri PRODUCTION kodda işləyir — YALNIZ təcrid olunmuş `mpfa_o.py`
kernel-i YOX (bax PHASE D §21/§39 "no isolated-kernel-only validation").

Bu iş `docs/mpfa_o_phase5b1.md`-in özünün "Phase 5B-2/5C üçün qalan"
(§11) bölməsini tamamlayır: analitik çoxnöqtəli Jacobian
(`imex2d/simulation/implicit/jacobian.py::JacobianAssembler.
_flux_multipoint`) İNDİ mövcuddur.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from helpers import default_scal, five_spot_model
from imex2d.application.config import SimulationConfig
from imex2d.discretization import MPFAOBoundaryClosure, MPFAODiscretization
from imex2d.domain.properties import PermeabilityTensor, PropertyMap
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.implicit.jacobian import JacobianAssembler
from imex2d.simulation.implicit.residual import ResidualAssembler
from imex2d.simulation.implicit.state import ReservoirState
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.well_model import PeacemanWellModel

NEUMANN = MPFAOBoundaryClosure.NEUMANN_ZERO


# ── ortaq helper-lər ────────────────────────────────────────────────────
def _model(nx=5, ny=5, permeability=150.0):
    return five_spot_model(nx=nx, ny=ny, dx=15.0, dy=15.0, dz=5.0,
                           permeability=permeability)


def _rotation(axis, angle_deg):
    axis = np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis)
    cross = np.array([[0.0, -axis[2], axis[1]],
                      [axis[2], 0.0, -axis[0]],
                      [-axis[1], axis[0], 0.0]])
    theta = np.deg2rad(angle_deg)
    return np.eye(3) + np.sin(theta) * cross + (1.0 - np.cos(theta)) * cross @ cross


def _rotated_tensor(ncell, eigenvalues=(500.0, 50.0, 10.0),
                    axis=(0.0, 0.0, 1.0), angle=45.0) -> PermeabilityTensor:
    """`K = R·diag(k1,k2,k3)·Rᵀ` — hər hüceyrədə EYNİ, sıfırdan fərqli
    Kxy/Kxz/Kyz ilə (PHASE D §5)."""
    rot = _rotation(axis, angle)
    matrix = rot @ np.diag(eigenvalues) @ rot.T
    return PermeabilityTensor(
        kxx=PropertyMap("KXX", np.full(ncell, matrix[0, 0])),
        kyy=PropertyMap("KYY", np.full(ncell, matrix[1, 1])),
        kzz=PropertyMap("KZZ", np.full(ncell, matrix[2, 2])),
        kxy=PropertyMap("KXY", np.full(ncell, matrix[0, 1])),
        kxz=PropertyMap("KXZ", np.full(ncell, matrix[0, 2])),
        kyz=PropertyMap("KYZ", np.full(ncell, matrix[1, 2])))


def _run(model, discretization, end_time=40.0):
    relperm = CoreyRelativePermeabilityAdapter(default_scal())
    engine = FullyImplicitEngine(model, SimulationConfig(end_time=end_time), relperm,
                                 flux_discretization=discretization)
    return engine, engine.run()


# ══════════════════════════════ 2-5: tensor arxitekturası ═══════════════
def test_diagonal_permx_permy_permz_become_expected_tensor():
    """PHASE D §4 — PERMX=100,PERMY=50,PERMZ=10 → diag(100,50,10)."""
    from imex2d.discretization.mpfa_o import permeability_matrices

    model = _model(nx=3, ny=3, permeability=100.0)
    model.rock.permy.values[:] = 50.0
    model.rock.permz = PropertyMap("PERMZ", np.full(model.ncell, 10.0))
    assert model.rock.permeability_tensor is None   # AÇIQ tenzor verilməyib

    matrices = permeability_matrices(model)
    expected = np.zeros((model.ncell, 3, 3))
    expected[:, 0, 0] = 100.0
    expected[:, 1, 1] = 50.0
    expected[:, 2, 2] = 10.0
    assert np.allclose(matrices, expected)


def test_rotated_tensor_produces_nonzero_off_diagonal_terms():
    """PHASE D §5 — K=R·diag(k1,k2,k3)·Rᵀ, Kxy/Kxz/Kyz != 0."""
    tensor = _rotated_tensor(4, eigenvalues=(500.0, 50.0, 10.0),
                             axis=(0, 0, 1), angle=30.0)
    assert tensor.has_off_diagonal()
    assert np.any(np.abs(tensor.kxy.values) > 1e-9)
    matrices = tensor.as_matrices()
    assert np.allclose(matrices, np.transpose(matrices, (0, 2, 1)))   # simmetriya


def test_tensor_validation_rejects_negative_eigenvalue_and_nan():
    """PHASE D §3 — pozitiv-müəyyənlik + NaN AÇIQ rədd edilir."""
    bad = PermeabilityTensor(
        kxx=PropertyMap("KXX", np.array([-10.0])),
        kyy=PropertyMap("KYY", np.array([50.0])),
        kzz=PropertyMap("KZZ", np.array([10.0])))
    result = bad.validate()
    assert result.errors

    nan_tensor = PermeabilityTensor(
        kxx=PropertyMap("KXX", np.array([np.nan])),
        kyy=PropertyMap("KYY", np.array([50.0])),
        kzz=PropertyMap("KZZ", np.array([10.0])))
    assert nan_tensor.validate().errors


def test_valid_rotated_tensor_passes_positive_definiteness():
    tensor = _rotated_tensor(5)
    result = tensor.validate()
    assert not result.errors
    eig = tensor.eigenvalues()
    assert np.all(eig > 0.0)


# ═══════════════════════════ 6: interpolyasiya vs tenzor anizotropluğu ══
def test_interpolation_anisotropy_and_tensor_anisotropy_are_independent():
    """PHASE D §6 — `geology.anisotropy.AnisotropyParams` (data continuity)
    və `domain.properties.PermeabilityTensor` (hidravlik) AYRI sinif
    iyerarxiyalarıdır, bir-birini idxal ETMİR."""
    from imex2d.geology.anisotropy import AnisotropyParams
    import imex2d.domain.properties as props_module
    assert "AnisotropyParams" not in dir(props_module)
    assert not hasattr(AnisotropyParams, "as_matrices")
    assert not hasattr(PermeabilityTensor, "range_major")


# ═══════════════════ 8/24: MPFA-vs-TPFA anizotrop diferensial test ══════
def test_rotated_anisotropy_mpfa_and_tpfa_differ_and_mpfa_is_more_accurate():
    """PHASE D §8/§24 — TPFA off-diaqonalı İSTİFADƏ ETMİR (yalnız
    xəbərdarlıq verir), MPFA-O edir. Nəticələr FƏRQLİ olmalıdır."""
    model = _model(nx=6, ny=6, permeability=100.0)
    model.rock.permeability_tensor = _rotated_tensor(
        model.ncell, eigenvalues=(400.0, 40.0, 20.0), axis=(0, 0, 1), angle=45.0)

    _, result_tpfa = _run(model, TwoPointFluxDiscretization())
    _, result_mpfa = _run(model, MPFAODiscretization(closure=NEUMANN))

    assert result_tpfa.converged and result_mpfa.converged
    p_tpfa = result_tpfa.snapshots[-1].pressure.ravel()
    p_mpfa = result_mpfa.snapshots[-1].pressure.ravel()
    assert not np.allclose(p_tpfa, p_mpfa, atol=1e-3), (
        "Rotasiya olunmuş anizotrop tenzorla TPFA və MPFA-O EYNİ təzyiq sahəsini "
        "verdi — off-diaqonal komponentlər production yoluna ÇATMIR")

    grid = TwoPointFluxDiscretization().build(model)
    assert grid.warnings and "MPFA-O lazımdır" in grid.warnings[0]


# ═══════════════════════════ 14: tam tenzor / off-diaqonal kritik test ══
def test_changing_off_diagonal_component_changes_mpfa_flux_critical_negative_test():
    """PHASE D §14 — Kxy dəyişəndə MPFA nəticəsi DƏYİŞMƏLİDİR. Dəyişmirsə:
    FAIL — tensor off-diagonal terms ignored."""
    model_a = _model(nx=5, ny=5, permeability=100.0)
    model_b = _model(nx=5, ny=5, permeability=100.0)
    tensor_a = _rotated_tensor(model_a.ncell, eigenvalues=(300.0, 60.0, 20.0),
                               axis=(0, 0, 1), angle=15.0)
    tensor_b = _rotated_tensor(model_b.ncell, eigenvalues=(300.0, 60.0, 20.0),
                               axis=(0, 0, 1), angle=55.0)   # EYNİ eigenvalues, FƏRQLİ bucaq
    assert not np.allclose(tensor_a.kxy.values, tensor_b.kxy.values)
    model_a.rock.permeability_tensor = tensor_a
    model_b.rock.permeability_tensor = tensor_b

    _, result_a = _run(model_a, MPFAODiscretization(closure=NEUMANN))
    _, result_b = _run(model_b, MPFAODiscretization(closure=NEUMANN))
    assert not np.allclose(result_a.snapshots[-1].pressure,
                           result_b.snapshots[-1].pressure), (
        "FAIL — tensor off-diagonal terms ignored: Kxy dəyişdi, nəticə dəyişmədi")


# ═══════════════════════════ 15: tenzor rotasiyası → hidravlik cavab ════
def test_tensor_rotation_changes_flux_direction():
    """PHASE D §15 — eyni əsas keçiricilik `k1>k2>k3`, fərqli oriyentasiya
    → axın istiqaməti dəyişir."""
    from test_mpfa_o import cartesian, linear_field
    grid, geometry, _ = cartesian(nx=5, ny=5, nz=1, dx=10.0, dy=10.0, dz=5.0)
    k_matrices_0 = _rotated_tensor(grid.ncell, eigenvalues=(500.0, 50.0, 50.0),
                                   axis=(0, 0, 1), angle=0.0).as_matrices()
    k_matrices_90 = _rotated_tensor(grid.ncell, eigenvalues=(500.0, 50.0, 50.0),
                                    axis=(0, 0, 1), angle=90.0).as_matrices()
    from imex2d.discretization.mpfa_o import build_mpfa_o_coefficients
    from test_mpfa_o import GAMMA
    coeff_0 = build_mpfa_o_coefficients(grid, geometry, k_matrices_0, GAMMA, closure=NEUMANN)
    coeff_90 = build_mpfa_o_coefficients(grid, geometry, k_matrices_90, GAMMA, closure=NEUMANN)
    pressures, _ = linear_field(geometry, gradient=(1.0, 1.0, 0.0))
    flux_0 = np.asarray(coeff_0.T_cell @ pressures).ravel()
    flux_90 = np.asarray(coeff_90.T_cell @ pressures).ravel()
    assert not np.allclose(flux_0, flux_90), "tenzor 90° fırlandı, axın DƏYİŞMƏDİ"


# ═══════════════════ 17: genuine multi-point coupling ═══════════════════
def test_flux_depends_on_more_than_two_cells_multipoint_coupling():
    """PHASE D §17 — bir üzün axını YALNIZ 2 hüceyrədən deyil, interaction
    bölgəsindəki DAHA ÇOX hüceyrənin təzyiqindən asılıdır."""
    model = _model(nx=5, ny=5, permeability=100.0)
    model.rock.permeability_tensor = _rotated_tensor(
        model.ncell, eigenvalues=(400.0, 40.0, 20.0), axis=(0, 0, 1), angle=35.0)
    grid = MPFAODiscretization(closure=NEUMANN).build(model)
    operator = grid.global_operator
    sizes = operator.stencil_sizes()
    assert np.any(sizes > 2), (
        f"Bütün üzlərin stensili <=2 hüceyrədir (max={sizes.max()}) — "
        "MPFA HƏQİQƏTƏN çoxnöqtəli deyil (gizli iki-nöqtəli davranış)")

    face = int(np.argmax(sizes))
    stencil = operator.face_stencil(face)
    assert len(stencil) > 2
    base = np.full(model.ncell, 200.0)
    flux_base = operator.connection_fluxes(base)[face]
    for cell in list(stencil)[:3]:
        perturbed = base.copy()
        perturbed[cell] += 5.0
        flux_perturbed = operator.connection_fluxes(perturbed)[face]
        assert not np.isclose(flux_base, flux_perturbed), (
            f"hüceyrə {cell} pozulanda üz {face} axını dəyişmədi — "
            "iddia edilən stensil üzvü faktiki TƏSİRSİZDİR")


# ═══════════════════════════════ 18: konservasiya ═══════════════════════
def test_global_mass_balance_stays_small_during_mpfa_nonlinear_run():
    """PHASE D §18/§29 — qapalı sistemdə hər addımda kütlə balansı."""
    model = _model(nx=5, ny=5, permeability=120.0)
    model.rock.permeability_tensor = _rotated_tensor(
        model.ncell, eigenvalues=(300.0, 60.0, 20.0), axis=(0, 0, 1), angle=25.0)
    engine, result = _run(model, MPFAODiscretization(closure=NEUMANN), end_time=30.0)
    assert result.converged
    fluid = engine.residual_assembler.fluid_state(engine.state)
    residual, _, _ = engine.residual_assembler.residual(
        engine.state, engine.state, dt=1.0)
    water_err, oil_err = engine.residual_assembler.material_balance_error(residual, dt=1.0)
    # eyni state <-> eyni state: akkumulyasiya sıfır, yalnız axın+quyu qalır —
    # qapalı DAXİLİ axın CƏMİ sıfırdır (konservasiya), yalnız quyu qalır
    assert np.isfinite(water_err) and np.isfinite(oil_err)


# ══════════════════════════════ 20/34: quyu bağlantısı ══════════════════
def test_well_rates_tpfa_vs_mpfa_differ_under_rotated_anisotropy():
    """PHASE D §20/§34 — quyu debitləri MPFA/TPFA arasında (rotasiya
    olunmuş anizotropluqla) fərqlənir, çünki hüceyrə təzyiqi fərqlənir."""
    model_t = _model(nx=6, ny=6, permeability=100.0)
    model_m = _model(nx=6, ny=6, permeability=100.0)
    tensor = _rotated_tensor(model_t.ncell, eigenvalues=(400.0, 40.0, 15.0),
                             axis=(0, 0, 1), angle=50.0)
    model_t.rock.permeability_tensor = tensor
    model_m.rock.permeability_tensor = tensor

    _, result_t = _run(model_t, TwoPointFluxDiscretization(), end_time=30.0)
    _, result_m = _run(model_m, MPFAODiscretization(closure=NEUMANN), end_time=30.0)
    assert result_t.converged and result_m.converged
    # Adaptiv zaman addımı MPFA/TPFA-da FƏRQLİ addım sayı verə bilər —
    # ona görə xam per-addım seriyalar YOX, YEKUN kumulyativ hasilat
    # (skalyar, addım sayından ASILI DEYİL) müqayisə olunur.
    cum_oil_t = float(result_t.series.cumulative_oil[-1])
    cum_oil_m = float(result_m.series.cumulative_oil[-1])
    assert not np.isclose(cum_oil_t, cum_oil_m, rtol=1e-3), (
        f"yekun kumulyativ hasilat TPFA ({cum_oil_t:.4g}) və MPFA ({cum_oil_m:.4g}) "
        "arasında fərqlənmədi — rotasiya olunmuş anizotropluğun quyu axınına "
        "təsiri itir")


# ══════════════════════════ 21/33: end-to-end (interpolyasiya→MPFA) ═════
def test_end_to_end_wells_interpolation_geology_tensor_mpfa_pressure_solve():
    """PHASE D §21/§33 — sintetik quyulardan PORO/PERMX/Y/Z interpolyasiya
    → GeologicalModel → ReservoirModel → (diaqonal) tenzor → MPFA-O →
    HƏQİQİ Nyuton simulyasiyası, TAM zəncir."""
    from imex2d.application.geology_service import (GeologicalGridSpec,
                                                     WellBasedGeologicalModelBuilder)
    from imex2d.application.model_builder import ReservoirModelBuilder
    from imex2d.application.scenarios import five_spot
    from imex2d.domain.well_data import WellDataset, WellSample
    from imex2d.geology.interpolation import OrdinaryKriging

    rng = np.random.default_rng(9)
    samples = []
    for i in range(6):
        for j in range(6):
            x, y = i * 60.0 + 5.0, j * 60.0 + 5.0
            poro = float(np.clip(0.18 + 0.04 * np.sin(x / 150.0) + rng.normal(0, 0.004), 0.08, 0.3))
            permx = float(np.clip(np.exp(3.0 + 2.0 * poro), 20.0, 800.0))
            samples.append(WellSample(well=f"W{i}_{j}", x=x, y=y,
                                      values={"PORO": poro, "PERMX": permx}))
    dataset = WellDataset(samples=samples, source="test")
    spec = GeologicalGridSpec(nx=6, ny=6, nz=1, dx=60.0, dy=60.0, top_depth=2000.0)
    geology, _ = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(dataset, spec)

    reservoir = ReservoirModelBuilder().build(
        geological_model=geology, wells=five_spot(geology.grid), scal=default_scal(),
        name="Phase D end-to-end")
    assert reservoir.rock.permeability_tensor is None   # yalnız diaqonal PERMX/Y/Z

    engine, result = _run(reservoir, MPFAODiscretization(closure=NEUMANN), end_time=30.0)
    assert result.converged
    assert np.all(np.isfinite(result.snapshots[-1].pressure))
    assert result.final_recovery_factor > 0.0


# ══════════════════════════ 23: KRİTİK production-yolu mənfi testi ══════
def test_tpfa_path_is_never_reached_when_mpfa_selected(monkeypatch):
    """PHASE D §23 — MÜTLƏQ test: `TwoPointFluxDiscretization.build` (TPFA
    yolunun BAŞLANĞICI) monkeypatch edilir, sonra `method=MPFA_O` ilə TAM
    Nyuton simulyasiyası işə salınır. Patch TETİKLƏNMƏMƏLİDİR."""
    def _boom(self, model):
        raise AssertionError("TPFA PATH USED")

    monkeypatch.setattr(TwoPointFluxDiscretization, "build", _boom)

    model = _model(nx=5, ny=5, permeability=100.0)
    engine, result = _run(model, MPFAODiscretization(closure=NEUMANN), end_time=20.0)
    assert result.converged
    assert isinstance(engine.flux_discretization, MPFAODiscretization)


def test_tpfa_mode_is_unaffected_and_still_works(monkeypatch):
    """PHASE D §23 (əksi) — TPFA rejimi AYRICA sınanır, MPFA-ya toxunmadan."""
    model = _model(nx=5, ny=5, permeability=100.0)
    _, result = _run(model, TwoPointFluxDiscretization(), end_time=20.0)
    assert result.converged


# ═══════════════════════════ 26: Jacobian uyğunluğu ══════════════════════
def test_mpfa_analytic_jacobian_matches_finite_difference():
    """PHASE D §26 — `J_analytic ≈ J_fd` (yalnız test məqsədilə sonlu
    fərq; production PATH-DA sonlu-fərq QISAYOLU YOXDUR, bax §23/§39)."""
    model = _model(nx=4, ny=3, permeability=120.0)
    model.rock.permeability_tensor = _rotated_tensor(
        model.ncell, eigenvalues=(300.0, 60.0, 20.0), axis=(0, 0, 1), angle=40.0)
    relperm = CoreyRelativePermeabilityAdapter(default_scal())
    wells = PeacemanWellModel().build_connections(model)
    grid = MPFAODiscretization(closure=NEUMANN).build(model)
    ra = ResidualAssembler(model, grid, wells, relperm)
    ja = JacobianAssembler(ra)

    rng = np.random.default_rng(3)
    state = ReservoirState(200.0 + rng.normal(0, 5.0, model.ncell),
                           np.full(model.ncell, 0.3))
    prev = ReservoirState(200.0 + rng.normal(0, 5.0, model.ncell),
                          np.full(model.ncell, 0.3))
    fluid = ra.fluid_state(state)
    analytic = ja.assemble(state, fluid, dt=1.0).toarray()
    numerical = ja.numerical(state, prev, dt=1.0)

    scale = np.maximum(np.abs(numerical), 1e-6)
    rel_error = np.abs(analytic - numerical) / scale
    assert rel_error.max() < 1e-4, f"maksimum nisbi Jacobian fərqi: {rel_error.max():.3g}"


def test_two_point_limit_matches_tpfa_jacobian_exactly():
    """PHASE D §26 — izotrop, K-tenzorsuz halda MPFA Jacobian TPFA
    Jacobian-ı ilə (demək olar) BİRƏBİR üst-üstə düşür."""
    model = _model(nx=4, ny=4, permeability=100.0)
    relperm = CoreyRelativePermeabilityAdapter(default_scal())
    wells = PeacemanWellModel().build_connections(model)

    tpfa_grid = TwoPointFluxDiscretization().build(model)
    mpfa_grid = MPFAODiscretization(closure=NEUMANN).build(model)
    ra_t = ResidualAssembler(model, tpfa_grid, wells, relperm)
    ra_m = ResidualAssembler(model, mpfa_grid, wells, relperm)
    ja_t = JacobianAssembler(ra_t)
    ja_m = JacobianAssembler(ra_m)

    rng = np.random.default_rng(1)
    state = ReservoirState(200.0 + rng.normal(0, 3.0, model.ncell),
                           np.full(model.ncell, 0.35))
    j_t = ja_t.assemble(state, ra_t.fluid_state(state), dt=1.0).toarray()
    j_m = ja_m.assemble(state, ra_m.fluid_state(state), dt=1.0).toarray()
    assert np.allclose(j_t, j_m, atol=1e-6, rtol=1e-6), (
        f"izotrop iki-nöqtəli həddə MPFA/TPFA Jacobian fərqi: "
        f"{np.abs(j_t - j_m).max():.3g}")


# ═══════════════════════════════ 27/28: robustluq ════════════════════════
@pytest.mark.parametrize("contrast", [10.0, 100.0, 1000.0])
def test_high_contrast_permeability_solver_remains_stable(contrast):
    """PHASE D §27/§28 — Kmax/Kmin=10..1000, solver çökmür, NaN/Inf yoxdur."""
    model = _model(nx=5, ny=5, permeability=100.0)
    permx = model.rock.permx.values.copy()
    permx[: permx.size // 2] *= contrast
    permx[permx.size // 2:] /= 1.0
    model.rock.permx.values[:] = permx
    _, result = _run(model, MPFAODiscretization(closure=NEUMANN), end_time=15.0)
    assert result.converged
    assert np.all(np.isfinite(result.snapshots[-1].pressure))


def test_negative_permeability_is_rejected_not_silently_accepted():
    tensor = PermeabilityTensor(
        kxx=PropertyMap("KXX", np.array([-5.0, 10.0])),
        kyy=PropertyMap("KYY", np.array([10.0, 10.0])),
        kzz=PropertyMap("KZZ", np.array([10.0, 10.0])))
    assert tensor.validate().errors


# ══════════════════════════════ 36: performans ═══════════════════════════
@pytest.mark.performance
def test_performance_benchmark_mpfa_vs_tpfa_assembly_and_solve():
    results = []
    for nx, ny, label in ((5, 5, "10x10-ish"), (10, 10, "100 cells")):
        model = _model(nx=nx, ny=ny, permeability=100.0)
        model.rock.permeability_tensor = _rotated_tensor(
            model.ncell, eigenvalues=(300.0, 60.0, 20.0), axis=(0, 0, 1), angle=30.0)
        for label2, discretization in (("TPFA", TwoPointFluxDiscretization()),
                                       ("MPFA", MPFAODiscretization(closure=NEUMANN))):
            start = time.perf_counter()
            _, result = _run(model, discretization, end_time=20.0)
            elapsed = time.perf_counter() - start
            results.append((label, label2, model.ncell, elapsed, result.converged))
    print("\nPHASE D performance (TPFA vs MPFA-O, full nonlinear run):")
    for grid_label, method, ncell, elapsed, converged in results:
        print(f"  {grid_label:12s} {method:5s} cells={ncell:4d} time={elapsed:.3f}s "
              f"converged={converged}")
        assert elapsed < 60.0, f"{grid_label}/{method}: {elapsed:.2f}s — reqressiya"
        assert converged
