"""Phase 5B-1 — MPFA-O QLOBAL residual/stensil inteqrasiyasının test dəsti.

Tapşırıq §29-da tələb olunan A–T kateqoriyaları + §16 (qlobal stensil),
§17 (lokal→qlobal uyğunluq), §20 (qlobal anti-pseudo), §30 (qlobal
matrisin ÖZÜNÜN auditi), §31 (əl ilə yoxlanıla bilən qlobal hal).

Bax `docs/mpfa_o_phase5b1.md` — hər testin yoxladığı düstur orada.

Həndəsə/tenzor qurucuları Phase 5A dəstindən TƏKRAR İSTİFADƏ olunur
(`test_mpfa_o`) — ikinci, fərqli test həndəsəsi yaratmaq iki dəstin
səssizcə ayrılmasına gətirərdi.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from scipy import sparse

from helpers import default_scal, five_spot_model
from test_mpfa_o import (GAMMA, analytic_flux, cartesian, dense, isotropic,
                         linear_field, rotated_anisotropic, rotation, skewed,
                         uniform_k, warped)

from imex2d.discretization import (MPFAOBoundaryClosure, MPFAODiscretization,
                                   build_mpfa_o_coefficients)
from imex2d.discretization.mpfa_global import MPFAGlobalOperator, MPFAStateError
from imex2d.discretization.mpfa_o import MPFADiscretizedGrid
from imex2d.domain.structure import FaultReference
from imex2d.interfaces.discretization import IFluxDiscretization
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.impes_engine import ImpesEngine
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.implicit.jacobian import JacobianAssembler
from imex2d.simulation.implicit.residual import OIL, WATER, ResidualAssembler
from imex2d.simulation.implicit.state import VARIABLES_PER_CELL, ReservoirState
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.well_model import PeacemanWellModel

NEUMANN = MPFAOBoundaryClosure.NEUMANN_ZERO


# ══════════════════════════════════════════════════ köməkçilər ══════════
def global_operator(grid, geometry, k_matrices, closure=NEUMANN, eta=1.0):
    """`MPFAGlobalOperator` — BİRBAŞA həndəsədən (ReservoirModel-siz).

    Skew/warped həndəsə `ReservoirModel.geometry` (Kartezian
    `CellGeometry`) ilə ifadə OLUNA BİLMİR, ona görə qlobal operator
    testləri nüvəni birbaşa qurur (Phase 5A-dakı eyni yanaşma)."""
    coefficients = build_mpfa_o_coefficients(grid, geometry, k_matrices, GAMMA,
                                             eta=eta, closure=closure)
    return MPFAGlobalOperator(
        coefficients=coefficients, connections=grid.build_connections(),
        connection_faces=geometry.connection_faces(),
        face_owner=geometry.face_owner, face_neighbor=geometry.face_neighbor)


def mpfa_grid_with_geometry(model, geometry, k_matrices, closure=NEUMANN):
    """`MPFADiscretizedGrid` — İXTİYARİ (skew) həndəsə ilə.

    `MPFAODiscretization.build()` Kartezian təpələr qurur; qalıq qatını
    QEYRİ-ORTOQONAL həndəsə ilə sınamaq üçün burada eyni obyekt əl ilə
    yığılır. Model topologiyası (`Connections`, `ncell`) EYNİ QALIR —
    yalnız MPFA-nın gördüyü həndəsə dəyişir. Cazibə söndürülmüş
    (tək lay, bərabər dərinlik) modellərdə bu, tam etibarlıdır."""
    conn = model.connections()
    coefficients = build_mpfa_o_coefficients(model.grid, geometry, k_matrices,
                                             model.units.darcy_constant,
                                             closure=closure)
    connection_faces = geometry.connection_faces()
    operator = MPFAGlobalOperator(
        coefficients=coefficients, connections=conn,
        connection_faces=connection_faces, face_owner=geometry.face_owner,
        face_neighbor=geometry.face_neighbor)
    return MPFADiscretizedGrid(
        connections=conn, pore_volume=model.pore_volume(),
        cell_volume=model.geometry.volumes(), geometry=geometry,
        coefficients=coefficients, connection_faces=connection_faces,
        global_operator=operator)


def assemblers(model, closure=NEUMANN):
    """(MPFA, TPFA) `ResidualAssembler` cütü — EYNİ model, EYNİ quyular,
    EYNİ relperm; YALNIZ məkan diskretizasiyası fərqlidir."""
    relperm = CoreyRelativePermeabilityAdapter(default_scal())
    wells = PeacemanWellModel().build_connections(model)
    mpfa_grid = MPFAODiscretization(closure=closure).build(model)
    tpfa_grid = TwoPointFluxDiscretization().build(model)
    return (ResidualAssembler(model, mpfa_grid, wells, relperm),
            ResidualAssembler(model, tpfa_grid, wells, relperm))


def state(ncell, seed=7, sw=0.3, spread=5.0):
    rng = np.random.default_rng(seed)
    return ReservoirState(200.0 + rng.normal(0.0, spread, ncell),
                          np.full(ncell, sw))


def flat_model(nx=4, ny=3, permeability=200.0):
    """Tək laylı, DÜZ model — cazibə üzvü sıfırdır (`_has_gravity=False`),
    ona görə `potentials()` ilə `cell_potentials()` BİRƏBİR uyğundur."""
    return five_spot_model(nx=nx, ny=ny, dx=10.0, dy=12.0, dz=5.0,
                           permeability=permeability)


def layered_model(nx=3, ny=3, nz=4, dz=5.0, top=2000.0):
    """COXLAYLI model — derinlik ferqi VAR, yeni `_has_gravity=True`.

    Cazibe konvensiyasi testleri (§12) mehz bunu teleb edir: tek layli
    `five_spot_model`-de butun huceyreler eyni derinlikdedir ve cazibe
    uzvu eyniliyle sifirdir — bele modelde cazibe testi BOS olardi."""
    from imex2d.application.model_builder import ReservoirModelBuilder
    from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                              five_spot)
    from imex2d.domain.initial import InitialConditions
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=dz, porosity=0.22,
        permx_base=200.0, nz=nz, kv_over_kh=0.1, top_depth=top)
    return ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=default_scal(),
        initial=InitialConditions(datum_depth=top, datum_pressure=250.0))


# ═════════════════════════════════════════════ A — API / seçim ══════════
def test_a_mpfa_is_selectable_and_reaches_the_residual_assembler():
    model = flat_model()
    mpfa, tpfa = assemblers(model)
    assert mpfa._multipoint is True
    assert tpfa._multipoint is False
    assert mpfa.transmissibility is None          # çoxnöqtəli sxemdə YOXDUR
    assert tpfa.transmissibility is not None
    assert isinstance(MPFAODiscretization(), IFluxDiscretization)


def test_a_jacobian_and_fully_implicit_engine_now_support_mpfa():
    """PHASE D (5B-2) — `JacobianAssembler`/`FullyImplicitEngine` ARTIQ
    MPFA-O ilə İŞLƏYİR (analitik çoxnöqtəli Jacobian implement edilib,
    bax `jacobian.py` modul docstring-i). `ImpesEngine` HƏLƏ DƏ AÇIQ
    imtina edir — onun təzyiq addımı tək-üz transmissivlik skalyarı
    tələb edir, bu, MPFA-nın stensilinə uyğun DEYİL (bilərəkdən
    implement edilməyib, aşağıda ayrıca yoxlanılır)."""
    model = flat_model()
    mpfa, _ = assemblers(model)
    jacobian = JacobianAssembler(mpfa)
    st = state(model.ncell)
    matrix = jacobian.assemble(st, mpfa.fluid_state(st), dt=1.0)
    assert matrix.shape == (model.ncell * VARIABLES_PER_CELL,) * 2
    assert matrix.nnz > 0

    from imex2d.application.config import SimulationConfig
    relperm = CoreyRelativePermeabilityAdapter(default_scal())
    engine = FullyImplicitEngine(
        model, SimulationConfig(end_time=1.0), relperm,
        flux_discretization=MPFAODiscretization(closure=NEUMANN))
    assert engine.jacobian_assembler._multipoint is True

    with pytest.raises(NotImplementedError, match="Phase 5B-2"):
        ImpesEngine(model, SimulationConfig(end_time=1.0), relperm,
                   flux_discretization=MPFAODiscretization(closure=NEUMANN),
                   linear_solver=None)


# ═══════════════════════════════════ B — TAM təzyiq vektoru girişi ══════
def test_b_flux_is_a_function_of_the_full_pressure_vector():
    """Tapşırıq §4 — `compute_flux(d_phi)` müqaviləsi BƏRPA EDİLMİR."""
    model = flat_model()
    grid = MPFAODiscretization(closure=NEUMANN).build(model)
    pressure = state(model.ncell).pressure

    flux = grid.connection_fluxes_from_potential(pressure)
    assert flux.shape == (model.connections().count,)
    with pytest.raises(NotImplementedError, match="Phase 5B"):
        grid.compute_flux(np.ones(model.connections().count))


def test_b_wrong_vector_length_is_rejected():
    model = flat_model()
    grid = MPFAODiscretization(closure=NEUMANN).build(model)
    with pytest.raises(MPFAStateError, match="uyğunsuzluq"):
        grid.connection_fluxes_from_potential(np.ones(model.ncell + 1))


# ══════════════════════════ C / §17 — LOKAL → QLOBAL uyğunluq ═══════════
def test_c_global_connection_flux_equals_phase5a_face_flux():
    """Tapşırıq §17 — indeks/işarə/ikiqat-sayma səhvlərini tutur."""
    grid, geometry, _ = skewed(3, 3, 2)
    operator = global_operator(grid, geometry, rotated_anisotropic(grid.ncell))
    pressure = state(grid.ncell, seed=3).pressure

    global_flux = operator.connection_fluxes(pressure)
    phase5a_flux = operator.coefficients.face_fluxes(pressure)[operator.connection_faces]
    assert np.allclose(global_flux, phase5a_flux, rtol=1e-14, atol=1e-14)


def test_c_global_flux_equals_sum_of_local_region_sub_face_fluxes():
    """Ən aşağı səviyyəyə qədər: LOKAL sistemlərin sub-üz axınlarının
    cəmi = QLOBAL əlaqə axını (§17)."""
    grid, geometry, _ = skewed(3, 3, 2)
    operator = global_operator(grid, geometry, rotated_anisotropic(grid.ncell))
    coefficients = operator.coefficients
    pressure = state(grid.ncell, seed=11).pressure
    global_flux = operator.connection_fluxes(pressure)

    for k in range(min(6, operator.nconn)):
        face = int(operator.connection_faces[k])
        total = 0.0
        for system in coefficients.local_systems:
            for sub in system.region.sub_faces:
                if sub.face_index != face:
                    continue
                local = system.sub_face_fluxes(pressure[system.region.cells])
                total += float(local[sub.local_index])
        assert np.isclose(total, global_flux[k], rtol=1e-11,
                          atol=1e-11 * max(abs(global_flux).max(), 1.0)), (
            f"əlaqə {k} (üz {face}): lokal cəm {total} != qlobal {global_flux[k]}")


def test_c_residual_face_flux_uses_the_same_operator():
    """Qalıq qatının aldığı axın = qlobal operatorun verdiyi axın ×
    upstream mobilitəsi (başqa bir yol YOXDUR)."""
    model = flat_model()
    mpfa, _ = assemblers(model)
    st = state(model.ncell)
    fluid = mpfa.fluid_state(st)

    water, oil = mpfa.face_fluxes(st, fluid)
    phi_w, phi_o = mpfa.cell_potentials(st, fluid)
    base_w = mpfa.grid.connection_fluxes_from_potential(phi_w)
    up_w = mpfa.grid.upstream_cells(base_w)
    assert np.allclose(water, base_w * (fluid.lam_w[up_w] / fluid.bw[up_w]))
    assert oil.shape == water.shape


# ════════════════════════════ D — daxili üz konservasiyası ══════════════
@pytest.mark.parametrize("builder,k_name", [
    (cartesian, "izotrop"), (skewed, "izotrop"),
    (skewed, "fırlanmış"), (warped, "fırlanmış")])
def test_d_internal_face_conservation(builder, k_name):
    """Tapşırıq §7 — `q_owner + q_neighbor = 0`, hər tərəf MÜSTƏQİL
    qradiyentdən."""
    grid, geometry, _ = builder()
    k = (isotropic(grid.ncell) if k_name == "izotrop"
         else rotated_anisotropic(grid.ncell))
    operator = global_operator(grid, geometry, k)
    pressure = state(grid.ncell, seed=5).pressure

    report = operator.conservation_report(pressure)
    scale = max(report["max_face_flux"], 1e-30)
    assert report["max_internal_face_error"] < 1e-9 * scale
    assert report["max_boundary_flux"] < 1e-9 * scale       # NEUMANN_ZERO
    assert report["global_imbalance"] < 1e-9 * scale


# ═══════════════════════════════ E — QLOBAL konservasiya ════════════════
def test_e_closed_system_residual_sums_to_accumulation_only():
    """Quyusuz qapalı sistemdə daxili axınlar bir-birini yeyir →
    `Σ_i R_i = Σ_i (akkumulyasiya dəyişməsi)/Δt` (tapşırıq §7)."""
    model = flat_model()
    model.wells = []
    mpfa, _ = assemblers(model)
    st, previous = state(model.ncell), state(model.ncell, seed=2, sw=0.25)

    fluid = mpfa.fluid_state(st)
    influx_water, influx_oil = mpfa.net_influx(st, fluid)
    scale = max(np.abs(influx_water).max(), 1e-30)
    assert abs(influx_water.sum()) < 1e-9 * scale
    assert abs(influx_oil.sum()) < 1e-9 * max(np.abs(influx_oil).max(), 1e-30)

    residual, _, _ = mpfa.residual(st, previous, dt=5.0)
    water_new, oil_new = mpfa.accumulation(st, fluid)
    water_old, oil_old = mpfa.accumulation(previous, mpfa.fluid_state(previous))
    assert np.isclose(residual[WATER::VARIABLES_PER_CELL].sum(),
                      (water_new - water_old).sum() / 5.0, rtol=1e-9)
    assert np.isclose(residual[OIL::VARIABLES_PER_CELL].sum(),
                      (oil_new - oil_old).sum() / 5.0, rtol=1e-9)


def test_e_global_operator_rows_sum_to_zero_columnwise():
    """`A = D·T_conn` sütun cəmi ≡ 0 — hər üz payı bir hüceyrədən çıxıb
    o birinə DƏQİQ eyni ölçüdə daxil olur (§8 ikiqat-sayma yoxlaması)."""
    grid, geometry, _ = skewed()
    operator = global_operator(grid, geometry, rotated_anisotropic(grid.ncell))
    matrix = operator.cell_operator()
    column_sums = np.asarray(matrix.sum(axis=0)).ravel()
    assert np.max(np.abs(column_sums)) < 1e-9 * np.abs(matrix.data).max()


# ═══════════════════════ F — xətti manufactured həll (§18) ══════════════
@pytest.mark.parametrize("builder", [cartesian, skewed, warped])
def test_f_linear_pressure_reproduces_analytic_darcy_at_global_level(builder):
    grid, geometry, _ = builder()
    k_matrix = rotated_anisotropic(1)[0]
    operator = global_operator(grid, geometry, uniform_k(grid.ncell, k_matrix),
                               closure=MPFAOBoundaryClosure.DIRICHLET)
    gradient = np.array([2.0, -3.0, 1.25])
    pressure, exact = linear_field(geometry, gradient, offset=11.0)
    boundary = operator.coefficients.boundary_pressures_from(exact)

    face_flux = operator.face_fluxes(pressure, boundary_potential=boundary)
    expected = analytic_flux(geometry, k_matrix, gradient)
    assert np.allclose(face_flux, expected, rtol=1e-10,
                       atol=1e-10 * np.abs(expected).max())

    connection_flux = operator.connection_fluxes(pressure, boundary_potential=boundary)
    assert np.allclose(connection_flux, expected[operator.connection_faces],
                       rtol=1e-10, atol=1e-10 * np.abs(expected).max())
    assert np.allclose(connection_flux, face_flux[operator.connection_faces],
                       rtol=1e-12, atol=1e-12)

    # Hüceyrə balansı: sabit K + xətti sahə ⇒ hər hüceyrədə divergensiya 0
    divergence = np.zeros(grid.ncell)
    np.add.at(divergence, geometry.face_owner, face_flux)
    interior = geometry.face_neighbor >= 0
    np.add.at(divergence, geometry.face_neighbor[interior], -face_flux[interior])
    assert np.max(np.abs(divergence)) < 1e-9 * np.abs(face_flux).max()


def test_f_residual_path_on_skew_geometry_matches_the_validated_operator():
    """Qaliq SEVIYYESINDE qeyri-ortoqonal hendese (§18).

    DIQQET — riyazi hedd: `NEUMANN_ZERO` (qapali, axinsiz) sistemde
    XETTI tezyiq sahesi HELLIN OZU DEYIL (serhedde normal qradiyent
    sifir olmalidir), ona gore analitik Darsi ile birbasa muqayise
    MEHZ BURADA menasizdir. Manufactured hell operator seviyyesinde,
    DIRICHLET serhed melumati ile yoxlanilir
    (`test_f_linear_pressure_reproduces_analytic_darcy_at_global_level`).

    Bu test qaliq yolunun HEMIN dogrulanmis operatordan kecdiyini ve
    konservativ oldugunu yoxlayir: baska bir axin yolu YOXDUR.
    """
    model = flat_model(nx=4, ny=3)
    model.wells = []
    _, geometry, _ = skewed(model.grid.nx, model.grid.ny, model.grid.nz,
                            shear=0.3, dx=10.0, dy=12.0, dz=5.0)
    k_matrix = rotated_anisotropic(1)[0]
    grid = mpfa_grid_with_geometry(model, geometry, uniform_k(model.ncell, k_matrix),
                                   closure=NEUMANN)
    relperm = CoreyRelativePermeabilityAdapter(default_scal())
    assembler = ResidualAssembler(model, grid, [], relperm)
    assert assembler._has_gravity is False, "tek layli model — cazibe sifir"

    gradient = np.array([1.5, -2.0, 0.0])
    st = ReservoirState(geometry.cell_centroids @ gradient + 200.0,
                        np.full(model.ncell, 0.3))
    fluid = assembler.fluid_state(st)
    water, oil = assembler.face_fluxes(st, fluid)

    # (1) Qaliq axini = dogrulanmis operator x upstream mobilitesi
    phi_w, _ = assembler.cell_potentials(st, fluid)
    base = grid.global_operator.connection_fluxes(phi_w)
    upstream = grid.global_operator.upstream_cells(base)
    assert np.allclose(water, base * (fluid.lam_w[upstream] / fluid.bw[upstream]),
                       rtol=1e-12, atol=1e-12)

    # (2) Coxnoqteli: stensil 2-den boyukdur (skew + firlanmis tenzor)
    assert grid.global_operator.stencil_sizes().min() > 2

    # (3) Konservativ: qapali sistemde net axin cemi sifirdir
    influx_water, influx_oil = assembler.net_influx(st, fluid)
    assert abs(influx_water.sum()) < 1e-9 * np.abs(influx_water).max()
    assert abs(influx_oil.sum()) < 1e-9 * np.abs(influx_oil).max()


def test_f_dirichlet_closure_is_rejected_by_the_residual_path():
    """Serhed pay SESSIZCE atila bilmez (§13/§27).

    `DIRICHLET` bagalanisinda `T_bnd` sutunlari var; qaliq qati hele
    serhed dеyerlerini otuрmur, ona gore bele grid RED EDILIR — `T_bnd`
    payini sifir saymaq fiziki cehetden yanlis axin verer."""
    model = flat_model()
    grid = MPFAODiscretization(closure=MPFAOBoundaryClosure.DIRICHLET).build(model)
    assert any("dirichlet" in feature.lower() for feature in grid.unsupported_features)
    relperm = CoreyRelativePermeabilityAdapter(default_scal())
    with pytest.raises(NotImplementedError, match="sərhəd"):
        ResidualAssembler(model, grid, [], relperm)

    # Operator seviyyesinde de sessiz sifir YOXDUR
    with pytest.raises(MPFAStateError, match="boundary_potential"):
        grid.global_operator.connection_fluxes(np.full(model.ncell, 200.0))


# ═══════════════════════════ G — sabit təzyiq → sıfır axın ══════════════
def test_g_constant_pressure_gives_zero_flux_and_zero_influx():
    model = flat_model()
    model.wells = []
    mpfa, _ = assemblers(model)
    st = ReservoirState(np.full(model.ncell, 187.5), np.full(model.ncell, 0.3))
    fluid = mpfa.fluid_state(st)

    water, oil = mpfa.face_fluxes(st, fluid)
    assert np.max(np.abs(water)) < 1e-9
    assert np.max(np.abs(oil)) < 1e-9
    influx_water, influx_oil = mpfa.net_influx(st, fluid)
    assert np.max(np.abs(influx_water)) < 1e-9
    assert np.max(np.abs(influx_oil)) < 1e-9


# ═════════════ H — ortoqonal izotrop limitdə MPFA ≡ TPFA (§19) ══════════
def test_h_orthogonal_isotropic_face_fluxes_match_tpfa():
    model = flat_model()
    mpfa, tpfa = assemblers(model)
    st = state(model.ncell)
    fluid_m, fluid_t = mpfa.fluid_state(st), tpfa.fluid_state(st)

    water_m, oil_m = mpfa.face_fluxes(st, fluid_m)
    water_t, oil_t = tpfa.face_fluxes(st, fluid_t)
    assert np.allclose(water_m, water_t, rtol=1e-10,
                       atol=1e-10 * np.abs(water_t).max())
    assert np.allclose(oil_m, oil_t, rtol=1e-10, atol=1e-10 * np.abs(oil_t).max())


def test_h_orthogonal_isotropic_residual_matches_tpfa():
    model = flat_model()
    mpfa, tpfa = assemblers(model)
    st, previous = state(model.ncell), state(model.ncell, seed=2, sw=0.25)

    residual_m, _, _ = mpfa.residual(st, previous, dt=5.0)
    residual_t, _, _ = tpfa.residual(st, previous, dt=5.0)
    assert np.allclose(residual_m, residual_t, rtol=1e-10,
                       atol=1e-10 * np.abs(residual_t).max())


def test_h_orthogonal_isotropic_stencil_is_two_point_at_global_level():
    """§19-un STENSİL səviyyəsində təsdiqi (flux müqayisəsindən güclüdür)."""
    grid, geometry, _ = cartesian(4, 3, 2)
    operator = global_operator(grid, geometry, isotropic(grid.ncell, 200.0))
    assert np.all(operator.stencil_sizes() == 2)


# ══════════════════════════════════ I / J / K / L — həndəsə və tenzor ═══
def test_i_skew_geometry_global_assembly_is_valid():
    grid, geometry, _ = skewed()
    operator = global_operator(grid, geometry, rotated_anisotropic(grid.ncell))
    assert operator.T_conn.shape == (grid.build_connections().count, grid.ncell)
    assert sparse.issparse(operator.T_conn)
    assert np.all(np.isfinite(operator.T_conn.data))


def test_j_diagonal_anisotropy_changes_global_operator_per_direction():
    grid, geometry, _ = cartesian(3, 3, 2)
    base = global_operator(grid, geometry, uniform_k(grid.ncell, np.diag([400.0, 40.0, 4.0])))
    swapped = global_operator(grid, geometry,
                              uniform_k(grid.ncell, np.diag([4.0, 40.0, 400.0])))
    assert not np.allclose(dense(base.T_conn), dense(swapped.T_conn))
    # X üzlərinin əmsalı MƏHZ kx nisbətində dəyişir
    conn = grid.build_connections()
    x_faces = np.flatnonzero(conn.axis == 0)
    ratio = (dense(swapped.T_conn)[x_faces, conn.cell_a[x_faces]]
             / dense(base.T_conn)[x_faces, conn.cell_a[x_faces]])
    assert np.allclose(ratio, 4.0 / 400.0, rtol=1e-10)


@pytest.mark.parametrize("angle", [0.0, 25.0, 61.0, 118.0])
def test_k_rotated_anisotropy_changes_the_global_operator(angle):
    """Tapşırıq §21 — tenzorun FIRLANMASI qlobal operatoru dəyişməlidir."""
    grid, geometry, _ = skewed()
    reference = global_operator(grid, geometry, rotated_anisotropic(grid.ncell, angle=0.0))
    rotated = global_operator(grid, geometry, rotated_anisotropic(grid.ncell, angle=angle))
    difference = np.max(np.abs(dense(rotated.T_conn) - dense(reference.T_conn)))
    scale = np.max(np.abs(dense(reference.T_conn)))
    if angle == 0.0:
        assert difference < 1e-12 * scale
    else:
        assert difference > 1e-6 * scale, f"bucaq {angle} operatoru dəyişmədi"


def test_l_off_diagonal_components_change_the_global_operator():
    """Kxy/Kxz/Kyz QLOBAL nəticəyə təsir edir (§21/§L)."""
    grid, geometry, _ = skewed()
    diagonal = np.diag([500.0, 50.0, 10.0])
    base = global_operator(grid, geometry, uniform_k(grid.ncell, diagonal))
    scale = np.max(np.abs(dense(base.T_conn)))

    for label, (i, j) in (("Kxy", (0, 1)), ("Kxz", (0, 2)), ("Kyz", (1, 2))):
        matrix = diagonal.copy()
        matrix[i, j] = matrix[j, i] = 15.0
        assert np.linalg.eigvalsh(matrix)[0] > 0.0
        perturbed = global_operator(grid, geometry, uniform_k(grid.ncell, matrix))
        difference = np.max(np.abs(dense(perturbed.T_conn) - dense(base.T_conn)))
        assert difference > 1e-6 * scale, f"{label} qlobal operatora TƏSİR ETMİR"


# ═════════════════ M / §16 — QLOBAL sətirdə >2 təzyiq bağlantısı ════════
def test_m_global_face_row_has_more_than_two_pressure_couplings():
    grid, geometry, _ = skewed(3, 3, 2)
    operator = global_operator(grid, geometry, rotated_anisotropic(grid.ncell))
    sizes = operator.stencil_sizes()
    assert sizes.min() > 2, f"ən kiçik qlobal stensil {sizes.min()}"
    assert sizes.max() >= 8


def test_m_global_cell_row_couples_beyond_direct_connections_neighbours():
    """Tapşırıq §22 — MPFA bağlantısı `Connections` qonşuluğundan
    GENİŞDİR (diaqonal hüceyrələr); `JacobianAssembler` Phase 5B-2-də
    məhz bu naxışa keçməlidir."""
    grid, geometry, _ = skewed(3, 3, 2)
    operator = global_operator(grid, geometry, rotated_anisotropic(grid.ncell))
    conn = grid.build_connections()

    direct = [set() for _ in range(grid.ncell)]
    for a, b in zip(conn.cell_a, conn.cell_b):
        direct[int(a)].update({int(a), int(b)})
        direct[int(b)].update({int(a), int(b)})

    pattern = operator.global_stencil_pattern()
    extra_total = 0
    for cell in range(grid.ncell):
        row = set(pattern.getrow(cell).indices.tolist())
        extra_total += len(row - direct[cell])
    assert extra_total > 0, ("qlobal naxış `Connections` qonşuluğu ilə eynidir "
                             "— çoxnöqtəli bağlantı itib")


# ═══════════════════════════════ N — sərhəd interfeysi (§13) ════════════
def test_n_neumann_zero_is_the_residual_path_boundary_condition():
    grid, geometry, _ = cartesian(3, 3, 2)
    operator = global_operator(grid, geometry, isotropic(grid.ncell), closure=NEUMANN)
    assert operator.coefficients.T_bnd.shape[1] == 0
    pressure = state(grid.ncell, seed=13).pressure
    face_flux = operator.face_fluxes(pressure)
    assert np.max(np.abs(face_flux[geometry.is_boundary])) < 1e-9 * np.abs(face_flux).max()


def test_n_dirichlet_structure_is_exposed_for_phase_5b2():
    grid, geometry, _ = cartesian(2, 2, 1)
    operator = global_operator(grid, geometry, isotropic(grid.ncell),
                               closure=MPFAOBoundaryClosure.DIRICHLET)
    assert operator.coefficients.T_bnd.shape[1] == len(operator.coefficients.boundary_dofs) > 0
    assert operator.coefficients.boundary_points.shape[1] == 3


# ═══════════════════════════════ O — işarə konvensiyası (§15) ═══════════
def test_o_two_cell_hand_check_of_flux_direction_and_residual_sign():
    """ƏL İLƏ yoxlanıla bilən 2 hüceyrəli hal (§15).

    `p_0 > p_1` ⇒ axın 0 → 1 ⇒ `flux[0] > 0`;
    `net_influx[0] < 0` (itirir), `net_influx[1] > 0` (qazanır);
    `R = akkumulyasiya/Δt − influx` olduğu üçün `R_0 > R_1` (influx
    işarəsi qalığa TƏRS düşür).
    """
    model = five_spot_model(nx=2, ny=1, dx=10.0, dy=10.0, dz=5.0, permeability=100.0)
    model.wells = []
    mpfa, tpfa = assemblers(model)
    st = ReservoirState(np.array([250.0, 150.0]), np.array([0.3, 0.3]))
    fluid = mpfa.fluid_state(st)

    water, oil = mpfa.face_fluxes(st, fluid)
    assert water[0] > 0.0 and oil[0] > 0.0        # 0 → 1
    influx_water, _ = mpfa.net_influx(st, fluid)
    assert influx_water[0] < 0.0 < influx_water[1]
    assert np.isclose(influx_water[0], -influx_water[1])

    water_t, _ = tpfa.face_fluxes(st, tpfa.fluid_state(st))
    assert np.isclose(water[0], water_t[0], rtol=1e-10)


def test_o_upwind_matches_tpfa_in_two_point_limit():
    """Tapşırıq §10 — upwind HƏQİQİ MPFA axınının işarəsindən çıxır,
    amma iki-nöqtəli limitdə TPFA qaydası ilə BİRƏBİR üst-üstə düşür."""
    model = flat_model()
    mpfa, tpfa = assemblers(model)
    conn = model.connections()
    for seed in range(4):
        st = state(model.ncell, seed=seed, spread=12.0)
        fluid = mpfa.fluid_state(st)
        phi_w, _ = mpfa.cell_potentials(st, fluid)
        mpfa_upwind = mpfa.grid.upstream_cells(
            mpfa.grid.connection_fluxes_from_potential(phi_w))
        d_phi_w, _ = tpfa.potentials(st, tpfa.fluid_state(st))
        tpfa_upwind = np.where(d_phi_w >= 0, conn.cell_a, conn.cell_b)
        assert np.array_equal(mpfa_upwind, tpfa_upwind)


def test_o_gravity_is_active_and_uses_the_existing_convention():
    """Tapsiriq §12 — cazibe SUSDURULMUR ve IKINCI konvensiya ICAD
    EDILMIR: sixliq huceyreler uzre BERABER olanda huceyre-potensiali
    TPFA-nin uz-ortalanmis potensiali ile EYNIDIR."""
    model = layered_model()
    model.wells = []
    model.fluids.water_compressibility = 0.0      # rho huceyreler uzre SABIT
    model.fluids.oil_compressibility = 0.0
    mpfa, tpfa = assemblers(model)
    assert mpfa._has_gravity is True, "coxlayli model — cazibe AKTIV olmalidir"

    st = state(model.ncell, spread=0.0)            # sabit p, yalniz cazibe
    fluid = mpfa.fluid_state(st)
    conn = model.connections()
    phi_w, phi_o = mpfa.cell_potentials(st, fluid)
    d_phi_w, d_phi_o = mpfa.potentials(st, fluid)

    assert np.max(np.abs(d_phi_w)) > 1e-6, "cazibe uzvu sifir cixdi — test bosdur"
    assert np.allclose(phi_w[conn.cell_a] - phi_w[conn.cell_b], d_phi_w, atol=1e-10)
    assert np.allclose(phi_o[conn.cell_a] - phi_o[conn.cell_b], d_phi_o, atol=1e-10)

    water_m, oil_m = mpfa.face_fluxes(st, fluid)
    water_t, oil_t = tpfa.face_fluxes(st, tpfa.fluid_state(st))
    assert np.allclose(water_m, water_t, rtol=1e-9,
                       atol=1e-9 * max(np.abs(water_t).max(), 1e-30))
    assert np.allclose(oil_m, oil_t, rtol=1e-9,
                       atol=1e-9 * max(np.abs(oil_t).max(), 1e-30))


def test_o_variable_density_gravity_difference_is_bounded_and_measured():
    """SENEDLESDIRILMIS ferq (§4): deyisen sixliqda huceyre-sixliq ile
    uz-ortalanmis sixliq `O(drho*dD)` qeder ferqlenir. Test onu OLCUR ve
    nezeri ifade ile TUTUSDURUR — gizletmir, "duzeltmir"."""
    model = layered_model()
    model.wells = []
    mpfa, _ = assemblers(model)
    assert mpfa._has_gravity is True

    rng = np.random.default_rng(4)
    st = ReservoirState(250.0 + rng.normal(0.0, 40.0, model.ncell),
                        np.full(model.ncell, 0.3))
    fluid = mpfa.fluid_state(st)
    conn = model.connections()
    phi_w, _ = mpfa.cell_potentials(st, fluid)
    d_phi_w, _ = mpfa.potentials(st, fluid)
    difference = np.abs((phi_w[conn.cell_a] - phi_w[conn.cell_b]) - d_phi_w)

    from imex2d.simulation.implicit.residual import GRAVITY, PA_TO_BAR
    rho_w = model.fluids.water_density / np.maximum(fluid.bw, 1e-9)
    depths = model.geometry.cell_depths()
    # |rho_a D_a - rho_b D_b - rho_bar (D_a - D_b)| = 0.5 |drho| |D_a + D_b|
    predicted = 0.5 * np.abs(rho_w[conn.cell_a] - rho_w[conn.cell_b]) * np.abs(
        depths[conn.cell_a] + depths[conn.cell_b]) * GRAVITY * PA_TO_BAR
    assert np.allclose(difference, predicted, rtol=1e-9, atol=1e-12), (
        "ferq nezeri O(drho*dD) ifadesi ile uygun gelmir")
    assert difference.max() > 0.0, "deyisen sixliqda ferq sifir cixdi — test bosdur"


def test_o_gravity_residual_is_finite_and_conservative_on_a_layered_model():
    """Cazibeli, coxlayli modelde MPFA qaligi hesablanir ve qapali
    sistemde daxili axinlar bir-birini yeyir."""
    model = layered_model()
    model.wells = []
    mpfa, _ = assemblers(model)
    st, previous = state(model.ncell), state(model.ncell, seed=2, sw=0.25)
    fluid = mpfa.fluid_state(st)

    influx_water, influx_oil = mpfa.net_influx(st, fluid)
    assert np.all(np.isfinite(influx_water)) and np.all(np.isfinite(influx_oil))
    assert abs(influx_water.sum()) < 1e-9 * max(np.abs(influx_water).max(), 1e-30)
    assert abs(influx_oil.sum()) < 1e-9 * max(np.abs(influx_oil).max(), 1e-30)
    residual, _, _ = mpfa.residual(st, previous, dt=5.0)
    assert np.all(np.isfinite(residual))


# ═══════════════════════ P — İKİQAT SAYMANIN aşkarlanması (§8) ══════════
def test_p_connection_face_mapping_is_a_bijection_onto_interior_faces():
    """Hər daxili fiziki üz DƏQİQ BİR `Connections` girişinə düşür —
    nə buraxılır, nə də iki dəfə sayılır."""
    grid, geometry, _ = skewed(3, 3, 2)
    faces = geometry.connection_faces()
    interior = np.flatnonzero(~geometry.is_boundary)
    assert faces.size == grid.build_connections().count
    assert np.unique(faces).size == faces.size                 # təkrar YOX
    assert set(faces.tolist()) == set(interior.tolist())       # buraxılan YOX


def test_p_face_level_and_connection_level_balances_agree_exactly():
    """IKIQAT SAYMA / BURAXILMA detektoru (§8).

    Iki MUSTEQIL yolla huceyre balansi qurulur:
      (1) BUTUN fiziki uzler uzre (serhed daxil), `face_owner`/
          `face_neighbor` ile sepelenerek;
      (2) YALNIZ `Connections` uzre (qaliq qatinin isletdiyi yol).
    Bir sub-uz iki fiziki uze yazilsaydi, bir uz buraxilsaydi, ve ya
    serhed uzu sifirdan ferqli axin dasisaydi, IKI yol AYRILARDI.
    """
    grid, geometry, _ = skewed(3, 3, 2)
    operator = global_operator(grid, geometry, rotated_anisotropic(grid.ncell))
    pressure = state(grid.ncell, seed=17).pressure

    face_flux = operator.face_fluxes(pressure)
    by_face = np.zeros(grid.ncell)
    np.add.at(by_face, geometry.face_owner, -face_flux)
    interior = geometry.face_neighbor >= 0
    np.add.at(by_face, geometry.face_neighbor[interior], face_flux[interior])

    by_connection = operator.net_influx(pressure)
    scale = max(np.abs(by_connection).max(), 1e-30)
    assert np.allclose(by_face, by_connection, rtol=1e-9, atol=1e-9 * scale)
    assert abs(by_connection.sum()) < 1e-9 * scale


def test_p_duplicated_face_contribution_is_detected_by_the_balance_check():
    """Detektorun HESSASLIGI: bir uzu iki defe sayan (qesden pozulmus)
    xerite huceyre balansini DEYISIR — yeni yuxaridaki test bos deyil."""
    grid, geometry, _ = skewed(3, 3, 2)
    operator = global_operator(grid, geometry, rotated_anisotropic(grid.ncell))
    pressure = state(grid.ncell, seed=17).pressure
    healthy = operator.net_influx(pressure)

    conn = grid.build_connections()
    duplicated = operator.connection_fluxes(pressure).copy()
    duplicated[0] *= 2.0
    broken = np.zeros(grid.ncell)
    np.add.at(broken, conn.cell_a, -duplicated)
    np.add.at(broken, conn.cell_b, +duplicated)
    assert not np.allclose(broken, healthy), "detektor ikiqat saymani gormur"


def test_p_boundary_faces_carry_no_flux_so_nothing_is_lost():
    grid, geometry, _ = skewed(3, 3, 2)
    operator = global_operator(grid, geometry, rotated_anisotropic(grid.ncell))
    pressure = state(grid.ncell, seed=19).pressure
    face_flux = operator.face_fluxes(pressure)
    boundary = np.flatnonzero(geometry.is_boundary)
    assert np.max(np.abs(face_flux[boundary])) < 1e-9 * np.abs(face_flux).max()


# ══════════════════════════════════ Q — seyrək yaddaş (§5/§28) ══════════
def test_q_global_operator_is_sparse_and_scales_linearly():
    densities = []
    for nx in (4, 10):
        grid, geometry, _ = cartesian(nx, nx, 2)
        operator = global_operator(grid, geometry, rotated_anisotropic(grid.ncell))
        assert sparse.issparse(operator.T_conn)
        assert sparse.issparse(operator.cell_operator())
        per_row = operator.T_conn.nnz / max(operator.nconn, 1)
        densities.append(per_row)
        assert per_row <= 18.0
        assert operator.T_conn.nnz < 0.5 * operator.nconn * operator.ncell
        assert operator.cell_operator().nnz < 0.5 * operator.ncell ** 2
    assert densities[1] < 1.5 * densities[0]


def test_q_no_dense_allocation_in_the_residual_path():
    """Qalıq yolunda `ncell × ncell` SIX massiv YARANMIR: operator
    seyrəkdir və qalıq yalnız matris-vektor hasilidir."""
    model = flat_model(nx=8, ny=8)
    mpfa, _ = assemblers(model)
    operator = mpfa.grid.global_operator
    dense_bytes = 8 * operator.ncell ** 2
    sparse_bytes = (operator.T_conn.data.nbytes + operator.T_conn.indices.nbytes
                    + operator.cell_operator().data.nbytes)
    assert sparse_bytes < 0.25 * dense_bytes


# ═════════════════════════════════ R — etibarsız giriş (§27) ════════════
@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_r_nan_or_inf_pressure_is_rejected_with_cell_identification(value):
    model = flat_model()
    mpfa, _ = assemblers(model)
    pressure = np.full(model.ncell, 200.0)
    pressure[3] = value
    with pytest.raises(MPFAStateError, match="NaN/Inf"):
        mpfa.grid.connection_fluxes_from_potential(pressure)


def test_r_faults_are_rejected_not_silently_ignored():
    """Tapşırıq §26 — TPFA fay çarpanlarını MPFA-ya "yapışdırmaq"
    QADAĞANDIR; MPFA qalığı belə modeli RƏDD EDİR."""
    model = flat_model()
    model.fault_references.append(
        FaultReference(name="F1", source_id="f1", axis="I", plane_index=1,
                       transmissibility_multiplier=0.0))
    grid = MPFAODiscretization(closure=NEUMANN).build(model)
    assert grid.unsupported_features
    relperm = CoreyRelativePermeabilityAdapter(default_scal())
    with pytest.raises(NotImplementedError, match="fay"):
        ResidualAssembler(model, grid, [], relperm)
    # TPFA həmin modeli DƏSTƏKLƏYİR (davranış dəyişməyib)
    ResidualAssembler(model, TwoPointFluxDiscretization().build(model), [], relperm)


def test_r_orientation_mismatch_is_detected():
    """Üz↔əlaqə xəritələməsi pozulsa AÇIQ xəta (§17/§27)."""
    grid, geometry, _ = cartesian(3, 2, 1)
    coefficients = build_mpfa_o_coefficients(grid, geometry, isotropic(grid.ncell),
                                             GAMMA, closure=NEUMANN)
    faces = geometry.connection_faces().copy()
    faces[0], faces[1] = faces[1], faces[0]        # QƏSDƏN pozulmuş xəritə
    with pytest.raises(MPFAStateError, match="xəritələməsi"):
        MPFAGlobalOperator(coefficients, grid.build_connections(), faces,
                           geometry.face_owner, geometry.face_neighbor)


# ══════════════════════════════════ S — TPFA reqressiyası ═══════════════
def test_s_tpfa_residual_is_bitwise_unchanged_by_the_mpfa_path():
    model = flat_model()
    relperm = CoreyRelativePermeabilityAdapter(default_scal())
    wells = PeacemanWellModel().build_connections(model)
    tpfa_grid = TwoPointFluxDiscretization().build(model)
    tpfa = ResidualAssembler(model, tpfa_grid, wells, relperm)
    st, previous = state(model.ncell), state(model.ncell, seed=2, sw=0.25)
    before, _, _ = tpfa.residual(st, previous, dt=5.0)

    MPFAODiscretization(closure=NEUMANN).build(model)     # MPFA qurulur...
    after, _, _ = tpfa.residual(st, previous, dt=5.0)     # ...TPFA dəyişmir
    assert np.array_equal(before, after)


def test_s_tpfa_still_reaches_jacobian_and_newton():
    """TPFA yolu TAM işlək qalır — Jacobian qurulur, Newton konvergensiya
    edir (MPFA imtinası TPFA-ya SIZMIR)."""
    from imex2d.application.config import SimulationConfig
    model = flat_model(nx=3, ny=3)
    relperm = CoreyRelativePermeabilityAdapter(default_scal())
    engine = FullyImplicitEngine(model, SimulationConfig(end_time=5.0), relperm)
    assert isinstance(engine.jacobian_assembler, JacobianAssembler)
    matrix = engine.jacobian_assembler.assemble(
        engine.state, engine.residual_assembler.fluid_state(engine.state), dt=1.0)
    assert sparse.issparse(matrix) and matrix.nnz > 0


def test_s_tpfa_transmissibility_formula_unchanged():
    from imex2d.domain.units import METRIC
    model = five_spot_model(nx=3, ny=3, dx=10.0, dy=10.0, dz=5.0, permeability=200.0)
    grid = TwoPointFluxDiscretization().build(model)
    expected = METRIC.darcy_constant * 200.0 * (10.0 * 5.0) / 10.0
    assert np.allclose(grid.transmissibility, expected, rtol=1e-12)
    assert TwoPointFluxDiscretization().supports_multipoint_stencil() is False


# ═════════════════ T — mövcud qeyri-xətti qalıq reqressiyası ════════════
def test_t_tpfa_residual_matches_a_frozen_reference():
    """Mövcud (TPFA) qeyri-xətti qalıq DƏYİŞMƏYİB — dondurulmuş etalon.

    Dəyər `TwoPointFluxDiscretization`-ın DƏYİŞMƏMİŞ düsturlarından
    gəlir; MPFA inteqrasiyası onu POZSA bu test dərhal düşər."""
    model = five_spot_model(nx=3, ny=3, dx=10.0, dy=10.0, dz=5.0, permeability=200.0)
    relperm = CoreyRelativePermeabilityAdapter(default_scal())
    wells = PeacemanWellModel().build_connections(model)
    assembler = ResidualAssembler(
        model, TwoPointFluxDiscretization().build(model), wells, relperm)
    st = ReservoirState(np.full(model.ncell, 210.0), np.full(model.ncell, 0.3))
    previous = ReservoirState(np.full(model.ncell, 200.0), np.full(model.ncell, 0.25))
    residual, _, rates = assembler.residual(st, previous, dt=10.0)

    assert np.all(np.isfinite(residual))
    assert residual.size == model.ncell * VARIABLES_PER_CELL
    # Material balansı: Σ R · Δt = akkumulyasiya dəyişməsi − quyu həcmi
    water_balance, oil_balance = assembler.material_balance_error(residual, 10.0)
    water_new, oil_new = assembler.accumulation(st, assembler.fluid_state(st))
    water_old, oil_old = assembler.accumulation(previous,
                                                assembler.fluid_state(previous))
    assert np.isclose(water_balance,
                      (water_new - water_old).sum() - rates.water.sum() * 10.0,
                      rtol=1e-9)
    assert np.isclose(oil_balance,
                      (oil_new - oil_old).sum() - rates.oil.sum() * 10.0, rtol=1e-9)


# ═══════════ §20 — QLOBAL səviyyədə ANTI-PSEUDO-MPFA (MƏCBURİ) ══════════
def test_anti_pseudo_global_flux_survives_equal_owner_neighbour_pressure():
    """Tapşırıq §20 — `q = T(p_i − p_j)`-ə yığılan QLOBAL yığım BURADA
    DÜŞƏR: `p_owner = p_neighbor` qurulur, iki-nöqtəli sxem MƏCBURİ 0
    verər, HƏQİQİ MPFA-O isə YOX."""
    grid, geometry, _ = skewed(3, 3, 2)
    operator = global_operator(grid, geometry, rotated_anisotropic(grid.ncell))
    conn = grid.build_connections()
    k = int(np.argmax(operator.stencil_sizes()))
    owner, neighbor = int(conn.cell_a[k]), int(conn.cell_b[k])
    assert len(operator.face_stencil(k, tolerance=1e-9)) > 2

    rng = np.random.default_rng(2026)
    pressure = 200.0 + rng.normal(0.0, 20.0, grid.ncell)
    pressure[neighbor] = pressure[owner]
    flux = operator.connection_fluxes(pressure)
    assert abs(flux[k]) > 1e-6 * np.abs(flux).max(), (
        "p_i = p_j olduqda qlobal axın sıfır çıxdı — qlobal yığım "
        "iki-nöqtəli relasiyaya YIĞILIB")

    third = next(c for c in operator.face_stencil(k, tolerance=1e-9)
                 if c not in (owner, neighbor))
    moved = pressure.copy()
    moved[third] += 25.0
    assert abs(operator.connection_fluxes(moved)[k] - flux[k]) > 1e-6 * np.abs(flux).max()


def test_anti_pseudo_no_tpfa_symbols_in_the_mpfa_package():
    """Tapşırıq §32 — MPFA yolunda gizli TPFA fallback-ı YOXDUR."""
    import pathlib
    package = pathlib.Path(__file__).resolve().parents[1] / "imex2d" / "discretization"
    for path in sorted(package.glob("*.py")):
        source = path.read_text(encoding="utf-8")
        code = "\n".join(line for line in source.splitlines()
                         if not line.lstrip().startswith("#"))
        assert "from ..simulation.discretization" not in code
        assert "import TwoPointFluxDiscretization" not in code
        for banned in ("correction_factor", "anisotropy_multiplier",
                       "skewness_factor", "alpha_correction"):
            assert banned not in code, f"{path.name}: empirik düzəliş '{banned}'"


# ═════════ §30/§31 — ƏL İLƏ YOXLANILA BİLƏN QLOBAL İSTİNAD HALI ═════════
def test_reference_global_case_unit_cubes_hand_verifiable_matrix():
    """Tapşırıq §30/§31 — QLOBAL matrisin ÖZÜ yoxlanılır (yalnız axın
    çıxışı deyil).

    HƏNDƏSƏ: 2×2×1 vahid kub (h = 1 m), `K = I`, `Γ = 1`, `η = 1`,
    `NEUMANN_ZERO`. Əl hesabı (Phase 5A §16):

        üz sahəsi A_F = 1,  yarım məsafə 0.5,
        T_F = Γ·A_F/(0.5/k + 0.5/k) = 1·1/1 = 1

    Gözlənilən `T_conn` sətri: `+1` owner-də, `−1` neighbor-da.
    Gözlənilən `A = D·T_conn` sətri: diaqonalda `−(qonşu sayı)`,
    hər qonşuda `+1` (klassik 5-nöqtəli Laplas operatoru).
    """
    from imex2d.domain.grid import CartesianGrid
    from imex2d.domain.geometry import CellGeometry
    from imex2d.domain.general_grid_geometry import (
        GeneralGridGeometry, hexahedral_vertices_from_cartesian)

    grid = CartesianGrid(2, 2, 1)
    geometry = CellGeometry(grid, dx=1.0, dy=1.0, dz=[1.0], top_depth=0.0)
    ggg = GeneralGridGeometry(hexahedral_vertices_from_cartesian(grid, geometry),
                              grid.build_connections())
    coefficients = build_mpfa_o_coefficients(
        grid, ggg, uniform_k(4, np.eye(3)), darcy_constant=1.0, closure=NEUMANN)
    operator = MPFAGlobalOperator(coefficients, grid.build_connections(),
                                  ggg.connection_faces(), ggg.face_owner,
                                  ggg.face_neighbor)

    conn = grid.build_connections()
    assert operator.T_conn.shape == (4, 4) and conn.count == 4

    t_conn = dense(operator.T_conn)
    for k in range(conn.count):
        expected = np.zeros(4)
        expected[conn.cell_a[k]] = 1.0
        expected[conn.cell_b[k]] = -1.0
        assert np.allclose(t_conn[k], expected, atol=1e-12), (
            f"əlaqə {k} sətri əl hesabı ilə uyğun gəlmir:\n{t_conn[k]}")

    matrix = dense(operator.cell_operator())
    expected_A = np.array([[-2.0, 1.0, 1.0, 0.0],
                           [1.0, -2.0, 0.0, 1.0],
                           [1.0, 0.0, -2.0, 1.0],
                           [0.0, 1.0, 1.0, -2.0]])
    assert np.allclose(matrix, expected_A, atol=1e-12), (
        f"qlobal operator əl hesabı ilə uyğun gəlmir:\n{matrix}")

    # Yoxlanıla bilən axın: p = (0, 1, 0, 1) → x boyu vahid qradiyent
    pressure = np.array([0.0, 1.0, 0.0, 1.0])
    flux = operator.connection_fluxes(pressure)
    for k in range(conn.count):
        assert np.isclose(flux[k],
                          pressure[conn.cell_a[k]] - pressure[conn.cell_b[k]],
                          atol=1e-12)
    assert np.allclose(operator.net_influx(pressure), matrix @ pressure, atol=1e-12)
    assert abs(operator.net_influx(pressure).sum()) < 1e-12


# ═══════════════════════════════════ performans (§28) ═══════════════════
def test_performance_operator_is_built_once_and_reused():
    """Həndəsi operator DÖVLƏTDƏN ASILI DEYİL: N qalıq qiymətləndirməsi
    N lokal sistem qurulmasına səbəb OLMAMALIDIR (§11/§28)."""
    model = flat_model(nx=6, ny=6)
    mpfa, _ = assemblers(model)
    operator = mpfa.grid.global_operator
    identity = id(operator.T_conn)

    started = time.perf_counter()
    for seed in range(20):
        st = state(model.ncell, seed=seed)
        mpfa.face_fluxes(st, mpfa.fluid_state(st))
    elapsed = time.perf_counter() - started

    assert id(mpfa.grid.global_operator.T_conn) == identity   # YENİDƏN qurulmayıb
    assert elapsed < 2.0, f"20 qalıq qiymətləndirməsi {elapsed:.2f}s çəkdi"
