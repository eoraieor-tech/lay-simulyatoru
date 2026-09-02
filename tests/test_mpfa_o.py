"""Phase 5A — TRUE MPFA-O riyazi nüvəsinin test dəsti.

Tapşırıq §28-də tələb olunan A–R kateqoriyaları + §29 (anti-pseudo-MPFA)
+ §30 (əl ilə yoxlanıla bilən istinad halı) + §31 (performans).

Bax `docs/mpfa_o_phase5a.md` — hər testin yoxladığı düstur orada.

TOLERANTLIQ SİYASƏTİ (tapşırıq §32): heç bir tolerantlıq "testi keçirmək
üçün" boşaldılmayıb. Xətti sahə testləri MAŞIN DƏQİQLİYİ səviyyəsindədir
(`rtol=1e-11`) — MPFA-O xətti təzyiq sahələrini DƏQİQ bərpa edir, ona
görə bundan zəif tolerantlıq gizli səhvi ört-basdır edərdi.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from helpers import five_spot_model
from imex2d.discretization import (MPFAOBoundaryClosure, MPFAOCoefficients,
                                   MPFAODiscretization, MPFAOLocalSystem,
                                   MPFAOSingularSystemError, MPFAOTensorError,
                                   build_interaction_regions,
                                   build_mpfa_o_coefficients,
                                   validate_permeability_matrices)
from imex2d.discretization.mpfa_o import permeability_matrices
from imex2d.discretization.mpfa_o_interaction import validate_interaction_regions
from imex2d.domain.general_grid_geometry import (GeneralGridGeometry,
                                                 hexahedral_vertices_from_cartesian)
from imex2d.domain.geometry import CellGeometry
from imex2d.domain.grid import CartesianGrid
from imex2d.domain.properties import PermeabilityTensor, PropertyMap
from imex2d.domain.units import METRIC
from imex2d.interfaces.discretization import IFluxDiscretization
from imex2d.simulation.discretization import (TwoPointFluxDiscretization,
                                              default_flux_discretization)

GAMMA = METRIC.darcy_constant


# ══════════════════════════════════════════════════ köməkçilər ══════════
def cartesian(nx=3, ny=3, nz=2, dx=10.0, dy=10.0, dz=5.0, top=0.0):
    grid = CartesianGrid(nx, ny, nz)
    geometry = CellGeometry(grid, dx=dx, dy=dy, dz=[dz] * nz, top_depth=top)
    vertices = hexahedral_vertices_from_cartesian(grid, geometry)
    return grid, GeneralGridGeometry(vertices, grid.build_connections()), vertices


def skewed(nx=3, ny=3, nz=2, shear=0.35, **kwargs):
    """QEYRİ-ORTOQONAL grid: x ← x + s·y + 0.2s·z, y ← y + 0.15s·z.

    Affin çevirmə olduğu üçün hüceyrələr paralelepiped qalır (etibarlı
    həndəsə), AMMA d_ij ARTIQ n_f-ə PARALEL DEYİL — TPFA-nın fərziyyəsi
    pozulur, MPFA-O isə dəqiq qalmalıdır."""
    grid, _, vertices = cartesian(nx, ny, nz, **kwargs)
    skew = vertices.copy()
    skew[:, :, 0] = vertices[:, :, 0] + shear * vertices[:, :, 1] \
        + 0.2 * shear * vertices[:, :, 2]
    skew[:, :, 1] = vertices[:, :, 1] + 0.15 * shear * vertices[:, :, 2]
    return grid, GeneralGridGeometry(skew, grid.build_connections()), skew


def warped(nx=3, ny=2, nz=2, amplitude=0.9, seed=0, **kwargs):
    """UYĞUN (conforming) qeyri-ortoqonal grid: koordinatlara HAMAR,
    yalnız MÖVQEDƏN asılı bir çevirmə tətbiq olunur.

    Təpələri müstəqil "titrətmək" (per-cell noise) grid-i QEYRİ-UYĞUN
    edərdi (paylaşılan təpə iki hüceyrədə fərqli yerə düşər) — MPFA-O
    bunu `validate_interaction_regions` ilə RƏDD EDİR, ona görə testlər
    də HƏQİQİ, uyğun bir deformasiyadan istifadə etməlidir."""
    grid, _, vertices = cartesian(nx, ny, nz, **kwargs)
    rng = np.random.default_rng(seed)
    phase = rng.uniform(0.0, 2.0 * np.pi, 9).reshape(3, 3)
    frequency = rng.uniform(0.02, 0.09, 9).reshape(3, 3)
    moved = vertices.copy()
    for axis in range(3):
        shift = np.zeros(vertices.shape[:2])
        for source in range(3):
            shift += np.sin(frequency[axis, source] * vertices[:, :, source]
                            + phase[axis, source])
        moved[:, :, axis] = vertices[:, :, axis] + amplitude * shift
    return grid, GeneralGridGeometry(moved, grid.build_connections()), moved


def rotation(axis, angle_deg):
    axis = np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis)
    cross = np.array([[0.0, -axis[2], axis[1]],
                      [axis[2], 0.0, -axis[0]],
                      [-axis[1], axis[0], 0.0]])
    theta = np.deg2rad(angle_deg)
    return np.eye(3) + np.sin(theta) * cross + (1.0 - np.cos(theta)) * cross @ cross


def uniform_k(ncell, matrix):
    return np.tile(np.asarray(matrix, float), (ncell, 1, 1))


def isotropic(ncell, k=100.0):
    return uniform_k(ncell, k * np.eye(3))


def rotated_anisotropic(ncell, eigenvalues=(500.0, 50.0, 10.0),
                        axis=(0.3, -0.5, 0.8), angle=37.0):
    """K = R diag(k1,k2,k3) Rᵀ — sıfırdan FƏRQLİ Kxy/Kxz/Kyz ilə."""
    rot = rotation(axis, angle)
    return uniform_k(ncell, rot @ np.diag(eigenvalues) @ rot.T)


def analytic_flux(geometry, k_matrix, gradient, gamma=GAMMA):
    """`q_F = Γ (A_F n_F) · (−K ∇p)` — ANALİTİK Darsi axını (§14/§16)."""
    velocity = -gamma * (np.asarray(k_matrix, float) @ np.asarray(gradient, float))
    return (geometry.face_areas[:, None] * geometry.face_normals) @ velocity


def linear_field(geometry, gradient, offset=0.0):
    gradient = np.asarray(gradient, float)
    pressures = geometry.cell_centroids @ gradient + offset
    return pressures, (lambda x: float(x @ gradient + offset))


def coefficients(grid, geometry, k_matrices, eta=1.0,
                 closure=MPFAOBoundaryClosure.DIRICHLET):
    return build_mpfa_o_coefficients(grid, geometry, k_matrices, GAMMA,
                                     eta=eta, closure=closure)


def dense(matrix):
    """`MPFAOCoefficients.T_cell/T_bnd` SEYRƏK matrisdir (bax §31 — sıx
    saxlama `O(N²)` yaddaş tutardı); testlərdə müqayisə üçün sıxlaşdırılır."""
    return np.asarray(matrix.toarray())


def cell_divergence(geometry, flux):
    """Hər hüceyrə üzrə net çıxan axın — konservasiya yoxlaması."""
    div = np.zeros(geometry.ncell)
    np.add.at(div, geometry.face_owner, flux)
    interior = geometry.face_neighbor >= 0
    np.add.at(div, geometry.face_neighbor[interior], -flux[interior])
    return div


# ═══════════════════════════════════════ A — API / interfeys uyğunluğu ══
def test_a_mpfa_implements_flux_discretization_interface():
    scheme = MPFAODiscretization()
    assert isinstance(scheme, IFluxDiscretization)
    assert scheme.supports_multipoint_stencil() is True
    # TPFA — DƏYİŞMƏMİŞ defolt davranış (geriyə uyğunluq, tapşırıq §23)
    assert TwoPointFluxDiscretization().supports_multipoint_stencil() is False
    assert isinstance(default_flux_discretization(), TwoPointFluxDiscretization)


def test_a_build_from_reservoir_model_produces_coefficients():
    model = five_spot_model(nx=3, ny=3, dx=10.0, dy=10.0, dz=5.0, permeability=200.0)
    discretized = MPFAODiscretization().build(model)

    assert discretized.connections.count == model.connections().count
    assert discretized.pore_volume.shape == (model.ncell,)
    assert discretized.cell_volume.shape == (model.ncell,)
    assert discretized.coefficients.n_cell == model.ncell
    assert discretized.coefficients.n_face == len(discretized.geometry.faces)
    assert discretized.connection_faces.shape == (model.connections().count,)


def test_a_compute_flux_from_delta_phi_is_explicitly_not_implemented():
    """Tapşırıq §1/§22: saxta TPFA uyğunluğu YARADILMIR — AÇIQ xəta."""
    model = five_spot_model(nx=2, ny=2, dx=10.0, dy=10.0, dz=5.0)
    discretized = MPFAODiscretization().build(model)
    with pytest.raises(NotImplementedError, match="Phase 5B"):
        discretized.compute_flux(np.ones(discretized.connections.count))


def test_a_eta_outside_unit_interval_is_rejected():
    grid, geometry, _ = cartesian(2, 2, 1)
    for bad in (0.0, -0.5, 1.5):
        with pytest.raises(ValueError, match="η"):
            build_interaction_regions(grid, geometry, eta=bad)


# ══════════════════════════════════════════ B — sabit təzyiq → sıfır axın
@pytest.mark.parametrize("name", ["ortoqonal", "skew", "anizotrop", "fırlanmış"])
def test_b_constant_pressure_gives_zero_flux(name):
    """Tapşırıq §17 — p=C ⇒ u=0, DÖRD konfiqurasiyanın HAMISINDA."""
    if name in ("ortoqonal", "anizotrop"):
        grid, geometry, _ = cartesian()
    else:
        grid, geometry, _ = skewed()
    ncell = grid.ncell
    k = {"ortoqonal": isotropic(ncell),
         "skew": isotropic(ncell),
         "anizotrop": uniform_k(ncell, np.diag([500.0, 50.0, 10.0])),
         "fırlanmış": rotated_anisotropic(ncell)}[name]

    coef = coefficients(grid, geometry, k)
    pressures = np.full(ncell, 137.5)
    boundary = coef.boundary_pressures_from(lambda x: 137.5)
    flux = coef.face_fluxes(pressures, boundary)

    assert np.max(np.abs(flux)) < 1e-9 * 137.5


# ═════════════════════════════════════ C — xətti təzyiq → analitik Darsi
@pytest.mark.parametrize("builder,k_name", [
    (cartesian, "izotrop"), (cartesian, "fırlanmış"),
    (skewed, "izotrop"), (skewed, "anizotrop"), (skewed, "fırlanmış"),
])
def test_c_linear_pressure_reproduces_analytic_darcy_flux(builder, k_name):
    """Tapşırıq §14/§16 — TPFA ilə DEYİL, ANALİTİK axınla müqayisə."""
    grid, geometry, _ = builder()
    k_matrix = {"izotrop": 100.0 * np.eye(3),
                "anizotrop": np.diag([500.0, 50.0, 10.0]),
                "fırlanmış": rotated_anisotropic(1)[0]}[k_name]
    k = uniform_k(grid.ncell, k_matrix)
    gradient = np.array([2.0, -3.0, 1.25])

    coef = coefficients(grid, geometry, k)
    pressures, exact = linear_field(geometry, gradient, offset=11.0)
    flux = coef.face_fluxes(pressures, coef.boundary_pressures_from(exact))
    expected = analytic_flux(geometry, k_matrix, gradient)

    assert np.allclose(flux, expected, rtol=1e-11, atol=1e-11 * np.abs(expected).max())


def test_c_linear_pressure_is_locally_divergence_free():
    """Xətti sahə + sabit K ⇒ hər hüceyrədə Σ(üz axınları) = 0."""
    grid, geometry, _ = skewed()
    k = rotated_anisotropic(grid.ncell)
    coef = coefficients(grid, geometry, k)
    pressures, exact = linear_field(geometry, [1.5, 0.75, -2.0])
    flux = coef.face_fluxes(pressures, coef.boundary_pressures_from(exact))

    divergence = cell_divergence(geometry, flux)
    assert np.max(np.abs(divergence)) < 1e-9 * np.abs(flux).max()


# ═══════════════════════════ D — izotrop ortoqonal limit: MPFA ≡ TPFA ══
def test_d_orthogonal_isotropic_limit_matches_tpfa_transmissibility():
    """Tapşırıq §13. MPFA DAXİLƏN TPFA ÇAĞIRMIR — bu, formulyasiyanın
    öz limitidir (bax `docs/mpfa_o_phase5a.md` §16)."""
    model = five_spot_model(nx=4, ny=3, dx=10.0, dy=12.0, dz=5.0, permeability=200.0)
    tpfa = TwoPointFluxDiscretization().build(model)
    mpfa = MPFAODiscretization().build(model)
    conn = model.connections()

    rng = np.random.default_rng(20260902)
    pressures = 200.0 + rng.normal(0.0, 5.0, model.ncell)
    # Sərhəd π-ləri: MPFA-nın DAXİLİ üz axınları sərhəd dəyərlərindən
    # ASILI OLMAMALIDIR (izotrop ortoqonal limitdə stensil 2 hüceyrəlidir),
    # ona görə ixtiyari (sıfır) sərhəd dəyəri verilir.
    boundary = np.zeros(len(mpfa.coefficients.boundary_dofs))
    mpfa_flux = mpfa.connection_fluxes(pressures, boundary)
    tpfa_flux = tpfa.transmissibility * (pressures[conn.cell_a] - pressures[conn.cell_b])

    assert np.allclose(mpfa_flux, tpfa_flux, rtol=1e-11,
                       atol=1e-11 * np.abs(tpfa_flux).max())


def test_d_orthogonal_isotropic_stencil_collapses_to_two_points():
    """§16-nın STENSİL səviyyəsində təsdiqi: K-ortoqonal Kartezian
    gridd-də MPFA-O əmsalları DƏQİQ iki hüceyrəlidir."""
    grid, geometry, _ = cartesian(4, 3, 2)
    coef = coefficients(grid, geometry, isotropic(grid.ncell, 200.0))
    interior = np.flatnonzero(~geometry.is_boundary)

    assert np.all(coef.stencil_sizes()[interior] == 2)
    for face in interior:
        stencil = coef.face_stencil(face, tolerance=1e-9)
        owner, neighbor = geometry.face_owner[face], geometry.face_neighbor[face]
        assert set(stencil) == {int(owner), int(neighbor)}
        assert stencil[int(owner)] > 0.0 and stencil[int(neighbor)] < 0.0
        assert np.isclose(stencil[int(owner)], -stencil[int(neighbor)], rtol=1e-12)


def test_d_boundary_faces_do_not_pollute_interior_transmissibility():
    """İzotrop ortoqonal limitdə DAXİLİ üzlərin `T_bnd` sətri sıfırdır —
    yəni daxili axın sərhəd dəyərlərindən ASILI DEYİL."""
    grid, geometry, _ = cartesian(3, 3, 2)
    coef = coefficients(grid, geometry, isotropic(grid.ncell))
    interior = np.flatnonzero(~geometry.is_boundary)
    assert np.max(np.abs(dense(coef.T_bnd)[interior])) < 1e-9


# ═════════════════════════════ E — skew həndəsə: etibarlı çoxnöqtəlilik
def test_e_skew_geometry_regions_are_valid_and_well_conditioned():
    grid, geometry, _ = skewed()
    regions = build_interaction_regions(grid, geometry)
    assert validate_interaction_regions(regions) == []
    assert len(regions) == (grid.nx + 1) * (grid.ny + 1) * (grid.nz + 1)

    interior = [r for r in regions if not r.is_boundary_region]
    assert interior, "daxili bölgə tapılmadı"
    for region in interior:
        assert region.n_cells == 8              # §2: daxili node → 8 hüceyrə
        assert region.n_sub_faces == 12         # §2: 8·3/2 = 12 sub-üz
        assert all(len(s) == 3 for s in region.cell_sub_faces)

    coef = coefficients(grid, geometry, rotated_anisotropic(grid.ncell))
    assert np.all(np.isfinite(coef.condition_numbers()))
    assert coef.ill_conditioned_regions() == []


def test_e_skew_geometry_is_genuinely_non_orthogonal():
    """Test ssenarisinin ÖZÜNÜN etibarlılığı: grid həqiqətən skewdirmi?"""
    _, geometry, _ = skewed()
    assert geometry.quality_metrics()["max_non_orthogonality_angle_deg"] > 10.0
    assert geometry.validate().ok


# ══════════════════════ F — diaqonal anizotropluq: tenzorun düzgün təsiri
def test_f_diagonal_anisotropy_scales_each_direction_independently():
    """`K = diag(kx,ky,kz)` — hər ox üzrə axın MƏHZ öz kx/ky/kz-i ilə."""
    grid, geometry, _ = cartesian(3, 3, 2)
    kx, ky, kz = 400.0, 40.0, 4.0
    coef = coefficients(grid, geometry, uniform_k(grid.ncell, np.diag([kx, ky, kz])))

    for axis, k_axis in enumerate((kx, ky, kz)):
        gradient = np.zeros(3)
        gradient[axis] = 1.0
        pressures, exact = linear_field(geometry, gradient)
        flux = coef.face_fluxes(pressures, coef.boundary_pressures_from(exact))
        expected = analytic_flux(geometry, np.diag([kx, ky, kz]), gradient)
        assert np.allclose(flux, expected, rtol=1e-11,
                           atol=1e-11 * np.abs(expected).max())
        # Yalnız bu oxun keçiriciliyi işə düşür
        normal_faces = np.flatnonzero(np.abs(geometry.face_normals[:, axis]) > 0.99)
        assert np.abs(flux[normal_faces]).max() > 0.0


# ══════════════════ G — fırlanmış anizotropluq: off-diaqonal HƏQİQƏTƏN işlənir
def test_g_off_diagonal_components_change_the_coefficients():
    """Tapşırıq §11/§E — Kxy/Kxz/Kyz DƏYİŞƏNDƏ əmsallar DƏYİŞMƏLİDİR.

    Əks halda implementasiya tenzoru səssizcə diaqonala yığır."""
    grid, geometry, _ = skewed()
    diagonal = np.diag([500.0, 50.0, 10.0])
    base = coefficients(grid, geometry, uniform_k(grid.ncell, diagonal))

    for label, (i, j) in (("Kxy", (0, 1)), ("Kxz", (0, 2)), ("Kyz", (1, 2))):
        perturbed_matrix = diagonal.copy()
        perturbed_matrix[i, j] = perturbed_matrix[j, i] = 15.0
        assert np.linalg.eigvalsh(perturbed_matrix)[0] > 0.0, f"{label} testi SPD deyil"
        perturbed = coefficients(grid, geometry, uniform_k(grid.ncell, perturbed_matrix))
        difference = np.max(np.abs(dense(perturbed.T_cell) - dense(base.T_cell)))
        scale = np.max(np.abs(dense(base.T_cell)))
        assert difference > 1e-6 * scale, (
            f"{label} əmsalları DƏYİŞMƏDİ — off-diaqonal komponent İSTİFADƏ "
            "OLUNMUR (tapşırıq §11 pozuntusu).")


def test_g_rotated_tensor_flux_follows_the_rotation():
    """Fırlanmış tenzor + xətti sahə → analitik `−K∇p` ilə üst-üstə."""
    grid, geometry, _ = skewed()
    for angle in (0.0, 15.0, 37.0, 72.0, 115.0):
        k_matrix = rotated_anisotropic(1, angle=angle)[0]
        assert np.max(np.abs(k_matrix - np.diag(np.diag(k_matrix)))) > 1.0 or angle == 0.0
        coef = coefficients(grid, geometry, uniform_k(grid.ncell, k_matrix))
        gradient = np.array([1.0, 2.0, -0.5])
        pressures, exact = linear_field(geometry, gradient)
        flux = coef.face_fluxes(pressures, coef.boundary_pressures_from(exact))
        expected = analytic_flux(geometry, k_matrix, gradient)
        assert np.allclose(flux, expected, rtol=1e-11,
                           atol=1e-11 * np.abs(expected).max()), f"bucaq {angle}"


def test_g_permeability_tensor_from_model_uses_all_six_components():
    """`PermeabilityTensor` → MPFA girişi: 6 komponentin HAMISI keçir."""
    model = five_spot_model(nx=2, ny=2, dx=10.0, dy=10.0, dz=5.0)
    n = model.ncell
    model.rock.permeability_tensor = PermeabilityTensor(
        kxx=PropertyMap.uniform("KXX", 300.0, n),
        kyy=PropertyMap.uniform("KYY", 200.0, n),
        kzz=PropertyMap.uniform("KZZ", 100.0, n),
        kxy=PropertyMap.uniform("KXY", 25.0, n),
        kxz=PropertyMap.uniform("KXZ", -15.0, n),
        kyz=PropertyMap.uniform("KYZ", 10.0, n))

    matrices = permeability_matrices(model)
    assert np.allclose(matrices[0], [[300.0, 25.0, -15.0],
                                     [25.0, 200.0, 10.0],
                                     [-15.0, 10.0, 100.0]])
    discretized = MPFAODiscretization().build(model)
    assert discretized.coefficients.n_cell == n


# ═════════════════════════════════ H — çoxnöqtəli stensil (>2 hüceyrə) ══
def test_h_non_orthogonal_full_tensor_produces_multi_point_stencil():
    """Tapşırıq §10 — |S_f| > 2 OLMALIDIR."""
    grid, geometry, _ = skewed(3, 3, 2)
    coef = coefficients(grid, geometry, rotated_anisotropic(grid.ncell))
    interior = np.flatnonzero(~geometry.is_boundary)
    sizes = coef.stencil_sizes()[interior]

    assert sizes.min() > 2, f"ən kiçik stensil {sizes.min()} — çoxnöqtəli DEYİL"
    assert sizes.max() >= 8


def test_h_interior_region_couples_up_to_eight_cells_per_sub_face():
    grid, geometry, _ = skewed(3, 3, 2)
    coef = coefficients(grid, geometry, rotated_anisotropic(grid.ncell))
    interior_regions = [s for s in coef.local_systems
                        if not s.region.is_boundary_region]
    assert interior_regions
    system = interior_regions[0]
    assert system.T_cell.shape == (12, 8)      # 12 sub-üz × 8 hüceyrə
    # Hər sub-üz sətri 2-dən çox hüceyrəyə toxunur
    touched = np.sum(np.abs(system.T_cell) > 1e-9 * np.abs(system.T_cell).max(), axis=1)
    assert touched.min() > 2


# ══════════════════════════════════════════ I — lokal konservasiya ══════
@pytest.mark.parametrize("builder", [cartesian, skewed])
def test_i_local_conservation_of_sub_face_fluxes(builder):
    """Tapşırıq §9/§15 — hər tərəf ÖZ qradiyentindən müstəqil hesablanır."""
    grid, geometry, _ = builder()
    coef = coefficients(grid, geometry, rotated_anisotropic(grid.ncell))
    rng = np.random.default_rng(5)
    pressures = 100.0 + rng.normal(0.0, 10.0, grid.ncell)
    boundary = rng.normal(100.0, 10.0, len(coef.boundary_dofs))

    residual = coef.conservation_residual(pressures, boundary)
    flux_scale = np.abs(coef.face_fluxes(pressures, boundary)).max()
    assert residual < 1e-9 * flux_scale, f"maks. konservasiya qalığı {residual:.3e}"


def test_i_face_flux_is_antisymmetric_between_owner_and_neighbour():
    """`q_{i,F} + q_{j,F} = 0` — İDENTİK, konstruksiyaca (§9)."""
    grid, geometry, _ = skewed()
    coef = coefficients(grid, geometry, rotated_anisotropic(grid.ncell))
    rng = np.random.default_rng(7)
    pressures = rng.normal(150.0, 8.0, grid.ncell)
    boundary = rng.normal(150.0, 8.0, len(coef.boundary_dofs))
    flux = coef.face_fluxes(pressures, boundary)

    for system in coef.local_systems:
        halves = system.half_fluxes(pressures[system.region.cells],
                                    _local_boundary(coef, system, boundary))
        for sub in system.region.sub_faces:
            if sub.is_boundary:
                continue
            pair = halves[(sub.owner, sub.local_index)] + halves[(sub.neighbor,
                                                                 sub.local_index)]
            assert abs(pair) < 1e-9 * max(1.0, np.abs(flux).max())


def _local_boundary(coef, system, boundary):
    if not system.known_boundary_sub_faces:
        return None
    return boundary[[coef._boundary_dof[(system.region.sub_faces[s].face_index,
                                         system.region.node_id)]
                     for s in system.known_boundary_sub_faces]]


# ═════════════════════════════════════════ J — sərhəd bölgələri ═════════
def test_j_boundary_regions_are_explicitly_represented():
    """Tapşırıq §20 — sərhəd bölgəsi AÇIQ fərqləndirilir, "hər bölgənin
    8 hüceyrəsi var" fərz EDİLMİR."""
    grid, geometry, _ = cartesian(3, 3, 2)
    regions = build_interaction_regions(grid, geometry)
    boundary = [r for r in regions if r.is_boundary_region]
    interior = [r for r in regions if not r.is_boundary_region]

    assert boundary and interior
    assert min(r.n_cells for r in boundary) == 1       # künc node → 1 hüceyrə
    assert all(r.n_cells == 8 for r in interior)
    assert all(len(r.boundary_sub_faces) > 0 for r in boundary)
    assert all(len(r.boundary_sub_faces) == 0 for r in interior)


def test_j_dirichlet_closure_exposes_boundary_coefficients():
    grid, geometry, _ = cartesian(2, 2, 1)
    coef = coefficients(grid, geometry, isotropic(grid.ncell))
    assert len(coef.boundary_dofs) == coef.T_bnd.shape[1] > 0
    assert coef.boundary_points.shape == (len(coef.boundary_dofs), 3)
    with pytest.raises(ValueError, match="boundary_pressures"):
        coef.face_fluxes(np.zeros(grid.ncell))


def test_j_neumann_zero_closure_gives_no_flow_boundaries():
    """`NEUMANN_ZERO` — AÇIQ seçilən FİZİKİ şərt (§10); sərhəd üzlərində
    axın DƏQİQ sıfır olmalıdır."""
    grid, geometry, _ = cartesian(3, 3, 2)
    coef = coefficients(grid, geometry, isotropic(grid.ncell),
                        closure=MPFAOBoundaryClosure.NEUMANN_ZERO)
    assert coef.T_bnd.shape[1] == 0
    rng = np.random.default_rng(11)
    pressures = 200.0 + rng.normal(0.0, 5.0, grid.ncell)
    flux = coef.face_fluxes(pressures)

    boundary_faces = np.flatnonzero(geometry.is_boundary)
    assert np.max(np.abs(flux[boundary_faces])) < 1e-9 * np.abs(flux).max()
    assert np.abs(flux[~geometry.is_boundary]).max() > 0.0


def test_j_single_layer_grid_is_a_confined_3d_layer_not_a_fake_2d_method():
    """Tapşırıq §26 — nz=1 AYRI 2D alqoritmi ilə DEYİL, məhdud 3D lay
    kimi həll olunur: üst/alt üzlər sərhəddir və bağlanış alır."""
    grid, geometry, _ = cartesian(3, 3, 1)
    coef = coefficients(grid, geometry, uniform_k(grid.ncell, np.diag([200.0, 200.0, 5.0])),
                        closure=MPFAOBoundaryClosure.NEUMANN_ZERO)
    vertical = np.flatnonzero(np.abs(geometry.face_normals[:, 2]) > 0.99)
    assert np.all(geometry.is_boundary[vertical])

    rng = np.random.default_rng(3)
    pressures = 180.0 + rng.normal(0.0, 4.0, grid.ncell)
    flux = coef.face_fluxes(pressures)
    assert np.max(np.abs(flux[vertical])) < 1e-9 * np.abs(flux).max()
    # Lateral axın YENƏ DƏ 3D formulyasiyadan gəlir və sıfırdan fərqlidir
    assert np.abs(flux[~geometry.is_boundary]).max() > 0.0


# ═══════════════════════════════ K — sinqulyar lokal sistem: AÇIQ xəta ══
def test_k_degenerate_flattened_cell_is_rejected_explicitly():
    """Tapşırıq §19/§32 — degenerativ həndəsə AÇIQ rədd edilir; səssiz
    requlyarizasiya və ya TPFA-ya keçid YOXDUR.

    Yastılanmış hüceyrədə sub-üzlərin sahəsi sıfırdır, ona görə
    `validate_interaction_regions` LOKAL SİSTEMƏ ÇATMAMIŞ dayandırır —
    bu, daha erkən və daha aydın diaqnostikadır.
    """
    grid, _, vertices = cartesian(2, 2, 1)
    flattened = vertices.copy()
    flattened[0, 4:, 2] = flattened[0, :4, 2]      # 0-cı hüceyrəni yastılaşdır
    geometry = GeneralGridGeometry(flattened, grid.build_connections())

    with pytest.raises(ValueError, match="degenerativ"):
        coefficients(grid, geometry, isotropic(grid.ncell))


def test_k_singular_local_system_raises_with_diagnostics():
    """`MPFAOSingularSystemError` yolunun BİRBAŞA yoxlanması.

    Etibarlı hekzahedral həndəsədə `D_(c,v)` praktiki olaraq heç vaxt
    sinqulyar olmur (paralelepiped üçün `det D ∝ həcm`), ona görə bu
    test MÜHAFİZƏNİ ÖZÜNÜ yoxlayır: 0-cı hüceyrənin mərkəzi QƏSDƏN öz
    3 kəsilməzlik nöqtəsinin MÜSTƏVİSİNƏ qoyulur (vahid kub bölgəsində
    həmin müstəvi `x+y+z = 2`-dir) → `D` sətirləri asılı olur.

    Gözlənilən: AÇIQ xəta + diaqnostika; SƏSSİZ ε-requlyarizasiya YOX.
    """
    grid = CartesianGrid(2, 2, 2)
    geometry = CellGeometry(grid, dx=1.0, dy=1.0, dz=[1.0, 1.0], top_depth=0.0)
    ggg = GeneralGridGeometry(hexahedral_vertices_from_cartesian(grid, geometry),
                              grid.build_connections())
    region = next(r for r in build_interaction_regions(grid, ggg)
                  if not r.is_boundary_region)

    centroids = ggg.cell_centroids.copy()
    centroids[0] = np.array([2.0, 2.0, 2.0]) / 3.0      # x+y+z = 2 müstəvisində
    with pytest.raises(MPFAOSingularSystemError) as error:
        MPFAOLocalSystem(region=region, cell_centroids=centroids,
                         k_matrices=uniform_k(grid.ncell, np.eye(3)),
                         darcy_constant=1.0)
    assert "sinqulyar" in str(error.value)
    assert error.value.diagnostics.get("region_id") == region.node_id
    assert error.value.diagnostics.get("cell") == 0


def test_k_diagnostics_expose_condition_rank_and_determinant():
    grid, geometry, _ = skewed()
    coef = coefficients(grid, geometry, rotated_anisotropic(grid.ncell))
    sample = coef.diagnostics[len(coef.diagnostics) // 2].as_dict()
    for key in ("region_id", "condition_number", "rank", "determinant",
                "singular", "ill_conditioned", "n_unknowns", "is_boundary_region",
                "max_sub_cell_condition"):
        assert key in sample
    assert not any(d.singular for d in coef.diagnostics)


def test_k_extremely_ill_conditioned_region_is_flagged_not_hidden():
    """Şərt ədədi həddi aşanda AÇIQ bayraq qalxır (hesablama davam edir)."""
    grid, geometry, _ = skewed(shear=0.9)
    coef = build_mpfa_o_coefficients(
        grid, geometry, rotated_anisotropic(grid.ncell, eigenvalues=(1e4, 1.0, 1e-2)),
        GAMMA, condition_warning_threshold=1.0)
    assert len(coef.ill_conditioned_regions()) > 0
    assert all(not d.singular for d in coef.diagnostics)


# ═══════════════════════════════════════ L/M — etibarsız tenzorun rəddi ══
def test_l_non_positive_definite_tensor_is_rejected():
    """Tapşırıq §18 — eigenvalue klipləmə/ε əlavəsi YOXDUR."""
    grid, geometry, _ = cartesian(2, 2, 1)
    # Diaqonal MÜSBƏT, amma λ_min < 0 — səthi yoxlama bunu buraxardı
    bad = np.array([[1.0, 5.0, 0.0], [5.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    assert np.all(np.diag(bad) > 0.0) and np.linalg.eigvalsh(bad)[0] < 0.0
    with pytest.raises(MPFAOTensorError, match="müsbət-müəyyən"):
        coefficients(grid, geometry, uniform_k(grid.ncell, bad))


def test_l_singular_tensor_is_rejected():
    grid, geometry, _ = cartesian(2, 2, 1)
    with pytest.raises(MPFAOTensorError, match="müsbət-müəyyən"):
        coefficients(grid, geometry, uniform_k(grid.ncell, np.diag([100.0, 100.0, 0.0])))


def test_l_asymmetric_tensor_is_rejected_not_symmetrised():
    grid, geometry, _ = cartesian(2, 2, 1)
    asymmetric = np.array([[100.0, 20.0, 0.0], [5.0, 100.0, 0.0], [0.0, 0.0, 100.0]])
    with pytest.raises(MPFAOTensorError, match="simmetrik"):
        coefficients(grid, geometry, uniform_k(grid.ncell, asymmetric))


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_m_nan_and_inf_permeability_are_rejected(value):
    grid, geometry, _ = cartesian(2, 2, 1)
    k = isotropic(grid.ncell)
    k[1, 0, 0] = value
    with pytest.raises(MPFAOTensorError, match="NaN/Inf"):
        coefficients(grid, geometry, k)


def test_m_validation_matches_permeability_tensor_authority():
    """`PermeabilityTensor.validate()` AVTORİTETDİR (tapşırıq §18) —
    MPFA yoxlaması onunla EYNİ nəticəni verməlidir."""
    n = 4
    tensor = PermeabilityTensor(
        kxx=PropertyMap.uniform("KXX", 1.0, n),
        kyy=PropertyMap.uniform("KYY", 1.0, n),
        kzz=PropertyMap.uniform("KZZ", 1.0, n),
        kxy=PropertyMap.uniform("KXY", 5.0, n))
    assert not tensor.validate().ok
    with pytest.raises(MPFAOTensorError):
        validate_permeability_matrices(tensor.as_matrices())

    good = PermeabilityTensor(
        kxx=PropertyMap.uniform("KXX", 100.0, n),
        kyy=PropertyMap.uniform("KYY", 80.0, n),
        kzz=PropertyMap.uniform("KZZ", 10.0, n),
        kxy=PropertyMap.uniform("KXY", 20.0, n))
    assert good.validate().ok
    validate_permeability_matrices(good.as_matrices())


# ══════════════════════════════════ N — fırlanma invariantlığı ══════════
def test_n_rotation_invariance_of_physical_flux():
    """Tapşırıq §12 — həndəsə + K + qradiyent BİRLİKDƏ fırlananda
    fiziki axın DƏYİŞMƏMƏLİDİR (koordinat fırlanması ≠ pertürbasiya)."""
    grid, geometry, vertices = skewed()
    k_matrix = rotated_anisotropic(1)[0]
    gradient = np.array([1.5, -2.0, 0.75])

    base = coefficients(grid, geometry, uniform_k(grid.ncell, k_matrix))
    pressures, exact = linear_field(geometry, gradient)
    flux = base.face_fluxes(pressures, base.boundary_pressures_from(exact))

    rot = rotation([0.2, 0.9, -0.35], 41.0)
    rotated_geometry = GeneralGridGeometry(vertices @ rot.T, grid.build_connections())
    rotated_k = uniform_k(grid.ncell, rot @ k_matrix @ rot.T)
    rotated_gradient = rot @ gradient
    rotated_coef = coefficients(grid, rotated_geometry, rotated_k)
    rotated_pressures, rotated_exact = linear_field(rotated_geometry, rotated_gradient)
    rotated_flux = rotated_coef.face_fluxes(
        rotated_pressures, rotated_coef.boundary_pressures_from(rotated_exact))

    assert np.allclose(flux, rotated_flux, rtol=1e-10,
                       atol=1e-10 * np.abs(flux).max())
    # Ssenarinin etibarlılığı: fırlanma HƏQİQƏTƏN həndəsəni dəyişib
    assert not np.allclose(geometry.cell_centroids, rotated_geometry.cell_centroids)


def test_n_arbitrary_perturbation_is_distinguished_from_rotation():
    """Fırlanma DEYİL, təsadüfi pertürbasiya axını DƏYİŞMƏLİDİR — əks
    halda test heç nə yoxlamırdı."""
    grid, geometry, vertices = skewed()
    k = rotated_anisotropic(grid.ncell)
    gradient = np.array([1.5, -2.0, 0.75])
    base = coefficients(grid, geometry, k)
    pressures, exact = linear_field(geometry, gradient)
    flux = base.face_fluxes(pressures, base.boundary_pressures_from(exact))

    _, perturbed, _ = warped(grid.nx, grid.ny, grid.nz, amplitude=0.7, seed=99)
    perturbed_coef = coefficients(grid, perturbed, k)
    perturbed_pressures, perturbed_exact = linear_field(perturbed, gradient)
    perturbed_flux = perturbed_coef.face_fluxes(
        perturbed_pressures, perturbed_coef.boundary_pressures_from(perturbed_exact))
    assert not np.allclose(flux, perturbed_flux, rtol=1e-6)


# ═══════════════════════════════ O/P — miqyas və sürüşdürmə ════════════
@pytest.mark.parametrize("scale", [0.25, 4.0, 100.0])
def test_o_dimensional_scaling_of_flux(scale):
    """Həndəsə `s` dəfə böyüyəndə: `A ~ s²`, `∇p ~ 1/s` ⇒ `q ~ s`."""
    grid, geometry, vertices = skewed()
    k = rotated_anisotropic(grid.ncell)
    gradient = np.array([1.0, -1.5, 0.5])

    base = coefficients(grid, geometry, k)
    pressures, exact = linear_field(geometry, gradient)
    flux = base.face_fluxes(pressures, base.boundary_pressures_from(exact))

    scaled_geometry = GeneralGridGeometry(vertices * scale, grid.build_connections())
    scaled_coef = coefficients(grid, scaled_geometry, k)
    # EYNİ fiziki təzyiq sahəsi: p(x) = g·x ⇒ böyüdülmüş gridd-də g/s
    scaled_pressures, scaled_exact = linear_field(scaled_geometry, gradient / scale)
    scaled_flux = scaled_coef.face_fluxes(
        scaled_pressures, scaled_coef.boundary_pressures_from(scaled_exact))

    assert np.allclose(scaled_flux, scale * flux, rtol=1e-10,
                       atol=1e-10 * np.abs(scale * flux).max())


def test_p_translation_does_not_change_flux():
    grid, geometry, vertices = skewed()
    k = rotated_anisotropic(grid.ncell)
    gradient = np.array([2.0, 1.0, -0.5])
    shift = np.array([1234.5, -678.9, 4321.0])

    base = coefficients(grid, geometry, k)
    pressures, exact = linear_field(geometry, gradient)
    flux = base.face_fluxes(pressures, base.boundary_pressures_from(exact))

    moved = GeneralGridGeometry(vertices + shift, grid.build_connections())
    moved_coef = coefficients(grid, moved, k)
    moved_pressures, moved_exact = linear_field(moved, gradient)
    moved_flux = moved_coef.face_fluxes(
        moved_pressures, moved_coef.boundary_pressures_from(moved_exact))

    assert np.allclose(flux, moved_flux, rtol=1e-10, atol=1e-10 * np.abs(flux).max())
    assert np.allclose(dense(moved_coef.T_cell), dense(base.T_cell), rtol=1e-10,
                       atol=1e-10 * np.abs(dense(base.T_cell)).max())


# ═════════════════════════════ Q — təsadüfi skew hallar (izahsız çökmə yox)
@pytest.mark.parametrize("seed", range(8))
def test_q_randomised_skew_cases_stay_exact_for_linear_fields(seed):
    rng = np.random.default_rng(1000 + seed)
    grid, geometry, _ = warped(3, 2, 2, amplitude=0.9, seed=seed)
    assert geometry.validate().ok, "test ssenarisi etibarsız həndəsə yaratdı"

    k_matrix = rotated_anisotropic(1, eigenvalues=tuple(rng.uniform(10.0, 800.0, 3)),
                                   axis=rng.normal(size=3),
                                   angle=float(rng.uniform(0.0, 180.0)))[0]
    coef = coefficients(grid, geometry, uniform_k(grid.ncell, k_matrix))
    gradient = rng.normal(size=3)
    pressures, exact = linear_field(geometry, gradient)
    flux = coef.face_fluxes(pressures, coef.boundary_pressures_from(exact))
    expected = analytic_flux(geometry, k_matrix, gradient)

    assert np.allclose(flux, expected, rtol=1e-9,
                       atol=1e-9 * np.abs(expected).max())
    assert coef.conservation_residual(
        pressures, coef.boundary_pressures_from(exact)) < 1e-9 * np.abs(flux).max()
    # Ssenarinin etibarlılığı: deformasiya HƏQİQİDİR (stensil çoxnöqtəlidir)
    interior = np.flatnonzero(~geometry.is_boundary)
    assert coef.stencil_sizes()[interior].min() > 2


# ════════════════════════════ R — TPFA reqressiyası (DƏYİŞMƏMƏLİDİR) ════
def test_r_tpfa_results_are_unchanged_by_mpfa_existence():
    """Tapşırıq §22 — TPFA-nın düsturu/nəticəsi TOXUNULMAZ."""
    model = five_spot_model(nx=3, ny=3, dx=10.0, dy=10.0, dz=5.0, permeability=200.0)
    grid = TwoPointFluxDiscretization().build(model)
    expected = METRIC.darcy_constant * 200.0 * (10.0 * 5.0) / 10.0
    assert np.allclose(grid.transmissibility, expected, rtol=1e-12)

    d_phi = np.linspace(-2.0, 3.0, grid.connections.count)
    assert np.allclose(grid.compute_flux(d_phi), grid.transmissibility * d_phi)


def test_r_mpfa_build_does_not_mutate_the_model():
    model = five_spot_model(nx=3, ny=3, dx=10.0, dy=10.0, dz=5.0, permeability=200.0)
    before = TwoPointFluxDiscretization().build(model).transmissibility.copy()
    permx_before = model.rock.permx.values.copy()

    MPFAODiscretization().build(model)

    after = TwoPointFluxDiscretization().build(model).transmissibility
    assert np.array_equal(before, after)
    assert np.array_equal(permx_before, model.rock.permx.values)


def test_r_mpfa_warns_about_unsupported_faults_instead_of_silently_ignoring():
    """Fay dəstəyi Phase 5D-dədir — SƏSSİZ yanlış nəticə QADAĞANDIR."""
    from imex2d.domain.structure import FaultReference
    model = five_spot_model(nx=3, ny=3, dx=10.0, dy=10.0, dz=5.0)
    model.fault_references.append(
        FaultReference(name="F1", source_id="f1", axis="I", plane_index=1,
                       transmissibility_multiplier=0.0))
    discretized = MPFAODiscretization().build(model)
    assert any("fay" in w.lower() for w in discretized.warnings)


# ══════════════════ §29 — ANTI-PSEUDO-MPFA (MƏCBURİ) ════════════════════
def test_anti_pseudo_mpfa_flux_is_not_a_two_point_relation():
    """Tapşırıq §29 — ƏN VACİB TEST.

    Qeyri-ortoqonal grid + fırlanmış tam tenzor. `q = T(p_i − p_j)`
    formasına YIĞILAN implementasiya bu testdən KEÇƏ BİLMƏZ: elə bir
    təzyiq sahəsi seçilir ki, `p_i = p_j` OLSUN — iki-nöqtəli sxem
    məcburi SIFIR verər, HƏQİQİ MPFA-O isə YOX.
    """
    grid, geometry, _ = skewed(3, 3, 2)
    coef = coefficients(grid, geometry, rotated_anisotropic(grid.ncell))
    interior = np.flatnonzero(~geometry.is_boundary)
    face = int(interior[np.argmax(coef.stencil_sizes()[interior])])
    owner = int(geometry.face_owner[face])
    neighbor = int(geometry.face_neighbor[face])

    stencil = coef.face_stencil(face, tolerance=1e-9)
    assert len(stencil) > 2, "stensil iki-nöqtəlidir — implementasiya MPFA-O DEYİL"

    rng = np.random.default_rng(2026)
    pressures = 100.0 + rng.normal(0.0, 20.0, grid.ncell)
    pressures[neighbor] = pressures[owner]          # ΔP = 0 QURULDU
    boundary = np.zeros(len(coef.boundary_dofs))
    flux = coef.face_fluxes(pressures, boundary)

    scale = np.abs(coef.face_fluxes(pressures, boundary)).max()
    assert abs(flux[face]) > 1e-6 * scale, (
        "p_i = p_j olduqda axın sıfır çıxdı — axın YALNIZ iki hüceyrədən "
        "asılıdır, yəni bu, TPFA-dır, MPFA-O deyil.")

    # Əlavə sübut: ÜÇÜNCÜ bir hüceyrənin təzyiqi TƏK BAŞINA axını dəyişir
    third = next(c for c in stencil if c not in (owner, neighbor))
    moved = pressures.copy()
    moved[third] += 25.0
    assert abs(coef.face_fluxes(moved, boundary)[face] - flux[face]) > 1e-6 * scale


def test_anti_pseudo_mpfa_no_hidden_tpfa_fallback_on_hard_geometry():
    """Çətin (çox skew + kəskin anizotrop) halda da stensil çoxnöqtəli
    QALIR — yəni heç bir gizli "TPFA-ya qayıt" yolu yoxdur."""
    grid, geometry, _ = skewed(3, 3, 2, shear=0.8)
    coef = coefficients(grid, geometry,
                        rotated_anisotropic(grid.ncell, eigenvalues=(2000.0, 20.0, 2.0),
                                            angle=63.0))
    interior = np.flatnonzero(~geometry.is_boundary)
    assert coef.stencil_sizes()[interior].min() > 2


# ═══════════ §30 — ƏL İLƏ YOXLANILA BİLƏN İSTİNAD HALI ══════════════════
def test_reference_case_unit_cubes_hand_verifiable_coefficients():
    """Tapşırıq §30 — başqa mühəndisin MÜSTƏQİL yoxlaya biləcəyi hal.

    HƏNDƏSƏ: 2×2×2 vahid kub (h = 1 m), mərkəzlər (±0.5, ±0.5, ±0.5).
    K       : I (izotrop, k = 1).
    Γ       : 1 (vahid Darsi sabiti).
    η       : 1 → kəsilməzlik nöqtəsi = üz mərkəzi.
    BÖLGƏ   : mərkəzi node (1,1,1) — 8 hüceyrə, 12 sub-üz, 12 naməlum π.

    ƏL HESABI (bax `docs/mpfa_o_phase5a.md` §16):
      sub-üz sahəsi   A_σ = (h/2)² = 0.25
      yarım məsafə    h_o = h_n = 0.5
      T_σ = Γ A_σ / (h_o/k + h_n/k) = 1·0.25/(0.5+0.5) = 0.25
      ⇒ T_cell[σ, owner] = +0.25,  T_cell[σ, neighbor] = −0.25, qalan 0

    Yoxlama: p = (1,0,0)·x sahəsində sub-üz axını
      q_σ = 0.25·(p_owner − p_neighbor) = 0.25·1.0 = 0.25
    """
    grid = CartesianGrid(2, 2, 2)
    geometry = CellGeometry(grid, dx=1.0, dy=1.0, dz=[1.0, 1.0], top_depth=0.0)
    vertices = hexahedral_vertices_from_cartesian(grid, geometry)
    ggg = GeneralGridGeometry(vertices, grid.build_connections())

    regions = build_interaction_regions(grid, ggg, eta=1.0)
    centre = [r for r in regions if not r.is_boundary_region]
    assert len(centre) == 1, "2×2×2 gridd-də DƏQİQ bir daxili node olmalıdır"
    region = centre[0]
    assert region.node_ijk == (1, 1, 1)
    assert region.n_cells == 8 and region.n_sub_faces == 12

    system = MPFAOLocalSystem(region=region, cell_centroids=ggg.cell_centroids,
                              k_matrices=uniform_k(8, np.eye(3)), darcy_constant=1.0)
    assert system.C.shape == (12, 12)
    assert system.T_cell.shape == (12, 8)
    assert system.T_bnd.shape[1] == 0            # daxili bölgə → sərhəd DOF yoxdur

    for s, sub in enumerate(region.sub_faces):
        assert np.isclose(sub.area, 0.25, rtol=1e-12)
        owner_col = region.cell_local[sub.owner]
        neighbor_col = region.cell_local[sub.neighbor]
        expected = np.zeros(8)
        expected[owner_col] = 0.25
        expected[neighbor_col] = -0.25
        assert np.allclose(system.T_cell[s], expected, atol=1e-12), (
            f"σ{s} əmsalları əl hesabı ilə uyğun gəlmir:\n{system.T_cell[s]}")

    pressures = np.array([ggg.cell_centroids[c][0] for c in region.cells])
    fluxes = system.sub_face_fluxes(pressures)
    for s, sub in enumerate(region.sub_faces):
        delta = ggg.cell_centroids[sub.owner][0] - ggg.cell_centroids[sub.neighbor][0]
        assert np.isclose(fluxes[s], 0.25 * delta, atol=1e-12)

    assert system.conservation_residual(pressures) < 1e-12
    assert "T_cell" in system.describe() and "D_(c,v)" in system.describe()


def test_reference_case_sub_face_area_vectors_sum_to_the_face_area_vector():
    """`Σ_q a_(F,u_q) = A_F n_F` — §2 tiling xassəsi (müstəvi üz üçün DƏQİQ)."""
    grid, geometry, _ = skewed(2, 2, 2)
    regions = build_interaction_regions(grid, geometry)
    accumulated = {}
    for region in regions:
        for sub in region.sub_faces:
            accumulated.setdefault(sub.face_index, np.zeros(3))
            accumulated[sub.face_index] += sub.area_vector

    for face, total in accumulated.items():
        expected = geometry.face_areas[face] * geometry.face_normals[face]
        assert np.allclose(total, expected, rtol=1e-10, atol=1e-10)


# ═══════════════════════ η parametri (§4/§16 sənədləşdirilmiş davranış) ══
def test_eta_half_stays_exact_but_widens_the_cartesian_stencil():
    """η=1/2 — xətti dəqiqlik SAXLANILIR, amma Kartezian stensili
    2-nöqtəli DEYİL. Bu, sənəddə (§4) AÇIQ yazılmış xassədir."""
    grid, geometry, _ = cartesian(3, 3, 2)
    k = isotropic(grid.ncell, 150.0)
    gradient = np.array([1.0, 2.0, -1.0])

    half = coefficients(grid, geometry, k, eta=0.5)
    pressures, exact = linear_field(geometry, gradient)
    flux = half.face_fluxes(pressures, half.boundary_pressures_from(exact))
    expected = analytic_flux(geometry, 150.0 * np.eye(3), gradient)
    assert np.allclose(flux, expected, rtol=1e-11, atol=1e-11 * np.abs(expected).max())

    interior = np.flatnonzero(~geometry.is_boundary)
    assert half.stencil_sizes()[interior].max() > 2
    full = coefficients(grid, geometry, k, eta=1.0)
    assert full.stencil_sizes()[interior].max() == 2


# ═══════════════════════════════════ performans (§31) ═══════════════════
def test_performance_scales_linearly_with_cell_count():
    """`O(N²)` qonşu axtarışı YOXDUR — bölgə sayı hüceyrə sayı ilə
    XƏTTİ artır və qurma vaxtı da təxminən xətti qalır."""
    timings = []
    for nx in (4, 8):
        grid, geometry, _ = cartesian(nx, nx, 2)
        k = isotropic(grid.ncell)
        started = time.perf_counter()
        coef = coefficients(grid, geometry, k)
        timings.append((grid.ncell, time.perf_counter() - started, len(coef.regions)))

    (n_small, t_small, r_small), (n_large, t_large, r_large) = timings
    assert r_large / r_small < 1.6 * (n_large / n_small)
    if t_small > 5e-3:                       # çox sürətli hallarda ölçmə küylüdür
        assert t_large / t_small < 4.0 * (n_large / n_small)


def test_performance_coefficient_storage_is_linear_not_quadratic():
    """`T_cell` SEYRƏKDİR: sətir başına sıfırdan fərqli əmsal sayı
    SABİT-məhduddur (§14: ≤18), ona görə yaddaş `O(N)`-dir.

    Sıx `(nface × ncell)` massiv `O(N²)` olardı və tapşırıq §31-i
    pozardı — bu test məhz həmin reqressiyanı tutur."""
    densities = []
    for nx in (4, 10):
        grid, geometry, _ = cartesian(nx, nx, 2)
        coef = coefficients(grid, geometry, rotated_anisotropic(grid.ncell))
        per_face = coef.T_cell.nnz / coef.n_face
        densities.append(per_face)
        assert per_face <= 18.0, f"sətir başına {per_face:.1f} əmsal — §14 həddini aşır"
        assert coef.T_cell.nnz < 0.5 * coef.n_face * coef.n_cell

    # Grid 6 dəfə böyüyəndə sətir sıxlığı DƏYİŞMİR (sabit stensil)
    assert densities[1] < 1.5 * densities[0]
