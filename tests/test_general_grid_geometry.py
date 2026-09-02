"""Phase 4 — "Grid-Level General Geometry Integration" test suite.

Audit tapşırığı §30-da tələb olunan A-T (20) test kateqoriyası, dəqiq
adlandırılmış ssenarilərlə + §29 performans reqressiyası.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from imex2d.domain.general_grid_geometry import (GeneralGridGeometry,
                                                  hexahedral_vertices_from_cartesian)
from imex2d.domain.geometry import CellGeometry
from imex2d.domain.grid import CartesianGrid
from imex2d.domain.polyhedral_geometry import Face, HexahedralCell


def _cartesian_ggg(nx=3, ny=2, nz=2, dx=10.0, dy=15.0, dz=None, top_depth=2000.0):
    dz = dz if dz is not None else [4.0, 6.0][:nz] or 5.0
    grid = CartesianGrid(nx, ny, nz)
    geometry = CellGeometry(grid, dx=dx, dy=dy, dz=dz, top_depth=top_depth)
    conn = grid.build_connections()
    vertices = hexahedral_vertices_from_cartesian(grid, geometry)
    return GeneralGridGeometry(vertices, conn), grid, geometry, conn, vertices


def _rotation_z(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


# ── A: Kartezian grid ────────────────────────────────────────────────────
def test_a_cartesian_grid_matches_cell_geometry_exactly():
    ggg, grid, geometry, conn, _ = _cartesian_ggg()
    assert np.allclose(ggg.cell_volumes, geometry.volumes())
    assert np.allclose(ggg.cell_centroids[:, 2], geometry.cell_depths())
    assert len(ggg.faces) - int((~ggg.is_boundary).sum()) == int(ggg.is_boundary.sum())
    assert int((~ggg.is_boundary).sum()) == conn.count   # HƏR əlaqə == 1 daxili üz
    assert ggg.validate().ok


# ── B: fırlanmış grid ─────────────────────────────────────────────────────
def test_b_rotated_grid_preserves_invariants():
    ggg, *_ , conn, vertices = _cartesian_ggg()
    rotation = _rotation_z(np.deg2rad(30.0))
    rotated = GeneralGridGeometry(vertices @ rotation.T, conn)

    assert np.allclose(rotated.cell_volumes, ggg.cell_volumes)
    assert np.allclose(rotated.cell_centroids, ggg.cell_centroids @ rotation.T)
    assert np.allclose(rotated.face_areas, ggg.face_areas)
    assert np.allclose(rotated.face_normals, ggg.face_normals @ rotation.T)
    assert np.array_equal(rotated.is_boundary, ggg.is_boundary)
    assert np.array_equal(rotated.face_owner, ggg.face_owner)
    assert np.array_equal(rotated.face_neighbor, ggg.face_neighbor)


# ── C: sürüşdürülmüş (translated) grid ───────────────────────────────────
def test_c_translated_grid_preserves_invariants():
    ggg, *_ , conn, vertices = _cartesian_ggg()
    t = np.array([100.0, -50.0, 2000.0])
    translated = GeneralGridGeometry(vertices + t, conn)

    assert np.allclose(translated.cell_volumes, ggg.cell_volumes)
    assert np.allclose(translated.face_areas, ggg.face_areas)
    assert np.allclose(translated.face_normals, ggg.face_normals)
    assert np.allclose(translated.cell_centroids, ggg.cell_centroids + t)
    assert np.array_equal(translated.face_owner, ggg.face_owner)
    assert np.array_equal(translated.face_neighbor, ggg.face_neighbor)


# ── D: miqyaslanmış (scaled) grid ─────────────────────────────────────────
@pytest.mark.parametrize("scale", [0.1, 1.0, 10.0])
def test_d_scaled_grid_follows_dimensional_analysis(scale):
    ggg, *_ , conn, vertices = _cartesian_ggg()
    scaled = GeneralGridGeometry(vertices * scale, conn)

    assert np.allclose(scaled.cell_volumes, ggg.cell_volumes * scale ** 3)
    assert np.allclose(scaled.face_areas, ggg.face_areas * scale ** 2)
    assert np.allclose(scaled.face_normals, ggg.face_normals)   # normal vahid vektor, DƏYİŞMİR
    assert np.allclose(scaled.cell_centroids, ggg.cell_centroids * scale)


# ── E: skewed/qeyri-ortoqonal grid ────────────────────────────────────────
def test_e_skewed_grid_valid_geometry_and_diagnostics():
    grid = CartesianGrid(3, 1, 1)
    geometry = CellGeometry(grid, dx=1.0, dy=1.0, dz=1.0, top_depth=0.0)
    conn = grid.build_connections()
    vertices = hexahedral_vertices_from_cartesian(grid, geometry)

    skewed = vertices.copy()
    top_mask = skewed[..., 2] > 0.5
    skewed[..., 0] += np.where(top_mask, 0.3, 0.0)   # bütün lövhəni Z-də sürüşdür
    ggg = GeneralGridGeometry(skewed, conn)

    assert ggg.validate().ok
    metrics = ggg.quality_metrics()
    assert metrics["max_non_orthogonality_angle_deg"] > 0.0   # artıq ORTOQONAL DEYİL
    assert np.isclose(metrics["max_non_orthogonality_angle_deg"],
                      np.degrees(np.arctan(0.3 / 1.0)), atol=1e-6)
    for gf in ggg.faces:
        if not gf.is_boundary:
            assert np.allclose(gf.normal_from(gf.owner), -gf.normal_from(gf.neighbor))


# ── F: daxili paylaşılan üz ────────────────────────────────────────────────
def test_f_interior_shared_face_is_single_physical_face():
    ggg, grid, *_ = _cartesian_ggg(nx=3, ny=1, nz=1)
    interior_faces = [gf for gf in ggg.faces if not gf.is_boundary]
    assert len(interior_faces) == 2   # 3 hüceyrə xətti sırada -> 2 daxili üz

    gf = interior_faces[0]
    assert np.allclose(gf.normal_from(gf.owner), -gf.normal_from(gf.neighbor))
    # TƏK Face obyekti — owner VƏ neighbor EYNİ sahə/mərkəzi paylaşır
    neighbor_local = {"X-": "X+", "X+": "X-", "Y-": "Y+", "Y+": "Y-",
                      "Z-": "Z+", "Z+": "Z-"}[gf.owner_local_name]
    neighbor_face = ggg.cells[gf.neighbor].faces()[neighbor_local]
    assert np.isclose(gf.face.area(), neighbor_face.area())
    assert np.allclose(gf.face.centroid(), neighbor_face.centroid())


# ── G: sərhəd üzləri ───────────────────────────────────────────────────────
def test_g_boundary_faces_have_single_owner_and_valid_geometry():
    ggg, grid, *_ = _cartesian_ggg(nx=2, ny=2, nz=1)
    boundary_faces = [gf for gf in ggg.faces if gf.is_boundary]
    assert len(boundary_faces) > 0
    for gf in boundary_faces:
        assert gf.neighbor is None
        assert gf.face.area() > 0.0
        assert np.all(np.isfinite(gf.face.centroid()))
        assert np.isclose(np.linalg.norm(gf.face.normal()), 1.0)
        assert ggg.is_boundary_face(gf.index)


# ── H: üz orientasiyası ───────────────────────────────────────────────────
def test_h_face_orientation_points_outward_from_owner():
    cell = HexahedralCell(hexahedral_vertices_from_cartesian(
        CartesianGrid(1, 1, 1), CellGeometry(CartesianGrid(1, 1, 1), 2.0, 2.0, 2.0))[0])
    centroid = cell.centroid()
    for face in cell.faces().values():
        direction = np.dot(face.normal(), face.centroid() - centroid)
        assert direction > 0.0   # qabarıq hüceyrə üçün HƏMİŞƏ doğrudur (bax audit §17)


# ── I: hüceyrə həcmi ───────────────────────────────────────────────────────
def test_i_cell_volume_matches_analytical_box_volume():
    ggg, grid, geometry, *_ = _cartesian_ggg(nx=2, ny=2, nz=2, dx=3.0, dy=4.0, dz=[5.0, 6.0])
    expected = geometry.volumes()
    assert np.allclose(ggg.cell_volumes, expected)


# ── J: üz sahəsi ───────────────────────────────────────────────────────────
def test_j_face_area_matches_analytical_value():
    ggg, grid, geometry, conn, _ = _cartesian_ggg(nx=2, ny=1, nz=1, dx=3.0, dy=4.0, dz=5.0)
    expected_areas = geometry.face_areas(conn)
    interior = [gf for gf in ggg.faces if not gf.is_boundary]
    assert len(interior) == 1
    assert np.isclose(interior[0].face.area(), expected_areas[0])
    assert np.isclose(interior[0].face.area(), 4.0 * 5.0)   # dy*dz, X-istiqamətli üz


# ── K: üz mərkəzi ──────────────────────────────────────────────────────────
def test_k_face_centroid_matches_analytical_value():
    ggg, grid, geometry, conn, _ = _cartesian_ggg(nx=2, ny=1, nz=1, dx=3.0, dy=4.0, dz=5.0,
                                                  top_depth=0.0)
    interior = [gf for gf in ggg.faces if not gf.is_boundary][0]
    expected = np.array([3.0, 2.0, 2.5])   # x=dx (hüceyrələr arası sərhəd), y=dy/2, z=dz/2
    assert np.allclose(interior.face.centroid(), expected)


# ── L: hüceyrə mərkəzi ─────────────────────────────────────────────────────
def test_l_cell_centroid_matches_analytical_value():
    ggg, grid, geometry, *_ = _cartesian_ggg(nx=2, ny=1, nz=1, dx=3.0, dy=4.0, dz=5.0,
                                             top_depth=100.0)
    assert np.allclose(ggg.cell_centroids[0], [1.5, 2.0, 102.5])
    assert np.allclose(ggg.cell_centroids[1], [4.5, 2.0, 102.5])


# ── M: qapanma münasibəti Σ A·n ≈ 0 ────────────────────────────────────────
def test_m_closure_relation_holds_near_machine_precision():
    ggg, *_ = _cartesian_ggg(nx=3, ny=2, nz=2)
    for cell in ggg.cells:
        residual = cell.closure_residual()
        assert np.linalg.norm(residual) < 1e-10

    # skewed hüceyrə üçün DƏ qapanma qorunur (affinə çevirmə planarlığı saxlayır)
    grid = CartesianGrid(1, 1, 1)
    geometry = CellGeometry(grid, dx=1.0, dy=1.0, dz=1.0)
    vertices = hexahedral_vertices_from_cartesian(grid, geometry)[0].copy()
    vertices[4:] += [0.4, 0.0, 0.0]
    skewed_cell = HexahedralCell(vertices)
    assert np.linalg.norm(skewed_cell.closure_residual()) < 1e-10


# ── N: qonşu xəritələməsi ──────────────────────────────────────────────────
def test_n_neighbor_mapping_matches_connections_topology():
    ggg, grid, geometry, conn, _ = _cartesian_ggg(nx=3, ny=2, nz=1)
    for k in range(conn.count):
        a, b = int(conn.cell_a[k]), int(conn.cell_b[k])
        assert b in ggg.neighbors(a)
        assert a in ggg.neighbors(b)
    # təcrid olunmuş (əlaqəsiz) hüceyrə — bax audit §25, ACTNUM hazırlığı
    isolated = GeneralGridGeometry(hexahedral_vertices_from_cartesian(
        CartesianGrid(1, 1, 1), CellGeometry(CartesianGrid(1, 1, 1), 1.0, 1.0, 1.0)),
        connections=None)
    assert isolated.neighbors(0) == []
    assert bool(isolated.is_boundary.all())


# ── O: hüceyrə-üz xəritələməsi ─────────────────────────────────────────────
def test_o_cell_face_mapping_is_consistent():
    ggg, *_ = _cartesian_ggg(nx=2, ny=2, nz=1)
    for cell in range(ggg.ncell):
        face_indices = ggg.cell_faces(cell)
        assert len(face_indices) == 6   # hər hekzahedral hüceyrənin DƏQİQ 6 üzü var
        for idx in face_indices:
            gf = ggg.faces[idx]
            assert cell in (gf.owner, gf.neighbor)


# ── P: etibarsız həndəsə ────────────────────────────────────────────────────
def test_p_invalid_degenerate_geometry_fails_validation():
    ggg, grid, geometry, conn, vertices = _cartesian_ggg(nx=3, ny=1, nz=1)
    bad = vertices.copy()
    bad[1, 4:] = bad[1, :4]   # hüceyrə 1-in yuxarı təpələrini aşağıya çökərt -> sıfır həcm
    bad_ggg = GeneralGridGeometry(bad, conn)
    result = bad_ggg.validate()
    assert not result.ok
    assert any("həcm" in e for e in result.errors)


# ── Q: NaN/Inf ──────────────────────────────────────────────────────────────
def test_q_nan_and_inf_geometry_fails_validation():
    """`np.errstate` YALNIZ konsolu təmizləyir — NaN-ın ÖZÜ `validate()`-də
    hələ də AÇIQ xəta kimi tutulur (bax bu testin assert-ləri)."""
    ggg, grid, geometry, conn, vertices = _cartesian_ggg(nx=2, ny=1, nz=1)
    with np.errstate(invalid="ignore"):
        nan_vertices = vertices.copy()
        nan_vertices[0, 0, 0] = np.nan
        assert not GeneralGridGeometry(nan_vertices, conn).validate().ok

        inf_vertices = vertices.copy()
        inf_vertices[0, 0, 0] = np.inf
        assert not GeneralGridGeometry(inf_vertices, conn).validate().ok


# ── R: degenerativ hüceyrə (təkrarlanan test, əlavə forma) ──────────────────
def test_r_zero_thickness_layer_is_rejected():
    grid = CartesianGrid(1, 1, 1)
    geometry = CellGeometry(grid, dx=1.0, dy=1.0, dz=1.0)
    vertices = hexahedral_vertices_from_cartesian(grid, geometry)
    vertices[0, 4:, 2] = vertices[0, :4, 2]   # dz -> 0
    result = GeneralGridGeometry(vertices, None).validate()
    assert not result.ok


# ── S: performans miqyaslanması ─────────────────────────────────────────────
def test_s_construction_scales_linearly_with_cell_count():
    def _time_build(n_per_axis):
        grid = CartesianGrid(n_per_axis, n_per_axis, 1)
        geometry = CellGeometry(grid, dx=10.0, dy=10.0, dz=5.0)
        conn = grid.build_connections()
        vertices = hexahedral_vertices_from_cartesian(grid, geometry)
        start = time.perf_counter()
        GeneralGridGeometry(vertices, conn)
        return time.perf_counter() - start

    _time_build(5)             # isinmə
    small = _time_build(15)    # 225 hüceyrə
    large = _time_build(60)    # 3600 hüceyrə (16x)
    ratio = large / max(small, 1e-9)
    # O(N) gözlənilir (~16x); O(N²) olsaydı ~256x — 50x həddi ikisini aydın ayırır
    assert ratio < 50, f"Qurma vaxtı N-dən SUPERXƏTTİ artır (nisbət={ratio:.1f})"


# ── T: TPFA reqressiyası ─────────────────────────────────────────────────────
def test_t_general_grid_geometry_does_not_affect_tpfa_regression():
    """Bu modul heç bir mövcud sinif tərəfindən ÇAĞIRILMIR (bax audit §28)
    — mövcud TPFA qızıl-etalon nəticəsi DƏYİŞMƏDƏN qalmalıdır."""
    import os
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from helpers import REFERENCE_FIVE_SPOT, default_scal, five_spot_model, make_service
    from imex2d.application.config import SimulationConfig

    scal = default_scal()
    model = five_spot_model(scal=scal)
    result = make_service(scal).run(model, SimulationConfig(end_time=1500.0))

    assert abs(result.ooip - REFERENCE_FIVE_SPOT["ooip"]) < 1.0
    assert result.steps == REFERENCE_FIVE_SPOT["steps"]
    assert result.converged
