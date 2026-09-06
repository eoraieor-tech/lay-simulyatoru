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


# ═══════════════════════════════════════════════════════════════════════════
# FAZA 4 AUDİTİ — əlavə, DAHA SƏRT invariant testləri
# ═══════════════════════════════════════════════════════════════════════════
# Yuxarıdakı A–T testləri əsasən Kartezian (və AFFİN çevrilmiş, yəni bütün
# üzləri MÜSTƏVİ qalan) həndəsəni yoxlayırdı. Aşağıdakılar auditin §2/§3/
# §4/§6/§7/§10 tələblərini — paylaşılan üzlərin TƏKLİYİ, ORİYENTASİYALI
# sahə-vektor qapanması, ACTNUM indeks xəritələməsi, HƏQİQƏTƏN əyri
# (non-planar) grid, pinch-out və ETİBARSIZ `Connections` — ölçür.

from imex2d.domain.grid import Connections                             # noqa: E402
from imex2d.domain.polyhedral_geometry import HEX_FACE_VERTEX_INDICES  # noqa: E402

_OPP = {"X-": "X+", "X+": "X-", "Y-": "Y+", "Y+": "Y-", "Z-": "Z+", "Z+": "Z-"}


def _oriented_area_vector(face) -> np.ndarray:
    """Üzün HƏQİQİ oriyentasiyalı sahə vektoru `Σ_tri ½·(b−a)×(c−a)`.

    Divergensiya teoremi (`∮ n dS`) məhz BUNU tələb edir. `Face.area() *
    Face.normal()` isə SKALYAR sahə × VAHİD normaldır — MÜSTƏVİ üzdə
    ikisi eynidir, ƏYRİ (warped) üzdə isə FƏRQLİDİR (bax
    `polyhedral_geometry.py` docstring-i, "Normal barədə").
    """
    total = np.zeros(3)
    for a, b, c in face._triangles():
        total += 0.5 * np.cross(b - a, c - a)
    return total


def _closure_residuals(ggg, oriented: bool = True) -> np.ndarray:
    """`(ncell, 3)` — hər hüceyrə üçün YIĞILMIŞ `Σ_f (±A_f n_f)`.

    `HexahedralCell.closure_residual()`-dan FƏRQİ: bu, hüceyrənin
    `GeneralGridGeometry`-də QEYD OLUNMUŞ üzləri üzərindən gedir və
    qonşudan görünən üz üçün İŞARƏNİ çevirir — yəni TOPOLOGİYANI da
    yoxlayır, təkcə hüceyrə həndəsəsini yox.
    """
    residuals = np.zeros((ggg.ncell, 3))
    for cell in range(ggg.ncell):
        for face_index in ggg.cell_faces(cell):
            gf = ggg.faces[face_index]
            sign = 1.0 if gf.owner == cell else -1.0
            contribution = (_oriented_area_vector(gf.face) if oriented
                            else gf.face.area() * gf.face.normal())
            residuals[cell] += sign * contribution
    return residuals


def _scenario(name):
    """Auditin §6-da tələb olunan qeyri-Kartezian ssenariləri —
    `(vertices, connections, grid)`. Təpə pozulmaları GRID-səviyyəlidir
    (koordinatın FUNKSİYASI), ona görə PAYLAŞILAN təpələr hər iki
    hüceyrədə EYNİ qalır — qonşu hüceyrələr arasında boşluq açılmır."""
    actnum = None
    if name == "actnum-inactive":
        actnum = np.ones(3 * 3 * 1, dtype=np.int8)
        actnum[4] = 0                                   # mərkəzi hüceyrə söndürülüb
    elif name == "actnum-checkerboard":
        actnum = np.array([1, 0, 1, 0, 1, 0, 1, 0, 1], dtype=np.int8)

    nx, ny, nz = (3, 3, 1) if actnum is not None else (3, 2, 2)
    dz = [1.0] * nz if nz > 1 else 1.0
    if name == "variable-thickness":
        nx, ny, nz, dz = 3, 1, 3, [1.0, 3.0, 0.5]

    grid = CartesianGrid(nx, ny, nz, actnum)
    geometry = CellGeometry(grid, dx=2.0, dy=3.0, dz=dz, top_depth=0.0)
    conn = grid.build_connections()
    v = hexahedral_vertices_from_cartesian(grid, geometry)

    if name == "dipping":
        v[..., 2] += 0.5 * v[..., 0]                    # x-ə görə maili laylar
    elif name == "skewed":
        v[..., 0] += 0.4 * v[..., 2]                    # z-ə görə sürüşmə (shear)
    elif name == "non-orthogonal":
        v[..., 0] += 0.3 * v[..., 1]
        v[..., 1] += 0.2 * v[..., 2]
    elif name == "warped":
        # z pozulması (x,y) üzrə BİLİNEAR DEYİL -> Z± üzləri qeyri-müstəvi
        # olur; `(1 + 0.35·z)` amili isə hüceyrənin TAVANI ilə DÖŞƏMƏSİNİ
        # FƏRQLİ əyir — əks halda ikisi eyni səthin sürüşməsi olardı və
        # `A·n` töhfələri bir-birini SÜNİ şəkildə ixtisar edərdi.
        # Pozulma yalnız KOORDİNATIN funksiyasıdır, ona görə PAYLAŞILAN
        # təpə hər iki hüceyrədə EYNİ yerə düşür (boşluq açılmır).
        v[..., 2] += (0.3 * np.sin(1.7 * v[..., 0]) * np.cos(1.1 * v[..., 1])
                      * (1.0 + 0.35 * v[..., 2]))
    elif name == "boundary-only":
        empty = np.zeros(0, dtype=int)
        conn = Connections(empty, empty, np.zeros(0, dtype=np.int8))
    return v, conn, grid


_SCENARIOS = ["uniform", "variable-thickness", "dipping", "warped", "skewed",
              "non-orthogonal", "actnum-inactive", "actnum-checkerboard",
              "boundary-only"]


# ── U: hər (hüceyrə, yerli üz) DƏQİQ bir dəfə — duplikat üz YOXDUR ────────
@pytest.mark.parametrize("scenario", _SCENARIOS)
def test_u_every_local_face_slot_is_realised_exactly_once(scenario):
    v, conn, grid = _scenario(scenario)
    ggg = GeneralGridGeometry(v, conn)

    seen: dict = {}
    for gf in ggg.faces:
        assert gf.owner != gf.neighbor, f"Üz {gf.index}: owner == neighbor"
        slots = [(gf.owner, gf.owner_local_name)]
        if gf.neighbor is not None:
            slots.append((gf.neighbor, _OPP[gf.owner_local_name]))
        for slot in slots:
            assert slot not in seen, (
                f"{slot} İKİ üzdə görünür: {seen.get(slot)} və {gf.index}")
            seen[slot] = gf.index

    # 6 slot × ncell, heç biri əskik/artıq deyil
    assert len(seen) == 6 * ggg.ncell
    # daxili üz sayı DƏQİQ əlaqə sayına bərabərdir (deduplikasiya işləyir)
    assert int((~ggg.is_boundary).sum()) == conn.count
    assert len(ggg.faces) == 6 * ggg.ncell - conn.count
    for cell in range(ggg.ncell):
        assert len(ggg.cell_faces(cell)) == 6


# ── V: A_owner = −A_neighbor (eyni üz, eyni sahə, ƏKS oriyentasiya) ───────
@pytest.mark.parametrize("scenario", _SCENARIOS)
def test_v_shared_face_is_identical_from_both_sides(scenario):
    v, conn, grid = _scenario(scenario)
    ggg = GeneralGridGeometry(v, conn)

    for gf in ggg.faces:
        if gf.is_boundary:
            continue
        neighbor_face = ggg.cells[gf.neighbor].faces()[_OPP[gf.owner_local_name]]
        scale = max(gf.face.area(), 1.0)
        # EYNİ üz: sahə, mərkəz eyni; normal TAM ƏKS
        assert abs(gf.face.area() - neighbor_face.area()) <= 1e-12 * scale
        assert np.allclose(gf.face.centroid(), neighbor_face.centroid(), atol=1e-12)
        assert np.allclose(gf.face.normal(), -neighbor_face.normal(), atol=1e-12)
        # A_owner = −A_neighbor (ORİYENTASİYALI sahə vektoru)
        assert np.allclose(_oriented_area_vector(gf.face),
                           -_oriented_area_vector(neighbor_face), atol=1e-12 * scale)
        assert np.allclose(gf.normal_from(gf.owner), -gf.normal_from(gf.neighbor))


# ── W: YIĞILMIŞ qapanma Σ(A_f n_f) ≈ 0 (oriyentasiyalı sahə vektoru) ──────
@pytest.mark.parametrize("scenario", _SCENARIOS)
def test_w_assembled_geometric_conservation_holds(scenario):
    """Tolerans əsaslandırması: qapanma cəmi ~`max(A_f)` tərtibli ədədlərin
    fərqidir, ona görə MÜTLƏQ deyil, NİSBİ hədd (max üz sahəsinə görə)
    götürülür; `1e-12` ikili üzən-nöqtə toplamasının (`~6·eps·A ≈ 1.3e-15·A`)
    təxminən 700 qat üstündədir — HƏQİQİ topoloji səhvi (nisbət ~1e-2)
    isə rahat tutur."""
    v, conn, grid = _scenario(scenario)
    ggg = GeneralGridGeometry(v, conn)

    assert np.all(ggg.cell_volumes > 0.0)               # V > 0
    scale = float(ggg.face_areas.max())
    residual = np.linalg.norm(_closure_residuals(ggg, oriented=True), axis=1)
    assert residual.max() <= 1e-12 * scale, (
        f"{scenario}: qapanma pozulub, max |Σ A·n| = {residual.max():.3e}")


# ── X: skalyar-sahə konvensiyası ƏYRİ üzdə qapanmanı POZUR (SƏNƏDLİ) ──────
def test_x_scalar_area_convention_breaks_closure_on_warped_faces():
    """`area()` və `oriented_area_vector()` AYRI kəmiyyətlərdir (FINDING-3).

    `Face.area()` — üçbucaq sahələrinin SKALYAR cəmi (əyri səthin HƏQİQİ
    sahəsi), `Face.normal()` — vahid normal. Qeyri-müstəvi üzdə
    `area()·normal()` ≠ `A⃗`, ona görə ONUNLA qapanma SIFIR VERMİR —
    HƏNDƏSƏ isə əslində QAPALIDIR (`oriented=True` ölçüsü maşın
    dəqiqliyindədir).

    Bu test həmin FƏRQİ ölçür və onun ölçülə bilən qaldığını təsdiqləyir;
    məhz buna görə istehsalat kodu (`closure_residual`,
    `GridFace.area_vector`, MPFA-O `a⃗_σ`) artıq `area()·normal()` DEYİL,
    `oriented_area_vector()` işlədir — bax `tests/test_oriented_area_
    vector.py`. `area()` yalnız səthin FİZİKİ sahəsi lazım olanda qalır.
    """
    v, conn, grid = _scenario("warped")
    ggg = GeneralGridGeometry(v, conn)

    non_planar = [gf for gf in ggg.faces if not gf.face.is_planar(1e-9)]
    assert non_planar, "ssenari HƏQİQƏTƏN qeyri-müstəvi üz yaratmadı"

    scale = float(ggg.face_areas.max())
    oriented = np.linalg.norm(_closure_residuals(ggg, oriented=True), axis=1).max()
    scalar = np.linalg.norm(_closure_residuals(ggg, oriented=False), axis=1).max()

    assert oriented <= 1e-12 * scale        # HƏQİQİ həndəsə QAPALIDIR
    assert scalar > 1e-6 * scale            # skalyar konvensiya İSƏ deyil

    # müstəvi üzlü (affin) griddə ikisi ÜST-ÜSTƏ DÜŞÜR — Kartezian reqressiyası
    v2, conn2, _ = _scenario("skewed")
    affine = GeneralGridGeometry(v2, conn2)
    assert all(gf.face.is_planar(1e-9) for gf in affine.faces)
    assert np.linalg.norm(_closure_residuals(affine, oriented=False),
                          axis=1).max() <= 1e-12 * float(affine.face_areas.max())


# ── Y: ACTNUM / indeks xəritələməsi ──────────────────────────────────────
@pytest.mark.parametrize("scenario", ["actnum-inactive", "actnum-checkerboard"])
def test_y_actnum_index_mapping_is_global_and_consistent(scenario):
    """`Connections` QLOBAL indekslərlə işləyir, `GeneralGridGeometry` də
    QLOBAL `(ncell,8,3)` təpələrlə qurulur — deməli `geometry cell index`
    ≡ `global cell index`. `active` indeks YALNIZ xətti sistem
    sərhədindədir (bax `grid.py` docstring-i) və həndəsəyə SIZMAMALIDIR."""
    v, conn, grid = _scenario(scenario)
    ggg = GeneralGridGeometry(v, conn)
    active = grid.active

    assert ggg.ncell == grid.ncell                      # QLOBAL, aktiv DEYİL
    assert ggg.ncell != grid.n_active                   # ssenaridə həqiqətən qeyri-aktiv var

    faces = ggg.connection_faces()
    assert faces.shape == (conn.count,)
    assert len(set(faces.tolist())) == conn.count       # BİYEKSİYA
    for k, face_index in enumerate(faces):
        gf = ggg.faces[int(face_index)]
        assert gf.owner == int(conn.cell_a[k])
        assert gf.neighbor == int(conn.cell_b[k])
        assert not gf.is_boundary
        assert ggg.face_owner[face_index] == gf.owner
        assert ggg.face_neighbor[face_index] == gf.neighbor
        # owner→neighbor istiqaməti: n · (c_j − c_i) > 0
        assert np.dot(ggg.face_normals[face_index], ggg.d_ij(int(face_index))) > 0.0
        # AKTİV↔QEYRİ-AKTİV daxili üz HEÇ VAXT qurulmamalıdır
        assert active.is_active(gf.owner) and active.is_active(gf.neighbor)

    # qeyri-aktiv hüceyrə axın topologiyasında İŞTİRAK ETMİR
    for cell in np.flatnonzero(active.actnum == 0):
        assert ggg.neighbors(int(cell)) == []
        assert all(ggg.is_boundary_face(f) for f in ggg.cell_faces(int(cell)))
    # daxili üzlərin heç biri qeyri-aktiv hüceyrəyə toxunmur
    interior = ~ggg.is_boundary
    assert bool(active.is_active(ggg.face_owner[interior]).all())
    assert bool(active.is_active(ggg.face_neighbor[interior]).all())


# ── Z: analitik affin etalon (qutu / maili prizma / sürüşmüş prizma) ──────
@pytest.mark.parametrize("label,matrix", [
    ("box", np.diag([2.0, 3.0, 5.0])),
    ("inclined prism", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.5, 0.0, 1.0]]),
    ("sheared prism", [[1.0, 0.3, 0.0], [0.0, 1.0, 0.4], [0.0, 0.0, 1.0]]),
    ("shear xz", [[1.0, 0.0, 0.7], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
])
def test_z_affine_reference_matches_closed_form_geometry(label, matrix):
    """Vahid kubun `S` affin çevrilməsi üçün həndəsə QAPALI ŞƏKİLDƏ
    məlumdur: `V = |det S|`, `centroid = S·(½,½,½)`, üz sahə-vektoru
    `A = cof(S)·e = det(S)·S⁻ᵀ·e` (Nanson düsturu). Bu, TƏTBİQDƏN TAM
    MÜSTƏQİL etalondur."""
    S = np.asarray(matrix, float)
    unit = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                     [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
    cell = HexahedralCell(unit @ S.T)

    assert np.isclose(cell.volume(), abs(np.linalg.det(S)), rtol=0, atol=1e-12)
    assert np.allclose(cell.centroid(), S @ np.full(3, 0.5), atol=1e-12)

    cofactor = np.linalg.det(S) * np.linalg.inv(S).T
    unit_normals = {"X+": [1, 0, 0], "X-": [-1, 0, 0], "Y+": [0, 1, 0],
                    "Y-": [0, -1, 0], "Z+": [0, 0, 1], "Z-": [0, 0, -1]}
    for name, face in cell.faces().items():
        expected = cofactor @ np.array(unit_normals[name], float)
        assert np.isclose(face.area(), np.linalg.norm(expected), atol=1e-12)
        assert np.allclose(face.normal(), expected / np.linalg.norm(expected),
                           atol=1e-12)
        assert np.allclose(_oriented_area_vector(face), expected, atol=1e-12)


# ── AA: pinch-out (sıfır qalınlıqlı lay) AÇIQ şəkildə bildirilir ──────────
def test_aa_pinch_out_layer_is_reported_not_silently_accepted():
    grid = CartesianGrid(2, 1, 3)
    geometry = CellGeometry(grid, dx=2.0, dy=3.0, dz=[1.0, 1.0, 1.0], top_depth=0.0)
    conn = grid.build_connections()
    v = hexahedral_vertices_from_cartesian(grid, geometry)
    middle = np.array([grid.ijk(c)[2] == 1 for c in range(grid.ncell)])
    v[middle, :, 2] = v[middle, 0, 2][:, None]          # orta lay -> sıfır qalınlıq

    ggg = GeneralGridGeometry(v, conn)
    assert np.any(ggg.cell_volumes <= 0.0)
    result = ggg.validate()
    assert not result.ok
    assert any("həcm" in message for message in result.errors)


# ── AB: ETİBARSIZ `Connections` SƏSSİZ qəbul edilmir ──────────────────────
def _linear_three_cells():
    grid = CartesianGrid(3, 1, 1)
    geometry = CellGeometry(grid, dx=2.0, dy=3.0, dz=1.0, top_depth=0.0)
    return grid, hexahedral_vertices_from_cartesian(grid, geometry)


def _conn(a, b, axis=0):
    return Connections(np.array(a, dtype=int), np.array(b, dtype=int),
                       np.full(len(a), axis, dtype=np.int8))


def test_ab_duplicate_connection_is_reported():
    """REQRESSİYA: eyni hüceyrə cütü İKİ dəfə verildikdə əvvəllər
    `validate()` `ok=True` qaytarırdı, HALBUKİ hüceyrə 7 üz alırdı,
    `_cell_local_face[(0,"X+")]` üzərinə YAZILIRDI və yığılmış
    `Σ A·n` = ±(bir üzün sahəsi) ≠ 0 olurdu — SƏSSİZ yanlış həndəsə."""
    grid, v = _linear_three_cells()
    ggg = GeneralGridGeometry(v, _conn([0, 0], [1, 1]))

    result = ggg.validate()
    assert not result.ok, "təkrarlanan əlaqə SƏSSİZ qəbul edildi"
    assert any("İKİ DƏFƏ" in message for message in result.errors), result.errors


def test_ab2_self_connection_is_reported():
    grid, v = _linear_three_cells()
    result = GeneralGridGeometry(v, _conn([1], [1])).validate()
    assert not result.ok
    assert any("owner == neighbor" in message for message in result.errors), result.errors


@pytest.mark.parametrize("a,b", [([0], [3]), ([3], [1]), ([0], [-1]), ([-2], [1])])
def test_ab3_out_of_range_connection_raises_clear_error(a, b):
    """REQRESSİYA: MƏNFİ indeks Python-un dövri indeksləməsi ilə SƏSSİZ
    olaraq SON hüceyrəyə bağlanırdı (`_neighbor_lists[-1]`), üstəlik
    `face_neighbor` massivində `-1` "sərhəd" sentineli ilə qarışırdı —
    yəni EYNİ üz həm daxili (qonşu siyahısında), həm də sərhəd
    (`is_boundary`) görünürdü."""
    grid, v = _linear_three_cells()
    with pytest.raises(ValueError, match="hüceyrə indeksi"):
        GeneralGridGeometry(v, _conn(a, b))


def test_ab4_inconsistent_axis_is_reported_by_validation():
    """Qonşu OLMAYAN (və ya yanlış oxda elan olunmuş) əlaqə — həndəsə
    uyğunsuzluğu `validate()`-də tutulur."""
    grid, v = _linear_three_cells()
    assert not GeneralGridGeometry(v, _conn([0], [1], axis=2)).validate().ok
    assert not GeneralGridGeometry(v, _conn([0], [2], axis=0)).validate().ok


# ── AC: yerli↔qlobal üz axtarışı hər iki tərəfdən EYNİ üzü verir ──────────
@pytest.mark.parametrize("scenario", _SCENARIOS)
def test_ac_face_index_round_trip_is_symmetric(scenario):
    v, conn, grid = _scenario(scenario)
    ggg = GeneralGridGeometry(v, conn)
    for cell in range(ggg.ncell):
        for name in HEX_FACE_VERTEX_INDICES:
            face_index = ggg.face_index(cell, name)
            gf = ggg.faces[face_index]
            assert cell in (gf.owner, gf.neighbor)
            if gf.is_boundary:
                assert gf.owner == cell and gf.owner_local_name == name
            else:
                other = gf.neighbor if gf.owner == cell else gf.owner
                assert ggg.face_index(other, _OPP[name]) == face_index


# ── AD: normal HƏMİŞƏ owner-dan KƏNARA, owner→neighbor istiqamətindədir ───
@pytest.mark.parametrize("scenario", _SCENARIOS)
def test_ad_all_normals_point_outward_from_their_owner(scenario):
    v, conn, grid = _scenario(scenario)
    ggg = GeneralGridGeometry(v, conn)
    for gf in ggg.faces:
        outward = np.dot(gf.face.normal(),
                         gf.face.centroid() - ggg.cell_centroids[gf.owner])
        assert outward > 0.0, f"Üz {gf.index} owner {gf.owner}-dan KƏNARA baxmır"
        if not gf.is_boundary:
            assert np.dot(ggg.face_normals[gf.index], ggg.d_ij(gf.index)) > 0.0
            assert np.dot(gf.normal_from(gf.neighbor), -ggg.d_ij(gf.index)) > 0.0


# ── AE: vektorlaşdırılmış təpə qurucusu skalyar etalona BƏRABƏRDİR ────────
def test_ae_vectorised_vertices_match_scalar_reference():
    grid = CartesianGrid(4, 3, 3)
    geometry = CellGeometry(grid, dx=7.0, dy=11.0, dz=[3.0, 5.0, 2.0],
                            top_depth=1234.0)
    vectorised = hexahedral_vertices_from_cartesian(grid, geometry)

    layer_top = np.concatenate(([0.0], np.cumsum(np.atleast_1d(geometry.dz))))
    reference = np.zeros((grid.ncell, 8, 3))
    for c in range(grid.ncell):
        i, j, k = grid.ijk(c)
        x0, x1 = i * geometry.dx, (i + 1) * geometry.dx
        y0, y1 = j * geometry.dy, (j + 1) * geometry.dy
        z0 = geometry.top_depth + layer_top[k]
        z1 = geometry.top_depth + layer_top[k + 1]
        reference[c] = [[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                        [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]]
    assert np.array_equal(vectorised, reference)

    # üz sırası DETERMİNİSTİKDİR — eyni giriş, BİT-BƏRABƏR eyni çıxış
    conn = grid.build_connections()
    first = GeneralGridGeometry(vectorised, conn)
    second = GeneralGridGeometry(vectorised, conn)
    assert np.array_equal(first.face_owner, second.face_owner)
    assert np.array_equal(first.face_neighbor, second.face_neighbor)
    assert np.array_equal(first.face_areas, second.face_areas)
    assert np.array_equal(first.face_normals, second.face_normals)
    assert np.array_equal(first.connection_faces(), second.connection_faces())


# ── AF: Kartezian analitik etalon (həcm, mərkəz, sahə, normal) ────────────
def test_af_cartesian_regression_against_closed_form_values():
    grid = CartesianGrid(3, 2, 2)
    geometry = CellGeometry(grid, dx=7.0, dy=11.0, dz=[3.0, 5.0], top_depth=1234.0)
    conn = grid.build_connections()
    ggg = GeneralGridGeometry(hexahedral_vertices_from_cartesian(grid, geometry),
                              conn)

    i, j, k = grid.ijk_array(np.arange(grid.ncell))
    layer_top = np.concatenate(([0.0], np.cumsum(geometry.dz)))
    expected_centroids = np.stack([
        (i + 0.5) * geometry.dx, (j + 0.5) * geometry.dy,
        geometry.top_depth + 0.5 * (layer_top[k] + layer_top[k + 1])], axis=-1)

    assert np.allclose(ggg.cell_volumes, geometry.volumes(), rtol=0, atol=1e-9)
    assert np.allclose(ggg.cell_centroids, expected_centroids, atol=1e-9)

    axis_unit = {0: [1.0, 0.0, 0.0], 1: [0.0, 1.0, 0.0], 2: [0.0, 0.0, 1.0]}
    for k_conn, face_index in enumerate(ggg.connection_faces()):
        axis = int(conn.axis[k_conn])
        gf = ggg.faces[int(face_index)]
        layer = grid.ijk(gf.owner)[2]
        expected_area = {0: geometry.dy * geometry.dz[layer],
                         1: geometry.dx * geometry.dz[layer],
                         2: geometry.dx * geometry.dy}[axis]
        assert np.isclose(gf.face.area(), expected_area, rtol=0, atol=1e-9)
        assert np.allclose(gf.face.normal(), axis_unit[axis], atol=1e-12)


# ── AG: degenerativ / patoloji hüceyrə həndəsəsi (audit §7) ───────────────
def _unit_cell_vertices() -> np.ndarray:
    return np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                     [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)


def _mutate(kind: str) -> np.ndarray:
    v = _unit_cell_vertices()
    if kind == "zero-volume":
        v[4:] = v[:4]                       # tavan döşəməyə çökür
    elif kind == "extremely-thin":
        v[4:, 2] = 1e-14                    # dz -> maşın sıfırı
    elif kind == "inverted":
        v[4:, 2] = -1.0                     # mənfi (tərs-yönümlü) həcm
    elif kind == "all-vertices-identical":
        v[:] = v[0]
    elif kind == "nan-vertex":
        v[0, 0] = np.nan
    elif kind == "inf-vertex":
        v[2, 1] = np.inf
    elif kind == "collapsed-face":
        v[1] = v[0]
        v[5] = v[4]                         # Y- üzü XƏTTƏ çökür (sahə = 0)
    elif kind == "collapsed-edge":
        v[1] = v[0]                         # YALNIZ bir təpə birləşir -> qanuni paz
    elif kind == "thin-but-finite":
        v[4:, 2] = 1e-6
    else:                                   # pragma: no cover
        raise AssertionError(kind)
    return v


@pytest.mark.parametrize("kind", ["zero-volume", "extremely-thin", "inverted",
                                  "all-vertices-identical", "nan-vertex",
                                  "inf-vertex", "collapsed-face"])
def test_ag_degenerate_cells_are_rejected_with_an_explicit_error(kind):
    """Audit §7 — bunların HEÇ BİRİ SƏSSİZ qəbul edilməməlidir."""
    with np.errstate(invalid="ignore"):
        result = GeneralGridGeometry(_mutate(kind)[None, ...], None).validate()
    assert not result.ok, f"{kind} SƏSSİZ qəbul edildi"
    assert result.errors, f"{kind} üçün AÇIQ xəta mesajı yoxdur"


@pytest.mark.parametrize("kind", ["collapsed-edge", "thin-but-finite"])
def test_ag2_legitimate_degenerate_hexahedra_are_accepted(kind):
    """SƏNƏDLİ HƏDD: bir təpənin qonşusuna çökməsi (paz/pinch-out kənarı)
    corner-point gridlərdə QANUNİDİR — həcm/üz sahələri hələ də müsbətdir,
    ona görə `validate()` bunu xəta SAYMIR. Eyni şəkildə `1e-12`-dən
    QALIN, amma çox nazik hüceyrə də qəbul edilir: hədd MÜTLƏQDİR
    (`volume <= 1e-12`), NİSBİ deyil — yəni SAHƏ miqyasında (dx~100 m)
    "nazik" hüceyrə üçün fiziki filtr DEYİL, yalnız ədədi sıfır
    mühafizəsidir. Bu test həmin həddi AÇIQ şəkildə sabitləyir."""
    ggg = GeneralGridGeometry(_mutate(kind)[None, ...], None)
    assert ggg.cell_volumes[0] > 0.0
    assert ggg.validate().ok
    # qanuni olsa da, QAPANMA hələ də qorunur
    assert np.linalg.norm(_closure_residuals(ggg, oriented=True)[0]) <= 1e-12
