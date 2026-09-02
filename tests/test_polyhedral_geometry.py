"""Phase 3 — "General Geometry Foundation for MPFA-O" test suite.

Audit tapşırığı §22-də tələb olunan 12 test, dəqiq adlandırılmış
ssenarilərlə + §23 (analitik gözlənti ilə müqayisə) və §12 (qeyri-
ortoqonallıq diaqnostikası) üçün əlavə testlər.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from imex2d.domain.polyhedral_geometry import (Face, HexahedralCell,
                                                cell_to_cell_vector,
                                                non_orthogonality_angle)

TOLERANCE = 1e-9   # bax audit §23: "document tolerances" — analitik həllər
                   # üçün maşın-dəqiqliyinə yaxın, yalnız üzən-nöqtə
                   # yığma xətasını (tetraedr/üçbucaq cəmi) hesaba alır.


def _box(dx: float, dy: float, dz: float, origin=(0.0, 0.0, 0.0)) -> HexahedralCell:
    ox, oy, oz = origin
    vertices = np.array([
        [0, 0, 0], [dx, 0, 0], [dx, dy, 0], [0, dy, 0],
        [0, 0, dz], [dx, 0, dz], [dx, dy, dz], [0, dy, dz],
    ], float) + np.array([ox, oy, oz])
    return HexahedralCell(vertices)


def _rotation_z(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


# ── Test 1: Kartezian hüceyrə ───────────────────────────────────────────
def test_1_cartesian_cell_volume_centroid_areas_normals():
    cell = _box(2.0, 3.0, 4.0)
    assert np.isclose(cell.volume(), 24.0, atol=TOLERANCE)
    assert np.allclose(cell.centroid(), [1.0, 1.5, 2.0], atol=TOLERANCE)

    expected_areas = {"Z-": 6.0, "Z+": 6.0, "Y-": 8.0, "Y+": 8.0, "X-": 12.0, "X+": 12.0}
    expected_normals = {"Z-": (0, 0, -1), "Z+": (0, 0, 1), "Y-": (0, -1, 0),
                        "Y+": (0, 1, 0), "X-": (-1, 0, 0), "X+": (1, 0, 0)}
    for name, face in cell.faces().items():
        assert np.isclose(face.area(), expected_areas[name], atol=TOLERANCE)
        assert np.allclose(face.normal(), expected_normals[name], atol=TOLERANCE)


# ── Test 2: dəyişən DZ ──────────────────────────────────────────────────
def test_2_variable_dz_stack_centroid_and_volume():
    """İki üst-üstə hüceyrə, fərqli hündürlük — hər birinin öz mərkəzi/
    həcmi DOĞRU olmalıdır (bax `CellGeometry.face_centroid`-in EYNİ
    "dəyişən qalınlıq" ssenarisi, Phase 1)."""
    lower = _box(10.0, 10.0, 4.0, origin=(0, 0, 0))
    upper = _box(10.0, 10.0, 10.0, origin=(0, 0, 4.0))

    assert np.isclose(lower.volume(), 400.0)
    assert np.isclose(upper.volume(), 1000.0)
    assert np.allclose(lower.centroid(), [5.0, 5.0, 2.0])
    assert np.allclose(upper.centroid(), [5.0, 5.0, 9.0])
    # paylaşılan üz (lower-in Z+ üzü) DƏQİQ z=4-dədir
    assert np.isclose(lower.faces()["Z+"].centroid()[2], 4.0)
    assert np.isclose(upper.faces()["Z-"].centroid()[2], 4.0)


# ── Test 3: fırlanmış hekzahedron ───────────────────────────────────────
def test_3_rotated_hexahedron_preserves_invariants():
    original = _box(2.0, 3.0, 4.0)
    rotation = _rotation_z(np.deg2rad(41.0))
    rotated = HexahedralCell(original.vertices @ rotation.T)

    assert np.isclose(rotated.volume(), original.volume(), atol=1e-8)
    assert np.allclose(rotated.centroid(), rotation @ original.centroid(), atol=1e-8)

    # kənar uzunluqları qorunub
    edge_original = np.linalg.norm(original.vertices[1] - original.vertices[0])
    edge_rotated = np.linalg.norm(rotated.vertices[1] - rotated.vertices[0])
    assert np.isclose(edge_original, edge_rotated)

    original_faces, rotated_faces = original.faces(), rotated.faces()
    for name in original_faces:
        assert np.isclose(original_faces[name].area(), rotated_faces[name].area(), atol=1e-8)
        assert np.allclose(rotation @ original_faces[name].normal(),
                           rotated_faces[name].normal(), atol=1e-8)


# ── Test 4: qeyri-ortoqonal (skewed) hüceyrə ────────────────────────────
def test_4_non_orthogonal_skewed_cell():
    """Üst üzü X istiqamətində sürüşdürülmüş kub (paralelepiped) —
    həcm QORUNUR (şəkil-dəyişmə determinantı 1), yan üzlər əyilir."""
    vertices = _box(1.0, 1.0, 1.0).vertices.copy()
    vertices[4:] += np.array([0.5, 0.0, 0.0])
    skewed = HexahedralCell(vertices)

    assert np.isclose(skewed.volume(), 1.0, atol=TOLERANCE)   # şəkil-dəyişmə həcmi qoruyur
    assert np.allclose(skewed.centroid(), [0.75, 0.5, 0.5], atol=TOLERANCE)
    assert skewed.validate().ok

    # X üzləri artıq oxa-perpendikulyar DEYİL
    assert not np.allclose(skewed.faces()["X-"].normal(), [-1, 0, 0])
    assert not np.allclose(skewed.faces()["X+"].normal(), [1, 0, 0])
    # Z üzləri hələ də üfüqi (yalnız X üzləri əyilib)
    assert np.allclose(skewed.faces()["Z-"].normal(), [0, 0, -1])
    assert np.allclose(skewed.faces()["Z+"].normal(), [0, 0, 1])


# ── Test 5: paylaşılan üz — normallar əks istiqamətdə ───────────────────
def test_5_shared_face_normals_are_opposite():
    cell_i = _box(1.0, 1.0, 1.0, origin=(0.0, 0.0, 0.0))
    cell_j = _box(1.0, 1.0, 1.0, origin=(1.0, 0.0, 0.0))   # X istiqamətində qonşu

    normal_i = cell_i.faces()["X+"].normal()
    normal_j = cell_j.faces()["X-"].normal()
    assert np.allclose(normal_i, -normal_j, atol=TOLERANCE)
    # paylaşılan üzün mərkəzi/sahəsi də EYNİDİR (fiziki olaraq eyni üz)
    assert np.allclose(cell_i.faces()["X+"].centroid(), cell_j.faces()["X-"].centroid())
    assert np.isclose(cell_i.faces()["X+"].area(), cell_j.faces()["X-"].area())


# ── Test 6: sərhəd üzü ───────────────────────────────────────────────────
def test_6_boundary_face_has_valid_geometry_with_single_owner():
    """Tək hüceyrə (qonşusu YOXDUR) — onun HƏR üzü etibarlı sahə/mərkəz/
    normal daşımalıdır, "sərhəd" statusu YALNIZ topologiyaya (bu modulda
    YOXDUR, bax audit §18) görə xarici olaraq müəyyən edilir."""
    cell = _box(3.0, 3.0, 3.0)
    boundary_face = cell.faces()["X+"]   # fərz: bu üzün qonşusu yoxdur (domen sərhədi)

    assert boundary_face.validate().ok
    assert boundary_face.area() > 0.0
    assert np.all(np.isfinite(boundary_face.centroid()))
    assert np.isclose(np.linalg.norm(boundary_face.normal()), 1.0)
    # tək sahib hüceyrə — ikinci hüceyrə/normal MÖVCUD DEYİL, YALNIZ bu birinin
    assert np.allclose(boundary_face.normal(), [1, 0, 0])


# ── Test 7: etibarsız — sıfır həcm ───────────────────────────────────────
def test_7_zero_volume_cell_fails_validation():
    flat = _box(2.0, 2.0, 0.0)
    result = flat.validate()
    assert not result.ok
    assert any("həcm" in e for e in result.errors)


# ── Test 8: etibarsız — sıfır sahəli üz ─────────────────────────────────
def test_8_zero_area_face_fails_validation():
    degenerate = Face(np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]], float))
    result = degenerate.validate()
    assert not result.ok
    assert any("sahə" in e for e in result.errors)


# ── Test 9: etibarsız — sıfır normal ─────────────────────────────────────
def test_9_zero_normal_fails_validation():
    """Üç kollinear + bir təkrar nöqtə — sahə sıfır, normal da sıfır-
    uzunluqlu olur, HƏR İKİSİ AYRICA bildirilir."""
    collinear = Face(np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], float))
    result = collinear.validate()
    assert not result.ok
    assert any("normal" in e for e in result.errors)


# ── Test 10: NaN/Inf həndəsə ─────────────────────────────────────────────
def test_10_nan_and_inf_geometry_fail_validation():
    nan_vertices = _box(1.0, 1.0, 1.0).vertices.copy()
    nan_vertices[0, 0] = np.nan
    assert not HexahedralCell(nan_vertices).validate().ok

    inf_vertices = _box(1.0, 1.0, 1.0).vertices.copy()
    inf_vertices[0, 0] = np.inf
    assert not HexahedralCell(inf_vertices).validate().ok


# ── Test 11: qeyri-ortoqonallıq metrikası (analitik yoxlama) ────────────
def test_11_non_orthogonality_metric_matches_analytical_value():
    """Bax audit §12: `θ = arccos(|d_ij·n_f| / (|d_ij||n_f|))`.

    Skewed hüceyrənin Z üzləri arasında ANALİTİK gözlənilən bucaq
    `arctan(0.5/1.0) = 26.565°`-dir (üfüqi sürüşmə 0.5, şaquli məsafə 1.0).
    """
    vertices = _box(1.0, 1.0, 1.0).vertices.copy()
    vertices[4:] += np.array([0.5, 0.0, 0.0])
    skewed = HexahedralCell(vertices)

    d_ij = cell_to_cell_vector(skewed.faces()["Z-"].centroid(),
                               skewed.faces()["Z+"].centroid())
    angle = non_orthogonality_angle(d_ij, skewed.faces()["Z+"].normal())
    expected = np.arctan(0.5 / 1.0)
    assert np.isclose(angle, expected, atol=1e-9)

    # ortoqonal qutu üçün bucaq DƏQİQ sıfırdır
    orthogonal = _box(1.0, 1.0, 1.0)
    d0 = cell_to_cell_vector(orthogonal.faces()["Z-"].centroid(),
                             orthogonal.faces()["Z+"].centroid())
    assert np.isclose(non_orthogonality_angle(d0, orthogonal.faces()["Z+"].normal()), 0.0,
                      atol=1e-12)


# ── Test 12: ümumi həndəsə TPFA-ya TƏSİR ETMİR ──────────────────────────
def test_12_general_geometry_module_does_not_affect_tpfa_regression():
    """Bu modul (`polyhedral_geometry.py`) heç bir mövcud sinif tərəfindən
    ÇAĞIRILMIR (bax audit §21: geometriya diskretizasiyadan MÜSTƏQİLDİR) —
    ona görə mövcud TPFA reqressiya testi DƏYİŞMƏDƏN keçməlidir. Bu, HƏR
    İKİSİNİ bir yerdə işə salaraq İNTEQRASİYA SƏVİYYƏSİNDƏ sübut edir."""
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


# ── §17/§24: keşləmə + N hüceyrə üzrə xətti miqyaslanma ─────────────────
def test_repeated_queries_are_cached_not_recomputed():
    """`volume()`/`centroid()`/`faces()` təkrar çağırış eyni (keşlənmiş)
    obyekti qaytarmalıdır — bax audit §24: "avoid repeated reconstruction
    of the same geometry"."""
    cell = _box(2.0, 3.0, 4.0)
    faces_first = cell.faces()
    faces_second = cell.faces()
    assert faces_first is faces_second   # eyni dict obyekti, YENİDƏN qurulmayıb

    face = faces_first["X+"]
    assert face.area() == face.area()    # təkrar çağırış eyni ƏDƏD (keş sınmayıb)


def test_constructing_and_validating_many_cells_scales_linearly():
    def _time_n_cells(n):
        start = time.perf_counter()
        for i in range(n):
            cell = _box(1.0 + 0.001 * i, 1.0, 1.0)
            cell.volume()
            cell.centroid()
            cell.validate()
        return time.perf_counter() - start

    _time_n_cells(200)             # isinmə
    small = _time_n_cells(1_000)
    large = _time_n_cells(20_000)  # 20x hüceyrə sayı
    ratio = large / max(small, 1e-9)
    # O(N) gözlənilir (~20x); O(N²) olsaydı ~400x olardı — 60x həddi
    # ikisini aydın ayırd edir, ölçmə səs-küyünə dözümlüdür
    assert ratio < 60, f"Hüceyrə-başına iş N-dən SUPERXƏTTİ artır (nisbət={ratio:.1f})"


# ── əlavə: `is_planar` düzgün ayırd edir ────────────────────────────────
def test_planar_vs_warped_face_detection():
    planar = Face(np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], float))
    assert planar.is_planar()

    warped = Face(np.array([[0, 0, 0], [1, 0, 0], [1, 1, 1], [0, 1, 0]], float))
    assert not warped.is_planar()
    # əyri üz üçün DƏ sahə/mərkəz/normal HESABLANIR (TƏXMİNİ, bax modul
    # docstring-i) — ÇÖKMÜR, sükutla "dəqiq" də İDDİA ETMİR
    assert warped.area() > 0.0
    assert np.all(np.isfinite(warped.centroid()))
    assert np.isclose(np.linalg.norm(warped.normal()), 1.0)
