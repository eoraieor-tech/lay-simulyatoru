"""FINDING-3 — ORİYENTASİYALI ÜZ SAHƏ VEKTORU (`oriented_area_vector`).

Bu fayl YALNIZ bir problemi əhatə edir: qeyri-müstəvi (warped, corner-
point) üzdə `Face.area() * Face.normal()` HƏQİQİ oriyentasiyalı sahə
vektoru DEYİL, ona görə divergensiya/qapanma və axın diskretizasiyası
üçün İSTİFADƏ EDİLƏ BİLMƏZ.

ÜÇ AYRI ANLAYIŞ (bax `polyhedral_geometry.py` docstring-i):

    area()                  SKALYAR   Σ ½‖(b−a)×(c−a)‖    (əyri səthin sahəsi)
    oriented_area_vector()  VEKTOR    Σ ½ (b−a)×(c−a)     (müstəvi proyeksiya)
    normal()                VAHİD     A⃗/‖A⃗‖

MÜSTƏVİ üzdə `area()·normal() == A⃗` (BİT-BƏRABƏR, aşağıda yoxlanılır);
ƏYRİ üzdə `‖A⃗‖ < area()` və YALNIZ `A⃗` ilə `Σ_f A⃗_f = 0` alınır.

§8-in tələb etdiyi 10 test kateqoriyası burada 1–10 kimi nömrələnib.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.discretization.mpfa_o_interaction import build_interaction_regions
from imex2d.domain.corner_point_geometry import quad_metrics
from imex2d.domain.general_grid_geometry import (GeneralGridGeometry,
                                                 hexahedral_vertices_from_cartesian)
from imex2d.domain.geometry import CellGeometry
from imex2d.domain.grid import CartesianGrid
from imex2d.domain.polyhedral_geometry import (HEX_FACE_VERTEX_INDICES, Face,
                                               HexahedralCell)

#: Vahid kubun təpələri, `HEX_FACE_VERTEX_INDICES` konvensiyasında.
UNIT_CUBE = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                      [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)

_OPP = {"X-": "X+", "X+": "X-", "Y-": "Y+", "Y+": "Y-", "Z-": "Z+", "Z+": "Z-"}


def _reference_area_vector(vertices: np.ndarray) -> np.ndarray:
    """TƏTBİQDƏN MÜSTƏQİL etalon: `A⃗ = Σ_tri ½·(b−a)×(c−a)`, mərkəz-fan
    üçbucaqları AÇIQ şəkildə burada qurulur (`Face`-in daxili
    `_triangles`-inə GÜVƏNMİR)."""
    vertices = np.asarray(vertices, float)
    apex = vertices.mean(axis=0)
    n = len(vertices)
    total = np.zeros(3)
    for i in range(n):
        a, b, c = apex, vertices[i], vertices[(i + 1) % n]
        total += 0.5 * np.cross(b - a, c - a)
    return total


def _warped_cell_vertices(amplitude: float = 0.5) -> np.ndarray:
    """Bir təpəsi z-də qaldırılmış kub — Z+ üzü QEYRİ-MÜSTƏVİ olur.
    Bu, FINDING-3-ün orijinal reproduksiya halıdır."""
    v = UNIT_CUBE.copy()
    v[6, 2] += amplitude
    return v


def _cell_scale(cell: HexahedralCell) -> float:
    return max(face.area() for face in cell.faces().values())


# ── 1: müstəvi Kartezian üz — köhnə və yeni ifadə BİT-BƏRABƏRDİR ──────────
def test_1_planar_cartesian_face_area_vector_is_bit_identical_to_old_expression():
    """REQRESSİYA MÜHAFİZƏSİ (§7): Kartezian üzdə `area()*normal()`
    HƏLƏ DƏ düzgündür, ona görə yeni ifadə ondan BİR BİT belə
    fərqlənməməlidir — `==`, `allclose` DEYİL."""
    grid = CartesianGrid(3, 2, 2)
    geometry = CellGeometry(grid, dx=7.0, dy=11.0, dz=[3.0, 5.0], top_depth=1234.0)
    vertices = hexahedral_vertices_from_cartesian(grid, geometry)

    for c in range(grid.ncell):
        for name, face in HexahedralCell(vertices[c]).faces().items():
            old = face.area() * face.normal()
            new = face.oriented_area_vector()
            assert np.array_equal(old, new), f"hüceyrə {c} üzü {name}: {old} != {new}"
            # analitik etalon da eynidir
            assert np.allclose(new, _reference_area_vector(face.vertices), atol=1e-12)
            assert face.is_planar(1e-12)
            assert np.isclose(np.linalg.norm(new), face.area(), rtol=0, atol=1e-12)


# ── 2: affin çevrilmiş hekzahedr — Nanson düsturu ─────────────────────────
@pytest.mark.parametrize("label,matrix", [
    ("identity", np.eye(3)),
    ("box", np.diag([2.0, 3.0, 5.0])),
    ("inclined prism", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.5, 0.0, 1.0]]),
    ("sheared prism", [[1.0, 0.3, 0.0], [0.0, 1.0, 0.4], [0.0, 0.0, 1.0]]),
    ("rotation+shear", [[0.8, -0.6, 0.2], [0.6, 0.8, 0.0], [0.0, 0.0, 1.3]]),
])
def test_2_affine_hexahedron_matches_nanson_formula(label, matrix):
    """Vahid kubun `S` affin təsviri üçün üz sahə-vektoru QAPALI
    ŞƏKİLDƏ məlumdur: `A⃗ = cof(S)·e = det(S)·S⁻ᵀ·e` (Nanson). Affin
    çevirmə müstəviliyi SAXLAYIR, ona görə burada `area()·normal()` də
    hələ düzgündür — yeni ifadə ONUNLA və ANALİTİK etalonla üst-üstə
    düşməlidir."""
    S = np.asarray(matrix, float)
    cell = HexahedralCell(UNIT_CUBE @ S.T)
    cofactor = np.linalg.det(S) * np.linalg.inv(S).T
    unit_normals = {"X+": [1, 0, 0], "X-": [-1, 0, 0], "Y+": [0, 1, 0],
                    "Y-": [0, -1, 0], "Z+": [0, 0, 1], "Z-": [0, 0, -1]}

    for name, face in cell.faces().items():
        expected = cofactor @ np.array(unit_normals[name], float)
        assert np.allclose(face.oriented_area_vector(), expected, atol=1e-12)
        assert np.allclose(face.area() * face.normal(), expected, atol=1e-12)
        assert np.isclose(face.area(), np.linalg.norm(expected), atol=1e-12)


# ── 3: sürüşdürülmüş (sheared) hüceyrə ────────────────────────────────────
def test_3_sheared_cell_area_vectors_are_exact_and_close():
    """Shear affin olduğu üçün BÜTÜN üzlər müstəvi qalır — sahə
    vektorları analitikdir və qapanma maşın dəqiqliyindədir."""
    S = np.array([[1.0, 0.0, 0.7], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    cell = HexahedralCell(UNIT_CUBE @ S.T)

    assert all(face.is_planar(1e-12) for face in cell.faces().values())
    for face in cell.faces().values():
        assert np.allclose(face.oriented_area_vector(),
                           _reference_area_vector(face.vertices), atol=1e-12)
        assert np.allclose(face.oriented_area_vector(),
                           face.area() * face.normal(), atol=1e-12)
    assert np.linalg.norm(cell.closure_residual()) <= 1e-14 * _cell_scale(cell)


# ── 4: warped / qeyri-müstəvi üz ──────────────────────────────────────────
def test_4_warped_face_area_vector_differs_from_scalar_area_times_normal():
    """ƏSAS HAL. Qeyri-müstəvi üzdə `‖A⃗‖ < area()` OLMALIDIR (skalyar
    sahə əyri səthin həqiqi sahəsidir, vektor isə onun proyeksiyası) və
    `A⃗` təpə sırasından qurulmuş ANALİTİK etalona bərabər olmalıdır."""
    cell = HexahedralCell(_warped_cell_vertices())
    top = cell.faces()["Z+"]

    assert not top.is_planar(1e-9), "Z+ üzü qeyri-müstəvi OLMALIDIR"
    area_vector = top.oriented_area_vector()
    assert np.allclose(area_vector, _reference_area_vector(top.vertices), atol=1e-14)

    # skalyar sahə STROQ olaraq vektorun uzunluğundan BÖYÜKDÜR
    assert np.linalg.norm(area_vector) < top.area()
    # və köhnə ifadə ilə fərq ÖLÇÜLƏ BİLƏNDİR (gizlədilmir)
    old = top.area() * top.normal()
    assert np.linalg.norm(old - area_vector) > 1e-3 * top.area()

    # müstəvi üzlərdə isə heç nə dəyişmir
    for name in ("X+", "X-", "Y+", "Y-", "Z-"):
        face = cell.faces()[name]
        assert face.is_planar(1e-12)
        assert np.allclose(face.oriented_area_vector(),
                           face.area() * face.normal(), atol=1e-14)


# ── 5: daxili üz — owner/neighbor ANTİSİMMETRİYASI ────────────────────────
@pytest.mark.parametrize("scenario", ["uniform", "skewed", "warped"])
def test_5_internal_face_area_vectors_are_antisymmetric(scenario):
    """`A⃗_owner = −A⃗_neighbor` — TƏK `Face` obyekti saxlandığı üçün bu
    KONSTRUKSİYA ilə təmin olunur, amma burada AÇIQ ölçülür. Həm də üz
    owner→neighbor istiqamətindədir."""
    ggg = _grid(scenario)
    interior = [gf for gf in ggg.faces if not gf.is_boundary]
    assert interior, f"{scenario} daxili üz yaratmadı"

    for gf in interior:
        owner_vector = gf.area_vector_from(gf.owner)
        neighbor_vector = gf.area_vector_from(gf.neighbor)
        assert np.array_equal(owner_vector, -neighbor_vector)
        assert np.allclose(owner_vector, ggg.face_area_vectors[gf.index], atol=0)
        # A⃗ · (c_neighbor − c_owner) > 0
        assert float(np.dot(owner_vector, ggg.d_ij(gf.index))) > 0.0
        # eyni fiziki üz qonşu hüceyrədən də EYNİ vektoru verir (əks işarə ilə)
        neighbor_face = ggg.cells[gf.neighbor].faces()[_OPP[gf.owner_local_name]]
        assert np.allclose(neighbor_face.oriented_area_vector(), -owner_vector,
                           atol=1e-12)


# ── 6: sərhəd üzü — owner-dan KƏNARA ──────────────────────────────────────
@pytest.mark.parametrize("scenario", ["uniform", "skewed", "warped"])
def test_6_boundary_face_area_vector_points_out_of_its_owner(scenario):
    ggg = _grid(scenario)
    boundary = [gf for gf in ggg.faces if gf.is_boundary]
    assert boundary

    for gf in boundary:
        area_vector = ggg.face_area_vectors[gf.index]
        outward = gf.face.centroid() - ggg.cell_centroids[gf.owner]
        assert float(np.dot(area_vector, outward)) > 0.0
        assert np.allclose(area_vector, gf.area_vector_from(gf.owner), atol=0)
        # vahid normal EYNİ istiqamətdədir
        assert float(np.dot(area_vector, ggg.face_normals[gf.index])) > 0.0
        with pytest.raises(ValueError):
            gf.area_vector_from(gf.owner + 10_000)


# ── 7: qapalı hüceyrə üçün Σ A⃗_f = 0 ──────────────────────────────────────
@pytest.mark.parametrize("label,vertices", [
    ("cube", UNIT_CUBE),
    ("sheared", UNIT_CUBE @ np.array([[1.0, 0.0, 0.7], [0, 1, 0], [0, 0, 1]]).T),
    ("warped-small", _warped_cell_vertices(0.1)),
    ("warped", _warped_cell_vertices(0.5)),
    ("warped-strong", _warped_cell_vertices(2.0)),
])
def test_7_closed_cell_area_vector_closure_is_machine_precision(label, vertices):
    """Divergensiya teoremi: QAPALI səth üçün `∮ n dS = 0`. Tolerans
    ƏSASLANDIRMASI — qalıq ~`max(A_f)` tərtibli 6 vektorun cəmidir, ona
    görə NİSBİ hədd götürülür; `1e-14` ikili toplamanın gözlənilən
    yığılmış xətasının (`~6·eps·A ≈ 1.3e-15·A`) cəmi bir neçə qatıdır,
    yəni KEÇMİŞ 1.15e-2 səhvini ötürməsi MÜMKÜN DEYİL."""
    cell = HexahedralCell(vertices)
    scale = _cell_scale(cell)

    residual = np.linalg.norm(cell.closure_residual())
    assert residual <= 1e-14 * scale, f"{label}: Σ A⃗ = {residual:.3e}"

    # etalon triangulyasiya ilə də eyni
    reference = sum(_reference_area_vector(f.vertices) for f in cell.faces().values())
    assert np.linalg.norm(reference) <= 1e-14 * scale


@pytest.mark.parametrize("scenario", ["uniform", "skewed", "warped"])
def test_7b_assembled_grid_closure_uses_shared_faces(scenario):
    """Grid səviyyəsində: hər hüceyrə üçün ÖZ üzlərinin owner-outward
    sahə vektorlarının cəmi sıfırdır (paylaşılan üzlərdə işarə çevrilir)."""
    ggg = _grid(scenario)
    scale = float(ggg.face_areas.max())
    for cell in range(ggg.ncell):
        total = np.zeros(3)
        for face_index in ggg.cell_faces(cell):
            total += ggg.faces[face_index].area_vector_from(cell)
        assert np.linalg.norm(total) <= 1e-14 * scale


# ── 8: MPFA-O inteqrasiyası ───────────────────────────────────────────────
@pytest.mark.parametrize("scenario", ["uniform", "skewed", "warped"])
def test_8_mpfa_o_sub_face_area_vectors_are_conservative(scenario):
    """MPFA-O sub-üzləri ana üzü TAM örtür, deməli onların sahə
    vektorlarının cəmi ana üzün sahə vektoruna BƏRABƏR olmalıdır — və
    hüceyrə üzrə cəm sıfır. `area()·normal()` ilə bu ƏYRİ üzdə
    POZULURDU (aşağıdakı test 10-a bax)."""
    ggg = _grid(scenario)
    grid = _grid_topology(scenario)
    regions = build_interaction_regions(grid, ggg, eta=1.0)
    scale = float(ggg.face_areas.max())

    per_face: dict = {}
    per_cell = np.zeros((ggg.ncell, 3))
    for region in regions:
        for sub in region.sub_faces:
            per_face.setdefault(sub.face_index, np.zeros(3))
            per_face[sub.face_index] += sub.area_vector
            for cell in sub.cells():
                per_cell[cell] += sub.outward_area_vector(cell)
            # işarə konvensiyası §6 — ana üzlə EYNİ tərəfə
            assert float(np.dot(sub.area_vector,
                                ggg.face_area_vectors[sub.face_index])) > 0.0
            # skalyar sahə ilə vektoru qarışdırmırıq
            assert np.linalg.norm(sub.area_vector) <= sub.area * (1.0 + 1e-12)

    for face_index, total in per_face.items():
        assert np.allclose(total, ggg.face_area_vectors[face_index], atol=1e-12), (
            f"sub-üzlər ana üz {face_index}-i bərpa etmir")
    assert np.abs(per_cell).max() <= 1e-12 * scale


def test_8b_mpfa_o_region_area_vector_sum_vanishes():
    """Daxili bölgə üçün `Σ_σ a⃗_σ` (hər hüceyrədən kənara) sıfırdır."""
    ggg = _grid("warped")
    regions = build_interaction_regions(_grid_topology("warped"), ggg, eta=1.0)
    scale = float(ggg.face_areas.max())
    interior = [r for r in regions if not r.is_boundary_region]
    assert interior, "daxili bölgə tapılmadı"
    for region in interior:
        total = np.zeros(3)
        for cell in region.cells:
            for local in region.cell_sub_faces[region.cell_local[cell]]:
                total += region.sub_faces[local].outward_area_vector(cell)
        assert np.linalg.norm(total) <= 1e-12 * scale


# ── 9: triangulyasiya uzlaşması (bir neçə MÜSTƏQİL yol) ───────────────────
@pytest.mark.parametrize("amplitude", [0.0, 0.1, 0.5, 2.0])
def test_9_area_vector_agrees_across_independent_implementations(amplitude):
    """Dörd MÜSTƏQİL yol EYNİ vektoru verməlidir:
      (a) `Face.oriented_area_vector()`;
      (b) bu faylın açıq etalon triangulyasiyası;
      (c) `corner_point_geometry.quad_metrics` (VEKTORLAŞDIRILMIŞ, ayrı kod);
      (d) dördbucaqlı üçün qapalı düstur `½·(C−A)×(D−B)`.
    """
    cell = HexahedralCell(_warped_cell_vertices(amplitude))
    quads = np.array([cell.vertices[list(idx)]
                      for idx in HEX_FACE_VERTEX_INDICES.values()])
    _, _, vectorised = quad_metrics(quads)

    for k, (name, face) in enumerate(cell.faces().items()):
        mine = face.oriented_area_vector()
        assert np.allclose(mine, _reference_area_vector(face.vertices), atol=1e-14)
        assert np.allclose(mine, vectorised[k], atol=1e-14)
        A, B, C, D = face.vertices
        assert np.allclose(mine, 0.5 * np.cross(C - A, D - B), atol=1e-14)


@pytest.mark.parametrize("shift", [0, 1, 2, 3])
def test_9b_area_vector_is_invariant_under_vertex_relabelling(shift):
    """Mərkəz-fan triangulyasiyası dövri sürüşməyə görə SİMMETRİKDİR
    (paylaşılan üzün owner və neighbor tərəfindən EYNİ görünməsinin
    səbəbi budur); sıranın TƏRSİNƏ çevrilməsi isə YALNIZ işarəni
    dəyişir."""
    face = HexahedralCell(_warped_cell_vertices()).faces()["Z+"]
    reference = face.oriented_area_vector()

    rotated = Face(np.roll(face.vertices, shift, axis=0))
    assert np.allclose(rotated.oriented_area_vector(), reference, atol=1e-14)

    reversed_face = Face(face.vertices[::-1])
    assert np.allclose(reversed_face.oriented_area_vector(), -reference, atol=1e-14)


# ── 10: KÖHNƏ davranışın reproduksiyası (1.15e-2) ─────────────────────────
def test_10_regression_old_scalar_convention_produced_percent_level_error():
    """ƏSAS TEST (§8): köhnə ifadə `area()·normal()` ilə qapanma
    NİSBİ səhvi PROSENT səviyyəsindədir; yeni `oriented_area_vector()`
    ilə maşın dəqiqliyindədir.

    Tolerans ARTIRILMIR — hər iki kəmiyyət EYNİ testdə, EYNİ həndəsədə
    ölçülür və aralarındakı ~13 tərtib fərq göstərilir."""
    cell = HexahedralCell(_warped_cell_vertices(0.5))
    scale = _cell_scale(cell)

    old = np.linalg.norm(sum(f.area() * f.normal() for f in cell.faces().values()))
    new = np.linalg.norm(cell.closure_residual())

    assert old / scale > 1e-3, f"köhnə səhv reproduksiya olunmadı ({old / scale:.3e})"
    assert new / scale <= 1e-14, f"yeni qapanma maşın dəqiqliyində DEYİL ({new / scale:.3e})"
    assert new < old * 1e-10


def test_10b_regression_grid_and_mpfa_level_error_is_removed():
    """Eyni müqayisə GRID və MPFA-O səviyyəsində — auditdə ölçülmüş
    1.15e-2 (grid) və 2.69e-3 (MPFA sub-üz) səhvləri."""
    ggg = _grid("warped")
    scale = float(ggg.face_areas.max())

    old_cell = np.zeros((ggg.ncell, 3))
    new_cell = np.zeros((ggg.ncell, 3))
    for cell in range(ggg.ncell):
        for face_index in ggg.cell_faces(cell):
            gf = ggg.faces[face_index]
            sign = 1.0 if gf.owner == cell else -1.0
            old_cell[cell] += sign * gf.face.area() * gf.face.normal()
            new_cell[cell] += gf.area_vector_from(cell)

    old_error = float(np.linalg.norm(old_cell, axis=1).max()) / scale
    new_error = float(np.linalg.norm(new_cell, axis=1).max()) / scale
    assert old_error > 1e-3, f"grid səviyyəsində köhnə səhv yoxdur ({old_error:.3e})"
    assert new_error <= 1e-14

    regions = build_interaction_regions(_grid_topology("warped"), ggg, eta=1.0)
    old_sub = np.zeros((ggg.ncell, 3))
    new_sub = np.zeros((ggg.ncell, 3))
    for region in regions:
        for sub in region.sub_faces:
            polygon = Face(sub.vertices)
            legacy = polygon.area() * polygon.normal()
            if float(np.dot(legacy, ggg.face_area_vectors[sub.face_index])) < 0.0:
                legacy = -legacy
            for cell in sub.cells():
                sign = 1.0 if cell == sub.owner else -1.0
                old_sub[cell] += sign * legacy
                new_sub[cell] += sub.outward_area_vector(cell)

    old_mpfa = float(np.linalg.norm(old_sub, axis=1).max()) / scale
    new_mpfa = float(np.linalg.norm(new_sub, axis=1).max()) / scale
    assert old_mpfa > 1e-4, f"MPFA səviyyəsində köhnə səhv yoxdur ({old_mpfa:.3e})"
    assert new_mpfa <= 1e-12


# ── ədədi möhkəmlik (§9) — dejenerativ hallar SƏSSİZ NaN vermir ───────────
@pytest.mark.parametrize("label,vertices", [
    ("collinear", [[0.0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]),
    ("repeated vertex", [[0.0, 0, 0], [0, 0, 0], [1, 1, 0], [0, 1, 0]]),
    ("all identical", [[1.0, 2, 3]] * 4),
    ("collapsed to a point", [[0.0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]),
])
def test_degenerate_faces_return_finite_vectors_not_nan(label, vertices):
    """Dejenerativ üz `NaN` DEYİL, SONLU vektor qaytarmalıdır; sıfır
    sahəli üz üçün `A⃗ = 0` yeganə düzgün cavabdır və `validate()` bunu
    AÇIQ xəta kimi bildirir (mövcud semantika POZULMUR)."""
    face = Face(np.asarray(vertices, float))
    area_vector = face.oriented_area_vector()

    assert np.all(np.isfinite(area_vector)), f"{label}: NaN/Inf sızdı"
    assert np.allclose(area_vector, _reference_area_vector(face.vertices), atol=1e-14)
    if face.area() <= 1e-12:
        assert np.allclose(area_vector, 0.0)
        assert not face.validate().ok          # SƏSSİZ qəbul edilmir


def test_repeated_vertex_wedge_keeps_closure():
    """Bir təpəsi qonşusuna çökmüş (paz/pinch-out) hekzahedr QANUNİDİR —
    qapanma orada da qorunmalıdır."""
    v = UNIT_CUBE.copy()
    v[1] = v[0]
    cell = HexahedralCell(v)
    assert cell.volume() > 0.0
    assert np.linalg.norm(cell.closure_residual()) <= 1e-14 * _cell_scale(cell)


def test_area_vector_does_not_mutate_the_decomposition_cache():
    """`oriented_area_vector()` KEŞLƏNMİŞ massivin ÖZÜNÜ qaytarmır —
    çağıran onu yerində dəyişsə, `normal()` korlanardı."""
    face = HexahedralCell(_warped_cell_vertices()).faces()["Z+"]
    first = face.oriented_area_vector()
    first += 1000.0
    assert np.allclose(face.oriented_area_vector(), _reference_area_vector(face.vertices))
    assert np.isclose(np.linalg.norm(face.normal()), 1.0)


# ── ssenari köməkçiləri ───────────────────────────────────────────────────
def _grid_topology(scenario: str) -> CartesianGrid:
    return CartesianGrid(3, 2, 2)


def _grid(scenario: str) -> GeneralGridGeometry:
    """`uniform` / `skewed` (affin, müstəvi üzlər) / `warped` (qeyri-
    müstəvi üzlər). Pozulma KOORDİNATIN funksiyasıdır, ona görə
    paylaşılan təpələr hər iki hüceyrədə EYNİ qalır."""
    grid = _grid_topology(scenario)
    geometry = CellGeometry(grid, dx=2.0, dy=3.0, dz=[1.0, 1.0], top_depth=0.0)
    conn = grid.build_connections()
    v = hexahedral_vertices_from_cartesian(grid, geometry)
    if scenario == "skewed":
        v[..., 0] += 0.4 * v[..., 2]
    elif scenario == "warped":
        v[..., 2] += (0.3 * np.sin(1.7 * v[..., 0]) * np.cos(1.1 * v[..., 1])
                      * (1.0 + 0.35 * v[..., 2]))
    return GeneralGridGeometry(v, conn)
