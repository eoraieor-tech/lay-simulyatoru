"""Phase 5D — HƏQİQİ corner-point (COORD/ZCORN) həndəsəsi.

Tapşırıq §1–§5: 8 təpəli rekonstruksiya, dəqiq çoxüzlü həcm, əyri üz
sahəsi, həqiqi normal, MPFA-O inteqrasiyası — və §3 GERİYƏ UYĞUNLUQ
(Kartezian model CPG-nin xüsusi halı kimi DƏYİŞMƏDƏN işləyir).
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

from imex2d.domain.corner_point_geometry import (CornerPointGeometry, cartesian_nodes,
                                                 corner_point_nodes, hex_metrics,
                                                 quad_metrics, unit_normals)
from imex2d.domain.general_grid_geometry import (GeneralGridGeometry,
                                                 hexahedral_vertices)
from imex2d.domain.geometry import CellGeometry
from imex2d.domain.grid import CartesianGrid
from imex2d.domain.polyhedral_geometry import HEX_FACE_VERTEX_INDICES, Face, HexahedralCell
from imex2d.io.grdecl import GrdeclError, read_grdecl
from imex2d.io.grdecl_import import GrdeclImporter
from imex2d.domain.diagnostics import DiagnosticReport

NX, NY, NZ = 4, 3, 2
DX, DY, DZ, TOP = 100.0, 80.0, 10.0, 2000.0


# ── deck qurucuları ──────────────────────────────────────────────────────
def _coord(nx=NX, ny=NY, dx=DX, dy=DY, top=TOP, base=3000.0,
           tilt_x=0.0, tilt_wave=0.0, y_sign=1.0) -> np.ndarray:
    """`(nx+1)·(ny+1)` pillar. `tilt_x` BÜTÜN pillarların dabanını X üzrə
    EYNİ qədər sürüşdürür (şaquli olmayan, amma paralel pillarlar);
    `tilt_wave` isə sürüşməni PİLLARDAN PİLLARA dəyişdirir — bu, üzləri
    həqiqətən qeyri-müstəvi edir (paralel pillarlarda üz müstəvi qalır,
    çünki x koordinatı z-in affin funksiyası olur). `y_sign=-1` sol-əlli
    deck yaradır."""
    values = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            x, y = i * dx, y_sign * j * dy
            shift_x = tilt_x + tilt_wave * np.sin(1.3 * i + 0.7 * j)
            shift_y = tilt_wave * np.cos(0.9 * i - 1.1 * j)
            values += [x, y, top, x + shift_x, y + shift_y, base]
    return np.array(values, float)


def _zcorn(nx=NX, ny=NY, nz=NZ, dz=DZ, top=TOP, dip_x=0.0, dip_y=0.0) -> np.ndarray:
    """`[k, üst/alt, j, j−/j+, i, i−/i+]`. `dip_x`/`dip_y` — lay meyli
    (hüceyrə başına dərinlik artımı)."""
    z = np.zeros((nz, 2, ny, 2, nx, 2))
    for k in range(nz):
        for plane in range(2):
            for j in range(ny):
                for cj in range(2):
                    for i in range(nx):
                        for ci in range(2):
                            column = top + dip_x * (i + ci) + dip_y * (j + cj)
                            z[k, plane, j, cj, i, ci] = column + (k + plane) * dz
    return z.ravel()


def _grid(nx=NX, ny=NY, nz=NZ) -> CartesianGrid:
    return CartesianGrid(nx, ny, nz)


def _box_geometry(grid=None) -> CellGeometry:
    grid = grid if grid is not None else _grid()
    return CellGeometry(grid=grid, dx=DX, dy=DY, dz=DZ, top_depth=TOP)


def _warped_nodes(seed=7, scale=4.0) -> np.ndarray:
    """Hər təpəsi MÜSTƏQİL sürüşdürülmüş hüceyrələr — üzlər əyridir, amma
    qonşular künc PAYLAŞMIR. YALNIZ "iki implementasiya eyni təpələrdə
    eyni ədədi verirmi" testləri üçün (mesh kimi mənası yoxdur)."""
    grid = _grid()
    nodes = cartesian_nodes(grid, _box_geometry(grid))
    rng = np.random.default_rng(seed)
    return nodes + rng.normal(0.0, scale, nodes.shape)


def _conformal_warped_zcorn(nx=NX, ny=NY, nz=NZ, dz=DZ, top=TOP, amplitude=6.0):
    """Künc dərinliyi QLOBAL künc indeksinin (`I=i+ci`, `J=j+cj`) hamar
    funksiyasıdır — ona görə qonşu hüceyrələr künc PAYLAŞIR (konformal
    mesh), amma üzlər müstəvi DEYİL."""
    z = np.zeros((nz, 2, ny, 2, nx, 2))
    for k in range(nz):
        for plane in range(2):
            for j in range(ny):
                for cj in range(2):
                    for i in range(nx):
                        for ci in range(2):
                            gi, gj = i + ci, j + cj
                            surface = top + amplitude * np.sin(1.7 * gi) * np.cos(1.1 * gj)
                            z[k, plane, j, cj, i, ci] = surface + (k + plane) * dz
    return z.ravel()


def _conformal_warped_geometry(tilt_wave=140.0):
    """Künclər PAYLAŞILAN (konformal), üzlər qeyri-müstəvi — həqiqi
    corner-point mesh.

    `tilt_wave=0` → pillarlar PARALEL (hamısı eyni qədər maili): X/Y
    üzləri müstəvi qalır (x, z-in affin funksiyası olur), amma əyri
    ZCORN səthinə görə Z üzləri qeyri-müstəvidir. Oblast bu halda
    XƏLƏNMİŞ PRİZMADIR, ona görə Kavalyeri prinsipi TAM tətbiq olunur.
    `tilt_wave>0` → pillar meyli pillardan pillara dəyişir, BÜTÜN üzlər
    qeyri-müstəvi olur (amma oblastın kənarı da əyrilir).
    """
    grid = _grid()
    geometry, _ = CornerPointGeometry.from_grdecl(
        grid, _coord(tilt_x=350.0, tilt_wave=tilt_wave), _conformal_warped_zcorn())
    return grid, geometry


# ═════════════════════════════════ §1 — 8 təpəli rekonstruksiya
def test_grdecl_box_reconstructs_exact_cartesian_cells():
    """Bərabər bloklu COORD/ZCORN → təpələr qutu təpələri ilə EYNİ."""
    grid = _grid()
    geometry, notes = CornerPointGeometry.from_grdecl(grid, _coord(), _zcorn())
    assert np.allclose(geometry.nodes, cartesian_nodes(grid, _box_geometry(grid)))
    assert notes == {"flipped_orientation": False, "degenerate_pillars": 0,
                     "negative_volume_cells": 0, "collapsed_cells": 0}


def test_node_ordering_follows_hex_face_convention():
    """`v0..v3` TAVAN (dayaz), `v4..v7` DABAN; ayaq izi
    `(i−,j−)→(i+,j−)→(i+,j+)→(i−,j+)`."""
    grid = _grid()
    geometry, _ = CornerPointGeometry.from_grdecl(grid, _coord(), _zcorn())
    first = geometry.nodes[0]
    assert np.allclose(first[0:4, 2], TOP)              # tavan müstəvisi
    assert np.allclose(first[4:8, 2], TOP + DZ)         # daban müstəvisi
    assert np.allclose(first[0], [0.0, 0.0, TOP])
    assert np.allclose(first[1], [DX, 0.0, TOP])
    assert np.allclose(first[2], [DX, DY, TOP])
    assert np.allclose(first[3], [0.0, DY, TOP])


def test_cell_index_matches_cartesian_grid_ordering():
    grid = _grid()
    geometry, _ = CornerPointGeometry.from_grdecl(grid, _coord(), _zcorn())
    for (i, j, k) in ((0, 0, 0), (3, 2, 1), (2, 1, 0)):
        cell = grid.index(i, j, k)
        nodes = geometry.cell_nodes(cell)
        assert np.isclose(nodes[0, 0], i * DX)
        assert np.isclose(nodes[0, 1], j * DY)
        assert np.isclose(nodes[0, 2], TOP + k * DZ)


def test_non_vertical_pillars_interpolate_x_along_the_pillar():
    """Pillar dabanı X üzrə sürüşübsə, künc X-i DƏRİNLİKLƏ dəyişməlidir."""
    grid = _grid()
    tilt = 500.0
    geometry, _ = CornerPointGeometry.from_grdecl(
        grid, _coord(tilt_x=tilt), _zcorn())
    nodes = geometry.cell_nodes(0)
    # t = (z − top) / (base − top); base−top = 1000, dz = 10 → t = 0.01
    assert np.isclose(nodes[0, 0], 0.0)                       # tavan: t = 0
    assert np.isclose(nodes[4, 0], tilt * (DZ / (3000.0 - TOP)))
    assert nodes[4, 0] > nodes[0, 0]


def test_degenerate_pillar_is_reported_not_silently_nan():
    grid = _grid()
    coord = _coord()
    coord[5] = coord[2]                                # ilk pillar: daban = tavan
    nodes, notes = corner_point_nodes(NX, NY, NZ, coord, _zcorn())
    assert notes["degenerate_pillars"] == 1
    assert np.all(np.isfinite(nodes))


def test_size_mismatch_is_an_explicit_error():
    with pytest.raises(ValueError, match="COORD ölçüsü"):
        corner_point_nodes(NX, NY, NZ, _coord()[:-6], _zcorn())
    with pytest.raises(ValueError, match="ZCORN ölçüsü"):
        corner_point_nodes(NX, NY, NZ, _coord(), _zcorn()[:-8])


# ═════════════════════════════════ §2 — dəqiq çoxüzlü həcm
def test_volume_is_exact_for_box_cells():
    grid = _grid()
    geometry, _ = CornerPointGeometry.from_grdecl(grid, _coord(), _zcorn())
    assert np.allclose(geometry.volumes(), DX * DY * DZ)


def test_shear_preserves_volume_but_not_the_box_formula():
    """Maili lay = XƏLƏMƏ (shear): həcm Kavalyeri prinsipinə görə
    DƏYİŞMİR, amma normal/mərkəz dəyişir."""
    grid = _grid()
    geometry, _ = CornerPointGeometry.from_grdecl(
        grid, _coord(), _zcorn(dip_x=30.0))
    assert np.allclose(geometry.volumes(), DX * DY * DZ)
    assert not np.allclose(geometry.cell_centroid()[:, 2],
                           _box_geometry(grid).cell_depths())


def test_vectorized_volume_matches_scalar_hexahedral_cell():
    """`hex_metrics` == `HexahedralCell.volume()/.centroid()` — ƏYRİ
    hüceyrələrdə də (iki kod yolu ayrılmasın)."""
    nodes = _warped_nodes()
    volumes, centroids = hex_metrics(nodes)
    for cell in range(nodes.shape[0]):
        reference = HexahedralCell(nodes[cell])
        assert np.isclose(volumes[cell], reference.volume(), rtol=1e-12)
        assert np.allclose(centroids[cell], reference.centroid(), rtol=1e-12)


def test_tetrahedral_decomposition_uses_cell_and_face_centres():
    """Tapşırıq §2 düsturunun BİRBAŞA yoxlanışı: `P_c = ⅛Σp_i` və üz
    mərkəzləri ilə 24 tetraedrin cəmi."""
    nodes = _warped_nodes()[3]
    apex = nodes.mean(axis=0)
    total = 0.0
    for indices in HEX_FACE_VERTEX_INDICES.values():
        quad = nodes[list(indices)]
        centre = quad.mean(axis=0)
        for corner in range(4):
            a, b = quad[corner], quad[(corner + 1) % 4]
            total += np.dot(centre - apex, np.cross(a - apex, b - apex)) / 6.0
    assert np.isclose(total, hex_metrics(nodes[None])[0][0], rtol=1e-12)


def test_cells_tile_the_domain_without_gaps_under_shear():
    """Həcmlərin CƏMİ bütöv blokun həcminə bərabərdir — hüceyrələr
    boşluqsuz döşənir (konservativlik ön şərti)."""
    grid = _grid()
    geometry, _ = CornerPointGeometry.from_grdecl(
        grid, _coord(tilt_x=400.0), _zcorn(dip_x=20.0, dip_y=15.0))
    assert np.isclose(geometry.volumes().sum(), NX * NY * NZ * DX * DY * DZ)


# ═════════════════════════════════ §3 — qeyri-müstəvi üz sahəsi
def test_vectorized_face_metrics_match_scalar_face():
    nodes = _warped_nodes()
    quads = nodes[:, list(HEX_FACE_VERTEX_INDICES["X+"]), :]
    area, centroid, area_vector = quad_metrics(quads)
    normal = unit_normals(area_vector)
    for cell in range(quads.shape[0]):
        reference = Face(quads[cell])
        assert np.isclose(area[cell], reference.area(), rtol=1e-12)
        assert np.allclose(centroid[cell], reference.centroid(), rtol=1e-12)
        assert np.allclose(normal[cell], reference.normal(), rtol=1e-12)


def test_warped_face_area_exceeds_its_vector_area():
    """ƏYRİ üzdə skalyar sahə (həqiqi səth) > ‖sahə-vektoru‖; MÜSTƏVİ
    üzdə ikisi ÜST-ÜSTƏ DÜŞÜR."""
    warped = np.array([[[0, 0, 0], [10, 0, 0], [10, 10, 6], [0, 10, 0]]], float)
    area, _, area_vector = quad_metrics(warped)
    assert area[0] > np.linalg.norm(area_vector[0]) + 1e-9

    planar = np.array([[[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]], float)
    area, _, area_vector = quad_metrics(planar)
    assert np.isclose(area[0], np.linalg.norm(area_vector[0]))
    assert np.isclose(area[0], 100.0)


def test_face_area_is_not_the_bounding_box_area_on_skewed_cells():
    grid = _grid()
    conn = grid.build_connections()
    geometry, _ = CornerPointGeometry.from_grdecl(
        grid, _coord(tilt_x=600.0), _zcorn())
    areas = geometry.face_areas(conn)
    box = _box_geometry(grid).face_areas(conn)
    assert not np.allclose(areas, box)
    assert np.all(areas > 0.0)


def test_shared_face_area_is_identical_from_both_cells():
    """ƏYRİ üzdə owner və neighbor EYNİ sahəni/mərkəzi/əks normalı
    görməlidir — `Face` mərkəz-fan parçalanmasının reqressiya testi
    (`v0`-dan sadə fan burada İKİ FƏRQLİ diaqonal seçib fərqli ədəd
    verərdi, yəni qonşu hüceyrələr arasında boşluq/örtüşmə yaranardı)."""
    grid, geometry = _conformal_warped_geometry()
    # Ön şərt: üzlər HƏQİQƏTƏN qeyri-müstəvidir, yoxsa test heç nə sınamır.
    quad = geometry.nodes[0][list(HEX_FACE_VERTEX_INDICES["X+"])]
    assert not Face(quad).is_planar()

    general = GeneralGridGeometry(geometry.nodes, grid.build_connections())
    result = general.validate()
    assert result.ok, result.errors


def test_warped_cells_still_tile_without_gaps():
    """Qeyri-müstəvi Z üzləri olan konformal meshdə həcmlərin CƏMİ bütöv
    oblastın həcminə bərabər qalır — hüceyrələr boşluqsuz/örtüşməsiz
    döşənir (mərkəz-fan parçalanmasının konservativlik nəticəsi).

    Paralel pillarlarla (`tilt_wave=0`) oblast XƏLƏNMİŞ PRİZMADIR, ona
    görə Kavalyeri prinsipi ilə həcm `ayaq izi × ümumi qalınlıq`-dır."""
    _, geometry = _conformal_warped_geometry(tilt_wave=0.0)
    quad = geometry.nodes[0][list(HEX_FACE_VERTEX_INDICES["Z+"])]
    assert not Face(quad).is_planar()               # ön şərt: mesh əyridir

    expected = (NX * DX) * (NY * DY) * (NZ * DZ)
    assert np.isclose(geometry.volumes().sum(), expected, rtol=1e-9)


# ═════════════════════════════════ §4 — həqiqi normal vektorları
def test_normal_equals_the_diagonal_cross_product_formula():
    """Tapşırıq §4: `n̂ = (C−A)×(D−B) / ‖·‖` — mərkəz-fan sahə-vektoru
    məhz bu düstura BƏRABƏRDİR (eynilik)."""
    quads = _warped_nodes()[:, list(HEX_FACE_VERTEX_INDICES["Z+"]), :]
    _, _, area_vector = quad_metrics(quads)
    a, b, c, d = quads[:, 0], quads[:, 1], quads[:, 2], quads[:, 3]
    assert np.allclose(area_vector, 0.5 * np.cross(c - a, d - b))


def test_normals_tilt_with_a_dipping_layer():
    """Maili layda Z-üzünün normalı `[0,0,1]` OLMAMALIDIR."""
    grid = _grid()
    conn = grid.build_connections()
    dip = 25.0
    geometry, _ = CornerPointGeometry.from_grdecl(grid, _coord(), _zcorn(dip_x=dip))
    normals = geometry.face_normal(conn)[conn.axis == 2]
    expected = np.array([-dip, 0.0, DX])
    expected = expected / np.linalg.norm(expected)
    assert np.allclose(normals, expected)
    assert not np.allclose(normals, [0.0, 0.0, 1.0])


def test_normals_point_from_cell_a_to_cell_b():
    grid = _grid()
    conn = grid.build_connections()
    geometry, _ = CornerPointGeometry.from_grdecl(
        grid, _coord(tilt_x=300.0), _zcorn(dip_x=20.0))
    normal = geometry.face_normal(conn)
    centroids = geometry.cell_centroid()
    d_ij = centroids[conn.cell_b] - centroids[conn.cell_a]
    assert np.all(np.einsum("fj,fj->f", d_ij, normal) > 0.0)


def test_face_area_vectors_close_each_cell():
    """Qapalı hüceyrə üçün SAHƏ-VEKTORLARININ cəmi sıfırdır: `Σ S⃗_f = 0`
    (divergensiya teoremi).

    DİQQƏT — bu, `HexahedralCell.closure_residual()`-un hesabladığı
    `Σ A_f·n̂_f`-dən FƏRQLİDİR: ƏYRİ üzdə skalyar sahə `A_f` sahə-vektorun
    uzunluğundan BÖYÜKDÜR (əyri səth müstəvi proyeksiyadan geniş olur),
    ona görə həqiqi eynilik VEKTOR cəmi üzərindədir. Konservativlik məhz
    bu eynilikdən asılıdır."""
    nodes = _warped_nodes()                       # ən sərt hal: hər üz əyri
    for cell in range(nodes.shape[0]):
        quads = np.array([nodes[cell][list(idx)]
                          for idx in HEX_FACE_VERTEX_INDICES.values()])
        _, _, area_vectors = quad_metrics(quads)
        scale = np.abs(hex_metrics(nodes[cell][None])[0][0]) ** (2.0 / 3.0)
        assert np.linalg.norm(area_vectors.sum(axis=0)) < 1e-9 * scale


def test_unit_normals_have_unit_length():
    quads = _warped_nodes()[:, list(HEX_FACE_VERTEX_INDICES["Y+"]), :]
    _, _, area_vector = quad_metrics(quads)
    assert np.allclose(np.linalg.norm(unit_normals(area_vector), axis=1), 1.0)


def test_degenerate_face_normal_is_zero_not_silently_normalised():
    assert np.allclose(unit_normals(np.zeros((3, 3))), 0.0)


# ═════════════════════════════════ §5 — MPFA-O inteqrasiyası
def test_hexahedral_vertices_passes_corner_point_nodes_through_unchanged():
    grid = _grid()
    geometry, _ = CornerPointGeometry.from_grdecl(
        grid, _coord(tilt_x=250.0), _zcorn(dip_x=18.0))
    assert hexahedral_vertices(grid, geometry) is geometry.nodes


def test_hexahedral_vertices_builds_box_nodes_for_plain_cell_geometry():
    grid = _grid()
    geometry = _box_geometry(grid)
    assert np.allclose(hexahedral_vertices(grid, geometry),
                       cartesian_nodes(grid, geometry))


def _reservoir_model_from_deck(text: str):
    """GRDECL → `GeologicalModel` → `ReservoirModel` — HƏQİQİ zəncir,
    saxta obyekt YOX."""
    from helpers import default_scal
    from imex2d.application.model_builder import ReservoirModelBuilder
    from imex2d.application.scenarios import five_spot

    geology, _ = _import(text)
    return ReservoirModelBuilder().build(
        geological_model=geology, wells=five_spot(geology.grid),
        scal=default_scal(), name="CPG test")


def test_mpfa_o_consumes_true_corner_point_geometry():
    """MPFA-O `build()` artıq həndəsəni Kartezian FƏRZ ETMİR: qeyri-
    ortoqonal modeldə qurduğu üz normalları ox vektorları DEYİL."""
    from imex2d.discretization import MPFAODiscretization
    from imex2d.discretization.mpfa_o_local_system import MPFAOBoundaryClosure

    model = _reservoir_model_from_deck(_corner_point_deck(nx=4, ny=4, nz=2, dip=20.0))
    assert isinstance(model.geometry, CornerPointGeometry)

    discretized = MPFAODiscretization(
        closure=MPFAOBoundaryClosure.NEUMANN_ZERO).build(model)

    # Təpələr həndəsədən BİRBAŞA gəlir — Kartezianlaşdırılmayıb.
    assert np.array_equal(
        np.array([cell.vertices for cell in discretized.geometry.cells]),
        model.geometry.nodes)
    normals = discretized.geometry.face_normals
    assert not np.all(np.abs(normals).max(axis=1) > 1.0 - 1e-9)   # əyilmiş normallar
    assert np.allclose(discretized.geometry.cell_volumes, model.geometry.volumes())
    assert discretized.geometry.validate().ok


def test_tpfa_uses_true_corner_point_areas_and_volumes():
    """TPFA də (dəyişmədən, EYNİ `CellGeometry` müqaviləsi ilə) həqiqi
    həndəsəni görür — `dx·dy·dz` qutusunu YOX."""
    from imex2d.simulation.discretization import TwoPointFluxDiscretization

    model = _reservoir_model_from_deck(_corner_point_deck(nx=4, ny=4, nz=2, dip=20.0))
    discretized = TwoPointFluxDiscretization().build(model)
    assert np.allclose(discretized.cell_volume, model.geometry.volumes())


def test_mpfa_o_on_cartesian_model_is_unchanged_by_the_dispatcher():
    """Kartezian model — köhnə yol (`hexahedral_vertices_from_cartesian`)
    ilə BİT-BƏ-BİT eyni təpələr."""
    from imex2d.domain.general_grid_geometry import hexahedral_vertices_from_cartesian

    grid = _grid()
    geometry = _box_geometry(grid)
    assert np.array_equal(hexahedral_vertices(grid, geometry),
                          hexahedral_vertices_from_cartesian(grid, geometry))


# ═════════════════════════════════ §3 (tapşırıq) — geriyə uyğunluq
@pytest.mark.parametrize("dz,top_map", [
    (DZ, None),
    ([8.0, 14.0], None),
    ([8.0, 14.0], np.linspace(1990.0, 2010.0, NX * NY)),
])
def test_cartesian_is_an_exact_special_case_of_corner_point(dz, top_map):
    """Tapşırıq §3 — Kartezian modellər "seamless" işləyir: HƏR
    kəmiyyət qutu düsturu ilə maşın dəqiqliyində üst-üstə düşür."""
    grid = _grid()
    box = CellGeometry(grid=grid, dx=DX, dy=DY, dz=dz, top_depth=TOP,
                       top_depth_map=top_map)
    cpg = CornerPointGeometry.from_cartesian(box)
    conn = grid.build_connections()

    assert np.allclose(cpg.volumes(), box.volumes())
    assert np.allclose(cpg.cell_depths(), box.cell_depths())
    assert np.allclose(cpg.cell_centroid(), box.cell_centroid())
    assert np.allclose(cpg.dz_per_cell(), box.dz_per_cell())
    assert np.allclose(cpg.face_areas(conn), box.face_areas(conn))
    assert np.allclose(cpg.face_normal(conn), box.face_normal(conn))
    assert np.allclose(cpg.face_centroid(conn), box.face_centroid(conn))
    half_a, half_b = cpg.face_half_distances(conn)
    box_a, box_b = box.face_half_distances(conn)
    assert np.allclose(half_a, box_a) and np.allclose(half_b, box_b)
    assert cpg.areal_extent() == box.areal_extent()
    assert cpg.validate() == []


def test_corner_point_geometry_is_a_cell_geometry():
    """Model zənciri (`GeologicalModel` → `ReservoirModel` → TPFA) tip
    yoxlaması etmədən işləsin deyə."""
    grid = _grid()
    geometry, _ = CornerPointGeometry.from_grdecl(grid, _coord(), _zcorn())
    assert isinstance(geometry, CellGeometry)


def test_from_cartesian_is_idempotent():
    grid = _grid()
    cpg = CornerPointGeometry.from_cartesian(_box_geometry(grid))
    assert CornerPointGeometry.from_cartesian(cpg) is cpg


def test_nodes_are_required():
    with pytest.raises(ValueError, match="nodes"):
        CornerPointGeometry(grid=_grid(), dx=DX, dy=DY, dz=DZ)


def test_nominal_scalars_are_derived_but_never_used_for_volume():
    """`dx`/`dy`/`dz` NOMİNALDIR — `xy_to_ij`/`depth_to_k` üçün qalır,
    amma həcmi TƏYİN ETMİR."""
    grid = _grid()
    geometry, _ = CornerPointGeometry.from_grdecl(
        grid, _coord(), _zcorn(dip_x=40.0))
    assert np.isclose(geometry.dx, DX) and np.isclose(geometry.dy, DY)
    assert len(geometry.approximation_notes()) == 3
    # Meylli layda hüceyrə mərkəzi ARTIQ `top + (k+½)dz` deyil.
    nominal = TOP + 0.5 * DZ
    assert not np.allclose(geometry.cell_depths(), nominal)


# ═════════════════════════════════ yönüm / dejenerativ hallar
def test_left_handed_deck_is_flipped_once_globally_and_reported():
    grid = _grid()
    geometry, notes = CornerPointGeometry.from_grdecl(
        grid, _coord(y_sign=-1.0), _zcorn())
    assert notes["flipped_orientation"] is True
    assert np.all(geometry.volumes() > 0.0)
    assert np.allclose(geometry.volumes(), DX * DY * DZ)


def test_orientation_fix_can_be_disabled():
    nodes, notes = corner_point_nodes(NX, NY, NZ, _coord(y_sign=-1.0), _zcorn(),
                                      fix_orientation=False)
    volumes, _ = hex_metrics(nodes)
    assert notes["flipped_orientation"] is False
    assert np.all(volumes < 0.0)
    assert notes["negative_volume_cells"] == volumes.size


def test_pinched_out_cells_are_counted_not_hidden():
    """ZCORN-da tavan == daban → sıfır həcmli lay (Eclipse pinch-out)."""
    z = _zcorn().reshape(NZ, 2, NY, 2, NX, 2)
    z[0, 1] = z[0, 0]                                  # 0-cı layın dabanı = tavanı
    _, notes = corner_point_nodes(NX, NY, NZ, _coord(), z.ravel())
    assert notes["collapsed_cells"] == NX * NY


def test_validate_flags_degenerate_geometry():
    grid = _grid()
    geometry, _ = CornerPointGeometry.from_grdecl(grid, _coord(), _zcorn())
    assert geometry.validate() == []

    collapsed = geometry.nodes.copy()
    collapsed[0, 4:8, 2] = collapsed[0, 0:4, 2]        # ilk hüceyrəni yastıla
    flat = CornerPointGeometry.from_nodes(grid, collapsed)
    assert flat.validate()                             # xəta bildirilir


# ═════════════════════════════════ diaqnostika
def test_non_orthogonality_is_zero_for_a_box_and_positive_when_skewed():
    grid = _grid()
    conn = grid.build_connections()
    box, _ = CornerPointGeometry.from_grdecl(grid, _coord(), _zcorn())
    assert box.quality_metrics(conn)["max_non_orthogonality_angle_deg"] < 1e-9

    dip = 25.0
    skewed, _ = CornerPointGeometry.from_grdecl(grid, _coord(), _zcorn(dip_x=dip))
    metrics = skewed.quality_metrics(conn)
    assert metrics["max_non_orthogonality_angle_deg"] > 1.0
    assert metrics["min_cell_volume"] > 0.0


# ═════════════════════════════════ GRDECL idxalı (uçdan-uca)
def _corner_point_deck(nx=3, ny=2, nz=2, dip=20.0) -> str:
    n = nx * ny * nz
    coord = _coord(nx=nx, ny=ny, dx=50.0, dy=40.0)
    zcorn = _zcorn(nx=nx, ny=ny, nz=nz, dz=10.0, dip_x=dip)
    return f"""-- corner-point nümunəsi
SPECGRID
  {nx} {ny} {nz} 1 F /

COORD
  {' '.join(f'{v:.4f}' for v in coord)} /

ZCORN
  {' '.join(f'{v:.4f}' for v in zcorn)} /

PORO
  {n}*0.2 /
PERMX
  {n}*150 /
"""


def _import(text: str):
    handle, path = tempfile.mkstemp(suffix=".GRDECL")
    with os.fdopen(handle, "w", encoding="utf-8") as file:
        file.write(text)
    try:
        report = DiagnosticReport()
        model = GrdeclImporter().build(read_grdecl(path), report)
        return model, report
    finally:
        os.unlink(path)


def test_grdecl_import_produces_true_corner_point_geometry():
    model, report = _import(_corner_point_deck())
    assert isinstance(model.geometry, CornerPointGeometry)
    assert model.geometry.nodes.shape == (model.grid.ncell, 8, 3)
    # Maili lay → normallar ox vektoru DEYİL.
    conn = model.grid.build_connections()
    normals = model.geometry.face_normal(conn)[conn.axis == 2]
    assert not np.allclose(normals, [0.0, 0.0, 1.0])


def test_grdecl_import_no_longer_warns_about_uniform_block_approximation():
    _, report = _import(_corner_point_deck())
    text = " ".join(str(entry) for entry in report.items)
    assert "APPROKSİMASİYA YOXDUR" in text
    assert "approksimasiya olundu" not in text


def test_grdecl_import_reports_nominal_scalar_limits_honestly():
    _, report = _import(_corner_point_deck())
    text = " ".join(str(entry) for entry in report.items)
    assert "nominal" in text and "xy_to_ij" in text


def test_grdecl_import_rejects_mismatched_corner_point_arrays():
    """COORD-dan bir pillar (6 dəyər) çıxarılır — grid ölçüləri və bütün
    xassə massivləri DÜZGÜN qalır, ona görə xəta məhz COORD/ZCORN
    uyğunsuzluğuna görə verilməlidir."""
    head, rest = _corner_point_deck().split("COORD\n  ", 1)
    values, tail = rest.split(" /", 1)
    truncated = " ".join(values.split()[:-6])
    with pytest.raises(GrdeclError, match="COORD/ZCORN"):
        _import(f"{head}COORD\n  {truncated} /{tail}")


def test_block_centred_deck_still_builds_plain_cell_geometry():
    """COORD/ZCORN OLMAYAN deck köhnə yolu saxlayır (geriyə uyğunluq)."""
    n = 3 * 2 * 2
    deck = f"""SPECGRID
  3 2 2 1 F /
DX
  {n}*50 /
DY
  {n}*40 /
DZ
  {n}*10 /
TOPS
  6*2000 /
PORO
  {n}*0.2 /
PERMX
  {n}*150 /
"""
    model, _ = _import(deck)
    assert type(model.geometry) is CellGeometry
    assert np.allclose(model.geometry.volumes(), 50.0 * 40.0 * 10.0)


# ═════════════════════════════════ serializasiya
def test_corner_point_nodes_survive_a_serialization_round_trip():
    """Layihəni saxlayıb açmaq HƏQİQİ həndəsəni İTİRMƏMƏLİDİR — `nodes`
    yazılmasaydı model sükutla bərabər bloka çevrilərdi."""
    from imex2d.application.serialization import ProjectSerializer

    geology, _ = _import(_corner_point_deck(nx=3, ny=2, nz=2, dip=20.0))
    serializer = ProjectSerializer()
    restored = serializer.geological_model_from_dict(
        serializer.geological_model_to_dict(geology))

    assert isinstance(restored.geometry, CornerPointGeometry)
    assert np.allclose(restored.geometry.nodes, geology.geometry.nodes)
    assert np.allclose(restored.geometry.volumes(), geology.geometry.volumes())


def test_plain_cartesian_geometry_round_trip_is_unchanged():
    """`nodes` OLMAYAN (köhnə) qeyd əvvəlki kimi `CellGeometry` qaytarır."""
    from imex2d.application.serialization import ProjectSerializer

    grid = _grid()
    payload = ProjectSerializer()._geometry_to_dict(_box_geometry(grid))
    assert "nodes" not in payload
    restored = ProjectSerializer()._geometry_from_dict(grid, payload)
    assert type(restored) is CellGeometry
