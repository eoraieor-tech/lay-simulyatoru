"""Peaceman quyu indeksi — HƏQİQİ corner-point wellblock həndəsəsi.

    Corner-Point Geometry
            ↓
    real local wellblock geometry   (V/h en kəsiyi + yerli kənar oxları)
            ↓
    effective Peaceman radius / dimensions
            ↓
    WI

Əvvəl Peaceman `geometry.cell_extents()` oxuyurdu — o isə OX-BOYU SƏRHƏD
QUTUSUDUR (`max(x)−min(x)`, `max(y)−min(y)`). Fırlanmış hüceyrədə qutu
hüceyrənin öz eninə deyil, DİAQONALINA yaxınlaşır; kəsilmiş (skewed)
hüceyrədə isə həm şişir, həm də həqiqi ayaq izi sahəsini itirir. Bu
modul həmin fərqi hər model tipi üçün ayrıca qıfıllayır.

Kartezian modelin nəticəsi maşın dəqiqliyində DƏYİŞMƏMƏLİDİR — bu, ən
sərt tələbdir və ilk testlərdə yoxlanılır.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.domain.corner_point_geometry import CornerPointGeometry
from imex2d.domain.geometry import CellGeometry
from imex2d.domain.properties import PropertyMap, RockProperties
from imex2d.domain.reservoir_model import ReservoirModel
from imex2d.domain.wells import (ControlMode, Perforation, Well, WellControl,
                                 WellType)
from imex2d.simulation.well_model import PeacemanWellModel
from imex2d.domain.units import METRIC

from test_local_cell_metrics import DX, DY, NX, NY, NZ, TOP, _cpg, _grid

DZ_LAYERS = [6.0, 21.0]
RADIUS = 0.1
WELL_I, WELL_J = 1, 1


# ══════════════════════════════════════════════════ köməkçilər ══════════
def _rock(ncell, kx=150.0, ky=150.0, kz=15.0) -> RockProperties:
    return RockProperties(
        porosity=PropertyMap.uniform("PORO", 0.24, ncell),
        permx=PropertyMap.uniform("PERMX", kx, ncell),
        permy=PropertyMap.uniform("PERMY", ky, ncell),
        permz=PropertyMap.uniform("PERMZ", kz, ncell))


def _connections(geometry, kx=150.0, ky=150.0, kz=15.0, direction="Z",
                 i=WELL_I, j=WELL_J, skin=0.0):
    grid = geometry.grid
    well = Well(name="P1", well_type=WellType.PRODUCER,
                control=WellControl(ControlMode.BHP, 150.0),
                perforations=[Perforation(i, j, k, True, skin, direction)
                              for k in range(grid.nz)],
                radius=RADIUS)
    model = ReservoirModel(name="peaceman-cpg", grid=grid, geometry=geometry,
                           rock=_rock(grid.ncell, kx, ky, kz), wells=[well])
    return PeacemanWellModel().build_connections(model)


def _well_index(geometry, **kwargs) -> np.ndarray:
    return np.array([c.well_index for c in _connections(geometry, **kwargs)])


def _peaceman(k1, k2, d1, d2, h, rw=RADIUS, skin=0.0) -> float:
    """Klassik anizotrop Peaceman — TESTİN ÖZ, müstəqil realizasiyası.

    Düstur `well_model.py`-dən KÖÇÜRÜLMÜR (o dəyişsəydi test də
    səssizcə dəyişərdi); burada mənbədən asılı olmayan istinad kimi
    yenidən yazılıb."""
    ratio = k2 / k1
    re = (0.28 * np.sqrt(np.sqrt(ratio) * d1 ** 2 + np.sqrt(1.0 / ratio) * d2 ** 2)
          / (ratio ** 0.25 + (1.0 / ratio) ** 0.25))
    return (METRIC.darcy_constant * 2.0 * np.pi * np.sqrt(k1 * k2) * h
            / (np.log(re / rw) + skin))


def _bounding_box_well_index(geometry, cell, kx=150.0, ky=150.0) -> float:
    """KÖHNƏ (düzəlişdən əvvəlki) hesab — ox-boyu sərhəd qutusu ilə.

    Testlərin "bounding-box-dan ASILI DEYİL" iddiasını ÖLÇƏ bilməsi
    üçün lazımdır: yeni nəticə bundan FƏRQLƏNMƏLİDİR."""
    extents = geometry.cell_extents()[cell]
    return _peaceman(kx, ky, float(extents[0]), float(extents[1]),
                     float(geometry.cell_thickness()[cell]))


def _cartesian(dz=None) -> CellGeometry:
    return CellGeometry(grid=_grid(), dx=DX, dy=DY,
                        dz=list(dz if dz is not None else DZ_LAYERS),
                        top_depth=TOP)


def _transformed(geometry: CornerPointGeometry, matrix) -> CornerPointGeometry:
    """Təpələrə 3×3 xətti çevirmə tətbiq edir — həndəsə DEFORMASİYA
    OLUNUR, hüceyrə sayı/topologiyası dəyişmir."""
    nodes = np.asarray(geometry.nodes, float) @ np.asarray(matrix, float).T
    return CornerPointGeometry.from_nodes(geometry.grid, nodes)


def _rotated(geometry, degrees) -> CornerPointGeometry:
    """XY müstəvisində FIRLANMA — həcm, ayaq izi sahəsi və hər kənarın
    uzunluğu DƏYİŞMİR, yalnız oxlar dönür."""
    angle = np.deg2rad(degrees)
    cos, sin = np.cos(angle), np.sin(angle)
    return _transformed(geometry, [[cos, -sin, 0.0], [sin, cos, 0.0],
                                   [0.0, 0.0, 1.0]])


def _sheared(geometry, factor) -> CornerPointGeometry:
    """X boyunca KƏSİLMƏ (skew) — həcm Kavalyeri prinsipinə görə
    dəyişmir, ox-boyu sərhəd qutusu isə BÖYÜYÜR."""
    return _transformed(geometry, [[1.0, factor, 0.0], [0.0, 1.0, 0.0],
                                   [0.0, 0.0, 1.0]])


# ══════════════════════════ 1 — Kartezian nəticə DƏYİŞMİR ═══════════════
def test_cartesian_well_index_matches_the_classic_hand_calculation():
    """Bərabər izotrop blok: `re = 0.28·√(dx²+dy²)/2`, `WI ∝ dz`."""
    values = _well_index(_cartesian())
    expected = [_peaceman(150.0, 150.0, DX, DY, dz) for dz in DZ_LAYERS]
    assert values == pytest.approx(expected, rel=1e-12)


def test_corner_point_reproduces_the_cartesian_answer_exactly():
    """Kartezian model CPG-nin XÜSUSİ HALIDIR — WI maşın dəqiqliyində
    eyni olmalıdır (yeni yol köhnəni əvəz etdiyi üçün bu, ən sərt
    reqressiya qıfılıdır)."""
    box = _cartesian()
    assert _well_index(CornerPointGeometry.from_cartesian(box)) == pytest.approx(
        _well_index(box), rel=1e-12)


def test_skin_still_reduces_the_well_index():
    """Peaceman düsturunun ÖZÜ dəyişməyib — yalnız həndəsi girişi."""
    geometry = CornerPointGeometry.from_cartesian(_cartesian())
    assert np.all(_well_index(geometry, skin=4.0) < _well_index(geometry))


# ══════════════════════ 2 — dəyişkən qalınlıq: WI lay-lay ═══════════════
def test_well_index_follows_each_layers_own_real_thickness():
    values = _well_index(_cpg(layer_dz=DZ_LAYERS))
    assert len(values) == NZ
    assert values[1] / values[0] == pytest.approx(DZ_LAYERS[1] / DZ_LAYERS[0],
                                                  rel=1e-12)


def test_variable_thickness_perforation_length_is_the_cells_own_k_extent():
    """`h` layın ORTALAMASI deyil, hüceyrənin öz K uzanmasıdır."""
    geometry = _cpg(layer_dz=DZ_LAYERS)
    cells = [geometry.grid.index(WELL_I, WELL_J, k) for k in range(NZ)]
    block = geometry.wellblock_geometry(cells)
    assert block.length == pytest.approx(DZ_LAYERS, rel=1e-12)


# ═══════════════ 3 — maili (dipping) hüceyrə: REAL yerli həndəsə ════════
def test_dipping_column_well_index_comes_from_the_real_local_geometry():
    """Maili pillar: quyu oxu ŞAQULİ DEYİL, hüceyrənin öz K oxudur.

    `h` həmin ox boyunca həqiqi uzunluqdur, en kəsiyi isə `V/h` — nə
    nominal `dx·dy`, nə də sərhəd qutusu."""
    geometry = _cpg(tilt_x=420.0, dip_x=18.0)
    cell = geometry.grid.index(WELL_I, WELL_J, 0)
    block = geometry.wellblock_geometry([cell])
    k1, k2 = block.directional_permeability([[150.0, 150.0, 15.0]])

    manual = _peaceman(float(k1[0]), float(k2[0]), float(block.d1[0]),
                       float(block.d2[0]), float(block.length[0]))
    assert _well_index(geometry)[0] == pytest.approx(manual, rel=1e-12)

    # ... və bu, ox-boyu sərhəd qutusundan gələn cavab DEYİL
    assert _well_index(geometry)[0] != pytest.approx(
        _bounding_box_well_index(geometry, cell), rel=1e-6)


def test_dipping_cell_perforation_length_is_longer_than_the_vertical_drop():
    """Maili pillarda quyu oxu boyunca yol ŞAQULİ enişdən UZUNDUR."""
    geometry = _cpg(tilt_x=420.0)
    cell = geometry.grid.index(WELL_I, WELL_J, 0)
    block = geometry.wellblock_geometry([cell])
    vertical = float(geometry.cell_thickness()[cell])
    assert float(block.length[0]) > vertical
    assert float(np.abs(block.well_axis[0, 0])) > 1e-3      # ox əyilib


# ════════ 4 — fırlanmış/kəsilmiş hüceyrə: sərhəd qutusundan ASILI DEYİL ═
def test_rotating_the_grid_does_not_change_the_well_index():
    """FİZİKA ARQUMENTİ: XY müstəvisində fırlanma hüceyrənin nə həcmini,
    nə ayaq izini, nə də quyuya nisbətini dəyişir — WI SABİT qalmalıdır.

    Sərhəd qutusu ilə bu MÜMKÜN DEYİL: 45°-də 100×80 m blokun qutusu
    ~127×127 m olur, yəni köhnə hesab eyni fiziki hüceyrəyə fərqli WI
    verirdi."""
    base = CornerPointGeometry.from_cartesian(_cartesian())
    rotated = _rotated(base, 45.0)
    cell = base.grid.index(WELL_I, WELL_J, 0)

    assert _well_index(rotated) == pytest.approx(_well_index(base), rel=1e-12)
    # sərhəd qutusu HƏQİQƏTƏN şişib — test boş yerə keçmir:
    # 45°-də qutunun hər tərəfi `(dx+dy)/√2` olur.
    box = rotated.cell_extents()[cell]
    assert box[0] == pytest.approx((DX + DY) / np.sqrt(2.0), rel=1e-9)
    assert box[1] == pytest.approx((DX + DY) / np.sqrt(2.0), rel=1e-9)
    assert _bounding_box_well_index(rotated, cell) < 0.97 * _well_index(base)[0]


def test_warped_cell_uses_the_real_footprint_area_not_the_bounding_box():
    """Kəsilmiş hüceyrə: `d1·d2` HƏQİQİ en kəsiyinə bərabərdir."""
    base = CornerPointGeometry.from_cartesian(_cartesian())
    sheared = _sheared(base, 0.6)
    cell = sheared.grid.index(WELL_I, WELL_J, 0)
    block = sheared.wellblock_geometry([cell])

    real_area = sheared.volumes()[cell] / float(block.length[0])
    assert float(block.cross_section()[0]) == pytest.approx(real_area, rel=1e-12)
    assert real_area == pytest.approx(DX * DY, rel=1e-12)   # kəsilmə sahəni saxlayır

    box = sheared.cell_extents()[cell]
    assert box[0] * box[1] > 1.4 * real_area                # qutu şişib
    assert _well_index(sheared)[0] != pytest.approx(
        _bounding_box_well_index(sheared, cell), rel=1e-6)


def test_warped_deck_well_index_matches_its_own_vertex_geometry():
    """Tam qeyri-müntəzəm deck (maili pillar + əyri səth): hər
    perforasiya öz təpələrindən hesablanan dəyəri verir."""
    geometry = _cpg(tilt_x=300.0, tilt_wave=120.0, warp=7.0,
                    layer_dz=DZ_LAYERS)
    cells = [geometry.grid.index(WELL_I, WELL_J, k) for k in range(NZ)]
    block = geometry.wellblock_geometry(cells)
    k1, k2 = block.directional_permeability(
        np.tile([150.0, 150.0, 15.0], (len(cells), 1)))
    manual = [_peaceman(k1[n], k2[n], block.d1[n], block.d2[n], block.length[n])
              for n in range(len(cells))]
    assert _well_index(geometry) == pytest.approx(manual, rel=1e-12)


# ══════════════════════════ 5 — anizotrop keçiricilik ═══════════════════
def test_anisotropic_permeability_changes_the_well_index():
    geometry = CornerPointGeometry.from_cartesian(_cartesian())
    isotropic = _well_index(geometry, kx=150.0, ky=150.0)
    anisotropic = _well_index(geometry, kx=450.0, ky=50.0)
    assert not np.allclose(isotropic, anisotropic)
    expected = [_peaceman(450.0, 50.0, DX, DY, dz) for dz in DZ_LAYERS]
    assert anisotropic == pytest.approx(expected, rel=1e-12)


def test_anisotropy_is_read_along_the_cells_own_axes_after_rotation():
    """90° fırlanmadan sonra hüceyrənin I oxu QLOBAL Y-ə baxır — yönlü
    keçiricilik ONUNLA BİRLİKDƏ dönməlidir: `K1 = Ky`, `K2 = Kx`.

    Keçiricilik sahəsi QLOBAL oxlarda verilir (PERMX qlobal X boyuncadır)
    və hüceyrə ilə birlikdə FIRLANMIR — ona görə bu, fiziki cəhətdən
    BAŞQA məsələdir və WI fırlanmamış blokdakı ilə eyni OLMAMALIDIR.
    Yoxlanılan budur: `K1`/`K2` yerli oxlara `uᵀ·K·u` ilə düzgün
    proyeksiya olunur və WI məhz həmin cütdən çıxır."""
    base = CornerPointGeometry.from_cartesian(_cartesian())
    rotated = _rotated(base, 90.0)
    cell = base.grid.index(WELL_I, WELL_J, 0)

    block = rotated.wellblock_geometry([cell])
    k1, k2 = block.directional_permeability([[450.0, 50.0, 15.0]])
    assert float(k1[0]) == pytest.approx(50.0, abs=1e-9)      # yerli I ‖ qlobal Y
    assert float(k2[0]) == pytest.approx(450.0, abs=1e-9)     # yerli J ‖ qlobal X
    assert float(block.d1[0]) == pytest.approx(DX, rel=1e-12)
    assert float(block.d2[0]) == pytest.approx(DY, rel=1e-12)

    values = _well_index(rotated, kx=450.0, ky=50.0)
    expected = [_peaceman(50.0, 450.0, DX, DY, dz) for dz in DZ_LAYERS]
    assert values == pytest.approx(expected, rel=1e-12)
    assert not np.allclose(values, _well_index(base, kx=450.0, ky=50.0))


def test_anisotropy_projects_onto_a_rotated_cell_axis():
    """45°-də hər iki yerli ox `Kx` və `Ky`-nin ORTALAMASINI görür
    (`uᵀ·K·u`, `u = (√2/2, √2/2, 0)`)."""
    rotated = _rotated(CornerPointGeometry.from_cartesian(_cartesian()), 45.0)
    block = rotated.wellblock_geometry([rotated.grid.index(WELL_I, WELL_J, 0)])
    k1, k2 = block.directional_permeability([[450.0, 50.0, 15.0]])
    assert float(k1[0]) == pytest.approx(250.0, rel=1e-12)
    assert float(k2[0]) == pytest.approx(250.0, rel=1e-12)


# ═════════════════ 6 — tamamlama: HƏR perforasiya ÖZ hüceyrəsi ══════════
def test_each_perforation_uses_its_own_cell_geometry():
    """Dəyişkən qalınlıq + maili pillar: heç iki perforasiya eyni
    həndəsəni PAYLAŞMIR və hər biri öz hüceyrəsinin dəyərini alır."""
    geometry = _cpg(layer_dz=DZ_LAYERS, tilt_x=380.0, warp=6.0)
    connections = _connections(geometry)
    assert [c.cell for c in connections] == [
        geometry.grid.index(WELL_I, WELL_J, k) for k in range(NZ)]

    block = geometry.wellblock_geometry([c.cell for c in connections])
    assert len(set(np.round(block.length, 9))) == NZ
    assert len(set(np.round([c.well_index for c in connections], 9))) == NZ


def test_two_wells_in_different_columns_see_different_geometry():
    """Eyni deck, fərqli sütun → fərqli yerli həndəsə → fərqli WI."""
    geometry = _cpg(tilt_x=380.0, tilt_wave=140.0, warp=6.0)
    left = _well_index(geometry, i=0, j=0)
    right = _well_index(geometry, i=NX - 1, j=NY - 1)
    assert not np.allclose(left, right)


def test_horizontal_completion_uses_the_i_axis_of_its_own_cell():
    """`direction="X"` — perforasiya hüceyrəni I oxu boyunca delir:
    `h = dx`, effektiv ölçülər isə `dy`/`dz`."""
    geometry = CornerPointGeometry.from_cartesian(_cartesian())
    cells = [geometry.grid.index(WELL_I, WELL_J, k) for k in range(NZ)]
    block = geometry.wellblock_geometry(cells, direction="X")

    assert block.length == pytest.approx([DX] * NZ, rel=1e-12)
    assert block.d1 == pytest.approx([DY] * NZ, rel=1e-12)
    assert block.d2 == pytest.approx(DZ_LAYERS, rel=1e-12)

    expected = [_peaceman(150.0, 15.0, DY, dz, DX) for dz in DZ_LAYERS]
    assert _well_index(geometry, direction="X") == pytest.approx(expected,
                                                                 rel=1e-12)


def test_completion_direction_changes_the_well_index():
    geometry = CornerPointGeometry.from_cartesian(_cartesian())
    assert not np.allclose(_well_index(geometry, direction="Z"),
                           _well_index(geometry, direction="X"))


def test_mixed_completion_directions_stay_aligned_with_their_perforations():
    """Bir quyuda müxtəlif istiqamətlər — `build_connections` onları
    QRUPLARLA hesablayır, nəticə sırası isə perforasiya sırası ilə
    QALMALIDIR (indeks sürüşməsi reqressiyası)."""
    geometry = CornerPointGeometry.from_cartesian(_cartesian())
    grid = geometry.grid
    well = Well(name="P-MIX", well_type=WellType.PRODUCER,
                control=WellControl(ControlMode.BHP, 150.0),
                perforations=[Perforation(WELL_I, WELL_J, 0, True, 0.0, "X"),
                              Perforation(WELL_I, WELL_J, 1, True, 0.0, "Z")],
                radius=RADIUS)
    model = ReservoirModel(name="mixed", grid=grid, geometry=geometry,
                           rock=_rock(grid.ncell), wells=[well])
    values = [c.well_index for c in PeacemanWellModel().build_connections(model)]

    assert values[0] == pytest.approx(
        _peaceman(150.0, 15.0, DY, DZ_LAYERS[0], DX), rel=1e-12)
    assert values[1] == pytest.approx(
        _peaceman(150.0, 150.0, DX, DY, DZ_LAYERS[1]), rel=1e-12)
