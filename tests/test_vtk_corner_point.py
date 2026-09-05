"""VTK görüntüsü — HƏQİQİ corner-point həndəsəsi.

    COORD + ZCORN
            ↓
    real 8 vertices / cell
            ↓
    VTK geometry

Əvvəl `VtkReservoirScene` HƏMİŞƏ nizamlı Kartezian şəbəkə qururdu
(`x = i·dx`, `y = j·dy`, `z = mərkəz ± dz/2`) — yəni corner-point deck
oxunsa da, ekranda onun maili pilları, dəyişkən qalınlığı və əyri
hüceyrələri GÖRÜNMÜRDÜ: `dx`/`dy` corner-point modeldə yalnız NOMİNAL
ortalamadır.

Bu modul iki şeyi ayrıca qıfıllayır:
  · Kartezian model — görüntü DƏYİŞMƏYİB (paylaşılan düyünlü
    `vtkStructuredGrid`, `tests/test_vtk_volume.py` onu ətraflı yoxlayır);
  · corner-point model — hər hüceyrə deck-dəki ÖZ 8 təpəsi ilə çəkilir.

Həqiqiliyin ölçüsü kimi VTK-nın ÖZ `vtkCellSizeFilter`-i işlədilir:
o, ekrana gedən həndəsədən həcmi hesablayır. Əgər görüntü nominal
qutulardan qurulsaydı, bu həcm domendəki dəqiq çoxüzlü həcmə
UYĞUN GƏLMƏZDİ.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.domain.corner_point_geometry import CornerPointGeometry
from imex2d.domain.geometry import CellGeometry
from imex2d.domain.properties import PropertyMap, RockProperties
from imex2d.domain.reservoir_model import ReservoirModel
from imex2d.rendering import vtk_volume

from test_local_cell_metrics import (DX, DY, DZ, NX, NY, NZ, TOP, _coord, _cpg,
                                     _grid, _zcorn)

vtk = pytest.importorskip("vtk", reason="VTK quraşdırılmayıb")


# ══════════════════════════════════════════════════ köməkçilər ══════════
def _model(geometry):
    grid = geometry.grid
    rock = RockProperties(
        porosity=PropertyMap.uniform("PORO", 0.24, grid.ncell),
        permx=PropertyMap.uniform("PERMX", 150.0, grid.ncell),
        permy=PropertyMap.uniform("PERMY", 150.0, grid.ncell))
    return ReservoirModel(name="vtk-cpg", grid=grid, geometry=geometry,
                          rock=rock, wells=[])


def _scene(geometry, **settings):
    return vtk_volume.VtkReservoirScene(
        _model(geometry), vtk_volume.VtkViewSettings(**settings))


def _cell_hexahedra(scene) -> np.ndarray:
    """`(ncell, 8, 3)` — VTK-nın SƏHNƏDƏ saxladığı təpələr, hüceyrə-hüceyrə.

    Səhnədən OXUNUR (giriş massivindən yox) — testin ölçdüyü şey məhz
    ekrana gedən həndəsədir."""
    grid = scene._grid
    cells = np.empty((grid.GetNumberOfCells(), 8, 3))
    for index in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(index)
        points = cell.GetPoints()
        for corner in range(8):
            cells[index, corner] = points.GetPoint(corner)
    return cells


def _vtk_cell_volumes(scene) -> np.ndarray:
    """VTK-nın ÖZ hesabladığı hüceyrə həcmləri — görüntü həndəsəsindən."""
    from vtkmodules.util import numpy_support

    sizer = vtk.vtkCellSizeFilter()
    sizer.SetInputData(scene._grid)
    sizer.ComputeVolumeOn()
    sizer.ComputeSumOff()
    sizer.Update()
    return numpy_support.vtk_to_numpy(
        sizer.GetOutput().GetCellData().GetArray("Volume"))


def _cartesian_geometry(dz=DZ) -> CellGeometry:
    return CellGeometry(grid=_grid(), dx=DX, dy=DY, dz=dz, top_depth=TOP)


# ═══════════════════ 1 — Kartezian görüntü DƏYİŞMİR ═════════════════════
def test_cartesian_model_still_uses_the_shared_node_structured_grid():
    """Bərabər bloklarda paylaşılan düyünlü struktur şəbəkə HƏM yığcam,
    HƏM tam dəqiqdir — yeni yol ona TOXUNMUR."""
    scene = _scene(_cartesian_geometry())
    assert scene.uses_corner_point_geometry() is False
    assert isinstance(scene._grid, vtk.vtkStructuredGrid)
    dimensions = [0, 0, 0]
    scene._grid.GetDimensions(dimensions)
    assert tuple(dimensions) == (NX + 1, NY + 1, NZ + 1)


def test_cartesian_node_coordinates_are_unchanged():
    """Nizamlı düyün şəbəkəsi (`i·dx`, `j·dy`) Kartezian modeldə
    HƏLƏ DƏ doğrudur və hərfi olaraq saxlanılıb."""
    scene = _scene(_cartesian_geometry())
    points = scene._cell_corner_points()
    x = np.unique(np.round(points[:, 0], 9))
    y = np.unique(np.round(points[:, 1], 9))
    assert x == pytest.approx(np.arange(NX + 1) * DX)
    assert y == pytest.approx(np.arange(NY + 1) * DY)


def test_cartesian_scene_renders_and_carries_cell_values():
    scene = _scene(_cartesian_geometry())
    scene.update_values(np.arange(float(scene.model.ncell)), label="test")
    scalars = scene._grid.GetCellData().GetScalars()
    assert scalars.GetNumberOfTuples() == scene.model.ncell


# ═════════════ 2 — corner-point: REAL 8 təpə / hüceyrə ══════════════════
def test_corner_point_model_switches_to_per_cell_hexahedra():
    scene = _scene(_cpg())
    assert scene.uses_corner_point_geometry() is True
    assert isinstance(scene._grid, vtk.vtkUnstructuredGrid)
    assert scene._grid.GetNumberOfCells() == NX * NY * NZ
    # təpələr PAYLAŞILMIR — hər hüceyrə öz 8 nöqtəsini daşıyır
    assert scene._grid.GetNumberOfPoints() == NX * NY * NZ * 8
    assert scene._grid.GetCellType(0) == vtk.VTK_HEXAHEDRON


def test_cell_order_still_matches_the_model_cell_indexing():
    """Hüceyrə sırası (`i` ən sürətli) dəyişməməlidir — `update_values()`
    skalyarları məhz bu sıra ilə bağlayır."""
    geometry = _cpg(dip_x=30.0, layer_dz=[6.0, 21.0])
    scene = _scene(geometry)
    centres = _cell_hexahedra(scene).mean(axis=1)
    reference = geometry.cell_centroid()
    assert centres[:, 0] == pytest.approx(reference[:, 0], abs=1e-9)
    assert centres[:, 1] == pytest.approx(reference[:, 1], abs=1e-9)
    assert centres[:, 2] == pytest.approx(-reference[:, 2], abs=1e-9)


# ════════════════ 3 — maili (dipping) corner-point grid ═════════════════
def test_dipping_grid_produces_tilted_cells_not_flat_layers():
    """Maili layda hüceyrənin tavanı I boyunca ENMƏLİDİR — nizamlı
    Kartezian şəbəkədə bütün tavan künc dərinlikləri EYNİ olardı."""
    dip = 30.0
    scene = _scene(_cpg(dip_x=dip))
    cells = _cell_hexahedra(scene)
    grid = _grid()

    for i in range(NX - 1):
        near = cells[grid.index(i, 0, 0)].mean(axis=0)
        far = cells[grid.index(i + 1, 0, 0)].mean(axis=0)
        # VTK-da Z yuxarı müsbətdir: dərinləşən lay AŞAĞI düşür
        assert far[2] < near[2]
        assert near[2] - far[2] == pytest.approx(dip, rel=1e-9)


def test_dipping_cell_top_face_is_not_horizontal():
    """Tək bir hüceyrənin ÖZÜ maildir — bu, "lay sürüşdürülüb" ilə
    "hüceyrə əyilib" arasındakı fərqdir."""
    scene = _scene(_cpg(dip_x=30.0))
    top_face = _cell_hexahedra(scene)[0][4:8]          # VTK-da yuxarı üz
    assert top_face[:, 2].max() - top_face[:, 2].min() > 1.0


def test_inclined_pillars_make_the_cell_lean_in_xy():
    """Maili pillar (`COORD`-un alt nöqtəsi sürüşüb): hüceyrənin tavanı
    və dabanı FƏRQLİ X ayaq izində olur. Nizamlı şəbəkədə bu MÜMKÜN
    DEYİL — orada X yalnız `i`-dən asılıdır."""
    scene = _scene(_cpg(tilt_x=400.0))
    cell = _cell_hexahedra(scene)[0]
    bottom_x = cell[0:4, 0].mean()
    top_x = cell[4:8, 0].mean()
    assert abs(bottom_x - top_x) > 1.0


# ═════════════════ 4 — dəyişkən qalınlıq: hər lay ÖZÜ ═══════════════════
def test_each_layer_keeps_its_own_real_thickness():
    layers = [6.0, 21.0]
    scene = _scene(_cpg(layer_dz=layers))
    cells = _cell_hexahedra(scene)
    grid = _grid()
    for k, thickness in enumerate(layers):
        cell = cells[grid.index(1, 1, k)]
        drawn = cell[4:8, 2].mean() - cell[0:4, 2].mean()      # yuxarı − aşağı
        assert drawn == pytest.approx(thickness, rel=1e-9)


def test_variable_thickness_is_not_smoothed_into_a_mean_layer():
    """Köhnə yol hər layı `dz` (ortalama) ilə çəkirdi — 6/21 m fərqi
    itirdi. İndi nisbət olduğu kimi görünür."""
    layers = [6.0, 21.0]
    scene = _scene(_cpg(layer_dz=layers))
    cells = _cell_hexahedra(scene)
    grid = _grid()
    drawn = []
    for k in range(NZ):
        cell = cells[grid.index(1, 1, k)]
        drawn.append(cell[4:8, 2].mean() - cell[0:4, 2].mean())
    assert drawn[1] / drawn[0] == pytest.approx(layers[1] / layers[0], rel=1e-9)
    assert not np.isclose(drawn[0], np.mean(layers))


def test_vertical_exaggeration_scales_only_depth_on_corner_point_cells():
    plain = _cell_hexahedra(_scene(_cpg(dip_x=25.0, layer_dz=[6.0, 21.0])))
    stretched = _cell_hexahedra(_scene(_cpg(dip_x=25.0, layer_dz=[6.0, 21.0]),
                                       vertical_exaggeration=4.0))
    assert stretched[..., 0] == pytest.approx(plain[..., 0])
    assert stretched[..., 1] == pytest.approx(plain[..., 1])
    assert stretched[..., 2] == pytest.approx(plain[..., 2] * 4.0)


# ═══════════════ 5 — əyri (warped) hüceyrə: 8 REAL təpə ═════════════════
def test_warped_cell_geometry_matches_its_eight_real_vertices():
    """Ən birbaşa iddia: VTK-dakı hər hüceyrə deck-in HƏMİN 8 təpəsidir
    (yalnız dərinlik işarəsi çevrilmiş halda)."""
    geometry = _cpg(tilt_x=300.0, tilt_wave=120.0, warp=7.0,
                    layer_dz=[6.0, 21.0])
    cells = _cell_hexahedra(_scene(geometry))
    nodes = np.asarray(geometry.nodes, float)

    # VTK sırası domen sırasının yarımlarını dəyişir (aşağı üz əvvəl)
    expected = np.empty_like(nodes)
    expected[:, 0:4] = nodes[:, 4:8]
    expected[:, 4:8] = nodes[:, 0:4]
    expected[..., 2] = -expected[..., 2]
    assert cells == pytest.approx(expected, abs=1e-9)


def test_warped_cells_have_positive_volume_in_the_scene():
    """Təpə sırası səhv olsaydı yakobian mənfi çıxar, normallar içəri
    baxar və işıqlandırma tərsinə düşərdi — VTK-nın öz həcmi bunu tutur."""
    scene = _scene(_cpg(tilt_x=300.0, tilt_wave=120.0, warp=7.0))
    assert np.all(_vtk_cell_volumes(scene) > 0.0)


@pytest.mark.parametrize("deck", [
    {},
    {"dip_x": 30.0, "dip_y": 12.0},
    {"layer_dz": [6.0, 21.0]},
    {"tilt_x": 300.0, "tilt_wave": 120.0, "warp": 7.0},
])
def test_scene_volume_equals_the_exact_corner_point_volume(deck):
    """VTK-nın ekrandakı həndəsədən hesabladığı həcm domendəki DƏQİQ
    çoxüzlü həcmə bərabərdir — yəni görüntü heç bir yerdə nominal
    qutuya qayıtmır."""
    geometry = _cpg(**deck)
    assert _vtk_cell_volumes(_scene(geometry)) == pytest.approx(
        geometry.volumes(), rel=1e-9)


def test_nominal_box_volume_would_have_been_visibly_wrong():
    """Testin boş yerə keçmədiyinin sübutu: nominal `dx·dy·dz` qutusu
    həmin deck üçün FƏRQLİ həcm verir."""
    geometry = _cpg(layer_dz=[6.0, 21.0], tilt_x=300.0, warp=7.0)
    nominal = (float(geometry.dx) * float(geometry.dy)
               * np.repeat(geometry.dz, NX * NY))
    assert not np.allclose(_vtk_cell_volumes(_scene(geometry)), nominal)


# ══════════ 6 — COORD/ZCORN dəyişəndə görüntü də DƏYİŞİR ════════════════
def test_changing_zcorn_changes_the_drawn_geometry():
    grid = _grid()
    coord = _coord(nx=grid.nx, ny=grid.ny)
    flat = _zcorn(nx=grid.nx, ny=grid.ny, nz=grid.nz)
    dipping = _zcorn(nx=grid.nx, ny=grid.ny, nz=grid.nz, dip_x=35.0)

    first, _ = CornerPointGeometry.from_grdecl(grid, coord, flat)
    second, _ = CornerPointGeometry.from_grdecl(grid, coord, dipping)
    drawn_first = _cell_hexahedra(_scene(first))
    drawn_second = _cell_hexahedra(_scene(second))

    assert not np.allclose(drawn_first[..., 2], drawn_second[..., 2])
    assert np.allclose(drawn_first[..., 0], drawn_second[..., 0])   # X toxunulmayıb


def test_changing_coord_changes_the_drawn_geometry():
    """`COORD` yalnız pillar OXLARINI dəyişir — ZCORN eynidir, amma
    hüceyrələrin X/Y ayaq izi sürüşməlidir."""
    grid = _grid()
    zcorn = _zcorn(nx=grid.nx, ny=grid.ny, nz=grid.nz)
    straight = _coord(nx=grid.nx, ny=grid.ny)
    tilted = _coord(nx=grid.nx, ny=grid.ny, tilt_x=400.0)

    first, _ = CornerPointGeometry.from_grdecl(grid, straight, zcorn)
    second, _ = CornerPointGeometry.from_grdecl(grid, tilted, zcorn)
    drawn_first = _cell_hexahedra(_scene(first))
    drawn_second = _cell_hexahedra(_scene(second))

    assert not np.allclose(drawn_first[..., 0], drawn_second[..., 0])
    assert np.allclose(drawn_first[..., 2], drawn_second[..., 2])   # Z toxunulmayıb


def test_scene_bounds_follow_the_deck():
    """Kamera/oxlar `GetBounds()`-a güvənir — o da real təpələrdən
    gəlməlidir."""
    geometry = _cpg(dip_x=40.0)
    bounds = _scene(geometry)._grid.GetBounds()
    nodes = np.asarray(geometry.nodes, float)
    assert bounds[0] == pytest.approx(nodes[..., 0].min(), abs=1e-9)
    assert bounds[1] == pytest.approx(nodes[..., 0].max(), abs=1e-9)
    assert bounds[4] == pytest.approx(-nodes[..., 2].max(), abs=1e-9)
    assert bounds[5] == pytest.approx(-nodes[..., 2].min(), abs=1e-9)


# ═════════════════ görüntü boru xəttinin qalanı işləyir ═════════════════
def test_values_filters_and_offscreen_render_work_on_corner_point_cells():
    """Kəsim həddi (`BlankCell`), rəng legendi və oflayn render — hamısı
    struktur şəbəkəyə bağlı deyil, ümumi `vtkDataSet` API-si üzərindədir."""
    geometry = _cpg(tilt_x=300.0, warp=7.0, layer_dz=[6.0, 21.0])
    scene = _scene(geometry, value_min=0.5, k_range=(0, 0))
    values = np.linspace(0.0, 1.0, geometry.grid.ncell)
    scene.update_values(values, label="test")
    assert scene._grid.GetCellGhostArray() is not None

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetSize(120, 90)
    window.AddRenderer(scene.renderer)
    scene.reset_camera("İzometrik")
    window.Render()

    capture = vtk.vtkWindowToImageFilter()
    capture.SetInput(window)
    capture.Update()
    assert capture.GetOutput().GetNumberOfPoints() > 0


def test_wells_sit_on_the_real_cell_centre_of_a_corner_point_column():
    """Lülə nominal `(i+0.5)·dx`-də deyil, hüceyrənin HƏQİQİ mərkəzində
    olmalıdır — maili pillarda ikisi onlarla metr fərqlənir."""
    from imex2d.domain.wells import (ControlMode, Perforation, Well,
                                     WellControl, WellType)

    geometry = _cpg(tilt_x=400.0)
    model = _model(geometry)
    well = Well(name="P1", well_type=WellType.PRODUCER,
                control=WellControl(ControlMode.BHP, 150.0),
                perforations=[Perforation(1, 1, k) for k in range(NZ)])
    object.__setattr__(model, "wells", [well])

    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(np.ones(model.ncell))
    assert scene._well_actors

    centroid = geometry.cell_centroid()[geometry.grid.index(1, 1, 0)]
    nominal_x = 1.5 * float(geometry.dx)
    # `abs=1e-3` — VTK poliqon nöqtələrini `float32`-də saxlayır, yəni
    # 150 m-lik koordinatda ~1e-5 m yuvarlaqlaşma qaçılmazdır.
    bore = scene._well_actors[0].GetCenter()          # aktyorun səhnədəki yeri
    assert bore[0] == pytest.approx(centroid[0], abs=1e-3)
    assert bore[1] == pytest.approx(centroid[1], abs=1e-3)
    # ... və nominal `(i+0.5)·dx` fərqli yerdədir (bu deck-də bir neçə
    # metr; dərin laylarda maili pillar fərqi onlarla metrə çatdırır)
    assert abs(centroid[0] - nominal_x) > 1.0
