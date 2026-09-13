"""B6-a — interaktiv kəsik müstəvisi (`vtkCutter` + `vtkImplicitPlaneWidget2`).

ÖLÇÜLMÜŞ FAKTLAR (`ISH_HESABATI.md` → Seans 21):
  * VTK 9.7.0-da ekransız render və `vtkImplicitPlaneWidget2` İŞLƏYİR —
    widget-i məntiqdən ayırmağa ehtiyac yoxdur, birbaşa sınanır;
  * `vtkCutter` hüceyrə skalyarlarını daşıyır → kəsik əsas həcmlə EYNİ
    rəng cədvəlini işlədə bilir;
  * `vtkCutter` HƏR kəsilən hüceyrəni İKİ ÜÇBUCAQ kimi verir (8×6×3-də
    X-kəsiyi 36 üçbucaq = 18 hüceyrə). Ona görə testlər poliqon sayını
    DEYİL, unikal hüceyrə (skalyar) sayını yoxlayır;
  * müstəvi dəqiq hüceyrə sərhədinə düşəndə unikal hüceyrə sayı dəyişmir —
    ilkin "üst-üstə poliqon" fərziyyəsi ölçmə ilə TƏKZİB olundu;
  * `vtkCutter` gizlədilmiş (blank) hüceyrələrə hörmət edir → K-filtri,
    kəsim həddi, status filtri kəsiyə AVTOMATİK tətbiq olunur;
  * `SetNormalToZAxis(1)` normal VEKTORUNU dəyişmir — vektor açıq verilir.
"""

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.rendering import vtk_volume
from imex2d.rendering.vtk_volume import (VtkViewSettings, slice_fraction,
                                         slice_plane)

# `test_vtk_volume.py`-dəki ilə eyni qoruma: Smart App Control VTK
# DLL-lərini bloklayanda `ImportError` bütün dəsti dayandırmamalıdır.
try:
    import vtk
    from vtkmodules.util import numpy_support
except ImportError as exc:                                   # pragma: no cover
    pytest.skip(f"VTK əlçatan deyil: {exc}", allow_module_level=True)

NX, NY, NZ = 8, 6, 3
BOUNDS = (0.0, 200.0, 0.0, 150.0, -2030.0, -2000.0)


def _model():
    geology = SyntheticGeologicalModelBuilder().build(
        nx=NX, ny=NY, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=NZ, top_depth=2000.0)
    return ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                         scal=default_scal())


def _scene(**settings):
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model, VtkViewSettings(**settings))
    scene.update_values(np.arange(model.ncell, dtype=float), "sınaq")
    return scene, model


def _cut_cells(scene) -> set:
    """Kəsilən UNİKAL hüceyrələr — hər hüceyrə öz indeksini dəyər kimi daşıyır."""
    scalars = scene.slice_output().GetCellData().GetScalars()
    return set(numpy_support.vtk_to_numpy(scalars).tolist())


# ═══════════════════════ saf həndəsə ═════════════════════════════════

@pytest.mark.parametrize("axis, index", [("X", 0), ("Y", 1), ("Z", 2)])
def test_plane_normal_points_along_the_chosen_axis(axis, index):
    _origin, normal = slice_plane(BOUNDS, axis, 0.5)
    expected = [0.0, 0.0, 0.0]
    expected[index] = 1.0
    assert normal == tuple(expected)


def test_plane_origin_follows_the_fraction():
    origin, _ = slice_plane(BOUNDS, "X", 0.25)
    assert origin[0] == pytest.approx(50.0)
    # digər iki koordinat modelin mərkəzindədir
    assert origin[1] == pytest.approx(75.0)
    assert origin[2] == pytest.approx(-2015.0)


def test_fraction_is_clamped_inside_the_model():
    low, _ = slice_plane(BOUNDS, "Y", -3.0)
    high, _ = slice_plane(BOUNDS, "Y", 7.0)
    assert low[1] == pytest.approx(BOUNDS[2])
    assert high[1] == pytest.approx(BOUNDS[3])


def test_invalid_axis_is_rejected():
    with pytest.raises(ValueError):
        slice_plane(BOUNDS, "W", 0.5)


@pytest.mark.parametrize("axis", ["X", "Y", "Z"])
def test_slice_fraction_inverts_slice_plane(axis):
    for fraction in (0.0, 0.1, 0.5, 0.9, 1.0):
        origin, _ = slice_plane(BOUNDS, axis, fraction)
        assert slice_fraction(BOUNDS, axis, origin) == pytest.approx(fraction)


# ═══════════════════════ səhnə — kəsik söndürülü ═════════════════════

def test_default_has_no_slice_and_the_volume_is_visible():
    """GERİYƏ UYĞUNLUQ: parametr verilməyəndə heç nə dəyişmir."""
    scene, _ = _scene()
    assert scene.slice_output() is None
    assert scene._actor.GetVisibility() == 1


# ═══════════════════════ səhnə — kəsik aktiv ═════════════════════════

@pytest.mark.parametrize("axis, expected", [
    ("X", NY * NZ),          # bir sütun təbəqəsi: 6 × 3
    ("Y", NX * NZ),          # 8 × 3
    ("Z", NX * NY),          # 8 × 6
])
def test_slice_cuts_exactly_one_layer_of_cells(axis, expected):
    scene, _ = _scene(slice_axis=axis, slice_position=0.43)
    assert len(_cut_cells(scene)) == expected


def test_each_cut_cell_is_drawn_as_two_triangles():
    """Ölçülmüş VTK davranışı — üçbucaqlaşdırma, TƏKRAR deyil."""
    scene, _ = _scene(slice_axis="X", slice_position=0.43)
    output = scene.slice_output()
    types = {output.GetCellType(i) for i in range(output.GetNumberOfCells())}
    assert types == {vtk.VTK_TRIANGLE}
    assert output.GetNumberOfCells() == 2 * len(_cut_cells(scene))


def test_plane_on_a_cell_boundary_still_cuts_one_layer():
    """X = 100 m dəqiq hüceyrə sərhədidir — iki qonşu sütun ikiqat sayılmamalıdır."""
    scene, _ = _scene(slice_axis="X", slice_position=0.5)
    assert len(_cut_cells(scene)) == NY * NZ


def test_slice_hides_the_volume_and_restores_it_when_off():
    scene, _ = _scene(slice_axis="X", slice_position=0.5)
    assert scene._actor.GetVisibility() == 0
    assert scene._slice_actor.GetVisibility() == 1

    scene.settings = VtkViewSettings(slice_axis=None)
    scene.update_values(np.arange(scene.model.ncell, dtype=float), "sınaq")
    assert scene._actor.GetVisibility() == 1
    assert scene.slice_output() is None


def test_slice_carries_the_cell_values_of_the_cut_cells():
    """Kəsiyin rəngi uydurma deyil — kəsilən hüceyrələrin öz dəyəridir."""
    scene, model = _scene(slice_axis="X", slice_position=0.43)
    assert _cut_cells(scene) <= set(np.arange(model.ncell, dtype=float))


def test_slice_uses_the_same_colour_table_as_the_volume():
    scene, _ = _scene(slice_axis="Y", slice_position=0.3)
    slice_mapper = scene._slice_actor.GetMapper()
    assert slice_mapper.GetLookupTable() is scene._mapper.GetLookupTable()
    assert slice_mapper.GetScalarRange() == scene._mapper.GetScalarRange()


def test_slice_respects_the_layer_filter():
    """Gizlədilmiş hüceyrələr kəsikdə də görünməməlidir (yalnız K = 1)."""
    scene, _ = _scene(slice_axis="X", slice_position=0.43, k_range=(0, 0))
    assert len(_cut_cells(scene)) == NY * 1


def test_moving_the_slice_changes_which_cells_are_cut():
    first, _ = _scene(slice_axis="X", slice_position=0.1)
    last, _ = _scene(slice_axis="X", slice_position=0.9)
    assert _cut_cells(first).isdisjoint(_cut_cells(last))


# ═══════════════════════ sürüklənən widget ═══════════════════════════

def _interactive(scene):
    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetSize(120, 120)
    window.AddRenderer(scene.renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    window.Render()
    return window, interactor


@pytest.mark.parametrize("axis, normal", [
    ("X", (1.0, 0.0, 0.0)), ("Y", (0.0, 1.0, 0.0)), ("Z", (0.0, 0.0, 1.0))])
def test_widget_attaches_offscreen_and_follows_the_axis(axis, normal):
    scene, _ = _scene(slice_axis=axis, slice_position=0.5)
    _window, interactor = _interactive(scene)
    scene.attach_slice_widget(interactor)
    assert scene._slice_widget is not None
    assert scene._slice_widget.GetEnabled() == 1
    assert scene._slice_widget.GetRepresentation().GetNormal() == \
        pytest.approx(normal)


def test_widget_normal_follows_an_axis_change():
    """Oxu Z-dən X-ə dəyişəndə widget də X-ə keçməlidir."""
    scene, _ = _scene(slice_axis="Z", slice_position=0.5)
    _window, interactor = _interactive(scene)
    scene.attach_slice_widget(interactor)
    scene.settings = VtkViewSettings(slice_axis="X", slice_position=0.5)
    scene.update_values(np.arange(scene.model.ncell, dtype=float), "sınaq")
    assert scene._slice_widget.GetRepresentation().GetNormal() == \
        pytest.approx((1.0, 0.0, 0.0))


def test_widget_is_off_when_the_slice_is_off():
    scene, _ = _scene()
    _window, interactor = _interactive(scene)
    scene.attach_slice_widget(interactor)
    assert scene._slice_widget.GetEnabled() == 0


def test_dragging_reports_the_new_fraction():
    """Sürükləmə → `on_moved(nisbət)` → interfeysin sürgüsü sinxron qalır."""
    scene, _ = _scene(slice_axis="X", slice_position=0.5)
    _window, interactor = _interactive(scene)
    reported = []
    scene.attach_slice_widget(interactor, reported.append)

    bounds = scene._grid.GetBounds()
    representation = scene._slice_widget.GetRepresentation()
    origin = list(representation.GetOrigin())
    origin[0] = bounds[0] + 0.8 * (bounds[1] - bounds[0])
    representation.SetOrigin(*origin)
    scene._slice_widget.InvokeEvent("InteractionEvent")

    assert reported and reported[-1] == pytest.approx(0.8, abs=1e-3)


def test_widget_is_skipped_without_an_interactor():
    scene, _ = _scene(slice_axis="X")
    scene.attach_slice_widget(None)
    assert scene._slice_widget is None


# ═══════════════════════ UI bağlantısı ═══════════════════════════════

def test_main_window_passes_slice_settings_to_the_scene():
    """Pəncərə QURULMADAN — `MainWindow()` testdə `pytest`-i çökdürür."""
    pytest.importorskip("PyQt5.QtWidgets")
    import inspect

    from imex2d.ui import main_window as MW
    source = inspect.getsource(MW.MainWindow)
    for fragment in ("volume_slice_axis", "volume_slice_position",
                     "slice_axis=self.volume_slice_axis.currentData()",
                     "attach_slice_widget", "_on_slice_dragged"):
        assert fragment in source, fragment
