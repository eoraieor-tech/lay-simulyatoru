"""B6-c — PNG kadr və animasiyalı GIF ixracı.

TAPILAN SƏHV (yol boyu): 3D tabındakı "Şəkli saxla…" HƏMİŞƏ matplotlib
fiqurunu yazırdı. VTK motoru aktiv olanda həmin fiqur gizli kanvasdadır —
saxlanılan şəkil istifadəçinin gördüyü VTK görüntüsü deyildi.

`MainWindow()` testdə `pytest`-i çökdürür (bax `test_playback_controls.py`),
ona görə interfeys hissəsi mənbə mətnindən yoxlanılır; ixracın özü isə
`rendering/animation_export.py`-də Qt-siz sınanır.
"""

import inspect

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.config import OutputConfig, SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.scal import GasCoreyParameters
from imex2d.rendering import animation_export as AE
from imex2d.rendering import renderers as R
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter

try:
    import vtk
    from imex2d.rendering import vtk_volume
    from PIL import Image
except ImportError as exc:                                   # pragma: no cover
    pytest.skip(f"VTK/Pillow əlçatan deyil: {exc}", allow_module_level=True)

WIDTH, HEIGHT = 160, 120


def _model(with_gas=False, nx=8, ny=6, nz=3):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=nz, top_depth=2000.0)
    return ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=default_scal(),
        gas_scal=GasCoreyParameters() if with_gas else None,
        pvt_table=build_pvt_table(pressure_min=1.0, pressure_max=400.0,
                                  n_points=40, bubble_point_bar=240.0,
                                  include_gas=with_gas))


def _offscreen(scene):
    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetSize(WIDTH, HEIGHT)
    window.AddRenderer(scene.renderer)
    scene.reset_camera()
    return window


def _frame(scene, window, values, caption=""):
    scene.update_values(values, "sınaq")
    scene.set_caption(caption)
    return AE.capture_render_window(window)


def _gif_frames(path):
    image = Image.open(path)
    count = 0
    durations = []
    try:
        while True:
            count += 1
            durations.append(image.info.get("duration"))
            image.seek(image.tell() + 1)
    except EOFError:
        pass
    return count, durations, image.info.get("loop")


# ═══════════════════════ VTK kadr tutma ══════════════════════════════

def test_vtk_capture_has_the_window_size_and_rgb_mode():
    scene = vtk_volume.VtkReservoirScene(_model())
    window = _offscreen(scene)
    frame = _frame(scene, window, np.arange(scene.model.ncell, dtype=float))
    assert frame.size == (WIDTH, HEIGHT)
    assert frame.mode == "RGB"


def test_different_values_give_different_frames():
    """Kadr cari dəyərləri göstərməlidir — köhnə buferi yox."""
    scene = vtk_volume.VtkReservoirScene(_model())
    window = _offscreen(scene)
    n = scene.model.ncell
    low = _frame(scene, window, np.zeros(n))
    high = _frame(scene, window, np.linspace(0.0, 1.0, n))
    assert not np.array_equal(np.asarray(low), np.asarray(high))


def test_capture_is_not_blank():
    scene = vtk_volume.VtkReservoirScene(_model())
    window = _offscreen(scene)
    pixels = np.asarray(_frame(scene, window,
                               np.linspace(0.0, 1.0, scene.model.ncell)))
    assert len(np.unique(pixels.reshape(-1, 3), axis=0)) > 10


# ═══════════════════════ kadr yazısı ═════════════════════════════════

def test_caption_is_set_and_shown():
    scene = vtk_volume.VtkReservoirScene(_model())
    scene.set_caption("t = 400 gün  ·  Sg")
    assert scene.caption() == "t = 400 gün  ·  Sg"
    assert scene._caption_actor.GetVisibility() == 1


def test_empty_caption_hides_the_text():
    scene = vtk_volume.VtkReservoirScene(_model())
    scene.set_caption("t = 1 gün")
    scene.set_caption("")
    assert scene.caption() == ""
    assert scene._caption_actor.GetVisibility() == 0


def test_caption_changes_the_captured_pixels():
    """Yazı səhnənin İÇİNDƏDİR — tutulan kadra həqiqətən düşür."""
    scene = vtk_volume.VtkReservoirScene(_model())
    window = _offscreen(scene)
    values = np.linspace(0.0, 1.0, scene.model.ncell)
    without = _frame(scene, window, values, caption="")
    with_text = _frame(scene, window, values, caption="t = 400 gün")
    assert not np.array_equal(np.asarray(without), np.asarray(with_text))


# ═══════════════════════ matplotlib kadr tutma ═══════════════════════

def test_figure_capture_matches_the_figure_size():
    from matplotlib.figure import Figure
    figure = Figure(figsize=(2.0, 1.5), dpi=80)
    figure.add_subplot().plot([0, 1], [0, 1])
    frame = AE.capture_figure(figure)
    assert frame.size == (160, 120)
    assert frame.mode == "RGB"


# ═══════════════════════ yazıcılar ═══════════════════════════════════

def _solid(colour, size=(WIDTH, HEIGHT)):
    return Image.new("RGB", size, colour)


def test_gif_round_trip_keeps_frames_duration_and_loop(tmp_path):
    frames = [_solid((i * 40, 80, 200 - i * 40)) for i in range(5)]
    path = AE.write_gif(frames, str(tmp_path / "a.gif"), frame_ms=140)
    count, durations, loop = _gif_frames(path)
    assert count == 5
    assert set(durations) == {140}
    assert loop == 0


def test_gif_frame_duration_has_a_floor(tmp_path):
    """<20 ms brauzerlərdə 100 ms kimi oynayır — ona görə aşağı salınmır."""
    frames = [_solid((0, 0, 0)), _solid((255, 255, 255))]
    path = AE.write_gif(frames, str(tmp_path / "b.gif"), frame_ms=5)
    _count, durations, _loop = _gif_frames(path)
    assert min(durations) >= AE.MIN_GIF_FRAME_MS


def test_gif_without_frames_is_rejected(tmp_path):
    with pytest.raises(ValueError):
        AE.write_gif([], str(tmp_path / "c.gif"), frame_ms=100)


def test_gif_with_mixed_frame_sizes_is_rejected(tmp_path):
    """Pəncərə ixrac zamanı dəyişsə səssizcə sıçrayan GIF yazılmamalıdır."""
    frames = [_solid((0, 0, 0)), _solid((0, 0, 0), size=(80, 60))]
    with pytest.raises(ValueError, match="ölçüsü"):
        AE.write_gif(frames, str(tmp_path / "d.gif"), frame_ms=100)


def test_png_is_written_and_readable(tmp_path):
    path = AE.write_png(_solid((10, 120, 30)), str(tmp_path / "kadr.png"))
    back = Image.open(path)
    assert back.size == (WIDTH, HEIGHT)
    assert back.convert("RGB").getpixel((5, 5)) == (10, 120, 30)


# ═══════════════════════ uc-uca: real simulyasiya → GIF ══════════════

def test_real_three_phase_run_exports_one_gif_frame_per_snapshot(tmp_path):
    """Hər snapshot bir kadr: kadr sayı = snapshot sayı, kadrlar fərqli."""
    model = _model(with_gas=True, nx=8, ny=6, nz=1)
    service = ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)
    result = service.run(model, SimulationConfig(
        end_time=400.0, output=OutputConfig(snapshot_count=6)))
    assert result.converged, result.message

    scene = vtk_volume.VtkReservoirScene(model)
    window = _offscreen(scene)
    frames = []
    for snapshot in result.snapshots:
        values, *_ = R.MapRenderer()._select_volume(
            model, R.SATURATION, snapshot)
        frames.append(_frame(scene, window, np.asarray(values).ravel(),
                             caption=f"t = {snapshot.time:.0f} gün"))

    path = AE.write_gif(frames, str(tmp_path / "suru.gif"), frame_ms=140)
    count, _durations, _loop = _gif_frames(path)
    assert count == len(result.snapshots)
    assert not np.array_equal(np.asarray(frames[0]), np.asarray(frames[-1]))


# ═══════════════════════ interfeys bağlantısı ═════════════════════════

def _source():
    pytest.importorskip("PyQt5.QtWidgets")
    from imex2d.ui import main_window as MW
    return inspect.getsource(MW.MainWindow)


def test_png_save_captures_the_active_engine():
    """SƏHVİN DÜZƏLİŞİ: VTK aktivdirsə kadr VTK pəncərəsindən tutulur."""
    source = _source()
    body = source.split("def save_volume_image", 1)[1].split("\n    def ", 1)[0]
    assert "_volume_uses_vtk()" in body
    assert "write_png" in body


def test_gif_export_is_wired_and_restores_the_slider():
    source = _source()
    assert "save_volume_animation" in source
    body = source.split("def save_volume_animation", 1)[1].split("\n    def ", 1)[0]
    assert "write_gif" in body
    assert "interval_ms(" in body, "kadr müddəti oynatma sürətindən gəlməlidir"
    assert "setValue(original)" in body, "sürgü sonda geri qaytarılmalıdır"


def test_caption_is_passed_to_the_vtk_scene():
    assert "set_caption(" in _source()
