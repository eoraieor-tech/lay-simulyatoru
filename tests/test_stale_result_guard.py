"""Seans 53 — X-1: nəticə ekranda ikən grid ölçüsü dəyişəndə çökmə.

ƏVVƏLKİ VƏZİYYƏT (Seans 52-də ölçüldü). Hesablamadan sonra NX/NY/NZ
dəyişdiriləndə, başqa ölçülü layihə açılanda və ya GRDECL idxal olunanda
proses dərhal dayanırdı (çıxış kodu `0xC0000409`, 5 ssenaridən 4-ü):

    main_window.rebuild_model → update_map → MapRenderer._select
    ValueError: cannot reshape array of size 1681 into shape (1,41,40)

İki səbəb: (1) köhnə nəticənin anı yeni gridin formasına salınırdı;
(2) `sys.excepthook` standart idi — PyQt5 slotdakı istisnada `qFatal`
çağırır. İstifadəçi mesaj görmürdü, loga heç nə yazılmırdı.

DÜZƏLİŞ. (1) `SimulationResult.fits_grid` — anlar yalnız öz gridinə
çəkilir (`MainWindow._snapshot_at`); (2) `ui/error_hook.py` — tutulmamış
istisna loga yazılır, istifadəçiyə göstərilir, proses DAVAM EDİR.

`MainWindow` burada qurulmur (VTK + kanvaslar pytest-i çökdürür —
bax `test_playback_controls.py`); pəncərə əlaqəsi mənbə mətnindən,
PyQt davranışı isə ayrıca prosesdə yoxlanılır.
"""

from __future__ import annotations

import inspect
import logging
import os
import subprocess
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest
from matplotlib.figure import Figure

from helpers import five_spot_model, make_service, short_config
from imex2d.application.project import Project
from imex2d.application.serialization import ProjectSerializer
from imex2d.rendering import renderers as R
from imex2d.simulation.results import SimulationResult, Snapshot
from imex2d.ui import error_hook

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _result(shape, with_shape=True):
    ncell = int(np.prod(shape))
    return SimulationResult(
        grid_shape=tuple(shape) if with_shape else (),
        snapshots=[Snapshot(time=1.0, pressure=np.full(ncell, 200.0),
                            water_saturation=np.full(ncell, 0.3))])


# ═══════════════════════ SimulationResult.fits_grid ══════════════════

def test_result_fits_its_own_grid():
    assert _result((1, 41, 41)).fits_grid((1, 41, 41))


@pytest.mark.parametrize("new_shape", [(1, 41, 40), (2, 41, 41), (1, 50, 50),
                                       (5, 25, 30)])
def test_result_does_not_fit_a_resized_grid(new_shape):
    """Seans 52-nin ssenariləri: NX, NZ, başqa layihə, GRDECL."""
    assert not _result((1, 41, 41)).fits_grid(new_shape)


def test_same_cell_count_but_different_shape_does_not_fit():
    """40×41 ↔ 41×40: `reshape` ÇÖKMƏZDİ, amma xəritə yanlış olardı."""
    assert not _result((1, 41, 40)).fits_grid((1, 40, 41))


def test_result_without_snapshots_fits_any_grid():
    """Yalnız sıralar (əyrilər) gridə bağlı deyil."""
    assert SimulationResult(grid_shape=(1, 41, 41)).fits_grid((3, 7, 9))


def test_result_without_recorded_shape_falls_back_to_cell_count():
    old = _result((1, 41, 41), with_shape=False)
    assert old.fits_grid((1, 41, 41))
    assert not old.fits_grid((1, 41, 40))


def test_engine_result_records_the_shape_of_its_grid():
    model = five_spot_model(nx=9, ny=7)
    result = make_service().run(model, short_config(end_time=30.0))
    assert result.snapshots
    assert result.fits_grid(model.grid.shape)
    assert not result.fits_grid((1, 7, 10))


def test_shape_survives_the_project_file(tmp_path):
    """Açılan layihənin nəticəsi də öz gridini tanımalıdır."""
    model = five_spot_model(nx=9, ny=7)
    result = make_service().run(model, short_config(end_time=30.0))
    project = Project("X-1")
    project.add_reservoir_model(model)
    run = project.new_run(model.name, short_config(end_time=30.0))
    run.result = result
    run.status = "FINISHED"
    path = str(tmp_path / "x1.imx")
    ProjectSerializer().save(project, path)
    back = ProjectSerializer().load(path).latest_run().result
    assert back.fits_grid(model.grid.shape)
    assert not back.fits_grid((1, 41, 41))


def test_mismatched_snapshot_is_exactly_what_crashed_the_map():
    """Kök səbəbin sənədi: uyğunsuz an rendererdə `ValueError` verir,
    `None` (statik xassə) isə çəkilir — `_snapshot_at` məhz bunu seçir."""
    model = five_spot_model(nx=9, ny=7)
    stale = _result((1, 9, 9)).snapshots[0]
    figure = Figure()
    axes = figure.subplots()
    with pytest.raises(ValueError):
        R.MapRenderer().draw(axes, figure, model, R.SATURATION, stale)
    R.MapRenderer().draw(axes, figure, model, R.SATURATION, None)


# ═══════════════════════ MainWindow əlaqəsi (mənbədən) ════════════════

def _source(name):
    pytest.importorskip("PyQt5.QtWidgets")
    from imex2d.ui import main_window as MW
    return inspect.getsource(getattr(MW.MainWindow, name))


@pytest.mark.parametrize("method", ["update_map", "update_volume"])
def test_views_take_snapshots_only_through_the_guard(method):
    source = _source(method)
    assert "self._snapshot_at(" in source
    assert "self.result.snapshots[" not in source


def test_guard_checks_the_grid_of_the_current_model():
    source = _source("_snapshot_at")
    assert "_result_fits_model" in source
    assert "fits_grid" in _source("_result_fits_model")


@pytest.mark.parametrize("method", ["export_snapshot", "save_volume_animation"])
def test_snapshot_exports_refuse_a_stale_result(method):
    assert "_result_fits_model" in _source(method)


def test_opened_project_sets_up_the_3d_time_slider_like_a_finished_run():
    """Seans 53-də tapıldı: açılan layihədə 3D sürgüsü qurulmurdu."""
    source = _source("open_project")
    for call in ("self.volume_time.setEnabled(True)",
                 "self.volume_time.setRange(", "self.volume_time.setValue("):
        assert call in source, call


def test_rebuild_tells_the_user_why_the_map_is_static():
    assert "STALE_RESULT_MESSAGE" in _source("rebuild_model")


# ═══════════════════════ qlobal xəta tutucusu ════════════════════════

class _Recorder(logging.Handler):
    def __init__(self):
        super().__init__()
        self.records = []

    def emit(self, record):
        self.records.append(record)


def _logger():
    logger = logging.getLogger("imex2d.test_error_hook")
    logger.propagate = False
    recorder = _Recorder()
    logger.handlers = [recorder]
    return logger, recorder


def _raise_and_hook(hook, exc):
    try:
        raise exc
    except BaseException:
        hook(*sys.exc_info())


def test_hook_logs_the_full_traceback_and_notifies_the_user():
    logger, recorder = _logger()
    shown = []
    hook = error_hook.make_hook(logger, shown.append)
    _raise_and_hook(hook, ValueError("cannot reshape array of size 1681"))

    assert len(recorder.records) == 1
    record = recorder.records[0]
    assert record.levelno == logging.CRITICAL
    assert record.exc_info and record.exc_info[0] is ValueError
    assert len(shown) == 1
    assert "ValueError" in shown[0] and "1681" in shown[0]
    assert error_hook.LOG_HINT in shown[0]


def test_failing_notification_is_logged_not_raised():
    logger, recorder = _logger()

    def broken(_text):
        raise RuntimeError("mesaj qutusu açılmadı")

    hook = error_hook.make_hook(logger, broken)
    _raise_and_hook(hook, KeyError("x"))          # istisna çölə çıxmamalıdır
    assert [r.levelno for r in recorder.records] == [logging.CRITICAL,
                                                     logging.ERROR]


def test_exception_while_notifying_does_not_recurse():
    """Mesaj göstərilərkən yeni istisna — yalnız loga, ikinci qutu yox."""
    logger, recorder = _logger()
    shown = []

    def nested(text):
        shown.append(text)
        _raise_and_hook(hook, ZeroDivisionError("iç xəta"))

    hook = error_hook.make_hook(logger, nested)
    _raise_and_hook(hook, ValueError("xarici xəta"))
    assert len(shown) == 1
    assert sum(r.levelno == logging.CRITICAL for r in recorder.records) == 2


def test_keyboard_interrupt_goes_to_the_default_hook(monkeypatch):
    logger, recorder = _logger()
    seen = []
    monkeypatch.setattr(sys, "__excepthook__", lambda *args: seen.append(args[0]))
    shown = []
    _raise_and_hook(error_hook.make_hook(logger, shown.append), KeyboardInterrupt())
    assert seen == [KeyboardInterrupt]
    assert not shown and not recorder.records


def test_install_replaces_and_returns_the_previous_hook(monkeypatch):
    logger, _ = _logger()
    sentinel = object()
    monkeypatch.setattr(sys, "excepthook", sentinel)
    previous = error_hook.install(logger)
    assert previous is sentinel
    assert callable(sys.excepthook) and sys.excepthook is not sentinel


def test_app_installs_the_hook():
    with open(os.path.join(ROOT, "app.py"), encoding="utf-8") as handle:
        source = handle.read()
    assert "error_hook.install(" in source


# ═══════════════════════ PyQt davranışı (ayrıca prosesdə) ═════════════

_SLOT_SCRIPT = textwrap.dedent("""
    import logging, sys
    sys.path.insert(0, {root!r})
    from PyQt5.QtCore import QTimer
    from PyQt5.QtWidgets import QApplication
    from imex2d.ui import error_hook
    app = QApplication(sys.argv)
    if {install}:
        error_hook.install(logging.getLogger("imex2d.subprocess"),
                           notify=lambda text: print("MESAJ:", text.splitlines()[2]))
    def slot():
        raise ValueError("cannot reshape array of size 1681 into shape (1,41,40)")
    QTimer.singleShot(0, slot)
    QTimer.singleShot(200, app.quit)
    app.exec_()
    print("PROSES SAĞ QALDI")
""")


def _run_slot_script(install: bool):
    pytest.importorskip("PyQt5.QtWidgets")
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen")
    return subprocess.run(
        [sys.executable, "-c", _SLOT_SCRIPT.format(root=ROOT, install=install)],
        capture_output=True, text=True, encoding="utf-8", errors="replace",
        env=env, timeout=120)


def test_slot_exception_no_longer_kills_the_process():
    done = _run_slot_script(install=True)
    assert done.returncode == 0, done.stderr[-2000:]
    assert "PROSES SAĞ QALDI" in done.stdout
    assert "MESAJ: ValueError" in done.stdout


def test_without_the_hook_pyqt_aborts():
    """Nəzarət sınağı: tutucusuz eyni skript prosesi DAYANDIRIR —
    yuxarıdakı sınağın mənalı olduğunun sübutu (Seans 52-dəki çökmə)."""
    done = _run_slot_script(install=False)
    assert done.returncode != 0
    assert "PROSES SAĞ QALDI" not in done.stdout
