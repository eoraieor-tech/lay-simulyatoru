"""Seans 46 — iş seansı məntiqi (`application/session.py`).

`MainWindow` pytest-də qurula bilmədiyi üçün (VTK) qərar məntiqi Qt-siz
modulda saxlanılır. Pəncərə ilə birlikdə real proqramda yoxlanıb:
saxla → yeni proses → aç → yenidən hesabla (bit-bit eyni), qəza →
bərpa, bərpadan imtina, saxlanmamış işlə bağlama (bax ISH_HESABATI →
Seans 46).
"""

from __future__ import annotations

import os

from imex2d.application import session


def test_push_recent_puts_newest_first_and_limits(tmp_path):
    paths = []
    for number in range(7):
        paths = session.push_recent(paths, str(tmp_path / f"p{number}.imx"))
    assert len(paths) == session.RECENT_LIMIT
    assert os.path.basename(paths[0]) == "p6.imx"
    assert os.path.basename(paths[-1]) == "p2.imx"


def test_push_recent_moves_existing_entry_to_top(tmp_path):
    a, b = str(tmp_path / "a.imx"), str(tmp_path / "b.imx")
    paths = session.push_recent(session.push_recent([], a), b)
    paths = session.push_recent(paths, a)
    assert [os.path.basename(p) for p in paths] == ["a.imx", "b.imx"]


def test_same_file_written_differently_is_one_entry(tmp_path):
    path = str(tmp_path / "Layihe.imx")
    other_spelling = os.path.join(str(tmp_path), ".", "Layihe.imx")
    paths = session.push_recent(session.push_recent([], path), other_spelling)
    assert len(paths) == 1


def test_remove_recent(tmp_path):
    a, b = str(tmp_path / "a.imx"), str(tmp_path / "b.imx")
    assert session.remove_recent([a, b], a) == [b]


def test_recovery_path_follows_data_dir_override(tmp_path, monkeypatch):
    monkeypatch.setenv(session.DATA_DIR_ENV, str(tmp_path))
    assert session.recovery_path() == str(tmp_path / session.RECOVERY_FILE_NAME)
    assert session.is_recovery_file(str(tmp_path / session.RECOVERY_FILE_NAME))
    assert not session.is_recovery_file(str(tmp_path / "layihe.imx"))
    assert not session.is_recovery_file(None)


def test_recovery_file_is_never_a_recent_project(tmp_path, monkeypatch):
    """Bərpa faylı istifadəçinin faylı deyil — `_remember_recent` onu atır;
    burada yalnız tanınmasının əsası yoxlanılır."""
    monkeypatch.setenv(session.DATA_DIR_ENV, str(tmp_path))
    assert session.is_recovery_file(session.recovery_path())


def test_signature_changes_only_when_work_changes():
    ui = {"panels": {"numerical": {"end_time": 1500.0}}}
    base = session.state_signature(ui, ["RUN-001"], [{"name": "W1"}])
    assert base == session.state_signature(
        {"panels": {"numerical": {"end_time": 1500.0}}}, ["RUN-001"],
        [{"name": "W1"}])
    assert base != session.state_signature(
        {"panels": {"numerical": {"end_time": 900.0}}}, ["RUN-001"],
        [{"name": "W1"}])
    assert base != session.state_signature(ui, ["RUN-001", "RUN-002"],
                                           [{"name": "W1"}])
    assert base != session.state_signature(ui, ["RUN-001"], [])
    # işə salınmaların sırası izi dəyişmir
    assert (session.state_signature(ui, ["B", "A"], [])
            == session.state_signature(ui, ["A", "B"], []))
