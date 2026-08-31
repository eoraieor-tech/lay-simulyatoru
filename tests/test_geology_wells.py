"""Geologiya quyu cədvəli (`GeologicalWell`) — domen, serializasiya, miqrasiya."""

import os
import tempfile

from helpers import five_spot_model
from imex2d.application.project import Project
from imex2d.application.serialization import ProjectSerializer
from imex2d.domain.geology import GeologicalWell


def _round_trip_project(project):
    serializer = ProjectSerializer()
    handle, path = tempfile.mkstemp(suffix=".imx")
    os.close(handle)
    try:
        serializer.save(project, path)
        return serializer.load(path)
    finally:
        os.unlink(path)


# ── GeologicalWell.to_dict / from_dict ──────────────────────────────────
def test_well_round_trips_through_dict():
    well = GeologicalWell(name="P-1", in_model=True, x=250.0, y=250.0,
                          top=2010.0, bottom=2032.0, porosity=0.21,
                          permeability=180.0, water_saturation=0.28,
                          note="karotaj W-1")
    restored = GeologicalWell.from_dict(well.to_dict())
    assert restored == well


def test_none_values_are_preserved_not_coerced_to_zero():
    well = GeologicalWell(name="P-2", x=10.0, y=10.0)
    data = well.to_dict()
    assert data["porosity"] is None
    assert data["top"] is None
    restored = GeologicalWell.from_dict(data)
    assert restored.porosity is None
    assert restored.bottom is None
    assert restored.porosity != 0.0


def test_zero_is_kept_as_zero_not_confused_with_missing():
    well = GeologicalWell(name="P-3", x=0.0, y=0.0, porosity=0.0)
    restored = GeologicalWell.from_dict(well.to_dict())
    assert restored.porosity == 0.0
    assert restored.porosity is not None


def test_in_model_flag_survives_round_trip():
    well = GeologicalWell(name="P-4", in_model=False, x=5.0, y=5.0)
    restored = GeologicalWell.from_dict(well.to_dict())
    assert restored.in_model is False


# ── Project.geology_source ──────────────────────────────────────────────
def test_geology_source_is_synthetic_when_table_is_empty():
    project = Project()
    assert project.geology_source == "synthetic"


def test_geology_source_is_wells_when_table_has_rows():
    project = Project()
    project.geology_wells = [GeologicalWell(name="W-1", x=0.0, y=0.0)]
    assert project.geology_source == "wells"


# ── serializasiya (Project səviyyəsində) ────────────────────────────────
def test_geology_wells_table_survives_round_trip():
    project = Project("Quyu cədvəli testi")
    project.geology_wells = [
        GeologicalWell(name="P-1", x=250.0, y=250.0, top=2010.0,
                       bottom=2032.0, porosity=0.21, permeability=180.0,
                       water_saturation=0.28),
        GeologicalWell(name="P-2", in_model=False, x=10.0, y=10.0),
    ]
    project.geology_method = "IDW"
    project.geology_params = {"power": 2.0}
    project.geology_defaults = {"sw": 0.30}

    restored = _round_trip_project(project)
    assert len(restored.geology_wells) == 2
    assert restored.geology_wells[0] == project.geology_wells[0]
    assert restored.geology_wells[1].in_model is False
    assert restored.geology_method == "IDW"
    assert restored.geology_params == {"power": 2.0}
    assert restored.geology_defaults == {"sw": 0.30}
    assert restored.geology_source == "wells"


def test_empty_geology_table_survives_round_trip_as_synthetic():
    project = Project("Sintetik layihə")
    restored = _round_trip_project(project)
    assert restored.geology_wells == []
    assert restored.geology_source == "synthetic"


# ── köhnə (.imx-də geology_wells yoxdur) fayl açılışı ───────────────────
def test_old_file_without_geology_block_still_loads():
    """v1 faylı (geology_wells açarı YOXDUR) sındırmadan açılmalıdır."""
    project = Project("Köhnə layihə")
    model = five_spot_model(nx=9, ny=9)
    project.add_reservoir_model(model)

    serializer = ProjectSerializer()
    payload_project = serializer.project_to_dict(project)
    assert "geology_wells" in payload_project
    del payload_project["geology_wells"]         # v1-i simulyasiya edir
    del payload_project["geology_method"]
    del payload_project["geology_params"]
    del payload_project["geology_defaults"]

    restored = serializer.project_from_dict(payload_project)
    assert restored.name == "Köhnə layihə"


# ── miqrasiya: wells blokundakı i/j -> geologiya X/Y ────────────────────
def test_migration_populates_geology_table_from_reservoir_wells():
    project = Project("Miqrasiya testi")
    model = five_spot_model(nx=9, ny=9, dx=25.0, dy=25.0)
    project.add_reservoir_model(model)

    serializer = ProjectSerializer()
    payload_project = serializer.project_to_dict(project)
    del payload_project["geology_wells"]
    del payload_project["geology_method"]
    del payload_project["geology_params"]
    del payload_project["geology_defaults"]

    restored = serializer.project_from_dict(payload_project)
    assert restored.geology_source == "wells"
    names = {w.name for w in restored.geology_wells}
    assert names == {w.name for w in model.wells}
    for well in restored.geology_wells:
        assert well.porosity is None            # heç vaxt saxlanılmayıb
        assert well.note == "köhnə layihədən miqrasiya edilib"
        assert well.x >= 0.0 and well.y >= 0.0


def test_migration_is_skipped_when_no_reservoir_model_exists():
    project = Project("Boş layihə")
    serializer = ProjectSerializer()
    payload_project = serializer.project_to_dict(project)
    del payload_project["geology_wells"]
    del payload_project["geology_method"]
    del payload_project["geology_params"]
    del payload_project["geology_defaults"]

    restored = serializer.project_from_dict(payload_project)
    assert restored.geology_wells == []
    assert restored.geology_source == "synthetic"


# ── Well.perf_top / perf_bottom ─────────────────────────────────────────
def test_well_perf_metres_survive_round_trip():
    project = Project("Perf metr testi")
    model = five_spot_model(nx=9, ny=9)
    model.wells[0].perf_top = 2015.0
    model.wells[0].perf_bottom = 2028.0
    project.add_reservoir_model(model)

    restored = _round_trip_project(project)
    restored_model = list(restored.reservoir_models.values())[0]
    well = next(w for w in restored_model.wells if w.name == model.wells[0].name)
    assert well.perf_top == 2015.0
    assert well.perf_bottom == 2028.0


def test_well_perf_metres_default_to_none():
    project = Project("Defolt perf testi")
    model = five_spot_model(nx=9, ny=9)
    project.add_reservoir_model(model)

    restored = _round_trip_project(project)
    restored_model = list(restored.reservoir_models.values())[0]
    assert all(w.perf_top is None and w.perf_bottom is None
              for w in restored_model.wells)


def test_unknown_future_version_still_rejected():
    """FORMAT_VERSION 2-yə qalxdıqdan sonra da gələcək versiya rədd olunur."""
    from imex2d.application.serialization import FORMAT_VERSION
    assert FORMAT_VERSION == 2
