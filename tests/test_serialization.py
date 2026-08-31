"""Layihə faylı (.imx) — B1 testləri."""

import gzip
import json
import os
import tempfile

import numpy as np

from helpers import default_scal, five_spot_model, make_service, short_config
from imex2d.application.config import OutputConfig, SimulationConfig
from imex2d.application.project import Project
from imex2d.application.serialization import (FORMAT_VERSION, ProjectFileError,
                                              ProjectSerializer)
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder
from imex2d.domain.wells import Perforation
from imex2d.simulation.pvt.correlations import build_pvt_table


def _rich_project(with_result=True, nz=1):
    """Bütün xüsusiyyətləri olan layihə: PVT, Pc, 3D, perforasiya, nəticə."""
    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=9, ny=9, dx=25.0, dy=25.0, dz=6.0, porosity=0.21,
        permx_base=180.0, nz=nz, top_depth=2000.0, heterogeneous=True,
        sigma=0.6, seed=11)
    model = five_spot_model(nx=9, ny=9, scal=scal)
    model.rock = model.rock
    model.pvt_table = build_pvt_table(bubble_point_bar=230.0)
    model.capillary_parameters.entry_pressure = 0.35
    model.wells[1].perforations = [Perforation(8, 8, 0)]

    project = Project("Serializasiya testi")
    project.add_geological_model(geology)
    project.add_reservoir_model(model)
    config = SimulationConfig(end_time=100.0,
                              output=OutputConfig(snapshot_count=4))
    run = project.new_run(model.name, config)
    if with_result:
        run.result = make_service(scal).run(model, config)
        run.status = "FINISHED"
    return project, model


def _round_trip(project, include_snapshots=True):
    serializer = ProjectSerializer()
    handle, path = tempfile.mkstemp(suffix=".imx")
    os.close(handle)
    try:
        serializer.save(project, path, include_snapshots)
        return serializer.load(path), os.path.getsize(path)
    finally:
        os.unlink(path)


# ── format ────────────────────────────────────────────────────────────
def test_saved_file_is_gzipped_json_with_version():
    project, _ = _rich_project(with_result=False)
    handle, path = tempfile.mkstemp(suffix=".imx")
    os.close(handle)
    try:
        ProjectSerializer().save(project, path)
        with gzip.open(path, "rt", encoding="utf-8") as file:
            payload = json.load(file)
        assert payload["version"] == FORMAT_VERSION
        assert "project" in payload and "saved_at" in payload
    finally:
        os.unlink(path)


def test_unknown_version_is_rejected():
    handle, path = tempfile.mkstemp(suffix=".imx")
    os.close(handle)
    try:
        with gzip.open(path, "wt", encoding="utf-8") as file:
            json.dump({"version": 999, "project": {}}, file)
        try:
            ProjectSerializer().load(path)
        except ProjectFileError:
            return
        raise AssertionError("Naməlum versiya qəbul edildi")
    finally:
        os.unlink(path)


def test_corrupt_file_raises_project_file_error():
    handle, path = tempfile.mkstemp(suffix=".imx")
    os.close(handle)
    try:
        with open(path, "w", encoding="utf-8") as file:
            file.write("bu JSON deyil")
        try:
            ProjectSerializer().load(path)
        except ProjectFileError:
            return
        raise AssertionError("Pozulmuş fayl qəbul edildi")
    finally:
        os.unlink(path)


# ── model bərpası ─────────────────────────────────────────────────────
def test_grid_and_geometry_survive_round_trip():
    project, model = _rich_project(with_result=False, nz=3)
    restored = _round_trip(project)[0].reservoir_models[model.name]
    assert restored.grid.nx == model.grid.nx
    assert restored.grid.nz == model.grid.nz
    assert restored.geometry.dx == model.geometry.dx
    assert restored.ncell == model.ncell


def test_variable_layer_thickness_survives_round_trip():
    """Hər təbəqənin öz DZ-i saxlanılmalıdır — tək orta ədədə sıxılmamalı."""
    from imex2d.application.model_builder import ReservoirModelBuilder
    from imex2d.application.scenarios import five_spot

    dz = [4.0, 6.0, 10.0]
    geology = SyntheticGeologicalModelBuilder().build(
        nx=5, ny=5, dx=25.0, dy=25.0, dz=dz, porosity=0.2,
        permx_base=150.0, nz=3, top_depth=2000.0)
    model = ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=default_scal(),
        name="Dəyişən DZ testi")

    project = Project("Dəyişən DZ layihəsi")
    project.add_geological_model(geology)
    project.add_reservoir_model(model)

    restored = _round_trip(project, include_snapshots=False)[0]
    restored_model = restored.reservoir_models[model.name]
    assert np.allclose(restored_model.geometry.dz, dz)
    assert not np.allclose(restored_model.geometry.dz, restored_model.geometry.dz[0])


def test_heterogeneous_properties_survive_round_trip():
    project, model = _rich_project(with_result=False)
    restored = _round_trip(project)[0].reservoir_models[model.name]
    assert np.allclose(restored.rock.permx.values, model.rock.permx.values)
    assert np.allclose(restored.rock.porosity.values, model.rock.porosity.values)


def test_wells_and_perforations_survive_round_trip():
    project, model = _rich_project(with_result=False)
    restored = _round_trip(project)[0].reservoir_models[model.name]
    assert [w.name for w in restored.wells] == [w.name for w in model.wells]
    for original, copy in zip(model.wells, restored.wells):
        assert original.well_type is copy.well_type
        assert original.control.mode is copy.control.mode
        assert original.control.target == copy.control.target
        assert [(p.i, p.j, p.k) for p in original.perforations] == \
               [(p.i, p.j, p.k) for p in copy.perforations]


def test_pvt_scal_capillary_and_initial_conditions_survive():
    project, model = _rich_project(with_result=False)
    restored = _round_trip(project)[0].reservoir_models[model.name]
    assert restored.pvt_table is not None
    assert restored.pvt_table.bubble_point == model.pvt_table.bubble_point
    assert np.allclose(restored.pvt_table.oil_fvf, model.pvt_table.oil_fvf)
    assert restored.capillary_parameters.entry_pressure == 0.35
    assert restored.scal_parameters.swc == model.scal_parameters.swc
    assert restored.initial_conditions.datum_pressure == \
           model.initial_conditions.datum_pressure


def test_restored_model_passes_validation():
    project, model = _rich_project(with_result=False)
    restored = _round_trip(project)[0].reservoir_models[model.name]
    assert restored.validate() == []


def test_restored_model_reproduces_the_same_result():
    """Ən vacib test: fayldan açılan model eyni rəqəmi verməlidir."""
    scal = default_scal()
    project, model = _rich_project(with_result=False)
    config = SimulationConfig(end_time=100.0,
                              output=OutputConfig(snapshot_count=3))
    original = make_service(scal).run(model, config)

    restored = _round_trip(project)[0].reservoir_models[model.name]
    repeated = make_service(scal).run(restored, config)

    assert abs(repeated.final_recovery_factor
               - original.final_recovery_factor) < 1e-9
    assert repeated.steps == original.steps


# ── nəticələr ─────────────────────────────────────────────────────────
def test_results_and_snapshots_survive_round_trip():
    project, model = _rich_project(with_result=True)
    restored = _round_trip(project)[0]
    run = restored.runs["RUN-001"]
    original = project.runs["RUN-001"].result

    assert run.status == "FINISHED"
    assert abs(run.result.ooip - original.ooip) < 1e-6
    assert run.result.steps == original.steps
    assert np.allclose(run.result.series.recovery_factor,
                       original.series.recovery_factor)
    assert len(run.result.snapshots) == len(original.snapshots)
    assert np.allclose(run.result.snapshots[-1].water_saturation,
                       original.snapshots[-1].water_saturation)
    assert run.result.snapshots[-1].water_saturation.shape == model.grid.shape


def test_snapshots_can_be_excluded_to_shrink_the_file():
    project, _ = _rich_project(with_result=True)
    with_snapshots = _round_trip(project, include_snapshots=True)
    without = _round_trip(project, include_snapshots=False)
    assert without[1] < with_snapshots[1]
    assert without[0].runs["RUN-001"].result.snapshots == []
    # zaman sıraları hər halda qalır
    assert without[0].runs["RUN-001"].result.series.time


def test_run_configuration_survives_round_trip():
    project, _ = _rich_project(with_result=True)
    restored = _round_trip(project)[0].runs["RUN-001"]
    original = project.runs["RUN-001"]
    assert restored.config.end_time == original.config.end_time
    assert restored.config.time_stepping.cfl_factor == \
           original.config.time_stepping.cfl_factor
    assert restored.config.output.snapshot_count == \
           original.config.output.snapshot_count


def test_project_counter_survives_so_new_runs_do_not_collide():
    project, model = _rich_project(with_result=True)
    restored = _round_trip(project)[0]
    new_run = restored.new_run(model.name, short_config())
    assert new_run.run_id not in {"RUN-001"}
    assert new_run.run_id == "RUN-002"
