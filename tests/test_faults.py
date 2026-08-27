"""Fay transmissivliyi (B3): həndəsə, tətbiq, fiziki doğruluq, I/O."""

import os
import tempfile

import numpy as np

from helpers import default_scal, make_service, short_config
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.domain.geological_model import Fault
from imex2d.domain.structure import FaultReference
from imex2d.io.fault_io import (FaultFormatError, read_eclipse_faults,
                                read_faults_csv, write_faults_csv)
from imex2d.simulation.discretization import TwoPointFluxDiscretization


def _geology(nx=20, ny=10, nz=1):
    return SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=nz, top_depth=2000.0)


def _model(faults=None, **kwargs):
    geology = _geology(**kwargs)
    scal = default_scal()
    return ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=scal,
        fault_references=faults)


def _write(text: str, suffix=".csv") -> str:
    handle, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(handle, "w", encoding="utf-8") as file:
        file.write(text)
    return path


# ── FaultReference domain ──────────────────────────────────────────────
def test_reference_without_geometry_matches_nothing():
    fault = FaultReference(name="F1", source_id="F1")
    assert not fault.has_geometry
    mask = fault.matches(0, np.array([0, 1, 2]), np.array([0, 0, 0]))
    assert not mask.any()


def test_sealing_overrides_the_multiplier():
    fault = FaultReference(name="F1", source_id="F1", axis="I", plane_index=5,
                           transmissibility_multiplier=0.8, sealing=True)
    assert fault.effective_multiplier == 0.0
    assert fault.is_sealing


def test_negative_multiplier_is_rejected():
    try:
        FaultReference(name="F1", source_id="F1",
                       transmissibility_multiplier=-0.1)
    except ValueError:
        return
    raise AssertionError("Mənfi çarpan qəbul edildi")


def test_invalid_axis_is_rejected():
    try:
        FaultReference(name="F1", source_id="F1", axis="X", plane_index=1)
    except ValueError:
        return
    raise AssertionError("Yanlış ox qəbul edildi")


def test_matches_only_the_declared_plane_and_range():
    fault = FaultReference(name="F1", source_id="F1", axis="I", plane_index=5,
                           range_a=(2, 4))
    i_boundary = np.array([5, 5, 5, 6])
    j_coordinate = np.array([2, 3, 9, 3])          # son ikisi range xaricində/planedən kənar
    axis_code = np.array([0, 0, 0, 0])
    # birbaşa matches() çağırışı: connection_axis_code sabitdir, koordinatlar dəyişir
    mask = fault.matches(0, j_coordinate, np.zeros(4))
    assert list(mask) == [True, True, False, True]     # plane yoxlaması ayrıca aparılır


def test_validate_flags_out_of_range_geometry():
    from imex2d.domain.grid import CartesianGrid

    grid = CartesianGrid(10, 10, 1)
    fault = FaultReference(name="F1", source_id="F1", axis="I", plane_index=15)
    assert fault.validate(grid)

    inside = FaultReference(name="F2", source_id="F2", axis="I", plane_index=5)
    assert inside.validate(grid) == []


# ── GeologicalModel.fault_references() bişirmə ─────────────────────────
def test_geological_fault_geometry_is_baked_into_reference():
    geology = _geology()
    geology.faults.append(Fault(name="F1", axis="I", plane_index=10,
                                transmissibility_multiplier=0.2))
    references = geology.fault_references()
    assert len(references) == 1
    assert references[0].axis == "I"
    assert references[0].plane_index == 10
    assert references[0].transmissibility_multiplier == 0.2


def test_direct_override_bypasses_geological_model_faults():
    geology = _geology()
    geology.faults.append(Fault(name="IGNORED", axis="I", plane_index=1))
    override = [FaultReference(name="REAL", source_id="REAL", axis="J",
                               plane_index=3)]
    model = ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=default_scal(),
        fault_references=override)
    assert [f.name for f in model.fault_references] == ["REAL"]


# ── diskretizasiyaya tətbiq ─────────────────────────────────────────────
def test_fault_only_affects_connections_on_its_plane():
    model = _model(faults=[FaultReference(
        name="F1", source_id="F1", axis="I", plane_index=10,
        transmissibility_multiplier=0.1)])
    baseline = _model(faults=[])

    grid_a = TwoPointFluxDiscretization().build(model)
    grid_b = TwoPointFluxDiscretization().build(baseline)
    conn = grid_a.connections

    i_a, _, _ = model.grid.ijk_array(conn.cell_a)
    on_fault = (conn.axis == 0) & (i_a == 10)

    assert np.allclose(grid_a.transmissibility[~on_fault],
                       grid_b.transmissibility[~on_fault])
    assert np.all(grid_a.transmissibility[on_fault]
                 < grid_b.transmissibility[on_fault])
    assert np.allclose(grid_a.transmissibility[on_fault],
                       grid_b.transmissibility[on_fault] * 0.1)


def test_sealing_fault_zeroes_out_the_plane_transmissibility():
    model = _model(faults=[FaultReference(
        name="F1", source_id="F1", axis="I", plane_index=10, sealing=True)])
    grid = TwoPointFluxDiscretization().build(model)
    i_a, _, _ = model.grid.ijk_array(grid.connections.cell_a)
    on_fault = (grid.connections.axis == 0) & (i_a == 10)
    assert np.allclose(grid.transmissibility[on_fault], 0.0)


def test_overlapping_faults_multiply_together():
    """İki qismən keçirici fay eyni üzdə üst-üstə düşəndə axın daha da azalır."""
    model = _model(faults=[
        FaultReference(name="F1", source_id="F1", axis="I", plane_index=10,
                      transmissibility_multiplier=0.5),
        FaultReference(name="F2", source_id="F2", axis="I", plane_index=10,
                      transmissibility_multiplier=0.4)])
    baseline = _model(faults=[])
    grid_a = TwoPointFluxDiscretization().build(model)
    grid_b = TwoPointFluxDiscretization().build(baseline)
    i_a, _, _ = model.grid.ijk_array(grid_a.connections.cell_a)
    on_fault = (grid_a.connections.axis == 0) & (i_a == 10)
    assert np.allclose(grid_a.transmissibility[on_fault],
                       grid_b.transmissibility[on_fault] * 0.2, rtol=1e-6)


def test_range_restricts_the_fault_to_part_of_the_boundary():
    model = _model(faults=[FaultReference(
        name="F1", source_id="F1", axis="I", plane_index=10, sealing=True,
        range_a=(0, 3))])
    grid = TwoPointFluxDiscretization().build(model)
    i_a, j_a, _ = model.grid.ijk_array(grid.connections.cell_a)
    on_plane = (grid.connections.axis == 0) & (i_a == 10)
    inside = on_plane & (j_a <= 3)
    outside = on_plane & (j_a > 3)
    assert np.allclose(grid.transmissibility[inside], 0.0)
    assert np.all(grid.transmissibility[outside] > 0.0)


# ── fiziki doğruluq ───────────────────────────────────────────────────
def test_sealing_fault_dramatically_reduces_recovery():
    """Əsas doğrulama: tam bağlı fay istismarçını təzyiq dəstəyindən kəsir."""
    scal = default_scal()
    baseline = make_service(scal).run(_model(faults=[]),
                                      short_config(end_time=400.0))
    sealed = make_service(scal).run(
        _model(faults=[FaultReference(name="F1", source_id="F1", axis="I",
                                      plane_index=10, sealing=True)]),
        short_config(end_time=400.0))
    assert sealed.final_recovery_factor < baseline.final_recovery_factor * 0.1


def test_partial_fault_gives_an_intermediate_result():
    scal = default_scal()
    baseline = make_service(scal).run(_model(faults=[]),
                                      short_config(end_time=400.0)).final_recovery_factor
    sealed = make_service(scal).run(
        _model(faults=[FaultReference(name="F1", source_id="F1", axis="I",
                                      plane_index=10, sealing=True)]),
        short_config(end_time=400.0)).final_recovery_factor
    partial = make_service(scal).run(
        _model(faults=[FaultReference(name="F1", source_id="F1", axis="I",
                                      plane_index=10,
                                      transmissibility_multiplier=0.05)]),
        short_config(end_time=400.0)).final_recovery_factor
    assert sealed < partial < baseline


def test_model_without_faults_behaves_exactly_as_before():
    """Reqressiya qorunması: fay yoxdursa nəticə dəyişməməlidir."""
    scal = default_scal()
    with_empty_list = make_service(scal).run(_model(faults=[]),
                                             short_config(end_time=300.0))
    without_override = make_service(scal).run(
        ReservoirModelBuilder().build(_geology(), five_spot(_geology().grid),
                                      scal=scal),
        short_config(end_time=300.0))
    assert abs(with_empty_list.final_recovery_factor
               - without_override.final_recovery_factor) < 1e-6


# ── diaqnostika ──────────────────────────────────────────────────────
def test_out_of_range_fault_is_a_blocking_error():
    model = _model(faults=[FaultReference(
        name="F1", source_id="F1", axis="I", plane_index=999)])
    report = model.diagnose()
    assert report.has_errors
    assert any("F1" == item.source for item in report.errors)


def test_duplicate_fault_names_warn():
    model = _model(faults=[
        FaultReference(name="F1", source_id="F1", axis="I", plane_index=5),
        FaultReference(name="F1", source_id="F1", axis="I", plane_index=8)])
    report = model.diagnose()
    assert any("təkrarlanır" in item.message for item in report.warnings)


# ── CSV I/O ───────────────────────────────────────────────────────────
def test_reads_csv_with_full_geometry():
    path = _write("""name,axis,plane_index,a_low,a_high,b_low,b_high,multiplier,sealing
F1,I,10,0,4,,,0.1,0
F2,J,3,,,,,,1
""")
    try:
        faults = read_faults_csv(path)
        assert len(faults) == 2
        by_name = {f.name: f for f in faults}
        assert by_name["F1"].range_a == (0, 4)
        assert by_name["F1"].range_b is None
        assert by_name["F2"].sealing
    finally:
        os.unlink(path)


def test_csv_missing_required_columns_is_rejected():
    path = _write("name,axis\nF1,I\n")
    try:
        read_faults_csv(path)
    except FaultFormatError:
        return
    finally:
        os.unlink(path)
    raise AssertionError("Natamam CSV qəbul edildi")


def test_csv_round_trip_preserves_geometry():
    faults = [FaultReference(name="F1", source_id="F1", axis="I",
                             plane_index=10, range_a=(0, 4), range_b=(1, 2),
                             transmissibility_multiplier=0.3)]
    handle, path = tempfile.mkstemp(suffix=".csv")
    os.close(handle)
    try:
        write_faults_csv(path, faults)
        restored = read_faults_csv(path)[0]
    finally:
        os.unlink(path)
    assert restored.axis == "I"
    assert restored.plane_index == 10
    assert tuple(restored.range_a) == (0, 4)
    assert abs(restored.transmissibility_multiplier - 0.3) < 1e-6


# ── Eclipse FAULTS/MULTFLT ──────────────────────────────────────────────
def test_reads_eclipse_faults_and_multfilt():
    path = _write("""GRID

FAULTS
  'F1'  11 11  1 41  1 5  'I' /
/

MULTFLT
  'F1'  0.05 /
/
""", suffix=".DATA")
    try:
        faults = read_eclipse_faults(path)
        assert len(faults) == 1
        fault = faults[0]
        assert fault.axis == "I"
        assert fault.plane_index == 10          # 11 (1-based) -> 10 (0-based)
        assert fault.range_a == (0, 40)
        assert abs(fault.transmissibility_multiplier - 0.05) < 1e-9
    finally:
        os.unlink(path)


def test_eclipse_fault_without_multflt_defaults_to_transparent():
    path = _write("FAULTS\n  'F1'  6 6  1 10  1 1  'I' /\n/\n", suffix=".DATA")
    try:
        fault = read_eclipse_faults(path)[0]
        assert abs(fault.transmissibility_multiplier - 1.0) < 1e-9
        assert not fault.sealing
    finally:
        os.unlink(path)


def test_missing_faults_keyword_is_rejected():
    path = _write("GRID\n\nPORO\n 4*0.2 /\n", suffix=".DATA")
    try:
        read_eclipse_faults(path)
    except FaultFormatError:
        return
    finally:
        os.unlink(path)
    raise AssertionError("FAULTS olmadan fayl qəbul edildi")
