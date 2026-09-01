"""Phase 4.1 §4 — sərt fasiya datasının hüceyrə xəritələnməsi + ziddiyyət."""

from __future__ import annotations

import pytest

from imex2d.domain.geometry import CellGeometry
from imex2d.domain.grid import CartesianGrid
from imex2d.domain.well_data import WellSample
from imex2d.geology.hard_data import (HardDataConflictError, detect_hard_data_conflicts,
                                      map_samples_to_cells, resolve_hard_data)


def _geometry(nx=5, ny=5, nz=2, dx=20.0, dy=20.0, dz=10.0, top=2000.0):
    grid = CartesianGrid(nx, ny, nz)
    return grid, CellGeometry(grid, dx, dy, dz, top_depth=top)


def test_no_conflict_when_samples_are_in_different_cells():
    grid, geometry = _geometry()
    samples = [
        WellSample(well="A", x=10.0, y=10.0, values={"FACIES": 0.0}, layer=0),
        WellSample(well="B", x=90.0, y=90.0, values={"FACIES": 1.0}, layer=0),
    ]
    assert detect_hard_data_conflicts(samples, "FACIES", grid, geometry) == []
    resolved = resolve_hard_data(samples, "FACIES", grid, geometry)
    assert len(resolved) == 2


def test_duplicate_same_code_in_same_cell_is_not_a_conflict():
    grid, geometry = _geometry()
    samples = [
        WellSample(well="A", x=10.0, y=10.0, values={"FACIES": 1.0}, layer=0),
        WellSample(well="A2", x=12.0, y=11.0, values={"FACIES": 1.0}, layer=0),  # eyni hüceyrə
    ]
    assert detect_hard_data_conflicts(samples, "FACIES", grid, geometry) == []
    resolved = resolve_hard_data(samples, "FACIES", grid, geometry)
    assert len(resolved) == 2   # DUPLİKAT REDDİ EDİLMİR, ikisi də saxlanılır


def test_conflicting_codes_in_same_cell_are_detected():
    grid, geometry = _geometry()
    samples = [
        WellSample(well="A", x=10.0, y=10.0, values={"FACIES": 0.0}, layer=0),
        WellSample(well="B", x=12.0, y=11.0, values={"FACIES": 1.0}, layer=0),  # eyni hüceyrə, FƏRQLİ kod
    ]
    conflicts = detect_hard_data_conflicts(samples, "FACIES", grid, geometry)
    assert len(conflicts) == 1
    assert set(conflicts[0].wells) == {"A", "B"}
    assert set(conflicts[0].codes) == {0, 1}


def test_resolve_raises_by_default_on_conflict():
    grid, geometry = _geometry()
    samples = [
        WellSample(well="A", x=10.0, y=10.0, values={"FACIES": 0.0}, layer=0),
        WellSample(well="B", x=12.0, y=11.0, values={"FACIES": 1.0}, layer=0),
    ]
    with pytest.raises(HardDataConflictError):
        resolve_hard_data(samples, "FACIES", grid, geometry)


def test_resolve_keep_first_and_keep_last_are_deterministic():
    grid, geometry = _geometry()
    samples = [
        WellSample(well="A", x=10.0, y=10.0, values={"FACIES": 0.0}, layer=0),
        WellSample(well="B", x=12.0, y=11.0, values={"FACIES": 1.0}, layer=0),
    ]
    first = resolve_hard_data(samples, "FACIES", grid, geometry, on_conflict="keep_first")
    last = resolve_hard_data(samples, "FACIES", grid, geometry, on_conflict="keep_last")
    assert len(first) == 1 and first[0].well == "A"
    assert len(last) == 1 and last[0].well == "B"


def test_resolve_majority_picks_the_more_common_code():
    grid, geometry = _geometry()
    samples = [
        WellSample(well="A", x=10.0, y=10.0, values={"FACIES": 1.0}, layer=0),
        WellSample(well="B", x=11.0, y=10.0, values={"FACIES": 1.0}, layer=0),
        WellSample(well="C", x=12.0, y=11.0, values={"FACIES": 0.0}, layer=0),
    ]
    resolved = resolve_hard_data(samples, "FACIES", grid, geometry, on_conflict="majority")
    assert len(resolved) == 1
    assert int(resolved[0].values["FACIES"]) == 1


def test_resolve_majority_raises_on_exact_tie():
    grid, geometry = _geometry()
    samples = [
        WellSample(well="A", x=10.0, y=10.0, values={"FACIES": 0.0}, layer=0),
        WellSample(well="B", x=12.0, y=11.0, values={"FACIES": 1.0}, layer=0),
    ]
    with pytest.raises(HardDataConflictError):
        resolve_hard_data(samples, "FACIES", grid, geometry, on_conflict="majority")


def test_unknown_conflict_strategy_rejected():
    grid, geometry = _geometry()
    samples = [WellSample(well="A", x=10.0, y=10.0, values={"FACIES": 0.0}, layer=0)]
    with pytest.raises(ValueError):
        resolve_hard_data(samples, "FACIES", grid, geometry, on_conflict="bogus")


def test_samples_without_layer_or_depth_are_unmapped_and_never_conflict():
    grid, geometry = _geometry()
    samples = [
        WellSample(well="A", x=10.0, y=10.0, values={"FACIES": 0.0}),           # layer/depth yoxdur
        WellSample(well="B", x=12.0, y=11.0, values={"FACIES": 1.0}),           # eyni areal hüceyrə
    ]
    mapping, samples_for = map_samples_to_cells(samples, "FACIES", grid, geometry)
    assert mapping == {}   # heç biri xəritələnmir (K qeyri-müəyyən)
    assert detect_hard_data_conflicts(samples, "FACIES", grid, geometry) == []
    resolved = resolve_hard_data(samples, "FACIES", grid, geometry)
    assert len(resolved) == 2   # ikisi də toxunulmadan saxlanılır


# ── Phase 5: kəsilməz rejim (tolerance-əsaslı) — PORO/PERMX ─────────────
def test_continuous_close_values_are_not_a_conflict_within_tolerance():
    grid, geometry = _geometry()
    samples = [
        WellSample(well="A", x=10.0, y=10.0, layer=0, values={"PORO": 0.201}),
        WellSample(well="B", x=12.0, y=11.0, layer=0, values={"PORO": 0.199}),
    ]
    conflicts = detect_hard_data_conflicts(samples, "PORO", grid, geometry, tolerance=0.01)
    assert conflicts == []
    resolved = resolve_hard_data(samples, "PORO", grid, geometry, tolerance=0.01)
    assert len(resolved) == 2


def test_continuous_large_difference_is_a_conflict():
    grid, geometry = _geometry()
    samples = [
        WellSample(well="A", x=10.0, y=10.0, layer=0, values={"PERMX": 100.0}),
        WellSample(well="B", x=12.0, y=11.0, layer=0, values={"PERMX": 900.0}),
    ]
    conflicts = detect_hard_data_conflicts(samples, "PERMX", grid, geometry, tolerance=50.0)
    assert len(conflicts) == 1
    with pytest.raises(HardDataConflictError):
        resolve_hard_data(samples, "PERMX", grid, geometry, tolerance=50.0)


def test_continuous_average_strategy_synthesizes_mean_value():
    grid, geometry = _geometry()
    samples = [
        WellSample(well="A", x=10.0, y=10.0, layer=0, values={"PORO": 0.20}),
        WellSample(well="B", x=12.0, y=11.0, layer=0, values={"PORO": 0.30}),
    ]
    resolved = resolve_hard_data(samples, "PORO", grid, geometry, on_conflict="average",
                                 tolerance=0.01)
    assert len(resolved) == 1
    assert resolved[0].values["PORO"] == pytest.approx(0.25)


def test_continuous_keep_first_and_keep_last_deterministic():
    grid, geometry = _geometry()
    samples = [
        WellSample(well="A", x=10.0, y=10.0, layer=0, values={"PERMX": 100.0}),
        WellSample(well="B", x=12.0, y=11.0, layer=0, values={"PERMX": 400.0}),
    ]
    first = resolve_hard_data(samples, "PERMX", grid, geometry, on_conflict="keep_first",
                              tolerance=10.0)
    last = resolve_hard_data(samples, "PERMX", grid, geometry, on_conflict="keep_last",
                             tolerance=10.0)
    assert first[0].well == "A" and last[0].well == "B"


def test_majority_rejected_in_continuous_mode_and_average_rejected_in_categorical_mode():
    grid, geometry = _geometry()
    samples = [WellSample(well="A", x=10.0, y=10.0, layer=0, values={"PORO": 0.2})]
    with pytest.raises(ValueError):
        resolve_hard_data(samples, "PORO", grid, geometry, on_conflict="majority", tolerance=0.01)
    with pytest.raises(ValueError):
        resolve_hard_data(samples, "PORO", grid, geometry, on_conflict="average")   # tolerance yoxdur


def test_permy_permz_duplicate_and_conflicting_observations():
    grid, geometry = _geometry()
    for prop in ("PERMY", "PERMZ"):
        duplicate = [
            WellSample(well="A", x=10.0, y=10.0, layer=0, values={prop: 50.0}),
            WellSample(well="B", x=11.0, y=10.0, layer=0, values={prop: 50.2}),
        ]
        assert detect_hard_data_conflicts(duplicate, prop, grid, geometry, tolerance=1.0) == []

        conflicting = [
            WellSample(well="A", x=10.0, y=10.0, layer=0, values={prop: 50.0}),
            WellSample(well="B", x=11.0, y=10.0, layer=0, values={prop: 500.0}),
        ]
        assert len(detect_hard_data_conflicts(conflicting, prop, grid, geometry,
                                              tolerance=1.0)) == 1


def test_depth_based_layer_lookup_maps_correctly():
    grid, geometry = _geometry(nz=2, dz=10.0, top=2000.0)
    samples = [
        WellSample(well="A", x=10.0, y=10.0, values={"FACIES": 0.0}, depth=2005.0),   # K=0
        WellSample(well="B", x=10.0, y=10.0, values={"FACIES": 1.0}, depth=2015.0),   # K=1
    ]
    mapping, _ = map_samples_to_cells(samples, "FACIES", grid, geometry)
    assert len(mapping) == 2   # fərqli laylar -> fərqli hüceyrələr, ziddiyyət YOXDUR
    assert detect_hard_data_conflicts(samples, "FACIES", grid, geometry) == []
