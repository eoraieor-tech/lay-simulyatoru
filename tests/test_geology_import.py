"""Quyu məlumatından geoloji model (B2) testləri."""

import os
import tempfile

import numpy as np

from helpers import default_scal
from imex2d.application.geology_service import (DEFAULT_RULES,
                                                GeologicalGridSpec,
                                                PropertyRule,
                                                WellBasedGeologicalModelBuilder)
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import five_spot
from imex2d.geology.interpolation import (InverseDistance, NearestNeighbour,
                                          OrdinaryKriging, interpolate_property)
from imex2d.geology.well_data_io import (WellDataFormatError, read_well_csv,
                                         write_example_csv)

POINTS = np.array([[0., 0.], [100., 0.], [0., 100.], [100., 100.], [50., 50.]])
VALUES = np.array([0.15, 0.25, 0.20, 0.30, 0.22])


def _write(text: str) -> str:
    handle, path = tempfile.mkstemp(suffix=".csv")
    with os.fdopen(handle, "w", encoding="utf-8") as file:
        file.write(text)
    return path


# ── CSV oxunması ──────────────────────────────────────────────────────
def test_reads_minimal_csv():
    path = _write("well,x,y,PORO,PERMX\nW-1,10,20,0.20,150\nW-2,90,80,0.25,400\n")
    try:
        dataset = read_well_csv(path)
        assert len(dataset) == 2
        assert dataset.well_names == ["W-1", "W-2"]
        assert dataset.property_names() == ["PERMX", "PORO"]
        assert not dataset.is_layered()
    finally:
        os.unlink(path)


def test_reads_layered_csv_and_converts_to_zero_based():
    path = _write("well,x,y,k,PORO\nW-1,10,20,1,0.20\nW-1,10,20,2,0.18\n")
    try:
        dataset = read_well_csv(path)
        assert dataset.is_layered()
        assert sorted(s.layer for s in dataset.samples) == [0, 1]
    finally:
        os.unlink(path)


def test_accepts_alternative_column_names():
    path = _write("Quyu;Easting;Northing;PORO;PERMX\n"
                  "W-1;10;20;0,20;150\nW-2;90;80;0,25;400\n")
    try:
        dataset = read_well_csv(path)
        assert len(dataset) == 2
        assert abs(dataset.samples[0].values["PORO"] - 0.20) < 1e-9
    finally:
        os.unlink(path)


def test_missing_coordinate_columns_is_rejected():
    path = _write("well,PORO\nW-1,0.2\n")
    try:
        read_well_csv(path)
    except WellDataFormatError:
        return
    finally:
        os.unlink(path)
    raise AssertionError("Koordinatsız fayl qəbul edildi")


def test_single_point_dataset_is_rejected():
    path = _write("well,x,y,PORO\nW-1,10,20,0.2\n")
    try:
        read_well_csv(path)
    except WellDataFormatError:
        return
    finally:
        os.unlink(path)
    raise AssertionError("Bir nöqtəli fayl qəbul edildi")


def test_example_file_is_readable():
    handle, path = tempfile.mkstemp(suffix=".csv")
    os.close(handle)
    try:
        write_example_csv(path, nx=21, ny=21, dx=25.0, dy=25.0, nz=2)
        dataset = read_well_csv(path)
        assert len(dataset.well_names) == 5
        assert dataset.is_layered()
        assert "PORO" in dataset.property_names()
    finally:
        os.unlink(path)


# ── interpolyasiya ────────────────────────────────────────────────────
def test_all_interpolators_reproduce_values_at_data_points():
    for interpolator in (NearestNeighbour(), InverseDistance(),
                         OrdinaryKriging()):
        result = interpolator.interpolate(POINTS, VALUES, POINTS)
        assert np.allclose(result, VALUES), interpolator.name


def test_interpolated_values_stay_within_data_range():
    targets = np.array([[25., 25.], [75., 25.], [50., 80.], [10., 90.]])
    for interpolator in (NearestNeighbour(), InverseDistance()):
        result = interpolator.interpolate(POINTS, VALUES, targets)
        assert result.min() >= VALUES.min() - 1e-9
        assert result.max() <= VALUES.max() + 1e-9


def test_higher_idw_power_moves_result_towards_nearest_point():
    target = np.array([[10., 10.]])
    smooth = InverseDistance(power=1.0).interpolate(POINTS, VALUES, target)[0]
    sharp = InverseDistance(power=6.0).interpolate(POINTS, VALUES, target)[0]
    nearest = NearestNeighbour().interpolate(POINTS, VALUES, target)[0]
    assert abs(sharp - nearest) < abs(smooth - nearest)


def test_idw_search_radius_falls_back_to_nearest_when_empty():
    target = np.array([[5000., 5000.]])
    result = InverseDistance(power=2.0, search_radius=10.0).interpolate(
        POINTS, VALUES, target)
    assert abs(result[0] - VALUES[3]) < 1e-9      # ən yaxın nöqtə (100,100)


def test_kriging_with_nugget_is_no_longer_exact():
    exact = OrdinaryKriging(nugget=0.0).interpolate(POINTS, VALUES, POINTS)
    smoothed = OrdinaryKriging(nugget=0.05).interpolate(POINTS, VALUES, POINTS)
    assert np.allclose(exact, VALUES)
    assert not np.allclose(smoothed, VALUES)


def test_single_data_point_returns_constant_field():
    result = OrdinaryKriging().interpolate(POINTS[:1], VALUES[:1],
                                           np.array([[10., 90.], [70., 20.]]))
    assert np.allclose(result, VALUES[0])


def test_log_transform_keeps_permeability_positive():
    permeability = np.array([5.0, 500.0, 50.0, 2000.0, 100.0])
    targets = np.array([[25., 25.], [80., 80.], [10., 90.]])
    result = interpolate_property(OrdinaryKriging(), POINTS, permeability,
                                  targets, log_transform=True, minimum=0.01)
    assert np.all(result > 0.0)
    assert result.max() <= permeability.max() * 1.5


def test_log_transform_rejects_non_positive_values():
    try:
        interpolate_property(InverseDistance(), POINTS,
                             np.array([1.0, 0.0, 1.0, 1.0, 1.0]),
                             POINTS, log_transform=True)
    except ValueError:
        return
    raise AssertionError("Sıfır dəyər log interpolyasiyada qəbul edildi")


def test_limits_are_applied():
    result = interpolate_property(InverseDistance(), POINTS, VALUES,
                                  POINTS, minimum=0.18, maximum=0.26)
    assert result.min() >= 0.18 - 1e-12
    assert result.max() <= 0.26 + 1e-12


# ── geoloji model qurulması ───────────────────────────────────────────
def _dataset(nz=2):
    handle, path = tempfile.mkstemp(suffix=".csv")
    os.close(handle)
    try:
        write_example_csv(path, nx=15, ny=15, dx=30.0, dy=30.0, nz=nz)
        return read_well_csv(path)
    finally:
        os.unlink(path)


def test_builder_creates_all_required_property_maps():
    spec = GeologicalGridSpec(nx=15, ny=15, nz=2, dx=30.0, dy=30.0, dz=5.0)
    model, report = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(
        _dataset(nz=2), spec)
    for key in ("PORO", "PERMX", "PERMY", "PERMZ", "NTG"):
        assert key in model.property_maps, key
        assert model.property_maps[key].values.size == model.grid.ncell
    assert model.validate() == []
    assert report.entries


def test_anisotropy_factors_are_applied_to_missing_permeability():
    spec = GeologicalGridSpec(nx=11, ny=11, nz=1, dx=30.0, dy=30.0)
    model, _ = WellBasedGeologicalModelBuilder(InverseDistance()).build(
        _dataset(nz=1), spec, ky_over_kx=0.5, kv_over_kh=0.05)
    permx = model.property_maps["PERMX"].values
    assert np.allclose(model.property_maps["PERMY"].values, permx * 0.5)
    assert np.allclose(model.property_maps["PERMZ"].values, permx * 0.05)


def test_unlayered_data_is_repeated_across_layers():
    path = _write("well,x,y,PORO,PERMX\n"
                  "W-1,50,50,0.20,150\nW-2,400,50,0.25,400\n"
                  "W-3,50,400,0.22,220\nW-4,400,400,0.18,90\n")
    try:
        dataset = read_well_csv(path)
    finally:
        os.unlink(path)
    spec = GeologicalGridSpec(nx=15, ny=15, nz=4, dx=30.0, dy=30.0)
    model, _ = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(
        dataset, spec)
    porosity = model.property_maps["PORO"].values.reshape(model.grid.shape)
    for k in range(1, 4):
        assert np.allclose(porosity[0], porosity[k])


def test_layered_data_produces_different_layers():
    spec = GeologicalGridSpec(nx=15, ny=15, nz=3, dx=30.0, dy=30.0)
    model, _ = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(
        _dataset(nz=3), spec)
    permeability = model.property_maps["PERMX"].values.reshape(model.grid.shape)
    assert not np.allclose(permeability[0], permeability[-1])
    assert permeability[0].mean() > permeability[-1].mean()


def test_built_model_feeds_the_reservoir_model_builder():
    """B2-nin əsl yoxlanışı: qurulan geoloji model zəncirin qalanına uyğun gəlir."""
    spec = GeologicalGridSpec(nx=13, ny=13, nz=2, dx=30.0, dy=30.0, dz=5.0,
                              top_depth=2000.0)
    geology, _ = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(
        _dataset(nz=2), spec)
    model = ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=default_scal())
    assert model.validate() == []
    assert model.ncell == 13 * 13 * 2
    assert model.rock.permz is not None
    assert model.source_geological_model == geology.name


def test_dipping_surface_is_created_from_spec():
    spec = GeologicalGridSpec(nx=11, ny=11, nz=1, dx=30.0, dy=30.0,
                              top_depth=2000.0, dip_x=4.0)
    model, _ = WellBasedGeologicalModelBuilder(NearestNeighbour()).build(
        _dataset(nz=1), spec)
    depths = model.geometry.cell_depths().reshape(model.grid.shape)[0]
    assert np.all(np.diff(depths, axis=1) > 0)


def test_custom_rule_overrides_default():
    rules = dict(DEFAULT_RULES)
    rules["PORO"] = PropertyRule("PORO", log_transform=False,
                                 minimum=0.19, maximum=0.21)
    spec = GeologicalGridSpec(nx=11, ny=11, nz=1, dx=30.0, dy=30.0)
    model, _ = WellBasedGeologicalModelBuilder(InverseDistance(), rules).build(
        _dataset(nz=1), spec)
    porosity = model.property_maps["PORO"].values
    assert porosity.min() >= 0.19 - 1e-12
    assert porosity.max() <= 0.21 + 1e-12


def test_empty_dataset_is_rejected_by_builder():
    from imex2d.domain.well_data import WellDataset
    try:
        WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(
            WellDataset(), GeologicalGridSpec())
    except ValueError:
        return
    raise AssertionError("Boş məlumat qəbul edildi")
