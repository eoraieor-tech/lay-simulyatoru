"""Domain qatının vahid testləri — hesablama mühərriki cəlb edilmir."""

import numpy as np

from helpers import default_scal, five_spot_model
from imex2d.domain.grid import CartesianGrid
from imex2d.domain.properties import PropertyMap, RockProperties
from imex2d.domain.scal import CoreyParameters
from imex2d.domain.structure import RegionSet
from imex2d.domain.wells import (ControlMode, Perforation, Well, WellControl,
                                 WellType)


def test_grid_cell_count_and_indexing():
    grid = CartesianGrid(5, 4, 3)
    assert grid.ncell == 60
    assert grid.index(0, 0, 0) == 0
    assert grid.index(4, 3, 2) == 59
    for cell in (0, 17, 42, 59):
        i, j, k = grid.ijk(cell)
        assert grid.index(i, j, k) == cell


def test_grid_connection_count_matches_structured_formula():
    grid = CartesianGrid(5, 4, 3)
    conn = grid.build_connections()
    expected = ((5 - 1) * 4 * 3) + (5 * (4 - 1) * 3) + (5 * 4 * (3 - 1))
    assert conn.count == expected
    assert conn.cell_a.size == conn.cell_b.size == conn.axis.size


def test_single_cell_grid_has_no_connections():
    conn = CartesianGrid(1, 1, 1).build_connections()
    assert conn.count == 0


def test_property_map_rejects_wrong_size():
    try:
        PropertyMap.from_array("PORO", np.zeros(10), ncell=20)
    except ValueError:
        return
    raise AssertionError("Yanlış ölçülü massiv qəbul edildi")


def test_rock_properties_validation_detects_zero_permeability():
    n = 4
    rock = RockProperties(
        porosity=PropertyMap.uniform("PORO", 0.2, n),
        permx=PropertyMap.from_array("PERMX", np.array([100., 0., 100., 100.]), n),
        permy=PropertyMap.uniform("PERMY", 100.0, n))
    assert rock.validate(), "Sıfır keçiricilik aşkarlanmadı"


def test_corey_endpoints_and_monotonicity():
    scal = CoreyParameters()
    assert abs(scal.krw(scal.swc)) < 1e-12
    assert abs(scal.kro(1.0 - scal.sor)) < 1e-12
    assert abs(scal.krw(1.0 - scal.sor) - scal.krw_end) < 1e-12
    assert abs(scal.kro(scal.swc) - scal.kro_end) < 1e-12
    sw = np.linspace(scal.swc, 1.0 - scal.sor, 50)
    assert np.all(np.diff(scal.krw(sw)) >= -1e-12), "krw monoton artmır"
    assert np.all(np.diff(scal.kro(sw)) <= 1e-12), "kro monoton azalmır"


def test_corey_validation_rejects_impossible_saturations():
    assert CoreyParameters(swc=0.6, sor=0.5).validate()


def test_region_set_defaults_to_single_region():
    regions = RegionSet.single(9)
    assert list(regions.ids) == [1]


def test_reservoir_model_requires_pressure_constrained_well():
    model = five_spot_model(nx=5, ny=5)
    model.wells = [Well("P1", WellType.PRODUCER,
                        WellControl(ControlMode.RATE, 50.0), [Perforation(0, 0, 0)])]
    issues = model.validate()
    assert any("BHP" in issue for issue in issues)


def test_reservoir_model_detects_perforation_outside_grid():
    model = five_spot_model(nx=5, ny=5)
    model.wells[0].perforations = [Perforation(99, 0, 0)]
    assert any("kənar" in issue for issue in model.validate())


def test_valid_model_reports_no_issues():
    assert five_spot_model(nx=9, ny=9).validate() == []


def test_pore_volume_matches_manual_calculation():
    model = five_spot_model(nx=4, ny=3, dx=10.0, dy=20.0, dz=5.0, porosity=0.25)
    expected = 4 * 3 * 10.0 * 20.0 * 5.0 * 0.25
    assert abs(model.pore_volume().sum() - expected) < 1e-6
