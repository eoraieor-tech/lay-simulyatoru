"""Domain qatının vahid testləri — hesablama mühərriki cəlb edilmir."""

import numpy as np

from helpers import default_scal, five_spot_model
from imex2d.domain.geometry import CellGeometry
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


# ── dəyişən təbəqə qalınlığı (per-layer DZ) ─────────────────────────────
def _variable_dz_geometry():
    """1x1x3 grid, təbəqə qalınlıqları 2/4/6 m — əl ilə yoxlanıla bilən ölçü."""
    grid = CartesianGrid(1, 1, 3)
    return grid, CellGeometry(grid, dx=10.0, dy=20.0, dz=[2.0, 4.0, 6.0],
                              top_depth=100.0)


def test_scalar_dz_is_normalised_to_a_per_layer_array():
    grid = CartesianGrid(2, 2, 3)
    geometry = CellGeometry(grid, dx=10.0, dy=10.0, dz=5.0)
    assert isinstance(geometry.dz, np.ndarray)
    assert geometry.dz.shape == (3,)
    assert np.all(geometry.dz == 5.0)


def test_wrong_length_dz_array_is_rejected():
    grid = CartesianGrid(2, 2, 3)
    try:
        CellGeometry(grid, dx=10.0, dy=10.0, dz=[1.0, 2.0])
    except ValueError:
        return
    raise AssertionError("NZ-ə uyğun olmayan dz massivi qəbul edildi")


def test_volumes_use_each_layers_own_thickness():
    grid, geometry = _variable_dz_geometry()
    expected = np.array([10.0 * 20.0 * 2.0, 10.0 * 20.0 * 4.0, 10.0 * 20.0 * 6.0])
    assert np.allclose(geometry.volumes(), expected)


def test_cell_depths_accumulate_variable_layer_thickness():
    grid, geometry = _variable_dz_geometry()
    # təbəqə tavanları: 0, 2, 6 (kumulyativ dz) -> mərkəzlər: 1, 4, 9
    expected = 100.0 + np.array([1.0, 4.0, 9.0])
    assert np.allclose(geometry.cell_depths(), expected)


def test_face_half_distances_differ_per_side_when_layers_differ():
    grid, geometry = _variable_dz_geometry()
    conn = grid.build_connections()
    half_a, half_b = geometry.face_half_distances(conn)
    # yeganə bağlantılar K istiqamətindədir (nx=ny=1): 0-1 və 1-2
    assert np.allclose(half_a, [1.0, 2.0])   # dz[0]/2, dz[1]/2
    assert np.allclose(half_b, [2.0, 3.0])   # dz[1]/2, dz[2]/2

    area = geometry.face_areas(conn)
    assert np.allclose(area, 10.0 * 20.0)    # K üzləri sadəcə dx*dy-dir
