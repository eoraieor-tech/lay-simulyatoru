"""`validate_wells` — hər yoxlama üçün ən azı bir test."""

from imex2d.domain.geology import GeologicalWell, ValidationIssue, validate_wells
from imex2d.domain.geometry import CellGeometry
from imex2d.domain.grid import CartesianGrid


def _geometry(nx=10, ny=10, dx=20.0, dy=20.0, nz=1):
    grid = CartesianGrid(nx, ny, nz)
    return CellGeometry(grid, dx, dy, dz=10.0, top_depth=2000.0)


def _levels(issues, level=None):
    if level is None:
        return [i.level for i in issues]
    return [i for i in issues if i.level == level]


def test_no_issues_for_a_clean_table():
    wells = [
        GeologicalWell(name="A", x=50.0, y=50.0, top=2000.0, bottom=2010.0,
                       porosity=0.2, permeability=100.0, water_saturation=0.3),
        GeologicalWell(name="B", x=150.0, y=150.0, top=2005.0, bottom=2015.0,
                       porosity=0.22, permeability=120.0, water_saturation=0.28),
        GeologicalWell(name="C", x=90.0, y=30.0, top=2001.0, bottom=2011.0,
                       porosity=0.19, permeability=95.0, water_saturation=0.31),
    ]
    issues = validate_wells(wells, _geometry(), method="Kriging (adi)")
    assert _levels(issues, "error") == []


def test_empty_name_is_error():
    issues = validate_wells([GeologicalWell(name="", x=0.0, y=0.0)])
    assert any(i.level == "error" and "boşdur" in i.message for i in issues)


def test_duplicate_name_is_error():
    wells = [GeologicalWell(name="A", x=0.0, y=0.0),
            GeologicalWell(name="A", x=10.0, y=10.0)]
    issues = validate_wells(wells)
    assert any(i.level == "error" and "təkrarlanır" in i.message for i in issues)


def test_out_of_bounds_xy_is_error():
    wells = [GeologicalWell(name="A", x=-5.0, y=50.0)]
    issues = validate_wells(wells, _geometry(nx=10, ny=10, dx=20.0, dy=20.0))
    assert any(i.level == "error" and "kənardadır" in i.message for i in issues)


def test_wells_in_same_cell_is_warning():
    wells = [GeologicalWell(name="A", x=15.0, y=15.0),
            GeologicalWell(name="B", x=18.0, y=12.0)]
    issues = validate_wells(wells, _geometry(nx=10, ny=10, dx=20.0, dy=20.0))
    assert any(i.level == "warning" and "eyni hüceyrədədir" in i.message for i in issues)


def test_porosity_out_of_range_is_error():
    wells = [GeologicalWell(name="A", x=10.0, y=10.0, porosity=1.4)]
    issues = validate_wells(wells)
    assert any(i.level == "error" and "φ" in i.message for i in issues)


def test_water_saturation_out_of_range_is_error():
    wells = [GeologicalWell(name="A", x=10.0, y=10.0, water_saturation=-0.1)]
    issues = validate_wells(wells)
    assert any(i.level == "error" and "Sw" in i.message for i in issues)


def test_non_positive_permeability_is_error():
    wells = [GeologicalWell(name="A", x=10.0, y=10.0, permeability=0.0)]
    issues = validate_wells(wells)
    assert any(i.level == "error" and "k = " in i.message for i in issues)


def test_bottom_shallower_than_top_is_error():
    wells = [GeologicalWell(name="A", x=10.0, y=10.0, top=2020.0, bottom=2010.0)]
    issues = validate_wells(wells)
    assert any(i.level == "error" and "lay altı" in i.message and "lay üstü" in i.message
              for i in issues)


def test_insufficient_wells_for_kriging_is_warning_not_blocking():
    """Xassə-üzrə çatışmazlıq BLOKLAMIR (yalnız həmin xassə buraxılır) —
    permeability 3 quyuda var (kifayətdir), porosity yalnız 2-də (azdır)."""
    wells = [
        GeologicalWell(name="A", x=10.0, y=10.0, porosity=0.2, permeability=100.0),
        GeologicalWell(name="B", x=50.0, y=50.0, porosity=0.25, permeability=110.0),
        GeologicalWell(name="C", x=90.0, y=30.0, permeability=90.0),
    ]
    issues = validate_wells(wells, method="Kriging (adi)")
    assert any(i.level == "warning" and "φ" in i.message for i in issues)
    assert _levels(issues, "error") == []


def test_no_property_has_enough_wells_is_error():
    wells = [GeologicalWell(name="A", x=10.0, y=10.0, porosity=0.2)]
    issues = validate_wells(wells, method="Kriging (adi)")
    assert any(i.level == "error" and "Heç bir xassə" in i.message for i in issues)


def test_in_model_without_reservoir_regime_is_warning():
    wells = [GeologicalWell(name="A", x=10.0, y=10.0, in_model=True)]
    issues = validate_wells(wells, reservoir_well_names=[])
    assert any(i.level == "warning" and "rejim" in i.message for i in issues)


def test_orphan_reservoir_well_is_warning():
    wells = [GeologicalWell(name="A", x=10.0, y=10.0, in_model=True)]
    issues = validate_wells(wells, reservoir_well_names=["A", "GHOST"])
    assert any(i.level == "warning" and "sahibsiz" in i.message for i in issues)


def test_well_on_cell_edge_is_warning():
    wells = [GeologicalWell(name="A", x=20.0, y=10.0)]     # tam hüceyrə sərhəddində
    issues = validate_wells(wells, _geometry(nx=10, ny=10, dx=20.0, dy=20.0))
    assert any(i.level == "warning" and "kənar" in i.message for i in issues)


def test_missing_fields_are_info():
    wells = [GeologicalWell(name="A", x=10.0, y=10.0)]      # heç bir xassə yoxdur
    issues = validate_wells(wells)
    assert any(i.level == "info" and "boş xanalar" in i.message for i in issues)


def test_validation_issue_carries_level_and_well():
    issue = ValidationIssue("error", "test mesajı", "A")
    assert issue.level == "error"
    assert issue.well == "A"


def test_no_grid_related_checks_when_geometry_is_none():
    """Grid hələ qurulmayıbsa sərhəd/hüceyrə yoxlamaları keçilir, xəta atılmır."""
    wells = [GeologicalWell(name="A", x=999999.0, y=999999.0)]
    issues = validate_wells(wells, geometry=None)
    assert not any("kənardadır" in i.message for i in issues)
