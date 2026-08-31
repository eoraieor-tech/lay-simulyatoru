"""`wells_to_dataset` — GeologicalWell cədvəli → WellDataset adapteri."""

from imex2d.application.geology_adapter import wells_to_dataset
from imex2d.domain.geology import GeologicalWell


def test_wells_with_full_data_produce_all_properties():
    wells = [
        GeologicalWell(name="A", x=0.0, y=0.0, porosity=0.2,
                       permeability=100.0, water_saturation=0.3),
        GeologicalWell(name="B", x=100.0, y=0.0, porosity=0.22,
                       permeability=120.0, water_saturation=0.28),
        GeologicalWell(name="C", x=0.0, y=100.0, porosity=0.19,
                       permeability=90.0, water_saturation=0.31),
    ]
    dataset, skipped = wells_to_dataset(wells, method="Kriging (adi)")
    assert set(dataset.property_names()) == {"PORO", "PERMX", "SW"}
    assert skipped == {}
    assert len(dataset) == 3


def test_missing_value_drops_only_that_property_not_the_row():
    """B2-nin ən böyük qüsuru düzəldilib: Sw yoxdursa sətir düşmür."""
    wells = [
        GeologicalWell(name="A", x=0.0, y=0.0, porosity=0.2,
                       permeability=100.0, water_saturation=0.3),
        GeologicalWell(name="B", x=100.0, y=0.0, porosity=0.22,
                       permeability=120.0),                    # Sw yoxdur
        GeologicalWell(name="C", x=0.0, y=100.0, porosity=0.19,
                       permeability=90.0, water_saturation=0.31),
    ]
    dataset, _ = wells_to_dataset(wells, method="Əks məsafə (IDW)")
    _, poro_values = dataset.points("PORO")
    _, sw_values = dataset.points("SW")
    assert poro_values.size == 3
    assert sw_values.size == 2
    assert "B" in [s.well for s in dataset.samples]


def test_property_with_too_few_wells_is_skipped_entirely_for_kriging():
    wells = [
        GeologicalWell(name="A", x=0.0, y=0.0, porosity=0.2, water_saturation=0.3),
        GeologicalWell(name="B", x=100.0, y=0.0, porosity=0.22, water_saturation=0.28),
        GeologicalWell(name="C", x=0.0, y=100.0, porosity=0.19),  # Sw yoxdur -> Sw=2 quyu
    ]
    dataset, skipped = wells_to_dataset(wells, method="Kriging (adi)")
    assert "SW" not in dataset.property_names()      # 2 < 3 (Kriging)
    assert "PORO" in dataset.property_names()          # 3 >= 3
    assert "SW" in skipped


def test_idw_accepts_two_wells_where_kriging_would_skip():
    wells = [
        GeologicalWell(name="A", x=0.0, y=0.0, porosity=0.2),
        GeologicalWell(name="B", x=100.0, y=0.0, porosity=0.22),
    ]
    dataset, skipped = wells_to_dataset(wells, method="Əks məsafə (IDW)")
    assert "PORO" in dataset.property_names()
    assert skipped == {}


def test_top_and_bottom_are_treated_as_properties():
    wells = [
        GeologicalWell(name="A", x=0.0, y=0.0, top=2000.0, bottom=2010.0),
        GeologicalWell(name="B", x=100.0, y=0.0, top=2005.0, bottom=2015.0),
    ]
    dataset, _ = wells_to_dataset(wells, method="Əks məsafə (IDW)")
    assert {"TOP", "BOTTOM"} <= set(dataset.property_names())


def test_empty_well_list_produces_empty_dataset():
    dataset, skipped = wells_to_dataset([], method="Kriging (adi)")
    assert len(dataset) == 0
    assert skipped == {}
