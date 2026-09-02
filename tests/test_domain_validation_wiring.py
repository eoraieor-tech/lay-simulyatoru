"""Phase 1 — mərkəzləşdirilmiş yoxlamanın domen obyektlərinə bağlanması.

Yeni `domain/validation.py` funksiyalarının özü `test_validation.py`-da
sınanıb — bura YALNIZ onların mövcud domen siniflərinə (CellGeometry,
RockProperties, PVTTable, SaturationTable, WellControl, GeologicalModel)
DOĞRU bağlandığını yoxlayır."""

from __future__ import annotations

import numpy as np

from imex2d.domain.geological_model import GeologicalModel
from imex2d.domain.geometry import CellGeometry
from imex2d.domain.grid import CartesianGrid
from imex2d.domain.properties import FluidProperties, PropertyMap, RockProperties
from imex2d.domain.pvt import PVTTable
from imex2d.domain.scal_tables import SaturationTable
from imex2d.domain.structure import RegionSet
from imex2d.domain.wells import ControlMode, WellControl


# ── PropertyMap vahid reyestri (Phase 1, giriş boru xətti) ──────────────
def test_property_map_accepts_valid_permeability_units():
    for unit in ("mD", "D", "m2"):
        PropertyMap.from_array("PERMX", [100.0, 200.0], 2, unit)


def test_property_map_rejects_permeability_with_pressure_unit():
    try:
        PropertyMap.from_array("PERMX", [100.0, 200.0], 2, "psi")
    except ValueError as exc:
        assert "PERMX" in str(exc) and "psi" in str(exc)
        return
    raise AssertionError("PERMX + 'psi' qəbul edildi")


def test_property_map_rejects_pressure_with_permeability_unit():
    try:
        PropertyMap.from_array("PRESSURE", [200.0, 210.0], 2, "mD")
    except ValueError as exc:
        assert "PRESSURE" in str(exc) and "mD" in str(exc)
        return
    raise AssertionError("PRESSURE + 'mD' qəbul edildi")


def test_property_map_missing_unit_never_rejected():
    """Boş `unit` (defolt) HEÇ VAXT rədd edilmir — mövcud kod bazasında
    PERMX/PORO HƏMİŞƏ ya boş, ya da 'mD' ilə qurulur (bax audit)."""
    PropertyMap.from_array("PERMX", [100.0, 200.0], 2)
    PropertyMap.from_array("PORO", [0.2, 0.3], 2)


def test_property_map_unregistered_name_accepts_any_unit_label():
    """Reyestrdə olmayan ad (məs. PORO) üçün `unit` sərbəst mətn olaraq
    qalır — bu davranış DƏYİŞMƏYİB."""
    PropertyMap.from_array("PORO", [0.2, 0.3], 2, "anything")


# ── CellGeometry.validate() — əvvəllər HEÇ YERDƏ yox idi (audit) ───────
def test_cell_geometry_validate_accepts_normal_grid():
    grid = CartesianGrid(4, 4, 2)
    geometry = CellGeometry(grid, dx=20.0, dy=20.0, dz=10.0)
    assert geometry.validate() == []


def test_cell_geometry_validate_rejects_zero_dx():
    grid = CartesianGrid(4, 4, 1)
    geometry = CellGeometry(grid, dx=0.0, dy=20.0, dz=10.0)
    assert geometry.validate()


def test_cell_geometry_validate_rejects_negative_thickness():
    grid = CartesianGrid(4, 4, 2)
    geometry = CellGeometry(grid, dx=20.0, dy=20.0, dz=[10.0, -5.0])
    assert geometry.validate()


# ── RockProperties ───────────────────────────────────────────────────────
def test_rock_properties_validate_rejects_porosity_at_or_above_one():
    n = 3
    rock = RockProperties(
        porosity=PropertyMap.from_array("PORO", np.array([0.2, 1.0, 0.3]), n),
        permx=PropertyMap.uniform("PERMX", 100.0, n),
        permy=PropertyMap.uniform("PERMY", 100.0, n))
    assert rock.validate()


def test_rock_properties_warnings_flag_extreme_permeability_without_rejecting():
    n = 2
    rock = RockProperties(
        porosity=PropertyMap.uniform("PORO", 0.2, n),
        permx=PropertyMap.from_array("PERMX", np.array([100.0, 50000.0]), n),
        permy=PropertyMap.uniform("PERMY", 100.0, n))
    assert rock.validate() == []
    assert rock.validate_warnings()


def test_rock_properties_validate_rejects_nan_porosity():
    """Reqressiya: `validate()` əvvəllər `values <= 0` kimi XAM müqayisə
    işlədirdi — `NaN <= 0` HƏMİŞƏ `False` olduğu üçün NaN PORO SƏSSİZCƏ
    keçib gedirdi (bax `properties.py::RockProperties.validate`)."""
    n = 3
    rock = RockProperties(
        porosity=PropertyMap.from_array("PORO", np.array([0.2, np.nan, 0.3]), n),
        permx=PropertyMap.uniform("PERMX", 100.0, n),
        permy=PropertyMap.uniform("PERMY", 100.0, n))
    assert rock.validate()


def test_rock_properties_validate_rejects_nan_and_inf_permeability():
    n = 3
    rock_nan = RockProperties(
        porosity=PropertyMap.uniform("PORO", 0.2, n),
        permx=PropertyMap.from_array("PERMX", np.array([100.0, np.nan, 100.0]), n),
        permy=PropertyMap.uniform("PERMY", 100.0, n))
    assert rock_nan.validate()

    rock_inf = RockProperties(
        porosity=PropertyMap.uniform("PORO", 0.2, n),
        permx=PropertyMap.uniform("PERMX", 100.0, n),
        permy=PropertyMap.from_array("PERMY", np.array([100.0, np.inf, 100.0]), n))
    assert rock_inf.validate()


# ── FluidProperties (əvvəllər HEÇ BİR yoxlama yox idi) ──────────────────
def test_fluid_properties_default_is_valid():
    assert FluidProperties().validate() == []


def test_fluid_properties_rejects_non_positive_viscosity():
    assert FluidProperties(oil_viscosity=0.0).validate()
    assert FluidProperties(water_viscosity=-1.0).validate()


def test_fluid_properties_rejects_non_positive_density():
    assert FluidProperties(oil_density=-10.0).validate()


def test_fluid_properties_warns_on_heavy_oil_without_rejecting():
    props = FluidProperties(oil_viscosity=90000.0)
    assert props.validate() == []
    assert props.validate_warnings()


# ── PVTTable ─────────────────────────────────────────────────────────────
def _sample_pvt(**overrides):
    n = 5
    defaults = dict(
        pressure=np.linspace(1.0, 300.0, n),
        oil_fvf=np.linspace(1.2, 1.1, n),
        oil_viscosity=np.linspace(3.0, 2.0, n),
        solution_gor=np.linspace(10.0, 80.0, n),
        water_fvf=np.full(n, 1.02),
        water_viscosity=np.full(n, 0.5))
    defaults.update(overrides)
    return PVTTable(**defaults)


def test_pvt_validate_detects_nan_explicitly():
    table = _sample_pvt(oil_fvf=np.array([1.2, np.nan, 1.15, 1.1, 1.05]))
    issues = table.validate()
    assert any("NaN" in issue for issue in issues)


def test_pvt_validate_warnings_flag_unusual_viscosity_without_rejecting():
    table = _sample_pvt(oil_viscosity=np.full(5, 80000.0))
    assert table.validate() == []
    assert table.validate_warnings()


def test_pvt_check_query_range_flags_out_of_bounds():
    table = _sample_pvt()
    warnings = table.check_query_range([0.5, 150.0, 500.0])
    assert len(warnings) == 2


# ── SaturationTable ──────────────────────────────────────────────────────
def test_saturation_table_validate_detects_nan_explicitly():
    table = SaturationTable(sw=np.array([0.2, 0.4, np.nan, 0.8]),
                            krw=np.array([0.0, 0.1, 0.2, 0.4]),
                            kro=np.array([0.8, 0.5, 0.2, 0.0]))
    issues = table.validate()
    assert any("NaN" in issue for issue in issues)


def test_saturation_table_check_query_range_flags_extrapolation():
    table = SaturationTable(sw=np.array([0.2, 0.5, 0.8]),
                            krw=np.array([0.0, 0.1, 0.4]),
                            kro=np.array([0.6, 0.2, 0.0]))
    warnings = table.check_query_range([0.1, 0.5, 0.95])
    assert len(warnings) == 2


# ── WellControl ──────────────────────────────────────────────────────────
def test_well_control_rejects_negative_rate():
    control = WellControl(ControlMode.RATE, -10.0)
    assert control.validate()


def test_well_control_accepts_normal_bhp_and_rate():
    assert WellControl(ControlMode.BHP, 200.0).validate() == []
    assert WellControl(ControlMode.RATE, 60.0).validate() == []


def test_well_control_rejects_non_positive_bhp():
    assert WellControl(ControlMode.BHP, 0.0).validate()
    assert WellControl(ControlMode.BHP, -50.0).validate()


def test_well_control_warns_on_zero_rate_without_rejecting():
    control = WellControl(ControlMode.RATE, 0.0)
    assert control.validate() == []
    assert control.validate_warnings()


# ── GeologicalModel: geometriya + PORO/PERM sərt yoxlaması bağlanıb ─────
def _model_with(poro_values, permx_values, dx=20.0, dy=20.0, dz=10.0):
    grid = CartesianGrid(2, 2, 1)
    geometry = CellGeometry(grid, dx=dx, dy=dy, dz=dz)
    model = GeologicalModel(name="test", grid=grid, geometry=geometry,
                            regions=RegionSet.single(grid.ncell))
    model.add_property(PropertyMap.from_array("PORO", poro_values, grid.ncell))
    model.add_property(PropertyMap.from_array("PERMX", permx_values, grid.ncell))
    return model


def test_geological_model_validate_accepts_normal_case():
    model = _model_with([0.2, 0.2, 0.2, 0.2], [100.0, 100.0, 100.0, 100.0])
    assert model.validate() == []


def test_geological_model_validate_rejects_impossible_porosity():
    model = _model_with([0.2, 1.5, 0.2, 0.2], [100.0, 100.0, 100.0, 100.0])
    assert model.validate()


def test_geological_model_validate_rejects_degenerate_geometry():
    model = _model_with([0.2] * 4, [100.0] * 4, dx=0.0)
    assert model.validate()


def test_geological_model_validate_warnings_flags_extreme_permeability():
    model = _model_with([0.2] * 4, [100.0, 100.0, 100.0, 90000.0])
    assert model.validate() == []
    assert model.validate_warnings()
