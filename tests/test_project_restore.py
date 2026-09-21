"""Seans 46 — layihə faylı hesablamanı TAM bərpa edir.

Seans 46-ya qədər `.imx` bunları səssizcə itirirdi (ölçüldü):

  * qaz SCAL (`gas_scal_parameters`) — açılan qazlı model defolt qaz
    SCAL ilə hesablanırdı: RF 64.75 % → 64.92 %, heç bir xəbərdarlıq yox;
  * `FluidProperties.gas_density` — əl siyahısında yox idi;
  * SWOF/SGOF/CSV cədvəlləri (`scal_tables`, `gas_scal_tables`);
  * nəticədə qaz sıraları, Sg anları, BHP/THP, idarə rejimi, OGIP.

ƏSAS MÜQAVİLƏ: saxlanılıb açılan model yenidən işlədiləndə nəticə
BİT-BİT eyni olmalıdır.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.project import Project
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.application.serialization import ProjectSerializer
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.properties import FluidProperties
from imex2d.domain.scal import GasCoreyParameters
from imex2d.domain.scal_tables import (GasSaturationTable, GasSaturationTableSet,
                                      SaturationTable, SaturationTableSet)
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.results import SimulationResult, Snapshot, TimeSeries
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter


def _gas_model(**extra):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=8, ny=8, dx=20.0, dy=20.0, dz=10.0, porosity=0.22,
        permx_base=150.0, top_depth=1200.0)
    return geology, ReservoirModelBuilder().build(
        geological_model=geology, wells=five_spot(geology.grid),
        scal=default_scal(), gas_scal=GasCoreyParameters(nog=3.3, sgc=0.07),
        pvt_table=build_pvt_table(pressure_min=1.0, pressure_max=400.0,
                                  n_points=40, bubble_point_bar=240.0,
                                  include_gas=True),
        name="bərpa sınağı", **extra)


def _service():
    return ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)


def _save_and_load(tmp_path, project) -> Project:
    path = str(tmp_path / "layihe.imx")
    ProjectSerializer().save(project, path)
    return ProjectSerializer().load(path)


def _project(geology, model) -> Project:
    project = Project("bərpa")
    project.add_geological_model(geology)
    project.add_reservoir_model(model)
    return project


# ═══════════════════════════ model ═══════════════════════════════════

def test_reopened_gas_model_reproduces_result_bit_for_bit(tmp_path):
    """Seans 46-dan əvvəl: 64.7536 % → 64.9197 % (qaz SCAL itirdi)."""
    geology, model = _gas_model()
    config = SimulationConfig(end_time=300.0)
    original = _service().run(model, config)

    back = _save_and_load(tmp_path, _project(geology, model))
    reopened = list(back.reservoir_models.values())[-1]
    again = _service().run(reopened, config)

    assert reopened.gas_scal_parameters == model.gas_scal_parameters
    assert again.series.recovery_factor == original.series.recovery_factor
    assert again.series.gas_oil_ratio == original.series.gas_oil_ratio


def test_all_fluid_fields_survive(tmp_path):
    geology, model = _gas_model()
    model.fluids = dataclasses.replace(FluidProperties(), gas_density=0.93,
                                       water_density=1030.0)
    back = _save_and_load(tmp_path, _project(geology, model))
    assert list(back.reservoir_models.values())[-1].fluids == model.fluids


def test_scal_tables_survive(tmp_path):
    geology, model = _gas_model()
    model.scal_tables = SaturationTableSet(tables={
        1: SaturationTable(sw=[0.2, 0.5, 0.8], krw=[0.0, 0.1, 0.4],
                           kro=[0.9, 0.3, 0.0], pc=[2.0, 0.5, 0.0], name="kern-1"),
        3: SaturationTable(sw=[0.25, 0.75], krw=[0.0, 0.3], kro=[0.8, 0.0])},
        default_region=3)
    model.gas_scal_tables = GasSaturationTableSet(tables={
        1: GasSaturationTable(sg=[0.0, 0.4, 0.7], krg=[0.0, 0.2, 0.8],
                              krog=[0.9, 0.2, 0.0])})
    back = list(_save_and_load(tmp_path, _project(geology, model))
                .reservoir_models.values())[-1]

    assert back.scal_tables.default_region == 3
    assert set(back.scal_tables.tables) == {1, 3}
    first = back.scal_tables.tables[1]
    np.testing.assert_array_equal(first.kro, [0.9, 0.3, 0.0])
    np.testing.assert_array_equal(first.pc, [2.0, 0.5, 0.0])
    assert first.name == "kern-1"
    assert back.scal_tables.tables[3].pc is None
    np.testing.assert_array_equal(back.gas_scal_tables.tables[1].krg,
                                  [0.0, 0.2, 0.8])


def test_two_phase_model_keeps_no_gas_scal(tmp_path):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=4, ny=4, dx=20.0, dy=20.0, dz=10.0, porosity=0.2, permx_base=100.0)
    model = ReservoirModelBuilder().build(
        geological_model=geology, wells=five_spot(geology.grid),
        scal=default_scal(), name="iki fazalı")
    back = list(_save_and_load(tmp_path, _project(geology, model))
                .reservoir_models.values())[-1]
    assert back.gas_scal_parameters is None
    assert back.scal_tables is None and back.gas_scal_tables is None


# ═══════════════════════════ nəticə ══════════════════════════════════

def _rich_result() -> SimulationResult:
    series = TimeSeries(time=[1.0, 2.0], oil_rate=[10.0, 9.0],
                        gas_rate=[500.0, 520.0], gas_injection_rate=[0.0, 0.0],
                        cumulative_gas=[500.0, 1020.0],
                        gas_oil_ratio=[50.0, 57.8])
    shape = (2, 1, 1)
    return SimulationResult(
        model_name="m", grid_shape=shape, series=series, ogip=1.5e6,
        snapshots=[Snapshot(time=1.0, pressure=np.full(shape, 200.0),
                            water_saturation=np.full(shape, 0.2),
                            gas_saturation=np.array([0.0, 0.1]).reshape(shape))],
        well_gas_rate={"P": [500.0, 520.0]},
        well_bhp={"P": [100.0, 101.0]},
        well_thp={"P": [12.0, float("nan")]},
        well_control_mode={"P": ["RATE", "BHP"]})


def test_result_gas_and_well_pressure_series_survive(tmp_path):
    geology, model = _gas_model()
    project = _project(geology, model)
    run = project.new_run(model.name, SimulationConfig(end_time=2.0))
    run.result = _rich_result()

    back = _save_and_load(tmp_path, project).latest_run().result
    assert back.series.gas_rate == [500.0, 520.0]
    assert back.series.gas_oil_ratio == [50.0, 57.8]
    assert back.series.cumulative_gas == [500.0, 1020.0]
    assert back.ogip == 1.5e6
    assert back.well_gas_rate == {"P": [500.0, 520.0]}
    assert back.well_bhp == {"P": [100.0, 101.0]}
    assert back.well_thp["P"][0] == 12.0 and math.isnan(back.well_thp["P"][1])
    assert back.well_control_mode == {"P": ["RATE", "BHP"]}
    np.testing.assert_array_equal(back.snapshots[0].gas_saturation.ravel(),
                                  [0.0, 0.1])


def test_old_result_without_new_keys_still_loads():
    """Seans 46-dan əvvəlki fayl: yeni açarlar yoxdur, xəta da yoxdur."""
    data = ProjectSerializer._result_to_dict(_rich_result())
    for key in ("ogip", "well_gas_rate", "well_bhp", "well_thp",
                "well_control_mode"):
        del data[key]
    del data["snapshots"][0]["gas_saturation"]
    data["series"] = {name: data["series"][name] for name in (
        "time", "oil_rate", "water_rate", "water_injection_rate",
        "cumulative_oil", "cumulative_water", "water_cut",
        "average_pressure", "recovery_factor")}
    data["series"]["obsolete_series"] = [1.0]      # naməlum açar atılır

    back = ProjectSerializer._result_from_dict(data)
    assert back.ogip == 0.0 and back.well_bhp == {}
    assert back.series.gas_rate == []
    assert back.snapshots[0].gas_saturation is None


def test_old_model_without_new_keys_still_loads():
    geology, model = _gas_model()
    data = ProjectSerializer().reservoir_model_to_dict(model)
    for key in ("gas_scal", "scal_tables", "gas_scal_tables"):
        del data[key]
    del data["fluids"]["gas_density"]
    back = ProjectSerializer().reservoir_model_from_dict(data)
    assert back.gas_scal_parameters is None
    assert back.fluids.gas_density == FluidProperties().gas_density
