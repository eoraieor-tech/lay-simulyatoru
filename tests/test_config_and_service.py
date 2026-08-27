"""Konfiqurasiya yoxlaması və servis qatının davranışı."""

from helpers import default_scal, five_spot_model, make_service, short_config
from imex2d.application.config import SimulationConfig, TimeSteppingConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.project import Project
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder, five_spot
from imex2d.application.simulation_service import (ModelValidationError,
                                                   SimulationService)
from imex2d.domain.wells import ControlMode, Perforation, Well, WellControl, WellType
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter


def test_config_rejects_negative_end_time():
    assert SimulationConfig(end_time=-1.0).validate()


def test_config_rejects_invalid_cfl():
    config = SimulationConfig(time_stepping=TimeSteppingConfig(cfl_factor=1.5))
    assert config.validate()


def test_default_config_is_valid():
    assert SimulationConfig().validate() == []


def test_service_raises_on_invalid_model():
    model = five_spot_model(nx=5, ny=5)
    model.wells = [Well("P1", WellType.PRODUCER,
                        WellControl(ControlMode.RATE, 10.0), [Perforation(0, 0, 0)])]
    try:
        make_service().run(model, short_config(end_time=5.0))
    except ModelValidationError as error:
        assert error.issues
        return
    raise AssertionError("Yanlış model qəbul edildi")


def test_pvt_provider_is_accepted_by_engine():
    """A1-dən sonra: PVT provider inject edilə bilər və mühərrik onu qəbul edir."""
    from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
    from imex2d.simulation.pvt.correlations import build_pvt_table

    scal = default_scal()
    service = SimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(scal),
        pvt_provider=BlackOilPVTProvider(build_pvt_table()))
    engine = service.create_engine(five_spot_model(nx=5, ny=5), short_config())
    assert engine.pvt is not None


def test_capillary_provider_is_accepted_by_engine():
    """A4-dən sonra: kapilyar provider mühərrikə inject edilə bilir."""
    from imex2d.domain.scal import CapillaryParameters
    from imex2d.simulation.capillary import BrooksCoreyCapillaryProvider

    scal = default_scal()
    capillary = BrooksCoreyCapillaryProvider(
        CapillaryParameters(entry_pressure=0.3), scal)
    service = SimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(scal),
        capillary_provider=capillary)
    engine = service.create_engine(five_spot_model(nx=5, ny=5), short_config())
    assert engine.capillary is not None


def test_builder_rejects_incomplete_geological_model():
    geology = SyntheticGeologicalModelBuilder().build(
        nx=5, ny=5, dx=10, dy=10, dz=5, porosity=0.2, permx_base=100.0)
    del geology.property_maps["PORO"]
    try:
        ReservoirModelBuilder().build(geology, five_spot(geology.grid))
    except ValueError:
        return
    raise AssertionError("Natamam geoloji model qəbul edildi")


def test_reservoir_model_records_source_geological_model():
    geology = SyntheticGeologicalModelBuilder().build(
        nx=5, ny=5, dx=10, dy=10, dz=5, porosity=0.2, permx_base=100.0,
        name="Geo-A")
    model = ReservoirModelBuilder().build(geology, five_spot(geology.grid))
    assert model.source_geological_model == "Geo-A"


def test_project_tracks_runs():
    project = Project("Test")
    model = five_spot_model(nx=5, ny=5)
    project.add_reservoir_model(model)
    run = project.new_run(model.name, short_config())
    assert run.run_id in project.runs
    assert project.latest_run().run_id == run.run_id
