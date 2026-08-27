"""Quyu perforasiya intervalı — qismən açılmış quyular."""

import numpy as np

from helpers import SKIP_SLOW, default_scal
from imex2d.application.config import OutputConfig, SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder
from imex2d.application.simulation_service import SimulationService
from imex2d.domain.initial import InitialConditions
from imex2d.domain.wells import (ControlMode, Perforation, Well, WellControl,
                                 WellType)
from imex2d.simulation.initialization.equilibrium import EquilibriumInitializationProvider
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.well_model import PeacemanWellModel


def _layered_model(producer_layers, nx=7, ny=7, nz=6, owc=2016.0, scal=None):
    """Alt təbəqələri su zonasında olan maili olmayan 3D model."""
    scal = scal or default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=30.0, dy=30.0, dz=4.0, porosity=0.22,
        permx_base=200.0, nz=nz, kv_over_kh=0.2, top_depth=2000.0)
    wells = [
        Well("INJ-1", WellType.INJECTOR, WellControl(ControlMode.BHP, 320.0),
             [Perforation(0, 0, k) for k in range(nz)]),
        Well("PROD-1", WellType.PRODUCER, WellControl(ControlMode.BHP, 150.0),
             [Perforation(nx - 1, ny - 1, k) for k in producer_layers]),
    ]
    return ReservoirModelBuilder().build(
        geology, wells, scal=scal,
        initial=InitialConditions(datum_depth=2000.0, datum_pressure=250.0,
                                  oil_water_contact=owc, use_equilibration=True))


def _service(scal):
    return SimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(scal),
        initialization_provider=EquilibriumInitializationProvider())


# ── bağlantılar ───────────────────────────────────────────────────────
def test_connection_count_matches_perforated_layers():
    model = _layered_model(producer_layers=[0, 1, 2], nz=6)
    connections = PeacemanWellModel().build_connections(model)
    producer = [c for c in connections if c.well_name == "PROD-1"]
    injector = [c for c in connections if c.well_name == "INJ-1"]
    assert len(producer) == 3
    assert len(injector) == 6


def test_closed_perforations_are_excluded():
    model = _layered_model(producer_layers=list(range(6)), nz=6)
    producer = next(w for w in model.wells if w.name == "PROD-1")
    for perforation in producer.perforations[3:]:
        perforation.open = False
    connections = PeacemanWellModel().build_connections(model)
    assert len([c for c in connections if c.well_name == "PROD-1"]) == 3


def test_perforated_cells_are_the_expected_layers():
    model = _layered_model(producer_layers=[2, 3], nz=6)
    connections = PeacemanWellModel().build_connections(model)
    layers = sorted(model.grid.ijk(c.cell)[2]
                    for c in connections if c.well_name == "PROD-1")
    assert layers == [2, 3]


def test_perforation_outside_grid_is_rejected():
    model = _layered_model(producer_layers=[0], nz=4)
    producer = next(w for w in model.wells if w.name == "PROD-1")
    producer.perforations = [Perforation(0, 0, 9)]
    assert any("kənar" in issue for issue in model.validate())


# ── fiziki nəticə ─────────────────────────────────────────────────────
def test_partial_completion_reduces_water_production():
    """Su zonasındakı təbəqələri bağlamaq sulaşmanı azaldır.

    Bu, EOR-dan əvvəlki ən sadə optimallaşdırmadır: quyu daha az neft
    verir, amma neft/su nisbəti xeyli yaxşılaşır — su emalı xərci
    hasilatın iqtisadi həddini müəyyən edən əsas amillərdəndir.
    """
    if SKIP_SLOW:
        return
    scal = default_scal()
    config = SimulationConfig(end_time=300.0,
                              output=OutputConfig(snapshot_count=4))

    full = _service(scal).run(_layered_model(list(range(6)), scal=scal), config)
    partial = _service(scal).run(_layered_model([0, 1, 2], scal=scal), config)

    assert partial.series.water_cut[-1] < full.series.water_cut[-1], \
        "Qismən perforasiya sulaşmanı azaltmadı"
    assert partial.series.cumulative_water[-1] < full.series.cumulative_water[-1]

    def oil_to_water(result):
        return result.series.cumulative_oil[-1] / max(
            result.series.cumulative_water[-1], 1e-9)

    assert oil_to_water(partial) > oil_to_water(full), \
        "Neft/su nisbəti yaxşılaşmadı"
    # daha az təbəqə -> daha az ümumi hasilat
    assert partial.series.cumulative_oil[-1] < full.series.cumulative_oil[-1]


def test_upper_layers_only_still_produces_oil():
    if SKIP_SLOW:
        return
    scal = default_scal()
    result = _service(scal).run(
        _layered_model([0, 1], scal=scal),
        SimulationConfig(end_time=150.0, output=OutputConfig(snapshot_count=3)))
    assert result.converged
    assert result.series.cumulative_oil[-1] > 0.0
