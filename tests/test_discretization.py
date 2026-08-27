"""Diskretizasiya və quyu modeli — düsturlar əl hesabı ilə yoxlanılır."""

import numpy as np

from helpers import five_spot_model
from imex2d.domain.units import METRIC
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.well_model import PeacemanWellModel


def test_homogeneous_transmissibility_matches_analytic_value():
    model = five_spot_model(nx=3, ny=3, dx=10.0, dy=10.0, dz=5.0,
                            permeability=200.0)
    grid = TwoPointFluxDiscretization().build(model)
    # Homogen halda: T = C * k * A / d,  A = dy*dz, d = dx
    expected = METRIC.darcy_constant * 200.0 * (10.0 * 5.0) / 10.0
    assert np.allclose(grid.transmissibility, expected, rtol=1e-12)


def test_transmissibility_uses_harmonic_mean():
    """Bir hüceyrənin keçiriciliyi kəskin aşağı olanda üz transmissivliyi
    harmonik ortaya uyğun gəlməlidir (arifmetik ortaya yox)."""
    model = five_spot_model(nx=2, ny=1, dx=10.0, dy=10.0, dz=5.0)
    model.rock.permx.values[:] = np.array([1000.0, 10.0])
    model.rock.permy.values[:] = model.rock.permx.values
    grid = TwoPointFluxDiscretization().build(model)

    k_harmonic = 2.0 / (1.0 / 1000.0 + 1.0 / 10.0)
    expected = METRIC.darcy_constant * k_harmonic * (10.0 * 5.0) / 10.0
    assert np.allclose(grid.transmissibility, expected, rtol=1e-10)
    arithmetic = METRIC.darcy_constant * 505.0 * (10.0 * 5.0) / 10.0
    assert grid.transmissibility[0] < arithmetic * 0.1


def test_pore_volume_from_discretization_matches_model():
    model = five_spot_model(nx=6, ny=5)
    grid = TwoPointFluxDiscretization().build(model)
    assert np.allclose(grid.pore_volume, model.pore_volume())


def test_peaceman_well_index_matches_hand_calculation():
    model = five_spot_model(nx=5, ny=5, dx=20.0, dy=20.0, dz=10.0,
                            permeability=150.0)
    connections = PeacemanWellModel().build_connections(model)
    assert len(connections) == 2

    # İzotrop kvadrat hüceyrə: re = 0.28 * sqrt(2*dx^2) / 2 = 0.198 * dx
    re = 0.28 * np.sqrt(2 * 20.0 ** 2) / 2.0
    expected = (METRIC.darcy_constant * 2.0 * np.pi * 150.0 * 10.0
                / np.log(re / 0.1))
    assert abs(connections[0].well_index - expected) / expected < 1e-10


def test_skin_reduces_well_index():
    model = five_spot_model(nx=5, ny=5)
    base = PeacemanWellModel().build_connections(model)[0].well_index
    for well in model.wells:
        well.perforations[0].skin = 5.0
    damaged = PeacemanWellModel().build_connections(model)[0].well_index
    assert damaged < base
