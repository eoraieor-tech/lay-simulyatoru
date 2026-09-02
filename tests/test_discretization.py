"""Diskretizasiya və quyu modeli — düsturlar əl hesabı ilə yoxlanılır."""

import numpy as np

from helpers import five_spot_model
from imex2d.domain.properties import PermeabilityTensor, PropertyMap
from imex2d.domain.units import METRIC
from imex2d.interfaces.discretization import IFluxDiscretization
from imex2d.simulation.discretization import (TwoPointFluxDiscretization,
                                              default_flux_discretization)
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


# ── MPFA arxitektura hazırlığı (Phase: Numerical Discretization
# Architecture Preparation) — TPFA-nın ÖZÜ dəyişməyib, bax §9/§12 ────────
def test_tpfa_implements_flux_discretization_interface():
    assert isinstance(TwoPointFluxDiscretization(), IFluxDiscretization)


def test_default_flux_discretization_is_tpfa():
    assert isinstance(default_flux_discretization(), TwoPointFluxDiscretization)


def test_discretized_grid_compute_flux_matches_transmissibility_times_delta():
    """`compute_flux` — gələcək MPFA-O-nun əvəz edəcəyi tək giriş nöqtəsi
    — TPFA üçün riyazi olaraq `T · ΔΦ`-dən FƏRQLƏNMƏMƏLİDİR."""
    model = five_spot_model(nx=3, ny=3, dx=10.0, dy=10.0, dz=5.0, permeability=200.0)
    grid = TwoPointFluxDiscretization().build(model)
    d_phi = np.linspace(-5.0, 5.0, grid.connections.count)
    assert np.allclose(grid.compute_flux(d_phi), grid.transmissibility * d_phi)


def test_flux_discretization_conserves_locally():
    """Audit §8: quyu/mənbə olmadıqda, daxili üz axınlarının ÜMUMİ CƏMİ
    bütün grid üzrə sıfır olmalıdır — hər üz `compute_flux`-dan TƏK bir
    ədəd alır və bu, əks işarə ilə iki tərəfə yazılır (bax
    `ResidualAssembler.net_influx`-un eyni pattern-i), ona görə daxili
    mübadilə HEÇ VAXT kütlə yaratmır/itirmir."""
    model = five_spot_model(nx=4, ny=4, dx=15.0, dy=15.0, dz=8.0, permeability=120.0)
    grid = TwoPointFluxDiscretization().build(model)
    rng = np.random.default_rng(0)
    pressure = rng.uniform(180.0, 220.0, size=model.ncell)
    conn = grid.connections
    d_phi = pressure[conn.cell_a] - pressure[conn.cell_b]
    flux = grid.compute_flux(d_phi)

    net = np.zeros(model.ncell)
    np.add.at(net, conn.cell_a, -flux)
    np.add.at(net, conn.cell_b, +flux)
    assert abs(net.sum()) < 1e-8, "Daxili axınların ümumi cəmi sıfır olmalıdır"


def test_full_tensor_permeability_triggers_explicit_warning_not_silent_scalarization():
    """Audit §5: TPFA tam tenzoru SƏSSİZCƏ diaqonala yumşaltmamalı —
    AÇIQ xəbərdarlıq verməli, VƏ öz nəticəsini DƏYİŞMƏMƏLİDİR (yalnız
    diaqonaldan istifadə etməyə davam etməlidir)."""
    baseline = five_spot_model(nx=3, ny=3, dx=10.0, dy=10.0, dz=5.0, permeability=200.0)
    baseline_grid = TwoPointFluxDiscretization().build(baseline)
    assert baseline_grid.warnings == []

    model = five_spot_model(nx=3, ny=3, dx=10.0, dy=10.0, dz=5.0, permeability=200.0)
    ncell = model.ncell
    model.rock.permeability_tensor = PermeabilityTensor(
        kxx=PropertyMap.uniform("KXX", 200.0, ncell),
        kyy=PropertyMap.uniform("KYY", 200.0, ncell),
        kzz=PropertyMap.uniform("KZZ", 20.0, ncell),
        kxy=PropertyMap.uniform("KXY", 50.0, ncell))

    grid = TwoPointFluxDiscretization().build(model)
    assert any("MPFA" in w for w in grid.warnings)
    # TPFA-nın öz nəticəsi dəyişməyib — tenzor İSTİFADƏ OLUNMUR, sadəcə bildirilir
    assert np.allclose(grid.transmissibility, baseline_grid.transmissibility)


def test_permeability_tensor_has_off_diagonal_detection():
    tensor_isotropic = PermeabilityTensor(
        kxx=PropertyMap.uniform("KXX", 100.0, 4),
        kyy=PropertyMap.uniform("KYY", 100.0, 4),
        kzz=PropertyMap.uniform("KZZ", 10.0, 4))
    assert not tensor_isotropic.has_off_diagonal()

    tensor_full = PermeabilityTensor(
        kxx=PropertyMap.uniform("KXX", 100.0, 4),
        kyy=PropertyMap.uniform("KYY", 100.0, 4),
        kzz=PropertyMap.uniform("KZZ", 10.0, 4),
        kxz=PropertyMap.uniform("KXZ", 5.0, 4))
    assert tensor_full.has_off_diagonal()
