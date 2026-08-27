"""Equilibration (A3) testləri."""

import dataclasses

import numpy as np

from helpers import default_scal, five_spot_model, make_service, short_config
from imex2d.application.simulation_service import SimulationService
from imex2d.domain.initial import InitialConditions
from imex2d.simulation.initialization.equilibrium import (GRAVITY, PA_TO_BAR,
                                                          EquilibriumInitializationProvider)
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter


def _dipping_model(nx=9, ny=9, top=2000.0, dip=5.0, owc=2050.0, scal=None):
    scal = scal or default_scal()
    model = five_spot_model(nx=nx, ny=ny, scal=scal)
    surface = np.add.outer(np.arange(ny) * 0.0, np.arange(nx) * dip) + top
    model.geometry = dataclasses.replace(model.geometry,
                                         top_depth=top,
                                         top_depth_map=surface)
    model.initial_conditions = InitialConditions(
        datum_depth=top, datum_pressure=250.0, oil_water_contact=owc,
        use_equilibration=True)
    return model


# ── həndəsə ───────────────────────────────────────────────────────────
def test_cell_depths_are_uniform_without_dip():
    model = five_spot_model(nx=5, ny=5)
    depths = model.geometry.cell_depths()
    assert np.ptp(depths) < 1e-12


def test_cell_depths_follow_dipping_surface():
    model = _dipping_model(nx=9, dip=5.0)
    depths = model.geometry.cell_depths().reshape(9, 9)
    assert np.all(np.diff(depths, axis=1) > 0), "Dərinlik X üzrə artmır"
    assert np.ptp(depths[:, 0]) < 1e-12, "Y üzrə maillik olmamalıdır"


# ── təzyiq profili ────────────────────────────────────────────────────
def test_pressure_equals_datum_pressure_at_datum_depth():
    model = _dipping_model(top=2000.0, dip=0.0, owc=9999.0)
    model.geometry = dataclasses.replace(model.geometry, dz=0.0001)
    state = EquilibriumInitializationProvider().initialize(model)
    assert abs(state.pressure.min() - 250.0) < 0.05


def test_pressure_increases_with_depth():
    model = _dipping_model(dip=5.0, owc=9999.0)
    state = EquilibriumInitializationProvider().initialize(model)
    depths = model.geometry.cell_depths()
    order = np.argsort(depths)
    assert np.all(np.diff(state.pressure[order]) >= -1e-9)


def test_oil_gradient_matches_hydrostatic_formula():
    """Neft zonasında qradiyent ρo·g·h/1e5 olmalıdır (ρo = ρ_səth/Bo)."""
    model = _dipping_model(nx=11, dip=10.0, owc=9999.0)
    state = EquilibriumInitializationProvider().initialize(model)
    depths = model.geometry.cell_depths()

    rho_reservoir = model.fluids.oil_density / model.fluids.oil_fvf
    expected_gradient = rho_reservoir * GRAVITY * PA_TO_BAR    # bar/m
    span = depths.max() - depths.min()
    actual_gradient = (state.pressure.max() - state.pressure.min()) / span
    assert abs(actual_gradient - expected_gradient) / expected_gradient < 1e-6


def test_water_zone_uses_water_density():
    """Kontaktdan aşağıda qradiyent daha dikdir (su neftdən ağırdır)."""
    model = _dipping_model(nx=21, dip=5.0, owc=2050.0)
    state = EquilibriumInitializationProvider().initialize(model)
    depths = model.geometry.cell_depths()
    oil_zone = depths < 2050.0
    water_zone = depths >= 2050.0
    assert oil_zone.any() and water_zone.any()

    def gradient(mask):
        d, p = depths[mask], state.pressure[mask]
        order = np.argsort(d)
        return np.polyfit(d[order], p[order], 1)[0]

    assert gradient(water_zone) > gradient(oil_zone)


# ── doyumluluq ────────────────────────────────────────────────────────
def test_saturation_splits_at_oil_water_contact():
    scal = default_scal()
    model = _dipping_model(nx=21, dip=5.0, owc=2050.0, scal=scal)
    state = EquilibriumInitializationProvider().initialize(model)
    depths = model.geometry.cell_depths()

    above = state.water_saturation[depths < 2050.0]
    below = state.water_saturation[depths >= 2050.0]
    assert np.allclose(above, scal.swc)
    assert np.allclose(below, 1.0 - scal.sor)


def test_missing_contact_puts_whole_model_in_oil_zone():
    scal = default_scal()
    model = _dipping_model(scal=scal)
    model.initial_conditions.oil_water_contact = None
    state = EquilibriumInitializationProvider().initialize(model)
    assert np.allclose(state.water_saturation, scal.swc)


def test_saturation_stays_inside_scal_limits():
    scal = default_scal()
    model = _dipping_model(nx=15, dip=8.0, owc=2040.0, scal=scal)
    sw = EquilibriumInitializationProvider().initialize(model).water_saturation
    assert sw.min() >= scal.swc - 1e-12
    assert sw.max() <= 1.0 - scal.sor + 1e-12


# ── PVT ilə birlikdə ──────────────────────────────────────────────────
def test_pvt_changes_reservoir_density_and_gradient():
    """Bo(p) > Bo_static olduğu üçün lay sıxlığı azalır, qradiyent yumşalır."""
    model = _dipping_model(nx=15, dip=10.0, owc=9999.0)
    table = build_pvt_table(bubble_point_bar=240.0)
    static_state = EquilibriumInitializationProvider().initialize(model)
    pvt_state = EquilibriumInitializationProvider(
        BlackOilPVTProvider(table)).initialize(model)
    assert np.ptp(pvt_state.pressure) < np.ptp(static_state.pressure)


# ── mühərriklə inteqrasiya ────────────────────────────────────────────
def test_engine_accepts_initialization_provider():
    scal = default_scal()
    model = _dipping_model(nx=15, dip=5.0, owc=2050.0, scal=scal)
    service = SimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(scal),
        initialization_provider=EquilibriumInitializationProvider())
    engine = service.create_engine(model, short_config(end_time=50.0))
    assert np.ptp(engine.pressure) > 0.0, "Təzyiq sahəsi bərabər qaldı"
    assert len(np.unique(engine.sw)) == 2, "Su zonası yaranmadı"


def test_water_zone_reduces_ooip():
    scal = default_scal()
    model = _dipping_model(nx=15, dip=5.0, owc=2040.0, scal=scal)
    service = SimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(scal),
        initialization_provider=EquilibriumInitializationProvider())
    equilibrated = service.create_engine(model, short_config(end_time=50.0))
    uniform = make_service(scal).create_engine(
        five_spot_model(nx=15, ny=15, scal=scal), short_config(end_time=50.0))
    assert equilibrated.original_oil_in_place() < uniform.original_oil_in_place()


def test_engine_runs_to_completion_with_equilibration():
    scal = default_scal()
    model = _dipping_model(nx=15, dip=4.0, owc=2045.0, scal=scal)
    service = SimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(scal),
        initialization_provider=EquilibriumInitializationProvider())
    result = service.run(model, short_config(end_time=200.0))
    assert result.converged


def test_uniform_initialization_used_when_provider_absent():
    """A3-ün zəmanəti: provider verilmədikdə köhnə yol işləyir."""
    scal = default_scal()
    engine = make_service(scal).create_engine(
        five_spot_model(nx=9, ny=9, scal=scal), short_config(end_time=10.0))
    assert engine.initialization is None
    assert np.ptp(engine.pressure) < 1e-12
    assert np.ptp(engine.sw) < 1e-12
