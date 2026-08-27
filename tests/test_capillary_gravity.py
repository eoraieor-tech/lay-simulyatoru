"""Kapilyar təzyiq və cazibə (A4) testləri."""

import dataclasses

import numpy as np

from helpers import default_scal, five_spot_model, make_service, short_config
from imex2d.application.simulation_service import SimulationService
from imex2d.domain.initial import InitialConditions
from imex2d.domain.scal import CapillaryParameters
from imex2d.simulation.capillary import BrooksCoreyCapillaryProvider
from imex2d.simulation.initialization.equilibrium import (GRAVITY, PA_TO_BAR,
                                                          EquilibriumInitializationProvider)
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from test_initialization import _dipping_model


def _capillary(scal=None, entry=0.3, lam=2.0, pc_max=4.0):
    scal = scal or default_scal()
    return BrooksCoreyCapillaryProvider(
        CapillaryParameters(entry_pressure=entry, lambda_exponent=lam,
                            max_pressure=pc_max), scal)


# ── Brooks-Corey modeli ───────────────────────────────────────────────
def test_capillary_pressure_decreases_with_water_saturation():
    scal = default_scal()
    provider = _capillary(scal)
    sw = np.linspace(scal.swc, 1.0 - scal.sor, 40)
    pc = provider.pcow(sw)
    assert np.all(np.diff(pc) <= 1e-12), "Pc doyumluluqla azalmır"


def test_capillary_pressure_equals_entry_pressure_at_maximum_saturation():
    scal = default_scal()
    provider = _capillary(scal, entry=0.3)
    assert abs(provider.pcow(1.0 - scal.sor) - 0.3) < 1e-9


def test_capillary_pressure_is_capped():
    provider = _capillary(entry=0.3, pc_max=2.0)
    assert provider.pcow(default_scal().swc) <= 2.0 + 1e-12


def test_inversion_round_trips():
    scal = default_scal()
    provider = _capillary(scal)
    sw = np.linspace(scal.swc + 0.02, 1.0 - scal.sor, 30)
    assert np.allclose(provider.saturation_from_pc(provider.pcow(sw)), sw, atol=1e-6)


def test_analytic_derivative_matches_numeric():
    scal = default_scal()
    provider = _capillary(scal)
    sw = np.linspace(scal.swc + 0.05, 1.0 - scal.sor - 0.01, 20)
    h = 1e-6
    numeric = (provider.pcow(sw + h) - provider.pcow(sw - h)) / (2 * h)
    analytic = provider.dpcow_dsw(sw)
    assert np.max(np.abs(numeric - analytic) / np.abs(numeric)) < 1e-4


def test_disabled_model_is_reported_as_disabled():
    assert not CapillaryParameters().enabled
    assert CapillaryParameters(entry_pressure=0.2).enabled


def test_invalid_capillary_parameters_rejected():
    assert CapillaryParameters(entry_pressure=1.0, max_pressure=0.5).validate()
    assert CapillaryParameters(lambda_exponent=-1.0).validate()


# ── equilibration keçid zonası ────────────────────────────────────────
def test_transition_zone_is_smooth_above_contact():
    scal = default_scal()
    model = _dipping_model(nx=21, ny=1, dip=4.0, owc=2050.0, scal=scal)
    state = EquilibriumInitializationProvider(
        capillary=_capillary(scal, entry=0.25)).initialize(model)
    depths = model.geometry.cell_depths()
    above = depths < 2050.0
    sw_above = state.water_saturation[above]
    assert len(np.unique(np.round(sw_above, 4))) > 3, "Keçid zonası yaranmadı"
    order = np.argsort(depths[above])
    assert np.all(np.diff(sw_above[order]) >= -1e-9), "Sw dərinliklə artmır"


def test_sharp_contact_without_capillary_provider():
    scal = default_scal()
    model = _dipping_model(nx=21, ny=1, dip=4.0, owc=2050.0, scal=scal)
    sw = EquilibriumInitializationProvider().initialize(model).water_saturation
    assert len(np.unique(np.round(sw, 6))) == 2


def test_transition_zone_height_matches_capillary_gravity_balance():
    """Pc_max-a uyğun sütun hündürlüyü: h = Pc·1e5 / (Δρ·g)."""
    scal = default_scal()
    provider = _capillary(scal, entry=0.25, pc_max=4.0)
    model = _dipping_model(nx=41, ny=1, dip=3.0, owc=2130.0, scal=scal)
    state = EquilibriumInitializationProvider(capillary=provider).initialize(model)
    depths = model.geometry.cell_depths()

    delta_rho = (model.fluids.water_density / model.fluids.water_fvf
                 - model.fluids.oil_density / model.fluids.oil_fvf)
    expected_height = provider.params.entry_pressure / (delta_rho * GRAVITY * PA_TO_BAR)

    # Pe-dən aşağı kapilyar təzyiqdə Sw maksimumdadır (kontaktdan bu hündürlüyə qədər)
    saturated = state.water_saturation >= (1.0 - scal.sor) - 1e-6
    top_of_saturated_zone = depths[saturated].min()
    actual_height = 2130.0 - top_of_saturated_zone
    assert abs(actual_height - expected_height) < 6.0


# ── cazibə axını ──────────────────────────────────────────────────────
def test_flat_model_has_no_gravity_term():
    engine = make_service(default_scal()).create_engine(
        five_spot_model(nx=9, ny=9), short_config(end_time=10.0))
    assert not engine._has_gravity


def test_dipping_model_activates_gravity_term():
    scal = default_scal()
    model = _dipping_model(nx=9, ny=9, dip=5.0, owc=9999.0, scal=scal)
    service = SimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(scal),
        initialization_provider=EquilibriumInitializationProvider())
    engine = service.create_engine(model, short_config(end_time=10.0))
    assert engine._has_gravity


def test_gravity_drives_water_downdip_at_oil_equilibrium():
    """Neft hidrostatik tarazlıqda olanda su hələ də aşağı sürüklənir.

    Φ_o ≈ const → dΦ_o ≈ 0, lakin ρw > ρo olduğu üçün dΦ_w > 0 qalır.
    Bu, cazibə seqreqasiyasının hərəkətverici qüvvəsidir.
    """
    scal = default_scal()
    model = _dipping_model(nx=21, ny=1, dip=6.0, owc=9999.0, scal=scal)
    service = SimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(scal),
        initialization_provider=EquilibriumInitializationProvider())
    engine = service.create_engine(model, short_config(end_time=10.0))
    engine.sw[:] = 0.5     # hər iki faza hərəkətli

    mu_w, mu_o, bw, bo, _ = engine._fluid_state(engine.pressure, engine.sw)
    lam_w, lam_o = engine._mobilities(engine.sw, mu_w, mu_o)
    d_phi_o, d_phi_w = engine._face_potential_terms(
        engine.pressure, engine.sw, lam_w, lam_o, bw, bo)

    depths = model.geometry.cell_depths()
    conn = model.connections()
    downdip = depths[conn.cell_b] > depths[conn.cell_a]
    assert downdip.any()

    assert np.allclose(d_phi_o[downdip], 0.0, atol=1e-3), \
        "Neft tarazlıqda deyil"
    assert np.all(d_phi_w[downdip] > 1e-4), \
        "Su üçün aşağıya doğru potensial fərqi yaranmadı"


def test_capillary_gradient_opposes_water_advance():
    """Kapilyar üzv su cəbhəsinin qarşısında əks istiqamətdə işləyir."""
    scal = default_scal()
    model = five_spot_model(nx=11, ny=1, scal=scal)
    capillary = _capillary(scal, entry=0.4)
    service = SimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(scal),
        capillary_provider=capillary)
    engine = service.create_engine(model, short_config(end_time=10.0))

    engine.pressure[:] = 250.0
    engine.sw[:] = scal.swc
    engine.sw[:4] = 0.6                      # sol tərəfdə su cəbhəsi

    mu_w, mu_o, bw, bo, _ = engine._fluid_state(engine.pressure, engine.sw)
    lam_w, lam_o = engine._mobilities(engine.sw, mu_w, mu_o)
    _, d_phi_w = engine._face_potential_terms(
        engine.pressure, engine.sw, lam_w, lam_o, bw, bo)

    conn = model.connections()
    front = (engine.sw[conn.cell_a] > engine.sw[conn.cell_b] + 0.05)
    assert front.any()
    # Pc quru tərəfdə daha yüksəkdir -> suyu irəli çəkir (dΦw > 0)
    assert np.all(d_phi_w[front] > 0.0)


def test_engine_runs_with_capillary_and_gravity_together():
    scal = default_scal()
    model = _dipping_model(nx=15, ny=15, dip=4.0, owc=2045.0, scal=scal)
    capillary = _capillary(scal, entry=0.25)
    service = SimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(scal),
        capillary_provider=capillary,
        initialization_provider=EquilibriumInitializationProvider(capillary=capillary))
    result = service.run(model, short_config(end_time=200.0))
    assert result.converged
    for snapshot in result.snapshots:
        sw = snapshot.water_saturation
        assert sw.min() >= scal.swc - 1e-9
        assert sw.max() <= 1.0 - scal.sor + 1e-9


def test_material_balance_holds_with_capillary_and_gravity():
    scal = default_scal()
    model = _dipping_model(nx=15, ny=15, dip=4.0, owc=2045.0, scal=scal)
    capillary = _capillary(scal, entry=0.25)
    service = SimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(scal),
        capillary_provider=capillary,
        initialization_provider=EquilibriumInitializationProvider(capillary=capillary))
    engine = service.create_engine(model, short_config(end_time=200.0))
    initial_water = float(np.sum(model.pore_volume() * engine.sw))
    result = engine.run()

    final_water = float(np.sum(model.pore_volume() * engine.sw))
    series = result.series
    injected = float(np.trapezoid(series.water_injection_rate, series.time))
    produced = float(np.trapezoid(series.water_rate, series.time)) * model.fluids.water_fvf
    error = abs((final_water - initial_water) - (injected - produced)) / max(injected, 1e-9)
    assert error < 0.005, f"Material balans xətası {error * 100:.3f} %"
