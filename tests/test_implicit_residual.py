"""Fully implicit qalıq vektoru (A6, mərhələ 1)."""

import numpy as np

from helpers import default_scal, five_spot_model, make_service, short_config
from imex2d.domain.initial import InitialConditions
from imex2d.domain.scal import CapillaryParameters
from imex2d.simulation.capillary import BrooksCoreyCapillaryProvider
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.implicit.residual import OIL, WATER, ResidualAssembler
from imex2d.simulation.implicit.state import (PRESSURE, VARIABLES_PER_CELL,
                                              WATER_SATURATION, ReservoirState,
                                              index_of)
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.well_model import PeacemanWellModel


def _assembler(model, scal, with_wells=True, pvt=None, capillary=None):
    grid = TwoPointFluxDiscretization().build(model)
    wells = PeacemanWellModel().build_connections(model) if with_wells else []
    return ResidualAssembler(model, grid, wells,
                             CoreyRelativePermeabilityAdapter(scal),
                             pvt=pvt, capillary=capillary)


def _uniform_state(model, pressure=None, sw=None):
    ic = model.initial_conditions
    return ReservoirState(
        np.full(model.ncell, ic.datum_pressure if pressure is None else pressure),
        np.full(model.ncell, ic.water_saturation if sw is None else sw))


# ── vəziyyət vektoru ──────────────────────────────────────────────────
def test_state_vector_round_trip():
    state = ReservoirState(np.linspace(200, 300, 7), np.linspace(0.2, 0.6, 7))
    restored = ReservoirState.from_vector(state.to_vector())
    assert np.allclose(restored.pressure, state.pressure)
    assert np.allclose(restored.water_saturation, state.water_saturation)


def test_state_vector_is_interleaved_by_cell():
    state = ReservoirState(np.array([10.0, 20.0]), np.array([0.3, 0.4]))
    assert list(state.to_vector()) == [10.0, 0.3, 20.0, 0.4]
    assert index_of(1, PRESSURE) == 2
    assert index_of(1, WATER_SATURATION) == 3


def test_newton_update_respects_saturation_limits():
    state = ReservoirState(np.array([250.0]), np.array([0.5]))
    delta = np.array([10.0, 5.0])                # nəhəng doyumluluq addımı
    updated = state.updated(delta, sw_min=0.2, sw_max=0.75)
    assert updated.water_saturation[0] == 0.75
    assert updated.pressure[0] == 260.0


def test_newton_update_chopping_limits_the_step():
    state = ReservoirState(np.array([250.0]), np.array([0.5]))
    updated = state.updated(np.array([100.0, 0.4]), 0.2, 0.9,
                            max_pressure_change=20.0,
                            max_saturation_change=0.1)
    assert updated.pressure[0] == 270.0
    assert abs(updated.water_saturation[0] - 0.6) < 1e-12


# ── qalığın fundamental xassələri ─────────────────────────────────────
def test_residual_is_zero_at_rest_without_wells():
    """Quyusuz, bərabər paylanmış vəziyyət tarazlıqdadır."""
    scal = default_scal()
    model = five_spot_model(nx=7, ny=7, scal=scal)
    assembler = _assembler(model, scal, with_wells=False)
    state = _uniform_state(model, sw=0.5)
    residual, _, _ = assembler.residual(state, state, dt=1.0)
    assert np.abs(residual).max() < 1e-12


def test_residual_length_is_two_per_cell():
    scal = default_scal()
    model = five_spot_model(nx=5, ny=5, scal=scal)
    residual, _, _ = _assembler(model, scal).residual(
        _uniform_state(model), _uniform_state(model), dt=1.0)
    assert residual.size == model.ncell * VARIABLES_PER_CELL


def test_residual_sum_equals_negative_well_rates():
    """Ən vacib xassə: daxili axınlar cəmdə bir-birini yeyir.

    Σ R_p = −Σ q_p. Əgər bu pozulsa, üzlər üzrə axın ya iki dəfə
    sayılır, ya da işarəsi səhvdir.
    """
    scal = default_scal()
    model = five_spot_model(nx=9, ny=9, scal=scal)
    assembler = _assembler(model, scal)
    state = _uniform_state(model)
    residual, _, rates = assembler.residual(state, state, dt=1.0)

    water, oil = assembler.material_balance_error(residual, dt=1.0)
    assert abs(water + rates.water.sum()) < 1e-6 * max(abs(water), 1.0)
    assert abs(oil + rates.oil.sum()) < 1e-6 * max(abs(oil), 1.0)


def test_flux_between_two_cells_is_antisymmetric():
    """Bir hüceyrədən çıxan digərinə daxil olur."""
    scal = default_scal()
    model = five_spot_model(nx=6, ny=1, scal=scal)
    assembler = _assembler(model, scal, with_wells=False)
    state = ReservoirState(np.linspace(260.0, 240.0, model.ncell),
                           np.full(model.ncell, 0.5))
    fluid = assembler.fluid_state(state)
    water, oil = assembler.net_influx(state, fluid)
    assert abs(water.sum()) < 1e-9
    assert abs(oil.sum()) < 1e-9


def test_accumulation_matches_pore_volume_times_saturation():
    scal = default_scal()
    model = five_spot_model(nx=5, ny=5, scal=scal)
    assembler = _assembler(model, scal)
    state = _uniform_state(model, sw=0.4)
    fluid = assembler.fluid_state(state)
    water, oil = assembler.accumulation(state, fluid)

    expected_water = model.pore_volume() * 0.4 / model.fluids.water_fvf
    expected_oil = model.pore_volume() * 0.6 / model.fluids.oil_fvf
    assert np.allclose(water, expected_water)
    assert np.allclose(oil, expected_oil)


def test_pressure_gradient_drives_flow_in_the_expected_direction():
    scal = default_scal()
    model = five_spot_model(nx=5, ny=1, scal=scal)
    assembler = _assembler(model, scal, with_wells=False)
    state = ReservoirState(np.linspace(260.0, 240.0, model.ncell),
                           np.full(model.ncell, 0.5))
    fluid = assembler.fluid_state(state)
    water, oil = assembler.net_influx(state, fluid)
    assert water[0] < 0 and water[-1] > 0        # yüksəkdən alçağa
    assert oil[0] < 0 and oil[-1] > 0


# ── IMPES ilə ardıcıllıq ──────────────────────────────────────────────
def test_residual_of_an_impes_step_shrinks_with_time_step():
    """İki sxem eyni fizikanı həll edir.

    IMPES doyumluluğu explicit yenilədiyi üçün onun nəticəsi implicit
    tənliyi tam ödəmir, lakin Δt kiçildikcə qalıq da azalmalıdır —
    bu, iki sxemin ardıcıllığını (consistency) təsdiqləyir.
    """
    scal = default_scal()
    model = five_spot_model(nx=9, ny=9, scal=scal)
    assembler = _assembler(model, scal)

    norms = []
    for dt in (1.0, 0.01):
        engine = make_service(scal).create_engine(model, short_config(end_time=10.0))
        previous = ReservoirState(engine.pressure.copy(), engine.sw.copy())
        pressure, lam_w, lam_o, lam_t, bw, bo, injection = engine._solve_pressure(dt)
        sw_new = engine._update_saturation(pressure, lam_w, lam_o, lam_t,
                                           bw, bo, injection, dt)[0]
        residual, _, _ = assembler.residual(ReservoirState(pressure, sw_new),
                                            previous, dt)
        norms.append(assembler.scaled_residual_norm(residual, dt))
    assert norms[1] < norms[0] * 0.5, f"Qalıq Δt ilə azalmadı: {norms}"


# ── provider-lərlə ────────────────────────────────────────────────────
def test_pvt_changes_accumulation_through_formation_volume_factor():
    scal = default_scal()
    model = five_spot_model(nx=5, ny=5, scal=scal)
    table = build_pvt_table(bubble_point_bar=240.0)
    state = _uniform_state(model)

    static = _assembler(model, scal)
    with_pvt = _assembler(model, scal, pvt=BlackOilPVTProvider(table))

    oil_static = static.accumulation(state, static.fluid_state(state))[1]
    oil_pvt = with_pvt.accumulation(state, with_pvt.fluid_state(state))[1]
    assert not np.allclose(oil_static, oil_pvt)
    assert np.all(oil_pvt < oil_static)      # Bo daha böyükdür -> səth həcmi az


def test_capillary_pressure_enters_the_water_potential():
    scal = default_scal()
    model = five_spot_model(nx=6, ny=1, scal=scal)
    capillary = BrooksCoreyCapillaryProvider(
        CapillaryParameters(entry_pressure=0.4), scal)

    state = ReservoirState(np.full(model.ncell, 250.0),
                           np.linspace(0.25, 0.7, model.ncell))
    plain = _assembler(model, scal, with_wells=False)
    with_pc = _assembler(model, scal, with_wells=False, capillary=capillary)

    plain_water, _ = plain.potentials(state, plain.fluid_state(state))
    pc_water, _ = with_pc.potentials(state, with_pc.fluid_state(state))
    assert not np.allclose(plain_water, pc_water)


def test_gravity_activates_only_when_depths_differ():
    import dataclasses
    scal = default_scal()
    flat = five_spot_model(nx=5, ny=5, scal=scal)
    assert not _assembler(flat, scal)._has_gravity

    dipping = five_spot_model(nx=5, ny=5, scal=scal)
    surface = np.add.outer(np.zeros(5), np.arange(5) * 5.0) + 2000.0
    dipping.geometry = dataclasses.replace(dipping.geometry,
                                           top_depth_map=surface)
    assert _assembler(dipping, scal)._has_gravity


# ── konvergensiya ölçüsü ──────────────────────────────────────────────
def test_scaled_norm_is_grid_size_independent():
    """Normallaşdırılmış qalıq eyni fiziki vəziyyətdə ölçüdən asılı olmamalıdır."""
    scal = default_scal()
    norms = []
    for size in (7, 15):
        model = five_spot_model(nx=size, ny=size, scal=scal)
        assembler = _assembler(model, scal, with_wells=False)
        state = _uniform_state(model, sw=0.5)
        perturbed = ReservoirState(state.pressure.copy(),
                                   state.water_saturation + 0.01)
        residual, _, _ = assembler.residual(perturbed, state, dt=1.0)
        norms.append(assembler.scaled_residual_norm(residual, dt=1.0))
    assert abs(norms[0] - norms[1]) / norms[0] < 1e-9


def test_material_balance_error_is_reported_in_surface_volume():
    scal = default_scal()
    model = five_spot_model(nx=7, ny=7, scal=scal)
    assembler = _assembler(model, scal)
    state = _uniform_state(model)
    residual, _, rates = assembler.residual(state, state, dt=2.0)
    water, _ = assembler.material_balance_error(residual, dt=2.0)
    assert abs(water - (-rates.water.sum() * 2.0)) < 1e-6
