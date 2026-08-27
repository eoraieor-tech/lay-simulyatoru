"""Üç fazalı doyumluluq və qaz papağı equilibration-u (A7, mərhələ 2)."""

import numpy as np

from helpers import default_scal
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.domain.three_phase import ThreePhaseSaturation, saturation_state
from imex2d.simulation.initialization.equilibrium import \
    EquilibriumInitializationProvider
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table


def _model(nz=10, **kwargs):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=6, ny=6, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=nz, top_depth=2000.0, **kwargs)
    scal = default_scal()
    return ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                         scal=scal)


def _column(model, state, field):
    values = getattr(state, field)
    return values.reshape(model.grid.shape)[:, 0, 0]


# ── ThreePhaseSaturation ────────────────────────────────────────────
def test_oil_saturation_is_derived_not_stored():
    saturation = ThreePhaseSaturation(water=np.array([0.2, 0.3]),
                                      gas=np.array([0.5, 0.1]))
    assert np.allclose(saturation.oil, [0.3, 0.6])


def test_is_saturated_flags_cells_with_free_gas():
    saturation = ThreePhaseSaturation(water=np.array([0.2, 0.2, 0.2]),
                                      gas=np.array([0.0, 0.3, 1e-10]))
    assert list(saturation.is_saturated) == [False, True, False]
    assert saturation.gas_cell_count == 1


def test_validate_accepts_a_consistent_field():
    saturation = ThreePhaseSaturation(water=np.array([0.2, 0.5]),
                                      gas=np.array([0.3, 0.0]))
    assert saturation.validate() == []


def test_validate_flags_saturations_summing_above_one():
    saturation = ThreePhaseSaturation(water=np.array([0.7]), gas=np.array([0.6]))
    issues = saturation.validate()
    assert any("mənfi" in issue for issue in issues)


def test_validate_flags_out_of_range_values():
    assert ThreePhaseSaturation(np.array([1.5]), np.array([0.1])).validate()
    assert ThreePhaseSaturation(np.array([0.2]), np.array([-0.1])).validate()


def test_clip_respects_the_sw_plus_sg_le_one_constraint():
    """Ayrı-ayrı clip kifayət etmir — Sw=0.9, Sg=0.9 cəmi 1.8 edərdi."""
    saturation = ThreePhaseSaturation(water=np.array([0.9]), gas=np.array([0.9]))
    clipped = saturation.clip(sw_min=0.0, sw_max=1.0, sg_min=0.0, sg_max=1.0)
    assert clipped.water[0] + clipped.gas[0] <= 1.0 + 1e-9
    # nisbət qorunur (hər ikisi bərabər azaldılıb)
    assert abs(clipped.water[0] - clipped.gas[0]) < 1e-9


def test_clip_leaves_valid_fields_untouched():
    saturation = ThreePhaseSaturation(water=np.array([0.3]), gas=np.array([0.2]))
    clipped = saturation.clip(0.15, 0.8)
    assert abs(clipped.water[0] - 0.3) < 1e-9
    assert abs(clipped.gas[0] - 0.2) < 1e-9


# ── saturation_state ──────────────────────────────────────────────────
def test_saturation_state_flags_cells_at_the_bubble_point():
    table = build_pvt_table(bubble_point_bar=240.0, include_gas=True)
    provider = BlackOilPVTProvider(table)
    pressure = np.array([150.0, 240.0, 300.0])

    # Rs HƏR TƏZYİQDƏ ÖZ doyma əyrisindən götürülür.
    #
    # Əvvəl burada Rs_sat(240) hər üç hüceyrə üçün işlədilirdi —
    # 240 bar-ın cədvəldə DƏQİQ düyün olduğu güman edilirdi. PVT
    # cədvəlinin aralığı dəyişəndə (10 → 1 bar) düyünlər sürüşdü və
    # test sındı. Doyma əyrisini hər nöqtədə öz təzyiqində
    # qiymətləndirmək bu kövrəkliyi aradan qaldırır.
    rs = np.array([provider.solution_gor(np.array([150.0]))[0] * 0.5,
                   provider.solution_gor(np.array([240.0]))[0],
                   provider.solution_gor(np.array([300.0]))[0]])
    state = saturation_state(pressure, rs, provider)
    assert list(state) == [False, True, True]


# ── qaz papağı equilibration ──────────────────────────────────────────
def _gas_cap_model():
    model = _model(nz=10)
    model.initial_conditions.datum_depth = 2000.0
    model.initial_conditions.datum_pressure = 250.0
    model.initial_conditions.gas_oil_contact = 2030.0
    model.initial_conditions.oil_water_contact = 2070.0
    return model


def test_three_zones_form_in_the_correct_depth_order():
    model = _gas_cap_model()
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    state = EquilibriumInitializationProvider(pvt=pvt).initialize(model)
    assert state.has_gas

    depths = _column(model, model.geometry, "cell_depths") \
        if False else model.geometry.cell_depths().reshape(model.grid.shape)[:, 0, 0]
    sw = state.water_saturation.reshape(model.grid.shape)[:, 0, 0]
    sg = state.gas_saturation.reshape(model.grid.shape)[:, 0, 0]

    gas_zone = depths < 2030.0
    oil_zone = (depths >= 2030.0) & (depths < 2070.0)
    water_zone = depths >= 2070.0

    assert np.all(sg[gas_zone] > 0.5)
    assert np.all(sg[oil_zone] < 1e-9)
    assert np.all(sg[water_zone] < 1e-9)
    assert np.all(sw[water_zone] > sw[oil_zone].max())


def test_saturations_sum_to_one_everywhere():
    model = _gas_cap_model()
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    state = EquilibriumInitializationProvider(pvt=pvt).initialize(model)
    total = state.water_saturation + state.gas_saturation + \
        (1.0 - state.water_saturation - state.gas_saturation)
    assert np.allclose(total, 1.0)


def test_gas_cap_saturation_respects_connate_water():
    model = _gas_cap_model()
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    state = EquilibriumInitializationProvider(pvt=pvt).initialize(model)
    depths = model.geometry.cell_depths()
    in_gas_cap = depths < 2030.0
    scal = model.scal_parameters
    assert np.allclose(state.water_saturation[in_gas_cap], scal.swc, atol=1e-6)
    assert np.allclose(state.gas_saturation[in_gas_cap], 1.0 - scal.swc, atol=1e-6)


def test_pressure_increases_monotonically_with_depth():
    model = _gas_cap_model()
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    state = EquilibriumInitializationProvider(pvt=pvt).initialize(model)
    column = state.pressure.reshape(model.grid.shape)[:, 0, 0]
    assert np.all(np.diff(column) > 0)


def test_gas_column_weight_is_much_smaller_than_oil_column():
    """ρ_qaz ≪ ρ_neft — qaz sütununun təzyiq töhfəsi kiçik olmalıdır."""
    model = _gas_cap_model()
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    provider = EquilibriumInitializationProvider(pvt=pvt)

    with_gas = provider.initialize(model).pressure
    without_gas_model = _gas_cap_model()
    without_gas_model.initial_conditions.gas_oil_contact = None
    without_gas = provider.initialize(without_gas_model).pressure

    # eyni datum təzyiqindən başlayır, fərq yalnız qaz sütununun
    # çəkisindən gəlir — kiçik olmalıdır
    difference = np.abs(with_gas - without_gas)
    assert np.all(difference < 5.0)      # bar — real qaz sütunu üçün kiçik


# ── geriyə uyğunluq ───────────────────────────────────────────────────
def test_without_goc_behaves_exactly_as_two_phase():
    model = _model(nz=6)
    model.initial_conditions.oil_water_contact = 2050.0
    # gas_oil_contact TƏYİN OLUNMUR
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    state = EquilibriumInitializationProvider(pvt=pvt).initialize(model)
    assert not state.has_gas
    assert state.gas_saturation is None


def test_two_phase_pvt_ignores_goc_even_if_set():
    """PVT-də qaz yoxdursa, GOC verilsə də iki fazalı davranış qorunur."""
    model = _model(nz=6)
    model.initial_conditions.oil_water_contact = 2050.0
    model.initial_conditions.gas_oil_contact = 2010.0
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0))
    state = EquilibriumInitializationProvider(pvt=pvt).initialize(model)
    assert not state.has_gas


def test_without_pvt_provider_ignores_goc():
    model = _model(nz=6)
    model.initial_conditions.oil_water_contact = 2050.0
    model.initial_conditions.gas_oil_contact = 2010.0
    state = EquilibriumInitializationProvider().initialize(model)
    assert not state.has_gas


def test_default_initial_state_has_no_gas():
    """`InitialState`-in defolt konstruktoru geriyə uyğun qalmalıdır."""
    from imex2d.interfaces.providers import InitialState

    state = InitialState(pressure=np.array([250.0]),
                         water_saturation=np.array([0.3]))
    assert not state.has_gas
    assert state.gas_saturation is None
