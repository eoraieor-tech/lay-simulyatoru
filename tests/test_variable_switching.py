"""Dəyişən keçid — üç fazalı primary dəyişənlər (A7, mərhələ 4)."""

import numpy as np

from imex2d.simulation.implicit.three_phase_state import (PRESSURE,
                                                          THIRD_VARIABLE,
                                                          VARIABLES_PER_CELL,
                                                          WATER_SATURATION,
                                                          ThreePhaseState,
                                                          index_of)
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table


def _pvt():
    table = build_pvt_table(bubble_point_bar=240.0, include_gas=True)
    return BlackOilPVTProvider(table)


def _state(pressure, sw, third, saturated):
    return ThreePhaseState(np.asarray(pressure, float), np.asarray(sw, float),
                           np.asarray(third, float),
                           np.asarray(saturated, bool))


# ── əsas doyumluluq törəmələri ─────────────────────────────────────────
def test_gas_saturation_is_zero_for_undersaturated_cells():
    state = _state([200.0], [0.3], [100.0], [False])
    assert state.gas_saturation[0] == 0.0


def test_gas_saturation_equals_third_variable_when_saturated():
    state = _state([200.0], [0.3], [0.15], [True])
    assert state.gas_saturation[0] == 0.15


def test_oil_saturation_completes_the_balance():
    state = _state([200.0], [0.3], [0.15], [True])
    assert abs(state.oil_saturation[0] - 0.55) < 1e-12


def test_solution_gor_uses_saturation_curve_when_saturated():
    pvt = _pvt()
    rs_sat = float(pvt.solution_gor(np.array([200.0]))[0])
    state = _state([200.0], [0.3], [0.1], [True])
    assert abs(state.solution_gor(pvt)[0] - rs_sat) < 1e-9


def test_solution_gor_uses_third_variable_when_undersaturated():
    state = _state([200.0], [0.3], [80.0], [False])
    assert abs(state.solution_gor(_pvt())[0] - 80.0) < 1e-9


# ── dəyişən keçid: doymamış -> doymuş ──────────────────────────────────
def test_cell_switches_to_saturated_when_rs_exceeds_the_curve():
    pvt = _pvt()
    rs_sat = float(pvt.solution_gor(np.array([200.0]))[0])
    state = _state([200.0], [0.3], [rs_sat * 1.1], [False])
    switched = state.switch_variables(pvt)
    assert switched.is_saturated[0]
    assert switched.third_variable[0] == 0.0     # Sg sərhəddə sıfırdan başlayır


def test_cell_stays_undersaturated_when_rs_is_below_the_curve():
    pvt = _pvt()
    rs_sat = float(pvt.solution_gor(np.array([200.0]))[0])
    state = _state([200.0], [0.3], [rs_sat * 0.5], [False])
    switched = state.switch_variables(pvt)
    assert not switched.is_saturated[0]
    assert abs(switched.third_variable[0] - rs_sat * 0.5) < 1e-9


# ── dəyişən keçid: doymuş -> doymamış ──────────────────────────────────
def test_cell_switches_to_undersaturated_when_sg_goes_negative():
    pvt = _pvt()
    rs_sat = float(pvt.solution_gor(np.array([200.0]))[0])
    state = _state([200.0], [0.3], [-0.02], [True])
    switched = state.switch_variables(pvt)
    assert not switched.is_saturated[0]
    # kəsilməzlik: Rs sərhəd dəyərinə (Rs_sat) qayıdır
    assert abs(switched.third_variable[0] - rs_sat) < 1e-9


def test_cell_stays_saturated_when_sg_is_still_positive():
    state = _state([200.0], [0.3], [0.05], [True])
    switched = state.switch_variables(_pvt())
    assert switched.is_saturated[0]
    assert abs(switched.third_variable[0] - 0.05) < 1e-9


# ── kəsilməzlik və kütlə balansı ─────────────────────────────────────
def test_switching_preserves_mass_balance():
    """Sw+So+Sg=1 keçiddən əvvəl və sonra qorunmalıdır."""
    pvt = _pvt()
    rs_sat = float(pvt.solution_gor(np.array([200.0]))[0])
    state = _state([200.0, 200.0], [0.3, 0.3],
                   [rs_sat * 1.2, -0.03], [False, True])
    switched = state.switch_variables(pvt)
    total = (switched.water_saturation + switched.oil_saturation
            + switched.gas_saturation)
    assert np.allclose(total, 1.0)


def test_boundary_crossing_is_physically_continuous():
    """Sg=0 ⟺ Rs=Rs_sat(p) — keçid anında fiziki vəziyyət dəyişməməlidir."""
    pvt = _pvt()
    rs_sat = float(pvt.solution_gor(np.array([200.0]))[0])

    # dəqiq sərhəddə: undersaturated təsviri ilə saturated təsviri
    # eyni Sg=0 vəziyyətini göstərməlidir
    at_boundary_under = _state([200.0], [0.3], [rs_sat], [False])
    at_boundary_over = _state([200.0], [0.3], [1e-9], [True])
    assert abs(at_boundary_under.gas_saturation[0]
              - at_boundary_over.gas_saturation[0]) < 1e-6


def test_multiple_cells_switch_independently():
    pvt = _pvt()
    rs_sat = pvt.solution_gor(np.array([200.0, 220.0, 180.0]))
    state = _state([200.0, 220.0, 180.0], [0.3, 0.25, 0.35],
                   [rs_sat[0] * 1.3, rs_sat[1] * 0.7, -0.01],
                   [False, False, True])
    switched = state.switch_variables(pvt)
    assert list(switched.is_saturated) == [True, False, False]


def test_switch_is_idempotent_once_stable():
    """Sərhədi keçməyən hüceyrələrdə təkrar çağırış heç nə dəyişməməlidir."""
    state = _state([200.0], [0.3], [0.1], [True])
    once = state.switch_variables(_pvt())
    twice = once.switch_variables(_pvt())
    assert once.is_saturated[0] == twice.is_saturated[0]
    assert abs(once.third_variable[0] - twice.third_variable[0]) < 1e-12


def test_switch_does_not_touch_pressure_or_water_saturation():
    state = _state([215.0], [0.42], [0.08], [True])
    state.third_variable[0] = -0.05          # doymuşdan keçəcək
    switched = state.switch_variables(_pvt())
    assert abs(switched.pressure[0] - 215.0) < 1e-12
    assert abs(switched.water_saturation[0] - 0.42) < 1e-12


# ── vektor kodlaşdırması ────────────────────────────────────────────
def test_vector_round_trip_preserves_all_three_variables():
    state = _state([200.0, 210.0], [0.3, 0.35], [0.1, 90.0], [True, False])
    vector = state.to_vector()
    restored = ThreePhaseState.from_vector(vector, state.is_saturated)
    assert np.allclose(restored.pressure, state.pressure)
    assert np.allclose(restored.water_saturation, state.water_saturation)
    assert np.allclose(restored.third_variable, state.third_variable)


def test_vector_layout_is_interleaved_per_cell():
    """[p0,Sw0,x0, p1,Sw1,x1, ...] — CPR-in blok strukturu ilə uyğun (A6)."""
    state = _state([100.0, 200.0], [0.2, 0.4], [50.0, 0.1], [False, True])
    vector = state.to_vector()
    assert vector[index_of(0, PRESSURE)] == 100.0
    assert vector[index_of(0, WATER_SATURATION)] == 0.2
    assert vector[index_of(0, THIRD_VARIABLE)] == 50.0
    assert vector[index_of(1, PRESSURE)] == 200.0
    assert len(vector) == 2 * VARIABLES_PER_CELL


# ── Nyuton addımının tətbiqi ────────────────────────────────────────
def test_updated_applies_delta_to_all_three_variables():
    state = _state([200.0], [0.3], [0.1], [True])
    delta = np.zeros(VARIABLES_PER_CELL)
    delta[PRESSURE] = 5.0
    delta[WATER_SATURATION] = 0.02
    delta[THIRD_VARIABLE] = 0.03
    updated = state.updated(delta, sw_min=0.0, sw_max=1.0)
    assert abs(updated.pressure[0] - 205.0) < 1e-9
    assert abs(updated.water_saturation[0] - 0.32) < 1e-9
    assert abs(updated.third_variable[0] - 0.13) < 1e-9


def test_updated_clips_pressure_and_saturation_changes():
    state = _state([200.0], [0.3], [0.1], [True])
    delta = np.zeros(VARIABLES_PER_CELL)
    delta[PRESSURE] = 100.0
    delta[WATER_SATURATION] = 0.5
    updated = state.updated(delta, sw_min=0.0, sw_max=1.0,
                            max_pressure_change=20.0,
                            max_saturation_change=0.1)
    assert abs(updated.pressure[0] - 220.0) < 1e-9
    assert abs(updated.water_saturation[0] - 0.4) < 1e-9


def test_updated_does_not_limit_rs_change_for_undersaturated_cells():
    """Rs vahidi (sm3/sm3) Sg-dən fərqlidir — doyumluluq həddi tətbiq olunmamalıdır."""
    state = _state([200.0], [0.3], [80.0], [False])
    delta = np.zeros(VARIABLES_PER_CELL)
    delta[THIRD_VARIABLE] = 40.0
    updated = state.updated(delta, sw_min=0.0, sw_max=1.0,
                            max_saturation_change=0.1)
    assert abs(updated.third_variable[0] - 120.0) < 1e-9


def test_updated_preserves_is_saturated_flag():
    """`updated()` yalnız kəmiyyəti dəyişir — vəziyyəti `switch_variables()` dəyişir."""
    state = _state([200.0], [0.3], [-0.05], [True])
    delta = np.zeros(VARIABLES_PER_CELL)
    updated = state.updated(delta, sw_min=0.0, sw_max=1.0)
    assert updated.is_saturated[0]        # hələ keçməyib
