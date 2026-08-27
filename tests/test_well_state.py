"""Quyu naməlumları — OPM tipli quyu modeli, MƏRHƏLƏ 1.

Bu mərhələdə yalnız VƏZİYYƏT var (naməlumlar, vektor yerləşməsi).
Tənliklər, Jakobian və Nyuton inteqrasiyası sonrakı mərhələlərdədir.
"""

import numpy as np

from helpers import default_scal
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.domain.wells import ControlMode
from imex2d.simulation.implicit.three_phase_state import ThreePhaseState
from imex2d.simulation.implicit.well_state import CoupledState, WellUnknowns
from imex2d.simulation.well_model import PeacemanWellModel


def _model(nx=4, ny=4):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    return ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                         scal=default_scal())


def _setup():
    model = _model()
    connections = PeacemanWellModel().build_connections(model)
    n = model.ncell
    reservoir = ThreePhaseState(np.full(n, 213.5), np.full(n, 0.35),
                                np.zeros(n), np.zeros(n, bool))
    wells = WellUnknowns.from_connections(connections, reservoir.pressure)
    return model, connections, reservoir, wells


# ── WellUnknowns ─────────────────────────────────────────────────────
def test_one_unknown_per_well_not_per_perforation():
    """Quyu BİR naməlumdur (BHP) — perforasiya sayından asılı deyil."""
    model, connections, reservoir, wells = _setup()
    distinct = {c.well_name for c in connections}
    assert wells.count == len(distinct)
    assert len(connections) >= wells.count


def test_bhp_controlled_wells_start_at_their_target():
    """BHP idarəli quyu üçün ən yaxşı ilkin qiymət hədəfin özüdür."""
    model, connections, reservoir, wells = _setup()
    for connection in connections:
        if connection.mode is ControlMode.BHP:
            assert abs(wells.bhp_of(connection.well_name)
                      - connection.target) < 1e-9


def test_rate_controlled_wells_start_at_the_cell_pressure():
    """RATE idarəli quyuda hədəf BHP yoxdur — hüceyrə təzyiqi işlədilir."""
    from imex2d.simulation.well_model import WellConnection

    model, connections, reservoir, _ = _setup()
    producer = next(c for c in connections if not c.is_injector)
    rate_connection = WellConnection(
        well_name=producer.well_name, cell=producer.cell,
        well_index=producer.well_index, is_injector=False,
        mode=ControlMode.RATE, target=-50.0)
    wells = WellUnknowns.from_connections([rate_connection],
                                          reservoir.pressure)
    assert abs(wells.bhp_of(producer.well_name)
              - reservoir.pressure[producer.cell]) < 1e-9


def test_explicit_initial_bhp_overrides_the_default():
    model, connections, reservoir, _ = _setup()
    wells = WellUnknowns.from_connections(
        connections, reservoir.pressure, initial_bhp={"PROD-1": 177.0})
    assert abs(wells.bhp_of("PROD-1") - 177.0) < 1e-9


def test_mismatched_names_and_values_are_rejected():
    try:
        WellUnknowns(["A", "B"], np.array([100.0]))
    except ValueError:
        return
    raise AssertionError("uyğunsuz ölçü qəbul edildi")


def test_copy_is_independent():
    model, connections, reservoir, wells = _setup()
    duplicate = wells.copy()
    duplicate.bhp[0] += 50.0
    assert wells.bhp[0] != duplicate.bhp[0]


# ── CoupledState: vektor yerləşməsi ─────────────────────────────────
def test_vector_size_is_reservoir_plus_wells():
    model, connections, reservoir, wells = _setup()
    state = CoupledState(reservoir, wells)
    assert state.size == model.ncell * 3 + wells.count


def test_wells_are_appended_after_the_reservoir_block():
    """Rezervuarın 3×3 blok strukturu TOXUNULMAMALIDIR — quyular sona
    əlavə olunur (CPR ön-şərtçisi bu quruluşa əsaslanır)."""
    model, connections, reservoir, wells = _setup()
    state = CoupledState(reservoir, wells)
    assert state.well_offset == model.ncell * 3
    vector = state.to_vector()
    assert np.allclose(vector[:state.well_offset], reservoir.to_vector())
    assert np.allclose(vector[state.well_offset:], wells.bhp)


def test_global_well_index_points_at_the_right_slot():
    model, connections, reservoir, wells = _setup()
    state = CoupledState(reservoir, wells)
    vector = state.to_vector()
    for name in wells.names:
        assert abs(vector[state.well_index(name)] - wells.bhp_of(name)) < 1e-9


def test_vector_round_trip_preserves_everything():
    model, connections, reservoir, wells = _setup()
    state = CoupledState(reservoir, wells)
    restored = CoupledState.from_vector(state.to_vector(),
                                        reservoir.is_saturated, wells.names)
    assert np.allclose(restored.to_vector(), state.to_vector())
    assert np.allclose(restored.reservoir.pressure, reservoir.pressure)
    assert np.allclose(restored.wells.bhp, wells.bhp)
    assert restored.wells.names == wells.names


def test_from_vector_rejects_an_inconsistent_length():
    model, connections, reservoir, wells = _setup()
    state = CoupledState(reservoir, wells)
    broken = np.zeros(state.size + 1)
    try:
        CoupledState.from_vector(broken, reservoir.is_saturated, wells.names)
    except ValueError:
        return
    raise AssertionError("uyğunsuz uzunluq qəbul edildi")


# ── Nyuton addımının tətbiqi ────────────────────────────────────────
def test_update_applies_reservoir_and_well_parts_separately():
    model, connections, reservoir, wells = _setup()
    state = CoupledState(reservoir, wells)
    delta = np.zeros(state.size)
    delta[0] = 4.0                                    # hüceyrə 0 təzyiqi
    delta[state.well_index("PROD-1")] = 9.0           # quyu BHP-si

    updated = state.updated(delta, 0.19, 0.76)
    assert abs(updated.reservoir.pressure[0] - (213.5 + 4.0)) < 1e-9
    assert abs(updated.wells.bhp_of("PROD-1")
              - (wells.bhp_of("PROD-1") + 9.0)) < 1e-9


def test_bhp_change_limit_is_applied():
    """Quyu təzyiqi bir iterasiyada həddindən çox sıçrasa, perforasiya
    debitləri qeyri-real dəyərlərə gedər."""
    model, connections, reservoir, wells = _setup()
    state = CoupledState(reservoir, wells)
    delta = np.zeros(state.size)
    delta[state.well_index("PROD-1")] = 500.0

    updated = state.updated(delta, 0.19, 0.76, max_bhp_change=15.0)
    change = updated.wells.bhp_of("PROD-1") - wells.bhp_of("PROD-1")
    assert abs(change - 15.0) < 1e-9


def test_reservoir_limits_still_apply_unchanged():
    """Mövcud Appleyard kəsməsi dəyişməməlidir."""
    model, connections, reservoir, wells = _setup()
    state = CoupledState(reservoir, wells)
    delta = np.zeros(state.size)
    delta[0] = 500.0
    updated = state.updated(delta, 0.19, 0.76, max_pressure_change=20.0)
    assert abs(updated.reservoir.pressure[0] - (213.5 + 20.0)) < 1e-9


def test_update_without_limits_applies_the_full_step():
    model, connections, reservoir, wells = _setup()
    state = CoupledState(reservoir, wells)
    delta = np.zeros(state.size)
    delta[state.well_index("INJ-1")] = 33.0
    updated = state.updated(delta, 0.19, 0.76)
    assert abs(updated.wells.bhp_of("INJ-1")
              - (wells.bhp_of("INJ-1") + 33.0)) < 1e-9


def test_copy_of_coupled_state_is_independent():
    model, connections, reservoir, wells = _setup()
    state = CoupledState(reservoir, wells)
    duplicate = state.copy()
    duplicate.wells.bhp[0] += 10.0
    duplicate.reservoir.pressure[0] += 10.0
    assert state.wells.bhp[0] != duplicate.wells.bhp[0]
    assert state.reservoir.pressure[0] != duplicate.reservoir.pressure[0]
