"""B7 (addım 1) — QAZ VURAN quyu.

TAPINTI: `WellControl.injected_phase` domendə əvvəldən var idi və layihə
faylında saxlanılırdı, lakin MÜHƏRRİK ONU OXUMURDU — bütün vurucular su
vururdu. SPE1 etalonu qaz vurur, ona görə bu boşluq bağlandı.

TƏHLÜKƏ (testlə kilidlənib): qaz fazası olmayan modeldə qaz vuran quyu
səssizcə SU vurardı — istifadəçinin seçdiyindən tamamilə fərqli hesab.
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder
from imex2d.application.serialization import ProjectSerializer
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.scal import CoreyParameters, GasCoreyParameters
from imex2d.domain.wells import (ControlMode, Perforation, Phase, Well,
                                 WellControl, WellType)
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.implicit.three_phase_residual import (
    ThreePhaseAccumulator, ThreePhaseFlux, ThreePhaseFluidState,
    ThreePhaseJacobianAssembler, ThreePhaseWellModel)
from imex2d.simulation.implicit.three_phase_state import ThreePhaseState
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.stone_relperm import StoneRelativePermeabilityProvider
from imex2d.simulation.well_model import PeacemanWellModel

WATER_OIL = CoreyParameters()
GAS = GasCoreyParameters()
ENDPOINT_WATER = 0.35


def _wells(injected=Phase.GAS, injector_mode=ControlMode.BHP,
           injector_target=320.0, nx=5, ny=5):
    return [
        Well("INJ", WellType.INJECTOR,
             WellControl(injector_mode, injector_target, injected_phase=injected),
             [Perforation(0, 0, 0)]),
        Well("PROD", WellType.PRODUCER, WellControl(ControlMode.BHP, 150.0),
             [Perforation(nx - 1, ny - 1, 0)]),
    ]


def _model(injected=Phase.GAS, include_gas=True, nx=5, ny=5,
           injector_mode=ControlMode.BHP, injector_target=320.0):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    return ReservoirModelBuilder().build(
        geology, _wells(injected, injector_mode, injector_target, nx, ny),
        scal=WATER_OIL, gas_scal=GAS if include_gas else None,
        pvt_table=build_pvt_table(pressure_min=1.0, pressure_max=400.0,
                                  n_points=40, bubble_point_bar=240.0,
                                  include_gas=include_gas),
        name="Qaz vurulması sınağı")


def _pieces(model):
    grid = TwoPointFluxDiscretization().build(model)
    connections = PeacemanWellModel().build_connections(model)
    relperm = StoneRelativePermeabilityProvider.from_corey(WATER_OIL, GAS)
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    well_model = ThreePhaseWellModel(model, connections, ENDPOINT_WATER,
                                     GAS.krg_end)
    return grid, connections, relperm, pvt, well_model


def _fluid(state, relperm, pvt):
    n = state.ncell
    sw, sg = state.water_saturation, state.gas_saturation
    return ThreePhaseFluidState(
        mu_w=np.full(n, 0.5), mu_o=pvt.oil_viscosity(state.pressure),
        mu_g=pvt.gas_viscosity(state.pressure), bw=np.full(n, 1.0),
        bo=pvt.oil_fvf(state.pressure), bg=pvt.gas_fvf(state.pressure),
        rs=state.solution_gor(pvt), krw=relperm.krw(sw),
        kro=relperm.kro_three_phase(sw, sg), krg=relperm.krg(sg))


def _state(model, pressure=200.0, sw=0.3, sg=0.1):
    n = model.ncell
    return ThreePhaseState(np.full(n, pressure), np.full(n, sw),
                           np.full(n, sg), np.ones(n, bool))


# ═══════════════════════ bağlantı səviyyəsi ══════════════════════════

def test_connection_carries_the_injected_phase():
    connections = PeacemanWellModel().build_connections(_model())
    injector = [c for c in connections if c.is_injector]
    assert injector and all(c.injected_phase is Phase.GAS for c in injector)
    producer = [c for c in connections if not c.is_injector]
    assert all(c.injected_phase is Phase.WATER for c in producer)


def test_water_injection_is_still_the_default():
    connections = PeacemanWellModel().build_connections(_model(Phase.WATER))
    assert all(c.injected_phase is Phase.WATER for c in connections)


# ═══════════════════════ debit hesabı ════════════════════════════════

def test_gas_injector_feeds_the_gas_equation_not_water():
    model = _model(Phase.GAS)
    grid, connections, relperm, pvt, well_model = _pieces(model)
    state = _state(model)
    rates = well_model.well_rates(state, _fluid(state, relperm, pvt))
    cell = next(c.cell for c in connections if c.is_injector)

    assert rates.gas[cell] > 0.0, "vurulan qaz qaz tənliyinə düşməlidir"
    assert rates.per_well_gas["INJ"] > 0.0
    assert rates.per_well_water["INJ"] == 0.0, "qaz vuran quyu su vurmamalıdır"


def test_water_injector_is_unchanged():
    """GERİYƏ UYĞUNLUQ: su vuran quyu əvvəlki kimi işləməlidir."""
    model = _model(Phase.WATER)
    grid, connections, relperm, pvt, well_model = _pieces(model)
    state = _state(model)
    rates = well_model.well_rates(state, _fluid(state, relperm, pvt))
    cell = next(c.cell for c in connections if c.is_injector)
    assert rates.water[cell] > 0.0
    assert rates.per_well_gas["INJ"] == 0.0


def test_gas_injection_rate_matches_the_peaceman_formula():
    model = _model(Phase.GAS)
    grid, connections, relperm, pvt, well_model = _pieces(model)
    state = _state(model)
    fluid = _fluid(state, relperm, pvt)
    rates = well_model.well_rates(state, fluid)

    connection = next(c for c in connections if c.is_injector)
    cell = connection.cell
    expected = (connection.well_index * GAS.krg_end / fluid.mu_g[cell]
                * (connection.target - state.pressure[cell]) / fluid.bg[cell])
    assert rates.gas[cell] == pytest.approx(expected, rel=1e-12)


def test_rate_controlled_gas_injector_uses_the_target_directly():
    model = _model(Phase.GAS, injector_mode=ControlMode.RATE,
                   injector_target=50000.0)
    grid, connections, relperm, pvt, well_model = _pieces(model)
    state = _state(model)
    fluid = _fluid(state, relperm, pvt)
    rates = well_model.well_rates(state, fluid)
    cell = next(c.cell for c in connections if c.is_injector)
    assert rates.gas[cell] == pytest.approx(50000.0 / fluid.bg[cell])


# ═══════════════════════ Jakobian — sonlu fərq ═══════════════════════

def _jacobian_error(model, columns=range(9)):
    grid, connections, relperm, pvt, well_model = _pieces(model)
    accumulator = ThreePhaseAccumulator(model, grid.pore_volume)
    flux = ThreePhaseFlux(model, grid)
    assembler = ThreePhaseJacobianAssembler(model, accumulator, flux,
                                            well_model, relperm, pvt)
    n = model.ncell
    state = _state(model)
    previous = state.copy()

    def residual(s):
        fluid, fluid_prev = _fluid(s, relperm, pvt), _fluid(previous, relperm, pvt)
        n_w, n_o, n_g = accumulator.accumulation(s, fluid)
        n_w0, n_o0, n_g0 = accumulator.accumulation(previous, fluid_prev)
        in_w, in_o, in_g = flux.net_influx(s, fluid)
        rates = well_model.well_rates(s, fluid)
        out = np.empty(n * 3)
        out[0::3] = (n_w - n_w0) - in_w - rates.water
        out[1::3] = (n_o - n_o0) - in_o - rates.oil
        out[2::3] = (n_g - n_g0) - in_g - rates.gas
        return out

    jacobian = assembler.assemble(state, _fluid(state, relperm, pvt), dt=1.0)
    vector = state.to_vector()
    worst = 0.0
    for column in columns:
        step = 1e-6 * max(1.0, abs(vector[column]))
        forward, backward = vector.copy(), vector.copy()
        forward[column] += step
        backward[column] -= step
        numeric = (residual(ThreePhaseState.from_vector(forward, state.is_saturated))
                   - residual(ThreePhaseState.from_vector(backward, state.is_saturated))
                   ) / (2 * step)
        analytic = np.asarray(jacobian[:, column].todense()).ravel()
        worst = max(worst, np.max(np.abs(numeric - analytic))
                    / max(1.0, np.max(np.abs(numeric))))
    return worst


def test_gas_injector_jacobian_matches_finite_difference():
    """Vurulan faza SƏTRİ dəyişdirir (su → 0-cı, qaz → 2-ci tənlik)."""
    assert _jacobian_error(_model(Phase.GAS)) < 1e-5


def test_rate_controlled_gas_injector_jacobian_matches_finite_difference():
    assert _jacobian_error(_model(Phase.GAS, injector_mode=ControlMode.RATE,
                                  injector_target=50000.0)) < 1e-5


def test_water_injector_jacobian_is_unchanged():
    """SU yolu TOXUNULMAYIB — ölçülmüş dəqiqliyi olduğu kimi qalır.

    ÖLÇÜLDÜ: su vurucusunun Jakobian xətası 1.34e-5-dir və bu, B7-dən
    ƏVVƏL də belə idi (vurucu mobilliyinin doyumluluq törəməsi nəzərə
    alınmır — köhnə sadələşdirmə). Qaz yolu isə 1e-5-dən kiçikdir, yəni
    yeni kod köhnəsindən PİS deyil.
    """
    water_error = _jacobian_error(_model(Phase.WATER))
    gas_error = _jacobian_error(_model(Phase.GAS))
    assert water_error < 1e-4, water_error
    assert gas_error <= water_error, (gas_error, water_error)


# ═══════════════════════ uc-uca simulyasiya ══════════════════════════

def _run(model, end_time=200.0):
    service = ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)
    return service.run(model, SimulationConfig(end_time=end_time))


def test_gas_injection_run_converges_and_builds_a_gas_zone():
    model = _model(Phase.GAS)
    result = _run(model)
    assert result.converged, result.message
    last = result.snapshots[-1]
    assert last.gas_saturation is not None
    injector_cell = 0
    assert last.gas_saturation[injector_cell] > 0.1, \
        "vurucu quyunun ətrafında qaz doymuşluğu artmalıdır"


def test_gas_injection_supports_pressure_more_than_no_injection():
    """Fiziki yoxlama: qaz vurulması lay təzyiqini SAXLAYIR."""
    with_gas = _run(_model(Phase.GAS))
    shut = _model(Phase.GAS)
    for well in shut.wells:
        if well.well_type is WellType.INJECTOR:
            well.active = False
    without = _run(shut)
    assert with_gas.converged and without.converged
    assert (with_gas.series.average_pressure[-1]
            > without.series.average_pressure[-1])


# ═══════════════════════ qoruyucular ═════════════════════════════════

def test_gas_injection_without_a_gas_phase_is_an_error():
    """TƏHLÜKƏNİN QARŞISI: iki fazalı mühərrik səssizcə SU vurardı."""
    model = _model(Phase.GAS, include_gas=False)
    messages = [d.message for d in model.diagnose().errors]
    assert any("qaz" in m.lower() for m in messages), messages


def test_water_injection_without_a_gas_phase_is_fine():
    model = _model(Phase.WATER, include_gas=False)
    assert not [d for d in model.diagnose().errors
                if "vur" in d.message.lower()]


def test_injected_phase_survives_the_project_round_trip():
    well = Well("INJ", WellType.INJECTOR,
                WellControl(ControlMode.RATE, 1000.0, injected_phase=Phase.GAS))
    restored = ProjectSerializer._well_from_dict(ProjectSerializer._well_to_dict(well))
    assert restored.control.injected_phase is Phase.GAS


def test_well_panel_offers_the_injected_phase_column():
    pytest.importorskip("PyQt5.QtWidgets")
    import inspect
    from imex2d.ui import panels
    source = inspect.getsource(panels)
    assert "Vurulan faza" in source and "COL_PHASE" in source
