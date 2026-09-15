"""B7 addım 3 — SƏTH debiti hədəfi (Eclipse ORAT / qaz vurucusunda RATE).

SPE1 istismarçıya neftin SƏTH debitini (20 000 STB/gün), vurucuya qazın SƏTH
debitini (100 MMscf/gün) verir. Mühərrikin RATE hədəfi isə LAY həcmidir.
DİZAYN (bax `simulation/well_constraints.py::SurfaceRateController`): qalıq
və Jakobian toxunulmur — səth hədəfi addımdan əvvəl yığılmış vəziyyətin
əmsalları ilə lay həcminə çevrilir, addımdan sonra əldə olunan səth debitinə
görə düzəldilir.
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers import default_scal
from test_eclipse_well_controls import _record
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder
from imex2d.application.serialization import ProjectSerializer
from imex2d.application.simulation_service import (ModelAwareSimulationService,
                                                   ModelValidationError)
from imex2d.domain.scal import GasCoreyParameters
from imex2d.domain.wells import (ControlMode, Perforation, Phase, RateBasis,
                                 Well, WellControl, WellType)
from imex2d.io.eclipse_export import EclipseDeckWriter
from imex2d.simulation.impes_engine import ImpesEngine
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.well_constraints import (SURFACE_RATE_TOLERANCE,
                                                 BhpLimitController,
                                                 SurfaceRateController)
from imex2d.simulation.well_model import PeacemanWellModel, WellConnection

OIL = 30.0          # sm³/gün
GAS = 20000.0       # sm³/gün


def _model(include_gas=False, injector=True, injected=Phase.WATER,
           injector_mode=ControlMode.BHP, injector_target=320.0,
           injector_basis=RateBasis.RESERVOIR, producer_target=OIL,
           producer_basis=RateBasis.SURFACE, producer_mode=ControlMode.RATE,
           limit=None, nz=1):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=7, ny=7, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=nz, top_depth=2000.0)
    wells = [Well("PROD", WellType.PRODUCER,
                  WellControl(producer_mode, producer_target, bhp_limit=limit,
                              rate_basis=producer_basis),
                  [Perforation(6, 6, k) for k in range(nz)])]
    if injector:
        wells.insert(0, Well("INJ", WellType.INJECTOR,
                             WellControl(injector_mode, injector_target,
                                         injected_phase=injected,
                                         rate_basis=injector_basis),
                             [Perforation(0, 0, k) for k in range(nz)]))
    return ReservoirModelBuilder().build(
        geology, wells, scal=default_scal(),
        gas_scal=GasCoreyParameters() if include_gas else None,
        pvt_table=build_pvt_table(pressure_min=1.0, pressure_max=400.0,
                                  n_points=40, bubble_point_bar=140.0,
                                  include_gas=include_gas),
        name="səth debiti sınağı")


def _service(engine=FullyImplicitEngine):
    return ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=engine)


def _relative_error(values, target):
    values = np.asarray(values, float)
    return float(np.max(np.abs(values - target) / target))


# ═══════════════════════ domain ═════════════════════════════════════

def test_rate_basis_defaults_to_reservoir_volume():
    assert WellControl(ControlMode.RATE, 50.0).rate_basis is RateBasis.RESERVOIR


def test_rate_basis_survives_the_imx_round_trip():
    well = Well("P", control=WellControl(ControlMode.RATE, 80.0,
                                         rate_basis=RateBasis.SURFACE))
    data = ProjectSerializer._well_to_dict(well)
    assert ProjectSerializer._well_from_dict(data).control.rate_basis is RateBasis.SURFACE
    del data["control"]["rate_basis"]          # B7 addım 3-dən əvvəlki fayl
    assert ProjectSerializer._well_from_dict(data).control.rate_basis is RateBasis.RESERVOIR


def test_connection_carries_the_rate_basis():
    connections = PeacemanWellModel().build_connections(_model())
    basis = {c.well_name: c.rate_basis for c in connections}
    assert basis == {"PROD": RateBasis.SURFACE, "INJ": RateBasis.RESERVOIR}


def test_reservoir_basis_leaves_the_controller_inactive():
    engine = _service().create_engine(_model(producer_basis=RateBasis.RESERVOIR),
                                      SimulationConfig(end_time=5.0))
    assert not engine.surface_rate.active


# ═══════════════════════ çevirmə ════════════════════════════════════

def test_prediction_delivers_the_surface_oil_rate_exactly_at_the_start_state():
    engine = _service().create_engine(_model(nz=2), SimulationConfig(end_time=5.0))
    engine._update_rate_shares()
    engine.surface_rate.predict(engine.state)
    assembler = engine.residual_assembler
    rates = assembler.well_rates(engine.state, assembler.fluid_state(engine.state))
    assert -rates.per_well_oil["PROD"] == pytest.approx(OIL, rel=1e-12)
    assert engine.surface_rate.reservoir_target["PROD"] > OIL, "Bo > 1 → lay həcmi böyükdür"


def test_prediction_is_exact_for_three_phase_oil_and_injected_gas():
    # İki RATE quyusu: təzyiq idarəsi qaydası üçün istismarçıda BHP limiti
    model = _model(include_gas=True, injected=Phase.GAS,
                   injector_mode=ControlMode.RATE, injector_target=GAS,
                   injector_basis=RateBasis.SURFACE, limit=50.0)
    engine = _service().create_engine(model, SimulationConfig(end_time=5.0))
    engine.surface_rate.predict(engine.state)
    well_model = engine.newton.well_model
    rates = well_model.well_rates(engine.state, engine.newton.build_fluid(engine.state))
    assert -rates.per_well_oil["PROD"] == pytest.approx(OIL, rel=1e-12)
    assert rates.per_well_gas["INJ"] == pytest.approx(GAS, rel=1e-12)


def test_limit_restores_the_converted_target_not_the_surface_number():
    """BHP limitindən RATE-ə qayıdan quyu SƏTH rəqəmini lay hədəfi kimi yazmamalıdır."""
    connection = WellConnection("W", 0, 2.0, False, ControlMode.RATE, OIL,
                                bhp_limit=100.0, rate_basis=RateBasis.SURFACE)
    surface = SurfaceRateController([connection], lambda state: [0.5])
    surface.predict(None)
    assert connection.target == pytest.approx(OIL / 0.5)

    limit = BhpLimitController([connection], lambda state: [1.0], surface.rate_target)
    state = type("S", (), {"pressure": np.array([0.0])})()
    limit.begin_step()
    state.pressure[0] = 90.0                      # tələb olunan BHP < 100 → limit
    assert limit.update(state) == ["W"]
    assert connection.mode is ControlMode.BHP
    limit.begin_step()
    state.pressure[0] = 500.0                     # bol təzyiq → RATE-ə qayıdış
    assert limit.update(state) == ["W"]
    assert connection.target == pytest.approx(OIL / 0.5)


# ═══════════════════════ uc-uca simulyasiya ═════════════════════════

def test_two_phase_run_holds_the_surface_oil_rate():
    engine = _service().create_engine(_model(nz=2), SimulationConfig(end_time=300.0))
    result = engine.run()
    assert result.converged, result.message
    error = _relative_error(result.well_oil_rate["PROD"], OIL)
    assert error <= SURFACE_RATE_TOLERANCE, error


def test_three_phase_run_holds_surface_oil_and_injected_gas():
    model = _model(include_gas=True, injected=Phase.GAS,
                   injector_mode=ControlMode.RATE, injector_target=GAS,
                   injector_basis=RateBasis.SURFACE, limit=50.0)
    engine = _service().create_engine(model, SimulationConfig(end_time=300.0))
    result = engine.run()
    assert result.converged, result.message
    assert _relative_error(result.series.gas_injection_rate, GAS) <= SURFACE_RATE_TOLERANCE
    modes = result.well_control_mode.get("PROD", ["RATE"] * result.steps)
    rates = [q for q, mode in zip(result.well_oil_rate["PROD"], modes) if mode == "RATE"]
    assert rates, "istismarçı ən azı bir addım RATE rejimində olmalıdır"
    assert _relative_error(rates, OIL) <= SURFACE_RATE_TOLERANCE


def test_depleting_surface_rate_producer_switches_to_its_limit():
    """SPE1 davranışı: neft debiti saxlanılır, lay tükənəndə BHP limitinə keçid."""
    model = _model(injector=False, producer_target=120.0, limit=180.0, nz=2)
    engine = _service().create_engine(model, SimulationConfig(end_time=200.0))
    result = engine.run()
    assert result.converged, result.message
    modes = result.well_control_mode["PROD"]
    assert modes[0] == "RATE" and "BHP" in modes
    first = modes.index("BHP")
    assert _relative_error(result.well_oil_rate["PROD"][:first], 120.0) <= SURFACE_RATE_TOLERANCE
    assert min(result.well_bhp["PROD"]) >= 180.0 - 1e-9


# ═══════════════════════ qoruyucular ════════════════════════════════

def test_impes_rejects_the_surface_basis_in_user_language():
    with pytest.raises(ModelValidationError, match="Səth debiti"):
        _service(ImpesEngine).create_engine(_model(), SimulationConfig(end_time=10.0))


def test_surface_basis_outside_rate_mode_is_reported():
    warnings = [d.message for d in
                _model(producer_mode=ControlMode.BHP, producer_target=150.0)
                .diagnose().warnings]
    assert any("səth" in m.lower() and "RATE" in m for m in warnings), warnings


def test_eclipse_export_writes_orat_and_surface_injection_rate():
    model = _model(injector_mode=ControlMode.RATE, injector_target=45.0,
                   injector_basis=RateBasis.SURFACE, limit=90.0)
    deck = EclipseDeckWriter().render(model)
    producer = _record(deck, "WCONPROD", "PROD")
    assert producer[2] == "ORAT" and float(producer[3]) == pytest.approx(OIL)
    assert producer[4:8] == ["*"] * 4 and float(producer[8]) == pytest.approx(90.0)
    injector = _record(deck, "WCONINJE", "INJ")
    assert injector[3] == "RATE" and float(injector[4]) == pytest.approx(45.0)
    assert injector[5] == "*"


def test_well_panel_offers_the_rate_basis_column():
    pytest.importorskip("PyQt5.QtWidgets")
    import inspect
    from imex2d.ui import panels
    source = inspect.getsource(panels.WellPanel)
    assert "COL_RATE_BASIS" in source and "rate_basis=" in source
