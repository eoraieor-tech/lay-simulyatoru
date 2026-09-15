"""B7 addım 2 — RATE rejimində BHP limiti və rejim keçidi.

SPE1-in istismarçısı debitlə idarə olunur, lakin BHP minimal həddən aşağı
düşə bilməz; lay tükəndikcə quyu BHP idarəsinə keçir. DİZAYN (bax
`simulation/well_constraints.py::BhpLimitController`): qalıq və Jakobian
toxunulmur — nəzarətçi yığılmış addımdan sonra bağlantının `mode`/`target`-ini
dəyişir və addım yenidən həll olunur. Histerezis və "addımda bir keçid"
qaydası rejim rəqsinin qarşısını alır.
"""

from __future__ import annotations

import math
from types import SimpleNamespace

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
from imex2d.domain.wells import ControlMode, Perforation, Well, WellControl, WellType
from imex2d.io.eclipse_export import EclipseDeckWriter
from imex2d.simulation.impes_engine import ImpesEngine
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.well_constraints import (RATE_RESTORE_MARGIN_BAR,
                                                 BhpLimitController)
from imex2d.simulation.well_model import PeacemanWellModel, WellConnection

RATE = 150.0
LIMIT = 180.0


def _model(limit=LIMIT, injector=False, include_gas=False, nz=2,
           producer_mode=ControlMode.RATE):
    """8×8×2, istismarçı RATE + BHP limiti. Vurucu olmayanda lay TÜKƏNİR
    (ölçüldü: iki fazalıda quyu t ≈ 21 gündə limitə keçir)."""
    geology = SyntheticGeologicalModelBuilder().build(
        nx=8, ny=8, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=nz, top_depth=2000.0)
    wells = [Well("PROD", WellType.PRODUCER,
                  WellControl(producer_mode, RATE, bhp_limit=limit),
                  [Perforation(7, 7, k) for k in range(nz)])]
    if injector:
        wells.insert(0, Well("INJ", WellType.INJECTOR,
                             WellControl(ControlMode.BHP, 255.0),
                             [Perforation(0, 0, k) for k in range(nz)]))
    return ReservoirModelBuilder().build(
        geology, wells, scal=default_scal(),
        gas_scal=GasCoreyParameters() if include_gas else None,
        pvt_table=build_pvt_table(pressure_min=1.0, pressure_max=400.0,
                                  n_points=40, bubble_point_bar=140.0,
                                  include_gas=include_gas),
        name="BHP limiti sınağı")


def _service(engine=FullyImplicitEngine):
    return ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=engine)


def _controller(pressures, mobility=(1.0, 1.0), well_index=(2.0, 2.0),
                injector=False, limit=LIMIT, target=RATE):
    """Saxta bağlantılarla nəzarətçi; `pressures` hüceyrə təzyiqləridir."""
    connections = [WellConnection("W", cell, wi, injector, ControlMode.RATE,
                                  target, bhp_limit=limit)
                   for cell, wi in enumerate(well_index)]
    controller = BhpLimitController(connections, lambda state: list(mobility))
    state = SimpleNamespace(pressure=np.asarray(pressures, float))
    return controller, connections, state


def _state_with_implied_bhp(value, mobility=(1.0, 1.0), well_index=(2.0, 2.0),
                            injector=False):
    """Hüceyrə təzyiqini elə seçir ki, tələb olunan BHP = `value` olsun."""
    weight = sum(w * m for w, m in zip(well_index, mobility))
    pressure = value - RATE / weight if injector else value + RATE / weight
    return SimpleNamespace(pressure=np.full(len(well_index), pressure))


# ═══════════════════════ domain ═════════════════════════════════════

def test_bhp_limit_defaults_to_none():
    assert WellControl(ControlMode.RATE, 50.0).bhp_limit is None


def test_non_positive_bhp_limit_is_a_validation_error():
    assert WellControl(ControlMode.RATE, 50.0, bhp_limit=-10.0).validate()
    assert not WellControl(ControlMode.RATE, 50.0, bhp_limit=120.0).validate()


def test_bhp_limit_survives_the_imx_round_trip():
    well = Well("P", control=WellControl(ControlMode.RATE, 80.0, bhp_limit=95.5))
    data = ProjectSerializer._well_to_dict(well)
    assert ProjectSerializer._well_from_dict(data).control.bhp_limit == pytest.approx(95.5)
    del data["control"]["bhp_limit"]           # B7-dən əvvəlki fayl
    assert ProjectSerializer._well_from_dict(data).control.bhp_limit is None


def test_connection_carries_the_limit_only_in_rate_mode():
    rate = PeacemanWellModel().build_connections(_model())
    assert all(c.bhp_limit == LIMIT for c in rate)
    bhp = PeacemanWellModel().build_connections(
        _model(producer_mode=ControlMode.BHP))
    assert all(c.bhp_limit is None for c in bhp)


# ═══════════════════════ nəzarətçi qaydaları ═════════════════════════

def test_implied_bhp_inverts_the_peaceman_relation():
    controller, _, state = _controller([200.0, 220.0], mobility=(1.0, 3.0))
    weight = 2.0 * 1.0 + 2.0 * 3.0
    expected = (2.0 * 200.0 + 6.0 * 220.0 - RATE) / weight
    assert controller.implied_bhp("W", state.pressure, [1.0, 3.0]) == pytest.approx(expected)

    injector, _, _ = _controller([200.0, 220.0], mobility=(1.0, 3.0), injector=True)
    assert injector.implied_bhp("W", state.pressure, [1.0, 3.0]) == pytest.approx(
        (2.0 * 200.0 + 6.0 * 220.0 + RATE) / weight)


def test_rate_producer_switches_to_its_limit_when_violated():
    controller, connections, _ = _controller([0.0, 0.0])
    controller.begin_step()
    assert controller.update(_state_with_implied_bhp(LIMIT + 5.0)) == []
    assert all(c.mode is ControlMode.RATE for c in connections)

    assert controller.update(_state_with_implied_bhp(LIMIT - 5.0)) == ["W"]
    assert controller.limited["W"]
    assert all(c.mode is ControlMode.BHP and c.target == LIMIT for c in connections)


def test_limited_well_returns_to_rate_only_beyond_the_margin():
    """HİSTEREZİS: sərhəddə qayıtmaq quyunu dərhal yenidən limitə salardı."""
    controller, connections, _ = _controller([0.0, 0.0])
    controller.begin_step()
    controller.update(_state_with_implied_bhp(LIMIT - 5.0))
    for implied in (LIMIT + 0.5, LIMIT + RATE_RESTORE_MARGIN_BAR - 0.1):
        controller.begin_step()
        assert controller.update(_state_with_implied_bhp(implied)) == []
        assert controller.limited["W"], implied

    controller.begin_step()
    assert controller.update(
        _state_with_implied_bhp(LIMIT + RATE_RESTORE_MARGIN_BAR + 0.5)) == ["W"]
    assert all(c.mode is ControlMode.RATE and c.target == RATE for c in connections)


def test_well_switches_at_most_once_per_step():
    """Addımın təkrar həlləri arasında aç-qapa rəqsi mümkün deyil."""
    controller, _, _ = _controller([0.0, 0.0])
    controller.begin_step()
    assert controller.update(_state_with_implied_bhp(LIMIT - 5.0)) == ["W"]
    assert controller.update(_state_with_implied_bhp(LIMIT + 50.0)) == []
    assert controller.limited["W"]
    assert controller.switches["W"] == 1


def test_well_that_cannot_flow_counts_as_a_violation():
    controller, connections, state = _controller([250.0, 250.0], mobility=(0.0, 0.0))
    controller.begin_step()
    assert controller.update(state) == ["W"]
    assert all(c.mode is ControlMode.BHP for c in connections)


def test_injector_limit_is_a_maximum():
    controller, connections, _ = _controller([0.0, 0.0], injector=True, limit=300.0)
    controller.begin_step()
    assert controller.update(_state_with_implied_bhp(290.0, injector=True)) == []
    assert controller.update(_state_with_implied_bhp(310.0, injector=True)) == ["W"]
    assert all(c.mode is ControlMode.BHP and c.target == 300.0 for c in connections)


def test_controller_is_inactive_without_limits():
    connections = PeacemanWellModel().build_connections(_model(limit=None, injector=True))
    assert not BhpLimitController(connections, lambda state: []).active


def test_switch_keeps_the_total_rate_continuous():
    """Keçid anında debit SIÇRAMIR: tələb olunan BHP-də BHP rejimi eyni
    cəmi debiti verir — tərs düstur qalığın öz mobilliyi ilə qurulub."""
    model = _model()
    engine = _service().create_engine(model, SimulationConfig(end_time=5.0))
    engine._update_rate_shares()
    assembler = engine.residual_assembler
    mobility = engine._connection_mobilities(engine.state)
    implied = engine.bhp_limit.implied_bhp("PROD", engine.state.pressure, mobility)

    def produced():
        fluid = assembler.fluid_state(engine.state)
        rates = assembler.well_rates(engine.state, fluid)
        return -sum(rates.water[c.cell] * fluid.bw[c.cell]
                    + rates.oil[c.cell] * fluid.bo[c.cell] for c in assembler.wells)

    assert produced() == pytest.approx(RATE, rel=1e-12)
    for connection in assembler.wells:
        connection.mode, connection.target = ControlMode.BHP, implied
    assert produced() == pytest.approx(RATE, rel=1e-9)


# ═══════════════════════ uc-uca simulyasiya ═════════════════════════

@pytest.mark.parametrize("include_gas", [False, True], ids=["iki-fazali", "uc-fazali"])
def test_depleting_rate_producer_switches_to_its_bhp_limit(include_gas):
    """ƏSAS TƏLƏB: lay tükənəndə quyu limitə keçir, BHP limitdən AŞAĞI düşmür,
    rejim bir dəfə dəyişir (rəqs yoxdur)."""
    model = _model(include_gas=include_gas)
    engine = _service().create_engine(model, SimulationConfig(end_time=200.0))
    result = engine.run()
    assert result.converged, result.message

    modes = result.well_control_mode["PROD"]
    bhp = np.asarray(result.well_bhp["PROD"])
    liquid = (np.asarray(result.well_oil_rate["PROD"])
              + np.asarray(result.well_water_rate["PROD"]))
    assert len(modes) == len(bhp) == result.steps
    assert modes[0] == "RATE" and "BHP" in modes
    assert sum(a != b for a, b in zip(modes, modes[1:])) == 1, "rejim rəqs edir"
    assert engine.bhp_limit.switches["PROD"] == 1
    assert engine.bhp_limit_resolves >= 1

    finite = bhp[np.isfinite(bhp)]
    assert finite.min() >= LIMIT - 1e-9, "BHP limitdən aşağı düşdü"
    first = modes.index("BHP")
    assert np.all(bhp[first:] == LIMIT)
    assert liquid[-1] < liquid[first - 1], "limitdə debit azalmalıdır"


def test_limit_that_is_never_reached_changes_nothing():
    """Vurucu təzyiqi saxlayanda limit toxunulmur — nəticə limitsizlə EYNİDİR."""
    with_limit = _service().run(_model(injector=True), SimulationConfig(end_time=200.0))
    without = _service().run(_model(limit=None, injector=True),
                             SimulationConfig(end_time=200.0))
    assert with_limit.converged and without.converged
    assert set(with_limit.well_control_mode["PROD"]) == {"RATE"}
    assert np.array_equal(with_limit.series.cumulative_oil, without.series.cumulative_oil)
    assert np.array_equal(with_limit.series.average_pressure,
                          without.series.average_pressure)


# ═══════════════════════ qoruyucular və diaqnostika ══════════════════

def test_impes_rejects_the_limit_in_user_language():
    with pytest.raises(ModelValidationError, match="BHP limiti"):
        _service(ImpesEngine).create_engine(_model(), SimulationConfig(end_time=10.0))


def test_rate_well_with_a_limit_satisfies_the_pressure_control_rule():
    """SPE1 kimi yalnız debitlə, lakin BHP həddi ilə idarə olunan model bloklanmır."""
    messages = [d.message for d in _model().diagnose().errors]
    assert not any("idarə olunmalıdır" in m for m in messages), messages
    messages = [d.message for d in _model(limit=None).diagnose().errors]
    assert any("idarə olunmalıdır" in m for m in messages), messages


def test_limit_outside_rate_mode_is_reported():
    """SƏSSİZ ATILMIR: BHP rejimində limit işləmir — xəbərdarlıq."""
    warnings = [d.message for d in
                _model(producer_mode=ControlMode.BHP).diagnose().warnings]
    assert any("yalnız RATE" in m for m in warnings), warnings


def test_producer_limit_above_reservoir_pressure_is_reported():
    warnings = [d.message for d in _model(limit=400.0).diagnose().warnings]
    assert any("minimal BHP" in m for m in warnings), warnings


def test_eclipse_export_writes_the_limit_into_the_bhp_column():
    geology_model = _model(injector=True)
    geology_model.wells[0].control = WellControl(ControlMode.RATE, 90.0,
                                                 bhp_limit=320.0)
    deck = EclipseDeckWriter().render(geology_model)
    assert float(_record(deck, "WCONPROD", "PROD")[8]) == pytest.approx(LIMIT)
    assert float(_record(deck, "WCONINJE", "INJ")[6]) == pytest.approx(320.0)


def test_well_panel_offers_the_bhp_limit_column():
    pytest.importorskip("PyQt5.QtWidgets")
    import inspect
    from imex2d.ui import panels
    source = inspect.getsource(panels.WellPanel)
    assert '"BHP limiti"' in source and "COL_BHP_LIMIT" in source
    assert "bhp_limit=bhp_limit" in source
