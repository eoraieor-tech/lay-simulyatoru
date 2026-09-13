"""B4-B — quyunu quyu başı təzyiqi (THP) ilə idarə etmək.

DİZAYN (bax `simulation/wellbore/thp_control.py`): THP quyusu mühərrik
üçün ADİ BHP bağlantısıdır; `ThpController` hər addımdan sonra son
debitlərlə traversi tərsinə həll edib BHP-ni yeniləyir (açıq birləşmə).

TƏHLÜKƏLİ BOŞLUQLAR (tətbiqdən əvvəl tapıldı, testlə kilidləndi):
  * `connection.mode is ControlMode.BHP` yoxlaması qalıqda, Jakobianda,
    IMPES-də var — yeni rejim bağlantıya ötürülsəydi SƏSSİZCƏ RATE kimi
    işləyərdi;
  * Eclipse ixracı BHP olmayan istismarçını LRAT kimi yazır — THP ədədi
    debit kimi deck-ə düşərdi.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.application.serialization import ProjectSerializer
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.properties import FluidProperties
from imex2d.domain.scal import GasCoreyParameters
from imex2d.domain.tubing import TubingGeometry
from imex2d.domain.wells import ControlMode, Well, WellControl, WellType
from imex2d.simulation.impes_engine import ImpesEngine
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.well_model import PeacemanWellModel
from imex2d.simulation.wellbore import WellStream, pressure_traverse
from imex2d.simulation.wellbore.thp_control import (MAX_BHP_CHANGE_BAR,
                                                     ThpController,
                                                     bhp_from_thp)

GRAVITY = 9.80665


def _fluids():
    return FluidProperties(water_fvf=1.0, oil_fvf=1.0)


def _model(thp=20.0, with_gas=False, tubing=True, top_depth=1200.0,
           producer_mode=ControlMode.THP):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=8, ny=8, dx=20.0, dy=20.0, dz=10.0, porosity=0.22,
        permx_base=150.0, top_depth=top_depth)
    wells = five_spot(geology.grid)
    for well in wells:
        if well.well_type is WellType.PRODUCER:
            well.control = WellControl(producer_mode, thp)
            if tubing:
                well.tubing = TubingGeometry(diameter=0.062, segments=20)
    return ReservoirModelBuilder().build(
        geological_model=geology, wells=wells, scal=default_scal(),
        gas_scal=GasCoreyParameters() if with_gas else None,
        pvt_table=build_pvt_table(pressure_min=1.0, pressure_max=400.0,
                                  n_points=40, bubble_point_bar=240.0,
                                  include_gas=with_gas),
        name="THP idarəsi sınağı")


def _service(engine=FullyImplicitEngine):
    return ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=engine)


def _producer(model):
    return next(w for w in model.wells if w.well_type is WellType.PRODUCER)


# ═══════════════════════ domain ═════════════════════════════════════

def test_thp_is_a_control_mode():
    assert ControlMode("THP") is ControlMode.THP


def test_thp_target_is_validated_as_a_pressure():
    assert WellControl(ControlMode.THP, -5.0).validate()
    assert not WellControl(ControlMode.THP, 20.0).validate()


def test_thp_survives_the_imx_round_trip():
    well = Well("P1", control=WellControl(ControlMode.THP, 18.5),
                tubing=TubingGeometry())
    restored = ProjectSerializer._well_from_dict(ProjectSerializer._well_to_dict(well))
    assert restored.control.mode is ControlMode.THP
    assert restored.control.target == pytest.approx(18.5)


# ═══════════════════════ tərs traverse ═══════════════════════════════

def test_bhp_from_thp_inverts_the_traverse():
    """ƏSAS RİYAZİ ŞƏRT: tapılan BHP geri traversdə eyni THP-ni verir."""
    tubing = TubingGeometry(diameter=0.062, segments=20)
    stream = WellStream(oil=80.0, water=20.0, water_density=1000.0)
    for thp in (5.0, 20.0, 60.0):
        bhp = bhp_from_thp(thp, 1500.0, tubing, stream, fluids=_fluids())
        back = pressure_traverse(bhp, 1500.0, tubing, stream, fluids=_fluids()).thp
        assert back == pytest.approx(thp, abs=0.02)


def test_higher_thp_needs_a_higher_bhp():
    tubing = TubingGeometry()
    stream = WellStream(oil=80.0, water=20.0)
    values = [bhp_from_thp(thp, 1500.0, tubing, stream, fluids=_fluids())
              for thp in (5.0, 20.0, 40.0, 80.0)]
    assert all(b > a for a, b in zip(values, values[1:])), values


def test_static_column_estimate_is_the_hydrostatic_head():
    """Cüzi debitdə BHP ≈ THP + ρ·g·h (sürtünmə yox)."""
    tubing = TubingGeometry(segments=40)
    stream = WellStream(water=1e-3, water_density=1000.0)
    bhp = bhp_from_thp(10.0, 1000.0, tubing, stream, fluids=_fluids())
    assert bhp == pytest.approx(10.0 + 1000.0 * GRAVITY * 1000.0 / 1e5, rel=1e-3)


def test_bhp_is_never_below_thp():
    stream = WellStream(oil=50.0)
    bhp = bhp_from_thp(30.0, 800.0, TubingGeometry(), stream, fluids=_fluids())
    assert bhp >= 30.0


def test_unreachable_thp_returns_nan():
    """Axtarış həddinə qədər heç bir BHP bu THP-ni vermirsə uydurma dəyər yox."""
    stream = WellStream(water=10.0, water_density=1000.0)
    bhp = bhp_from_thp(500.0, 3000.0, TubingGeometry(), stream,
                       fluids=_fluids(), bhp_limit=300.0)
    assert math.isnan(bhp)


# ═══════════════════════ bağlantı səviyyəsi ══════════════════════════

def test_thp_well_is_built_as_a_bhp_connection():
    """TƏHLÜKƏLİ BOŞLUĞUN QARŞISI: qalıq yalnız BHP tanıyır.

    THP bağlantıya ötürülsəydi `mode is BHP` yoxlamaları onu RATE kimi
    işlədərdi — yəni THP ədədi debit hədəfinə çevrilərdi.
    """
    model = _model(thp=20.0)
    connections = PeacemanWellModel().build_connections(model)
    producer = [c for c in connections if not c.is_injector]
    assert producer
    for connection in producer:
        assert connection.mode is ControlMode.BHP
        assert connection.thp_target == pytest.approx(20.0)
    injector = [c for c in connections if c.is_injector]
    assert all(c.thp_target is None for c in injector)


def test_controller_initialises_bhp_above_thp():
    model = _model(thp=20.0)
    connections = PeacemanWellModel().build_connections(model)
    controller = ThpController(model, connections, fluids=model.fluids)
    assert controller.active
    controller.initialize()
    name = _producer(model).name
    assert controller.bhp[name] > 20.0 + 50.0, "1200 m sütun BHP-ni xeyli qaldırmalıdır"
    for connection in connections:
        if connection.well_name == name:
            assert connection.target == pytest.approx(controller.bhp[name])


def test_controller_limits_the_change_per_step():
    """Açıq birləşmənin rəqsə qarşı qoruyucusu — addım başına maksimal dəyişmə."""
    model = _model(thp=20.0)
    connections = PeacemanWellModel().build_connections(model)
    controller = ThpController(model, connections, fluids=model.fluids,
                               relaxation=1.0)
    controller.initialize()
    name = _producer(model).name
    before = controller.bhp[name]
    controller.update({name: -5000.0}, {name: -5000.0})   # nəhəng debit
    assert abs(controller.bhp[name] - before) <= MAX_BHP_CHANGE_BAR + 1e-9


def test_no_thp_wells_means_an_inactive_controller():
    model = _model(producer_mode=ControlMode.BHP, thp=150.0)
    controller = ThpController(model, PeacemanWellModel().build_connections(model))
    assert not controller.active


# ═══════════════════════ uc-uca: simulyasiya ═════════════════════════

def test_two_phase_run_with_thp_control_converges_and_tracks_thp():
    """ƏSAS TƏLƏB: THP ilə idarə olunan quyu işləyir və hədəfi izləyir."""
    model = _model(thp=20.0)
    result = _service().run(model, SimulationConfig(end_time=400.0))
    assert result.converged, result.message

    name = _producer(model).name
    bhp = result.well_bhp[name]
    thp = result.well_thp[name]
    assert len(bhp) == len(result.well_oil_rate[name])
    assert all(b >= 20.0 for b in bhp), "BHP THP-dən aşağı düşməməlidir"
    assert len(set(np.round(bhp, 6))) > 1, "BHP addım-addım yenilənməlidir"

    finite = [v for v in thp[len(thp) // 2:] if math.isfinite(v)]
    assert finite, "ikinci yarıda THP hesablanmalıdır"
    assert np.median(np.abs(np.asarray(finite) - 20.0)) < 3.0, \
        f"THP hədəfə yaxın qalmalıdır: {finite[-5:]}"


def test_three_phase_run_with_thp_control_converges():
    model = _model(thp=20.0, with_gas=True)
    result = _service().run(model, SimulationConfig(end_time=300.0))
    assert result.converged, result.message
    assert result.well_bhp[_producer(model).name]


def test_higher_thp_target_produces_less_oil():
    """Fiziki monotonluq: quyu başında böyük əks-təzyiq → az hasilat."""
    low = _service().run(_model(thp=10.0), SimulationConfig(end_time=300.0))
    high = _service().run(_model(thp=60.0), SimulationConfig(end_time=300.0))
    assert low.converged and high.converged
    assert high.series.cumulative_oil[-1] < low.series.cumulative_oil[-1]


def test_bhp_controlled_runs_are_unchanged():
    """GERİYƏ UYĞUNLUQ: BHP modellərində nəzarətçi heç nə etmir."""
    model = _model(producer_mode=ControlMode.BHP, thp=150.0)
    engine = _service().create_engine(model, SimulationConfig(end_time=100.0))
    assert not engine.thp_control.active


# ═══════════════════════ imtinalar və diaqnostika ════════════════════

def test_impes_rejects_thp_in_user_language():
    model = _model(thp=20.0)
    with pytest.raises(Exception, match="THP"):
        _service(ImpesEngine).create_engine(model, SimulationConfig(end_time=100.0))


def test_thp_without_tubing_is_a_diagnostic_error():
    model = _model(thp=20.0, tubing=False)
    messages = [d.message for d in model.diagnose().errors]
    assert any("lülə" in m for m in messages), messages


def test_thp_only_producer_satisfies_the_pressure_control_rule():
    """"Ən azı bir BHP quyusu" şərti THP quyusunu da saymalıdır."""
    model = _model(thp=20.0)
    messages = [d.message for d in model.diagnose().errors]
    assert not any("idarə olunmalıdır" in m for m in messages), messages


def test_eclipse_export_refuses_thp_instead_of_writing_a_rate():
    """TƏHLÜKƏLİ BOŞLUĞUN QARŞISI: THP səssizcə LRAT kimi yazılmamalıdır."""
    import inspect
    from imex2d.io import eclipse_export
    source = inspect.getsource(eclipse_export)
    block = source.split("for well in producers:", 1)[1][:600]
    assert "ControlMode.THP" in block and "raise ValueError" in block


def test_ui_offers_thp_mode():
    pytest.importorskip("PyQt5.QtWidgets")
    import inspect
    from imex2d.ui import panels
    assert 'mode_box.addItems(["BHP", "RATE", "THP"])' in inspect.getsource(panels)
