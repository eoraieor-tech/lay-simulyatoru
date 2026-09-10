"""B4b — ilkin həll olmuş qaz (Rs).

PROBLEM (bax `ISH_HESABATI.md` → Seans 6). Üç fazalı mühərrik qaz
papağı (GOC) verilmədikdə nefti "ölü" başladırdı: Rs = 0. Nəticə
zəncirvari idi — OGIP = 0, təzyiq doyma təzyiqindən aşağı düşsə belə
qaz AYRILMIRDI (ayrılacaq həll olmuş qaz yox idi), GOR = 0. Yəni üç
fazalı mühərrik faktiki olaraq iki fazalı işləyirdi.

Bu, v69-un sildiyi bir şey DEYİL — heç vaxt olmayıb. A7-nin açıq
sənədləşdirilmiş mühafizəkar defoltu idi.

HƏLL: `InitialConditions.solution_gor`. `None` (defolt) olanda dəyər
PVT cədvəlindən çıxarılır: `Rs = Rs_sat(min(P, Pb))` — sənaye
standartı (Eclipse `EQUIL`/`RSVD` ilə eyni məntiq).
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.application.serialization import ProjectSerializer
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.initial import InitialConditions
from imex2d.domain.scal import GasCoreyParameters
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter


def _model(bubble_point=200.0, initial=None, nx=8, ny=8):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=20.0, dy=20.0, dz=10.0,
        porosity=0.22, permx_base=150.0)
    table = build_pvt_table(pressure_min=1.0, pressure_max=400.0, n_points=40,
                            bubble_point_bar=bubble_point, include_gas=True)
    return ReservoirModelBuilder().build(
        geological_model=geology, wells=five_spot(geology.grid),
        scal=default_scal(), gas_scal=GasCoreyParameters(),
        pvt_table=table, initial=initial, name="B4b sınağı")


def _service():
    return ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(),
        engine_factory=FullyImplicitEngine)


def _engine(model, end_time=400.0):
    return _service().create_engine(model, SimulationConfig(end_time=end_time))


# ══════════════════════════ ilkin vəziyyət ══════════════════════════════
def test_initial_rs_is_no_longer_zero():
    """ƏSAS TEST — əvvəl burada dəqiq 0 var idi ("ölü neft")."""
    engine = _engine(_model())
    rs = engine.state.third_variable
    assert np.all(rs > 0.0), "neft həll olmuş qazla başlamalıdır"


def test_initial_rs_matches_the_pvt_table_at_cell_pressure():
    """Rs = Rs_sat(min(P, Pb)) — sənaye standartı."""
    model = _model(bubble_point=200.0)
    engine = _engine(model)
    pressure = engine.state.pressure
    expected = engine.pvt.solution_gor(
        np.minimum(pressure, model.pvt_table.bubble_point))
    assert np.allclose(engine.state.third_variable, expected)


def test_cells_start_undersaturated_without_a_gas_cap():
    """Həll olmuş qaz var, SƏRBƏST qaz yoxdur — Sg = 0."""
    engine = _engine(_model())
    assert not engine.state.is_saturated.any()
    assert np.allclose(engine.state.gas_saturation, 0.0)


def test_explicit_override_is_used_as_given():
    initial = InitialConditions(datum_pressure=250.0, water_saturation=0.2,
                                solution_gor=42.0)
    engine = _engine(_model(initial=initial))
    assert np.allclose(engine.state.third_variable, 42.0)


def test_higher_bubble_point_means_more_dissolved_gas():
    values = [float(np.mean(_engine(_model(pb)).state.third_variable))
              for pb in (100.0, 200.0, 300.0)]
    assert values == sorted(values), values


# ══════════════════════════ ehtiyat və hasilat ══════════════════════════
def test_original_gas_in_place_is_positive():
    """Əvvəl OGIP = 0 idi — yerdə qaz "yox" idi."""
    engine = _engine(_model())
    assert engine.original_gas_in_place() > 0.0


def test_gas_is_liberated_when_pressure_drops_below_the_bubble_point():
    """M4-ün ƏSL TƏLƏBİ: P < Psat → qaz ayrılır, GOR qalxır.

    Pb = 200 bar, istismarçı BHP = 150 bar → quyu ətrafı hüceyrələr
    doyma təzyiqindən aşağı düşür və sərbəst qaz yaranır.
    """
    result = _service().run(_model(bubble_point=200.0),
                            SimulationConfig(end_time=400.0))
    assert result.converged, result.message

    peak_gas = max(float(np.max(s.gas_saturation)) for s in result.snapshots
                   if s.gas_saturation is not None)
    assert peak_gas > 1e-3, "sərbəst qaz yaranmalıdır"

    initial_rs = float(np.mean(_engine(_model(200.0)).state.third_variable))
    assert result.series.gas_oil_ratio[-1] > initial_rs, \
        "ayrılan qaz GOR-u həll olmuş qaz nisbətindən yuxarı qaldırmalıdır"


def test_no_free_gas_when_pressure_stays_above_the_bubble_point():
    """Fiziki nəzarət: Pb = 100 bar, ən aşağı təzyiq ~150 bar →
    heç bir hüceyrə doyma təzyiqinə çatmır, qaz AYRILMAMALIDIR."""
    result = _service().run(_model(bubble_point=100.0),
                            SimulationConfig(end_time=400.0))
    assert result.converged, result.message
    peak_gas = max(float(np.max(s.gas_saturation)) for s in result.snapshots
                   if s.gas_saturation is not None)
    assert peak_gas < 1e-6


# ══════════════════════════ `.imx` faylı ════════════════════════════════
def test_solution_gor_survives_the_imx_round_trip(tmp_path):
    from imex2d.application.project import Project

    project = Project("B4b")
    model = _model(initial=InitialConditions(datum_pressure=250.0,
                                             solution_gor=37.5))
    project.add_reservoir_model(model)
    path = str(tmp_path / "b4b.imx")
    ProjectSerializer().save(project, path)

    reopened = ProjectSerializer().load(path)
    restored = list(reopened.reservoir_models.values())[0]
    assert restored.initial_conditions.solution_gor == pytest.approx(37.5)


def test_old_imx_without_the_key_derives_from_pvt():
    """GERİYƏ UYĞUNLUQ — B4b-dən əvvəlki fayllarda açar yoxdur."""
    serializer = ProjectSerializer()
    data = {"datum_depth": 0.0, "datum_pressure": 250.0,
            "water_saturation": 0.2, "use_equilibration": False}
    from imex2d.application.serialization import _known_fields
    restored = InitialConditions(**_known_fields(InitialConditions, data))
    assert restored.solution_gor is None      # → PVT-dən çıxarılır


# ══════════════════════════ UI bağlantısı ═══════════════════════════════
@pytest.fixture(scope="module")
def qt_app():
    QtWidgets = pytest.importorskip("PyQt5.QtWidgets")
    yield QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def test_pvt_panel_defaults_to_deriving_from_the_table(qt_app):
    from imex2d.ui.panels import PvtPanel

    panel = PvtPanel()
    assert panel.initial_solution_gor() is None
    assert not panel.solution_gor.isEnabled()


def test_pvt_panel_manual_value_is_passed_through(qt_app):
    from imex2d.ui.panels import PvtPanel

    panel = PvtPanel()
    panel.manual_rs.setChecked(True)
    panel.solution_gor.setValue(75.0)
    assert panel.solution_gor.isEnabled()
    assert panel.initial_solution_gor() == pytest.approx(75.0)


def test_numerical_panel_forwards_the_value_into_initial_conditions(qt_app):
    from imex2d.ui.panels import NumericalPanel

    panel = NumericalPanel()
    assert panel.initial_conditions().solution_gor is None

    panel.set_solution_gor(64.0)
    assert panel.initial_conditions().solution_gor == pytest.approx(64.0)

