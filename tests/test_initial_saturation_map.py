"""İlkin doyumluluğun geologiya `SW` xəritəsindən qurulması (A1).

`use_saturation_map` bayrağı SÖNDÜRÜLÜ olanda hər şey əvvəlki kimi
qalmalıdır — testlərin yarısı məhz bunu qoruyur.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from helpers import default_scal, make_service, short_config
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder
from imex2d.application.serialization import ProjectSerializer
from imex2d.application.simulation_service import (ModelAwareSimulationService,
                                                   ModelValidationError)
from imex2d.domain.diagnostics import Severity
from imex2d.domain.initial import InitialConditions
from imex2d.domain.properties import PropertyMap
from imex2d.domain.wells import (ControlMode, Perforation, Well, WellControl,
                                 WellType)
from imex2d.history.parameters import (ModelModifier, ParameterSet,
                                       standard_parameters)
from imex2d.simulation.impes_engine import ImpesEngine
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.initialization.equilibrium import (
    EquilibriumInitializationProvider)
from imex2d.simulation.initialization.saturation_map import (
    SaturationMapInitializationProvider, SaturationMapOverride,
    count_outside_scal_limits, water_saturation_from_map)


def _wells(nx: int, nz: int = 1):
    perfs_inj = [Perforation(0, 0, k) for k in range(nz)]
    perfs_prod = [Perforation(nx - 1, 0, k) for k in range(nz)]
    return [
        Well("INJ", WellType.INJECTOR,
             WellControl(ControlMode.BHP, 320.0), perfs_inj),
        Well("PROD", WellType.PRODUCER,
             WellControl(ControlMode.BHP, 150.0), perfs_prod),
    ]


def make_model(nx=8, nz=1, sw_map=None, top_depth=2000.0, **initial):
    """Kiçik model; `sw_map` verilibsə `property_maps["SW"]` kimi qoşulur."""
    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=1, nz=nz, dx=20.0, dy=100.0, dz=10.0,
        porosity=0.22, permx_base=200.0, top_depth=top_depth)
    conditions = InitialConditions(datum_depth=top_depth, datum_pressure=250.0,
                                   water_saturation=0.20, **initial)
    model = ReservoirModelBuilder().build(
        geological_model=geology, wells=_wells(nx, nz), scal=scal,
        initial=conditions, name="SW xəritəsi testi")
    if sw_map is not None:
        values = np.asarray(sw_map, float)
        if values.ndim == 0:
            values = np.full(model.ncell, float(values))
        model.property_maps["SW"] = PropertyMap("SW", values)
    return model


def two_zone_map(ncell: int, low=0.25, high=0.45) -> np.ndarray:
    """Yarısı `low`, yarısı `high` — skalyar ilə qarışdırıla bilməyən xəritə."""
    values = np.full(ncell, float(low))
    values[ncell // 2:] = float(high)
    return values


# ══════════════════════════════════════════ saf funksiya: oxuma və yoxlama

def test_map_is_read_as_is():
    model = make_model(sw_map=two_zone_map(8))
    sw = water_saturation_from_map(model)
    assert np.allclose(sw, two_zone_map(8))


def test_missing_map_raises():
    model = make_model()
    with pytest.raises(ValueError, match="SW"):
        water_saturation_from_map(model)


def test_wrong_size_raises():
    model = make_model(nx=8)
    model.property_maps["SW"] = PropertyMap("SW", np.full(5, 0.3))
    with pytest.raises(ValueError, match="ölçüsü"):
        water_saturation_from_map(model)


def test_nan_in_active_cell_raises():
    values = two_zone_map(8)
    values[3] = np.nan
    model = make_model(sw_map=values)
    with pytest.raises(ValueError, match="NaN"):
        water_saturation_from_map(model)


def test_values_outside_unit_interval_raise():
    values = two_zone_map(8)
    values[2] = 1.4
    model = make_model(sw_map=values)
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        water_saturation_from_map(model)


def test_inactive_cells_may_hold_nan_and_get_the_scalar():
    values = two_zone_map(8)
    values[7] = np.nan
    actnum = np.ones(8, dtype=np.int8)
    actnum[7] = 0
    model = make_model(sw_map=values)
    model.grid = model.grid.with_actnum(actnum)

    sw = water_saturation_from_map(model)
    assert sw[7] == pytest.approx(model.initial_conditions.water_saturation)
    assert np.allclose(sw[:7], values[:7])


def test_clipped_cell_count_is_reported_but_not_applied():
    scal = default_scal()
    values = two_zone_map(8)
    values[0] = scal.swc / 2.0             # Swc-dən kiçik
    values[-1] = 1.0 - scal.sor / 2.0      # 1−Sor-dan böyük
    model = make_model(sw_map=values)

    sw = water_saturation_from_map(model)
    assert sw[0] == pytest.approx(values[0])     # kəsmə BURADA edilmir
    assert sw[-1] == pytest.approx(values[-1])
    assert count_outside_scal_limits(sw, model) == 2


# ═══════════════════════════════════════════════ provider seçimi (4 hal)

def test_provider_selection_follows_the_priority_table():
    values = two_zone_map(8)
    select = ModelAwareSimulationService._initialization

    assert select(make_model(sw_map=values)) is None
    assert isinstance(select(make_model(sw_map=values, use_saturation_map=True)),
                      SaturationMapInitializationProvider)
    assert isinstance(select(make_model(sw_map=values, use_equilibration=True,
                                        oil_water_contact=2020.0)),
                      EquilibriumInitializationProvider)

    both = select(make_model(sw_map=values, use_equilibration=True,
                             use_saturation_map=True, oil_water_contact=2020.0))
    assert isinstance(both, SaturationMapOverride)
    assert isinstance(both.inner, EquilibriumInitializationProvider)


def test_requested_but_unreadable_map_is_a_validation_error():
    model = make_model(use_saturation_map=True)     # xəritə ümumiyyətlə yoxdur
    with pytest.raises(ModelValidationError, match="SW"):
        ModelAwareSimulationService._initialization(model)


def test_engine_construction_fails_when_map_is_broken():
    values = two_zone_map(8)
    values[4] = np.nan
    model = make_model(sw_map=values, use_saturation_map=True)
    with pytest.raises(ModelValidationError, match="NaN"):
        make_service().create_engine(model, short_config())


# ════════════════════════════════════════════ mühərriklər xəritəyə hörmət edir

ENGINES = [ImpesEngine, FullyImplicitEngine]


def engine_saturation(model, engine_factory):
    engine = make_service().with_engine(engine_factory).create_engine(
        model, short_config())
    return np.asarray(engine.sw if engine_factory is ImpesEngine
                      else engine.state.water_saturation, float)


@pytest.mark.parametrize("engine_factory", ENGINES)
def test_engine_uses_the_map_when_the_flag_is_on(engine_factory):
    values = two_zone_map(8)
    model = make_model(sw_map=values, use_saturation_map=True)
    assert np.allclose(engine_saturation(model, engine_factory), values)


@pytest.mark.parametrize("engine_factory", ENGINES)
def test_map_is_ignored_while_the_flag_is_off(engine_factory):
    values = two_zone_map(8)
    model = make_model(sw_map=values)
    scalar = model.initial_conditions.water_saturation
    assert np.allclose(engine_saturation(model, engine_factory), scalar)


@pytest.mark.parametrize("engine_factory", ENGINES)
def test_map_values_are_clipped_to_scal_limits_by_the_engine(engine_factory):
    scal = default_scal()
    values = two_zone_map(8)
    values[0] = scal.swc / 2.0
    values[-1] = 1.0 - scal.sor / 2.0
    model = make_model(sw_map=values, use_saturation_map=True)

    sw = engine_saturation(model, engine_factory)
    assert sw[0] == pytest.approx(scal.swc)
    assert sw[-1] == pytest.approx(1.0 - scal.sor)
    assert np.allclose(sw[1:-1], values[1:-1])


@pytest.mark.parametrize("engine_factory", ENGINES)
def test_clipped_cell_count_reaches_the_log(engine_factory, caplog):
    """İstifadəçi xəritəsinin nə qədərinin dəyişdiyini JURNALDAN görməlidir."""
    scal = default_scal()
    values = two_zone_map(8)
    values[0] = scal.swc / 2.0
    values[-1] = 1.0 - scal.sor / 2.0
    model = make_model(sw_map=values, use_saturation_map=True)

    with caplog.at_level(logging.INFO, logger="imex2d"):
        engine_saturation(model, engine_factory)
    assert "kənarda 2 hüceyrə" in caplog.text


@pytest.mark.parametrize("engine_factory", ENGINES)
def test_equilibration_logs_that_the_transition_zone_is_replaced(engine_factory, caplog):
    model = make_model(nx=8, nz=3, sw_map=two_zone_map(24),
                       use_saturation_map=True, use_equilibration=True,
                       oil_water_contact=2100.0)
    with caplog.at_level(logging.INFO, logger="imex2d"):
        engine_saturation(model, engine_factory)
    assert "ƏVƏZ OLUNDU" in caplog.text


@pytest.mark.parametrize("engine_factory", ENGINES)
def test_equilibration_keeps_hydrostatic_pressure_with_the_map(engine_factory):
    values = two_zone_map(24)          # nx=8, nz=3
    model = make_model(nx=8, nz=3, sw_map=values, use_saturation_map=True,
                       use_equilibration=True, oil_water_contact=2100.0)
    engine = make_service().with_engine(engine_factory).create_engine(
        model, short_config())
    pressure = np.asarray(engine.pressure if engine_factory is ImpesEngine
                          else engine.state.pressure, float)
    sw = np.asarray(engine.sw if engine_factory is ImpesEngine
                    else engine.state.water_saturation, float)

    depths = model.geometry.cell_depths()
    deepest, shallowest = int(np.argmax(depths)), int(np.argmin(depths))
    assert pressure[deepest] > pressure[shallowest]     # hidrostatik
    assert np.allclose(sw, values)                      # kontaktdan ASILI DEYİL


# ═══════════════════════════════════════════════════════════════ OOIP

@pytest.mark.parametrize("engine_factory", ENGINES)
def test_ooip_follows_the_map(engine_factory):
    values = two_zone_map(8)
    with_map = make_model(sw_map=values, use_saturation_map=True)
    without_map = make_model(sw_map=values)

    service = make_service().with_engine(engine_factory)
    mapped = service.create_engine(with_map, short_config()).original_oil_in_place()
    scalar = service.create_engine(without_map, short_config()).original_oil_in_place()
    assert abs(mapped - scalar) > 1.0

    expected = float(np.sum(with_map.pore_volume() * (1.0 - values)
                            / with_map.fluids.oil_fvf))
    assert mapped == pytest.approx(expected, rel=1e-6)


def test_recovery_factor_reacts_to_the_flag():
    """İstifadəçinin GUI-də görəcəyi son nəticə: bayraq açılanda OOIP və
    RF dəyişir, bağlananda ƏVVƏLKİ dəyərə qayıdır."""
    values = two_zone_map(8)
    service, config = make_service(), short_config()

    off = service.run(make_model(sw_map=values), config)
    on = service.run(make_model(sw_map=values, use_saturation_map=True), config)
    off_again = service.run(make_model(sw_map=values), config)

    assert on.ooip != pytest.approx(off.ooip)
    assert on.final_recovery_factor != pytest.approx(off.final_recovery_factor)
    assert off_again.final_recovery_factor == pytest.approx(off.final_recovery_factor)


# ══════════════════════════════════════════════════════════ diaqnostika

def test_unused_map_produces_a_warning():
    report = make_model(sw_map=two_zone_map(8)).diagnose()
    assert any("SW xəritəsi modeldə var" in message
               for message in report.messages(Severity.WARNING))


def test_requested_but_missing_map_is_a_diagnostic_error():
    messages = make_model(use_saturation_map=True).validate()
    assert any("SW xəritəsi yoxdur" in message for message in messages)


def test_used_map_reports_statistics_as_info():
    model = make_model(sw_map=two_zone_map(8), use_saturation_map=True)
    report = model.diagnose()
    assert not report.has_errors
    info = "\n".join(str(item) for item in report.of(Severity.INFO))
    assert "min 0.250" in info and "maks 0.450" in info


def test_model_without_a_map_stays_silent():
    report = make_model().diagnose()
    assert all("SW" not in message for message in report.messages())


# ═════════════════════════════════════════════════════════ serializasiya

def test_flag_survives_a_round_trip():
    serializer = ProjectSerializer()
    model = make_model(sw_map=two_zone_map(8), use_saturation_map=True)
    restored = serializer.reservoir_model_from_dict(
        serializer.reservoir_model_to_dict(model))

    assert restored.initial_conditions.use_saturation_map is True
    assert np.allclose(restored.property_maps["SW"].values, two_zone_map(8))


def test_old_file_without_the_key_loads_with_the_flag_off():
    serializer = ProjectSerializer()
    payload = serializer.reservoir_model_to_dict(make_model())
    payload["initial_conditions"].pop("use_saturation_map")

    restored = serializer.reservoir_model_from_dict(payload)
    assert restored.initial_conditions.use_saturation_map is False


# ═══════════════════════════════════════════════ uyğunlaşdırma (history)

def apply_swc(model, value=0.35):
    parameters = ParameterSet(standard_parameters(model))
    values = parameters.initial_values.copy()
    values[parameters.names.index("SWC")] = value
    return ModelModifier(model, parameters).apply(values)


def test_scalar_saturation_is_left_alone_in_map_mode():
    """SCAL parametri dəyişəndə skalyar Sw TƏNZİMLƏNMİR — o, nəticəyə
    təsir etmir və optimallaşdırıcı onu "uyğunlaşdırmamalıdır"."""
    model = make_model(sw_map=two_zone_map(8), use_saturation_map=True)
    model.initial_conditions.water_saturation = 0.20
    assert apply_swc(model).initial_conditions.water_saturation == pytest.approx(0.20)


def test_scalar_saturation_is_still_reconciled_without_the_map():
    model = make_model()
    model.initial_conditions.water_saturation = 0.20
    assert apply_swc(model).initial_conditions.water_saturation == pytest.approx(0.35)


# ═══════════════════════════════════════════════════════════════════ UI

@pytest.fixture(scope="module")
def qapp():
    QtWidgets = pytest.importorskip("PyQt5.QtWidgets")
    yield QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def numerical_panel():
    from imex2d.ui.panels import NumericalPanel
    return NumericalPanel()


def test_panel_checkbox_reaches_initial_conditions(qapp):
    panel = numerical_panel()
    assert panel.initial_conditions().use_saturation_map is False

    panel.set_saturation_map({"min": 0.25, "mean": 0.35, "max": 0.45})
    panel.use_saturation_map.setChecked(True)
    assert panel.initial_conditions().use_saturation_map is True


def test_panel_disables_the_scalar_field_in_map_mode(qapp):
    panel = numerical_panel()
    panel.set_saturation_map({"min": 0.25, "mean": 0.35, "max": 0.45})
    assert panel.initial_sw.isEnabled()

    panel.use_saturation_map.setChecked(True)
    assert not panel.initial_sw.isEnabled()
    panel.use_saturation_map.setChecked(False)
    assert panel.initial_sw.isEnabled()


def test_panel_without_a_map_cannot_turn_the_option_on(qapp):
    panel = numerical_panel()
    panel.set_saturation_map(None)
    assert not panel.use_saturation_map.isEnabled()
    assert "interpolyasiya" in panel.use_saturation_map.toolTip()

    panel.set_saturation_map({"min": 0.25, "mean": 0.35, "max": 0.45})
    assert panel.use_saturation_map.isEnabled()
    assert "0.250" in panel.saturation_map_info.text()


def test_panel_keeps_an_already_checked_option_reachable(qapp):
    """Fayldan xəritəsiz açılan model istifadəçini tələyə salmamalıdır."""
    panel = numerical_panel()
    panel.set_saturation_map({"min": 0.25, "mean": 0.35, "max": 0.45})
    panel.use_saturation_map.setChecked(True)

    panel.set_saturation_map(None)
    assert panel.use_saturation_map.isEnabled()      # geri söndürülə bilər


def test_main_window_map_statistics_ignore_nan():
    from imex2d.ui.main_window import _saturation_map_stats
    assert _saturation_map_stats({}) is None
    values = np.array([0.25, np.nan, 0.45])
    stats = _saturation_map_stats({"SW": PropertyMap("SW", values)})
    assert stats == pytest.approx({"min": 0.25, "mean": 0.35, "max": 0.45})
