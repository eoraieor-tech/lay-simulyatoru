"""Seans 45 — günlük göstəricilər (`reporting/daily.py`).

Mühərrik adaptiv addımlarla gedir; günlük cədvəl addımlardan qurulur.
Testlər üç müqaviləni kilidləyir:

  1. Debit addım boyu SABİTDİR → günün debiti = günün həcmi / uzunluğu,
     günlük həcmlərin cəmi son kumulyativə BƏRABƏRDİR.
  2. Addımlar tam günlərə düşəndə cədvəl xam nəticə ilə EYNİDİR.
  3. Nəticədə olmayan sıra sütun yaratmır; THP-nin `nan`-ı sıfıra dönmür.

Əlavə olaraq: mühərriklər vurucular üzrə sıranı yazır və onların cəmi
sahə vurmasına bərabərdir (Seans 45-ə qədər vurucu üzrə heç nə yox idi).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from helpers import five_spot_model, make_service, short_config
from imex2d.application.serialization import ProjectSerializer
from imex2d.reporting import daily, results_export
from imex2d.simulation.results import SimulationResult, TimeSeries
from test_results_export import _run


def _synthetic(times, oil, water=None, pressure=None, bhp=None, thp=None,
               ooip=1000.0) -> SimulationResult:
    """Əl ilə qurulmuş nəticə — mühərrikin qaydası ilə (kumulyativ = Σ q·Δt)."""
    times = list(map(float, times))
    oil = list(map(float, oil))
    water = list(map(float, water if water is not None else [0.0] * len(oil)))
    steps = np.diff(np.concatenate(([0.0], times)))
    cumulative_oil = list(np.cumsum(np.array(oil) * steps))
    series = TimeSeries(
        time=times, oil_rate=oil, water_rate=water,
        cumulative_oil=cumulative_oil,
        average_pressure=list(pressure if pressure is not None
                              else [200.0] * len(oil)),
        recovery_factor=[c / ooip * 100.0 for c in cumulative_oil])
    result = SimulationResult(series=series, ooip=ooip)
    result.well_oil_rate = {"P1": list(oil)}
    result.well_water_rate = {"P1": list(water)}
    if bhp is not None:
        result.well_bhp = {"P1": list(bhp)}
    if thp is not None:
        result.well_thp = {"P1": list(thp)}
    return result


# ═══════════════════════════════ günlər ══════════════════════════════

def test_day_edges_whole_and_fractional_end():
    assert list(daily.day_edges(3.0)) == [1.0, 2.0, 3.0]
    assert list(daily.day_edges(2.5)) == [1.0, 2.0, 2.5]
    # yuvarlaqlaşdırma xətası əlavə "gün" yaratmır
    assert list(daily.day_edges(3.0 - 1e-12)) == [1.0, 2.0, 3.0]
    assert len(daily.day_edges(0.0)) == 0


def test_every_simulated_day_has_a_row():
    result = _synthetic([0.25, 7.0, 30.0], [5.0, 6.0, 7.0])
    table = daily.daily_table(result)
    assert len(table) == 30
    assert list(table.days) == list(range(1, 31))


# ═════════════════════════ debit və kumulyativ ═══════════════════════

def test_day_split_by_step_boundary_is_time_weighted():
    """Addım sərhədi günün ortasındadır → iki debitin zamanla çəkili ortası."""
    result = _synthetic([0.5, 2.5, 3.0], [10.0, 20.0, 40.0])
    field = daily.daily_table(result).field_columns

    assert field[daily.OIL] == pytest.approx([15.0, 20.0, 30.0], abs=1e-12)
    assert field[daily.CUM_OIL] == pytest.approx([15.0, 35.0, 65.0], abs=1e-12)
    assert field[daily.RF] == pytest.approx([1.5, 3.5, 6.5], abs=1e-12)


def test_daily_volumes_sum_to_engine_cumulative():
    """Əsas müqavilə real qaçışda: Σ günlük həcm = mühərrikin kumulyativi."""
    result = _run(with_gas=False)
    table = daily.daily_table(result)
    lengths = np.diff(np.concatenate(([0.0], table.days)))
    field = table.field_columns

    assert np.sum(field[daily.OIL] * lengths) == pytest.approx(
        result.series.cumulative_oil[-1], rel=1e-12)
    assert field[daily.CUM_OIL][-1] == pytest.approx(
        result.series.cumulative_oil[-1], rel=1e-12)
    assert field[daily.RF][-1] == pytest.approx(
        result.series.recovery_factor[-1], rel=1e-12)
    well = table.wells["PROD-1"]
    assert np.sum(well[daily.OIL] * lengths) == pytest.approx(
        sum(q * dt for q, dt in zip(
            result.well_oil_rate["PROD-1"],
            np.diff(np.concatenate(([0.0], result.series.time))))), rel=1e-12)


def test_whole_day_steps_reproduce_raw_result():
    """Addımlar tam günlərdədirsə, heç nə interpolyasiya olunmur."""
    result = _synthetic([1.0, 2.0, 3.0], [10.0, 12.0, 9.0],
                        water=[1.0, 3.0, 9.0], pressure=[250.0, 240.0, 235.0],
                        bhp=[100.0, 110.0, 120.0])
    table = daily.daily_table(result)
    field = table.field_columns

    assert field[daily.OIL] == pytest.approx(result.series.oil_rate, abs=1e-12)
    assert field[daily.PRESSURE] == pytest.approx(
        result.series.average_pressure, abs=1e-12)
    assert field[daily.RF] == pytest.approx(result.series.recovery_factor,
                                            abs=1e-12)
    assert field[daily.WATER_CUT] == pytest.approx([100 / 11, 20.0, 50.0],
                                                   abs=1e-12)
    assert table.wells["P1"][daily.BHP] == pytest.approx([100.0, 110.0, 120.0])


def test_fractional_last_day_uses_its_own_length():
    """Qaçış 2.5-ci gündə bitib: son "gün" 0.5 gündür, debit ona bölünür."""
    result = _synthetic([2.5], [8.0])
    table = daily.daily_table(result)
    assert list(table.days) == [1.0, 2.0, 2.5]
    assert table.field_columns[daily.OIL] == pytest.approx([8.0, 8.0, 8.0])
    assert table.field_columns[daily.CUM_OIL][-1] == pytest.approx(20.0)


# ═══════════════════════ təzyiq, BHP, THP ════════════════════════════

def test_pressure_is_interpolated_between_step_ends():
    result = _synthetic([0.5, 2.5, 3.0], [1.0, 1.0, 1.0],
                        pressure=[200.0, 180.0, 170.0])
    pressure = daily.daily_table(result).field_columns[daily.PRESSURE]
    assert pressure == pytest.approx([195.0, 185.0, 170.0], abs=1e-12)
    assert "interp" in daily.PRESSURE, "təxmini olduğu sütun adında görünməlidir"


def test_bhp_taken_from_step_covering_day_end():
    """Gün (d-1, d]: sonu 1.0 olan addım 1-ci günə aiddir, 2-ciyə yox."""
    result = _synthetic([0.5, 1.0, 2.7, 3.0], [1.0] * 4,
                        bhp=[90.0, 91.0, 92.0, 93.0])
    bhp = daily.daily_table(result).wells["P1"][daily.BHP]
    assert list(bhp) == [91.0, 92.0, 93.0]


def test_thp_nan_is_kept_not_zeroed():
    result = _synthetic([1.0, 2.0], [1.0, 0.0], thp=[15.0, float("nan")])
    thp = daily.daily_table(result).wells["P1"][daily.THP]
    assert thp[0] == 15.0
    assert math.isnan(thp[1])


# ═════════════════════════ olmayan sıralar ═══════════════════════════

def test_missing_series_create_no_columns():
    result = _synthetic([1.0, 2.0], [1.0, 2.0])
    table = daily.daily_table(result)
    for absent in (daily.GAS, daily.GOR, daily.GAS_INJ, daily.WATER_INJ):
        assert absent not in table.field_columns
    assert daily.BHP not in table.wells["P1"]


def test_series_with_wrong_length_is_skipped():
    """Uyğunsuz uzunluqlu sıra uydurma dəyərlə doldurulmur."""
    result = _synthetic([1.0, 2.0], [1.0, 2.0])
    result.well_bhp = {"P1": [100.0]}
    assert daily.BHP not in daily.daily_table(result).wells["P1"]


def test_empty_result_gives_empty_table():
    table = daily.daily_table(SimulationResult())
    assert len(table) == 0
    assert table.targets() == [daily.FIELD]


def test_row_and_day_lookup():
    table = daily.daily_table(_synthetic([0.5, 2.5, 3.0], [10.0, 20.0, 40.0]))
    assert table.index_of_day(2) == 1
    assert table.index_of_day(1.4) == 1          # (1, 2] → 2-ci gün
    assert table.index_of_day(999) == 2
    row = table.row(0, "P1")
    assert row[daily.DAY] == 1.0 and row[daily.OIL] == pytest.approx(15.0)


# ═══════════════════ mühərriklər: vurucu üzrə sıra ═══════════════════

@pytest.mark.parametrize("with_gas", [False, True],
                         ids=["iki-fazalı FIM", "üç fazalı"])
def test_fim_engines_record_injector_rates(with_gas):
    result = _run(with_gas=with_gas)
    assert set(result.well_water_injection_rate) == {"INJ-1"}
    injected = np.sum([result.well_water_injection_rate[name]
                       for name in result.well_water_injection_rate], axis=0)
    assert injected == pytest.approx(result.series.water_injection_rate,
                                     rel=1e-12)
    table = daily.daily_table(result)
    assert daily.WATER_INJ in table.wells["INJ-1"]
    assert daily.OIL not in table.wells["INJ-1"]


def test_impes_engine_records_injector_rates():
    model = five_spot_model(nx=9, ny=9)
    result = make_service().run(model, short_config(end_time=60.0))
    assert result.converged, result.message
    assert set(result.well_water_injection_rate) == {"INJ-1"}
    assert result.well_water_injection_rate["INJ-1"] == pytest.approx(
        result.series.water_injection_rate, rel=1e-12)


# ═══════════════════════ ixrac və saxlama ════════════════════════════

def test_daily_csv_round_trip(tmp_path):
    result = _run(with_gas=False)
    table = daily.daily_table(result)
    path = str(tmp_path / "gunluk.csv")
    results_export.write_daily_csv(result, path)

    back = results_export.read_csv(path)
    assert back[daily.DAY] == list(table.days)
    assert back[daily.OIL] == list(table.field_columns[daily.OIL])
    assert back[f"PROD-1: {daily.OIL}"] == list(table.wells["PROD-1"][daily.OIL])
    assert back[f"INJ-1: {daily.WATER_INJ}"] == list(
        table.wells["INJ-1"][daily.WATER_INJ])


def test_project_file_keeps_injector_rates():
    result = _run(with_gas=False)
    data = ProjectSerializer._result_to_dict(result, include_snapshots=False)
    back = ProjectSerializer._result_from_dict(data)
    assert back.well_water_injection_rate == result.well_water_injection_rate

    # köhnə layihə faylında açar yoxdur — xəta yox, boş lüğət
    del data["well_water_injection_rate"], data["well_gas_injection_rate"]
    old = ProjectSerializer._result_from_dict(data)
    assert old.well_water_injection_rate == {}


# ═══════════════════════ UI modeli və qrafik ═════════════════════════

def test_table_model_shows_selected_target_and_blanks_nan():
    from PyQt5.QtCore import Qt
    from imex2d.ui.daily_view import BLANK, DailyTableModel, format_value

    result = _synthetic([1.0, 2.0], [10.0, 0.0], thp=[15.0, float("nan")])
    model = DailyTableModel()
    model.set_table(daily.daily_table(result), "P1")
    assert model.target == "P1"
    assert model.rowCount() == 2
    assert model.headers[0] == daily.DAY and daily.THP in model.headers
    thp_column = model.headers.index(daily.THP)
    assert model.data(model.index(1, thp_column), Qt.DisplayRole) == BLANK
    assert model.data(model.index(0, 0), Qt.DisplayRole) == "1"

    # olmayan obyekt → yataq
    model.set_table(daily.daily_table(result), "YOXDUR")
    assert model.target == daily.FIELD
    assert format_value(daily.CUM_OIL, 12345.67) == "12 345.7"


def test_daily_renderer_draws_steps_and_marker():
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
    from imex2d.rendering.renderers import DailyRenderer

    result = _run(with_gas=True)
    table = daily.daily_table(result)
    figure = Figure()
    axes = figure.subplots(1, 2)
    renderer = DailyRenderer()
    for _ in range(2):                    # təkrar çəkmə ikinci oxları yığmır
        renderer.draw(figure, axes, table.days, table.columns("PROD-1"),
                      "PROD-1", selected_day=100.0)
    rate_ax, cumulative_ax = axes
    labels = [line.get_label() for line in rate_ax.get_lines()]
    assert "Neft" in labels and "Su" in labels
    assert len(figure.axes) == 3, "qaz üçün TƏK ikinci ox olmalıdır"
    assert any(list(line.get_xdata()) == [100.0, 100.0]
               for line in cumulative_ax.get_lines())
