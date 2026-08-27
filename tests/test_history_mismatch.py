"""History matching: müşahidə məlumatı və uyğunsuzluq ölçüsü (C5/1)."""

import os
import tempfile

import numpy as np

from helpers import default_scal, make_service, short_config
from imex2d.domain.observations import (ObservationSet, ObservedQuantity,
                                        ObservedSeries)
from imex2d.history.mismatch import (DEFAULT_WEIGHTS, MismatchCalculator,
                                     MismatchReport, SeriesMismatch)
from imex2d.history.observation_io import (ObservationFormatError,
                                           read_observations_csv,
                                           write_observations_csv)
from test_implicit_newton import _rate_controlled


def _truth(nx=11, end_time=500.0):
    scal = default_scal()
    model = _rate_controlled(nx=nx, scal=scal)
    return model, make_service(scal).run(
        model, short_config(end_time=end_time, snapshots=4))


def _observations_from(result, points=15, noise=0.0, seed=1):
    times = np.asarray(result.series.time, float)
    sample_times = np.linspace(times[1], times[-1] * 0.97, points)
    rng = np.random.default_rng(seed)

    def sample(values):
        interpolated = np.interp(sample_times, times, np.asarray(values, float))
        if noise:
            interpolated = interpolated * (1 + rng.normal(0, noise,
                                                          interpolated.size))
        return interpolated

    return ObservationSet(series=[
        ObservedSeries("", ObservedQuantity.OIL_RATE, sample_times,
                       sample(result.series.oil_rate)),
        ObservedSeries("", ObservedQuantity.CUMULATIVE_OIL, sample_times,
                       sample(result.series.cumulative_oil)),
        ObservedSeries("", ObservedQuantity.AVERAGE_PRESSURE, sample_times,
                       sample(result.series.average_pressure)),
    ])


def _write(text: str) -> str:
    handle, path = tempfile.mkstemp(suffix=".csv")
    with os.fdopen(handle, "w", encoding="utf-8") as file:
        file.write(text)
    return path


# ── domain ────────────────────────────────────────────────────────────
def test_series_validation_catches_common_errors():
    good = ObservedSeries("W-1", ObservedQuantity.OIL_RATE,
                          [0.0, 10.0, 20.0], [5.0, 4.0, 3.0])
    assert good.validate() == []

    assert ObservedSeries("W-1", ObservedQuantity.OIL_RATE,
                          [0.0, 10.0], [1.0]).validate()
    assert ObservedSeries("W-1", ObservedQuantity.OIL_RATE,
                          [10.0, 5.0], [1.0, 2.0]).validate()
    assert ObservedSeries("W-1", ObservedQuantity.OIL_RATE,
                          [0.0], [1.0]).validate()


def test_field_level_series_get_a_readable_label():
    field = ObservedSeries("", ObservedQuantity.OIL_RATE,
                           [0.0, 1.0], [1.0, 2.0])
    well = ObservedSeries("PROD-1", ObservedQuantity.OIL_RATE,
                          [0.0, 1.0], [1.0, 2.0])
    assert field.label.startswith("Yataq")
    assert well.label.startswith("PROD-1")


def test_observation_set_summarises_content():
    _, result = _truth()
    observations = _observations_from(result)
    summary = observations.summary()
    assert summary["sıra"] == 3
    assert summary["nöqtə"] == 45


# ── uyğunsuzluq ───────────────────────────────────────────────────────
def test_perfect_match_gives_near_zero_mismatch():
    """Model öz nəticəsi ilə müqayisə olunanda uyğunsuzluq sıfıra yaxındır.

    Tam sıfır olmur, çünki müşahidə vaxtları model addımlarına düşmür
    və interpolyasiya kiçik xəta verir.
    """
    _, result = _truth()
    report = MismatchCalculator().evaluate(result, _observations_from(result))
    assert report.total < 0.02, report.as_text()
    assert not report.skipped


def test_wrong_model_gives_a_much_larger_mismatch():
    """Ölçü həssas olmalıdır — səhv model aydın seçilməlidir."""
    scal = default_scal()
    _, truth = _truth(nx=11)
    observations = _observations_from(truth)

    wrong = _rate_controlled(nx=11, scal=scal)
    wrong.rock.permx.values[:] *= 0.3
    wrong.rock.permy.values[:] *= 0.3
    wrong_result = make_service(scal).run(wrong, short_config(end_time=500.0,
                                                             snapshots=4))

    good_score = MismatchCalculator().evaluate(truth, observations).total
    bad_score = MismatchCalculator().evaluate(wrong_result, observations).total
    assert bad_score > good_score * 5, (good_score, bad_score)


def test_mismatch_is_dimensionless_and_scale_invariant():
    """Vahid dəyişəndə NRMSE dəyişməməlidir.

    Bu vacibdir: təzyiq barla, debit m³/günlə ölçülür. Xam SSE-ləri
    toplasaydıq, böyük ədədli kəmiyyət yekun ölçüyə hakim olardı.
    """
    time = np.linspace(0, 100, 11)
    observed = np.linspace(10.0, 20.0, 11)
    simulated = observed + 0.5

    base = SeriesMismatch("a", ObservedQuantity.OIL_RATE, "", time,
                          observed, simulated)
    scaled = SeriesMismatch("b", ObservedQuantity.OIL_RATE, "", time,
                            observed * 1000.0, simulated * 1000.0)
    assert abs(base.nrmse - scaled.nrmse) < 1e-12
    assert scaled.rmse > base.rmse * 100


def test_bias_reveals_systematic_error():
    time = np.linspace(0, 10, 6)
    observed = np.full(6, 100.0)
    high = SeriesMismatch("a", ObservedQuantity.OIL_RATE, "", time,
                          observed, observed + 5.0)
    low = SeriesMismatch("b", ObservedQuantity.OIL_RATE, "", time,
                         observed, observed - 5.0)
    assert high.bias > 0 and low.bias < 0
    assert abs(high.rmse - low.rmse) < 1e-12      # RMSE işarəni gizlədir


def test_constant_observation_does_not_divide_by_zero():
    time = np.linspace(0, 10, 5)
    constant = np.full(5, 250.0)
    mismatch = SeriesMismatch("a", ObservedQuantity.AVERAGE_PRESSURE, "",
                              time, constant, constant + 2.5)
    assert np.isfinite(mismatch.nrmse)
    assert mismatch.nrmse > 0


def test_observations_outside_the_simulated_period_are_skipped():
    """Susmaq təhlükəlidir — atlanan sıra hesabatda görünməlidir."""
    _, result = _truth(end_time=300.0)
    late = ObservationSet(series=[
        ObservedSeries("", ObservedQuantity.OIL_RATE,
                       [900.0, 1000.0], [10.0, 8.0])])
    report = MismatchCalculator().evaluate(result, late)
    assert report.skipped
    assert not report.series


def test_unknown_well_is_skipped_not_silently_ignored():
    _, result = _truth()
    observations = ObservationSet(series=[
        ObservedSeries("YOXDUR-1", ObservedQuantity.OIL_RATE,
                       [50.0, 100.0], [10.0, 9.0])])
    report = MismatchCalculator().evaluate(result, observations)
    assert report.skipped == [observations.series[0].label]


def test_per_well_series_use_well_curves():
    scal = default_scal()
    model = _rate_controlled(nx=9, scal=scal)
    result = make_service(scal).run(model, short_config(end_time=200.0))
    producer = next(name for name in result.well_oil_rate)

    times = np.asarray(result.series.time, float)
    sample_times = np.linspace(times[1], times[-1] * 0.9, 8)
    observations = ObservationSet(series=[
        ObservedSeries(producer, ObservedQuantity.OIL_RATE, sample_times,
                       np.interp(sample_times, times,
                                 result.well_oil_rate[producer]))])
    report = MismatchCalculator().evaluate(result, observations)
    assert report.series and report.series[0].nrmse < 0.02


def test_weights_change_the_total_but_not_the_individual_scores():
    _, result = _truth()
    observations = _observations_from(result, noise=0.05, seed=3)
    plain = MismatchCalculator().evaluate(result, observations)
    weighted = MismatchCalculator(
        {ObservedQuantity.AVERAGE_PRESSURE: 20.0}).evaluate(result,
                                                            observations)
    assert abs(plain.total - weighted.total) > 1e-9
    plain_scores = {item.label: item.nrmse for item in plain.series}
    for item in weighted.series:
        assert abs(plain_scores[item.label] - item.nrmse) < 1e-12


def test_report_identifies_the_worst_series():
    _, result = _truth()
    observations = _observations_from(result, noise=0.05, seed=5)
    report = MismatchCalculator().evaluate(result, observations)
    assert report.worst is not None
    assert report.worst.nrmse == max(item.nrmse for item in report.series)


def test_empty_result_reports_everything_as_skipped():
    from imex2d.simulation.results import SimulationResult

    _, result = _truth()
    report = MismatchCalculator().evaluate(SimulationResult(),
                                           _observations_from(result))
    assert not report.series
    assert len(report.skipped) == 3
    assert report.total == float("inf")


# ── CSV ───────────────────────────────────────────────────────────────
def test_reads_long_format_csv():
    path = _write("""time,well,quantity,value
30,PROD-1,OIL_RATE,142.5
60,PROD-1,OIL_RATE,138.0
30,PROD-1,WATER_RATE,3.1
60,PROD-1,WATER_RATE,7.4
30,,AVERAGE_PRESSURE,247.8
60,,AVERAGE_PRESSURE,244.1
""")
    try:
        observations = read_observations_csv(path)
        assert len(observations) == 3
        assert observations.wells == ["PROD-1"]
        oil = observations.get("PROD-1", ObservedQuantity.OIL_RATE)
        assert list(oil.time) == [30.0, 60.0]
    finally:
        os.unlink(path)


def test_eclipse_style_names_are_recognised():
    path = _write("""time,well,quantity,value
10,W1,WOPR,100
20,W1,WOPR,95
10,,FPR,250
20,,FPR,248
""")
    try:
        observations = read_observations_csv(path)
        assert observations.get("W1", ObservedQuantity.OIL_RATE)
        assert observations.get("", ObservedQuantity.AVERAGE_PRESSURE)
    finally:
        os.unlink(path)


def test_rows_are_sorted_by_time():
    path = _write("""time,well,quantity,value
90,W1,OIL_RATE,80
30,W1,OIL_RATE,120
60,W1,OIL_RATE,100
""")
    try:
        series = read_observations_csv(path).series[0]
        assert list(series.time) == [30.0, 60.0, 90.0]
        assert list(series.values) == [120.0, 100.0, 80.0]
    finally:
        os.unlink(path)


def test_missing_columns_are_rejected():
    path = _write("time,value\n10,100\n20,90\n")
    try:
        read_observations_csv(path)
    except ObservationFormatError:
        return
    finally:
        os.unlink(path)
    raise AssertionError("Natamam CSV qəbul edildi")


def test_unknown_quantity_is_reported():
    path = _write("time,well,quantity,value\n10,W1,NAMELUM,5\n20,W1,NAMELUM,6\n")
    try:
        read_observations_csv(path)
    except ObservationFormatError as error:
        assert "NAMELUM" in str(error)
        return
    finally:
        os.unlink(path)
    raise AssertionError("Naməlum kəmiyyət qəbul edildi")


def test_csv_round_trip_preserves_values():
    _, result = _truth()
    observations = _observations_from(result)
    handle, path = tempfile.mkstemp(suffix=".csv")
    os.close(handle)
    try:
        write_observations_csv(path, observations)
        restored = read_observations_csv(path)
    finally:
        os.unlink(path)

    assert len(restored) == len(observations)
    for original in observations.series:
        copy = restored.get(original.well, original.quantity)
        assert copy is not None, original.label
        assert np.allclose(copy.time, original.time, atol=1e-3)
        assert np.allclose(copy.values, original.values, rtol=1e-4)


def test_default_weights_favour_cumulative_over_daily_rate():
    """Gündəlik debitdə səs-küy çox olur, toplam isə ehtiyatı göstərir."""
    assert (DEFAULT_WEIGHTS[ObservedQuantity.CUMULATIVE_OIL]
            > DEFAULT_WEIGHTS[ObservedQuantity.OIL_RATE])
