"""Həssaslıq analizi (C6): Tornado və yerli elastiklik."""

import numpy as np

from helpers import default_scal, make_service, short_config
from imex2d.history.parameters import ParameterSet, standard_parameters
from imex2d.history.sensitivity import (OUTPUT_METRICS, SensitivityAnalyzer,
                                        SensitivityReport)
from test_implicit_newton import _rate_controlled


def _model(nx=9):
    return _rate_controlled(nx=nx, scal=default_scal())


def _analyzer(model, names, end_time=1200.0):
    parameters = ParameterSet([item for item in standard_parameters(model)
                               if item.name in names])
    config = short_config(end_time=end_time, snapshots=4)
    return SensitivityAnalyzer(model, parameters, make_service(default_scal()),
                               config), parameters


# ── tornado ───────────────────────────────────────────────────────────
def test_tornado_ranks_sor_above_permeability_for_recovery_factor():
    """SOR — RF-in nəzəri yuxarı həddi, ən böyük hərəkətverici olmalıdır."""
    model = _model()
    analyzer, _ = _analyzer(model, ("SOR", "PERM_MULT"))
    report = analyzer.run_tornado(metric="RF (%)")
    ranked = {item.name: item.swing for item in report.sorted_by_swing()}
    assert ranked["SOR"] > ranked["PERM_MULT"]


def test_breakthrough_time_is_nearly_insensitive_to_permeability_under_rate_control():
    """Buckley-Leverett nəticəsi: RATE idarəetməsində cəbhə sürəti
    vurulan həcmdən asılıdır, mütləq keçiricilikdən demək olar yox."""
    model = _model()
    analyzer, _ = _analyzer(model, ("PERM_MULT", "SOR"))
    report = analyzer.run_tornado(metric="Su gəlişi vaxtı (gün)")
    items = {item.name: item for item in report.items}
    assert items["PERM_MULT"].swing < items["SOR"].swing / 10.0


def test_tornado_baseline_matches_unmodified_model():
    model = _model()
    analyzer, parameters = _analyzer(model, ("SOR",))
    report = analyzer.run_tornado(metric="RF (%)")

    from imex2d.history.parameters import ModelModifier
    direct = make_service(default_scal()).run(
        ModelModifier(model, parameters).apply(parameters.initial_values),
        short_config(end_time=1200.0, snapshots=4))
    assert abs(report.baseline_output - direct.final_recovery_factor) < 1e-9


def test_higher_sor_reduces_recovery_factor():
    """SOR artanda RF azalmalıdır — istiqamət `direction_reversed` ilə uyğun."""
    model = _model()
    analyzer, _ = _analyzer(model, ("SOR",))
    report = analyzer.run_tornado(metric="RF (%)")
    item = report.items[0]
    assert item.direction_reversed
    assert item.low_output > item.high_output


def test_tornado_sorts_from_largest_swing_to_smallest():
    model = _model()
    analyzer, _ = _analyzer(model, ("SOR", "KRW_END", "PERM_MULT"))
    report = analyzer.run_tornado(metric="RF (%)")
    swings = [item.swing for item in report.sorted_by_swing()]
    assert swings == sorted(swings, reverse=True)


def test_all_output_metrics_can_be_computed():
    model = _model()
    analyzer, _ = _analyzer(model, ("SOR",))
    for metric in OUTPUT_METRICS:
        report = analyzer.run_tornado(metric=metric)
        assert np.isfinite(report.baseline_output)
        assert len(report.items) == 1


def test_progress_callback_can_stop_the_scan():
    model = _model()
    analyzer, _ = _analyzer(model, ("SOR", "KRW_END", "PORO_MULT"))

    calls = []

    def stop_after_two(done, total):
        calls.append(done)
        return done < 2

    report = analyzer.run_tornado(metric="RF (%)", progress=stop_after_two)
    assert len(report.items) < 3
    assert calls[-1] >= 2


def test_failed_bound_is_reported_not_silently_dropped():
    """Hədd yığılmasa, baza dəyəri ilə əvəzlənir və işarələnir — susmaq yox."""
    model = _model()
    analyzer, parameters = _analyzer(model, ("SOR",))
    definition = parameters.definitions[0]
    definition.maximum = 0.499999999          # hərəkətli interval demək olar boş

    report = analyzer.run_tornado(metric="RF (%)")
    item = report.items[0]
    if item.failed_high:
        assert item.high_output == report.baseline_output
        assert report.failures >= 1


def test_raises_when_baseline_itself_does_not_converge():
    from imex2d.history.parameters import ParameterDefinition

    model = _model()
    broken = ParameterSet([ParameterDefinition(
        "BROKEN", lambda m, v: setattr(m.scal_parameters, "sor", 0.999),
        minimum=0.9, maximum=0.999, initial=0.999)])
    analyzer = SensitivityAnalyzer(
        model, broken, make_service(default_scal()),
        short_config(end_time=1200.0, snapshots=4))
    try:
        analyzer.run_tornado(metric="RF (%)")
    except RuntimeError:
        return
    raise AssertionError("Yığılmayan baza səssizcə qəbul edildi")


# ── yerli elastiklik ──────────────────────────────────────────────────
def test_local_sensitivity_uses_a_small_step_around_baseline():
    model = _model()
    analyzer, parameters = _analyzer(model, ("SOR",))
    report = analyzer.run_local(metric="RF (%)", step_fraction=0.05)
    item = report.items[0]

    span = report.items[0].high_value - report.items[0].low_value
    definition = parameters.definitions[0]
    assert span < (definition.maximum - definition.minimum) * 0.3


def test_local_sensitivity_swing_grows_with_step_size():
    model = _model()
    analyzer, _ = _analyzer(model, ("SOR",))
    small = analyzer.run_local(metric="RF (%)", step_fraction=0.03).items[0]
    large = analyzer.run_local(metric="RF (%)", step_fraction=0.20).items[0]
    assert large.swing > small.swing


def test_tornado_and_local_can_disagree_in_ranking():
    """Geniş hədli, lakin lokal az təsirli parametr fərqli sıralana bilər.

    Bu, iki üsulun eyni şeyi ölçmədiyini göstərir — testin özü
    fərqi TƏLƏB ETMİR, sadəcə hər ikisinin işlədiyini yoxlayır.
    """
    model = _model()
    analyzer, _ = _analyzer(model, ("SOR", "KRW_END", "PORO_MULT", "SWC"))
    tornado = analyzer.run_tornado(metric="RF (%)")
    local = analyzer.run_local(metric="RF (%)", step_fraction=0.1)
    assert len(tornado.items) == len(local.items) == 4
    assert all(np.isfinite(item.swing) for item in tornado.items)
    assert all(np.isfinite(item.swing) for item in local.items)


def test_report_text_lists_every_parameter():
    model = _model()
    analyzer, _ = _analyzer(model, ("SOR", "PERM_MULT"))
    text = analyzer.run_tornado(metric="RF (%)").as_text()
    assert "SOR" in text and "PERM_MULT" in text
