"""Avtomatik uyğunlaşdırma (C5/3).

Əsas yoxlama üsulu — TWIN EXPERIMENT: modelin parametrləri məlum
dəyərlərə qoyulur, nəticə "müşahidə" kimi işlədilir, sonra
optimallaşdırıcının həmin dəyərləri bərpa edib-etmədiyi yoxlanılır.
"""

import numpy as np

from helpers import default_scal, make_service, short_config
from imex2d.domain.observations import (ObservationSet, ObservedQuantity,
                                        ObservedSeries)
from imex2d.history.mismatch import MismatchCalculator
from imex2d.history.optimizer import (FAILURE_PENALTY, Evaluation,
                                      HistoryMatchingService, MatchResult)
from imex2d.history.parameters import (ModelModifier, ParameterSet,
                                       standard_parameters)
from test_implicit_newton import _rate_controlled

END_TIME = 1000.0
TRUE_PERMEABILITY = 2.5
TRUE_RESIDUAL_OIL = 0.32


def _base_model(nx=9):
    return _rate_controlled(nx=nx, scal=default_scal())


def _config(end_time=END_TIME):
    return short_config(end_time=end_time, snapshots=4)


def _synthetic_truth(model, include_pressure=True, noise=0.0, seed=1,
                     points=12):
    """Məlum parametrlərlə "həqiqət" qurur və ondan müşahidə çıxarır."""
    scal = default_scal()
    parameters = ParameterSet(standard_parameters(model))
    values = parameters.initial_values.copy()
    values[parameters.names.index("PERM_MULT")] = TRUE_PERMEABILITY
    values[parameters.names.index("SOR")] = TRUE_RESIDUAL_OIL

    truth = make_service(scal).run(
        ModelModifier(model, parameters).apply(values), _config())

    times = np.asarray(truth.series.time, float)
    sample_times = np.linspace(times[1], times[-1] * 0.95, points)
    rng = np.random.default_rng(seed)

    def sample(values_):
        interpolated = np.interp(sample_times, times,
                                 np.asarray(values_, float))
        if noise:
            interpolated = interpolated * (1 + rng.normal(0, noise,
                                                          interpolated.size))
        return np.maximum(interpolated, 0.0)

    series = [
        ObservedSeries("", ObservedQuantity.OIL_RATE, sample_times,
                       sample(truth.series.oil_rate)),
        ObservedSeries("", ObservedQuantity.CUMULATIVE_OIL, sample_times,
                       sample(truth.series.cumulative_oil)),
        ObservedSeries("", ObservedQuantity.WATER_CUT, sample_times,
                       sample(truth.series.water_cut)),
    ]
    if include_pressure:
        series.append(ObservedSeries(
            "", ObservedQuantity.AVERAGE_PRESSURE, sample_times,
            sample(truth.series.average_pressure)))
    return truth, ObservationSet(series=series)


def _two_parameter_set(model):
    return ParameterSet([item for item in standard_parameters(model)
                         if item.name in ("PERM_MULT", "SOR")])


def _service(model, observations, parameters=None):
    return HistoryMatchingService(
        model, parameters or _two_parameter_set(model), observations,
        make_service(default_scal()), _config())


# ── qiymətləndirmə ────────────────────────────────────────────────────
def test_evaluation_returns_the_mismatch():
    model = _base_model()
    _, observations = _synthetic_truth(model)
    service = _service(model, observations)
    parameters = service.parameters

    at_truth = service.evaluate(parameters.to_unit(
        np.array([TRUE_PERMEABILITY, TRUE_RESIDUAL_OIL])))
    away = service.evaluate(parameters.to_unit(np.array([0.3, 0.1])))
    assert at_truth < 1e-3
    assert away > at_truth * 10


def test_repeated_points_are_cached():
    """Hər qiymətləndirmə bir simulyasiyadır — təkrar hesablama israfdır."""
    model = _base_model()
    _, observations = _synthetic_truth(model)
    service = _service(model, observations)
    point = np.array([0.5, 0.5])
    service.evaluate(point)
    service.evaluate(point)
    assert len(service.history) == 1


def test_invalid_parameters_return_a_penalty_not_an_exception():
    """Optimallaşdırıcı hədlərin kənarını sınayır — istisna axtarışı dayandırardı."""
    model = _base_model()
    _, observations = _synthetic_truth(model)
    parameters = _two_parameter_set(model)
    service = HistoryMatchingService(model, parameters, observations,
                                     make_service(default_scal()), _config())

    class _Broken:
        def run(self, *args, **kwargs):
            raise RuntimeError("simulyasiya çökdü")

    service.service = _Broken()
    score = service.evaluate(np.array([0.5, 0.5]))
    assert score == FAILURE_PENALTY
    assert service.history[-1].succeeded is False
    assert "çökdü" in service.history[-1].message


def test_penalty_is_finite():
    """Sonsuzluq Nelder-Mead simpleksini və DE seçimini pozur."""
    assert np.isfinite(FAILURE_PENALTY)
    assert FAILURE_PENALTY > 1e3


def test_evaluation_is_recorded_with_physical_values():
    model = _base_model()
    _, observations = _synthetic_truth(model)
    service = _service(model, observations)
    service.evaluate(np.array([1.0, 0.0]))
    record = service.history[-1]
    assert isinstance(record, Evaluation)
    assert record.values[0] > 1.0          # PERM_MULT yuxarı həddə
    assert record.seconds > 0


# ── twin experiment ───────────────────────────────────────────────────
def test_nelder_mead_recovers_the_true_parameters():
    """Ən vacib test: gizlədilmiş parametrlər bərpa olunmalıdır."""
    model = _base_model()
    _, observations = _synthetic_truth(model)
    result = _service(model, observations).run("Nelder-Mead",
                                               max_evaluations=45)
    recovered = result.as_dict()
    assert abs(recovered["PERM_MULT"] - TRUE_PERMEABILITY) < 0.3, recovered
    assert abs(recovered["SOR"] - TRUE_RESIDUAL_OIL) < 0.03, recovered
    assert result.best_mismatch < result.initial_mismatch / 10.0


def test_powell_also_converges():
    model = _base_model()
    _, observations = _synthetic_truth(model)
    result = _service(model, observations).run("Powell", max_evaluations=45)
    assert abs(result.as_dict()["PERM_MULT"] - TRUE_PERMEABILITY) < 0.5
    assert result.improvement > 90.0


def test_differential_evolution_runs_within_its_budget():
    model = _base_model(nx=7)
    _, observations = _synthetic_truth(model)
    result = _service(model, observations).run("Differential Evolution",
                                               max_evaluations=40, seed=3)
    assert result.evaluations <= 120
    assert result.best_mismatch < result.initial_mismatch


# ── nəticə obyekti ────────────────────────────────────────────────────
def test_result_reports_improvement_and_history():
    model = _base_model()
    _, observations = _synthetic_truth(model)
    result = _service(model, observations).run("Nelder-Mead",
                                               max_evaluations=25)
    assert isinstance(result, MatchResult)
    assert result.evaluations > 1
    assert 0.0 <= result.improvement <= 100.0
    assert result.best_report is not None
    assert result.best_result is not None
    assert "PERM_MULT" in result.summary()


def test_convergence_curve_is_monotonic():
    """Əyri o ana qədərki ƏN YAXŞI nəticəni göstərir — pisləşə bilməz."""
    model = _base_model()
    _, observations = _synthetic_truth(model)
    curve = _service(model, observations).run(
        "Nelder-Mead", max_evaluations=25).convergence_curve
    assert np.all(np.diff(curve) <= 1e-12)


def test_best_values_stay_within_the_bounds():
    model = _base_model()
    _, observations = _synthetic_truth(model)
    result = _service(model, observations).run("Nelder-Mead",
                                               max_evaluations=25)
    for definition, value in zip(result.parameters.definitions,
                                 result.best_values):
        assert definition.minimum - 1e-9 <= value <= definition.maximum + 1e-9


def test_search_never_worsens_the_starting_point():
    model = _base_model()
    _, observations = _synthetic_truth(model)
    result = _service(model, observations).run("Nelder-Mead",
                                               max_evaluations=20)
    assert result.best_mismatch <= result.initial_mismatch + 1e-12


# ── idarəetmə ─────────────────────────────────────────────────────────
def test_progress_callback_can_stop_the_search():
    model = _base_model()
    _, observations = _synthetic_truth(model)
    service = _service(model, observations)

    def stop_after_five(evaluation):
        return evaluation.iteration < 5

    result = service.run("Nelder-Mead", max_evaluations=60,
                         progress=stop_after_five)
    assert result.stopped_early
    assert result.evaluations < 30


def test_cancelled_search_still_returns_the_best_so_far():
    model = _base_model()
    _, observations = _synthetic_truth(model)
    service = _service(model, observations)

    def stop_immediately(evaluation):
        return evaluation.iteration < 3

    result = service.run("Nelder-Mead", max_evaluations=40,
                         progress=stop_immediately)
    assert np.isfinite(result.best_mismatch)
    assert len(result.best_values) == len(service.parameters)


def test_initial_simplex_moves_every_parameter():
    """Defolt simpleks nisbi addım işlədir — sıfıra yaxın parametrlər
    üçün praktiki olaraq hərəkətsiz qalır."""
    start = np.array([0.0, 0.5, 1.0])
    simplex = HistoryMatchingService._initial_simplex(start, step=0.15)
    assert simplex.shape == (4, 3)
    for index in range(3):
        assert abs(simplex[index + 1, index] - start[index]) > 0.1
        assert 0.0 <= simplex[index + 1, index] <= 1.0


# ── təyin olunma qabiliyyəti ──────────────────────────────────────────
def test_pressure_observations_sharpen_the_permeability_signal():
    """Müşahidə dəsti parametri təyin edə bilməyə bilər.

    RATE ilə idarə olunan vurucuda neft debiti vurulan həcmlə
    müəyyənləşir və keçiricilikdən demək olar asılı deyil. Yalnız
    təzyiq keçiriciliyə həssasdır — o olmadan xəta funksiyası yastı
    qalır və səs-küy siqnalı basdırır.
    """
    scal = default_scal()
    model = _base_model()
    calculator = MismatchCalculator()
    parameters = ParameterSet(standard_parameters(model))
    modifier = ModelModifier(model, parameters)

    _, without = _synthetic_truth(model, include_pressure=False)
    _, with_pressure = _synthetic_truth(model, include_pressure=True)

    def penalty_at(observations, multiplier):
        values = parameters.initial_values.copy()
        values[parameters.names.index("PERM_MULT")] = multiplier
        values[parameters.names.index("SOR")] = TRUE_RESIDUAL_OIL
        result = make_service(scal).run(modifier.apply(values), _config())
        return calculator.evaluate(result, observations).total

    # Siqnal AŞAĞI keçiricilik tərəfində kəskinləşir: orada təzyiq
    # düşməsi böyükdür və orta lay təzyiqi həqiqətdən aydın ayrılır.
    # Yuxarı tərəfdə təzyiq onsuz da yastıdır, ona görə fərq azdır.
    assert penalty_at(with_pressure, 1.0) > penalty_at(without, 1.0) * 10

    # hər iki halda həqiqi dəyər minimumdur
    for observations in (without, with_pressure):
        assert penalty_at(observations, TRUE_PERMEABILITY) < 1e-3
