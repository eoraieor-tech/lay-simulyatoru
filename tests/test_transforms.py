"""B8 №8-12 — xassə-fəzası çevirmələri: loq, logit, normal-score, bias
düzəlişi, dəqiq 0/1 halları, sabit sahə.

Bu testlərin mərkəzi iddiası: çevirmələr RİYAZİ olaraq düzgündür və
`inverse(forward(z)) == z` MAŞIN DƏQİQLİYİ ilə ödənir — yəni sərt
datanın dəqiq honor edilməsi çevirmə səbəbindən POZULMUR.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.geology.transforms import (IDENTITY_TRANSFORM, BackTransform,
                                       LogitTransform, LogTransform,
                                       NormalScoreValueTransform, TransformError,
                                       ValueTransform, VarianceKind,
                                       apply_back_transform)


# ── 8. loq çevirməsi ──────────────────────────────────────────────────
def test_log_transform_round_trips_exactly():
    values = np.array([1e-3, 0.5, 1.0, 42.0, 5_000.0, 1e6])
    transform = LogTransform()
    assert np.allclose(transform.inverse(transform.forward(values)), values,
                       rtol=1e-12, atol=0.0)


def test_log_transform_is_the_natural_logarithm():
    values = np.array([1.0, np.e, 100.0])
    assert np.allclose(LogTransform().forward(values), np.log(values), atol=0.0)


def test_log_transform_offset_allows_zero_and_round_trips():
    values = np.array([0.0, 1.0, 10.0])
    transform = LogTransform(offset=1.0)
    assert np.all(np.isfinite(transform.forward(values)))
    assert np.allclose(transform.inverse(transform.forward(values)), values, atol=1e-12)


@pytest.mark.parametrize("bad", [
    np.array([1.0, 0.0, 2.0]),        # sıfır
    np.array([1.0, -5.0, 2.0]),       # mənfi
])
def test_log_transform_rejects_non_positive_instead_of_inventing_values(bad):
    """Sıfır/mənfi keçiricilik SƏSSİZCƏ kiçik müsbət ədədə çevrilmir (B1.2)."""
    with pytest.raises(TransformError, match="müsbət"):
        LogTransform().forward(bad)


def test_log_transform_rejects_non_finite():
    with pytest.raises(TransformError):
        LogTransform().forward(np.array([1.0, np.nan]))
    with pytest.raises(TransformError):
        LogTransform().forward(np.array([1.0, np.inf]))


# ── 9. loq-normal bias düzəlişi ───────────────────────────────────────
def test_lognormal_back_transform_modes_are_mathematically_distinct():
    """`median`, `mean`, `mean_ok` ÜÇ FƏRQLİ kəmiyyətdir (B1.3)."""
    transform = LogTransform()
    y = np.array([np.log(100.0)])
    sigma2 = np.array([0.5])
    mu = np.array([0.1])

    median = transform.inverse(y, sigma2, mode=BackTransform.MEDIAN)
    mean = transform.inverse(y, sigma2, mode=BackTransform.MEAN)
    mean_ok = transform.inverse(y, sigma2, mu, mode=BackTransform.MEAN_OK)

    assert median[0] == pytest.approx(100.0)
    assert mean[0] == pytest.approx(100.0 * np.exp(0.25))
    assert mean_ok[0] == pytest.approx(100.0 * np.exp(0.25 - 0.1))
    assert median[0] < mean_ok[0] < mean[0]


def test_lognormal_mean_exceeds_median_by_the_right_factor():
    """`E[K]/median(K) = exp(σ²/2)` — loq-normal paylanmanın tərifi."""
    transform = LogTransform()
    y = np.array([2.0, 2.0, 2.0])
    for sigma2 in (0.1, 1.0, 3.0):
        variance = np.full(3, sigma2)
        ratio = (transform.inverse(y, variance, mode=BackTransform.MEAN)
                 / transform.inverse(y, variance, mode=BackTransform.MEDIAN))
        assert np.allclose(ratio, np.exp(sigma2 / 2.0))


def test_mean_ok_requires_the_lagrange_multiplier():
    with pytest.raises(TransformError, match="MEAN_OK"):
        LogTransform().inverse(np.array([1.0]), np.array([0.5]),
                               mode=BackTransform.MEAN_OK)


def test_lognormal_variance_matches_the_closed_form():
    """`Var[K] = exp(2ŷ + σ²)(exp(σ²) − 1)` — yaxınlaşma DEYİL."""
    transform = LogTransform()
    y = np.array([np.log(50.0), np.log(500.0)])
    sigma2 = np.array([0.2, 0.8])
    variance, kind = transform.inverse_variance(y, sigma2)
    expected = np.exp(2 * y + sigma2) * (np.exp(sigma2) - 1.0)
    assert kind is VarianceKind.EXACT
    assert np.allclose(variance, expected, rtol=1e-10)


def test_lognormal_variance_matches_monte_carlo():
    """Qapalı düsturu MÜSTƏQİL yolla (Monte-Karlo) yoxlayırıq."""
    rng = np.random.default_rng(7)
    y_hat, sigma2 = 1.5, 0.4
    sample = np.exp(rng.normal(y_hat, np.sqrt(sigma2), size=400_000))
    analytic, _ = LogTransform().inverse_variance(np.array([y_hat]),
                                                  np.array([sigma2]))
    assert analytic[0] == pytest.approx(float(np.var(sample)), rel=0.03)


def test_zero_variance_makes_all_back_transform_modes_agree():
    transform = LogTransform()
    y = np.array([1.234])
    zero = np.array([0.0])
    assert transform.inverse(y, zero, zero, mode=BackTransform.MEAN)[0] == \
        pytest.approx(transform.inverse(y, zero, mode=BackTransform.MEDIAN)[0])
    assert transform.inverse(y, zero, zero, mode=BackTransform.MEAN_OK)[0] == \
        pytest.approx(np.exp(1.234))


# ── 10. logit çevirməsi ───────────────────────────────────────────────
def test_logit_round_trips_exactly_including_the_bounds():
    """`eps` sıxması GERİ AÇILIR — hədlərdə də dəqiq bərpa (B1.4)."""
    transform = LogitTransform(lower=0.0, upper=1.0, eps=1e-4)
    values = np.array([0.0, 1e-9, 0.01, 0.5, 0.73, 0.999, 1.0])
    assert np.allclose(transform.inverse(transform.forward(values)), values,
                       atol=1e-12)


def test_logit_inverse_can_never_leave_the_bounds():
    """RİYAZİ zəmanət: hansı kriginq nəticəsi gəlirsə gəlsin, geri
    çevirmə `[0,1]`-dən çıxa BİLMİR."""
    transform = LogitTransform()
    extreme = np.array([-1e6, -50.0, -3.0, 0.0, 3.0, 50.0, 1e6])
    result = transform.inverse(extreme)
    assert np.all(result >= 0.0) and np.all(result <= 1.0)
    assert np.all(np.isfinite(result))


def test_logit_is_monotone_increasing():
    transform = LogitTransform()
    values = np.linspace(0.0, 1.0, 50)
    assert np.all(np.diff(transform.forward(values)) > 0.0)


def test_logit_maps_the_midpoint_to_zero():
    assert LogitTransform().forward(np.array([0.5]))[0] == pytest.approx(0.0)


def test_logit_supports_arbitrary_bounds():
    transform = LogitTransform(lower=0.2, upper=0.8)
    values = np.array([0.2, 0.35, 0.5, 0.8])
    assert np.allclose(transform.inverse(transform.forward(values)), values, atol=1e-12)
    assert np.all(transform.inverse(np.array([-99.0, 99.0])) >= 0.2 - 1e-12)
    assert np.all(transform.inverse(np.array([-99.0, 99.0])) <= 0.8 + 1e-12)


def test_logit_rejects_values_outside_the_bounds():
    with pytest.raises(TransformError, match="aralığında"):
        LogitTransform().forward(np.array([0.5, 1.4]))
    with pytest.raises(TransformError):
        LogitTransform().forward(np.array([-0.1, 0.5]))


@pytest.mark.parametrize("kwargs", [
    {"lower": 1.0, "upper": 0.0},
    {"eps": 0.0},
    {"eps": 0.5},
    {"lower": np.nan},
])
def test_invalid_logit_configuration_raises(kwargs):
    base = {"lower": 0.0, "upper": 1.0, "eps": 1e-4}
    base.update(kwargs)
    with pytest.raises(TransformError):
        LogitTransform(**base)


def test_logit_variance_uses_the_documented_delta_method():
    """`Var[z] ≈ (dz/dy)²·σ²_y`, `dz/dy = span·(1−2ε)·p̃(1−p̃)`."""
    transform = LogitTransform()
    y = np.array([0.0, 2.0])
    sigma2 = np.array([0.25, 0.25])
    variance, kind = transform.inverse_variance(y, sigma2)
    assert kind is VarianceKind.DELTA
    squeezed = 1.0 / (1.0 + np.exp(-y))
    derivative = 1.0 * (1.0 - 2.0 * transform.eps) * squeezed * (1.0 - squeezed)
    assert np.allclose(variance, derivative ** 2 * sigma2)
    # mərkəzdə (p=0.5) qeyri-müəyyənlik ƏN BÖYÜK, kənarda kiçilir
    assert variance[0] > variance[1]


def test_logit_delta_variance_approximates_monte_carlo_for_small_sigma():
    rng = np.random.default_rng(11)
    transform = LogitTransform()
    y_hat, sigma2 = 0.4, 0.05
    sample = transform.inverse(rng.normal(y_hat, np.sqrt(sigma2), size=300_000))
    analytic, _ = transform.inverse_variance(np.array([y_hat]), np.array([sigma2]))
    assert analytic[0] == pytest.approx(float(np.var(sample)), rel=0.05)


# ── 11. dəqiq 0/1 halları ─────────────────────────────────────────────
def test_exact_zero_and_one_become_finite_and_symmetric():
    transform = LogitTransform(eps=1e-4)
    at_zero = float(transform.forward(np.array([0.0]))[0])
    at_one = float(transform.forward(np.array([1.0]))[0])
    assert np.isfinite(at_zero) and np.isfinite(at_one)
    assert at_zero == pytest.approx(np.log(1e-4 / (1 - 1e-4)))
    assert at_one == pytest.approx(-at_zero)


@pytest.mark.parametrize("eps", [1e-2, 1e-3, 1e-4, 1e-6])
def test_eps_controls_only_how_far_the_bounds_reach_not_the_round_trip(eps):
    """`eps` sənədləşdirilmiş TƏSİRƏ malikdir: yalnız hədlərin logit
    fəzasında nə qədər uzağa düşdüyünü dəyişir; geri çevirmə HƏR `eps`
    üçün dəqiq qalır."""
    transform = LogitTransform(eps=eps)
    values = np.array([0.0, 0.25, 0.5, 1.0])
    assert np.allclose(transform.inverse(transform.forward(values)), values, atol=1e-12)
    extreme = abs(float(transform.forward(np.array([0.0]))[0]))
    assert extreme == pytest.approx(abs(np.log(eps / (1 - eps))))


def test_smaller_eps_pushes_the_bounds_further_out():
    coarse = abs(float(LogitTransform(eps=1e-2).forward(np.array([0.0]))[0]))
    fine = abs(float(LogitTransform(eps=1e-6).forward(np.array([0.0]))[0]))
    assert fine > coarse


# ── 12. sabit sahə ────────────────────────────────────────────────────
def test_constant_field_survives_every_transform():
    constant = np.full(20, 0.42)
    for transform in (IDENTITY_TRANSFORM, LogTransform(),
                      LogitTransform(), NormalScoreValueTransform().fit(constant)):
        forward = transform.forward(constant)
        assert np.all(np.isfinite(forward))
        assert np.allclose(transform.inverse(forward), constant, atol=1e-9)


def test_normal_score_of_a_constant_field_is_degenerate_but_safe():
    constant = np.full(15, 3.5)
    transform = NormalScoreValueTransform().fit(constant)
    assert np.allclose(transform.forward(constant), 0.0)
    assert np.allclose(transform.inverse(np.array([-2.0, 0.0, 2.0])), 3.5)
    variance, kind = transform.inverse_variance(np.array([0.0]), np.array([1.0]))
    assert variance[0] == 0.0 and kind is VarianceKind.EXACT


def test_constant_permeability_round_trips_in_log_space():
    constant = np.full(10, 250.0)
    transform = LogTransform()
    assert np.allclose(transform.inverse(transform.forward(constant)), constant)


# ── normal-score ──────────────────────────────────────────────────────
def test_normal_score_is_data_dependent_and_requires_fit():
    transform = NormalScoreValueTransform()
    assert transform.data_dependent is True
    with pytest.raises(TransformError, match="fit"):
        transform.forward(np.array([1.0]))


def test_normal_score_round_trips_on_its_own_data():
    rng = np.random.default_rng(5)
    values = np.exp(rng.normal(2.0, 1.0, size=200))
    transform = NormalScoreValueTransform().fit(values)
    assert np.allclose(transform.inverse(transform.forward(values)), values, atol=1e-9)


def test_normal_score_output_is_approximately_standard_normal():
    rng = np.random.default_rng(6)
    values = np.exp(rng.normal(0.0, 2.0, size=2000))
    scores = NormalScoreValueTransform().fit(values).forward(values)
    assert abs(float(np.mean(scores))) < 0.05
    assert float(np.std(scores)) == pytest.approx(1.0, rel=0.05)


def test_normal_score_fit_returns_a_new_object_leaving_the_prototype_unfitted():
    """Sızma auditinin əsası: `fit()` YENİ obyekt qaytarır, prototip
    dəyişmir — qatlar arasında statistika sızmır (B3.1)."""
    prototype = NormalScoreValueTransform()
    fitted = prototype.fit(np.array([1.0, 2.0, 3.0, 4.0]))
    assert fitted is not prototype
    assert prototype.table is None and fitted.table is not None


# ── eynilik + köməkçi ─────────────────────────────────────────────────
def test_identity_transform_changes_nothing():
    values = np.array([-3.0, 0.0, 7.5])
    assert IDENTITY_TRANSFORM.is_identity
    assert np.array_equal(IDENTITY_TRANSFORM.forward(values), values)
    assert np.array_equal(IDENTITY_TRANSFORM.inverse(values), values)
    variance, kind = IDENTITY_TRANSFORM.inverse_variance(values, np.ones(3))
    assert kind is VarianceKind.IDENTITY and np.array_equal(variance, np.ones(3))


def test_base_class_is_usable_directly_for_backward_compatibility():
    """Phase A `ValueTransform()`-u eynilik kimi işlədirdi — davam edir."""
    transform = ValueTransform()
    values = np.array([1.0, 2.0])
    assert np.array_equal(transform.inverse(transform.forward(values)), values)


def test_apply_back_transform_reports_the_variance_kind():
    estimate = np.array([np.log(10.0)])
    variance = np.array([0.3])
    values, back_variance, kind = apply_back_transform(
        LogTransform(), estimate, variance, np.zeros(1), BackTransform.MEAN)
    assert kind is VarianceKind.EXACT
    assert values[0] == pytest.approx(10.0 * np.exp(0.15))
    assert back_variance[0] > 0.0


def test_apply_back_transform_without_variance_marks_it_undefined():
    values, variance, kind = apply_back_transform(
        LogTransform(), np.array([0.0]), None, None, BackTransform.MEDIAN)
    assert kind is VarianceKind.UNDEFINED
    assert np.isnan(variance[0]) and values[0] == pytest.approx(1.0)
