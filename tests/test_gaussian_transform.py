"""Phase 5 §3 — normal-score (Gauss) çevirməsi."""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.geology.gaussian_transform import NormalScoreTransform


def test_forward_then_inverse_reproduces_original_values_at_data_points():
    rng = np.random.default_rng(0)
    values = rng.lognormal(mean=2.0, sigma=0.5, size=50)
    transform = NormalScoreTransform.fit(values)
    gaussian = transform.forward(values)
    back = transform.inverse(gaussian)
    assert np.allclose(back, values, atol=1e-9)


def test_forward_transform_is_approximately_standard_normal_distributed():
    rng = np.random.default_rng(1)
    values = rng.lognormal(mean=1.0, sigma=1.0, size=500)
    transform = NormalScoreTransform.fit(values)
    gaussian = transform.forward(values)
    assert abs(np.mean(gaussian)) < 0.1
    assert abs(np.std(gaussian) - 1.0) < 0.1


def test_reproducible_deterministic_fit():
    values = np.array([1.0, 5.0, 2.0, 8.0, 3.0, 7.0])
    t1 = NormalScoreTransform.fit(values)
    t2 = NormalScoreTransform.fit(values)
    assert np.array_equal(t1.sorted_gaussian, t2.sorted_gaussian)


def test_tied_values_get_identical_gaussian_score():
    values = np.array([1.0, 2.0, 2.0, 2.0, 5.0])
    transform = NormalScoreTransform.fit(values)
    scores = transform.forward(np.array([2.0]))
    tied_scores = transform.forward(np.array([2.0, 2.0, 2.0]))
    assert np.allclose(tied_scores, scores[0])


def test_constant_property_is_handled_without_fabricating_variability():
    values = np.full(10, 0.22)
    transform = NormalScoreTransform.fit(values)
    assert transform.is_constant
    gaussian = transform.forward(values)
    assert np.all(gaussian == 0.0)
    back = transform.inverse(np.array([-2.0, 0.0, 3.0]))
    assert np.all(back == 0.22)   # sabit dəyər QORUNUR, UYDURULMUŞ dəyişkənlik YOXDUR


def test_small_dataset_does_not_crash():
    for n in (1, 2, 3):
        values = np.linspace(0.1, 0.3, n)
        transform = NormalScoreTransform.fit(values)
        back = transform.inverse(transform.forward(values))
        assert np.allclose(back, values, atol=1e-9)


def test_extreme_tail_queries_clamp_to_boundary_not_extrapolated():
    values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    transform = NormalScoreTransform.fit(values)
    gaussian = transform.forward(values)
    far_low = transform.inverse(np.array([gaussian.min() - 10.0]))
    far_high = transform.inverse(np.array([gaussian.max() + 10.0]))
    assert far_low[0] == values.min()
    assert far_high[0] == values.max()


def test_rejects_empty_and_nan_input():
    with pytest.raises(ValueError):
        NormalScoreTransform.fit(np.array([]))
    with pytest.raises(ValueError):
        NormalScoreTransform.fit(np.array([1.0, np.nan, 2.0]))
