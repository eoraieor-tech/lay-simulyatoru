"""Phase 5 §4 — xassə paylanması təhlili (SGS-dən əvvəl)."""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.geology.distribution_analysis import (log_transform_is_justified,
                                                   summarize_distribution)


def test_summary_matches_numpy_reference():
    rng = np.random.default_rng(0)
    values = rng.normal(0.2, 0.03, size=200)
    summary = summarize_distribution(values)
    assert summary.n == 200
    assert summary.mean == pytest.approx(float(np.mean(values)))
    assert summary.std == pytest.approx(float(np.std(values)))
    assert summary.p50 == pytest.approx(float(np.median(values)))
    assert summary.minimum == pytest.approx(float(values.min()))
    assert summary.maximum == pytest.approx(float(values.max()))


def test_summary_rejects_empty_and_nan():
    with pytest.raises(ValueError):
        summarize_distribution([])
    with pytest.raises(ValueError):
        summarize_distribution([1.0, np.nan])


def test_skewness_is_none_for_constant_or_tiny_sample():
    assert summarize_distribution([1.0, 1.0, 1.0]).skewness is None
    assert summarize_distribution([1.0, 2.0]).skewness is None


# ── log-fəza uyğunluğu: DATA-ƏSASLI, kor-koranə DEYİL ──────────────────
def test_log_transform_is_justified_for_lognormal_permeability():
    rng = np.random.default_rng(1)
    permx = rng.lognormal(mean=4.0, sigma=1.2, size=300)   # tipik keçiricilik forması
    assert log_transform_is_justified(permx)


def test_log_transform_not_justified_for_already_symmetric_data():
    rng = np.random.default_rng(2)
    poro = rng.normal(0.20, 0.03, size=300)   # simmetrik, artıq az çarpıq
    assert not log_transform_is_justified(poro)


def test_log_transform_rejects_non_positive_values():
    values = np.array([-1.0, 2.0, 3.0, 4.0])
    assert not log_transform_is_justified(values)


def test_log_transform_false_for_constant_values():
    assert not log_transform_is_justified(np.full(10, 150.0))
