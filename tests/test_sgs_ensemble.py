"""B8 №34-42 — SGS finalizasiyası: kondisiyalaşdırma, paylanma, variogram,
ansambl, seed təkrarlanması (GATE B10/B11).

`test_sgs.py` (Phase 5) tək realizasiyanın mövcud davranışını qoruyur —
bu fayl B5-in ƏLAVƏ tələblərini yoxlayır.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.geology.sgs import PropertyVariogramParams, simulate_sgs
from imex2d.geology.sgs_ensemble import (SGSEnsemble, ensemble_statistics,
                                         simulate_sgs_ensemble, validate_ensemble,
                                         validate_realization)
from imex2d.geology.variogram import VariogramParameters, fit_variogram_from_data


def _gaussian_field(points, ranges, seed):
    rng = np.random.default_rng(seed)
    n = points.shape[0]
    scaled = points[:, :len(ranges)] / np.asarray(ranges, float)[None, :]
    diff = scaled[:, None, :] - scaled[None, :, :]
    cov = np.exp(-3.0 * np.sqrt(np.sum(diff * diff, axis=-1))) + 1e-8 * np.eye(n)
    return np.linalg.cholesky(cov) @ rng.standard_normal(n)


def _grid(n=24, high=800.0):
    axis = np.linspace(0.0, high, n)
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


def _case(n_wells=30, seed=1, high=800.0, ranges=(220.0, 220.0)):
    """Sərt data + hədəflər (sərt data mövqeləri hədəflərə DAXİLDİR ki,
    kondisiyalaşdırma yoxlana bilsin)."""
    rng = np.random.default_rng(seed)
    points = rng.uniform(0.0, high, size=(n_wells, 2))
    values = 0.20 + 0.03 * _gaussian_field(points, ranges, seed + 1)
    targets = np.vstack([points, _grid(24, high)])
    return points, values, targets


VARIOGRAM = PropertyVariogramParams(model="spherical", nugget=0.0, range_=220.0)


# ── 34. sərt data kondisiyalaşdırması ─────────────────────────────────
def test_hard_data_is_reproduced_exactly_in_every_realization():
    """GATE B10 — yalnız qlobal histoqram uyğunluğu KİFAYƏT DEYİL."""
    points, values, targets = _case(seed=2)
    ensemble = simulate_sgs_ensemble(5, points, values, targets,
                                     variogram=VARIOGRAM, base_seed=11,
                                     max_neighbors=16)
    for realization in ensemble.realizations:
        assert np.allclose(realization.values[:len(values)], values, atol=1e-12)
        assert np.all(realization.hard_data_mask[:len(values)])


def test_validation_report_confirms_hard_data_conditioning():
    points, values, targets = _case(seed=3)
    realization = simulate_sgs(points, values, targets, variogram=VARIOGRAM,
                               seed=5, max_neighbors=16)
    report = validate_realization(realization, points, values, targets)
    assert report.hard_data_honored
    assert report.hard_data_max_error == pytest.approx(0.0, abs=1e-12)
    assert report.n_hard_data == len(values)


def test_hard_data_violation_would_be_detected():
    """Doğrulama funksiyası HƏQİQƏTƏN yoxlayır — süni pozuntu tutulur."""
    points, values, targets = _case(seed=4)
    realization = simulate_sgs(points, values, targets, variogram=VARIOGRAM,
                               seed=6, max_neighbors=16)
    realization.values[0] += 0.05                      # süni pozuntu
    report = validate_realization(realization, points, values, targets)
    assert not report.hard_data_honored
    assert report.hard_data_max_error == pytest.approx(0.05, rel=1e-6)
    assert any("POZULUB" in w for w in report.warnings)


# ── 35. histoqram (marjinal paylanma) ─────────────────────────────────
def test_realization_reproduces_the_conditioning_distribution():
    points, values, targets = _case(n_wells=40, seed=7)
    realization = simulate_sgs(points, values, targets, variogram=VARIOGRAM,
                               seed=8, max_neighbors=16)
    report = validate_realization(realization, points, values, targets)
    assert report.ks_statistic < 0.25, f"KS D={report.ks_statistic:.3f}"
    assert abs(report.mean_shift) < 0.8
    assert 0.6 < report.std_ratio < 1.6


def test_ensemble_mean_is_close_to_the_data_mean():
    points, values, targets = _case(n_wells=40, seed=9)
    ensemble = simulate_sgs_ensemble(8, points, values, targets,
                                     variogram=VARIOGRAM, base_seed=12,
                                     max_neighbors=16)
    assert float(np.mean(ensemble.mean)) == pytest.approx(float(np.mean(values)),
                                                          abs=0.02)


# ── 36. variogram təkrarlanması ───────────────────────────────────────
def test_realization_reproduces_the_target_spatial_continuity():
    points, values, targets = _case(n_wells=45, seed=10, ranges=(220.0, 220.0))
    realization = simulate_sgs(points, values, targets, variogram=VARIOGRAM,
                               seed=13, max_neighbors=20)
    target = VariogramParameters("spherical", 0.0, float(np.var(values)), 220.0,
                                 0.0, 0)
    report = validate_realization(realization, points, values, targets,
                                  target_variogram=target)
    assert np.isfinite(report.realized_range)
    assert np.isfinite(report.range_ratio)
    assert 0.6 <= report.range_ratio <= 1.7, (
        f"realizə/hədəf radius nisbəti {report.range_ratio:.2f}")


def test_shorter_target_range_produces_a_rougher_realization():
    """Variogram REAL təsir edir: kiçik radius → daha kobud sahə."""
    points, values, targets = _case(n_wells=35, seed=11)
    smooth = simulate_sgs(points, values, targets, seed=14, max_neighbors=16,
                          variogram=PropertyVariogramParams(model="spherical",
                                                            nugget=0.0, range_=500.0))
    rough = simulate_sgs(points, values, targets, seed=14, max_neighbors=16,
                         variogram=PropertyVariogramParams(model="spherical",
                                                           nugget=0.0, range_=60.0))
    grid = targets[len(values):]
    smooth_range = fit_variogram_from_data(grid, smooth.values[len(values):],
                                           model="spherical").range_
    rough_range = fit_variogram_from_data(grid, rough.values[len(values):],
                                          model="spherical").range_
    assert smooth_range > rough_range


# ── 37. anizotrop davamlılıq ──────────────────────────────────────────
def test_anisotropic_variogram_produces_directional_continuity():
    """major ox boyunca davamlılıq minor oxdan BÖYÜK olmalıdır."""
    from imex2d.geology.variogram import experimental_variogram
    points, values, targets = _case(n_wells=35, seed=12, high=1000.0)
    anisotropic = PropertyVariogramParams(model="spherical", nugget=0.0,
                                          range_=700.0, range_minor=90.0,
                                          azimuth_deg=90.0)
    realization = simulate_sgs(points, values, targets, variogram=anisotropic,
                               seed=15, max_neighbors=20)
    grid = targets[len(values):]
    simulated = realization.values[len(values):]

    along = experimental_variogram(grid, simulated, n_lags=8, max_lag=400.0,
                                   azimuth_deg=90.0, azimuth_tolerance_deg=20.0)
    across = experimental_variogram(grid, simulated, n_lags=8, max_lag=400.0,
                                    azimuth_deg=0.0, azimuth_tolerance_deg=20.0)
    mid = 4
    assert along.gamma[mid] < across.gamma[mid], (
        "major ox boyunca γ daha kiçik (davamlılıq daha çox) olmalıdır")


# ── 38. çoxlu realizasiya ─────────────────────────────────────────────
def test_ensemble_produces_the_requested_number_of_distinct_realizations():
    points, values, targets = _case(seed=13)
    ensemble = simulate_sgs_ensemble(6, points, values, targets,
                                     variogram=VARIOGRAM, base_seed=20,
                                     max_neighbors=16)
    assert isinstance(ensemble, SGSEnsemble)
    assert ensemble.n_realizations == 6
    assert ensemble.n_cells == targets.shape[0]
    assert [r.realization_id for r in ensemble.realizations] == list(range(6))
    simulated = ensemble.values[:, len(values):]
    for i in range(6):
        for j in range(i + 1, 6):
            assert not np.allclose(simulated[i], simulated[j]), (
                "realizasiyalar MÜSTƏQİL olmalıdır")


def test_ensemble_rejects_a_non_positive_realization_count():
    points, values, targets = _case(seed=14)
    with pytest.raises(ValueError, match="n_realizations"):
        simulate_sgs_ensemble(0, points, values, targets, variogram=VARIOGRAM)


def test_small_ensembles_warn_about_rough_quantiles():
    points, values, targets = _case(seed=15)
    ensemble = simulate_sgs_ensemble(3, points, values, targets,
                                     variogram=VARIOGRAM, base_seed=21,
                                     max_neighbors=16)
    assert any("KOBUDDUR" in w for w in ensemble.warnings)


# ── 39. seed təkrarlanması ────────────────────────────────────────────
def test_the_same_seed_reproduces_the_ensemble_bitwise():
    """GATE B11."""
    points, values, targets = _case(seed=16)
    kwargs = dict(variogram=VARIOGRAM, max_neighbors=16)
    a = simulate_sgs_ensemble(4, points, values, targets, base_seed=33, **kwargs)
    b = simulate_sgs_ensemble(4, points, values, targets, base_seed=33, **kwargs)
    assert np.array_equal(a.values, b.values)


def test_different_seeds_produce_different_ensembles():
    points, values, targets = _case(seed=17)
    kwargs = dict(variogram=VARIOGRAM, max_neighbors=16)
    a = simulate_sgs_ensemble(4, points, values, targets, base_seed=33, **kwargs)
    b = simulate_sgs_ensemble(4, points, values, targets, base_seed=77, **kwargs)
    assert not np.allclose(a.values[:, len(values):], b.values[:, len(values):])


def test_ensemble_seed_convention_matches_the_phase5_helper():
    """`run_realizations_sgs` ilə EYNİ seed konvensiyası — geriyə uyğunluq."""
    from imex2d.geology.sgs import run_realizations_sgs
    points, values, targets = _case(seed=18)
    kwargs = dict(variogram=VARIOGRAM, max_neighbors=16)
    ensemble = simulate_sgs_ensemble(3, points, values, targets, base_seed=5, **kwargs)
    legacy = run_realizations_sgs(3, points, values, targets, seed=5, **kwargs)
    for new, old in zip(ensemble.realizations, legacy):
        assert new.seed == old.seed
        assert np.array_equal(new.values, old.values)


# ── 40. ansambl statistikası ──────────────────────────────────────────
def test_ensemble_statistics_are_ordered_and_consistent():
    points, values, targets = _case(n_wells=30, seed=19)
    ensemble = simulate_sgs_ensemble(12, points, values, targets,
                                     variogram=VARIOGRAM, base_seed=44,
                                     max_neighbors=16)
    assert np.all(ensemble.p10 <= ensemble.p50 + 1e-12)
    assert np.all(ensemble.p50 <= ensemble.p90 + 1e-12)
    assert np.all(ensemble.variance >= -1e-15)
    assert np.allclose(ensemble.std ** 2, ensemble.variance, atol=1e-12)
    assert ensemble.mean.shape == (targets.shape[0],)


def test_ensemble_spread_is_zero_at_hard_data_and_positive_away_from_it():
    """Kondisiyalaşdırılmış nöqtədə bütün realizasiyalar EYNİDİR."""
    points, values, targets = _case(n_wells=25, seed=20)
    ensemble = simulate_sgs_ensemble(8, points, values, targets,
                                     variogram=VARIOGRAM, base_seed=45,
                                     max_neighbors=16)
    hard = ensemble.hard_data_mask
    assert np.allclose(ensemble.std[hard], 0.0, atol=1e-12)
    assert float(np.mean(ensemble.std[~hard])) > 0.0


def test_ensemble_statistics_helper_returns_named_quantiles():
    points, values, targets = _case(seed=21)
    ensemble = simulate_sgs_ensemble(5, points, values, targets,
                                     variogram=VARIOGRAM, base_seed=46,
                                     max_neighbors=16)
    stats = ensemble_statistics(ensemble, quantiles=(10.0, 50.0, 90.0))
    assert set(stats) == {"mean", "variance", "std", "p10", "p50", "p90"}
    assert np.allclose(stats["p50"], ensemble.p50)


def test_quantiles_are_never_called_confidence_intervals():
    """B5.7 — terminologiya qəsdən dəqiqdir."""
    import imex2d.geology.sgs_ensemble as module
    assert "etibar intervalı" in module.__doc__
    assert "confidence interval" not in module.SGSEnsemble.__doc__.lower()


def test_ensemble_exposes_grid_ready_arrays():
    points, values, targets = _case(seed=22)
    ensemble = simulate_sgs_ensemble(5, points, values, targets,
                                     variogram=VARIOGRAM, base_seed=47,
                                     max_neighbors=16)
    grids = ensemble.as_grids()
    assert set(grids) == {"mean", "variance", "std", "p10", "p50", "p90",
                          "hard_data"}
    for array in grids.values():
        assert array.shape == (targets.shape[0],)


def test_validate_ensemble_reports_every_realization():
    points, values, targets = _case(seed=23)
    ensemble = simulate_sgs_ensemble(4, points, values, targets,
                                     variogram=VARIOGRAM, base_seed=48,
                                     max_neighbors=16)
    reports = validate_ensemble(ensemble, points, values, targets)
    assert len(reports) == 4
    assert all(r.hard_data_honored for r in reports)
    assert all(isinstance(r.as_dict(), dict) for r in reports)
    assert all("SGS doğrulaması" in r.as_text() for r in reports)


# ── 41. sabit sahə ────────────────────────────────────────────────────
def test_constant_property_is_propagated_without_invented_variability():
    points, values, targets = _case(seed=24)
    constant = np.full(values.size, 0.21)
    ensemble = simulate_sgs_ensemble(4, points, constant, targets,
                                     variogram=VARIOGRAM, base_seed=49,
                                     max_neighbors=16)
    assert np.allclose(ensemble.values, 0.21, atol=1e-12)
    assert np.allclose(ensemble.std, 0.0, atol=1e-12)
    assert any("SABİTDİR" in w for w in ensemble.warnings)


# ── 42. seyrək məlumat ────────────────────────────────────────────────
def test_sparse_data_still_produces_finite_conditioned_realizations():
    rng = np.random.default_rng(25)
    points = rng.uniform(0.0, 800.0, size=(4, 2))
    values = np.array([0.18, 0.22, 0.20, 0.24])
    targets = np.vstack([points, _grid(12, 800.0)])
    ensemble = simulate_sgs_ensemble(3, points, values, targets, base_seed=50,
                                     max_neighbors=8)
    assert np.all(np.isfinite(ensemble.values))
    assert np.allclose(ensemble.values[:, :4], values, atol=1e-12)


def test_single_hard_data_point_is_handled_without_crashing():
    points = np.array([[400.0, 400.0]])
    values = np.array([0.2])
    targets = np.vstack([points, _grid(10, 800.0)])
    ensemble = simulate_sgs_ensemble(2, points, values, targets, base_seed=51,
                                     max_neighbors=4)
    assert np.all(np.isfinite(ensemble.values))
    assert np.allclose(ensemble.values[:, 0], 0.2)


def test_sparse_data_produces_larger_ensemble_spread_than_dense_data():
    """B9 Case 7 — məlumat seyrəkləşdikcə qeyri-müəyyənlik ARTIR."""
    rng = np.random.default_rng(26)
    dense_points = rng.uniform(0.0, 800.0, size=(60, 2))
    dense_values = 0.20 + 0.03 * _gaussian_field(dense_points, (220.0, 220.0), 27)
    sparse_points = dense_points[:6]
    sparse_values = dense_values[:6]
    grid = _grid(16, 800.0)

    dense = simulate_sgs_ensemble(6, dense_points, dense_values,
                                  np.vstack([dense_points, grid]),
                                  variogram=VARIOGRAM, base_seed=60,
                                  max_neighbors=16)
    sparse = simulate_sgs_ensemble(6, sparse_points, sparse_values,
                                   np.vstack([sparse_points, grid]),
                                   variogram=VARIOGRAM, base_seed=60,
                                   max_neighbors=16)
    dense_spread = float(np.mean(dense.std[len(dense_values):]))
    sparse_spread = float(np.mean(sparse.std[len(sparse_values):]))
    assert sparse_spread > dense_spread
