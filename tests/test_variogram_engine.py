"""A3 — peşəkar variogram özəyi: deneysel γ(h), lag tolerantlığı, minimal
cüt sayı, istiqamətli/şaquli variogramlar, model ailəsi, fit, model
müqayisəsi və PARAMETR DOĞRULAMASI.

`test_variogram.py` (Phase 3) mövcud davranışı qoruyur — burada YENİ
(A3) imkanlar yoxlanılır.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.geology.variogram import (EFFECTIVE_RANGE_FACTOR, KNOWN_MODELS,
                                      MODEL_EXPONENTIAL, MODEL_GAUSSIAN,
                                      MODEL_SPHERICAL, VariogramValidationError,
                                      compare_variogram_models, covariance,
                                      directional_variograms, effective_range,
                                      experimental_variogram, fit_variogram,
                                      fit_variogram_from_data, gaussian,
                                      repair_variogram_parameters, semivariance,
                                      spherical, validate_variogram_parameters,
                                      vertical_variogram)


def _gaussian_random_field(n=260, seed=0, range_x=250.0, range_y=250.0, high=1000.0):
    """Verilmiş korrelyasiya radiusları ilə (anizotrop ola bilən) sahə.

    Sahə kovariasiya matrisinin Xolesskiy parçalanması ilə qurulur —
    "təxminən korrelyasiyalı görünən" səs-küy DEYİL, radiusları MƏLUM
    olan həqiqi qauss sahəsi."""
    rng = np.random.default_rng(seed)
    points = rng.uniform(0.0, high, size=(n, 2))
    dx = (points[:, 0][:, None] - points[:, 0][None, :]) / range_x
    dy = (points[:, 1][:, None] - points[:, 1][None, :]) / range_y
    cov = np.exp(-3.0 * np.sqrt(dx ** 2 + dy ** 2)) + 1e-8 * np.eye(n)
    values = np.linalg.cholesky(cov) @ rng.standard_normal(n)
    return points, values


# ── 21. deneysel variogram ────────────────────────────────────────────
def test_experimental_variogram_matches_manual_definition():
    """γ(h) = 1/(2N) Σ (Zᵢ − Zⱼ)² — bir binin dəyəri ƏLLƏ yoxlanılır."""
    rng = np.random.default_rng(4)
    points = rng.uniform(0.0, 400.0, size=(50, 2))
    values = rng.standard_normal(50)
    exp = experimental_variogram(points, values, n_lags=8, max_lag=300.0)

    iu, ju = np.triu_indices(50, k=1)
    h = np.linalg.norm(points[iu] - points[ju], axis=1)
    dv = values[iu] - values[ju]
    width = 300.0 / 8
    # `max_lag`-dan UZAQ cütlər tamamilə ATILIR (sonuncu binə yığılmır)
    used = h <= 300.0
    bin_index = np.minimum((h / width).astype(int), 7)
    bin_index = np.where(used, bin_index, -1)
    assert int(exp.counts.sum()) == int(used.sum()) == exp.n_pairs_total
    assert exp.n_pairs_total < h.size, "uzaq cütlər həqiqətən atılmalıdır"
    for b in range(8):
        in_bin = bin_index == b
        if exp.counts[b] == 0:
            continue
        assert exp.counts[b] == int(in_bin.sum())
        assert exp.gamma[b] == pytest.approx(0.5 * np.mean(dv[in_bin] ** 2), rel=1e-9)


def test_experimental_variogram_rises_towards_the_sill():
    points, values = _gaussian_random_field(n=200, seed=1, range_x=200.0, range_y=200.0)
    exp = experimental_variogram(points, values, n_lags=12, max_lag=600.0).valid()
    assert exp.gamma[0] < exp.gamma[-1]
    assert exp.gamma[-1] == pytest.approx(np.var(values), rel=0.6)


def test_pair_counts_are_conserved_in_contiguous_bins():
    points, values = _gaussian_random_field(n=80, seed=2)
    exp = experimental_variogram(points, values, n_lags=10, max_lag=1e9)
    assert int(exp.counts.sum()) == exp.n_pairs_total == 80 * 79 // 2


# ── 22. lag tolerantlığı ──────────────────────────────────────────────
def test_lag_tolerance_overlapping_bins_share_pairs():
    """`lag_tolerance > lag_width/2` — binlər ÜST-ÜSTƏ düşür, ona görə
    cütlərin ümumi sayı ümumi cüt sayından ÇOX olur."""
    points, values = _gaussian_random_field(n=90, seed=3)
    width_only = experimental_variogram(points, values, n_lags=10, max_lag=800.0)
    overlapping = experimental_variogram(points, values, n_lags=10, max_lag=800.0,
                                         lag_tolerance=80.0)   # yarım-en 80 > 40
    assert overlapping.counts.sum() > width_only.counts.sum()
    assert overlapping.lag_tolerance == 80.0
    assert np.allclose(overlapping.lags, width_only.lags)


def test_narrow_lag_tolerance_leaves_gaps_and_drops_pairs():
    points, values = _gaussian_random_field(n=90, seed=3)
    baseline = experimental_variogram(points, values, n_lags=10, max_lag=800.0)
    narrow = experimental_variogram(points, values, n_lags=10, max_lag=800.0,
                                    lag_tolerance=5.0)
    assert narrow.counts.sum() < baseline.counts.sum()


def test_zero_lag_tolerance_is_rejected():
    points, values = _gaussian_random_field(n=40, seed=3)
    with pytest.raises(ValueError, match="lag_tolerance"):
        experimental_variogram(points, values, lag_tolerance=0.0)


# ── 23. minimal cüt sayı ──────────────────────────────────────────────
def test_min_pairs_removes_statistically_meaningless_bins():
    points, values = _gaussian_random_field(n=40, seed=5, high=500.0)
    loose = experimental_variogram(points, values, n_lags=40, max_lag=700.0)
    strict = experimental_variogram(points, values, n_lags=40, max_lag=700.0,
                                    min_pairs=25)
    assert strict.n_bins < loose.valid().n_bins
    assert np.all(strict.counts >= 25)
    assert any("min_pairs" in w for w in strict.warnings)


def test_valid_min_pairs_filter_is_explicit():
    points, values = _gaussian_random_field(n=60, seed=6)
    exp = experimental_variogram(points, values, n_lags=20, max_lag=900.0)
    assert exp.valid(min_pairs=1).n_bins >= exp.valid(min_pairs=50).n_bins
    assert np.all(exp.valid(min_pairs=50).counts >= 50)


def test_fit_raises_when_min_pairs_leaves_too_few_bins():
    points, values = _gaussian_random_field(n=30, seed=7)
    exp = experimental_variogram(points, values, n_lags=12, max_lag=900.0)
    with pytest.raises(ValueError, match="dolu lag-bin"):
        fit_variogram(exp, model=MODEL_SPHERICAL, min_pairs=10_000)


# ── 24-26. model ailəsi ───────────────────────────────────────────────
@pytest.mark.parametrize("model,func", [
    (MODEL_SPHERICAL, spherical),
    (MODEL_EXPONENTIAL, None),
    (MODEL_GAUSSIAN, gaussian),
])
def test_models_are_monotone_start_at_nugget_and_reach_practical_sill(model, func):
    h = np.linspace(0.0, 400.0, 401)
    nugget, sill, range_ = 0.2, 1.0, 150.0
    gamma = semivariance(h, nugget, sill, range_, model)
    assert gamma[0] == pytest.approx(nugget)
    assert np.all(np.diff(gamma) >= -1e-12), "γ(h) monoton artan olmalıdır"
    # praktiki radius konvensiyası: γ(range_) ≈ 0.95·(nugget+sill)
    at_range = float(semivariance(np.array([range_]), nugget, sill, range_, model)[0])
    assert at_range >= 0.94 * (nugget + sill)
    assert gamma[-1] <= nugget + sill + 1e-12


def test_spherical_reaches_the_sill_exactly_at_the_range():
    value = float(semivariance(np.array([150.0]), 0.0, 1.0, 150.0, MODEL_SPHERICAL)[0])
    assert value == pytest.approx(1.0)
    beyond = float(semivariance(np.array([400.0]), 0.0, 1.0, 150.0, MODEL_SPHERICAL)[0])
    assert beyond == pytest.approx(1.0)


def test_gaussian_is_flatter_than_exponential_near_the_origin():
    """Qauss modeli başlanğıcda parabolik (çox hamar), eksponensial isə
    xəttidir — bu, iki modeli fərqləndirən əsas əlamətdir."""
    h = np.array([5.0])
    g = float(semivariance(h, 0.0, 1.0, 200.0, MODEL_GAUSSIAN)[0])
    e = float(semivariance(h, 0.0, 1.0, 200.0, MODEL_EXPONENTIAL)[0])
    assert g < e


def test_covariance_and_semivariance_are_consistent():
    """`C(h) = C(0) − γ(h)` (h > 0), `C(0) = nugget + sill` (A1.6).

    Başlanğıcdakı nugget SIÇRAYIŞI qəsdən ayrıca yoxlanılır: `C(0)`
    tam dispersiyadır, `C(0⁺)` isə ondan nugget qədər aşağıdır."""
    h = np.linspace(0.0, 500.0, 51)
    for model in KNOWN_MODELS:
        gamma = semivariance(h, 0.1, 0.9, 200.0, model)
        cov = covariance(h, 0.1, 0.9, 200.0, model)
        assert np.allclose(gamma[1:] + cov[1:], 1.0, atol=1e-12)
        assert cov[0] == pytest.approx(1.0)                      # C(0) = nugget+sill
        assert float(covariance(np.array([1e-6]), 0.1, 0.9, 200.0, model)[0])             == pytest.approx(0.9, abs=1e-6)                      # C(0+) = sill
        assert cov[-1] <= cov[0]
        assert cov[-1] == pytest.approx(0.0, abs=0.02)


def test_effective_range_is_documented_and_distinct_from_practical_range():
    assert effective_range(300.0, MODEL_SPHERICAL) == pytest.approx(300.0)
    assert effective_range(300.0, MODEL_EXPONENTIAL) == pytest.approx(100.0)
    assert effective_range(300.0, MODEL_GAUSSIAN) == pytest.approx(300.0 / np.sqrt(3.0))
    assert set(EFFECTIVE_RANGE_FACTOR) == set(KNOWN_MODELS)


def test_unknown_model_name_is_rejected_everywhere():
    for call in (lambda: semivariance(np.array([1.0]), 0.0, 1.0, 1.0, "linear"),
                 lambda: covariance(np.array([1.0]), 0.0, 1.0, 1.0, "linear"),
                 lambda: effective_range(1.0, "linear")):
        with pytest.raises(ValueError, match="Naməlum variogram modeli"):
            call()


# ── 27. parametr fitting ──────────────────────────────────────────────
@pytest.mark.parametrize("model", KNOWN_MODELS)
def test_fit_recovers_synthetic_parameters(model):
    """Modelin ÖZÜNDƏN qurulmuş "deneysel" nöqtələr → fit onları
    geri tapmalıdır."""
    lags = np.linspace(20.0, 400.0, 20)
    truth = semivariance(lags, 0.15, 1.0, 220.0, model)
    from imex2d.geology.variogram import ExperimentalVariogram
    exp = ExperimentalVariogram(lags, truth, np.full(20, 100), lag_width=20.0)
    fit = fit_variogram(exp, model=model)
    assert fit.model == model
    assert fit.range_ == pytest.approx(220.0, rel=0.15)
    assert fit.sill == pytest.approx(1.0, rel=0.25)
    assert fit.nugget == pytest.approx(0.15, abs=0.12)
    assert fit.weighted_rmse < 0.02


def test_fit_from_real_field_gives_a_sensible_range():
    points, values = _gaussian_random_field(n=260, seed=11, range_x=250.0,
                                            range_y=250.0, high=1200.0)
    fit = fit_variogram_from_data(points, values, n_lags=14, max_lag=700.0,
                                  model=MODEL_EXPONENTIAL)
    assert 80.0 < fit.range_ < 900.0
    assert fit.sill > 0.0 and fit.nugget >= 0.0
    assert fit.n_pairs > 0 and fit.n_lags_used >= 3


def test_fit_records_the_lag_and_direction_configuration():
    points, values = _gaussian_random_field(n=120, seed=12)
    exp = experimental_variogram(points, values, n_lags=10, max_lag=600.0,
                                 azimuth_deg=45.0, lag_tolerance=45.0)
    fit = fit_variogram(exp, model=MODEL_SPHERICAL)
    assert fit.azimuth_deg == 45.0
    assert fit.lag_tolerance == 45.0
    assert fit.lag_width == pytest.approx(60.0)
    assert fit.vertical is False


def test_weighted_fitting_favours_well_populated_lags():
    """Çəki `sqrt(cüt sayı)`-dır: bir binin cüt sayı 1000 dəfə artanda
    fit ona daha çox "qulaq asmalıdır"."""
    from imex2d.geology.variogram import ExperimentalVariogram
    lags = np.linspace(20.0, 400.0, 12)
    gamma = semivariance(lags, 0.0, 1.0, 200.0, MODEL_SPHERICAL)
    gamma[-1] = 3.0                                    # kənar (outlier) bin
    light = np.full(12, 100); light[-1] = 100
    heavy = np.full(12, 100); heavy[-1] = 100_000

    fit_light = fit_variogram(ExperimentalVariogram(lags, gamma, light),
                              model=MODEL_SPHERICAL)
    fit_heavy = fit_variogram(ExperimentalVariogram(lags, gamma, heavy),
                              model=MODEL_SPHERICAL)
    assert fit_heavy.sill > fit_light.sill


# ── 28. model müqayisəsi/seçimi ───────────────────────────────────────
def test_compare_variogram_models_returns_every_candidate_with_a_score():
    points, values = _gaussian_random_field(n=200, seed=13)
    exp = experimental_variogram(points, values, n_lags=14, max_lag=800.0)
    results = compare_variogram_models(exp)
    assert set(results) == set(KNOWN_MODELS)
    for model, params in results.items():
        assert params.model == model
        assert np.isfinite(params.weighted_rmse)


@pytest.mark.parametrize("model", KNOWN_MODELS)
def test_auto_selects_the_model_the_data_was_generated_from(model):
    """Sintetik γ məhz bir modeldən qurulub — `model="auto"` onu ədədi
    meyarla (çəkili RMSE) TAPMALIDIR, sabit favorit YOXDUR."""
    from imex2d.geology.variogram import ExperimentalVariogram
    lags = np.linspace(10.0, 400.0, 30)
    gamma = semivariance(lags, 0.0, 1.0, 200.0, model)
    exp = ExperimentalVariogram(lags, gamma, np.full(30, 200), lag_width=13.0)
    assert fit_variogram(exp, model="auto").model == model


def test_auto_choice_is_not_hard_coded_to_spherical():
    from imex2d.geology.variogram import ExperimentalVariogram
    lags = np.linspace(10.0, 400.0, 30)
    chosen = set()
    for model in KNOWN_MODELS:
        exp = ExperimentalVariogram(lags, semivariance(lags, 0.05, 1.0, 180.0, model),
                                    np.full(30, 150), lag_width=13.0)
        chosen.add(fit_variogram(exp, model="auto").model)
    assert len(chosen) == 3


# ── 29. istiqamətli variogram ─────────────────────────────────────────
def test_directional_variogram_detects_stronger_continuity_along_the_major_axis():
    """Major ox X (azimut 90°) boyunca 6 dəfə uzundur: 90° istiqamətdə
    γ DAHA YAVAŞ artmalıdır (davamlılıq daha çoxdur)."""
    points, values = _gaussian_random_field(n=300, seed=21, range_x=600.0,
                                            range_y=100.0, high=1000.0)
    along = experimental_variogram(points, values, n_lags=10, max_lag=350.0,
                                   azimuth_deg=90.0, azimuth_tolerance_deg=20.0)
    across = experimental_variogram(points, values, n_lags=10, max_lag=350.0,
                                    azimuth_deg=0.0, azimuth_tolerance_deg=20.0)
    mid = 5
    assert along.gamma[mid] < across.gamma[mid]
    assert along.azimuth_deg == 90.0 and across.azimuth_deg == 0.0


def test_directional_variograms_helper_returns_one_per_azimuth():
    points, values = _gaussian_random_field(n=200, seed=22)
    results = directional_variograms(points, values, azimuths=(0.0, 45.0, 90.0, 135.0),
                                     n_lags=8, max_lag=500.0)
    assert set(results) == {0.0, 45.0, 90.0, 135.0}
    for azimuth, exp in results.items():
        assert exp.azimuth_deg == azimuth
        assert exp.n_pairs_total > 0


def test_azimuth_tolerance_controls_how_many_pairs_are_kept():
    points, values = _gaussian_random_field(n=150, seed=23)
    narrow = experimental_variogram(points, values, azimuth_deg=30.0,
                                    azimuth_tolerance_deg=5.0)
    wide = experimental_variogram(points, values, azimuth_deg=30.0,
                                  azimuth_tolerance_deg=60.0)
    assert narrow.n_pairs_total < wide.n_pairs_total


def test_bandwidth_limits_the_lateral_spread_of_the_search_cone():
    points, values = _gaussian_random_field(n=200, seed=24, high=1000.0)
    wide = experimental_variogram(points, values, azimuth_deg=0.0,
                                  azimuth_tolerance_deg=45.0)
    narrow = experimental_variogram(points, values, azimuth_deg=0.0,
                                    azimuth_tolerance_deg=45.0, bandwidth=50.0)
    assert narrow.n_pairs_total < wide.n_pairs_total


def test_empty_direction_raises_instead_of_returning_empty_bins():
    points = np.array([[0., 0.], [100., 0.], [200., 0.], [300., 0.]])
    values = np.array([1.0, 2.0, 3.0, 4.0])
    with pytest.raises(ValueError, match="cüt tapılmadı"):
        experimental_variogram(points, values, azimuth_deg=0.0,
                               azimuth_tolerance_deg=5.0)


# ── 30. şaquli / 3D istiqamətli variogram ─────────────────────────────
def _layered_well_data(n_wells=12, n_depths=14, seed=31):
    """Quyu-tipli 3D data: güclü üfüqi davamlılıq, zəif şaquli."""
    rng = np.random.default_rng(seed)
    wells = rng.uniform(0.0, 1000.0, size=(n_wells, 2))
    depths = np.linspace(2000.0, 2130.0, n_depths)
    points, values = [], []
    for x, y in wells:
        trend = 0.20 + 0.00002 * x
        for depth in depths:
            points.append((x, y, depth))
            values.append(trend + 0.05 * np.sin(depth / 4.0) + rng.normal(0.0, 0.002))
    return np.asarray(points), np.asarray(values)


def test_vertical_variogram_uses_only_same_well_pairs():
    points, values = _layered_well_data()
    exp = vertical_variogram(points, values, horizontal_tolerance=1e-6, n_lags=8,
                             max_lag=200.0)      # dərinlik aralığı 130 m
    assert exp.vertical is True
    # 12 quyu × C(14,2) = 12 × 91 şaquli cüt
    assert exp.n_pairs_total == 12 * 14 * 13 // 2
    assert exp.gamma[exp.counts > 0].max() > 0.0


def test_vertical_variogram_raises_when_no_repeated_wells():
    rng = np.random.default_rng(41)
    points = np.column_stack([rng.uniform(0, 500, size=(20, 2)),
                              rng.uniform(2000, 2100, 20)])
    values = rng.standard_normal(20)
    with pytest.raises(ValueError, match="Şaquli variogram"):
        vertical_variogram(points, values, horizontal_tolerance=1e-9)


def test_horizontal_and_vertical_variograms_differ_for_layered_data():
    points, values = _layered_well_data()
    horizontal = experimental_variogram(points, values, dip_deg=0.0,
                                        dip_tolerance_deg=10.0, n_lags=8,
                                        max_lag=600.0)
    vertical = vertical_variogram(points, values, horizontal_tolerance=1e-6,
                                  n_lags=8, max_lag=100.0)
    assert horizontal.dip_deg == 0.0
    assert horizontal.n_pairs_total > 0 and vertical.n_pairs_total > 0
    assert not np.allclose(horizontal.gamma[:4], vertical.gamma[:4])


def test_dip_filter_selects_only_pairs_within_the_dip_window():
    points, values = _layered_well_data()
    steep = experimental_variogram(points, values, dip_deg=90.0, dip_tolerance_deg=5.0,
                                   n_lags=8, max_lag=140.0)
    flat = experimental_variogram(points, values, dip_deg=0.0, dip_tolerance_deg=5.0,
                                  n_lags=8, max_lag=1400.0)
    assert steep.n_pairs_total == 12 * 14 * 13 // 2      # yalnız eyni quyunun cütləri
    assert flat.n_pairs_total > 0
    assert steep.n_pairs_total != flat.n_pairs_total


def test_dip_requires_a_z_column():
    rng = np.random.default_rng(51)
    points = rng.uniform(0.0, 100.0, size=(20, 2))
    with pytest.raises(ValueError, match="Z sütunu yoxdur"):
        experimental_variogram(points, rng.standard_normal(20), dip_deg=30.0)


# ── 31. etibarsız parametrlərin rədd edilməsi ─────────────────────────
@pytest.mark.parametrize("nugget,sill,range_", [
    (-0.1, 1.0, 100.0),        # mənfi nugget
    (0.0, 0.0, 100.0),         # sıfır sill
    (0.0, -1.0, 100.0),        # mənfi sill
    (0.0, 1.0, 0.0),           # sıfır radius
    (0.0, 1.0, -50.0),         # mənfi radius
    (np.nan, 1.0, 100.0),      # qeyri-sonlu
    (0.0, np.inf, 100.0),
    (0.0, 1.0, np.nan),
])
def test_invalid_parameters_are_rejected_before_reaching_the_solver(nugget, sill, range_):
    with pytest.raises(VariogramValidationError):
        validate_variogram_parameters(MODEL_SPHERICAL, nugget, sill, range_)


def test_unknown_model_is_rejected_by_the_validator():
    with pytest.raises(VariogramValidationError, match="Naməlum"):
        validate_variogram_parameters("linear", 0.0, 1.0, 100.0)


def test_valid_parameters_pass_but_dominant_nugget_warns():
    assert validate_variogram_parameters(MODEL_SPHERICAL, 0.1, 1.0, 100.0) == []
    warnings = validate_variogram_parameters(MODEL_SPHERICAL, 9.5, 0.5, 100.0)
    assert warnings and "nugget" in warnings[0]


def test_repair_is_explicit_and_reports_every_change():
    nugget, sill, range_, warnings = repair_variogram_parameters(
        MODEL_SPHERICAL, -1.0, 0.0, -5.0)
    assert nugget == 0.0 and sill > 0.0 and range_ > 0.0
    assert len(warnings) == 3


def test_kriging_rejects_an_invalid_model_name_at_construction():
    from imex2d.geology.interpolation import OrdinaryKriging
    with pytest.raises(ValueError, match="Naməlum variogram modeli"):
        OrdinaryKriging(model="linear")


def test_fit_result_exposes_total_sill_and_helper_methods():
    from imex2d.geology.variogram import ExperimentalVariogram
    lags = np.linspace(10.0, 300.0, 15)
    exp = ExperimentalVariogram(lags, semivariance(lags, 0.2, 1.0, 150.0,
                                                   MODEL_SPHERICAL),
                                np.full(15, 80), lag_width=20.0)
    fit = fit_variogram(exp, model=MODEL_SPHERICAL)
    assert fit.total_sill == pytest.approx(fit.nugget + fit.sill)
    assert fit.semivariance(0.0) == pytest.approx(fit.nugget)
    assert fit.covariance(0.0) == pytest.approx(fit.total_sill)
    assert fit.validate() == [] or isinstance(fit.validate(), list)
