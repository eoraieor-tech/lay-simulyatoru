"""B8 №25-33 — MƏLUMAT KEYFİYYƏTİ boru xətti (GATE B8/B9).

İki mərkəzi iddia:

1. QC variogramdan/interpolyasiyadan ƏVVƏL işləyir və heç bir qərar
   SƏSSİZ deyil — hər şey `DataQualityReport`-dadır.
2. FİZİKİ CƏHƏTDƏN ETİBARSIZ ≠ STATİSTİK KƏNAR-DƏYƏR. Birincisi
   çıxarılır, ikincisi DEFOLT olaraq YALNIZ işarələnir.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.geology.data_quality import (DataQualityError, QCConfig, detect_outliers,
                                         run_quality_control)
from imex2d.geology.property_config import (DuplicatePolicy, OutlierMethod,
                                            resolve_strategy)
from imex2d.geology.property_interpolation import interpolate_property_field


def _points(n, seed=0, high=1000.0):
    rng = np.random.default_rng(seed)
    return rng.uniform(0.0, high, size=(n, 2))


def _poro(n, seed=0):
    rng = np.random.default_rng(seed + 500)
    return 0.15 + 0.05 * rng.random(n)


# ── 25. NaN ───────────────────────────────────────────────────────────
def test_nan_values_are_removed_and_reported():
    points, values = _points(10, 1), _poro(10, 1)
    values[3] = np.nan
    values[7] = np.nan
    result = run_quality_control(points, values, resolve_strategy("PORO"))
    assert result.report.n_non_finite == 2
    assert result.report.n_valid == 8
    assert result.values.size == 8
    assert np.all(np.isfinite(result.values))
    assert any(f.kind == "non_finite" and f.action == "removed"
               for f in result.report.findings)


def test_nan_policy_raise_is_available_and_explicit():
    points, values = _points(8, 2), _poro(8, 2)
    values[1] = np.nan
    with pytest.raises(DataQualityError, match="NaN"):
        run_quality_control(points, values, resolve_strategy("PORO"),
                            QCConfig(non_finite_policy="raise"))


# ── 26. Inf ───────────────────────────────────────────────────────────
def test_positive_and_negative_infinity_are_both_removed():
    points, values = _points(10, 3), _poro(10, 3)
    values[2] = np.inf
    values[5] = -np.inf
    result = run_quality_control(points, values, resolve_strategy("PORO"))
    assert result.report.n_non_finite == 2
    assert np.all(np.isfinite(result.values))


# ── 27. etibarsız koordinatlar ────────────────────────────────────────
def test_non_finite_coordinates_are_removed_not_repaired():
    points, values = _points(8, 4), _poro(8, 4)
    points[2, 0] = np.nan
    points[5, 1] = np.inf
    result = run_quality_control(points, values, resolve_strategy("PORO"))
    assert result.report.n_invalid_coordinates == 2
    assert result.report.n_valid == 6
    assert np.all(np.isfinite(result.points))
    finding = next(f for f in result.report.findings if f.kind == "coordinate")
    assert finding.severity == "error"
    assert "TƏXMİN EDİLMİR" in finding.detail


def test_mismatched_lengths_raise_immediately():
    with pytest.raises(DataQualityError, match="uyğun gəlmir"):
        run_quality_control(_points(6, 5), _poro(4, 5), resolve_strategy("PORO"))


@pytest.mark.parametrize("bad", [np.zeros((4, 1)), np.zeros((4, 5))])
def test_wrong_coordinate_dimensionality_is_rejected(bad):
    with pytest.raises(DataQualityError, match="n,2"):
        run_quality_control(bad, np.zeros(4), resolve_strategy("PORO"))


def test_two_dimensional_input_is_padded_to_three():
    result = run_quality_control(_points(6, 6), _poro(6, 6), resolve_strategy("PORO"))
    assert result.points.shape[1] == 3
    assert np.all(result.points[:, 2] == 0.0)


# ── 28. dublikatlar ───────────────────────────────────────────────────
def test_identical_duplicates_are_merged_and_counted():
    points = np.array([[0., 0.], [0., 0.], [100., 0.], [0., 100.], [100., 100.]])
    values = np.array([0.20, 0.20, 0.25, 0.22, 0.24])
    result = run_quality_control(points, values, resolve_strategy("PORO"))
    assert result.report.n_duplicate_locations == 1
    assert result.report.n_duplicate_observations == 2
    assert result.report.n_conflicting_duplicates == 0
    assert result.values.size == 4


@pytest.mark.parametrize("policy,expected", [
    (DuplicatePolicy.MEAN, 0.30),
    (DuplicatePolicy.MEDIAN, 0.30),
    (DuplicatePolicy.KEEP_FIRST, 0.20),
    (DuplicatePolicy.KEEP_LAST, 0.40),
])
def test_duplicate_policies_are_explicit_and_deterministic(policy, expected):
    points = np.array([[0., 0.], [0., 0.], [100., 0.], [0., 100.], [100., 100.]])
    values = np.array([0.20, 0.40, 0.25, 0.22, 0.24])
    result = run_quality_control(
        points, values, resolve_strategy("PORO").derive(duplicate_policy=policy))
    assert result.values[0] == pytest.approx(expected)
    assert result.report.duplicate_policy == policy.value


def test_duplicate_policy_keep_separate_does_not_merge():
    points = np.array([[0., 0.], [0., 0.], [100., 0.], [0., 100.], [100., 100.]])
    values = np.array([0.20, 0.40, 0.25, 0.22, 0.24])
    result = run_quality_control(
        points, values,
        resolve_strategy("PORO").derive(duplicate_policy=DuplicatePolicy.KEEP_SEPARATE))
    assert result.values.size == 5
    assert any(f.action == "kept" for f in result.report.findings)


def test_duplicate_policy_raise_refuses_to_guess():
    points = np.array([[0., 0.], [0., 0.], [100., 0.], [0., 100.]])
    values = np.array([0.20, 0.40, 0.25, 0.22])
    with pytest.raises(DataQualityError, match="AÇIQ strategiya"):
        run_quality_control(
            points, values,
            resolve_strategy("PORO").derive(duplicate_policy=DuplicatePolicy.RAISE))


# ── 29. ZİDDİYYƏTLİ dublikatlar ───────────────────────────────────────
def test_conflicting_duplicates_are_flagged_separately_from_plain_ones():
    points = np.array([[0., 0.], [0., 0.], [50., 50.], [50., 50.],
                       [100., 0.], [0., 100.]])
    values = np.array([0.20, 0.20,       # eyni → ziddiyyət DEYİL
                       0.10, 0.40,       # fərqli → ZİDDİYYƏT
                       0.25, 0.22])
    result = run_quality_control(points, values, resolve_strategy("PORO"))
    assert result.report.n_duplicate_locations == 2
    assert result.report.n_conflicting_duplicates == 1
    finding = next(f for f in result.report.findings if f.kind == "duplicate")
    assert finding.severity == "warning"
    assert result.values[1] == pytest.approx(0.25)     # (0.10+0.40)/2


def test_duplicate_groups_are_reported_with_raw_indices_for_audit():
    points = np.array([[0., 0.], [100., 0.], [0., 0.], [0., 100.]])
    values = np.array([0.20, 0.25, 0.40, 0.22])
    result = run_quality_control(points, values, resolve_strategy("PORO"))
    assert len(result.duplicate_groups) == 1
    assert set(result.duplicate_groups[0]) == {0, 2}


def test_categorical_majority_policy_resolves_duplicate_codes():
    points = np.array([[0., 0.], [0., 0.], [0., 0.], [100., 0.], [0., 100.]])
    codes = np.array([1.0, 1.0, 2.0, 2.0, 1.0])
    result = run_quality_control(points, codes, resolve_strategy("FACIES"))
    assert result.values[0] == pytest.approx(1.0)      # 2 səs vs 1


def test_categorical_majority_refuses_a_tie():
    points = np.array([[0., 0.], [0., 0.], [100., 0.], [0., 100.]])
    codes = np.array([1.0, 2.0, 2.0, 1.0])
    with pytest.raises(DataQualityError, match="BƏRABƏRDİR"):
        run_quality_control(points, codes, resolve_strategy("FACIES"))


# ── 30. fiziki cəhətdən etibarsız hədlər ──────────────────────────────
def test_physically_invalid_porosity_is_removed_with_a_reason():
    points, values = _points(10, 7), _poro(10, 7)
    values[1] = -0.05          # mənfi məsaməlik — FİZİKİ olaraq mümkünsüz
    values[6] = 1.4            # 140% məsaməlik
    result = run_quality_control(points, values, resolve_strategy("PORO"))
    assert result.report.n_bound_violations == 2
    assert result.report.n_valid == 8
    finding = next(f for f in result.report.findings if f.kind == "bounds")
    assert finding.severity == "error"
    assert "STATİSTİK kənar-dəyər DEYİL" in finding.detail


def test_saturation_outside_zero_one_is_rejected():
    points = _points(8, 8)
    values = np.array([0.1, 0.4, 1.6, 0.5, -0.2, 0.7, 0.3, 0.9])
    result = run_quality_control(points, values, resolve_strategy("SW"))
    assert result.report.n_bound_violations == 2
    assert np.all(result.values >= 0.0) and np.all(result.values <= 1.0)


def test_non_positive_permeability_is_rejected_for_the_log_path():
    points = _points(8, 9)
    values = np.array([10.0, 50.0, 0.0, 200.0, -3.0, 80.0, 15.0, 400.0])
    result = run_quality_control(points, values, resolve_strategy("PERMX"))
    assert result.report.n_bound_violations == 2
    assert np.all(result.values > 0.0)


def test_invalid_facies_code_outside_the_declared_set_is_rejected():
    points = _points(6, 10)
    codes = np.array([1.0, 2.0, 9.0, 1.0, 2.0, 2.0])
    strategy = resolve_strategy("FACIES").derive(categories=(1, 2))
    result = run_quality_control(points, codes, strategy)
    assert result.report.n_bound_violations == 1
    assert set(np.unique(result.values).tolist()) == {1.0, 2.0}


def test_bounds_policy_keep_preserves_invalid_values_with_a_warning():
    points, values = _points(8, 11), _poro(8, 11)
    values[2] = 1.7
    result = run_quality_control(points, values, resolve_strategy("PORO"),
                                 QCConfig(invalid_bounds_policy="keep"))
    assert result.report.n_valid == 8
    finding = next(f for f in result.report.findings if f.kind == "bounds")
    assert finding.action == "kept"


def test_bounds_policy_raise_is_available():
    points, values = _points(8, 12), _poro(8, 12)
    values[2] = -1.0
    with pytest.raises(DataQualityError, match="etibarsız"):
        run_quality_control(points, values, resolve_strategy("PORO"),
                            QCConfig(invalid_bounds_policy="raise"))


# ── 31. NADİR amma ETİBARLI dəyərlər ──────────────────────────────────
def test_a_rare_but_physically_valid_value_is_never_removed_by_default():
    """5000 mD keçiricilik ÇAT ZONASI ola bilər — fiziki olaraq etibarlıdır,
    ona görə DEFOLT olaraq SİLİNMİR, yalnız işarələnir (B4.4/B4.5)."""
    points = _points(30, 13)
    values = np.full(30, 50.0)
    values[10] = 5000.0
    result = run_quality_control(points, values, resolve_strategy("PERMX"))
    assert result.report.n_bound_violations == 0        # fiziki cəhətdən etibarlı
    assert result.report.n_outlier_candidates >= 1      # amma kənar-dəyər namizədi
    assert result.report.n_outliers_removed == 0        # SİLİNMƏYİB
    assert result.values.size == 30
    assert 5000.0 in result.values


def test_extreme_valid_value_still_reaches_the_interpolator():
    points = _points(25, 14)
    values = np.full(25, 100.0)
    values[5] = 9000.0
    result = interpolate_property_field(points, values, points,
                                        property_name="PERMX")
    assert result.estimate[5] == pytest.approx(9000.0, rel=1e-6)


# ── 32. kənar-dəyər aşkarlanması ──────────────────────────────────────
def test_mad_outlier_detection_finds_the_planted_outlier():
    values = np.concatenate([np.full(29, 10.0) + np.linspace(-0.5, 0.5, 29),
                             [500.0]])
    flagged, score = detect_outliers(values, OutlierMethod.MAD, threshold=3.5)
    assert flagged[-1] and int(np.sum(flagged)) == 1
    assert score[-1] > score[:-1].max()


def test_iqr_outlier_detection_finds_both_tails():
    values = np.concatenate([np.linspace(9.0, 11.0, 28), [-40.0], [60.0]])
    flagged, _ = detect_outliers(values, OutlierMethod.IQR, threshold=1.5)
    assert flagged[-1] and flagged[-2]


def test_spatial_outlier_detection_finds_a_locally_anomalous_point():
    """Qlobal paylanmada NORMAL, amma QONŞULARINDAN kəskin fərqlənən nöqtə —
    geostatistikada ən mənalı kənar-dəyər növü (B4.5)."""
    axis = np.linspace(0.0, 300.0, 8)
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    points = np.column_stack([xx.ravel(), yy.ravel()])
    values = 0.1 + 0.002 * points[:, 0]           # hamar trend
    rogue = 30
    values[rogue] = values.max()                  # qlobal olaraq normal…
    flagged, _ = detect_outliers(values, OutlierMethod.SPATIAL, threshold=3.5,
                                 points=np.column_stack([points, np.zeros(64)]))
    global_flagged, _ = detect_outliers(values, OutlierMethod.MAD, threshold=3.5)
    assert flagged[rogue], "…amma yerli olaraq anomaldır"
    assert not global_flagged[rogue], "qlobal MAD bunu görə BİLMİR"


def test_outlier_detection_is_disabled_cleanly():
    values = np.concatenate([np.full(20, 1.0), [1000.0]])
    flagged, _ = detect_outliers(values, OutlierMethod.NONE, threshold=3.5)
    assert not np.any(flagged)


def test_constant_data_produces_no_outliers():
    """MAD = 0 olanda "hər şey kənardır" nəticəsi VERİLMİR."""
    flagged, _ = detect_outliers(np.full(30, 7.0), OutlierMethod.MAD, 3.5)
    assert not np.any(flagged)


def test_outliers_are_removed_only_with_an_explicit_opt_in():
    points = _points(30, 15)
    values = np.full(30, 20.0) + np.linspace(-1.0, 1.0, 30)
    values[7] = 900.0
    kept = run_quality_control(points, values, resolve_strategy("PERMX"))
    removed = run_quality_control(
        points, values, resolve_strategy("PERMX").derive(remove_outliers=True))
    assert kept.report.n_outliers_removed == 0 and kept.values.size == 30
    assert removed.report.n_outliers_removed >= 1
    assert removed.values.size < 30
    assert any("AÇIQ `remove_outliers=True`" in f.detail
               for f in removed.report.findings if f.kind == "outlier")


def test_outlier_mask_is_exposed_for_audit_even_when_nothing_is_removed():
    points = _points(30, 16)
    values = np.full(30, 5.0) + np.linspace(-0.2, 0.2, 30)
    values[11] = 400.0
    report = run_quality_control(points, values, resolve_strategy("PERMX")).report
    assert report.outlier_mask is not None
    assert bool(report.outlier_mask[11])
    assert int(report.outlier_mask.sum()) == report.n_outlier_candidates


# ── 33. QC hesabatı ───────────────────────────────────────────────────
def test_report_is_structured_and_machine_readable():
    points, values = _points(12, 17), _poro(12, 17)
    values[0] = np.nan
    values[1] = 1.5
    points[2] = points[3]
    report = run_quality_control(points, values, resolve_strategy("PORO")).report

    data = report.as_dict()
    for key in ("property", "n_input", "n_valid", "n_removed", "n_non_finite",
                "n_invalid_coordinates", "n_duplicate_locations",
                "n_bound_violations", "n_outlier_candidates", "duplicate_policy",
                "outlier_method", "findings", "warnings"):
        assert key in data
    assert data["n_input"] == 12
    assert data["n_valid"] == report.n_valid
    assert isinstance(data["findings"], list)
    assert all(isinstance(f, dict) for f in data["findings"])


def test_report_text_lists_every_stage():
    points, values = _points(10, 18), _poro(10, 18)
    values[0] = np.inf
    text = run_quality_control(points, values, resolve_strategy("PORO")).report.as_text()
    for fragment in ("Məlumat keyfiyyəti", "qeyri-sonlu", "dublikat", "kənar-dəyər"):
        assert fragment in text


def test_kept_mask_maps_back_to_the_raw_input_rows():
    points, values = _points(10, 19), _poro(10, 19)
    values[4] = np.nan
    result = run_quality_control(points, values, resolve_strategy("PORO"))
    assert result.report.kept_mask.shape == (10,)
    assert not result.report.kept_mask[4]
    assert int(result.report.kept_mask.sum()) == result.values.size


def test_qc_runs_before_interpolation_and_is_attached_to_the_estimate():
    """GATE B8 — QC interpolyasiyadan ƏVVƏL, nəticə obyektində görünür."""
    points, values = _points(20, 20), _poro(20, 20)
    values[3] = np.nan
    values[9] = 2.0
    result = interpolate_property_field(points, values, _points(10, 21),
                                        property_name="PORO")
    assert result.quality is not None
    assert result.quality.n_input == 20
    assert result.quality.n_valid == 18
    assert any("QC-də çıxarıldı" in w for w in result.warnings)


def test_qc_can_be_disabled_explicitly():
    points, values = _points(15, 22), _poro(15, 22)
    result = interpolate_property_field(points, values, _points(8, 23),
                                        property_name="PORO", run_qc=False)
    assert result.quality is None


def test_all_observations_invalid_produces_nan_and_says_so():
    points = _points(6, 24)
    values = np.full(6, np.nan)
    result = interpolate_property_field(points, values, _points(5, 25),
                                        property_name="PORO")
    assert np.all(np.isnan(result.estimate))
    assert result.quality.n_valid == 0
    assert any("etibarlı" in w for w in result.warnings)


def test_qc_warns_when_too_few_points_remain_for_a_variogram():
    points = _points(6, 26)
    values = _poro(6, 26)
    values[:3] = np.nan
    report = run_quality_control(points, values, resolve_strategy("PORO")).report
    assert report.n_valid == 3
    assert any("variogram" in w for w in report.warnings)
