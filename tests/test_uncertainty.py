"""B8 №13-17 — QEYRİ-MÜƏYYƏNLİK SİSTEMİ (GATE B5).

Mərkəzi iddia: qaytarılan qeyri-müəyyənlik ƏSL kriginq/simulyasiya
qeyri-müəyyənliyidir — `1/(1+məsafə)` kimi uydurma "bal" DEYİL, və hər
kəmiyyətin MƏNASI (çevrilmiş fəza? geri çevrilmiş? diaqnostika?) AÇIQ
ayrılır (B2.1).
"""

from __future__ import annotations

import numpy as np

from imex2d.geology.property_config import UncertaintyKind, resolve_strategy
from imex2d.geology.property_interpolation import (CONFIDENCE_VARIANCE_HIGH,
                                                   Confidence, classify_confidence,
                                                   compute_data_density,
                                                   interpolate_categorical_field,
                                                   interpolate_property_field)
from imex2d.geology.spatial_search import (SUPPORT_BOUNDARY, SUPPORT_EXTRAPOLATED,
                                           SUPPORT_WEAK, SUPPORT_WELL)
from imex2d.geology.transforms import VarianceKind


def _cluster(n=40, seed=1, high=600.0):
    rng = np.random.default_rng(seed)
    return rng.uniform(0.0, high, size=(n, 2))


# ── 13. kriginq variansı ──────────────────────────────────────────────
def test_variance_is_zero_at_data_and_grows_away_from_it():
    points = _cluster(30, seed=2)
    values = 0.18 + 0.05 * np.sin(points[:, 0] / 130.0)
    targets = np.vstack([points[:3], [[300.0, 300.0]], [[5000.0, 5000.0]]])
    result = interpolate_property_field(points, values, targets,
                                        property_name="PORO")
    assert np.allclose(result.variance[:3], 0.0, atol=1e-9)
    assert result.variance[3] < result.variance[4]


def test_variance_is_never_negative():
    points = _cluster(50, seed=3)
    values = np.exp(2.0 + np.cos(points[:, 1] / 90.0))
    targets = _cluster(80, seed=4, high=900.0)
    result = interpolate_property_field(points, values, targets,
                                        property_name="PERMX")
    finite = np.isfinite(result.variance)
    assert np.all(result.variance[finite] >= -1e-12)
    assert np.all(result.std[finite] >= 0.0)


def test_idw_reports_no_variance_instead_of_inventing_one():
    """GATE B5 — variansı olmayan üsul üçün ədəd UYDURULMUR."""
    from imex2d.geology.property_config import InterpolationKind
    points = _cluster(25, seed=5)
    values = 0.2 + 0.03 * points[:, 0] / 600.0
    result = interpolate_property_field(
        points, values, _cluster(20, seed=6),
        strategy=resolve_strategy("PORO").derive(
            interpolation=InterpolationKind.IDW))
    assert np.all(np.isnan(result.variance))
    assert result.variance_kind is VarianceKind.UNDEFINED
    assert any("UYDURULMUR" in w for w in result.warnings)


# ── 14. çevrilmiş fəza variansı ───────────────────────────────────────
def test_transformed_and_back_transformed_variance_are_distinct_quantities():
    """B2.1 — ikisi QARIŞDIRILMIR və hansının hansı olduğu bildirilir."""
    rng = np.random.default_rng(7)
    points = _cluster(40, seed=7)
    values = np.exp(rng.normal(4.0, 1.2, size=40))
    result = interpolate_property_field(points, values, _cluster(30, seed=8),
                                        property_name="PERMX")

    assert result.uncertainty_kind is UncertaintyKind.BACK_TRANSFORMED
    assert result.variance_kind is VarianceKind.EXACT
    away = result.transformed_variance > 1e-9
    assert np.any(away)
    # loq fəzasında varians kiçik ədəddir (ln miqyası), fiziki fəzada
    # isə mD² miqyasında — QARIŞDIRILA BİLMƏZ
    assert np.all(result.variance[away] > result.transformed_variance[away])


def test_back_transformed_variance_matches_the_closed_form():
    rng = np.random.default_rng(9)
    points = _cluster(35, seed=9)
    values = np.exp(rng.normal(3.5, 1.0, size=35))
    result = interpolate_property_field(points, values, _cluster(25, seed=10),
                                        property_name="PERMX")
    y, s2 = result.transformed_estimate, result.transformed_variance
    expected = np.exp(2 * y + s2) * (np.exp(s2) - 1.0)
    assert np.allclose(result.variance, expected, rtol=1e-9)


def test_identity_transform_leaves_variance_untouched():
    points = _cluster(30, seed=11)
    values = 0.2 + 0.04 * np.sin(points[:, 0] / 110.0)
    result = interpolate_property_field(points, values, _cluster(20, seed=12),
                                        property_name="PORO")
    assert result.variance_kind is VarianceKind.IDENTITY
    assert np.allclose(result.variance, result.transformed_variance, atol=0.0)
    assert result.uncertainty_kind is UncertaintyKind.KRIGING_VARIANCE


def test_bounded_variance_is_marked_as_a_delta_approximation():
    """Logit-normalın variansı qapalı formada YOXDUR — DELTA kimi bildirilir."""
    points = _cluster(30, seed=13)
    values = np.clip(0.4 + 0.2 * np.sin(points[:, 1] / 140.0), 0.0, 1.0)
    result = interpolate_property_field(points, values, _cluster(20, seed=14),
                                        property_name="SW")
    assert result.variance_kind is VarianceKind.DELTA


# ── 15. dəstək təsnifatı ──────────────────────────────────────────────
def test_confidence_classes_are_derived_from_real_diagnostics():
    support = np.array([SUPPORT_WELL, SUPPORT_WELL, SUPPORT_BOUNDARY,
                        SUPPORT_WEAK, SUPPORT_EXTRAPOLATED], dtype=object)
    variance = np.array([0.05, 0.90, 0.05, 0.05, 0.05])
    counts = np.array([12, 12, 12, 12, 0])
    labels = classify_confidence(support, variance, total_sill=1.0,
                                 neighbor_count=counts)
    assert labels[0] == Confidence.HIGH.value        # sıx + kiçik varians
    assert labels[1] == Confidence.LOW.value         # varians böyükdür
    assert labels[2] == Confidence.LOW.value         # bir tərəfdə məlumat
    assert labels[3] == Confidence.LOW.value
    assert labels[4] == Confidence.EXTRAPOLATED.value


def test_confidence_threshold_is_the_documented_relative_variance():
    support = np.full(3, SUPPORT_WELL, dtype=object)
    counts = np.full(3, 10)
    variance = np.array([CONFIDENCE_VARIANCE_HIGH * 0.5,
                         CONFIDENCE_VARIANCE_HIGH * 1.5, 0.9])
    labels = classify_confidence(support, variance, 1.0, counts)
    assert labels[0] == Confidence.HIGH.value
    assert labels[1] == Confidence.MEDIUM.value
    assert labels[2] == Confidence.LOW.value


def test_confidence_is_documented_as_interpretation_not_probability():
    """GATE B5 — sinif docstring-i bunun kalibrlənmiş ehtimal OLMADIĞINI
    AÇIQ yazır; test sənədin yerində olduğunu qoruyur."""
    assert "KALİBRLƏNMİŞ EHTİMAL DEYİL" in Confidence.__doc__


def test_support_classification_reaches_the_estimate_object():
    points = _cluster(40, seed=15)
    values = 0.2 + 0.05 * np.cos(points[:, 0] / 100.0)
    targets = np.vstack([[[300.0, 300.0]], [[9000.0, 9000.0]]])
    result = interpolate_property_field(points, values, targets,
                                        property_name="PORO")
    assert result.support[0] == SUPPORT_WELL
    assert result.support[1] == SUPPORT_EXTRAPOLATED


# ── 16. ekstrapolyasiya bayrağı ───────────────────────────────────────
def test_extrapolation_flag_is_geometric_not_a_neighbour_count():
    points = _cluster(50, seed=16)
    values = 0.2 + 0.05 * np.sin(points[:, 0] / 120.0)
    inside = np.array([[300.0, 300.0]])
    far = np.array([[20000.0, 20000.0]])
    a = interpolate_property_field(points, values, inside, property_name="PORO")
    b = interpolate_property_field(points, values, far, property_name="PORO")
    assert a.neighbor_count[0] == b.neighbor_count[0]      # EYNİ qonşu sayı
    assert not a.extrapolated[0] and b.extrapolated[0]     # FƏRQLİ təsnifat


def test_extrapolated_cells_are_labelled_extrapolated_in_confidence():
    points = _cluster(30, seed=17)
    values = 0.2 + 0.02 * points[:, 0] / 600.0
    result = interpolate_property_field(points, values,
                                        np.array([[50000.0, 50000.0]]),
                                        property_name="PORO")
    assert result.confidence[0] == Confidence.EXTRAPOLATED.value


# ── data sıxlığı — MÜSTƏQİL diaqnostika ───────────────────────────────
def test_data_density_is_a_count_not_a_probability():
    points = np.array([[0., 0.], [10., 0.], [20., 0.], [1000., 1000.]])
    targets = np.array([[10.0, 0.0], [1000.0, 1000.0]])
    density = compute_data_density(np.column_stack([points, np.zeros(4)]),
                                   np.column_stack([targets, np.zeros(2)]),
                                   radius=50.0)
    assert density[0] == 3 and density[1] == 1
    assert density.dtype.kind in "iu"


def test_data_density_reaches_the_estimate_and_differs_from_variance():
    rng = np.random.default_rng(18)
    dense = rng.uniform(0.0, 200.0, size=(40, 2))
    sparse = rng.uniform(800.0, 1000.0, size=(4, 2))
    points = np.vstack([dense, sparse])
    values = 0.2 + 0.01 * points[:, 0] / 1000.0
    targets = np.array([[100.0, 100.0], [900.0, 900.0]])
    result = interpolate_property_field(points, values, targets,
                                        property_name="PORO")
    assert result.data_density[0] > result.data_density[1]


# ── 17. standartlaşdırılmış xəta ──────────────────────────────────────
def test_standardized_error_is_computed_from_real_kriging_variance():
    """B2.2 — `e = (z − ẑ)/σ`; kalibrləmə metrikləri hesablanır."""
    from imex2d.geology.cross_validation import (ValidationDesign, ValidationKind,
                                                 cross_validate_property)
    rng = np.random.default_rng(19)
    points = rng.uniform(0.0, 1000.0, size=(60, 2))
    values = 0.2 + 0.04 * np.sin(points[:, 0] / 200.0) + 0.005 * rng.standard_normal(60)
    metrics = cross_validate_property(
        points, values, resolve_strategy("PORO"),
        ValidationDesign(kind=ValidationKind.LEAVE_ONE_OUT))

    assert metrics.n_with_variance > 0
    assert np.isfinite(metrics.mean_standardized_error)
    assert np.isfinite(metrics.variance_standardized_error)
    assert 0.0 <= metrics.coverage_68 <= 1.0
    assert 0.0 <= metrics.coverage_95 <= 1.0
    assert metrics.coverage_95 >= metrics.coverage_68
    assert np.isfinite(metrics.calibration_error)


def test_standardized_error_detects_an_underestimated_variance():
    """Variansı SÜNİ olaraq kiçildəndə `var(e)` 1-dən BÖYÜK olmalıdır —
    metrik həqiqətən kalibrləməni ölçür, formal ədəd deyil."""
    from imex2d.geology.cross_validation import (ValidationDesign,
                                                 cross_validate_property)
    rng = np.random.default_rng(20)
    points = rng.uniform(0.0, 800.0, size=(50, 2))
    values = 0.2 + 0.05 * np.sin(points[:, 0] / 150.0) + 0.01 * rng.standard_normal(50)

    honest = cross_validate_property(points, values, resolve_strategy("PORO"),
                                     ValidationDesign())
    shrunk = cross_validate_property(
        points, values,
        resolve_strategy("PORO").derive(variogram_model="spherical"),
        ValidationDesign(), kriging_overrides={"sill": 1e-8, "range_": 500.0})
    assert shrunk.variance_standardized_error > honest.variance_standardized_error


def test_metrics_report_when_variance_is_unavailable_instead_of_faking_it():
    from imex2d.geology.cross_validation import (ValidationDesign,
                                                 cross_validate_property)
    from imex2d.geology.property_config import InterpolationKind
    rng = np.random.default_rng(21)
    points = rng.uniform(0.0, 500.0, size=(30, 2))
    values = 0.2 + 0.03 * rng.random(30)
    metrics = cross_validate_property(
        points, values,
        resolve_strategy("PORO").derive(interpolation=InterpolationKind.IDW),
        ValidationDesign())
    assert metrics.n_with_variance == 0
    assert np.isnan(metrics.variance_standardized_error)
    assert any("UYDURULMUR" in w for w in metrics.warnings)


# ── B2.4 xəritə dəstəyi ───────────────────────────────────────────────
def test_estimate_exposes_grid_ready_arrays():
    points = _cluster(30, seed=22)
    values = np.exp(3.0 + np.sin(points[:, 0] / 100.0))
    result = interpolate_property_field(points, values, _cluster(24, seed=23),
                                        property_name="PERMX")
    grids = result.as_grids()
    for key in ("estimate", "variance", "std", "transformed_variance",
                "nearest_distance", "neighbor_count", "data_density",
                "extrapolated", "confidence_rank"):
        assert key in grids
        assert grids[key].shape == (24,)
        assert grids[key].dtype.kind == "f"


def test_categorical_estimate_exposes_probability_grids():
    points = _cluster(30, seed=24)
    codes = (points[:, 0] > 300).astype(int)
    result = interpolate_categorical_field(points, codes, _cluster(20, seed=25),
                                           property_name="FACIES")
    grids = result.as_grids()
    assert "probability_0" in grids and "probability_1" in grids
    assert "entropy" in grids and "normalized_entropy" in grids
    assert np.allclose(grids["probability_0"] + grids["probability_1"], 1.0)


def test_summary_text_mentions_the_uncertainty_kind():
    points = _cluster(25, seed=26)
    values = np.exp(3.0 + 0.5 * np.cos(points[:, 1] / 80.0))
    result = interpolate_property_field(points, values, _cluster(15, seed=27),
                                        property_name="PERMX")
    text = result.summary()
    assert "varians" in text and "exact" in text
