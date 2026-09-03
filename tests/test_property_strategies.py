"""B8 №1-7 — hər rezervuar xassəsinin ÖZ statistik rejimi (GATE B1-B4).

Bu faylın mərkəzi iddiası: **PORO ≠ PERMX ≠ SW ≠ NTG ≠ FACIES**.
Eyni nöqtələr və eyni hədəflərlə fərqli xassələr FƏRQLİ yollardan keçir
və nəticələr həmin xassənin fiziki/statistik təbiətinə uyğun olur.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.geology.property_config import (DEFAULT_STRATEGIES, BackTransform,
                                            BoundPolicy, InterpolationKind,
                                            PropertyConfigError, VariableType,
                                            normal_score_strategy, resolve_strategy,
                                            strategy_table)
from imex2d.geology.property_interpolation import (CategoricalEstimate, Confidence,
                                                   PropertyEstimate,
                                                   interpolate_by_name,
                                                   interpolate_categorical_field,
                                                   interpolate_property_field)
from imex2d.geology.transforms import VarianceKind


def _wells(n=40, seed=1, high=1000.0):
    rng = np.random.default_rng(seed)
    return rng.uniform(0.0, high, size=(n, 2))


def _grid(n=15, high=1000.0):
    axis = np.linspace(50.0, high - 50.0, n)
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


# ── 1. PORO ───────────────────────────────────────────────────────────
def test_porosity_uses_untransformed_kriging_and_stays_in_bounds():
    points = _wells(seed=11)
    values = 0.12 + 0.10 * np.sin(points[:, 0] / 260.0) ** 2
    result = interpolate_property_field(points, values, _grid(), property_name="PORO")

    assert isinstance(result, PropertyEstimate)
    assert result.strategy.variable_type is VariableType.CONTINUOUS
    assert result.strategy.transform.is_identity
    assert result.variance_kind is VarianceKind.IDENTITY
    assert np.all(result.estimate >= 0.0) and np.all(result.estimate <= 1.0)
    # çevirmə yoxdursa çevrilmiş fəza = orijinal fəza
    assert np.allclose(result.transformed_estimate, result.raw_estimate, atol=1e-12)


def test_porosity_honours_hard_data_exactly():
    points = _wells(20, seed=12)
    values = 0.15 + 0.08 * np.cos(points[:, 1] / 180.0)
    result = interpolate_property_field(points, values, points, property_name="PORO")
    assert np.allclose(result.estimate, values, atol=1e-9)
    assert np.allclose(result.variance[np.isfinite(result.variance)], 0.0, atol=1e-9)


def test_porosity_bound_clipping_is_counted_not_silent():
    """Hədd kəsilməsi HƏMİŞƏ sayılır və xam dəyər saxlanılır (B1.1)."""
    strategy = resolve_strategy("PORO").derive(output_bounds=(0.0, 0.16))
    points = _wells(25, seed=13)
    values = 0.10 + 0.15 * (points[:, 0] / 1000.0)
    result = interpolate_property_field(points, values, _grid(), strategy=strategy)
    assert np.any(result.bound_adjusted)
    assert np.all(result.estimate <= 0.16 + 1e-12)
    assert np.any(result.raw_estimate > 0.16)          # xam dəyər DƏYİŞMƏYİB
    assert any("hədləri" in w for w in result.warnings)


# ── 2-4. PERMX / PERMY / PERMZ ────────────────────────────────────────
@pytest.mark.parametrize("name", ["PERMX", "PERMY", "PERMZ"])
def test_permeability_defaults_to_log_space(name):
    """GATE B2 — keçiricilik DEFOLT olaraq loq fəzasında."""
    strategy = resolve_strategy(name)
    assert strategy.variable_type is VariableType.LOGNORMAL
    assert strategy.transform.name == "log"
    assert strategy.back_transform is BackTransform.MEDIAN


@pytest.mark.parametrize("name", ["PERMX", "PERMY", "PERMZ"])
def test_permeability_interpolation_stays_positive_and_uses_log_space(name):
    rng = np.random.default_rng(21)
    points = _wells(50, seed=21)
    values = np.exp(4.0 + 2.0 * np.sin(points[:, 0] / 300.0) + 0.3 * rng.standard_normal(50))
    result = interpolate_property_field(points, values, _grid(), property_name=name)

    assert np.all(result.estimate > 0.0), "keçiricilik heç vaxt ≤ 0 ola bilməz"
    assert result.variance_kind is VarianceKind.EXACT
    # çevrilmiş fəza LOQ fəzasıdır — qiymətlər ln(K) miqyasındadır
    assert np.allclose(np.exp(result.transformed_estimate), result.raw_estimate,
                       rtol=1e-9)
    assert np.all(np.isfinite(result.transformed_variance))


def test_permeability_log_kriging_differs_from_raw_kriging():
    """Loq fəzası REAL fərq yaradır — sadəcə konfiqurasiya sahəsi deyil."""
    from imex2d.geology.transforms import IDENTITY_TRANSFORM
    rng = np.random.default_rng(22)
    points = _wells(45, seed=22)
    values = np.exp(rng.normal(4.0, 1.6, size=45))
    targets = _grid(10)

    log_space = interpolate_property_field(points, values, targets,
                                           property_name="PERMX")
    raw = interpolate_property_field(
        points, values, targets,
        strategy=resolve_strategy("PERMX").derive(
            transform=IDENTITY_TRANSFORM, variable_type=VariableType.CONTINUOUS,
            legacy_log_transform=False))
    assert not np.allclose(log_space.estimate, raw.estimate, rtol=1e-3)


def test_permeability_rejects_non_positive_input_through_qc():
    """Sıfır keçiricilik SƏSSİZCƏ müsbət ədədə çevrilmir — QC çıxarır."""
    points = _wells(20, seed=23)
    values = np.exp(np.linspace(2.0, 6.0, 20))
    values[3] = 0.0
    values[7] = -5.0
    result = interpolate_property_field(points, values, _grid(6), property_name="PERMX")
    assert result.quality.n_bound_violations == 2
    assert result.quality.n_valid == 18
    assert np.all(result.estimate > 0.0)


def test_lognormal_mean_back_transform_exceeds_median_everywhere():
    """B1.3 — `MEAN` modu sistematik olaraq medianın ÜSTÜNDƏDİR."""
    rng = np.random.default_rng(24)
    points = _wells(30, seed=24)
    values = np.exp(rng.normal(3.0, 1.0, size=30))
    targets = _grid(8)

    median = interpolate_property_field(points, values, targets,
                                        property_name="PERMX")
    mean = interpolate_property_field(
        points, values, targets,
        strategy=resolve_strategy("PERMX").derive(back_transform=BackTransform.MEAN))
    mean_ok = interpolate_property_field(
        points, values, targets,
        strategy=resolve_strategy("PERMX").derive(back_transform=BackTransform.MEAN_OK))

    away = ~np.isclose(median.transformed_variance, 0.0)
    assert np.all(mean.raw_estimate[away] > median.raw_estimate[away])
    assert not np.allclose(mean.raw_estimate, mean_ok.raw_estimate)


# ── 5-6. SW / NTG ─────────────────────────────────────────────────────
@pytest.mark.parametrize("name", ["SW", "NTG"])
def test_bounded_properties_use_logit_and_can_never_leave_zero_one(name):
    """GATE B3 — hədli xassələr RİYAZİ olaraq `[0,1]`-dən çıxa bilmir."""
    strategy = resolve_strategy(name)
    assert strategy.variable_type is VariableType.BOUNDED
    assert strategy.transform.name == "logit"

    rng = np.random.default_rng(31)
    points = _wells(35, seed=31)
    values = np.clip(0.5 + 0.45 * np.sin(points[:, 0] / 200.0)
                     + 0.02 * rng.standard_normal(35), 0.001, 0.999)
    result = interpolate_property_field(points, values, _grid(), property_name=name)

    assert np.all(result.raw_estimate >= 0.0) and np.all(result.raw_estimate <= 1.0)
    assert not np.any(result.bound_adjusted), (
        "logit geri çevirməsi hədləri poza bilmir → kəsməyə EHTİYAC olmamalıdır")
    assert result.variance_kind is VarianceKind.DELTA


@pytest.mark.parametrize("name", ["SW", "NTG"])
def test_bounded_properties_honour_hard_data_including_exact_bounds(name):
    points = _wells(18, seed=32)
    values = np.linspace(0.0, 1.0, 18)          # DƏQİQ 0 və 1 daxil
    result = interpolate_property_field(points, values, points, property_name=name)
    assert np.allclose(result.estimate, values, atol=1e-8)


def test_unbounded_kriging_would_violate_bounds_but_logit_does_not():
    """Müqayisə: EYNİ data ilə çevirməsiz kriginq hədləri POZUR, logit
    isə POZMUR — çevirmənin faydası ÖLÇÜLÜR (nəzəri iddia deyil).

    Qurğu kəskin keçidli (step) doyma sahəsidir: dəyərlər hər iki hədə
    yaxındır, gradient isə diktdir — adi kriginqin mənfi çəkiləri məhz
    belə həndəsədə qiyməti aralıqdan çıxarır."""
    from imex2d.geology.transforms import IDENTITY_TRANSFORM
    rng = np.random.default_rng(4)
    points = rng.uniform(0.0, 500.0, size=(30, 2))
    values = np.clip(0.02 + 0.96 / (1.0 + np.exp(-(points[:, 0] - 250.0) / 12.0)),
                     0.0, 1.0)
    axis = np.linspace(0.0, 500.0, 25)
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    targets = np.column_stack([xx.ravel(), yy.ravel()])

    plain = interpolate_property_field(
        points, values, targets,
        strategy=resolve_strategy("SW").derive(
            transform=IDENTITY_TRANSFORM, variable_type=VariableType.CONTINUOUS,
            bound_policy=BoundPolicy.FLAG, variogram_model="spherical"))
    logit = interpolate_property_field(
        points, values, targets,
        strategy=resolve_strategy("SW").derive(
            bound_policy=BoundPolicy.FLAG, variogram_model="spherical"))

    assert np.any((plain.raw_estimate < 0.0) | (plain.raw_estimate > 1.0)), (
        "çevirməsiz kriginq bu həndəsədə hədləri POZMALIDIR")
    assert np.all(logit.raw_estimate >= 0.0) and np.all(logit.raw_estimate <= 1.0)
    assert not np.any(logit.bound_adjusted)


def test_ntg_and_sw_share_the_bounded_treatment_but_are_separate_objects():
    assert resolve_strategy("NTG") is not resolve_strategy("SW")
    assert resolve_strategy("NTG").name == "NTG"
    assert resolve_strategy("SW").name == "SW"


# ── 7. FACIES ─────────────────────────────────────────────────────────
def test_facies_is_categorical_and_never_continuous():
    """GATE B4 — kateqorik kod kəsilməz kriginqdən KEÇMİR."""
    strategy = resolve_strategy("FACIES")
    assert strategy.variable_type is VariableType.CATEGORICAL
    assert strategy.interpolation is InterpolationKind.INDICATOR

    points = _wells(30, seed=41)
    codes = (points[:, 0] > 500).astype(int) + 2 * (points[:, 1] > 500).astype(int)
    with pytest.raises(ValueError, match="KATEQORİK"):
        interpolate_property_field(points, codes, _grid(), property_name="FACIES")


def test_facies_returns_valid_probability_distributions():
    points = _wells(40, seed=42)
    codes = (points[:, 0] > 500).astype(int) + 2 * (points[:, 1] > 500).astype(int)
    result = interpolate_categorical_field(points, codes, _grid(),
                                           property_name="FACIES")

    assert isinstance(result, CategoricalEstimate)
    assert np.all(result.probabilities >= 0.0) and np.all(result.probabilities <= 1.0)
    assert np.allclose(result.probabilities.sum(axis=1), 1.0, atol=1e-12)
    assert set(result.most_probable.tolist()).issubset(set(result.categories.tolist()))
    # ən ehtimallı kod HƏQİQƏTƏN ən böyük ehtimallı sütundur
    assert np.array_equal(result.most_probable,
                          result.categories[np.argmax(result.probabilities, axis=1)])


def test_facies_never_produces_fractional_category_codes():
    points = _wells(30, seed=43)
    codes = np.where(points[:, 0] > 500, 3, 1)
    result = interpolate_categorical_field(points, codes, _grid(),
                                           property_name="FACIES")
    assert set(np.unique(result.most_probable).tolist()).issubset({1, 3})
    assert result.most_probable.dtype.kind in "iu"


def test_facies_entropy_is_low_near_data_and_high_at_boundaries():
    """Entropiya ƏSL qeyri-müəyyənlik ölçüsüdür: quyunun üstündə ~0,
    iki fasiya arasındakı sərhəddə yüksək."""
    points = np.array([[100., 500.], [150., 500.], [850., 500.], [900., 500.],
                       [100., 200.], [900., 200.], [100., 800.], [900., 800.]])
    codes = np.array([1, 1, 2, 2, 1, 2, 1, 2])
    targets = np.array([[100., 500.], [500., 500.], [900., 500.]])
    result = interpolate_categorical_field(points, codes, targets,
                                           property_name="FACIES")
    assert result.entropy[0] < result.entropy[1]
    assert result.entropy[2] < result.entropy[1]
    assert np.all(result.normalized_entropy >= 0.0)
    assert np.all(result.normalized_entropy <= 1.0 + 1e-12)


def test_facies_probability_corrections_are_counted():
    points = _wells(25, seed=44)
    codes = (points[:, 1] > 400).astype(int)
    result = interpolate_categorical_field(points, codes, _grid(),
                                           property_name="FACIES")
    assert result.n_probability_corrections >= 0
    if result.n_probability_corrections:
        assert any("normallaşdırıldı" in w for w in result.warnings)


def test_single_category_is_certain():
    points = _wells(15, seed=45)
    codes = np.full(15, 2)
    result = interpolate_categorical_field(points, codes, _grid(5),
                                           property_name="FACIES")
    assert np.all(result.most_probable == 2)
    assert np.allclose(result.entropy, 0.0)
    assert np.allclose(result.max_probability, 1.0)


def test_facies_probability_of_helper_matches_the_matrix():
    points = _wells(30, seed=46)
    codes = (points[:, 0] > 500).astype(int)
    result = interpolate_categorical_field(points, codes, _grid(6),
                                           property_name="FACIES")
    assert np.array_equal(result.probability_of(1), result.probabilities[:, 1])


# ── strategiya reyestri (B1.7) ────────────────────────────────────────
def test_registry_covers_every_required_property():
    for name in ("PORO", "PERMX", "PERMY", "PERMZ", "SW", "NTG", "FACIES"):
        assert name in DEFAULT_STRATEGIES
        assert resolve_strategy(name).name == name


def test_unknown_property_gets_a_neutral_continuous_strategy():
    strategy = resolve_strategy("SOME_NEW_LOG")
    assert strategy.variable_type is VariableType.CONTINUOUS
    assert strategy.transform.is_identity
    assert strategy.bound_policy is BoundPolicy.NONE


def test_unknown_categorical_name_still_routes_to_indicator():
    """`property_types` reyestrindəki kateqorik ad kəsilməz yola DÜŞMÜR."""
    assert resolve_strategy("LITHOLOGY").is_categorical
    assert resolve_strategy("ROCKTYPE").interpolation is InterpolationKind.INDICATOR


def test_strategies_are_immutable_and_derive_returns_a_copy():
    base = resolve_strategy("PERMX")
    derived = base.derive(variogram_model="spherical")
    assert derived is not base
    assert base.variogram_model == "auto" and derived.variogram_model == "spherical"
    with pytest.raises(Exception):
        base.variogram_model = "gaussian"          # frozen dataclass


@pytest.mark.parametrize("kwargs,match", [
    ({"variogram_model": "linear"}, "variogram"),
    ({"min_neighbors": 0}, "min_neighbors"),
    ({"max_neighbors": 1, "min_neighbors": 4}, "max_neighbors"),
    ({"honor_hard_data": "sometimes"}, "honor_hard_data"),
])
def test_invalid_strategy_configuration_is_rejected(kwargs, match):
    with pytest.raises(PropertyConfigError, match=match):
        resolve_strategy("PORO").derive(**kwargs)


def test_lognormal_strategy_refuses_silent_identity_transform():
    """GATE B2 sərhədi: loq-normal xassəni sükutla xam fəzaya keçirmək OLMAZ."""
    from imex2d.geology.transforms import IDENTITY_TRANSFORM
    with pytest.raises(PropertyConfigError, match="GATE B2"):
        resolve_strategy("PERMX").derive(transform=IDENTITY_TRANSFORM,
                                         legacy_log_transform=False)


def test_categorical_strategy_refuses_continuous_interpolation():
    with pytest.raises(PropertyConfigError, match="GATE B4"):
        resolve_strategy("FACIES").derive(interpolation=InterpolationKind.KRIGING)


def test_normal_score_strategy_keeps_bounds_but_swaps_the_transform():
    base = resolve_strategy("PERMX")
    gaussian = normal_score_strategy(base)
    assert gaussian.transform.name == "normal_score"
    assert gaussian.physical_bounds == base.physical_bounds
    assert gaussian.output_bounds == base.output_bounds


def test_normal_score_strategy_rejects_categorical_properties():
    with pytest.raises(PropertyConfigError, match="SIS"):
        normal_score_strategy(resolve_strategy("FACIES"))


def test_strategy_table_lists_the_required_properties():
    table = strategy_table()
    for name in ("PORO", "PERMX", "SW", "NTG", "FACIES"):
        assert name in table
    assert "lognormal" in table and "bounded" in table and "categorical" in table


def test_interpolate_by_name_routes_categorical_and_continuous_differently():
    points = _wells(25, seed=51)
    continuous = interpolate_by_name(points, 0.2 + 0.05 * np.sin(points[:, 0] / 100.0),
                                     _grid(6), "PORO")
    categorical = interpolate_by_name(points, (points[:, 0] > 500).astype(int),
                                      _grid(6), "FACIES")
    assert isinstance(continuous, PropertyEstimate)
    assert isinstance(categorical, CategoricalEstimate)


def test_every_property_produces_a_confidence_label():
    points = _wells(30, seed=52)
    for name, values in (("PORO", 0.2 + 0.04 * np.sin(points[:, 0] / 200.0)),
                         ("PERMX", np.exp(4 + np.cos(points[:, 1] / 150.0))),
                         ("SW", np.clip(0.4 + 0.2 * np.sin(points[:, 0] / 170.0), 0, 1))):
        result = interpolate_property_field(points, values, _grid(6),
                                            property_name=name)
        labels = set(np.asarray(result.confidence).astype(str).tolist())
        assert labels.issubset({c.value for c in Confidence})
