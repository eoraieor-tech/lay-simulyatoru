"""A1 — istehsal səviyyəli Kriging özəyi: riyaziyyat + ədədi dayanıqlıq.

Bu fayl `OrdinaryKriging`-in YENİ (A1) müqaviləsini yoxlayır:
`krige()` → `KrigingResult` (qiymət, varians, qonşu sayı, ən yaxın
məsafə, ekstrapolyasiya bayrağı, solver statusu), yerli sistemlərin
ədədi dayanıqlığı (təkil/pis şərtlənmiş matris, dublikat koordinat,
NaN/±inf) və yansızlıq şərti `Σwᵢ = 1`.

Mövcud (Phase 2-5) davranış `test_kriging_3d_anisotropy.py` /
`test_kriging_variance.py` ilə qorunur — burada onlar TƏKRARLANMIR.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.geology.interpolation import (SOLVER_DIRECT, SOLVER_EXACT, SOLVER_NONE,
                                          UNBIASED_TOLERANCE, KrigingResult,
                                          OrdinaryKriging)
from imex2d.geology.spatial_search import NeighborhoodConfig

POINTS = np.array([[0., 0.], [100., 0.], [0., 100.], [100., 100.], [50., 50.]])
VALUES = np.array([0.15, 0.25, 0.20, 0.30, 0.22])


def _grid(n=12, high=100.0):
    axis = np.linspace(0.0, high, n)
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


def _scattered(n, seed, high=1000.0, ndim=2):
    rng = np.random.default_rng(seed)
    return rng.uniform(0.0, high, size=(n, ndim))


# ── 1. xətti sahənin bərpası ───────────────────────────────────────────
def test_linear_field_is_reproduced_within_tolerance():
    """`Z = 2X + 3Y + 10` — radius domendən böyük olanda variogram
    başlanğıcda demək olar XƏTTİDİR, ona görə adi kriging xətti trendi
    kiçik xəta ilə bərpa etməlidir (A8.1/A9 Test 1)."""
    points = _scattered(120, seed=11, high=1000.0)
    values = 2.0 * points[:, 0] + 3.0 * points[:, 1] + 10.0
    targets = _scattered(60, seed=12, high=800.0) + 100.0

    kriging = OrdinaryKriging(range_=5000.0, sill=1.0, nugget=0.0, max_neighbors=16)
    estimate = kriging.interpolate(points, values, targets)
    truth = 2.0 * targets[:, 0] + 3.0 * targets[:, 1] + 10.0
    relative = np.abs(estimate - truth) / (truth.max() - truth.min())
    assert np.all(np.isfinite(estimate))
    assert relative.max() < 0.02, f"maksimum nisbi xəta {relative.max():.4f}"


def test_linear_field_in_3d_is_reproduced():
    """`Z = 2X + 3Y + 4Z + 10` (A9 Test 1 — 3D forması)."""
    points = _scattered(150, seed=21, high=500.0, ndim=3)
    values = 2.0 * points[:, 0] + 3.0 * points[:, 1] + 4.0 * points[:, 2] + 10.0
    targets = _scattered(40, seed=22, high=400.0, ndim=3) + 50.0

    kriging = OrdinaryKriging(range_=5000.0, range_v=5000.0, sill=1.0, nugget=0.0,
                              max_neighbors=20)
    estimate = kriging.interpolate(points, values, targets)
    truth = 2.0 * targets[:, 0] + 3.0 * targets[:, 1] + 4.0 * targets[:, 2] + 10.0
    relative = np.abs(estimate - truth) / (truth.max() - truth.min())
    assert relative.max() < 0.05, f"maksimum nisbi xəta {relative.max():.4f}"


# ── 2. sabit sahə ──────────────────────────────────────────────────────
def test_constant_field_is_reproduced_exactly():
    """Bütün dəyərlər eynidirsə `Σwᵢ=1` şərti nəticəni DƏYİŞMƏZ saxlayır."""
    points = _scattered(60, seed=3)
    values = np.full(60, 0.23)
    targets = _scattered(30, seed=4)
    for kriging in (OrdinaryKriging(range_=300.0, sill=0.01, nugget=0.0),
                    OrdinaryKriging(range_=300.0, sill=0.01, nugget=0.002,
                                    max_neighbors=8),
                    OrdinaryKriging(range_=300.0, sill=0.01, nugget=0.0,
                                    search_radius=400.0, max_neighbors=12)):
        estimate = kriging.interpolate(points, values, targets)
        assert np.allclose(estimate[np.isfinite(estimate)], 0.23, atol=1e-9)


# ── 3. sərt datanın DƏQİQ honor edilməsi ───────────────────────────────
def test_exact_hard_data_is_honoured_global_and_local():
    for kriging in (OrdinaryKriging(nugget=0.0),
                    OrdinaryKriging(nugget=0.0, max_neighbors=3),
                    OrdinaryKriging(nugget=0.0, search_radius=200.0, max_neighbors=4),
                    OrdinaryKriging(nugget=0.0, max_neighbors=4, sectors=4)):
        result = kriging.krige(POINTS, VALUES, POINTS)
        assert np.allclose(result.estimate, VALUES, atol=1e-12)
        assert np.allclose(result.variance, 0.0, atol=1e-12)
        assert np.all(result.solver.astype(str) == SOLVER_EXACT)


def test_hard_data_policy_never_and_always():
    """Siyasət AÇIQDIR: nugget > 0 olanda defolt honor ETMİR (ölçmə
    səhvi süzülür), `always` isə məcbur edir."""
    filtered = OrdinaryKriging(nugget=0.05, range_=200.0, sill=0.01)
    honoured = OrdinaryKriging(nugget=0.05, range_=200.0, sill=0.01,
                               honor_hard_data="always")
    never = OrdinaryKriging(nugget=0.0, range_=200.0, sill=0.01,
                            honor_hard_data="never")
    assert not np.allclose(filtered.interpolate(POINTS, VALUES, POINTS), VALUES,
                           atol=1e-6)
    assert np.allclose(honoured.interpolate(POINTS, VALUES, POINTS), VALUES, atol=1e-12)
    # never: nugget=0 olsa belə sistemin öz nəticəsi qalır (praktikada
    # yenə çox yaxındır, amma DƏQİQ bərpa MƏCBUR EDİLMİR)
    assert never.interpolate(POINTS, VALUES, POINTS).shape == VALUES.shape


def test_conflicting_duplicates_are_averaged_deterministically():
    points = np.vstack([POINTS, [[50., 50.]]])
    values = np.append(VALUES, 0.40)          # (50,50)-də 0.22 və 0.40
    kriging = OrdinaryKriging(nugget=0.0)
    result = kriging.krige(points, values, np.array([[50., 50.]]))
    assert result.estimate[0] == pytest.approx(0.31)      # (0.22+0.40)/2
    assert any("ziddiyyətli" in w for w in result.warnings)


# ── 4-5. sıfır məsafə / dublikat koordinat ────────────────────────────
def test_zero_distance_target_matches_data_point():
    kriging = OrdinaryKriging(range_=200.0, sill=0.01, nugget=0.0, max_neighbors=4)
    result = kriging.krige(POINTS, VALUES, np.array([[50., 50.], [50. + 1e-15, 50.]]))
    assert np.allclose(result.estimate, 0.22, atol=1e-9)
    assert np.allclose(result.variance, 0.0, atol=1e-9)


def test_duplicate_coordinates_do_not_break_the_solver():
    """Eyni koordinatda EYNİ dəyərli təkrarlar səssizcə birləşdirilir —
    sistem təkil olmur, xəbərdarlıq da yaranmır (ziddiyyət yoxdur)."""
    points = np.vstack([POINTS, POINTS[:3]])
    values = np.append(VALUES, VALUES[:3])
    kriging = OrdinaryKriging(range_=200.0, sill=0.01, nugget=0.0, max_neighbors=6)
    result = kriging.krige(points, values, np.array([[25., 25.], [0., 0.]]))
    assert np.all(np.isfinite(result.estimate))
    assert result.estimate[1] == pytest.approx(0.15)
    assert not any("ziddiyyətli" in w for w in result.warnings)


# ── 6-7. deqenerativ həndəsə / pis şərtlənmiş sistem ──────────────────
def test_collinear_points_are_solved_deterministically():
    """Bütün nöqtələr bir xətt üzərində — həndəsə deqenerativdir, amma
    kriging sistemi (variogram matrisi) hələ də həll edilə bilər."""
    points = np.column_stack([np.linspace(0.0, 500.0, 21), np.zeros(21)])
    values = np.sin(points[:, 0] / 120.0)
    kriging = OrdinaryKriging(range_=250.0, sill=0.5, nugget=0.0, max_neighbors=8)
    first = kriging.krige(points, values, np.array([[123.0, 0.0], [123.0, 40.0]]))
    second = kriging.krige(points, values, np.array([[123.0, 0.0], [123.0, 40.0]]))
    assert np.all(np.isfinite(first.estimate))
    assert np.array_equal(first.estimate, second.estimate)


def test_nearly_coincident_points_trigger_documented_fallback_not_nan():
    """1e-12 məsafədə (dublikat kimi görünməyən, amma matrisi demək olar
    təkil edən) nöqtələr — solver ehtiyat yolu ilə SONLU nəticə verir."""
    eps = 1e-12
    points = np.array([[0., 0.], [eps, 0.], [0., eps], [100., 0.], [0., 100.]])
    values = np.array([0.10, 0.10 + 1e-9, 0.10 - 1e-9, 0.30, 0.20])
    kriging = OrdinaryKriging(range_=200.0, sill=0.01, nugget=0.0, max_neighbors=5)
    result = kriging.krige(points, values, np.array([[10., 10.]]))
    assert np.isfinite(result.estimate[0])
    assert result.solver[0] in ("direct", "jitter", "lstsq", "renormalized",
                                "idw_fallback")


def test_extreme_anisotropy_ratio_stays_finite():
    """1:10⁶ üfüqi anizotropluq — matris çox pis şərtlənir, nəticə yenə
    SONLU olmalıdır (ehtiyat yolları işə düşür, NaN qaytarılmır)."""
    points = _scattered(40, seed=7, high=500.0)
    values = points[:, 0] * 0.001
    kriging = OrdinaryKriging(range_=500.0, range_minor=5e-4, azimuth_deg=30.0,
                              sill=0.01, nugget=0.0, max_neighbors=10)
    result = kriging.krige(points, values, _scattered(20, seed=8, high=500.0))
    assert np.all(np.isfinite(result.estimate))


# ── 8. qonşu çatışmazlığı ─────────────────────────────────────────────
def test_insufficient_neighbors_returns_nan_not_invented_value():
    kriging = OrdinaryKriging(nugget=0.0, search_radius=5.0, min_neighbors=3)
    result = kriging.krige(POINTS, VALUES, np.array([[500., 500.], [50., 50.]]))
    assert np.isnan(result.estimate[0]) and np.isnan(result.variance[0])
    assert result.solver[0] == SOLVER_NONE
    assert result.neighbor_count[0] == 0
    assert result.extrapolated[0]


def test_min_neighbors_is_enforced_even_when_some_points_are_in_range():
    """Radiusda 1 nöqtə var, amma `min_neighbors=3` — dəyər UYDURULMUR."""
    kriging = OrdinaryKriging(nugget=0.0, search_radius=20.0, min_neighbors=3)
    result = kriging.krige(POINTS, VALUES, np.array([[10., 0.]]))
    assert np.isnan(result.estimate[0])


# ── 9. NaN / ±inf ─────────────────────────────────────────────────────
def test_non_finite_hard_data_is_dropped_with_explicit_warning():
    points = np.vstack([POINTS, [[np.nan, 0.]], [[np.inf, 10.]]])
    values = np.append(VALUES, [0.5, 0.6])
    kriging = OrdinaryKriging(range_=200.0, sill=0.01, nugget=0.0)
    result = kriging.krige(points, values, np.array([[25., 25.]]))
    assert np.isfinite(result.estimate[0])
    assert any("ÇIXARILDI" in w for w in result.warnings)


def test_non_finite_value_is_dropped_too():
    values = VALUES.copy()
    values[2] = np.nan
    kriging = OrdinaryKriging(range_=200.0, sill=0.01, nugget=0.0)
    result = kriging.krige(POINTS, values, np.array([[25., 25.]]))
    assert np.isfinite(result.estimate[0])


def test_drop_non_finite_false_raises_instead_of_silently_filtering():
    points = np.vstack([POINTS, [[np.nan, 0.]]])
    values = np.append(VALUES, 0.5)
    kriging = OrdinaryKriging(range_=200.0, sill=0.01, drop_non_finite=False)
    with pytest.raises(ValueError, match="NaN"):
        kriging.interpolate(points, values, np.array([[25., 25.]]))


def test_non_finite_target_yields_nan_without_breaking_others():
    kriging = OrdinaryKriging(range_=200.0, sill=0.01, nugget=0.0)
    targets = np.array([[25., 25.], [np.nan, 10.], [75., 75.]])
    result = kriging.krige(POINTS, VALUES, targets)
    assert np.isfinite(result.estimate[0]) and np.isfinite(result.estimate[2])
    assert np.isnan(result.estimate[1])


def test_all_hard_data_non_finite_returns_nan_and_says_so():
    points = np.array([[np.nan, 0.], [np.inf, 1.]])
    values = np.array([0.1, 0.2])
    result = OrdinaryKriging(range_=100.0, sill=0.01).krige(
        points, values, np.array([[0., 0.]]))
    assert np.isnan(result.estimate[0])
    assert any("etibarlı sərt data" in w for w in result.warnings)


def test_mismatched_lengths_raise():
    with pytest.raises(ValueError, match="uzunluğu uyğun gəlmir"):
        OrdinaryKriging().interpolate(POINTS, VALUES[:3], POINTS)


# ── 10. çəkilər cəmi = 1 ──────────────────────────────────────────────
def test_kriging_weights_sum_to_one():
    """Yansızlıq şərti sistemin İÇİNDƏ qurulub; onu birbaşa yoxlamaq üçün
    SABİT sahədən istifadə edirik: `Σwᵢ=1` ⇔ sabit sahə dəqiq bərpa olunur
    (yuxarıdakı test) — burada isə solver statusunun heç bir hədəfdə
    "yenidən normallanmış"a düşmədiyini təsdiqləyirik, yəni orijinal həll
    `UNBIASED_TOLERANCE` daxilində Σw=1 verib."""
    points = _scattered(80, seed=31, high=600.0)
    values = np.sin(points[:, 0] / 90.0) + np.cos(points[:, 1] / 70.0)
    targets = _scattered(50, seed=32, high=600.0)
    result = OrdinaryKriging(range_=250.0, sill=1.0, nugget=0.01,
                             max_neighbors=12).krige(points, values, targets)
    assert np.all(result.solver.astype(str) == SOLVER_DIRECT)
    assert UNBIASED_TOLERANCE > 0.0


def test_estimate_is_a_convex_like_combination_of_neighbour_values():
    """Σwᵢ=1 olduğundan qiymət heç vaxt qonşu dəyərlərin diapazonundan
    ÇOX uzaqda ola bilməz — kəskin ekstrapolyasiyanın olmaması yoxlanılır."""
    points = _grid(8, 700.0)
    values = 0.1 + 0.0002 * points[:, 0]
    targets = _scattered(40, seed=41, high=700.0)
    result = OrdinaryKriging(range_=300.0, sill=0.01, nugget=0.0,
                             max_neighbors=9).krige(points, values, targets)
    span = values.max() - values.min()
    assert np.all(result.estimate >= values.min() - 0.25 * span)
    assert np.all(result.estimate <= values.max() + 0.25 * span)


# ── 11. kriging variansı ──────────────────────────────────────────────
def test_variance_is_non_negative_and_zero_at_data():
    points = _scattered(50, seed=51, high=500.0)
    values = np.sin(points[:, 0] / 80.0)
    targets = np.vstack([points[:5], _scattered(30, seed=52, high=500.0)])
    result = OrdinaryKriging(range_=200.0, sill=1.0, nugget=0.0,
                             max_neighbors=10).krige(points, values, targets)
    assert np.all(result.variance[np.isfinite(result.variance)] >= -1e-12)
    assert np.allclose(result.variance[:5], 0.0, atol=1e-9)


def test_variance_grows_away_from_data_and_saturates_below_total_sill():
    """Uzaqda kriging variansı a-priori sillə YAXINLAŞIR (adi kriging-də
    Laqranj vuruğu səbəbindən onu bir qədər KEÇƏ bilər, ona görə hədd
    `2·(nugget+sill)` kimi qoyulur — sonsuz böyümə OLMAMALIDIR)."""
    kriging = OrdinaryKriging(range_=100.0, sill=0.02, nugget=0.001, max_neighbors=5)
    result = kriging.krige(POINTS, VALUES,
                           np.array([[50., 50.], [150., 150.], [4000., 4000.]]))
    assert result.variance[0] < result.variance[1] < result.variance[2]
    assert result.variance[2] <= 2.0 * (0.02 + 0.001)


def test_variance_matches_manual_lagrange_formula():
    """`σ² = Σ wᵢ·γ(xᵢ,x₀) + μ` — sistemi ƏLLƏ qurub müqayisə edirik
    (A1.6: variogram konvensiyası ilə tam uyğunluq)."""
    from imex2d.geology.variogram import spherical
    points = POINTS[:4]
    values = VALUES[:4]
    target = np.array([[30.0, 40.0]])
    range_, sill, nugget = 150.0, 0.02, 0.0

    n = points.shape[0]
    left = np.ones((n + 1, n + 1))
    diff = points[:, None, :] - points[None, :, :]
    hh = np.sqrt(np.sum(diff * diff, axis=-1))
    left[:n, :n] = np.where(hh <= 1e-12, 0.0, spherical(hh, nugget, sill, range_))
    left[n, n] = 0.0
    right = np.ones(n + 1)
    h0 = np.linalg.norm(points - target, axis=1)
    right[:n] = spherical(h0, nugget, sill, range_)
    solution = np.linalg.solve(left, right)
    expected_estimate = float(np.dot(solution[:n], values))
    expected_variance = float(np.dot(solution[:n], right[:n]) + solution[n])

    result = OrdinaryKriging(range_=range_, sill=sill, nugget=nugget).krige(
        points, values, target)
    assert result.estimate[0] == pytest.approx(expected_estimate, rel=1e-10)
    assert result.variance[0] == pytest.approx(expected_variance, rel=1e-10)


def test_kriging_matrix_is_symmetric_by_construction():
    """`Γ[i,j] = γ(hᵢⱼ) = γ(hⱼᵢ)` — simmetriya variogramın özündən gəlir."""
    kriging = OrdinaryKriging(range_=200.0, sill=0.01, nugget=0.003)
    diff = POINTS[:, None, :] - POINTS[None, :, :]
    h = np.sqrt(np.sum(diff * diff, axis=-1))
    gamma = kriging._variogram(h, 200.0, 0.01, 0.003, "spherical")
    assert np.allclose(gamma, gamma.T, atol=0.0)
    assert np.allclose(np.diag(gamma), 0.0, atol=0.0)


# ── 12. determinizm ───────────────────────────────────────────────────
def test_repeated_calls_are_bitwise_identical():
    points = _scattered(90, seed=61, high=800.0)
    values = np.cos(points[:, 1] / 60.0)
    targets = _scattered(70, seed=62, high=800.0)
    kriging = OrdinaryKriging(range_=250.0, sill=1.0, nugget=0.0, max_neighbors=12,
                              sectors=4)
    first = kriging.interpolate(points, values, targets)
    second = kriging.interpolate(points, values, targets)
    assert np.array_equal(first, second)


def test_batched_and_single_target_calls_agree():
    """Hədəflər BİRLİKDƏ verilsə də, TƏK-TƏK verilsə də nəticə eynidir —
    toplu (batched) həll heç nəyi dəyişmir."""
    points = _scattered(40, seed=71, high=400.0)
    values = np.sin(points[:, 0] / 50.0)
    targets = _scattered(9, seed=72, high=400.0)
    kriging = OrdinaryKriging(range_=180.0, sill=1.0, nugget=0.0, max_neighbors=7)
    batched = kriging.interpolate(points, values, targets)
    single = np.array([kriging.interpolate(points, values, targets[i:i + 1])[0]
                       for i in range(targets.shape[0])])
    assert np.allclose(batched, single, atol=1e-12)


# ── nəticə obyekti + geriyə uyğunluq ──────────────────────────────────
def test_kriging_result_exposes_required_fields():
    result = OrdinaryKriging(range_=200.0, sill=0.01, nugget=0.0,
                             max_neighbors=3).krige(POINTS, VALUES,
                                                    np.array([[25., 25.]]))
    assert isinstance(result, KrigingResult)
    assert result.neighbor_count[0] == 3
    assert np.isfinite(result.nearest_distance[0])
    assert result.support[0] in ("well_supported", "boundary", "weak", "extrapolated")
    assert result.local is True
    assert np.asarray(result).shape == (1,)
    assert len(result) == 1
    assert "Kriging" in result.summary()


def test_interpolate_and_krige_agree():
    kriging = OrdinaryKriging(range_=200.0, sill=0.01, nugget=0.01, max_neighbors=4)
    targets = _grid(5, 100.0)
    plain = kriging.interpolate(POINTS, VALUES, targets)
    with_variance, variance = kriging.interpolate_with_variance(POINTS, VALUES, targets)
    result = kriging.krige(POINTS, VALUES, targets)
    assert np.allclose(plain, with_variance, atol=0.0, equal_nan=True)
    assert np.allclose(plain, result.estimate, atol=0.0, equal_nan=True)
    assert np.allclose(variance, result.variance, atol=0.0, equal_nan=True)


# ── Gate 3: qlobal sıx sistem defolt istehsal yolu DEYİL ──────────────
def test_large_dataset_switches_to_local_system_automatically():
    points = _scattered(400, seed=81, high=2000.0)
    values = np.sin(points[:, 0] / 200.0)
    targets = _scattered(50, seed=82, high=2000.0)
    kriging = OrdinaryKriging(range_=400.0, sill=1.0, nugget=0.0)
    result = kriging.krige(points, values, targets)
    assert result.local is True
    assert result.neighbor_count.max() <= kriging.auto_local_max_neighbors
    assert any("auto_local_threshold" in w for w in result.warnings)


def test_small_dataset_still_uses_the_global_system():
    result = OrdinaryKriging(range_=200.0, sill=0.01, nugget=0.0).krige(
        POINTS, VALUES, np.array([[25., 25.]]))
    assert result.local is False


def test_explicit_neighborhood_config_overrides_shortcuts():
    config = NeighborhoodConfig(min_neighbors=2, max_neighbors=3, search_radius=90.0,
                                support_range=200.0)
    kriging = OrdinaryKriging(range_=200.0, sill=0.01, nugget=0.0, neighborhood=config)
    result = kriging.krige(POINTS, VALUES, np.array([[50., 50.]]))
    assert result.local is True
    assert result.neighbor_count[0] == 3
