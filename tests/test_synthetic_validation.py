"""A9 — nəzarət altındakı SİNTETİK doğrulama.

Beş ssenari, hər biri MƏLUM həqiqətlə:

    Test 1  xətti sahə            Z = 2X + 3Y + 4Z + 10
    Test 2  izotrop qauss sahəsi  məlum korrelyasiya radiusu
    Test 3  anizotrop sahə        major ≫ minor ≫ şaquli
    Test 4  dönmüş anizotropluq   baş istiqamət 45°
    Test 5  laylı rezervuar       güclü üfüqi, zəif şaquli davamlılıq

Sahələr kovariasiya matrisinin Xolesskiy parçalanması ilə qurulur —
"korrelyasiyalı görünən" səs-küy DEYİL, radiusları DƏQİQ məlum olan
qauss sahəsidir; ona görə variogram fitinin/interpolyasiyanın nə qədər
düz işlədiyi ÖLÇÜLƏ bilir.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.geology.interpolation import OrdinaryKriging
from imex2d.geology.variogram import (MODEL_EXPONENTIAL, detect_anisotropy,
                                      experimental_variogram, fit_variogram,
                                      fit_variogram_from_data, vertical_variogram)


def _gaussian_field(points: np.ndarray, ranges, seed: int) -> np.ndarray:
    """`exp(−3·d_ani)` kovariasiyalı qauss sahəsi — layihənin EKSPONENSİAL
    variogram modeli ilə eyni ailədən (`γ = sill·(1 − exp(−3h/a))`),
    ona görə fit edilən radius nəzəri radiusla müqayisə edilə bilər."""
    rng = np.random.default_rng(seed)
    n = points.shape[0]
    scaled = points / np.asarray(ranges, float)[None, :points.shape[1]]
    diff = scaled[:, None, :] - scaled[None, :, :]
    distance = np.sqrt(np.sum(diff * diff, axis=-1))
    cov = np.exp(-3.0 * distance) + 1e-8 * np.eye(n)
    return np.linalg.cholesky(cov) @ rng.standard_normal(n)


def _rotate(points_xy: np.ndarray, degrees: float) -> np.ndarray:
    """(X,Y)-i saat əqrəbinin ƏKSİNƏ `degrees` qədər döndərir."""
    theta = np.radians(degrees)
    c, s = np.cos(theta), np.sin(theta)
    return points_xy @ np.array([[c, s], [-s, c]])


def _loo_rmse(kriging, points, values):
    errors = []
    mask = np.ones(points.shape[0], dtype=bool)
    for i in range(points.shape[0]):
        mask[:] = True
        mask[i] = False
        predicted = kriging.interpolate(points[mask], values[mask], points[i:i + 1])[0]
        if np.isfinite(predicted):
            errors.append(predicted - values[i])
    return float(np.sqrt(np.mean(np.square(errors))))


# ══ Test 1 — xətti sahə ═══════════════════════════════════════════════
def test_1_linear_field_is_interpolated_and_honoured():
    """`Z = 2X + 3Y + 4Z + 10`. Adi kriging (radius domendən böyük)
    xətti trendi kiçik nisbi xəta ilə bərpa etməli, sərt data
    nöqtələrini isə DƏQİQ saxlamalıdır."""
    rng = np.random.default_rng(101)
    points = np.column_stack([rng.uniform(0.0, 1000.0, size=(200, 2)),
                              rng.uniform(0.0, 300.0, 200)])
    values = 2 * points[:, 0] + 3 * points[:, 1] + 4 * points[:, 2] + 10.0
    targets = np.column_stack([rng.uniform(100.0, 900.0, size=(120, 2)),
                               rng.uniform(50.0, 250.0, 120)])
    truth = 2 * targets[:, 0] + 3 * targets[:, 1] + 4 * targets[:, 2] + 10.0

    kriging = OrdinaryKriging(range_=1e5, range_v=1e5, sill=1.0, nugget=0.0,
                              max_neighbors=24)
    result = kriging.krige(points, values, targets)
    relative = np.abs(result.estimate - truth) / (truth.max() - truth.min())
    assert relative.max() < 0.02, f"maksimum nisbi xəta {relative.max():.4f}"
    assert np.allclose(kriging.interpolate(points, values, points), values, atol=1e-6)


def test_1_linear_field_variance_is_smallest_where_data_is_densest():
    rng = np.random.default_rng(102)
    points = rng.uniform(0.0, 500.0, size=(60, 2))
    values = 2 * points[:, 0] + 3 * points[:, 1] + 10.0
    kriging = OrdinaryKriging(range_=600.0, sill=1e5, nugget=0.0, max_neighbors=12)
    result = kriging.krige(points, values, np.array([[250., 250.], [1200., 1200.]]))
    assert result.variance[0] < result.variance[1]
    assert result.support[0] == "well_supported"
    assert result.extrapolated[1]


# ══ Test 2 — izotrop qauss sahəsi ═════════════════════════════════════
def test_2_isotropic_field_variogram_recovers_the_known_range():
    """Sahə `range = 250` ilə qurulub — fit edilən praktiki radius eyni
    tərtibdə olmalıdır (dəqiq bərabərlik gözlənilmir: sonlu nümunə +
    `max_lag` kəsimi radiusu bir qədər aşağı çəkir, bu, geostatistikanın
    məlum davranışıdır və GİZLƏDİLMİR)."""
    rng = np.random.default_rng(201)
    points = rng.uniform(0.0, 1500.0, size=(350, 2))
    values = _gaussian_field(points, (250.0, 250.0), seed=202)

    fit = fit_variogram_from_data(points, values, n_lags=14, max_lag=900.0,
                                  model=MODEL_EXPONENTIAL)
    assert 60.0 < fit.range_ < 900.0
    assert fit.sill == pytest.approx(np.var(values), rel=0.8)
    assert fit.nugget < 0.5 * fit.total_sill, "qurulmuş sahədə nugget kiçik olmalıdır"


def test_2_isotropic_field_variogram_is_direction_independent():
    """İzotrop sahədə istiqamətli radiuslar bir-birinə yaxın olmalıdır."""
    rng = np.random.default_rng(203)
    points = rng.uniform(0.0, 1500.0, size=(400, 2))
    values = _gaussian_field(points, (300.0, 300.0), seed=204)
    detection = detect_anisotropy(points, values, n_directions=6, n_lags=8)
    ranges = np.array(list(detection.directional_ranges.values()))
    assert ranges.min() > 0.0
    assert ranges.max() / ranges.min() < 3.0, "izotrop sahədə kəskin anizotropluq olmamalı"


def test_2_isotropic_kriging_is_better_than_nearest_neighbour():
    from imex2d.geology.interpolation import NearestNeighbour
    rng = np.random.default_rng(205)
    points = rng.uniform(0.0, 1000.0, size=(120, 2))
    values = _gaussian_field(points, (300.0, 300.0), seed=206)
    kriging = OrdinaryKriging(range_=300.0, sill=1.0, nugget=0.0,
                              model=MODEL_EXPONENTIAL, max_neighbors=16)
    assert _loo_rmse(kriging, points, values) < _loo_rmse(NearestNeighbour(),
                                                          points, values)


def test_2_kriging_variance_grows_with_distance_from_the_data_cloud():
    rng = np.random.default_rng(207)
    points = rng.uniform(0.0, 600.0, size=(80, 2))
    values = _gaussian_field(points, (200.0, 200.0), seed=208)
    kriging = OrdinaryKriging(range_=200.0, sill=1.0, nugget=0.0,
                              model=MODEL_EXPONENTIAL, max_neighbors=12)
    targets = np.array([[300., 300.], [700., 300.], [1200., 300.]])
    variance = kriging.krige(points, values, targets).variance
    assert variance[0] < variance[1] < variance[2]


# ══ Test 3 — anizotrop sahə (major ≫ minor ≫ şaquli) ═════════════════
def _anisotropic_3d_field(n=260, seed=301, ranges=(900.0, 150.0, 20.0),
                          depth_span=25.0):
    """major(X) ≫ minor(Y) ≫ şaquli(Z).

    `depth_span` qəsdən şaquli radiusdan bir qədər KİÇİKDİR: əks halda
    nöqtələr şaquli olaraq TAM dekorrelyasiya olur və (X,Y) müstəvisinə
    proyeksiya edilmiş variogram sırf nugget görünür — yəni üfüqi
    anizotropluq PRİNSİPCƏ ölçülə bilməz olur. Bu, kodun yox,
    NÜMUNƏLƏMƏNİN həddidir və qəsdən qeyd olunur."""
    rng = np.random.default_rng(seed)
    points = np.column_stack([rng.uniform(0.0, 1500.0, size=(n, 2)),
                              rng.uniform(0.0, depth_span, n)])
    return points, _gaussian_field(points, ranges, seed=seed + 1)


def test_3_anisotropic_field_directional_ranges_are_ordered_major_over_minor():
    """major (X, azimut 90°) ≫ minor (Y, azimut 0°) — istiqamətli
    variogramlar bunu GÖRMƏLİDİR."""
    points, values = _anisotropic_3d_field(n=320, seed=311)
    along = experimental_variogram(points[:, :2], values, n_lags=10, max_lag=500.0,
                                   azimuth_deg=90.0, azimuth_tolerance_deg=20.0)
    across = experimental_variogram(points[:, :2], values, n_lags=10, max_lag=500.0,
                                    azimuth_deg=0.0, azimuth_tolerance_deg=20.0)
    range_along = fit_variogram(along, model=MODEL_EXPONENTIAL).range_
    range_across = fit_variogram(across, model=MODEL_EXPONENTIAL).range_
    assert range_along > 1.5 * range_across


def test_3_anisotropic_kriging_uses_the_directional_continuity():
    """major ≫ minor ≫ şaquli həndəsəsini VERƏN model, izotrop modeldən
    daha yaxşı proqnoz verməlidir."""
    points, values = _anisotropic_3d_field(n=200, seed=321)
    isotropic = OrdinaryKriging(range_=900.0, sill=1.0, nugget=0.0,
                                model=MODEL_EXPONENTIAL, max_neighbors=20)
    anisotropic = OrdinaryKriging(range_=900.0, range_minor=150.0, range_v=20.0,
                                  azimuth_deg=90.0, sill=1.0, nugget=0.0,
                                  model=MODEL_EXPONENTIAL, max_neighbors=20)
    assert _loo_rmse(anisotropic, points, values) < _loo_rmse(isotropic, points, values)


def test_3_detection_recovers_the_major_direction():
    points, values = _anisotropic_3d_field(n=400, seed=331)
    detection = detect_anisotropy(points[:, :2], values, n_directions=12, n_lags=10)
    assert detection.reliable
    delta = abs(detection.azimuth_deg - 90.0)
    assert min(delta, 180.0 - delta) <= 30.0
    assert detection.ratio < 0.8, "minor/major nisbəti 1-dən açıq şəkildə kiçik olmalıdır"


# ══ Test 4 — DÖNMÜŞ anizotropluq ═════════════════════════════════════
def _rotated_field(angle_deg: float, n=250, seed=401):
    """Baş davamlılıq istiqaməti `angle_deg` (şimaldan saat əqrəbi ilə)
    olan sahə: nöqtələr ox-uyğun sahədə qurulub SONRA döndərilir."""
    rng = np.random.default_rng(seed)
    base = rng.uniform(-700.0, 700.0, size=(n, 2))
    values = _gaussian_field(base, (800.0, 90.0), seed=seed + 1)   # major = +X
    # +X (azimut 90°) → `angle_deg` azimutuna aparan dönmə
    return _rotate(base, 90.0 - angle_deg), values


def test_4_rotated_field_is_detected_at_the_right_azimuth():
    for angle in (45.0, 135.0):
        points, values = _rotated_field(angle, n=400, seed=411)
        detection = detect_anisotropy(points, values, n_directions=12, n_lags=8)
        assert detection.reliable
        delta = abs(detection.azimuth_deg - angle)
        assert min(delta, 180.0 - delta) <= 25.0, (
            f"gözlənilən {angle}°, aşkarlanan {detection.azimuth_deg}°")


def test_4_rotation_actually_changes_the_interpolation_geometry():
    """45°-yə dönmüş sahədə DOĞRU azimut (45°) yanlış azimutdan (135°)
    daha yaxşı olmalıdır — dönmə həndəsəyə REAL təsir edir."""
    points, values = _rotated_field(45.0, n=180, seed=421)

    def rmse(azimuth):
        return _loo_rmse(OrdinaryKriging(range_=800.0, range_minor=90.0,
                                         azimuth_deg=azimuth, sill=1.0, nugget=0.0,
                                         model=MODEL_EXPONENTIAL, max_neighbors=20),
                         points, values)

    assert rmse(45.0) < rmse(135.0)


def test_4_rotating_data_and_model_together_leaves_the_result_invariant():
    """Məlumatı və modeli EYNİ bucaqda döndərsək, nəticə DƏYİŞMƏMƏLİDİR
    — transformasiyanın həqiqətən dönmə olduğunun sübutu."""
    rng = np.random.default_rng(431)
    points = rng.uniform(-500.0, 500.0, size=(60, 2))
    values = _gaussian_field(points, (600.0, 80.0), seed=432)
    targets = rng.uniform(-400.0, 400.0, size=(30, 2))

    base = OrdinaryKriging(range_=600.0, range_minor=80.0, azimuth_deg=90.0,
                           sill=1.0, nugget=0.0, model=MODEL_EXPONENTIAL)
    turned = OrdinaryKriging(range_=600.0, range_minor=80.0, azimuth_deg=60.0,
                             sill=1.0, nugget=0.0, model=MODEL_EXPONENTIAL)
    # `_rotate(·, θ)` azimutu `a → a − θ` aparır, ona görə modelin azimutunu
    # 90°-dən 60°-yə endirmək üçün məlumat +30° döndərilir.
    rotated_points = _rotate(points, 30.0)
    rotated_targets = _rotate(targets, 30.0)
    assert np.allclose(base.interpolate(points, values, targets),
                       turned.interpolate(rotated_points, values, rotated_targets),
                       atol=1e-8)


# ══ Test 5 — laylı rezervuar ═════════════════════════════════════════
def _layered_reservoir(n_wells=16, n_layers=10, seed=501):
    """Güclü üfüqi davamlılıq (radius ≈ domen), çox zəif şaquli
    (radius ≈ bir lay qalınlığı): hər layın öz səviyyəsi var."""
    rng = np.random.default_rng(seed)
    wells = rng.uniform(0.0, 1000.0, size=(n_wells, 2))
    layer_depth = 2000.0 + 10.0 * np.arange(n_layers)
    layer_level = np.array([0.25, 0.05, 0.30, 0.08, 0.28, 0.06, 0.32, 0.04, 0.27, 0.07])
    points, values = [], []
    for x, y in wells:
        for k in range(n_layers):
            points.append((x, y, layer_depth[k]))
            values.append(layer_level[k] + 0.00002 * x + rng.normal(0.0, 0.002))
    return np.asarray(points), np.asarray(values), layer_depth, layer_level


def test_5_vertical_variogram_shows_much_shorter_continuity_than_horizontal():
    points, values, _, _ = _layered_reservoir()
    # `vertical_tolerance` OLMADAN "üfüqi" variogram əslində üfüqi OLMUR:
    # 500 m aralıqda ±5° dip pəncərəsi ±44 m şaquli fərqə, yəni bir neçə
    # LAYA icazə verir (bax `experimental_variogram` docstring-i).
    horizontal = experimental_variogram(points, values, dip_deg=0.0,
                                        dip_tolerance_deg=5.0, vertical_tolerance=1e-6,
                                        n_lags=8, max_lag=800.0)
    vertical = vertical_variogram(points, values, horizontal_tolerance=1e-6,
                                  n_lags=8, max_lag=95.0)
    # üfüqi: eyni layda dəyərlər demək olar eynidir → γ kiçik
    # şaquli: qonşu lay tamam başqa səviyyədədir → γ böyük
    assert vertical.gamma[vertical.counts > 0].max() > \
        5.0 * horizontal.gamma[horizontal.counts > 0].max()

    leaky = experimental_variogram(points, values, dip_deg=0.0, dip_tolerance_deg=5.0,
                                   n_lags=8, max_lag=800.0)
    assert leaky.gamma[leaky.counts > 0].max() > \
        10.0 * horizontal.gamma[horizontal.counts > 0].max(), (
        "şaquli bant eni olmadan üfüqi variogram laylararası fərqi udur")


def test_5_vertically_distant_points_do_not_dominate():
    """Şaquli radius bir lay qalınlığı olanda hədəf ÖZ layının
    dəyərini almalıdır — qonşu layın kəskin fərqli dəyəri onu
    "çəkməməlidir"."""
    points, values, layer_depth, layer_level = _layered_reservoir()
    target = np.array([[500.0, 500.0, layer_depth[2]]])

    layered = OrdinaryKriging(range_=2000.0, range_v=6.0, sill=0.02, nugget=0.0,
                              model=MODEL_EXPONENTIAL, max_neighbors=24)
    mixed = OrdinaryKriging(range_=2000.0, range_v=2000.0, sill=0.02, nugget=0.0,
                            model=MODEL_EXPONENTIAL, max_neighbors=24)
    layered_value = layered.interpolate(points, values, target)[0]
    mixed_value = mixed.interpolate(points, values, target)[0]

    assert abs(layered_value - layer_level[2]) < 0.03, (
        f"laylı model öz layının səviyyəsini ({layer_level[2]}) verməlidir, "
        f"aldı {layered_value:.4f}")
    assert abs(mixed_value - layer_level[2]) > abs(layered_value - layer_level[2])


def test_5_vertical_cut_removes_other_layers_entirely():
    """`max_vertical_distance` XAM |ΔZ| ilə kəsir — qonşuluğa yalnız
    öz layının nöqtələri düşür (A2.6)."""
    from imex2d.geology.spatial_search import NeighborhoodConfig, NeighborhoodSelector
    points, values, layer_depth, _ = _layered_reservoir()
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(
        max_neighbors=40, max_vertical_distance=5.0, support_range=2000.0))
    result = selector.select(np.array([[500.0, 500.0, layer_depth[4]]]))
    assert result.count == 16, "yalnız 16 quyunun bu laydakı nöqtəsi qalmalıdır"
    assert np.allclose(points[result.indices, 2], layer_depth[4])


def test_5_layered_kriging_is_more_accurate_than_isotropic_3d_kriging():
    points, values, _, _ = _layered_reservoir(seed=511)
    layered = OrdinaryKriging(range_=2000.0, range_v=6.0, sill=0.02, nugget=0.0,
                              model=MODEL_EXPONENTIAL, max_neighbors=24)
    isotropic = OrdinaryKriging(range_=2000.0, range_v=2000.0, sill=0.02, nugget=0.0,
                                model=MODEL_EXPONENTIAL, max_neighbors=24)
    assert _loo_rmse(layered, points, values) < _loo_rmse(isotropic, points, values)


# ══ boru xəttinin bütövlüyü (A5) ═════════════════════════════════════
def test_full_pipeline_auto_fit_and_auto_anisotropy_runs_end_to_end():
    """XAM data → deneysel variogram → fit → anizotropluq aşkarlanması →
    transformasiya → indeks → qonşuluq → yerli kriging → varians."""
    points, values = _rotated_field(60.0, n=300, seed=601)
    targets = np.random.default_rng(602).uniform(-500.0, 500.0, size=(200, 2))
    kriging = OrdinaryKriging(auto_fit=True, model="auto",
                              auto_detect_anisotropy=True, max_neighbors=16,
                              sectors=4)
    result = kriging.krige(points, values, targets)

    assert result.local is True
    assert result.fit is not None and result.fit.model in ("spherical", "exponential",
                                                           "gaussian")
    assert result.anisotropy.range_major > 0.0
    assert result.anisotropy.horizontal_ratio < 1.0, "anizotropluq aşkarlanmalıdır"
    delta = abs(result.anisotropy.azimuth_deg - 60.0)
    assert min(delta, 180.0 - delta) <= 30.0
    assert np.all(np.isfinite(result.estimate))
    assert np.all(result.variance >= -1e-12)
    assert np.all(result.neighbor_count[np.isfinite(result.estimate)] > 0)
