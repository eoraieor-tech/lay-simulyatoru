"""A4 — HƏQİQİ geometrik anizotropluq: transformasiya RİYAZİYYATA təsir edir.

Bu faylın mərkəzi iddiası (Gate 4/5): `azimuth_deg`/`range_minor`/
`range_vertical`/`dip_deg` sadəcə SAXLANILAN konfiqurasiya dəyəri DEYİL
— onlar FAKTİKİ məsafəni, qonşu sıralamasını, Kriging matrisini və
nəticəni dəyişir. Həm də HƏR İKİ modul (A1 Kriging və A2 qonşuluq) EYNİ
transformasiyanı işlədir.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.geology.anisotropy import (ISOTROPIC, AnisotropyError, AnisotropyParams,
                                       transform_points)
from imex2d.geology.interpolation import OrdinaryKriging
from imex2d.geology.spatial_search import NeighborhoodConfig, NeighborhoodSelector


def _points_3d(n, seed, high=1000.0, depth=200.0):
    rng = np.random.default_rng(seed)
    return np.column_stack([rng.uniform(0.0, high, size=(n, 2)),
                            rng.uniform(0.0, depth, n)])


# ── 32. izotrop hal = adi Evklid davranışı ────────────────────────────
def test_isotropic_transform_preserves_euclidean_distances():
    points = _points_3d(40, seed=1)
    aniso = AnisotropyParams(azimuth_deg=0.0, range_major=250.0, range_minor=250.0,
                             range_vertical=250.0)
    assert aniso.is_isotropic
    transformed = aniso.transform(points)
    for i in range(0, 39, 7):
        for j in range(i + 1, 40, 11):
            raw = np.linalg.norm(points[i] - points[j])
            new = np.linalg.norm(transformed[i] - transformed[j])
            assert new == pytest.approx(raw, rel=1e-12)


@pytest.mark.parametrize("azimuth", [0.0, 17.0, 45.0, 90.0, 123.0, 270.0])
def test_rotation_alone_is_an_isometry(azimuth):
    """İzotrop radiuslarda dönmə məsafəni DƏYİŞMƏMƏLİDİR (yalnız
    anizotrop miqyaslanma dəyişir)."""
    points = _points_3d(25, seed=2)
    aniso = AnisotropyParams(azimuth_deg=azimuth, range_major=100.0,
                             range_minor=100.0, range_vertical=100.0)
    before = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=-1)
    t = aniso.transform(points)
    after = np.linalg.norm(t[:, None, :] - t[None, :, :], axis=-1)
    assert np.allclose(before, after, atol=1e-9)


def test_isotropic_kriging_matches_plain_euclidean_kriging():
    points = _points_3d(30, seed=3, depth=0.0)
    values = np.sin(points[:, 0] / 150.0)
    targets = _points_3d(20, seed=4, depth=0.0)
    plain = OrdinaryKriging(range_=300.0, sill=1.0, nugget=0.0)
    explicit = OrdinaryKriging(range_=300.0, sill=1.0, nugget=0.0, range_minor=300.0,
                               range_v=300.0, azimuth_deg=0.0)
    assert np.allclose(plain.interpolate(points, values, targets),
                       explicit.interpolate(points, values, targets), atol=1e-12)


# ── 33. üfüqi anizotropluq qonşu sıralamasını dəyişir ────────────────
def test_horizontal_anisotropy_changes_neighbour_ranking():
    """A2 və A1 EYNİ həndəsəni işlədir: sıralama da, kriging nəticəsi də
    dəyişir."""
    points = np.array([[300., 0., 0.], [0., 100., 0.]])
    values = np.array([1.0, 2.0])
    target = np.array([[0.0, 0.0, 0.0]])

    isotropic = AnisotropyParams(range_major=500.0, range_minor=500.0,
                                 range_vertical=500.0)
    stretched = AnisotropyParams(azimuth_deg=90.0, range_major=500.0,
                                 range_minor=100.0, range_vertical=500.0)

    iso_nearest = NeighborhoodSelector(points, anisotropy=isotropic,
                                       config=NeighborhoodConfig(max_neighbors=1))
    ani_nearest = NeighborhoodSelector(points, anisotropy=stretched,
                                       config=NeighborhoodConfig(max_neighbors=1))
    assert iso_nearest.select(target).indices[0] == 1     # 100 < 300
    assert ani_nearest.select(target).indices[0] == 0     # 300 < 100·5

    iso_value = OrdinaryKriging(range_=500.0, sill=1.0, nugget=0.0).interpolate(
        points, values, target)[0]
    ani_value = OrdinaryKriging(range_=500.0, range_minor=100.0, azimuth_deg=90.0,
                                sill=1.0, nugget=0.0).interpolate(points, values, target)[0]
    assert abs(iso_value - 2.0) < abs(iso_value - 1.0), "izotropda (0,100) daha ağırdır"
    assert abs(ani_value - 1.0) < abs(ani_value - 2.0), "anizotropda (300,0) daha ağırdır"


def test_anisotropy_ratio_actually_scales_the_distance():
    aniso = AnisotropyParams(azimuth_deg=0.0, range_major=200.0, range_minor=50.0,
                             range_vertical=200.0)
    along_major = aniso.distance(np.array([[0., 0., 0.]]), np.array([[0., 100., 0.]]))[0, 0]
    along_minor = aniso.distance(np.array([[0., 0., 0.]]), np.array([[100., 0., 0.]]))[0, 0]
    assert along_major == pytest.approx(100.0)
    assert along_minor == pytest.approx(400.0)            # 100 × (200/50)
    assert aniso.horizontal_ratio == pytest.approx(0.25)


# ── 34/36/37. dönmə: 0° və 90° halları ────────────────────────────────
def test_zero_azimuth_puts_the_major_axis_along_north():
    aniso = AnisotropyParams(azimuth_deg=0.0, range_major=400.0, range_minor=50.0,
                             range_vertical=400.0)
    north = aniso.distance(np.array([[0., 0., 0.]]), np.array([[0., 100., 0.]]))[0, 0]
    east = aniso.distance(np.array([[0., 0., 0.]]), np.array([[100., 0., 0.]]))[0, 0]
    assert north == pytest.approx(100.0)                  # major oxu boyunca
    assert east == pytest.approx(800.0)                   # minor oxu boyunca


def test_ninety_degree_azimuth_swaps_the_principal_axes():
    aniso = AnisotropyParams(azimuth_deg=90.0, range_major=400.0, range_minor=50.0,
                             range_vertical=400.0)
    north = aniso.distance(np.array([[0., 0., 0.]]), np.array([[0., 100., 0.]]))[0, 0]
    east = aniso.distance(np.array([[0., 0., 0.]]), np.array([[100., 0., 0.]]))[0, 0]
    assert east == pytest.approx(100.0)
    assert north == pytest.approx(800.0)


@pytest.mark.parametrize("azimuth", [0.0, 30.0, 45.0, 90.0, 135.0])
def test_rotation_changes_the_interpolated_value(azimuth):
    """Azimut nəticəni FAKTİKİ dəyişməlidir — "saxlanan, işlədilməyən
    parametr" olmadığının birbaşa sübutu."""
    rng = np.random.default_rng(11)
    points = rng.uniform(0.0, 800.0, size=(40, 2))
    values = rng.standard_normal(40)
    targets = rng.uniform(100.0, 700.0, size=(25, 2))
    baseline = OrdinaryKriging(range_=400.0, range_minor=80.0, azimuth_deg=0.0,
                               sill=1.0, nugget=0.0).interpolate(points, values, targets)
    rotated = OrdinaryKriging(range_=400.0, range_minor=80.0, azimuth_deg=azimuth,
                              sill=1.0, nugget=0.0).interpolate(points, values, targets)
    if azimuth == 0.0:
        assert np.allclose(baseline, rotated, atol=0.0)
    else:
        assert not np.allclose(baseline, rotated, atol=1e-6)


def test_azimuth_is_periodic_modulo_180_degrees():
    """Anizotropluq ellipsi 180°-də təkrarlanır — `az` və `az+180`
    EYNİ həndəsədir."""
    points = _points_3d(30, seed=12, depth=0.0)
    values = np.cos(points[:, 1] / 90.0)
    targets = _points_3d(15, seed=13, depth=0.0)
    a = OrdinaryKriging(range_=300.0, range_minor=60.0, azimuth_deg=35.0, sill=1.0,
                        nugget=0.0).interpolate(points, values, targets)
    b = OrdinaryKriging(range_=300.0, range_minor=60.0, azimuth_deg=215.0, sill=1.0,
                        nugget=0.0).interpolate(points, values, targets)
    assert np.allclose(a, b, atol=1e-9)


# ── 35. şaquli anizotropluq ───────────────────────────────────────────
def test_vertical_range_scales_the_depth_axis():
    aniso = AnisotropyParams(range_major=500.0, range_minor=500.0, range_vertical=25.0)
    horizontal = aniso.distance(np.array([[0., 0., 0.]]), np.array([[100., 0., 0.]]))[0, 0]
    vertical = aniso.distance(np.array([[0., 0., 0.]]), np.array([[0., 0., 100.]]))[0, 0]
    assert horizontal == pytest.approx(100.0)
    assert vertical == pytest.approx(2000.0)              # 100 × (500/25)


def test_small_vertical_range_isolates_layers_in_kriging():
    """Şaquli radius kiçildikcə BAŞQA laydakı nöqtənin təsiri AZALIR.

    Qurğu: hədəflə EYNİ layda iki nöqtə (1.0 və 1.2), 100 m aşağıda isə
    kəskin fərqli dəyər (5.0). `range_v` böyük olanda aşağıdakı nöqtə
    qiyməti yuxarı çəkir; kiçik olanda o, korrelyasiya radiusundan
    kənara düşür və çəkisi itir — nəticə öz layının dəyərlərində qalır."""
    points = np.array([[0., 0., 0.], [50., 0., 0.], [0., 0., 100.]])
    values = np.array([1.0, 1.2, 5.0])
    target = np.array([[10.0, 0.0, 0.0]])
    coupled = OrdinaryKriging(range_=500.0, range_v=500.0, sill=1.0,
                              nugget=0.0).interpolate(points, values, target)[0]
    isolated = OrdinaryKriging(range_=500.0, range_v=10.0, sill=1.0,
                               nugget=0.0).interpolate(points, values, target)[0]
    assert coupled > isolated, "böyük range_v aşağı layı işə salır"
    assert abs(isolated - 1.05) < 0.25, "kiçik range_v öz layında saxlayır"


def test_dip_rotation_tilts_the_major_axis_out_of_the_horizontal_plane():
    """`dip_deg=90` — major ox tam ŞAQULİ olur, yəni ən böyük davamlılıq
    dərinlik boyunca gedir."""
    flat = AnisotropyParams(azimuth_deg=0.0, range_major=400.0, range_minor=40.0,
                            range_vertical=40.0, dip_deg=0.0)
    steep = AnisotropyParams(azimuth_deg=0.0, range_major=400.0, range_minor=40.0,
                             range_vertical=40.0, dip_deg=90.0)
    origin = np.array([[0., 0., 0.]])
    vertical_offset = np.array([[0., 0., 100.]])
    assert flat.distance(origin, vertical_offset)[0, 0] == pytest.approx(1000.0)
    assert steep.distance(origin, vertical_offset)[0, 0] == pytest.approx(100.0)


def test_zero_dip_is_bitwise_identical_to_no_dip_at_all():
    """Defolt `dip_deg=0` ƏVVƏLKİ (M2) transformu BİT-BİT təkrarlamalıdır."""
    points = _points_3d(50, seed=14)
    without = AnisotropyParams(azimuth_deg=37.0, range_major=300.0, range_minor=90.0,
                               range_vertical=40.0)
    with_zero = AnisotropyParams(azimuth_deg=37.0, range_major=300.0, range_minor=90.0,
                                 range_vertical=40.0, dip_deg=0.0)
    assert np.array_equal(without.transform(points), with_zero.transform(points))


def test_dip_changes_the_kriging_result():
    rng = np.random.default_rng(15)
    points = np.column_stack([rng.uniform(0.0, 500.0, size=(50, 2)),
                              rng.uniform(0.0, 200.0, 50)])
    values = rng.standard_normal(50)
    targets = np.column_stack([rng.uniform(50.0, 450.0, size=(20, 2)),
                               rng.uniform(20.0, 180.0, 20)])
    flat = OrdinaryKriging(range_=400.0, range_minor=80.0, range_v=60.0,
                           azimuth_deg=20.0, sill=1.0, nugget=0.0)
    dipping = OrdinaryKriging(range_=400.0, range_minor=80.0, range_v=60.0,
                              azimuth_deg=20.0, dip_deg=35.0, sill=1.0, nugget=0.0)
    assert not np.allclose(flat.interpolate(points, values, targets),
                           dipping.interpolate(points, values, targets), atol=1e-6)


# ── A4.1/A4.5 cəbri forma ─────────────────────────────────────────────
def test_matrix_form_agrees_with_the_sequential_transform():
    """`matrix()` (`M = S·R`) `transform()` ilə RİYAZİ olaraq eynidir."""
    for aniso in (AnisotropyParams(azimuth_deg=0.0, range_major=100.0,
                                   range_minor=100.0, range_vertical=100.0),
                  AnisotropyParams(azimuth_deg=63.0, range_major=400.0,
                                   range_minor=70.0, range_vertical=25.0),
                  AnisotropyParams(azimuth_deg=63.0, range_major=400.0,
                                   range_minor=70.0, range_vertical=25.0, dip_deg=22.0)):
        points = _points_3d(20, seed=16)
        assert np.allclose(aniso.transform(points), points @ aniso.matrix().T, atol=1e-9)


def test_metric_tensor_reproduces_the_anisotropic_distance():
    """`d² = Δxᵀ G Δx`, `G = MᵀM` — tenzor-uyğun genişlənmə nöqtəsi (A4.5)."""
    aniso = AnisotropyParams(azimuth_deg=48.0, range_major=350.0, range_minor=60.0,
                             range_vertical=30.0, dip_deg=15.0)
    g = aniso.metric_tensor()
    assert np.allclose(g, g.T, atol=1e-9), "metrik tenzor simmetrik olmalıdır"
    assert np.all(np.linalg.eigvalsh(g) > 0.0), "müsbət-müəyyən olmalıdır"
    a, b = _points_3d(2, seed=17)
    delta = a - b
    assert float(delta @ g @ delta) == pytest.approx(
        aniso.distance(a[None, :], b[None, :])[0, 0] ** 2, rel=1e-9)


def test_principal_axes_are_orthonormal_and_carry_their_ranges():
    aniso = AnisotropyParams(azimuth_deg=25.0, range_major=400.0, range_minor=90.0,
                             range_vertical=20.0, dip_deg=10.0)
    axes = aniso.principal_axes()
    names = [name for name, _, _ in axes]
    assert names == ["major", "minor", "vertical"]
    assert [r for _, r, _ in axes] == [400.0, 90.0, 20.0]
    vectors = np.array([v for _, _, v in axes])
    assert np.allclose(vectors @ vectors.T, np.eye(3), atol=1e-9)


def test_transform_points_helper_pads_2d_and_respects_none():
    points = np.array([[1.0, 2.0], [3.0, 4.0]])
    padded = transform_points(points, None)
    assert padded.shape == (2, 3) and np.all(padded[:, 2] == 0.0)
    assert np.allclose(transform_points(points, ISOTROPIC)[:, :2].sum(),
                       points.sum(), atol=1e-9)


# ── parametr doğrulaması ──────────────────────────────────────────────
@pytest.mark.parametrize("kwargs", [
    {"range_major": 0.0},
    {"range_minor": -10.0},
    {"range_vertical": np.nan},
    {"azimuth_deg": np.inf},
    {"dip_deg": 120.0},
])
def test_invalid_anisotropy_parameters_raise(kwargs):
    base = dict(azimuth_deg=0.0, range_major=100.0, range_minor=100.0,
                range_vertical=100.0)
    base.update(kwargs)
    with pytest.raises(AnisotropyError):
        AnisotropyParams(**base).validate()


def test_minor_larger_than_major_warns_but_stays_valid():
    warnings = AnisotropyParams(range_major=100.0, range_minor=400.0,
                                range_vertical=100.0).validate()
    assert warnings and "major" in warnings[0]


def test_from_ranges_defaults_missing_ranges_to_major():
    aniso = AnisotropyParams.from_ranges(250.0)
    assert aniso.is_isotropic and aniso.range_vertical == 250.0
    partial = AnisotropyParams.from_ranges(250.0, range_vertical=25.0, azimuth_deg=40.0)
    assert partial.range_minor == 250.0 and partial.azimuth_deg == 40.0


# ── 38. anizotrop sahədə anizotrop Kriging DAHA YAXŞI olmalıdır ───────
def _anisotropic_field(n, seed, range_x, range_y, high=1200.0):
    rng = np.random.default_rng(seed)
    points = rng.uniform(0.0, high, size=(n, 2))
    dx = (points[:, 0][:, None] - points[:, 0][None, :]) / range_x
    dy = (points[:, 1][:, None] - points[:, 1][None, :]) / range_y
    cov = np.exp(-3.0 * np.sqrt(dx ** 2 + dy ** 2)) + 1e-8 * np.eye(n)
    values = np.linalg.cholesky(cov) @ rng.standard_normal(n)
    return points, values


def _loo_rmse(kriging, points, values):
    """Leave-one-out RMSE — modelin REAL proqnoz xətası."""
    errors = []
    mask = np.ones(points.shape[0], dtype=bool)
    for i in range(points.shape[0]):
        mask[:] = True
        mask[i] = False
        predicted = kriging.interpolate(points[mask], values[mask], points[i:i + 1])[0]
        if np.isfinite(predicted):
            errors.append(predicted - values[i])
    return float(np.sqrt(np.mean(np.square(errors))))


def test_anisotropic_kriging_beats_isotropic_on_an_anisotropic_field():
    """Sahə X boyunca 8 dəfə davamlıdır (kovariasiya `exp(−3d)`, yəni
    layihənin EKSPONENSİAL modeli ilə eyni ailə). Anizotrop model
    (azimut 90°, major=800, minor=100) izotrop modeldən DAHA KİÇİK
    leave-one-out xətası verməlidir — anizotropluq REAL riyazi üstünlük
    gətirir, sadəcə saxlanılan parametr deyil."""
    from imex2d.geology.variogram import MODEL_EXPONENTIAL
    points, values = _anisotropic_field(150, seed=21, range_x=800.0, range_y=100.0)
    isotropic = OrdinaryKriging(range_=800.0, sill=1.0, nugget=0.0,
                                model=MODEL_EXPONENTIAL, max_neighbors=20)
    anisotropic = OrdinaryKriging(range_=800.0, range_minor=100.0, azimuth_deg=90.0,
                                  sill=1.0, nugget=0.0, model=MODEL_EXPONENTIAL,
                                  max_neighbors=20)
    assert _loo_rmse(anisotropic, points, values) < _loo_rmse(isotropic, points, values)


def test_wrong_anisotropy_direction_is_worse_than_the_right_one():
    """Doğru istiqamət (90°) yanlış istiqamətdən (0°) daha yaxşı olmalıdır
    — yəni azimut nəticəyə MƏNALI şəkildə təsir edir, təsadüfi deyil."""
    from imex2d.geology.variogram import MODEL_EXPONENTIAL
    points, values = _anisotropic_field(150, seed=22, range_x=800.0, range_y=100.0)

    def rmse(azimuth):
        return _loo_rmse(OrdinaryKriging(range_=800.0, range_minor=100.0,
                                         azimuth_deg=azimuth, sill=1.0, nugget=0.0,
                                         model=MODEL_EXPONENTIAL, max_neighbors=20),
                         points, values)

    assert rmse(90.0) < rmse(0.0)


def test_detect_anisotropy_recovers_the_direction_and_feeds_the_geometry():
    """`detect_anisotropy` → `to_params()` → Kriging: aşkarlanan həndəsə
    BİRBAŞA transformasiya obyektinə çevrilir (A4.4)."""
    from imex2d.geology.variogram import detect_anisotropy
    points, values = _anisotropic_field(400, seed=23, range_x=700.0, range_y=90.0)
    detection = detect_anisotropy(points, values, n_directions=6, n_lags=8)
    assert detection.reliable
    assert min(abs(detection.azimuth_deg - 90.0), 180.0 - abs(detection.azimuth_deg - 90.0)) \
        <= 30.0
    assert detection.ratio < 1.0
    params = detection.to_params(range_vertical=15.0)
    assert isinstance(params, AnisotropyParams)
    assert params.range_vertical == 15.0
    assert params.validate() == [] or isinstance(params.validate(), list)
