"""Phase 2/3 — deneysel variogram, model fitting, anizotropluq aşkarlanması.

Sintetik sınaqlar: MƏLUM (ground-truth) parametrlərlə anizotrop/izotrop
Gauss təsadüfi sahəsi qurulur (kovaryans matrisinin Cholesky dekompozisiyası
ilə — bu, `variogram.py`-dəki fit alqoritmindən TAMAMILƏ MÜSTƏQİL yoldur),
sonra fit/aşkarlanma funksiyalarının bu məlum parametrləri (təxmini)
bərpa edib-etmədiyi yoxlanılır. Statistik təbiətli testlər (sonlu təsadüfi
nümunə) olduğu üçün sərt bərabərlik deyil, GENİŞ tolerantlıq işlədilir və
sabit seed istifadə olunur.
"""

from __future__ import annotations

import numpy as np

from imex2d.geology.cross_validation import CrossValidationResult
from imex2d.geology.interpolation import OrdinaryKriging
from imex2d.geology.variogram import (KNOWN_MODELS, MODEL_EXPONENTIAL, MODEL_GAUSSIAN,
                                      MODEL_SPHERICAL, AnisotropyParams, VariogramParameters,
                                      detect_anisotropy, experimental_variogram,
                                      exponential, fit_variogram, fit_variogram_from_data,
                                      gaussian, select_best_variogram_model, spherical)


def _dmat(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    diff = a[:, None, :] - b[None, :, :]
    return np.sqrt(np.sum(diff ** 2, axis=-1))


def _gaussian_random_field(rng, points_xyz, model_func, nugget, sill, range_, mean=100.0):
    """Bilinən variogram modelinə uyğun kovaryanslı Gauss sahəsi.

    `C(h) = sill - γ(h)` (nugget-siz hissə) + diaqonala nugget-dispersiya —
    standart geostatistik "sill-bounded" model üçün kovaryans-variogram
    əlaqəsi. Nəticə fit funksiyalarının test etdiyi koddan asılı DEYİL.
    """
    d = _dmat(points_xyz, points_xyz)
    gamma = model_func(d, 0.0, sill, range_)
    cov = (sill - gamma) + np.eye(points_xyz.shape[0]) * (nugget + 1e-9)
    chol = np.linalg.cholesky(cov)
    return mean + chol @ rng.standard_normal(points_xyz.shape[0])


# ── deneysel variogram ──────────────────────────────────────────────────
def test_experimental_variogram_rejects_too_few_points():
    try:
        experimental_variogram(np.array([[0., 0.], [1., 1.], [2., 2.]]),
                               np.array([1.0, 2.0, 3.0]))
    except ValueError as exc:
        assert "4 nöqtə" in str(exc)
        return
    raise AssertionError("3 nöqtə ilə ValueError gözlənilirdi")


def test_experimental_variogram_gamma_increases_towards_sill():
    """Sferik sahədə deneysel gamma lag artdıqca (range-ə qədər) artmalı,
    böyük lag-larda ~sill ətrafında sabitləşməlidir."""
    rng = np.random.default_rng(1)
    xy = rng.uniform(0, 300, size=(150, 2))
    points = np.column_stack([xy, np.zeros(150)])
    values = _gaussian_random_field(rng, points, spherical, nugget=0.0, sill=4.0, range_=80.0)
    exp = experimental_variogram(xy, values, n_lags=15, max_lag=250.0)
    valid = exp.valid()
    assert valid.lags.size >= 8
    # ilk bin son bindən açıq-aydın kiçik olmalıdır (artan trend, monoton deyil amma)
    assert valid.gamma[0] < valid.gamma[-1]
    # uzaq lag-larda sill-ə yaxın olmalı (± geniş tolerantlıq, sonlu nümunə)
    assert 1.5 < valid.gamma[-3:].mean() < 8.0


# ── model fitting ────────────────────────────────────────────────────────
def test_fit_variogram_recovers_known_spherical_parameters():
    rng = np.random.default_rng(2)
    xy = rng.uniform(0, 300, size=(200, 2))
    points = np.column_stack([xy, np.zeros(200)])
    true_nugget, true_sill, true_range = 0.3, 5.0, 90.0
    values = _gaussian_random_field(rng, points, spherical, true_nugget, true_sill,
                                    true_range)
    fit = fit_variogram_from_data(xy, values, model=MODEL_SPHERICAL)
    assert isinstance(fit, VariogramParameters)
    assert 40.0 < fit.range_ < 160.0, f"range_ bərpa edilmədi: {fit.range_}"
    assert 2.5 < fit.sill < 8.5, f"sill bərpa edilmədi: {fit.sill}"
    assert fit.nugget >= 0.0


def test_fit_variogram_recovers_known_exponential_parameters():
    rng = np.random.default_rng(3)
    xy = rng.uniform(0, 300, size=(200, 2))
    points = np.column_stack([xy, np.zeros(200)])
    values = _gaussian_random_field(rng, points, exponential, 0.0, 4.0, 70.0)
    fit = fit_variogram_from_data(xy, values, model=MODEL_EXPONENTIAL)
    assert 30.0 < fit.range_ < 140.0, f"range_ bərpa edilmədi: {fit.range_}"
    assert 2.0 < fit.sill < 7.0, f"sill bərpa edilmədi: {fit.sill}"


def test_fit_variogram_auto_selects_from_known_models():
    rng = np.random.default_rng(4)
    xy = rng.uniform(0, 300, size=(200, 2))
    points = np.column_stack([xy, np.zeros(200)])
    values = _gaussian_random_field(rng, points, gaussian, 0.0, 3.0, 60.0)
    fit = fit_variogram_from_data(xy, values, model="auto")
    assert fit.model in KNOWN_MODELS
    assert 20.0 < fit.range_ < 120.0


def test_fit_variogram_raises_on_too_few_bins():
    # kvadratın 4 küncü: cütlərin məsafələri yalnız 2 fərqli dəyər alır
    # (tərəf və diaqonal) -> çox az sayda dolu lag-bin
    xy = np.array([[0., 0.], [10., 0.], [0., 10.], [10., 10.]])
    values = np.array([1.0, 2.0, 3.0, 4.0])
    try:
        fit_variogram_from_data(xy, values, n_lags=12)
    except ValueError as exc:
        assert "dolu lag-bin" in str(exc)
        return
    raise AssertionError("Nöqtə/lag azlığında ValueError gözlənilirdi")


def test_fit_variogram_unknown_model_name_raises():
    exp = experimental_variogram(
        np.column_stack([np.linspace(0, 100, 20), np.zeros(20)]),
        np.linspace(1, 2, 20))
    try:
        fit_variogram(exp, model="linear")
    except ValueError as exc:
        assert "Naməlum variogram modeli" in str(exc)
        return
    raise AssertionError("naməlum model adı ilə ValueError gözlənilirdi")


# ── model seçimi (cross-validation əsaslı) ──────────────────────────────
def test_select_best_variogram_model_returns_valid_structure():
    rng = np.random.default_rng(6)
    xy = rng.uniform(0, 300, size=(60, 2))
    points = np.column_stack([xy, np.zeros(60)])
    values = _gaussian_random_field(rng, points, spherical, 0.2, 3.0, 80.0)
    best_model, results = select_best_variogram_model(xy, values)
    assert best_model in results
    assert best_model in KNOWN_MODELS
    for model, (fit, cv) in results.items():
        assert isinstance(fit, VariogramParameters)
        assert isinstance(cv, CrossValidationResult)
        assert np.isfinite(cv.rmse)
    # seçilən model həqiqətən ən kiçik CV RMSE-ə malikdir (seçim məntiqinin
    # özünü qoruyan reqressiya yoxlaması)
    assert results[best_model][1].rmse == min(r.rmse for _, r in results.values())


# ── tam 3D anizotropluq transformu ──────────────────────────────────────
def test_anisotropy_transform_reduces_to_z_only_scaling_when_isotropic_horizontal():
    """azimuth=0, major=minor olanda transform (X,Y)-i sadəcə yerdəyişdirir
    — Evklid məsafəsinə təsir etmir (əvvəlki M2 davranışı qorunur)."""
    aniso = AnisotropyParams(azimuth_deg=0.0, range_major=100.0, range_minor=100.0,
                             range_vertical=20.0)
    points = np.array([[10., 20., 5.], [30., -5., 15.]])
    transformed = aniso.transform(points)
    old_style = points.copy()
    old_style[:, 2] *= 100.0 / 20.0
    assert np.allclose(np.sort(np.abs(transformed[:, :2]), axis=1),
                       np.sort(np.abs(old_style[:, :2]), axis=1))
    assert np.allclose(transformed[:, 2], old_style[:, 2])
    assert np.allclose(_dmat(transformed, transformed), _dmat(old_style, old_style))


def test_anisotropy_transform_stretches_minor_axis_more_than_major():
    """azimuth=0 -> major ox = Y. Eyni fiziki məsafədə (50) Y-boyu (major)
    nöqtə X-boyu (minor) nöqtədən daha 'yaxın' (kiçik transformasiya olunmuş
    məsafə) olmalıdır, çünki minor radius kiçikdir."""
    aniso = AnisotropyParams(azimuth_deg=0.0, range_major=200.0, range_minor=20.0,
                             range_vertical=200.0)
    origin = np.array([[0., 0., 0.]])
    along_major = np.array([[0., 50., 0.]])   # Y-boyu
    along_minor = np.array([[50., 0., 0.]])   # X-boyu
    d_major = _dmat(aniso.transform(origin), aniso.transform(along_major))[0, 0]
    d_minor = _dmat(aniso.transform(origin), aniso.transform(along_minor))[0, 0]
    assert d_major == 50.0
    assert d_minor == 500.0
    assert d_minor > d_major


def test_kriging_with_azimuth_weights_major_axis_neighbour_more():
    """Güclü azimutlu anizotropluqla, hədəfə major-ox istiqamətində olan
    qonşu minor-ox istiqamətindəki (eyni fiziki məsafəli) qonşudan DAHA ÇOX
    təsir etməlidir — nəticə major qonşunun dəyərinə daha yaxın olmalıdır."""
    points = np.array([[0., 50.], [50., 0.], [-200., -200.]])   # 3-cü: qlobal lövbər
    values = np.array([10.0, 20.0, 15.0])
    kriging = OrdinaryKriging(range_=200.0, range_minor=20.0, azimuth_deg=0.0,
                              range_v=1e9, nugget=0.0)
    result = kriging.interpolate(points, values, np.array([[0., 0.]]))[0]
    assert abs(result - 10.0) < abs(result - 20.0), (
        f"major-ox qonşusuna (dəyər 10) daha yaxın olmalı idi, alındı: {result}")


# ── anizotropluq aşkarlanması ────────────────────────────────────────────
def test_detect_anisotropy_recovers_known_azimuth_and_ratio():
    rng = np.random.default_rng(7)
    xy = rng.uniform(0, 400, size=(180, 2))
    true_azimuth, true_major, true_minor = 30.0, 150.0, 40.0
    aniso = AnisotropyParams(azimuth_deg=true_azimuth, range_major=true_major,
                             range_minor=true_minor, range_vertical=true_major)
    points3 = np.column_stack([xy, np.zeros(180)])
    transformed = aniso.transform(points3)
    d = _dmat(transformed, transformed)
    gamma = spherical(d, 0.0, 4.0, true_major)
    cov = (4.0 - gamma) + np.eye(180) * 1e-6
    values = 100.0 + np.linalg.cholesky(cov) @ rng.standard_normal(180)

    result = detect_anisotropy(xy, values, n_directions=6)
    assert result.reliable, f"aşkarlanma etibarsız oldu: {result.warnings}"
    az_error = min(abs(result.azimuth_deg - true_azimuth),
                   180.0 - abs(result.azimuth_deg - true_azimuth))
    assert az_error < 35.0, (
        f"azimut bərpa edilmədi: gözlənilən {true_azimuth}, alınan {result.azimuth_deg}")
    assert result.ratio < 0.85, f"anizotropluq nisbəti aşkarlanmadı: {result.ratio}"


def test_detect_anisotropy_reports_unreliable_with_sparse_wells():
    """Tipik reservoir ssenarisi: 5-6 quyu. İstiqamətli variogram üçün
    kifayət qədər cüt YOXDUR — nəticə UYDURULMAMALI, açıq bildirilməlidir."""
    points = np.array([[0., 0.], [100., 20.], [40., 130.], [180., 90.], [60., 60.]])
    values = np.array([10.0, 12.0, 9.0, 14.0, 11.0])
    result = detect_anisotropy(points, values)
    assert result.reliable is False
    assert result.warnings, "az quyu ilə xəbərdarlıq mətni gözlənilirdi"
    assert result.ratio == 1.0   # izotrop defolta qayıtmalıdır, UYDURULMUR


# ── OrdinaryKriging inteqrasiyası: auto_fit / auto_detect_anisotropy ────
def test_kriging_auto_fit_is_exact_interpolator_at_data_points():
    rng = np.random.default_rng(8)
    xy = rng.uniform(0, 300, size=(30, 2))
    values = rng.uniform(5, 15, size=30)
    kriging = OrdinaryKriging(auto_fit=True, nugget=0.0)
    result = kriging.interpolate(xy, values, xy)
    assert np.allclose(result, values, atol=1e-6)
    assert kriging.last_fit_ is not None
    assert kriging.last_fit_.range_ > 0


def test_kriging_auto_fit_falls_back_with_warning_when_too_few_points():
    points = np.array([[0., 0.], [10., 0.], [0., 10.], [10., 10.]])
    values = np.array([1.0, 2.0, 3.0, 4.0])
    kriging = OrdinaryKriging(auto_fit=True, nugget=0.0)
    result = kriging.interpolate(points, values, np.array([[5., 5.]]))
    assert np.isfinite(result[0])
    assert kriging.last_fit_ is None
    assert any("domen/3" in w for w in kriging.last_warnings_)


def test_kriging_default_still_uses_domain_third_heuristic_without_auto_fit():
    """`auto_fit=False` (defolt) olanda DAVRANIŞ dəyişməməlidir — heç bir
    fit cəhdi edilmir, `last_fit_` həmişə `None` qalır."""
    rng = np.random.default_rng(9)
    xy = rng.uniform(0, 300, size=(30, 2))
    values = rng.uniform(5, 15, size=30)
    kriging = OrdinaryKriging(nugget=0.05)
    kriging.interpolate(xy, values, xy[:1])
    assert kriging.last_fit_ is None
    assert kriging.last_anisotropy_ is None


def test_model_auto_without_auto_fit_raises():
    try:
        OrdinaryKriging(model="auto", auto_fit=False)
    except ValueError as exc:
        assert "auto_fit=True" in str(exc)
        return
    raise AssertionError("model='auto' + auto_fit=False ValueError verməli idi")


def test_unknown_model_name_raises_at_construction():
    try:
        OrdinaryKriging(model="linear")
    except ValueError as exc:
        assert "Naməlum variogram modeli" in str(exc)
        return
    raise AssertionError("naməlum model adı ilə ValueError gözlənilirdi")
