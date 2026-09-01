"""Deneysel (empirical) variogram, model fitting, istiqamətli anizotropluq.

`interpolation.py`-dəki köhnə `OrdinaryKriging` sferik variogramı SABİT
parametrlərlə (və ya `range = domen/3`, `sill = var(dəyərlər)` evristikası
ilə) işlədirdi — bu modul əsl geostatistik iş axınını əlavə edir:

    məlumat → deneysel variogram (lag-binlənmiş yarım-dəyişkənlik)
            → model fit (sferik/eksponensial/qauss, çəkili ən kiçik kvadrat)
            → (istəyə görə) istiqamətli fit → anizotropluq aşkarlanması
            → `OrdinaryKriging`-ə ötürülən parametrlər

Heç bir funksiya səssizcə "kifayət qədər nöqtə yoxdur" halını gizlətmir —
ya `ValueError` atılır (fit mümkün deyil), ya da nəticədə açıq
`warnings`/`reliable=False` bayrağı olur (bax `AnisotropyDetectionResult`).
Bu, tapşırığın "heç vaxt ədədi problemi gizlətmə" qaydasına uyğundur.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import least_squares

# ── variogram modelləri ("praktiki radius" konvensiyası: hər üç modeldə
#    γ(range_) ≈ 0.95·sill, ona görə range_ modellər arasında müqayisə
#    edilə bilər) ────────────────────────────────────────────────────────
MODEL_SPHERICAL = "spherical"
MODEL_EXPONENTIAL = "exponential"
MODEL_GAUSSIAN = "gaussian"
KNOWN_MODELS: Tuple[str, ...] = (MODEL_SPHERICAL, MODEL_EXPONENTIAL, MODEL_GAUSSIAN)


def spherical(h: np.ndarray, nugget: float, sill: float, range_: float) -> np.ndarray:
    ratio = np.clip(h / max(range_, 1e-9), 0.0, 1.0)
    return nugget + sill * (1.5 * ratio - 0.5 * ratio ** 3)


def exponential(h: np.ndarray, nugget: float, sill: float, range_: float) -> np.ndarray:
    return nugget + sill * (1.0 - np.exp(-3.0 * h / max(range_, 1e-9)))


def gaussian(h: np.ndarray, nugget: float, sill: float, range_: float) -> np.ndarray:
    return nugget + sill * (1.0 - np.exp(-3.0 * (h / max(range_, 1e-9)) ** 2))


MODEL_FUNCS: Dict[str, Callable[[np.ndarray, float, float, float], np.ndarray]] = {
    MODEL_SPHERICAL: spherical,
    MODEL_EXPONENTIAL: exponential,
    MODEL_GAUSSIAN: gaussian,
}


# ── deneysel variogram ──────────────────────────────────────────────────
@dataclass
class ExperimentalVariogram:
    lags: np.ndarray
    gamma: np.ndarray
    counts: np.ndarray
    azimuth_deg: Optional[float] = None
    n_pairs_total: int = 0

    def valid(self) -> "ExperimentalVariogram":
        """Yalnız cüt sayı > 0 olan (dolu) lag-binləri qaytarır."""
        mask = self.counts > 0
        return ExperimentalVariogram(self.lags[mask], self.gamma[mask], self.counts[mask],
                                     self.azimuth_deg, self.n_pairs_total)


def experimental_variogram(points: np.ndarray, values: np.ndarray, n_lags: int = 12,
                           max_lag: Optional[float] = None,
                           azimuth_deg: Optional[float] = None,
                           azimuth_tolerance_deg: float = 22.5,
                           vertical: bool = False,
                           horizontal_tolerance: Optional[float] = None
                           ) -> ExperimentalVariogram:
    """Cüt-cüt lag məsafələrindən binlənmiş yarım-dəyişkənlik hesablayır.

    `azimuth_deg` verilməyəndə (defolt) OMNİDİREKSİONAL — bütün
    istiqamətlər eyni binə düşür. Verilsə, yalnız (X,Y) müstəvisində həmin
    azimutdan (±`azimuth_tolerance_deg`) kənara çıxmayan cütlər saxlanılır
    (Z fərqi nəzərə alınmır — üfüqi istiqamətlilik üçündür).

    `vertical=True` — şaquli variogram: yalnız üfüqi məsafəsi
    `horizontal_tolerance`-dan kiçik olan cütlər (eyni/yaxın quyu, fərqli
    dərinlik) saxlanılır, lag isə |Δz|-dir. Uyğun cüt tapılmazsa `ValueError`
    atılır — sükutla boş nəticə YOX.
    """
    points = np.asarray(points, float)
    values = np.asarray(values, float)
    n = points.shape[0]
    if n < 4:
        raise ValueError(f"Deneysel variogram üçün ən azı 4 nöqtə lazımdır (tapıldı: {n}).")
    if points.ndim != 2 or points.shape[1] not in (2, 3):
        raise ValueError(f"Nöqtələr (n,2) və ya (n,3) olmalıdır, alındı: {points.shape}")

    has_z = points.shape[1] == 3
    xy = points[:, :2]
    z = points[:, 2] if has_z else np.zeros(n)

    iu, ju = np.triu_indices(n, k=1)
    dxy = xy[iu] - xy[ju]
    dz = z[iu] - z[ju]
    dv = values[iu] - values[ju]
    horiz_dist = np.sqrt(np.sum(dxy ** 2, axis=-1))

    if vertical:
        tol = horizontal_tolerance if horizontal_tolerance is not None else 1e-6
        keep = horiz_dist <= tol
        h = np.abs(dz[keep])
        dv = dv[keep]
        if h.size == 0:
            raise ValueError(
                "Şaquli variogram üçün üfüqi məsafəsi tolerantlıq daxilində olan cüt "
                "tapılmadı (eyni X,Y-də fərqli dərinlikli nöqtə yoxdur).")
    else:
        h = horiz_dist
        if azimuth_deg is not None:
            ang = np.degrees(np.arctan2(dxy[:, 0], dxy[:, 1])) % 180.0
            target = azimuth_deg % 180.0
            delta = np.minimum(np.abs(ang - target), 180.0 - np.abs(ang - target))
            keep = delta <= azimuth_tolerance_deg
            h = h[keep]
            dv = dv[keep]
        if has_z:
            dz_keep = dz if azimuth_deg is None else dz[keep]
            h = np.sqrt(h ** 2 + dz_keep ** 2)
        if h.size == 0:
            raise ValueError(
                f"Azimut {azimuth_deg}° (±{azimuth_tolerance_deg}°) daxilində heç bir "
                "cüt tapılmadı — bu istiqamətdə nöqtə sıxlığı kifayət deyil.")

    n_pairs_total = int(h.size)
    if max_lag is None:
        max_lag = float(np.percentile(h, 75)) if h.size else 1.0
    max_lag = max(max_lag, 1e-9)
    lag_width = max_lag / n_lags
    bin_idx = np.minimum((h / lag_width).astype(int), n_lags - 1)

    lags = np.empty(n_lags)
    gamma = np.zeros(n_lags)
    counts = np.zeros(n_lags, dtype=int)
    for b in range(n_lags):
        in_bin = bin_idx == b
        lags[b] = (b + 0.5) * lag_width
        c = int(np.sum(in_bin))
        counts[b] = c
        if c > 0:
            gamma[b] = 0.5 * np.mean(dv[in_bin] ** 2)

    return ExperimentalVariogram(lags, gamma, counts, azimuth_deg, n_pairs_total)


# ── model fitting ────────────────────────────────────────────────────────
@dataclass
class VariogramParameters:
    model: str
    nugget: float
    sill: float
    range_: float
    weighted_rmse: float
    n_pairs: int
    warnings: List[str] = field(default_factory=list)


def _fit_single_model(lags: np.ndarray, gamma: np.ndarray, counts: np.ndarray,
                      model: str, fix_nugget: Optional[float]) -> VariogramParameters:
    func = MODEL_FUNCS[model]
    weights = np.sqrt(counts.astype(float))
    sill0 = float(gamma.max()) if gamma.max() > 0 else 1.0
    # 95%-lik sill həddinə çatan ən kiçik lag — praktiki radiusun ilkin təxmini
    above = np.where(gamma >= 0.95 * sill0)[0]
    range0 = float(lags[above[0]]) if above.size else float(lags[-1])
    range0 = max(range0, float(lags[0]) if lags.size else 1e-6, 1e-6)
    nugget0 = float(np.clip(gamma[0] - (gamma[1] - gamma[0]) if gamma.size > 1 else gamma[0],
                            0.0, sill0)) if fix_nugget is None else float(fix_nugget)

    if fix_nugget is None:
        def resid(p):
            nugget, sill, range_ = p
            return weights * (func(lags, nugget, sill, range_) - gamma)
        x0 = [nugget0, sill0, range0]
        lo = [0.0, 1e-12, 1e-6]
        hi = [sill0 * 4.0 + 1e-9, sill0 * 4.0 + 1e-9, max(lags[-1] * 10.0, range0 * 10.0)]
    else:
        def resid(p):
            sill, range_ = p
            return weights * (func(lags, fix_nugget, sill, range_) - gamma)
        x0 = [sill0, range0]
        lo = [1e-12, 1e-6]
        hi = [sill0 * 4.0 + 1e-9, max(lags[-1] * 10.0, range0 * 10.0)]

    result = least_squares(resid, x0, bounds=(lo, hi))
    if fix_nugget is None:
        nugget, sill, range_ = (float(v) for v in result.x)
    else:
        nugget = float(fix_nugget)
        sill, range_ = (float(v) for v in result.x)

    fitted = func(lags, nugget, sill, range_)
    rmse = float(np.sqrt(np.average((fitted - gamma) ** 2, weights=weights ** 2))) \
        if np.any(weights > 0) else float(np.sqrt(np.mean((fitted - gamma) ** 2)))
    return VariogramParameters(model=model, nugget=nugget, sill=max(sill, 1e-12),
                               range_=max(range_, 1e-6), weighted_rmse=rmse,
                               n_pairs=int(counts.sum()))


def fit_variogram(experimental: ExperimentalVariogram, model: str = "auto",
                  fix_nugget: Optional[float] = None) -> VariogramParameters:
    """Deneysel variograma nugget/sill/range fit edir.

    `model="auto"` — sferik/eksponensial/qauss ÜÇÜNÜN də fit edilib, çəkili
    RMSE-si ən kiçik olanı seçilir (bu, "ən hamar görünən" seçimi DEYİL —
    ədədi uyğunluğa əsaslanır). Konkret model adı verilsə, yalnız o fit
    edilir.
    """
    exp = experimental.valid()
    if exp.lags.size < 3:
        raise ValueError(
            f"Variogram model fit üçün ən azı 3 dolu lag-bin lazımdır "
            f"(tapıldı: {exp.lags.size}). Nöqtə sayı/sıxlığı kifayət deyil.")

    candidates = KNOWN_MODELS if model == "auto" else (model,)
    if model != "auto" and model not in KNOWN_MODELS:
        raise ValueError(f"Naməlum variogram modeli: {model!r}. Dəstəklənən: {KNOWN_MODELS}")

    best: Optional[VariogramParameters] = None
    for m in candidates:
        candidate = _fit_single_model(exp.lags, exp.gamma, exp.counts, m, fix_nugget)
        if best is None or candidate.weighted_rmse < best.weighted_rmse:
            best = candidate
    assert best is not None
    return best


def fit_variogram_from_data(points: np.ndarray, values: np.ndarray, n_lags: int = 12,
                            max_lag: Optional[float] = None, model: str = "auto",
                            fix_nugget: Optional[float] = None) -> VariogramParameters:
    """`experimental_variogram` + `fit_variogram` bir addımda (omnidireksional)."""
    exp = experimental_variogram(points, values, n_lags=n_lags, max_lag=max_lag)
    return fit_variogram(exp, model=model, fix_nugget=fix_nugget)


def select_best_variogram_model(points: np.ndarray, values: np.ndarray,
                                candidates: Tuple[str, ...] = KNOWN_MODELS,
                                cv_method: str = "loo", k: int = 5, seed: int = 42,
                                log_transform: bool = False, n_lags: int = 12,
                                max_lag: Optional[float] = None
                                ) -> Tuple[str, Dict[str, Tuple[VariogramParameters, object]]]:
    """Hər namizəd model üçün fit + REAL cross-validation dəqiqliyi.

    Seçim empirik variogramla ən kiçik SSE-yə görə DEYİL, faktiki Kriging
    ilə leave-one-out/k-fold proqnoz xətasına (RMSE) görə edilir — bu,
    "vizual ən hamar model" seçimindən fərqli olaraq, doğrulama
    performansına əsaslanan seçimdir (bax tapşırıq Phase 20).
    """
    from .cross_validation import k_fold, leave_one_out  # dövri idxalın qarşısı
    from .interpolation import OrdinaryKriging

    exp = experimental_variogram(points, values, n_lags=n_lags, max_lag=max_lag)
    results: Dict[str, Tuple[VariogramParameters, object]] = {}
    for m in candidates:
        try:
            fit = fit_variogram(exp, model=m)
        except ValueError:
            continue
        interpolator = OrdinaryKriging(model=m, nugget=fit.nugget, sill=fit.sill,
                                       range_=fit.range_)
        runner = k_fold if cv_method == "k-fold" else leave_one_out
        kwargs = {"k": k, "seed": seed} if cv_method == "k-fold" else {}
        cv = runner(interpolator, points, values, log_transform=log_transform, **kwargs)
        results[m] = (fit, cv)

    if not results:
        raise ValueError(
            "Heç bir model üçün fit/CV alınmadı — nöqtə sayı model seçimi üçün kifayət deyil.")
    best_model = min(results, key=lambda m: results[m][1].rmse)
    return best_model, results


# ── tam 3D anizotropluq (azimut + major/minor + şaquli) ─────────────────
@dataclass
class AnisotropyParams:
    """Geometrik anizotropluq transformu.

    `azimuth_deg` — major oxun istiqaməti, +Y (Şimal) oxundan SAAT
    ƏQRƏBİ istiqamətində ölçülür (0°=Şimal, 90°=Şərq) — geoloji
    proqramlarda adət olunan konvensiya. `range_major`/`range_minor` —
    üfüqi müstəvidə uyğun radiuslar, `range_vertical` — şaquli radius.
    """
    azimuth_deg: float = 0.0
    range_major: float = 1.0
    range_minor: float = 1.0
    range_vertical: float = 1.0

    def transform(self, points_xyz: np.ndarray) -> np.ndarray:
        """Nöqtələri anizotrop fəzadan izotrop (vahid `range_major`
        radiuslu) fəzaya çevirir — nəticədə adi Evklid məsafəsi birbaşa
        `range_major` ilə istifadə oluna bilər.

        `range_minor == range_major` və `azimuth_deg == 0` olanda
        transform (X,Y)-i sadəcə YERDƏYİŞDİRİR (major=Y, minor=X) —
        Evklid normuna təsir etmir, ona görə əvvəlki "yalnız Z miqyaslama"
        davranışı ilə ƏDƏD-ƏDƏD eynidir (bax `interpolation.py` test
        `test_default_parameters_reproduce_pre_m2_2d_behaviour`).
        """
        pts = np.asarray(points_xyz, float)
        x, y = pts[:, 0], pts[:, 1]
        z = pts[:, 2] if pts.shape[1] > 2 else np.zeros(pts.shape[0])
        theta = np.radians(self.azimuth_deg)
        c, s = np.cos(theta), np.sin(theta)
        major = x * s + y * c
        minor = x * c - y * s
        ref = max(self.range_major, 1e-9)
        return np.column_stack([
            major,
            minor * (ref / max(self.range_minor, 1e-9)),
            z * (ref / max(self.range_vertical, 1e-9)),
        ])


@dataclass
class AnisotropyDetectionResult:
    azimuth_deg: float
    range_major: float
    range_minor: float
    ratio: float
    directional_ranges: Dict[float, float]
    reliable: bool
    warnings: List[str] = field(default_factory=list)


def detect_anisotropy(points: np.ndarray, values: np.ndarray, n_directions: int = 6,
                      n_lags: int = 8, max_lag: Optional[float] = None,
                      azimuth_tolerance_deg: Optional[float] = None,
                      min_pairs_per_direction: int = 10) -> AnisotropyDetectionResult:
    """`n_directions` üfüqi istiqamətdə sferik variogram fit edir, ən böyük
    radiuslu istiqaməti major ox kimi seçir.

    Quyu sayı az olanda (tipik reservoir modelində 5-20 quyu) istiqamətli
    binlərdə kifayət qədər cüt olmaya bilər — bu vəziyyət GİZLƏDİLMİR:
    `reliable=False` və izahlı `warnings` qaytarılır, izotrop defoltla
    (major=minor=omnidireksional radius) geri qayıdılır.
    """
    points = np.asarray(points, float)
    if azimuth_tolerance_deg is None:
        azimuth_tolerance_deg = 90.0 / n_directions

    azimuths = np.linspace(0.0, 180.0, n_directions, endpoint=False)
    directional_ranges: Dict[float, float] = {}
    for az in azimuths:
        try:
            exp = experimental_variogram(points, values, n_lags=n_lags, max_lag=max_lag,
                                         azimuth_deg=float(az),
                                         azimuth_tolerance_deg=azimuth_tolerance_deg)
        except ValueError:
            continue
        if exp.n_pairs_total < min_pairs_per_direction:
            continue
        try:
            fit = fit_variogram(exp, model=MODEL_SPHERICAL)
        except ValueError:
            continue
        directional_ranges[float(az)] = fit.range_

    if len(directional_ranges) < 2:
        warnings = [
            "İstiqamətli variogram üçün kifayət qədər cüt yoxdur (nöqtə sayı/sıxlığı "
            "az) — anizotropluq etibarlı aşkarlanmadı, izotrop defolt işlədilir."]
        try:
            omni = fit_variogram_from_data(points, values, n_lags=n_lags, max_lag=max_lag,
                                           model=MODEL_SPHERICAL)
            range_omni = omni.range_
        except ValueError as exc:
            warnings.append(str(exc))
            span = float(np.sqrt(np.sum((points[:, :2].max(axis=0)
                                         - points[:, :2].min(axis=0)) ** 2)))
            range_omni = max(span / 3.0, 1e-6)
        return AnisotropyDetectionResult(azimuth_deg=0.0, range_major=range_omni,
                                         range_minor=range_omni, ratio=1.0,
                                         directional_ranges=directional_ranges,
                                         reliable=False, warnings=warnings)

    best_az = max(directional_ranges, key=directional_ranges.get)
    range_major = directional_ranges[best_az]
    minor_target = (best_az + 90.0) % 180.0
    closest_minor_az = min(
        directional_ranges,
        key=lambda a: min(abs(a - minor_target), 180.0 - abs(a - minor_target)))
    range_minor = directional_ranges[closest_minor_az]
    ratio = range_minor / range_major if range_major > 0 else 1.0
    return AnisotropyDetectionResult(azimuth_deg=best_az, range_major=range_major,
                                     range_minor=range_minor, ratio=ratio,
                                     directional_ranges=directional_ranges,
                                     reliable=True, warnings=[])
