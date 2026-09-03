"""Deneysel (empirical) variogram, model fitting, istiqamətli anizotropluq (A3).

`interpolation.py`-dəki köhnə `OrdinaryKriging` sferik variogramı SABİT
parametrlərlə (və ya `range = domen/3`, `sill = var(dəyərlər)` evristikası
ilə) işlədirdi — bu modul əsl geostatistik iş axınını əlavə edir:

    məlumat → deneysel variogram (lag-binlənmiş yarım-dəyişkənlik)
            → model fit (sferik/eksponensial/qauss, çəkili ən kiçik kvadrat)
            → model müqayisəsi/seçimi (ədədi meyarla)
            → (istəyə görə) istiqamətli fit → anizotropluq aşkarlanması
            → PARAMETR DOĞRULAMASI (etibarsız model Kriging-ə ÇATMIR)
            → `OrdinaryKriging`-ə ötürülən parametrlər

Heç bir funksiya səssizcə "kifayət qədər nöqtə yoxdur" halını gizlətmir —
ya `ValueError` atılır (fit mümkün deyil), ya da nəticədə açıq
`warnings`/`reliable=False` bayrağı olur (bax `AnisotropyDetectionResult`).
Bu, tapşırığın "heç vaxt ədədi problemi gizlətmə" qaydasına uyğundur.

KONVENSİYA (A1.6/A3.4) — bu modul YALNIZ YARIM-DƏYİŞKƏNLİK (semivariogram)
`γ(h)` istehsal edir::

    γ(h) = nugget + c·g(h/a),      g(0)=0, g(∞)=1
    C(h) = (nugget + c) − γ(h)     (kovariasiya, `covariance()`)
    C(0) = nugget + c = "total sill"

`sill` PARAMETRİ QURULUŞLU (structured) hissədir — `c`, yəni nugget
DAXİL DEYİL. `MODEL_FUNCS[m](h, nugget, sill, range_)` düsturuna baxın:
`nugget + sill·g(...)`. Bu, layihənin ƏVVƏLKİ konvensiyasıdır və
`OrdinaryKriging` də məhz bu funksiyaları çağırır, ona görə Kriging
matrisi ilə variogram arasında konvensiya UYĞUNSUZLUĞU MÜMKÜN DEYİL.

RADİUS KONVENSİYASI (A3.4) — hər üç modeldə `range_` PRAKTİKİ radiusdur:
`γ(range_) ≈ 0.95·(nugget+sill)` (sferikdə DƏQİQ `nugget+sill`,
eksponensial/qaussda `−3` əmsalı ilə 95%). Beləliklə `range_` modellər
arasında birbaşa müqayisə edilə bilər — eksponensialın "effektiv
radiusu" (`a_e = range_/3`) ilə QARIŞDIRILMAMALIDIR.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import least_squares

# `AnisotropyParams` HƏNDƏSƏSİ `geology/anisotropy.py`-yə köçüb (A4 — vahid
# transformasiya yolu); burada GERİYƏ-UYĞUNLUQ üçün yenidən eksport olunur,
# çünki `sgs.py`/`facies.py`/`spatial_search.py` və mövcud testlər onu bu
# moduldan idxal edir.
from .anisotropy import AnisotropyParams  # noqa: F401

# ── variogram modelləri ("praktiki radius" konvensiyası: hər üç modeldə
#    γ(range_) ≈ 0.95·sill, ona görə range_ modellər arasında müqayisə
#    edilə bilər) ────────────────────────────────────────────────────────
MODEL_SPHERICAL = "spherical"
MODEL_EXPONENTIAL = "exponential"
MODEL_GAUSSIAN = "gaussian"
KNOWN_MODELS: Tuple[str, ...] = (MODEL_SPHERICAL, MODEL_EXPONENTIAL, MODEL_GAUSSIAN)

#: `range_` (praktiki radius) → modelin öz "effektiv" (e-qatlama) radiusu.
#: Yalnız SƏNƏDLƏŞDİRMƏ/hesabat üçündür — hesablamada `range_` işlədilir.
EFFECTIVE_RANGE_FACTOR: Dict[str, float] = {
    MODEL_SPHERICAL: 1.0,       # sferikdə praktiki = həqiqi radius
    MODEL_EXPONENTIAL: 1.0 / 3.0,
    MODEL_GAUSSIAN: 1.0 / np.sqrt(3.0),
}


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


def semivariance(h: np.ndarray, nugget: float, sill: float, range_: float,
                 model: str = MODEL_SPHERICAL) -> np.ndarray:
    """`γ(h)` — modul docstring-indəki konvensiya ilə (A1.6 vahid giriş)."""
    if model not in MODEL_FUNCS:
        raise ValueError(f"Naməlum variogram modeli: {model!r}. Dəstəklənən: {KNOWN_MODELS}")
    return MODEL_FUNCS[model](np.asarray(h, float), nugget, sill, range_)


def covariance(h: np.ndarray, nugget: float, sill: float, range_: float,
               model: str = MODEL_SPHERICAL) -> np.ndarray:
    """`C(h) = C(0) − γ(h)`, `C(0) = nugget + sill`.

    NUGGET-in başlanğıcdakı SIÇRAYIŞI burada DÜZGÜN işlənir: ciddi
    mənada `γ(0) = 0` (nöqtənin özü ilə fərqi sıfırdır), `γ(0⁺) =
    nugget`. `MODEL_FUNCS[...]` isə `h=0`-da `nugget` qaytarır (davamlı
    ifadə), ona görə burada başlanğıc AÇIQ şəkildə sıfırlanır — nəticədə::

        C(0)  = nugget + sill        (tam dispersiya)
        C(0⁺) = sill                 (nugget qədər sıçrayış)
        C(∞)  = 0

    Bu, `interpolation.OrdinaryKriging._variogram(zero_at_origin=True)`-
    in məlumat-məlumat matrisində etdiyi ilə EYNİ qaydadır — iki modul
    arasında işarə/konvensiya uyğunsuzluğu ola bilməz (A1.6).

    Adi Kriging sistemi `γ` ilə qurulur (bax `interpolation.py`) — bu
    funksiya yalnız kovariasiya dilində yoxlama/hesabat üçündür."""
    h = np.asarray(h, float)
    total_sill = float(nugget) + float(sill)
    gamma = np.where(h <= 1e-12, 0.0,
                     semivariance(h, nugget, sill, range_, model))
    return total_sill - gamma


def effective_range(range_: float, model: str) -> float:
    """Praktiki radiusdan modelin e-qatlama ("effektiv") radiusu."""
    if model not in EFFECTIVE_RANGE_FACTOR:
        raise ValueError(f"Naməlum variogram modeli: {model!r}. Dəstəklənən: {KNOWN_MODELS}")
    return float(range_) * EFFECTIVE_RANGE_FACTOR[model]


# ── deneysel variogram ──────────────────────────────────────────────────
@dataclass
class ExperimentalVariogram:
    """Binlənmiş deneysel yarım-dəyişkənlik + binin necə qurulduğu.

    `lags` bin MƏRKƏZLƏRİ, `gamma` bin ortası, `counts` bindəki CÜT sayı.
    `lag_width`/`lag_tolerance` binin həndəsəsini (A3.1) saxlayır ki,
    hesabatda "hansı konfiqurasiya ilə" sualı cavabsız qalmasın.
    """
    lags: np.ndarray
    gamma: np.ndarray
    counts: np.ndarray
    azimuth_deg: Optional[float] = None
    n_pairs_total: int = 0
    lag_width: float = 0.0
    lag_tolerance: Optional[float] = None
    dip_deg: Optional[float] = None
    vertical: bool = False
    max_lag: float = 0.0
    warnings: List[str] = field(default_factory=list)

    def valid(self, min_pairs: int = 1) -> "ExperimentalVariogram":
        """Yalnız cüt sayı `min_pairs`-dan AZ OLMAYAN binləri qaytarır.

        `min_pairs` (A3.1) — 1-2 cütdən hesablanan `γ` statistik olaraq
        MƏNASIZDIR (dispersiyası sillin özü qədərdir); praktikada 30
        tövsiyə olunur, defolt 1 isə ƏVVƏLKİ davranışı (counts > 0)
        birəbir saxlayır."""
        threshold = max(int(min_pairs), 1)
        mask = self.counts >= threshold
        return ExperimentalVariogram(
            self.lags[mask], self.gamma[mask], self.counts[mask], self.azimuth_deg,
            self.n_pairs_total, self.lag_width, self.lag_tolerance, self.dip_deg,
            self.vertical, self.max_lag, list(self.warnings))

    @property
    def n_bins(self) -> int:
        return int(self.lags.size)


#: Deneysel variogramda işlədilən MAKSİMUM cüt sayı. Bundan çox cüt
#: olanda DETERMİNİSTİK təsadüfi alt-nümunə götürülür (bax
#: `_pair_arrays`). 2·10⁶ hədd `n ≈ 2000` nöqtəyə uyğundur — bütün mövcud
#: testlər (ən böyüyü 400 nöqtə = 80 min cüt) bu həddin ALTINDADIR, ona
#: görə onların nəticəsi DƏQİQ və DƏYİŞMƏZ qalır.
MAX_VARIOGRAM_PAIRS = 2_000_000

#: Alt-nümunənin seed-i — SABİT, ona görə nəticə təkrarlana biləndir.
_PAIR_SUBSAMPLE_SEED = 20240917


def _pair_arrays(points: np.ndarray, values: np.ndarray,
                 max_pairs: Optional[int] = MAX_VARIOGRAM_PAIRS):
    """Yuxarı üçbucaq cüt indeksləri üçün `(dxy, dz, dv, horiz_dist, has_z)`.

    `max_pairs` — cüt sayı bu həddi keçəndə DETERMİNİSTİK (sabit seed-li)
    təsadüfi alt-nümunə götürülür. Səbəb: cüt sayı `n(n−1)/2`, yəni 5000
    nöqtə üçün 12.5 milyon cüt — ölçülüb ki, bu, 858 MB pik yaddaş və
    1.6 s vaxt deməkdir, halbuki variogram üçün bir neçə yüz min cüt
    onsuz da statistik olaraq kifayətdir (hər lag-bində minlərlə cüt).

    Alt-nümunə BƏRABƏR ehtimallıdır, ona görə lag paylanmasını TƏHRİF
    ETMİR; `max_pairs=None` tam (dəqiq) hesablamanı bərpa edir."""
    n = points.shape[0]
    has_z = points.shape[1] == 3
    xy = points[:, :2]
    z = points[:, 2] if has_z else np.zeros(n)
    iu, ju = np.triu_indices(n, k=1)
    if max_pairs is not None and iu.size > max_pairs:
        rng = np.random.default_rng(_PAIR_SUBSAMPLE_SEED)
        keep = rng.choice(iu.size, size=int(max_pairs), replace=False)
        keep.sort()
        iu, ju = iu[keep], ju[keep]
    dxy = xy[iu] - xy[ju]
    dz = z[iu] - z[ju]
    dv = values[iu] - values[ju]
    return dxy, dz, dv, np.sqrt(np.sum(dxy ** 2, axis=-1)), has_z


def _azimuth_mask(dxy: np.ndarray, azimuth_deg: float, tolerance_deg: float) -> np.ndarray:
    """+Y-dən saat əqrəbi ilə ölçülmüş azimut süzgəci (0°=Şimal)."""
    ang = np.degrees(np.arctan2(dxy[:, 0], dxy[:, 1])) % 180.0
    target = azimuth_deg % 180.0
    delta = np.minimum(np.abs(ang - target), 180.0 - np.abs(ang - target))
    return delta <= tolerance_deg


def _dip_mask(horiz: np.ndarray, dz: np.ndarray, dip_deg: float,
              tolerance_deg: float) -> np.ndarray:
    """Cütün üfüqi müstəvidən çıxma bucağı `dip_deg` ± `tolerance_deg`
    aralığındadırmı (A3.2 — 3D istiqamətli variogram).

    Bucaq işarəsizdir (cüt istiqaməti simmetrikdir): `|arctan(Δz/Δh)|`.
    Tam şaquli cüt (Δh=0) 90° sayılır."""
    with np.errstate(divide="ignore", invalid="ignore"):
        ang = np.degrees(np.arctan2(np.abs(dz), horiz))
    ang = np.where(np.isfinite(ang), ang, 90.0)
    return np.abs(ang - abs(dip_deg)) <= tolerance_deg


def experimental_variogram(points: np.ndarray, values: np.ndarray, n_lags: int = 12,
                           max_lag: Optional[float] = None,
                           azimuth_deg: Optional[float] = None,
                           azimuth_tolerance_deg: float = 22.5,
                           vertical: bool = False,
                           horizontal_tolerance: Optional[float] = None,
                           lag_tolerance: Optional[float] = None,
                           dip_deg: Optional[float] = None,
                           dip_tolerance_deg: float = 22.5,
                           min_pairs: int = 1,
                           bandwidth: Optional[float] = None,
                           vertical_tolerance: Optional[float] = None,
                           max_pairs: Optional[int] = MAX_VARIOGRAM_PAIRS
                           ) -> ExperimentalVariogram:
    """Cüt-cüt lag məsafələrindən binlənmiş yarım-dəyişkənlik hesablayır::

        γ(h) = 1/(2·N(h)) · Σ [Z(xᵢ) − Z(xᵢ+h)]²

    `azimuth_deg` verilməyəndə (defolt) OMNİDİREKSİONAL — bütün
    istiqamətlər eyni binə düşür. Verilsə, yalnız (X,Y) müstəvisində həmin
    azimutdan (±`azimuth_tolerance_deg`) kənara çıxmayan cütlər saxlanılır.

    `dip_deg` (A3.2/A3.3) — 3D istiqamətli variogram: cütün üfüqi
    müstəvidən çıxma bucağı `dip_deg ± dip_tolerance_deg` olmalıdır
    (`dip_deg=0` → üfüqi variogram, `dip_deg=90` → şaquli). `azimuth_deg`
    ilə BİRLİKDƏ işlədilə bilər.

    `lag_tolerance` (A3.1) — bin YARIM-ENİ. Verilməyəndə binlər bitişikdir
    (`tol = lag_width/2`, əvvəlki davranış); verilsə binlər ÜST-ÜSTƏ düşə
    (tol > width/2) və ya ARALIQ buraxa (tol < width/2) bilər — hər cüt
    ÖZ mərkəzinə ən yaxın binə YOX, tolerantlığa düşən HƏR binə qatılır
    (standart geostatistik lag-tolerantlıq davranışı).

    `bandwidth` — istiqamətli axtarışda maksimum yanal (ox-dan) sapma
    (uzunluq vahidində). Uzaq lag-larda bucaq tolerantlığı çox geniş
    "yelpik" verir; bandwidth bunu məhdudlaşdırır. `None` = məhdudiyyət yox.

    `vertical_tolerance` (A3.3, ŞAQULİ BANT ENİ) — ÜFÜQİ/istiqamətli
    variogramda cütün |ΔZ|-si bu həddi keçməməlidir. Laylı rezervuarda
    bu OLMADAN "üfüqi" variogram əslində üfüqi OLMUR: 500 m aralıqda
    ±5° dip pəncərəsi ±44 m şaquli fərqə icazə verir, yəni bir neçə
    LAYI qarışdırır. `None` = məhdudiyyət yox (əvvəlki davranış).

    `max_lag` — bundan UZAQ cütlər ÜMUMİYYƏTLƏ nəzərə alınmır (sonuncu
    binə yığılmır). Verilməyəndə cüt məsafələrinin 75%-lik kvantili
    işlədilir — variogramın etibarlı hissəsi adətən domenin yarısına
    qədərdir.

    `min_pairs` — bu qədər cütü olmayan binlər NƏTİCƏYƏ DAXİL EDİLMİR
    (mənasız `γ` nöqtəsi istehsal edilmir, A3.1).

    `max_pairs` — çox böyük nöqtə çoxluğunda cüt alt-nümunəsi (bax
    `_pair_arrays`); `None` tam hesablamadır.

    `vertical=True` — şaquli variogram: yalnız üfüqi məsafəsi
    `horizontal_tolerance`-dan kiçik olan cütlər (eyni/yaxın quyu, fərqli
    dərinlik) saxlanılır, lag isə |Δz|-dir. Uyğun cüt tapılmazsa `ValueError`
    atılır — sükutla boş nəticə YOX.
    """
    points = np.asarray(points, float)
    values = np.asarray(values, float)
    if points.ndim != 2 or points.shape[1] not in (2, 3):
        raise ValueError(f"Nöqtələr (n,2) və ya (n,3) olmalıdır, alındı: {points.shape}")
    n = points.shape[0]
    if n < 4:
        raise ValueError(f"Deneysel variogram üçün ən azı 4 nöqtə lazımdır (tapıldı: {n}).")
    if values.shape[0] != n:
        raise ValueError(
            f"points ({n}) və values ({values.shape[0]}) uzunluğu uyğun gəlmir.")
    if n_lags < 1:
        raise ValueError(f"n_lags ≥ 1 olmalıdır, alındı: {n_lags}")

    dxy, dz, dv, horiz_dist, has_z = _pair_arrays(points, values, max_pairs)
    warnings: List[str] = []

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
        keep = np.ones(horiz_dist.shape, dtype=bool)
        if azimuth_deg is not None:
            keep &= _azimuth_mask(dxy, float(azimuth_deg), float(azimuth_tolerance_deg))
            if bandwidth is not None:
                # ox-dan yanal sapma: |d| · sin(Δbucaq) — birbaşa hesablanır
                ang = np.degrees(np.arctan2(dxy[:, 0], dxy[:, 1])) % 180.0
                target = float(azimuth_deg) % 180.0
                delta = np.minimum(np.abs(ang - target), 180.0 - np.abs(ang - target))
                keep &= horiz_dist * np.sin(np.radians(delta)) <= float(bandwidth)
        if dip_deg is not None:
            if not has_z:
                raise ValueError(
                    "dip_deg yalnız (n,3) [X,Y,Z] nöqtələrlə mənalıdır — Z sütunu yoxdur.")
            keep &= _dip_mask(horiz_dist, dz, float(dip_deg), float(dip_tolerance_deg))
        if vertical_tolerance is not None:
            if not has_z:
                raise ValueError(
                    "vertical_tolerance yalnız (n,3) [X,Y,Z] nöqtələrlə mənalıdır "
                    "— Z sütunu yoxdur.")
            keep &= np.abs(dz) <= float(vertical_tolerance)
        h = np.sqrt(horiz_dist ** 2 + dz ** 2) if has_z else horiz_dist.copy()
        h = h[keep]
        dv = dv[keep]
        if h.size == 0:
            raise ValueError(
                f"Azimut {azimuth_deg}° (±{azimuth_tolerance_deg}°) / dip {dip_deg}° "
                "daxilində heç bir cüt tapılmadı — bu istiqamətdə nöqtə sıxlığı "
                "kifayət deyil.")

    if max_lag is None:
        max_lag = float(np.percentile(h, 75)) if h.size else 1.0
    max_lag = max(float(max_lag), 1e-9)
    lag_width = max_lag / n_lags
    half = lag_width / 2.0 if lag_tolerance is None else float(lag_tolerance)
    if half <= 0.0:
        raise ValueError(f"lag_tolerance müsbət olmalıdır, alındı: {lag_tolerance!r}")

    # `max_lag`-dan UZAQ cütlər ATILIR — sonuncu binə YIĞILMIR.
    # (Əvvəllər `min(bin, n_lags-1)` ilə bütün uzaq cütlər sonuncu binə
    # düşürdü: bin öz mərkəzindən qat-qat böyük məsafələri təmsil edir və
    # cüt sayı — yəni fit ÇƏKİSİ — nəhəng olur. Nəticədə fit süni uzun
    # radiusa çəkilirdi; bu, istiqamətli radiusların müqayisəsini
    # tamamilə etibarsız edən REAL qüsur idi.)
    within = h <= max_lag
    h, dv = h[within], dv[within]
    if h.size == 0:
        raise ValueError(
            f"max_lag={max_lag:.4g} daxilində heç bir cüt qalmadı — hədd çox kiçikdir.")
    n_pairs_total = int(h.size)

    lags = (np.arange(n_lags) + 0.5) * lag_width
    gamma = np.zeros(n_lags)
    counts = np.zeros(n_lags, dtype=int)

    if lag_tolerance is None:
        # bitişik binlər — sürətli, tək keçidli yol
        bin_idx = np.minimum((h / lag_width).astype(int), n_lags - 1)
        for b in range(n_lags):
            in_bin = bin_idx == b
            c = int(np.sum(in_bin))
            counts[b] = c
            if c > 0:
                gamma[b] = 0.5 * np.mean(dv[in_bin] ** 2)
    else:
        # açıq tolerantlıq — bir cüt bir neçə binə düşə bilər
        for b in range(n_lags):
            in_bin = np.abs(h - lags[b]) <= half
            c = int(np.sum(in_bin))
            counts[b] = c
            if c > 0:
                gamma[b] = 0.5 * np.mean(dv[in_bin] ** 2)

    result = ExperimentalVariogram(
        lags, gamma, counts, azimuth_deg, n_pairs_total, lag_width,
        lag_tolerance, dip_deg, vertical, max_lag, warnings)

    if min_pairs > 1:
        dropped = int(np.sum((counts > 0) & (counts < min_pairs)))
        result = result.valid(min_pairs=min_pairs)
        if dropped:
            result.warnings.append(
                f"{dropped} lag-bin `min_pairs={min_pairs}` həddindən az cütə malik "
                "olduğu üçün ÇIXARILDI (statistik mənasız γ nöqtəsi istehsal edilmir).")
    return result


def directional_variograms(points: np.ndarray, values: np.ndarray,
                           azimuths: Sequence[float] = (0.0, 45.0, 90.0, 135.0),
                           **kwargs) -> Dict[float, ExperimentalVariogram]:
    """Bir neçə azimut üçün deneysel variogram (A3.2).

    Cüt tapılmayan istiqamət NƏTİCƏDƏ OLMUR (səssiz boş bin YOX) —
    çağıran hansı istiqamətlərin hesablandığını açarlardan görür.
    `azimuth_tolerance_deg` verilməyibsə `90/len(azimuths)` (bitişik,
    üst-üstə düşməyən sektorlar) işlədilir."""
    kwargs.setdefault("azimuth_tolerance_deg", 90.0 / max(len(azimuths), 1))
    result: Dict[float, ExperimentalVariogram] = {}
    for az in azimuths:
        try:
            result[float(az)] = experimental_variogram(
                points, values, azimuth_deg=float(az), **kwargs)
        except ValueError:
            continue
    return result


def vertical_variogram(points: np.ndarray, values: np.ndarray,
                       horizontal_tolerance: Optional[float] = None,
                       **kwargs) -> ExperimentalVariogram:
    """Şaquli (laylararası) variogram — `experimental_variogram(vertical=True)`
    üçün adlandırılmış qısayol (A3.3)."""
    return experimental_variogram(points, values, vertical=True,
                                  horizontal_tolerance=horizontal_tolerance, **kwargs)


# ── model fitting ────────────────────────────────────────────────────────
class VariogramValidationError(ValueError):
    """Riyazi olaraq ETİBARSIZ variogram parametrləri (A3.7).

    Belə model Kriging matrisini müsbət-müəyyən olmayan edər (mənfi sill,
    mənfi radius və s.), ona görə solver-ə ÇATMADAN AÇIQ atılır."""


@dataclass
class VariogramParameters:
    model: str
    nugget: float
    sill: float
    range_: float
    weighted_rmse: float
    n_pairs: int
    warnings: List[str] = field(default_factory=list)
    #: fit-in hansı konfiqurasiya ilə alındığı (A3.6 — nəticə saxlanılır)
    n_lags_used: int = 0
    azimuth_deg: Optional[float] = None
    dip_deg: Optional[float] = None
    lag_width: float = 0.0
    lag_tolerance: Optional[float] = None
    vertical: bool = False

    @property
    def total_sill(self) -> float:
        """`C(0) = nugget + sill` (bax modul docstring-i)."""
        return float(self.nugget) + float(self.sill)

    @property
    def effective_range(self) -> float:
        """Modelin e-qatlama radiusu — `range_` praktiki radiusdur."""
        return effective_range(self.range_, self.model)

    def semivariance(self, h) -> np.ndarray:
        return semivariance(h, self.nugget, self.sill, self.range_, self.model)

    def covariance(self, h) -> np.ndarray:
        return covariance(h, self.nugget, self.sill, self.range_, self.model)

    def validate(self, strict: bool = True) -> List[str]:
        """Bax `validate_variogram_parameters` — metod formasında."""
        return validate_variogram_parameters(self.model, self.nugget, self.sill,
                                             self.range_, strict=strict)


#: QAUSS modeli üçün MİNİMUM nugget (sillə nisbətdə) — bax
#: `stabilizing_nugget()`. GSLIB-in klassik tövsiyəsi (~1%).
GAUSSIAN_MIN_NUGGET_RATIO = 0.01


def stabilizing_nugget(model: str, nugget: float, sill: float,
                       ratio: float = GAUSSIAN_MIN_NUGGET_RATIO
                       ) -> Tuple[float, Optional[str]]:
    """QAUSS modelinin Kriging sistemini SABİTLƏŞDİRƏN minimum nugget.

    PROBLEM (bu repozitoriyada ÖLÇÜLÜB, nəzəri iddia deyil): qauss
    variogramı başlanğıcda PARABOLİKDİR (`γ ∝ h²`), ona görə bir-birinə
    yaxın nöqtələrin kovariasiya matrisi demək olar TƏKİLDİR. Sistem
    "uğurla" həll olunur (`Σwᵢ = 1` ödənir, solver `direct` bildirir),
    AMMA çəkilər böyük müsbət/mənfi qiymətlərlə OSSİLYASİYA edir və
    qiymət məlumat diapazonundan kənara çıxır.

    ÖLÇÜLMÜŞ nümunə (60 nöqtə, `range=557`, `sill=2.6e-5`, 24 qonşu;
    məlumat aralığı `[0.1508, 0.2100]`)::

        nugget/sill = 0       → qiymət aralığı [−0.300, +1.126]   ← YARARSIZ
        nugget/sill = 1e-6    → [−0.089, 0.794]
        nugget/sill = 1e-4    → [ 0.028, 0.295]
        nugget/sill = 1e-2    → [ 0.147, 0.225]                   ← sabit
        (eksponensial/sferik model HƏR nugget üçün sabitdir)

    HƏLL: qauss modeli seçiləndə nugget sillin ən azı `ratio` hissəsinə
    qaldırılır. Bu, standart geostatistik təcrübədir (Deutsch & Journel,
    GSLIB: qauss modeli HƏMİŞƏ kiçik nugget-lə işlədilir) və AÇIQ
    xəbərdarlıqla bildirilir — səssiz düzəliş DEYİL.

    Qaytarır: `(nugget, xəbərdarlıq|None)`. `ratio = 0` bu davranışı
    tamamilə söndürür (istifadəçinin AÇIQ seçimi)."""
    if model != MODEL_GAUSSIAN or ratio <= 0.0:
        return float(nugget), None
    minimum = float(ratio) * float(sill)
    if nugget >= minimum:
        return float(nugget), None
    return minimum, (
        f"Qauss variogramı üçün nugget {nugget:.4g} → {minimum:.4g} qaldırıldı "
        f"(sillin {ratio * 100:.3g}%-i). Səbəb: qauss modeli başlanğıcda "
        "parabolikdir, sıfıra yaxın nugget-lə Kriging matrisi demək olar təkil "
        "olur və çəkilər ossilyasiya edərək qiyməti məlumat diapazonundan "
        "kənara çıxarır (ölçülüb). `gaussian_min_nugget_ratio=0` ilə söndürülə bilər.")


def validate_variogram_parameters(model: str, nugget: float, sill: float, range_: float,
                                  strict: bool = True) -> List[str]:
    """Kriging üçün RİYAZİ ETİBARLILIQ yoxlaması (A3.7).

    Atılır (`VariogramValidationError`):
        * naməlum model adı
        * qeyri-sonlu (NaN/±inf) parametr
        * mənfi nugget
        * müsbət olmayan sill (`sill ≤ 0` → struktur yoxdur, matris
          deqenerativ olur)
        * müsbət olmayan radius

    Xəbərdarlıq (siyahı kimi qaytarılır, ATILMIR — model hələ etibarlıdır):
        * nugget ümumi sillin >90%-i (demək olar tam təsadüfi sahə)
        * `strict=False` olanda yuxarıdakı bəzi hallar üçün izah

    Heç bir parametr SƏSSİZ DÜZƏLDİLMİR — düzəlişi `repair_variogram_
    parameters()` AÇIQ şəkildə edir."""
    warnings: List[str] = []
    if model not in MODEL_FUNCS:
        raise VariogramValidationError(
            f"Naməlum variogram modeli: {model!r}. Dəstəklənən: {KNOWN_MODELS}")
    values = {"nugget": nugget, "sill": sill, "range_": range_}
    for name, value in values.items():
        if not np.isfinite(value):
            raise VariogramValidationError(
                f"Variogram parametri '{name}' sonlu olmalıdır, alındı: {value!r}")
    if nugget < 0.0:
        raise VariogramValidationError(f"nugget mənfi ola bilməz, alındı: {nugget!r}")
    if sill <= 0.0:
        raise VariogramValidationError(
            f"sill (quruluşlu hissə) müsbət olmalıdır, alındı: {sill!r} — "
            "sıfır sill məkan strukturunu tamamilə silir, Kriging matrisi deqenerativ olur.")
    if range_ <= 0.0:
        raise VariogramValidationError(f"range_ müsbət olmalıdır, alındı: {range_!r}")

    total = nugget + sill
    if total > 0.0 and nugget / total > 0.9:
        warnings.append(
            f"nugget ümumi sillin {100.0 * nugget / total:.0f}%-idir — məkan "
            "korrelyasiyası demək olar yoxdur, Kriging nəticəsi qlobal ortaya yaxın olacaq.")
    if not strict:
        return warnings
    return warnings


def repair_variogram_parameters(model: str, nugget: float, sill: float, range_: float
                                ) -> Tuple[float, float, float, List[str]]:
    """Etibarsız parametrləri AÇIQ (xəbərdarlıqlı) şəkildə düzəldir.

    Yalnız çağıran "fit nəticəsini nə olursa olsun işlət" dedikdə
    istifadə olunur — defolt yol `validate_variogram_parameters()`-in
    XƏTA ATMASIDIR. Qaytarır: `(nugget, sill, range_, warnings)`."""
    warnings: List[str] = []
    if model not in MODEL_FUNCS:
        raise VariogramValidationError(
            f"Naməlum variogram modeli: {model!r}. Dəstəklənən: {KNOWN_MODELS}")
    if not np.isfinite(nugget) or nugget < 0.0:
        warnings.append(f"Etibarsız nugget ({nugget!r}) → 0.0 ilə əvəz edildi.")
        nugget = 0.0
    if not np.isfinite(sill) or sill <= 0.0:
        warnings.append(f"Etibarsız sill ({sill!r}) → 1e-12 ilə əvəz edildi.")
        sill = 1e-12
    if not np.isfinite(range_) or range_ <= 0.0:
        warnings.append(f"Etibarsız range_ ({range_!r}) → 1e-6 ilə əvəz edildi.")
        range_ = 1e-6
    return float(nugget), float(sill), float(range_), warnings


def _fit_single_model(lags: np.ndarray, gamma: np.ndarray, counts: np.ndarray,
                      model: str, fix_nugget: Optional[float],
                      fix_sill: Optional[float] = None) -> VariogramParameters:
    """Çəkili (çəki = √cüt sayı) ən kiçik kvadratlarla bir model fiti.

    `fix_nugget`/`fix_sill` — həmin parametr SABİT saxlanılır, yalnız
    qalanları optimallaşdırılır. `fix_sill` xüsusilə İSTİQAMƏTLİ fit
    üçün vacibdir (bax `detect_anisotropy`): GEOMETRİK anizotropluqda
    sill istiqamətdən ASILI DEYİL, yalnız radius dəyişir — sill də
    sərbəst buraxılsa, sill/radius bir-birini əvəz edə bilir və fit
    zəif müəyyən (ill-posed) olur.
    """
    func = MODEL_FUNCS[model]
    weights = np.sqrt(counts.astype(float))
    sill0 = float(fix_sill) if fix_sill is not None else (
        float(gamma.max()) if gamma.max() > 0 else 1.0)
    scale = max(sill0, 1e-12)
    # 95%-lik sill həddinə çatan ən kiçik lag — praktiki radiusun ilkin təxmini
    above = np.where(gamma >= 0.95 * scale)[0]
    range0 = float(lags[above[0]]) if above.size else float(lags[-1])
    range0 = max(range0, float(lags[0]) if lags.size else 1e-6, 1e-6)
    nugget0 = (float(fix_nugget) if fix_nugget is not None else
               float(np.clip(gamma[0] - (gamma[1] - gamma[0]) if gamma.size > 1
                             else gamma[0], 0.0, scale)))

    range_hi = max(lags[-1] * 10.0, range0 * 10.0)
    x0: List[float] = []
    lo: List[float] = []
    hi: List[float] = []
    if fix_nugget is None:
        x0.append(nugget0); lo.append(0.0); hi.append(scale * 4.0 + 1e-9)
    if fix_sill is None:
        x0.append(sill0); lo.append(1e-12); hi.append(scale * 4.0 + 1e-9)
    x0.append(range0); lo.append(1e-6); hi.append(range_hi)

    def unpack(p):
        cursor = 0
        if fix_nugget is None:
            nugget = float(p[cursor]); cursor += 1
        else:
            nugget = float(fix_nugget)
        if fix_sill is None:
            sill = float(p[cursor]); cursor += 1
        else:
            sill = float(fix_sill)
        return nugget, sill, float(p[cursor])

    def resid(p):
        nugget, sill, range_ = unpack(p)
        return weights * (func(lags, nugget, sill, range_) - gamma)

    result = least_squares(resid, x0, bounds=(lo, hi))
    nugget, sill, range_ = unpack(result.x)

    fitted = func(lags, nugget, sill, range_)
    rmse = float(np.sqrt(np.average((fitted - gamma) ** 2, weights=weights ** 2)))         if np.any(weights > 0) else float(np.sqrt(np.mean((fitted - gamma) ** 2)))
    return VariogramParameters(model=model, nugget=nugget, sill=max(sill, 1e-12),
                               range_=max(range_, 1e-6), weighted_rmse=rmse,
                               n_pairs=int(counts.sum()), n_lags_used=int(lags.size))


def _stamp_configuration(params: VariogramParameters,
                         exp: ExperimentalVariogram) -> VariogramParameters:
    """Fit nəticəsinə DENEYSEL konfiqurasiyanı yazır (A3.6)."""
    params.azimuth_deg = exp.azimuth_deg
    params.dip_deg = exp.dip_deg
    params.lag_width = exp.lag_width
    params.lag_tolerance = exp.lag_tolerance
    params.vertical = exp.vertical
    return params


def compare_variogram_models(experimental: ExperimentalVariogram,
                             models: Sequence[str] = KNOWN_MODELS,
                             fix_nugget: Optional[float] = None,
                             min_pairs: int = 1,
                             fix_sill: Optional[float] = None
                             ) -> Dict[str, VariogramParameters]:
    """Hər namizəd model üçün AYRICA fit + çəkili RMSE (A3.5).

    Nəticə lüğəti çağırana bütün namizədləri (təkcə qalibi yox) verir ki,
    model seçimi ŞƏFFAF olsun və hesabatda göstərilə bilsin. Fit alınmayan
    model nəticədə OLMUR."""
    exp = experimental.valid(min_pairs=min_pairs)
    if exp.lags.size < 3:
        raise ValueError(
            f"Variogram model fit üçün ən azı 3 dolu lag-bin lazımdır "
            f"(tapıldı: {exp.lags.size}, min_pairs={min_pairs}). "
            "Nöqtə sayı/sıxlığı kifayət deyil.")
    unknown = [m for m in models if m not in KNOWN_MODELS]
    if unknown:
        raise ValueError(f"Naməlum variogram modeli: {unknown}. Dəstəklənən: {KNOWN_MODELS}")

    results: Dict[str, VariogramParameters] = {}
    for m in models:
        try:
            candidate = _fit_single_model(exp.lags, exp.gamma, exp.counts, m,
                                          fix_nugget, fix_sill)
        except (ValueError, np.linalg.LinAlgError):
            continue
        try:
            candidate.warnings.extend(candidate.validate())
        except VariogramValidationError:
            # etibarsız fit NAMİZƏD SİYAHISINA BURAXILMIR (A3.7)
            continue
        results[m] = _stamp_configuration(candidate, exp)
    return results


def fit_variogram(experimental: ExperimentalVariogram, model: str = "auto",
                  fix_nugget: Optional[float] = None,
                  min_pairs: int = 1,
                  fix_sill: Optional[float] = None) -> VariogramParameters:
    """Deneysel variograma nugget/sill/range fit edir.

    `model="auto"` — sferik/eksponensial/qauss ÜÇÜNÜN də fit edilib, çəkili
    RMSE-si ən kiçik olanı seçilir (bu, "ən hamar görünən" seçimi DEYİL —
    ədədi uyğunluğa əsaslanır). Konkret model adı verilsə, yalnız o fit
    edilir. Nəticə HƏMİŞƏ `validate_variogram_parameters()`-dən keçir —
    etibarsız model Kriging-ə ÇATA BİLMİR (A3.7).
    """
    if model != "auto" and model not in KNOWN_MODELS:
        raise ValueError(f"Naməlum variogram modeli: {model!r}. Dəstəklənən: {KNOWN_MODELS}")
    candidates = KNOWN_MODELS if model == "auto" else (model,)
    results = compare_variogram_models(experimental, candidates, fix_nugget,
                                       min_pairs, fix_sill)
    if not results:
        raise ValueError(
            f"Heç bir namizəd model ({', '.join(candidates)}) üçün ETİBARLI fit "
            "alınmadı — deneysel variogram struktursuzdur və ya nöqtə sayı azdır.")
    return min(results.values(), key=lambda p: p.weighted_rmse)


def fit_variogram_from_data(points: np.ndarray, values: np.ndarray, n_lags: int = 12,
                            max_lag: Optional[float] = None, model: str = "auto",
                            fix_nugget: Optional[float] = None,
                            min_pairs: int = 1,
                            azimuth_deg: Optional[float] = None,
                            lag_tolerance: Optional[float] = None
                            ) -> VariogramParameters:
    """`experimental_variogram` + `fit_variogram` bir addımda.

    `azimuth_deg` verilməyəndə omnidireksional (əvvəlki davranış)."""
    exp = experimental_variogram(points, values, n_lags=n_lags, max_lag=max_lag,
                                 azimuth_deg=azimuth_deg, lag_tolerance=lag_tolerance)
    return fit_variogram(exp, model=model, fix_nugget=fix_nugget, min_pairs=min_pairs)


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


# ── anizotropluq aşkarlanması (həndəsə `geology/anisotropy.py`-dədir) ───
@dataclass
class AnisotropyDetectionResult:
    azimuth_deg: float
    range_major: float
    range_minor: float
    ratio: float
    directional_ranges: Dict[float, float]
    reliable: bool
    warnings: List[str] = field(default_factory=list)
    range_vertical: Optional[float] = None

    def to_params(self, range_vertical: Optional[float] = None) -> AnisotropyParams:
        """Aşkarlanma nəticəsindən HƏNDƏSƏ obyekti (A4.4).

        Şaquli radius aşkarlanmayıbsa (`None`) `range_major` işlədilir —
        yəni şaquli anizotropluq TƏTBİQ EDİLMİR (uydurulmuş nisbət YOX)."""
        vertical = (range_vertical if range_vertical is not None
                    else (self.range_vertical if self.range_vertical is not None
                          else self.range_major))
        return AnisotropyParams(azimuth_deg=self.azimuth_deg,
                                range_major=self.range_major,
                                range_minor=self.range_minor,
                                range_vertical=vertical)


def detect_anisotropy(points: np.ndarray, values: np.ndarray, n_directions: int = 6,
                      n_lags: int = 8, max_lag: Optional[float] = None,
                      azimuth_tolerance_deg: Optional[float] = None,
                      min_pairs_per_direction: int = 10,
                      detect_vertical: bool = False,
                      vertical_horizontal_tolerance: Optional[float] = None,
                      bandwidth: Optional[float] = None
                      ) -> AnisotropyDetectionResult:
    """`n_directions` üfüqi istiqamətdə radius təxmin edib ən böyüyünü
    major ox kimi seçir (A3.2/A4.4).

    METODİKA (bu, sadə "hər istiqaməti ayrıca fit et" yanaşmasından
    fərqlidir və QƏSDƏN belədir):

    1. Əvvəlcə OMNİDİREKSİONAL fit aparılır və oradan ÜMUMİ SİLL
       (`nugget + sill`) götürülür. İstiqamətli fitlərdə bu sill SABİT,
       nugget isə SIFIR saxlanılır — yeganə sərbəst parametr RADİUSDUR.
       Səbəb: GEOMETRİK anizotropluqda ümumi dispersiya istiqamətdən
       ASILI DEYİL, yalnız davamlılıq radiusu dəyişir. Sill (və ya
       nugget) da sərbəst buraxılsa, onlar radiusu kompensasiya edir və
       fit zəif müəyyən (ill-posed) olur — praktikada bu, izotrop sahədə
       belə 100 dəfə fərqli "radiuslar", dönmüş sahədə isə 45-75° azimut
       xətası verirdi (ölçülüb). Bu qayda ilə üç fərqli dönmə bucağı
       (45°/60°/135°) DƏQİQ bərpa olunur.

       DİQQƏT: buradakı istiqamətli radius MÜQAYİSƏ ÜÇÜN ölçüdür
       (nugget=0 fərziyyəsi ilə), tam model deyil — nugget/sill
       parçalanması omnidireksional (və ya istifadəçinin öz) fitindən
       gəlir.
    2. Bütün istiqamətlər EYNİ `max_lag`/lag-eni ilə binlənir (omni
       hesablamasından götürülür) — əks halda binlər müqayisə edilə
       bilməz.
    3. Müşahidə HÜDUDU: `max_lag`-dan uzağı məlumatda YOXDUR, ona görə
       istiqamətli radius `max_lag`-a KƏSİLİR və bu, `warnings`-də açıq
       bildirilir — "radius 9600 m" kimi uydurma ədəd qaytarılmır.

    Quyu sayı az olanda (tipik reservoir modelində 5-20 quyu) istiqamətli
    binlərdə kifayət qədər cüt olmaya bilər — bu vəziyyət GİZLƏDİLMİR:
    `reliable=False` və izahlı `warnings` qaytarılır, izotrop defoltla
    (major=minor=omnidireksional radius) geri qayıdılır.

    `detect_vertical=True` (A3.3) — əlavə olaraq ŞAQULİ variogram da fit
    edilir; alınmazsa `range_vertical=None` qalır (uydurulmur).
    """
    points = np.asarray(points, float)
    if azimuth_tolerance_deg is None:
        azimuth_tolerance_deg = 90.0 / n_directions

    warnings: List[str] = []
    omni: Optional[VariogramParameters] = None
    common_max_lag = max_lag
    try:
        omni_exp = experimental_variogram(points, values, n_lags=n_lags, max_lag=max_lag,
                                          bandwidth=bandwidth)
        common_max_lag = omni_exp.max_lag
        omni = fit_variogram(omni_exp, model=MODEL_SPHERICAL)
    except ValueError as exc:
        warnings.append(f"Omnidireksional fit alınmadı ({exc}) — sill sabitlənmədi.")

    range_vertical: Optional[float] = None
    if detect_vertical:
        try:
            exp_v = vertical_variogram(
                points, values, horizontal_tolerance=vertical_horizontal_tolerance,
                n_lags=n_lags)
            range_vertical = fit_variogram(exp_v, model=MODEL_SPHERICAL).range_
        except ValueError as exc:
            warnings.append(f"Şaquli radius aşkarlanmadı: {exc}")

    azimuths = np.linspace(0.0, 180.0, n_directions, endpoint=False)
    directional_ranges: Dict[float, float] = {}
    clamped: List[float] = []
    for az in azimuths:
        try:
            exp = experimental_variogram(points, values, n_lags=n_lags,
                                         max_lag=common_max_lag,
                                         azimuth_deg=float(az),
                                         azimuth_tolerance_deg=azimuth_tolerance_deg,
                                         bandwidth=bandwidth)
        except ValueError:
            continue
        if exp.n_pairs_total < min_pairs_per_direction:
            continue
        try:
            fit = fit_variogram(
                exp, model=MODEL_SPHERICAL,
                fix_nugget=None if omni is None else 0.0,
                fix_sill=None if omni is None else omni.total_sill)
        except ValueError:
            continue
        limit = common_max_lag if common_max_lag else fit.range_
        if fit.range_ > limit:
            clamped.append(float(az))
        directional_ranges[float(az)] = min(fit.range_, limit)

    if clamped:
        warnings.append(
            f"{len(clamped)} istiqamətdə ({', '.join(f'{a:.0f}°' for a in clamped)}) "
            f"variogram müşahidə həddi (max_lag={common_max_lag:.4g}) daxilində sillə "
            "çatmadı — radius bu hədlə MƏHDUDLAŞDIRILDI (uydurma ekstrapolyasiya yox).")

    if len(directional_ranges) < 2:
        warnings.insert(0,
            "İstiqamətli variogram üçün kifayət qədər cüt yoxdur (nöqtə sayı/sıxlığı "
            "az) — anizotropluq etibarlı aşkarlanmadı, izotrop defolt işlədilir.")
        if omni is not None:
            range_omni = omni.range_
        else:
            span = float(np.sqrt(np.sum((points[:, :2].max(axis=0)
                                         - points[:, :2].min(axis=0)) ** 2)))
            range_omni = max(span / 3.0, 1e-6)
        return AnisotropyDetectionResult(azimuth_deg=0.0, range_major=range_omni,
                                         range_minor=range_omni, ratio=1.0,
                                         directional_ranges=directional_ranges,
                                         reliable=False, warnings=warnings,
                                         range_vertical=range_vertical)

    spread = np.asarray(list(directional_ranges.values()), float)
    if spread.max() > 0.0 and (spread.max() - spread.min()) <= 1e-9 * spread.max():
        warnings.append(
            "Bütün istiqamətlərdə radius EYNİ alındı — bu lag ayırdetməsi ilə "
            "anizotropluq HƏLL EDİLƏ BİLMİR (major ox seçimi mənasızdır). "
            "n_lags-i artırın, max_lag-i kiçildin və ya `bandwidth` verin.")

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
                                     reliable=True, warnings=warnings,
                                     range_vertical=range_vertical)
