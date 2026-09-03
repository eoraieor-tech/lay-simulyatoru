"""Məkan interpolyasiyası — IPropertyInterpolator implementasiyaları (A1/A5).

Üç üsul, artan mürəkkəblik sırası ilə:

    NearestNeighbour  — ən yaxın quyunun dəyəri. Sürətli, pilləli.
    InverseDistance   — 1/d^p çəkiləri. Hamar, amma quyular arasında
                        həmişə orta dəyərə meyl edir ("öküz gözü" effekti).
    OrdinaryKriging   — variogram əsaslı, statistik optimal xətti qiymət.
                        Məkan korrelyasiyasını nəzərə alır.

BORU XƏTTİ (A5) — `OrdinaryKriging` bir çağırışda bunları icra edir::

    XAM SƏRT DATA
      → DOĞRULAMA (`_prepare`: NaN/±inf, uzunluq, dublikat siyasəti)
      → KOORDİNAT NORMALLAŞDIRMASI ((n,2) → (n,3))
      → DENEYSEL VARİOGRAM + MODEL FİT (`geology/variogram.py`, auto_fit)
      → PARAMETR DOĞRULAMASI (`validate_variogram_parameters`)
      → ANİZOTROP TRANSFORMASİYA (`geology/anisotropy.py`)
      → MƏKAN İNDEKSİ + QONŞULUQ SEÇİMİ (`geology/spatial_search.py`)
      → YERLİ KRİGİNG (kiçik, ədədi dayanıqlı sistemlər)
      → QİYMƏT + KRİGİNG VARİANSI + DİAQNOSTİKA (`KrigingResult`)
      → İNTERPOLYASİYA EDİLMİŞ ŞƏBƏKƏ

Həndəsə (anizotrop məsafə) YALNIZ `geology/anisotropy.py`-dən gəlir,
qonşuluq seçimi YALNIZ `geology/spatial_search.py`-dən — yəni "A2 bir
həndəsə, A1 başqa həndəsə işlədir" uyğunsuzluğu struktur olaraq MÜMKÜN
DEYİL (A4.3/Gate 5).

Keçiricilik üçün `log_transform=True` tövsiyə olunur: keçiricilik
log-normal paylanır, ona görə interpolyasiya ln(k) fəzasında aparılır
(bax `interpolate_property` və `ValueTransform` — A6 genişlənmə nöqtəsi).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

from ..interfaces.interpolation import IPropertyInterpolator
from .anisotropy import AnisotropyParams
from .spatial_search import (SUPPORT_EXTRAPOLATED, BatchNeighborhood,
                             NeighborhoodConfig, NeighborhoodSelector)
from .transforms import IDENTITY_TRANSFORM, LOG_TRANSFORM, LogTransform, ValueTransform
from .variogram import (GAUSSIAN_MIN_NUGGET_RATIO, MODEL_FUNCS, MODEL_SPHERICAL,
                        AnisotropyDetectionResult, VariogramParameters,
                        detect_anisotropy, fit_variogram_from_data, stabilizing_nugget,
                        validate_variogram_parameters)

# ── solver statusları (A1.5: hər ehtiyat yolu AÇIQ görünür) ────────────
SOLVER_DIRECT = "direct"              #: LAPACK `solve` — normal yol
SOLVER_JITTER = "jitter"              #: diaqonal requlyarlaşdırma tələb olundu
SOLVER_LSTSQ = "lstsq"                #: minimal-norma (psevdo-tərs) həlli
SOLVER_RENORMALIZED = "renormalized"  #: çəkilər Σw=1-ə yenidən normallandı
SOLVER_IDW = "idw_fallback"           #: sistem tamamilə həll edilmədi → 1/d²
SOLVER_EXACT = "exact_hard_data"      #: hədəf sərt data nöqtəsi ilə üst-üstə
SOLVER_NONE = "no_neighbors"          #: qonşu yoxdur → NaN
SOLVER_SINGLE_VALUE = "single_value"  #: yeganə sərt data nöqtəsi

#: `Σ wᵢ = 1` (yansızlıq) şərtinin qəbul edilən sapması. Bundan böyük
#: sapma sistemin faktiki olaraq həll OLUNMADIĞINI göstərir.
UNBIASED_TOLERANCE = 1e-6

#: Diaqonal requlyarlaşdırma (jitter) cəhdləri — matrisin izinə (trace)
#: NİSBİ əmsallar. Kiçikdən böyüyə sınanır, İLK uğurlu dayandırır.
JITTER_FACTORS: Tuple[float, ...] = (1e-12, 1e-10, 1e-8, 1e-6, 1e-4)

#: `_max_pairwise_distance` tam (n,n) matris qurmaqdan qaçmağa keçdiyi
#: hədd — bundan aşağıda ƏVVƏLKİ (bit-bit eyni) yol işlədilir.
_FULL_MATRIX_LIMIT = 1500

#: Yerli sistemlərin toplu (batched) həllində bir dəstədə saxlanılan
#: MAKSIMUM matris elementi — yaddaş tavanı (≈ 16 MB double).
_BATCH_ELEMENTS = 2_000_000


def _distance_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """(n,2) və (m,2) arasında Evklid məsafələri -> (n,m)."""
    diff = a[:, None, :] - b[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=-1))


def _max_pairwise_distance(points: np.ndarray) -> float:
    """Ən böyük cüt-cüt məsafə — TAM (n,n) matris QURMADAN (A7/Gate 10).

    Kiçik çoxluqlarda (`n ≤ _FULL_MATRIX_LIMIT`) əvvəlki yol (tam matris)
    saxlanılır ki, `range = domen/3` evristikası BİT-BİT eyni qalsın.
    Böyük çoxluqda əvvəlcə QABARIQ ÖRTÜK (convex hull) təpələri arasında
    maksimum axtarılır — diametr HƏMİŞƏ örtük təpələrində realizə olunur,
    ona görə nəticə TƏQRİBİ DEYİL, DƏQİQDİR. Örtük deqenerativ (kollinear
    nöqtələr) olanda parçalı (chunked) dəqiq hesablamaya keçilir —
    yaddaş `O(n)`, nəticə yenə dəqiq.
    """
    n = points.shape[0]
    if n < 2:
        return 0.0
    if n <= _FULL_MATRIX_LIMIT:
        return float(_distance_matrix(points, points).max())

    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(points)
        vertices = points[hull.vertices]
        if vertices.shape[0] >= 2:
            return float(_distance_matrix(vertices, vertices).max())
    except Exception:   # QhullError daxil — deqenerativ həndəsə
        pass

    chunk = max(1, _BATCH_ELEMENTS // n)
    best = 0.0
    for start in range(0, n, chunk):
        block = points[start:start + chunk]
        best = max(best, float(_distance_matrix(block, points).max()))
    return best


class NearestNeighbour(IPropertyInterpolator):
    name = "Ən yaxın qonşu"

    def interpolate(self, points, values, targets) -> np.ndarray:
        points = np.asarray(points, float).reshape(-1, 2)
        values = np.asarray(values, float)
        targets = np.asarray(targets, float).reshape(-1, 2)
        distances = _distance_matrix(targets, points)
        return values[np.argmin(distances, axis=1)]


@dataclass
class InverseDistance(IPropertyInterpolator):
    """IDW: w_i = 1 / d_i^p."""
    power: float = 2.0
    search_radius: Optional[float] = None
    name: str = "Əks məsafə (IDW)"

    def interpolate(self, points, values, targets) -> np.ndarray:
        points = np.asarray(points, float).reshape(-1, 2)
        values = np.asarray(values, float)
        targets = np.asarray(targets, float).reshape(-1, 2)

        distances = _distance_matrix(targets, points)
        exact = distances < 1e-9
        distances = np.maximum(distances, 1e-9)

        weights = 1.0 / distances ** self.power
        if self.search_radius:
            weights = np.where(distances <= self.search_radius, weights, 0.0)
            empty = weights.sum(axis=1) <= 0.0
            if np.any(empty):
                # radius daxilində quyu yoxdursa, ən yaxınına qayıdırıq
                nearest = np.argmin(distances[empty], axis=1)
                weights[empty] = 0.0
                weights[np.where(empty)[0], nearest] = 1.0

        result = (weights * values[None, :]).sum(axis=1) / weights.sum(axis=1)

        # quyu nöqtəsində dəqiq dəyər qaytarılır
        hit_rows, hit_columns = np.where(exact)
        result[hit_rows] = values[hit_columns]
        return result


@dataclass
class KrigingResult:
    """Bir `krige()` çağırışının TAM nəticəsi (A1.5).

    Bütün massivlər `(m,)` — hədəf sayı qədər. Heç bir sahə UYDURULMUR:
    qonşu tapılmayan hədəfdə `estimate`/`variance` NaN-dır, `solver`
    `SOLVER_NONE` göstərir.
    """

    estimate: np.ndarray
    variance: np.ndarray
    neighbor_count: np.ndarray          #: yerli sistemə daxil olan nöqtə sayı
    nearest_distance: np.ndarray        #: ANİZOTROP fəzada ən yaxın qonşu məsafəsi
    extrapolated: np.ndarray            #: bool — `support == "extrapolated"`
    support: np.ndarray                 #: `spatial_search.SUPPORT_*`
    neighborhood_status: np.ndarray     #: `spatial_search.STATUS_*`
    solver: np.ndarray                  #: `SOLVER_*`
    #: Laqranj vuruğu `μ` — yansızlıq şərtinin qiyməti. Loq-normal geri
    #: çevirmənin ADİ KRİGİNQ üçün düzgün orta düsturu buna ehtiyac duyur
    #: (`exp(ŷ + σ²/2 − μ)`, bax `transforms.BackTransform.MEAN_OK`).
    lagrange: np.ndarray
    anisotropy: AnisotropyParams
    model: str
    range_: float
    sill: float
    nugget: float
    local: bool                         #: yerli (True) yoxsa qlobal (False) sistem
    fit: Optional[VariogramParameters] = None
    warnings: List[str] = field(default_factory=list)

    def __array__(self, dtype=None, copy=None):
        """`np.asarray(result)` → qiymətlər. Mövcud, massiv gözləyən
        çağıranlar üçün geriyə-uyğunluq körpüsü (A12)."""
        if dtype is None:
            return np.array(self.estimate, copy=bool(copy)) if copy else self.estimate
        return self.estimate.astype(dtype, copy=True)

    def __len__(self) -> int:
        return int(self.estimate.size)

    @property
    def standard_deviation(self) -> np.ndarray:
        """`sqrt(variance)` — mənfi (ədədi) varians sıfıra kəsilir."""
        return np.sqrt(np.clip(self.variance, 0.0, None))

    def summary(self) -> str:
        """Qısa mətn hesabatı — UI/PDF üçün."""
        finite = np.isfinite(self.estimate)
        lines = [
            f"Kriging: model={self.model} range={self.range_:.4g} "
            f"sill={self.sill:.4g} nugget={self.nugget:.4g} "
            f"({'yerli' if self.local else 'qlobal'} sistem)",
            f"  hədəf: {self.estimate.size}, etibarlı: {int(finite.sum())}, "
            f"ekstrapolyasiya: {int(np.sum(self.extrapolated))}",
        ]
        if finite.any():
            lines.append(
                f"  qonşu sayı: min {int(self.neighbor_count[finite].min())} "
                f"orta {self.neighbor_count[finite].mean():.1f} "
                f"maks {int(self.neighbor_count[finite].max())}")
        unique, counts = np.unique(self.solver.astype(str), return_counts=True)
        lines.append("  solver: " + ", ".join(f"{u}={c}" for u, c in zip(unique, counts)))
        for message in self.warnings:
            lines.append(f"  ⚠ {message}")
        return "\n".join(lines)


def _check_weights(solution: np.ndarray, n: int) -> bool:
    """Həllin ETİBARLILIĞI: sonlu + yansızlıq şərti `Σ wᵢ = 1`."""
    if not np.all(np.isfinite(solution)):
        return False
    total = solution[:n].sum(axis=0)
    return bool(np.all(np.abs(total - 1.0) <= UNBIASED_TOLERANCE))


def _solve_single_robust(left: np.ndarray, right: np.ndarray, distances: np.ndarray,
                         total_sill: float) -> Tuple[np.ndarray, str]:
    """Bir yerli Kriging sisteminin ÇOX PİLLƏLİ dayanıqlı həlli (A1.3).

    Ardıcıllıq (hər pillə AÇIQ, statusla bildirilir):

    1. `np.linalg.solve` — LAPACK `gesv` (LU + qismən pivotlama). Kor-koranə
       `inv()` İSTİFADƏ EDİLMİR.
    2. **Jitter** — məlumat blokunun diaqonalına matrisin izinə nisbi
       kiçik `ε` əlavə edilir (Tixonov requlyarlaşdırması). Bu, riyazi
       cəhətdən SONSUZ KİÇİK nugget əlavə etməyə bərabərdir — dublikat/
       demək olar üst-üstə düşən nöqtələrin yaratdığı təkliyi (singularity)
       aradan qaldırır. `JITTER_FACTORS` kiçikdən böyüyə sınanır, İLK
       uğurlu dayandırır (lazım olandan artıq requlyarlaşdırma YOX).
    3. **`lstsq`** — minimal-norma (psevdo-tərs) həlli. Yalnız 1-2 uğursuz
       olanda; nəticə yansızlıq şərtini pozarsa çəkilər AÇIQ şəkildə
       yenidən normallanır (`SOLVER_RENORMALIZED`).
    4. **IDW ehtiyatı** — sistem ümumiyyətlə həll edilmirsə `1/d²` çəkiləri.
       Bu, Kriging DEYİL: varians a-priori sill (`nugget+sill`) kimi
       qaytarılır, çünki həqiqi kriging variansı mövcud deyil. Status
       `SOLVER_IDW` ilə HƏMİŞƏ görünür.

    Qaytarır: `(solution (n+1,), status)`.
    """
    n = left.shape[0] - 1

    try:
        solution = np.linalg.solve(left, right)
        if _check_weights(solution, n):
            return solution, SOLVER_DIRECT
    except np.linalg.LinAlgError:
        pass

    scale = float(np.trace(left[:n, :n])) / max(n, 1)
    if not np.isfinite(scale) or scale <= 0.0:
        scale = 1.0
    for factor in JITTER_FACTORS:
        jittered = left.copy()
        jittered[:n, :n] += np.eye(n) * (factor * scale)
        try:
            solution = np.linalg.solve(jittered, right)
        except np.linalg.LinAlgError:
            continue
        if _check_weights(solution, n):
            return solution, SOLVER_JITTER

    solution, *_ = np.linalg.lstsq(left, right, rcond=None)
    if np.all(np.isfinite(solution)):
        if _check_weights(solution, n):
            return solution, SOLVER_LSTSQ
        total = solution[:n].sum()
        if np.isfinite(total) and abs(total) > 1e-12:
            solution[:n] /= total
            return solution, SOLVER_RENORMALIZED

    safe = np.maximum(np.asarray(distances, float), 1e-12)
    weights = 1.0 / safe ** 2
    weights = weights / weights.sum()
    solution = np.empty(n + 1)
    solution[:n] = weights
    # Laqranj vuruğu elə seçilir ki, `σ² = Σ w·γ + μ` a-priori sill versin
    solution[n] = total_sill - float(np.dot(weights, right[:n]))
    return solution, SOLVER_IDW


def _solve_batch(left: np.ndarray, right: np.ndarray, distances: np.ndarray,
                 total_sill: float) -> Tuple[np.ndarray, np.ndarray]:
    """Eyni ölçülü yerli sistemlərin TOPLU həlli (A7 vektorlaşdırma).

    `left` (g, n+1, n+1), `right` (g, n+1), `distances` (g, n).
    Əvvəlcə HAMISI birlikdə `np.linalg.solve` ilə həll olunur (LAPACK
    yığcam çağırışı); YALNIZ yoxlamadan keçməyən sistemlər tək-tək
    `_solve_single_robust`-a göndərilir — yəni ehtiyat yolları
    performansı normal halda HEÇ CÜR yavaşlatmır.
    """
    g, size = left.shape[0], left.shape[1]
    n = size - 1
    status = np.full(g, SOLVER_DIRECT, dtype=object)
    try:
        solution = np.linalg.solve(left, right[:, :, None])[:, :, 0]
    except np.linalg.LinAlgError:
        solution = np.full((g, size), np.nan)

    finite = np.all(np.isfinite(solution), axis=1)
    unbiased = np.zeros(g, dtype=bool)
    unbiased[finite] = np.abs(solution[finite, :n].sum(axis=1) - 1.0) <= UNBIASED_TOLERANCE
    for row in np.where(~(finite & unbiased))[0]:
        solution[row], status[row] = _solve_single_robust(
            left[row], right[row], distances[row], total_sill)
    return solution, status


@dataclass
class OrdinaryKriging(IPropertyInterpolator):
    """Adi kriging — çoxlu variogram modeli, avtomatik fit, tam 3D
    anizotropluq (azimut + dip + major/minor/şaquli), peşəkar yerli
    qonşuluq seçimi və ədədi dayanıqlı solver ilə (A1).

        γ(h) = c0 + c·[1.5·(h/a) − 0.5·(h/a)³],  h < a      (sferik, defolt)
        γ(h) = c0 + c,                            h ≥ a

    Eksponensial/qauss modelləri də var (bax `geology/variogram.py`).
    c0 — nugget, c — QURULUŞLU sill, a — major üfüqi radius (`range_`).

    **Sistem (A1.1)** — N qonşu üçün::

        [Γ  1] [w]   [γ₀]
        [1ᵀ 0] [μ] = [ 1]

    `Γ[i,j] = γ(hᵢⱼ)` (diaqonal 0), `γ₀[i] = γ(hᵢ₀)`, `Σwᵢ = 1` MƏCBURDUR
    (Laqranj vuruğu `μ`). Qiymət `ẑ = Σ wᵢ zᵢ`, varians `σ² = Σ wᵢγ₀ᵢ + μ`
    — hər ikisi EYNİ həll edilmiş sistemdən oxunur (A1.6 konvensiya
    vahidliyi).

    **Parametr mənbəyi** — `range_`/`sill` birbaşa verilə bilər; verilməyib
    `auto_fit=True` olsa `geology/variogram.py`-dəki deneysel variogram +
    çəkili ən kiçik kvadrat fit işə düşür (bax `fit_variogram_from_data`);
    fit mümkün olmasa (nöqtə azdır) AÇIQ xəbərdarlıqla köhnə `domen/3`
    evristikasına geri qayıdılır. `auto_fit=False` (defolt) olanda dəyişən
    YOXDUR — köhnə `domen/3`/`var(dəyərlər)` evristikası birbaşa işlədilir,
    tamamilə əvvəlki davranış.

    **3D/tam anizotropluq (A4)** — `points`/`targets` (n,2) [X,Y] və ya (n,3)
    [X,Y,Z] ola bilər. `range_v` (şaquli radius) verilməyibsə `range_`-ə
    bərabər qəbul edilir. `azimuth_deg`/`range_minor` verilməyəndə (defolt)
    üfüqi müstəvi İZOTROPDUR — yalnız Z miqyaslanır (əvvəlki M2 davranışı).
    `azimuth_deg`/`dip_deg` verilsə (və ya `auto_detect_anisotropy=True`
    etibarlı nəticə tapsa) tam geometrik anizotropluq transformu işə düşür
    (`geology/anisotropy.AnisotropyParams`). Transform HƏM Kriging
    matrisinə, HƏM DƏ qonşuluq axtarışına eyni obyektlə tətbiq olunur.

    **Yerli axtarış (moving neighbourhood, A1.2/A2)** — `search_radius`,
    `max_neighbors`, `sectors` və ya açıq `neighborhood` veriləndə hər
    hədəf üçün YALNIZ seçilmiş qonşularla AYRICA yerli sistem qurulur.
    `min_neighbors`-dan az uyğun qonşu olan hədəf üçün NaN qaytarılır —
    dəyər UYDURULMUR.

    **Avtomatik yerliləşmə (Gate 3)** — heç bir yerli parametr
    verilməyəndə də nöqtə sayı `auto_local_threshold`-u KEÇƏRSƏ sistem
    AVTOMATİK yerli rejimə keçir (`auto_local_max_neighbors` qonşu ilə),
    çünki minlərlə nöqtə üçün qlobal sıx sistem nə praktik, nə də
    geostatistik olaraq düzgündür. Bu həddən AŞAĞIDA (tipik quyu
    çoxluğu) nəticə köhnə QLOBAL kriging ilə BİRƏBİR eynidir — bu yol
    reqressiya testləri ilə qorunur.
    """
    range_: Optional[float] = None
    range_v: Optional[float] = None
    sill: Optional[float] = None
    nugget: float = 0.0
    model: str = MODEL_SPHERICAL
    auto_fit: bool = False
    auto_fit_nugget: bool = False
    azimuth_deg: Optional[float] = None
    range_minor: Optional[float] = None
    dip_deg: float = 0.0
    auto_detect_anisotropy: bool = False
    search_radius: Optional[float] = None
    min_neighbors: int = 1
    max_neighbors: Optional[int] = None
    #: istiqamətli balanslaşdırma (A2.5): 0=söndürülmüş, 4=kvadrant, …
    sectors: int = 0
    vertical_sectors: bool = False
    max_per_sector: Optional[int] = None
    #: XAM |ΔZ| həddi (A2.6) — geoloji olaraq uzaq layı tamamilə kəsir
    max_vertical_distance: Optional[float] = None
    #: seyrək məlumat ehtiyat yolları (A2.7)
    max_radius_expansions: int = 0
    radius_expansion_factor: float = 2.0
    max_search_radius: Optional[float] = None
    allow_knn_fallback: bool = False
    allow_global_fallback: bool = False
    #: tam `NeighborhoodConfig` — verilsə yuxarıdakı qısayolları ƏVƏZ EDİR
    neighborhood: Optional[NeighborhoodConfig] = None
    #: Gate 3: bu qədər nöqtədən sonra qlobal sistem AVTOMATİK yerliləşir
    auto_local_threshold: int = 100
    auto_local_max_neighbors: int = 40
    #: sərt data siyasəti (A1.4): "auto" = nugget==0 olanda dəqiq honor
    honor_hard_data: str = "auto"
    #: qeyri-sonlu (NaN/±inf) giriş nöqtələri çıxarılsınmı (A1.3)
    drop_non_finite: bool = True
    #: QAUSS variogramı üçün minimum nugget (sillə nisbətdə) — Kriging
    #: sistemini sabitləşdirir, bax `variogram.stabilizing_nugget()`.
    #: `0.0` bu davranışı söndürür.
    gaussian_min_nugget_ratio: float = GAUSSIAN_MIN_NUGGET_RATIO
    name: str = "Kriging (adi)"

    supports_z = True

    #: sonuncu `interpolate()` çağırışının introspeksiya/hesabat nəticələri
    #: (istəyə görə UI/PDF hesabatda göstərilə bilər) — bax `_parameters`.
    last_fit_: Optional[VariogramParameters] = field(default=None, repr=False, compare=False)
    last_anisotropy_: Optional[AnisotropyDetectionResult] = field(
        default=None, repr=False, compare=False)
    last_warnings_: list = field(default_factory=list, repr=False, compare=False)

    def __post_init__(self):
        if self.model == "auto" and not self.auto_fit:
            raise ValueError("model='auto' üçün auto_fit=True lazımdır.")
        if self.model != "auto" and self.model not in MODEL_FUNCS:
            raise ValueError(
                f"Naməlum variogram modeli: {self.model!r}. "
                f"Dəstəklənən: {tuple(MODEL_FUNCS) + ('auto',)}")
        if self.honor_hard_data not in ("auto", "always", "never"):
            raise ValueError(
                f"honor_hard_data 'auto'/'always'/'never' olmalıdır, "
                f"alındı: {self.honor_hard_data!r}")
        if self.min_neighbors < 1:
            raise ValueError(f"min_neighbors ≥ 1 olmalıdır, alındı: {self.min_neighbors}")
        if self.auto_local_threshold < 1:
            raise ValueError(
                f"auto_local_threshold ≥ 1 olmalıdır, alındı: {self.auto_local_threshold}")

    # ── variogram qiymətləndirməsi ────────────────────────────────────
    def _variogram(self, h: np.ndarray, range_: float, sill: float, nugget: float,
                   model: str, zero_at_origin: bool = True) -> np.ndarray:
        """`zero_at_origin` yalnız məlumat-məlumat matrisi üçün doğrudur:
        nöqtənin özü ilə arasındakı variogram sıfırdır. Sağ tərəfdə
        (məlumat-hədəf) isə nugget saxlanılır — məhz buna görə nugget
        sıfırdan böyük olanda kriging DƏQİQ interpolyator olmaqdan
        çıxır və ölçmə səhvini süzür.
        """
        func = MODEL_FUNCS[model]
        gamma = func(h, nugget, sill, range_)
        if zero_at_origin:
            return np.where(h <= 1e-12, 0.0, gamma)
        return gamma

    # ── giriş doğrulaması / normallaşdırma ────────────────────────────
    @staticmethod
    def _as_points(arr) -> np.ndarray:
        """(n,2) -> Z=0 sütunu ilə (n,3); (n,3) olduğu kimi qalır."""
        arr = np.asarray(arr, float)
        if arr.ndim != 2:
            raise ValueError("Nöqtələr 2D massiv olmalıdır: (n,2) və ya (n,3)")
        if arr.shape[1] == 2:
            return np.column_stack([arr, np.zeros(arr.shape[0])])
        if arr.shape[1] == 3:
            return arr
        raise ValueError(f"Nöqtələr (n,2) və ya (n,3) olmalıdır, alındı: {arr.shape}")

    def _prepare(self, points, values, targets):
        """XAM giriş → `(points3, values, targets3, target_ok, warnings)`.

        * uzunluq uyğunluğu YOXLANILIR (səssiz kəsmə YOX);
        * qeyri-sonlu (NaN/±inf) KOORDİNAT və ya DƏYƏR daşıyan sərt data
          nöqtələri `drop_non_finite=True` (defolt) olanda ÇIXARILIR və
          xəbərdarlıq yazılır — belə sətir Kriging matrisini bütövlüklə
          NaN edərdi (A1.3);
        * qeyri-sonlu HƏDƏF koordinatı üçün `target_ok=False` — nəticə
          NaN olur, sistem qurulmur;
        * ziddiyyətli dublikatlar `_dedupe_conflicting_points` ilə
          deterministik ORTALANIR (A1.4 dublikat siyasəti).
        """
        points3 = self._as_points(points)
        values = np.asarray(values, float).ravel()
        targets3 = self._as_points(targets)
        if values.shape[0] != points3.shape[0]:
            raise ValueError(
                f"points ({points3.shape[0]}) və values ({values.shape[0]}) "
                "uzunluğu uyğun gəlmir.")

        warnings: List[str] = []
        finite = np.all(np.isfinite(points3), axis=1) & np.isfinite(values)
        if not np.all(finite):
            if not self.drop_non_finite:
                raise ValueError(
                    f"{int(np.sum(~finite))} sərt data nöqtəsində NaN/sonsuz "
                    "koordinat və ya dəyər var (drop_non_finite=False).")
            warnings.append(
                f"{int(np.sum(~finite))} sərt data nöqtəsi NaN/sonsuz koordinat və ya "
                "dəyər daşıdığı üçün ÇIXARILDI — Kriging sistemi onlarla qurula bilməz.")
            points3, values = points3[finite], values[finite]

        target_ok = np.all(np.isfinite(targets3), axis=1)
        if not np.all(target_ok):
            warnings.append(
                f"{int(np.sum(~target_ok))} hədəf koordinatı NaN/sonsuzdur — "
                "onlar üçün NaN qaytarılır (dəyər uydurulmur).")

        points3, values, dup_warnings = self._dedupe_conflicting_points(points3, values)
        return points3, values, targets3, target_ok, warnings + dup_warnings

    def _dedupe_conflicting_points(self, points3: np.ndarray, values: np.ndarray):
        """Tam üst-üstə düşən (məsafə < 1e-9) giriş nöqtələrini araşdırır.

        Əvvəlki sükut davranış: `_solve_global`-da tam-üst-üstə-düşən
        hədəflər üçün dəyər `result[rows] = values[columns]` fancy-
        indexing ilə bərpa edilirdi — bir neçə eyni-koordinatlı GİRİŞ
        nöqtəsi ziddiyyətli dəyərlə verildikdə, hansının "qazanacağı"
        NumPy-ın sənədləşdirilməmiş sırasına görə həll olunurdu.

        İndi: ziddiyyətli dublikatlar DETERMİNİSTİK olaraq ORTALANIR və
        `last_warnings_`-ə yazılır; dəyərləri praktik eyni olan
        dublikatlar səssizcə (xəbərdarlıqsız) birləşdirilir. Bu, A1.4-ün
        tələb etdiyi AÇIQ dublikat siyasətidir: sərt datanın DƏQİQ
        honor edilməsi bu ortalanmış (təkil) dəyərə görədir.

        Bu funksiya HƏR `interpolate()`/`interpolate_with_variance()`
        çağırışında (o cümlədən SGS/SIS-in artan-nöqtə simulyasiya
        dövründə hər addımda, bax `sgs.py`/`facies.py`) işə düşür, ona
        görə İKİ PİLLƏLİDİR: əvvəlcə ucuz `O(n log n)` lexsort-əsaslı
        "heç bir dublikat yoxdurmu?" yoxlaması (adi, dublikatsız hal —
        demək olar bütün çağırışlar) sürətlə keçir; yalnız faktiki
        dublikat tapılan NADIR halda daha ətraflı (yenə `O(n log n)`,
        `np.unique` hash-əsaslı) qruplaşdırma işə düşür. Tam cüt-cüt
        (n,n) məsafə matrisi HEÇ VAXT hesablanmır — bu, əvvəlki versiyada
        `test_facies_performance.py`-ni O(n³)-ə qədər yavaşladan real
        performans reqressiyası idi (performans auditində tapılıb və
        canlı təsdiqlənib), buna görə iki pilləli struktur seçildi.
        """
        n = points3.shape[0]
        if n < 2:
            return points3, values, []

        order = np.lexsort((points3[:, 2], points3[:, 1], points3[:, 0]))
        adjacent_close = np.all(
            np.abs(np.diff(points3[order], axis=0)) < 1e-9, axis=1)
        if not np.any(adjacent_close):
            return points3, values, []

        rounded = np.round(np.ascontiguousarray(points3), 9)
        view = rounded.view([("", rounded.dtype)] * rounded.shape[1])
        _, first_idx, inverse, counts = np.unique(
            view, return_index=True, return_inverse=True, return_counts=True)
        inverse = inverse.reshape(-1)

        keep = np.zeros(n, dtype=bool)
        keep[first_idx] = True
        resolved = values.copy()
        warnings: list = []
        for group_id in np.where(counts > 1)[0]:
            idx = np.where(inverse == group_id)[0]
            group_values = values[idx]
            if not np.allclose(group_values, group_values[0], atol=1e-9, rtol=0.0):
                xy = points3[idx[0], :2]
                warnings.append(
                    f"Eyni koordinatda ({xy[0]:.4g}, {xy[1]:.4g}) {idx.size} "
                    f"ziddiyyətli sərt-data nöqtəsi tapıldı (dəyərlər: "
                    f"{group_values.tolist()}) — deterministik ORTA dəyər "
                    "istifadə olundu (əvvəlki sükut 'son yazı qazanır' "
                    "davranışı əvəzinə).")
            resolved[first_idx[group_id]] = float(np.mean(group_values))

        return points3[keep], resolved[keep], warnings

    # ── variogram/anizotropluq parametrləri ───────────────────────────
    def _parameters(self, points_xy: np.ndarray, values: np.ndarray):
        """`points_xy` — YALNIZ X,Y (Z auto-range təxminini korlamasın).

        Qaytarır: `(range_h, range_v, range_minor, azimuth_deg, sill,
        nugget, model, model_nugget)` — hamısı bu çağırış üçün İSTİFADƏ
        OLUNACAQ faktiki dəyərlər (fit/aşkarlanma nəticələri daxil).
        Nəticə HƏMİŞƏ `validate_variogram_parameters()`-dən keçir (A3.7)
        — etibarsız model solver-ə ÇATA BİLMİR.

        `model_nugget` — SABİTLƏŞDİRMƏDƏN ƏVVƏLKİ nugget, yəni MODELİN
        ÖZÜNÜN nugget effekti. Fərq VACİBDİR: `stabilizing_nugget()`
        qauss modelinə ƏDƏDİ requlyarlaşdırma kimi kiçik nugget əlavə edir,
        bu isə "ölçmədə səhv var" DEMƏK DEYİL. `honor_hard_data="auto"`
        qərarı məhz `model_nugget`-ə baxır — əks halda sırf ədədi düzəliş
        sərt datanın DƏQİQ honor edilməsini SÜKUTLA söndürərdi (tutulmuş
        həqiqi səhv, bax `tests/test_property_strategies.py::
        test_porosity_honours_hard_data_exactly`).
        """
        warnings: list = []
        fit: Optional[VariogramParameters] = None
        model = self.model if self.model != "auto" else MODEL_SPHERICAL

        if self.range_ is not None:
            range_h = float(self.range_)
        elif self.auto_fit:
            try:
                fit = fit_variogram_from_data(points_xy, values, model=self.model)
                range_h = fit.range_
                model = fit.model
            except ValueError as exc:
                warnings.append(
                    f"Avtomatik variogram fit alınmadı ({exc}) — ehtiyat evristika "
                    "(domen/3) işlədildi.")
                range_h = max(_max_pairwise_distance(points_xy) / 3.0, 1e-6)
        else:
            range_h = max(_max_pairwise_distance(points_xy) / 3.0, 1e-6)

        range_v = float(self.range_v) if self.range_v is not None else range_h
        sill = (max(fit.sill, 1e-12) if (fit is not None and self.sill is None)
               else max(float(self.sill) if self.sill is not None else float(np.var(values)),
                        1e-12))
        nugget = fit.nugget if (fit is not None and self.auto_fit_nugget) else self.nugget
        range_minor = float(self.range_minor) if self.range_minor is not None else range_h
        azimuth_deg = float(self.azimuth_deg) if self.azimuth_deg is not None else 0.0

        anisotropy: Optional[AnisotropyDetectionResult] = None
        if self.auto_detect_anisotropy and self.azimuth_deg is None:
            try:
                anisotropy = detect_anisotropy(points_xy, values)
                if anisotropy.reliable:
                    azimuth_deg = anisotropy.azimuth_deg
                    if self.range_ is None:
                        range_h = anisotropy.range_major
                    if self.range_minor is None:
                        range_minor = anisotropy.range_minor
                else:
                    warnings.extend(anisotropy.warnings)
            except ValueError as exc:
                warnings.append(f"Anizotropluq aşkarlanması alınmadı: {exc}")

        # Qauss modelinin ədədi sabitliyi (ölçülmüş qüsur — bax
        # `variogram.stabilizing_nugget` docstring-i). `model_nugget`
        # SABİTLƏŞDİRMƏDƏN ƏVVƏLKİ dəyəri saxlayır (bax docstring).
        model_nugget = float(nugget)
        nugget, nugget_warning = stabilizing_nugget(
            model, nugget, sill, self.gaussian_min_nugget_ratio)
        if nugget_warning:
            warnings.append(nugget_warning)

        # A3.7 — etibarsız parametr solver-ə ÇATMIR
        warnings.extend(validate_variogram_parameters(model, nugget, sill, range_h))

        self.last_fit_ = fit
        self.last_anisotropy_ = anisotropy
        self.last_warnings_ = warnings
        return (range_h, range_v, range_minor, azimuth_deg, sill, nugget, model,
                model_nugget)

    def _build_anisotropy(self, range_h, range_v, range_minor, azimuth_deg
                          ) -> Tuple[AnisotropyParams, List[str]]:
        params = AnisotropyParams(azimuth_deg=azimuth_deg, range_major=range_h,
                                  range_minor=range_minor, range_vertical=range_v,
                                  dip_deg=float(self.dip_deg))
        return params, params.validate()

    def _neighborhood_config(self, n_points: int, range_h: float
                             ) -> Tuple[Optional[NeighborhoodConfig], List[str]]:
        """Yerli rejimin konfiqurasiyası — `None` = QLOBAL sistem.

        Qlobal yol YALNIZ heç bir yerli parametr verilməyəndə VƏ nöqtə
        sayı `auto_local_threshold`-dan çox olmayanda seçilir (Gate 3).
        """
        if self.neighborhood is not None:
            config = self.neighborhood
            config.validate()
            return config, []

        explicit = (self.search_radius is not None or self.max_neighbors is not None
                    or self.sectors or self.max_vertical_distance is not None
                    or self.allow_knn_fallback or self.allow_global_fallback
                    or self.max_radius_expansions)
        warnings: List[str] = []
        max_neighbors = self.max_neighbors
        if not explicit:
            if n_points <= self.auto_local_threshold:
                return None, []
            max_neighbors = self.auto_local_max_neighbors
            warnings.append(
                f"Sərt data nöqtəsi sayı ({n_points}) `auto_local_threshold`"
                f"={self.auto_local_threshold} həddini keçdi — qlobal sıx sistem "
                f"əvəzinə AVTOMATİK yerli kriging işlədildi "
                f"(max_neighbors={max_neighbors}).")

        config = NeighborhoodConfig(
            min_neighbors=self.min_neighbors,
            max_neighbors=max_neighbors,
            search_radius=self.search_radius,
            max_search_radius=self.max_search_radius,
            radius_expansion_factor=self.radius_expansion_factor,
            max_radius_expansions=self.max_radius_expansions,
            sectors=self.sectors,
            vertical_sectors=self.vertical_sectors,
            max_per_sector=self.max_per_sector,
            max_vertical_distance=self.max_vertical_distance,
            allow_knn_fallback=self.allow_knn_fallback,
            allow_global_fallback=self.allow_global_fallback,
            support_range=range_h)
        config.validate()
        return config, warnings

    # ── ictimai giriş nöqtələri ───────────────────────────────────────
    def interpolate(self, points, values, targets) -> np.ndarray:
        """Yalnız qiymətlər `(m,)`. Diaqnostika üçün `krige()` işlədin."""
        return self._run(points, values, targets, want_variance=False,
                         want_diagnostics=False)[0]

    def interpolate_with_variance(self, points, values, targets):
        """`interpolate()`-in EYNİSİ, ƏLAVƏ olaraq kriging VARİANSINI da
        qaytarır (Phase 5/SGS: şərti Gauss paylanması `N(estimate,
        variance)` üçün tələb olunur).

        Qaytarır: `(estimate, variance)`, hər ikisi `(m,)` massiv.
        """
        estimate, variance, _ = self._run(points, values, targets, want_variance=True,
                                          want_diagnostics=False)
        return estimate, variance

    def krige(self, points, values, targets) -> KrigingResult:
        """TAM nəticə obyekti (A1.5): qiymət, varians, qonşu sayı, ən
        yaxın məsafə, ekstrapolyasiya bayrağı, dəstək təsnifatı, hər
        hədəf üçün qonşuluq və solver statusu, işlədilmiş variogram/
        anizotropluq parametrləri və bütün xəbərdarlıqlar."""
        _, _, result = self._run(points, values, targets, want_variance=True,
                                 want_diagnostics=True)
        return result

    # ── boru xəttinin özəyi ───────────────────────────────────────────
    def _run(self, points, values, targets, want_variance: bool,
             want_diagnostics: bool):
        """A5 boru xətti — bax modul docstring-i.

        Qaytarır: `(estimate, variance|None, KrigingResult|None)`.
        `want_diagnostics=False` olanda dəstək təsnifatı/status massivləri
        QURULMUR — SGS/SIS hər hüceyrədə bu metodu çağırır, ona görə
        diaqnostika ƏLAVƏ xərc kimi yalnız tələb olunanda ödənilir.
        """
        points3, values, targets3, target_ok, warnings = self._prepare(
            points, values, targets)
        m = targets3.shape[0]
        n = points3.shape[0]

        if n == 0:
            self.last_warnings_ = warnings + ["Heç bir etibarlı sərt data nöqtəsi qalmadı."]
            estimate = np.full(m, np.nan)
            variance = np.full(m, np.nan)
            if not want_diagnostics:
                return estimate, variance, None
            return estimate, variance, self._empty_diagnostics(
                m, estimate, variance, self.last_warnings_)

        (range_h, range_v, range_minor, azimuth_deg, sill, nugget, model,
         model_nugget) = self._parameters(points3[:, :2], values)
        warnings = warnings + self.last_warnings_
        anisotropy, aniso_warnings = self._build_anisotropy(
            range_h, range_v, range_minor, azimuth_deg)
        warnings += aniso_warnings

        if n == 1 and not want_variance and not want_diagnostics:
            # tək nöqtə + yalnız qiymət: sistem qurulmadan da nəticə
            # eynidir (çəki = 1) — ƏVVƏLKİ qısayol, ədəd-ədəd saxlanılır.
            self.last_warnings_ = warnings
            return np.where(target_ok, values[0], np.nan), None, None

        scaled_points = anisotropy.transform(points3)
        scaled_targets = anisotropy.transform(targets3)

        config, config_warnings = self._neighborhood_config(n, range_h)
        warnings += config_warnings
        self.last_warnings_ = warnings

        if config is None:
            estimate, variance, solver, lagrange = self._solve_global_path(
                scaled_points, values, scaled_targets, target_ok, range_h, sill,
                nugget, model, model_nugget)
            neighbor_count = np.where(target_ok, n, 0)
            batch: Optional[BatchNeighborhood] = None
            nearest_index, nearest_distance = self._nearest_hard_data(
                scaled_points, scaled_targets, model_nugget)
        else:
            (estimate, variance, solver, neighbor_count, batch,
             nearest_index, nearest_distance, lagrange) = self._solve_local_path(
                scaled_points, values, scaled_targets, target_ok, range_h, sill,
                nugget, model, config, points3, targets3)

        estimate, variance, solver = self._honor_hard_data(
            values, nearest_index, nearest_distance, target_ok, model_nugget,
            estimate, variance, solver)

        if not want_diagnostics:
            return estimate, variance, None

        result = self._build_diagnostics(
            estimate, variance, solver, lagrange, neighbor_count, batch,
            scaled_points, scaled_targets, target_ok, anisotropy, model, range_h,
            sill, nugget, config is not None, warnings)
        return estimate, variance, result

    # ── qlobal (kiçik məlumat çoxluğu) yol ────────────────────────────
    def _solve_global_path(self, points, values, targets, target_ok, range_h, sill,
                           nugget, model, model_nugget=None):
        """Bütün nöqtələrlə TƏK sistem — kiçik quyu çoxluğu üçün.

        Sistem BİR DƏFƏ qurulub `m` sağ tərəflə həll olunur (`O(n³ + n²m)`),
        yəni hər hüceyrə üçün YENİDƏN həll YOXDUR. `auto_local_threshold`-
        dan böyük çoxluqda bu yol ÇAĞIRILMIR (bax `_neighborhood_config`).
        """
        estimate, variance, lagrange = self._solve_global(
            points, values, targets, range_h, sill, nugget, model,
            return_variance=True, return_lagrange=True, model_nugget=model_nugget)
        solver = np.full(targets.shape[0], SOLVER_DIRECT, dtype=object)
        if points.shape[0] == 1:
            solver[:] = SOLVER_SINGLE_VALUE
        if not np.all(target_ok):
            estimate = np.where(target_ok, estimate, np.nan)
            variance = np.where(target_ok, variance, np.nan)
            lagrange = np.where(target_ok, lagrange, np.nan)
            solver[~target_ok] = SOLVER_NONE
        return estimate, variance, solver, lagrange

    def _solve_global(self, points, values, targets, range_h, sill, nugget, model,
                      return_variance: bool = False, return_lagrange: bool = False,
                      model_nugget: Optional[float] = None):
        """Bütün nöqtələrlə TƏK kriging sistemi.

        `return_variance=True` — kriging variansı `σ²(x0) = Σ w_i·
        γ(x_i,x0) + μ` ARTIQ HƏLL EDİLMİŞ sistemdən (`weights`, `right`,
        Laqranj vuruğu `solution[n,:]`) hesablanır — YENİDƏN HƏLL YOXDUR,
        sadəcə eyni nəticədən ƏLAVƏ oxuma (A1.6 vahid konvensiya).
        """
        n = points.shape[0]
        # sol tərəf: variogram matrisi + Laqranj sətri/sütunu
        left = np.ones((n + 1, n + 1))
        left[:n, :n] = self._variogram(_distance_matrix(points, points), range_h,
                                       sill, nugget, model, zero_at_origin=True)
        left[n, n] = 0.0

        right = np.ones((n + 1, targets.shape[0]))
        effective_nugget = nugget if model_nugget is None else model_nugget
        right[:n, :] = self._variogram(
            _distance_matrix(points, targets), range_h, sill, nugget, model,
            zero_at_origin=effective_nugget <= 0.0)

        try:
            solution = np.linalg.solve(left, right)
        except np.linalg.LinAlgError:
            solution = np.linalg.lstsq(left, right, rcond=None)[0]

        weights = solution[:n, :]
        result = (weights * values[:, None]).sum(axis=0)

        # MODELİN nugget-i sıfırdırsa kriging dəqiq interpolyatordur —
        # nöqtələri bərpa edirik. `model_nugget` verilməyəndə `nugget`-in
        # özü işlədilir (Phase A davranışı, `_solve_local` üçün).
        exact_rows = None
        if (nugget if model_nugget is None else model_nugget) <= 0.0:
            hits = _distance_matrix(targets, points) < 1e-9
            rows, columns = np.where(hits)
            result[rows] = values[columns]
            exact_rows = rows

        if not return_variance:
            return result

        lagrange = solution[n, :].copy()
        variance = np.sum(weights * right[:n, :], axis=0) + solution[n, :]
        variance = np.clip(variance, 0.0, None)
        if exact_rows is not None:
            variance[exact_rows] = 0.0
        if return_lagrange:
            return result, variance, lagrange
        return result, variance

    # ── yerli yol (istehsal yolu) ─────────────────────────────────────
    def _solve_local_path(self, points, values, targets, target_ok, range_h, sill,
                          nugget, model, config: NeighborhoodConfig,
                          raw_points, raw_targets):
        """Hər hədəf üçün AYRICA yerli sistem (A1.2).

        Qonşuluq `spatial_search.NeighborhoodSelector` ilə seçilir — EYNİ
        anizotrop transformasiya fəzasında (A4.3). Seçim TOPLU aparılır
        (`select_batch`: sadə k-ən-yaxın halında tək cKDTree çağırışı),
        sonra eyni ölçülü sistemlər qruplaşdırılıb TOPLU həll olunur
        (`_solve_batch`) — hədəf başına nə Python axtarış dövrəsi, nə də
        ayrıca LAPACK çağırışı var (A7).

        `NeighborhoodSelector` artıq TRANSFORMASİYA EDİLMİŞ nöqtələr
        üzərində qurulur (`anisotropy=None`), çünki `points`/`targets`
        bu metoda daxil olanda transformdan KEÇMİŞDİR — ikiqat
        transformasiya riyazi səhv olardı.
        """
        m = targets.shape[0]
        estimate = np.full(m, np.nan)
        variance = np.full(m, np.nan)
        lagrange = np.full(m, np.nan)
        solver = np.full(m, SOLVER_NONE, dtype=object)

        # A2.6 — XAM |ΔZ| kəsiyi XAM Z-yə görə ölçülməlidir, ona görə
        # selektora xam Z sütunu ötürülür (transform edilmiş Z DEYİL).
        selector = NeighborhoodSelector(points, anisotropy=None, config=config,
                                        index="kdtree" if self._use_tree(m, points)
                                        else "brute")
        selector.set_raw_vertical(raw_points[:, 2])

        if config.max_vertical_distance is not None:
            # şaquli kəsik hədəfin XAM Z-sini tələb edir → sətir-sətir yol
            batch = selector._select_batch_by_rows_with_depth(targets, raw_targets[:, 2])
        else:
            batch = selector.select_batch(targets)

        counts = batch.counts.copy()
        counts[~target_ok] = 0
        neighbor_count = counts

        nearest_index = np.where(counts > 0, batch.indices[:, 0], -1)
        nearest_distance = np.where(counts > 0, batch.distances[:, 0], np.inf)

        total_sill = float(nugget) + float(sill)
        for size in np.unique(counts[counts > 0]):
            rows_array = np.where(counts == size)[0]
            index_matrix = batch.indices[rows_array, :size]
            chunk = max(1, _BATCH_ELEMENTS // max(int(size) * int(size), 1))
            for start in range(0, rows_array.size, chunk):
                block = rows_array[start:start + chunk]
                idx = index_matrix[start:start + chunk]
                est, var, status, mu = self._solve_group(
                    points, values, targets[block], idx, range_h, sill, nugget,
                    model, total_sill)
                estimate[block] = est
                variance[block] = var
                solver[block] = status
                lagrange[block] = mu

        return (estimate, variance, solver, neighbor_count, batch,
                nearest_index, nearest_distance, lagrange)

    @staticmethod
    def _use_tree(n_targets: int, points: np.ndarray) -> bool:
        """cKDTree qurmaq sərfəlidirmi.

        Ağac qurmaq `O(n log n)`, hər sorğu `O(log n)`; kobud güc isə
        qurmasızdır, amma hər sorğu `O(n)`. Hədəf sayı `log₂ n`-dən
        azdırsa qurma xərci ödənmir — məhz SGS/SIS-in hər hüceyrədə TƏK
        hədəflə çağırdığı hal (bax `sgs.py`/`facies.py`), orada ağac
        qurmaq REAL yavaşlama olardı.
        """
        return n_targets >= max(4, math.log2(max(points.shape[0], 2)))

    def _solve_group(self, points, values, targets, index_matrix, range_h, sill,
                     nugget, model, total_sill):
        """Eyni qonşu SAYINA malik hədəflər dəstəsi üçün toplu həll.

        `index_matrix` (g, k) — hər sətir bir hədəfin qonşu indeksləri.
        """
        g, k = index_matrix.shape
        neighbors = points[index_matrix]                        # (g,k,3)
        diff = neighbors[:, :, None, :] - neighbors[:, None, :, :]
        pair_distance = np.sqrt(np.sum(diff * diff, axis=-1))   # (g,k,k)
        target_distance = np.sqrt(
            np.sum((neighbors - targets[:, None, :]) ** 2, axis=-1))   # (g,k)

        left = np.ones((g, k + 1, k + 1))
        left[:, :k, :k] = self._variogram(pair_distance, range_h, sill, nugget,
                                          model, zero_at_origin=True)
        left[:, k, k] = 0.0
        right = np.ones((g, k + 1))
        right[:, :k] = self._variogram(target_distance, range_h, sill, nugget,
                                       model, zero_at_origin=nugget <= 0.0)

        solution, status = _solve_batch(left, right, target_distance, total_sill)
        weights = solution[:, :k]
        estimate = np.sum(weights * values[index_matrix], axis=1)
        lagrange = solution[:, k]
        variance = np.sum(weights * right[:, :k], axis=1) + lagrange
        return estimate, np.clip(variance, 0.0, None), status, lagrange

    # ── sərt datanın DƏQİQ honor edilməsi (A1.4) ──────────────────────
    def _nearest_hard_data(self, points, targets, nugget):
        """Qlobal yolda hər hədəf üçün ən yaxın sərt data (indeks, məsafə).

        Yalnız honor siyasəti FAKTİKİ olaraq lazım olanda hesablanır.
        `nugget == 0` halında `_solve_global` dəqiq honor-u ARTIQ edib,
        amma statusu `SOLVER_EXACT` kimi işarələmək üçün indekslər yenə
        lazımdır — hesablama `(m,n)` matrisdir, sistemin özünün onsuz da
        qurduğu sağ tərəflə eyni ölçüdə, yəni ƏLAVƏ asimptotik xərc yox.
        Honor söndürülübsə heç nə hesablanmır.
        """
        m = targets.shape[0]
        if self.honor_hard_data == "never" or points.shape[0] == 0 or m == 0:
            return np.full(m, -1, dtype=int), np.full(m, np.inf)
        if self.honor_hard_data == "auto" and nugget > 0.0:
            return np.full(m, -1, dtype=int), np.full(m, np.inf)
        distances = _distance_matrix(targets, points)
        index = np.argmin(distances, axis=1)
        return index.astype(int), distances[np.arange(m), index]

    def _honor_hard_data(self, values, nearest_index, nearest_distance, target_ok,
                         nugget, estimate, variance, solver):
        """Hədəf sərt data nöqtəsi ilə ÜST-ÜSTƏ düşəndə dəyəri DƏQİQ qaytarır.

        SİYASƏT (`honor_hard_data`):
            "auto" (defolt) — YALNIZ `nugget == 0` olanda. Nugget > 0
                Kriging-in ölçmə səhvini SÜZMƏSİ deməkdir; belə halda
                ölçülmüş dəyərə "geri qaytarmaq" modelin öz fərziyyəsini
                pozardı, ona görə honor EDİLMİR (əvvəlki davranış).
            "always" — nugget-dən ASILI OLMAYARAQ dəqiq honor.
            "never"  — heç vaxt; yalnız sistemin öz nəticəsi.

        Üst-üstə düşən dəyər `_dedupe_conflicting_points`-in verdiyi
        (ziddiyyət halında ORTALANMIŞ) təkil dəyərdir — AÇIQ dublikat
        siyasəti budur.

        Yerli yolda `nearest_index`/`nearest_distance` qonşuluq seçiminin
        ARTIQ hesabladığı nəticədir — ƏLAVƏ məsafə axtarışı YOXDUR (A7).
        """
        if self.honor_hard_data == "never":
            return estimate, variance, solver
        if self.honor_hard_data == "auto" and nugget > 0.0:
            return estimate, variance, solver

        hit = (nearest_index >= 0) & (nearest_distance < 1e-9) & target_ok
        if not np.any(hit):
            return estimate, variance, solver
        estimate = estimate.copy()
        variance = variance.copy()
        estimate[hit] = values[nearest_index[hit]]
        variance[hit] = 0.0
        solver = solver.copy()
        solver[hit] = SOLVER_EXACT
        return estimate, variance, solver

    # ── diaqnostika ───────────────────────────────────────────────────
    def _empty_diagnostics(self, m, estimate, variance, warnings) -> KrigingResult:
        return KrigingResult(
            estimate=estimate, variance=variance,
            neighbor_count=np.zeros(m, dtype=int),
            nearest_distance=np.full(m, np.inf),
            extrapolated=np.ones(m, dtype=bool),
            support=np.full(m, SUPPORT_EXTRAPOLATED, dtype=object),
            neighborhood_status=np.full(m, "empty", dtype=object),
            solver=np.full(m, SOLVER_NONE, dtype=object),
            lagrange=np.full(m, np.nan),
            anisotropy=AnisotropyParams(), model=self.model if self.model != "auto"
            else MODEL_SPHERICAL, range_=float("nan"), sill=float("nan"),
            nugget=float(self.nugget), local=False, warnings=list(warnings))

    def _build_diagnostics(self, estimate, variance, solver, lagrange, neighbor_count,
                           batch, points, targets, target_ok, anisotropy,
                           model, range_h, sill, nugget, local, warnings) -> KrigingResult:
        """Hər hədəf üçün dəstək təsnifatı + status massivləri.

        Qlobal yolda qonşuluq obyekti yoxdur, ona görə təsnifat üçün
        AYRICA (yalnız burada, `krige()` çağırılanda) `NeighborhoodSelector`
        qurulur — `interpolate()` bu xərci ÖDƏMİR."""
        if batch is None:
            config = NeighborhoodConfig(support_range=range_h)
            selector = NeighborhoodSelector(points, anisotropy=None, config=config)
            batch = selector.select_batch(targets)

        nearest = np.where(batch.counts > 0, batch.distances[:, 0], np.inf)
        support = np.asarray(batch.support, dtype=object).copy()
        status = np.asarray(batch.status, dtype=object).copy()
        if not np.all(target_ok):
            nearest[~target_ok] = np.inf
            support[~target_ok] = SUPPORT_EXTRAPOLATED
            status[~target_ok] = "empty"
        for message in batch.warnings:
            if message not in warnings:
                warnings.append(message)

        return KrigingResult(
            estimate=estimate, variance=variance,
            neighbor_count=np.asarray(neighbor_count, dtype=int),
            nearest_distance=nearest,
            extrapolated=(support == SUPPORT_EXTRAPOLATED),
            support=support, neighborhood_status=status, solver=solver,
            lagrange=np.asarray(lagrange, float),
            anisotropy=anisotropy, model=model, range_=float(range_h),
            sill=float(sill), nugget=float(nugget), local=local,
            fit=self.last_fit_, warnings=list(warnings))


INTERPOLATORS = {
    "Ən yaxın qonşu": NearestNeighbour,
    "Əks məsafə (IDW)": InverseDistance,
    "Kriging (adi)": OrdinaryKriging,
}


# ── xassə-yönlü çevirmələr (A6 genişlənmə nöqtəsi → Phase B) ──────────
# Bu siniflər Phase A-da BURADA yaradılmışdı; Phase B onları
# `geology/transforms.py`-yə köçürdü və genişləndirdi (logit, normal-score,
# geri-çevirmə modları, varians köçürməsi). İKİ implementasiya SAXLANMIR
# — burada YALNIZ geriyə-uyğun ad bağlantısı var, ona görə mövcud
# `from .interpolation import ValueTransform/LogTransform` idxalları
# işləməyə davam edir və EYNİ obyektləri qaytarır.
_LEGACY_TRANSFORM_EXPORTS = (ValueTransform, LogTransform, IDENTITY_TRANSFORM,
                             LOG_TRANSFORM)


def interpolate_property(interpolator: IPropertyInterpolator,
                         points: np.ndarray, values: np.ndarray,
                         targets: np.ndarray,
                         log_transform: bool = False,
                         minimum: Optional[float] = None,
                         maximum: Optional[float] = None,
                         transform: Optional[ValueTransform] = None) -> np.ndarray:
    """İnterpolyasiya + dəyər çevirməsi + hədlərin tətbiqi.

    Log çevirmə keçiricilik üçündür: ln(k) fəzasında interpolyasiya
    həm mənfi dəyərin qarşısını alır, həm də fiziki cəhətdən daha
    doğrudur, çünki keçiricilik log-normal paylanır.

    `transform` (A6) — açıq `ValueTransform`; verilməyəndə `log_transform`
    bayrağı `LOG_TRANSFORM`/`IDENTITY_TRANSFORM` seçir (mövcud
    çağırışlar DƏYİŞMİR). İkisi birlikdə verilə bilməz.
    """
    if transform is not None and log_transform:
        raise ValueError("`transform` və `log_transform=True` birlikdə verilə bilməz.")
    if transform is None:
        transform = LOG_TRANSFORM if log_transform else IDENTITY_TRANSFORM

    values = np.asarray(values, float)
    if isinstance(transform, ValueTransform) and type(transform) is ValueTransform:
        result = interpolator.interpolate(points, values, targets)
    else:
        result = transform.inverse(
            interpolator.interpolate(points, transform.forward(values), targets))

    if minimum is not None:
        result = np.maximum(result, minimum)
    if maximum is not None:
        result = np.minimum(result, maximum)
    return result
