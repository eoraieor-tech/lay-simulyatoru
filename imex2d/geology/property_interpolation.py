"""XASSƏ-YÖNLÜ İNTERPOLYASİYA + QEYRİ-MÜƏYYƏNLİK (B1 + B2).

Bu modul Phase A-nın ÜZƏRİNDƏ oturur və heç bir riyazi mühərriki
TƏKRARLAMIR:

    məsafə/anizotropluq   → `anisotropy.AnisotropyParams`
    qonşuluq              → `spatial_search.NeighborhoodSelector`
    variogram             → `variogram.py`
    kriginq + varians     → `interpolation.OrdinaryKriging.krige()`

Buraya ƏLAVƏ olunan yalnız budur:

    XASSƏNİN STATİSTİK TƏBİƏTİ (`property_config.PropertyStrategy`)
      → məlumat keyfiyyəti  (`data_quality.run_quality_control`)
      → çevirmə             (`transforms`)
      → Phase A kriginqi ÇEVRİLMİŞ fəzada
      → geri çevirmə + variansın köçürülməsi
      → fiziki hədlərin AÇIQ tətbiqi (kəsilən hüceyrə SAYILIR)
      → QEYRİ-MÜƏYYƏNLİK SEMANTİKASI + dəstək təsnifatı

Kateqorik xassələr (FACIES) TAMAM AYRI yoldan gedir: indikator kriginq
→ kateqoriya ehtimalları → ən ehtimallı kod + entropiya. Kateqorik kod
HEÇ VAXT kəsilməz kriginqdən keçmir (GATE B4).

QEYRİ-MÜƏYYƏNLİK SEMANTİKASI (B2.1) — nəyin nə olduğu QARIŞDIRILMIR:

    `transformed_variance`  çevrilmiş fəzada ƏSL kriginq variansı σ²_y
    `variance`              orijinal fəzaya köçürülmüş varians;
                            `variance_kind` onun DƏQİQ (loq-normal) yoxsa
                            DELTA (logit/normal-score) olduğunu deyir
    `nearest_distance`      anizotrop fəzada ən yaxın sərt data məsafəsi
    `data_density`          korrelyasiya radiusu daxilindəki nöqtə sayı —
                            DİAQNOSTİKADIR, ehtimal DEYİL
    `support`               Phase A həndəsi dəstək təsnifatı
    `confidence`            HIGH/MEDIUM/LOW/EXTRAPOLATED — İNTERPRETASİYA
                            kateqoriyasıdır, KALİBRLƏNMİŞ EHTİMAL DEYİL
                            (B2.3 — bu, docstring-də və sahə adında yazılır)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional

import numpy as np

from .anisotropy import AnisotropyParams
from .data_quality import DataQualityReport, QCConfig, run_quality_control
from .interpolation import InverseDistance, KrigingResult, NearestNeighbour, OrdinaryKriging
from .property_config import (BoundPolicy, InterpolationKind, PropertyStrategy,
                              UncertaintyKind, VariableType, resolve_strategy)
from .spatial_search import (SUPPORT_BOUNDARY, SUPPORT_EXTRAPOLATED, SUPPORT_WEAK,
                             SUPPORT_WELL)
from .transforms import TransformError, VarianceKind


class Confidence(str, Enum):
    """İNTERPRETASİYA kateqoriyaları (B2.3).

    QƏTİ QEYD: bunlar KALİBRLƏNMİŞ EHTİMAL DEYİL. "HIGH" 95% demək
    deyil — bu, "sıx məlumatla əhatələnmiş, kiçik nisbi kriginq variansı"
    deməkdir. Statistik ölçü lazım olanda `variance`/`std` və
    `cross_validation`-un standartlaşdırılmış xəta metrikləri işlədilir.
    """

    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    EXTRAPOLATED = "extrapolated"


#: `confidence` təsnifatının hədləri — nisbi kriginq variansı
#: `σ²/(nugget+sill)` üzrə. Sənədləşdirilmiş, sabit, deterministik.
CONFIDENCE_VARIANCE_HIGH = 0.30
CONFIDENCE_VARIANCE_MEDIUM = 0.70


@dataclass
class PropertyEstimate:
    """Kəsilməz xassə üçün TAM nəticə (B2).

    Bütün massivlər `(m,)` — hədəf sayı qədər.
    """

    property_name: str
    estimate: np.ndarray                #: FİZİKİ fəzada, hədlər tətbiq olunmuş
    variance: np.ndarray                #: fiziki fəzada (bax `variance_kind`)
    std: np.ndarray                     #: `sqrt(variance)`
    raw_estimate: np.ndarray            #: hədlər tətbiq olunmadan ƏVVƏLKİ dəyər
    transformed_estimate: np.ndarray    #: çevrilmiş fəzada kriginq qiyməti
    transformed_variance: np.ndarray    #: çevrilmiş fəzada ƏSL kriginq variansı
    variance_kind: VarianceKind
    neighbor_count: np.ndarray
    nearest_distance: np.ndarray
    data_density: np.ndarray            #: korrelyasiya radiusundakı nöqtə sayı
    support: np.ndarray                 #: Phase A `SUPPORT_*`
    extrapolated: np.ndarray
    confidence: np.ndarray              #: `Confidence` — bax sinif docstring-i
    solver: np.ndarray
    bound_adjusted: np.ndarray          #: hədlərə görə dəyişdirilən hüceyrələr
    strategy: PropertyStrategy
    kriging: Optional[KrigingResult] = None
    quality: Optional[DataQualityReport] = None
    warnings: List[str] = field(default_factory=list)

    def __array__(self, dtype=None, copy=None):
        """`np.asarray(result)` → qiymətlər (geriyə-uyğunluq körpüsü)."""
        if dtype is None:
            return (np.array(self.estimate, copy=True) if copy else self.estimate)
        return self.estimate.astype(dtype, copy=True)

    def __len__(self) -> int:
        return int(self.estimate.size)

    @property
    def uncertainty_kind(self) -> UncertaintyKind:
        return self.strategy.uncertainty_kind

    def as_grids(self) -> Dict[str, np.ndarray]:
        """Xəritə kimi göstərilə bilən massivlər (B2.4).

        UI/hesabat qatına heç bir yeni asılılıq gətirmədən qeyri-
        müəyyənlik xəritəsi qurmağa imkan verir."""
        return {
            "estimate": self.estimate,
            "variance": self.variance,
            "std": self.std,
            "transformed_variance": self.transformed_variance,
            "nearest_distance": self.nearest_distance,
            "neighbor_count": self.neighbor_count.astype(float),
            "data_density": self.data_density.astype(float),
            "extrapolated": self.extrapolated.astype(float),
            "confidence_rank": np.asarray(
                [_CONFIDENCE_RANK[str(c)] for c in self.confidence], dtype=float),
        }

    def summary(self) -> str:
        finite = np.isfinite(self.estimate)
        lines = [f"{self.property_name}: {self.strategy.describe()}",
                 f"  hədəf {self.estimate.size} · etibarlı {int(finite.sum())} · "
                 f"ekstrapolyasiya {int(np.sum(self.extrapolated))} · "
                 f"hədd düzəlişi {int(np.sum(self.bound_adjusted))}"]
        if finite.any():
            lines.append(
                f"  dəyər: min {np.nanmin(self.estimate):.5g} "
                f"orta {np.nanmean(self.estimate):.5g} "
                f"maks {np.nanmax(self.estimate):.5g}")
            lines.append(f"  varians ({self.variance_kind.value}): "
                         f"orta {np.nanmean(self.variance):.4g}")
        unique, counts = np.unique(self.confidence.astype(str), return_counts=True)
        lines.append("  etimad: " + ", ".join(f"{u}={c}" for u, c in zip(unique, counts)))
        lines.extend(f"  ⚠ {w}" for w in self.warnings)
        return "\n".join(lines)


_CONFIDENCE_RANK = {Confidence.HIGH.value: 3.0, Confidence.MEDIUM.value: 2.0,
                    Confidence.LOW.value: 1.0, Confidence.EXTRAPOLATED.value: 0.0}


@dataclass
class CategoricalEstimate:
    """Kateqorik xassə üçün nəticə (B1.6).

    `probabilities` (m, k) — SÜTUNLAR `categories` sırasındadır; hər
    sətir `[0,1]`-dədir və cəmi 1-dir (düzəliş edilibsə sayılır).
    """

    property_name: str
    categories: np.ndarray              #: (k,) kateqoriya kodları
    probabilities: np.ndarray           #: (m, k)
    most_probable: np.ndarray           #: (m,) kod
    entropy: np.ndarray                 #: (m,) Şennon entropiyası (nat)
    max_probability: np.ndarray         #: (m,)
    neighbor_count: np.ndarray
    nearest_distance: np.ndarray
    support: np.ndarray
    extrapolated: np.ndarray
    confidence: np.ndarray
    n_probability_corrections: int      #: [0,1]-ə kəsilib normallaşdırılan hüceyrə
    strategy: PropertyStrategy
    quality: Optional[DataQualityReport] = None
    warnings: List[str] = field(default_factory=list)

    def __len__(self) -> int:
        return int(self.most_probable.size)

    @property
    def normalized_entropy(self) -> np.ndarray:
        """`entropy / ln(k)` ∈ [0,1] — kateqoriya sayından asılı olmayan
        qeyri-müəyyənlik ölçüsü (1 = tam qeyri-müəyyən)."""
        k = self.categories.size
        if k < 2:
            return np.zeros_like(self.entropy)
        return self.entropy / np.log(k)

    def probability_of(self, category: int) -> np.ndarray:
        idx = int(np.where(self.categories == category)[0][0])
        return self.probabilities[:, idx]

    def as_grids(self) -> Dict[str, np.ndarray]:
        grids = {
            "most_probable": self.most_probable.astype(float),
            "max_probability": self.max_probability,
            "entropy": self.entropy,
            "normalized_entropy": self.normalized_entropy,
            "nearest_distance": self.nearest_distance,
            "neighbor_count": self.neighbor_count.astype(float),
            "extrapolated": self.extrapolated.astype(float),
        }
        for index, code in enumerate(self.categories):
            grids[f"probability_{int(code)}"] = self.probabilities[:, index]
        return grids

    def summary(self) -> str:
        lines = [f"{self.property_name}: kateqorik ({self.categories.size} kod)",
                 f"  hədəf {len(self)} · ehtimal düzəlişi "
                 f"{self.n_probability_corrections}"]
        unique, counts = np.unique(self.most_probable, return_counts=True)
        lines.append("  proqnoz: " + ", ".join(
            f"{int(u)}={c}" for u, c in zip(unique, counts)))
        lines.append(f"  orta normallaşdırılmış entropiya: "
                     f"{float(np.mean(self.normalized_entropy)):.3f}")
        lines.extend(f"  ⚠ {w}" for w in self.warnings)
        return "\n".join(lines)


# ── etimad təsnifatı ──────────────────────────────────────────────────
def classify_confidence(support: np.ndarray, variance: np.ndarray,
                        total_sill: float, neighbor_count: np.ndarray) -> np.ndarray:
    """Phase A həndəsi dəstəyini + NİSBİ kriginq variansını birləşdirir.

    Qaydalar (deterministik, sıra ilə):

        dəstək EXTRAPOLATED və ya qonşu yoxdur ... EXTRAPOLATED
        dəstək BOUNDARY/WEAK ..................... LOW
        σ²/(nugget+sill) > 0.70 .................. LOW
        σ²/(nugget+sill) > 0.30 .................. MEDIUM
        qalan .................................... HIGH

    Nisbi varians işlədilir (mütləq deyil), çünki `σ²`-nin miqyası
    xassədən xassəyə dəyişir; `nugget+sill`-ə bölmək onu ölçüsüz edir və
    "a-priori dispersiyanın neçə faizi izah olunmayıb" kimi oxunur."""
    support = np.asarray(support, dtype=object)
    variance = np.asarray(variance, float)
    neighbor_count = np.asarray(neighbor_count, int)
    scale = total_sill if (np.isfinite(total_sill) and total_sill > 0.0) else 1.0
    ratio = np.where(np.isfinite(variance), variance / scale, np.inf)

    result = np.full(support.shape, Confidence.HIGH.value, dtype=object)
    result[ratio > CONFIDENCE_VARIANCE_HIGH] = Confidence.MEDIUM.value
    result[ratio > CONFIDENCE_VARIANCE_MEDIUM] = Confidence.LOW.value
    weak = np.isin(support.astype(str), [SUPPORT_BOUNDARY, SUPPORT_WEAK])
    result[weak] = Confidence.LOW.value
    lost = (support.astype(str) == SUPPORT_EXTRAPOLATED) | (neighbor_count == 0)
    result[lost] = Confidence.EXTRAPOLATED.value
    return result


def compute_data_density(points: np.ndarray, targets: np.ndarray,
                         radius: float,
                         anisotropy: Optional[AnisotropyParams] = None) -> np.ndarray:
    """Hər hədəfin ətrafında `radius` daxilindəki sərt data sayı.

    DİAQNOSTİKADIR (B2.1 №4) — ehtimal və ya varians DEYİL. Kriginq
    variansından fərqli, MÜSTƏQİL məlumat verir: varians həndəsəyə,
    sıxlıq isə say/örtüyə baxır. cKDTree ilə `O((n+m) log n)`."""
    from scipy.spatial import cKDTree
    from .anisotropy import transform_points

    if points.shape[0] == 0 or targets.shape[0] == 0:
        return np.zeros(targets.shape[0], dtype=int)
    tree = cKDTree(transform_points(points, anisotropy))
    query = transform_points(targets, anisotropy)
    finite = np.all(np.isfinite(query), axis=1)
    counts = np.zeros(targets.shape[0], dtype=int)
    if np.any(finite):
        counts[finite] = np.asarray(
            tree.query_ball_point(query[finite], r=float(max(radius, 1e-12)),
                                  return_length=True), dtype=int)
    return counts


# ── kəsilməz xassə yolu ───────────────────────────────────────────────
def _build_interpolator(strategy: PropertyStrategy) -> object:
    if strategy.interpolation is InterpolationKind.NEAREST:
        return NearestNeighbour()
    if strategy.interpolation is InterpolationKind.IDW:
        return InverseDistance(search_radius=strategy.search_radius)
    return OrdinaryKriging(
        model=strategy.variogram_model if strategy.variogram_model != "auto" else "auto",
        auto_fit=strategy.variogram_model == "auto",
        max_neighbors=strategy.max_neighbors,
        min_neighbors=strategy.min_neighbors,
        search_radius=strategy.search_radius,
        sectors=strategy.sectors,
        honor_hard_data=strategy.honor_hard_data,
        azimuth_deg=(strategy.anisotropy.azimuth_deg if strategy.anisotropy else None),
        range_minor=(strategy.anisotropy.range_minor if strategy.anisotropy else None),
        range_=(strategy.anisotropy.range_major if strategy.anisotropy else None),
        range_v=(strategy.anisotropy.range_vertical if strategy.anisotropy else None),
        dip_deg=(strategy.anisotropy.dip_deg if strategy.anisotropy else 0.0))


def interpolate_property_field(points, values, targets,
                               strategy: Optional[PropertyStrategy] = None,
                               property_name: Optional[str] = None,
                               qc: Optional[QCConfig] = None,
                               run_qc: bool = True,
                               kriging_overrides: Optional[Dict[str, object]] = None
                               ) -> PropertyEstimate:
    """Bir KƏSİLMƏZ xassənin tam boru xətti (B6).

    QC → çevirmə → Phase A kriginqi → geri çevirmə → hədlər → etimad.

    `strategy` verilməyəndə `property_name` reyestrdən həll edilir.
    `kriging_overrides` Phase A interpolyatorunun sahələrini birbaşa
    dəyişməyə imkan verir (məs. çarpaz-doğrulamada namizəd variogram
    modelini sabitləmək) — strategiya obyektinə toxunmadan.
    """
    if strategy is None:
        if property_name is None:
            raise ValueError("`strategy` və ya `property_name` verilməlidir.")
        strategy = resolve_strategy(property_name)
    if strategy.is_categorical:
        raise ValueError(
            f"'{strategy.name}' KATEQORİKDİR — `interpolate_categorical_field()` "
            "işlədin (kəsilməz kriginq kateqorik koda tətbiq EDİLMİR, GATE B4).")

    targets = np.atleast_2d(np.asarray(targets, float))
    if targets.shape[1] == 2:
        targets = np.column_stack([targets, np.zeros(targets.shape[0])])
    m = targets.shape[0]
    warnings: List[str] = []

    # ── 1. məlumat keyfiyyəti ─────────────────────────────────────────
    if run_qc:
        qc_result = run_quality_control(points, values, strategy, qc)
        clean_points, clean_values = qc_result.points, qc_result.values
        report = qc_result.report
        warnings.extend(report.warnings)
    else:
        clean_points = np.atleast_2d(np.asarray(points, float))
        if clean_points.shape[1] == 2:
            clean_points = np.column_stack([clean_points,
                                            np.zeros(clean_points.shape[0])])
        clean_values = np.asarray(values, float).ravel()
        report = None

    if clean_values.size == 0:
        return _empty_estimate(strategy, m, report,
                               warnings + ["Etibarlı sərt data yoxdur — nəticə NaN."])

    # ── 2. çevirmə (datadan asılıdırsa TƏMİZ datadan fit edilir) ──────
    transform = strategy.fit_transform(clean_values)
    try:
        transformed_values = transform.forward(clean_values)
    except TransformError as exc:
        raise TransformError(
            f"'{strategy.name}': çevirmə tətbiq edilə bilmədi — {exc}") from exc

    # ── 3. Phase A kriginqi ÇEVRİLMİŞ fəzada ──────────────────────────
    interpolator = _build_interpolator(strategy)
    for key, value in (kriging_overrides or {}).items():
        setattr(interpolator, key, value)

    if isinstance(interpolator, OrdinaryKriging):
        kriging = interpolator.krige(clean_points, transformed_values, targets)
        transformed_estimate = kriging.estimate
        transformed_variance = kriging.variance
        lagrange = kriging.lagrange
        support = kriging.support
        neighbor_count = kriging.neighbor_count
        nearest_distance = kriging.nearest_distance
        extrapolated = kriging.extrapolated
        solver = kriging.solver
        total_sill = float(kriging.nugget + kriging.sill)
        range_major = float(kriging.range_)
        anisotropy = kriging.anisotropy
        warnings.extend(kriging.warnings)
    else:
        # IDW/NearestNeighbour — Phase A-nın sadə üsulları: VARİANS YOXDUR
        # və UYDURULMUR (B2: "do not fake uncertainty").
        kriging = None
        transformed_estimate = np.asarray(
            interpolator.interpolate(clean_points[:, :2], transformed_values,
                                     targets[:, :2]), float)
        transformed_variance = np.full(m, np.nan)
        lagrange = np.full(m, np.nan)
        support = np.full(m, SUPPORT_WELL, dtype=object)
        neighbor_count = np.full(m, clean_values.size, dtype=int)
        nearest_distance = np.full(m, np.nan)
        extrapolated = np.zeros(m, dtype=bool)
        solver = np.full(m, strategy.interpolation.value, dtype=object)
        total_sill = float(np.var(transformed_values))
        span = float(np.ptp(clean_points[:, :2])) if clean_points.shape[0] > 1 else 1.0
        range_major = max(span / 3.0, 1e-6)
        anisotropy = strategy.anisotropy or AnisotropyParams()
        warnings.append(
            f"'{strategy.interpolation.value}' üsulunun kriginq variansı YOXDUR — "
            "qeyri-müəyyənlik UYDURULMUR, `variance` NaN qalır.")

    # ── 4. geri çevirmə + variansın köçürülməsi ───────────────────────
    raw_estimate = transform.inverse(
        transformed_estimate,
        variance=None if np.all(~np.isfinite(transformed_variance))
        else transformed_variance,
        lagrange=lagrange, mode=strategy.back_transform)
    if np.all(~np.isfinite(transformed_variance)):
        back_variance, variance_kind = np.full(m, np.nan), VarianceKind.UNDEFINED
    else:
        back_variance, variance_kind = transform.inverse_variance(
            transformed_estimate, transformed_variance)

    # ── 5. fiziki hədlər — AÇIQ, sayılan ──────────────────────────────
    bounded, violated = strategy.apply_output_bounds(raw_estimate)
    if np.any(violated):
        action = ("kəsildi (clip)" if strategy.bound_policy is BoundPolicy.CLIP
                  else "YALNIZ işarələndi")
        warnings.append(
            f"'{strategy.name}': {int(np.sum(violated))}/{m} hüceyrədə nəticə fiziki "
            f"hədləri {strategy.output_bounds} keçdi — {action}. Səbəb adətən "
            "ekstrapolyasiya və ya çox böyük kriginq variansıdır; xam dəyər "
            "`raw_estimate`-də saxlanılır.")

    # ── 6. diaqnostika + etimad ───────────────────────────────────────
    density = compute_data_density(clean_points, targets, range_major, anisotropy)
    confidence = classify_confidence(support, transformed_variance, total_sill,
                                     neighbor_count)

    return PropertyEstimate(
        property_name=strategy.name, estimate=bounded, variance=back_variance,
        std=np.sqrt(np.clip(back_variance, 0.0, None)), raw_estimate=raw_estimate,
        transformed_estimate=transformed_estimate,
        transformed_variance=transformed_variance, variance_kind=variance_kind,
        neighbor_count=np.asarray(neighbor_count, int),
        nearest_distance=np.asarray(nearest_distance, float),
        data_density=density, support=np.asarray(support, dtype=object),
        extrapolated=np.asarray(extrapolated, bool), confidence=confidence,
        solver=np.asarray(solver, dtype=object), bound_adjusted=violated,
        strategy=strategy, kriging=kriging, quality=report, warnings=warnings)


def _empty_estimate(strategy: PropertyStrategy, m: int,
                    report: Optional[DataQualityReport],
                    warnings: List[str]) -> PropertyEstimate:
    nan = np.full(m, np.nan)
    return PropertyEstimate(
        property_name=strategy.name, estimate=nan.copy(), variance=nan.copy(),
        std=nan.copy(), raw_estimate=nan.copy(), transformed_estimate=nan.copy(),
        transformed_variance=nan.copy(), variance_kind=VarianceKind.UNDEFINED,
        neighbor_count=np.zeros(m, dtype=int), nearest_distance=np.full(m, np.inf),
        data_density=np.zeros(m, dtype=int),
        support=np.full(m, SUPPORT_EXTRAPOLATED, dtype=object),
        extrapolated=np.ones(m, dtype=bool),
        confidence=np.full(m, Confidence.EXTRAPOLATED.value, dtype=object),
        solver=np.full(m, "no_neighbors", dtype=object),
        bound_adjusted=np.zeros(m, dtype=bool), strategy=strategy,
        quality=report, warnings=warnings)


# ── kateqorik xassə yolu (B1.6) ───────────────────────────────────────
def interpolate_categorical_field(points, codes, targets,
                                  strategy: Optional[PropertyStrategy] = None,
                                  property_name: Optional[str] = None,
                                  qc: Optional[QCConfig] = None,
                                  run_qc: bool = True,
                                  kriging_overrides: Optional[Dict[str, object]] = None
                                  ) -> CategoricalEstimate:
    """İNDİKATOR kriginq → kateqoriya ehtimalları (B1.6).

    Hər kateqoriya `k` üçün indikator `I_k(x) ∈ {0,1}` qurulur və Phase
    A-nın adi kriginqi ilə interpolyasiya olunur. Xam indikator kriginq
    nəticəsi `[0,1]`-dən kənara çıxa VƏ cəmi 1 OLMAYA bilər (bu, indikator
    kriginqin MƏLUM və qaçılmaz xüsusiyyətidir, gizlədilmir):

        1. hər ehtimal `[0, 1]`-ə kəsilir;
        2. sətir cəmi 1-ə normallaşdırılır;
        3. DÜZƏLİŞ EDİLƏN hüceyrələr SAYILIR (`n_probability_corrections`).

    Nəticədə `Σ P(k|x) = 1` və `0 ≤ P(k|x) ≤ 1` TƏMİN EDİLİR. Bütün
    ehtimallar sıfıra kəsilirsə (heç bir məlumat yoxdur) müşahidə olunan
    QLOBAL nisbətlərə keçilir və bu, xəbərdarlıqla bildirilir.
    """
    if strategy is None:
        if property_name is None:
            raise ValueError("`strategy` və ya `property_name` verilməlidir.")
        strategy = resolve_strategy(property_name)
    if not strategy.is_categorical:
        raise ValueError(f"'{strategy.name}' kateqorik deyil — "
                         "`interpolate_property_field()` işlədin.")

    targets = np.atleast_2d(np.asarray(targets, float))
    if targets.shape[1] == 2:
        targets = np.column_stack([targets, np.zeros(targets.shape[0])])
    m = targets.shape[0]
    warnings: List[str] = []

    if run_qc:
        qc_result = run_quality_control(points, codes, strategy, qc)
        clean_points, clean_codes = qc_result.points, qc_result.values
        report = qc_result.report
        warnings.extend(report.warnings)
    else:
        clean_points = np.atleast_2d(np.asarray(points, float))
        if clean_points.shape[1] == 2:
            clean_points = np.column_stack([clean_points,
                                            np.zeros(clean_points.shape[0])])
        clean_codes = np.asarray(codes, float).ravel()
        report = None

    categories = (np.asarray(strategy.categories, int) if strategy.categories
                  else np.unique(np.rint(clean_codes).astype(int)))
    categories = np.asarray(sorted(set(categories.tolist())), dtype=int)
    if categories.size == 0 or clean_codes.size == 0:
        raise ValueError(f"'{strategy.name}': heç bir etibarlı kateqoriya kodu yoxdur.")

    integer_codes = np.rint(clean_codes).astype(int)
    global_proportions = np.array(
        [float(np.mean(integer_codes == code)) for code in categories])

    if categories.size == 1:
        probabilities = np.ones((m, 1))
        most_probable = np.full(m, int(categories[0]))
        entropy = np.zeros(m)
        return CategoricalEstimate(
            property_name=strategy.name, categories=categories,
            probabilities=probabilities, most_probable=most_probable, entropy=entropy,
            max_probability=np.ones(m), neighbor_count=np.full(m, clean_codes.size, int),
            nearest_distance=np.zeros(m), support=np.full(m, SUPPORT_WELL, dtype=object),
            extrapolated=np.zeros(m, dtype=bool),
            confidence=np.full(m, Confidence.HIGH.value, dtype=object),
            n_probability_corrections=0, strategy=strategy, quality=report,
            warnings=warnings + ["Yalnız bir kateqoriya var — qeyri-müəyyənlik yoxdur."])

    raw = np.empty((m, categories.size))
    support = None
    neighbor_count = np.zeros(m, dtype=int)
    nearest_distance = np.full(m, np.inf)
    extrapolated = np.zeros(m, dtype=bool)
    for index, code in enumerate(categories):
        indicator = (integer_codes == code).astype(float)
        interpolator = _build_interpolator(strategy.derive(
            interpolation=InterpolationKind.KRIGING,
            variable_type=VariableType.CONTINUOUS))
        for key, value in (kriging_overrides or {}).items():
            setattr(interpolator, key, value)
        result = interpolator.krige(clean_points, indicator, targets)
        raw[:, index] = result.estimate
        if support is None:
            support = np.asarray(result.support, dtype=object)
            neighbor_count = np.asarray(result.neighbor_count, int)
            nearest_distance = np.asarray(result.nearest_distance, float)
            extrapolated = np.asarray(result.extrapolated, bool)
        warnings.extend(w for w in result.warnings if w not in warnings)

    clipped = np.clip(np.nan_to_num(raw, nan=0.0), 0.0, 1.0)
    row_sum = clipped.sum(axis=1)
    needs_fix = (np.abs(row_sum - 1.0) > 1e-9) | ~np.all(np.isfinite(raw), axis=1)
    empty = row_sum <= 1e-12
    if np.any(empty):
        clipped[empty] = global_proportions
        row_sum = clipped.sum(axis=1)
        warnings.append(
            f"{int(np.sum(empty))} hüceyrədə bütün indikator ehtimalları sıfıra "
            "kəsildi — müşahidə olunan QLOBAL fasiya nisbətlərinə keçildi "
            "(məkan məlumatı yoxdur, bu hüceyrələr üçün proqnoz məkanca "
            "əsaslandırılmayıb).")
    probabilities = clipped / row_sum[:, None]

    most_probable = categories[np.argmax(probabilities, axis=1)]
    max_probability = probabilities.max(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        logs = np.where(probabilities > 0.0, np.log(probabilities), 0.0)
    entropy = -np.sum(probabilities * logs, axis=1)

    n_corrections = int(np.sum(needs_fix))
    if n_corrections:
        warnings.append(
            f"{n_corrections}/{m} hüceyrədə xam indikator ehtimalları [0,1] "
            "aralığından kənarda idi və/və ya cəmi 1 deyildi — kəsilib "
            "normallaşdırıldı (indikator kriginqin məlum xüsusiyyəti; "
            "düzəliş SAYILIR, gizlədilmir).")

    normalized_entropy = entropy / np.log(categories.size)
    confidence = np.full(m, Confidence.HIGH.value, dtype=object)
    confidence[normalized_entropy > 0.5] = Confidence.MEDIUM.value
    confidence[normalized_entropy > 0.85] = Confidence.LOW.value
    confidence[extrapolated] = Confidence.EXTRAPOLATED.value

    return CategoricalEstimate(
        property_name=strategy.name, categories=categories, probabilities=probabilities,
        most_probable=most_probable, entropy=entropy, max_probability=max_probability,
        neighbor_count=neighbor_count, nearest_distance=nearest_distance,
        support=support if support is not None else np.full(m, SUPPORT_WELL, dtype=object),
        extrapolated=extrapolated, confidence=confidence,
        n_probability_corrections=n_corrections, strategy=strategy, quality=report,
        warnings=warnings)


def interpolate_by_name(points, values, targets, property_name: str, **kwargs):
    """Adına görə düzgün yolu seçən körpü: kateqorik → indikator,
    kəsilməz → çevirməli kriginq. Çağıran tipi bilmək məcburiyyətində
    deyil, amma nəticə tipi FƏRQLİDİR (bilərəkdən — kateqorik nəticəni
    kəsilməz kimi göstərmək məhz qadağan olunan şeydir)."""
    strategy = resolve_strategy(property_name)
    if strategy.is_categorical:
        return interpolate_categorical_field(points, values, targets,
                                             strategy=strategy, **kwargs)
    return interpolate_property_field(points, values, targets, strategy=strategy,
                                      **kwargs)
