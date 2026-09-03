"""MƏLUMAT KEYFİYYƏTİ — variogramdan və interpolyasiyadan ƏVVƏL (B4).

Boru xətti (GATE B8: QC variogramdan ƏVVƏL gəlir):

    XAM DATA
      → KOORDİNAT DOĞRULAMASI     (sonluluq, ölçü, uyğunluq)
      → SONSUZ/ÇATIŞMAYAN DƏYƏR   (açıq siyasət)
      → DUBLİKAT AŞKARLANMASI     (açıq, deterministik həll)
      → FİZİKİ ETİBARLILIQ        (xassəyə-xas hədlər)
      → KƏNAR-DƏYƏR DİAQNOSTİKASI (işarələ — SİLMƏ)
      → TƏMİZ SƏRT DATA + HESABAT

ÜÇ FƏRQLİ PROBLEMİ QARIŞDIRMIRIQ (tapşırıq B4.4):

    1. FİZİKİ CƏHƏTDƏN ETİBARSIZ  — mənfi məsaməlik, SW = 1.4.
       Bu, ölçmə/emal səhvidir; defolt olaraq ÇIXARILIR (səbəb yazılır).
    2. STATİSTİK KƏNAR-DƏYƏR       — 5000 mD keçiricilik.
       Bu, FİZİKİ OLARAQ MÜMKÜNDÜR; defolt olaraq YALNIZ İŞARƏLƏNİR,
       SİLİNMİR (`remove_outliers=True` açıq şəkildə verilməyibsə).
    3. NADİR AMMA ETİBARLI         — çat zonası, yüksək keçiricilik.
       Bu, məlumatın ÖZÜDÜR; heç bir avtomatik qərar verilmir.

Heç bir addım SƏSSİZ deyil: hər qərar `DataQualityReport`-a düşür.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from .anisotropy import AnisotropyParams
from .property_config import DuplicatePolicy, OutlierMethod, PropertyStrategy
from .spatial_search import NeighborhoodConfig, NeighborhoodSelector

#: Robust z üçün normal paylanma sabiti: `0.6745 = Φ⁻¹(0.75)`, MAD-ı
#: standart kənarlanmaya çevirir (normal data üçün `MAD/0.6745 ≈ σ`).
_MAD_SCALE = 0.6745

#: `E|X − μ| = σ·√(2/π) ≈ 0.7979·σ` — ORTA MÜTLƏQ KƏNARLANMANI standart
#: kənarlanmaya çevirən sabit (MAD sıfır olanda ehtiyat miqyas).
_MEAN_AD_SCALE = 0.7978845608028654


def _robust_scale(deviations: np.ndarray) -> float:
    """Robust `σ̂` — İKİ PİLLƏLİ (MAD, sonra orta mütləq kənarlanma).

    NİYƏ İKİ PİLLƏ (ölçülmüş qüsur): MAD kənar dəyərlərə qarşı ən güclü
    miqyasdır, AMMA nümunənin yarıdan çoxu EYNİ dəyərdirsə `MAD = 0` olur
    və bölmə mümkün deyil. Bu, nadir kənar hal DEYİL — məhz TƏMİZ məlumatda
    baş verir:

      · 29 nöqtədə `K = 50 mD`, bir nöqtədə `5000 mD` → MAD = 0;
      · müntəzəm şəbəkədə xətti trend üçün yerli qalıqlar DƏQİQ sıfırdır
        (simmetrik qonşuluğun medianı mərkəzi dəyərə bərabərdir) → MAD = 0.

    Yalnız MAD-a baxan versiya bu hallarda "kənar-dəyər yoxdur" deyirdi —
    yəni aşkarlanma məhz ideal şəraitdə SÖNÜRDÜ (tutulub, bax
    `tests/test_data_quality.py`).

    Ona görə `MAD = 0` olanda ORTA MÜTLƏQ KƏNARLANMAYA keçilir — o, yalnız
    BÜTÜN dəyərlər eyni olanda sıfırdır, yəni "həqiqətən kənar-dəyər yoxdur"
    halını düzgün saxlayır. Hər iki miqyas normal paylanma üçün `σ`-ya
    uyğunlaşdırılıb, ona görə `threshold` eyni mənanı saxlayır."""
    centred = np.abs(deviations - np.median(deviations))
    mad = float(np.median(centred))
    if mad > 1e-12:
        return mad / _MAD_SCALE
    mean_ad = float(np.mean(centred))
    if mean_ad > 1e-12:
        return mean_ad / _MEAN_AD_SCALE
    return 0.0


class DataQualityError(ValueError):
    """QC siyasəti pozulub və siyasət "xəta at"dır — səssiz keçid YOX."""


@dataclass(frozen=True)
class QCConfig:
    """QC davranışının AÇIQ konfiqurasiyası.

    `strategy`-dən gələn defoltlar (dublikat/kənar-dəyər siyasəti) burada
    əvəz edilə bilər; verilməyən sahələr strategiyadan götürülür."""

    non_finite_policy: str = "drop"          #: "drop" | "raise"
    invalid_bounds_policy: str = "drop"      #: "drop" | "raise" | "keep"
    duplicate_policy: Optional[DuplicatePolicy] = None
    duplicate_tolerance: Optional[float] = None
    outlier_method: Optional[OutlierMethod] = None
    outlier_threshold: Optional[float] = None
    remove_outliers: Optional[bool] = None
    #: məkan kənar-dəyər analizi üçün qonşu sayı
    spatial_neighbors: int = 8
    #: `SPATIAL` metodunda anizotrop həndəsə (Phase A ilə EYNİ obyekt)
    anisotropy: Optional[AnisotropyParams] = None

    def resolved(self, strategy: PropertyStrategy) -> "QCConfig":
        """Verilməyən sahələri strategiyadan doldurulmuş nüsxə."""
        return QCConfig(
            non_finite_policy=self.non_finite_policy,
            invalid_bounds_policy=self.invalid_bounds_policy,
            duplicate_policy=(self.duplicate_policy
                              if self.duplicate_policy is not None
                              else strategy.duplicate_policy),
            duplicate_tolerance=(self.duplicate_tolerance
                                 if self.duplicate_tolerance is not None
                                 else strategy.duplicate_tolerance),
            outlier_method=(self.outlier_method if self.outlier_method is not None
                            else strategy.outlier_method),
            outlier_threshold=(self.outlier_threshold
                               if self.outlier_threshold is not None
                               else strategy.outlier_threshold),
            remove_outliers=(self.remove_outliers if self.remove_outliers is not None
                             else strategy.remove_outliers),
            spatial_neighbors=self.spatial_neighbors,
            anisotropy=self.anisotropy)


@dataclass
class QCFinding:
    """Bir QC tapıntısı — NƏ tapıldı, HARADA, NƏ EDİLDİ."""

    kind: str                       #: "non_finite", "duplicate", "bounds", "outlier"…
    severity: str                   #: "info" | "warning" | "error"
    count: int
    action: str                     #: "removed" | "merged" | "flagged" | "kept"
    detail: str
    indices: Tuple[int, ...] = ()   #: XAM giriş massivindəki sətir indeksləri

    def __str__(self) -> str:
        return f"[{self.severity}] {self.kind}: {self.count} → {self.action}. {self.detail}"


@dataclass
class DataQualityReport:
    """Struktur QC hesabatı (B4.6) — maşınla oxuna bilən."""

    property_name: str
    n_input: int = 0
    n_valid: int = 0
    n_removed: int = 0
    n_non_finite: int = 0
    n_invalid_coordinates: int = 0
    n_duplicate_locations: int = 0
    n_duplicate_observations: int = 0
    n_conflicting_duplicates: int = 0
    n_bound_violations: int = 0
    n_outlier_candidates: int = 0
    n_outliers_removed: int = 0
    duplicate_policy: str = ""
    outlier_method: str = ""
    findings: List[QCFinding] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    #: XAM indeks → saxlanıldımı (True) / çıxarıldımı (False)
    kept_mask: Optional[np.ndarray] = None
    #: XAM indeks → kənar-dəyər namizədidirmi (silinməsindən ASILI OLMAYARAQ)
    outlier_mask: Optional[np.ndarray] = None

    @property
    def ok(self) -> bool:
        """Ən azı bir etibarlı müşahidə qaldımı."""
        return self.n_valid > 0

    def add(self, finding: QCFinding) -> None:
        self.findings.append(finding)

    def as_dict(self) -> Dict[str, object]:
        return {
            "property": self.property_name, "n_input": self.n_input,
            "n_valid": self.n_valid, "n_removed": self.n_removed,
            "n_non_finite": self.n_non_finite,
            "n_invalid_coordinates": self.n_invalid_coordinates,
            "n_duplicate_locations": self.n_duplicate_locations,
            "n_duplicate_observations": self.n_duplicate_observations,
            "n_conflicting_duplicates": self.n_conflicting_duplicates,
            "n_bound_violations": self.n_bound_violations,
            "n_outlier_candidates": self.n_outlier_candidates,
            "n_outliers_removed": self.n_outliers_removed,
            "duplicate_policy": self.duplicate_policy,
            "outlier_method": self.outlier_method,
            "findings": [f.__dict__ for f in self.findings],
            "warnings": list(self.warnings),
        }

    def as_text(self) -> str:
        lines = [
            f"Məlumat keyfiyyəti — {self.property_name}",
            f"  giriş {self.n_input} → etibarlı {self.n_valid} "
            f"(çıxarılan {self.n_removed})",
            f"  qeyri-sonlu {self.n_non_finite} · etibarsız koordinat "
            f"{self.n_invalid_coordinates} · hədd pozması {self.n_bound_violations}",
            f"  dublikat mövqe {self.n_duplicate_locations} "
            f"(ziddiyyətli {self.n_conflicting_duplicates}) — siyasət: "
            f"{self.duplicate_policy}",
            f"  kənar-dəyər namizədi {self.n_outlier_candidates} "
            f"(silinən {self.n_outliers_removed}) — metod: {self.outlier_method}",
        ]
        lines.extend(f"  {finding}" for finding in self.findings)
        lines.extend(f"  ⚠ {message}" for message in self.warnings)
        return "\n".join(lines)


@dataclass
class QCResult:
    """QC-dən sonra TƏMİZ sərt data + hesabat."""

    points: np.ndarray
    values: np.ndarray
    report: DataQualityReport
    #: birləşdirilmiş dublikatların XAM indeks qrupları (audit üçün)
    duplicate_groups: Tuple[Tuple[int, ...], ...] = ()


# ── addımlar ──────────────────────────────────────────────────────────
def _as_points(points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, float)
    if points.ndim == 1:
        points = points.reshape(1, -1)
    if points.ndim != 2 or points.shape[1] not in (2, 3):
        raise DataQualityError(
            f"Koordinatlar (n,2) və ya (n,3) olmalıdır, alındı: {points.shape}")
    if points.shape[1] == 2:
        points = np.column_stack([points, np.zeros(points.shape[0])])
    return points


def detect_outliers(values: np.ndarray, method: OutlierMethod, threshold: float,
                    points: Optional[np.ndarray] = None,
                    n_neighbors: int = 8,
                    anisotropy: Optional[AnisotropyParams] = None
                    ) -> Tuple[np.ndarray, np.ndarray]:
    """`(maska, bal)` — kənar-dəyər NAMİZƏDLƏRİ (B4.5). HEÇ NƏ SİLMİR.

    * `MAD` — robust z: `|z − median| / σ̂`, `σ̂` İKİ PİLLƏLİ robust
      miqyasdır (bax `_robust_scale`). Ortadan/standart kənarlanmadan
      fərqli olaraq kənar dəyərlərin ÖZÜNDƏN təsirlənmir. YALNIZ bütün
      dəyərlər eyni olanda heç nə işarələnmir.
    * `IQR` — Tukey çəpərləri `Q1 − k·IQR`, `Q3 + k·IQR` (`k = threshold`,
      ənənəvi 1.5; burada defolt daha mühafizəkardır).
    * `SPATIAL` — MƏKAN kənarı: hər nöqtə ÖZ qonşuluğunun medianı ilə
      müqayisə edilir (qonşuluq Phase A-nın `NeighborhoodSelector`-u ilə,
      EYNİ anizotrop həndəsədə). Qlobal paylanmada normal, amma
      qonşularından kəskin fərqlənən nöqtəni tapır — geostatistikada ən
      mənalı kənar-dəyər növü.
    """
    values = np.asarray(values, float)
    n = values.size
    score = np.zeros(n)
    if method is OutlierMethod.NONE or n < 4:
        return np.zeros(n, dtype=bool), score

    if method is OutlierMethod.MAD:
        scale = _robust_scale(values)
        if scale <= 0.0:
            return np.zeros(n, dtype=bool), score      # bütün dəyərlər eynidir
        score = np.abs(values - np.median(values)) / scale
        return score > threshold, score

    if method is OutlierMethod.IQR:
        q1, q3 = np.percentile(values, [25.0, 75.0])
        iqr = float(q3 - q1)
        if iqr <= 1e-12:
            return np.zeros(n, dtype=bool), score
        low, high = q1 - threshold * iqr, q3 + threshold * iqr
        score = np.maximum((low - values) / iqr, (values - high) / iqr)
        return (values < low) | (values > high), np.maximum(score, 0.0)

    if method is OutlierMethod.SPATIAL:
        if points is None:
            raise DataQualityError(
                "SPATIAL kənar-dəyər analizi üçün koordinatlar lazımdır.")
        selector = NeighborhoodSelector(
            points, anisotropy=anisotropy,
            config=NeighborhoodConfig(max_neighbors=min(n_neighbors + 1, n),
                                      min_neighbors=1))
        batch = selector.select_batch(points)
        residual = np.zeros(n)
        for row in range(n):
            neighbours = batch.neighbours(row)
            neighbours = neighbours[neighbours != row]      # ÖZÜNÜ çıxar
            if neighbours.size < 2:
                continue
            residual[row] = values[row] - float(np.median(values[neighbours]))
        scale = _robust_scale(residual)
        if scale <= 0.0:
            return np.zeros(n, dtype=bool), score      # bütün qalıqlar eynidir
        score = np.abs(residual - np.median(residual)) / scale
        return score > threshold, score

    raise DataQualityError(f"Naməlum kənar-dəyər metodu: {method!r}")


def _resolve_duplicates(points: np.ndarray, values: np.ndarray, indices: np.ndarray,
                        policy: DuplicatePolicy, tolerance: float,
                        strategy: PropertyStrategy, report: DataQualityReport
                        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                   Tuple[Tuple[int, ...], ...]]:
    """Eyni koordinatlı müşahidələri AÇIQ siyasətlə həll edir (B4.2).

    Qaytarır: `(points, values, xam_indekslər, dublikat_qrupları)`.
    Heç bir müşahidə səssizcə ATILMIR — atılan/birləşdirilən hər qrup
    hesabata düşür.
    """
    n = points.shape[0]
    if n < 2:
        return points, values, indices, ()

    # `tolerance` daxilində eyni sayılan koordinatları qruplaşdır:
    # yuvarlaqlaşdırma + `np.unique` (O(n log n), tam (n,n) matris YOX).
    decimals = max(0, int(round(-np.log10(max(tolerance, 1e-15)))))
    rounded = np.round(np.ascontiguousarray(points), decimals)
    view = rounded.view([("", rounded.dtype)] * rounded.shape[1])
    _, first_idx, inverse, counts = np.unique(
        view, return_index=True, return_inverse=True, return_counts=True)
    inverse = inverse.reshape(-1)
    duplicate_groups_ids = np.where(counts > 1)[0]
    if duplicate_groups_ids.size == 0:
        return points, values, indices, ()

    groups: List[Tuple[int, ...]] = []
    conflicting = 0
    n_dupe_observations = 0
    for group_id in duplicate_groups_ids:
        member = np.where(inverse == group_id)[0]
        groups.append(tuple(int(indices[i]) for i in member))
        n_dupe_observations += int(member.size)
        group_values = values[member]
        spread = float(np.max(group_values) - np.min(group_values))
        if spread > 1e-9:
            conflicting += 1

    report.n_duplicate_locations = int(duplicate_groups_ids.size)
    report.n_duplicate_observations = n_dupe_observations
    report.n_conflicting_duplicates = conflicting

    if policy is DuplicatePolicy.RAISE:
        raise DataQualityError(
            f"'{strategy.name}': {duplicate_groups_ids.size} mövqedə dublikat müşahidə "
            f"var ({conflicting} ziddiyyətli). `duplicate_policy` ilə AÇIQ strategiya "
            "seçin (mean/median/keep_first/keep_last/keep_separate/majority).")

    if policy is DuplicatePolicy.KEEP_SEPARATE:
        report.add(QCFinding(
            "duplicate", "warning" if conflicting else "info",
            duplicate_groups_ids.size, "kept",
            "Dublikatlar BİRLƏŞDİRİLMƏDİ (keep_separate) — Kriging matrisi təkil "
            "ola bilər; Phase A solver-i jitter/lstsq ehtiyat yolu ilə həll edir.",
            tuple(int(i) for group in groups for i in group)))
        return points, values, indices, tuple(groups)

    keep = np.ones(n, dtype=bool)
    resolved = values.copy()
    for group_id in duplicate_groups_ids:
        member = np.where(inverse == group_id)[0]
        group_values = values[member]
        target = member[0]
        if policy is DuplicatePolicy.MEAN:
            resolved[target] = float(np.mean(group_values))
        elif policy is DuplicatePolicy.MEDIAN:
            resolved[target] = float(np.median(group_values))
        elif policy is DuplicatePolicy.KEEP_FIRST:
            resolved[target] = group_values[0]
        elif policy is DuplicatePolicy.KEEP_LAST:
            resolved[target] = group_values[-1]
        elif policy is DuplicatePolicy.MAJORITY:
            codes, code_counts = np.unique(np.rint(group_values).astype(int),
                                           return_counts=True)
            winners = codes[code_counts == code_counts.max()]
            if winners.size > 1:
                raise DataQualityError(
                    f"'{strategy.name}': dublikat qrupda səs BƏRABƏRDİR "
                    f"({dict(zip(codes.tolist(), code_counts.tolist()))}) — 'majority' "
                    "həll edə bilmir; keep_first/keep_last seçin və ya datanı düzəldin.")
            resolved[target] = float(winners[0])
        keep[member] = False
        keep[target] = True

    report.add(QCFinding(
        "duplicate", "warning" if conflicting else "info",
        duplicate_groups_ids.size, "merged",
        f"{n_dupe_observations} müşahidə {duplicate_groups_ids.size} mövqedə "
        f"'{policy.value}' siyasəti ilə birləşdirildi ({conflicting} qrupda dəyərlər "
        "ziddiyyətli idi).",
        tuple(int(i) for group in groups for i in group)))
    return points[keep], resolved[keep], indices[keep], tuple(groups)


def run_quality_control(points, values, strategy: PropertyStrategy,
                        config: Optional[QCConfig] = None) -> QCResult:
    """B4 boru xəttinin TAM icrası — variogramdan ƏVVƏL çağırılır.

    Sıra qəsdən belədir: koordinat → sonluluq → dublikat → fiziki hədd →
    kənar-dəyər. Kənar-dəyər ƏN SONDA gəlir, çünki o, qalan (təmiz)
    paylanmaya nəzərən hesablanmalıdır — etibarsız dəyərlər hələ
    içindəykən robust statistika da təhrif olunardı.
    """
    config = (config or QCConfig()).resolved(strategy)
    report = DataQualityReport(
        property_name=strategy.name,
        duplicate_policy=config.duplicate_policy.value,
        outlier_method=config.outlier_method.value)

    points = _as_points(points)
    values = np.asarray(values, float).ravel()
    if values.shape[0] != points.shape[0]:
        raise DataQualityError(
            f"'{strategy.name}': koordinat ({points.shape[0]}) və dəyər "
            f"({values.shape[0]}) sayı uyğun gəlmir.")
    report.n_input = int(values.size)
    indices = np.arange(values.size)

    # ── 1. koordinat doğrulaması (B4.1) ───────────────────────────────
    finite_points = np.all(np.isfinite(points), axis=1)
    if not np.all(finite_points):
        count = int(np.sum(~finite_points))
        report.n_invalid_coordinates = count
        report.add(QCFinding(
            "coordinate", "error", count, "removed",
            "NaN/sonsuz koordinat — nöqtə məkanda yerləşdirilə bilmir, "
            "koordinat TƏXMİN EDİLMİR.",
            tuple(indices[~finite_points].tolist())))
        points, values, indices = (points[finite_points], values[finite_points],
                                   indices[finite_points])

    # ── 2. qeyri-sonlu DƏYƏR (B4.3) ───────────────────────────────────
    finite_values = np.isfinite(values)
    if not np.all(finite_values):
        count = int(np.sum(~finite_values))
        report.n_non_finite = count
        if config.non_finite_policy == "raise":
            raise DataQualityError(
                f"'{strategy.name}': {count} müşahidədə NaN/sonsuz dəyər var "
                "(non_finite_policy='raise').")
        report.add(QCFinding(
            "non_finite", "warning", count, "removed",
            "NaN/±sonsuz dəyər — interpolyasiya sistemini bütövlükdə NaN edərdi.",
            tuple(indices[~finite_values].tolist())))
        points, values, indices = (points[finite_values], values[finite_values],
                                   indices[finite_values])

    # ── 3. dublikatlar (B4.2) ─────────────────────────────────────────
    points, values, indices, groups = _resolve_duplicates(
        points, values, indices, config.duplicate_policy,
        float(config.duplicate_tolerance), strategy, report)

    # ── 4. fiziki etibarlılıq (B4.4) ──────────────────────────────────
    invalid = strategy.invalid_value_mask(values)
    if np.any(invalid):
        count = int(np.sum(invalid))
        report.n_bound_violations = count
        lo, hi = strategy.physical_bounds
        detail = (f"Fiziki etibarlılıq aralığından kənar dəyər "
                  f"(hədlər: {lo}, {hi}; növ: {strategy.variable_type.value}). "
                  "Bu, STATİSTİK kənar-dəyər DEYİL — ölçmə/emal səhvidir.")
        if config.invalid_bounds_policy == "raise":
            raise DataQualityError(f"'{strategy.name}': {count} etibarsız dəyər. {detail}")
        if config.invalid_bounds_policy == "keep":
            report.add(QCFinding("bounds", "warning", count, "kept",
                                 detail + " Siyasət 'keep' — SAXLANILDI; çevirmə "
                                 "mərhələsində xəta verə bilər.",
                                 tuple(indices[invalid].tolist())))
        else:
            report.add(QCFinding("bounds", "error", count, "removed", detail,
                                 tuple(indices[invalid].tolist())))
            keep = ~invalid
            points, values, indices = points[keep], values[keep], indices[keep]

    # ── 5. kənar-dəyər DİAQNOSTİKASI (B4.5) — SİLMİR ──────────────────
    outlier_raw_indices: Tuple[int, ...] = ()
    if values.size and config.outlier_method is not OutlierMethod.NONE:
        flagged, score = detect_outliers(
            values, config.outlier_method, float(config.outlier_threshold),
            points=points, n_neighbors=config.spatial_neighbors,
            anisotropy=config.anisotropy)
        if np.any(flagged):
            count = int(np.sum(flagged))
            report.n_outlier_candidates = count
            outlier_raw_indices = tuple(indices[flagged].tolist())
            worst = float(np.max(score[flagged]))
            if config.remove_outliers:
                report.n_outliers_removed = count
                report.add(QCFinding(
                    "outlier", "warning", count, "removed",
                    f"'{config.outlier_method.value}' metodu ilə hədd "
                    f"{config.outlier_threshold:g} aşıldı (ən böyük bal {worst:.2f}). "
                    "AÇIQ `remove_outliers=True` seçimi ilə çıxarıldı.",
                    outlier_raw_indices))
                keep = ~flagged
                points, values, indices = points[keep], values[keep], indices[keep]
            else:
                report.add(QCFinding(
                    "outlier", "info", count, "flagged",
                    f"'{config.outlier_method.value}' metodu ilə hədd "
                    f"{config.outlier_threshold:g} aşıldı (ən böyük bal {worst:.2f}). "
                    "DEFOLT: SİLİNMİR — fiziki cəhətdən mümkün, nadir dəyər ola "
                    "bilər (çat zonası və s.). Silmək üçün `remove_outliers=True`.",
                    outlier_raw_indices))

    # ── hesabat ───────────────────────────────────────────────────────
    report.n_valid = int(values.size)
    report.n_removed = report.n_input - report.n_valid
    kept_mask = np.zeros(report.n_input, dtype=bool)
    kept_mask[indices] = True
    report.kept_mask = kept_mask
    outlier_mask = np.zeros(report.n_input, dtype=bool)
    if outlier_raw_indices:
        outlier_mask[list(outlier_raw_indices)] = True
    report.outlier_mask = outlier_mask

    if report.n_valid == 0:
        report.warnings.append(
            f"'{strategy.name}': QC-dən sonra heç bir etibarlı müşahidə qalmadı — "
            "interpolyasiya mümkün deyil.")
    elif report.n_valid < 4:
        report.warnings.append(
            f"'{strategy.name}': yalnız {report.n_valid} etibarlı müşahidə qaldı — "
            "deneysel variogram (≥4 nöqtə) qurula bilməyəcək, evristik radiusa "
            "keçiləcək.")
    if report.n_removed:
        report.warnings.append(
            f"{report.n_removed}/{report.n_input} müşahidə QC-də çıxarıldı — "
            "səbəblər yuxarıdakı tapıntılardadır.")

    return QCResult(points=points, values=values, report=report,
                    duplicate_groups=groups)


def quality_control_many(points, value_map: Dict[str, np.ndarray],
                         strategies: Dict[str, PropertyStrategy],
                         config: Optional[QCConfig] = None
                         ) -> Dict[str, QCResult]:
    """Bir neçə xassə üçün QC — hər biri ÖZ strategiyası ilə.

    Xassələr AYRI-AYRI təmizlənir: PERMX-də etibarsız bir ölçmə PORO
    müşahidəsini SİLMİR (sətir yox, DƏYƏR düşür) — bu, `geology_adapter.py`
    -dəki mövcud "xassə-üzrə davranış" qaydası ilə eynidir."""
    return {name: run_quality_control(points, values,
                                      strategies.get(name.upper()) or strategies[name],
                                      config)
            for name, values in value_map.items()}
