"""Cross-validation — Kriging/IDW proqnozunun REAL dəqiqliyini ölçür.

Heç bir model "100% dəqiq" deyil. Bu modul mövcud quyu nöqtələrindən
bəzilərini müvəqqəti gizlədib qalanları ilə proqnozlaşdırır, sonra
gizlədilmiş nöqtənin ƏSL (ölçülmüş) dəyəri ilə müqayisə edir. Nəticə pis
çıxa bilər (az nöqtə ilə R² mənfi ola bilər) — bu GİZLƏDİLMİR, olduğu
kimi göstərilir (bax `CrossValidationResult.warnings`).

İki üsul:
    leave_one_out — hər nöqtə növbə ilə gizlədilir (n dəfə yenidən fit).
    k_fold        — nöqtələr `k` təsadüfi qrupa bölünür, hər qrup bir
                    dəfə gizlədilib qalanlarla proqnozlaşdırılır.

Metriklər: MAE, RMSE, R² (ss_tot≈0 olanda təyin edilmir, NaN),
MAPE (sıfıra-yaxın ölçülmüş dəyərlər hesaba qatılmır, sayı bildirilir).
`compute_log_metrics=True` olanda (keçiricilik üçün tövsiyə olunur, bax
`geology_service.DEFAULT_RULES`-də `log_transform=True`) eyni metriklər
log(dəyər) fəzasında da hesablanır — Kriging özü də log fəzasında
işlədiyi üçün bu, statistik cəhətdən daha uyğundur, amma mD-də oxunan
xəta da fiziki mənalı olduğu üçün İKİSİ DƏ göstərilir.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..interfaces.interpolation import IPropertyInterpolator
from .interpolation import interpolate_property
from .variogram import KNOWN_MODELS


@dataclass
class CrossValidationResult:
    method: str
    n_points: int
    folds: int
    mae: float
    rmse: float
    r2: float
    mape: Optional[float]
    mape_excluded: int
    mae_log: Optional[float] = None
    rmse_log: Optional[float] = None
    r2_log: Optional[float] = None
    warnings: List[str] = field(default_factory=list)

    def as_text(self, label: str = "") -> str:
        head = f"{label} " if label else ""
        r2_text = f"{self.r2:.4f}" if np.isfinite(self.r2) else "təyin edilmədi"
        mape_text = f"{self.mape:.2f}%" if self.mape is not None else "təyin edilmədi"
        lines = [
            f"{head}({self.method}, n={self.n_points}, {self.folds} qat):",
            f"  RMSE = {self.rmse:.4g}",
            f"  MAE  = {self.mae:.4g}",
            f"  R²   = {r2_text}",
            f"  MAPE = {mape_text}",
        ]
        if self.rmse_log is not None:
            r2_log_text = f"{self.r2_log:.4f}" if np.isfinite(self.r2_log) else "təyin edilmədi"
            lines.append(
                f"  (log fəzasında) RMSE = {self.rmse_log:.4g}  "
                f"MAE = {self.mae_log:.4g}  R² = {r2_log_text}")
        for warning in self.warnings:
            lines.append(f"  ⚠ {warning}")
        return "\n".join(lines)


def _summarize(method: str, folds: int, actual: np.ndarray, predicted: np.ndarray,
              compute_log_metrics: bool) -> CrossValidationResult:
    actual = np.asarray(actual, float)
    predicted = np.asarray(predicted, float)
    mask = np.isfinite(predicted)
    dropped = int(predicted.size - np.sum(mask))
    warnings: List[str] = []
    if dropped:
        warnings.append(
            f"{dropped} nöqtə üçün proqnoz NaN oldu (axtarış radiusu/qonşu limiti "
            "ilə uyğun məlumat tapılmadı) — bu nöqtələr metrikdən çıxarıldı.")
    actual, predicted = actual[mask], predicted[mask]
    if actual.size == 0:
        return CrossValidationResult(
            method, 0, folds, float("nan"), float("nan"), float("nan"), None, 0,
            warnings=warnings + ["Heç bir etibarlı proqnoz qalmadı."])

    residual = predicted - actual
    mae = float(np.mean(np.abs(residual)))
    rmse = float(np.sqrt(np.mean(residual ** 2)))
    ss_res = float(np.sum(residual ** 2))
    ss_tot = float(np.sum((actual - actual.mean()) ** 2))
    if ss_tot < 1e-12:
        r2 = float("nan")
        warnings.append("R² təyin edilmədi (ölçülmüş dəyərlər demək olar sabitdir).")
    else:
        r2 = 1.0 - ss_res / ss_tot

    safe = np.abs(actual) > 1e-9
    excluded = int(np.sum(~safe))
    if np.any(safe):
        mape = float(np.mean(np.abs(residual[safe] / actual[safe]))) * 100.0
    else:
        mape = None
        warnings.append("MAPE təyin edilmədi (bütün ölçülmüş dəyərlər sıfıra çox yaxındır).")
    if excluded:
        warnings.append(
            f"MAPE {excluded} sıfıra-yaxın nöqtəni nəzərə almadı (sıfıra bölmənin "
            "qarşısı alınıb).")

    result = CrossValidationResult(method, int(actual.size), folds, mae, rmse, r2,
                                   mape, excluded, warnings=warnings)
    if compute_log_metrics:
        positive = (actual > 0) & (predicted > 0)
        if np.sum(positive) >= 2:
            log_a, log_p = np.log(actual[positive]), np.log(predicted[positive])
            log_res = log_p - log_a
            result.mae_log = float(np.mean(np.abs(log_res)))
            result.rmse_log = float(np.sqrt(np.mean(log_res ** 2)))
            ss_res_l = float(np.sum(log_res ** 2))
            ss_tot_l = float(np.sum((log_a - log_a.mean()) ** 2))
            result.r2_log = float("nan") if ss_tot_l < 1e-12 else 1.0 - ss_res_l / ss_tot_l
        else:
            result.warnings.append(
                "log-fəza metrikləri üçün kifayət qədər müsbət cüt yoxdur.")
    return result


def leave_one_out(interpolator: IPropertyInterpolator, points: np.ndarray,
                  values: np.ndarray, log_transform: bool = False,
                  compute_log_metrics: bool = False) -> CrossValidationResult:
    """Hər nöqtəni növbə ilə gizlədib qalan nöqtələrlə proqnozlaşdırır."""
    points = np.asarray(points, float)
    values = np.asarray(values, float)
    n = values.size
    if n < 3:
        raise ValueError(f"Cross-validation üçün ən azı 3 nöqtə lazımdır (tapıldı: {n}).")
    predicted = np.empty(n)
    mask = np.ones(n, dtype=bool)
    for i in range(n):
        mask[:] = True
        mask[i] = False
        predicted[i] = interpolate_property(
            interpolator, points[mask], values[mask], points[i:i + 1],
            log_transform=log_transform)[0]
    return _summarize("leave-one-out", n, values, predicted, compute_log_metrics)


def k_fold(interpolator: IPropertyInterpolator, points: np.ndarray, values: np.ndarray,
          k: int = 5, seed: int = 42, log_transform: bool = False,
          compute_log_metrics: bool = False) -> CrossValidationResult:
    """`k` təsadüfi qatlı spatial cross-validation.

    Nöqtə sayı `k+1`-dən azdırsa, etibarsız (bəlkə boş) qatlar yaratmaq
    əvəzinə açıq xəbərdarlıqla leave-one-out-a keçir — dəyər UYDURULMUR.
    """
    points = np.asarray(points, float)
    values = np.asarray(values, float)
    n = values.size
    if n < k + 1:
        result = leave_one_out(interpolator, points, values, log_transform,
                               compute_log_metrics)
        result.warnings.insert(
            0, f"{k}-fold üçün nöqtə azdır (n={n} < {k + 1}) — leave-one-out işlədildi.")
        return result

    rng = np.random.default_rng(seed)
    order = rng.permutation(n)
    folds = np.array_split(order, k)
    predicted = np.full(n, np.nan)
    for fold in folds:
        if fold.size == 0:
            continue
        mask = np.ones(n, dtype=bool)
        mask[fold] = False
        predicted[fold] = interpolate_property(
            interpolator, points[mask], values[mask], points[fold],
            log_transform=log_transform)
    return _summarize(f"{k}-fold", k, values, predicted, compute_log_metrics)


# ══════════════════════════════════════════════════════════════════════
# B3 — ÇARPAZ-DOĞRULAMA MODEL SEÇİMİ SİSTEMİ KİMİ
# ══════════════════════════════════════════════════════════════════════
"""Aşağıdakı hissə çarpaz-doğrulamanı "metrik hesabla"dan MODEL SEÇİMİ
sisteminə çevirir. Yuxarıdakı `leave_one_out`/`k_fold` funksiyaları
DƏYİŞMƏYİB (mövcud çağıranlar üçün) — bura YENİ, xassə-yönlü qatdır.

SIZMA (data leakage) QAYDASI — B3.1, bu qatın ƏSAS zəmanəti:

    Gizlədilmiş `i` nöqtəsi HEÇ NƏYƏ təsir etməməlidir:
      · variogram fitinə          · qonşuluq seçiminə
      · Kriging çəkilərinə        · çevirmə statistikasına
      · kənar-dəyər qərarlarına   · model seçiminə

    Bu, TƏSADÜFƏN deyil, MEMARLIQLA təmin edilir: hər qat üçün
    `property_interpolation.interpolate_property_field()` YALNIZ TƏLİM
    alt-çoxluğu ilə çağırılır. Həmin funksiya QC-ni, çevirməni və
    variogram fitini ÖZÜ, ötürülən nöqtələrdən qurur — yəni test nöqtəsi
    heç bir mərhələyə fiziki olaraq DAXİL OLMUR. `test_model_selection.py`
    bunu ölçü ilə sübut edir (gizlədilmiş nöqtəni dəyişmək proqnozu
    DƏYİŞMİR).
"""


class ValidationKind(str, Enum):
    """Doğrulama dizaynı (B3.3)."""

    LEAVE_ONE_OUT = "loo"
    RANDOM_KFOLD = "random_kfold"
    SPATIAL_BLOCK = "spatial_block"


@dataclass(frozen=True)
class ValidationDesign:
    """Qatların NECƏ qurulduğu.

    `SPATIAL_BLOCK` niyə vacibdir (B3.3): təsadüfi bölgü ilə gizlədilmiş
    nöqtənin qonşuları demək olar həmişə TƏLİM çoxluğunda qalır, ona görə
    metrik "interpolyasiya" performansını ölçür. Məkan blokları BÜTÖV
    bölgələri gizlədir — bu, modelin MƏLUMATSIZ sahədə davranışını ölçür
    və ekstrapolyasiya keyfiyyəti barədə yeganə dürüst göstəricidir.
    Təsadüfi bölgüyə əsaslanıb "məkan ümumiləşdirməsi yaxşıdır" DEMƏK
    OLMAZ — ona görə hər nəticə hansı dizaynla alındığını daşıyır.
    """

    kind: ValidationKind = ValidationKind.LEAVE_ONE_OUT
    k: int = 5
    seed: int = 42
    #: `SPATIAL_BLOCK` üçün blok ölçüsü (uzunluq vahidi). `None` olanda
    #: domen `k`-ya bölünəcək şəkildə avtomatik hesablanır.
    block_size: Optional[float] = None

    def describe(self) -> str:
        if self.kind is ValidationKind.LEAVE_ONE_OUT:
            return "leave-one-out"
        if self.kind is ValidationKind.RANDOM_KFOLD:
            return f"{self.k}-fold (təsadüfi, seed={self.seed})"
        size = "avtomatik" if self.block_size is None else f"{self.block_size:g}"
        return f"{self.k}-fold (MƏKAN blokları, blok={size})"


def build_folds(points: np.ndarray, design: ValidationDesign
                ) -> List[Tuple[np.ndarray, np.ndarray]]:
    """`[(təlim_indeksləri, test_indeksləri), …]` — DETERMİNİSTİK.

    Məkan bloklarında: domen bərabər kvadrat bloklara bölünür, hər blok
    bütövlükdə BİR qata düşür (blok id-sinə görə növbəli paylama), ona
    görə eyni bloka düşən qonşular heç vaxt təlim/test arasında BÖLÜNMÜR.
    """
    points = np.atleast_2d(np.asarray(points, float))
    n = points.shape[0]
    if n < 2:
        raise ValueError(f"Doğrulama üçün ən azı 2 nöqtə lazımdır (tapıldı: {n}).")
    indices = np.arange(n)

    if design.kind is ValidationKind.LEAVE_ONE_OUT:
        return [(np.delete(indices, i), np.array([i])) for i in range(n)]

    if design.kind is ValidationKind.RANDOM_KFOLD:
        k = max(2, min(design.k, n))
        order = np.random.default_rng(design.seed).permutation(n)
        folds = np.array_split(order, k)
        return [(np.setdiff1d(indices, fold, assume_unique=False), np.sort(fold))
                for fold in folds if fold.size]

    # ── məkan blokları ────────────────────────────────────────────────
    xy = points[:, :2]
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    span = np.maximum(hi - lo, 1e-12)
    k = max(2, min(design.k, n))
    if design.block_size is not None:
        block = float(design.block_size)
    else:
        # domeni təxminən `2k` bloka bölürük ki, hər qata bir neçə AYRI
        # blok düşsün (tək böyük blok qatı təmsilçi olmazdı)
        block = float(np.max(span) / max(int(np.ceil(np.sqrt(2 * k))), 1))
    block = max(block, 1e-9)
    cell = np.floor((xy - lo) / block).astype(int)
    block_ids = cell[:, 0] * 100003 + cell[:, 1]      # deterministik qarışdırma
    unique_blocks = np.unique(block_ids)
    assignment = {int(b): i % k for i, b in enumerate(unique_blocks)}
    fold_of_point = np.array([assignment[int(b)] for b in block_ids])

    folds = []
    for fold_id in range(k):
        test = indices[fold_of_point == fold_id]
        if test.size == 0:
            continue
        folds.append((indices[fold_of_point != fold_id], test))
    return folds


@dataclass
class ContinuousCVMetrics:
    """Kəsilməz xassə üçün doğrulama metrikləri (B3.2)."""

    n: int
    rmse: float
    mae: float
    bias: float                       #: orta xəta (ME) — işarəli
    r2: float
    #: STANDARTLAŞDIRILMIŞ xəta `e = (z − ẑ)/σ` statistikası (B2.2).
    #: İdeal: orta ≈ 0, dispersiya ≈ 1. Dispersiya ≫ 1 → varians AZ
    #: qiymətləndirilib; ≪ 1 → ÇOX qiymətləndirilib.
    mean_standardized_error: float = float("nan")
    variance_standardized_error: float = float("nan")
    coverage_68: float = float("nan")   #: |e| ≤ 1 payı (ideal 0.683)
    coverage_95: float = float("nan")   #: |e| ≤ 1.96 payı (ideal 0.954)
    n_with_variance: int = 0
    #: hədli xassələr üçün: proqnozun fiziki hədləri pozma sayı
    bound_violations: int = 0
    #: ÇEVRİLMİŞ fəzada RMSE (çevirmə varsa) — model müqayisəsi üçün
    transformed_rmse: float = float("nan")
    n_failed: int = 0                  #: NaN proqnoz sayı
    warnings: List[str] = field(default_factory=list)

    @property
    def calibration_error(self) -> float:
        """`|var(e) − 1|` — 0 = mükəmməl kalibrlənmiş qeyri-müəyyənlik."""
        if not np.isfinite(self.variance_standardized_error):
            return float("nan")
        return abs(self.variance_standardized_error - 1.0)

    def as_dict(self) -> Dict[str, float]:
        return {"n": self.n, "rmse": self.rmse, "mae": self.mae, "bias": self.bias,
                "r2": self.r2, "mean_std_error": self.mean_standardized_error,
                "var_std_error": self.variance_standardized_error,
                "coverage_68": self.coverage_68, "coverage_95": self.coverage_95,
                "bound_violations": self.bound_violations,
                "transformed_rmse": self.transformed_rmse, "n_failed": self.n_failed}

    def as_text(self) -> str:
        r2 = f"{self.r2:.4f}" if np.isfinite(self.r2) else "təyin edilmədi"
        lines = [f"n={self.n}  RMSE={self.rmse:.5g}  MAE={self.mae:.5g}  "
                 f"bias={self.bias:+.4g}  R²={r2}"]
        if self.n_with_variance:
            lines.append(
                f"  standartlaşdırılmış xəta: orta={self.mean_standardized_error:+.3f} "
                f"(ideal 0) dispersiya={self.variance_standardized_error:.3f} (ideal 1) "
                f"| örtük 68%={self.coverage_68:.3f} 95%={self.coverage_95:.3f}")
        if self.bound_violations:
            lines.append(f"  fiziki hədd pozması: {self.bound_violations}")
        if self.n_failed:
            lines.append(f"  proqnoz alınmadı: {self.n_failed}")
        return "\n".join(lines)


@dataclass
class CategoricalCVMetrics:
    """Kateqorik xassə üçün doğrulama metrikləri (B3.2)."""

    n: int
    accuracy: float
    log_loss: float          #: −(1/n)·Σ ln P(doğru kateqoriya) — kiçik yaxşıdır
    brier_score: float       #: çoxsinifli Brier: (1/n)·Σ Σ (p_k − 1{k=doğru})²
    categories: Tuple[int, ...] = ()
    confusion: Optional[np.ndarray] = None   #: (k,k) sətir=doğru, sütun=proqnoz
    n_failed: int = 0
    warnings: List[str] = field(default_factory=list)

    def as_dict(self) -> Dict[str, object]:
        return {"n": self.n, "accuracy": self.accuracy, "log_loss": self.log_loss,
                "brier_score": self.brier_score, "categories": list(self.categories),
                "confusion": None if self.confusion is None else self.confusion.tolist(),
                "n_failed": self.n_failed}

    def as_text(self) -> str:
        lines = [f"n={self.n}  dəqiqlik={self.accuracy:.4f}  "
                 f"log-loss={self.log_loss:.4f}  Brier={self.brier_score:.4f}"]
        if self.confusion is not None:
            lines.append("  qarışıqlıq matrisi (sətir=doğru, sütun=proqnoz): "
                         + str(self.confusion.tolist()))
        return "\n".join(lines)


def cross_validate_property(points, values, strategy, design: Optional[ValidationDesign] = None,
                            qc=None, kriging_overrides: Optional[Dict[str, object]] = None):
    """SIZMASIZ çarpaz-doğrulama — kəsilməz və kateqorik xassələr üçün.

    Hər qatda YALNIZ təlim alt-çoxluğu ilə tam boru xətti (QC → çevirmə →
    variogram → kriginq) yenidən qurulur; test nöqtələri heç bir mərhələyə
    daxil olmur (B3.1).
    """
    from .property_interpolation import (interpolate_categorical_field,
                                         interpolate_property_field)

    design = design or ValidationDesign()
    points = np.atleast_2d(np.asarray(points, float))
    values = np.asarray(values, float).ravel()
    if points.shape[0] != values.size:
        raise ValueError("points və values uzunluğu uyğun gəlmir.")
    folds = build_folds(points, design)

    if strategy.is_categorical:
        return _cross_validate_categorical(points, values, strategy, folds, qc,
                                           kriging_overrides,
                                           interpolate_categorical_field)
    return _cross_validate_continuous(points, values, strategy, folds, qc,
                                      kriging_overrides, interpolate_property_field)


def _cross_validate_continuous(points, values, strategy, folds, qc, overrides, runner):
    n = values.size
    predicted = np.full(n, np.nan)
    variance = np.full(n, np.nan)
    transformed_predicted = np.full(n, np.nan)
    transformed_actual = np.full(n, np.nan)
    bound_violations = 0
    warnings: List[str] = []

    for train_idx, test_idx in folds:
        if train_idx.size < 2:
            continue
        estimate = runner(points[train_idx], values[train_idx], points[test_idx],
                          strategy=strategy, qc=qc, kriging_overrides=overrides)
        predicted[test_idx] = estimate.raw_estimate
        variance[test_idx] = estimate.transformed_variance
        transformed_predicted[test_idx] = estimate.transformed_estimate
        bound_violations += int(np.sum(estimate.bound_adjusted))
        # gizlədilmiş nöqtənin ÇEVRİLMİŞ həqiqi dəyəri — çevirmə YALNIZ
        # TƏLİM datasından fit edilib (sızma yoxdur)
        try:
            transform = strategy.fit_transform(values[train_idx])
            transformed_actual[test_idx] = transform.forward(values[test_idx])
        except Exception:      # noqa: BLE001 — çevrilə bilməyən test dəyəri
            pass

    valid = np.isfinite(predicted)
    n_failed = int(np.sum(~valid))
    if n_failed:
        warnings.append(
            f"{n_failed} nöqtə üçün proqnoz alınmadı (qonşu tapılmadı və ya "
            "çevirmə mümkün olmadı) — metriklərdən çıxarıldı.")
    if not np.any(valid):
        return ContinuousCVMetrics(
            n=0, rmse=float("nan"), mae=float("nan"), bias=float("nan"),
            r2=float("nan"), n_failed=n_failed,
            warnings=warnings + ["Heç bir etibarlı proqnoz qalmadı."])

    residual = predicted[valid] - values[valid]
    rmse = float(np.sqrt(np.mean(residual ** 2)))
    mae = float(np.mean(np.abs(residual)))
    bias = float(np.mean(residual))
    ss_tot = float(np.sum((values[valid] - values[valid].mean()) ** 2))
    r2 = float("nan") if ss_tot < 1e-12 else 1.0 - float(np.sum(residual ** 2)) / ss_tot

    metrics = ContinuousCVMetrics(
        n=int(np.sum(valid)), rmse=rmse, mae=mae, bias=bias, r2=r2,
        bound_violations=bound_violations, n_failed=n_failed, warnings=warnings)

    # ── standartlaşdırılmış xəta (B2.2) — ÇEVRİLMİŞ fəzada ────────────
    # Varians məhz orada ƏSL kriginq variansıdır; geri çevrilmiş varians
    # (delta metodu) kalibrləmə üçün uyğun deyil.
    usable = (valid & np.isfinite(variance) & (variance > 1e-15)
              & np.isfinite(transformed_actual) & np.isfinite(transformed_predicted))
    if np.any(usable):
        standardized = ((transformed_actual[usable] - transformed_predicted[usable])
                        / np.sqrt(variance[usable]))
        metrics.n_with_variance = int(np.sum(usable))
        metrics.mean_standardized_error = float(np.mean(standardized))
        metrics.variance_standardized_error = float(np.var(standardized, ddof=0))
        metrics.coverage_68 = float(np.mean(np.abs(standardized) <= 1.0))
        metrics.coverage_95 = float(np.mean(np.abs(standardized) <= 1.96))
        metrics.transformed_rmse = float(np.sqrt(np.mean(
            (transformed_predicted[usable] - transformed_actual[usable]) ** 2)))
    else:
        warnings.append(
            "Kriginq variansı mövcud deyil — standartlaşdırılmış xəta/kalibrləmə "
            "metrikləri hesablanmadı (UYDURULMUR).")
    return metrics


def _cross_validate_categorical(points, values, strategy, folds, qc, overrides, runner):
    n = values.size
    codes = np.rint(values).astype(int)
    categories = np.asarray(sorted(set(codes.tolist())), dtype=int)
    k = categories.size
    probabilities = np.full((n, k), np.nan)
    predicted = np.full(n, -1, dtype=int)
    warnings: List[str] = []

    for train_idx, test_idx in folds:
        if train_idx.size < 2:
            continue
        train_categories = np.unique(codes[train_idx])
        estimate = runner(points[train_idx], values[train_idx], points[test_idx],
                          strategy=strategy.derive(
                              categories=tuple(int(c) for c in train_categories)),
                          qc=qc, kriging_overrides=overrides)
        # təlim qatında olmayan kateqoriya üçün ehtimal 0 qalır — bu,
        # doğrudur: model onu HEÇ GÖRMƏYİB, ehtimal UYDURULMUR
        column = {int(c): i for i, c in enumerate(categories)}
        block = np.zeros((test_idx.size, k))
        for local_index, code in enumerate(estimate.categories):
            block[:, column[int(code)]] = estimate.probabilities[:, local_index]
        probabilities[test_idx] = block
        predicted[test_idx] = estimate.most_probable

    valid = predicted >= 0
    n_failed = int(np.sum(~valid))
    if not np.any(valid):
        return CategoricalCVMetrics(0, float("nan"), float("nan"), float("nan"),
                                    tuple(categories.tolist()), None, n_failed,
                                    ["Heç bir etibarlı proqnoz qalmadı."])

    truth = codes[valid]
    accuracy = float(np.mean(predicted[valid] == truth))
    column = {int(c): i for i, c in enumerate(categories)}
    truth_index = np.array([column[int(c)] for c in truth])
    prob = np.clip(np.nan_to_num(probabilities[valid], nan=0.0), 0.0, 1.0)
    row_sum = prob.sum(axis=1, keepdims=True)
    prob = np.where(row_sum > 0, prob / np.maximum(row_sum, 1e-15), 1.0 / k)

    #: log-loss üçün sıfır ehtimalı kəsmə həddi — `ln(0) = −∞` sonsuz cəza
    #: verərdi; `1e-12` sənədləşdirilmiş, bütün namizədlərə EYNİ tətbiq olunan
    #: hədddir, ona görə müqayisəni təhrif etmir.
    truth_probability = prob[np.arange(truth.size), truth_index]
    log_loss = float(-np.mean(np.log(np.maximum(truth_probability, 1e-12))))
    onehot = np.zeros_like(prob)
    onehot[np.arange(truth.size), truth_index] = 1.0
    brier = float(np.mean(np.sum((prob - onehot) ** 2, axis=1)))

    confusion = np.zeros((k, k), dtype=int)
    for actual, guess in zip(truth, predicted[valid]):
        confusion[column[int(actual)], column[int(guess)]] += 1

    if n_failed:
        warnings.append(f"{n_failed} nöqtə üçün proqnoz alınmadı.")
    return CategoricalCVMetrics(
        n=int(np.sum(valid)), accuracy=accuracy, log_loss=log_loss, brier_score=brier,
        categories=tuple(int(c) for c in categories), confusion=confusion,
        n_failed=n_failed, warnings=warnings)


# ── model seçimi (B3.4/B3.5/B3.6) ─────────────────────────────────────
@dataclass
class ModelCandidate:
    """Bir namizəd model = etiket + TAM xassə strategiyası."""

    label: str
    strategy: object                      #: `property_config.PropertyStrategy`
    kriging_overrides: Dict[str, object] = field(default_factory=dict)


#: Çox-metrikli sıralamanın DEFOLT çəkiləri (B3.5).
#:
#: Hər meyar [0, ∞) aralığında "cərimə" kimi normallaşdırılır (0 = ən
#: yaxşı), sonra çəkili cəmlənir. Çəkilər SƏBƏBİ İLƏ belədir:
#:
#:   accuracy 0.50 — proqnoz xətası əsas məqsəddir, amma TƏK meyar deyil;
#:   calibration 0.25 — qeyri-müəyyənlik yanlışdırsa, aşağı RMSE də
#:                      aldadıcıdır (SGS/qərar analizi ona güvənir);
#:   bias 0.15 — sistematik sürüşmə həcm hesablamasını birbaşa təhrif edir;
#:   validity 0.07 — fiziki hədd pozması modelin uyğunsuzluğunun əlamətidir;
#:   stability 0.03 — alınmayan proqnozlar (NaN) az, amma sıfır olmayan cəza.
#:
#: Cəmi 1.0. Çağıran öz çəkilərini verə bilər — dəyərlər SABİT KODLANMIR.
DEFAULT_SELECTION_WEIGHTS: Dict[str, float] = {
    "accuracy": 0.50, "calibration": 0.25, "bias": 0.15,
    "validity": 0.07, "stability": 0.03,
}


@dataclass
class CandidateResult:
    candidate: ModelCandidate
    metrics: object                      #: `ContinuousCVMetrics`/`CategoricalCVMetrics`
    score: float = float("nan")
    penalties: Dict[str, float] = field(default_factory=dict)
    error: Optional[str] = None

    @property
    def ok(self) -> bool:
        return self.error is None and self.metrics is not None


@dataclass
class ModelSelectionReport:
    """Struktur, maşınla oxuna bilən model-seçimi hesabatı (B3.6)."""

    property_name: str
    design: ValidationDesign
    n_samples: int
    results: List[CandidateResult]
    selected: Optional[CandidateResult]
    weights: Dict[str, float]
    n_excluded: int = 0
    warnings: List[str] = field(default_factory=list)

    @property
    def ranking(self) -> List[CandidateResult]:
        """Ən yaxşıdan ən pisə — bal, sonra etiket (deterministik)."""
        usable = [r for r in self.results if r.ok and np.isfinite(r.score)]
        return sorted(usable, key=lambda r: (r.score, r.candidate.label))

    def as_dict(self) -> Dict[str, object]:
        return {
            "property": self.property_name,
            "design": {"kind": self.design.kind.value, "k": self.design.k,
                       "seed": self.design.seed, "description": self.design.describe()},
            "n_samples": self.n_samples, "n_excluded": self.n_excluded,
            "weights": dict(self.weights),
            "selected": None if self.selected is None else self.selected.candidate.label,
            "candidates": [
                {"label": r.candidate.label, "score": r.score,
                 "penalties": dict(r.penalties), "error": r.error,
                 "metrics": None if r.metrics is None else r.metrics.as_dict()}
                for r in self.results],
            "warnings": list(self.warnings),
        }

    def as_text(self) -> str:
        lines = [f"Model seçimi — {self.property_name}",
                 f"  dizayn: {self.design.describe()} · nümunə: {self.n_samples}",
                 "  çəkilər: " + ", ".join(f"{k}={v:g}" for k, v in self.weights.items())]
        for rank, result in enumerate(self.ranking, start=1):
            marker = "★" if result is self.selected else " "
            lines.append(f"  {marker}{rank}. {result.candidate.label:<34} "
                         f"bal={result.score:.4f}")
            lines.append("      " + result.metrics.as_text().replace("\n", "\n      "))
        for result in self.results:
            if not result.ok:
                lines.append(f"   ✗ {result.candidate.label}: {result.error}")
        lines.extend(f"  ⚠ {w}" for w in self.warnings)
        return "\n".join(lines)


def _relative_penalty(value: float, best: float) -> float:
    """`value/best − 1` — ən yaxşı namizədə nəzərən NİSBİ pisləşmə."""
    if not np.isfinite(value) or not np.isfinite(best) or best <= 0.0:
        return 0.0 if value == best else 1.0
    return max(float(value) / float(best) - 1.0, 0.0)


def score_candidates(results: List[CandidateResult],
                     weights: Optional[Dict[str, float]] = None) -> None:
    """Namizədləri ÇOX-METRİKLİ, deterministik balla qiymətləndirir (B3.5).

    Bal `0` = bütün meyarlarda ən yaxşı. Aşağı bal daha yaxşıdır.
    Bütün cərimələr ən yaxşı namizədə NƏZƏRƏN nisbidir, ona görə
    metriklərin miqyası (mD vs kəsr) nəticəyə təsir etmir.
    """
    weights = dict(weights or DEFAULT_SELECTION_WEIGHTS)
    usable = [r for r in results if r.ok]
    if not usable:
        return

    categorical = isinstance(usable[0].metrics, CategoricalCVMetrics)
    if categorical:
        best_loss = min(r.metrics.log_loss for r in usable
                        if np.isfinite(r.metrics.log_loss))
        best_error = min(1.0 - r.metrics.accuracy for r in usable)
        for result in usable:
            m = result.metrics
            penalties = {
                "accuracy": _relative_penalty(1.0 - m.accuracy, max(best_error, 1e-9)),
                "calibration": _relative_penalty(m.log_loss, best_loss),
                "bias": float(m.brier_score),
                "validity": 0.0,
                "stability": float(m.n_failed) / max(m.n + m.n_failed, 1),
            }
            result.penalties = penalties
            result.score = float(sum(weights.get(k, 0.0) * v
                                     for k, v in penalties.items()))
        return

    finite_rmse = [r.metrics.rmse for r in usable if np.isfinite(r.metrics.rmse)]
    if not finite_rmse:
        return
    best_rmse = min(finite_rmse)

    # KALİBRLƏMƏ MƏLUMATI OLMAYAN NAMİZƏD (məs. IDW — kriginq variansı
    # yoxdur) NEYTRAL 0 cərimə ALMIR. Səbəb (ölçülmüş nəticə): neytral 0
    # ilə IDW balın 25%-ini PULSUZ qazanırdı və PORO ssenarisində RMSE-si
    # DAHA PİS olmasına baxmayaraq (0.02684 vs 0.02594) qalib gəlirdi —
    # yəni qeyri-müəyyənliyi ÜMUMİYYƏTLƏ verə bilməyən üsul, verə bilən
    # üsuldan üstün tutulurdu. Bu, Phase B-nin bütün məntiqinə ziddir.
    #
    # ƏVƏZİNƏ: qiymətləndirilə bilməyən meyar üçün MÜŞAHİDƏ OLUNAN ƏN PİS
    # cərimə verilir (çox-meyarlı qərar analizində standart mühafizəkar
    # doldurma). Yəni "yoxlaya bilmirik" → "kredit vermirik". Heç bir
    # namizəddə kalibrləmə yoxdursa, hamısı 0 alır (müqayisə pozulmur).
    calibrations = [m.calibration_error for m in (r.metrics for r in usable)
                    if np.isfinite(m.calibration_error)]
    worst_calibration = max(calibrations) if calibrations else 0.0

    for result in usable:
        m = result.metrics
        calibration = m.calibration_error
        penalties = {
            "accuracy": _relative_penalty(m.rmse, best_rmse),
            # kalibrləmə: `|var(e) − 1|`; ölçülə bilmirsə ən pis müşahidə
            "calibration": (float(calibration) if np.isfinite(calibration)
                            else float(worst_calibration)),
            "bias": (abs(m.bias) / m.rmse if np.isfinite(m.bias) and m.rmse > 0 else 0.0),
            "validity": float(m.bound_violations) / max(m.n, 1),
            "stability": float(m.n_failed) / max(m.n + m.n_failed, 1),
        }
        result.penalties = penalties
        result.score = float(sum(weights.get(k, 0.0) * v for k, v in penalties.items()))


def select_property_model(points, values, candidates: Sequence[ModelCandidate],
                          property_name: str = "",
                          design: Optional[ValidationDesign] = None,
                          weights: Optional[Dict[str, float]] = None,
                          qc=None) -> ModelSelectionReport:
    """Namizədləri SIZMASIZ çarpaz-doğrulama ilə müqayisə edib seçir (B3.4).

    "Adi kriginq həmişə ən yaxşıdır" kimi SABİT qərar YOXDUR — qalib
    yalnız doğrulama metriklərindən çıxır. Bərabər ballarda etiket üzrə
    əlifba sırası ilə deterministik seçim edilir (təkrarlana bilən).
    """
    design = design or ValidationDesign()
    values = np.asarray(values, float).ravel()
    results: List[CandidateResult] = []
    warnings: List[str] = []

    for candidate in candidates:
        try:
            metrics = cross_validate_property(
                points, values, candidate.strategy, design, qc=qc,
                kriging_overrides=candidate.kriging_overrides or None)
            results.append(CandidateResult(candidate=candidate, metrics=metrics))
        except Exception as exc:            # noqa: BLE001 — namizəd uğursuzluğu
            results.append(CandidateResult(candidate=candidate, metrics=None,
                                           error=f"{type(exc).__name__}: {exc}"))

    score_candidates(results, weights)
    report = ModelSelectionReport(
        property_name=property_name or getattr(candidates[0].strategy, "name", ""),
        design=design, n_samples=int(values.size), results=results, selected=None,
        weights=dict(weights or DEFAULT_SELECTION_WEIGHTS),
        n_excluded=sum(1 for r in results if not r.ok), warnings=warnings)
    ranking = report.ranking
    report.selected = ranking[0] if ranking else None
    if report.selected is None:
        report.warnings.append(
            "Heç bir namizəd doğrulanmadı — model seçilə bilmədi.")
    if report.n_excluded:
        report.warnings.append(
            f"{report.n_excluded} namizəd xəta ilə kənarda qaldı (səbəblər yuxarıda).")
    return report


def default_candidates(property_name: str,
                       include_idw: bool = True) -> List[ModelCandidate]:
    """Bir xassə üçün standart namizəd dəsti (B3.4).

    Kəsilməz xassələr: IDW + üç variogram modeli ilə kriginq + (loq-normal
    xassələr üçün) xam-fəza müqayisəsi + hədli xassələr üçün çevirməsiz
    müqayisə. Kateqorik xassələr: üç variogram modeli ilə indikator kriginq.

    Xam/çevirməsiz variant QƏSDƏN daxil edilir: "loq fəzası daha yaxşıdır"
    iddiası da DATA ilə YOXLANILMALIDIR, fərz edilməməlidir."""
    from .property_config import (BackTransform, InterpolationKind, VariableType,
                                  resolve_strategy)
    from .transforms import IDENTITY_TRANSFORM

    base = resolve_strategy(property_name)
    candidates: List[ModelCandidate] = []

    if base.is_categorical:
        for model in KNOWN_MODELS:
            candidates.append(ModelCandidate(
                label=f"indikator kriginq + {model}",
                strategy=base.derive(variogram_model=model)))
        return candidates

    for model in KNOWN_MODELS:
        candidates.append(ModelCandidate(
            label=f"kriginq + {model} ({base.transform.name})",
            strategy=base.derive(variogram_model=model)))

    if not base.transform.is_identity:
        # çevirməsiz (xam fəza) müqayisə — hipotezi yoxlamaq üçün
        raw = base.derive(transform=IDENTITY_TRANSFORM,
                          back_transform=BackTransform.MEDIAN,
                          variable_type=VariableType.CONTINUOUS,
                          variogram_model="spherical", legacy_log_transform=False)
        candidates.append(ModelCandidate(label="kriginq + spherical (xam fəza)",
                                         strategy=raw))
    if include_idw:
        candidates.append(ModelCandidate(
            label="IDW (çevirməli)",
            strategy=base.derive(interpolation=InterpolationKind.IDW)))
    return candidates
