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
from typing import List, Optional

import numpy as np

from ..interfaces.interpolation import IPropertyInterpolator
from .interpolation import interpolate_property


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
