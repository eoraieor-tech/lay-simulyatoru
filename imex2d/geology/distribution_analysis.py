"""Xassə paylanması təhlili — SGS-dən ƏVVƏL (Phase 5 §4).

Heç bir xassə "kor-koranə" Gauss/log-çevirməyə göndərilmir — bu modul
əvvəlcə TƏSVİRİ statistikanı (mean/median/var/std/quantiles/skewness)
hesablayır və LOG fəzasının uyğun olub-olmadığını (perkeabilite üçün)
DATA ilə (çarpıqlıq müqayisəsi) qərar verir, sabit ehtimal ilə DEYİL.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.stats import skew


@dataclass
class DistributionSummary:
    n: int
    mean: float
    median: float
    variance: float
    std: float
    minimum: float
    maximum: float
    p10: float
    p50: float
    p90: float
    skewness: Optional[float]

    def as_dict(self) -> dict:
        return {"n": self.n, "mean": self.mean, "median": self.median,
                "variance": self.variance, "std": self.std, "min": self.minimum,
                "max": self.maximum, "p10": self.p10, "p50": self.p50, "p90": self.p90,
                "skewness": self.skewness}


def summarize_distribution(values) -> DistributionSummary:
    values = np.asarray(values, float).ravel()
    if values.size == 0:
        raise ValueError("Paylanma təhlili boş massiv üçün aparıla bilməz.")
    if np.any(~np.isfinite(values)):
        raise ValueError("Paylanma təhlili NaN/sonsuz dəyər qəbul etmir.")
    skewness = float(skew(values)) if values.size >= 3 and np.ptp(values) > 1e-12 else None
    return DistributionSummary(
        n=int(values.size), mean=float(np.mean(values)), median=float(np.median(values)),
        variance=float(np.var(values, ddof=0)), std=float(np.std(values, ddof=0)),
        minimum=float(values.min()), maximum=float(values.max()),
        p10=float(np.percentile(values, 10)), p50=float(np.percentile(values, 50)),
        p90=float(np.percentile(values, 90)), skewness=skewness)


def log_transform_is_justified(values, skew_threshold: float = 1.0) -> bool:
    """Log-fəzanın uyğun olub-olmadığını DATA ilə qərarlaşdırır (kor-
    koranə HƏR keçiricilik sütununa tətbiq etmək ƏVƏZİNƏ, tapşırıq §4):
    dəyərlərin ÖZ çarpıqlığı `log(dəyər)`-in çarpıqlığından DAHA BÖYÜKSƏ
    (mütləq qiymətcə), log fəzası statistik cəhətdən daha simmetrikdir —
    üstünlük verilir. Bütün dəyərlər müsbət olmalıdır (əks halda log
    mümkün deyil, `False` qaytarılır)."""
    values = np.asarray(values, float).ravel()
    if values.size < 3 or np.any(values <= 0) or np.any(~np.isfinite(values)):
        return False
    if np.ptp(values) < 1e-12:
        return False   # sabit dəyər — çevirmə lazımsızdır
    direct_skew = abs(float(skew(values)))
    log_skew = abs(float(skew(np.log(values))))
    return log_skew < direct_skew
