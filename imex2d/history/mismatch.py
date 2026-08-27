"""Uyğunsuzluq ölçüsü — history matching-in "xəta funksiyası".

Model müşahidəni nə qədər yaxşı təkrarlayır? Bu sual bir rəqəmə
endirilməlidir ki, müxtəlif variantlar müqayisə oluna bilsin.

    SSE   = Σ w·(hesablanmış − ölçülmüş)²        xam kvadrat cəm
    RMSE  = √(SSE / n)                            ölçü vahidində
    NRMSE = RMSE / (müşahidənin diapazonu)        ÖLÇÜSÜZ

Yekun uyğunsuzluq NRMSE-lərin çəkili ortasıdır. Ölçüsüz olması
vacibdir: təzyiq barla, debit m³/günlə ölçülür — xam SSE-ləri toplamaq
böyük ədədli kəmiyyəti üstün edərdi.

ZAMAN UYĞUNLAŞDIRMASI
Simulyator öz addımlarında nəticə verir, müşahidələr isə başqa
tarixlərdə ölçülüb. Model nəticəsi müşahidə vaxtlarına interpolyasiya
olunur — əksinə yox, çünki müşahidə həqiqətdir və dəyişdirilməməlidir.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from ..domain.observations import (ObservationSet, ObservedQuantity,
                                   ObservedSeries)
from ..simulation.results import SimulationResult

# kəmiyyət -> yekun uyğunsuzluqdakı çəkisi
DEFAULT_WEIGHTS: Dict[ObservedQuantity, float] = {
    ObservedQuantity.OIL_RATE: 1.0,
    ObservedQuantity.WATER_RATE: 1.0,
    ObservedQuantity.WATER_CUT: 1.0,
    ObservedQuantity.WATER_INJECTION: 0.5,
    ObservedQuantity.BOTTOM_HOLE_PRESSURE: 0.7,
    ObservedQuantity.AVERAGE_PRESSURE: 0.7,
    ObservedQuantity.CUMULATIVE_OIL: 1.5,
}
"""Kumulyativ neft daha ağırdır: gündəlik debitdə səs-küy çox olur,
toplam isə ehtiyatın nə qədər çıxarıldığını birbaşa göstərir."""


@dataclass
class SeriesMismatch:
    """Bir zaman sırası üzrə nəticə."""
    label: str
    quantity: ObservedQuantity
    well: str
    time: np.ndarray
    observed: np.ndarray
    simulated: np.ndarray
    weight: float = 1.0

    @property
    def residuals(self) -> np.ndarray:
        return self.simulated - self.observed

    @property
    def sse(self) -> float:
        return float(np.sum(self.residuals ** 2))

    @property
    def rmse(self) -> float:
        return float(np.sqrt(np.mean(self.residuals ** 2)))

    @property
    def scale(self) -> float:
        """Normallaşdırma miqyası — müşahidənin diapazonu."""
        span = float(np.ptp(self.observed))
        if span > 1e-12:
            return span
        magnitude = float(np.mean(np.abs(self.observed)))
        return magnitude if magnitude > 1e-12 else 1.0

    @property
    def nrmse(self) -> float:
        return self.rmse / self.scale

    @property
    def bias(self) -> float:
        """Sistematik meyl: müsbət = model çox verir."""
        return float(np.mean(self.residuals))

    @property
    def correlation(self) -> float:
        if self.observed.size < 2:
            return 0.0
        if np.ptp(self.observed) < 1e-12 or np.ptp(self.simulated) < 1e-12:
            return 0.0
        return float(np.corrcoef(self.observed, self.simulated)[0, 1])


@dataclass
class MismatchReport:
    series: List[SeriesMismatch] = field(default_factory=list)
    skipped: List[str] = field(default_factory=list)
    """Modeldə qarşılığı olmayan müşahidələr — susmaq təhlükəlidir."""

    @property
    def total(self) -> float:
        """Yekun ölçüsüz uyğunsuzluq (çəkili orta NRMSE)."""
        if not self.series:
            return float("inf")
        weights = np.array([item.weight for item in self.series])
        values = np.array([item.nrmse for item in self.series])
        return float(np.sum(weights * values) / max(np.sum(weights), 1e-30))

    @property
    def worst(self) -> Optional[SeriesMismatch]:
        return max(self.series, key=lambda item: item.nrmse, default=None)

    def as_rows(self) -> List[tuple]:
        return [(item.label, item.rmse, item.nrmse, item.bias,
                 item.correlation, item.weight)
                for item in sorted(self.series, key=lambda x: -x.nrmse)]

    def as_text(self) -> str:
        lines = [f"Yekun uyğunsuzluq: {self.total:.4f}  "
                 f"({len(self.series)} sıra)"]
        for label, rmse, nrmse, bias, correlation, _ in self.as_rows():
            lines.append(f"  {label:<28} NRMSE {nrmse:6.3f}  "
                         f"RMSE {rmse:10.3f}  meyl {bias:+9.3f}  "
                         f"r {correlation:5.2f}")
        for name in self.skipped:
            lines.append(f"  {name:<28} ATLANDI (modeldə qarşılığı yoxdur)")
        return "\n".join(lines)


class MismatchCalculator:
    """Simulyasiya nəticəsini müşahidə ilə tutuşdurur."""

    def __init__(self, weights: Optional[Dict[ObservedQuantity, float]] = None):
        self.weights = dict(DEFAULT_WEIGHTS)
        if weights:
            self.weights.update(weights)

    def evaluate(self, result: SimulationResult,
                 observations: ObservationSet) -> MismatchReport:
        report = MismatchReport()
        if not result.series.time:
            report.skipped = [item.label for item in observations.series]
            return report

        model_time = np.asarray(result.series.time, float)
        for observed in observations.series:
            simulated = self._model_curve(result, observed)
            if simulated is None:
                report.skipped.append(observed.label)
                continue
            inside = ((observed.time >= model_time[0] - 1e-9)
                      & (observed.time <= model_time[-1] + 1e-9))
            if not np.any(inside):
                report.skipped.append(observed.label)
                continue

            report.series.append(SeriesMismatch(
                label=observed.label,
                quantity=observed.quantity,
                well=observed.well,
                time=observed.time[inside],
                observed=observed.values[inside],
                simulated=np.interp(observed.time[inside], model_time,
                                    simulated),
                weight=self.weights.get(observed.quantity, 1.0)))
        return report

    # ─────────────────────────────────────────── model əyrisinin seçimi
    @staticmethod
    def _model_curve(result: SimulationResult,
                     observed: ObservedSeries) -> Optional[np.ndarray]:
        series = result.series
        quantity = observed.quantity

        if quantity is ObservedQuantity.AVERAGE_PRESSURE:
            return np.asarray(series.average_pressure, float)

        # yataq səviyyəsi: quyu adı verilməyibsə və ya "FIELD" ise
        field_level = observed.well in ("", "FIELD", "FIELD-TOTAL", None)

        if quantity is ObservedQuantity.OIL_RATE:
            if field_level:
                return np.asarray(series.oil_rate, float)
            values = result.well_oil_rate.get(observed.well)
            return np.asarray(values, float) if values else None

        if quantity is ObservedQuantity.WATER_RATE:
            if field_level:
                return np.asarray(series.water_rate, float)
            values = result.well_water_rate.get(observed.well)
            return np.asarray(values, float) if values else None

        if quantity is ObservedQuantity.WATER_INJECTION:
            return np.asarray(series.water_injection_rate, float)

        if quantity is ObservedQuantity.CUMULATIVE_OIL:
            return np.asarray(series.cumulative_oil, float)

        if quantity is ObservedQuantity.WATER_CUT:
            if field_level:
                return np.asarray(series.water_cut, float)
            oil = result.well_oil_rate.get(observed.well)
            water = result.well_water_rate.get(observed.well)
            if not oil or not water:
                return None
            oil = np.asarray(oil, float)
            water = np.asarray(water, float)
            return water / np.maximum(oil + water, 1e-12) * 100.0

        if quantity is ObservedQuantity.BOTTOM_HOLE_PRESSURE:
            return None      # quyudibi təzyiq hələ nəticədə saxlanılmır
        return None
