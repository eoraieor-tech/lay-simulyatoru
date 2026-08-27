"""Müşahidə məlumatı — history matching-in giriş nöqtəsi.

Bu, yataqdan ölçülmüş faktiki hasilat tarixçəsidir: quyu debitləri,
quyudibi təzyiqlər, sulaşma. Model bunları təkrarlaya bilirsə,
proqnozuna güvənmək olar.

Vahidlər modelin öz sistemi ilə eynidir (METRIC): m3/gün, bar, gün.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional

import numpy as np


class ObservedQuantity(Enum):
    """Ölçülə bilən kəmiyyətlər."""
    OIL_RATE = "neft debiti"
    WATER_RATE = "su debiti"
    WATER_INJECTION = "su vurma"
    WATER_CUT = "sulaşma"
    BOTTOM_HOLE_PRESSURE = "quyudibi təzyiq"
    AVERAGE_PRESSURE = "orta lay təzyiqi"
    CUMULATIVE_OIL = "kumulyativ neft"

    @property
    def unit(self) -> str:
        return {
            ObservedQuantity.OIL_RATE: "m³/gün",
            ObservedQuantity.WATER_RATE: "m³/gün",
            ObservedQuantity.WATER_INJECTION: "m³/gün",
            ObservedQuantity.WATER_CUT: "%",
            ObservedQuantity.BOTTOM_HOLE_PRESSURE: "bar",
            ObservedQuantity.AVERAGE_PRESSURE: "bar",
            ObservedQuantity.CUMULATIVE_OIL: "m³",
        }[self]

    @property
    def is_field_level(self) -> bool:
        """Quyuya deyil, bütün yataq üçün ölçülən kəmiyyət."""
        return self in (ObservedQuantity.AVERAGE_PRESSURE,)


@dataclass
class ObservedSeries:
    """Bir quyunun bir kəmiyyəti üzrə zaman sırası."""
    well: str
    quantity: ObservedQuantity
    time: np.ndarray            # gün
    values: np.ndarray
    uncertainty: Optional[np.ndarray] = None
    """Ölçmə xətası (σ). Verilməyibsə, çəki hesablanarkən dəyərin
    özündən nisbi qiymət götürülür."""

    def __post_init__(self):
        self.time = np.asarray(self.time, dtype=float).ravel()
        self.values = np.asarray(self.values, dtype=float).ravel()
        if self.uncertainty is not None:
            self.uncertainty = np.asarray(self.uncertainty, float).ravel()

    def __len__(self) -> int:
        return int(self.time.size)

    @property
    def key(self) -> tuple:
        return (self.well, self.quantity)

    FIELD_LABEL = "Yataq"

    @property
    def label(self) -> str:
        if self.quantity.is_field_level or not self.well:
            return f"{self.FIELD_LABEL} · {self.quantity.value}"
        return f"{self.well} · {self.quantity.value}"

    def validate(self) -> list:
        issues = []
        if self.time.size != self.values.size:
            issues.append(f"{self.label}: zaman və dəyər sayı fərqlidir.")
        if self.time.size < 2:
            issues.append(f"{self.label}: ən azı iki nöqtə lazımdır.")
        if np.any(np.diff(self.time) <= 0):
            issues.append(f"{self.label}: zaman artan sıralı olmalıdır.")
        if np.any(~np.isfinite(self.values)):
            issues.append(f"{self.label}: etibarsız dəyər var.")
        return issues


@dataclass
class ObservationSet:
    """Bütün müşahidələr — bir yataq tarixçəsi."""
    series: List[ObservedSeries] = field(default_factory=list)
    source: str = ""

    def __len__(self) -> int:
        return len(self.series)

    @property
    def wells(self) -> List[str]:
        seen = []
        for item in self.series:
            if item.well and item.well not in seen:
                seen.append(item.well)
        return seen

    @property
    def quantities(self) -> List[ObservedQuantity]:
        seen = []
        for item in self.series:
            if item.quantity not in seen:
                seen.append(item.quantity)
        return seen

    @property
    def time_span(self) -> tuple:
        if not self.series:
            return (0.0, 0.0)
        starts = [item.time.min() for item in self.series]
        ends = [item.time.max() for item in self.series]
        return (float(min(starts)), float(max(ends)))

    def get(self, well: str, quantity: ObservedQuantity
            ) -> Optional[ObservedSeries]:
        for item in self.series:
            if item.well == well and item.quantity is quantity:
                return item
        return None

    def validate(self) -> list:
        issues = []
        if not self.series:
            issues.append("Müşahidə məlumatı boşdur.")
        for item in self.series:
            issues.extend(item.validate())
        return issues

    def summary(self) -> Dict[str, object]:
        start, end = self.time_span
        return {
            "sıra": len(self.series),
            "quyu": len(self.wells),
            "kəmiyyət": ", ".join(q.value for q in self.quantities) or "—",
            "müddət": f"{start:.0f} – {end:.0f} gün",
            "nöqtə": sum(len(item) for item in self.series),
        }
