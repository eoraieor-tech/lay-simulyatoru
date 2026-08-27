"""Quyu məlumatı — interpolyasiyanın giriş nöqtəsi.

Bu, karotaj interpretasiyasının NƏTİCƏSİDİR: hər quyuda hesablanmış
petrofiziki parametrlər (məsaməlilik, keçiricilik, gillilik...).
Loqların özü burada saxlanılmır — yalnız interpretasiya olunmuş dəyərlər.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

import numpy as np


@dataclass
class WellSample:
    """Bir quyuda (və istəyə görə bir təbəqədə) ölçülmüş dəyərlər."""
    well: str
    x: float
    y: float
    values: Dict[str, float] = field(default_factory=dict)
    layer: Optional[int] = None      # None -> bütün təbəqələrə aiddir
    depth: Optional[float] = None


@dataclass
class WellDataset:
    samples: List[WellSample] = field(default_factory=list)
    source: str = ""

    def __len__(self) -> int:
        return len(self.samples)

    @property
    def well_names(self) -> List[str]:
        seen = []
        for sample in self.samples:
            if sample.well not in seen:
                seen.append(sample.well)
        return seen

    def property_names(self) -> List[str]:
        names: Set[str] = set()
        for sample in self.samples:
            names.update(sample.values)
        return sorted(names)

    @property
    def layers(self) -> List[Optional[int]]:
        return sorted({s.layer for s in self.samples},
                      key=lambda k: (k is None, k))

    def is_layered(self) -> bool:
        return any(s.layer is not None for s in self.samples)

    def points(self, prop: str, layer: Optional[int] = None):
        """(koordinatlar, dəyərlər) — verilmiş xassə və təbəqə üçün."""
        coordinates, values = [], []
        for sample in self.samples:
            if prop not in sample.values:
                continue
            if layer is not None and sample.layer is not None and sample.layer != layer:
                continue
            coordinates.append((sample.x, sample.y))
            values.append(sample.values[prop])
        return np.asarray(coordinates, float).reshape(-1, 2), np.asarray(values, float)

    def bounds(self):
        if not self.samples:
            return (0.0, 0.0, 0.0, 0.0)
        xs = [s.x for s in self.samples]
        ys = [s.y for s in self.samples]
        return (min(xs), max(xs), min(ys), max(ys))

    def validate(self) -> list:
        issues = []
        if not self.samples:
            issues.append("Quyu məlumatı boşdur.")
            return issues
        if not self.property_names():
            issues.append("Heç bir xassə sütunu tapılmadı.")
        for prop in self.property_names():
            _, values = self.points(prop)
            if values.size < 2:
                issues.append(f"'{prop}' üçün ən azı iki nöqtə lazımdır "
                              f"(tapıldı: {values.size}).")
            if np.any(~np.isfinite(values)):
                issues.append(f"'{prop}' sütununda etibarsız dəyər var.")
        return issues

    def summary(self) -> dict:
        x_min, x_max, y_min, y_max = self.bounds()
        return {"quyu": len(self.well_names), "nöqtə": len(self.samples),
                "xassə": ", ".join(self.property_names()) or "—",
                "təbəqəli": self.is_layered(),
                "sahə": f"X {x_min:.0f}–{x_max:.0f}, Y {y_min:.0f}–{y_max:.0f} m"}
