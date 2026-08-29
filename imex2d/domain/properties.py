"""Xassə xəritələri və süxur/flüid xassələri."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class PropertyMap:
    """Adlandırılmış hüceyrə massivi. Hər xassə öz adını və vahidini daşıyır."""
    name: str
    values: np.ndarray
    unit: str = ""

    @classmethod
    def uniform(cls, name: str, value: float, ncell: int, unit: str = "") -> "PropertyMap":
        return cls(name, np.full(ncell, float(value)), unit)

    @classmethod
    def from_array(cls, name: str, arr, ncell: int, unit: str = "") -> "PropertyMap":
        a = np.asarray(arr, dtype=float)
        if a.ndim == 0:
            return cls.uniform(name, float(a), ncell, unit)
        if a.size != ncell:
            raise ValueError(f"{name}: {a.size} dəyər, gözlənilən {ncell}")
        return cls(name, a.ravel().copy(), unit)

    def as_grid(self, shape) -> np.ndarray:
        return self.values.reshape(shape)

    def stats(self) -> dict:
        return {"min": float(self.values.min()), "max": float(self.values.max()),
                "mean": float(self.values.mean())}


@dataclass
class RockProperties:
    """Statik süxur xassələri — geoloji modeldən gəlir."""
    porosity: PropertyMap
    permx: PropertyMap
    permy: PropertyMap
    permz: Optional[PropertyMap] = None
    net_to_gross: Optional[PropertyMap] = None
    compressibility: float = 4.5e-5

    def validate(self) -> list:
        issues = []
        if np.any(self.porosity.values <= 0):
            issues.append("Məsaməlilik sıfır və ya mənfi hüceyrələr var.")
        if np.any(self.permx.values <= 0) or np.any(self.permy.values <= 0):
            issues.append("Keçiricilik sıfır və ya mənfi hüceyrələr var.")
        return issues


@dataclass
class FluidProperties:
    """PLACEHOLDER — sabit flüid xassələri.

    Bu sinif müvəqqətidir. Real black-oil davranışı IPVTProvider vasitəsilə
    gələcək; PVT provider inject edilmədikdə mühərrik bu dəyərləri oxuyur.
    """
    water_viscosity: float = 0.5
    oil_viscosity: float = 3.0
    water_fvf: float = 1.0
    oil_fvf: float = 1.15
    water_compressibility: float = 4.4e-5
    oil_compressibility: float = 1.4e-4
    water_density: float = 1010.0   # səth sıxlığı, kg/m3
    oil_density: float = 850.0      # səth sıxlığı, kg/m3
