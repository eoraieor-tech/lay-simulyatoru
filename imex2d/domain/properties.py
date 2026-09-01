"""Xassə xəritələri və süxur/flüid xassələri."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

import numpy as np

from .unit_conversions import known_units
from .validation import (validate_compressibility, validate_density,
                         validate_permeability, validate_porosity, validate_viscosity)

#: Xassə adı -> gözlənilən kəmiyyət növü (bax `unit_conversions.py`).
#: YALNIZ real vahid-qarışıqlığı riski olan xassələr üçün (Phase 1
#: audit: keçiricilik mD/Darcy/m² arasında, təzyiq bar/psi arasında
#: real çaşqınlıq mənbəyidir). Reyestrdə OLMAYAN ad (PORO, NTG, SW,
#: REGION_ID...) üçün `unit` sərbəst mətn olaraq qalır — DƏYİŞMİR.
PROPERTY_QUANTITY = {
    "PERMX": "permeability", "PERMY": "permeability", "PERMZ": "permeability",
    "PRESSURE": "pressure",
}


@dataclass
class PropertyMap:
    """Adlandırılmış hüceyrə massivi. Hər xassə öz adını və vahidini daşıyır.

    `unit` boş deyilsə (istifadəçi/idxal AÇIQ vahid göstərib) VƏ `name`
    `PROPERTY_QUANTITY`-də qeydiyyatdan keçibsə, vahid bu xassənin
    kəmiyyət növünə (məs. keçiricilik üçün mD/D/m²) uyğun OLMALIDIR —
    əks halda `ValueError` (məs. `PropertyMap("PERMX", ..., "psi")`).
    Boş `unit` (defolt) HEÇ VAXT rədd edilmir — bu, "vahid göstərilməyib"
    halıdır, mövcud mühərrik vahidi kimi qəbul edilir (bax `UNITS.md`).
    """
    name: str
    values: np.ndarray
    unit: str = ""

    def __post_init__(self):
        quantity = PROPERTY_QUANTITY.get(self.name)
        if quantity is not None and self.unit:
            allowed = known_units(quantity)
            if self.unit not in allowed:
                raise ValueError(
                    f"{self.name}: vahid {self.unit!r} bu xassə üçün etibarsızdır "
                    f"(gözlənilən kəmiyyət: {quantity}, dəstəklənən vahidlər: {allowed}).")

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
        """Sərt fiziki xətalar. `validate_warnings()` — qeyri-adi (amma
        mümkün) diapazon xəbərdarlıqları üçün, bax Phase 1 hesabatı."""
        issues = []
        if np.any(self.porosity.values <= 0):
            issues.append("Məsaməlilik sıfır və ya mənfi hüceyrələr var.")
        if np.any(self.porosity.values >= 1.0):
            issues.append("Məsaməlilik >= 1.0 olan hüceyrələr var (fiziki cəhətdən qeyri-mümkün).")
        if np.any(self.permx.values <= 0) or np.any(self.permy.values <= 0):
            issues.append("Keçiricilik sıfır və ya mənfi hüceyrələr var.")
        return issues

    def validate_warnings(self) -> list:
        """Rədd edilməyən, amma qeyri-adi diapazon xəbərdarlıqları."""
        warnings = validate_porosity(self.porosity.values, "PORO").warnings
        warnings += validate_permeability(self.permx.values, "PERMX").warnings
        warnings += validate_permeability(self.permy.values, "PERMY").warnings
        if self.permz is not None:
            warnings += validate_permeability(self.permz.values, "PERMZ").warnings
        return warnings


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

    def validate(self) -> list:
        """Sərt fiziki xətalar (əvvəllər bu sinifdə HEÇ BİR yoxlama yox
        idi — bax Phase 1 audit)."""
        issues = []
        issues += validate_viscosity(self.water_viscosity, "su lözlüyü").errors
        issues += validate_viscosity(self.oil_viscosity, "neft lözlüyü").errors
        issues += validate_density(self.water_density, "su sıxlığı").errors
        issues += validate_density(self.oil_density, "neft sıxlığı").errors
        issues += validate_compressibility(self.water_compressibility, "su sıxılması").errors
        issues += validate_compressibility(self.oil_compressibility, "neft sıxılması").errors
        if self.water_fvf <= 0 or self.oil_fvf <= 0:
            issues.append("Formasiya həcm əmsalı (Bw/Bo) müsbət olmalıdır.")
        return issues

    def validate_warnings(self) -> list:
        warnings = []
        warnings += validate_viscosity(self.water_viscosity, "su lözlüyü").warnings
        warnings += validate_viscosity(self.oil_viscosity, "neft lözlüyü").warnings
        warnings += validate_density(self.water_density, "su sıxlığı").warnings
        warnings += validate_density(self.oil_density, "neft sıxlığı").warnings
        warnings += validate_compressibility(self.water_compressibility, "su sıxılması").warnings
        warnings += validate_compressibility(self.oil_compressibility, "neft sıxılması").warnings
        return warnings
