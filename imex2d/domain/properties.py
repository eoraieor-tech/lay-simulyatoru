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
class PermeabilityTensor:
    """Tam simmetrik permeabilite tenzoru — gələcək MPFA-O üçün HAZIRLIQ
    (bax audit tapşırığı §5, `imex2d/simulation/discretization.py` modul
    docstring-i).

        K = [[Kxx, Kxy, Kxz],
             [Kxy, Kyy, Kyz],
             [Kxz, Kyz, Kzz]]     (simmetrik fərz edilir — standart fiziki qəbul)

    BU FAZADA HEÇ BİR HƏLLEDİCİ (TPFA) bunu İSTİFADƏ ETMİR — `RockProperties.
    permx/permy/permz` (diaqonal) YEGANƏ TPFA-nın oxuduğu mənbədir, DƏYİŞMİR.
    Bu sinif YALNIZ off-diaqonal (Kxy/Kxz/Kyz) anizotropluq məlumatını
    İTİRMƏDƏN daşımaq üçündür ki, gələcək MPFA-O onu birbaşa istifadə edə
    bilsin — TPFA-nın bunu SƏSSİZCƏ diaqonala "yumşaltması" (scalarize)
    QADAĞANDIR (bax `has_off_diagonal`/`TwoPointFluxDiscretization.build`
    xəbərdarlığı).
    """
    kxx: PropertyMap
    kyy: PropertyMap
    kzz: PropertyMap
    kxy: Optional[PropertyMap] = None
    kxz: Optional[PropertyMap] = None
    kyz: Optional[PropertyMap] = None

    def has_off_diagonal(self, tol: float = 1e-12) -> bool:
        """TPFA-nın DÜZGÜN HƏLL EDƏ BİLMƏDİYİ off-diaqonal komponent varmı
        (yəni real, sıfırdan fərqli anizotropluq bucağı)."""
        for component in (self.kxy, self.kxz, self.kyz):
            if component is not None and np.any(np.abs(component.values) > tol):
                return True
        return False


@dataclass
class RockProperties:
    """Statik süxur xassələri — geoloji modeldən gəlir."""
    porosity: PropertyMap
    permx: PropertyMap
    permy: PropertyMap
    permz: Optional[PropertyMap] = None
    net_to_gross: Optional[PropertyMap] = None
    compressibility: float = 4.5e-5
    #: HƏLƏ HEÇ BİR HƏLLEDİCİ TƏRƏFİNDƏN İSTİFADƏ OLUNMUR (bax
    #: `PermeabilityTensor` docstring-i) — yalnız gələcək MPFA-O üçün
    #: opt-in verilənlər daşıyıcısı. `None` (defolt) — mövcud bütün
    #: modellər ÜÇÜN DAVRANIŞ TAM EYNİDİR.
    permeability_tensor: Optional[PermeabilityTensor] = None

    def validate(self) -> list:
        """Sərt fiziki xətalar. `validate_warnings()` — qeyri-adi (amma
        mümkün) diapazon xəbərdarlıqları üçün, bax Phase 1 hesabatı.

        HƏQİQİ SƏHV (bax tapşırıq: NaN/inf idarəetməsi auditi): əvvəllər
        bu metod `self.porosity.values <= 0` kimi XAM müqayisələr
        işlədirdi — NaN dəyər üçün `NaN <= 0` HƏMİŞƏ `False` qaytarır,
        ona görə NaN/sonsuz PORO/PERMX/PERMY SƏSSİZCƏ bu yoxlamadan
        keçib `ReservoirModel.validate()` (bax `reservoir_model.py`)
        vasitəsilə simulyasiyaya buraxıla bilərdi. İndi eyni fayldan
        onsuz da idxal edilən (`validate_warnings()`-in artıq işlətdiyi)
        `validate_porosity`/`validate_permeability` istifadə olunur —
        bunlar NaN/sonsuzu AÇIQ xəta kimi tuturlar (bax `validation.py`
        `_finite_issue`)."""
        issues = []
        issues += validate_porosity(self.porosity.values, "PORO").errors
        issues += validate_permeability(self.permx.values, "PERMX").errors
        issues += validate_permeability(self.permy.values, "PERMY").errors
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
