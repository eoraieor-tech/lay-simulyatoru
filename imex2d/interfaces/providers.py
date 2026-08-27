"""Gələcək modullar üçün İNTERFEYSLƏR.

Bu fayl YALNIZ müqavilə (contract) təyin edir. Heç bir implementasiya
yoxdur və bu qəsdəndir — tapşırıq PVT, kapilyar təzyiq və
initialization modullarının yazılmasını qadağan edir.

Mühərrik bu interfeyslərdən asılıdır, konkret sinifdən yox
(Dependency Inversion Principle). Provider verilmədikdə mühərrik
modelin statik dəyərləri ilə işləyir — yəni hazırkı davranış qorunur.
"""

from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Optional

import numpy as np


class IPVTProvider(ABC):
    """Təzyiqdən asılı flüid xassələri. HƏLƏ İMPLEMENTASİYA EDİLMƏYİB."""

    @abstractmethod
    def oil_fvf(self, pressure: np.ndarray, region: Optional[np.ndarray] = None) -> np.ndarray:
        """Bo(p) — neftin formasiya həcm əmsalı."""

    @abstractmethod
    def oil_viscosity(self, pressure: np.ndarray, region: Optional[np.ndarray] = None) -> np.ndarray:
        """μo(p), cP."""

    @abstractmethod
    def water_fvf(self, pressure: np.ndarray, region: Optional[np.ndarray] = None) -> np.ndarray:
        """Bw(p)."""

    @abstractmethod
    def water_viscosity(self, pressure: np.ndarray, region: Optional[np.ndarray] = None) -> np.ndarray:
        """μw(p), cP."""

    @abstractmethod
    def total_compressibility(self, pressure: np.ndarray, sw: np.ndarray,
                              region: Optional[np.ndarray] = None) -> np.ndarray:
        """ct(p, Sw), 1/bar."""

    def bubble_point(self, region: Optional[np.ndarray] = None):
        """Pb — qaz fazası əlavə olunanda tələb olunacaq."""
        raise NotImplementedError

    def has_gas_phase(self, region: Optional[np.ndarray] = None) -> bool:
        """A7: cədvəldə Bg/μg varsa üç fazalı hesablama aktivləşir.

        Defolt `False` — köhnə iki fazalı provider-lər (məs. sabit
        dəyərli `StaticFluidPropertiesProvider`) heç nə etmədən bu
        interfeysi tətbiq edir.
        """
        return False

    def gas_fvf(self, pressure: np.ndarray,
               region: Optional[np.ndarray] = None) -> np.ndarray:
        """Bg(p) — YALNIZ `has_gas_phase()` True olanda çağırılır."""
        raise NotImplementedError

    def gas_viscosity(self, pressure: np.ndarray,
                      region: Optional[np.ndarray] = None) -> np.ndarray:
        """μg(p), cP — YALNIZ `has_gas_phase()` True olanda çağırılır."""
        raise NotImplementedError

    def solution_gor(self, pressure: np.ndarray,
                     region: Optional[np.ndarray] = None) -> np.ndarray:
        """Rs(p) — doyma təzyiqinə qədər artan, ondan sonra sabit.

        Bu, mühərrikin doyma nöqtəsini keçən hüceyrələrdə sərbəst qaz
        fazasının nə vaxt yarandığını təyin etməsi üçün tələb olunur
        (dəyişən keçid — bax `A7_PLAN.md`).
        """
        raise NotImplementedError


class IRelativePermeabilityProvider(ABC):
    """Nisbi keçiricilik. Region-a görə fərqli cədvəl qaytara bilər."""

    @abstractmethod
    def krw(self, sw: np.ndarray, region: Optional[np.ndarray] = None) -> np.ndarray: ...

    @abstractmethod
    def kro(self, sw: np.ndarray, region: Optional[np.ndarray] = None) -> np.ndarray: ...

    @abstractmethod
    def saturation_limits(self, region: Optional[int] = None) -> tuple:
        """(Swc, 1 - Sor) — doyumluluğun kəsilmə hədləri."""

    @abstractmethod
    def endpoint_water_mobility(self, water_viscosity: float,
                                region: Optional[int] = None) -> float:
        """Vurucu quyu üçün krw_end / μw."""

    @abstractmethod
    def max_fractional_flow_derivative(self, water_viscosity: float,
                                       oil_viscosity: float,
                                       region: Optional[int] = None) -> float:
        """max |dfw/dSw| — CFL limiti üçün."""

    def has_gas_phase(self, region: Optional[int] = None) -> bool:
        """A7: üç fazalı (Stone) provider `True` qaytarır.

        Defolt `False` — mövcud iki fazalı provider-lər (Corey
        adapteri, B4 cədvəl provider-i) heç nə etmədən bu interfeysi
        tətbiq edir.
        """
        return False

    def krg(self, sg: np.ndarray, region: Optional[np.ndarray] = None) -> np.ndarray:
        """krg(Sg) — YALNIZ `has_gas_phase()` True olanda çağırılır."""
        raise NotImplementedError

    def kro_three_phase(self, sw: np.ndarray, sg: np.ndarray,
                        region: Optional[np.ndarray] = None) -> np.ndarray:
        """kro(Sw, Sg) — Stone modeli, hər iki doyumluluqdan asılı.

        İki fazalı `kro(sw)`-dan fərqli olaraq üç fazalı sistemdə neft
        həm suyun, həm qazın "sıxışdırmasına" məruz qalır.
        """
        raise NotImplementedError

    def gas_saturation_limits(self, region: Optional[int] = None) -> tuple:
        """(Sgc, 1 − Swc − Sorg) — qazın hərəkətli doyumluluq intervalı."""
        raise NotImplementedError


class ICapillaryPressureProvider(ABC):
    """Kapilyar təzyiq. HƏLƏ İMPLEMENTASİYA EDİLMƏYİB."""

    @abstractmethod
    def pcow(self, sw: np.ndarray, region: Optional[np.ndarray] = None) -> np.ndarray:
        """Pc = Po - Pw, bar."""

    @abstractmethod
    def dpcow_dsw(self, sw: np.ndarray, region: Optional[np.ndarray] = None) -> np.ndarray:
        """Törəmə — implicit sxemdə Jakobian üçün lazımdır."""


class InitialState:
    """Initialization provider-in nəticəsi.

    `gas_saturation` YALNIZ üç fazalı modeldə (A7) verilir. `None`
    defoltu iki fazalı davranışı saxlayır — köhnə provider-lər və
    onları çağıran bütün kod (mühərriklər, testlər) toxunulmadan
    işləyir; `So = 1 − Sw` fərziyyəsi qorunur.
    """

    def __init__(self, pressure: np.ndarray, water_saturation: np.ndarray,
                 gas_saturation: Optional[np.ndarray] = None):
        self.pressure = pressure
        self.water_saturation = water_saturation
        self.gas_saturation = (None if gas_saturation is None
                               else np.asarray(gas_saturation, float))

    @property
    def has_gas(self) -> bool:
        return self.gas_saturation is not None


class IInitializationProvider(ABC):
    """Equilibration. HƏLƏ İMPLEMENTASİYA EDİLMƏYİB."""

    @abstractmethod
    def initialize(self, model) -> InitialState:
        """ReservoirModel-dən ilkin təzyiq və doyumluluq sahələri qurur."""
