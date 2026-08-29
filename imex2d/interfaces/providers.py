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
        """Pb — doyma təzyiqi. Qaz modelləşdirilmir, lakin diaqnostika
        (quyudibi təzyiqin Pb-dən aşağı düşməsi) bunu işlədir."""
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


class ICapillaryPressureProvider(ABC):
    """Kapilyar təzyiq. HƏLƏ İMPLEMENTASİYA EDİLMƏYİB."""

    @abstractmethod
    def pcow(self, sw: np.ndarray, region: Optional[np.ndarray] = None) -> np.ndarray:
        """Pc = Po - Pw, bar."""

    @abstractmethod
    def dpcow_dsw(self, sw: np.ndarray, region: Optional[np.ndarray] = None) -> np.ndarray:
        """Törəmə — implicit sxemdə Jakobian üçün lazımdır."""


class InitialState:
    """Initialization provider-in nəticəsi — iki fazalı (So = 1 − Sw)."""

    def __init__(self, pressure: np.ndarray, water_saturation: np.ndarray):
        self.pressure = pressure
        self.water_saturation = water_saturation


class IInitializationProvider(ABC):
    """Equilibration. HƏLƏ İMPLEMENTASİYA EDİLMƏYİB."""

    @abstractmethod
    def initialize(self, model) -> InitialState:
        """ReservoirModel-dən ilkin təzyiq və doyumluluq sahələri qurur."""
