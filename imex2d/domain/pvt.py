"""PVT cədvəli — MƏLUMAT strukturu.

Mühərrik bu sinfə birbaşa müraciət etmir; yalnız IPVTProvider
interfeysi vasitəsilə əlçatandır. Cədvəl həm laboratoriya hesabatından,
həm korrelyasiyadan gələ bilər — modelin bundan xəbəri yoxdur.

Vahidlər (METRIC): p [bar], Bo/Bw [rm3/sm3], Rs [sm3/sm3], μ [cP].
"""

from __future__ import annotations
from dataclasses import dataclass

import numpy as np

from .unit_conversions import convert, to_engine_units
from .validation import (check_extrapolation_range, validate_compressibility,
                         validate_viscosity)


@dataclass
class PVTTable:
    """Təzyiqdən asılı flüid xassələri cədvəli.

    Bütün massivlər eyni uzunluqda və `pressure` artan sıralı olmalıdır.
    """
    pressure: np.ndarray
    oil_fvf: np.ndarray
    oil_viscosity: np.ndarray
    solution_gor: np.ndarray
    water_fvf: np.ndarray
    water_viscosity: np.ndarray
    bubble_point: float = 0.0
    rock_compressibility: float = 4.5e-5
    source: str = "manual"

    def __post_init__(self):
        arrays = ("pressure", "oil_fvf", "oil_viscosity", "solution_gor",
                  "water_fvf", "water_viscosity")
        for name in arrays:
            setattr(self, name, np.asarray(getattr(self, name), dtype=float).ravel())

    @property
    def size(self) -> int:
        return int(self.pressure.size)

    @classmethod
    def from_values(cls, pressure, oil_fvf, oil_viscosity, solution_gor,
                    water_fvf, water_viscosity, bubble_point: float = 0.0,
                    rock_compressibility: float = 4.5e-5, source: str = "manual",
                    pressure_unit: str = "bar", viscosity_unit: str = "cP",
                    solution_gor_unit: str = "sm3/sm3",
                    rock_compressibility_unit: str = "bar") -> "PVTTable":
        """Xarici vahiddə (məs. FIELD: psi/cP/scf-stb) verilmiş PVT
        məlumatını mühərrik vahidinə (bar/cP/sm3sm3/1-bar) çevirib
        `PVTTable` qurur — idxal GİRİŞ sərhədidir.

        `oil_fvf`/`water_fvf` (Bo/Bw) ÇEVRİLMİR — bunlar ÖLÇÜSÜZ nisbətdir
        (rezervuar həcmi / səth həcmi, hər iki həcm EYNİ fiziki vahiddə
        ölçüldüyü üçün nisbət vahiddən asılı deyil, bax UNITS.md). Defolt
        parametrlər mühərrik vahidləridir, ona görə heç bir arqument
        DƏYİŞDİRİLMƏSƏ, nəticə birbaşa `PVTTable(...)` çağırışı ilə
        ƏDƏD-ƏDƏD eynidir (defolt vahid == hədəf vahid -> `convert` no-op).
        """
        pressure = to_engine_units(np.asarray(pressure, float), pressure_unit, "pressure")
        oil_viscosity = to_engine_units(np.asarray(oil_viscosity, float),
                                        viscosity_unit, "viscosity")
        water_viscosity = to_engine_units(np.asarray(water_viscosity, float),
                                          viscosity_unit, "viscosity")
        solution_gor = convert(np.asarray(solution_gor, float), solution_gor_unit,
                               "sm3/sm3", "solution_gor")
        bubble_point = to_engine_units(float(bubble_point), pressure_unit, "pressure")
        rock_compressibility = to_engine_units(float(rock_compressibility),
                                               rock_compressibility_unit, "compressibility")
        return cls(pressure=pressure, oil_fvf=np.asarray(oil_fvf, float),
                   oil_viscosity=oil_viscosity, solution_gor=solution_gor,
                   water_fvf=np.asarray(water_fvf, float), water_viscosity=water_viscosity,
                   bubble_point=bubble_point, rock_compressibility=rock_compressibility,
                   source=source)

    def validate(self) -> list:
        issues = []
        n = self.size
        if n < 2:
            issues.append("PVT cədvəlində ən azı iki sətir olmalıdır.")
            return issues
        for name in ("oil_fvf", "oil_viscosity", "solution_gor",
                     "water_fvf", "water_viscosity"):
            if getattr(self, name).size != n:
                issues.append(f"PVT: '{name}' sütununun uzunluğu təzyiq sütunundan fərqlidir.")
        # NaN/sonsuz: `np.diff(...) <= 0`/`<= 0` müqayisələri NaN üçün HƏMİŞƏ
        # False qaytarır (IEEE 754), ona görə NaN aşağıdakı yoxlamalardan
        # SƏSSİZCƏ keçə bilərdi — burada AYRICA, aydın mesajla tutulur.
        for name in ("pressure", "oil_fvf", "oil_viscosity", "solution_gor",
                     "water_fvf", "water_viscosity"):
            column = getattr(self, name)
            if column.size == n and np.any(~np.isfinite(column)):
                issues.append(f"PVT: '{name}' sütununda NaN/sonsuz dəyər var.")
        if np.any(np.diff(self.pressure) <= 0):
            issues.append("PVT: təzyiq sütunu artan sıralı olmalıdır (təkrarlanan dəyər də daxil).")
        if np.any(self.oil_fvf <= 0) or np.any(self.water_fvf <= 0):
            issues.append("PVT: formasiya həcm əmsalı müsbət olmalıdır.")
        if np.any(self.oil_viscosity <= 0) or np.any(self.water_viscosity <= 0):
            issues.append("PVT: lözlük müsbət olmalıdır.")
        if np.any(np.diff(self.solution_gor) < -1e-9):
            issues.append("PVT: Rs təzyiqlə azalmamalıdır.")
        return issues

    def validate_warnings(self) -> list:
        """Rədd edilməyən, amma qeyri-adi diapazon xəbərdarlıqları —
        `validate()`-in ÜSTÜNƏ (onu ƏVƏZ ETMİR). Bax `domain/validation.py`."""
        if self.size < 2:
            return []
        warnings = list(validate_viscosity(self.oil_viscosity, "PVT: neft lözlüyü").warnings)
        warnings += validate_viscosity(self.water_viscosity, "PVT: su lözlüyü").warnings
        warnings += validate_compressibility(self.rock_compressibility,
                                             "PVT: süxur sıxılması").warnings
        return warnings

    def check_query_range(self, pressures) -> list:
        """`pressures` cədvəlin [min, maks] diapazonundan kənara çıxırsa
        bildirir — mövcud `BlackOilPVTProvider` sərhədə KƏSİR (clamp),
        bu ekstrapolyasiyanın SƏSSİZ olmasının qarşısını alır."""
        if self.size < 2:
            return []
        return check_extrapolation_range(pressures, float(self.pressure[0]),
                                         float(self.pressure[-1]), "PVT sorğusu")

    def compressibility(self, key: str) -> np.ndarray:
        """c = -(1/B)·dB/dp — cədvəldən ədədi törəmə ilə.

        Bu, yeni fiziki model deyil, sıxılmanın tərifidir.
        """
        b = getattr(self, key)
        with np.errstate(divide="ignore", invalid="ignore"):
            c = -np.gradient(b, self.pressure) / np.maximum(b, 1e-30)
        return np.nan_to_num(np.abs(c), nan=0.0, posinf=0.0, neginf=0.0)

    def as_rows(self):
        """UI cədvəli üçün sətir-sətir görünüş."""
        return list(zip(self.pressure, self.oil_fvf, self.oil_viscosity,
                        self.solution_gor, self.water_fvf, self.water_viscosity))
