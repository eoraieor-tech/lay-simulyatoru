"""PVT cədvəli — MƏLUMAT strukturu.

Mühərrik bu sinfə birbaşa müraciət etmir; yalnız IPVTProvider
interfeysi vasitəsilə əlçatandır. Cədvəl həm laboratoriya hesabatından,
həm korrelyasiyadan gələ bilər — modelin bundan xəbəri yoxdur.

Vahidlər (METRIC): p [bar], Bo/Bw [rm3/sm3], Rs [sm3/sm3], μ [cP].
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

import numpy as np


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
    gas_fvf: Optional[np.ndarray] = None
    gas_viscosity: Optional[np.ndarray] = None
    """Bg(p), μg(p) — YALNIZ üç fazalı (A7) modeldə tələb olunur.

    `None` — köhnə iki fazalı cədvəllər üçün defolt; belə cədvəl
    `has_gas_phase` == False qaytarır və mühərrik onu əvvəlki kimi
    iki fazalı (su-neft) işlədir. Geriyə uyğunluq üçün heç bir mövcud
    layihə faylı və ya test pozulmur.
    """

    def __post_init__(self):
        arrays = ("pressure", "oil_fvf", "oil_viscosity", "solution_gor",
                  "water_fvf", "water_viscosity")
        for name in arrays:
            setattr(self, name, np.asarray(getattr(self, name), dtype=float).ravel())
        for name in ("gas_fvf", "gas_viscosity"):
            value = getattr(self, name)
            if value is not None:
                setattr(self, name, np.asarray(value, dtype=float).ravel())

    @property
    def has_gas_phase(self) -> bool:
        return self.gas_fvf is not None and self.gas_viscosity is not None

    @property
    def size(self) -> int:
        return int(self.pressure.size)

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
        if np.any(np.diff(self.pressure) <= 0):
            issues.append("PVT: təzyiq sütunu artan sıralı olmalıdır.")
        if np.any(self.oil_fvf <= 0) or np.any(self.water_fvf <= 0):
            issues.append("PVT: formasiya həcm əmsalı müsbət olmalıdır.")
        if np.any(self.oil_viscosity <= 0) or np.any(self.water_viscosity <= 0):
            issues.append("PVT: lözlük müsbət olmalıdır.")
        if np.any(np.diff(self.solution_gor) < -1e-9):
            issues.append("PVT: Rs təzyiqlə azalmamalıdır.")
        if self.has_gas_phase:
            for name in ("gas_fvf", "gas_viscosity"):
                values = getattr(self, name)
                if values.size != n:
                    issues.append(f"PVT: '{name}' sütununun uzunluğu "
                                  f"təzyiq sütunundan fərqlidir.")
            if np.any(self.gas_fvf <= 0):
                issues.append("PVT: Bg müsbət olmalıdır.")
            if np.any(self.gas_viscosity <= 0):
                issues.append("PVT: μg müsbət olmalıdır.")
            if np.any(np.diff(self.gas_fvf) > 1e-9):
                issues.append("PVT: Bg təzyiqlə artmamalıdır (qaz sıxılır).")
        return issues

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
