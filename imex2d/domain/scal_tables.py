"""SCAL cədvəlləri — laboratoriya ölçmələrindən nisbi keçiricilik.

Corey düsturu analitik və hamardır, lakin real kern məlumatı belə
davranmır: əyrilər asimmetrik olur, son nöqtələr ayrıca ölçülür,
bəzən orta hissədə əyilmə görünür. Laboratoriya cədvəli birbaşa
işlədilməlidir.

Format Eclipse `SWOF` ilə eynidir:

    Sw      krw       kro       Pc
    0.20    0.000     0.800     0.00
    0.30    0.015     0.520     0.00
    ...

REGION ANLAYIŞI
Bir yataqda litologiya dəyişir: qumdaşı, əhəngdaşı, gilli zona. Hər
birinin öz `kr` əyrisi var. `SATNUM` massivi hər hüceyrəni bir regiona
bağlayır (GRDECL-dən oxunur, bax `ECLIPSE_IO.md`).

MONOTONLUQ tələb olunur: `krw` artan, `kro` azalan. Pozulsa
diskretizasiya qeyri-stabil olur — ona görə yoxlama sərtdir.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np


@dataclass
class SaturationTable:
    """Bir region üçün su-neft nisbi keçiricilik cədvəli."""
    sw: np.ndarray
    krw: np.ndarray
    kro: np.ndarray
    pc: Optional[np.ndarray] = None
    name: str = ""

    def __post_init__(self):
        self.sw = np.asarray(self.sw, dtype=float).ravel()
        self.krw = np.asarray(self.krw, dtype=float).ravel()
        self.kro = np.asarray(self.kro, dtype=float).ravel()
        if self.pc is not None:
            self.pc = np.asarray(self.pc, dtype=float).ravel()

    # ─────────────────────────────────────────── son nöqtələr
    @property
    def swc(self) -> float:
        """Bağlı su: krw sıfırdan çıxan ilk nöqtə."""
        moving = np.nonzero(self.krw > 0.0)[0]
        return float(self.sw[0] if moving.size == 0
                     else self.sw[max(moving[0] - 1, 0)])

    @property
    def sor(self) -> float:
        """Qalıq neft: kro sıfıra çatan nöqtədən sonrası."""
        moving = np.nonzero(self.kro > 0.0)[0]
        if moving.size == 0:
            return float(1.0 - self.sw[-1])
        index = min(moving[-1] + 1, self.sw.size - 1)
        return float(1.0 - self.sw[index])

    @property
    def krw_end(self) -> float:
        return float(self.krw.max())

    @property
    def kro_end(self) -> float:
        return float(self.kro.max())

    @property
    def has_capillary(self) -> bool:
        return self.pc is not None and bool(np.any(np.abs(self.pc) > 1e-12))

    # ─────────────────────────────────────────── interpolyasiya
    def interpolate_krw(self, sw) -> np.ndarray:
        return np.interp(sw, self.sw, self.krw)

    def interpolate_kro(self, sw) -> np.ndarray:
        return np.interp(sw, self.sw, self.kro)

    def interpolate_pc(self, sw) -> np.ndarray:
        sw = np.asarray(sw, float)
        if self.pc is None:
            return np.zeros_like(sw)
        return np.interp(sw, self.sw, self.pc)

    def slope(self, values: np.ndarray, sw) -> np.ndarray:
        """Parçalı xətti cədvəlin DƏQİQ törəməsi — interval meyli.

        Hamar törəmə (`np.gradient`) daha "gözəl" görünür, lakin
        cədvəlin özü ilə uyğun gəlmir və Nyuton iterasiyasında
        Jakobianı sonlu fərqdən uzaqlaşdırır (bax `A6_PLAN.md`).
        """
        sw = np.atleast_1d(np.asarray(sw, float))
        if self.sw.size < 2:
            return np.zeros_like(sw)
        slopes = np.diff(values) / np.diff(self.sw)
        index = np.clip(np.searchsorted(self.sw, sw, side="right") - 1,
                        0, self.sw.size - 2)
        result = slopes[index]
        outside = (sw < self.sw[0]) | (sw > self.sw[-1])
        return np.where(outside, 0.0, result)

    # ─────────────────────────────────────────── yoxlama
    def validate(self) -> List[str]:
        label = self.name or "cədvəl"
        issues = []
        if not (self.sw.size == self.krw.size == self.kro.size):
            issues.append(f"{label}: sütun uzunluqları fərqlidir.")
            return issues
        if self.sw.size < 2:
            issues.append(f"{label}: ən azı iki sətir lazımdır.")
            return issues
        if self.pc is not None and self.pc.size != self.sw.size:
            issues.append(f"{label}: Pc sütununun uzunluğu uyğun deyil.")
        if np.any(np.diff(self.sw) <= 0):
            issues.append(f"{label}: Sw artan sıralı olmalıdır.")
        if self.sw[0] < -1e-9 or self.sw[-1] > 1.0 + 1e-9:
            issues.append(f"{label}: Sw [0, 1] intervalından kənardadır.")
        if np.any(self.krw < -1e-12) or np.any(self.kro < -1e-12):
            issues.append(f"{label}: nisbi keçiricilik mənfi ola bilməz.")
        if np.any(self.krw > 1.0 + 1e-9) or np.any(self.kro > 1.0 + 1e-9):
            issues.append(f"{label}: nisbi keçiricilik 1-dən böyükdür.")
        if np.any(np.diff(self.krw) < -1e-9):
            issues.append(f"{label}: krw azalır — monoton artan olmalıdır.")
        if np.any(np.diff(self.kro) > 1e-9):
            issues.append(f"{label}: kro artır — monoton azalan olmalıdır.")
        if self.swc >= 1.0 - self.sor:
            issues.append(f"{label}: hərəkətli doyumluluq intervalı boşdur.")
        return issues

    def summary(self) -> str:
        return (f"{self.name or 'SCAL'}: {self.sw.size} sətir, "
                f"Swc {self.swc:.3f}, Sor {self.sor:.3f}, "
                f"krw_end {self.krw_end:.3f}, kro_end {self.kro_end:.3f}"
                f"{', Pc var' if self.has_capillary else ''}")

    @classmethod
    def from_corey(cls, parameters, points: int = 21,
                   name: str = "Corey") -> "SaturationTable":
        """Corey parametrlərindən cədvəl — köhnə modellərin körpüsü."""
        sw = np.linspace(parameters.swc, 1.0 - parameters.sor, points)
        return cls(sw=sw,
                   krw=np.asarray(parameters.krw(sw), float),
                   kro=np.asarray(parameters.kro(sw), float),
                   name=name)


@dataclass
class SaturationTableSet:
    """Region -> cədvəl uyğunluğu.

    Region nömrələri `SATNUM` ilə eynidir (1-dən başlayır). Hüceyrənin
    regionu tapılmasa, defolt cədvəl işlədilir — susmaq təhlükəlidir,
    ona görə bu hal ayrıca sayılır.
    """
    tables: Dict[int, SaturationTable] = field(default_factory=dict)
    default_region: int = 1

    def __len__(self) -> int:
        return len(self.tables)

    @property
    def regions(self) -> List[int]:
        return sorted(self.tables)

    def get(self, region: Optional[int] = None) -> SaturationTable:
        if region is not None and region in self.tables:
            return self.tables[region]
        if self.default_region in self.tables:
            return self.tables[self.default_region]
        if not self.tables:
            raise ValueError("SCAL cədvəli yoxdur.")
        return self.tables[self.regions[0]]

    def add(self, region: int, table: SaturationTable) -> None:
        self.tables[int(region)] = table

    def validate(self) -> List[str]:
        if not self.tables:
            return ["Heç bir SCAL cədvəli yüklənməyib."]
        issues = []
        for region, table in sorted(self.tables.items()):
            for message in table.validate():
                issues.append(f"Region {region}: {message}")
        return issues

    def summary(self) -> str:
        return "\n".join(f"  Region {region}: {table.summary()}"
                         for region, table in sorted(self.tables.items()))
