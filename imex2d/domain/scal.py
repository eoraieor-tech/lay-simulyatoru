"""SCAL parametrləri — MƏLUMAT strukturu.

Diqqət: buradakı Corey düsturları YENİ deyil, mövcud koddan köçürülüb.
Hesablama mühərriki bu sinfə birbaşa müraciət ETMİR — yalnız
IRelativePermeabilityProvider interfeysi vasitəsilə əlçatandır.
"""

from __future__ import annotations
from dataclasses import dataclass

import numpy as np


@dataclass
class CapillaryParameters:
    """Brooks-Corey kapilyar təzyiq modeli.

        Pc(Sw) = Pe · Sw_norm^(-1/λ)

    `entry_pressure = 0` verildikdə model söndürülür və provider
    qurulmur — yəni köhnə (kapilyarsız) davranış qalır.
    """
    entry_pressure: float = 0.0     # Pe, bar
    lambda_exponent: float = 2.0    # λ — məsamə ölçüsü paylanması indeksi
    max_pressure: float = 5.0       # Pc-nin yuxarı kəsilmə həddi, bar

    @property
    def enabled(self) -> bool:
        return self.entry_pressure > 0.0

    def validate(self) -> list:
        issues = []
        if self.entry_pressure < 0.0:
            issues.append("Kapilyar giriş təzyiqi mənfi ola bilməz.")
        if self.lambda_exponent <= 0.0:
            issues.append("Brooks-Corey λ müsbət olmalıdır.")
        if self.enabled and self.max_pressure <= self.entry_pressure:
            issues.append("Pc yuxarı həddi giriş təzyiqindən böyük olmalıdır.")
        return issues


@dataclass
class CoreyParameters:
    swc: float = 0.20
    sor: float = 0.25
    krw_end: float = 0.35
    kro_end: float = 0.90
    nw: float = 2.5
    no: float = 2.0

    def normalized(self, sw) -> np.ndarray:
        den = max(1.0 - self.swc - self.sor, 1e-8)
        return np.clip((np.asarray(sw, float) - self.swc) / den, 0.0, 1.0)

    def krw(self, sw) -> np.ndarray:
        return self.krw_end * self.normalized(sw) ** self.nw

    def kro(self, sw) -> np.ndarray:
        return self.kro_end * (1.0 - self.normalized(sw)) ** self.no

    def validate(self) -> list:
        issues = []
        if self.swc + self.sor >= 0.95:
            issues.append("Swc + Sor ≥ 0.95 — hərəkətli doyumluluq intervalı yoxdur.")
        if not (0 < self.krw_end <= 1) or not (0 < self.kro_end <= 1):
            issues.append("Nisbi keçiricilik son nöqtələri (0, 1] intervalında olmalıdır.")
        return issues
