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


@dataclass
class GasCoreyParameters:
    """Qaz-neft SCAL (A7) — su-neft `CoreyParameters`-in yoldaşıdır.

    Su-neft sistemində `Sw` artdıqca `kro` azalır; qaz-neft sistemində
    isə `Sg` artdıqca eyni şəkildə azalır — buna görə `krog` (qaz-neft
    sistemində neftin keçiriciliyi) `CoreyParameters.kro`-nun `Sg`
    versiyasıdır. `kro_end` burada TƏKRARLANMIR — Stone II düsturunda
    su-neft əyrisinin `kro_end`-i işlədilir, çünki hər iki əyri eyni
    fiziki son nöqtəyə (Swc-də, Sg=0-da) istinad edir. Fərqli dəyər
    versək iki əyri uyğunsuz olardı.
    """
    sgc: float = 0.05        # bağlı (hərəkətsiz) qaz doyumluluğu
    sorg: float = 0.10       # qaza qarşı qalıq neft doyumluluğu
    krg_end: float = 0.80
    ng: float = 2.0          # qazın Corey göstəricisi
    nog: float = 2.0         # qaz-neft sistemində neftin Corey göstəricisi

    def normalized(self, sg, swc: float) -> np.ndarray:
        den = max(1.0 - swc - self.sgc - self.sorg, 1e-8)
        return np.clip((np.asarray(sg, float) - self.sgc) / den, 0.0, 1.0)

    def krg(self, sg, swc: float) -> np.ndarray:
        return self.krg_end * self.normalized(sg, swc) ** self.ng

    def krog(self, sg, swc: float, kro_end: float) -> np.ndarray:
        """Qaz-neft sistemində neftin keçiriciliyi — Sw = Swc-də."""
        return kro_end * (1.0 - self.normalized(sg, swc)) ** self.nog

    def krg_derivative(self, sg, swc: float) -> np.ndarray:
        """dkrg/dSg — analitik (su-neft `krw_derivative` ilə eyni üslub)."""
        span = max(1.0 - swc - self.sgc - self.sorg, 1e-8)
        sn = self.normalized(sg, swc)
        raw = (np.asarray(sg, float) - self.sgc) / span
        inside = (raw > 0.0) & (raw < 1.0)
        result = np.zeros_like(sn)
        clipped = np.clip(sn, 1e-12, 1.0)
        result[inside] = (self.krg_end * self.ng * clipped[inside] ** (self.ng - 1.0)
                          / span)
        return result

    def krog_derivative(self, sg, swc: float, kro_end: float) -> np.ndarray:
        """dkrog/dSg — analitik (mənfi işarəli)."""
        span = max(1.0 - swc - self.sgc - self.sorg, 1e-8)
        sn = self.normalized(sg, swc)
        raw = (np.asarray(sg, float) - self.sgc) / span
        inside = (raw > 0.0) & (raw < 1.0)
        result = np.zeros_like(sn)
        clipped = np.clip(1.0 - sn, 1e-12, 1.0)
        result[inside] = (-kro_end * self.nog * clipped[inside] ** (self.nog - 1.0)
                          / span)
        return result

    def validate(self, swc: float) -> list:
        issues = []
        if swc + self.sgc + self.sorg >= 0.95:
            issues.append("Swc + Sgc + Sorg ≥ 0.95 — hərəkətli qaz "
                          "doyumluluğu intervalı yoxdur.")
        if not (0 < self.krg_end <= 1):
            issues.append("krg_end (0, 1] intervalında olmalıdır.")
        if self.sgc < 0 or self.sorg < 0:
            issues.append("Sgc və Sorg mənfi ola bilməz.")
        return issues
