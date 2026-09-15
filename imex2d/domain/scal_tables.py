"""SCAL cədvəlləri — laboratoriya ölçmələrindən nisbi keçiricilik.

Corey düsturu analitik və hamardır, lakin real kern məlumatı belə
davranmır: əyrilər asimmetrik olur, son nöqtələr ayrıca ölçülür,
bəzən orta hissədə əyilmə görünür. Laboratoriya cədvəli birbaşa
işlədilməlidir.

Su-neft üçün format Eclipse `SWOF`, qaz-neft üçün `SGOF` ilə eynidir
(`GasSaturationTable`, G4 — SPE1 etalonu bu cədvəlləri tələb edir):

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

from .validation import check_extrapolation_range


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
        # NaN/sonsuz: `<=`/`np.diff` müqayisələri NaN üçün HƏMİŞƏ False
        # qaytarır — aşağıdakı monotonluq/həd yoxlamalarından SƏSSİZCƏ
        # keçməsin deyə burada AYRICA, aydın mesajla tutulur.
        for name in ("sw", "krw", "kro"):
            column = getattr(self, name)
            if np.any(~np.isfinite(column)):
                issues.append(f"{label}: '{name}' sütununda NaN/sonsuz dəyər var.")
        if self.pc is not None and self.pc.size == self.sw.size and np.any(~np.isfinite(self.pc)):
            issues.append(f"{label}: 'pc' sütununda NaN/sonsuz dəyər var.")
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

    def check_query_range(self, sw_values) -> List[str]:
        """`sw_values` cədvəlin ölçdüyü [min Sw, maks Sw] intervalından
        kənara çıxırsa bildirir — `interpolate_krw`/`kro`/`pc` (`np.interp`)
        sərhədə KƏSİR (clamp), bu ekstrapolyasiyanın SƏSSİZ olmasının
        qarşısını alır (bax `domain/validation.check_extrapolation_range`)."""
        if self.sw.size < 2:
            return []
        return check_extrapolation_range(sw_values, float(self.sw[0]), float(self.sw[-1]),
                                         f"{self.name or 'SCAL'}: Sw sorğusu")

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


@dataclass
class GasSaturationTable:
    """Bir region üçün QAZ-NEFT nisbi keçiricilik cədvəli — G4.

    Eclipse `SGOF` sütunları: `Sg`, `krg`, `krog`, `Pcog`.

    NİYƏ AYRI SİNİF: su-neft cədvəlində arqument `Sw`-dir və `kro`
    AZALIR; burada arqument `Sg`-dir və `krg` ARTIR, `krog` AZALIR.
    Eyni sinfə sığışdırmaq sütun adlarını yalan edərdi.

    MÜQAVİLƏ: bu sinif `GasCoreyParameters`-in YERİNƏ keçir —
    `StoneRelativePermeabilityProvider` yalnız `krg(sg, swc)`,
    `krog(sg, swc, kro_end)`, onların törəmələri, `sgc`/`sorg`/`krg_end`
    və `validate(swc)` çağırır. `swc`/`kro_end` arqumentləri QƏBUL
    EDİLİR, lakin İŞLƏDİLMİR: cədvəl onsuz da Swc-də ölçülüb və öz
    son nöqtəsini daşıyır (Corey düsturunda isə onlar hesaba girir).
    """

    sg: np.ndarray
    krg: np.ndarray
    krog: np.ndarray
    pcog: Optional[np.ndarray] = None
    name: str = ""

    def __post_init__(self):
        self.sg = np.asarray(self.sg, dtype=float).ravel()
        self.krg = np.asarray(self.krg, dtype=float).ravel()
        self.krog = np.asarray(self.krog, dtype=float).ravel()
        if self.pcog is not None:
            self.pcog = np.asarray(self.pcog, dtype=float).ravel()

    # ─────────────────────────────────────────── son nöqtələr
    @property
    def sgc(self) -> float:
        """Bağlı (hərəkətsiz) qaz: krg sıfırdan çıxan son nöqtədən əvvəlki."""
        moving = np.nonzero(self.krg > 0.0)[0]
        return float(self.sg[0] if moving.size == 0
                     else self.sg[max(moving[0] - 1, 0)])

    @property
    def sorg(self) -> float:
        """Qaza qarşı qalıq neft: krog sıfıra çatandan sonrası.

        `Sg + So + Swc = 1` olduğundan cədvəldəki ən böyük axan `Sg`-dən
        çıxarılır. Cədvəldə `Swc` açıq yazılmır, ona görə bu, ancaq
        `Sorg`-un cədvəldən görünən hissəsidir — Stone düsturunda
        yalnız `gas_saturation_limits` üçün işlədilir.
        """
        moving = np.nonzero(self.krog > 0.0)[0]
        if moving.size == 0:
            return 0.0
        index = min(moving[-1] + 1, self.sg.size - 1)
        return float(max(1.0 - self.sg[index], 0.0))

    @property
    def krg_end(self) -> float:
        return float(self.krg.max())

    @property
    def krog_end(self) -> float:
        return float(self.krog.max())

    @property
    def has_capillary(self) -> bool:
        return self.pcog is not None and bool(np.any(np.abs(self.pcog) > 1e-12))

    # ─────────────────────────────────────────── interpolyasiya
    def interpolate_krg(self, sg) -> np.ndarray:
        return np.interp(sg, self.sg, self.krg)

    def interpolate_krog(self, sg) -> np.ndarray:
        return np.interp(sg, self.sg, self.krog)

    def interpolate_pcog(self, sg) -> np.ndarray:
        sg = np.asarray(sg, float)
        if self.pcog is None:
            return np.zeros_like(sg)
        return np.interp(sg, self.sg, self.pcog)

    def slope(self, values: np.ndarray, sg) -> np.ndarray:
        """Parçalı xətti cədvəlin DƏQİQ törəməsi — `SaturationTable.slope`
        ilə eyni qayda (Jakobian qalıqla uyğun qalsın deyə)."""
        sg = np.atleast_1d(np.asarray(sg, float))
        if self.sg.size < 2:
            return np.zeros_like(sg)
        slopes = np.diff(values) / np.diff(self.sg)
        index = np.clip(np.searchsorted(self.sg, sg, side="right") - 1,
                        0, self.sg.size - 2)
        result = slopes[index]
        outside = (sg < self.sg[0]) | (sg > self.sg[-1])
        return np.where(outside, 0.0, result)

    # ──────────────────── `GasCoreyParameters` müqaviləsi (duck-typing)
    def krg_at(self, sg, swc=None) -> np.ndarray:
        return self.interpolate_krg(sg)

    def krog_at(self, sg, swc=None, kro_end=None) -> np.ndarray:
        return self.interpolate_krog(sg)

    # ─────────────────────────────────────────── yoxlama
    def validate(self, swc: float = 0.0) -> List[str]:
        label = self.name or "qaz cədvəli"
        issues = []
        if not (self.sg.size == self.krg.size == self.krog.size):
            issues.append(f"{label}: sütun uzunluqları fərqlidir.")
            return issues
        if self.sg.size < 2:
            issues.append(f"{label}: ən azı iki sətir lazımdır.")
            return issues
        if self.pcog is not None and self.pcog.size != self.sg.size:
            issues.append(f"{label}: Pcog sütununun uzunluğu uyğun deyil.")
        # NaN/sonsuz AYRICA tutulur — `np.diff` müqayisələri NaN üçün
        # həmişə False qaytarır və səhv SƏSSİZCƏ keçərdi (su-neft
        # cədvəlindəki eyni dərs).
        for name in ("sg", "krg", "krog"):
            column = getattr(self, name)
            if np.any(~np.isfinite(column)):
                issues.append(f"{label}: '{name}' sütununda NaN/sonsuz dəyər var.")
        if self.pcog is not None and self.pcog.size == self.sg.size \
                and np.any(~np.isfinite(self.pcog)):
            issues.append(f"{label}: 'pcog' sütununda NaN/sonsuz dəyər var.")
        if np.any(np.diff(self.sg) <= 0):
            issues.append(f"{label}: Sg artan sıralı olmalıdır.")
        if self.sg[0] < -1e-9 or self.sg[-1] > 1.0 + 1e-9:
            issues.append(f"{label}: Sg [0, 1] intervalından kənardadır.")
        if np.any(self.krg < -1e-12) or np.any(self.krog < -1e-12):
            issues.append(f"{label}: nisbi keçiricilik mənfi ola bilməz.")
        if np.any(self.krg > 1.0 + 1e-9) or np.any(self.krog > 1.0 + 1e-9):
            issues.append(f"{label}: nisbi keçiricilik 1-dən böyükdür.")
        if np.any(np.diff(self.krg) < -1e-9):
            issues.append(f"{label}: krg azalır — Sg artdıqca artmalıdır.")
        if np.any(np.diff(self.krog) > 1e-9):
            issues.append(f"{label}: krog artır — Sg artdıqca azalmalıdır.")
        return issues

    def check_query_range(self, sg_values) -> List[str]:
        if self.sg.size < 2:
            return []
        return check_extrapolation_range(
            sg_values, float(self.sg[0]), float(self.sg[-1]),
            f"{self.name or 'SGOF'}: Sg sorğusu")

    def summary(self) -> str:
        return (f"{self.name or 'SGOF'}: {self.sg.size} sətir, "
                f"Sgc {self.sgc:.3f}, Sorg {self.sorg:.3f}, "
                f"krg_end {self.krg_end:.3f}, krog_end {self.krog_end:.3f}"
                f"{', Pcog var' if self.has_capillary else ''}")

    @classmethod
    def from_corey(cls, parameters, swc: float, kro_end: float,
                   points: int = 21, name: str = "Corey (qaz)"):
        """`GasCoreyParameters`-dən cədvəl — köhnə modellərin körpüsü."""
        sg = np.linspace(parameters.sgc, 1.0 - swc - parameters.sorg, points)
        return cls(sg=sg,
                   krg=np.asarray(parameters.krg(sg, swc), float),
                   krog=np.asarray(parameters.krog(sg, swc, kro_end), float),
                   name=name)


@dataclass
class GasSaturationTableSet:
    """Region → qaz-neft cədvəli (su-neft `SaturationTableSet`-in yoldaşı).

    ⏳ REGION MƏHDUDİYYƏTİ: `StoneRelativePermeabilityProvider` qaz
    əyrisini region arqumenti OLMADAN çağırır (`self.gas.krg(sg, swc)`),
    ona görə hazırda yalnız DEFOLT region işlədilir. Regionlu qaz
    əyriləri Stone müqaviləsinin genişləndirilməsini tələb edir.
    """

    tables: Dict[int, GasSaturationTable] = field(default_factory=dict)
    default_region: int = 1

    def __len__(self) -> int:
        return len(self.tables)

    @property
    def regions(self) -> List[int]:
        return sorted(self.tables)

    def get(self, region: Optional[int] = None) -> GasSaturationTable:
        if region is not None and region in self.tables:
            return self.tables[region]
        if self.default_region in self.tables:
            return self.tables[self.default_region]
        if not self.tables:
            raise ValueError("Qaz-neft SCAL cədvəli yoxdur.")
        return self.tables[self.regions[0]]

    def add(self, region: int, table: GasSaturationTable) -> None:
        self.tables[int(region)] = table

    def validate(self, swc: float = 0.0) -> List[str]:
        if not self.tables:
            return ["Heç bir qaz-neft SCAL cədvəli yüklənməyib."]
        issues = []
        for region, table in sorted(self.tables.items()):
            for message in table.validate(swc):
                issues.append(f"Region {region}: {message}")
        return issues

    def summary(self) -> str:
        return "\n".join(f"  Region {region}: {table.summary()}"
                         for region, table in sorted(self.tables.items()))

    # ──────────────────── `GasCoreyParameters` müqaviləsi (duck-typing)
    #
    # `StoneRelativePermeabilityProvider` bu obyekti Corey parametrləri
    # kimi işlədir. `swc`/`kro_end` qəbul edilir, lakin cədvəl onları
    # işlətmir (bax `GasSaturationTable` sənədi).
    @property
    def sgc(self) -> float:
        return self.get().sgc

    @property
    def sorg(self) -> float:
        return self.get().sorg

    @property
    def krg_end(self) -> float:
        return self.get().krg_end

    def krg(self, sg, swc: float = 0.0) -> np.ndarray:
        return self.get().interpolate_krg(sg)

    def krog(self, sg, swc: float = 0.0, kro_end: float = 1.0) -> np.ndarray:
        return self.get().interpolate_krog(sg)

    def krg_derivative(self, sg, swc: float = 0.0) -> np.ndarray:
        table = self.get()
        return table.slope(table.krg, sg)

    def krog_derivative(self, sg, swc: float = 0.0,
                        kro_end: float = 1.0) -> np.ndarray:
        table = self.get()
        return table.slope(table.krog, sg)
