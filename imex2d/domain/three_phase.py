"""Üç fazalı doyumluluq — A7, mərhələ 2.

İki fazalı modeldə Sw + So = 1 və So sərbəst dəyişən deyil (Sw-dan
hesablanır). Üç fazalı modeldə üçüncü fazanın (qaz) əlavə olunması
tək bir asılılıq gətirir:

    Sw + So + Sg = 1

Yəni İKİ sərbəst doyumluluq var (adətən Sw və Sg), üçüncüsü (So)
onlardan hesablanır — Sw-dan So-nun hesablandığı iki fazalı modellə
eyni məntiq, sadəcə bir dəyişən artıq.

DOYMA VƏZİYYƏTİ (saturated / undersaturated)

Hüceyrə iki rejimdən birindədir:

    doymamış (undersaturated)   Sg = 0,  bütün qaz neftdə həll olub,
                                 Rs < Rs_sat(p) — neft "ac" qalıb
    doymuş (saturated)          Sg > 0,  sərbəst qaz fazası mövcuddur,
                                 Rs = Rs_sat(p) — neft artıq daha çox
                                 qaz həll edə bilməz

Bu tənlik qrupunun (mərhələ 4-6) əsas mürəkkəbliyi buradan gəlir:
hüceyrə bir addımda bir rejimdən digərinə keçə bilər (dəyişən keçid).
Bu modul yalnız VƏZİYYƏTİ TƏSVİR edir — keçid məntiqi mərhələ 4-dədir.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional

import numpy as np


@dataclass
class ThreePhaseSaturation:
    """Bir zaman qatının üç fazalı doyumluluq sahəsi."""
    water: np.ndarray
    gas: np.ndarray

    def __post_init__(self):
        self.water = np.asarray(self.water, dtype=float).ravel()
        self.gas = np.asarray(self.gas, dtype=float).ravel()

    @property
    def oil(self) -> np.ndarray:
        """So = 1 − Sw − Sg — sərbəst dəyişən deyil, həmişə hesablanır."""
        return 1.0 - self.water - self.gas

    @property
    def is_saturated(self) -> np.ndarray:
        """Sərbəst qaz olan hüceyrələr (bool massiv)."""
        return self.gas > 1e-9

    @property
    def gas_cell_count(self) -> int:
        return int(np.sum(self.is_saturated))

    def validate(self, tolerance: float = 1e-6) -> List[str]:
        issues = []
        if self.water.shape != self.gas.shape:
            issues.append("Sw və Sg massivlərinin ölçüsü fərqlidir.")
            return issues
        if np.any(self.water < -tolerance) or np.any(self.water > 1.0 + tolerance):
            issues.append("Sw [0, 1] intervalından kənardadır.")
        if np.any(self.gas < -tolerance) or np.any(self.gas > 1.0 + tolerance):
            issues.append("Sg [0, 1] intervalından kənardadır.")
        oil = self.oil
        if np.any(oil < -tolerance):
            issues.append("Sw + Sg > 1 — So mənfi çıxır "
                          f"(ən pis hüceyrədə {float(oil.min()):.4f}).")
        total = self.water + self.oil + self.gas
        if np.any(np.abs(total - 1.0) > tolerance):
            issues.append("Sw + So + Sg ≠ 1 — kütlə balansı pozulub.")
        return issues

    def clip(self, sw_min: float, sw_max: float, sg_min: float = 0.0,
            sg_max: float = 1.0) -> "ThreePhaseSaturation":
        """Hədlər daxilinə salır — Sw+Sg ≤ 1 məhdudiyyəti də qorunur.

        Sadə ayrı-ayrı `clip` kifayət etmir: Sw=0.9, Sg=0.9 ayrılıqda
        keçərli olsa da, cəmi 1-dən böyükdür. Ona görə əvvəlcə hər
        biri öz həddinə, sonra CƏM 1-dən böyükdürsə nisbi miqyaslanır.
        """
        water = np.clip(self.water, sw_min, sw_max)
        gas = np.clip(self.gas, sg_min, sg_max)
        total = water + gas
        overflow = total > 1.0
        if np.any(overflow):
            scale = np.where(overflow, 1.0 / np.maximum(total, 1e-12), 1.0)
            water = water * scale
            gas = gas * scale
        return ThreePhaseSaturation(water, gas)


def saturation_state(pressure: np.ndarray, solution_gor: np.ndarray,
                     pvt) -> np.ndarray:
    """Hər hüceyrə doymuşdur, ya doymamış — cari Rs-i cədvəllə tutuşdurur.

    `pvt` — `has_gas_phase()` True olan `IPVTProvider`. Qaytarılan
    bool massivdə `True` = doymuş (sərbəst qaz var).

    Bir hüceyrənin Rs-i cədvəlin doyma əyrisindəki `Rs_sat(p)`-ə
    demək olar bərabərdirsə (fərq < tolerans), hüceyrə doyma
    sərhəddindədir — praktikada bu, mərhələ 4-də Nyuton addımının
    nəticəsində yaranan sərhəd halıdır.
    """
    pressure = np.asarray(pressure, float)
    solution_gor = np.asarray(solution_gor, float)
    rs_saturated = np.asarray(pvt.solution_gor(pressure), float)
    return solution_gor >= rs_saturated - 1e-9
