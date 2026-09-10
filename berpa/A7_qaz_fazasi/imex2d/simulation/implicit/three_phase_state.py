"""Üç fazalı primary dəyişənlər — A7, mərhələ 4.

İki fazalı sxemdə (`state.py`) vəziyyət sadədir: hər hüceyrə üçün
(p, Sw) — Nyuton bunları birbaşa yeniləyir.

Üç fazalı sxemdə ÜÇÜNCÜ dəyişənin MƏNASI hüceyrənin doyma
vəziyyətindən asılıdır — bu, sənaye simulyatorlarının (Eclipse, IMEX)
standart "dəyişən keçid" (variable switching) üsuludur:

    doymuş hüceyrə (sərbəst qaz var)     3-cü dəyişən = Sg
    doymamış hüceyrə (bütün qaz həll olub) 3-cü dəyişən = Rs

Niyə belə: doymamış hüceyrədə Sg həmişə 0-dır — onu primary dəyişən
etmək mənasız olardı (Nyuton sabit 0-ı "həll edərdi"). Doymuş
hüceyrədə isə Rs həmişə Rs_sat(p)-ə bərabərdir (neft artıq daha çox
qaz həll edə bilməz) — onu izləməyə ehtiyac yoxdur, ÇÜNKİ TƏZYİQDƏN
birbaşa çıxarıla bilər.

DƏYİŞƏN KEÇİD

Hər Nyuton addımından sonra hüceyrə bir vəziyyətdən digərinə keçə
bilər:

    doymamış → doymuş   Rs yeni təzyiqdə Rs_sat(p)-i keçəndə
                         (neft "doydu", artıq qaz sərbəst qalır)
    doymuş → doymamış   Sg mənfiyə düşəndə
                         (bütün sərbəst qaz yenidən həll oldu)

Keçid ANINDA kəmiyyət KƏSİLMƏZ olmalıdır — sərhəddə hər iki təsvir
eyni fiziki vəziyyəti göstərməlidir (Sg=0 ⟺ Rs=Rs_sat(p)). Bu modul
yalnız VƏZİYYƏTİ və KEÇİD MƏNTİQİNİ təsvir edir; qalıq tənlikləri və
Jakobian (mərhələ 5) bunun üzərində qurulacaq.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

import numpy as np

PRESSURE = 0
WATER_SATURATION = 1
THIRD_VARIABLE = 2       # Sg (doymuş) və ya Rs (doymamış)
VARIABLES_PER_CELL = 3


@dataclass
class ThreePhaseState:
    """Bir zaman qatının üç fazalı vəziyyəti.

    `is_saturated` Nyuton vektorunun HİSSƏSİ DEYİL — kəsilməz dəyişən
    deyil, diskret bayraqdır. `switch_variables()` onu yeniləyir,
    `to_vector()`/`from_vector()` yalnız kəsilməz üç dəyişənlə işləyir.
    """
    pressure: np.ndarray
    water_saturation: np.ndarray
    third_variable: np.ndarray
    is_saturated: np.ndarray

    def __post_init__(self):
        self.pressure = np.asarray(self.pressure, float)
        self.water_saturation = np.asarray(self.water_saturation, float)
        self.third_variable = np.asarray(self.third_variable, float)
        self.is_saturated = np.asarray(self.is_saturated, bool)

    @property
    def ncell(self) -> int:
        return int(self.pressure.size)

    @property
    def gas_saturation(self) -> np.ndarray:
        """Sg — doymamış hüceyrələrdə həmişə 0 (sərbəst qaz yoxdur)."""
        return np.where(self.is_saturated, self.third_variable, 0.0)

    @property
    def oil_saturation(self) -> np.ndarray:
        return 1.0 - self.water_saturation - self.gas_saturation

    def solution_gor(self, pvt) -> np.ndarray:
        """Rs — doymuş hüceyrələrdə Rs_sat(p) (izlənmir, hesablanır)."""
        rs_saturated = np.asarray(pvt.solution_gor(self.pressure), float)
        return np.where(self.is_saturated, rs_saturated, self.third_variable)

    def copy(self) -> "ThreePhaseState":
        return ThreePhaseState(self.pressure.copy(),
                               self.water_saturation.copy(),
                               self.third_variable.copy(),
                               self.is_saturated.copy())

    # ─────────────────────────────────────────── Nyuton vektoru
    def to_vector(self) -> np.ndarray:
        vector = np.empty(self.ncell * VARIABLES_PER_CELL)
        vector[PRESSURE::VARIABLES_PER_CELL] = self.pressure
        vector[WATER_SATURATION::VARIABLES_PER_CELL] = self.water_saturation
        vector[THIRD_VARIABLE::VARIABLES_PER_CELL] = self.third_variable
        return vector

    @classmethod
    def from_vector(cls, vector: np.ndarray,
                    is_saturated: np.ndarray) -> "ThreePhaseState":
        vector = np.asarray(vector, float)
        return cls(vector[PRESSURE::VARIABLES_PER_CELL].copy(),
                   vector[WATER_SATURATION::VARIABLES_PER_CELL].copy(),
                   vector[THIRD_VARIABLE::VARIABLES_PER_CELL].copy(),
                   is_saturated)

    def updated(self, delta: np.ndarray, sw_min: float, sw_max: float,
               max_pressure_change: Optional[float] = None,
               max_saturation_change: Optional[float] = None
               ) -> "ThreePhaseState":
        """Nyuton addımını tətbiq edir (Appleyard chopping, `state.py` ilə eyni üsul).

        3-cü dəyişənin özü burada MƏHDUDLAŞDIRILMIR (Sg üçün [0,1],
        Rs üçün fərqli miqyas ola bilər) — o, `switch_variables()`-də
        həll olunur, çünki hədd YALNız cari vəziyyət bilinəndə (Sg
        yoxsa Rs) mənalıdır.
        """
        dp = delta[PRESSURE::VARIABLES_PER_CELL]
        dsw = delta[WATER_SATURATION::VARIABLES_PER_CELL]
        dthird = delta[THIRD_VARIABLE::VARIABLES_PER_CELL]

        if max_pressure_change:
            dp = np.clip(dp, -max_pressure_change, max_pressure_change)
        if max_saturation_change:
            dsw = np.clip(dsw, -max_saturation_change, max_saturation_change)
            # Sg üçün eyni fiziki miqyasda addım məhdudlaşdırılır;
            # Rs dəyişəni tamam fərqli vahiddədir (sm3/sm3), ona görə
            # yalnız doymuş (Sg) hüceyrələrdə tətbiq olunur.
            dthird = np.where(self.is_saturated,
                              np.clip(dthird, -max_saturation_change,
                                     max_saturation_change), dthird)

        return ThreePhaseState(
            self.pressure + dp,
            np.clip(self.water_saturation + dsw, sw_min, sw_max),
            self.third_variable + dthird,
            self.is_saturated.copy())

    # ─────────────────────────────────────────────── dəyişən keçid
    def switch_variables(self, pvt) -> "ThreePhaseState":
        """Doyma sərhədini keçən hüceyrələri yeni təsvirə köçürür.

        `pvt` — `has_gas_phase()` True olan `IPVTProvider`; `Rs_sat(p)`
        üçün lazımdır.

        Qaytarılan vəziyyətdə HƏR HÜCEYRƏNİN fiziki mənası (Sw, Sg, So)
        giriş ilə EYNİDİR — dəyişən yalnız NECƏ TƏSVİR OLUNDUĞUdur.
        """
        rs_saturated = np.asarray(pvt.solution_gor(self.pressure), float)
        third = self.third_variable.copy()
        saturated = self.is_saturated.copy()

        # doymamış -> doymuş: Rs (cari 3-cü dəyişən) doyma əyrisini keçib
        becomes_saturated = (~self.is_saturated) & (third > rs_saturated)
        third[becomes_saturated] = 0.0          # Sg sərhəddə sıfırdan başlayır
        saturated[becomes_saturated] = True

        # doymuş -> doymamış: Sg (cari 3-cü dəyişən) mənfiyə düşüb
        becomes_undersaturated = self.is_saturated & (third < 0.0)
        # kəsilməzlik: keçid anında Rs = Rs_sat(p) (sərhəd dəyəri)
        third[becomes_undersaturated] = rs_saturated[becomes_undersaturated]
        saturated[becomes_undersaturated] = False

        return ThreePhaseState(self.pressure.copy(),
                               self.water_saturation.copy(),
                               third, saturated)


def index_of(cell: int, variable: int) -> int:
    return cell * VARIABLES_PER_CELL + variable
