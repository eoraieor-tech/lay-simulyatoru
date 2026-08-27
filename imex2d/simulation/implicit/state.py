"""Primary dəyişənlər — fully implicit sxemin vəziyyət vektoru.

IMPES-də təzyiq və doyumluluq ayrı-ayrı həll olunurdu. Fully implicit
sxemdə hər ikisi eyni Nyuton sistemində iştirak edir, ona görə onlar
tək vektorda birləşdirilir:

    x = [p_0, Sw_0, p_1, Sw_1, …, p_{N-1}, Sw_{N-1}]

Hüceyrə üzrə növbələşdirmə (interleaved) qəsdən seçilib: Jakobian
2×2 blok-diaqonal struktura düşür və CPR ön-şərtçisi (A6-nın son
mərhələsi) məhz bu strukturu tələb edir.
"""

from __future__ import annotations
from dataclasses import dataclass

import numpy as np

PRESSURE = 0
WATER_SATURATION = 1
VARIABLES_PER_CELL = 2


@dataclass
class ReservoirState:
    """Bir zaman qatının vəziyyəti."""
    pressure: np.ndarray
    water_saturation: np.ndarray

    def copy(self) -> "ReservoirState":
        return ReservoirState(self.pressure.copy(), self.water_saturation.copy())

    @property
    def ncell(self) -> int:
        return int(self.pressure.size)

    def to_vector(self) -> np.ndarray:
        vector = np.empty(self.ncell * VARIABLES_PER_CELL)
        vector[PRESSURE::VARIABLES_PER_CELL] = self.pressure
        vector[WATER_SATURATION::VARIABLES_PER_CELL] = self.water_saturation
        return vector

    @classmethod
    def from_vector(cls, vector: np.ndarray) -> "ReservoirState":
        vector = np.asarray(vector, float)
        return cls(vector[PRESSURE::VARIABLES_PER_CELL].copy(),
                   vector[WATER_SATURATION::VARIABLES_PER_CELL].copy())

    def updated(self, delta: np.ndarray, sw_min: float, sw_max: float,
                max_pressure_change: float = None,
                max_saturation_change: float = None) -> "ReservoirState":
        """Nyuton addımını tətbiq edir və dəyişikliyi məhdudlaşdırır.

        Məhdudlaşdırma (Appleyard chopping) Nyutonun uzaq başlanğıcdan
        divergensiya etməsinin qarşısını alır: bir iterasiyada doyumluluq
        bütün intervalı keçə bilməz.
        """
        dp = delta[PRESSURE::VARIABLES_PER_CELL]
        dsw = delta[WATER_SATURATION::VARIABLES_PER_CELL]

        if max_pressure_change:
            dp = np.clip(dp, -max_pressure_change, max_pressure_change)
        if max_saturation_change:
            dsw = np.clip(dsw, -max_saturation_change, max_saturation_change)

        return ReservoirState(
            self.pressure + dp,
            np.clip(self.water_saturation + dsw, sw_min, sw_max))


def index_of(cell: int, variable: int) -> int:
    return cell * VARIABLES_PER_CELL + variable
