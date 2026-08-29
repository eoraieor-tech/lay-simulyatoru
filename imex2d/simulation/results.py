"""Simulyasiya nəticələri — SAF MƏLUMAT, çəkmə (rendering) yoxdur."""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np


@dataclass
class TimeSeries:
    time: List[float] = field(default_factory=list)
    oil_rate: List[float] = field(default_factory=list)
    water_rate: List[float] = field(default_factory=list)
    water_injection_rate: List[float] = field(default_factory=list)
    cumulative_oil: List[float] = field(default_factory=list)
    cumulative_water: List[float] = field(default_factory=list)
    water_cut: List[float] = field(default_factory=list)
    average_pressure: List[float] = field(default_factory=list)
    recovery_factor: List[float] = field(default_factory=list)


@dataclass
class Snapshot:
    time: float
    pressure: np.ndarray
    water_saturation: np.ndarray


@dataclass
class SimulationResult:
    model_name: str = ""
    grid_shape: tuple = ()
    series: TimeSeries = field(default_factory=TimeSeries)
    snapshots: List[Snapshot] = field(default_factory=list)
    well_oil_rate: Dict[str, List[float]] = field(default_factory=dict)
    well_water_rate: Dict[str, List[float]] = field(default_factory=dict)
    ooip: float = 0.0
    steps: int = 0
    converged: bool = True
    message: str = ""

    @property
    def final_recovery_factor(self) -> float:
        return self.series.recovery_factor[-1] if self.series.recovery_factor else 0.0

    @property
    def breakthrough_time(self):
        for t, wc in zip(self.series.time, self.series.water_cut):
            if wc > 1.0:
                return t
        return None
