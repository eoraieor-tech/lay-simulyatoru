"""KONFİQURASİYA OBYEKTLƏRİ — bütün sabitlər buradadır.

Əvvəllər bu dəyərlər üç yerdə səpələnmişdi: mühərrikin daxilində
(rtol=1e-8, drop_tol=1e-4, 25 addım, 1.15 artım əmsalı), UI spinbox
default-larında və nüvənin dataclass sahələrində. İndi tək mənbə.
"""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class LinearSolverConfig:
    """Xətti həlledicinin parametrləri.

    `preconditioner_refresh_steps` = 50 ölçmə ilə seçilib: 25-də ILU
    faktorlaşdırması çox tez-tez təkrarlanır, 400-də isə ön-şərtçi
    köhnəlir və KQ iterasiyaları artır.

    `ilu_drop_tolerance` böyüdülməməlidir — zəif ILU KQ-nin iterasiya
    sayını kəskin artırır (ölçmədə 1e-3 keçidi 30 dəfə yavaşlatdı).
    """
    tolerance: float = 1e-8
    max_iterations: int = 200
    preconditioner_refresh_steps: int = 50
    ilu_drop_tolerance: float = 1e-4
    ilu_fill_factor: float = 10.0
    fallback_to_direct: bool = True


@dataclass
class TimeSteppingConfig:
    initial_dt: float = 1.0
    max_dt: float = 20.0
    min_dt: float = 1e-4
    cfl_factor: float = 0.45
    growth_factor: float = 1.15
    max_steps: int = 20000


@dataclass
class OutputConfig:
    snapshot_count: int = 60
    record_well_rates: bool = True
    progress_every_n_steps: int = 5


@dataclass
class SimulationConfig:
    """Bir simulyasiya işə salınmasının tam konfiqurasiyası."""
    end_time: float = 1500.0
    time_stepping: TimeSteppingConfig = field(default_factory=TimeSteppingConfig)
    linear_solver: LinearSolverConfig = field(default_factory=LinearSolverConfig)
    output: OutputConfig = field(default_factory=OutputConfig)

    def validate(self) -> list:
        issues = []
        if self.end_time <= 0:
            issues.append("Simulyasiya müddəti müsbət olmalıdır.")
        ts = self.time_stepping
        if not (0 < ts.cfl_factor <= 1.0):
            issues.append("CFL əmsalı (0, 1] intervalında olmalıdır.")
        if ts.max_dt <= ts.min_dt:
            issues.append("Maksimal Δt minimaldan böyük olmalıdır.")
        return issues
