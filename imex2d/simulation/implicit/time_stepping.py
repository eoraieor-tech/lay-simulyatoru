"""Adaptiv zaman addımı — A6, mərhələ 4.

Strategiya (sənaye simulyatorlarında standart):

    Nyuton yığılmadı  -> Δt-ni kəs, addımı TƏKRARLA
    az iterasiya      -> Δt-ni böyüt
    çox iterasiya     -> Δt-ni kiçilt (təkrar yox, növbəti addım üçün)

Mərhələ 3-də ölçüldü ki, Nyutonun yığılması `Δt·q/PV` nisbətindən
asılıdır və hədd təxminən 10-15-dir. Adaptiv nəzarət bu həddi
avtomatik tapır — istifadəçi Δt seçməli deyil.

Əlavə məhdudiyyət: doyumluluq dəyişikliyi. Bir addımda cəbhənin çox
irəliləməsi Nyutonu çətinləşdirir və nəticənin dəqiqliyini azaldır.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from ...logging_setup import get_logger
from .newton import NewtonResult, NewtonSolver, NewtonStatus
from .state import ReservoirState

LOG = get_logger(__name__)


@dataclass
class AdaptiveTimeStepConfig:
    initial_dt: float = 1.0
    min_dt: float = 1e-3
    max_dt: float = 365.0

    target_iterations: int = 6
    """Bu qədər iterasiyada yığılma "ideal" sayılır."""

    growth_factor: float = 1.5
    cut_factor: float = 0.5

    soft_failure_cnv_tolerance: Optional[float] = None
    """İSTƏYƏ BAĞLI (opt-in) — A6-nın iki fazalı davranışını dəyişmir.

    `None` (defolt) — köhnə davranış: min Δt-də yığılmasa, tam
    dayan (`dt=0.0` qaytarılır).

    Verilibsə — min Δt-də DƏ yığılmasa, son sınağın CNV-si bu
    həddən aşağıdırsa, nəticə XƏBƏRDARLIQLA qəbul edilir (tam
    dayanmaq əvəzinə). Sənayedə tanınan "yumşaq uğursuzluq" (soft
    failure) təcrübəsidir: CNV fiziki cəhətdən kifayət qədər
    kiçikdirsə (mb. material balans da yoxlanılır), simulyasiyanın
    tam dayanması ƏSASSIZDIR — bir neçə hüceyrə çətin nöqtədə
    (məs. quyunun öz BHP sərhədinə çox yaxın olması) "asılıb" qala
    bilər, halbuki qlobal nəticə fiziki cəhətdən etibarlıdır.

    A7-nin üç fazalı mühərriki bunu aktivləşdirir — bax
    `three_phase_engine.py`.
    """
    max_growth_per_step: float = 2.0

    max_saturation_change: float = 0.2
    """Bir addımda hüceyrədə icazə verilən maksimal ΔSw.

    Yığılma meyarı deyil — DƏQİQLİK meyarıdır. Nyuton böyük ΔSw ilə
    də yığıla bilər, lakin cəbhə bir addımda çox irəliləyəndə həll
    fiziki cəhətdən kobudlaşır.
    """

    max_repeats: int = 10
    """Bir addımı bu qədər dəfə kəsdikdən sonra imtina edilir."""


@dataclass
class TimeStepRecord:
    time: float
    dt: float
    iterations: int
    repeats: int
    max_saturation_change: float
    converged: bool


class AdaptiveTimeStepper:
    """Nyuton həlledicisini adaptiv Δt ilə idarə edir."""

    def __init__(self, newton: NewtonSolver,
                 config: Optional[AdaptiveTimeStepConfig] = None):
        self.newton = newton
        self.config = config or AdaptiveTimeStepConfig()
        self.dt = self.config.initial_dt
        self.history: List[TimeStepRecord] = []

    # ═══════════════════════════════════════════════════ bir addım
    def advance(self, state: ReservoirState, time: float,
                remaining: float) -> tuple:
        """Bir zaman addımı atır (lazım gələrsə təkrarlayaraq).

        Qaytarır: (yeni_vəziyyət, atılan_dt, NewtonResult) və ya
        uğursuzluqda (state, 0.0, son_nəticə).
        """
        config = self.config
        dt = min(self.dt, config.max_dt, remaining)

        # MƏRHƏLƏ 6 (A7): quyu debitinə görə əlavə hədd — bax
        # `CoupledNewtonSolver.max_stable_dt()`-in sənədləşməsi.
        #
        # `hasattr` yoxlaması qəsdən işlədilir: A6-nın iki fazalı
        # `NewtonSolver`-i bu metodu TANIMIR, ona görə köhnə mühərrik
        # HEÇ BİR DƏYİŞİKLİK olmadan işləməyə davam edir — yalnız
        # quyu debiti həddini bilən həlledicilər (indi: üç fazalı
        # birləşmiş model) bundan faydalanır.
        well_limit = getattr(self.newton, "max_stable_dt", None)
        if well_limit is not None:
            dt = min(dt, well_limit(state))

        for repeat in range(config.max_repeats + 1):
            result = self.newton.solve(state, dt)
            change = self._saturation_change(state, result.state)

            if result.converged and change <= config.max_saturation_change:
                self.history.append(TimeStepRecord(
                    time, dt, result.iterations, repeat, change, True))
                self.dt = self._next_dt(dt, result.iterations, change)
                return result.state, dt, result

            reason = ("yığılmadı: " + result.status.value if not result.converged
                      else f"ΔSw = {change:.3f} həddi keçdi")
            if dt <= config.min_dt * (1.0 + 1e-9):
                soft_tolerance = config.soft_failure_cnv_tolerance
                history = getattr(result, "history", None)
                last_cnv = history[-1] if history else float("inf")
                if soft_tolerance is not None and history and last_cnv < soft_tolerance:
                    LOG.warning(
                        "t = %.2f gün: minimal Δt-də tam yığılmadı, lakin "
                        "CNV=%.2e yumşaq həddin (%.2e) altındadır — "
                        "XƏBƏRDARLIQLA qəbul edilir.", time, last_cnv,
                        soft_tolerance)
                    self.history.append(TimeStepRecord(
                        time, dt, result.iterations, repeat, change, True))
                    self.dt = config.min_dt
                    return result.state, dt, result
                LOG.error("t = %.2f gün: minimal Δt-də də %s", time, reason)
                self.history.append(TimeStepRecord(
                    time, dt, result.iterations, repeat, change, False))
                return state, 0.0, result

            dt = max(dt * config.cut_factor, config.min_dt)
            LOG.debug("t = %.2f gün: %s -> Δt = %.4f gün", time, reason, dt)

        return state, 0.0, result

    # ═══════════════════════════════════════════ növbəti addımın ölçüsü
    def _next_dt(self, dt: float, iterations: int, change: float) -> float:
        """İterasiya sayı və doyumluluq dəyişikliyinə görə tənzimləmə."""
        config = self.config

        if iterations <= config.target_iterations:
            factor = config.growth_factor
        elif iterations <= config.target_iterations * 1.5:
            factor = 1.0
        else:
            factor = config.cut_factor

        # doyumluluq dəyişikliyi hədə yaxınlaşırsa, böyüməni ləngit
        if change > 0.0:
            allowed = config.max_saturation_change / change
            factor = min(factor, max(allowed * 0.9, config.cut_factor))

        factor = min(factor, config.max_growth_per_step)
        return float(np.clip(dt * factor, config.min_dt, config.max_dt))

    @staticmethod
    def _saturation_change(before: ReservoirState, after: ReservoirState) -> float:
        return float(np.abs(after.water_saturation
                            - before.water_saturation).max())

    # ═══════════════════════════════════════════════════════ statistika
    def summary(self) -> dict:
        if not self.history:
            return {}
        steps = [record for record in self.history if record.converged]
        return {
            "addım": len(steps),
            "təkrar": sum(record.repeats for record in self.history),
            "orta Δt": float(np.mean([r.dt for r in steps])) if steps else 0.0,
            "maks Δt": float(max(r.dt for r in steps)) if steps else 0.0,
            "orta iterasiya": float(np.mean([r.iterations for r in steps]))
            if steps else 0.0,
        }
