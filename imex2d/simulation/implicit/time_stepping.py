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
from .newton import NewtonSolver
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

    Defolt olaraq söndürülüb; qəbul edən mühərrik bunu özü
    aktivləşdirir.
    """
    soft_failure_mb_tolerance: Optional[float] = None
    """`soft_failure_cnv_tolerance` ilə BİRLİKDƏ tələb olunur.

    CNV yalnız ƏN PİS HÜCEYRənin yerli qalığını ölçür — qlobal kütlə
    balansını yoxlamır. Bu hədd olmadan "yumşaq uğursuzluq" lokal
    cəhətdən yaxşı görünən, lakin qlobal kütləni itirən/yaradan bir
    vəziyyəti də qəbul edə bilərdi. İkisi BİRLİKDƏ tələb olunur.
    """
    max_consecutive_soft_failures: int = 20
    """TAPILAN SƏHV-in nəticəsi: bəzi hüceyrələrdə çətinlik KEÇİCİ
    deyil, DAVAMLI olanda (məs. PVT-nin doyma nöqtəsi ətrafında) hər
    addım YENİDƏN minimal Δt-yə düşüb yumşaq qəbul edilir — simulyasiya
    "dayanmır", lakin praktiki olaraq İRƏLİLƏMİR (min_dt qədər addımlarla
    sonsuza addımlayır). Bu hədd BU HALI aşkarlayır: ardıcıl yumşaq
    qəbulların sayı bunu keçsə, davam etmək əvəzinə TƏMİZ dayandırılır —
    "sonsuz sürünmə" "sonsuz gözləmə"dən daha pisdir, çünki nəticəsiz CPU
    yeyir və istifadəçini aldadır (`converged=True` amma faktiki t
    irəliləmir)."""
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
    soft_failure: bool = False
    """Tam yığılmadan (bax `soft_failure_cnv_tolerance`) XƏBƏRDARLIQLA
    qəbul edilib — nəticə fiziki cəhətdən etibarlıdır, lakin sərt
    yığılma meyarına tam çatmayıb."""


class AdaptiveTimeStepper:
    """Nyuton həlledicisini adaptiv Δt ilə idarə edir."""

    def __init__(self, newton: NewtonSolver,
                 config: Optional[AdaptiveTimeStepConfig] = None):
        self.newton = newton
        self.config = config or AdaptiveTimeStepConfig()
        self.dt = self.config.initial_dt
        self.history: List[TimeStepRecord] = []
        self._consecutive_soft_failures = 0

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
                self._consecutive_soft_failures = 0
                return result.state, dt, result

            reason = ("yığılmadı: " + result.status.value if not result.converged
                      else f"ΔSw = {change:.3f} həddi keçdi")
            if dt <= config.min_dt * (1.0 + 1e-9):
                soft_tolerance = config.soft_failure_cnv_tolerance
                mb_tolerance = config.soft_failure_mb_tolerance
                history = getattr(result, "cnv_history", None)
                last_cnv = history[-1] if history else float("inf")
                mb_water, mb_oil = getattr(
                    result, "material_balance", (float("inf"), float("inf")))
                soft_ok = (soft_tolerance is not None and mb_tolerance is not None
                          and history and last_cnv < soft_tolerance
                          and max(mb_water, mb_oil) < mb_tolerance
                          and (self._consecutive_soft_failures
                               < config.max_consecutive_soft_failures))
                if soft_ok:
                    self._consecutive_soft_failures += 1
                    LOG.warning(
                        "t = %.2f gün: minimal Δt-də tam yığılmadı, lakin "
                        "CNV=%.2e və MB=%.2e yumşaq hədlərin altındadır — "
                        "XƏBƏRDARLIQLA qəbul edilir (%d/%d ardıcıl).", time,
                        last_cnv, max(mb_water, mb_oil),
                        self._consecutive_soft_failures,
                        config.max_consecutive_soft_failures)
                    self.history.append(TimeStepRecord(
                        time, dt, result.iterations, repeat, change, True,
                        soft_failure=True))
                    self.dt = config.min_dt
                    return result.state, dt, result
                if (soft_tolerance is not None
                        and self._consecutive_soft_failures
                        >= config.max_consecutive_soft_failures):
                    LOG.error(
                        "t = %.2f gün: %d ardıcıl yumşaq qəbuldan sonra DA "
                        "irəliləmə yoxdur — model bu nöqtədə TIXANIB, "
                        "dayandırılır.", time, self._consecutive_soft_failures)
                else:
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
            "yumşaq qəbul": sum(1 for r in steps if r.soft_failure),
        }
