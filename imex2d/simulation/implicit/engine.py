"""FullyImplicitEngine — ISimulationEngine implementasiyası.

IMPES mühərriki ilə eyni interfeys, eyni ReservoirModel, eyni
SimulationResult. Fərq yalnız həll üsulundadır:

    IMPES            təzyiq implicit, doyumluluq explicit, CFL məhdud
    FullyImplicit    hər ikisi implicit, Nyuton, adaptiv Δt

Hər ikisi saxlanılır, çünki kiçik modellərdə IMPES daha sürətlidir
(hər addımı ucuzdur, Nyuton iterasiyası yoxdur). Seçim istifadəçinindir.
"""

from __future__ import annotations
from typing import Optional

import numpy as np

from ...application.config import SimulationConfig
from ...domain.reservoir_model import ReservoirModel
from ...interfaces.providers import (ICapillaryPressureProvider,
                                     IInitializationProvider, IPVTProvider,
                                     IRelativePermeabilityProvider)
from ...interfaces.services import (IProgressReporter, ISimulationEngine,
                                    NullProgressReporter)
from ...interfaces.discretization import IFluxDiscretization
from ...logging_setup import get_logger
from ..discretization import default_flux_discretization
from ..results import SimulationResult, Snapshot
from ..well_model import PeacemanWellModel
from .jacobian import JacobianAssembler
from .linear import NewtonLinearSolver
from .newton import NewtonConfig, NewtonSolver
from .residual import ResidualAssembler
from .state import ReservoirState
from .time_stepping import AdaptiveTimeStepConfig, AdaptiveTimeStepper

LOG = get_logger(__name__)

#: PHASE D (5B-2): `FullyImplicitEngine` ARTIQ MPFA-O ilə işləyir (bax
#: `jacobian.py` modul docstring-i, "PHASE 5B-2" bölməsi) — burada əvvəllər
#: mövcud olan `_reject_multipoint_engine(...)` çağırışı (AÇIQ imtina)
#: SİLİNİB. `ImpesEngine` HƏLƏ DƏ öz AYRICA `_reject_multipoint_impes`-i
#: ilə (bax `impes_engine.py`) rədd edir — onun təzyiq addımı tək-üz
#: transmissivlik skalyarına əsaslanır, MPFA-nın çoxnöqtəli `T_conn`-u ilə
#: RİYAZİ CƏHƏTDƏN uyğun deyil (saxta "orta transmissivlik" uydurmaq
#: QADAĞANDIR) — IMPES+MPFA-O HƏLƏ DƏ implement EDİLMƏYİB, gizlədilmir.


class FullyImplicitEngine(ISimulationEngine):

    def __init__(self,
                 model: ReservoirModel,
                 config: SimulationConfig,
                 relperm: IRelativePermeabilityProvider,
                 linear_solver=None,
                 pvt: Optional[IPVTProvider] = None,
                 capillary: Optional[ICapillaryPressureProvider] = None,
                 initialization: Optional[IInitializationProvider] = None,
                 newton_config: Optional[NewtonConfig] = None,
                 time_step_config: Optional[AdaptiveTimeStepConfig] = None,
                 flux_discretization: Optional[IFluxDiscretization] = None):
        self.model = model
        self.config = config
        self.relperm = relperm
        self.pvt = pvt
        self.capillary = capillary
        self.initialization = initialization

        #: DEFOLT = TPFA (bax audit tapşırığı §4). `flux_discretization`
        #: AÇIQ verilməyibsə mövcud davranış BİRƏBİR eynidir. PHASE D:
        #: `MPFAODiscretization()` ötürülsə, ResidualAssembler/Jacobian
        #: AVTOMATİK çoxnöqtəli yola keçir (bax `jacobian.py` "PHASE 5B-2").
        self.flux_discretization = flux_discretization or default_flux_discretization()
        grid = self.flux_discretization.build(model)
        wells = PeacemanWellModel().build_connections(model)
        self.residual_assembler = ResidualAssembler(
            model, grid, wells, relperm, pvt=pvt, capillary=capillary)
        self.jacobian_assembler = JacobianAssembler(self.residual_assembler)
        self.newton = NewtonSolver(
            self.residual_assembler, self.jacobian_assembler,
            NewtonLinearSolver(), newton_config)
        self.time_stepper = AdaptiveTimeStepper(
            self.newton, time_step_config or self._time_config(config))

        self._producers = sorted({c.well_name for c in wells
                                  if not c.is_injector})
        self.state = self._initial_state()
        # IMPES mühərriki ilə eyni atributlar — testlər və UI üçün
        self.pressure = self.state.pressure
        self.sw = self.state.water_saturation

    # ─────────────────────────────────────────────────────── qurulma
    @staticmethod
    def _time_config(config: SimulationConfig) -> AdaptiveTimeStepConfig:
        """Ümumi konfiqurasiyadan adaptiv addım parametrləri.

        `max_dt` istifadəçinin sorduğu kimi hörmət edilir. ƏVVƏLLƏR
        (TAPILAN SƏHV) bura süni minimum (30 gün) tətbiq olunurdu —
        istifadəçi 0.5 və ya 2 gün desə də, mühərrik səssizcə 30 günə
        keçirdi. Bu, `max_dt`-nin nəticəyə təsirini yoxlayan sınaqları
        çaşdırırdı, çünki 30-dan kiçik HƏR DƏYƏR eyni davranışı verirdi.

        `soft_failure_*` — son-çarə təhlükəsizlik toru: minimal Δt-də
        DƏ tam yığılmasa, YALNIZ HƏM CNV, HƏM DƏ qlobal kütlə balansı
        kifayət qədər kiçikdirsə, addım tam dayanmaq əvəzinə
        XƏBƏRDARLIQLA qəbul edilir (bax `AdaptiveTimeStepConfig`).
        """
        stepping = config.time_stepping
        return AdaptiveTimeStepConfig(
            initial_dt=stepping.initial_dt,
            min_dt=stepping.min_dt,
            max_dt=stepping.max_dt,
            growth_factor=stepping.growth_factor + 0.35,
            soft_failure_cnv_tolerance=1e-2,
            soft_failure_mb_tolerance=1e-4)

    def _initial_state(self) -> ReservoirState:
        ic = self.model.initial_conditions
        if self.initialization is not None:
            initial = self.initialization.initialize(self.model)
            state = ReservoirState(np.asarray(initial.pressure, float).copy(),
                                   np.asarray(initial.water_saturation,
                                              float).copy())
        else:
            state = ReservoirState(
                np.full(self.model.ncell, ic.datum_pressure),
                np.full(self.model.ncell, ic.water_saturation))
        sw_min, sw_max = self.relperm.saturation_limits()
        state.water_saturation = np.clip(state.water_saturation, sw_min, sw_max)
        return state

    def original_oil_in_place(self) -> float:
        fluid = self.residual_assembler.fluid_state(self.state)
        return float(self.residual_assembler.accumulation(self.state, fluid)[1].sum())

    # ══════════════════════════════════════════════════════════ run
    def run(self, reporter: Optional[IProgressReporter] = None) -> SimulationResult:
        reporter = reporter or NullProgressReporter()
        config = self.config
        output = config.output

        result = SimulationResult(model_name=self.model.name,
                                  grid_shape=self.model.grid.shape)
        result.ooip = self.original_oil_in_place()
        result.well_oil_rate = {name: [] for name in self._producers}
        result.well_water_rate = {name: [] for name in self._producers}
        series = result.series

        snapshot_interval = max(config.end_time / max(output.snapshot_count, 1),
                                1e-9)
        next_snapshot = 0.0
        self._record_snapshot(result, 0.0)
        next_snapshot += snapshot_interval

        time = 0.0
        steps = 0
        cumulative_oil = cumulative_water = 0.0

        while time < config.end_time - 1e-9:
            new_state, dt, newton_result = self.time_stepper.advance(
                self.state, time, config.end_time - time)

            if dt <= 0.0:
                result.converged = False
                result.message = (f"t = {time:.1f} gün: zaman addımı "
                                  f"minimal həddə də yığılmadı "
                                  f"({newton_result.status.value}).")
                break

            self.state = new_state
            self.pressure = self.state.pressure
            self.sw = self.state.water_saturation
            time += dt
            steps += 1

            rates = newton_result.rates
            oil_rate = float(-min(rates.oil.sum(), 0.0))
            water_rate = float(-min(rates.water[rates.water < 0].sum(), 0.0))
            injection = float(max(rates.water[rates.water > 0].sum(), 0.0))
            cumulative_oil += oil_rate * dt
            cumulative_water += water_rate * dt

            series.time.append(time)
            series.oil_rate.append(oil_rate)
            series.water_rate.append(water_rate)
            series.water_injection_rate.append(injection)
            series.cumulative_oil.append(cumulative_oil)
            series.cumulative_water.append(cumulative_water)
            series.water_cut.append(water_rate / max(oil_rate + water_rate, 1e-12)
                                    * 100.0)
            series.average_pressure.append(float(np.mean(self.pressure)))
            series.recovery_factor.append(cumulative_oil
                                          / max(result.ooip, 1e-12) * 100.0)
            if output.record_well_rates:
                for name in self._producers:
                    result.well_oil_rate[name].append(
                        float(-rates.per_well_oil.get(name, 0.0)))
                    result.well_water_rate[name].append(
                        float(-rates.per_well_water.get(name, 0.0)))

            if time >= next_snapshot - 1e-9:
                self._record_snapshot(result, time)
                next_snapshot += snapshot_interval

            if steps % max(output.progress_every_n_steps, 1) == 0:
                message = (f"t = {time:8.1f} gün | RF = "
                           f"{series.recovery_factor[-1]:5.2f} % | "
                           f"dt = {dt:6.2f} | Nyuton {newton_result.iterations}")
                if not reporter.report(time / config.end_time * 100.0, message):
                    result.message = "İstifadəçi tərəfindən dayandırıldı."
                    break

        if result.snapshots and result.snapshots[-1].time < time - 1e-9:
            self._record_snapshot(result, time)
        result.steps = steps
        if not result.message:
            statistics = self.time_stepper.summary()
            soft = statistics.get("yumşaq qəbul", 0)
            soft_note = (f" DİQQƏT: {soft} addım tam yığılmadan "
                        f"XƏBƏRDARLIQLA qəbul edildi — log-a baxın."
                        if soft else "")
            result.message = (
                f"Tamamlandı: {steps} addım, t = {time:.1f} gün "
                f"(orta Δt = {statistics.get('orta Δt', 0.0):.1f} gün, "
                f"orta {statistics.get('orta iterasiya', 0.0):.1f} Nyuton "
                f"iterasiyası, {statistics.get('təkrar', 0)} təkrar)."
                f"{soft_note}")
        LOG.info("%s  RF = %.2f %%", result.message,
                 result.final_recovery_factor)
        return result

    def _record_snapshot(self, result: SimulationResult, time: float) -> None:
        shape = self.model.grid.shape
        result.snapshots.append(Snapshot(
            time=time,
            pressure=self.state.pressure.reshape(shape).copy(),
            water_saturation=self.state.water_saturation.reshape(shape).copy()))
