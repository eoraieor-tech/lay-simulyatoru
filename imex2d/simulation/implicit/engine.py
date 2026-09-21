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
from ..well_constraints import (MAX_SURFACE_ITERATIONS, SURFACE_RATE_TOLERANCE,
                                BhpLimitController, SurfaceRateController,
                                assign_rate_shares, needs_rate_allocation)
from ..well_model import PeacemanWellModel
from ..wellbore.thp_control import (MAX_OUTER_ITERATIONS,
                                    OUTER_TOLERANCE_BAR, ThpController)
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
        self._injectors = sorted({c.well_name for c in wells
                                  if c.is_injector})
        #: Çox perforasiyalı RATE quyusu varmı (Seans 27) — yoxdursa pay
        #: heç vaxt yenilənmir və qaçış əvvəlki kimidir
        self._rate_allocation = needs_rate_allocation(wells)
        # B4-B: THP quyuları — bağlantı hədəfi addım-addım yenilənir
        self.thp_control = ThpController(model, wells, pvt=pvt,
                                         fluids=model.fluids)
        # B7 addım 3: SƏTH bazalı RATE hədəfi — addımlar arasında lay həcminə
        # çevrilir (bax `_surface_rate_loop`)
        self.surface_rate = SurfaceRateController(wells, self._surface_factors)
        #: səth debiti düzəlişinə görə addımın ƏLAVƏ həllərinin sayı
        self.surface_rate_resolves = 0
        # B7 addım 2: RATE quyularının BHP həddi — rejim addımlar arasında
        # dəyişir (bax `_bhp_limit_loop`)
        self.bhp_limit = BhpLimitController(wells, self._connection_mobilities,
                                            self.surface_rate.rate_target)
        #: BHP limiti keçidlərinə görə addımın ƏLAVƏ həllərinin sayı
        self.bhp_limit_resolves = 0

        #: Yarı-implicit THP dövrəsində edilən ƏLAVƏ həllərin sayı
        #: (diaqnostika üçün — sıfır olması dövrənin işə düşmədiyini bildirir)
        self.thp_outer_iterations = 0
        self.state = self._initial_state()
        # THP nəzarətçisi ilkin BHP-ni LAY TƏZYİQİNİ bilərək qoyur
        # (bax `ThpController.initialize`) — ona görə vəziyyət
        # qurulandan SONRA çağırılır.
        self.thp_control.initialize(self.state.pressure)
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
        result.well_water_injection_rate = {name: [] for name in self._injectors}
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
            self._update_rate_shares()
            if self.surface_rate.active:
                self.surface_rate.predict(self.state)
            new_state, dt, newton_result = self.time_stepper.advance(
                self.state, time, config.end_time - time)

            if dt <= 0.0:
                result.converged = False
                result.message = (f"t = {time:.1f} gün: zaman addımı "
                                  f"minimal həddə də yığılmadı "
                                  f"({newton_result.status.value}).")
                break

            # Yarı-implicit THP: BHP bu addımın öz debitləri ilə yenilənir
            # (bax `_thp_outer_loop`). `self.state` HƏLƏ köhnə vəziyyətdir —
            # təkrar məhz ondan başlamalıdır.
            new_state, newton_result = self._thp_outer_loop(
                new_state, dt, newton_result)
            new_state, newton_result = self._surface_rate_loop(
                new_state, dt, newton_result)
            new_state, newton_result = self._bhp_limit_loop(
                new_state, dt, newton_result)

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
            # orta təzyiq YALNIZ aktiv hüceyrələr üzrə — qeyri-aktiv
            # hüceyrənin təzyiqi naməlum deyil, sadəcə saxlanılan
            # ilkin dəyərdir (bax `domain/grid.py` KONVENSİYA)
            series.average_pressure.append(float(np.mean(
                self.model.active_values(self.pressure))))
            series.recovery_factor.append(cumulative_oil
                                          / max(result.ooip, 1e-12) * 100.0)
            if output.record_well_rates:
                for name in self._producers:
                    result.well_oil_rate[name].append(
                        float(-rates.per_well_oil.get(name, 0.0)))
                    result.well_water_rate[name].append(
                        float(-rates.per_well_water.get(name, 0.0)))
                for name in self._injectors:
                    result.well_water_injection_rate[name].append(
                        float(max(rates.per_well_water.get(name, 0.0), 0.0)))

            # B4-B: bu addımda İŞLƏDİLƏN BHP qeyd olunur (yenilənməsi
            # artıq `_thp_outer_loop`-da addımın öz debitləri ilə olub).
            if self.thp_control.active and output.record_well_rates:
                self.thp_control.record(result)
            if self.bhp_limit.active and output.record_well_rates:
                self.bhp_limit.record(result)

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

    def _update_rate_shares(self) -> None:
        """RATE hədəfinin perforasiyalara payı — addımın ƏVVƏLİNDƏ, yığılmış
        vəziyyətin λ-sı ilə; Nyuton daxilində sabit qalır (bax
        `well_constraints.py`)."""
        if not self._rate_allocation:
            return
        assign_rate_shares(self.residual_assembler.wells,
                           self._connection_mobilities(self.state))

    def _connection_mobilities(self, state) -> list:
        """Bağlantıların lay həcmi mobillikləri — RATE payı və BHP limiti üçün."""
        assembler = self.residual_assembler
        return assembler.connection_mobilities(assembler.fluid_state(state))

    def _surface_factors(self, state) -> list:
        """Lay həcmi → səth debiti əmsalları (bax `SurfaceRateController`)."""
        assembler = self.residual_assembler
        return assembler.connection_surface_factors(assembler.fluid_state(state))

    def _surface_rate_loop(self, new_state, dt, newton_result):
        """SƏTH debiti hədəfi — addım DAXİLİNDƏ düzəliş (B7 addım 3).

        Əldə olunan səth debiti hədəfdən `SURFACE_RATE_TOLERANCE`-dan çox
        fərqlənirsə lay hədəfi miqyaslanır və addım eyni Δt ilə yenidən həll
        olunur. Təkrar yığılmasa əvvəlki həll saxlanılır; düzəldilmiş hədəf
        növbəti addımın proqnozu ilə onsuz da yenilənir.
        """
        controller = self.surface_rate
        if not controller.active:
            return new_state, newton_result
        for _ in range(MAX_SURFACE_ITERATIONS):
            if controller.correct(newton_result.rates) <= SURFACE_RATE_TOLERANCE:
                break
            retry = self.time_stepper.resolve_step(self.state, dt)
            if retry is None:
                LOG.warning("Səth debiti: düzəldilmiş hədəflə addım yenidən "
                            "yığılmadı — əvvəlki həll saxlanıldı.")
                break
            new_state, dt, newton_result = retry
            self.surface_rate_resolves += 1
        return new_state, newton_result

    def _bhp_limit_loop(self, new_state, dt, newton_result):
        """RATE quyularının BHP həddi — addım DAXİLİNDƏ rejim keçidi (B7 addım 2).

        Yığılmış həll yoxlanılır; bir quyu rejim dəyişibsə addım eyni Δt ilə
        yenidən həll olunur (`resolve_step` — tarixçə şişmir). Hər quyu bir
        addımda ən çox bir dəfə keçdiyi üçün dövr ən çox `quyu sayı + 1` dəfə
        fırlanır və SONUNCU qiymətləndirmə həmişə qəbul olunan həll üzərindədir.

        Təkrar yığılmasa əvvəlki həll saxlanılır; yeni rejim isə növbəti addım
        üçün qüvvədə qalır (limit növbəti addımda artıq pozulmur).
        """
        controller = self.bhp_limit
        if not controller.active:
            return new_state, newton_result
        controller.begin_step()
        for _ in range(controller.well_count + 1):
            if not controller.update(new_state):
                break
            retry = self.time_stepper.resolve_step(self.state, dt)
            if retry is None:
                LOG.warning("BHP limiti: rejim keçidindən sonra addım yenidən "
                            "yığılmadı — əvvəlki həll saxlanıldı, yeni rejim "
                            "növbəti addımdan tətbiq olunur.")
                break
            new_state, dt, newton_result = retry
            self.bhp_limit_resolves += 1
        return new_state, newton_result

    def _thp_outer_loop(self, new_state, dt, newton_result):
        """Yarı-implicit THP dövrəsi — addım DAXİLİNDƏ təkrarlama (Seans 25).

        Addım həll olunandan sonra BHP həmin addımın ÖZ debitləri ilə
        yenilənir; dəyişmə `OUTER_TOLERANCE_BAR`-dan böyükdürsə addım
        eyni Δt ilə YENİDƏN həll olunur. Beləliklə açıq birləşmənin bir
        addımlıq gecikməsi praktiki olaraq aradan qalxır — QALIQ və
        JAKOBİAN isə TOXUNULMUR (tam implicit birləşmə ⏳ hələ lazım deyil).

        Təkrar yığılmasa əvvəlki (artıq qəbul olunmuş) həll saxlanılır —
        yəni bu dövrə qaçışı HEÇ VAXT pisləşdirə bilməz.
        """
        if not self.thp_control.active:
            return new_state, newton_result
        for _ in range(MAX_OUTER_ITERATIONS):
            rates = newton_result.rates
            moved = self.thp_control.update(
                rates.per_well_oil, rates.per_well_water,
                getattr(rates, "per_well_gas", None),
                pressure=new_state.pressure)
            if moved <= OUTER_TOLERANCE_BAR:
                break
            retry = self.time_stepper.resolve_step(self.state, dt)
            if retry is None:
                break
            new_state, dt, newton_result = retry
            self.thp_outer_iterations += 1
        return new_state, newton_result

    def _record_snapshot(self, result: SimulationResult, time: float) -> None:
        shape = self.model.grid.shape
        result.snapshots.append(Snapshot(
            time=time,
            pressure=self.state.pressure.reshape(shape).copy(),
            water_saturation=self.state.water_saturation.reshape(shape).copy()))
