"""ThreePhaseSimulationEngine — A7, mərhələ 6d (son hissə).

`ISimulationEngine` implementasiyası — `FullyImplicitEngine`-i
güzgüləyir, eyni interfeys, eyni `SimulationResult`. Fərq: üç fazalı
Nyuton döngəsi (`ThreePhaseNewtonSolver`) və qaz seriyaları
(`gas_rate`, `cumulative_gas`, `gas_oil_ratio`).

BU MÜHƏRRİK YALNIZ PVT PROVIDER `has_gas_phase() == True` OLANDA
İŞLƏDİLMƏLİDİR. Yoxdursa, `FullyImplicitEngine` (iki fazalı) işlədilir
— seçim `application` qatında edilir, bax `simulation_service.py`.

`AdaptiveTimeStepper` (A6) DƏYİŞDİRİLMƏDƏN yenidən işlədilir —
`ThreePhaseNewtonSolver.solve(state, dt)` eyni imzaya malikdir və
`ThreePhaseState.water_saturation` mövcuddur, ona görə adaptiv Δt
məntiqi (növbəti addımın ölçüsü, doyumluluq dəyişikliyi limiti)
heç bir dəyişiklik olmadan üç fazaya tətbiq olunur.
"""

from __future__ import annotations
from typing import Optional

import numpy as np

from ...application.config import SimulationConfig
from ...domain.reservoir_model import ReservoirModel
from ...interfaces.providers import IInitializationProvider, IPVTProvider
from ...interfaces.services import (IProgressReporter, ISimulationEngine,
                                    NullProgressReporter)
from ...logging_setup import get_logger
from ..discretization import TwoPointFluxDiscretization
from ..results import SimulationResult, Snapshot
from .newton import NewtonConfig
from .three_phase_newton import ThreePhaseNewtonSolver
from ..wellbore.thp_control import (MAX_OUTER_ITERATIONS,
                                    OUTER_TOLERANCE_BAR, ThpController)
from ..well_constraints import (MAX_SURFACE_ITERATIONS, SURFACE_RATE_TOLERANCE,
                                BhpLimitController, SurfaceRateController,
                                assign_rate_shares, needs_rate_allocation)
from .three_phase_state import ThreePhaseState
from .time_stepping import AdaptiveTimeStepConfig, AdaptiveTimeStepper

LOG = get_logger(__name__)


class ThreePhaseSimulationEngine(ISimulationEngine):

    def __init__(self, model: ReservoirModel, config: SimulationConfig,
                relperm, pvt: IPVTProvider,
                linear_solver=None, capillary=None,
                initialization: Optional[IInitializationProvider] = None,
                newton_config: Optional[NewtonConfig] = None,
                time_step_config: Optional[AdaptiveTimeStepConfig] = None,
                flux_discretization=None):
        """`linear_solver`, `capillary` — A6-dakı mühərriklərlə eyni
        çağırış imzasını qorumaq üçün qəbul edilir
        (`SimulationService.create_engine()` bütün mühərrikləri eyni
        açar sözlərlə qurur), LAKİN `linear_solver` FAKTİKİ İŞLƏDİLMİR
        — `SimulationService`-in defolt `ScipyCgIluSolver`-i A6-nın
        2×2 CPR/ILU blok strukturu üçün tənzimlənib, 3×3 struktura
        uyğun deyil (CPR-in 3×3 genişlənməsi qəsdən təxirə salınıb,
        bax `A7_PLAN.md`). Bu mühərrik həmişə öz `NewtonLinearSolver`
        nüsxəsini (kiçik-orta modellərdə birbaşa `splu`) işlədir.
        `capillary` (su-neft Pc) B5-b-dən ETİBARƏN İŞLƏDİLİR — Nyutona
        ötürülür. Qaz-neft Pc-si isə MODELDƏN qurulur
        (`gas_capillary_parameters`): mühərrik imzası bütün mühərriklər
        üçün ortaqdır və yalnız üç fazalı mühərrikin oxuduğu provider
        ora əlavə edilmir.
        """
        if not pvt.has_gas_phase():
            raise ValueError(
                "ThreePhaseSimulationEngine yalnız qaz xassələri olan PVT "
                "cədvəli ilə işləyir (build_pvt_table(..., include_gas=True)).")
        # B1 (MPFA-O seçimi) `SimulationService.create_engine()`-ə
        # `flux_discretization` açar sözü əlavə etdi və bütün mühərriklər
        # onu EYNİ imza ilə alır. Üç fazalı qalıq/Jakobian hələ YALNIZ
        # TPFA üçün yazılıb (çoxnöqtəli stensil üç fazalı Jakobiana
        # qoşulmayıb), ona görə burada AÇIQ RƏDD edilir — səssizcə
        # TPFA-ya keçmək istifadəçinin seçdiyindən BAŞQA bir hesab
        # aparmaq olardı.
        if flux_discretization is not None and getattr(
                flux_discretization, "supports_multipoint_stencil",
                lambda: False)():
            raise NotImplementedError(
                "Üç fazalı (qazlı) mühərrik çoxnöqtəli diskretizasiya "
                "(MPFA-O) ilə HƏLƏ İŞLƏMİR — TPFA seçin.")
        self.model = model
        self.config = config
        self.relperm = relperm
        self.pvt = pvt
        self.capillary = capillary
        self.initialization = initialization

        self.grid = TwoPointFluxDiscretization().build(model)
        self.newton = ThreePhaseNewtonSolver(
            model, relperm, pvt, self._linear_solver(),
            self.grid, config=newton_config,
            capillary=capillary,
            gas_capillary=self._gas_capillary(model, relperm))
        self.time_stepper = AdaptiveTimeStepper(
            self.newton, time_step_config or self._time_config(config))

        self._producers = sorted({c.well_name for c in self.newton.well_model.wells
                                  if not c.is_injector})
        self._injectors = sorted({c.well_name for c in self.newton.well_model.wells
                                  if c.is_injector})
        #: Çox perforasiyalı RATE quyusu varmı (Seans 27)
        self._rate_allocation = needs_rate_allocation(self.newton.well_model.wells)
        # B4-B: THP quyuları — bağlantı hədəfi addım-addım yenilənir
        self.thp_control = ThpController(model, self.newton.well_model.wells,
                                         pvt=pvt, fluids=model.fluids)
        # B7 addım 3: SƏTH bazalı RATE hədəfi (bax `_surface_rate_loop`)
        self.surface_rate = SurfaceRateController(self.newton.well_model.wells,
                                                  self._surface_factors)
        #: səth debiti düzəlişinə görə addımın ƏLAVƏ həllərinin sayı
        self.surface_rate_resolves = 0
        # B7 addım 2: RATE quyularının BHP həddi (bax `_bhp_limit_loop`)
        self.bhp_limit = BhpLimitController(self.newton.well_model.wells,
                                            self._connection_mobilities,
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
        # IMPES/FullyImplicit ilə eyni atributlar — UI/testlər üçün
        self.pressure = self.state.pressure
        self.sw = self.state.water_saturation
        self.sg = self.state.gas_saturation

    @staticmethod
    def _linear_solver():
        from .linear import NewtonLinearSolver
        return NewtonLinearSolver()

    @staticmethod
    def _gas_capillary(model, relperm):
        """Modeldəki Pcog parametrlərindən provider — söndürülübsə `None`.

        Sgc/Sorg və Swc RELPERM provider-indən götürülür ki, kapilyar
        əyri ilə nisbi keçiricilik əyrisi EYNİ son nöqtələrə istinad
        etsin (fərqli dəyərlər iki uyğunsuz əyri yaradardı).
        """
        params = getattr(model, "gas_capillary_parameters", None)
        if params is None or not params.enabled:
            return None
        from ...domain.scal import GasCoreyParameters
        from ..capillary import BrooksCoreyGasCapillaryProvider
        gas_scal = (getattr(relperm, "gas", None)
                    or getattr(model, "gas_scal_parameters", None)
                    or GasCoreyParameters())
        swc = float(getattr(relperm, "swc", model.scal_parameters.swc))
        return BrooksCoreyGasCapillaryProvider(params, gas_scal, swc)

    @staticmethod
    def _time_config(config: SimulationConfig) -> AdaptiveTimeStepConfig:
        """Adaptiv addım parametrləri — `max_dt` OLDUĞU KİMİ ötürülür.

        ƏVVƏL burada `max(stepping.max_dt, 30.0)` yazılırdı: istifadəçi
        interfeysdə 20 gün desə də, mühərrik səssizcə 30-a keçirdi.
        Ölçüldü — 5 və 20 gün TAM EYNİ nəticə verirdi (31 addım,
        maks Δt 30.00), çünki hər ikisi eyni həddə qaldırılırdı.

        Bu, iki fazalı mühərrikdə ARTIQ tapılıb düzəldilmiş səhvin
        eynisidir (bax `engine.py::_time_config` sənədi). Üç fazalı
        yol A7 bərpasında (B2) git tarixçəsindən qaytarıldığı üçün
        köhnə kodu özü ilə geri gətirmişdi.

        ⏳ Qeyd: iki fazalı variant `soft_failure_*` toleranslarını da
        verir, bu isə vermir. Bu, AYRI fərqdir — yığılma davranışına
        toxunur, ona görə bu düzəlişdə TOXUNULMADI.
        """
        stepping = config.time_stepping
        return AdaptiveTimeStepConfig(
            initial_dt=stepping.initial_dt, min_dt=stepping.min_dt,
            max_dt=stepping.max_dt,
            growth_factor=stepping.growth_factor + 0.35)

    # ─────────────────────────────────────────────────────────── ilkin
    def _initial_state(self) -> ThreePhaseState:
        ic = self.model.initial_conditions
        n = self.model.ncell

        if self.initialization is not None:
            initial = self.initialization.initialize(self.model)
            pressure = np.asarray(initial.pressure, float).copy()
            water = np.asarray(initial.water_saturation, float).copy()
            if initial.has_gas:
                gas = np.asarray(initial.gas_saturation, float)
                is_saturated = gas > 1e-9
                # Qaz papağının ALTINDA sərbəst qaz yoxdur, amma neft
                # ÖLÜ DEYİL — orada da həll olmuş qaz var (B4b).
                third = np.where(is_saturated, gas,
                                 self._initial_solution_gor(pressure, ic))
            else:
                is_saturated = np.zeros(n, dtype=bool)
                third = self._initial_solution_gor(pressure, ic)
        else:
            pressure = np.full(n, ic.datum_pressure)
            water = np.full(n, ic.water_saturation)
            is_saturated = np.zeros(n, dtype=bool)
            third = self._initial_solution_gor(pressure, ic)

        sw_min, sw_max = self.relperm.saturation_limits()
        water = np.clip(water, sw_min, sw_max)
        return ThreePhaseState(pressure, water, third, is_saturated)

    def _initial_solution_gor(self, pressure: np.ndarray,
                              ic) -> np.ndarray:
        """İLKİN HƏLL OLMUŞ QAZ (Rs) — doymamış hüceyrələr üçün (B4b).

        ƏVVƏL BURADA NƏ VAR İDİ: `np.zeros(n)`, yəni Rs = 0 — "ölü
        neft". Şərh belə əsaslandırırdı: "domain modelində ayrıca ilkin
        Rs sahəsi yoxdur; bu, ən mühafizəkar seçimdir — xəyali qaz
        yaratmır."

        NİYƏ DƏYİŞDİ. Ölçüldü (bax `ISH_HESABATI.md` → Seans 6):
        Rs = 0 olanda OGIP = 0 çıxır və təzyiq doyma təzyiqindən aşağı
        düşsə BELƏ qaz ayrılmır — ayrılacaq həll olmuş qaz YOXDUR.
        Yəni üç fazalı mühərrik faktiki olaraq iki fazalı işləyirdi və
        sahibkarın "P < Psat olduqda qazın ayrılması" tələbi (M4)
        ÖDƏNMİRDİ.

        İNDİ: `InitialConditions.solution_gor` sahəsi əlavə olundu.
        `None` (defolt) olanda dəyər PVT cədvəlindən çıxarılır:

            Rs = Rs_sat(min(P_hüceyrə, Pb))

        Cədvəldə Rs onsuz da Pb-dən yuxarı SABİTDİR (korrelyasiya belə
        qurur), ona görə `pvt.solution_gor(P)` düz həmin dəyəri verir —
        ayrıca `min()` lazım deyil, lakin AÇIQ yazılır ki, qeyri-standart
        (məs. Eclipse-dən idxal olunmuş) cədvəldə də düzgün işləsin.

        Bu, sənaye standartıdır: neft öz doyma təzyiqinə uyğun qədər
        qaz saxlayır (Eclipse `EQUIL`/`RSVD` ilə eyni məntiq).
        """
        override = getattr(ic, "solution_gor", None)
        if override is not None:
            return np.full(pressure.size, float(override))

        bubble_point = getattr(self.model.pvt_table, "bubble_point", 0.0) or 0.0
        reference = (np.minimum(pressure, bubble_point) if bubble_point > 0.0
                     else pressure)
        return np.asarray(self.pvt.solution_gor(reference), float).copy()

    def original_oil_in_place(self) -> float:
        fluid = self.newton.build_fluid(self.state)
        _, oil, _ = self.newton.accumulator.accumulation(self.state, fluid)
        return float(oil.sum())

    def original_gas_in_place(self) -> float:
        fluid = self.newton.build_fluid(self.state)
        _, _, gas = self.newton.accumulator.accumulation(self.state, fluid)
        return float(gas.sum())

    # ══════════════════════════════════════════════════════════════ run
    def run(self, reporter: Optional[IProgressReporter] = None) -> SimulationResult:
        """`ISTİFADƏÇİYƏ VERİLƏN GARANTİ`: bu metod HEÇ VAXT istisna
        atmır — hər zaman `SimulationResult` qaytarır (`converged`
        True ya da False). A7 sınaq statusundadır və gözlənilməz
        vəziyyətlər (NaN, sinqulyar Jakobian) yarana bilər — bunlar
        HEÇ VAXT tam proqram çökməsinə səbəb olmamalıdır (istifadəçi
        bildirişi: qaz fazası aktivkən real çökmə müşahidə edilib,
        bax A7_PLAN.md). `ThreePhaseNewtonSolver.solve()` artıq öz
        səviyyəsində qorunur (bax onun sənədləşməsi) — bu, əlavə,
        SON mühafizə təbəqəsidir: akkumulyasiya/snapshot/seriya
        yazma kimi digər hissələrdə gözlənilməz bir şey baş versə
        belə, yenə TƏHLÜKƏSİZ nəticə qaytarılır.
        """
        try:
            return self._run_inner(reporter)
        except Exception as error:
            LOG.exception("Üç fazalı simulyasiya gözlənilməz istisna ilə "
                          "dayandı (%s) — TƏHLÜKƏSİZ nəticə qaytarılır",
                          type(error).__name__)
            result = SimulationResult(model_name=self.model.name,
                                      grid_shape=self.model.grid.shape)
            result.converged = False
            result.message = (f"Gözlənilməz xəta ({type(error).__name__}): "
                             f"{error}")
            return result

    def _run_inner(self, reporter: Optional[IProgressReporter] = None
                   ) -> SimulationResult:
        reporter = reporter or NullProgressReporter()
        config = self.config
        output = config.output

        result = SimulationResult(model_name=self.model.name,
                                  grid_shape=self.model.grid.shape)
        result.ooip = self.original_oil_in_place()
        result.ogip = self.original_gas_in_place()
        result.well_oil_rate = {name: [] for name in self._producers}
        result.well_water_rate = {name: [] for name in self._producers}
        result.well_gas_rate = {name: [] for name in self._producers}
        series = result.series

        snapshot_interval = max(config.end_time / max(output.snapshot_count, 1),
                                1e-9)
        next_snapshot = 0.0
        self._record_snapshot(result, 0.0)
        next_snapshot += snapshot_interval

        time = 0.0
        steps = 0
        cumulative_oil = cumulative_water = cumulative_gas = 0.0

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
            self.sg = self.state.gas_saturation
            time += dt
            steps += 1

            rates = newton_result.rates
            oil_rate = float(-min(rates.oil.sum(), 0.0))
            water_rate = float(-min(rates.water[rates.water < 0].sum(), 0.0))
            # QAZ QUYULAR ÜZRƏ ayrılır. Əvvəl `-min(rates.gas.sum(), 0)` idi:
            # hüceyrə cəmi VURULAN (+) və HASİL OLUNAN (−) qazı qarışdırırdı,
            # vurucu hasilatdan güclü olanda cəm müsbət çıxır və `min` onu
            # SIFIRLAYIRDI — ölçüldü: qaz vuran modeldə 1500 gün boyu sahə qaz
            # debiti və GOR = 0, halbuki istismarçının öz GOR-u 109+ idi.
            gas_rate = float(max(0.0, -sum(rates.per_well_gas.get(name, 0.0)
                                           for name in self._producers)))
            gas_injection = float(max(0.0, sum(rates.per_well_gas.get(name, 0.0)
                                               for name in self._injectors)))
            injection = float(max(rates.water[rates.water > 0].sum(), 0.0))
            cumulative_oil += oil_rate * dt
            cumulative_water += water_rate * dt
            cumulative_gas += gas_rate * dt

            series.time.append(time)
            series.oil_rate.append(oil_rate)
            series.water_rate.append(water_rate)
            series.gas_rate.append(gas_rate)
            series.water_injection_rate.append(injection)
            series.gas_injection_rate.append(gas_injection)
            series.cumulative_oil.append(cumulative_oil)
            series.cumulative_water.append(cumulative_water)
            series.cumulative_gas.append(cumulative_gas)
            series.water_cut.append(water_rate / max(oil_rate + water_rate, 1e-12)
                                    * 100.0)
            series.gas_oil_ratio.append(gas_rate / max(oil_rate, 1e-12))
            series.average_pressure.append(float(np.mean(self.pressure)))
            series.recovery_factor.append(cumulative_oil
                                          / max(result.ooip, 1e-12) * 100.0)
            if output.record_well_rates:
                for name in self._producers:
                    result.well_oil_rate[name].append(
                        float(-rates.per_well_oil.get(name, 0.0)))
                    result.well_water_rate[name].append(
                        float(-rates.per_well_water.get(name, 0.0)))
                    result.well_gas_rate[name].append(
                        float(-rates.per_well_gas.get(name, 0.0)))

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
                          f"dt = {dt:6.2f} | Nyuton {newton_result.iterations} | "
                          f"GOR = {series.gas_oil_ratio[-1]:.1f}")
                if not reporter.report(time / config.end_time * 100.0, message):
                    result.message = "İstifadəçi tərəfindən dayandırıldı."
                    break

        result.steps = steps
        # Uğurla bitəndə YEKUN MESAJ — iki fazalı mühərriklə eyni format
        # (`implicit/engine.py`). Bərpa olunan A7 kodunda bu blok yox idi:
        # qaz aktiv olanda istifadəçi status sətrində və jurnalda BOŞ
        # mesaj görürdü (ölçülüb, B2). Fizika DƏYİŞMİR.
        if not result.message:
            statistics = self.time_stepper.summary()
            result.message = (
                f"Tamamlandı (3 fazalı): {steps} addım, t = {time:.1f} gün "
                f"(orta Δt = {statistics.get('orta Δt', 0.0):.1f} gün, "
                f"orta {statistics.get('orta iterasiya', 0.0):.1f} Nyuton "
                f"iterasiyası, {statistics.get('təkrar', 0)} təkrar).")
        LOG.info("%s  RF = %.2f %%", result.message,
                 result.final_recovery_factor)
        return result

    def _update_rate_shares(self) -> None:
        """RATE hədəfinin perforasiyalara payı — addımın ƏVVƏLİNDƏ, yığılmış
        vəziyyətin λ-sı ilə (bax `well_constraints.py`)."""
        if not self._rate_allocation:
            return
        assign_rate_shares(self.newton.well_model.wells,
                           self._connection_mobilities(self.state))

    def _connection_mobilities(self, state) -> list:
        """Bağlantıların lay həcmi mobillikləri — RATE payı və BHP limiti üçün."""
        return self.newton.well_model.connection_mobilities(
            self.newton.build_fluid(state))

    def _surface_factors(self, state) -> list:
        """Lay həcmi → səth debiti əmsalları (bax `SurfaceRateController`)."""
        return self.newton.well_model.connection_surface_factors(
            self.newton.build_fluid(state))

    def _surface_rate_loop(self, new_state, dt, newton_result):
        """SƏTH debiti hədəfi — addım DAXİLİNDƏ düzəliş (B7 addım 3).

        `FullyImplicitEngine._surface_rate_loop` ilə eynidir (bax onun sənədi).
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

        `FullyImplicitEngine._bhp_limit_loop` ilə eynidir (bax onun sənədi).
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
        result.snapshots.append(Snapshot(
            time=time, pressure=self.state.pressure.copy(),
            water_saturation=self.state.water_saturation.copy(),
            gas_saturation=self.state.gas_saturation.copy()))
