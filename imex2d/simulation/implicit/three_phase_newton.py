"""Üç fazalı Nyuton-Rafson döngəsi — A7, mərhələ 6d.

A6-dakı `NewtonSolver`-i (bax `newton.py`) güzgüləyir — eyni döngə
strukturu, eyni iki konvergensiya meyarı (CNV+MB), eyni Appleyard
kəsilməsi (chopping). YEGANƏ prinsipial fərq: hər addımdan sonra
`state.switch_variables(pvt)` çağırılır — bu, A7-nin əsas
mürəkkəbliyidir (bax `three_phase_state.py`): hüceyrə bir addımda
doymuş↔doymamış vəziyyət arasında keçə bilər, ona görə 3-cü
dəyişənin MƏNASI (Sg yoxsa Rs) hər iterasiyadan sonra yenidən
qiymətləndirilməlidir.

    x⁰ = xⁿ
    təkrarla:
        fluid = flüid_xassələri(x^k)
        R(x^k) hesabla
        yığılıbmı? -> switch_variables(x) sonra qaytar
        J(x^k) yığ
        δ = J⁻¹·(−R)
        x^{k+1} = chop(x^k + δ)              ← is_saturated DƏYİŞMİR
    yığılanda: switch_variables(x_yığılmış)   ← YALNIZ BURADA

DƏYİŞƏN KEÇİD NİYƏ YALNIZ SONDA

İlk versiyada `switch_variables()` HƏR iterasiyadan sonra çağırılırdı.
Ölçmə göstərdi ki, bu, OSİLYASİYAYA səbəb olur: sərhəddəki hüceyrə bir
iterasiyada doymuşa keçir, residual onu geri itələyir, növbəti
iterasiyada doymamışa qayıdır — CNV azalmaq əvəzinə tərəddüd edir və
Nyuton divergensiya edir (6-cı addımdan sonra müşahidə olundu).

Sənaye standartı (Aziz & Settari): hüceyrənin doyma vəziyyəti
(`is_saturated`) bir Nyuton addımının BÜTÜN iterasiyaları ərzində
SABİT saxlanılır — yalnız 3-cü dəyişənin (Sg və ya Rs) ƏDƏDİ
QİYMƏTİ yenilənir. Vəziyyət yalnız addım YIĞILDIQDAN SONRA
qiymətləndirilir və növbəti addım üçün istifadə olunur. Bu, "bir
addımda bir keçid" prinsipidir və osilyasiyanın qarşısını alır.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional

import numpy as np

from ...logging_setup import get_logger
from .newton import NewtonConfig, NewtonStatus
from .three_phase_residual import (ThreePhaseAccumulator, ThreePhaseFlux,
                                   ThreePhaseFluidState,
                                   ThreePhaseJacobianAssembler,
                                   ThreePhaseWellModel)
from .three_phase_state import ThreePhaseState

LOG = get_logger(__name__)


@dataclass
class ThreePhaseNewtonResult:
    status: NewtonStatus
    state: ThreePhaseState
    iterations: int
    history: List[float] = field(default_factory=list)
    material_balance: tuple = (0.0, 0.0, 0.0)
    fluid: Optional[ThreePhaseFluidState] = None
    rates: Optional[object] = None
    linear_iterations: int = 0

    @property
    def converged(self) -> bool:
        return self.status is NewtonStatus.CONVERGED


class ThreePhaseNewtonSolver:
    """Bir zaman addımını üç fazalı implicit həll edir.

    `NewtonConfig` A6-dan olduğu kimi işlədilir — sahələri (tolerans,
    addım hədləri) fazadan asılı deyil, ölçüsüzdür.
    """

    def __init__(self, model, relperm, pvt, linear_solver, grid,
                config: Optional[NewtonConfig] = None,
                endpoint_water_mobility: float = 0.35,
                capillary=None, gas_capillary=None):
        """`capillary` — su-neft Pc provider-i (`pcow`, `dpcow_dsw`);
        `gas_capillary` — qaz-neft Pc provider-i (`pcog`, `dpcog_dsg`), B5-b.
        İkisi də `None` ola bilər — onda kapilyar hədd YOXDUR (köhnə davranış).
        """
        self.model = model
        self.relperm = relperm
        self.pvt = pvt
        self.capillary = capillary
        self.gas_capillary = gas_capillary
        self.linear_solver = linear_solver
        self.grid = grid

        self.accumulator = ThreePhaseAccumulator(model, grid.pore_volume)
        self.flux = ThreePhaseFlux(model, grid)
        from ..well_model import PeacemanWellModel
        wells = PeacemanWellModel().build_connections(model)
        # Qaz vurulmasında (B7) son nöqtə kimi krg_end işlədilir —
        # relperm provider-i veribsə ondan, yoxsa mühafizəkar defolt.
        endpoint_gas = float(getattr(getattr(relperm, "gas", None),
                                     "krg_end", 0.8))
        self.well_model = ThreePhaseWellModel(model, wells,
                                              endpoint_water_mobility,
                                              endpoint_gas)
        self.jacobian = ThreePhaseJacobianAssembler(
            model, self.accumulator, self.flux, self.well_model, relperm, pvt)

        self.config = config or self._default_config()
        sw_min, sw_max = relperm.saturation_limits()
        relaxation = self.config.saturation_bound_relaxation
        self.sw_min = sw_min - relaxation
        self.sw_max = sw_max + relaxation
        self.physical_sw_min = sw_min
        self.physical_sw_max = sw_max
        self.sg_min, self.sg_max = relperm.gas_saturation_limits()

    @staticmethod
    def _default_config() -> NewtonConfig:
        """A6-nın iki fazalı defoltundan bir qədər boşaldılmış.

        Ölçmə göstərdi: üç fazalı sistemdə CNV mükəmməl yığılsa da
        (~10⁻⁵), material balans (mb_water) A6-nın sərt `1e-7`
        həddinə YAVAŞ-YAVAŞ yaxınlaşır — 3-cü tənliyin (qaz) əlavə
        etdiyi qoşulma bunu bir qədər ləngidir. 12 iterasiya bəzən
        kifayət etmir, bu da adaptiv addımlayıcını lazımsız yerə Δt-ni
        kəsməyə məcbur edir (real fiziki qeyri-sabitlik OLMADAN).

        Tolerans BOŞALDILMIR (`1e-7` qalır — dəqiqlikdən güzəştə
        getmirik), yalnız BÜDCƏ artırılır: 12-dən 25-ə. Ölçülmüş
        nəticə: əvvəllər 12-də uğursuz olan addım 19-da yığılırdı.
        """
        return NewtonConfig(max_iterations=25)

    # ═══════════════════════════════════════════════ flüid xassələri
    def build_fluid(self, state: ThreePhaseState) -> ThreePhaseFluidState:
        pressure = state.pressure
        sw, sg = state.water_saturation, state.gas_saturation
        rs = state.solution_gor(self.pvt)
        pc, dpc_dsw, pcog, dpcog_dsg = self._capillary_terms(state)
        bo, bo_p, bo_rs = self._oil_fvf(state, rs)
        mu_o, mu_o_p, mu_o_rs = self._oil_viscosity(state, rs)
        return ThreePhaseFluidState(
            mu_w=np.full(state.ncell, self.model.fluids.water_viscosity)
                if self.pvt is None else self._water_viscosity(pressure),
            mu_o=mu_o,
            mu_g=self.pvt.gas_viscosity(pressure),
            bw=self.pvt.water_fvf(pressure), bo=bo,
            bg=self.pvt.gas_fvf(pressure), rs=rs,
            krw=self.relperm.krw(sw),
            kro=self.relperm.kro_three_phase(sw, sg),
            krg=self.relperm.krg(sg), pc=pc,
            bo_p=bo_p, bo_rs=bo_rs,
            mu_o_p=mu_o_p, mu_o_rs=mu_o_rs,
            dpc_dsw=dpc_dsw, pcog=pcog, dpcog_dsg=dpcog_dsg)

    def _capillary_terms(self, state: ThreePhaseState):
        """(Pcow, ∂Pcow/∂Sw, Pcog, ∂Pcog/∂3-cü) — provider yoxdursa `None`.

        B5-b-dən ƏVVƏL burada `pc = None` SABİT yazılmışdı: servis Pc
        provider-ini mühərrikə ötürürdü, mühərrik isə onu Nyutona
        vermirdi — üç fazalı rejimdə kapilyar təzyiq səssizcə atılırdı.

        Pcog doymamış hüceyrədə Sg = 0-da oxunur və törəməsi SIFIRLANIR:
        orada 3-cü dəyişən Rs-dir, Pcog isə Rs-dən asılı deyil.
        """
        pc = dpc_dsw = pcog = dpcog_dsg = None
        sw = state.water_saturation
        if self.capillary is not None:
            pc = np.asarray(self.capillary.pcow(sw), float)
            dpc_dsw = np.asarray(self.capillary.dpcow_dsw(sw), float)
        if self.gas_capillary is not None:
            sg = state.gas_saturation
            pcog = np.asarray(self.gas_capillary.pcog(sg), float)
            dpcog_dsg = np.where(state.is_saturated,
                                 np.asarray(self.gas_capillary.dpcog_dsg(sg), float),
                                 0.0)
        return pc, dpc_dsw, pcog, dpcog_dsg

    def _oil_fvf(self, state: ThreePhaseState, rs: np.ndarray):
        """Bo və törəmələri — DOYMA VƏZİYYƏTİNƏ GÖRƏ (B3-B).

        Doymuş hüceyrə  → `Bo_sat(p)`, 3-cü dəyişən Sg-dir, `∂Bo/∂Sg = 0`.
        Doymamış hüceyrə → `Bo(p, Rs)` doymamış qoldan; `∂Bo/∂Rs ≠ 0`.

        Əvvəl HƏR hüceyrədə `oil_fvf(p)` (doymuş qol) işlədilirdi —
        doymamış hüceyrədə bu, termodinamik olaraq yanlışdır və
        Jakobianı kilidləyirdi. Tam izah: `BlackOilPVTProvider.
        oil_fvf_undersaturated`.

        Provider doymamış qolu dəstəkləmirsə (öz IPVTProvider
        implementasiyası) köhnə davranışa qayıdılır — mühərrik
        müqaviləsi genişlənmir.
        """
        pressure = state.pressure
        bo_sat = self.pvt.oil_fvf(pressure)
        branch = getattr(self.pvt, "oil_fvf_undersaturated", None)
        if branch is None:
            return bo_sat, None, None

        saturated = state.is_saturated
        bo_u = branch(pressure, rs)
        du_dp, du_drs = self.pvt.oil_fvf_undersaturated_derivatives(pressure, rs)
        sat_dp = self.pvt.oil_fvf_derivative(pressure)

        bo = np.where(saturated, bo_sat, bo_u)
        bo_p = np.where(saturated, sat_dp, du_dp)
        bo_rs = np.where(saturated, 0.0, du_drs)
        return bo, bo_p, bo_rs

    def _oil_viscosity(self, state: ThreePhaseState, rs: np.ndarray):
        """μo və törəmələri — DOYMA VƏZİYYƏTİNƏ GÖRƏ (G2b).

        `_oil_fvf`-in (B3-B) güzgüsüdür:
          doymuş hüceyrə  → `μo_sat(p)`, `∂μo/∂Sg = 0`;
          doymamış hüceyrə → `μo(p, Rs)` doymamış qoldan, `∂μo/∂Rs ≠ 0`.

        NİYƏ VACİBDİR: doymuş qoldan oxumaq doymamış neftin təzyiqlə
        QATILAŞMASINI itirir (SPE1: 0.51 → 0.74 cP, 45 %) və neft
        mobilliyini olduğundan YÜKSƏK göstərir.

        Provider doymamış qolu dəstəkləmirsə köhnə davranışa qayıdılır —
        mühərrik müqaviləsi genişlənmir.
        """
        pressure = state.pressure
        mu_sat = self.pvt.oil_viscosity(pressure)
        branch = getattr(self.pvt, "oil_viscosity_undersaturated", None)
        if branch is None:
            return mu_sat, None, None

        saturated = state.is_saturated
        mu_u = branch(pressure, rs)
        du_dp, du_drs = self.pvt.oil_viscosity_undersaturated_derivatives(
            pressure, rs)
        sat_dp = self.pvt.oil_viscosity_derivative(pressure)

        mu_o = np.where(saturated, mu_sat, mu_u)
        mu_o_p = np.where(saturated, sat_dp, du_dp)
        mu_o_rs = np.where(saturated, 0.0, du_drs)
        return mu_o, mu_o_p, mu_o_rs

    def _water_viscosity(self, pressure):
        analytic = getattr(self.pvt, "water_viscosity", None)
        if analytic is not None:
            return analytic(pressure)
        return np.full(pressure.shape, self.model.fluids.water_viscosity)

    def compute_residual(self, state: ThreePhaseState,
                         previous: ThreePhaseState,
                         previous_fluid: ThreePhaseFluidState, dt: float,
                         reference_upstream=None):
        fluid = self.build_fluid(state)
        n_w, n_o, n_g = self.accumulator.accumulation(state, fluid)
        n_w0, n_o0, n_g0 = self.accumulator.accumulation(previous, previous_fluid)
        influx_w, influx_o, influx_g = self.flux.net_influx(
            state, fluid, reference_upstream)
        # BHP aktiv/bağlı qərarı ADDIMIN ƏVVƏLKİ (yığılmış) təzyiqi ilə
        # verilir — bax `ThreePhaseWellModel.well_rates()`-in sənədləşməsi.
        # Cari iterasiyanın təzyiqi ilə versək, drawdown işarəsi
        # iterasiyalar arasında dəyişəndə Nyuton osilyasiya edir
        # (ölçülüb: RATE rejimində bu problem yoxdur, BHP-də var idi).
        rates = self.well_model.well_rates(state, fluid, previous.pressure)

        residual = np.empty(state.ncell * 3)
        residual[0::3] = (n_w - n_w0) / dt - influx_w - rates.water
        residual[1::3] = (n_o - n_o0) / dt - influx_o - rates.oil
        residual[2::3] = (n_g - n_g0) / dt - influx_g - rates.gas
        return residual, fluid, rates

    # ═══════════════════════════════════════════════ konvergensiya
    def convergence_measures(self, residual_vector: np.ndarray,
                             state: ThreePhaseState,
                             fluid: ThreePhaseFluidState, dt: float):
        """(CNV, MB_su, MB_neft, MB_qaz) — dördü də ölçüsüzdür."""
        pore_volume = self.accumulator.pore_volume_at(state.pressure)
        scale = np.maximum(pore_volume / dt, 1e-30)

        water = residual_vector[0::3]
        oil = residual_vector[1::3]
        gas = residual_vector[2::3]
        cnv = float(max((np.abs(water) / scale).max(),
                        (np.abs(oil) / scale).max(),
                        (np.abs(gas) / scale).max()))

        water_in_place, oil_in_place, gas_in_place = self.accumulator.accumulation(
            state, fluid)
        mb_water = float(abs(water.sum()) * dt
                         / max(water_in_place.sum(), 1e-30))
        mb_oil = float(abs(oil.sum()) * dt / max(oil_in_place.sum(), 1e-30))

        # Qaz üçün NİSBİ meyar TƏHLÜKƏLİDİR: Rs=0 (ölü neft, sərbəst
        # qaz yoxdur) ssenarisində gas_in_place.sum() demək olar
        # sıfırdır — `max(...,1e-30)`-a bölmə hətta çox kiçik ədədi
        # səs-küyü ASTRONOMİK mb_gas dəyərinə çevirir və yığılmanı
        # RİYAZİ CƏHƏTDƏN MÜMKÜNSÜZ edir (ölçülüb: bu, quyu BHP
        # sərhədi yaxınlığındakı "izahsız" yığılmama probleminin əsl
        # kökü idi — dt ölçüsündən ASILI OLMADAN baş verirdi, çünki
        # problem addımlamada deyil, bu bölmədə idi).
        #
        # Həll: gas_in_place ƏHƏMİYYƏTSİZ olanda (ümumi məsamə
        # həcminin çox kiçik bir hissəsi) mb_gas-ı avtomatik təmin
        # olunmuş sayırıq — CNV artıq bu hüceyrələrdə qazın əhəmiyyətsiz
        # olduğunu göstərir, əlavə nisbi yoxlamaya ehtiyac yoxdur.
        total_pore_volume = float(pore_volume.sum())
        gas_reference = gas_in_place.sum()
        if gas_reference < 1e-6 * max(total_pore_volume, 1e-30):
            mb_gas = 0.0
        else:
            mb_gas = float(abs(gas.sum()) * dt / gas_reference)
        return cnv, mb_water, mb_oil, mb_gas

    def _is_converged(self, cnv, mb_water, mb_oil, mb_gas) -> bool:
        tolerance = self.config.material_balance_tolerance
        return (cnv < self.config.cnv_tolerance and mb_water < tolerance
                and mb_oil < tolerance and mb_gas < tolerance)

    def _clamp_to_physical(self, state: ThreePhaseState) -> ThreePhaseState:
        """Ədədi zolaq yalnız həll müddətindədir — nəticə fiziki hədlərdədir."""
        water = np.clip(state.water_saturation, self.physical_sw_min,
                        self.physical_sw_max)
        gas = np.where(state.is_saturated,
                       np.clip(state.gas_saturation, self.sg_min, self.sg_max),
                       0.0)
        third = np.where(state.is_saturated, gas, state.third_variable)
        return ThreePhaseState(state.pressure, water, third, state.is_saturated)

    def _damped_update(self, state: ThreePhaseState, delta: np.ndarray,
                       previous: ThreePhaseState,
                       previous_fluid: ThreePhaseFluidState, dt: float,
                       current_cnv: float, pressure_limit, saturation_limit
                       ) -> ThreePhaseState:
        """Söndürülmüş (damped) Nyuton addımı — sənaye standartı.

        Tam addım CNV-ni PİSLƏŞDİRİRSƏ (adətən quyu BHP sərhədinə çox
        yaxın olan nöqtələrdə — well-switching diskretliyi Jakobianın
        yerli xəttiliyini pozur), addım ölçüsü YARIYA bölünüb yenidən
        sınanır. Bu, kök səbəbi diaqnoz etmədən Nyutonun qeyri-xətti
        çətin nöqtələrdə "sürüşməsini" dayandıran ümumi (heç bir
        modelə xas olmayan) bir texnikadır.

        Tam addımın özü hər zaman İLK sınanır (əksər addımlarda söndürmə
        lazım deyil, performans itkisi yoxdur) — yalnız pisləşmə halında
        geri çəkilir.
        """
        scale = 1.0
        trial = state
        for _ in range(6):          # 1, 0.5, 0.25, ... 1/32
            trial = state.updated(
                delta * scale, self.sw_min, self.sw_max,
                max_pressure_change=pressure_limit,
                max_saturation_change=saturation_limit)
            trial_fluid = self.build_fluid(trial)
            trial_residual, _, _ = self.compute_residual(
                trial, previous, previous_fluid, dt)
            trial_cnv, _, _, _ = self.convergence_measures(
                trial_residual, trial, trial_fluid, dt)
            # kiçik artımlara dözümlü oluruq (Nyuton adətən yerli-qeyri-
            # monoton ola bilər) — yalnız KƏSKİN pisləşmə söndürməni tələb edir
            if trial_cnv < current_cnv * 3.0 or scale < 0.05:
                return trial
            scale *= 0.5
        return trial

    # ═════════════════════════════════════════════════════════ döngə
    def solve(self, previous: ThreePhaseState, dt: float,
             initial_guess: Optional[ThreePhaseState] = None
             ) -> ThreePhaseNewtonResult:
        """`ISTİFADƏÇİYƏ VERİLƏN GARANTİ`: bu metod HEÇ VAXT istisna
        atmır — nə qədər qeyri-adi vəziyyət (NaN, sinqulyar Jakobian,
        platformaya xas scipy xətası) yaransa da, `ThreePhaseNewtonResult`
        (uğurlu ya uğursuz statusla) qaytarır. Kənar (`_solve_inner`)
        döngə daxilində gözlənilməyən bir şey baş versə, ən xaricdəki
        `except Exception` onu tutub `LINEAR_SOLVER_FAILED` kimi
        qaytarır — real çökmə əvəzinə (bax A7_PLAN.md, istifadəçi
        bildirişi: qaz fazası aktivkən tam proqram çökməsi).
        """
        try:
            return self._solve_inner(previous, dt, initial_guess)
        except Exception as error:
            LOG.exception("Üç fazalı Nyuton gözlənilməz istisna ilə "
                          "dayandı (%s) — TƏHLÜKƏSİZ uğursuzluqla qaytarılır",
                          type(error).__name__)
            fallback_state = (initial_guess or previous).copy()
            return ThreePhaseNewtonResult(
                NewtonStatus.LINEAR_SOLVER_FAILED, fallback_state, 0, [],
                (0.0, 0.0, 0.0), None, None, 0)

    def _solve_inner(self, previous: ThreePhaseState, dt: float,
                     initial_guess: Optional[ThreePhaseState] = None
                     ) -> ThreePhaseNewtonResult:
        config = self.config
        state = (initial_guess or previous).copy()
        previous_fluid = self.build_fluid(previous)
        # Axının upstream seçimi ADDIMIN ƏVVƏLKİ vəziyyətinə görə BİR
        # DƏFƏ dondurulur (bax `ThreePhaseFlux.face_fluxes`-in
        # sənədləşməsi) — potensial fərqi sıfıra yaxın olan üzlərdə
        # bu seçimin iterasiyalar arasında dəyişməsi PERİODİK
        # OSİLYASİYAYA səbəb olurdu (ölçülüb: CNV dövr-3 dövrədə
        # ilişib qalırdı, heç bir dt ölçüsü kömək etmirdi).
        reference_upstream = self.flux.upstream_masks(previous, previous_fluid)

        history: List[float] = []
        first_cnv = None
        linear_iterations = 0
        fluid = rates = None
        material_balance = (0.0, 0.0, 0.0)

        for iteration in range(config.max_iterations + 1):
            residual_vector, fluid, rates = self.compute_residual(
                state, previous, previous_fluid, dt, reference_upstream)
            cnv, mb_w, mb_o, mb_g = self.convergence_measures(
                residual_vector, state, fluid, dt)
            history.append(cnv)
            material_balance = (mb_w, mb_o, mb_g)

            if self._is_converged(cnv, mb_w, mb_o, mb_g):
                converged_state = self._clamp_to_physical(state).switch_variables(
                    self.pvt)
                return ThreePhaseNewtonResult(
                    NewtonStatus.CONVERGED, converged_state,
                    iteration, history, material_balance, fluid, rates,
                    linear_iterations)

            if first_cnv is None:
                first_cnv = max(cnv, 1e-30)
            elif cnv > first_cnv * config.divergence_factor:
                LOG.debug("Üç fazalı Nyuton divergensiya etdi: CNV %.3e -> %.3e",
                         first_cnv, cnv)
                return ThreePhaseNewtonResult(
                    NewtonStatus.DIVERGED, state, iteration, history,
                    material_balance, fluid, rates, linear_iterations)

            if iteration == config.max_iterations:
                break

            # NaN/Inf qorunması — QALIQ ÖZÜ artıq zədəlidirsə (güclü
            # osilyasiyadan sonra baş verə bilər), Jakobian yığımına və
            # xətti həllediciyə HEÇ VERİLMİR. NaN/Inf-li seyrək matrisi
            # native SuperLU/UMFPACK kitabxanasına ötürmək platformadan
            # asılı davranışa (bəzi sistemlərdə səssiz çökmə) səbəb ola
            # bilər — Python-un öz istisna tutma mexanizmi bunu HƏMİŞƏ
            # tuta bilmir. Əvvəlcədən yoxlamaq buna ehtiyacı aradan
            # qaldırır.
            if not np.all(np.isfinite(residual_vector)):
                LOG.debug("Qalıq NaN/Inf ehtiva edir — xətti həll keçilir")
                return ThreePhaseNewtonResult(
                    NewtonStatus.LINEAR_SOLVER_FAILED, state, iteration,
                    history, material_balance, fluid, rates, linear_iterations)

            try:
                matrix = self.jacobian.assemble(
                    state, fluid, dt, previous.pressure, reference_upstream)
                if not np.all(np.isfinite(matrix.data)):
                    raise ValueError("Jakobian matrisi NaN/Inf ehtiva edir")
                delta = self.linear_solver.solve(matrix, -residual_vector)
                linear_iterations += self.linear_solver.last_iterations
            except Exception as error:
                # QƏSDƏN GENİŞ TUTULUR: scipy-nin sparse həllediciləri
                # platformadan asılı olaraq müxtəlif istisna növləri
                # ata bilər (RuntimeError, LinAlgError, MemoryError və
                # s.) — yalnız FloatingPointError/ValueError kifayət
                # etmirdi (istifadəçi tərəfindən real çökmə bildirilib,
                # bax A7_PLAN.md). Nyuton üçün "xətti həll uğursuz oldu"
                # HƏMİŞƏ təhlükəsiz, gözlənilən bir haldır — proses heç
                # vaxt bundan ÇÖKMƏMƏLİDİR.
                LOG.debug("Xətti həlledici uğursuz (%s): %s",
                         type(error).__name__, error)
                return ThreePhaseNewtonResult(
                    NewtonStatus.LINEAR_SOLVER_FAILED, state, iteration,
                    history, material_balance, fluid, rates, linear_iterations)

            if not np.all(np.isfinite(delta)):
                return ThreePhaseNewtonResult(
                    NewtonStatus.LINEAR_SOLVER_FAILED, state, iteration,
                    history, material_balance, fluid, rates, linear_iterations)

            pressure_limit, saturation_limit = config.limits_at(iteration)
            state = self._damped_update(
                state, delta, previous, previous_fluid, dt, cnv,
                pressure_limit, saturation_limit)
            # Dəyişən keçid BURADA ÇAĞIRILMIR (bax sinif sənədləşməsi) —
            # `is_saturated` bu addımın bütün iterasiyaları boyu sabit
            # qalır, yalnız yığılmadan SONRA yenilənir.

        return ThreePhaseNewtonResult(
            NewtonStatus.MAX_ITERATIONS, state, config.max_iterations,
            history, material_balance, fluid, rates, linear_iterations)
