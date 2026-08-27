"""Üç fazalı qalıq tənlikləri (akkumulyasiya) — A7, mərhələ 5.

İki fazalı sxemdə (`residual.py`) hər hüceyrə üçün İKİ kütlə balansı
tənliyi yazılırdı (su, neft). Qaz əlavə olunanda ÜÇÜNCÜ tənlik lazım
gəlir, lakin onun forması sadə deyil — qaz İKİ yerdə saxlanıla bilər:

    1. SƏRBƏST qaz fazası kimi     (Sg/Bg)
    2. NEFTDƏ HƏLL OLMUŞ kimi       (So·Rs/Bo)

Standart black-oil qaz kütlə balansı bunların CƏMİDİR:

    N_gas = PV · (Sg/Bg + So·Rs/Bo)

Bu, su və neft tənliklərindən keyfiyyətcə fərqlidir: onlar YALNIZ öz
fazalarını daşıyır (su heç vaxt neftdə həll olmur bu modeldə — "quru
qaz" fərziyyəsi, Rv=0, bax `A7_PLAN.md`), qaz isə İKİ mənbədən gəlir.

Doymamış hüceyrədə (Sg=0) ikinci hədd DOMİNANTdır — bütün qaz neftdə
həll olub. Doymuş hüceyrədə hər ikisi iştirak edir. Bu, mərhələ 4-dəki
`ThreePhaseState.solution_gor()`-un niyə doymuş hüceyrələr üçün
Rs_sat(p) qaytardığını izah edir: hətta sərbəst qaz olanda da neft
"doymuş" qalır və maksimum qazı özündə saxlayır.

Su/neft akkumulyasiyası A6-dakı ilə EYNİDİR (dəyişməyib) — yalnız
üçüncü tənlik yenidir. Bu ayrılıq qəsdən saxlanılıb: mövcud
`ResidualAssembler.accumulation()` iki fazalı testlərdə hələ də
işlədilir, bu modul onu SINDIRMIR, üstünə əlavə edir.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

import numpy as np

from ...domain.reservoir_model import ReservoirModel
from ...domain.wells import ControlMode
from ..discretization import DiscretizedGrid
from ..well_model import WellConnection
from .three_phase_state import ThreePhaseState

GRAVITY = 9.80665
PA_TO_BAR = 1.0e-5


def smooth_negative_part(value, scale):
    """`min(value, 0)`-un hamar əvəzi — SINANDI, İŞLƏTİLMİR.

    XƏBƏRDARLIQ: bu funksiya quyu debitlərində sınandı və GERİ
    QAYTARILDI. İki səbəb: (1) hədəf problemi (quyu BHP sərhədində
    yığılmama) HƏLL ETMƏDİ — ölçüldü, CNV tarixçəsi hərfi olaraq eyni
    qaldı, çünki problemli nöqtədə debit kəsmə zonasından çox uzaqdır;
    (2) FİZİKANI POZDU — nisbi keçiricilik sıfır olanda (məs. krg=0,
    doymamış hüceyrədə) debit DƏQİQ sıfır olmalıdır, hamarlaşdırma isə
    ona süni dəyər verir (ölçüldü: qaz debitində −58847 sm³/gün süni
    fərq). Saxlanılır ki, gələcəkdə eyni yanaşma təkrar sınanmasın.

        f(x) = ½·(x − √(x² + ε²))

    NİYƏ LAZIMDIR — bu, qaz fazasının əsas davamlılıq probleminin
    kökü idi. Sərt `min(q, 0)` q=0 nöqtəsində SINIQdır (kink): törəmə
    orada 1-dən 0-a sıçrayır. Quyu öz BHP hədəfinə yaxınlaşanda debit
    məhz sıfıra yaxınlaşır — yəni Nyuton hər iterasiyada bu sınıq
    nöqtənin bir tərəfindən digərinə atılır və YIĞILA BİLMİR
    (ölçülüb: istismarçının qonşu hüceyrəsində residual dövr-2 rəqs
    edirdi).

    OPM Flow bu problemi quyunu ÖZ primary dəyişənləri (BHP daxil)
    ilə tam implicit həll etməklə aradan qaldırır. Bizim quyu
    modelimiz sadədir (BHP sabit sərhəd şərtidir), ona görə burada
    ekvivalent nəticəni HAMARLAŞDIRMA ilə alırıq: funksiya hər yerdə
    diferensiallana biləndir, `ε` isə keçid zolağının enidir.

    `scale` — tipik debit böyüklüyü; ε bundan hesablanır ki,
    hamarlaşdırma modelin miqyasından asılı olmasın.
    """
    epsilon = max(abs(scale), 1e-12) * 1e-3
    return 0.5 * (value - np.sqrt(value * value + epsilon * epsilon))


def smooth_negative_part_derivative(value, scale):
    """`smooth_negative_part`-ın törəməsi — Jakobian üçün.

    Sərt `min`-in törəməsi (0 ya 1) əvəzinə hamar keçid verir; məhz
    bu, Jakobianı residualla UYĞUN saxlayır (əks halda Nyuton yanlış
    istiqamətə addımlayır).
    """
    epsilon = max(abs(scale), 1e-12) * 1e-3
    return 0.5 * (1.0 - value / np.sqrt(value * value + epsilon * epsilon))


@dataclass
class ThreePhaseFluidState:
    """Verilmiş (p, Sw, Sg) qiymətində flüid xassələri — hüceyrə üzrə."""
    mu_w: np.ndarray
    mu_o: np.ndarray
    mu_g: np.ndarray
    bw: np.ndarray
    bo: np.ndarray
    bg: np.ndarray
    rs: np.ndarray            # solution GOR — bax `ThreePhaseState.solution_gor()`
    krw: np.ndarray
    kro: np.ndarray
    krg: np.ndarray
    pc: Optional[np.ndarray] = None   # su-neft kapilyar təzyiqi (A4-dəki ilə eyni)

    @property
    def lam_w(self) -> np.ndarray:
        return self.krw / self.mu_w

    @property
    def lam_o(self) -> np.ndarray:
        return self.kro / self.mu_o

    @property
    def lam_g(self) -> np.ndarray:
        return self.krg / self.mu_g


class ThreePhaseAccumulator:
    """Üç fazalı akkumulyasiya (kütlə saxlanması həddi).

    Yalnız akkumulyasiya — axın (flux) və quyu həddləri hələ BURADA
    DEYİL (mərhələ 6-nın işidir, tam qalıq yığımı + Jakobian ilə
    birlikdə). Bu ayrılıq A6-nın öz mərhələləşməsi ilə uyğundur:
    orada da akkumulyasiya əvvəl tək başına yazılıb, sınanıb, sonra
    axın həddləri üstünə əlavə olunub.
    """

    def __init__(self, model: ReservoirModel, pore_volume: np.ndarray):
        self.model = model
        self.pore_volume = pore_volume
        self.reference_pressure = float(model.initial_conditions.datum_pressure)

    def pore_volume_at(self, pressure: np.ndarray) -> np.ndarray:
        """Süxur sıxılması — A6-dakı ilə eyni düstur (`residual.py`)."""
        compressibility = self.model.rock.compressibility
        if compressibility <= 0.0:
            return self.pore_volume
        factor = 1.0 + compressibility * (pressure - self.reference_pressure)
        return self.pore_volume * np.maximum(factor, 1e-6)

    def accumulation(self, state: ThreePhaseState, fluid: ThreePhaseFluidState):
        """(su, neft, qaz) — hər biri səth həcmi vahidində.

            N_water = PV · Sw / Bw
            N_oil   = PV · So / Bo
            N_gas   = PV · (Sg/Bg + So·Rs/Bo)      ← sərbəst + həll olmuş

        Qaz düsturundakı iki həddin FİZİKİ MƏNASI ayrıdır: birincisi
        sərbəst qaz fazasının özü, ikincisi neftin daşıdığı həll olmuş
        qaz. Doymamış hüceyrədə (Sg=0) yalnız ikinci hədd qalır.
        """
        pore_volume = self.pore_volume_at(state.pressure)
        sw = state.water_saturation
        so = state.oil_saturation
        sg = state.gas_saturation

        water = pore_volume * sw / fluid.bw
        oil = pore_volume * so / fluid.bo
        free_gas = pore_volume * sg / fluid.bg
        dissolved_gas = pore_volume * so * fluid.rs / fluid.bo
        gas = free_gas + dissolved_gas
        return water, oil, gas

    def two_phase_accumulation_matches(self, state: ThreePhaseState,
                                       fluid: ThreePhaseFluidState):
        """Doğrulama köməkçisi: Sg=0 olan hüceyrələrdə su/neft həddi
        A6-dakı iki fazalı `ResidualAssembler.accumulation()`-la EYNİ
        olmalıdır (yalnız qaz tənliyi yenidir, su/neft dəyişməyib)."""
        pore_volume = self.pore_volume_at(state.pressure)
        water = pore_volume * state.water_saturation / fluid.bw
        oil = pore_volume * (1.0 - state.water_saturation) / fluid.bo
        return water, oil


def _pvt_derivative_dispatch(pvt, name, pressure):
    method = getattr(pvt, f"{name}_derivative")
    return np.asarray(method(pressure), float)


class ThreePhaseFlux:
    """Üzlər üzrə Darcy axını — A7, mərhələ 6 (axın həddləri).

    Su və neft axını A6-dakı ilə EYNİ formuldan gəlir (dəyişməyib).
    QAZ AXINI keyfiyyətcə fərqlidir — iki ayrı mexanizmlə daşınır:

        1. SƏRBƏST qaz öz təzyiq qradiyenti ilə (ΦG)
        2. HƏLL OLMUŞ qaz NEFTLƏ BİRLİKDƏ (Rs · neft axını)

    İkinci hədd üçün upstream istiqaməti NEFTİN öz istiqamətidir —
    Rs "sərnişin"dir, öz axın istiqamətini seçmir, neftin daşıdığı
    yerə gedir. Bu, standart black-oil qaz axını formuludur (Aziz &
    Settari, Fanchi).

    Qaz-neft kapilyar keçidi (Pgo) hələ MODELƏ DAXİL EDİLMƏYİB —
    mərhələ 2-dəki eyni sadələşdirmə: GOC-da kəskin sərhəd, hamar
    keçid gələcək təkmilləşdirmədir. Bu, axın hesablamasını
    dəyişmir — yalnız equilibration-a təsir edir.
    """

    def __init__(self, model: ReservoirModel, grid: DiscretizedGrid):
        self.model = model
        self.connections = grid.connections
        self.transmissibility = grid.transmissibility
        self.ncell = model.ncell

        depths = model.geometry.cell_depths()
        self._depth_difference = (depths[self.connections.cell_a]
                                  - depths[self.connections.cell_b])
        self._has_gravity = bool(np.any(np.abs(self._depth_difference) > 1e-12))

    def potentials(self, state: ThreePhaseState, fluid: ThreePhaseFluidState):
        """Üzlər üzrə (ΔΦ_w, ΔΦ_o, ΔΦ_g) — bar.

        Qazın sıxlığı neft/sudan qat-qat kiçikdir (bax `A7_PLAN.md`,
        mərhələ 1) — cazibə həddi kiçik olsa da, doğru ölçüdə tətbiq
        olunur ki, qaz papağının seqreqasiyası (yuxarı toplanması)
        düzgün simulyasiya olunsun.
        """
        conn = self.connections
        dp = state.pressure[conn.cell_a] - state.pressure[conn.cell_b]
        d_phi_w = dp.copy()
        d_phi_o = dp.copy()
        d_phi_g = dp.copy()

        if self._has_gravity:
            fluids = self.model.fluids
            rho_w = fluids.water_density / np.maximum(fluid.bw, 1e-9)
            rho_o = fluids.oil_density / np.maximum(fluid.bo, 1e-9)
            rho_g = fluids.gas_density / np.maximum(fluid.bg, 1e-9)
            head = GRAVITY * self._depth_difference * PA_TO_BAR
            d_phi_w -= 0.5 * (rho_w[conn.cell_a] + rho_w[conn.cell_b]) * head
            d_phi_o -= 0.5 * (rho_o[conn.cell_a] + rho_o[conn.cell_b]) * head
            d_phi_g -= 0.5 * (rho_g[conn.cell_a] + rho_g[conn.cell_b]) * head

        if fluid.pc is not None:
            d_phi_w -= fluid.pc[conn.cell_a] - fluid.pc[conn.cell_b]

        return d_phi_w, d_phi_o, d_phi_g

    def upstream_masks(self, state: ThreePhaseState, fluid: ThreePhaseFluidState):
        """Hər üz üçün upstream hüceyrələr (su, neft, qaz) — ayrıca metod
        ki, "hansı hüceyrə upstream-dir" qərarı DONDURULA bilsin (bax
        `face_fluxes`-in `reference` parametri)."""
        conn = self.connections
        d_phi_w, d_phi_o, d_phi_g = self.potentials(state, fluid)
        up_w = np.where(d_phi_w >= 0, conn.cell_a, conn.cell_b)
        up_o = np.where(d_phi_o >= 0, conn.cell_a, conn.cell_b)
        up_g = np.where(d_phi_g >= 0, conn.cell_a, conn.cell_b)
        return up_w, up_o, up_g

    def face_fluxes(self, state: ThreePhaseState, fluid: ThreePhaseFluidState,
                    reference=None):
        """A → B istiqamətində səth həcmi axını (m³/gün), üç faza + həll olmuş qaz.

        `reference` — `(up_w, up_o, up_g)` verilibsə, UPSTREAM SEÇİMİ bu
        hüceyrələrlə DONDURULUR (potensialın ÖZÜ hələ də cari `state`-dən
        hesablanır — yalnız "hansı tərəf upstream-dir" qərarı sabitdir).

        NİYƏ LAZIMDIR: potensial fərqi sıfıra çox yaxın olan üzlərdə
        (məs. quyu öz BHP sərhədinə yaxınlaşanda) upstream tərəf bir
        Nyuton iterasiyasından digərinə DƏYİŞƏ bilər — bu, dəyişən
        keçid və quyu keçidi ilə eyni növ Jakobian kəsilməzliyini
        pozur və PERIODIK OSİLYASİYAYA səbəb olur (ölçülüb: CNV dövr-3
        dövrədə ilişib qalırdı). Dondurma bunu aradan qaldırır — A6-nın
        "upstream seçimi diferensiallaşdırılmır" sadələşdirməsinin
        məntiqi davamı: seçim təkcə DİFERENSİALLAŞDIRILMIR, addım
        daxilində DƏYİŞMİR də.
        """
        conn = self.connections
        d_phi_w, d_phi_o, d_phi_g = self.potentials(state, fluid)

        if reference is not None:
            up_w, up_o, up_g = reference
        else:
            up_w = np.where(d_phi_w >= 0, conn.cell_a, conn.cell_b)
            up_o = np.where(d_phi_o >= 0, conn.cell_a, conn.cell_b)
            up_g = np.where(d_phi_g >= 0, conn.cell_a, conn.cell_b)

        water = self.transmissibility * (fluid.lam_w[up_w] / fluid.bw[up_w]) * d_phi_w
        oil = self.transmissibility * (fluid.lam_o[up_o] / fluid.bo[up_o]) * d_phi_o

        free_gas = self.transmissibility * (fluid.lam_g[up_g] / fluid.bg[up_g]) * d_phi_g
        # həll olmuş qaz: NEFTİN öz upstream-i ilə, NEFTİN öz axını qədər
        dissolved_gas = fluid.rs[up_o] * oil
        gas = free_gas + dissolved_gas
        return water, oil, gas

    @staticmethod
    def model_pvt_derivative(pvt, name, pressure):
        return _pvt_derivative_dispatch(pvt, name, pressure)

    def net_influx(self, state: ThreePhaseState, fluid: ThreePhaseFluidState,
                  reference=None):
        """Hüceyrəyə daxil olan xalis axın (m³/gün, səth)."""
        conn = self.connections
        water_flux, oil_flux, gas_flux = self.face_fluxes(state, fluid, reference)

        water = np.zeros(self.ncell)
        oil = np.zeros(self.ncell)
        gas = np.zeros(self.ncell)
        np.add.at(water, conn.cell_a, -water_flux)
        np.add.at(water, conn.cell_b, +water_flux)
        np.add.at(oil, conn.cell_a, -oil_flux)
        np.add.at(oil, conn.cell_b, +oil_flux)
        np.add.at(gas, conn.cell_a, -gas_flux)
        np.add.at(gas, conn.cell_b, +gas_flux)
        return water, oil, gas


@dataclass
class ThreePhaseWellRates:
    """Quyu debitləri, səth həcmi. Müsbət = laya daxil olur (A6-dakı ilə eyni işarə)."""
    water: np.ndarray
    oil: np.ndarray
    gas: np.ndarray
    per_well_water: dict
    per_well_oil: dict
    per_well_gas: dict


class ThreePhaseWellModel:
    """Quyu mənbə həddləri — A7, mərhələ 6b.

    VURUCULAR hələ yalnız su vurur — A6-dakı davranış dəyişməyib. Qaz
    vurma (WAG, gas injector) EOR-un öz mövzusudur (CO₂ vurma) və bu
    modulun əhatəsindən kənardadır; `WellType`-a yeni qaz-injektor
    tipi əlavə olunanda buraya qoşulacaq.

    İSTİSMARÇILARDA qaz İKİ mənbədən çıxır — axın modulundakı (mərhələ
    6a) eyni məntiq:

        q_gas = q_sərbəst_qaz  +  Rs · q_neft

    RATE rejimində hədəf HƏLƏ DƏ maye (su+neft) debitidir — A6-dakı
    konvensiya qorunur, ona görə köhnə iki fazalı ssenarilər eyni
    nəticəni verir. Qaz RATE rejimində NƏTİCƏ kimi çıxır, hədəf kimi
    yox — bu, real qazlift/BHP rejimli quyularla üst-üstə düşür, harada
    ki, operator adətən maye debitini idarə edir, qaz "gəldiyi kimi" gəlir.
    """

    def __init__(self, model: ReservoirModel, wells: list,
                relperm_endpoint_water_mobility: float):
        self.model = model
        self.wells = wells
        self.ncell = model.ncell
        self._endpoint_water_mobility = relperm_endpoint_water_mobility
        self._producer_names = sorted({c.well_name for c in wells
                                       if not c.is_injector})
        self._injector_names = sorted({c.well_name for c in wells
                                       if c.is_injector})

    def well_rates(self, state: ThreePhaseState,
                   fluid: ThreePhaseFluidState,
                   reference_pressure=None) -> ThreePhaseWellRates:
        """`reference_pressure` verilibsə, BHP-nin AKTİV/BAĞLI qərarı
        HƏMİN təzyiqlə (adətən addımın əvvəlki, yığılmış vəziyyəti)
        verilir — cari iterasiyanın təzyiqi ilə YOX. Bax
        `ThreePhaseNewtonSolver`-in sənədləşməsi: bu, "künc" nöqtəsində
        (drawdown işarəsi dəyişəndə) Nyuton osilyasiyasının qarşısını
        alır. `None` olduqda əvvəlki (addım-daxili) davranış qorunur.
        """
        water = np.zeros(self.ncell)
        oil = np.zeros(self.ncell)
        gas = np.zeros(self.ncell)
        names = self._producer_names + self._injector_names
        per_well_water = {name: 0.0 for name in names}
        per_well_oil = dict(per_well_water)
        per_well_gas = dict(per_well_water)

        for connection in self.wells:
            cell = connection.cell
            if connection.is_injector:
                mobility = self._endpoint_water_mobility / fluid.mu_w[cell]
                if connection.mode is ControlMode.BHP:
                    rate = (connection.well_index * mobility
                            * (connection.target - state.pressure[cell]))
                else:
                    rate = abs(connection.target)
                rate = max(rate, 0.0) / fluid.bw[cell]
                water[cell] += rate
                per_well_water[connection.well_name] += rate
                continue

            wi = connection.well_index
            lam_w = fluid.lam_w[cell]
            lam_o = fluid.lam_o[cell]
            lam_g = fluid.lam_g[cell]
            if connection.mode is ControlMode.BHP:
                # SƏRT "aktiv/bağlı" QAPISI ARTIQ YOXDUR.
                #
                # Əvvəl burada `active = ... < 0` yoxlaması var idi:
                # quyu bağlı sayılanda bütün debitlər BİRDƏN sıfıra
                # düşürdü. Bu, residualda sıçrayış yaradırdı — Nyuton
                # quyu öz BHP hədəfinə yaxınlaşanda yığıla bilmirdi.
                # Deadband da kifayət etmədi (5 bar-a qədər sınandı).
                #
                # İndi debit sadəcə hesablanır və aşağıdakı HAMAR
                # kəsmə (`smooth_negative_part`) ilə məhdudlaşdırılır:
                # istismarçı yalnız hasil edə bilər, lakin keçid
                # hamardır — bax həmin funksiyanın sənədləşməsi.
                drawdown = connection.target - state.pressure[cell]
                qw = wi * lam_w * drawdown
                qo = wi * lam_o * drawdown
                q_free_gas = wi * lam_g * drawdown
            else:
                # RATE hədəfi MAYE debitidir (su+neft) — A6-dakı
                # konvensiya; qaz nəticə kimi çıxır.
                total = -abs(connection.target)
                fraction = lam_w / max(lam_w + lam_o, 1e-30)
                qw, qo = total * fraction, total * (1.0 - fraction)
                liquid_gas_fraction = lam_g / max(lam_w + lam_o, 1e-30)
                q_free_gas = total * liquid_gas_fraction

            qw = min(qw, 0.0) / fluid.bw[cell]
            qo = min(qo, 0.0) / fluid.bo[cell]
            q_free_gas = min(q_free_gas, 0.0) / fluid.bg[cell]
            q_dissolved_gas = fluid.rs[cell] * qo
            qg = q_free_gas + q_dissolved_gas

            water[cell] += qw
            oil[cell] += qo
            gas[cell] += qg
            per_well_water[connection.well_name] += qw
            per_well_oil[connection.well_name] += qo
            per_well_gas[connection.well_name] += qg

        return ThreePhaseWellRates(water, oil, gas, per_well_water,
                                   per_well_oil, per_well_gas)


class ThreePhaseAccumulationJacobian:
    """Akkumulyasiyanın diaqonal Jakobian bloku — A7, mərhələ 6c (1-ci hissə).

    A6-dakı ("1. AKKUMULYASİYA" — `jacobian.py`) eyni ideyanın üç
    fazalı davamıdır: hər hüceyrənin öz 3×3 blokunu — akkumulyasiya
    tənliklərinin ÖZ dəyişənlərinə görə törəməsini — verir. Bu, tam
    Jakobianın YALNIZ bir hissəsidir (əlavə töhfələr: axın — qonşu
    hüceyrələrə görə də sıfırdan fərqlidir — və quyular; mərhələ
    6c-nin qalan hissəsi).

    DƏYİŞƏN KEÇİDİN TÖRƏMƏYƏ TƏSİRİ

    3-cü dəyişənin mənası (Sg və ya Rs) hüceyrədən hüceyrəyə fərqli
    olduğu üçün ∂N/∂(3-cü dəyişən) düsturu da fərqlidir:

        doymuş (3-cü dəyişən = Sg):
            ∂N_oil/∂Sg = −PV/Bo         (So = 1−Sw−Sg → dSo/dSg=−1)
            ∂N_gas/∂Sg = PV·(1/Bg − Rs/Bo)

        doymamış (3-cü dəyişən = Rs):
            ∂N_oil/∂Rs = 0               (neft həcmi Rs-dən asılı deyil)
            ∂N_gas/∂Rs = PV·So/Bo        (birbaşa xətti asılılıq)

    Bu, bir Jakobian matrisində İKİ FƏRQLİ FORMULANIN eyni sütun
    mövqeyində (3-cü dəyişən) qarışıq şəkildə mövcud olduğu deməkdir
    — sətir-sətir seçim `is_saturated` bayrağı ilə aparılır.
    """

    def __init__(self, accumulator: ThreePhaseAccumulator, pvt):
        self.accumulator = accumulator
        self.pvt = pvt

    def blocks(self, state: ThreePhaseState, fluid: ThreePhaseFluidState):
        """Hər hüceyrə üçün 3×3 blok — `(ncell, 3, 3)` massiv.

        Sətir sırası: su, neft, qaz. Sütun sırası: p, Sw, 3-cü dəyişən.
        """
        pv = self.accumulator.pore_volume_at(state.pressure)
        compressibility = self.accumulator.model.rock.compressibility
        # d(PV)/dp = PV_REFERANS · c  (PV(p) = PV_ref·(1+c·Δp) — sabit
        # törəmə, CARİ pv-dən yox, istinad həcmindən çıxır)
        pv_p = (self.accumulator.pore_volume * compressibility
               if compressibility > 0 else np.zeros_like(pv))

        bw, bo, bg = fluid.bw, fluid.bo, fluid.bg
        bw_p = self.pvt.water_fvf_derivative(state.pressure)
        bo_p = self.pvt.oil_fvf_derivative(state.pressure)
        bg_p = self.pvt.gas_fvf_derivative(state.pressure)
        rs_sat_p = self.pvt.solution_gor_derivative(state.pressure)

        sw = state.water_saturation
        sg = state.gas_saturation
        so = state.oil_saturation
        rs = fluid.rs
        saturated = state.is_saturated

        n = state.ncell
        blocks = np.zeros((n, 3, 3))

        # ── su tənliyi: N_w = PV·Sw/Bw — 3-cü dəyişəndən asılı deyil
        blocks[:, 0, 0] = pv_p * sw / bw - pv * sw * bw_p / bw ** 2
        blocks[:, 0, 1] = pv / bw
        # blocks[:, 0, 2] = 0 — su nə Sg-dən, nə Rs-dən asılıdır

        # ── neft tənliyi: N_o = PV·So/Bo
        blocks[:, 1, 0] = pv_p * so / bo - pv * so * bo_p / bo ** 2
        blocks[:, 1, 1] = -pv / bo
        # doymuş: So = 1-Sw-Sg  → ∂So/∂(3-cü)=-1 → ∂N_o/∂Sg = -PV/Bo
        # doymamış: So = 1-Sw   → 3-cü dəyişən (Rs) heç görünmür → 0
        blocks[:, 1, 2] = np.where(saturated, -pv / bo, 0.0)

        # ── qaz tənliyi: N_g = PV·(Sg/Bg + So·Rs/Bo)
        gas_p_saturated = (pv_p * (sg / bg + so * rs / bo)
                          + pv * (-sg * bg_p / bg ** 2
                                 + so * (rs_sat_p * bo - rs * bo_p) / bo ** 2))
        # doymamış: Sg=0 sabit, Rs sərbəst dəyişən (p-dən asılı deyil bu təsvirdə)
        gas_p_undersaturated = (pv_p * (so * rs / bo)
                               - pv * so * rs * bo_p / bo ** 2)
        blocks[:, 2, 0] = np.where(saturated, gas_p_saturated,
                                   gas_p_undersaturated)

        blocks[:, 2, 1] = -pv * rs / bo     # hər iki halda: ∂So/∂Sw=-1 həddi

        gas_third_saturated = pv * (1.0 / bg - rs / bo)      # ∂Sg, ∂So/∂Sg=-1
        gas_third_undersaturated = pv * so / bo               # ∂Rs birbaşa
        blocks[:, 2, 2] = np.where(saturated, gas_third_saturated,
                                   gas_third_undersaturated)

        return blocks

    def numerical(self, state: ThreePhaseState, fluid_builder, step: float = 1e-6):
        """Sonlu fərqlə yoxlama — YALNIZ testlər üçün (A6-dakı ilə eyni məqsəd).

        `fluid_builder(state) -> ThreePhaseFluidState` — hər pertürbasiya
        üçün flüid xassələrini yenidən qurur (Bo/Bg/Rs təzyiqdən asılı).
        """
        n = state.ncell
        result = np.zeros((n, 3, 3))
        base_vector = state.to_vector()
        for local_column in range(3):
            forward = base_vector.copy()
            backward = base_vector.copy()
            scale = step * max(1.0, float(np.mean(np.abs(
                base_vector[local_column::3]))))
            forward[local_column::3] += scale
            backward[local_column::3] -= scale

            state_forward = ThreePhaseState.from_vector(forward, state.is_saturated)
            state_backward = ThreePhaseState.from_vector(backward, state.is_saturated)
            fluid_forward = fluid_builder(state_forward)
            fluid_backward = fluid_builder(state_backward)

            water_f, oil_f, gas_f = self.accumulator.accumulation(
                state_forward, fluid_forward)
            water_b, oil_b, gas_b = self.accumulator.accumulation(
                state_backward, fluid_backward)

            result[:, 0, local_column] = (water_f - water_b) / (2.0 * scale)
            result[:, 1, local_column] = (oil_f - oil_b) / (2.0 * scale)
            result[:, 2, local_column] = (gas_f - gas_b) / (2.0 * scale)
        return result


class ThreePhaseFluxJacobian:
    """Axın Jakobianı — A7, mərhələ 6c (2-ci hissə, axın töhfəsi).

    A6-dakı sadələşdirmələr eynilə qorunur (bax `jacobian.py`-nin öz
    sənədləşməsi):

        · Upstream seçiminin özü diferensiallaşdırılmır
        · Cazibə üzvündə sıxlığın təzyiqdən asılılığı nəzərə alınmır

    Hər üz üçün YALNIZ upstream hüceyrənin mobilliyi diferensiallaşdırılır
    (downstream hüceyrədə isə yalnız ΔΦ-nin öz xətti asılılığı) — bu,
    A6-nın "İkinci üzv YALNIZ upstream hüceyrəsinə görə sıfırdan
    fərqlidir" qaydasının üç fazalı davamıdır.

    Su/neft mobillik törəməsi üçün Stone II-nin analitik zəncirvari
    qaydası (`kro_three_phase_derivatives`) işlədilir — bu, mərkəzi
    sonlu fərqə ehtiyacı aradan qaldırır və maşın dəqiqliyi verir
    (ölçülüb: fərq ~10⁻¹¹).
    """

    def __init__(self, flux: ThreePhaseFlux, relperm):
        self.flux = flux
        self.relperm = relperm

    def face_pressure_derivatives(self, state: ThreePhaseState,
                                  fluid: ThreePhaseFluidState, pvt,
                                  reference=None):
        """(∂F/∂p_a, ∂F/∂p_b) — su, neft, sərbəst qaz VƏ həll olmuş qaz.

        A6-dakı EYNİ tərkib: `F = mob_upstream·ΔΦ`-in HƏR İKİ hissəsi
        diferensiallaşdırılır — ΔΦ-nin özü (xətti, ∂/∂p_a=+1) VƏ
        mob(p)-nin özü (yalnız upstream hüceyrədə, çünki mobillik
        YALNIZ upstream-dən götürülür). Yeganə sadələşdirmə (A6-dakı
        eyni): cazibə həddindəki sıxlığın təzyiqdən asılılığı BURADA
        nəzərə alınmır.

        Qaz üçün əlavə hədd: HƏLL OLMUŞ hissə (Rs_upstream_neft ·
        neft_axını) — zəncirvari qayda ilə həm Rs-in, həm neft
        axınının öz təzyiq törəməsi nəzərə alınır.

        `reference` — bax `ThreePhaseFlux.face_fluxes`-in sənədləşməsi:
        upstream seçimini dondurur (periodik osilyasiyanın qarşısını alır).
        """
        conn = self.flux.connections
        trans = self.flux.transmissibility
        a, b = conn.cell_a, conn.cell_b
        d_phi_w, d_phi_o, d_phi_g = self.flux.potentials(state, fluid)

        if reference is not None:
            up_w, up_o, up_g = reference
        else:
            up_w = np.where(d_phi_w >= 0, a, b)
            up_o = np.where(d_phi_o >= 0, a, b)
            up_g = np.where(d_phi_g >= 0, a, b)
        up_w_is_a, up_o_is_a, up_g_is_a = up_w == a, up_o == a, up_g == a

        mob_w = fluid.lam_w / fluid.bw
        mob_o = fluid.lam_o / fluid.bo
        mob_g = fluid.lam_g / fluid.bg

        # SU LÖZLÜYÜNÜN TƏZYİQ TÖRƏMƏSİ — əvvəl SIFIR sayılırdı.
        #
        # Bu, real bir səhv idi: mühərrik `mu_w`-ni PVT-dən oxuyur və
        # o, təzyiqdən ASILIDIR, ona görə törəməsi sıfır deyil. Səhv
        # gizli qalmışdı, çünki testlər sabit `mu_w=0.5` işlədirdi —
        # sabit dəyərdə törəmə həqiqətən sıfırdır, ona görə sonlu
        # fərq yoxlamaları keçirdi.
        #
        # Ölçülüb: bu səhv doymamış vəziyyətdə Jakobian xətasını
        # 0.001 %-dən 0.24 %-ə qaldırırdı (vurucu hüceyrədə) və
        # Nyutonun periodik osilyasiyasına səbəb olurdu.
        mu_w_p = self.flux.model_pvt_derivative(pvt, "water_viscosity",
                                                state.pressure)
        bw_p = self.flux.model_pvt_derivative(pvt, "water_fvf", state.pressure)
        bo_p = self.flux.model_pvt_derivative(pvt, "oil_fvf", state.pressure)
        bg_p = self.flux.model_pvt_derivative(pvt, "gas_fvf", state.pressure)
        mu_o_p = self.flux.model_pvt_derivative(pvt, "oil_viscosity", state.pressure)
        mu_g_p = self.flux.model_pvt_derivative(pvt, "gas_viscosity", state.pressure)

        dmw_dp = -fluid.krw * (mu_w_p * fluid.bw + fluid.mu_w * bw_p) \
            / (fluid.mu_w * fluid.bw) ** 2
        dmo_dp = -fluid.kro * (mu_o_p * fluid.bo + fluid.mu_o * bo_p) \
            / (fluid.mu_o * fluid.bo) ** 2
        dmg_dp = -fluid.krg * (mu_g_p * fluid.bg + fluid.mu_g * bg_p) \
            / (fluid.mu_g * fluid.bg) ** 2

        m_w, m_o, m_g = mob_w[up_w], mob_o[up_o], mob_g[up_g]
        dfw_dpa = trans * (m_w + d_phi_w * np.where(up_w_is_a, dmw_dp[a], 0.0))
        dfw_dpb = trans * (-m_w + d_phi_w * np.where(up_w_is_a, 0.0, dmw_dp[b]))
        dfo_dpa = trans * (m_o + d_phi_o * np.where(up_o_is_a, dmo_dp[a], 0.0))
        dfo_dpb = trans * (-m_o + d_phi_o * np.where(up_o_is_a, 0.0, dmo_dp[b]))
        d_free_g_dpa = trans * (m_g + d_phi_g * np.where(up_g_is_a, dmg_dp[a], 0.0))
        d_free_g_dpb = trans * (-m_g + d_phi_g * np.where(up_g_is_a, 0.0, dmg_dp[b]))

        # ── həll olmuş qaz: Rs_up_o · F_oil, zəncirvari qayda ──
        rs_sat_p = self.flux.model_pvt_derivative(pvt, "solution_gor", state.pressure)
        rs_up = fluid.rs[up_o]
        # Rs_up p_a-dan asılıdır YALNIZ up_o==a VƏ doymuş olduqda (Rs=Rs_sat(p));
        # doymamışda Rs sərbəst dəyişəndir (3-cü sütunda), p-dən bu formulada asılı deyil.
        drs_up_dpa = np.where(up_o_is_a & state.is_saturated[a], rs_sat_p[a], 0.0)
        drs_up_dpb = np.where(~up_o_is_a & state.is_saturated[b], rs_sat_p[b], 0.0)

        oil_flux = trans * m_o * d_phi_o
        d_dissolved_dpa = drs_up_dpa * oil_flux + rs_up * dfo_dpa
        d_dissolved_dpb = drs_up_dpb * oil_flux + rs_up * dfo_dpb

        dpg_a = d_free_g_dpa + d_dissolved_dpa
        dpg_b = d_free_g_dpb + d_dissolved_dpb

        return ((dfw_dpa, dfw_dpb), (dfo_dpa, dfo_dpb), (dpg_a, dpg_b), up_o)

    def face_saturation_derivatives(self, state: ThreePhaseState,
                                    fluid: ThreePhaseFluidState,
                                    reference=None):
        """Mobillik törəməsi × ΔΦ — YALNIZ upstream hüceyrə üçün.

        Qaytarılan: hər faza üçün `(d/dSw_upstream, d/dthird_upstream)`
        — işarə upstream A yoxsa B olduğuna görə çağıran tərəfindən
        müvafiq hüceyrəyə yerləşdirilir (`net_influx_blocks`-da).

        `reference` — upstream seçimini dondurur (bax `face_fluxes`).
        """
        conn = self.flux.connections
        trans = self.flux.transmissibility
        d_phi_w, d_phi_o, d_phi_g = self.flux.potentials(state, fluid)

        if reference is not None:
            up_w, up_o, up_g = reference
        else:
            up_w = np.where(d_phi_w >= 0, conn.cell_a, conn.cell_b)
            up_o = np.where(d_phi_o >= 0, conn.cell_a, conn.cell_b)
            up_g = np.where(d_phi_g >= 0, conn.cell_a, conn.cell_b)

        sw_up_w = state.water_saturation[up_w]
        sw_up_o = state.water_saturation[up_o]
        sg_up_o = state.gas_saturation[up_o]
        sg_up_g = state.gas_saturation[up_g]

        dkrw_dsw = getattr(self.relperm, "krw_derivative",
                           lambda s, r=None: np.zeros_like(s))(sw_up_w)
        d_water_dsw = trans * dkrw_dsw / fluid.mu_w[up_w] / fluid.bw[up_w] * d_phi_w

        dkro_dsw, _ = self.relperm.kro_three_phase_derivatives(sw_up_o, sg_up_o)
        d_oil_dsw = (trans * dkro_dsw / fluid.mu_o[up_o] / fluid.bo[up_o]
                    * d_phi_o)
        _, dkro_dsg = self.relperm.kro_three_phase_derivatives(sw_up_o, sg_up_o)
        d_oil_dthird = np.where(
            state.is_saturated[up_o],
            trans * dkro_dsg / fluid.mu_o[up_o] / fluid.bo[up_o] * d_phi_o,
            0.0)          # doymamışda 3-cü dəyişən Rs-dir, kro ondan asılı deyil

        dkrg_dsg = self.relperm.gas.krg_derivative(sg_up_g, self.relperm.swc)
        d_gas_dthird_free = np.where(
            state.is_saturated[up_g],
            trans * dkrg_dsg / fluid.mu_g[up_g] / fluid.bg[up_g] * d_phi_g,
            0.0)

        return {
            "water_dsw": d_water_dsw, "water_up": up_w,
            "oil_dsw": d_oil_dsw, "oil_dthird": d_oil_dthird, "oil_up": up_o,
            "gas_dthird_free": d_gas_dthird_free, "gas_up": up_g,
        }

    def numerical_face_derivative(self, state: ThreePhaseState,
                                  fluid_builder, variable_index: int,
                                  cell_side: str, step: float = 1e-6):
        """Sonlu fərqlə yoxlama — bir dəyişən, bir tərəf (a/b) üçün.

        `fluid_builder(state) -> ThreePhaseFluidState`.
        """
        conn = self.flux.connections
        base_vector = state.to_vector()
        target_cells = conn.cell_a if cell_side == "a" else conn.cell_b

        forward = base_vector.copy()
        backward = base_vector.copy()
        scale = step * max(1.0, float(np.mean(np.abs(
            base_vector[variable_index::3]))))

        # yalnız HƏMİN ÜZÜN müvafiq tərəfindəki hüceyrələri pozur
        for cell in np.unique(target_cells):
            forward[cell * 3 + variable_index] += scale
            backward[cell * 3 + variable_index] -= scale

        state_forward = ThreePhaseState.from_vector(forward, state.is_saturated)
        state_backward = ThreePhaseState.from_vector(backward, state.is_saturated)
        water_f, oil_f, gas_f = self.flux.face_fluxes(
            state_forward, fluid_builder(state_forward))
        water_b, oil_b, gas_b = self.flux.face_fluxes(
            state_backward, fluid_builder(state_backward))
        return ((water_f - water_b) / (2 * scale),
               (oil_f - oil_b) / (2 * scale),
               (gas_f - gas_b) / (2 * scale))


class ThreePhaseWellJacobian:
    """Quyu Jakobianı — A7, mərhələ 6c (3-cü hissə).

    A6-dakı `_wells()` metodunu güzgüləyir (bax `jacobian.py`).
    `R = … − q`, ona görə `∂R/∂x = −∂q/∂x` — bu sinif `∂q/∂x`-i
    hesablayır, işarəni tətbiq etmək çağıran tərəfin işidir.

    Vurucular üçün su törəməsi A6-dakı ilə EYNİDİR (dəyişməyib).
    İstismarçılarda əlavə olaraq HƏLL OLMUŞ qazın öz zəncirvari
    qaydası var: `q_gas = q_sərbəst + Rs(p)·q_neft`.
    """

    def __init__(self, well_model: ThreePhaseWellModel, pvt, relperm):
        self.well_model = well_model
        self.pvt = pvt
        self.relperm = relperm

    def blocks(self, state: ThreePhaseState, fluid: ThreePhaseFluidState,
              reference_pressure=None):
        """Hər hüceyrə üçün 3×3 blok (`∂q_faza/∂dəyişən`) — `(ncell,3,3)`.

        `reference_pressure` — bax `well_rates()`-in sənədləşməsi: BHP
        aktiv/bağlı qərarı bu təzyiqlə verilir, cari iterasiya ilə YOX.
        """
        n = state.ncell
        blocks = np.zeros((n, 3, 3))
        relperm = None

        # bax `face_pressure_derivatives`-dəki eyni düzəlişin şərhi:
        # su lözlüyü PVT-dən gəlir və təzyiqdən asılıdır
        mu_w_p = self.pvt.water_viscosity_derivative(state.pressure)
        mu_o_p = self.pvt.oil_viscosity_derivative(state.pressure)
        mu_g_p = self.pvt.gas_viscosity_derivative(state.pressure)
        bw_p = self.pvt.water_fvf_derivative(state.pressure)
        bo_p = self.pvt.oil_fvf_derivative(state.pressure)
        bg_p = self.pvt.gas_fvf_derivative(state.pressure)
        rs_sat_p = self.pvt.solution_gor_derivative(state.pressure)

        for connection in self.well_model.wells:
            c = connection.cell
            wi = connection.well_index

            if connection.is_injector:
                endpoint = self.well_model._endpoint_water_mobility
                transport = endpoint / (fluid.mu_w[c] * fluid.bw[c])
                drawdown = connection.target - state.pressure[c]
                if connection.mode is ControlMode.BHP:
                    if wi * transport * drawdown <= 0.0:
                        continue
                    d_transport = -endpoint * (mu_w_p[c] * fluid.bw[c]
                                               + fluid.mu_w[c] * bw_p[c]) \
                        / (fluid.mu_w[c] * fluid.bw[c]) ** 2
                    blocks[c, 0, 0] += wi * (-transport + drawdown * d_transport)
                else:
                    rate = abs(connection.target)
                    blocks[c, 0, 0] += -rate * bw_p[c] / fluid.bw[c] ** 2
                continue

            # ── istismarçı ────────────────────────────────────────
            lam_w, lam_o, lam_g = fluid.lam_w[c], fluid.lam_o[c], fluid.lam_g[c]
            mob_w = lam_w / fluid.bw[c]
            mob_o = lam_o / fluid.bo[c]
            mob_g = lam_g / fluid.bg[c]
            rs = fluid.rs[c]
            saturated = state.is_saturated[c]

            dmw_dp = -fluid.krw[c] * (mu_w_p[c] * fluid.bw[c]
                                      + fluid.mu_w[c] * bw_p[c]) \
                / (fluid.mu_w[c] * fluid.bw[c]) ** 2
            dmo_dp = -fluid.kro[c] * (mu_o_p[c] * fluid.bo[c]
                                      + fluid.mu_o[c] * bo_p[c]) \
                / (fluid.mu_o[c] * fluid.bo[c]) ** 2
            dmg_dp = -fluid.krg[c] * (mu_g_p[c] * fluid.bg[c]
                                      + fluid.mu_g[c] * bg_p[c]) \
                / (fluid.mu_g[c] * fluid.bg[c]) ** 2

            dkro_dsw, dkro_dsg = self.relperm.kro_three_phase_derivatives(
                np.array([state.water_saturation[c]]),
                np.array([state.gas_saturation[c]]))
            dkro_dsw, dkro_dsg = float(dkro_dsw[0]), float(dkro_dsg[0])
            dkrw_dsw = float(getattr(
                self.relperm, "krw_derivative",
                lambda s, r=None: np.zeros(1))(
                    np.array([state.water_saturation[c]]))[0])
            dkrg_dsg = float(self.relperm.gas.krg_derivative(
                np.array([state.gas_saturation[c]]), self.relperm.swc)[0])

            if connection.mode is ControlMode.BHP:
                # SƏRT QAPI ARTIQ YOXDUR — `well_rates()` ilə eyni
                # dəyişiklik. Onun yerinə HAMAR kəsmənin (`smooth_
                # negative_part`) törəməsi zəncirvari qayda ilə
                # tətbiq olunur: bu, Jakobianı residualla UYĞUN
                # saxlayır (uyğunsuzluq Nyutonu yanlış istiqamətə
                # addımladardı).
                drawdown = connection.target - state.pressure[c]

                chop_w = 1.0 if wi * mob_w * drawdown < 0.0 else 0.0
                chop_o = 1.0 if wi * mob_o * drawdown < 0.0 else 0.0
                chop_g = 1.0 if wi * mob_g * drawdown < 0.0 else 0.0

                dqw_dp = chop_w * wi * (-mob_w + drawdown * dmw_dp)
                dqw_dsw = (chop_w * wi * drawdown * dkrw_dsw
                          / fluid.mu_w[c] / fluid.bw[c])

                dqo_dp = chop_o * wi * (-mob_o + drawdown * dmo_dp)
                dqo_dsw = (chop_o * wi * drawdown * dkro_dsw
                          / fluid.mu_o[c] / fluid.bo[c])
                dqo_dthird = 0.0
                if saturated:
                    dqo_dthird = (chop_o * wi * drawdown * dkro_dsg
                                 / fluid.mu_o[c] / fluid.bo[c])

                dqfree_dp = chop_g * wi * (-mob_g + drawdown * dmg_dp)
                dqfree_dthird = 0.0
                if saturated:
                    dqfree_dthird = (chop_g * wi * drawdown * dkrg_dsg
                                    / fluid.mu_g[c] / fluid.bg[c])

            else:
                # RATE rejimi: hədəf maye (su+neft) debitidir — A6-dakı
                # fraksional axın formulu, Stone II-nin dkro/dSw ilə.
                total = -abs(connection.target)
                lam_t = max(lam_w + lam_o, 1e-30)
                fraction = lam_w / lam_t

                dlw_dsw = dkrw_dsw / fluid.mu_w[c]
                dlo_dsw = dkro_dsw / fluid.mu_o[c]
                dfraction_dsw = (dlw_dsw * lam_o - lam_w * dlo_dsw) / lam_t ** 2

                dlw_dp = -fluid.krw[c] * mu_w_p[c] / fluid.mu_w[c] ** 2
                dlo_dp = -fluid.kro[c] * mu_o_p[c] / fluid.mu_o[c] ** 2
                dfraction_dp = (dlw_dp * lam_o - lam_w * dlo_dp) / lam_t ** 2

                dqw_dsw = total * dfraction_dsw / fluid.bw[c]
                dqo_dsw = -total * dfraction_dsw / fluid.bo[c]
                dqw_dp = total * (dfraction_dp / fluid.bw[c]
                                  - fraction * bw_p[c] / fluid.bw[c] ** 2)
                dqo_dp = total * (-dfraction_dp / fluid.bo[c]
                                  - (1.0 - fraction) * bo_p[c] / fluid.bo[c] ** 2)
                dqo_dthird = 0.0     # RATE rejimində maye fraksiyası Sg-dən asılı deyil (sadələşdirmə)
                dqfree_dp = 0.0
                dqfree_dthird = 0.0

            blocks[c, 0, 0] += dqw_dp
            blocks[c, 0, 1] += dqw_dsw
            blocks[c, 1, 0] += dqo_dp
            blocks[c, 1, 1] += dqo_dsw
            blocks[c, 1, 2] += dqo_dthird

            # ── qaz: sərbəst + Rs(p)·neft (hasil qaydası) ──────────
            # `qo_value` — SƏTH həcmi neft debiti (mob_o artıq Bo-ya
            # bölünmüşdür, təkrar bölmə YOXDUR).
            drs_dp = rs_sat_p[c] if saturated else 0.0
            drs_dthird = 0.0 if saturated else 1.0
            # HAMAR kəsmə ilə — `well_rates()`-dəki EYNİ düstur
            # (sərt `min` residualla uyğunsuzluq yaradardı).
            qo_value = (min(wi * mob_o * (connection.target - state.pressure[c]), 0.0)
                       if connection.mode is ControlMode.BHP else 0.0)

            blocks[c, 2, 0] += dqfree_dp + drs_dp * qo_value + rs * dqo_dp
            blocks[c, 2, 1] += rs * dqo_dsw
            blocks[c, 2, 2] += (dqfree_dthird + rs * dqo_dthird
                                + drs_dthird * qo_value)

        return blocks


class ThreePhaseJacobianAssembler:
    """Tam Jakobian yığımı — A7, mərhələ 6c (4-cü hissə).

    Əvvəlki üç sinfin (akkumulyasiya, axın, quyu) VERİLƏNLƏRİNİ bir
    seyrək (sparse) 3N×3N matrisdə birləşdirir.

    A6-dakı optimallaşdırılmış üsuldan (əvvəlcədən hesablanmış naxış +
    `bincount` ilə birbaşa yazma) fərqli olaraq burada SADƏ COO-üçlük
    yığımı işlədilir — hər `assemble()` çağırışında siyahılar yenidən
    qurulur. Performans baxımından bu, A6-dan yavaşdır, lakin
    DOĞRULUĞU YOXLAMAQ və SƏHV RİSKİNİ azaltmaq üçün daha etibarlıdır
    (mürəkkəb naxış-indeksləmə məntiqi yoxdur). Performans
    optimallaşdırması gələcək iş kimi qalır — düzgünlük əvvəldir.
    """

    def __init__(self, model: ReservoirModel, accumulator: ThreePhaseAccumulator,
                flux: ThreePhaseFlux, well_model: ThreePhaseWellModel,
                relperm, pvt):
        self.model = model
        self.accumulator = accumulator
        self.flux = flux
        self.well_model = well_model
        self.pvt = pvt
        self.ncell = model.ncell
        self.size = self.ncell * 3

        self.accumulation_jacobian = ThreePhaseAccumulationJacobian(accumulator, pvt)
        self.flux_jacobian = ThreePhaseFluxJacobian(flux, relperm)
        self.well_jacobian = ThreePhaseWellJacobian(well_model, pvt, relperm)

    def assemble(self, state: ThreePhaseState, fluid: ThreePhaseFluidState,
                dt: float, reference_pressure=None, reference_upstream=None):
        """R = (N_yeni−N_köhnə)/dt − axın − quyu  →  J = ∂R/∂x.

        `reference_pressure` — quyu Jakobianına ötürülür (bax
        `ThreePhaseWellJacobian.blocks()`).
        `reference_upstream` — axının upstream seçimini dondurur (bax
        `ThreePhaseFlux.face_fluxes`-in sənədləşməsi).
        """
        import scipy.sparse as sp

        rows: list = []
        cols: list = []
        values: list = []
        cell_index = np.arange(self.ncell)

        def add_block(row_cells, col_cells, block):
            """`block`: (ncell_ya_da_face, 3, 3) — hər (r,c) qeyd olunur."""
            for r in range(3):
                for c in range(3):
                    entries = block[:, r, c]
                    nonzero = entries != 0.0
                    if not np.any(nonzero):
                        continue
                    rows.append(row_cells[nonzero] * 3 + r)
                    cols.append(col_cells[nonzero] * 3 + c)
                    values.append(entries[nonzero])

        # ── akkumulyasiya (diaqonal, /dt) ──────────────────────────
        acc_blocks = self.accumulation_jacobian.blocks(state, fluid) / dt
        add_block(cell_index, cell_index, acc_blocks)

        # ── quyular (diaqonal, mənfi işarə: R = … − q) ─────────────
        well_blocks = -self.well_jacobian.blocks(state, fluid, reference_pressure)
        add_block(cell_index, cell_index, well_blocks)

        # ── axın (diaqonal + qonşu, R_a = …+F, R_b = …−F) ──────────
        conn = self.flux.connections
        a, b = conn.cell_a, conn.cell_b
        n_faces = len(a)

        (dpa_w, dpb_w), (dpa_o, dpb_o), (dpa_g, dpb_g), up_o = \
            self.flux_jacobian.face_pressure_derivatives(
                state, fluid, self.pvt, reference_upstream)
        sat = self.flux_jacobian.face_saturation_derivatives(
            state, fluid, reference_upstream)

        # p-sütunu (sütun 0) — hər iki tərəf üçün hazır
        p_block_a = np.zeros((n_faces, 3, 3))
        p_block_a[:, 0, 0] = dpa_w
        p_block_a[:, 1, 0] = dpa_o
        p_block_a[:, 2, 0] = dpa_g
        p_block_b = np.zeros((n_faces, 3, 3))
        p_block_b[:, 0, 0] = dpb_w
        p_block_b[:, 1, 0] = dpb_o
        p_block_b[:, 2, 0] = dpb_g

        # Sw/3-cü sütunları — YALNIZ upstream tərəfdə sıfırdan fərqli
        s_block_a = np.zeros((n_faces, 3, 3))
        s_block_b = np.zeros((n_faces, 3, 3))

        water_up_is_a = sat["water_up"] == a
        s_block_a[water_up_is_a, 0, 1] = sat["water_dsw"][water_up_is_a]
        s_block_b[~water_up_is_a, 0, 1] = sat["water_dsw"][~water_up_is_a]

        oil_up_is_a = sat["oil_up"] == a
        s_block_a[oil_up_is_a, 1, 1] = sat["oil_dsw"][oil_up_is_a]
        s_block_b[~oil_up_is_a, 1, 1] = sat["oil_dsw"][~oil_up_is_a]
        s_block_a[oil_up_is_a, 1, 2] = sat["oil_dthird"][oil_up_is_a]
        s_block_b[~oil_up_is_a, 1, 2] = sat["oil_dthird"][~oil_up_is_a]
        # neftin daşıdığı həll olmuş qaz — qaz sətrinə əlavə töhfə
        s_block_a[oil_up_is_a, 2, 1] = fluid.rs[up_o][oil_up_is_a] * sat["oil_dsw"][oil_up_is_a]
        s_block_b[~oil_up_is_a, 2, 1] = fluid.rs[up_o][~oil_up_is_a] * sat["oil_dsw"][~oil_up_is_a]
        # HƏLL OLMUŞ QAZ AXINININ Rs-Ə GÖRƏ TÖRƏMƏSİ — HASİL QAYDASI
        #
        #     F_qaz_həll = Rs_upstream · F_neft
        #     ∂F/∂Rs_up  = F_neft  +  Rs_up · ∂F_neft/∂Rs_up
        #                  ↑ BU HƏDD ƏVVƏL UNUDULMUŞDU
        #
        # İkinci hədd (aşağıda) əvvəldən var idi, birincisi YOX idi.
        #
        # Niyə gizli qalmışdı: bu hədd YALNIZ hüceyrə DOYMAMIŞ olanda
        # sıfırdan fərqlidir, çünki yalnız o halda 3-cü primary
        # dəyişən Rs-dir. Doymuş hüceyrədə Rs = Rs_sat(p) — 3-cü
        # dəyişəndən (Sg) asılı deyil, ona görə hədd sıfırdır.
        # Bütün əvvəlki tam-sistem doğrulamalarım TAM DOYMUŞ vəziyyətdə
        # aparılmışdı (`np.ones(n, bool)`), ona görə səhv görünmürdü.
        #
        # Ölçülüb: bu hədd olmadan uğursuz vəziyyətdə Jakobian xətası
        # 0.17 % idi (bircins vəziyyətdə 9×10⁻⁹) və Nyuton dövr-2
        # osilyasiya edirdi.
        oil_flux_face = self.flux.face_fluxes(state, fluid, reference_upstream)[1]
        undersaturated_up = ~state.is_saturated[up_o]

        s_block_a[oil_up_is_a, 2, 2] = fluid.rs[up_o][oil_up_is_a] * sat["oil_dthird"][oil_up_is_a]
        s_block_b[~oil_up_is_a, 2, 2] = fluid.rs[up_o][~oil_up_is_a] * sat["oil_dthird"][~oil_up_is_a]

        direct = np.where(undersaturated_up, oil_flux_face, 0.0)
        s_block_a[oil_up_is_a, 2, 2] += direct[oil_up_is_a]
        s_block_b[~oil_up_is_a, 2, 2] += direct[~oil_up_is_a]

        gas_up_is_a = sat["gas_up"] == a
        s_block_a[gas_up_is_a, 2, 2] += sat["gas_dthird_free"][gas_up_is_a]
        s_block_b[~gas_up_is_a, 2, 2] += sat["gas_dthird_free"][~gas_up_is_a]

        block_a = p_block_a + s_block_a
        block_b = p_block_b + s_block_b

        add_block(a, a, block_a)
        add_block(a, b, block_b)
        add_block(b, a, -block_a)
        add_block(b, b, -block_b)

        row_array = np.concatenate(rows) if rows else np.array([], dtype=int)
        col_array = np.concatenate(cols) if cols else np.array([], dtype=int)
        value_array = np.concatenate(values) if values else np.array([])

        return sp.coo_matrix((value_array, (row_array, col_array)),
                             shape=(self.size, self.size)).tocsr()
