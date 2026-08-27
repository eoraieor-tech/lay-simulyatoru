"""Quyu tənlikləri — OPM tipli standart quyu modeli, MƏRHƏLƏ 2.

MƏRHƏLƏ 1-DƏ nə var idi: quyu BHP-si naməlum dəyişən kimi, vektor
yerləşməsi. BU MƏRHƏLƏDƏ: həmin naməlumdan debitlərin hesablanması və
quyunun öz tənliyi.

ƏSAS FƏRQ — KƏSMƏ YOXDUR

Köhnə modeldə (`ThreePhaseWellModel`) debit belə hesablanırdı:

    q = WI · λ · (p_hədəf − p_hüceyrə)
    q = min(q, 0)            ← SƏRT KƏSMƏ (istismarçı üçün)

`min` sıfır nöqtəsində SINIQdır: törəməsi 0-dan 1-ə sıçrayır. Quyu öz
hədəfinə yaxınlaşanda debit məhz sıfıra yaxınlaşır — Nyuton hər
iterasiyada bu sınıq nöqtənin bir tərəfindən digərinə atılır.
(Bax `A7_PLAN.md`: hamarlaşdırma da sınandı, işləmədi və fizikanı
pozdu.)

Yeni modeldə BHP naməlumdur, debit isə sadəcə:

    q_α = WI · λ_α · (p_bhp − p_hüceyrə)

İşarə TƏBİİ olaraq çıxır: `p_bhp < p_hüceyrə` → q < 0 (hasilat),
`p_bhp > p_hüceyrə` → q > 0 (vurma). Heç bir kəsmə yoxdur, funksiya
hər yerdə hamardır. Bu, eyni zamanda ÇARPAZ AXINI (cross-flow) təbii
modelləşdirir — OPM də məhz belə edir.

Quyu isə öz tənliyi ilə idarə olunur:

    BHP idarəsində:   R_ctrl = p_bhp − p_hədəf
    RATE idarəsində:  R_ctrl = Σ q_maye − q_hədəf

İdarəetmə rejimi dəyişəndə DEBİT sıçramır — sadəcə hansı tənliyin
işlədiyi dəyişir.

MOBİLLİYİN "UPSTREAM" SEÇİMİ

Axın istiqamətindən asılıdır (OPM-dəki ilə eyni):

    laya DAXİL olur (vurma)  → vurulan fazanın mobilliyi
                                (su vurucusunda yalnız su)
    laydan ÇIXIR (hasilat)   → hüceyrənin öz mobillikləri
                                (bütün fazalar öz nisbətində)
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from ...domain.wells import ControlMode
from .three_phase_residual import ThreePhaseFluidState
from .three_phase_state import ThreePhaseState
from .well_state import WellUnknowns


@dataclass
class PerforationRates:
    """Perforasiya debitləri, səth həcmi.

    İşarə konvensiyası köhnə modellə EYNİDİR: müsbət = laya daxil
    olur, mənfi = laydan çıxır.
    """
    water: np.ndarray
    oil: np.ndarray
    gas: np.ndarray
    per_well_water: Dict[str, float] = field(default_factory=dict)
    per_well_oil: Dict[str, float] = field(default_factory=dict)
    per_well_gas: Dict[str, float] = field(default_factory=dict)


class StandardWellModel:
    """OPM tipli quyu modeli — BHP naməlum dəyişəndir.

    Bu sinif KÖHNƏ `ThreePhaseWellModel`-i ƏVƏZ ETMİR — mərhələli
    keçid üçün onun yanında yaşayır. Mühərrikə qoşulma mərhələ 4-də
    olacaq, ondan əvvəl mövcud kod tam toxunulmaz qalır.
    """

    def __init__(self, connections: Sequence, ncell: int,
                 endpoint_water_mobility: float = 0.35):
        self.connections = list(connections)
        self.ncell = ncell
        self._endpoint_water_mobility = endpoint_water_mobility
        self.names: List[str] = []
        for connection in self.connections:
            if connection.well_name not in self.names:
                self.names.append(connection.well_name)
        self.shut: Dict[str, bool] = {name: False for name in self.names}
        """Bağlı quyular — bax `update_shut_wells()`."""

    def update_shut_wells(self, state: ThreePhaseState,
                          fluid: ThreePhaseFluidState,
                          wells: WellUnknowns) -> bool:
        """Hasilat edə bilməyən quyunu BAĞLAYIR.

        FİZİKİ SƏBƏB: BHP idarəli istismarçının hüceyrə təzyiqi öz
        hədəfindən AŞAĞI düşəndə düstur `q > 0` verir — yəni quyu laya
        VURMAĞA başlayır. Tək perforasiyalı istismarçı üçün bu, absurddur
        (ölçüldü: neft 268 m³/gün laya vurulurdu). Real quyu belə halda
        sadəcə DAYANIR.

        OPM eyni prinsipi işlədir: BHP idarəsi mümkün olmayanda quyu
        idarəetməsi dəyişir (bağlanır və ya başqa rejimə keçir).

        VACİB — VAXTLAMA: bu qərar Nyuton İTERASİYALARI ARASINDA
        DEYİL, yalnız zaman addımının ƏVVƏLİNDƏ verilir. Əks halda
        quyu iterasiyadan-iterasiyaya açılıb-bağlanar və Nyuton
        osilyasiya edər — bu, məhz köhnə modeldəki `min(q,0)`
        problemidir (bax `A7_PLAN.md`). Bir addım daxilində quyunun
        vəziyyəti SABİTDİR.

        Qaytarır: vəziyyət dəyişibsə `True`.
        """
        changed = False
        for connection in self.connections:
            name = connection.well_name
            if connection.is_injector or connection.mode is not ControlMode.BHP:
                continue
            bhp = wells.bhp_of(name)
            would_inject = bhp > state.pressure[connection.cell]
            if would_inject != self.shut.get(name, False):
                self.shut[name] = would_inject
                changed = True
        return changed

    # ── perforasiya debitləri ──────────────────────────────────────
    def perforation_rates(self, state: ThreePhaseState,
                          fluid: ThreePhaseFluidState,
                          wells: WellUnknowns) -> PerforationRates:
        """Hər perforasiya üçün üç fazalı debit — KƏSMƏSİZ.

        Debitin işarəsi `p_bhp − p_hüceyrə` fərqindən təbii çıxır;
        heç bir `min`/`max` yoxdur, ona görə funksiya hamardır.
        """
        water = np.zeros(self.ncell)
        oil = np.zeros(self.ncell)
        gas = np.zeros(self.ncell)
        per_water = {name: 0.0 for name in self.names}
        per_oil = dict(per_water)
        per_gas = dict(per_water)

        for connection in self.connections:
            if self.shut.get(connection.well_name, False):
                continue                      # bağlı quyu debit vermir
            cell = connection.cell
            bhp = wells.bhp_of(connection.well_name)
            drawdown = bhp - state.pressure[cell]
            transmissibility = connection.well_index

            if drawdown > 0.0:
                # laya DAXİL olur — vurulan fazanın mobilliyi
                # (hazırda yalnız su vurucusu dəstəklənir; qaz/CO₂
                # vurma EOR-un öz mövzusudur)
                mobility_w = (self._endpoint_water_mobility
                             / fluid.mu_w[cell]) if connection.is_injector \
                    else fluid.lam_w[cell]
                mobility_o = 0.0 if connection.is_injector else fluid.lam_o[cell]
                mobility_g = 0.0 if connection.is_injector else fluid.lam_g[cell]
            else:
                # laydan ÇIXIR — hüceyrənin öz mobillikləri
                mobility_w = fluid.lam_w[cell]
                mobility_o = fluid.lam_o[cell]
                mobility_g = fluid.lam_g[cell]

            qw = transmissibility * mobility_w * drawdown / fluid.bw[cell]
            qo = transmissibility * mobility_o * drawdown / fluid.bo[cell]
            q_free_gas = (transmissibility * mobility_g * drawdown
                          / fluid.bg[cell])
            # həll olmuş qaz neftlə birlikdə hərəkət edir (bax 6a/6b)
            qg = q_free_gas + fluid.rs[cell] * qo

            water[cell] += qw
            oil[cell] += qo
            gas[cell] += qg
            per_water[connection.well_name] += qw
            per_oil[connection.well_name] += qo
            per_gas[connection.well_name] += qg

        return PerforationRates(water, oil, gas, per_water, per_oil, per_gas)

    # ── quyu idarəetmə tənlikləri ──────────────────────────────────
    def control_residuals(self, rates: PerforationRates,
                          wells: WellUnknowns) -> np.ndarray:
        """Hər quyu üçün bir qalıq — vektorun quyu hissəsi.

            BHP idarəsində:   R = p_bhp − p_hədəf
            RATE idarəsində:  R = Σ q_maye − q_hədəf

        RATE hədəfi MAYE (su+neft) debitidir — A6-dakı konvensiya
        (bax mərhələ 6b): qaz nəticə kimi çıxır, hədəf kimi yox.

        MİQYAS QEYDİ: BHP qalığı bar vahidindədir, RATE qalığı isə
        m³/gün — eyni vektorda çox fərqli böyüklüklər. Xətti həlledici
        üçün bu, pis şərtlənmə deməkdir, ona görə RATE qalığı quyu
        indeksinə bölünüb təzyiq miqyasına gətirilir.
        """
        residuals = np.zeros(wells.count)
        target_by_well: Dict[str, tuple] = {}
        scale_by_well: Dict[str, float] = {}

        for connection in self.connections:
            name = connection.well_name
            if name not in target_by_well:
                target_by_well[name] = (connection.mode, connection.target)
                scale_by_well[name] = 0.0
            scale_by_well[name] += connection.well_index

        for position, name in enumerate(wells.names):
            mode, target = target_by_well[name]
            if self.shut.get(name, False):
                # Bağlı quyunun BHP-si sərbəstdir — onu lay təzyiqinə
                # bərabər saxlayırıq ki, sistem təkil olmasın və quyu
                # yenidən açılmağa hazır qalsın.
                residuals[position] = 0.0
                continue
            if mode is ControlMode.BHP:
                residuals[position] = wells.bhp[position] - target
            else:
                liquid = rates.per_well_water[name] + rates.per_well_oil[name]
                scale = max(scale_by_well[name], 1e-12)
                residuals[position] = (liquid - _signed_rate_target(
                    target, name, self.connections)) / scale
        return residuals

    def is_bhp_controlled(self, name: str) -> bool:
        for connection in self.connections:
            if connection.well_name == name:
                return connection.mode is ControlMode.BHP
        raise KeyError(name)


def _signed_rate_target(target: float, name: str, connections) -> float:
    """RATE hədəfini işarə konvensiyamıza gətirir.

    İstifadəçi hədəfi MÜSBƏT böyüklük kimi verir (məs. "50 m³/gün
    hasilat"). Daxili konvensiyada isə hasilat MƏNFİdir.
    """
    is_injector = any(c.well_name == name and c.is_injector
                      for c in connections)
    magnitude = abs(target)
    return magnitude if is_injector else -magnitude


@dataclass
class WellJacobianBlocks:
    """Quyu tənliklərinin törəmələri — dörd hissə.

    Birləşmiş sistem belə görünür (quyular vektorun sonundadır):

        ┌─────────────┬──────────┐  ┌───┐   ┌───┐
        │  rezervuar  │ R↔Q      │  │ δx│   │ R │
        │   (3N×3N)   │ (3N×W)   │  │   │ = │   │
        ├─────────────┼──────────┤  ├───┤   ├───┤
        │  Q↔R (W×3N) │ Q (W×W)  │  │δbhp│   │Rc │
        └─────────────┴──────────┘  └───┘   └───┘

    · `rate_wrt_reservoir` — perforasiya debitlərinin hüceyrə
      dəyişənlərinə görə törəməsi (rezervuar sətirlərinə DİAQONAL
      töhfə, çünki quyu yalnız öz hüceyrəsinə təsir edir)
    · `rate_wrt_bhp`       — sağ yuxarı blok (R↔Q)
    · `control_wrt_reservoir` — sol aşağı blok (Q↔R)
    · `control_wrt_bhp`    — sağ aşağı blok (Q)
    """
    rate_wrt_reservoir: np.ndarray      # (ncell, 3, 3)
    rate_wrt_bhp: Dict[int, np.ndarray] # hüceyrə -> (3,) hər quyu üçün
    rate_bhp_owner: Dict[int, int]      # hüceyrə -> quyunun yerli nömrəsi
    control_wrt_reservoir: Dict[int, Dict[int, np.ndarray]]  # quyu -> hüceyrə -> (3,)
    control_wrt_bhp: np.ndarray         # (W,)


class StandardWellJacobian:
    """Quyu tənliklərinin analitik törəmələri — MƏRHƏLƏ 3.

    Debit düsturu (bax mərhələ 2):

        q_α = WI · λ_α(Sw, x) / B_α(p) · (p_bhp − p_hüceyrə)

    Törəmələr:

        ∂q/∂p_bhp     = WI · λ/B                      ← YENİ bağlantı
        ∂q/∂p_hüceyrə = WI · [−λ/B + Δp · d(λ/B)/dp]
        ∂q/∂Sw        = WI · (dkr/dSw)/(μ·B) · Δp
        ∂q/∂x         = WI · (dkr/dx)/(μ·B) · Δp

    `∂q/∂p_bhp` köhnə modeldə YOX idi — BHP orada sabit idi. Məhz bu
    bağlantı quyunu sistemin bir hissəsinə çevirir və kəsilməzliyi
    təmin edir.

    QEYD — mobilliyin keçidi: axın istiqaməti dəyişəndə hansı
    mobilliyin işlədildiyi dəyişir (bax `perforation_rates`).
    Δp = 0 nöqtəsində debit hər iki halda sıfırdır (funksiya kəsilməz),
    lakin TÖRƏMƏ sıçrayır. Bu, OPM-də də belədir və `min`-dəki
    sıçrayışdan qat-qat zəifdir: orada FUNKSİYANIN özü sınırdı.
    """

    def __init__(self, well_model: StandardWellModel, pvt, relperm):
        self.well_model = well_model
        self.pvt = pvt
        self.relperm = relperm

    def blocks(self, state: ThreePhaseState, fluid: ThreePhaseFluidState,
               wells: WellUnknowns) -> WellJacobianBlocks:
        ncell = self.well_model.ncell
        rate_wrt_reservoir = np.zeros((ncell, 3, 3))
        rate_wrt_bhp: Dict[int, np.ndarray] = {}
        rate_bhp_owner: Dict[int, int] = {}
        control_wrt_reservoir: Dict[int, Dict[int, np.ndarray]] = {
            index: {} for index in range(wells.count)}
        control_wrt_bhp = np.zeros(wells.count)

        pressure = state.pressure
        mu_w_p = self.pvt.water_viscosity_derivative(pressure)
        mu_o_p = self.pvt.oil_viscosity_derivative(pressure)
        mu_g_p = self.pvt.gas_viscosity_derivative(pressure)
        bw_p = self.pvt.water_fvf_derivative(pressure)
        bo_p = self.pvt.oil_fvf_derivative(pressure)
        bg_p = self.pvt.gas_fvf_derivative(pressure)
        rs_sat_p = self.pvt.solution_gor_derivative(pressure)

        scale_by_well = self._control_scales(wells)

        for connection in self.well_model.connections:
            if self.well_model.shut.get(connection.well_name, False):
                # bağlı quyu: debit sıfır (törəmələr də), BHP tənliyi
                # isə sadəcə diaqonal 1 (təkillik olmasın)
                control_wrt_bhp[wells.index_of(connection.well_name)] = 1.0
                continue
            cell = connection.cell
            well_position = wells.index_of(connection.well_name)
            transmissibility = connection.well_index
            bhp = wells.bhp[well_position]
            drawdown = bhp - pressure[cell]
            injecting = drawdown > 0.0 and connection.is_injector

            saturated = bool(state.is_saturated[cell])
            derivatives = self._relperm_derivatives(state, cell)

            if injecting:
                # vurulan faza: mobillik doyumluluqdan ASILI DEYİL
                # (son nöqtə mobilliyi), ona görə Sw/x törəmələri sıfır
                kr_w = self.well_model._endpoint_water_mobility
                dkrw_dsw = 0.0
                kr_o = kr_g = 0.0
                dkro_dsw = dkro_dx = dkrg_dx = 0.0
            else:
                kr_w, kr_o, kr_g = fluid.krw[cell], fluid.kro[cell], fluid.krg[cell]
                dkrw_dsw = derivatives["krw_dsw"]
                dkro_dsw = derivatives["kro_dsw"]
                dkro_dx = derivatives["kro_dx"] if saturated else 0.0
                dkrg_dx = derivatives["krg_dx"] if saturated else 0.0

            transport_w = kr_w / (fluid.mu_w[cell] * fluid.bw[cell])
            transport_o = kr_o / (fluid.mu_o[cell] * fluid.bo[cell])
            transport_g = kr_g / (fluid.mu_g[cell] * fluid.bg[cell])

            d_transport_w = -kr_w * (mu_w_p[cell] * fluid.bw[cell]
                                     + fluid.mu_w[cell] * bw_p[cell]) \
                / (fluid.mu_w[cell] * fluid.bw[cell]) ** 2
            d_transport_o = -kr_o * (mu_o_p[cell] * fluid.bo[cell]
                                     + fluid.mu_o[cell] * bo_p[cell]) \
                / (fluid.mu_o[cell] * fluid.bo[cell]) ** 2
            d_transport_g = -kr_g * (mu_g_p[cell] * fluid.bg[cell]
                                     + fluid.mu_g[cell] * bg_p[cell]) \
                / (fluid.mu_g[cell] * fluid.bg[cell]) ** 2

            # ── hüceyrə dəyişənlərinə görə ──────────────────────────
            dqw_dp = transmissibility * (-transport_w + drawdown * d_transport_w)
            dqo_dp = transmissibility * (-transport_o + drawdown * d_transport_o)
            dqfree_dp = transmissibility * (-transport_g
                                            + drawdown * d_transport_g)

            dqw_dsw = (transmissibility * dkrw_dsw * drawdown
                      / (fluid.mu_w[cell] * fluid.bw[cell]))
            dqo_dsw = (transmissibility * dkro_dsw * drawdown
                      / (fluid.mu_o[cell] * fluid.bo[cell]))
            dqo_dx = (transmissibility * dkro_dx * drawdown
                     / (fluid.mu_o[cell] * fluid.bo[cell]))
            dqfree_dx = (transmissibility * dkrg_dx * drawdown
                        / (fluid.mu_g[cell] * fluid.bg[cell]))

            # qaz = sərbəst + Rs·neft (hasil qaydası, bax 6c/3)
            rs = fluid.rs[cell]
            qo_value = transmissibility * transport_o * drawdown
            drs_dp = rs_sat_p[cell] if saturated else 0.0
            drs_dx = 0.0 if saturated else 1.0

            rate_wrt_reservoir[cell, 0, 0] += dqw_dp
            rate_wrt_reservoir[cell, 0, 1] += dqw_dsw
            rate_wrt_reservoir[cell, 1, 0] += dqo_dp
            rate_wrt_reservoir[cell, 1, 1] += dqo_dsw
            rate_wrt_reservoir[cell, 1, 2] += dqo_dx
            rate_wrt_reservoir[cell, 2, 0] += (dqfree_dp + drs_dp * qo_value
                                               + rs * dqo_dp)
            rate_wrt_reservoir[cell, 2, 1] += rs * dqo_dsw
            rate_wrt_reservoir[cell, 2, 2] += (dqfree_dx + rs * dqo_dx
                                               + drs_dx * qo_value)

            # ── BHP-yə görə (YENİ bağlantı) ─────────────────────────
            dqw_dbhp = transmissibility * transport_w
            dqo_dbhp = transmissibility * transport_o
            dqg_dbhp = transmissibility * transport_g + rs * dqo_dbhp
            column = np.array([dqw_dbhp, dqo_dbhp, dqg_dbhp])
            if cell in rate_wrt_bhp:
                rate_wrt_bhp[cell] = rate_wrt_bhp[cell] + column
            else:
                rate_wrt_bhp[cell] = column
                rate_bhp_owner[cell] = well_position

            # ── idarəetmə tənliyinin törəmələri ─────────────────────
            if self.well_model.is_bhp_controlled(connection.well_name):
                control_wrt_bhp[well_position] = 1.0
            else:
                scale = scale_by_well[connection.well_name]
                control_wrt_bhp[well_position] += (dqw_dbhp + dqo_dbhp) / scale
                existing = control_wrt_reservoir[well_position].get(
                    cell, np.zeros(3))
                control_wrt_reservoir[well_position][cell] = existing + np.array([
                    (dqw_dp + dqo_dp) / scale,
                    (dqw_dsw + dqo_dsw) / scale,
                    dqo_dx / scale])

        return WellJacobianBlocks(rate_wrt_reservoir, rate_wrt_bhp,
                                  rate_bhp_owner, control_wrt_reservoir,
                                  control_wrt_bhp)

    def _control_scales(self, wells: WellUnknowns) -> Dict[str, float]:
        scales = {name: 0.0 for name in wells.names}
        for connection in self.well_model.connections:
            scales[connection.well_name] += connection.well_index
        return {name: max(value, 1e-12) for name, value in scales.items()}

    def _relperm_derivatives(self, state: ThreePhaseState, cell: int) -> Dict:
        sw = np.array([state.water_saturation[cell]])
        sg = np.array([state.gas_saturation[cell]])
        kro_dsw, kro_dsg = self.relperm.kro_three_phase_derivatives(sw, sg)
        krw_dsw = self.relperm.krw_derivative(sw)
        krg_dsg = self.relperm.gas.krg_derivative(sg, self.relperm.swc)
        return {"krw_dsw": float(krw_dsw[0]),
                "kro_dsw": float(kro_dsw[0]),
                "kro_dx": float(kro_dsg[0]),
                "krg_dx": float(krg_dsg[0])}
