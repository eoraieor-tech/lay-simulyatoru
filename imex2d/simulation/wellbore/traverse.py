"""Çoxseqmentli təzyiq traversi — quyu dibindən quyu başına.

TƏNLİK (şaquli, aşağı istiqamət müsbət):

    dp/dz = [ ρ_qrav·g  +  f·ρ_ns·v_m²/(2·d) ] / 1e5        [bar/m]

Yuxarı gedərkən HƏR İKİ hədd təzyiqi AZALDIR: qravitasiya sütunun
çəkisidir, sürtünmə isə axın istiqamətinə (yuxarı) əks işləyir. Ona görə

    THP = BHP − Σ (Δp_qravitasiya + Δp_sürtünmə)

SÜRƏTLƏNMƏ HƏDDİ (`ρ·v·dv/dz`) DAXİL DEYİL. O, yalnız quyu başına yaxın,
qaz sürətlə genişlənəndə nəzərə çarpır və adətən ümumi itkinin < 1 %-idir.
⏳ V2-yə saxlanılır — sənəddə açıq göstərilir, uydurulmur.

NİYƏ ÇOX SEQMENT. Qaz yuxarı qalxdıqca təzyiq düşür və qaz genişlənir:
`Bg` qat-qat böyüyür, ona görə ρ_m azalır, v_m isə artır. Tək seqmentdə
(orta təzyiqdə) hesablamaq yalnız qazsız quyuda düzgündür. B4b-dən sonra
modeldə REAL sərbəst qaz var (maks Sg 0.11-ə çatır), ona görə lülə
seqmentlərə bölünür və hər seqmentdə PVT YERLİ təzyiqdə yenidən oxunur.

Hər seqmentdə təzyiq özü-özündən asılıdır (xassələr `p_orta`-dan
asılıdır, `p_orta` isə Δp-dən). Sabit-nöqtə iterasiyası ilə həll olunur —
bir neçə dövrdə yığılır, çünki Δp seqment boyunca kiçikdir.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from ...domain.tubing import TubingGeometry
from .friction import chen_friction_factor, reynolds
from .holdup import IHoldupCorrelation, NoSlipHoldup

#: Ağırlıq təcili, m/s².
GRAVITY = 9.80665

#: 1 bar = 1e5 Pa.
PA_PER_BAR = 1.0e5

#: m³/gün → m³/san.
SECONDS_PER_DAY = 86400.0

#: Bundan aşağı ümumi debit "quyu axmır" sayılır — axan traverse orada
#: təyin olunmayıb (dayanmış quyunun quyu başı təzyiqi AYRI hesabatdır).
MIN_RATE_M3_DAY = 1e-9

#: Traversin dayandığı aşağı hədd, bar. Bundan aşağı düşmək quyunun
#: səthə axa bilməməsi deməkdir (bax `pressure_traverse`).
MIN_PRESSURE_BAR = 1e-3

#: Seqment daxilindəki sabit-nöqtə iterasiyası.
MAX_SEGMENT_ITERATIONS = 12
SEGMENT_TOLERANCE_BAR = 1e-7


@dataclass
class WellStream:
    """Bir quyunun SƏTH debitləri (m³/gün) və səth sıxlıqları (kg/m³).

    İşarə konvensiyası: istismarçıda hamısı MÜSBƏT (hasilat). Bu,
    `SimulationResult.well_oil_rate`-in saxladığı formadır — mühərrik
    daxilində mənfi olan dəyər orada artıq işarəsi dəyişdirilib.

    `gas` — ÜMUMİ səth qazı (sərbəst + neftdən ayrılan). Lülədə axan
    kütlə bu dəyərlə hesablanır; sərbəst qazın YERİNDƏ həcmi isə
    `Rs`-i çıxmaqla alınır.
    """

    oil: float = 0.0
    water: float = 0.0
    gas: float = 0.0
    oil_density: float = 850.0
    water_density: float = 1010.0
    gas_density: float = 0.85

    @property
    def total_liquid(self) -> float:
        return self.oil + self.water

    @property
    def mass_rate(self) -> float:
        """Ümumi kütlə axını, kq/gün — lülə boyunca SABİTDİR.

        Həll olmuş qazın kütləsi `gas`-ın içindədir (səthdə ayrılır),
        ona görə burada ayrıca sayılmır — ikiqat sayma olardı.
        """
        return (self.oil * self.oil_density
                + self.water * self.water_density
                + self.gas * self.gas_density)

    def is_flowing(self) -> bool:
        return (abs(self.oil) + abs(self.water) + abs(self.gas)) > MIN_RATE_M3_DAY


@dataclass
class TraverseSegment:
    """Bir seqmentin diaqnostikası — testlər və jurnal üçün."""
    top_depth: float
    bottom_depth: float
    pressure_top: float
    pressure_bottom: float
    mixture_density: float
    mixture_velocity: float
    liquid_holdup: float
    friction_factor: float
    reynolds: float
    gravity_drop: float
    friction_drop: float


@dataclass
class TraverseResult:
    """Traversin nəticəsi."""
    thp: float
    bhp: float
    length: float
    gravity_drop: float
    friction_drop: float
    converged: bool = True
    message: str = ""
    segments: List[TraverseSegment] = field(default_factory=list)

    @property
    def total_drop(self) -> float:
        return self.gravity_drop + self.friction_drop


class _FluidLookup:
    """PVT-yə vahid müraciət — provider varsa ondan, yoxsa sabitlərdən.

    Mühərrikin qalan hissəsi kimi: `IPVTProvider` inject olunmayanda
    `FluidProperties` placeholder dəyərləri işlədilir (bax
    `domain/properties.py::FluidProperties`).
    """

    def __init__(self, pvt=None, fluids=None):
        self.pvt = pvt
        self.fluids = fluids
        self.has_gas = bool(pvt is not None and pvt.has_gas_phase())

    def at(self, pressure: float) -> dict:
        p = np.array([float(pressure)], dtype=float)
        if self.pvt is None:
            f = self.fluids
            return dict(bo=float(f.oil_fvf), bw=float(f.water_fvf), bg=1.0, rs=0.0,
                        mu_o=float(f.oil_viscosity), mu_w=float(f.water_viscosity),
                        mu_g=1e-2)
        out = dict(
            bo=float(self.pvt.oil_fvf(p)[0]),
            bw=float(self.pvt.water_fvf(p)[0]),
            mu_o=float(self.pvt.oil_viscosity(p)[0]),
            mu_w=float(self.pvt.water_viscosity(p)[0]),
            bg=1.0, rs=0.0, mu_g=1e-2)
        if self.has_gas:
            out["bg"] = float(self.pvt.gas_fvf(p)[0])
            out["rs"] = float(self.pvt.solution_gor(p)[0])
            out["mu_g"] = float(self.pvt.gas_viscosity(p)[0])
        return out

    def oil_fvf_at(self, pressure: float, rs_actual: float,
                   rs_saturated: float, props: dict) -> float:
        """Bo — neftin HƏQİQİ həll olmuş qazına görə.

        `rs_actual < rs_saturated` olduqda neft DOYMAMIŞDIR və Bo
        cədvəlin doymuş qolundan OXUNA BİLMƏZ (bax
        `BlackOilPVTProvider.oil_fvf_undersaturated` — B3-B-də düzəldilən
        eyni səhv). Provider doymamış qolu dəstəkləmirsə doymuş qiymət
        qalır və bu, mühafizəkar sadələşdirmə kimi sənədləşir.
        """
        if rs_actual >= rs_saturated - 1e-12:
            return props["bo"]
        branch = getattr(self.pvt, "oil_fvf_undersaturated", None)
        if branch is None:
            return props["bo"]
        return float(branch(np.array([float(pressure)]),
                            np.array([float(rs_actual)]))[0])


def _segment_gradient(pressure: float, stream: WellStream,
                      tubing: TubingGeometry, fluid: _FluidLookup,
                      holdup: IHoldupCorrelation) -> dict:
    """Verilmiş təzyiqdə qradiyentin komponentləri, bar/m."""
    props = fluid.at(pressure)
    area = tubing.area

    # ── neftin HƏQİQİ həll olmuş qazı ───────────────────────────────
    # Axının hasilat GOR-u neftin nə qədər qaz daşıdığını söyləyir.
    # `Rs_sat(p)` yalnız YUXARI HƏDDDİR: ondan artığı sərbəst qazdır,
    # azı isə neftin DOYMAMIŞ olması deməkdir.
    #
    # Bo-nu buna görə seçmək VACİBDİR (ölçüldü — B3-B-dəki səhvin
    # eynisi): doymuş cədvəldən oxunan Bo, doymamış neftə HƏDDİNDƏN
    # ARTIQ şişmə verir, sütun yüngülləşir və THP qaz artdıqca
    # QEYRİ-MONOTON çıxır.
    rs_saturated = props["rs"]
    producing_gor = (stream.gas / stream.oil
                     if stream.oil > MIN_RATE_M3_DAY else float("inf"))
    rs_actual = min(producing_gor, rs_saturated)
    bo = fluid.oil_fvf_at(pressure, rs_actual, rs_saturated, props)

    # ── yerində (in-situ) həcm debitləri, m³/gün ────────────────────
    q_oil = stream.oil * bo
    q_water = stream.water * props["bw"]
    free_gas_surface = max(stream.gas - rs_actual * stream.oil, 0.0)
    q_gas = free_gas_surface * props["bg"]

    q_liquid = q_oil + q_water
    q_total = q_liquid + q_gas
    if q_total <= MIN_RATE_M3_DAY:
        return dict(gravity=0.0, friction=0.0, density=0.0, velocity=0.0,
                    holdup=1.0, friction_factor=0.0, reynolds=0.0)

    # ── faza sıxlıqları YERİNDƏ, kq/m³ ──────────────────────────────
    # Mayenin kütləsinə həll olmuş qaz DA daxildir (o, lay şəraitində
    # neftin içindədir və yalnız səthdə ayrılır).
    dissolved_gas_surface = stream.gas - free_gas_surface
    liquid_mass = (stream.oil * stream.oil_density
                   + stream.water * stream.water_density
                   + dissolved_gas_surface * stream.gas_density)
    rho_liquid = liquid_mass / q_liquid if q_liquid > MIN_RATE_M3_DAY else 0.0
    gas_mass = free_gas_surface * stream.gas_density
    rho_gas = gas_mass / q_gas if q_gas > MIN_RATE_M3_DAY else 0.0

    # ── sürətlər, m/s ───────────────────────────────────────────────
    v_liquid = q_liquid / (SECONDS_PER_DAY * area)
    v_gas = q_gas / (SECONDS_PER_DAY * area)
    v_mixture = v_liquid + v_gas

    lambda_liquid = q_liquid / q_total
    h_liquid = holdup.liquid_holdup(lambda_liquid, v_liquid, v_gas,
                                    tubing.diameter, 0.0)

    # QRAVİTASİYA sürüşməli sıxlıqla, SÜRTÜNMƏ sürüşməsiz sıxlıqla —
    # sənaye standartı budur. Sürüşməsiz halda ikisi üst-üstə düşür.
    rho_gravity = rho_liquid * h_liquid + rho_gas * (1.0 - h_liquid)
    rho_no_slip = rho_liquid * lambda_liquid + rho_gas * (1.0 - lambda_liquid)

    mu_liquid = ((q_oil * props["mu_o"] + q_water * props["mu_w"]) / q_liquid
                 if q_liquid > MIN_RATE_M3_DAY else props["mu_o"])
    mu_no_slip = mu_liquid * lambda_liquid + props["mu_g"] * (1.0 - lambda_liquid)

    re = reynolds(rho_no_slip, v_mixture, tubing.diameter, mu_no_slip)
    f = chen_friction_factor(re, tubing.relative_roughness)

    gravity = rho_gravity * GRAVITY / PA_PER_BAR
    friction = (f * rho_no_slip * v_mixture * v_mixture
                / (2.0 * tubing.diameter) / PA_PER_BAR)
    return dict(gravity=gravity, friction=friction, density=rho_gravity,
                velocity=v_mixture, holdup=h_liquid, friction_factor=f,
                reynolds=re)


def pressure_traverse(bhp: float, perforation_depth: float,
                      tubing: TubingGeometry, stream: WellStream,
                      pvt=None, fluids=None,
                      holdup: Optional[IHoldupCorrelation] = None,
                      keep_segments: bool = False) -> TraverseResult:
    """Quyu dibi təzyiqindən quyu başı təzyiqini hesablayır.

    `bhp` — quyu dibi təzyiqi, bar (perforasiya dərinliyində).
    `perforation_depth` — perforasiyanın dərinliyi, m.
    `keep_segments` — diaqnostika üçün hər seqmenti saxlamaq.

    Qaytarır: `TraverseResult`. Quyu axmırsa (`stream.is_flowing()` yalan),
    `thp = nan` və `converged=False` qaytarılır — axan traverse dayanmış
    quyu üçün TƏYİN OLUNMAYIB, uydurma dəyər verilmir.
    """
    holdup = holdup or NoSlipHoldup()
    length = tubing.length_to(perforation_depth)

    if not stream.is_flowing():
        return TraverseResult(thp=float("nan"), bhp=float(bhp), length=length,
                              gravity_drop=0.0, friction_drop=0.0,
                              converged=False,
                              message="quyu axmır — axan traverse təyin olunmayıb")
    if length <= 0.0:
        return TraverseResult(thp=float(bhp), bhp=float(bhp), length=length,
                              gravity_drop=0.0, friction_drop=0.0,
                              message="perforasiya quyu başından aşağıda deyil")

    fluid = _FluidLookup(pvt, fluids)
    n = max(int(tubing.segments), 1)
    dz = length / n

    pressure = float(bhp)
    depth = float(perforation_depth)
    gravity_total = 0.0
    friction_total = 0.0
    segments: List[TraverseSegment] = []
    converged = True
    message = ""

    for _ in range(n):
        # Seqment ORTASINDAKI təzyiqdə xassələr — Δp-dən asılıdır, ona
        # görə sabit-nöqtə ilə həll olunur.
        drop = 0.0
        grad = None
        for _iteration in range(MAX_SEGMENT_ITERATIONS):
            midpoint = max(pressure - drop / 2.0, 1e-3)
            grad = _segment_gradient(midpoint, stream, tubing, fluid, holdup)
            new_drop = (grad["gravity"] + grad["friction"]) * dz
            if abs(new_drop - drop) < SEGMENT_TOLERANCE_BAR:
                drop = new_drop
                break
            drop = new_drop
        else:
            converged = False
            message = "seqment daxilində sabit-nöqtə iterasiyası yığılmadı"

        gravity_drop = grad["gravity"] * dz
        friction_drop = grad["friction"] * dz
        top_pressure = pressure - gravity_drop - friction_drop

        if top_pressure <= MIN_PRESSURE_BAR:
            # Quyu SƏTHƏ AXA BİLMİR: sütunun çəkisi quyu dibi təzyiqini
            # üstələyir. Bu, ədədi qüsur DEYİL, fiziki nəticədir —
            # belə quyu süni qaldırma (qazlift, nasos) tələb edir.
            #
            # Sıfıra "qısaldılmış" dəyər qaytarmaq YANLIŞ siqnal olardı
            # (qrafikdə 0 bar real ölçmə kimi görünərdi), ona görə `nan`
            # qaytarılır və səbəb mesajda yazılır.
            return TraverseResult(
                thp=float("nan"), bhp=float(bhp), length=length,
                gravity_drop=gravity_total + gravity_drop,
                friction_drop=friction_total + friction_drop,
                converged=False, segments=segments,
                message=("quyu səthə axa bilmir — lülədəki təzyiq itkisi "
                         f"({gravity_total + gravity_drop + friction_total + friction_drop:.1f} bar) "
                         f"quyu dibi təzyiqini ({bhp:.1f} bar) üstələyir; "
                         "süni qaldırma lazımdır"))

        if keep_segments:
            segments.append(TraverseSegment(
                top_depth=depth - dz, bottom_depth=depth,
                pressure_top=top_pressure, pressure_bottom=pressure,
                mixture_density=grad["density"], mixture_velocity=grad["velocity"],
                liquid_holdup=grad["holdup"], friction_factor=grad["friction_factor"],
                reynolds=grad["reynolds"], gravity_drop=gravity_drop,
                friction_drop=friction_drop))

        gravity_total += gravity_drop
        friction_total += friction_drop
        pressure = top_pressure
        depth -= dz

    return TraverseResult(thp=pressure, bhp=float(bhp), length=length,
                          gravity_drop=gravity_total, friction_drop=friction_total,
                          converged=converged, message=message, segments=segments)
