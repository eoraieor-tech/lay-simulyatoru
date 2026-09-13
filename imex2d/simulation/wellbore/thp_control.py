"""Quyunu quyu başı təzyiqi (THP) ilə idarə etmək — B4-B.

AÇIQ (EXPLICIT) BİRLƏŞMƏ. Nyuton həlledicisi THP-ni tanımır — yalnız
quyu dibi təzyiqini (BHP). Ona görə THP quyusu bağlantı səviyyəsində
ADİ BHP bağlantısı kimi qurulur (`well_model.build_connections`), bu
modul isə hər zaman addımından sonra son debitlərlə traversi TƏRSİNƏ
həll edib növbəti addımın BHP-sini yeniləyir:

    THP (istifadəçi)  +  son addımın debitləri  →  BHP (növbəti addım)

NİYƏ QALIQ/JAKOBİAN TOXUNULMUR. Mühərrikdə `connection.mode is BHP`
yoxlaması qalıqda, Jakobianda, üç fazalı qalıqda, IMPES-də və daha bir
neçə yerdə var. Yeni rejimi bağlantıya ötürmək onların HAMISINI səssizcə
RATE budağına göndərərdi. THP quyusunu BHP bağlantısı kimi saxlamaq bu
riski tamamilə aradan qaldırır.

AÇIQ BİRLƏŞMƏNİN BİLİNƏN RİSKİ — gecikmə. BHP bir addım əvvəlki debitlərə
görə hesablanır; debit sürətlə dəyişəndə BHP rəqs edə bilər. Qarşısı iki
yolla alınır: relaksasiya (`RELAXATION`) və addım başına maksimal dəyişmə
(`MAX_BHP_CHANGE_BAR`). Tam implicit THP birləşməsi ⏳ sonraya.

V1 MƏHDUDİYYƏTLƏRİ (B4-A ilə eyni): yalnız istismarçılar, tam şaquli
lülə, sürüşmə yoxdur.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional

from ...domain.wells import ControlMode, WellType
from ...logging_setup import get_logger
from .holdup import IHoldupCorrelation, NoSlipHoldup
from .traverse import WellStream, pressure_traverse

LOG = get_logger(__name__)

#: Axınsız (statik) təxmin üçün cüzi neft debiti, m³/gün. Traverse
#: axmayan quyu üçün `nan` qaytarır (Q-11), ona görə ilk addımda və quyu
#: dayananda sütun praktik olaraq sürtünməsiz neft sütunu kimi götürülür.
STATIC_PROBE_RATE = 1e-3

#: İkiqat bölmənin dəqiqliyi, bar.
BISECTION_TOLERANCE_BAR = 1e-3
MAX_BISECTION_ITERATIONS = 80

#: BHP axtarışının yuxarı həddi, bar — bundan böyük BHP fiziki deyil.
BHP_SEARCH_LIMIT_BAR = 2000.0

#: Açıq birləşmə relaksasiyası: yeni BHP = köhnə + ω·(hədəf − köhnə).
RELAXATION = 0.5

#: Bir addımda BHP-nin maksimal dəyişməsi, bar — rəqsin qarşısını alır.
MAX_BHP_CHANGE_BAR = 25.0


def bhp_from_thp(thp: float, perforation_depth: float, tubing, stream: WellStream,
                 pvt=None, fluids=None,
                 holdup: Optional[IHoldupCorrelation] = None,
                 bhp_limit: float = BHP_SEARCH_LIMIT_BAR) -> float:
    """Verilmiş THP üçün BHP — `pressure_traverse`-in TƏRSİ.

    Sabit axın üçün THP(BHP) artan funksiyadır (BHP böyüdükcə quyu başına
    daha çox təzyiq çatır), ona görə ikiqat bölmə ilə dəqiq tapılır.
    Aşağı hədd `THP`-nin özüdür: lülədəki təzyiq itkisi mənfi ola bilməz.

    Qaytarır: BHP (bar) və ya `nan` — `bhp_limit`-ə qədər heç bir BHP
    istənilən THP-ni vermirsə (quyu belə THP-də axa bilmir).
    """
    target = float(thp)

    def thp_at(bhp: float) -> float:
        return pressure_traverse(bhp, perforation_depth, tubing, stream,
                                 pvt=pvt, fluids=fluids, holdup=holdup).thp

    if target >= bhp_limit:
        return float("nan")
    low = target
    high = min(max(target + 1.0, 2.0 * target), bhp_limit)
    value = thp_at(high)
    while (not math.isfinite(value) or value < target) and high < bhp_limit:
        high = min(1.5 * high + 10.0, bhp_limit)
        value = thp_at(high)
    if not math.isfinite(value) or value < target:
        return float("nan")

    for _ in range(MAX_BISECTION_ITERATIONS):
        middle = 0.5 * (low + high)
        value = thp_at(middle)
        if not math.isfinite(value) or value < target:
            low = middle
        else:
            high = middle
        if high - low < BISECTION_TOLERANCE_BAR:
            break
    return high


class ThpController:
    """THP quyularının BHP hədəflərini addım-addım yeniləyir.

    `connections` — mühərrikin qalığının OXUDUĞU bağlantı obyektləri.
    Hədəf yerində (`connection.target`) dəyişdirilir, yəni qalıq və
    Jakobian növbəti qiymətləndirmədə yeni BHP-ni görür.
    """

    def __init__(self, model, connections, pvt=None, fluids=None,
                 holdup: Optional[IHoldupCorrelation] = None,
                 relaxation: float = RELAXATION,
                 max_change: float = MAX_BHP_CHANGE_BAR):
        self.model = model
        self.pvt = pvt
        self.fluids = fluids if fluids is not None else model.fluids
        self.holdup = holdup or NoSlipHoldup()
        self.relaxation = float(relaxation)
        self.max_change = float(max_change)

        wells = {well.name: well for well in model.active_wells()}
        self._connections: Dict[str, List] = {}
        for connection in connections:
            if getattr(connection, "thp_target", None) is None:
                continue
            well = wells.get(connection.well_name)
            if well is None or well.well_type is not WellType.PRODUCER:
                continue
            self._connections.setdefault(connection.well_name, []).append(connection)
        self._wells = {name: wells[name] for name in self._connections}

        self._depths: Dict[str, float] = {}
        if self._connections:
            from .hydraulics import WellboreHydraulics
            self._depths = WellboreHydraulics._perforation_depths(model)

        #: Cari (növbəti addımda işlədiləcək) BHP, bar.
        self.bhp: Dict[str, float] = {}

    @property
    def active(self) -> bool:
        """Heç bir THP quyusu yoxdursa nəzarətçi heç nə etmir."""
        return bool(self._connections)

    # ─────────────────────────────────────────────────────── köməkçilər
    def _stream(self, oil: float, water: float, gas: float) -> WellStream:
        return WellStream(oil=max(oil, 0.0), water=max(water, 0.0),
                          gas=max(gas, 0.0),
                          oil_density=self.fluids.oil_density,
                          water_density=self.fluids.water_density,
                          gas_density=self.fluids.gas_density)

    def _ready(self, name: str) -> bool:
        well = self._wells.get(name)
        return (well is not None and getattr(well, "tubing", None) is not None
                and name in self._depths)

    def _solve(self, name: str, stream: WellStream) -> float:
        well = self._wells[name]
        return bhp_from_thp(float(well.control.target), self._depths[name],
                            well.tubing, stream, pvt=self.pvt,
                            fluids=self.fluids, holdup=self.holdup)

    def _apply(self, name: str, bhp: float) -> None:
        # BHP heç vaxt THP-dən aşağı ola bilməz — lülədə itki mənfi deyil.
        bhp = max(float(bhp), float(self._wells[name].control.target))
        self.bhp[name] = bhp
        for connection in self._connections[name]:
            connection.target = bhp

    # ─────────────────────────────────────────────────────── dövr
    def initialize(self) -> None:
        """İlk addımdan ƏVVƏL: statik (axınsız) sütunla BHP təxmini.

        Debit hələ məlum deyil — ona görə sütun cüzi neft axını ilə
        hesablanır. Lülə həndəsəsi olmayan quyu atlanılır: diaqnostika
        onu model qurularkən XƏTA kimi göstərir.
        """
        for name in self._connections:
            if not self._ready(name):
                LOG.warning("THP quyusu '%s' üçün lülə həndəsəsi və ya "
                            "perforasiya dərinliyi yoxdur — BHP təyin edilmədi",
                            name)
                continue
            bhp = self._solve(name, self._stream(STATIC_PROBE_RATE, 0.0, 0.0))
            if not math.isfinite(bhp):
                LOG.warning("THP quyusu '%s': statik sütun belə THP = %.1f bar "
                            "verə bilmir — BHP = THP götürüldü", name,
                            float(self._wells[name].control.target))
                bhp = float(self._wells[name].control.target)
            self._apply(name, bhp)

    def update(self, per_well_oil, per_well_water, per_well_gas=None) -> None:
        """Qəbul olunmuş addımdan SONRA: növbəti addımın BHP-si.

        Debit lüğətləri mühərrikin daxili işarəsindədir — hasilat MƏNFİ.
        """
        per_well_gas = per_well_gas or {}
        for name in self._connections:
            if name not in self.bhp:
                continue
            stream = self._stream(-float(per_well_oil.get(name, 0.0)),
                                  -float(per_well_water.get(name, 0.0)),
                                  -float(per_well_gas.get(name, 0.0)))
            if not stream.is_flowing():
                stream = self._stream(STATIC_PROBE_RATE, 0.0, 0.0)
            target = self._solve(name, stream)
            if not math.isfinite(target):
                continue                       # cari BHP saxlanılır
            current = self.bhp[name]
            proposal = current + self.relaxation * (target - current)
            change = max(-self.max_change, min(self.max_change, proposal - current))
            self._apply(name, current + change)

    def record(self, result) -> None:
        """Bu addımda İŞLƏDİLƏN BHP-ni nəticəyə yazır (`update`-dən ƏVVƏL)."""
        for name, value in self.bhp.items():
            result.well_bhp.setdefault(name, []).append(float(value))
