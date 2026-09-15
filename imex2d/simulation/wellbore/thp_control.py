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

#: Addım daxilindəki təkrarın dayanma həddi, bar (yarı-implicit dövr).
OUTER_TOLERANCE_BAR = 1.0

#: Addım daxilində maksimal təkrar sayı. 2 kifayətdir: ölçüldü, üçüncü
#: təkrar BHP-ni 1 bar-dan az dəyişir, hesab xərci isə addımı ikiqat edir.
MAX_OUTER_ITERATIONS = 2

#: Quyunun yenidən açılması üçün tələb olunan ehtiyat, bar. Sərhəddə
#: (sıfır ehtiyatla) açmaq quyunu dərhal yenidən bağlanmağa məcbur edərdi
#: — yəni aç/bağla rəqsi. Ehtiyat bu histerezisi yaradır.
REOPEN_MARGIN_BAR = 2.0

#: Quyu HƏLƏ heç axmayıbsa BHP lay təzyiqindən bu qədər aşağı qoyulur.
#:
#: NİYƏ (ölçüldü, Seans 25): axın olmayanda sütun DURĞUN NEFT kimi
#: hesablanır — bu, ən AĞIR sütundur və dərin quyuda tələb olunan BHP-ni
#: lay təzyiqindən yuxarı atır. Həmin qiymətə baxıb quyunu bağlasaydıq,
#: real şəraitdə işləyən quyu (qaz ayrılan kimi sütun yüngülləşir) heç
#: vaxt işə düşə bilməzdi — ölçüldü: PROD-1 t = 0-da səhvən bağlanırdı.
STARTUP_DRAWDOWN_BAR = 5.0

#: J hesablanarkən etibarlı sayılan minimal drawdown, bar.
MIN_DRAWDOWN_BAR = 0.5

#: Nodal həllin dəqiqliyi (bar) və maksimal ikiqat bölmə sayı.
NODAL_TOLERANCE_BAR = 0.05
MAX_NODAL_ITERATIONS = 40


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
        #: Quyunun SON AXAN tərkibi — axın dayananda sütunu bununla
        #: hesablayırıq. Bax `update()`-dəki BİSTABİLLİK izahı.
        self._last_stream: Dict[str, WellStream] = {}
        #: Hədəf lay təzyiqindən yuxarı düşəndə bir dəfə xəbərdarlıq
        self._warned_shut: Dict[str, bool] = {}
        #: Məhsuldarlıq əmsalı J (IPR mailliyi), m³/gün/bar — hər addımda
        #: son həqiqi iş nöqtəsindən yenilənir.
        self._productivity: Dict[str, float] = {}
        #: AVTOMATİK BAĞLANMA (auto shut-in) — quyu adı → bağlıdırmı
        self.shut: Dict[str, bool] = {name: False for name in self._connections}
        #: Bağlananda sıfırlanan quyu indeksləri burada saxlanılır
        self._well_index: Dict[str, List[float]] = {
            name: [float(c.well_index) for c in connections]
            for name, connections in self._connections.items()}

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
    def initialize(self, pressure=None) -> None:
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
            reference = self._reservoir_pressure(name, pressure)
            if reference is not None:
                # Quyu İŞƏ DÜŞƏ BİLMƏLİDİR: durğun sütunun tələb etdiyi
                # BHP lay təzyiqindən yuxarıdırsa, onu sərhədin bir qədər
                # ALTINA qoyuruq ki, ilk axın başlasın və nəzarətçi həqiqi
                # (qazlı, yüngül) tərkibi öyrənsin.
                bhp = min(bhp, reference - STARTUP_DRAWDOWN_BAR)
            self._apply(name, bhp)

    def update(self, per_well_oil, per_well_water, per_well_gas=None,
               pressure=None) -> float:
        """Addımdan sonra BHP-ni yeniləyir; ƏN BÖYÜK dəyişməni qaytarır (bar).

        Debit lüğətləri mühərrikin daxili işarəsindədir — hasilat MƏNFİ.
        `pressure` — hüceyrə təzyiqləri massivi (verilsə, BHP lay
        təzyiqindən yuxarı qaldırılmır).

        NODAL ANALİZ (Seans 25) — BU METODUN ƏSASI.

        Əvvəl BHP "keçən addımın debiti" ilə hesablanırdı. Bu, prinsipcə
        səhv idi: tələb olunan BHP DEBİTDƏN güclü asılıdır (ölçüldü —
        5 m³/gün-də 71 bar, 500-də 145 bar, 1000-də 245 bar), debitin özü
        isə BHP-dən asılıdır. Nəticədə nəzarətçi VLP əyrisi boyunca
        sıçrayırdı və "quyu axa bilmir" kimi YANLIŞ qərarlar verirdi.

        İndi iş nöqtəsi hasilat mühəndisliyinin standart üsulu ilə —
        IPR və VLP əyrilərinin KƏSİŞMƏSİ ilə tapılır:

            IPR:  q(BHP) = J · (p_lay − BHP)        (düz xətt, J son addımdan)
            VLP:  THP(BHP, q(BHP))                  (çoxseqmentli traverse)
            həll: THP(BHP, q(BHP)) = THP_hədəf      (BHP üzrə ikiqat bölmə)

        BHP artdıqca debit AZALIR (IPR) və sürtünmə itkisi düşür, yəni
        alınan THP MONOTON artır — ona görə kəsişmə YEGANƏDİR və ikiqat
        bölmə etibarlıdır. Hesab xərci əvvəlki ilə eynidir (~40 traverse).
        """
        per_well_gas = per_well_gas or {}
        largest = 0.0
        for name in self._connections:
            if name not in self.bhp:
                continue
            stream = self._stream(-float(per_well_oil.get(name, 0.0)),
                                  -float(per_well_water.get(name, 0.0)),
                                  -float(per_well_gas.get(name, 0.0)))
            flowing = stream.is_flowing()
            reference = self._reservoir_pressure(name, pressure)
            if flowing:
                self._last_stream[name] = stream
                self._update_productivity(name, stream, reference)

            target = self._operating_point(name, reference)
            if self._update_shut_state(name, target, reference):
                continue          # quyu bağlıdır — BHP-ni dartmağa ehtiyac yox
            if not math.isfinite(target):
                continue                       # cari BHP saxlanılır
            current = self.bhp[name]
            # Nodal həll MÜMKÜN OLMAYANDA (lay təzyiqi və ya J yoxdur) hədəf
            # köhnə üsulla — son debitlərlə — hesablanır və həmin qiymət
            # axmayan quyuda BHP-ni səhvən YUXARI ata bilər (ölçülmüş
            # bistabillik). Belə halda qoruyucu qayda qüvvədə qalır.
            nodal = (reference is not None
                     and self._productivity.get(name) is not None)
            if not flowing and target > current and not nodal:
                continue
            proposal = current + self.relaxation * (target - current)
            change = max(-self.max_change, min(self.max_change, proposal - current))
            self._apply(name, current + change)
            largest = max(largest, abs(change))
        return largest

    def _update_productivity(self, name: str, stream: WellStream, reference) -> None:
        """Quyunun məhsuldarlıq əmsalı J = q_maye / (p_lay − BHP), m³/gün/bar.

        Bu, IPR düz xəttinin mailliyidir və HƏR ADDIMDA yenilənir — yəni
        lay tükəndikcə əyri də yenilənir. Drawdown sıfıra yaxındırsa
        (ədədi cəhətdən etibarsız) əvvəlki qiymət saxlanılır.
        """
        if reference is None:
            return
        drawdown = reference - self.bhp.get(name, reference)
        if drawdown <= MIN_DRAWDOWN_BAR:
            return
        liquid = stream.oil + stream.water
        if liquid <= 0.0:
            return
        self._productivity[name] = liquid / drawdown

    def _operating_point(self, name: str, reference):
        """İŞ NÖQTƏSİ — IPR ∩ VLP (bax `update`-in sənədləşməsi).

        Qaytarır: BHP (bar), və ya `nan` — quyu bu THP ilə heç bir debitdə
        səthə axa bilmirsə (məhz BU, "bağlanmalıdır" siqnalıdır).

        Məlumat çatmayanda (lay təzyiqi və ya J hələ yoxdur — ilk addım)
        köhnə üsula qayıdır: son axan tərkiblə sadə tərs traverse.
        """
        ratios = self._last_stream.get(name)
        productivity = self._productivity.get(name)
        if reference is None or productivity is None or ratios is None:
            fallback = self._solve(
                name, ratios or self._stream(STATIC_PROBE_RATE, 0.0, 0.0))
            if reference is not None and math.isfinite(fallback):
                fallback = min(fallback, reference - STARTUP_DRAWDOWN_BAR)
            return fallback

        well = self._wells[name]
        thp_target = float(well.control.target)
        depth, tubing = self._depths[name], well.tubing
        liquid_reference = max(ratios.oil + ratios.water, 1e-12)

        def achieved_thp(bhp: float) -> float:
            """Bu BHP-də quyu başına ÇATAN təzyiq (IPR debiti ilə)."""
            liquid = productivity * max(reference - bhp, 0.0)
            scale = liquid / liquid_reference
            stream = self._stream(ratios.oil * scale, ratios.water * scale,
                                  ratios.gas * scale)
            if not stream.is_flowing():
                return float("-inf")
            return pressure_traverse(bhp, depth, tubing, stream, pvt=self.pvt,
                                     fluids=self.fluids, holdup=self.holdup).thp

        low, high = thp_target, max(reference - 1e-3, thp_target)
        best = achieved_thp(high)
        if not math.isfinite(best) or best < thp_target:
            return float("nan")          # ən əlverişli nöqtədə belə axmır
        for _ in range(MAX_NODAL_ITERATIONS):
            middle = 0.5 * (low + high)
            value = achieved_thp(middle)
            if not math.isfinite(value) or value < thp_target:
                low = middle
            else:
                high = middle
            if high - low < NODAL_TOLERANCE_BAR:
                break
        return high

    def _reservoir_pressure(self, name: str, pressure):
        """Quyunun perforasiyalarındakı ƏN YÜKSƏK lay təzyiqi, bar."""
        if pressure is None:
            return None
        cells = [c.cell for c in self._connections[name]]
        if not cells:
            return None
        return float(max(float(pressure[c]) for c in cells))

    def _update_shut_state(self, name: str, target: float, reference) -> bool:
        """AVTOMATİK BAĞLANMA — quyu axa bilmirsə tənlikdən çıxarılır.

        NİYƏ LAZIMDIR (ölçüldü, Seans 25): tələb olunan BHP lay təzyiqini
        keçəndə quyunu sadəcə lay təzyiqində SAXLAMAQ sabit görünsə də,
        ədədi cəhətdən ən pis vəziyyətdir — drawdown sıfıra yaxındır, quyu
        həddi ilə "sürünür" və Nyuton Δt-ni 0.003 günə qədər kəsirdi.

        Kommersiya simulyatorlarının davranışı: belə quyu BAĞLANIR.
        Burada bağlamaq = bağlantıların quyu indeksini SIFIRLAMAQ, yəni
        hasilat DƏQİQ sıfır olur və quyu həddi Jakobiandan tamamilə
        çıxır (qalıq/Jakobian kodu TOXUNULMUR — onlar sadəcə WI = 0 görür).

        Yenidən açılma: tələb olunan BHP lay təzyiqindən `REOPEN_MARGIN_BAR`
        qədər aşağı düşəndə (məsələn vurucu quyu lay təzyiqini qaldıranda).
        Ehtiyat aç/bağla rəqsinin qarşısını alır.

        Qaytarır: quyu BAĞLI qaldısa `True` (BHP yenilənmir).
        """
        if reference is None or name not in self._last_stream:
            return False       # məlumat azdır — bağlamaq üçün əsas yoxdur
        can_flow = math.isfinite(target)
        if self.shut.get(name):
            if can_flow and target < reference - REOPEN_MARGIN_BAR:
                self._set_shut(name, False)
                LOG.info("THP quyusu '%s': şərait düzəldi — quyu yenidən "
                         "AÇILDI (iş nöqtəsi %.1f bar, lay təzyiqi %.1f bar).",
                         name, target, reference)
                return False
            return True
        if not can_flow:
            self._set_shut(name, True)
            self.bhp[name] = reference
            LOG.warning("THP quyusu '%s': HEÇ BİR debitdə THP = %.1f bar-a "
                        "çatmır (lay təzyiqi %.1f bar) — AVTOMATİK BAĞLANDI "
                        "(süni qaldırma lazımdır).", name,
                        float(self._wells[name].control.target), reference)
            return True
        return False

    def _set_shut(self, name: str, shut: bool) -> None:
        """Bağlantıların quyu indeksini sıfırlayır / bərpa edir."""
        self.shut[name] = shut
        for connection, original in zip(self._connections[name],
                                        self._well_index[name]):
            connection.well_index = 0.0 if shut else original

    def record(self, result) -> None:
        """Bu addımda İŞLƏDİLƏN BHP-ni nəticəyə yazır (`update`-dən ƏVVƏL)."""
        for name, value in self.bhp.items():
            result.well_bhp.setdefault(name, []).append(float(value))
