"""Quyu məhdudiyyətləri — Nyuton həlləri ARASINDA tətbiq olunan qaydalar (B7).

Bu moduldakı hər şey `ThpController` (Q-15, Q-17) ilə eyni fəlsəfədədir:
QALIQ və JAKOBİAN yeni rejim tanımır. Qayda bağlantı obyektinin sahəsini
addımlar arasında yerində dəyişir, qalıq isə sadəcə yeni ədədi görür.

RATE HƏDƏFİNİN PERFORASİYALARA BÖLÜNMƏSİ (Seans 27)

ƏVVƏLKİ SƏHV (ölçüldü): qalıq, Jakobian və IMPES RATE hədəfini HƏR
perforasiyaya TAM yazırdı — `total = -abs(connection.target)` hər bağlantı
üçün ayrıca. Nəticədə 3 perforasiyalı quyu hədəfin 3 QATINI hasil edirdi
(50 m³/gün hədəf → 150 m³/gün). Heç bir xəbərdarlıq yox idi.

İNDİ hədəf quyunun bağlantıları arasında Peaceman nisbətində bölünür:

    pay_c = WI_c · λ_c / Σ_k WI_k · λ_k,        Σ pay = 1

Bu, BHP rejimində debitin təbii paylanmasıdır (`q_c = WI_c·λ_c·Δp`,
bütün perforasiyalarda eyni BHP).

NİYƏ λ ƏVVƏLKİ ADDIMDANDIR. Pay cari iterasiyanın λ-sı ilə hesablansaydı,
bir perforasiyanın debiti quyunun BÜTÜN digər hüceyrələrinin doymuşluğundan
asılı olardı — Jakobianda hüceyrələr arası yeni törəmələr (yeni seyrəklik
strukturu) yaranardı. Pay addımın ƏVVƏLİNDƏ (yığılmış vəziyyətdən) verilir və
Nyuton daxilində SABİTDİR, ona görə Jakobianda RATE törəmələri sadəcə paya
vurulur — struktur dəyişmir. Quyunun CƏMİ debiti dəqiqdir; yalnız təbəqələr
arası paylanma bir addım gecikir.

GERİYƏ UYĞUNLUQ: tək perforasiyalı quyuda pay DƏQİQ 1.0-dır (`x/x`), yəni
`target * 1.0` — nəticə bit-bit eyni.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, List, Optional, Sequence

from ..domain.wells import ControlMode, Phase, RateBasis
from ..logging_setup import get_logger

LOG = get_logger(__name__)

#: RATE rejiminə QAYITMAQ üçün ehtiyat, bar. Sərhəddə qayıtmaq quyunu dərhal
#: yenidən limitə salardı — rejim rəqsi. THP-dəki `REOPEN_MARGIN_BAR` ilə
#: eyni rol və eyni qiymət.
RATE_RESTORE_MARGIN_BAR = 2.0


def _group_by_well(connections: Sequence) -> Dict[str, List[int]]:
    groups: Dict[str, List[int]] = {}
    for position, connection in enumerate(connections):
        groups.setdefault(connection.well_name, []).append(position)
    return groups


def needs_rate_allocation(connections: Sequence) -> bool:
    """Bir neçə perforasiyalı RATE quyusu varmı — yoxdursa mühərrik payı
    heç vaxt yeniləmir (əlavə hesab xərci yoxdur)."""
    connections = list(connections)
    return any(len(positions) > 1
               and connections[positions[0]].mode is ControlMode.RATE
               for positions in _group_by_well(connections).values())


def assign_rate_shares(connections: Sequence,
                       mobility: Optional[Sequence[float]] = None) -> None:
    """Hər bağlantının `rate_share`-ini yerində yazır.

    `mobility` — bağlantılarla EYNİ sırada LAY HƏCMİ mobilliyi (istismarçı
    `λw + λo`, vurucu vurulan fazanın son nöqtə mobilliyi). `None` olanda
    pay yalnız WI ilə verilir (mühərrik ilk addımda λ ilə yeniləyənə qədər).

    Bütün çəkilər sıfırdırsa (məs. hər perforasiyada λ = 0), WI nisbətinə,
    o da sıfırdırsa bərabər paya qayıdılır — cəm HƏMİŞƏ 1 qalır.
    """
    connections = list(connections)
    for positions in _group_by_well(connections).values():
        if len(positions) == 1:
            connections[positions[0]].rate_share = 1.0
            continue
        weights = [_weight(connections[p],
                           None if mobility is None else mobility[p])
                   for p in positions]
        total = sum(weights)
        if not total > 0.0:
            weights = [_weight(connections[p], None) for p in positions]
            total = sum(weights)
        for position, weight in zip(positions, weights):
            connections[position].rate_share = (
                weight / total if total > 0.0 else 1.0 / len(positions))


def _weight(connection, mobility: Optional[float]) -> float:
    value = float(connection.well_index)
    if mobility is not None:
        value *= float(mobility)
    return value if math.isfinite(value) and value > 0.0 else 0.0


class BhpLimitController:
    """RATE quyusunun BHP həddi — rejim keçidi addımlar ARASINDA (B7 addım 2).

    SPE1-in istismarçısı debitlə idarə olunur, lakin BHP minimal həddən aşağı
    düşə bilməz; lay tükəndikcə quyu debiti saxlaya bilmir və BHP idarəsinə
    keçir. Qalıq və Jakobian bu qaydanı TANIMIR — nəzarətçi bağlantının
    `mode`/`target`-ini yerində dəyişir (Q-15-dəki THP yanaşması).

    HƏDƏF DEBİTİ VERƏN BHP — Peaceman-ın tərsi, lay həcmində:

        istismarçı:  Σ WI·λ·(p − BHP) = q   →   BHP = (Σ WI·λ·p − q) / Σ WI·λ
        vurucu:      Σ WI·λ·(BHP − p) = q   →   BHP = (Σ WI·λ·p + q) / Σ WI·λ

    λ — `well_rates`-in BHP budağındakı EYNİ mobillik (`connection_mobilities`),
    ona görə bu BHP-yə keçəndə quyunun cəmi debiti sıçramır.

    QAYDALAR (hamısı yığılmış addımdan SONRA yoxlanılır, Nyuton daxilində yox):

    * RATE quyusu: tələb olunan BHP həddi pozursa → `mode = BHP`,
      `target = hədd`; mobillik sıfırdırsa (quyu axa bilmir) da keçir;
    * limitdə olan quyu: hədəf debit həddən `RATE_RESTORE_MARGIN_BAR` qədər
      ƏLVERİŞLİ BHP ilə əldə olunursa → RATE-ə qayıdır (histerezis);
    * hər quyu bir addımda ən çox BİR dəfə rejim dəyişir (`begin_step`) —
      addımın təkrarlanan həlli arasında aç-qapa mümkün deyil.

    Bilinən sadələşdirmə: BHP rejimində perforasiya debiti `min(q, 0)` ilə
    kəsilir (çarpaz axın yoxdur), tərs düstur isə kəsməni nəzərə almır —
    təbəqələr arasında təzyiq fərqi böyük olanda qiymət təqribidir.
    """

    def __init__(self, connections: Sequence, mobility: Callable,
                 rate_target: Optional[Callable] = None):
        """`mobility(state)` — bağlantılarla EYNİ sırada lay həcmi mobillikləri.

        `rate_target(ad)` — quyunun CARİ lay həcmi RATE hədəfi (`None` →
        qurulma anındakı `connection.target`). SƏTH bazalı quyuda hədəf
        addım-addım dəyişir (`SurfaceRateController.rate_target`); onsuz
        RATE-ə qayıdan quyu səth rəqəmini lay hədəfi kimi bərpa edərdi.
        """
        self._connections = list(connections)
        self._mobility = mobility
        self._rate_lookup = rate_target
        self._positions: Dict[str, List[int]] = {}
        for position, connection in enumerate(self._connections):
            if (connection.mode is ControlMode.RATE
                    and getattr(connection, "bhp_limit", None) is not None):
                self._positions.setdefault(connection.well_name, []).append(position)
        first = {name: self._connections[positions[0]]
                 for name, positions in self._positions.items()}
        self._rate_target = {name: abs(float(c.target)) for name, c in first.items()}
        self._limit = {name: float(c.bhp_limit) for name, c in first.items()}
        self._injector = {name: bool(c.is_injector) for name, c in first.items()}

        #: quyu adı → hazırda BHP limitindədirmi
        self.limited: Dict[str, bool] = {name: False for name in self._positions}
        #: son qiymətləndirmədə quyunun BHP-si, bar (RATE-də tələb olunan, limitdə hədd)
        self.bhp: Dict[str, float] = {}
        #: son qiymətləndirilən həllin HANSI rejimlə alındığı ("RATE" / "BHP")
        self.mode_used: Dict[str, str] = {}
        #: qaçış boyu rejim keçidlərinin sayı (rəqs diaqnostikası)
        self.switches: Dict[str, int] = {name: 0 for name in self._positions}
        self._switched_this_step: set = set()

    @property
    def active(self) -> bool:
        """BHP limitli RATE quyusu yoxdursa nəzarətçi heç nə etmir."""
        return bool(self._positions)

    @property
    def well_count(self) -> int:
        return len(self._positions)

    def begin_step(self) -> None:
        """Yeni zaman addımı — hər quyuya yenidən bir keçid hüququ verilir."""
        self._switched_this_step = set()

    def _rate(self, name: str) -> float:
        """Quyunun cari lay həcmi RATE hədəfi."""
        if self._rate_lookup is not None:
            value = self._rate_lookup(name)
            if value is not None:
                return float(value)
        return self._rate_target[name]

    def implied_bhp(self, name: str, pressure, mobility) -> float:
        """Hədəf debiti verən BHP, bar; mobillik sıfırdırsa `nan`."""
        weight = weighted = 0.0
        for position in self._positions[name]:
            connection = self._connections[position]
            value = float(connection.well_index) * float(mobility[position])
            weight += value
            weighted += value * float(pressure[connection.cell])
        if not weight > 0.0:
            return float("nan")
        rate = self._rate(name)
        if self._injector[name]:
            return (weighted + rate) / weight
        return (weighted - rate) / weight

    def update(self, state) -> List[str]:
        """Yığılmış həlli qiymətləndirir, lazım olsa rejimi dəyişir.

        Qaytarır: rejimi DƏYİŞƏN quyuların adları — boşdursa addım olduğu
        kimi qəbul olunur, əks halda mühərrik addımı yenidən həll etməlidir.
        """
        mobility = self._mobility(state)
        switched: List[str] = []
        for name in self._positions:
            limit = self._limit[name]
            implied = self.implied_bhp(name, state.pressure, mobility)
            limited = self.limited[name]
            self.bhp[name] = limit if limited else implied
            self.mode_used[name] = "BHP" if limited else "RATE"
            if name in self._switched_this_step:
                continue
            injector = self._injector[name]
            if limited:
                restore = math.isfinite(implied) and (
                    implied < limit - RATE_RESTORE_MARGIN_BAR if injector
                    else implied > limit + RATE_RESTORE_MARGIN_BAR)
                if restore:
                    self._set_limited(name, False, implied)
                    switched.append(name)
            else:
                violated = not math.isfinite(implied) or (
                    implied > limit if injector else implied < limit)
                if violated:
                    self._set_limited(name, True, implied)
                    switched.append(name)
        return switched

    def _set_limited(self, name: str, limited: bool, implied: float) -> None:
        self.limited[name] = limited
        self.switches[name] += 1
        self._switched_this_step.add(name)
        for position in self._positions[name]:
            connection = self._connections[position]
            connection.mode = ControlMode.BHP if limited else ControlMode.RATE
            connection.target = self._limit[name] if limited else self._rate(name)
        LOG.info("Quyu '%s': %s (hədəf debit üçün tələb olunan BHP %.1f bar, "
                 "limit %.1f bar).", name,
                 "BHP LİMİTİNƏ keçdi" if limited else "RATE rejiminə QAYITDI",
                 implied, self._limit[name])

    def record(self, result) -> None:
        """Bu addımın BHP-si və rejimi nəticəyə yazılır."""
        for name, value in self.bhp.items():
            result.well_bhp.setdefault(name, []).append(float(value))
            result.well_control_mode.setdefault(name, []).append(
                self.mode_used[name])


#: Səth debitinin hədəfdən icazə verilən nisbi sapması — bundan kiçik
#: sapmada addım təkrarlanmadan qəbul olunur.
SURFACE_RATE_TOLERANCE = 1e-3

#: Addım daxilində səth debiti düzəlişinin maksimal təkrar sayı.
MAX_SURFACE_ITERATIONS = 3


class SurfaceRateController:
    """SƏTH debiti hədəfi — addımlar arasında lay həcmi hədəfinə çevrilir (B7 addım 3).

    Qalıq RATE hədəfini LAY həcmində tanıyır (istismarçıda maye, vurucuda
    vurulan faza). SPE1 isə neftin SƏTH debitini (`ORAT`) və qazın SƏTH vurma
    debitini verir. Qalığa yeni hədd əlavə etmək əvəzinə (Q-15/Q-21 yanaşması)
    nəzarətçi `connection.target`-i addım-addım yenidən hesablayır.

    RATE budağının öz düsturları ilə, L lay həcmi hədəfində:

        istismarçı:  q_neft,səth = L · Σ_c pay_c · (1 − f_c) / Bo_c,   f = λw/(λw+λo)
        vurucu:      q_səth      = L · Σ_c pay_c / B_c

    * **Proqnoz** (`predict`, addımdan əvvəl): yığılmış vəziyyətin əmsalları ilə
      `L = q_səth / Σ pay·k` — B və f addım boyu az dəyişdiyi üçün adətən
      kifayətdir.
    * **Düzəliş** (`correct`, addımdan sonra): əldə olunan səth debiti hədəfdən
      `SURFACE_RATE_TOLERANCE`-dan çox fərqlənirsə `L ← L · hədəf / əldə olunan`
      və mühərrik addımı yenidən həll edir.

    BHP limitində olan (`mode = BHP`) quyunun hədəfinə toxunulmur.
    """

    def __init__(self, connections: Sequence, surface_factors: Callable):
        """`surface_factors(state)` — bağlantılarla EYNİ sırada `k_c` əmsalları."""
        self._connections = list(connections)
        self._factors = surface_factors
        self._positions: Dict[str, List[int]] = {}
        for position, connection in enumerate(self._connections):
            if (connection.mode is ControlMode.RATE
                    and getattr(connection, "rate_basis", RateBasis.RESERVOIR)
                    is RateBasis.SURFACE):
                self._positions.setdefault(connection.well_name, []).append(position)
        first = {name: self._connections[positions[0]]
                 for name, positions in self._positions.items()}
        #: istifadəçinin səth hədəfi (müsbət böyüklük)
        self.surface_target = {name: abs(float(c.target)) for name, c in first.items()}
        self._injector = {name: bool(c.is_injector) for name, c in first.items()}
        self._phase = {name: c.injected_phase for name, c in first.items()}
        #: qalığın gördüyü CARİ lay həcmi hədəfi
        self.reservoir_target: Dict[str, float] = {}

    @property
    def active(self) -> bool:
        return bool(self._positions)

    def rate_target(self, name: str) -> Optional[float]:
        """`BhpLimitController` üçün — səth bazalı olmayan quyuda `None`."""
        return self.reservoir_target.get(name)

    def predict(self, state) -> None:
        """Addımdan ƏVVƏL: yığılmış vəziyyətin əmsalları ilə lay hədəfi."""
        factors = self._factors(state)
        for name, positions in self._positions.items():
            total = sum(float(self._connections[p].rate_share) * float(factors[p])
                        for p in positions)
            if math.isfinite(total) and total > 0.0:
                self._set(name, self.surface_target[name] / total)

    def correct(self, rates) -> float:
        """Addımdan SONRA: əldə olunan səth debitinə görə düzəliş.

        Qaytarır: RATE rejimindəki quyular üzrə ƏN BÖYÜK nisbi sapma.
        """
        worst = 0.0
        for name, positions in self._positions.items():
            if self._connections[positions[0]].mode is not ControlMode.RATE:
                continue                          # BHP limitindədir
            achieved = self.achieved(name, rates)
            target = self.surface_target[name]
            if not (achieved > 0.0 and target > 0.0):
                continue
            error = abs(achieved - target) / target
            worst = max(worst, error)
            if error > SURFACE_RATE_TOLERANCE and name in self.reservoir_target:
                self._set(name, self.reservoir_target[name] * target / achieved)
        return worst

    def achieved(self, name: str, rates) -> float:
        """Həll olunmuş addımda quyunun səth debiti (müsbət böyüklük)."""
        if not self._injector[name]:
            return -float(rates.per_well_oil.get(name, 0.0))
        if self._phase[name] is Phase.GAS:
            return float((getattr(rates, "per_well_gas", None) or {}).get(name, 0.0))
        return float(rates.per_well_water.get(name, 0.0))

    def _set(self, name: str, value: float) -> None:
        self.reservoir_target[name] = float(value)
        for position in self._positions[name]:
            connection = self._connections[position]
            if connection.mode is ControlMode.RATE:
                connection.target = float(value)
