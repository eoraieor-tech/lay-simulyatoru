"""Günlük göstəricilər — Seans 45. **Qt-dən ASILI DEYİL.**

NİYƏ LAZIMDIR. Mühərrik günbəgün yox, ADAPTİV addımlarla gedir
(nümunə qaçışı: 1500 gün 88 addımda, orta Δt 17 gün). `TimeSeries`-də
yalnız addımların sonu var — "250-ci gün" adlı sətir yoxdur. Sahibkar
hər günün göstəricilərini (yataq və hər quyu üzrə) görmək istədi.

GÜNLÜK DƏYƏR HARADAN GƏLİR (mühərrik TOXUNULMUR):

    Debitlər     — implicit Eyler addım boyu debiti SABİT götürür, yəni
                   debit addım daxilində pilləli funksiyadır. Günün
                   debiti = həmin gündə çıxarılan HƏCM / günün uzunluğu.
                   Gün addımın içindədirsə, bu, addımın öz debitidir;
                   addım sərhədi günün ortasına düşürsə, iki addımın
                   zamanla çəkili ortasıdır. Mühərrikin hesabladığından
                   KƏNARA ÇIXMIR.
    Kumulyativ,  — debit sabit olduğu üçün addım boyu XƏTTİ artır, yəni
    RF             xətti interpolyasiya DƏQİQDİR. Günlük həcmlərin cəmi
                   son kumulyativə bərabərdir (testlə yoxlanır).
    Su kəsri,    — günün həcmlərindən hesablanır (debitlərin nisbəti).
    GOR
    Orta təzyiq  — yalnız addımın SONUNDA məlumdur; aralıq günlər qonşu
                   iki dəyər arasında xətti interpolyasiyadır → TƏXMİNİ.
                   Sütun adında bu açıq yazılır. İlk addımın sonundan
                   əvvəlki günlər ilk dəyəri alır (t = 0 dəyəri seriyada yoxdur).
    BHP / THP    — addım boyu sabit tətbiq olunur; günün SONUNU örtən
                   addımın dəyəri götürülür. THP-nin `nan`-ı (quyu axmır)
                   olduğu kimi saxlanılır — sıfıra çevrilmir.

Günlər: 1, 2, …, ⌊T⌋ və T tam deyilsə sonda T-nin özü (qısa son gün).
T — faktiki çatılan son an (`series.time[-1]`), qaçış dayandırılıbsa da.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..simulation.results import SimulationResult

#: Tam gün kimi qəbul olunan sapma — `1499.9999999` 1500-cü gündür.
_DAY_TOLERANCE = 1e-9

#: Adlar — başlıqda vahid MÜTLƏQDİR (bax `UNITS.md`), CSV ixracı ilə eyni üslub.
DAY = "Gün"
OIL = "q_neft [m³/gün]"
WATER = "q_su [m³/gün]"
GAS = "q_qaz [m³/gün]"
WATER_INJ = "q_vurulan_su [m³/gün]"
GAS_INJ = "q_vurulan_qaz [m³/gün]"
CUM_OIL = "kum_neft [m³]"
CUM_WATER = "kum_su [m³]"
CUM_GAS = "kum_qaz [m³]"
CUM_WATER_INJ = "kum_vurulan_su [m³]"
CUM_GAS_INJ = "kum_vurulan_qaz [m³]"
WATER_CUT = "su_kəsri [%]"
GOR = "GOR [sm³/sm³]"
RF = "RF [%]"
PRESSURE = "orta_P [bar] (interp.)"
BHP = "BHP [bar]"
THP = "THP [bar]"

#: Yataq sətri üçün seçici adı (UI və CSV-də quyu adları ilə yanaşı).
FIELD = "Yataq"


@dataclass
class DailyTable:
    """Günlük cədvəl — sütun adı → dəyərlər (hamısı `days` uzunluğunda).

    Sütun sırası lüğətin sırasıdır. Nəticədə OLMAYAN kəmiyyət (məs. iki
    fazalı qaçışda qaz) sütun YARATMIR — `results_export`-dakı qayda.
    """
    days: np.ndarray = field(default_factory=lambda: np.zeros(0))
    field_columns: Dict[str, np.ndarray] = field(default_factory=dict)
    wells: Dict[str, Dict[str, np.ndarray]] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.days)

    def targets(self) -> List[str]:
        """Seçilə bilən obyektlər: əvvəl yataq, sonra quyular (əlifba ilə)."""
        return [FIELD] + list(self.wells)

    def columns(self, target: str = FIELD) -> Dict[str, np.ndarray]:
        """`FIELD` üçün yataq sütunları, quyu adı üçün həmin quyununku."""
        return self.field_columns if target == FIELD else self.wells[target]

    def row(self, index: int, target: str = FIELD) -> Dict[str, float]:
        """Bir günün bütün göstəriciləri — `{DAY: gün, sütun: dəyər, …}`."""
        values = {DAY: float(self.days[index])}
        values.update({name: float(column[index])
                       for name, column in self.columns(target).items()})
        return values

    def index_of_day(self, day: float) -> int:
        """`day`-i örtən sətir (gün `(d-1, d]` aralığıdır)."""
        if len(self.days) == 0:
            raise IndexError("Günlük cədvəl boşdur.")
        index = int(np.searchsorted(self.days, day - _DAY_TOLERANCE))
        return min(max(index, 0), len(self.days) - 1)


# ═══════════════════════════════ zaman ═══════════════════════════════

def day_edges(end_time: float) -> np.ndarray:
    """Günlərin SONLARI: 1, 2, …, ⌊T⌋ (+ T, tam deyilsə)."""
    if end_time <= _DAY_TOLERANCE:
        return np.zeros(0)
    whole = int(np.floor(end_time + _DAY_TOLERANCE))
    edges = np.arange(1, whole + 1, dtype=float)
    if end_time - whole > _DAY_TOLERANCE:
        edges = np.append(edges, float(end_time))
    return edges


def _cumulative_at(times: np.ndarray, rates: np.ndarray,
                   edges: np.ndarray) -> np.ndarray:
    """Pilləli debitin inteqralı günlərin sonunda — DƏQİQ (debit addım
    boyu sabitdir, inteqral isə xəttidir)."""
    steps = np.diff(np.concatenate(([0.0], times)))
    cumulative = np.concatenate(([0.0], np.cumsum(rates * steps)))
    return np.interp(edges, np.concatenate(([0.0], times)), cumulative)


def _covering_step(times: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Hər günün SONUNU örtən addımın indeksi: `t[i-1] < son ≤ t[i]`."""
    index = np.searchsorted(times, edges - _DAY_TOLERANCE, side="left")
    return np.clip(index, 0, len(times) - 1)


def _series(values: Optional[Sequence], length: int) -> Optional[np.ndarray]:
    """Uzunluğu addım sayına UYĞUN olan sıra, yoxsa `None`.

    Boş və ya uyğunsuz sıra (məs. RATE quyusunun boş `well_bhp`-si)
    uydurma dəyərlə doldurulmur — sadəcə sütun yaranmır.
    """
    if not values or len(values) != length:
        return None
    return np.asarray(values, dtype=float)


def _ratio(numerator: np.ndarray, denominator: np.ndarray,
           scale: float = 1.0) -> np.ndarray:
    """Mühərriklərin öz qaydası: məxrəc `1e-12`-dən kiçik olmur (0/0 → 0)."""
    return numerator / np.maximum(denominator, 1e-12) * scale


# ═══════════════════════════════ qurucu ══════════════════════════════

class _Builder:
    """Bir nəticənin günlük cədvəlini sütun-sütun qurur."""

    def __init__(self, result: SimulationResult):
        self.result = result
        self.times = np.asarray(result.series.time, dtype=float)
        self.count = len(self.times)
        self.edges = day_edges(float(self.times[-1])) if self.count else np.zeros(0)
        self.lengths = np.diff(np.concatenate(([0.0], self.edges)))
        self.cover = (_covering_step(self.times, self.edges)
                      if self.count else np.zeros(0, dtype=int))

    # -------------------------------------------------------- köməkçilər
    def _volume_pair(self, rates: Optional[np.ndarray]
                     ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """(günlük debit, günün sonundakı kumulyativ) — sıra yoxdursa `None`."""
        if rates is None:
            return None
        cumulative = _cumulative_at(self.times, rates, self.edges)
        volume = np.diff(np.concatenate(([0.0], cumulative)))
        return volume / self.lengths, cumulative

    def _rates_block(self, source: Dict[str, Optional[np.ndarray]]
                     ) -> Dict[str, np.ndarray]:
        """Debit + kumulyativ + su kəsri + GOR — yataq və quyu üçün EYNİ qayda."""
        columns: Dict[str, np.ndarray] = {}
        pairs = {name: self._volume_pair(source.get(name))
                 for name in (OIL, WATER, GAS, WATER_INJ, GAS_INJ)}
        cumulative_name = {OIL: CUM_OIL, WATER: CUM_WATER, GAS: CUM_GAS,
                           WATER_INJ: CUM_WATER_INJ, GAS_INJ: CUM_GAS_INJ}

        for name in (OIL, WATER, GAS, WATER_INJ, GAS_INJ):
            if pairs[name] is not None:
                columns[name] = pairs[name][0]
        for name in (OIL, WATER, GAS, WATER_INJ, GAS_INJ):
            if pairs[name] is not None:
                columns[cumulative_name[name]] = pairs[name][1]

        oil, water, gas = pairs[OIL], pairs[WATER], pairs[GAS]
        if oil is not None and water is not None:
            columns[WATER_CUT] = _ratio(water[0], oil[0] + water[0], 100.0)
        if oil is not None and gas is not None:
            columns[GOR] = _ratio(gas[0], oil[0])
        return columns

    def _at_day_end(self, values: Optional[np.ndarray]) -> Optional[np.ndarray]:
        return None if values is None else values[self.cover]

    # -------------------------------------------------------------- yataq
    def field_columns(self) -> Dict[str, np.ndarray]:
        series = self.result.series
        get = lambda name: _series(getattr(series, name, None), self.count)
        columns = self._rates_block({
            OIL: get("oil_rate"), WATER: get("water_rate"),
            GAS: get("gas_rate"), WATER_INJ: get("water_injection_rate"),
            GAS_INJ: get("gas_injection_rate")})

        recovery = get("recovery_factor")
        if recovery is not None:
            # RF kumulyativ neftlə mütənasibdir → addım boyu xətti, t=0-da 0
            columns[RF] = np.interp(self.edges,
                                    np.concatenate(([0.0], self.times)),
                                    np.concatenate(([0.0], recovery)))
        pressure = get("average_pressure")
        if pressure is not None:
            columns[PRESSURE] = np.interp(self.edges, self.times, pressure)
        return columns

    # ------------------------------------------------------------- quyular
    def well_names(self) -> List[str]:
        names = set()
        for attribute in ("well_oil_rate", "well_water_rate", "well_gas_rate",
                          "well_water_injection_rate",
                          "well_gas_injection_rate", "well_bhp", "well_thp"):
            for name, values in (getattr(self.result, attribute, {}) or {}).items():
                if _series(values, self.count) is not None:
                    names.add(name)
        return sorted(names)

    def well_columns(self, name: str) -> Dict[str, np.ndarray]:
        get = lambda attribute: _series(
            (getattr(self.result, attribute, {}) or {}).get(name), self.count)
        columns = self._rates_block({
            OIL: get("well_oil_rate"), WATER: get("well_water_rate"),
            GAS: get("well_gas_rate"),
            WATER_INJ: get("well_water_injection_rate"),
            GAS_INJ: get("well_gas_injection_rate")})
        for header, attribute in ((BHP, "well_bhp"), (THP, "well_thp")):
            values = self._at_day_end(get(attribute))
            if values is not None:
                columns[header] = values
        return columns


def daily_table(result: SimulationResult) -> DailyTable:
    """Nəticədən günlük cədvəl. Nəticə boşdursa, boş cədvəl."""
    if not result.series.time:
        return DailyTable()
    builder = _Builder(result)
    table = DailyTable(days=builder.edges.copy(),
                       field_columns=builder.field_columns())
    for name in builder.well_names():
        columns = builder.well_columns(name)
        if columns:
            table.wells[name] = columns
    return table
