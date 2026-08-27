"""Quyu naməlumları — OPM tipli "standart quyu modeli", MƏRHƏLƏ 1.

NİYƏ LAZIMDIR

Hazırkı quyu modelimizdə BHP **sabit sərhəd şərtidir**: debiti ondan
hesablayırıq. Bu, sadədir, lakin quyu öz hədəfinə yaxınlaşanda güclü
qeyri-xəttilik yaradır — ölçüldü ki, Nyuton istiqaməti boyunca
residual yalnız ~17 % azala bilir və həll yığılmır (bax `A7_PLAN.md`,
v59-v60).

OPM Flow (və bütün sənaye simulyatorları) fərqli edir: quyunun BHP-si
**naməlum dəyişəndir** (primary variable) və ona ayrıca **idarəetmə
tənliyi** yazılır:

    BHP idarəsində:   R_ctrl = p_bhp − p_hədəf = 0
    RATE idarəsində:  R_ctrl = Σ q_perf − q_hədəf = 0

Beləliklə idarəetmə rejimi dəyişəndə DEBİT sıçramır — sadəcə HANSI
tənliyin işlədiyi dəyişir. Bu, kəsilməzliyi qoruyur.

BU MƏRHƏLƏDƏ NƏ VAR

Yalnız VƏZİYYƏTİN ÖZÜ: quyu naməlumları, onların Nyuton vektorunda
necə yerləşdiyi və oradan geri oxunması. Tənliklər (mərhələ 2),
Jakobian (mərhələ 3) və Nyuton inteqrasiyası (mərhələ 4) sonrakı
mərhələlərdədir. Bu modul heç bir mövcud kodu DƏYİŞMİR — yanında
yaşayır.

VEKTOR YERLƏŞMƏSİ

    [ p₀, Sw₀, x₀,  p₁, Sw₁, x₁,  …,  p_{N−1}, Sw_{N−1}, x_{N−1},
      bhp₀, bhp₁, …, bhp_{W−1} ]
      └─────────── rezervuar (3N) ───────────┘ └─── quyular (W) ───┘

Quyu naməlumları SONA əlavə olunur, rezervuar blokunun 3×3 strukturu
TOXUNULMUR. Bu, qəsdən belədir: CPR ön-şərtçisi (A6) rezervuar
blokunun müntəzəm strukturuna əsaslanır, quyuları sona qoymaq onu
pozmur. OPM də eyni yerləşməni işlədir (quyu blokunu Schur tamamlayıcı
ilə ayırır).
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np

from .state import VARIABLES_PER_CELL, ReservoirState


@dataclass
class WellUnknowns:
    """Hər quyu üçün bir naməlum: quyudibi təzyiq (BHP), bar.

    `names` sırası vektor sırasını təyin edir və bütün mərhələlərdə
    SABİT qalmalıdır — indekslər buna görə hesablanır.
    """
    names: List[str]
    bhp: np.ndarray

    def __post_init__(self):
        self.bhp = np.asarray(self.bhp, dtype=float).ravel()
        if len(self.names) != self.bhp.size:
            raise ValueError(
                f"quyu sayı uyğunsuzdur: {len(self.names)} ad, "
                f"{self.bhp.size} BHP dəyəri")

    @property
    def count(self) -> int:
        return len(self.names)

    def index_of(self, name: str) -> int:
        """Quyunun vektordakı YERLİ nömrəsi (0-dan, quyular arasında)."""
        return self.names.index(name)

    def bhp_of(self, name: str) -> float:
        return float(self.bhp[self.index_of(name)])

    def copy(self) -> "WellUnknowns":
        return WellUnknowns(list(self.names), self.bhp.copy())

    @classmethod
    def from_connections(cls, connections: Sequence,
                         reservoir_pressure: np.ndarray,
                         initial_bhp: Optional[dict] = None) -> "WellUnknowns":
        """Quyu bağlantılarından ilkin BHP dəyərləri qurur.

        İLKİN QİYMƏT SEÇİMİ: hədəf BHP verilibsə o, əks halda quyunun
        ilk perforasiyasındakı hüceyrə təzyiqi. Bu, Nyutonun ilk
        addımını qısaldır — quyu təzyiqi onsuz da lay təzyiqinə yaxın
        olur. RATE idarəli quyularda hədəf BHP olmadığı üçün hüceyrə
        təzyiqi yeganə ağıllı başlanğıcdır.
        """
        from ...domain.wells import ControlMode

        names: List[str] = []
        values: List[float] = []
        seen = set()
        for connection in connections:
            if connection.well_name in seen:
                continue
            seen.add(connection.well_name)
            names.append(connection.well_name)

            if initial_bhp and connection.well_name in initial_bhp:
                values.append(float(initial_bhp[connection.well_name]))
            elif connection.mode is ControlMode.BHP:
                values.append(float(connection.target))
            else:
                values.append(float(reservoir_pressure[connection.cell]))
        return cls(names, np.array(values, dtype=float))


@dataclass
class CoupledState:
    """Rezervuar + quyu naməlumları BİRLİKDƏ.

    Nyuton bu birləşmiş vektor üzərində işləyəcək (mərhələ 4).
    Rezervuar hissəsi mövcud `ReservoirState`-dir — heç nə
    dəyişmir, yalnız yanına quyu naməlumları əlavə olunur.
    """
    reservoir: ReservoirState
    wells: WellUnknowns

    @property
    def size(self) -> int:
        return self.reservoir.ncell * VARIABLES_PER_CELL + self.wells.count

    @property
    def well_offset(self) -> int:
        """Quyu naməlumlarının vektorda başladığı indeks."""
        return self.reservoir.ncell * VARIABLES_PER_CELL

    def well_index(self, name: str) -> int:
        """Quyunun QLOBAL vektor indeksi."""
        return self.well_offset + self.wells.index_of(name)

    @property
    def water_saturation(self) -> np.ndarray:
        """`AdaptiveTimeStepper._saturation_change()` üçün — o, ümumi
        (generic) dizayn edilib və `state.water_saturation`-ı birbaşa
        oxuyur (bax `time_stepping.py`). Bu, ReservoirState-in özünə
        BƏRABƏR davranışdır — CoupledState-i ondan fərqləndirməmək
        üçün əlavə olunub."""
        return self.reservoir.water_saturation

    def copy(self) -> "CoupledState":
        return CoupledState(self.reservoir.copy(), self.wells.copy())

    # ─────────────────────────────────────────── vektor çevirmələri
    def to_vector(self) -> np.ndarray:
        vector = np.empty(self.size)
        vector[:self.well_offset] = self.reservoir.to_vector()
        vector[self.well_offset:] = self.wells.bhp
        return vector

    @classmethod
    def from_vector(cls, vector: np.ndarray,
                    well_names: List[str]) -> "CoupledState":
        vector = np.asarray(vector, float)
        well_count = len(well_names)
        offset = vector.size - well_count
        if offset % VARIABLES_PER_CELL != 0:
            raise ValueError(
                "vektorun uzunluğu uyğun deyil: rezervuar hissəsi "
                f"({offset}) 2-ə bölünmür")
        reservoir = ReservoirState.from_vector(vector[:offset])
        wells = WellUnknowns(list(well_names), vector[offset:].copy())
        return cls(reservoir, wells)

    def updated(self, delta: np.ndarray, sw_min: float, sw_max: float,
               max_pressure_change: Optional[float] = None,
               max_saturation_change: Optional[float] = None,
               max_bhp_change: Optional[float] = None) -> "CoupledState":
        """Nyuton addımını tətbiq edir — rezervuar VƏ quyular.

        Rezervuar hissəsi mövcud `ReservoirState.updated()`-ə ötürülür
        (Appleyard kəsməsi orada, dəyişməyib). Quyu BHP-si üçün ayrıca
        hədd var: quyu təzyiqi bir iterasiyada həddindən çox sıçrasa,
        perforasiya debitləri qeyri-real dəyərlərə gedər.
        """
        reservoir_delta = delta[:self.well_offset]
        well_delta = delta[self.well_offset:]

        reservoir = self.reservoir.updated(
            reservoir_delta, sw_min, sw_max,
            max_pressure_change=max_pressure_change,
            max_saturation_change=max_saturation_change)

        if max_bhp_change:
            well_delta = np.clip(well_delta, -max_bhp_change, max_bhp_change)
        wells = WellUnknowns(list(self.wells.names),
                             self.wells.bhp + well_delta)
        return CoupledState(reservoir, wells)
