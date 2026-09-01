"""Vahid çevirmə qatı — mərkəzləşdirilmiş, SI-pivotlu.

**Niyə bura, `units.py`-a yox**: `domain/units.py`-dəki `UnitSystem`
(`METRIC`/`FIELD`) yalnız `darcy_constant`-ı və təsviri etiket sətirlərini
saxlayır — heç bir çevirmə funksiyası YOXDUR (audit təsdiqləyib: `FIELD`
heç yerdə istehlak edilmir). Bu fayl onu ƏVƏZ ETMİR (geriyə uyğunluq üçün
`units.py` toxunulmaz saxlanılıb), yalnız ÇATIŞMAYAN çevirmə qatını əlavə
edir.

**Mühərrikin daxili "kanonik" vahid sistemi SI DEYİL** — bar/m/mD/cP/
m³/gün/kg/m³ ("neft-sənayesi metrik", Eclipse-in METRIC vahid dəstinə
uyğun). Bu, BİLƏRƏKDƏN belədir: `discretization.py`-dəki Darsi sabiti
(`0.008527`), Nyuton qalığı (`residual.py`), quyu modeli
(`standard_well.py`) və PVT (`pvt.py`) bu sistemdə YAZILIB, sınanıb və
672 testlə qorunur. Bu modulları təmiz SI-ya (Pa, m³/s, Pa·s) köçürmək
bütöv Nyuton/Jakobian nüvəsini yenidən yazmaq deməkdir — "mövcud ədədi
davranışı SƏSSİZCƏ dəyişmə" qadağasına birbaşa ZİDD olardı. Ona görə bu
modul İKİ şey təklif edir:

    1. Həqiqi SI baza vahidlərinə (Pa, m, m², Pa·s, m³, m³/s, kg/m³, K)
       PİVOT edən, istənilən dəstəklənən vahid CÜTÜ arasında dəqiq
       çevirmə (`convert`, `to_si`, `from_si`) — elmi cəhətdən sərt,
       mühərrikdən müstəqil, ayrıca sınanan.
    2. Mühərrikin GERÇƏK işlədiyi vahidlərə (bar/m/mD/cP/m³/gün/kg/m³)
       çevirmək üçün rahatlıq funksiyaları (`to_engine_*`) — istənilən
       xarici vahiddə verilmiş məlumatı GİRİŞ sərhədində (CSV/GRDECL/UI)
       mühərrikin gözlədiyi ədədə çevirmək üçün.

Çevirmə YALNIZ GİRİŞ→DAXİLİ və DAXİLİ→ÇIXIŞ sərhədlərində aparılmalıdır
— mühərrik özü (residual.py, discretization.py, standard_well.py və s.)
bu modulu İSTEHLAK ETMİR və indiki mərhələdə DƏYİŞDİRİLMİR.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple

# ── fiziki sabitlər (dublikatların qarşısını almaq üçün mərkəzləşdirilib) ─
STANDARD_GRAVITY_M_S2 = 9.80665

# ── SI-pivot faktorları: `dəyər [vahid] * faktor = dəyər [SI baza vahidi]` ─
PRESSURE_TO_PA: Dict[str, float] = {
    "Pa": 1.0,
    "kPa": 1.0e3,
    "MPa": 1.0e6,
    "bar": 1.0e5,
    "psi": 6894.757293168361,
}

LENGTH_TO_M: Dict[str, float] = {
    "m": 1.0,
    "ft": 0.3048,
    "cm": 0.01,
}

AREA_TO_M2: Dict[str, float] = {
    "m2": 1.0,
    "ft2": 0.09290304,
    "acre": 4046.8564224,
}

#: Darsi tərifi: 1 D = 9.869232667160128e-13 m² (neft mühəndisliyinin
#: standart tərifi — 1 santipuazlıq maye 1 atm/sm təzyiq qradiyenti ilə
#: 1 sm²-lik en kəsikdən 1 sm/san sürətlə axanda keçiricilik 1 Darsi-dir).
PERMEABILITY_TO_M2: Dict[str, float] = {
    "m2": 1.0,
    "D": 9.869232667160128e-13,
    "mD": 9.869232667160128e-16,
}

VISCOSITY_TO_PAS: Dict[str, float] = {
    "Pa.s": 1.0,
    "cP": 1.0e-3,
}

#: 1 bbl = 42 ABŞ qallonu = 0.158987294928 m³ (dəqiq). `stb` FİZİKİ olaraq
#: EYNİ həcm vahididir (42 qallon) — fərq VAHİDDƏ deyil, ŞƏRTDƏDİR
#: (səth/stok-tank şəraiti ↔ `bbl`/`rb` adətən layın öz şəraitini bildirir).
#: Bu modul yalnız HƏCM ədədini çevirir; "hansı şəraitdə" sualı (Bo/Bw ilə
#: bağlıdır) mühərrikin PVT qatının işidir, bu modulun deyil.
VOLUME_TO_M3: Dict[str, float] = {
    "m3": 1.0,
    "ft3": 0.028316846592,
    "bbl": 0.158987294928,
    "stb": 0.158987294928,
    "rb": 0.158987294928,
}

#: Sürət = həcm/vaxt, SI pivotu m³/san (gerçək SI), mühərrik isə m³/gün
#: işlədir (bax `to_engine_rate`).
_SECONDS_PER_DAY = 86400.0
RATE_TO_M3_PER_S: Dict[str, float] = {
    "m3/s": 1.0,
    "m3/day": VOLUME_TO_M3["m3"] / _SECONDS_PER_DAY,
    "bbl/day": VOLUME_TO_M3["bbl"] / _SECONDS_PER_DAY,
    "stb/day": VOLUME_TO_M3["stb"] / _SECONDS_PER_DAY,
}

#: 1 lb = 0.45359237 kg (dəqiq), 1 ft³ = 0.028316846592 m³ (dəqiq).
DENSITY_TO_KG_M3: Dict[str, float] = {
    "kg/m3": 1.0,
    "lb/ft3": 0.45359237 / 0.028316846592,
}

#: Həll olmuş qaz-neft nisbəti (Rs) — "sm3/sm3" (mühərrik/METRIC) pivot.
#: `simulation/pvt/correlations.py`-dəki `SM3M3_TO_SCFSTB=5.61458` sabiti
#: ilə EYNİ faktor (o fayl TOXUNULMAYIB, bax UNITS.md — bura YALNIZ
#: idxal sərhədində istifadə üçün ayrıca mərkəzləşdirilib).
RS_TO_SM3_SM3: Dict[str, float] = {
    "sm3/sm3": 1.0,
    "scf/stb": 1.0 / 5.61458,
}

_SCALAR_TABLES: Dict[str, Dict[str, float]] = {
    "pressure": PRESSURE_TO_PA,
    "length": LENGTH_TO_M,
    "area": AREA_TO_M2,
    "permeability": PERMEABILITY_TO_M2,
    "viscosity": VISCOSITY_TO_PAS,
    "volume": VOLUME_TO_M3,
    "rate": RATE_TO_M3_PER_S,
    "density": DENSITY_TO_KG_M3,
    "solution_gor": RS_TO_SM3_SM3,
}

#: Mühərrikin faktiki gözlədiyi ("kanonik daxili") vahid — bax modul
#: docstring-i. `compressibility` üçün ayrıca `_engine_pressure` sahəsi
#: yoxdur, çünki sıxılma "1/[təzyiq vahidi]" formasındadır (bax aşağı).
ENGINE_UNITS: Dict[str, str] = {
    "pressure": "bar",
    "length": "m",
    "area": "m2",
    "permeability": "mD",
    "viscosity": "cP",
    "volume": "m3",
    "rate": "m3/day",
    "density": "kg/m3",
    #: sıxılma = 1/[təzyiq] — "bar" burada TƏZYİQ vahididir
    #: (`convert_compressibility`-yə ötürülür), "1/bar" yox.
    "compressibility": "bar",
    "solution_gor": "sm3/sm3",
}


@dataclass(frozen=True)
class Quantity:
    """(dəyər, vahid, kəmiyyət növü) üçlüyü — giriş sərhədində (CSV/
    GRDECL/UI/PVT) hər fiziki kəmiyyəti daşımaq üçün rahatlıq örtüyü.

        Quantity(3000.0, "psi", "pressure").to_engine()  # -> 206.8427 (bar)

    Bu, `to_engine_units()`-in nazik obyekt-yönümlü örtüyüdür — məcburi
    DEYİL (bütün mövcud kod birbaşa funksiyaları çağırır), yalnız
    "dəyər+vahid+növ" üçlüyünü BİR yerdə saxlamaq faydalı olanda üçündür.
    """
    value: float
    unit: str
    quantity: str

    def to_engine(self) -> float:
        return to_engine_units(self.value, self.unit, self.quantity)

    def to_si(self) -> float:
        return to_si(self.value, self.unit, self.quantity)


def known_units(quantity: str) -> Tuple[str, ...]:
    """`quantity` üçün dəstəklənən vahid adlarını qaytarır."""
    if quantity == "temperature":
        return ("K", "C", "F")
    if quantity == "compressibility":
        return tuple(PRESSURE_TO_PA)
    return tuple(_SCALAR_TABLES[quantity])


def convert(value: float, from_unit: str, to_unit: str, quantity: str) -> float:
    """`quantity` növü üzrə `from_unit`-dən `to_unit`-ə çevirir.

    Dəstəklənən `quantity`: pressure, length, area, permeability,
    viscosity, volume, rate, density (sadə miqyas faktoru), və
    temperature/compressibility (aşağıdakı xüsusi funksiyalarla, çünki
    xətti miqyas DEYİL / TƏRS miqyasdır).
    """
    if quantity == "temperature":
        return convert_temperature(value, from_unit, to_unit)
    if quantity == "compressibility":
        return convert_compressibility(value, from_unit, to_unit)
    table = _SCALAR_TABLES.get(quantity)
    if table is None:
        raise ValueError(
            f"Naməlum kəmiyyət növü: {quantity!r}. Dəstəklənən: "
            f"{tuple(_SCALAR_TABLES) + ('temperature', 'compressibility')}")
    if from_unit not in table:
        raise ValueError(f"Naməlum {quantity} vahidi: {from_unit!r}. Dəstəklənən: {tuple(table)}")
    if to_unit not in table:
        raise ValueError(f"Naməlum {quantity} vahidi: {to_unit!r}. Dəstəklənən: {tuple(table)}")
    if from_unit == to_unit:
        return value   # eyni vahid -> DƏQİQ no-op (üzən nöqtə itkisi belə YOXDUR)
    si_value = value * table[from_unit]
    return si_value / table[to_unit]


def to_si(value: float, unit: str, quantity: str) -> float:
    """Dəyəri müvafiq SI baza vahidinə çevirir (Pa/m/m²/Pa·s/m³/m³ san/kg/m³/K)."""
    if quantity == "temperature":
        return convert_temperature(value, unit, "K")
    table = _SCALAR_TABLES[quantity]
    return value * table[unit]


def from_si(value_si: float, unit: str, quantity: str) -> float:
    """SI baza vahidindəki dəyəri hədəf vahidə çevirir."""
    if quantity == "temperature":
        return convert_temperature(value_si, "K", unit)
    table = _SCALAR_TABLES[quantity]
    return value_si / table[unit]


def to_engine_units(value: float, from_unit: str, quantity: str) -> float:
    """Mühərrikin GERÇƏK işlətdiyi vahidə çevirir (bax `ENGINE_UNITS`)."""
    return convert(value, from_unit, ENGINE_UNITS[quantity], quantity)


def from_engine_units(value: float, to_unit: str, quantity: str) -> float:
    """Mühərrikin vahidindən istənilən xarici vahidə çevirir."""
    return convert(value, ENGINE_UNITS[quantity], to_unit, quantity)


# ── temperatur: xətti DEYİL (ofset var) — ayrıca işlədilir ──────────────
def convert_temperature(value: float, from_unit: str, to_unit: str) -> float:
    if from_unit == to_unit and from_unit in ("K", "C", "F"):
        return value   # eyni vahid -> DƏQİQ no-op
    if from_unit == "K":
        kelvin = value
    elif from_unit == "C":
        kelvin = value + 273.15
    elif from_unit == "F":
        kelvin = (value - 32.0) * 5.0 / 9.0 + 273.15
    else:
        raise ValueError(f"Naməlum temperatur vahidi: {from_unit!r}. Dəstəklənən: K, C, F")

    if to_unit == "K":
        return kelvin
    if to_unit == "C":
        return kelvin - 273.15
    if to_unit == "F":
        return (kelvin - 273.15) * 9.0 / 5.0 + 32.0
    raise ValueError(f"Naməlum temperatur vahidi: {to_unit!r}. Dəstəklənən: K, C, F")


# ── sıxılma: TƏRS kəmiyyətdir (1/[təzyiq]) — miqyas TƏRS istiqamətdə ────
def convert_compressibility(value: float, from_unit: str, to_unit: str) -> float:
    """`from_unit`/`to_unit` — TƏZYİQ vahidləridir (Pa, kPa, MPa, bar, psi).

    c = -(1/V)(dV/dP) tərifindən: `∂P_hədəf/∂P_mənbə = factor_mənbə/
    factor_hədəf` (zəncir qaydası), ona görə `c_hədəf = c_mənbə ·
    (∂P_mənbə/∂P_hədəf) = c_mənbə · factor_hədəf / factor_mənbə` — bu,
    ADİ kəmiyyət çevirməsinin (`value·factor_mənbə/factor_hədəf`) TAM
    TƏRSİDİR, çünki sıxılma təzyiqin TƏRSİDİR.

    Yoxlama: `co [1/psi] * BAR_TO_PSI(14.5037744) = co [1/bar]` — bu,
    `simulation/pvt/correlations.py`-dəki mövcud `vazquez_beggs_...`
    çevirməsi ilə EYNİ nəticəni verir (əl ilə yoxlanılıb, bax
    `test_compressibility_matches_existing_correlations_conversion`).
    """
    if from_unit not in PRESSURE_TO_PA or to_unit not in PRESSURE_TO_PA:
        raise ValueError(
            f"Sıxılma vahidləri təzyiq vahidi olmalıdır (Pa/kPa/MPa/bar/psi), "
            f"alındı: {from_unit!r}, {to_unit!r}")
    if from_unit == to_unit:
        return value
    return value * PRESSURE_TO_PA[to_unit] / PRESSURE_TO_PA[from_unit]


# ── açıq adlandırılmış rahatlıq funksiyaları (tapşırıqda tələb olunan
#    konkret cütlər) — hamısı yuxarıdakı `convert()`-in nazik örtüyüdür ──
def psi_to_pa(value: float) -> float: return convert(value, "psi", "Pa", "pressure")
def pa_to_psi(value: float) -> float: return convert(value, "Pa", "psi", "pressure")
def bar_to_pa(value: float) -> float: return convert(value, "bar", "Pa", "pressure")
def pa_to_bar(value: float) -> float: return convert(value, "Pa", "bar", "pressure")
def psi_to_bar(value: float) -> float: return convert(value, "psi", "bar", "pressure")
def bar_to_psi(value: float) -> float: return convert(value, "bar", "psi", "pressure")
def ft_to_m(value: float) -> float: return convert(value, "ft", "m", "length")
def m_to_ft(value: float) -> float: return convert(value, "m", "ft", "length")
def md_to_m2(value: float) -> float: return convert(value, "mD", "m2", "permeability")
def m2_to_md(value: float) -> float: return convert(value, "m2", "mD", "permeability")
def darcy_to_m2(value: float) -> float: return convert(value, "D", "m2", "permeability")
def cp_to_pas(value: float) -> float: return convert(value, "cP", "Pa.s", "viscosity")
def pas_to_cp(value: float) -> float: return convert(value, "Pa.s", "cP", "viscosity")
def bbl_to_m3(value: float) -> float: return convert(value, "bbl", "m3", "volume")
def m3_to_bbl(value: float) -> float: return convert(value, "m3", "bbl", "volume")
def stb_per_day_to_m3_per_day(value: float) -> float:
    return convert(value, "stb/day", "m3/day", "rate")
def m3_per_day_to_stb_per_day(value: float) -> float:
    return convert(value, "m3/day", "stb/day", "rate")
