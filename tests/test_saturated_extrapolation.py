"""G7 — doymuş qolun cədvəldən yuxarı XƏTTİ uzadılması (Q-32).

OPM/Eclipse doymuş PVTO cədvəlini SON SEQMENTİN meyli ilə uzadır. Qayda
MƏNBƏDƏN oxundu (Seans 40):

* `opm-common/opm/material/fluidsystems/blackoilpvt/LiveOilPvt.hpp:512` —
  `saturatedGasDissolutionFactor` cədvəli `extrapolate=true` ilə çağırılır;
* `opm-common/opm/material/common/Tabulated1DFunction.hpp:266-283` — bu
  bayraq qoyulanda son seqmentin xətti düsturu cədvəldən KƏNARDA da tətbiq
  olunur ("extended beyond its range by straight lines").

NİYƏ YALNIZ DECK YOLUNDA: korrelyasiya cədvəlində Pb-dən yuxarı Rs platosu
HƏQİQİ fizikadır (neftin tərkibi sabitdir), deck-də isə sadəcə cədvəlin
məlumatı bitir.

ÖLÇÜLMÜŞ TƏSİRİ (SPE1CASE2, `ISH_HESABATI.md` → Seans 41): istismarçının
BHP-si 1034-cü gündə −60.7 % → −1.2 %, FGOR +389 % → +38.7 %.
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

from imex2d.domain.unit_conversions import convert, to_engine_units
from imex2d.io.pvt_io import read_deck_pvt
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table

#: SPE1CASE2-nin PVT quruluşu — doymuş qol 5014.7 psia-da BİTİR.
DECK = """RUNSPEC
PROPS

PVTW
   4017.55  1.038   3.22E-6   0.318   0.0 /

PVDG
    14.7     166.666     0.0080
    400.0    5.4924      0.0096
    4014.7   0.5312      0.0126
    9014.7   0.2369      0.0140
/

PVTO
    0.0010   14.7      1.0620   1.0400 /
    0.0905   400.0     1.1500   0.9750 /
    1.2700   4014.7    1.6950   0.5100
             9014.7    1.5790   0.7400 /
    1.6180   5014.7    1.8270   0.4490
             9014.7    1.7370   0.6310 /
/

DENSITY
   53.66 64.49 0.0533 /
"""

#: Deck-in son İKİ doymuş düyünü — uzantının meyli bunlardan gəlir.
P_LOW = to_engine_units(4014.7, "psi", "pressure")
P_TOP = to_engine_units(5014.7, "psi", "pressure")
RS_LOW = convert(1.2700, "Mscf/stb", "sm3/sm3", "solution_gor")
RS_TOP = convert(1.6180, "Mscf/stb", "sm3/sm3", "solution_gor")
BO_LOW, BO_TOP = 1.6950, 1.8270
MU_LOW, MU_TOP = 0.5100, 0.4490
DRS_DP = (RS_TOP - RS_LOW) / (P_TOP - P_LOW)
DBO_DP = (BO_TOP - BO_LOW) / (P_TOP - P_LOW)
DMU_DP = (MU_TOP - MU_LOW) / (P_TOP - P_LOW)


def _deck_path(text: str = DECK) -> str:
    handle, path = tempfile.mkstemp(suffix=".DATA")
    with os.fdopen(handle, "w", encoding="utf-8") as stream:
        stream.write(text)
    return path


@pytest.fixture
def provider():
    """İSTEHSAL YOLU: cədvəl + qollar (G3), yəni uzantı işə düşür."""
    path = _deck_path()
    try:
        deck = read_deck_pvt(path)
    finally:
        os.unlink(path)
    return BlackOilPVTProvider(deck.to_pvt_table(), oil_branches=deck.oil.branches)


@pytest.fixture
def correlation_provider():
    """Qolsuz provider — uzantı QURULMAMALIDIR."""
    return BlackOilPVTProvider(build_pvt_table(
        pressure_min=1.0, pressure_max=700.0, n_points=40,
        bubble_point_bar=345.75, include_gas=True))


# ═════════════════ uzantının özü ══════════════════════════════════════

def test_extension_is_built_from_the_last_saturated_segment(provider):
    extension = provider._sat_extension
    assert extension is not None
    assert extension["pressure"] == pytest.approx(P_TOP, rel=1e-9)
    assert extension["rs"] == pytest.approx(RS_TOP, rel=1e-9)
    assert extension["solution_gor"] == pytest.approx(DRS_DP, rel=1e-9)
    assert extension["oil_fvf"] == pytest.approx(DBO_DP, rel=1e-9)
    assert extension["oil_viscosity"] == pytest.approx(DMU_DP, rel=1e-9)


@pytest.mark.parametrize("delta_bar", [0.0, 10.0, 50.0, 150.0])
def test_solution_gor_extends_linearly(provider, delta_bar):
    value = float(provider.solution_gor(np.array([P_TOP + delta_bar]))[0])
    assert value == pytest.approx(RS_TOP + DRS_DP * delta_bar, rel=1e-10)


def test_measured_value_matches_the_hand_calculation(provider):
    """6150 psia-da Rs_sat = 2.0131 Mscf/STB (Seans 40-da əl ilə yoxlandı)."""
    bar = to_engine_units(6150.0, "psi", "pressure")
    value = float(provider.solution_gor(np.array([bar]))[0])
    field = convert(value, "sm3/sm3", "Mscf/stb", "solution_gor")
    assert field == pytest.approx(2.0131, abs=5e-4)


def test_oil_fvf_and_viscosity_follow_the_saturated_trend(provider):
    delta = 100.0
    assert float(provider.oil_fvf(np.array([P_TOP + delta]))[0]) == \
        pytest.approx(BO_TOP + DBO_DP * delta, rel=1e-9)
    assert float(provider.oil_viscosity(np.array([P_TOP + delta]))[0]) == \
        pytest.approx(MU_TOP + DMU_DP * delta, rel=1e-9)


def test_inside_the_table_nothing_changes(provider):
    """Cədvəl daxilində qiymətlər deck-in öz düyünləridir."""
    assert float(provider.solution_gor(np.array([P_LOW]))[0]) == \
        pytest.approx(RS_LOW, rel=1e-9)
    assert float(provider.oil_fvf(np.array([P_LOW]))[0]) == \
        pytest.approx(BO_LOW, rel=1e-9)
    assert float(provider.oil_viscosity(np.array([P_TOP]))[0]) == \
        pytest.approx(MU_TOP, rel=1e-9)


def test_correlation_table_keeps_the_plateau(correlation_provider):
    """Qollar verilmirsə uzantı yoxdur — korrelyasiya yolu toxunulmamışdır."""
    assert correlation_provider._sat_extension is None
    above = np.array([400.0, 600.0])
    values = correlation_provider.solution_gor(above)
    assert values[0] == pytest.approx(values[1], rel=1e-12)


# ═════════════════ törəmələr — sonlu fərqlə ═══════════════════════════

@pytest.mark.parametrize("name", ["solution_gor", "oil_fvf", "oil_viscosity"])
def test_derivatives_match_finite_difference_above_the_top(provider, name):
    """Analitik meyl uzantı zonasında da qalıqla uyğundur.

    NİYƏ VACİB: uzantı qoyulub törəmə sıfır qalsaydı, Nyuton YANLIŞ
    Jakobianla işləyərdi — səssiz səhvin ən pis sinfi.
    """
    value = getattr(provider, name)
    derivative = getattr(provider, f"{name}_derivative")
    point, step = P_TOP + 60.0, 1e-4
    numeric = (float(value(np.array([point + step]))[0])
               - float(value(np.array([point - step]))[0])) / (2 * step)
    assert float(derivative(np.array([point]))[0]) == \
        pytest.approx(numeric, rel=1e-6)


def test_saturation_pressure_is_the_exact_inverse(provider):
    """Pb(Rs_sat(p)) = p — uzantı zonasında da (keçid kəsilməz qalır)."""
    for point in (P_TOP + 5.0, P_TOP + 60.0, P_TOP + 150.0):
        rs = provider.solution_gor(np.array([point]))
        assert float(provider.saturation_pressure(rs)[0]) == \
            pytest.approx(point, rel=1e-10)


def test_saturation_pressure_slope_matches_finite_difference(provider):
    """dPb/dRs uzantıda `1/(dRs/dp)`-dir — əvvəl orada SIFIR idi (TB-3)."""
    rs, step = RS_TOP + 20.0, 1e-6
    numeric = (float(provider.saturation_pressure(np.array([rs + step]))[0])
               - float(provider.saturation_pressure(np.array([rs - step]))[0])
               ) / (2 * step)
    analytic = float(provider._saturation_pressure_slope(np.array([rs]))[0])
    assert analytic == pytest.approx(numeric, rel=1e-6)
    assert analytic == pytest.approx(1.0 / DRS_DP, rel=1e-10)


# ═════════════════ qoruyucu ═══════════════════════════════════════════

def test_extension_never_returns_a_non_positive_viscosity(provider, caplog):
    """Xətti uzantı kifayət qədər yüksək təzyiqdə μo-nu sıfırdan keçirərdi.

    SPE1-də təzyiq 7600 psia-ya qədər qalxır və hədd İŞƏ DÜŞMÜR; bu test
    onun HƏQİQƏTƏN işlədiyini uzaq nöqtədə yoxlayır.
    """
    far = P_TOP + abs(MU_TOP / DMU_DP) * 2.0
    with caplog.at_level("WARNING"):
        value = float(provider.oil_viscosity(np.array([far]))[0])
    assert value > 0.0
    assert any("G7" in record.message for record in caplog.records)
