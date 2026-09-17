"""Q-34 — doymamış qol deck düyünləri arasında XƏTTİdir.

OPM cədvəl düyünləri arasında xətti interpolyasiya edir
(`Tabulated1DFunction.hpp:282`: `y0 + (y1−y0)·(x−x0)/(x1−x0)`). Bizdə
əvvəl ÜSTƏL qanun var idi (Q-28/Q-30) — ÖLÇÜLDÜ: μo-da 1.97 %, Bo-da
0.06 % fərq (5500–6500 psia aralığında).

QURULUŞ: lövbər (Pb-dəki doymuş qiymət) DƏYİŞMİR — Seans 36-da qazanılan
doymuş/doymamış keçid kəsilməzliyi qorunur; yalnız qolun FORMASI dəyişir:

    Bo(p, Rs)  = Bo_sat(Pb) + slope_Bo(Rs)·(p − Pb)
    μo(p, Rs)  = μo_sat(Pb) + slope_μ(Rs)·(p − Pb)

İki düyünlü qolda (SPE1 belədir) bu, OPM-in interpolyasiyası ilə
EYNİDİR. Korrelyasiya yolunda (deck qolları verilmir) üstəl qanun QALIR.
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

#: Rs = 1.27 qolunun deck düyünləri (FIELD)
P_LOW, P_HIGH = 4014.7, 9014.7
BO_LOW, BO_HIGH = 1.6950, 1.5790
MU_LOW, MU_HIGH = 0.5100, 0.7400
RS = convert(1.2700, "Mscf/stb", "sm3/sm3", "solution_gor")


def _deck_path(text: str = DECK) -> str:
    handle, path = tempfile.mkstemp(suffix=".DATA")
    with os.fdopen(handle, "w", encoding="utf-8") as stream:
        stream.write(text)
    return path


@pytest.fixture
def provider():
    path = _deck_path()
    try:
        deck = read_deck_pvt(path)
    finally:
        os.unlink(path)
    return BlackOilPVTProvider(deck.to_pvt_table(), oil_branches=deck.oil.branches)


@pytest.fixture
def correlation_provider():
    return BlackOilPVTProvider(build_pvt_table(
        pressure_min=1.0, pressure_max=700.0, n_points=40,
        bubble_point_bar=345.75, include_gas=True))


def _bar(psia):
    return np.array([to_engine_units(psia, "psi", "pressure")])


# ═════════════════ qiymətlər — deck ilə DƏQİQ ═════════════════════════

@pytest.mark.parametrize("psia", [4014.7, 4500.0, 4800.0, 5500.0,
                                  6500.0, 7500.0, 9014.7])
def test_values_match_the_deck_linear_interpolation(provider, psia):
    """ƏSAS XASSƏ: qol deck düyünləri arasında məhz xəttidir."""
    rs = np.array([RS])
    bo_expected = np.interp(psia, [P_LOW, P_HIGH], [BO_LOW, BO_HIGH])
    mu_expected = np.interp(psia, [P_LOW, P_HIGH], [MU_LOW, MU_HIGH])
    assert float(provider.oil_fvf_undersaturated(_bar(psia), rs)[0]) == \
        pytest.approx(bo_expected, rel=1e-9)
    assert float(provider.oil_viscosity_undersaturated(_bar(psia), rs)[0]) == \
        pytest.approx(mu_expected, rel=1e-9)


def test_anchor_at_the_bubble_point_is_unchanged(provider):
    """Pb-də doymuş qiymətlə üst-üstə düşür — keçid KƏSİLMƏZDİR."""
    rs = np.array([RS])
    pb = provider.saturation_pressure(rs)
    assert float(provider.oil_fvf_undersaturated(pb, rs)[0]) == \
        pytest.approx(BO_LOW, rel=1e-9)
    assert float(provider.oil_viscosity_undersaturated(pb, rs)[0]) == \
        pytest.approx(MU_LOW, rel=1e-9)


def test_correlation_path_still_uses_the_power_law(correlation_provider):
    """Qollar verilmirsə davranış DƏYİŞMİR (üstəl qanun).

    Yoxlama: üstəl qanun xətti DEYİL, ona görə orta nöqtədəki qiymət
    uclardan çəkilən düzün üstündə/altında olur.
    """
    assert correlation_provider._branch_rs is None
    rs = np.array([correlation_provider.table.solution_gor[-1]])
    low, high = np.array([350.0]), np.array([600.0])
    middle = np.array([475.0])
    values = [float(correlation_provider.oil_viscosity_undersaturated(p, rs)[0])
              for p in (low, middle, high)]
    linear_middle = 0.5 * (values[0] + values[2])
    assert values[1] != pytest.approx(linear_middle, rel=1e-6)


# ═════════════════ törəmələr ══════════════════════════════════════════

def test_pressure_derivative_is_the_branch_slope(provider):
    """∂/∂p qolun öz meylidir — sabit və dəqiq."""
    rs = np.array([RS])
    span = to_engine_units(P_HIGH, "psi", "pressure") - \
        to_engine_units(P_LOW, "psi", "pressure")
    expected_bo = (BO_HIGH - BO_LOW) / span
    expected_mu = (MU_HIGH - MU_LOW) / span
    for psia in (4500.0, 6000.0, 8000.0):
        d_bo = provider.oil_fvf_undersaturated_derivatives(_bar(psia), rs)[0]
        d_mu = provider.oil_viscosity_undersaturated_derivatives(_bar(psia), rs)[0]
        assert float(d_bo[0]) == pytest.approx(expected_bo, rel=1e-9)
        assert float(d_mu[0]) == pytest.approx(expected_mu, rel=1e-9)


@pytest.mark.parametrize("rs_field", [1.10, 1.40, 1.55])
def test_rs_derivative_matches_finite_difference_between_nodes(provider, rs_field):
    """∂/∂Rs qolLAR ARASINDA dəqiqdir.

    DÜYÜNÜN ÖZÜNDƏ yoxlanılmır: parçalı xətti modeldə orada törəmə
    birqiymətli deyil (mərkəzi fərq iki fərqli meyli ortalayır) — bu,
    cədvəl əsaslı bütün törəmələrimiz üçün eynidir.
    """
    rs = np.array([convert(rs_field, "Mscf/stb", "sm3/sm3", "solution_gor")])
    point, step = _bar(6000.0), 1e-3
    for name in ("oil_fvf_undersaturated", "oil_viscosity_undersaturated"):
        value = getattr(provider, name)
        derivative = getattr(provider, f"{name}_derivatives")
        numeric = float((value(point, rs + step) - value(point, rs - step))[0]
                        / (2 * step))
        assert float(derivative(point, rs)[1][0]) == \
            pytest.approx(numeric, rel=1e-6)


# ═════════════════ fiziki yoxlama ═════════════════════════════════════

def test_branch_with_rising_bo_is_rejected():
    """Doymamış Bo təzyiqlə ARTARSA — açıq xəta (səssiz qəbul yox)."""
    # Rs = 1.27 qolunun doymamış Bo-su 1.5790 → 1.9790 (yəni təzyiqlə
    # ARTIR). Yalnız bu rəqəm dəyişir, deck-in formatı toxunulmur.
    assert DECK.count("1.5790") == 1
    broken = DECK.replace("1.5790", "1.9790")
    path = _deck_path(broken)
    try:
        deck = read_deck_pvt(path)
    finally:
        os.unlink(path)
    with pytest.raises(ValueError, match="fiziki deyil"):
        BlackOilPVTProvider(deck.to_pvt_table(), oil_branches=deck.oil.branches)
