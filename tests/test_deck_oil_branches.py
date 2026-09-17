"""G3 (1-ci hissə) — doymamış Bo və μo HƏR PVTO qolunun öz parametrləri ilə.

ÖLÇÜLDÜ (Seans 37, SPE1CASE2.DATA): iki doymamış qolun parametrləri FƏRQLİDİR —

    Rs 1.270 Mscf/STB:  c_o = 2.056e-4 1/bar,  n = 0.460
    Rs 1.618 Mscf/STB:  c_o = 1.832e-4 1/bar,  n = 0.580

Tək doymamış hissəli `PVTTable` bunu daşıya bilmir: `reference_rs = 1.27` ilə
Rs 1.27-də kəsilir (Rs 1.618 qolunda Bo −9 %, μo +17 %), defolt qolla isə c_o
ehtiyat qiymətə düşür (Bo −80 %). Həll: `BlackOilPVTProvider(...,
oil_branches=...)` — c_o(Rs) və n(Rs) qollar arasında Rs üzrə interpolyasiya.

BU ADDIMDA MÜHƏRRİK TOXUNULMUR. Yeni törəmə hədləri (dc_o/dRs, dn/dRs)
provider-in `*_derivatives` metodlarındadır; mühərrik səviyyəsində yoxlama
növbəti commit-dir (G3b).
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from test_pvt_io import DECK, _deck
from imex2d.io.pvt_io import OilBranch, read_deck_pvt
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider


def _deck_pvt():
    return read_deck_pvt(_deck(DECK))


def _provider(deck=None):
    deck = deck or _deck_pvt()
    return BlackOilPVTProvider(deck.to_pvt_table(), oil_branches=deck.oil.branches)


def _undersaturated_branches(deck):
    return [branch for branch in deck.oil.branches if branch.has_undersaturated]


def _between_branches(deck):
    low, high = _undersaturated_branches(deck)
    rs = np.array([low.solution_gor + 0.25 * (high.solution_gor - low.solution_gor),
                   0.5 * (low.solution_gor + high.solution_gor),
                   high.solution_gor - 5.0])
    pressure = np.array([high.bubble_point + 20.0, high.bubble_point + 150.0,
                         high.pressure[-1] - 10.0])
    return pressure, rs


# ═══════════════════════ deck-in təkrarlanması ═══════════════════════

def test_every_deck_row_is_reproduced_exactly():
    """İki sətirli qolda doymamış sətir DƏQİQ çıxmalıdır — hər iki qolda."""
    deck = _deck_pvt()
    provider = _provider(deck)
    for branch in _undersaturated_branches(deck):
        rs = np.full(branch.pressure.size, branch.solution_gor)
        assert provider.oil_fvf_undersaturated(branch.pressure, rs) == pytest.approx(
            branch.formation_volume_factor, rel=1e-12)
        assert provider.oil_viscosity_undersaturated(branch.pressure, rs) == \
            pytest.approx(branch.viscosity, rel=1e-12)


def test_branch_parameters_are_the_measured_ones():
    provider = _provider()
    assert provider._branch_co == pytest.approx([2.0564e-4, 1.8317e-4], rel=1e-3)
    assert provider._branch_n == pytest.approx([0.4602, 0.5802], rel=1e-3)


def test_parameters_are_interpolated_linearly_between_branches():
    provider = _provider()
    low, high = provider._branch_rs
    middle = np.array([0.5 * (low + high)])
    co, co_rs = provider._compressibility(middle)
    n, n_rs = provider._viscosity_exponent(middle)
    assert co[0] == pytest.approx(provider._branch_co.mean(), rel=1e-12)
    assert n[0] == pytest.approx(provider._branch_n.mean(), rel=1e-12)
    assert co_rs[0] == pytest.approx(np.diff(provider._branch_co)[0] / (high - low))
    assert n_rs[0] == pytest.approx(np.diff(provider._branch_n)[0] / (high - low))


def test_parameters_are_held_constant_outside_the_branches():
    """Ekstrapolyasiya YOX — ən yaxın doymamış qolun parametri, törəmə sıfır."""
    provider = _provider()
    rs = np.array([provider._branch_rs[0] - 50.0, provider._branch_rs[-1] + 10.0])
    co, co_rs = provider._compressibility(rs)
    n, n_rs = provider._viscosity_exponent(rs)
    assert co == pytest.approx([provider._branch_co[0], provider._branch_co[-1]])
    assert n == pytest.approx([provider._branch_n[0], provider._branch_n[-1]])
    assert np.all(co_rs == 0.0) and np.all(n_rs == 0.0)


def test_undersaturated_branch_meets_the_saturated_curve_at_the_heads():
    """Kəsilməzlik: Rs = Rs_sat(Pb) olanda doymamış qol doymuş başla üst-üstədir."""
    deck = _deck_pvt()
    provider = _provider(deck)
    rs_heads, pb_heads, bo_heads, mu_heads = deck.oil.saturated
    assert provider.oil_fvf_undersaturated(pb_heads, rs_heads) == pytest.approx(
        bo_heads, rel=1e-12)
    assert provider.oil_viscosity_undersaturated(pb_heads, rs_heads) == \
        pytest.approx(mu_heads, rel=1e-12)


# ═══════════════════════ niyə lazımdır — ölçülmüş xəta ═════════════════

def test_single_branch_table_misses_the_second_branch():
    """KÖHNƏ YOL: `reference_rs = 1.27` ilə Rs 1.618 qolu kəskin səhvdir."""
    deck = _deck_pvt()
    low, high = _undersaturated_branches(deck)
    old = BlackOilPVTProvider(deck.to_pvt_table(reference_rs=low.solution_gor))
    rs = np.full(high.pressure.size, high.solution_gor)
    bo_error = old.oil_fvf_undersaturated(high.pressure, rs) / \
        high.formation_volume_factor - 1.0
    mu_error = old.oil_viscosity_undersaturated(high.pressure, rs) / high.viscosity - 1.0
    assert np.max(np.abs(bo_error)) > 0.05, bo_error
    assert np.max(np.abs(mu_error)) > 0.15, mu_error


# ═══════════════════════ törəmələr — sonlu fərq ══════════════════════

@pytest.mark.parametrize("name", ["oil_fvf", "oil_viscosity"])
def test_derivatives_match_finite_difference_between_branches(name):
    deck = _deck_pvt()
    provider = _provider(deck)
    value = getattr(provider, f"{name}_undersaturated")
    derivatives = getattr(provider, f"{name}_undersaturated_derivatives")
    pressure, rs = _between_branches(deck)
    d_dp, d_drs = derivatives(pressure, rs)
    h = 1e-4
    numeric_p = (value(pressure + h, rs) - value(pressure - h, rs)) / (2 * h)
    numeric_rs = (value(pressure, rs + h) - value(pressure, rs - h)) / (2 * h)
    assert d_dp == pytest.approx(numeric_p, rel=1e-6)
    assert d_drs == pytest.approx(numeric_rs, rel=1e-6)


@pytest.mark.parametrize("name, slope",
                         [("oil_fvf", "_branch_bo_slope_rs"),
                          ("oil_viscosity", "_branch_mu_slope_rs")])
def test_the_new_parameter_slope_term_is_necessary(name, slope):
    """ZƏRURİLİK: qol parametrinin Rs üzrə meyli söndürüləndə ∂/∂Rs pozulur.

    Q-34-DƏN SONRA (Seans 43) deck yolunda qol XƏTTİdir, yəni
    qiymətləndirmədə `c_o`/`n` YOX, `_branch_bo_slope`/`_branch_mu_slope`
    işlədilir. Testin məqsədi dəyişmir — yoxlanılan hədd həmin
    parametrlərin Rs üzrə meylidir (`..._slope_rs`). `c_o`/`n` isə
    korrelyasiya yolunda və diaqnostikada qalır.
    """
    deck = _deck_pvt()
    provider = _provider(deck)
    value = getattr(provider, f"{name}_undersaturated")
    pressure, rs = _between_branches(deck)
    h = 1e-4
    numeric_rs = (value(pressure, rs + h) - value(pressure, rs - h)) / (2 * h)

    setattr(provider, slope, np.zeros_like(getattr(provider, slope)))
    _, without_term = getattr(provider, f"{name}_undersaturated_derivatives")(pressure, rs)
    relative = np.max(np.abs(without_term - numeric_rs) / np.abs(numeric_rs))
    assert relative > 1e-2, relative


# ═══════════════════════ açıq xətalar və xəbərdarlıqlar ════════════════

def test_truncated_table_warns_that_solution_gor_will_be_capped():
    deck = _deck_pvt()
    low = _undersaturated_branches(deck)[0]
    logger = logging.getLogger("imex2d.simulation.pvt.black_oil")
    records = []

    class Collector(logging.Handler):
        def emit(self, record):
            records.append(record.getMessage())

    handler = Collector(level=logging.WARNING)
    logger.addHandler(handler)
    try:
        BlackOilPVTProvider(deck.to_pvt_table(reference_rs=low.solution_gor),
                            oil_branches=deck.oil.branches)
    finally:
        logger.removeHandler(handler)
    assert any("KƏSİLƏCƏK" in message for message in records), records


def test_branches_from_another_table_are_rejected():
    deck = _deck_pvt()
    foreign = [OilBranch(branch.solution_gor * 1.3, branch.pressure,
                         branch.formation_volume_factor, branch.viscosity)
               for branch in deck.oil.branches]
    with pytest.raises(ValueError, match="uyğun"):
        BlackOilPVTProvider(deck.to_pvt_table(), oil_branches=foreign)


def test_branches_without_undersaturated_rows_are_rejected():
    deck = _deck_pvt()
    heads = [OilBranch(branch.solution_gor, branch.pressure[:1],
                       branch.formation_volume_factor[:1], branch.viscosity[:1])
             for branch in deck.oil.branches]
    with pytest.raises(ValueError, match="doymamış sətir yoxdur"):
        BlackOilPVTProvider(deck.to_pvt_table(), oil_branches=heads)


def test_non_physical_branch_is_rejected():
    deck = _deck_pvt()
    branches = list(deck.oil.branches)
    top = branches[-1]
    branches[-1] = OilBranch(top.solution_gor, top.pressure,
                             top.formation_volume_factor[::-1].copy(),
                             top.viscosity)
    with pytest.raises(ValueError, match="fiziki deyil"):
        BlackOilPVTProvider(deck.to_pvt_table(), oil_branches=branches)


def test_without_branches_the_single_parameter_path_is_kept():
    """GERİYƏ UYĞUNLUQ: `oil_branches` verilməyəndə yeni hədlər YOXDUR."""
    deck = _deck_pvt()
    low = _undersaturated_branches(deck)[0]
    provider = BlackOilPVTProvider(deck.to_pvt_table(reference_rs=low.solution_gor))
    assert provider._branch_rs is None
    assert provider._compressibility(np.array([200.0]))[1] is None
    assert provider._viscosity_exponent(np.array([200.0]))[1] is None
