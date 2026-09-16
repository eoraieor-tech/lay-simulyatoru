"""G1 — Eclipse deck-dən PVT oxuyucuları (PVTW, PVDG, PVTO).

SPE1CASE2 PVT-ni üç açar sözlə verir və onlar FƏRQLİ təzyiq şəbəkələrindədir;
`PVTTable` isə TƏK şəbəkəlidir. Ona görə oxuma (itkisiz) və birləşdirmə
(təqribi) AYRI addımlardır.

ÖLÇÜLMÜŞ TƏHLÜKƏ (testlə kilidlənib): deck-də hər Rs qolunda çox vaxt cəmi
BİR doymamış sətir olur. `BlackOilPVTProvider` doymamış sıxılmanı ən azı iki
düyündən fit edir — belə halda o, səssizcə ehtiyat qiymətə düşürdü (ölçüldü:
həqiqi 2.06e-4 yerinə 3.12e-3 1/bar, 15 dəfə böyük). İndi xəbərdarlıq verilir.
"""

from __future__ import annotations

import logging
import os
import tempfile

import numpy as np
import pytest

from imex2d.domain.unit_conversions import convert, to_engine_units
from imex2d.io.pvt_io import (DeckPvt, PvtFormatError, read_deck_pvt, read_pvdg,
                              read_pvto, read_pvtw)
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider

DECK = """RUNSPEC
PROPS

PVTW
-- Pref     Bw      Cw        muw     Cmu
   4017.55  1.038   3.22E-6   0.318   0.0 /

PVDG
--  P        Bg          mug
    14.7     166.666     0.0080
    400.0    5.4924      0.0096
    4014.7   0.5312      0.0126
    9014.7   0.2369      0.0140
/

PVTO
--  Rs       Pb        Bo       muo
    0.0010   14.7      1.0620   1.0400 /
    0.0905   400.0     1.1500   0.9750 /
    1.2700   4014.7    1.6950   0.5100
             9014.7    1.5790   0.7400 /
    1.6180   5014.7    1.8270   0.4490
             9014.7    1.7260   0.6050 /
/

DENSITY
   53.66 64.49 0.0533 /
"""

#: SPE1-in öz rəqəmləri, mühərrik vahidlərinə çevrilmiş
RS_INITIAL = convert(1.2700, "Mscf/stb", "sm3/sm3", "solution_gor")
PB_INITIAL = to_engine_units(4014.7, "psi", "pressure")


def _deck(text: str = DECK) -> str:
    handle, path = tempfile.mkstemp(suffix=".DATA")
    with os.fdopen(handle, "w", encoding="utf-8") as file:
        file.write(text)
    return path


# ═══════════════════════════════ PVTW ════════════════════════════════

def test_pvtw_converts_field_units():
    water = read_pvtw(_deck())
    assert water.reference_pressure == pytest.approx(
        to_engine_units(4017.55, "psi", "pressure"))
    assert water.formation_volume_factor == pytest.approx(1.038), "Bw ölçüsüzdür"
    assert water.compressibility == pytest.approx(
        to_engine_units(3.22e-6, "psi", "compressibility"))
    assert water.viscosity == pytest.approx(0.318)


def test_water_fvf_decreases_with_pressure():
    water = read_pvtw(_deck())
    pressure = np.array([100.0, 300.0, 500.0])
    values = water.water_fvf(pressure)
    assert np.all(np.diff(values) < 0.0)


def test_water_viscosity_is_constant_when_viscosibility_is_zero():
    """SPE1-də `Cmu = 0` — μw təzyiqdən asılı deyil."""
    water = read_pvtw(_deck())
    values = water.water_viscosity(np.array([50.0, 400.0]))
    assert values[0] == pytest.approx(values[1]) == pytest.approx(0.318)


# ═══════════════════════════════ PVDG ════════════════════════════════

def test_pvdg_reads_all_rows_and_converts_bg():
    gas = read_pvdg(_deck())
    assert gas.pressure.size == 4
    assert gas.pressure[-1] == pytest.approx(
        to_engine_units(9014.7, "psi", "pressure"))
    assert gas.formation_volume_factor[-1] == pytest.approx(
        convert(0.2369, "rb/Mscf", "m3/sm3", "gas_fvf"))


def test_pvdg_bg_decreases_with_pressure():
    gas = read_pvdg(_deck())
    assert np.all(np.diff(gas.formation_volume_factor) < 0.0)


def test_pvdg_rejects_increasing_bg():
    bad = DECK.replace("    9014.7   0.2369      0.0140",
                       "    9014.7   9.9990      0.0140")
    with pytest.raises(PvtFormatError, match="Bg"):
        read_pvdg(_deck(bad))


# ═══════════════════════════════ PVTO ════════════════════════════════

def test_pvto_splits_branches_by_solution_gor():
    oil = read_pvto(_deck())
    assert len(oil) == 4
    assert [branch.pressure.size for branch in oil.branches] == [1, 1, 2, 2]


def test_pvto_converts_rs_and_bubble_point():
    oil = read_pvto(_deck())
    branch = oil.branch_for(RS_INITIAL)
    assert branch.solution_gor == pytest.approx(RS_INITIAL)
    assert branch.bubble_point == pytest.approx(PB_INITIAL)
    assert branch.has_undersaturated


def test_pvto_undersaturated_rows_belong_to_the_same_branch():
    """Davam sətri YENİ qol açmamalıdır — Rs eyni qalır."""
    branch = read_pvto(_deck()).branch_for(RS_INITIAL)
    assert branch.pressure.size == 2
    assert branch.pressure[1] > branch.pressure[0]
    assert branch.formation_volume_factor[1] < branch.formation_volume_factor[0]
    assert branch.viscosity[1] > branch.viscosity[0], "doymamış neft qatılaşır"


def test_pvto_rejects_non_increasing_rs():
    bad = DECK.replace("    1.6180   5014.7", "    0.5000   5014.7")
    with pytest.raises(PvtFormatError, match="Rs"):
        read_pvto(_deck(bad))


def test_missing_keyword_is_an_explicit_error():
    with pytest.raises(PvtFormatError, match="PVTO"):
        read_pvto(_deck("PVTW\n 4017.55 1.038 3.22E-6 0.318 0.0 /\n"))


def test_reader_stops_at_the_next_keyword():
    """`DENSITY` sətirləri PVTO-ya qarışmamalıdır."""
    oil = read_pvto(_deck())
    assert all(branch.solution_gor < 300.0 for branch in oil.branches)


# ═══════════════════════════ birləşdirmə ═════════════════════════════

def test_merged_table_is_valid_and_accepted_by_the_provider():
    table = read_deck_pvt(_deck()).to_pvt_table(reference_rs=RS_INITIAL)
    assert table.validate() == []
    assert table.has_gas_phase
    provider = BlackOilPVTProvider(table)
    assert provider.oil_fvf(200.0) > 1.0


def test_reference_rs_selects_the_branch():
    deck = read_deck_pvt(_deck())
    assert deck.to_pvt_table(reference_rs=RS_INITIAL).bubble_point == pytest.approx(
        PB_INITIAL)
    highest = deck.oil.branches[-1].bubble_point
    assert deck.to_pvt_table().bubble_point == pytest.approx(highest)


def test_solution_gor_is_a_plateau_above_the_bubble_point():
    table = read_deck_pvt(_deck()).to_pvt_table(reference_rs=RS_INITIAL)
    above = table.pressure > table.bubble_point
    assert np.any(above)
    assert np.allclose(table.solution_gor[above], RS_INITIAL)


def test_grid_is_the_union_of_the_deck_pressures():
    deck = read_deck_pvt(_deck())
    table = deck.to_pvt_table(reference_rs=RS_INITIAL)
    for pressure in deck.gas.pressure:
        assert np.any(np.isclose(table.pressure, pressure))
    for branch in deck.oil.branches:
        assert np.any(np.isclose(table.pressure, branch.bubble_point))


def test_all_branches_survive_in_the_container():
    """İTKİSİZ oxuma: cədvələ köçürülməyən qollar `DeckPvt`-də qalır."""
    deck = read_deck_pvt(_deck())
    deck.to_pvt_table(reference_rs=RS_INITIAL)
    assert len(deck.oil) == 4
    assert sum(1 for branch in deck.oil.branches if branch.has_undersaturated) == 2


def test_thin_undersaturated_branch_warns_instead_of_failing_silently():
    """ÖLÇÜLMÜŞ TƏHLÜKƏ: bir doymamış düyün c_o fit-ini mümkünsüz edir.

    Provider o zaman səssizcə ehtiyat qiymətə düşür (15 dəfə böyük) — indi
    birləşdirmə bundan AÇIQ xəbərdarlıq verir.
    """
    deck = read_deck_pvt(_deck())
    logger = logging.getLogger("imex2d.io.pvt_io")
    records = []

    class Collector(logging.Handler):
        def emit(self, record):
            records.append(record.getMessage())

    handler = Collector(level=logging.WARNING)
    logger.addHandler(handler)
    try:
        table = deck.to_pvt_table()          # ən böyük Rs qolu — bir düyün
    finally:
        logger.removeHandler(handler)

    above = int(np.count_nonzero(table.pressure > table.bubble_point))
    assert above < 2, "sınaq deck-i məhz bu halı yaratmalıdır"
    assert any("doymamış sıxılma" in message for message in records), records


def test_measured_conversions_match_the_deck_numbers():
    """SPE1-in öz rəqəmləri — çevirmə əmsalları kilidlənir."""
    assert RS_INITIAL == pytest.approx(226.1971, rel=1e-5)
    assert PB_INITIAL == pytest.approx(276.8038, rel=1e-5)
    assert convert(0.0093, "rb/Mscf", "m3/sm3", "gas_fvf") == pytest.approx(
        5.2216e-5, rel=1e-4)


def test_deck_branches_give_a_consistent_compressibility():
    """G3 üçün ÖLÇMƏ: iki Rs qolunun doymamış sıxılması praktik olaraq eynidir."""
    deck = read_deck_pvt(_deck())
    values = []
    for branch in deck.oil.branches:
        if not branch.has_undersaturated:
            continue
        bo = branch.formation_volume_factor
        values.append(float(np.log(bo[0] / bo[1])
                            / (branch.pressure[1] - branch.pressure[0])))
    assert len(values) == 2
    assert abs(values[0] - values[1]) / values[0] < 0.01, values
