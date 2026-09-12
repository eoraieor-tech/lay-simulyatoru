"""B5-a — nəticələrin CSV / JSON ixracı.

ƏSAS MÜQAVİLƏ: ixrac olunan fayl yenidən oxunanda **eyni ədədləri**
verir. Qalan testlər dörd konkret tələni kilidləyir:

  1. `nan` BOŞ xana kimi yazılır, `0` kimi YOX — `well_thp`-də quyunun
     səthə axa bilmədiyi addımlar `nan`-dır (bax `wellbore/traverse.py`);
     orada 0 bar yazmaq real ölçmə kimi oxunardı.
  2. İki fazalı nəticədə qaz sıraları BOŞDUR — onlar sütun yaratmamalıdır.
  3. `well_bhp`/`well_thp` yalnız bəzi quyularda var (BHP rejimli,
     lüləsi verilmiş istismarçılar) — yoxluq `KeyError` verməməlidir.
  4. Vahid HƏR sütun başlığındadır (`UNITS.md` prinsipi).
"""

from __future__ import annotations

import json
import math

import pytest

from helpers import default_scal
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.scal import GasCoreyParameters
from imex2d.domain.tubing import TubingGeometry
from imex2d.domain.wells import WellType
from imex2d.reporting import results_export
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.results import SimulationResult, TimeSeries
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter


def _run(with_gas=True, with_tubing=True, top_depth=1200.0, end_time=300.0):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=8, ny=8, dx=20.0, dy=20.0, dz=10.0, porosity=0.22,
        permx_base=150.0, top_depth=top_depth)
    wells = five_spot(geology.grid)
    if with_tubing:
        for well in wells:
            if well.well_type is WellType.PRODUCER:
                well.tubing = TubingGeometry(diameter=0.062, segments=20)
    model = ReservoirModelBuilder().build(
        geological_model=geology, wells=wells, scal=default_scal(),
        gas_scal=GasCoreyParameters() if with_gas else None,
        pvt_table=build_pvt_table(pressure_min=1.0, pressure_max=400.0,
                                  n_points=40, bubble_point_bar=240.0,
                                  include_gas=with_gas),
        name="ixrac sınağı")
    service = ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)
    result = service.run(model, SimulationConfig(end_time=end_time))
    assert result.converged, result.message
    return result


# ═══════════════════════ əsas müqavilə ═══════════════════════════════

def test_csv_round_trip_preserves_every_number(tmp_path):
    """İXRACIN ƏSAS ŞƏRTİ — geri oxunanda eyni ədədlər."""
    result = _run()
    path = str(tmp_path / "neticeler.csv")
    results_export.write_csv(result, path)

    back = results_export.read_csv(path)
    assert back, "fayl boş qayıtdı"

    time_column = "t [gün]"
    assert time_column in back
    assert len(back[time_column]) == len(result.series.time)
    for written, original in zip(back[time_column], result.series.time):
        assert written == pytest.approx(original, rel=0, abs=0)

    rf = back["RF [%]"]
    for written, original in zip(rf, result.series.recovery_factor):
        assert written == pytest.approx(original, rel=0, abs=0)


def test_json_round_trip_preserves_every_number(tmp_path):
    result = _run()
    path = str(tmp_path / "neticeler.json")
    results_export.write_json(result, path)

    with open(path, encoding="utf-8") as handle:
        payload = json.load(handle)

    assert payload["series"]["time"] == list(result.series.time)
    assert payload["series"]["recovery_factor"] == \
        list(result.series.recovery_factor)
    assert payload["metadata"]["converged"] is True
    assert payload["metadata"]["ooip"] == pytest.approx(result.ooip)


def test_write_picks_the_format_from_the_extension(tmp_path):
    result = _run()
    csv_path = str(tmp_path / "a.csv")
    json_path = str(tmp_path / "a.json")
    results_export.write(result, csv_path)
    results_export.write(result, json_path)

    assert results_export.read_csv(csv_path)
    with open(json_path, encoding="utf-8") as handle:
        assert "metadata" in json.load(handle)


# ═══════════════════════ tələ 1 — `nan` ══════════════════════════════

def test_nan_is_written_as_an_empty_cell_not_zero(tmp_path):
    """Quyunun axmadığı addım BOŞ qalmalıdır, `0` OLMAMALIDIR."""
    result = SimulationResult(model_name="nan sınağı")
    result.series = TimeSeries(time=[0.0, 1.0, 2.0],
                               recovery_factor=[0.0, 1.0, 2.0])
    result.well_thp["P1"] = [10.0, float("nan"), 30.0]

    path = str(tmp_path / "nan.csv")
    results_export.write_csv(result, path)

    raw = open(path, encoding=results_export.ENCODING).read().splitlines()
    assert raw[2].endswith(","), f"boş xana gözlənilirdi: {raw[2]!r}"
    assert "0.0" not in raw[2].split(",")[-1]

    back = results_export.read_csv(path)["P1: THP [bar]"]
    assert back[0] == pytest.approx(10.0)
    assert math.isnan(back[1]), "nan gedər-gələrdə qorunmalıdır"
    assert back[2] == pytest.approx(30.0)


def test_json_writes_nan_as_null(tmp_path):
    """JSON-da `NaN` standart deyil — `null` yazılmalıdır."""
    result = SimulationResult()
    result.series = TimeSeries(time=[0.0, 1.0])
    result.well_thp["P1"] = [5.0, float("nan")]

    path = str(tmp_path / "nan.json")
    results_export.write_json(result, path)

    text = open(path, encoding="utf-8").read()
    assert "NaN" not in text
    payload = json.loads(text)
    assert payload["wells"]["P1"]["well_thp"] == [5.0, None]


def test_real_run_with_unreachable_surface_exports_blanks(tmp_path):
    """Real qaçış: 5000 m-də quyu axmır → bütün THP xanaları boş."""
    result = _run(top_depth=5000.0, end_time=200.0)
    assert result.well_thp, "sınaq THP sırası tələb edir"

    path = str(tmp_path / "derin.csv")
    results_export.write_csv(result, path)
    back = results_export.read_csv(path)

    for name in result.well_thp:
        column = f"{name}: THP [bar]"
        assert all(math.isnan(v) for v in back[column])


# ═══════════════════════ tələ 2 — qazsız nəticə ══════════════════════

def test_two_phase_result_has_no_gas_columns(tmp_path):
    """Qaz sıraları boşdursa sütun YARADILMAMALIDIR."""
    result = _run(with_gas=False)
    path = str(tmp_path / "iki_fazali.csv")
    results_export.write_csv(result, path)

    headers = list(results_export.read_csv(path))
    assert not any("qaz" in h.lower() or "GOR" in h for h in headers), headers
    assert "t [gün]" in headers and "RF [%]" in headers


def test_three_phase_result_does_have_gas_columns(tmp_path):
    """NƏZARƏT: qaz aktivdirsə sütunlar MÖVCUD olmalıdır."""
    result = _run(with_gas=True)
    path = str(tmp_path / "uc_fazali.csv")
    results_export.write_csv(result, path)

    headers = list(results_export.read_csv(path))
    assert "GOR [sm³/sm³]" in headers
    assert "q_qaz [m³/gün]" in headers


# ═══════════════════════ tələ 3 — natamam quyu sıraları ══════════════

def test_wells_without_thp_do_not_break_the_export(tmp_path):
    """Lülə verilməyəndə `well_thp` boşdur — ixrac sınmamalıdır."""
    result = _run(with_tubing=False)
    assert not result.well_thp

    path = str(tmp_path / "lulesiz.csv")
    results_export.write_csv(result, path)
    headers = list(results_export.read_csv(path))

    assert not any("THP" in h for h in headers)
    assert any("q_neft" in h and ":" in h for h in headers), \
        "quyu debitləri yenə də olmalıdır"


def test_missing_well_key_is_not_an_error(tmp_path):
    """Bir quyuda BHP var, THP yox — `.get()` naxışı qorunmalıdır."""
    result = SimulationResult()
    result.series = TimeSeries(time=[0.0, 1.0])
    result.well_bhp["P1"] = [150.0, 150.0]
    result.well_oil_rate["P2"] = [3.0, 2.5]

    path = str(tmp_path / "qarishiq.csv")
    results_export.write_csv(result, path)
    headers = list(results_export.read_csv(path))

    assert "P1: BHP [bar]" in headers
    assert "P2: q_neft [m³/gün]" in headers
    assert not any("P1: THP" in h for h in headers)


def test_empty_result_does_not_crash(tmp_path):
    """Heç bir sıra yoxdursa da fayl yaranmalıdır."""
    path = str(tmp_path / "bos.csv")
    results_export.write_csv(SimulationResult(), path)
    assert open(path, encoding=results_export.ENCODING).read().strip()


def test_short_column_is_padded_with_blanks(tmp_path):
    """Qısa sütun uydurma dəyərlə DOLDURULMUR, boş qalır."""
    result = SimulationResult()
    result.series = TimeSeries(time=[0.0, 1.0, 2.0])
    result.well_bhp["P1"] = [150.0]          # qəsdən qısa

    path = str(tmp_path / "qisa.csv")
    results_export.write_csv(result, path)
    back = results_export.read_csv(path)

    assert len(back["P1: BHP [bar]"]) == 3
    assert back["P1: BHP [bar]"][0] == pytest.approx(150.0)
    assert all(math.isnan(v) for v in back["P1: BHP [bar]"][1:])


# ═══════════════════════ tələ 4 — vahidlər ═══════════════════════════

def test_every_header_carries_its_unit(tmp_path):
    """`UNITS.md` prinsipi: vahid başlıqda olmalıdır."""
    result = _run()
    path = str(tmp_path / "vahidler.csv")
    results_export.write_csv(result, path)

    for header in results_export.read_csv(path):
        assert header.endswith("]"), f"vahidsiz başlıq: {header!r}"
        unit = header[header.rindex("[") + 1:-1].strip()
        assert unit, f"boş vahid: {header!r}"
        # Vahid VERGÜLLƏ ayrılmamalıdır — vergül CSV-nin ayırıcısıdır,
        # başlığa düşsə hər başlıq dırnağa alınır.
        assert "," not in header, f"başlıqda vergül: {header!r}"


def test_csv_opens_with_a_bom_for_excel(tmp_path):
    """Windows Excel BOM-suz UTF-8-i tanımır — 'ə/ş/ğ' korlanır."""
    result = _run()
    path = str(tmp_path / "bom.csv")
    results_export.write_csv(result, path)

    assert open(path, "rb").read(3) == b"\xef\xbb\xbf"
    # BOM gedər-gələri POZMAMALIDIR
    assert "t [gün]" in results_export.read_csv(path)
