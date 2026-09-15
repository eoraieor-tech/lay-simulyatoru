"""G4 — qaz-neft SCAL CƏDVƏLİ (Eclipse `SGOF`) və onun oxuyucusu.

SPE1CASE2 deck-i qaz əyrisini SGOF cədvəli ilə verir; bizdə isə qaz əyrisi
yalnız Corey düsturundan gəlirdi (`GasCoreyParameters`).

MÜQAVİLƏ: cədvəl sinfi `GasCoreyParameters`-in YERİNƏ keçir —
`StoneRelativePermeabilityProvider` ondan yalnız `krg`, `krog`, törəmələri,
`sgc`/`sorg`/`krg_end` və `validate(swc)` istəyir.

YOL BOYU TAPILAN SƏSSİZ SƏHV (testlə kilidlənib): üç fazalı mühərrik
qurularkən su-neft üçün HƏMİŞƏ Corey işlədilirdi — modeldə SWOF cədvəli
olsa belə, o, səssizcə atılırdı.
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder, five_spot
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.scal import CoreyParameters, GasCoreyParameters
from imex2d.domain.scal_tables import (GasSaturationTable, GasSaturationTableSet,
                                       SaturationTable, SaturationTableSet)
from imex2d.io.scal_io import ScalFormatError, read_sgof
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.scal_tables_provider import TableRelativePermeabilityProvider
from imex2d.simulation.stone_relperm import StoneRelativePermeabilityProvider

WATER_OIL = CoreyParameters()
GAS = GasCoreyParameters()

SGOF_DECK = """-- SPE1 üslubunda qaz-neft cədvəli
SGOF
-- Sg      krg       krog      Pcog
   0.00    0.000     1.000     0.0
   0.10    0.020     0.600     0.0
   0.20    0.070     0.330     0.0
   0.40    0.250     0.100     0.0
   0.60    0.560     0.000     0.0
   0.80    1.000     0.000     0.0
/
"""


def _write(text: str, suffix=".DATA") -> str:
    handle, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(handle, "w", encoding="utf-8") as file:
        file.write(text)
    return path


def _gas_table(points=41) -> GasSaturationTable:
    return GasSaturationTable.from_corey(GAS, WATER_OIL.swc, WATER_OIL.kro_end,
                                         points=points)


def _gas_set(points=41) -> GasSaturationTableSet:
    tables = GasSaturationTableSet()
    tables.add(1, _gas_table(points))
    return tables


def _water_set() -> SaturationTableSet:
    tables = SaturationTableSet()
    tables.add(1, SaturationTable.from_corey(WATER_OIL, points=41))
    return tables


def _model(gas_tables=None, water_tables=None):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=4, ny=4, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    return ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=WATER_OIL, gas_scal=GAS,
        scal_tables=water_tables, gas_scal_tables=gas_tables,
        pvt_table=build_pvt_table(pressure_min=1.0, pressure_max=400.0,
                                  n_points=40, bubble_point_bar=240.0,
                                  include_gas=True),
        name="G4 sınağı")


def _engine(model):
    service = ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)
    return service.create_engine(model, SimulationConfig(end_time=10.0))


# ═══════════════════════ cədvəl sinfi ════════════════════════════════

def test_table_interpolates_between_rows():
    table = GasSaturationTable(sg=[0.0, 0.5], krg=[0.0, 1.0], krog=[1.0, 0.0])
    assert table.interpolate_krg(0.25) == pytest.approx(0.5)
    assert table.interpolate_krog(0.25) == pytest.approx(0.5)


def test_table_endpoints_come_from_the_columns():
    table = GasSaturationTable(sg=[0.0, 0.05, 0.4, 0.8],
                               krg=[0.0, 0.0, 0.3, 0.9],
                               krog=[1.0, 0.8, 0.1, 0.0])
    assert table.sgc == pytest.approx(0.05), "krg sıfırdan çıxan son nöqtə"
    assert table.krg_end == pytest.approx(0.9)
    assert table.krog_end == pytest.approx(1.0)


def test_slope_is_the_exact_piecewise_derivative():
    """Jakobian qalıqla uyğun qalsın deyə `np.gradient` DEYİL, interval meyli."""
    table = GasSaturationTable(sg=[0.0, 0.2, 0.6], krg=[0.0, 0.1, 0.5],
                               krog=[1.0, 0.7, 0.0])
    assert table.slope(table.krg, 0.1) == pytest.approx(0.5)
    assert table.slope(table.krg, 0.4) == pytest.approx(1.0)
    assert table.slope(table.krog, 0.4) == pytest.approx(-1.75)


def test_slope_is_zero_outside_the_table():
    table = GasSaturationTable(sg=[0.1, 0.5], krg=[0.0, 1.0], krog=[1.0, 0.0])
    assert table.slope(table.krg, 0.05) == pytest.approx(0.0)
    assert table.slope(table.krg, 0.9) == pytest.approx(0.0)


def test_table_reproduces_corey_at_its_own_nodes():
    """Cədvəl Corey-dən qurulubsa, ÖZ düyünlərində eyni qiyməti verməlidir."""
    table = _gas_table(points=21)
    nodes = table.sg
    assert np.allclose(table.interpolate_krg(nodes),
                       GAS.krg(nodes, WATER_OIL.swc), atol=1e-12)
    assert np.allclose(table.interpolate_krog(nodes),
                       GAS.krog(nodes, WATER_OIL.swc, WATER_OIL.kro_end), atol=1e-12)


def test_pcog_column_is_kept():
    """SGOF-un 4-cü sütunu itirilmir (Pcog-a bağlanması ⏳ ayrı işdir)."""
    table = GasSaturationTable(sg=[0.0, 0.5], krg=[0.0, 1.0], krog=[1.0, 0.0],
                               pcog=[0.0, 1.5])
    assert table.has_capillary
    assert table.interpolate_pcog(0.25) == pytest.approx(0.75)


# ═══════════════════════ yoxlama (validate) ══════════════════════════

def test_non_monotone_krg_is_rejected():
    table = GasSaturationTable(sg=[0.0, 0.3, 0.6], krg=[0.0, 0.5, 0.2],
                               krog=[1.0, 0.5, 0.0])
    assert any("krg azalır" in message for message in table.validate())


def test_increasing_krog_is_rejected():
    table = GasSaturationTable(sg=[0.0, 0.3], krg=[0.0, 0.5], krog=[0.2, 0.6])
    assert any("krog artır" in message for message in table.validate())


def test_nan_is_caught_explicitly():
    table = GasSaturationTable(sg=[0.0, 0.3], krg=[0.0, float("nan")],
                               krog=[1.0, 0.0])
    assert any("NaN" in message for message in table.validate())


def test_unsorted_saturation_is_rejected():
    table = GasSaturationTable(sg=[0.3, 0.1], krg=[0.0, 0.5], krog=[1.0, 0.0])
    assert any("artan sıralı" in message for message in table.validate())


# ═══════════════════════ SGOF oxuyucusu ══════════════════════════════

def test_read_sgof_parses_columns_and_comments():
    tables = read_sgof(_write(SGOF_DECK))
    table = tables.get()
    assert len(tables) == 1
    assert table.sg.size == 6
    assert table.krg[-1] == pytest.approx(1.0)
    assert table.krog[0] == pytest.approx(1.0)
    assert table.has_capillary is False, "sıfır Pcog sütunu 'var' sayılmamalıdır"


def test_read_sgof_handles_two_regions():
    deck = SGOF_DECK + """   0.00    0.000     1.000     0.0
   0.50    0.400     0.000     0.0
/
"""
    tables = read_sgof(_write(deck))
    assert len(tables) == 2
    assert tables.get(2).krg_end == pytest.approx(0.4)


def test_read_sgof_rejects_a_file_without_the_keyword():
    with pytest.raises(ScalFormatError, match="SGOF"):
        read_sgof(_write("SWOF\n 0.2 0.0 0.8 0.0\n/\n"))


def test_read_sgof_rejects_an_invalid_table():
    bad = """SGOF
   0.00    0.000     1.000
   0.30    0.500     0.200
   0.60    0.200     0.000
/
"""
    with pytest.raises(ScalFormatError, match="krg"):
        read_sgof(_write(bad))


def test_read_sgof_stops_at_the_next_keyword():
    deck = SGOF_DECK + "\nPVDG\n  14.7 200.0 0.0125 /\n"
    tables = read_sgof(_write(deck))
    assert tables.get().sg.size == 6


# ═══════════════════════ Stone provider ilə ══════════════════════════

def test_stone_accepts_the_table_in_place_of_corey():
    provider = StoneRelativePermeabilityProvider(
        CoreyRelativePermeabilityAdapter(WATER_OIL), _gas_set(),
        WATER_OIL.swc, WATER_OIL.kro_end)
    sg = np.array([0.1, 0.25])
    assert np.allclose(provider.krg(sg), _gas_set().krg(sg))
    assert provider.gas_saturation_limits()[0] == pytest.approx(_gas_set().sgc)


def test_table_and_corey_agree_within_interpolation_error():
    """Eyni Corey parametrlərindən qurulan cədvəl ona yaxın olmalıdır."""
    tables = _gas_set(points=101)
    sg = np.linspace(0.06, 0.5, 25)
    assert np.allclose(tables.krg(sg, WATER_OIL.swc),
                       GAS.krg(sg, WATER_OIL.swc), atol=2e-3)
    assert np.allclose(tables.krog(sg, WATER_OIL.swc, WATER_OIL.kro_end),
                       GAS.krog(sg, WATER_OIL.swc, WATER_OIL.kro_end), atol=2e-3)


# ═══════════════════════ mühərriklə bağlantı ═════════════════════════

def test_engine_uses_the_gas_table_when_given():
    engine = _engine(_model(gas_tables=_gas_set()))
    assert isinstance(engine.relperm.gas, GasSaturationTableSet)


def test_engine_falls_back_to_corey_without_a_table():
    engine = _engine(_model())
    assert isinstance(engine.relperm.gas, GasCoreyParameters)


def test_engine_no_longer_drops_the_water_oil_table():
    """SƏSSİZ SƏHVİN QARŞISI: üç fazalı yol SWOF cədvəlini ATIRDI."""
    engine = _engine(_model(gas_tables=_gas_set(), water_tables=_water_set()))
    assert isinstance(engine.relperm.water_oil, TableRelativePermeabilityProvider)


def test_run_with_gas_table_converges():
    engine = _engine(_model(gas_tables=_gas_set()))
    result = engine.run()
    assert result.converged, result.message


def test_scal_source_panel_offers_sgof():
    pytest.importorskip("PyQt5.QtWidgets")
    import inspect
    from imex2d.ui import panels
    source = inspect.getsource(panels.ScalSourcePanel)
    assert "SGOF" in source and "gas_tables_enabled" in source
