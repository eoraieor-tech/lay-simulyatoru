"""Phase 1 (giriş boru xətti) — PVT/SCAL üçün açıq vahid metadatası.

`PVTTable.from_values()` (bax `domain/pvt.py`) idxal sərhədidir: xarici
vahiddə (FIELD: psi/cP/scf-stb) verilmiş dəyərləri mühərrik vahidinə
(bar/cP/sm3sm3) çevirir. Bo/Bw ÖLÇÜSÜZ nisbətdir — ÇEVRİLMİR (bax
`from_values` docstring-i, UNITS.md).
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder, five_spot
from imex2d.domain import unit_conversions as uc
from imex2d.domain.pvt import PVTTable
from imex2d.domain.scal_tables import SaturationTable
from imex2d.domain.wells import ControlMode, Perforation, Well, WellControl, WellType
from imex2d.simulation.pvt.correlations import build_pvt_table


def _field_values():
    pressure_psi = np.array([100.0, 1500.0, 3000.0, 4500.0])
    return dict(
        pressure=pressure_psi,
        oil_fvf=np.array([1.05, 1.20, 1.30, 1.25]),
        oil_viscosity=np.array([5.0, 2.0, 1.2, 1.5]),          # cP
        solution_gor=np.array([50.0, 400.0, 800.0, 800.0]),    # scf/stb
        water_fvf=np.full(4, 1.01),
        water_viscosity=np.full(4, 0.6),                        # cP
        bubble_point=3000.0,
        rock_compressibility=3.1e-6,                             # 1/psi
    )


# ── test #7/#8: PVT FIELD vs METRIC ──────────────────────────────────────
def test_pvt_imported_in_field_units_converts_pressure_to_engine_bar():
    values = _field_values()
    table = PVTTable.from_values(**values, pressure_unit="psi", viscosity_unit="cP",
                                 solution_gor_unit="scf/stb",
                                 rock_compressibility_unit="psi")
    assert np.allclose(table.pressure, uc.psi_to_bar(values["pressure"]))
    assert table.validate() == []


def test_pvt_imported_in_field_units_does_not_convert_dimensionless_fvf():
    """Bo/Bw ÖLÇÜSÜZDÜR — FIELD vahidində idxal edilsə də DƏYİŞMİR."""
    values = _field_values()
    table = PVTTable.from_values(**values, pressure_unit="psi",
                                 solution_gor_unit="scf/stb",
                                 rock_compressibility_unit="psi")
    assert np.array_equal(table.oil_fvf, values["oil_fvf"])
    assert np.array_equal(table.water_fvf, values["water_fvf"])


def test_pvt_imported_in_field_units_converts_solution_gor_to_sm3sm3():
    values = _field_values()
    table = PVTTable.from_values(**values, pressure_unit="psi",
                                 solution_gor_unit="scf/stb",
                                 rock_compressibility_unit="psi")
    expected = values["solution_gor"] / 5.61458
    assert np.allclose(table.solution_gor, expected, rtol=1e-6)


def test_pvt_metric_and_field_inputs_produce_the_same_physical_table():
    """Test #8: eyni fiziki PVT-ni METRIC (bar/sm3sm3) və FIELD (psi/
    scf-stb) kimi ifadə edib eyni `PVTTable`-ı almalıyıq (dəyirmi-
    səyahət tolerantlığı daxilində)."""
    field = _field_values()
    table_field = PVTTable.from_values(**field, pressure_unit="psi",
                                       solution_gor_unit="scf/stb",
                                       rock_compressibility_unit="psi")

    metric = dict(field)
    metric["pressure"] = uc.convert(field["pressure"], "psi", "bar", "pressure")
    metric["solution_gor"] = uc.convert(field["solution_gor"], "scf/stb", "sm3/sm3",
                                        "solution_gor")
    metric["bubble_point"] = uc.psi_to_bar(field["bubble_point"])
    metric["rock_compressibility"] = uc.convert_compressibility(
        field["rock_compressibility"], "psi", "bar")
    table_metric = PVTTable.from_values(**metric)   # bütün defolt (bar/cP/sm3sm3)

    assert np.allclose(table_field.pressure, table_metric.pressure, rtol=1e-9)
    assert np.allclose(table_field.solution_gor, table_metric.solution_gor, rtol=1e-9)
    assert table_field.bubble_point == pytest.approx(table_metric.bubble_point, rel=1e-9)
    assert table_field.rock_compressibility == pytest.approx(
        table_metric.rock_compressibility, rel=1e-9)


def test_pvt_from_values_default_units_reproduce_plain_constructor():
    """Vahid arqumentləri DƏYİŞDİRİLMƏSƏ, `from_values` == birbaşa
    `PVTTable(...)` (geriyə uyğunluq)."""
    values = _field_values()
    values["pressure"] = np.array([1.0, 100.0, 200.0, 300.0])   # artıq bar
    values["rock_compressibility"] = 4.5e-5                      # artıq 1/bar
    values["solution_gor"] = np.array([10.0, 50.0, 80.0, 80.0])  # artıq sm3/sm3
    via_from_values = PVTTable.from_values(**values)
    direct = PVTTable(pressure=values["pressure"], oil_fvf=values["oil_fvf"],
                      oil_viscosity=values["oil_viscosity"],
                      solution_gor=values["solution_gor"],
                      water_fvf=values["water_fvf"],
                      water_viscosity=values["water_viscosity"],
                      bubble_point=values["bubble_point"],
                      rock_compressibility=values["rock_compressibility"])
    assert np.array_equal(via_from_values.pressure, direct.pressure)
    assert via_from_values.rock_compressibility == direct.rock_compressibility


# ── temperatur vahidi (PVT cədvəl generatoru) ────────────────────────────
def test_build_pvt_table_temperature_unit_matches_celsius_default():
    table_c = build_pvt_table(temperature_c=70.0)
    table_f = build_pvt_table(temperature_c=uc.convert_temperature(70.0, "C", "F"),
                              temperature_unit="F")
    assert np.allclose(table_c.oil_viscosity, table_f.oil_viscosity, rtol=1e-12)
    assert np.allclose(table_c.oil_fvf, table_f.oil_fvf, rtol=1e-12)


# ── test #9: SCAL Pc vahidi ──────────────────────────────────────────────
# ── extrapolyasiya darvazası: model qurularkən BİR DƏFƏ (Newton həlqəsində DEYİL) ─
def _small_model(pvt_table=None, bhp_target=200.0):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=3, ny=3, dx=20.0, dy=20.0, dz=10.0, porosity=0.2, permx_base=100.0)
    wells = [Well("PROD", WellType.PRODUCER, WellControl(ControlMode.BHP, bhp_target),
                  [Perforation(1, 1, 0)])]
    return ReservoirModelBuilder().build(geology, wells, pvt_table=pvt_table)


def test_pvt_range_gate_passes_silently_when_bhp_within_table():
    table = build_pvt_table(pressure_min=1.0, pressure_max=400.0)
    model = _small_model(pvt_table=table, bhp_target=200.0)
    report = model.diagnose()
    assert not any("PVT" in d.source for d in report.errors)
    assert not any("PVT" in d.source for d in report.warnings)


def test_pvt_range_gate_warns_on_mild_extrapolation_but_does_not_block():
    table = build_pvt_table(pressure_min=1.0, pressure_max=400.0)
    model = _small_model(pvt_table=table, bhp_target=420.0)   # yüngül kənar
    report = model.diagnose()
    assert not any("PVT" in d.source for d in report.errors)
    assert any("PVT" in d.source for d in report.warnings)
    assert model.validate() == []   # yüngül kənarlaşma BLOKLAMIR


def test_pvt_range_gate_hard_errors_on_severe_extrapolation():
    """BHP=3000 (məs. psi sanılıb bar kimi verilmiş dəyər) 1-400 bar-lıq
    PVT cədvəlindən HƏDDİNDƏN ARTIQ kənardadır — bu, `model.validate()`-i
    BLOKLAMALIDIR (vahid qarışıqlığı əlaməti)."""
    table = build_pvt_table(pressure_min=1.0, pressure_max=400.0)
    model = _small_model(pvt_table=table, bhp_target=3000.0)
    report = model.diagnose()
    assert any("PVT" in d.source for d in report.errors)
    assert model.validate() != []


def test_scal_table_check_query_range_uses_dimensionless_saturation():
    """Sw/kr ÖLÇÜSÜZDÜR — çevrilmir; Pc-nin öz vahidi (adətən bar,
    mühərrikin təzyiq vahidi) `check_query_range`-ə TƏSİR ETMİR, çünki
    sorğu Sw üzərindədir, Pc üzərində deyil."""
    table = SaturationTable(sw=np.array([0.2, 0.5, 0.8]),
                            krw=np.array([0.0, 0.1, 0.4]),
                            kro=np.array([0.6, 0.2, 0.0]),
                            pc=np.array([0.5, 0.1, 0.0]))
    assert table.check_query_range([0.2, 0.5, 0.8]) == []
    assert len(table.check_query_range([0.05, 0.95])) == 2
