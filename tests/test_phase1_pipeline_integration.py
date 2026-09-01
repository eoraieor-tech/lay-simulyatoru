"""Phase 1 (giriş boru xətti) — tapşırıqda sadalanan 14 inteqrasiya testi.

Bu fayl HƏR BİRİNİ AÇIQ, tək-tək təmsil edir (bəziləri digər test
fayllarında daha ətraflı sınanıb — bura referans üçün konsolidasiya
edir, tam əhatə üçün deyil).
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest
from helpers import default_scal, make_service, one_dimensional_model, short_config

from imex2d.domain import unit_conversions as uc
from imex2d.domain.pvt import PVTTable
from imex2d.geology.well_data_io import read_well_csv
from imex2d.io.grdecl import read_grdecl
from imex2d.io.grdecl_import import GrdeclImporter
from imex2d.domain.diagnostics import DiagnosticReport
from imex2d.domain.properties import PropertyMap


# 1. 3000 psi -> mühərrik bar
def test_1_psi_to_engine_bar():
    assert uc.psi_to_bar(3000.0) == pytest.approx(206.8427, rel=1e-6)


# 2. 100 ft -> mühərrik m
def test_2_ft_to_engine_m():
    assert uc.ft_to_m(100.0) == pytest.approx(30.48, rel=1e-9)


# 3. 1000 mD -> mühərrik mD (no-op, artıq mühərrik vahidi)
def test_3_md_to_engine_md_is_noop():
    assert uc.to_engine_units(1000.0, "mD", "permeability") == 1000.0


# 4. 1 Darsi -> mühərrik mD
def test_4_darcy_to_engine_md():
    assert uc.to_engine_units(1.0, "D", "permeability") == pytest.approx(1000.0, rel=1e-12)


# 5. 1 cP -> mühərrik cP (no-op)
def test_5_cp_to_engine_cp_is_noop():
    assert uc.to_engine_units(1.0, "cP", "viscosity") == 1.0


# 6. 1000 stb/day -> mühərrik m3/day
def test_6_stb_per_day_to_engine_m3_per_day():
    assert uc.stb_per_day_to_m3_per_day(1000.0) == pytest.approx(158.987294928, rel=1e-9)


# 7. PVT FIELD vahidində idxal -> mühərrik vahidi
def test_7_pvt_imported_in_field_units_converts_to_engine():
    table = PVTTable.from_values(
        pressure=[100.0, 4000.0], oil_fvf=[1.05, 1.30], oil_viscosity=[5.0, 1.2],
        solution_gor=[50.0, 800.0], water_fvf=[1.01, 1.01], water_viscosity=[0.6, 0.6],
        pressure_unit="psi", solution_gor_unit="scf/stb")
    assert table.pressure[0] == pytest.approx(uc.psi_to_bar(100.0), rel=1e-9)
    assert table.validate() == []


# 8. PVT METRIC vahidində idxal -> EYNİ fiziki nəticə
def test_8_pvt_metric_and_field_are_physically_equivalent():
    field = PVTTable.from_values(
        pressure=[100.0, 4000.0], oil_fvf=[1.05, 1.30], oil_viscosity=[5.0, 1.2],
        solution_gor=[50.0, 800.0], water_fvf=[1.01, 1.01], water_viscosity=[0.6, 0.6],
        pressure_unit="psi", solution_gor_unit="scf/stb")
    metric = PVTTable.from_values(
        pressure=uc.convert(np.array([100.0, 4000.0]), "psi", "bar", "pressure"),
        oil_fvf=[1.05, 1.30], oil_viscosity=[5.0, 1.2],
        solution_gor=uc.convert(np.array([50.0, 800.0]), "scf/stb", "sm3/sm3", "solution_gor"),
        water_fvf=[1.01, 1.01], water_viscosity=[0.6, 0.6])
    assert np.allclose(field.pressure, metric.pressure, rtol=1e-9)
    assert np.allclose(field.solution_gor, metric.solution_gor, rtol=1e-9)


# 9. SCAL təzyiq (Pc) çevirməsi — Pc mühərrik təzyiq vahidindədir (bar);
# Sw/kr ÖLÇÜSÜZDÜR, çevrilmir (bax `domain/scal_tables.py` docstring-i).
def test_9_scal_pc_is_a_pressure_quantity_saturation_is_dimensionless():
    from imex2d.domain.scal_tables import SaturationTable
    pc_psi = np.array([10.0, 5.0, 0.0])
    pc_bar = uc.convert(pc_psi, "psi", "bar", "pressure")
    assert pc_bar[0] == pytest.approx(uc.psi_to_bar(10.0), rel=1e-9)

    # Sw/kr [0,1] ÖLÇÜSÜZ nisbətdir — heç bir vahid çevirməsi tətbiq
    # edilmir, "vahid" olmadan da fiziki mənası tamdır.
    table_from_psi = SaturationTable(sw=np.array([0.2, 0.5, 0.8]),
                                     krw=np.array([0.0, 0.1, 0.4]),
                                     kro=np.array([0.6, 0.2, 0.0]),
                                     pc=pc_bar)   # Pc artıq mühərrik vahidində (bar)
    assert table_from_psi.validate() == []
    assert table_from_psi.interpolate_pc(0.5) == pytest.approx(pc_bar[1])


# 10. PropertyMap etibarsız vahid kombinasiyası rədd edilir
def test_10_property_map_rejects_invalid_unit_combination():
    with pytest.raises(ValueError):
        PropertyMap.from_array("PERMX", [100.0, 200.0], 2, "psi")


# 11. Vahid göstərilməyəndə xəbərdarlıq (CSV)
def test_11_missing_unit_produces_warning_not_silent_guess():
    handle, path = tempfile.mkstemp(suffix=".csv")
    with os.fdopen(handle, "w", encoding="utf-8") as f:
        f.write("well,x,y,PORO,PERMX\nW-1,10,20,0.2,150\nW-2,90,80,0.25,400\n")
    try:
        dataset = read_well_csv(path)
        assert dataset.samples[0].values["PERMX"] == 150.0   # dəyər DƏYİŞMİR
        assert any("PERMX" in w for w in dataset.warnings)   # amma bildirilir
    finally:
        os.unlink(path)


# 12. GRDECL FIELD sükutla METRIC olmur
def test_12_grdecl_field_never_silently_becomes_metric():
    handle, path = tempfile.mkstemp(suffix=".GRDECL")
    with os.fdopen(handle, "w", encoding="utf-8") as f:
        f.write("""RUNSPEC
DIMENS
  2 2 1 /
FIELD
GRID
DX
  4*100 /
DY
  4*100 /
DZ
  4*10 /
PORO
  4*0.2 /
PERMX
  4*150 /
""")
    try:
        report = DiagnosticReport()
        model = GrdeclImporter().build(read_grdecl(path, report), report)
        # DX=100 FT -> 30.48 m, DEYİL 100 m (sükutla METRIC sayılsaydı belə olardı)
        assert model.geometry.dx == pytest.approx(uc.ft_to_m(100.0), rel=1e-9)
        assert any("FIELD" in w.message for w in report.warnings)
    finally:
        os.unlink(path)


# 13. Mövcud (vahidsiz/legacy) giriş ƏDƏDİ olaraq DƏYİŞMƏZ qalır
def test_13_legacy_unitless_inputs_remain_numerically_unchanged():
    handle, path = tempfile.mkstemp(suffix=".csv")
    with os.fdopen(handle, "w", encoding="utf-8") as f:
        f.write("well,x,y,PORO,PERMX\nW-1,10,20,0.20,150\nW-2,90,80,0.25,400\n")
    try:
        dataset = read_well_csv(path)
        assert dataset.samples[0].x == 10.0 and dataset.samples[0].y == 20.0
        assert dataset.samples[0].values["PERMX"] == 150.0
        assert dataset.samples[1].values["PERMX"] == 400.0
    finally:
        os.unlink(path)


# 14. Tam simulyasiya: SI/field/metrik-ekvivalent giriş -> ekvivalent nəticə
def test_14_full_simulation_field_and_metric_equivalent_inputs_agree():
    scal = default_scal()

    def run(dx, permeability, injection_rate):
        model = one_dimensional_model(nx=20, dx=dx, permeability=permeability,
                                      injection_rate=injection_rate, scal=scal)
        return make_service(scal).run(model, short_config(end_time=60.0, snapshots=4))

    result_metric = run(8.0, 200.0, 60.0)

    dx_ft = uc.m_to_ft(8.0)
    rate_stbday = uc.m3_per_day_to_stb_per_day(60.0)
    result_field_roundtrip = run(uc.to_engine_units(dx_ft, "ft", "length"), 200.0,
                                 uc.to_engine_units(rate_stbday, "stb/day", "rate"))

    assert abs(result_metric.ooip - result_field_roundtrip.ooip) / result_metric.ooip < 1e-9
    assert abs(result_metric.final_recovery_factor
              - result_field_roundtrip.final_recovery_factor) < 1e-6
    assert result_metric.steps == result_field_roundtrip.steps
