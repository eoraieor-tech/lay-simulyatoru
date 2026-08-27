"""Eclipse deck I/O (B5): GRDECL oxuma və `.DATA` yazma."""

import os
import tempfile

import numpy as np

from helpers import default_scal, make_service, short_config
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.domain.diagnostics import DiagnosticReport
from imex2d.io.eclipse_export import EclipseDeckWriter, _compress
from imex2d.io.grdecl import GrdeclError, read_grdecl
from imex2d.io.grdecl_import import GrdeclImporter
from imex2d.simulation.pvt.correlations import build_pvt_table


def _write(text: str, suffix=".GRDECL") -> str:
    handle, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(handle, "w", encoding="utf-8") as file:
        file.write(text)
    return path


def _minimal_deck(nx=3, ny=2, nz=2) -> str:
    n = nx * ny * nz
    return f"""-- nümunə
SPECGRID
  {nx} {ny} {nz} 1 F /

DX
  {n}*50 /
DY
  {n}*40 /
DZ
  {n}*10 /
TOPS
  {nx * ny}*2000 /
PORO
  {n}*0.2 /
PERMX
  {n}*150 /
"""


def _model(nx=10, ny=8, nz=3, **kwargs):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=30.0, dz=6.0, porosity=0.21,
        permx_base=180.0, nz=nz, top_depth=2050.0, **kwargs)
    return ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                         scal=default_scal())


# ── oxuma ─────────────────────────────────────────────────────────────
def test_reads_dimensions_and_arrays():
    path = _write(_minimal_deck())
    try:
        deck = read_grdecl(path)
        assert deck.dimensions == (3, 2, 2)
        assert deck.ncell == 12
        for keyword in ("DX", "DY", "DZ", "PORO", "PERMX"):
            assert deck.has(keyword), keyword
        assert np.allclose(deck.get("PORO"), 0.2)
    finally:
        os.unlink(path)


def test_repeat_syntax_is_expanded():
    path = _write("""SPECGRID
  2 2 1 1 F /
PORO
  2*0.2 0.3 0.25 /
PERMX
  100 3*200 /
""")
    try:
        deck = read_grdecl(path)
        assert list(deck.get("PORO")) == [0.2, 0.2, 0.3, 0.25]
        assert list(deck.get("PERMX")) == [100.0, 200.0, 200.0, 200.0]
    finally:
        os.unlink(path)


def test_section_headers_do_not_swallow_the_next_array():
    """`GRID` kimi başlıqlar `/` tələb etmir.

    Onları adi açar söz saymaq incə, lakin ağır səhvdir: oxucu növbəti
    `/` işarəsinə qədər hər şeyi udur və ilk massiv (adətən DX) itir.
    Model səssizcə defolt ölçülərlə qurulur — nəticələr yanlış olur.
    """
    path = _write("""RUNSPEC
DIMENS
  2 2 1 /
OIL
WATER
METRIC

GRID

DX
  4*50 /
DY
  4*50 /
DZ
  4*10 /
PORO
  4*0.2 /
PERMX
  4*100 /
INIT
""")
    try:
        deck = read_grdecl(path)
        assert deck.has("DX"), "DX bölmə başlığı tərəfindən udulub"
        assert np.allclose(deck.get("DX"), 50.0)
    finally:
        os.unlink(path)


def test_comments_are_ignored():
    path = _write("""-- başlıq
SPECGRID
  2 2 1 1 F /   -- ölçülər
PORO
  4*0.2 /       -- məsaməlilik
PERMX
  4*100 /
""")
    try:
        assert read_grdecl(path).dimensions == (2, 2, 1)
    finally:
        os.unlink(path)


def test_missing_dimensions_is_rejected():
    path = _write("PORO\n  4*0.2 /\n")
    try:
        read_grdecl(path)
    except GrdeclError:
        return
    finally:
        os.unlink(path)
    raise AssertionError("Ölçüsüz fayl qəbul edildi")


def test_wrong_array_length_is_rejected():
    path = _write("SPECGRID\n  4 4 1 1 F /\nPORO\n  5*0.2 /\n")
    try:
        read_grdecl(path)
    except GrdeclError as error:
        assert "PORO" in str(error)
        return
    finally:
        os.unlink(path)
    raise AssertionError("Yanlış uzunluqlu massiv qəbul edildi")


def test_include_produces_a_warning():
    path = _write("""SPECGRID
  2 2 1 1 F /
INCLUDE
  'başqa.grdecl' /
PORO
  4*0.2 /
PERMX
  4*100 /
""")
    try:
        report = DiagnosticReport()
        read_grdecl(path, report)
        assert any("INCLUDE" in w.message for w in report.warnings)
    finally:
        os.unlink(path)


# ── model qurulması ───────────────────────────────────────────────────
def test_builds_a_valid_geological_model():
    path = _write(_minimal_deck(nx=5, ny=4, nz=3))
    try:
        model = GrdeclImporter().build(read_grdecl(path))
        assert model.grid.nx == 5 and model.grid.nz == 3
        assert model.geometry.dx == 50.0 and model.geometry.dz == 10.0
        assert model.validate() == []
        for key in ("PORO", "PERMX", "PERMY", "PERMZ"):
            assert key in model.property_maps
    finally:
        os.unlink(path)


def test_missing_permeability_is_rejected():
    path = _write("SPECGRID\n  2 2 1 1 F /\nPORO\n  4*0.2 /\n")
    try:
        GrdeclImporter().build(read_grdecl(path))
    except GrdeclError:
        return
    finally:
        os.unlink(path)
    raise AssertionError("PERMX olmadan model quruldu")


def test_variable_cell_size_falls_back_to_the_mean_with_a_warning():
    path = _write("""SPECGRID
  2 2 1 1 F /
DX
  40 60 40 60 /
DY
  4*50 /
DZ
  4*10 /
PORO
  4*0.2 /
PERMX
  4*100 /
""")
    try:
        report = DiagnosticReport()
        model = GrdeclImporter().build(read_grdecl(path, report), report)
        assert abs(model.geometry.dx - 50.0) < 1e-9
        assert any("DX" in w.message for w in report.warnings)
    finally:
        os.unlink(path)


def test_inactive_cells_produce_a_warning():
    """ACTNUM dəstəklənmir — susmaq həcm səhvinə aparır."""
    path = _write("""SPECGRID
  2 2 1 1 F /
DX
  4*50 /
DY
  4*50 /
DZ
  4*10 /
PORO
  4*0.2 /
PERMX
  4*100 /
ACTNUM
  1 1 0 1 /
""")
    try:
        report = DiagnosticReport()
        GrdeclImporter().build(read_grdecl(path, report), report)
        assert any("ACTNUM" in w.message for w in report.warnings)
    finally:
        os.unlink(path)


def test_satnum_becomes_regions():
    path = _write("""SPECGRID
  2 2 1 1 F /
DX
  4*50 /
DY
  4*50 /
DZ
  4*10 /
PORO
  4*0.2 /
PERMX
  4*100 /
SATNUM
  1 1 2 2 /
""")
    try:
        model = GrdeclImporter().build(read_grdecl(path))
        assert sorted(model.regions.ids) == [1, 2]
    finally:
        os.unlink(path)


def test_varying_tops_become_a_dipping_surface():
    path = _write("""SPECGRID
  2 2 1 1 F /
DX
  4*50 /
DY
  4*50 /
DZ
  4*10 /
TOPS
  2000 2010 2005 2015 /
PORO
  4*0.2 /
PERMX
  4*100 /
""")
    try:
        model = GrdeclImporter().build(read_grdecl(path))
        depths = model.geometry.cell_depths()
        assert np.ptp(depths) > 0, "Maili səth qurulmadı"
    finally:
        os.unlink(path)


# ── yazma ─────────────────────────────────────────────────────────────
def test_compression_produces_repeat_tokens():
    assert _compress(np.array([0.2] * 5 + [0.3, 0.3, 0.4])) == \
        ["5*0.2", "2*0.3", "0.4"]


def test_deck_contains_all_required_sections():
    text = EclipseDeckWriter().render(_model())
    for section in ("RUNSPEC", "GRID", "PROPS", "SOLUTION", "SUMMARY",
                    "SCHEDULE", "END"):
        assert f"\n{section}" in text or text.startswith(section), section


def test_deck_lists_every_well_and_perforation():
    model = _model(nz=4)
    text = EclipseDeckWriter().render(model)
    for well in model.active_wells():
        assert f"'{well.name}'" in text
    perforations = sum(len(w.open_perforations())
                       for w in model.active_wells())
    compdat = text.split("COMPDAT")[1].split("/\n\n")[0]
    assert compdat.count("OPEN") == perforations


def test_deck_includes_pvt_table_when_present():
    model = _model()
    model.pvt_table = build_pvt_table(bubble_point_bar=150.0)
    text = EclipseDeckWriter().render(model)
    assert "PVDO" in text
    assert "PVCDO" not in text


def test_deck_uses_constant_pvt_without_a_table():
    text = EclipseDeckWriter().render(_model())
    assert "PVCDO" in text
    assert "PVDO" not in text


def test_deck_writes_capillary_pressure_into_swof():
    model = _model()
    model.capillary_parameters.entry_pressure = 0.4
    rows = EclipseDeckWriter().render(model).split("SWOF")[1].split("/")[0]
    values = [line.split() for line in rows.strip().split("\n")
              if line.strip() and not line.strip().startswith("--")]
    assert all(float(row[3]) > 0 for row in values)


def test_written_file_can_be_read_back():
    model = _model()
    handle, path = tempfile.mkstemp(suffix=".DATA")
    os.close(handle)
    try:
        EclipseDeckWriter().write(model, path)
        deck = read_grdecl(path)
        assert deck.dimensions == (model.grid.nx, model.grid.ny,
                                   model.grid.nz)
    finally:
        os.unlink(path)


# ── tam dövrə ─────────────────────────────────────────────────────────
def test_round_trip_preserves_geometry_and_properties():
    model = _model(heterogeneous=True, sigma=0.6, seed=9)
    handle, path = tempfile.mkstemp(suffix=".DATA")
    os.close(handle)
    try:
        EclipseDeckWriter().write(model, path)
        report = DiagnosticReport()
        restored = GrdeclImporter().build(read_grdecl(path, report), report)
    finally:
        os.unlink(path)

    assert restored.geometry.dx == model.geometry.dx
    assert restored.geometry.dy == model.geometry.dy
    assert restored.geometry.dz == model.geometry.dz
    assert np.allclose(restored.geometry.cell_depths(),
                       model.geometry.cell_depths(), atol=1e-3)
    assert np.allclose(restored.property_maps["PORO"].values,
                       model.rock.porosity.values, atol=1e-5)
    assert np.allclose(restored.property_maps["PERMX"].values,
                       model.rock.permx.values, rtol=1e-4)


def test_round_trip_reproduces_the_simulation_result():
    """B5-in əsl yoxlanışı: fayldan geri qurulan model eyni cavabı verir."""
    scal = default_scal()
    model = _model(heterogeneous=True, sigma=0.6, seed=9)
    handle, path = tempfile.mkstemp(suffix=".DATA")
    os.close(handle)
    try:
        EclipseDeckWriter().write(model, path)
        restored_geology = GrdeclImporter().build(read_grdecl(path))
    finally:
        os.unlink(path)

    rebuilt = ReservoirModelBuilder().build(
        restored_geology, five_spot(restored_geology.grid), scal=scal)

    original = make_service(scal).run(model, short_config(end_time=200.0))
    repeated = make_service(scal).run(rebuilt, short_config(end_time=200.0))

    assert abs(original.ooip - repeated.ooip) / original.ooip < 1e-4
    assert abs(original.final_recovery_factor
               - repeated.final_recovery_factor) < 0.01
