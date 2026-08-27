"""OPM Flow nəticələrinin idxalı (strateji dönüş — bax jurnal qeydi)."""

import numpy as np
import pytest

resdata = pytest.importorskip("resdata", reason="resdata quraşdırılmayıb")

from resdata.grid import Grid
from resdata.resfile import ResdataKW, openFortIO
from resdata import ResDataType

from imex2d.io.opm_import import (OpmImportError, build_display_model,
                                  load_opm_case)


def _write_synthetic_case(tmp_path, nx=4, ny=3, nz=2, steps=3,
                          include_gas=True):
    root = str(tmp_path / "CASE")
    grid = Grid.create_rectangular((nx, ny, nz), (25.0, 25.0, 10.0))
    grid.save_EGRID(root + ".EGRID")

    ncell = nx * ny * nz
    with openFortIO(root + ".UNRST", mode=2) as handle:
        for step in range(steps):
            doub = ResdataKW("DOUBHEAD", 1, ResDataType.RD_DOUBLE)
            doub[0] = float(step * 10)
            doub.fwrite(handle)

            pressure = ResdataKW("PRESSURE", ncell, ResDataType.RD_FLOAT)
            water = ResdataKW("SWAT", ncell, ResDataType.RD_FLOAT)
            for i in range(ncell):
                pressure[i] = 200.0 - step * 2.0 + i * 0.1
                water[i] = 0.30 + step * 0.02
            pressure.fwrite(handle)
            water.fwrite(handle)

            if include_gas:
                gas = ResdataKW("SGAS", ncell, ResDataType.RD_FLOAT)
                for i in range(ncell):
                    gas[i] = 0.05 + step * 0.01
                gas.fwrite(handle)
    return root


# ── uğurlu yükləmə ──────────────────────────────────────────────────
def test_grid_geometry_is_read_correctly(tmp_path):
    root = _write_synthetic_case(tmp_path)
    case = load_opm_case(root, name="test")
    assert case.geometry.nx == 4
    assert case.geometry.ny == 3
    assert case.geometry.nz == 2
    assert abs(case.geometry.dx - 25.0) < 1e-6
    assert abs(case.geometry.dy - 25.0) < 1e-6
    assert abs(case.geometry.dz - 10.0) < 1e-6


def test_all_report_steps_are_read(tmp_path):
    root = _write_synthetic_case(tmp_path, steps=3)
    case = load_opm_case(root)
    assert len(case.snapshots) == 3


def test_snapshot_times_are_read_from_doubhead(tmp_path):
    root = _write_synthetic_case(tmp_path, steps=3)
    case = load_opm_case(root)
    assert [s.time for s in case.snapshots] == [0.0, 10.0, 20.0]


def test_pressure_and_saturation_values_are_correct(tmp_path):
    root = _write_synthetic_case(tmp_path, nx=4, ny=3, nz=2)
    case = load_opm_case(root)
    first = case.snapshots[0]
    assert abs(first.pressure[0] - 200.0) < 1e-3
    assert abs(first.water_saturation[0] - 0.30) < 1e-3
    assert abs(first.gas_saturation[0] - 0.05) < 1e-3


def test_gas_values_evolve_across_steps(tmp_path):
    root = _write_synthetic_case(tmp_path, steps=3)
    case = load_opm_case(root)
    gas_at_step = [s.gas_saturation[0] for s in case.snapshots]
    assert gas_at_step == pytest.approx([0.05, 0.06, 0.07])


# ── qazsız (iki fazalı) OPM halı ────────────────────────────────────
def test_missing_sgas_produces_a_clear_warning_not_a_crash(tmp_path):
    root = _write_synthetic_case(tmp_path, include_gas=False)
    case = load_opm_case(root)
    assert case.snapshots[0].gas_saturation is None
    assert any("SGAS" in w for w in case.warnings)


# ── xəta hallari ─────────────────────────────────────────────────────
def test_missing_files_raise_a_clear_opm_import_error(tmp_path):
    with pytest.raises(OpmImportError):
        load_opm_case(str(tmp_path / "NOEXIST"))


def test_error_message_mentions_the_missing_file(tmp_path):
    try:
        load_opm_case(str(tmp_path / "NOEXIST"))
        assert False, "istisna gözlənilirdi"
    except OpmImportError as error:
        assert "EGRID" in str(error) or "NOEXIST" in str(error)


# ── öz renderer-imizə uyğunluq ──────────────────────────────────────
def test_build_display_model_produces_a_compatible_reservoir_model(tmp_path):
    root = _write_synthetic_case(tmp_path)
    case = load_opm_case(root)
    model = build_display_model(case)
    assert model.grid.shape == (2, 3, 4)          # (nz, ny, nx)
    assert model.ncell == 24


def test_display_model_renders_without_crashing_through_our_own_renderer(tmp_path):
    """Ən vacib inteqrasiya testi: OPM məlumatı bizim öz VolumeRenderer-imizlə
    çəkilə bilirmi (OPM-in öz görüntüləyicisini əvəz etmək məqsədi)."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    from imex2d.rendering.volume import VolumeRenderer

    root = _write_synthetic_case(tmp_path)
    case = load_opm_case(root)
    model = build_display_model(case)

    figure = Figure()
    axes = figure.add_subplot(projection="3d")
    VolumeRenderer().draw(axes, figure, model,
                          case.snapshots[-1].pressure,
                          show_wells=False, show_faults=False)
    assert len(axes.collections) > 0


def test_display_model_can_render_saturation_too(tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    from imex2d.rendering.volume import VolumeRenderer

    root = _write_synthetic_case(tmp_path)
    case = load_opm_case(root)
    model = build_display_model(case)

    figure = Figure()
    axes = figure.add_subplot(projection="3d")
    VolumeRenderer().draw(axes, figure, model,
                          case.snapshots[-1].water_saturation,
                          show_wells=False, show_faults=False)
    assert len(axes.collections) > 0


# ── qeyri-aktiv hüceyrələr üzrə xəbərdarlıq ─────────────────────────
def test_inactive_cells_produce_a_warning(tmp_path):
    root = str(tmp_path / "CASE2")
    nx, ny, nz = 3, 3, 1
    actnum = [1] * (nx * ny)
    actnum[0] = 0
    grid = Grid.create_rectangular((nx, ny, nz), (25.0, 25.0, 10.0),
                                   actnum=actnum)
    grid.save_EGRID(root + ".EGRID")
    ncell = nx * ny * nz
    with openFortIO(root + ".UNRST", mode=2) as handle:
        doub = ResdataKW("DOUBHEAD", 1, ResDataType.RD_DOUBLE)
        doub[0] = 0.0
        doub.fwrite(handle)
        pressure = ResdataKW("PRESSURE", grid.get_num_active(),
                             ResDataType.RD_FLOAT)
        for i in range(grid.get_num_active()):
            pressure[i] = 200.0
        pressure.fwrite(handle)
    case = load_opm_case(root)
    assert any("aktiv" in w for w in case.warnings)
