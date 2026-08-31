"""`xy_to_ij` / `depth_to_k` — metr → indeks çevrilməsi."""

import numpy as np

from imex2d.domain.geometry import CellGeometry, depth_to_k, xy_to_ij
from imex2d.domain.grid import CartesianGrid


def _flat_geometry(nx=10, ny=8, nz=5, dx=20.0, dy=25.0, dz=4.0, top=2000.0):
    grid = CartesianGrid(nx, ny, nz)
    return grid, CellGeometry(grid, dx, dy, dz, top_depth=top)


# ── xy_to_ij ─────────────────────────────────────────────────────────────
def test_xy_to_ij_at_grid_centre():
    grid, geometry = _flat_geometry(nx=10, ny=8, dx=20.0, dy=25.0)
    i, j = xy_to_ij(105.0, 110.0, geometry)          # cell (5, 4)
    assert (i, j) == (5, 4)


def test_xy_to_ij_at_lower_boundary():
    grid, geometry = _flat_geometry(nx=10, ny=8, dx=20.0, dy=25.0)
    i, j = xy_to_ij(0.0, 0.0, geometry)
    assert (i, j) == (0, 0)


def test_xy_to_ij_at_exact_upper_boundary_does_not_overflow():
    """x == x_max halında i = nx çıxmamalıdır (tapşırığın açıq qaydası)."""
    grid, geometry = _flat_geometry(nx=10, ny=8, dx=20.0, dy=25.0)
    x_max, y_max = geometry.areal_extent()
    i, j = xy_to_ij(x_max, y_max, geometry)
    assert i == grid.nx - 1
    assert j == grid.ny - 1


def test_xy_to_ij_clamps_values_beyond_grid():
    grid, geometry = _flat_geometry(nx=10, ny=8, dx=20.0, dy=25.0)
    i, j = xy_to_ij(-50.0, 100000.0, geometry)
    assert i == 0
    assert j == grid.ny - 1


def test_xy_to_ij_just_inside_each_cell():
    grid, geometry = _flat_geometry(nx=4, ny=4, dx=10.0, dy=10.0)
    for idx in range(grid.nx):
        i, _ = xy_to_ij(idx * 10.0 + 0.01, 5.0, geometry)
        assert i == idx


# ── depth_to_k: düz lay ──────────────────────────────────────────────────
def test_depth_to_k_flat_layer_first_and_last():
    grid, geometry = _flat_geometry(nx=5, ny=5, nz=5, dz=4.0, top=2000.0)
    assert depth_to_k(50.0, 50.0, 2000.5, geometry) == 0
    assert depth_to_k(50.0, 50.0, 2019.9, geometry) == 4    # son təbəqə: 2016-2020


def test_depth_to_k_flat_layer_middle():
    grid, geometry = _flat_geometry(nx=5, ny=5, nz=5, dz=4.0, top=2000.0)
    assert depth_to_k(50.0, 50.0, 2010.0, geometry) == 2     # 2008-2012 aralığı


def test_depth_to_k_returns_none_above_grid():
    grid, geometry = _flat_geometry(nx=5, ny=5, nz=5, dz=4.0, top=2000.0)
    assert depth_to_k(50.0, 50.0, 1500.0, geometry) is None


def test_depth_to_k_returns_none_below_grid():
    grid, geometry = _flat_geometry(nx=5, ny=5, nz=5, dz=4.0, top=2000.0)
    assert depth_to_k(50.0, 50.0, 5000.0, geometry) is None


def test_depth_to_k_respects_variable_layer_thickness():
    grid = CartesianGrid(3, 3, 3)
    geometry = CellGeometry(grid, 10.0, 10.0, dz=[2.0, 8.0, 2.0], top_depth=2000.0)
    assert depth_to_k(15.0, 15.0, 2001.0, geometry) == 0     # 2000-2002
    assert depth_to_k(15.0, 15.0, 2005.0, geometry) == 1     # 2002-2010
    assert depth_to_k(15.0, 15.0, 2011.0, geometry) == 2     # 2010-2012


# ── depth_to_k: maili lay ─────────────────────────────────────────────────
def test_depth_to_k_dipping_layer_uses_column_specific_top():
    """Sütun dərinliyi maili olanda eyni mütləq dərinlik fərqli k verir."""
    nx, ny, nz = 5, 1, 3
    grid = CartesianGrid(nx, ny, nz)
    dip_per_cell = 10.0
    surface = 2000.0 + dip_per_cell * np.arange(nx)     # (nx,) → (ny=1, nx)
    surface = surface.reshape(ny, nx)
    geometry = CellGeometry(grid, 20.0, 20.0, dz=5.0,
                            top_depth=2000.0, top_depth_map=surface)

    # sütun 0: tavan 2000 → k0 = [2000,2005)
    assert depth_to_k(5.0, 5.0, 2002.0, geometry) == 0
    # sütun 4 (ən dərin): tavan 2040 → eyni mütləq dərinlik (2002) artıq
    # grid-dən yuxarıdadır (None), amma 2042 hələ k0-dadır
    assert depth_to_k(85.0, 5.0, 2002.0, geometry) is None
    assert depth_to_k(85.0, 5.0, 2042.0, geometry) == 0


def test_depth_to_k_column_lookup_uses_correct_ij_not_flat_index():
    """Sütun tapılması `xy_to_ij` ilə eyni məntiqi işlətməlidir."""
    nx, ny, nz = 3, 3, 2
    grid = CartesianGrid(nx, ny, nz)
    surface = np.array([[2000.0, 2000.0, 2000.0],
                        [2000.0, 2000.0, 2000.0],
                        [2000.0, 2000.0, 2100.0]])   # yalnız (i=2,j=2) dərin
    geometry = CellGeometry(grid, 10.0, 10.0, dz=5.0,
                            top_depth=2000.0, top_depth_map=surface)
    # (i=2, j=2) sütununda tavan 2100 -> 2002 grid-dən kənar
    assert depth_to_k(25.0, 25.0, 2002.0, geometry) is None
    # (i=0, j=0) sütununda tavan 2000 -> 2002 birinci təbəqə
    assert depth_to_k(5.0, 5.0, 2002.0, geometry) == 0
