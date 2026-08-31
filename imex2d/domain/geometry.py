"""Hüceyrə həndəsəsi — həcmlər, üz sahələri, mərkəzlər arası məsafələr.

Simulyator yalnız bu interfeysə güvənir. Corner-point və ya qeyri-struktur
grid gələndə burada yeni sinif yazılır, hesablama nüvəsi dəyişmir.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Union, Sequence

import numpy as np

from .grid import CartesianGrid, Connections


@dataclass(frozen=True)
class CellGeometry:
    """Kartezian bloklar üçün həndəsə.

    `dz` hər K-təbəqəsinin qalınlığıdır: tək `float` (bütün təbəqələr
    eyni) və ya `nz` uzunluqlu ardıcıllıq (hər təbəqə ayrı) qəbul edir.
    Daxildə həmişə `np.ndarray` şəklində (`nz`,) saxlanılır.
    """
    grid: CartesianGrid
    dx: float
    dy: float
    dz: Union[float, Sequence[float]]
    top_depth: float = 0.0
    top_depth_map: Optional[np.ndarray] = None

    def __post_init__(self):
        arr = np.asarray(self.dz, dtype=float)
        if arr.ndim == 0:
            arr = np.full(self.grid.nz, float(arr))
        elif arr.size != self.grid.nz:
            raise ValueError(
                f"dz: {arr.size} dəyər, gözlənilən {self.grid.nz} (NZ)")
        else:
            arr = arr.ravel().copy()
        object.__setattr__(self, "dz", arr)

    def dz_per_cell(self) -> np.ndarray:
        """Hər hüceyrənin öz təbəqəsinin qalınlığı (uzunluq = ncell)."""
        return np.repeat(self.dz, self.grid.nx * self.grid.ny)

    def volumes(self) -> np.ndarray:
        return self.dx * self.dy * self.dz_per_cell()

    def cell_depths(self) -> np.ndarray:
        """Hər hüceyrənin mərkəz dərinliyi, m.

        `top_depth_map` verilibsə lay maili/qırışıqlı ola bilər; verilməyibsə
        sabit `top_depth` işlədilir. Bu, həndəsə məlumatıdır — initialization
        provider-i buradan oxuyur, öz dərinlik modelini qurmur.

        Təbəqələr fərqli qalınlıqda ola bildiyi üçün mərkəz dərinliyi
        kumulyativ təbəqə tavanlarından hesablanır, sadə `(k+0.5)*dz`
        yox.
        """
        grid = self.grid
        k = np.repeat(np.arange(grid.nz), grid.nx * grid.ny)
        layer_top_offset = np.concatenate(([0.0], np.cumsum(self.dz)[:-1]))
        layer_centre_offset = layer_top_offset + self.dz * 0.5
        centre_offset = np.repeat(layer_centre_offset, grid.nx * grid.ny)
        if self.top_depth_map is None:
            top = np.full(grid.ncell, self.top_depth)
        else:
            areal = np.asarray(self.top_depth_map, float).ravel()
            if areal.size == grid.ncell:
                top = areal
            elif areal.size == grid.nx * grid.ny:
                top = np.tile(areal, grid.nz)
            else:
                raise ValueError("top_depth_map ölçüsü grid ilə uyğun gəlmir")
        return top + centre_offset

    def face_areas(self, conn: Connections) -> np.ndarray:
        dz_cell = self.dz_per_cell()
        area = np.empty(conn.count)
        m0, m1, m2 = conn.axis == 0, conn.axis == 1, conn.axis == 2
        area[m0] = self.dy * dz_cell[conn.cell_a[m0]]
        area[m1] = self.dx * dz_cell[conn.cell_a[m1]]
        area[m2] = self.dx * self.dy
        return area

    def face_half_distances(self, conn: Connections) -> tuple:
        """Hər üzün hər tərəfindən mərkəzə qədər yarım-məsafə.

        `(half_a, half_b)` qaytarır — K istiqamətində qonşu təbəqələrin
        qalınlığı fərqli ola bildiyi üçün iki tərəf ayrı hesablanır;
        I/J istiqamətində dz-dən asılı olmadığı üçün iki tərəf eynidir.
        """
        dz_cell = self.dz_per_cell()
        half_a = np.empty(conn.count)
        half_b = np.empty(conn.count)
        m0, m1, m2 = conn.axis == 0, conn.axis == 1, conn.axis == 2
        half_a[m0] = half_b[m0] = self.dx * 0.5
        half_a[m1] = half_b[m1] = self.dy * 0.5
        half_a[m2] = dz_cell[conn.cell_a[m2]] * 0.5
        half_b[m2] = dz_cell[conn.cell_b[m2]] * 0.5
        return half_a, half_b

    def areal_extent(self) -> tuple:
        return (self.grid.nx * self.dx, self.grid.ny * self.dy)
