"""Hüceyrə həndəsəsi — həcmlər, üz sahələri, mərkəzlər arası məsafələr.

Simulyator yalnız bu interfeysə güvənir. Corner-point və ya qeyri-struktur
grid gələndə burada yeni sinif yazılır, hesablama nüvəsi dəyişmir.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

import numpy as np

from .grid import CartesianGrid, Connections


@dataclass(frozen=True)
class CellGeometry:
    """Bərabər ölçülü kartezian bloklar üçün həndəsə."""
    grid: CartesianGrid
    dx: float
    dy: float
    dz: float
    top_depth: float = 0.0
    top_depth_map: Optional[np.ndarray] = None

    @property
    def cell_volume(self) -> float:
        return self.dx * self.dy * self.dz

    def volumes(self) -> np.ndarray:
        return np.full(self.grid.ncell, self.cell_volume)

    def cell_depths(self) -> np.ndarray:
        """Hər hüceyrənin mərkəz dərinliyi, m.

        `top_depth_map` verilibsə lay maili/qırışıqlı ola bilər; verilməyibsə
        sabit `top_depth` işlədilir. Bu, həndəsə məlumatıdır — initialization
        provider-i buradan oxuyur, öz dərinlik modelini qurmur.
        """
        grid = self.grid
        k = np.repeat(np.arange(grid.nz), grid.nx * grid.ny)
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
        return top + (k + 0.5) * self.dz

    def face_areas(self, conn: Connections) -> np.ndarray:
        area = np.empty(conn.count)
        area[conn.axis == 0] = self.dy * self.dz
        area[conn.axis == 1] = self.dx * self.dz
        area[conn.axis == 2] = self.dx * self.dy
        return area

    def face_half_distances(self, conn: Connections) -> np.ndarray:
        d = np.empty(conn.count)
        d[conn.axis == 0] = self.dx * 0.5
        d[conn.axis == 1] = self.dy * 0.5
        d[conn.axis == 2] = self.dz * 0.5
        return d

    def areal_extent(self) -> tuple:
        return (self.grid.nx * self.dx, self.grid.ny * self.dy)
