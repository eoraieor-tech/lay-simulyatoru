"""Grid topologiyası — YALNIZ hüceyrələrin bir-birinə necə bağlandığı.

Ölçülər (həcm, sahə, məsafə) burada deyil, CellGeometry-dədir.
Bu ayrılıq gələcəkdə corner-point grid-ə keçidi mümkün edir: topologiya
qalır, yalnız geometriya dəyişir.
"""

from __future__ import annotations
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Connections:
    """Hüceyrələr arası üzlərin siyahısı (qonşuluq qrafı)."""
    cell_a: np.ndarray
    cell_b: np.ndarray
    axis: np.ndarray

    @property
    def count(self) -> int:
        return int(self.cell_a.size)


@dataclass(frozen=True)
class CartesianGrid:
    """Struktur (kartezian) grid topologiyası."""
    nx: int
    ny: int
    nz: int = 1

    @property
    def ncell(self) -> int:
        return self.nx * self.ny * self.nz

    @property
    def shape(self) -> tuple:
        return (self.nz, self.ny, self.nx)

    def index(self, i: int, j: int, k: int = 0) -> int:
        return (k * self.ny + j) * self.nx + i

    def ijk(self, c: int) -> tuple:
        i = c % self.nx
        j = (c // self.nx) % self.ny
        k = c // (self.nx * self.ny)
        return i, j, k

    def ijk_array(self, cells: np.ndarray) -> tuple:
        """`ijk`-nin vektor forması — massiv üçün bir dəfə hesablanır."""
        cells = np.asarray(cells, dtype=np.int64)
        i = cells % self.nx
        j = (cells // self.nx) % self.ny
        k = cells // (self.nx * self.ny)
        return i, j, k

    def build_connections(self) -> Connections:
        idx = np.arange(self.ncell).reshape(self.shape)
        a_list, b_list, ax_list = [], [], []
        if self.nx > 1:
            a_list.append(idx[:, :, :-1].ravel())
            b_list.append(idx[:, :, 1:].ravel())
            ax_list.append(np.zeros(idx[:, :, :-1].size, dtype=np.int8))
        if self.ny > 1:
            a_list.append(idx[:, :-1, :].ravel())
            b_list.append(idx[:, 1:, :].ravel())
            ax_list.append(np.ones(idx[:, :-1, :].size, dtype=np.int8))
        if self.nz > 1:
            a_list.append(idx[:-1, :, :].ravel())
            b_list.append(idx[1:, :, :].ravel())
            ax_list.append(np.full(idx[:-1, :, :].size, 2, dtype=np.int8))
        if not a_list:
            e = np.zeros(0, dtype=int)
            return Connections(e, e, np.zeros(0, dtype=np.int8))
        return Connections(np.concatenate(a_list), np.concatenate(b_list),
                           np.concatenate(ax_list))
