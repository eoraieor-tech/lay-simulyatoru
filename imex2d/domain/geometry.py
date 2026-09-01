"""Hüceyrə həndəsəsi — həcmlər, üz sahələri, mərkəzlər arası məsafələr.

Simulyator yalnız bu interfeysə güvənir. Corner-point və ya qeyri-struktur
grid gələndə burada yeni sinif yazılır, hesablama nüvəsi dəyişmir.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Union, Sequence

import numpy as np

from .grid import CartesianGrid, Connections
from .validation import validate_cell_volumes, validate_grid_dimensions, validate_thickness


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

    def validate(self) -> list:
        """Dejenerativ həndəsəni aşkarlayır — sıfır/mənfi ölçü, sıfır/
        mənfi hüceyrə həcmi (audit: bu yoxlama əvvəllər HEÇ YERDƏ yox
        idi, bax `GEOSTATISTICS.md`-dən sonrakı Phase 1 hesabatı).

        YENİ, AYRICA metoddur — `__post_init__`-ə ƏLAVƏ EDİLMƏYİB ki,
        mövcud konstruksiya yolları (672 test) DƏYİŞMƏSİN. Çağıran
        (məs. `GeologicalModel.validate()`) bunu İSTƏYƏ görə çağırır.
        """
        issues = []
        grid_result = validate_grid_dimensions(self.grid.nx, self.grid.ny, self.grid.nz,
                                               self.dx, self.dy)
        issues.extend(grid_result.errors)
        thickness_result = validate_thickness(self.dz, label="DZ")
        issues.extend(thickness_result.errors)
        if not thickness_result.errors and not grid_result.errors:
            volume_result = validate_cell_volumes(self.volumes(), label="hüceyrə həcmi")
            issues.extend(volume_result.errors)
        return issues


def xy_to_ij(x: float, y: float, geometry: CellGeometry) -> tuple:
    """Metr koordinatını hüceyrə indeksinə çevirir.

    Origin (0, 0) qəbul edilir — bu kod bazasında koordinat sistemi
    həmişə grid-in aşağı-sol küncündən başlayır (bax `LOCAL`
    `coordinate_system`), ayrıca origin sahəsi yoxdur.

    `x == x_max` sərhədində `i = nx` çıxmasın deyə yuxarı hədd `nx - 1`-ə
    kəsilir (tapşırıqda tələb olunan qayda); aşağı hədddə də simmetrik
    olaraq `0`-a kəsilir ki, grid-dən kənar mənfi X/Y mənfi indeksə yox,
    ən yaxın kənar hüceyrəyə düşsün (özü `validate_wells`-də ayrıca
    "sərhəddən kənar" xətası kimi bildirilir — bu funksiya heç vaxt
    xəta atmır, yalnız ən yaxın hüceyrəni qaytarır).
    """
    grid = geometry.grid
    i = int((x - 0.0) / geometry.dx) if geometry.dx > 0 else 0
    j = int((y - 0.0) / geometry.dy) if geometry.dy > 0 else 0
    i = min(max(i, 0), grid.nx - 1)
    j = min(max(j, 0), grid.ny - 1)
    return i, j


def depth_to_k(x: float, y: float, depth: float, geometry: CellGeometry):
    """Dərinliyi (m) verilmiş `(x, y)` sütununda təbəqə indeksinə çevirir.

    Lay üfüqi deyil — hər `(i, j)` sütununun tavan dərinliyi
    `top_depth_map`-dan (varsa) və ya sabit `top_depth`-dən götürülür,
    sonra kumulyativ `dz` ilə təbəqə sərhədləri qurulur. Dərinlik grid-in
    diapazonundan kənardadırsa `None` qaytarır — çağıran bunu "grid
    qurulduqdan sonra" mesajı kimi göstərməlidir, xəta atılmır.
    """
    i, j = xy_to_ij(x, y, geometry)
    grid = geometry.grid
    if geometry.top_depth_map is None:
        column_top = geometry.top_depth
    else:
        areal = np.asarray(geometry.top_depth_map, float).ravel()
        if areal.size == grid.nx * grid.ny:
            column_top = float(areal.reshape(grid.ny, grid.nx)[j, i])
        elif areal.size == grid.ncell:
            column_top = float(areal.reshape(grid.shape)[0, j, i])
        else:
            column_top = geometry.top_depth

    edges = np.concatenate(([0.0], np.cumsum(geometry.dz)))
    if depth < column_top + edges[0] or depth > column_top + edges[-1]:
        return None
    k = int(np.searchsorted(edges, depth - column_top, side="right")) - 1
    return int(min(max(k, 0), grid.nz - 1))
