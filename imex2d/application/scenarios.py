"""Hazır quyu sxemləri və sintetik geoloji model qurucusu.

Bunlar əvvəl UI kodunun içində idi (five_spot_wells, perm_field metodu).
Biznes məntiqi olduğu üçün UI-dən çıxarılıb.
"""

from __future__ import annotations
from typing import List

import numpy as np

from ..domain.geological_model import GeologicalModel
from ..domain.geometry import CellGeometry
from ..domain.grid import CartesianGrid
from ..domain.properties import PropertyMap
from ..domain.structure import RegionSet
from ..domain.wells import (ControlMode, Perforation, Well, WellControl,
                            WellType)


def _vertical_well(name, grid: CartesianGrid, i, j, well_type, mode, target,
                   radius=0.1, skin=0.0) -> Well:
    """Bütün təbəqələrdə perforasiya olunmuş şaquli quyu.

    2D-də (nz = 1) nəticə əvvəlki ilə eynidir — bir perforasiya.
    """
    return Well(name=name, well_type=well_type,
                control=WellControl(mode, target),
                perforations=[Perforation(i, j, k, True, skin)
                              for k in range(grid.nz)],
                radius=radius)


def five_spot(grid: CartesianGrid, p_inj=320.0, p_prod=150.0) -> List[Well]:
    return [
        _vertical_well("INJ-1", grid, 0, 0, WellType.INJECTOR,
                       ControlMode.BHP, p_inj),
        _vertical_well("PROD-1", grid, grid.nx - 1, grid.ny - 1,
                       WellType.PRODUCER, ControlMode.BHP, p_prod),
    ]


def line_drive(grid: CartesianGrid, p_inj=320.0, p_prod=150.0) -> List[Well]:
    j = grid.ny // 2
    return [
        _vertical_well("INJ-1", grid, 0, j, WellType.INJECTOR,
                       ControlMode.BHP, p_inj),
        _vertical_well("PROD-1", grid, grid.nx - 1, j, WellType.PRODUCER,
                       ControlMode.BHP, p_prod),
    ]


def central_injector(grid: CartesianGrid, p_inj=320.0, p_prod=150.0) -> List[Well]:
    j, i = grid.ny // 2, grid.nx // 2
    return [
        _vertical_well("INJ-1", grid, i, j, WellType.INJECTOR,
                       ControlMode.BHP, p_inj),
        _vertical_well("PROD-1", grid, 0, j, WellType.PRODUCER,
                       ControlMode.BHP, p_prod),
        _vertical_well("PROD-2", grid, grid.nx - 1, j, WellType.PRODUCER,
                       ControlMode.BHP, p_prod),
    ]


WELL_PATTERNS = {
    "Five-spot (1/4)": five_spot,
    "Xətti sıxışdırma (line drive)": line_drive,
    "İki hasilat + mərkəzi vurucu": central_injector,
}


def lognormal_permeability(grid: CartesianGrid, base: float, sigma: float,
                           seed: int, smoothing: int = 4) -> np.ndarray:
    """Sintetik heterogen keçiricilik sahəsi, (nz, ny, nx).

    Hər təbəqə üçün ayrı realizasiya qurulur (seed təbəqə nömrəsi ilə
    dəyişir) — bu, şaquli heterogenliyi yaradır və 3D-də süpürmənin
    təbəqələr üzrə fərqlənməsini görməyə imkan verir.
    """
    layers = []
    for k in range(grid.nz):
        rng = np.random.default_rng(seed + k * 1000)
        f = rng.normal(0.0, 1.0, (grid.ny, grid.nx))
        for _ in range(smoothing):
            f = (f + np.roll(f, 1, 0) + np.roll(f, -1, 0)
                 + np.roll(f, 1, 1) + np.roll(f, -1, 1)) / 5.0
        f /= max(f.std(), 1e-9)
        layers.append(np.clip(base * np.exp(sigma * f), 0.05, 1e5))
    return np.stack(layers)


def dipping_surface(grid: CartesianGrid, top_depth: float,
                    dip_x: float, dip_y: float):
    """Maili lay tavanı: hər hüceyrədə dərinlik, m.

    Yalnız həndəsə generatorudur — mailliyi olmayan halda None qaytarır
    ki, sabit `top_depth` yolu (köhnə davranış) işləsin.
    """
    if abs(dip_x) < 1e-12 and abs(dip_y) < 1e-12:
        return None
    i = np.arange(grid.nx)
    j = np.arange(grid.ny)
    jj, ii = np.meshgrid(j, i, indexing="ij")
    return top_depth + ii * dip_x + jj * dip_y


class SyntheticGeologicalModelBuilder:
    """UI-dən gələn sadə parametrlərdən geoloji model qurur."""

    def build(self, nx, ny, dx, dy, dz, porosity, permx_base,
              ky_over_kx=1.0, heterogeneous=False, sigma=0.5, seed=7,
              top_depth=0.0, dip_x=0.0, dip_y=0.0, nz=1, kv_over_kh=0.1,
              name="Sintetik geoloji model") -> GeologicalModel:
        grid = CartesianGrid(nx, ny, nz)
        geometry = CellGeometry(grid, dx, dy, dz, top_depth=top_depth,
                                top_depth_map=dipping_surface(grid, top_depth,
                                                              dip_x, dip_y))
        n = grid.ncell

        if heterogeneous:
            kx = lognormal_permeability(grid, permx_base, sigma, seed).ravel()
        else:
            kx = np.full(n, float(permx_base))

        model = GeologicalModel(name=name, grid=grid, geometry=geometry,
                                regions=RegionSet.single(n))
        model.add_property(PropertyMap.from_array("PORO", porosity, n))
        model.add_property(PropertyMap.from_array("PERMX", kx, n, "mD"))
        model.add_property(PropertyMap.from_array("PERMY", kx * ky_over_kx, n, "mD"))
        model.add_property(PropertyMap.from_array("PERMZ", kx * kv_over_kh, n, "mD"))
        return model
