"""Quyu modeli — Peaceman bağlantı əmsalı.

Kod mövcud nüvədən köçürülüb. Fərq: quyu məlumatı (domain.wells.Well)
ilə hesablama (bu fayl) ayrılıb. Well artıq özünün hüceyrə indeksini
və ya WI-ni bilmir — bu, simulyasiya təfərrüatıdır.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List

import numpy as np

from ..domain.reservoir_model import ReservoirModel
from ..domain.wells import ControlMode, WellType


@dataclass
class WellConnection:
    """Bir perforasiyanın grid ilə əlaqəsi."""
    well_name: str
    cell: int
    well_index: float
    is_injector: bool
    mode: ControlMode
    target: float


class PeacemanWellModel:
    """WI = 2π·C·√(Kx·Ky)·h / (ln(re/rw) + S)"""

    def build_connections(self, model: ReservoirModel) -> List[WellConnection]:
        out: List[WellConnection] = []
        kx_all = model.rock.permx.values
        ky_all = model.rock.permy.values
        dx, dy = model.geometry.dx, model.geometry.dy
        dz_all = model.geometry.dz_per_cell()
        c_darcy = model.units.darcy_constant

        for well in model.active_wells():
            for perf in well.open_perforations():
                cell = model.grid.index(perf.i, perf.j, perf.k)
                kx, ky = kx_all[cell], ky_all[cell]
                dz = dz_all[cell]
                wi = self._well_index(kx, ky, dx, dy, dz, well.radius,
                                      perf.skin, c_darcy)
                out.append(WellConnection(
                    well_name=well.name,
                    cell=cell,
                    well_index=wi,
                    is_injector=well.well_type is WellType.INJECTOR,
                    mode=well.control.mode,
                    target=well.control.target,
                ))
        return out

    @staticmethod
    def _well_index(kx, ky, dx, dy, dz, rw, skin, c_darcy) -> float:
        kh = np.sqrt(kx * ky)
        r1 = np.sqrt(ky / kx) * dx ** 2
        r2 = np.sqrt(kx / ky) * dy ** 2
        re = 0.28 * np.sqrt(r1 + r2) / ((ky / kx) ** 0.25 + (kx / ky) ** 0.25)
        return float(c_darcy * 2.0 * np.pi * kh * dz /
                     (np.log(max(re / rw, 1.01)) + skin))
