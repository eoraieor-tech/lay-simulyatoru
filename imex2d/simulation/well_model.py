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
from ..logging_setup import get_logger

LOG = get_logger(__name__)


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
        """Açıq perforasiyalardan bağlantı siyahısı.

        ACTNUM (tapşırıq §4): QEYRİ-AKTİV hüceyrəyə (ACTNUM = 0) düşən
        perforasiya bu siyahıya SALINMIR — yəni onun WI-si effektiv
        SIFIRDIR. Səbəb sırf fizikidir: həmin hüceyrənin məsamə həcmi 0,
        qonşuluq bağlantısı yoxdur və xətti sistemdə naməlumu yoxdur —
        ora quyu mənbəyi yazmaq maddəni MODELDƏN KƏNARA vurmaq
        (hasilatda isə YOXDAN maye çıxarmaq) demək olardı.

        Bağlantı SƏSSİZCƏ atılmır: `journal`-a yazılır və
        `ReservoirModel._check_wells` eyni vəziyyəti diaqnostika
        hesabatında (xəbərdarlıq / bütün perforasiyalar qeyri-aktivdirsə
        XƏTA) göstərir.
        """
        out: List[WellConnection] = []
        kx_all = model.rock.permx.values
        ky_all = model.rock.permy.values
        dx, dy = model.geometry.dx, model.geometry.dy
        dz_all = model.geometry.dz_per_cell()
        c_darcy = model.units.darcy_constant
        actnum = model.grid.active.actnum if model.grid.has_inactive_cells else None

        for well in model.active_wells():
            for perf in well.open_perforations():
                cell = model.grid.index(perf.i, perf.j, perf.k)
                if actnum is not None and actnum[cell] <= 0:
                    LOG.warning(
                        "%s: perforasiya (i=%d, j=%d, k=%d) qeyri-aktiv "
                        "hüceyrədədir (ACTNUM = 0) — söndürüldü (WI = 0).",
                        well.name, perf.i, perf.j, perf.k + 1)
                    continue
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
