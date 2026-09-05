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
    """WI = 2π·C·√(K1·K2)·h / (ln(re/rw) + S)

    `K1`/`K2`/`d1`/`d2`/`h` — hüceyrənin HƏQİQİ wellblock həndəsəsindən
    (bax `domain/geometry.py::WellblockGeometry` və
    `CornerPointGeometry.wellblock_geometry`), ox-boyu sərhəd qutusundan
    DEYİL. Kartezian blokda onlar `Kx`/`Ky`/`dx`/`dy`/`dz_k`-ya
    eyniliklə bərabərdir, yəni klassik nəticə DƏYİŞMİR.
    """

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
        c_darcy = model.units.darcy_constant
        actnum = model.grid.active.actnum if model.grid.has_inactive_cells else None

        # Perforasiyalar ƏVVƏLCƏ toplanır, həndəsə isə TOPLU hesablanır:
        # `wellblock_geometry()` vektorlaşdırılıb, hər perforasiya üçün
        # ayrı çağırmaq eyni nəticəni verir, sadəcə israfdır.
        perforations = []                    # (well, perf, cell, direction)
        for well in model.active_wells():
            for perf in well.open_perforations():
                cell = model.grid.index(perf.i, perf.j, perf.k)
                if actnum is not None and actnum[cell] <= 0:
                    LOG.warning(
                        "%s: perforasiya (i=%d, j=%d, k=%d) qeyri-aktiv "
                        "hüceyrədədir (ACTNUM = 0) — söndürüldü (WI = 0).",
                        well.name, perf.i, perf.j, perf.k + 1)
                    continue
                perforations.append((well, perf, cell,
                                     getattr(perf, "direction", "Z")))
        if not perforations:
            return out

        metrics = self._wellblock_metrics(model, perforations)
        for slot, (well, perf, cell, _direction) in enumerate(perforations):
            wi = self._well_index(metrics["k1"][slot], metrics["k2"][slot],
                                  metrics["d1"][slot], metrics["d2"][slot],
                                  metrics["h"][slot], well.radius,
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
    def _wellblock_metrics(model: ReservoirModel, perforations) -> dict:
        """`{d1, d2, h, k1, k2}` — hər perforasiya üçün HƏQİQİ hüceyrə
        həndəsəsindən çıxarılmış Peaceman girişləri.

        NİYƏ `cell_extents()` DEYİL (Phase 3 düzəlişi): o, ox-boyu
        SƏRHƏD QUTUSUDUR (`max(x)−min(x)`), yəni fırlanmış/kəsilmiş
        corner-point hüceyrəsində həqiqi wellblock enindən böyükdür —
        `r_e` şişir, WI süni azalır. `wellblock_geometry()` isə quyu
        oxuna perpendikulyar müstəvidə HƏQİQİ en kəsiyi (`V/h`) və
        yerli kənar istiqamətlərini işlədir; Kartezian blokda ikisi
        maşın dəqiqliyində eynidir (bax `WellblockGeometry`).

        Keçiricilik də həmin YERLİ istiqamətlərə proyeksiya olunur
        (`uᵀ·K·u`), ona görə anizotrop Peaceman məntiqi olduğu kimi
        qalır — sadəcə `Kx`/`Ky` əvəzinə ox-uyğun `K1`/`K2` alır.
        Kartezian halda `K1 = Kx`, `K2 = Ky`.
        """
        rock = model.rock
        kx = rock.permx.values
        ky = rock.permy.values
        # PERMZ verilməyəndə `Kz = Kx` — bu kod bazasının hər yerində
        # işlənən konvensiya (bax `simulation/discretization.py`,
        # `discretization/mpfa_o.py`). Şaquli quyuda `Kz` onsuz da
        # heç bir çəki almır (yerli oxlar üfüqidir).
        kz = rock.permz.values if rock.permz is not None else kx
        k_diagonal = np.column_stack([kx, ky, kz])

        count = len(perforations)
        metrics = {name: np.empty(count) for name in ("d1", "d2", "h", "k1", "k2")}
        cells = np.array([entry[2] for entry in perforations], dtype=int)
        directions = [entry[3] for entry in perforations]

        for direction in dict.fromkeys(directions):          # sıra sabit
            selected = np.array([d == direction for d in directions])
            block = model.geometry.wellblock_geometry(cells[selected], direction)
            k1, k2 = block.directional_permeability(k_diagonal[cells[selected]])
            metrics["d1"][selected] = block.d1
            metrics["d2"][selected] = block.d2
            metrics["h"][selected] = block.length
            metrics["k1"][selected] = k1
            metrics["k2"][selected] = k2
        return metrics

    @staticmethod
    def _well_index(k1, k2, d1, d2, h, rw, skin, c_darcy) -> float:
        """Anizotrop Peaceman — DÜSTUR DƏYİŞMƏYİB.

        `k1`/`k2` quyu oxuna perpendikulyar iki YERLİ istiqamətdəki
        keçiricilik, `d1`/`d2` həmin istiqamətlərdəki effektiv wellblock
        ölçüləri, `h` isə perforasiyanın hüceyrə içindəki uzunluğudur.
        """
        kh = np.sqrt(k1 * k2)
        ratio = k2 / k1 if k1 > 0.0 else 1.0
        r1 = np.sqrt(ratio) * d1 ** 2
        r2 = np.sqrt(1.0 / ratio) * d2 ** 2
        re = 0.28 * np.sqrt(r1 + r2) / (ratio ** 0.25 + (1.0 / ratio) ** 0.25)
        return float(c_darcy * 2.0 * np.pi * kh * h /
                     (np.log(max(re / rw, 1.01)) + skin))
