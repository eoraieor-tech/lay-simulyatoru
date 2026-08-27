"""Diskretizasiya — ReservoirModel-dən transmissivlik hesablayır.

Bu kod mövcud nüvədən köçürülüb, riyaziyyat dəyişməyib (TPFA, harmonik
orta). Fərq: əvvəl bu hesablama simulyator sinfinin içində idi və
grid ölçülərini birbaşa oxuyurdu; indi model + həndəsə interfeysindən
işləyir, ona görə corner-point həndəsəsi gələndə dəyişmir.
"""

from __future__ import annotations
from dataclasses import dataclass

import numpy as np

from ..domain.grid import Connections
from ..domain.reservoir_model import ReservoirModel


@dataclass(frozen=True)
class DiscretizedGrid:
    connections: Connections
    transmissibility: np.ndarray
    pore_volume: np.ndarray
    cell_volume: np.ndarray


class TwoPointFluxDiscretization:
    """İki nöqtəli axın approksimasiyası (TPFA)."""

    def build(self, model: ReservoirModel) -> DiscretizedGrid:
        conn = model.connections()
        geom = model.geometry
        area = geom.face_areas(conn)
        half = geom.face_half_distances(conn)
        perm = self._directional_permeability(model, conn)

        k_a = np.maximum(perm[0], 1e-9)
        k_b = np.maximum(perm[1], 1e-9)
        trans = model.units.darcy_constant * area / (half / k_a + half / k_b)
        trans = self._apply_fault_multipliers(model, conn, trans)

        return DiscretizedGrid(
            connections=conn,
            transmissibility=trans,
            pore_volume=model.pore_volume(),
            cell_volume=geom.volumes(),
        )

    @staticmethod
    def _directional_permeability(model: ReservoirModel, conn: Connections):
        kx = model.rock.permx.values
        ky = model.rock.permy.values
        kz = model.rock.permz.values if model.rock.permz is not None else kx
        per_axis = (kx, ky, kz)
        k_a = np.empty(conn.count)
        k_b = np.empty(conn.count)
        for axis in (0, 1, 2):
            m = conn.axis == axis
            if np.any(m):
                k_a[m] = per_axis[axis][conn.cell_a[m]]
                k_b[m] = per_axis[axis][conn.cell_b[m]]
        return k_a, k_b

    @staticmethod
    def _apply_fault_multipliers(model: ReservoirModel, conn: Connections,
                                 trans: np.ndarray) -> np.ndarray:
        """Fay transmissivlik çarpanları (B3).

        Hər fay bir müstəvidir: yalnız həmin müstəvidəki üzlərin
        transmissivliyi çarpanla dəyişir, qalan bütün üzlər toxunulmaz
        qalır. Bir üz birdən çox faya düşərsə, çarpanlar VURULUR —
        iki qismən keçirici fay eyni yerdə üst-üstə düşəndə axının
        daha da azalması fiziki cəhətdən doğrudur.
        """
        faults = [f for f in model.fault_references if f.has_geometry]
        if not faults:
            return trans

        grid = model.grid
        i_a, j_a, k_a = grid.ijk_array(conn.cell_a)
        coordinates = {0: (j_a, k_a), 1: (i_a, k_a), 2: (i_a, j_a)}
        axis_code = {"I": 0, "J": 1, "K": 2}
        boundary = {0: i_a, 1: j_a, 2: k_a}

        trans = trans.copy()
        for fault in faults:
            code = axis_code[fault.axis.upper()]
            coordinate_a, coordinate_b = coordinates[code]
            mask = ((conn.axis == code)
                    & (boundary[code] == fault.plane_index)
                    & fault.matches(code, coordinate_a, coordinate_b))
            if np.any(mask):
                trans[mask] *= fault.effective_multiplier
        return trans
