"""Grid həndəsəsi interfeysi — Kartezian VƏ gələcək corner-point üçün
ORTAQ müqavilə (bax audit tapşırığı §6).

Hazırda YEGANƏ implementasiya `imex2d.domain.geometry.CellGeometry`-dir
(struktur/bərabər-blok Kartezian). Bu ABC HEÇ NƏYİ DƏYİŞMİR — `CellGeometry`
artıq bu metodların hamısını (bəziləri bu fazada ƏLAVƏ EDİLİB, bax
`geometry.py`) doğru tətbiq edir. Məqsəd: gələcək corner-point/qeyri-struktur
həndəsə sinfi YAZILANDA, hansı müqaviləyə əməl etməli olduğunu AÇIQ etmək —
`TwoPointFluxDiscretization` və gələcək `MPFAODiscretization` YALNIZ bu
metodlardan istifadə etməlidir, `CellGeometry`-nin daxili strukturundan
(`dx`/`dy`/`dz`) BİRBAŞA YOX (MPFA-O həndəsə-müstəqil qala bilsin deyə).

Native corner-point (COORD/ZCORN-dən BİRBAŞA) HƏLƏ İMPLEMENTASİYA
EDİLMİR — bax `ECLIPSE_IO.md`.
"""

from __future__ import annotations
from abc import ABC, abstractmethod

import numpy as np

from ..domain.grid import Connections


class IGridGeometry(ABC):
    """Hüceyrə/üz həndəsəsi — mərkəzlər, həcmlər, sahələr, normallar."""

    @abstractmethod
    def cell_centroid(self) -> np.ndarray:
        """(ncell, 3) — hər hüceyrənin mərkəzi [X, Y, Z], metr."""

    @abstractmethod
    def volumes(self) -> np.ndarray:
        """(ncell,) — hüceyrə həcmi, m³."""

    @abstractmethod
    def face_centroid(self, conn: Connections) -> np.ndarray:
        """(nface, 3) — hər üzün mərkəzi [X, Y, Z], metr."""

    @abstractmethod
    def face_areas(self, conn: Connections) -> np.ndarray:
        """(nface,) — üz sahəsi, m²."""

    @abstractmethod
    def face_normal(self, conn: Connections) -> np.ndarray:
        """(nface, 3) — vahid normal vektor, `cell_a`-dan `cell_b`-yə."""
