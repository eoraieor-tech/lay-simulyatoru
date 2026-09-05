"""Grid həndəsəsi interfeysi — Kartezian VƏ gələcək corner-point üçün
ORTAQ müqavilə (bax audit tapşırığı §6).

İKİ implementasiya var:

  * `imex2d.domain.geometry.CellGeometry` — struktur/bərabər-blok
    Kartezian (ilk və hələ də defolt yol);
  * `imex2d.domain.corner_point_geometry.CornerPointGeometry` — HƏQİQİ
    corner-point (COORD/ZCORN → hüceyrə başına 8 təpə), Phase 5D. O,
    `CellGeometry`-nin ALT SİNFİDİR: eyni metodları eyni imzalarla,
    amma qutu düsturu yerinə dəqiq çoxüzlü riyaziyyatla tətbiq edir.

Bu ABC HEÇ NƏYİ DƏYİŞMİR — hər iki sinif metodların hamısını doğru
tətbiq edir. Məqsəd: hansı müqaviləyə əməl edilməli olduğunu AÇIQ etmək
— `TwoPointFluxDiscretization` və `MPFAODiscretization` YALNIZ bu
metodlardan istifadə etməlidir, `CellGeometry`-nin daxili strukturundan
(`dx`/`dy`/`dz`) BİRBAŞA YOX. Məhz buna görə corner-point dəstəyi
diskretizasiya kodunda BİR SƏTİR belə dəyişiklik tələb etmədi (bax
`ARCHITECTURE.md` §5.1/§5.18).

Phase 3 ("General Geometry Foundation for MPFA-O") ƏLAVƏSİ: bu ABC
DƏYİŞDİRİLMƏDİ (audit §4 nəticəsi — "improve if necessary", NECESSARY
tapılmadı). Səbəb: bu, VEKTORLAŞDIRILMIŞ (bütün grid TƏK massivlərlə)
GRID-səviyyəli müqavilədir; `imex2d.domain.polyhedral_geometry.
HexahedralCell`/`Face` isə HÜCEYRƏ-BAŞINA (per-cell), potensial qeyri-
ortoqonal ÜMUMİ həndəsə NÜVƏSİDİR — bir səviyyə AŞAĞIDA yerləşir, bu
ABC-ni ƏVƏZ ETMİR. `CornerPointGeometry` (Phase 5D) məhz bunu edir —
eyni riyaziyyatı VEKTORLAŞDIRILMIŞ formada tətbiq edib bu ABC-nin
müqaviləsini doldurur; yəni `polyhedral_geometry.py` bir `IGridGeometry`
İMPLEMENTASİYASININ DAXİLİ ALƏTİDİR, ƏLAVƏ/ALTERNATİV İNTERFEYS YOX.
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
