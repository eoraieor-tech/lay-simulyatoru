"""Ümumi (potensial qeyri-ortoqonal) çoxüzlü hüceyrə həndəsəsi — MPFA-O
üçün geometriya NÜVƏSİ (Phase 3: "General Geometry Foundation for MPFA-O").

BU MODUL SAF HƏNDƏSƏDİR — heç bir transmissivlik, mobilite, təzyiq,
Darcy axını, Nyuton/Jakobian anlayışı YOXDUR (bax audit §21: "Geometry
must NOT contain transmissibility/mobility/.../Jacobian"). Diskretizasiya
qatı (`simulation/discretization.py`) bu modulu İSTEHLAK EDƏ BİLƏR,
əksinə YOX.

Niyə `domain/geometry.py`-ə ƏLAVƏ, ONUN ƏVƏZİNƏ YOX: `CellGeometry`
VEKTORLAŞDIRILMIŞ, bütün grid üçün TEK dəfəyə (numpy massivləri, struktur
Kartezian fərziyyəsi ilə) işləyir — bu, TPFA üçün son dərəcə səmərəlidir
və DƏYİŞDİRİLMİR (bax audit §14/§26: "existing CellGeometry must continue
to work... all existing results must remain unchanged"). Bu modul isə
HÜCEYRƏ-BAŞINA (per-cell) ÜMUMİ (potensial qeyri-ortoqonal, gələcəkdə
corner-point) həndəsəni təmsil edir — gələcək bir "ÜMUMI grid həndəsəsi"
sinfi (məs. `CornerPointGeometry`, HƏLƏ YAZILMAYIB) bu NÜVƏNİ hüceyrə-
hüceyrə çağırıb nəticələri `IGridGeometry`-nin vektorlaşdırılmış
müqaviləsinə (bax `interfaces/geometry.py`) "yığa" bilər — beləliklə
TPFA/MPFA-O bu iki tətbiqin HANSI olduğunu BİLMİR, YALNIZ `IGridGeometry`-ə
güvənir.

Qat diaqramı (bax audit §21, `ARCHITECTURE.md` §5.14):

    Geometry NÜVƏSİ (bu fayl: HexahedralCell/Face — saf riyaziyyat)
        ↓
    Grid Geometry (CellGeometry — Kartezian, VEKTORLAŞDIRILMIŞ; gələcək
                   CornerPointGeometry bu nüvəni İSTİFADƏ EDƏ BİLƏR)
        ↓
    Topology (`domain/grid.py::Connections` — HANSI hüceyrələr bağlıdır,
              HARADA olduqları İLƏ QARIŞDIRILMIR, bax audit §18)
        ↓
    Discretization (TPFA indiki, MPFA-O gələcək)
        ↓
    Flow

HEXAHEDRAL FƏRZİYYƏLƏR (bax audit §16 — "document assumptions"):
  - Hər üz DÖRDBUCAQLIDIR (4 təpə) və TAM MÜSTƏVİ OLMAYA BİLƏR (əyri/
    "warped" üz) — sahə/mərkəz/normal HƏR ÜZÜ İKİ ÜÇBUCAĞA bölərək
    (fan-triangulyasiya, `(v0,v1,v2)` və `(v0,v2,v3)`) hesablanır. Bu,
    HƏQİQİ əyri səthin YALNIZ TƏXMİNİDİR — tam müstəvi üzlər üçün DƏQİQ,
    əyri üzlər üçün TƏXMİNİ nəticə verir (xəta üzün əyriliyi ilə mütənasib
    böyüyür). Bu, silinməz bir riyazi məhdudiyyətdir, SÜKUTLA
    GİZLƏDİLMİR (bax `Face.is_planar`).
  - Həcm/mərkəz: tetraedr-parçalanması (bütün 8 təpənin ortası olan
    daxili istinad nöqtəsi `p0`-dan hər üzün üçbucaqlarına qədər) —
    standart, ÜMUMİ (qabarıq, YAXUD yüngül qeyri-qabarıq) çoxüzlülər
    üçün RİYAZİ CƏHƏTDƏN DÜZGÜN üsul (bax `HexahedralCell.volume`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from .validation import ValidationResult

#: Standart hüceyrə təpə sırası (VTK_HEXAHEDRON / Eclipse corner-point
#: konvensiyasına uyğun): aşağı üz (Z-) düz sırayla, sonra yuxarı üz (Z+)
#: EYNİ sırayla, birbaşa üstündə.
#:
#:     v3────v2      v7────v6
#:     │ alt │       │ üst │
#:     v0────v1      v4────v5
#:
#: Hər üzün təpə indeksləri elə seçilib ki, sıra ilə gedəndə (sağ əl
#: qaydası) normal HƏMİŞƏ HÜCEYRƏDƏN KƏNARA (outward) baxsın — bax modul
#: səviyyəli `test`/audit doğrulaması, `tests/test_polyhedral_geometry.py`.
HEX_FACE_VERTEX_INDICES = {
    "Z-": (0, 3, 2, 1),
    "Z+": (4, 5, 6, 7),
    "Y-": (0, 1, 5, 4),
    "Y+": (2, 3, 7, 6),
    "X-": (3, 0, 4, 7),
    "X+": (1, 2, 6, 5),
}


def _triangle_area_centroid_normal(a: np.ndarray, b: np.ndarray, c: np.ndarray
                                   ) -> Tuple[float, np.ndarray, np.ndarray]:
    """Bir üçbucağın sahəsi, mərkəzi, VAHİD normalı (sağ əl qaydası,
    `(b-a) x (c-a)`)."""
    cross = np.cross(b - a, c - a)
    norm = float(np.linalg.norm(cross))
    area = 0.5 * norm
    centroid = (a + b + c) / 3.0
    normal = cross / norm if norm > 1e-30 else np.zeros(3)
    return area, centroid, normal


@dataclass
class Face:
    """Bir çoxbucaqlı üz — ƏN AZI 3 təpə, sıra ilə gedəndə (sağ əl
    qaydası) normal HÜCEYRƏDƏN KƏNARA baxmalıdır (bu, ÇAĞIRANIN
    məsuliyyətidir — bax `HexahedralCell.faces`).
    """
    vertices: np.ndarray   # (k, 3), k >= 3

    def __post_init__(self):
        self.vertices = np.asarray(self.vertices, float)
        if self.vertices.ndim != 2 or self.vertices.shape[1] != 3 or self.vertices.shape[0] < 3:
            raise ValueError(
                f"Face.vertices (k,3) olmalıdır, k>=3, alındı {self.vertices.shape}")
        self._decomposition_cache: Optional[Tuple[float, np.ndarray, np.ndarray]] = None

    def _triangles(self):
        """Fan-triangulyasiya: `(v0,v1,v2), (v0,v2,v3), ...` — bax modul
        docstring-i, "HEXAHEDRAL FƏRZİYYƏLƏR"."""
        v0 = self.vertices[0]
        for i in range(1, self.vertices.shape[0] - 1):
            yield v0, self.vertices[i], self.vertices[i + 1]

    def _decompose(self) -> Tuple[float, np.ndarray, np.ndarray]:
        """`area()`/`centroid()`/`normal()` ARASINDA paylaşılan TƏK
        üçbucaqlaşdırma keçidi — bax audit §17/§24: "avoid repeated...
        expensive polygon calculations". Nəticə keşlənir (`vertices`
        konstruksiyadan sonra DƏYİŞMİR, ona görə keş etibarlıdır)."""
        if self._decomposition_cache is not None:
            return self._decomposition_cache
        total_area = 0.0
        weighted_centroid = np.zeros(3)
        weighted_normal = np.zeros(3)
        for a, b, c in self._triangles():
            area, centroid, normal = _triangle_area_centroid_normal(a, b, c)
            total_area += area
            weighted_centroid += area * centroid
            weighted_normal += area * normal
        self._decomposition_cache = (total_area, weighted_centroid, weighted_normal)
        return self._decomposition_cache

    def area(self) -> float:
        total_area, _, _ = self._decompose()
        return float(total_area)

    def centroid(self) -> np.ndarray:
        """Sahə-çəkili üçbucaq mərkəzlərinin ortalaması — sadə təpə-
        ortalamasından FƏRQLİ olaraq, qeyri-bərabər üçbucaqlarda DÜZGÜN
        (bax audit §10: "Do not assume face centroid is always the
        midpoint... geometric coordinate")."""
        total_area, weighted_centroid, _ = self._decompose()
        if total_area <= 1e-30:
            return self.vertices.mean(axis=0)
        return weighted_centroid / total_area

    def normal(self) -> np.ndarray:
        """Sahə-çəkili üçbucaq normallarının cəmi, VAHİD vektora
        normallaşdırılıb. Tam müstəvi üz üçün bu, üçbucaqlaşdırmadan
        ASILI OLMAYAN, DƏQİQ normaldır; əyri (warped) üz üçün İKİ
        üçbucağın ORTALAMA normalıdır (bax `is_planar`)."""
        _, _, weighted_normal = self._decompose()
        norm = float(np.linalg.norm(weighted_normal))
        return weighted_normal / norm if norm > 1e-30 else np.zeros(3)

    def is_planar(self, tol: float = 1e-9) -> bool:
        """Bütün üçbucaqların normalları (demək olar) EYNİDİRSƏ üz
        müstəvidir — əks halda `area()`/`centroid()` YALNIZ TƏXMİNİDİR
        (bax modul docstring-i)."""
        normals = [n for *_ , n in
                  (_triangle_area_centroid_normal(a, b, c) for a, b, c in self._triangles())]
        if len(normals) < 2:
            return True
        reference = normals[0]
        return all(np.linalg.norm(n - reference) < tol for n in normals[1:])

    def validate(self, label: str = "Üz") -> ValidationResult:
        result = ValidationResult()
        if not np.all(np.isfinite(self.vertices)):
            result.errors.append(f"{label}: təpə koordinatlarında NaN/sonsuz var.")
            return result
        area = self.area()
        if not np.isfinite(area):
            result.errors.append(f"{label}: sahə NaN/sonsuzdur.")
        elif area <= 1e-12:
            result.errors.append(f"{label}: sahə sıfır/mənfi ({area:.3g}) — degenerativ üz.")
        normal = self.normal()
        if not np.all(np.isfinite(normal)):
            result.errors.append(f"{label}: normal NaN/sonsuzdur.")
        elif np.linalg.norm(normal) < 1e-9:
            result.errors.append(f"{label}: normal sıfır-uzunluqludur (degenerativ üz).")
        return result


def _signed_tet_volume_centroid(p0: np.ndarray, a: np.ndarray, b: np.ndarray,
                                c: np.ndarray) -> Tuple[float, np.ndarray]:
    """İşarəli tetraedr həcmi (`p0,a,b,c`) + mərkəzi. Həcm düsturu:
    `V = (1/6)·(a-p0)·[(b-p0)×(c-p0)]` — standart, dəqiq."""
    volume = float(np.dot(a - p0, np.cross(b - p0, c - p0))) / 6.0
    centroid = (p0 + a + b + c) / 4.0
    return volume, centroid


@dataclass
class HexahedralCell:
    """8-təpəli ümumi hüceyrə (ortoqonal OLMAYA da bilər) — bax modul
    docstring-i, `HEX_FACE_VERTEX_INDICES` üçün təpə/üz konvensiyası.

    `vertices` (8,3) — sıra `HEX_FACE_VERTEX_INDICES`-ə uyğun OLMALIDIR
    (bax modul-səviyyəli ASCII diaqram).
    """
    vertices: np.ndarray

    def __post_init__(self):
        self.vertices = np.asarray(self.vertices, float)
        if self.vertices.shape != (8, 3):
            raise ValueError(f"HexahedralCell.vertices (8,3) olmalıdır, alındı "
                             f"{self.vertices.shape}")
        self._faces_cache: Optional["dict[str, Face]"] = None
        self._volume_centroid_cache: Optional[Tuple[float, np.ndarray]] = None

    def faces(self) -> "dict[str, Face]":
        """6 üz, HƏR BİRİ öz sırası ilə HÜCEYRƏDƏN KƏNARA baxan normal
        verəcək şəkildə (bax `HEX_FACE_VERTEX_INDICES`). Nəticə keşlənir
        (bax audit §17/§24: "avoid repeated reconstruction of the same
        geometry") — `vertices` konstruksiyadan sonra DƏYİŞMİR."""
        if self._faces_cache is None:
            self._faces_cache = {name: Face(self.vertices[list(idx)])
                                 for name, idx in HEX_FACE_VERTEX_INDICES.items()}
        return self._faces_cache

    def _reference_point(self) -> np.ndarray:
        """Tetraedr-parçalanması üçün DAXİLİ istinad nöqtəsi — 8 təpənin
        sadə ortalaması. Qabarıq (convex) hüceyrə üçün HƏMİŞƏ daxildədir."""
        return self.vertices.mean(axis=0)

    def _tetrahedra(self):
        p0 = self._reference_point()
        for face in self.faces().values():
            for a, b, c in face._triangles():
                yield _signed_tet_volume_centroid(p0, a, b, c)

    def _volume_and_centroid(self) -> Tuple[float, np.ndarray]:
        """`volume()`/`centroid()` ARASINDA paylaşılan TƏK tetraedr-
        parçalanması keçidi (bax `Face._decompose`-un EYNİ məntiqi)."""
        if self._volume_centroid_cache is not None:
            return self._volume_centroid_cache
        total_volume = 0.0
        weighted = np.zeros(3)
        for volume, centroid in self._tetrahedra():
            total_volume += volume
            weighted += volume * centroid
        self._volume_centroid_cache = (total_volume, weighted)
        return self._volume_centroid_cache

    def volume(self) -> float:
        """Bax modul docstring-i, "HEXAHEDRAL FƏRZİYYƏLƏR" — tetraedr-
        parçalanması cəmi, standart divergensiya-teoremi ekvivalenti."""
        total_volume, _ = self._volume_and_centroid()
        return float(total_volume)

    def centroid(self) -> np.ndarray:
        """Həcm-çəkili tetraedr mərkəzlərinin ortalaması — SADƏ təpə-
        ortalamasından FƏRQLİ (bax audit §9: "Do not assume (cx,cy,cz)
        = (dx/2,dy/2,dz/2) except for the actual Cartesian
        implementation" — bu, QEYRİ-bərabər paylanmış təpələr üçün
        DÜZGÜN ağırlıqlı mərkəzi verir)."""
        total_volume, weighted = self._volume_and_centroid()
        if abs(total_volume) <= 1e-30:
            return self.vertices.mean(axis=0)
        return weighted / total_volume

    def validate(self, label: str = "Hüceyrə") -> ValidationResult:
        result = ValidationResult()
        if not np.all(np.isfinite(self.vertices)):
            result.errors.append(f"{label}: təpə koordinatlarında NaN/sonsuz var.")
            return result
        volume = self.volume()
        if not np.isfinite(volume):
            result.errors.append(f"{label}: həcm NaN/sonsuzdur.")
        elif volume <= 1e-12:
            result.errors.append(
                f"{label}: həcm sıfır/mənfi ({volume:.6g}) — degenerativ və ya "
                "tərs-yönümlü (inverted) hüceyrə.")
        for name, face in self.faces().items():
            face_result = face.validate(label=f"{label} üzü {name}")
            result.errors.extend(face_result.errors)
            result.warnings.extend(face_result.warnings)
        return result


# ── qeyri-ortoqonallıq diaqnostikası (bax audit §12) — YALNIZ ─────────────
# DİAQNOSTİKA, TPFA-nı DƏYİŞDİRMİR/DÜZƏLTMİR.
def cell_to_cell_vector(centroid_i: np.ndarray, centroid_j: np.ndarray) -> np.ndarray:
    """`d_ij = c_j - c_i`."""
    return np.asarray(centroid_j, float) - np.asarray(centroid_i, float)


def non_orthogonality_angle(d_ij: np.ndarray, face_normal: np.ndarray) -> float:
    """`θ = arccos(|d_ij · n_f| / (|d_ij|·|n_f|))`, radian.

    `θ = 0` — tam ortoqonal (TPFA üçün İDEAL, d_ij normala paraleldir).
    `θ` böyüdükcə (MPFA-O TƏLƏB edən) qeyri-ortoqonallıq artır. YALNIZ
    DİAQNOSTİKADIR — heç bir düzəliş/kompensasiya TƏTBİQ EDİLMİR (bax
    audit §12: "Do not use them to modify TPFA yet").
    """
    d_ij = np.asarray(d_ij, float)
    face_normal = np.asarray(face_normal, float)
    d_norm = np.linalg.norm(d_ij)
    n_norm = np.linalg.norm(face_normal)
    if d_norm < 1e-30 or n_norm < 1e-30:
        return float("nan")
    cos_theta = abs(np.dot(d_ij, face_normal)) / (d_norm * n_norm)
    cos_theta = min(1.0, max(-1.0, cos_theta))   # üzən-nöqtə səhvi 1-i CÜZİ aşa bilər
    return float(np.arccos(cos_theta))
