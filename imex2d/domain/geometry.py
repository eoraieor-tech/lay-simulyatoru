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

    # ─────────────────────────────────────────────────────────────
    # MPFA/geometriya-abstraksiya HAZIRLIĞI (audit tapşırığı §6) —
    # `interfaces/geometry.py::IGridGeometry` müqaviləsinə uyğun,
    # AMMA bu sinif həmin ABC-ni İRSƏN ALMIR (bax izahat orada: bu
    # kod bazasında domain dataclass-ları interfeys-inject edilən
    # strategiyalar deyil, ona görə mövcud təbəqələşmə pozulmur —
    # yalnız METOD ADLARI/İMZALARI uyğunlaşdırılıb). Bu üç metod
    # TAMAMİLƏ ƏLAVƏDİR — heç bir mövcud metodun davranışını
    # DƏYİŞMİR, `TwoPointFluxDiscretization` bunlardan HƏLƏ İSTİFADƏ
    # ETMİR (TPFA öz köhnə `face_areas`/`face_half_distances` yolunu
    # saxlayır) — yalnız gələcək MPFA-O üçün lazım olan əlavə
    # həndəsi məlumatı əvvəlcədən təmin edir.
    def cell_centroid(self) -> np.ndarray:
        """(ncell, 3) — hər hüceyrənin mərkəzi [X, Y, Z], metr."""
        grid = self.grid
        i, j, _k = grid.ijk_array(np.arange(grid.ncell))
        x = (i.astype(float) + 0.5) * self.dx
        y = (j.astype(float) + 0.5) * self.dy
        z = self.cell_depths()
        return np.column_stack([x, y, z])

    def face_normal(self, conn: Connections) -> np.ndarray:
        """(nface, 3) — vahid normal, `cell_a`-dan `cell_b`-yə.

        `CartesianGrid.build_connections`-da `cell_a` HƏMİŞƏ aşağı
        indeksdir, ona görə normal HƏMİŞƏ müsbət ox istiqamətindədir.
        """
        normal = np.zeros((conn.count, 3))
        normal[conn.axis == 0, 0] = 1.0
        normal[conn.axis == 1, 1] = 1.0
        normal[conn.axis == 2, 2] = 1.0
        return normal

    def face_centroid(self, conn: Connections) -> np.ndarray:
        """(nface, 3) — üz mərkəzi.

        Oxa-perpendikulyar müstəvidə `cell_a`/`cell_b`-nin mərkəzləri
        EYNİDİR (struktur grid) — ona görə yalnız bağlayıcı ox üzrə
        `cell_a` mərkəzindən yarım-məsafə (`face_half_distances`) qədər
        sürüşdürülür. K istiqamətində qonşu təbəqələrin qalınlığı fərqli
        ola bildiyi üçün bu, sadə `(centroid_a+centroid_b)/2`
        ortalamasından DAHA DƏQİQDİR (üzün əsl mövqeyini verir).
        """
        centroid_a = self.cell_centroid()[conn.cell_a]
        half_a, _half_b = self.face_half_distances(conn)
        normal = self.face_normal(conn)
        return centroid_a + half_a[:, None] * normal

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


def column_top_depth(x: float, y: float, geometry: CellGeometry) -> float:
    """`(x, y)` sütununun TAVAN dərinliyi (m).

    `depth_to_k`-nın içindən ÇIXARILIB (davranış EYNİ) ki, `layer_edges`
    və quyu-interval → K-lay uyğunlaşdırması (bax
    `geology/layer_availability.py`) EYNİ həndəsə mənbəyini işlətsin,
    öz dərinlik modelini QURMASIN.
    """
    i, j = xy_to_ij(x, y, geometry)
    grid = geometry.grid
    if geometry.top_depth_map is None:
        return float(geometry.top_depth)
    areal = np.asarray(geometry.top_depth_map, float).ravel()
    if areal.size == grid.nx * grid.ny:
        return float(areal.reshape(grid.ny, grid.nx)[j, i])
    if areal.size == grid.ncell:
        return float(areal.reshape(grid.shape)[0, j, i])
    return float(geometry.top_depth)


def layer_edges(x: float, y: float, geometry: CellGeometry) -> np.ndarray:
    """`(x, y)` sütununda təbəqə sərhədləri — `(nz + 1,)` MÜTLƏQ dərinlik.

    `edges[k]` k-cı təbəqənin tavanı, `edges[k+1]` dabanıdır.
    """
    return column_top_depth(x, y, geometry) + np.concatenate(
        ([0.0], np.cumsum(np.asarray(geometry.dz, float))))


#: `interval_layers` üçün nisbi tolerans — layın qalınlığına görə (sərhəd
#: ÜSTÜNDƏKİ interval qonşu laya "yapışmasın", tapşırıq §23.6).
_OVERLAP_TOLERANCE = 1e-9


def interval_layers(x: float, y: float, top: float, bottom: float,
                    geometry: CellGeometry) -> list:
    """`[top, bottom]` dərinlik intervalının KƏSDİYİ K-təbəqələri (0-əsaslı).

    Kəsişmə HƏQİQİ ÖRTMƏ uzunluğu ilə hesablanır (təkcə uc nöqtələrin
    `depth_to_k`-sı ilə deyil) ki, kənar hallar DÜZGÜN işlənsin:

      · interval grid-dən TAM kənardadır      → `[]`
      · yalnız BİR layı kəsir                 → tək element
      · lay SƏRHƏDİ üzərindədir               → qonşu lay DAXİL EDİLMİR
        (sıfır qalınlıqlı örtmə "kəsir" sayılmır)
      · grid-dən yuxarı/aşağı çıxır           → yalnız örtülən hissə

    `bottom <= top` olanda `ValueError` — bu, cədvəl xətasıdır
    (`geology.validate_wells` onu ayrıca "error" kimi göstərir), SƏSSİZCƏ
    boş siyahıya çevrilmir.

    QEYD: bu, SAF HƏNDƏSƏDİR — "quyu intervalı = məlumat var" MƏNASI
    YOXDUR (bax `domain/data_availability.py` modul docstring-i).
    """
    if bottom <= top:
        raise ValueError(
            f"İnterval etibarsızdır: alt hədd ({bottom:g}) <= üst hədd ({top:g}).")
    edges = layer_edges(x, y, geometry)
    thickness = np.asarray(geometry.dz, float)
    layers = []
    for k in range(geometry.grid.nz):
        overlap = min(bottom, edges[k + 1]) - max(top, edges[k])
        if overlap > _OVERLAP_TOLERANCE * max(float(thickness[k]), 1.0):
            layers.append(k)
    return layers


def depth_to_k(x: float, y: float, depth: float, geometry: CellGeometry):
    """Dərinliyi (m) verilmiş `(x, y)` sütununda təbəqə indeksinə çevirir.

    Lay üfüqi deyil — hər `(i, j)` sütununun tavan dərinliyi
    `top_depth_map`-dan (varsa) və ya sabit `top_depth`-dən götürülür,
    sonra kumulyativ `dz` ilə təbəqə sərhədləri qurulur. Dərinlik grid-in
    diapazonundan kənardadırsa `None` qaytarır — çağıran bunu "grid
    qurulduqdan sonra" mesajı kimi göstərməlidir, xəta atılmır.
    """
    grid = geometry.grid
    edges = layer_edges(x, y, geometry)
    if depth < edges[0] or depth > edges[-1]:
        return None
    k = int(np.searchsorted(edges - edges[0], depth - edges[0], side="right")) - 1
    return int(min(max(k, 0), grid.nz - 1))
