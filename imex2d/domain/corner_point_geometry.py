"""HƏQİQİ corner-point (COORD/ZCORN) həndəsəsi — 8 təpəli hüceyrələr,
DƏQİQ həcm, ƏYRİ üz sahəsi, HƏQİQİ normal.

NƏ DƏYİŞİR
==========
Əvvəl `.GRDECL` faylının `COORD`/`ZCORN` açar sözləri OXUNUR, amma
dərhal SKALYAR `DX`/`DY`/`DZ` ortalamasına ÇEVRİLİRDİ (bax
`io/grdecl_import.py::_from_corner_point`-un köhnə forması) — yəni fay,
maili lay, əyri (skewed) hüceyrə, qeyri-konformal mesh İTİRİLİRDİ:

    həcm  = dx·dy·dz                     (qutu)
    sahə  = dy·dz / dx·dz / dx·dy        (oxa-perpendikulyar qutu üzü)
    normal = [1,0,0] / [0,1,0] / [0,0,1] (SABİT ox vektorları)

Bu modul həmin çevirməni ƏVƏZ EDİR: `COORD` (pillar oxları) və `ZCORN`
(künc dərinlikləri) hər hüceyrə üçün AÇIQ 8 təpə koordinatına açılır,
bütün ölçülər həmin təpələrdən HESABLANIR.

GERİYƏ UYĞUNLUQ (tapşırıq §3 — "simple Cartesian models still work")
====================================================================
`CornerPointGeometry` `CellGeometry`-DƏN İRSƏN GƏLİR. Yəni modeli
daşıyan bütün mövcud zəncir (`GeologicalModel` → `ReservoirModelBuilder`
→ `ReservoirModel` → TPFA → hesabat/görüntü) HEÇ BİR DƏYİŞİKLİK OLMADAN
işləyir — `volumes()`, `cell_depths()`, `face_areas(conn)`,
`face_half_distances(conn)`, `cell_centroid()`, `face_normal(conn)`,
`face_centroid(conn)` EYNİ İMZALARLA qalır, sadəcə QUTU düsturu yerinə
HƏQİQİ həndəsəni qaytarır.

İrs olunan `dx`/`dy`/`dz`/`top_depth`/`top_depth_map` sahələri SAXLANILIR
(TƏMSİLEDİCİ ortalamalar kimi) — çünki modul-səviyyəli köməkçi
funksiyalar (`geometry.xy_to_ij`, `layer_edges`, `depth_to_k` — quyu
yerləşdirmə və dərinlik→lay uyğunlaşdırması) hələ də onlardan istifadə
edir. BU SAHƏLƏR ARTIQ HEÇ BİR HƏCM/SAHƏ/NORMAL HESABINDA İŞTİRAK
ETMİR; onlar YALNIZ həmin köməkçilərin TƏXMİNİ (nominal grid ölçüsü)
girişidir və bu, `CornerPointGeometry.approximation_notes()`-da AÇIQ
bildirilir — sükutla gizlədilmir.

Kartezian model CPG-nin XÜSUSİ HALIDIR: `CornerPointGeometry.
from_cartesian(cell_geometry)` istənilən `CellGeometry`-ni eyni
təpə-əsaslı təmsilə çevirir və nəticələr qutu düsturu ilə MAŞIN
DƏQİQLİYİNDƏ üst-üstə düşür (bax `tests/test_corner_point_geometry.py`).

RİYAZİYYAT
==========
Bütün düsturlar `polyhedral_geometry.py`-dəki `Face`/`HexahedralCell`
nüvəsi ilə EYNİDİR — bu modul onların VEKTORLAŞDIRILMIŞ (bütün grid
üçün tək numpy keçidi, hüceyrə üzərində Python dövrü YOX) formasıdır.
İkisinin uyğunluğu testlə qıfıllanıb.

  · Təpə rekonstruksiyası — pillar boyunca XƏTTİ interpolyasiya:
        t = (z − z_top) / (z_bot − z_top)
        (x, y) = (x_top, y_top) + t·((x_bot, y_bot) − (x_top, y_top))
    Hüceyrənin dörd pillar-küncü birlikdə BİLİNEAR ayaq izi verir.

  · Üz (A,B,C,D) — MƏRKƏZ-fan triangulyasiyası (`c0 = ¼(A+B+C+D)`):
        S      = Σ ½‖(v_i − c0) × (v_{i+1} − c0)‖         (skalyar sahə)
        S⃗      = Σ ½ (v_i − c0) × (v_{i+1} − c0) ≡ ½(C−A)×(D−B)
        n̂      = S⃗ / ‖S⃗‖                                  (vahid normal)
    `S⃗`-nin diaqonal düsturuna BƏRABƏRLİYİ eynilikdir (bax
    `polyhedral_geometry.py` modul docstring-i, "Normal barədə") — yəni
    bu, klassik `(C−A)×(D−B)` normalıdır, sadəcə owner/neighbor
    arasında SİMMETRİK yolla alınıb.

  · Həcm — hüceyrə mərkəzi `P_c = ⅛Σp_i` və ÜZ MƏRKƏZLƏRİ ilə
    tetraedr parçalanması (6 üz × 4 kənar = 24 tetraedr):
        V_tet = ⅙ (c0 − P_c) · [(v_i − P_c) × (v_{i+1} − P_c)]
        V     = Σ V_tet
    ƏYRİ üzlərdə də qonşu hüceyrələr EYNİ üz parçalanmasını görür, ona
    görə hüceyrələr boşluqsuz/üst-üstə düşmədən DÖŞƏNİR.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from .geometry import CellGeometry, WellblockGeometry, completion_axis
from .grid import CartesianGrid, Connections
from .polyhedral_geometry import HEX_FACE_VERTEX_INDICES
from .validation import validate_cell_volumes, validate_grid_dimensions

#: Sıfıra bölmə mühafizəsi — koordinatlar metrdədir, bu hədd fiziki
#: cəhətdən mənalı heç bir uzunluqdan kiçikdir.
_EPS = 1e-12

#: `Connections.axis` (0=X, 1=Y, 2=Z) -> `cell_a`-nın həmin əlaqəni
#: daşıyan MÜSBƏT yerli üzü. `CartesianGrid.build_connections`-da
#: `cell_a` HƏMİŞƏ aşağı indeksdir, ona görə paylaşılan fiziki üz
#: HƏMİŞƏ `cell_a`-nın "+" üzüdür (və `cell_b`-nin "−" üzüdür).
_AXIS_POSITIVE_FACE = {0: "X+", 1: "Y+", 2: "Z+"}

#: `(b, c)` = (j-istiqamətli künc, i-istiqamətli künc) -> üz-daxili yerli
#: təpə indeksi. Bax `HEX_FACE_VERTEX_INDICES` diaqramı: v0=(i−,j−),
#: v1=(i+,j−), v2=(i+,j+), v3=(i−,j+).
_INPLANE_LOCAL_INDEX = {(0, 0): 0, (0, 1): 1, (1, 1): 2, (1, 0): 3}

#: Ayaq izinin sarğısını (winding) TƏRSİNƏ çevirən permutasiya — sol-əlli
#: (left-handed) deck-lərdə mənfi həcmi düzəltmək üçün, bax
#: `corner_point_nodes(..., fix_orientation=True)`.
_REVERSE_WINDING = (0, 3, 2, 1, 4, 7, 6, 5)


# ═══════════════════════════════════════════ vektorlaşdırılmış nüvə
def quad_metrics(quads: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Dördbucaqlı (ümumiyyətlə istənilən k-bucaqlı) üzlərin sahəsi,
    mərkəzi və SAHƏ-VEKTORU — mərkəz-fan triangulyasiyası ilə.

    `quads` — `(..., k, 3)`; qaytarır `(area (...,), centroid (...,3),
    area_vector (...,3))`.

    `area` üçbucaq sahələrinin SKALYAR cəmidir (əyri səthin həqiqi
    sahəsi), `area_vector` isə VEKTOR cəmidir (`≡ ½(C−A)×(D−B)`);
    müstəvi üzdə `‖area_vector‖ == area`, əyri üzdə KİÇİKDİR — bu fərq
    üzün əyriliyinin birbaşa ölçüsüdür.
    """
    quads = np.asarray(quads, float)
    centre = quads.mean(axis=-2)                       # (..., 3)
    tail = quads - centre[..., None, :]                # (..., k, 3)
    head = np.roll(tail, -1, axis=-2)
    cross = np.cross(tail, head)                       # (..., k, 3)

    tri_area = 0.5 * np.linalg.norm(cross, axis=-1)    # (..., k)
    tri_centroid = (centre[..., None, :] + quads + np.roll(quads, -1, axis=-2)) / 3.0

    area = tri_area.sum(axis=-1)
    area_vector = 0.5 * cross.sum(axis=-2)
    weighted = (tri_area[..., None] * tri_centroid).sum(axis=-2)
    # Dejenerativ (sıfır sahəli) üzdə mərkəz təpə-ortalamasına düşür —
    # NaN yaymaqdansa həndəsi cəhətdən mənalı yeganə cavab budur.
    safe = np.where(area > _EPS, area, 1.0)
    centroid = np.where((area > _EPS)[..., None], weighted / safe[..., None], centre)
    return area, centroid, area_vector


def _point_in_quads(point: np.ndarray, quads: np.ndarray) -> np.ndarray:
    """`(n,)` bool — 2D `point` hər dördbucağın (`(n, 4, 2)`) İÇİNDƏDİRMİ.

    Hər dördbucaq İKİ üçbucağa (`v0,v1,v2` və `v0,v2,v3`) bölünür və
    nöqtə barisentrik işarə testindən keçirilir. Bu, sarğı istiqamətindən
    (saat əqrəbi ilə / əksinə) ASILI DEYİL — hər iki halda işləyir, ona
    görə sol-əlli deck-lərdə də düzgün nəticə verir. Sərhəd üzərindəki
    nöqtə İÇƏRİDƏ sayılır (`>= 0`), yəni qonşu sütunların ortaq kənarına
    düşən quyu "heç bir hüceyrədə deyil" vəziyyətinə düşmür.
    """
    def _cross2(u, v):
        """2D skalyar vektor hasili. `np.cross` NumPy 2.0-da 2 ölçülü
        vektorlar üçün SİLİNİB, ona görə açıq yazılır."""
        return u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]

    def _in_triangle(a, b, c):
        # İşarəli sahələr; hamısı eyni işarədədirsə (və ya sıfırdırsa)
        # nöqtə üçbucağın içindədir.
        d1 = _cross2(b - a, point - a)
        d2 = _cross2(c - b, point - b)
        d3 = _cross2(a - c, point - c)
        non_negative = (d1 >= 0) & (d2 >= 0) & (d3 >= 0)
        non_positive = (d1 <= 0) & (d2 <= 0) & (d3 <= 0)
        return non_negative | non_positive

    quads = np.asarray(quads, float)
    v0, v1, v2, v3 = quads[:, 0], quads[:, 1], quads[:, 2], quads[:, 3]
    return _in_triangle(v0, v1, v2) | _in_triangle(v0, v2, v3)


def unit_normals(area_vectors: np.ndarray) -> np.ndarray:
    """Sahə-vektorlarını VAHİD normala çevirir; dejenerativ üzdə sıfır
    vektor (SÜKUTLA "1"-ə normallaşdırılmır — bax `validate`)."""
    area_vectors = np.asarray(area_vectors, float)
    norm = np.linalg.norm(area_vectors, axis=-1, keepdims=True)
    return np.divide(area_vectors, norm, out=np.zeros_like(area_vectors),
                     where=norm > _EPS)


def hex_metrics(nodes: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Hekzahedral hüceyrələrin DƏQİQ həcmi və HƏCM-ÇƏKİLİ mərkəzi.

    `nodes` — `(ncell, 8, 3)`, `HEX_FACE_VERTEX_INDICES` konvensiyası ilə.
    Qaytarır `(volumes (ncell,), centroids (ncell, 3))`.

    Parçalanma: hüceyrə mərkəzi `P_c = ⅛Σp_i` → hər üzün mərkəz-fan
    üçbucaqları → 6×4 = 24 tetraedr (bax modul docstring-i,
    "RİYAZİYYAT"). `HexahedralCell.volume()`/`.centroid()` ilə EYNİ
    ədədi verir, sadəcə hüceyrə üzərində Python dövrü OLMADAN.
    """
    nodes = np.asarray(nodes, float)
    if nodes.ndim != 3 or nodes.shape[1:] != (8, 3):
        raise ValueError(f"nodes (ncell,8,3) olmalıdır, alındı {nodes.shape}")

    apex = nodes.mean(axis=1)                          # (ncell, 3) = P_c
    total_volume = np.zeros(nodes.shape[0])
    weighted = np.zeros((nodes.shape[0], 3))

    for indices in HEX_FACE_VERTEX_INDICES.values():
        quad = nodes[:, list(indices), :]              # (ncell, 4, 3)
        face_centre = quad.mean(axis=1)                # (ncell, 3) = c0
        v_i = quad - apex[:, None, :]                  # (ncell, 4, 3)
        v_next = np.roll(quad, -1, axis=1) - apex[:, None, :]
        # V_tet = ⅙ (c0 − P_c) · [(v_i − P_c) × (v_{i+1} − P_c)]
        triple = np.einsum("cj,ckj->ck", face_centre - apex, np.cross(v_i, v_next))
        tet_volume = triple / 6.0                      # (ncell, 4)
        tet_centroid = (apex[:, None, :] + face_centre[:, None, :]
                        + quad + np.roll(quad, -1, axis=1)) / 4.0
        total_volume += tet_volume.sum(axis=1)
        weighted += np.einsum("ck,ckj->cj", tet_volume, tet_centroid)

    safe = np.where(np.abs(total_volume) > _EPS, total_volume, 1.0)
    centroid = np.where(np.abs(total_volume)[:, None] > _EPS,
                        weighted / safe[:, None], apex)
    return total_volume, centroid


# ═══════════════════════════════════════════ COORD/ZCORN -> 8 təpə
def corner_point_nodes(nx: int, ny: int, nz: int, coord: np.ndarray,
                       zcorn: np.ndarray, fix_orientation: bool = True
                       ) -> Tuple[np.ndarray, Dict[str, object]]:
    """`COORD`/`ZCORN` → `(ncell, 8, 3)` AÇIQ təpə koordinatları.

    Girişlərin Eclipse yerləşimi
    ----------------------------
    `COORD`  — `6·(nx+1)·(ny+1)`, j-major: hər pillar üçün
               `(x_top, y_top, z_top, x_bot, y_bot, z_bot)`.
    `ZCORN`  — `8·nx·ny·nz`, `[k, üst/alt, j, j−/j+, i, i−/i+]` sırası.

    Çıxışın təpə sırası `HEX_FACE_VERTEX_INDICES` konvensiyasıdır:
    `v0..v3` = TAVAN müstəvisi (kiçik z = dayaz), `v4..v7` = DABAN
    müstəvisi, hər ikisi `(i−,j−) → (i+,j−) → (i+,j+) → (i−,j+)`
    sırası ilə. Hüceyrə indeksi `(k·ny + j)·nx + i` — `CartesianGrid`
    ilə EYNİ (bax `CartesianGrid.index`).

    `fix_orientation` — deck sol-əlli (left-handed) koordinat sistemində
    yazılıbsa (məs. j şimala yox, cənuba artırsa) BÜTÜN həcmlər mənfi
    çıxar. Bu halda ayaq izinin sarğısı BİR DƏFƏ, BÜTÜN grid üçün
    tərsinə çevrilir — hüceyrə-hüceyrə "mütləq qiymət almaq" DEYİL (o,
    həqiqətən tərs-yönümlü tək hüceyrələri GİZLƏDƏRDİ). Nə edildiyi
    qaytarılan hesabatda (`notes`) AÇIQ yazılır.

    Qaytarır `(nodes, notes)` — `notes` diaqnostika lüğətidir:
    `flipped_orientation`, `degenerate_pillars`, `negative_volume_cells`,
    `collapsed_cells` (sıfıra yaxın həcmli — Eclipse-in "pinch-out"
    layları; bunlar XƏTA DEYİL, ÇAĞIRAN qərar verir).
    """
    nx, ny, nz = int(nx), int(ny), int(nz)
    coord = np.asarray(coord, float).ravel()
    zcorn = np.asarray(zcorn, float).ravel()

    expected_coord = 6 * (nx + 1) * (ny + 1)
    expected_zcorn = 8 * nx * ny * nz
    if coord.size != expected_coord:
        raise ValueError(f"COORD ölçüsü {coord.size}, gözlənilən "
                         f"{expected_coord} (6·(nx+1)·(ny+1)).")
    if zcorn.size != expected_zcorn:
        raise ValueError(f"ZCORN ölçüsü {zcorn.size}, gözlənilən "
                         f"{expected_zcorn} (8·nx·ny·nz).")

    pillars = coord.reshape(ny + 1, nx + 1, 6)
    top_xy = pillars[:, :, 0:2]                        # (ny+1, nx+1, 2)
    top_z = pillars[:, :, 2]
    bot_xy = pillars[:, :, 3:5]
    bot_z = pillars[:, :, 5]
    span = bot_z - top_z                               # pillar şaquli uzanması

    # Dejenerativ pillar (sıfır şaquli uzanma): `t` təyin olunmur. Belə
    # pillarda x/y onsuz da dərinlikdən ASILI DEYİL, ona görə TAVAN
    # nöqtəsi götürülür (t = 0) — bu, uydurma DEYİL, həmin pillar üçün
    # yeganə mövcud koordinatdır. Sayı hesabatda bildirilir.
    degenerate = np.abs(span) <= _EPS
    safe_span = np.where(degenerate, 1.0, span)

    # ZCORN -> [k, üst/alt, j, j−/j+, i, i−/i+]
    z = zcorn.reshape(nz, 2, ny, 2, nx, 2)

    nodes = np.empty((nz, ny, nx, 8, 3), dtype=float)
    for (corner_j, corner_i), local in _INPLANE_LOCAL_INDEX.items():
        # Bu ayaq-izi küncünə uyğun pillar altmassivi — (ny, nx)
        sl = (slice(corner_j, corner_j + ny), slice(corner_i, corner_i + nx))
        p_top_xy = top_xy[sl]                          # (ny, nx, 2)
        p_bot_xy = bot_xy[sl]
        p_top_z = top_z[sl]                            # (ny, nx)
        p_span = safe_span[sl]

        for plane in (0, 1):                           # 0 = tavan, 1 = daban
            depth = z[:, plane, :, corner_j, :, corner_i]         # (nz, ny, nx)
            t = (depth - p_top_z[None, :, :]) / p_span[None, :, :]
            t = np.where(degenerate[sl][None, :, :], 0.0, t)
            xy = (p_top_xy[None, :, :, :]
                  + t[..., None] * (p_bot_xy - p_top_xy)[None, :, :, :])
            node = 4 * plane + local
            nodes[:, :, :, node, 0:2] = xy
            nodes[:, :, :, node, 2] = depth

    nodes = nodes.reshape(nx * ny * nz, 8, 3)

    notes: Dict[str, object] = {
        "flipped_orientation": False,
        "degenerate_pillars": int(np.count_nonzero(degenerate)),
    }

    volumes, _ = hex_metrics(nodes)
    if fix_orientation and volumes.size and np.median(volumes) < 0.0:
        # BÜTÜN grid sistematik olaraq tərsdir (sol-əlli deck) — bir
        # dəfəlik, QLOBAL sarğı düzəlişi.
        nodes = nodes[:, list(_REVERSE_WINDING), :]
        volumes = -volumes
        notes["flipped_orientation"] = True

    notes["negative_volume_cells"] = int(np.count_nonzero(volumes < -_EPS))
    notes["collapsed_cells"] = int(np.count_nonzero(np.abs(volumes) <= _EPS))
    return nodes, notes


def cartesian_nodes(grid: CartesianGrid, geometry: CellGeometry) -> np.ndarray:
    """Kartezian `CellGeometry` → `(ncell, 8, 3)` təpələr.

    Kartezian grid CPG-nin XÜSUSİ HALIDIR (tapşırıq §3) — bu funksiya
    həmin xüsusi halı ÜMUMİ təmsilə çevirir, beləliklə MPFA-O və
    `CornerPointGeometry` üçün TƏK bir kod yolu qalır.

    `hexahedral_vertices_from_cartesian` (bax `general_grid_geometry.py`)
    ilə EYNİ nəticəni verir — o, tarixən MPFA-O testləri üçün yazılmış
    körpüdür və olduğu kimi saxlanılır; bu isə `top_depth_map`-i də
    (maili lay tavanı) NƏZƏRƏ ALIR, ona görə CPG konstruktoru bunu
    işlədir.
    """
    i, j, k = grid.ijk_array(np.arange(grid.ncell))
    thickness = np.asarray(geometry.dz, float)
    layer_top = np.concatenate(([0.0], np.cumsum(thickness)))

    if geometry.top_depth_map is None:
        column_top = np.full(grid.ncell, float(geometry.top_depth))
    else:
        areal = np.asarray(geometry.top_depth_map, float).ravel()
        if areal.size == grid.ncell:
            column_top = areal
        elif areal.size == grid.nx * grid.ny:
            column_top = np.tile(areal, grid.nz)
        else:
            raise ValueError("top_depth_map ölçüsü grid ilə uyğun gəlmir")

    x0, x1 = i * geometry.dx, (i + 1) * geometry.dx
    y0, y1 = j * geometry.dy, (j + 1) * geometry.dy
    z0 = column_top + layer_top[k]
    z1 = column_top + layer_top[k + 1]
    corners = (
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    )
    return np.stack([np.stack(corner, axis=-1) for corner in corners], axis=1)


# ═══════════════════════════════════════════════════════ həndəsə sinfi
@dataclass(frozen=True)
class CornerPointGeometry(CellGeometry):
    """Corner-point grid həndəsəsi — `CellGeometry` müqaviləsi, HƏQİQİ
    hesablamalar (bax modul docstring-i).

    `nodes` — `(ncell, 8, 3)`, `HEX_FACE_VERTEX_INDICES` konvensiyası.
    İrs olunan `dx`/`dy`/`dz`/`top_depth`/`top_depth_map` YALNIZ nominal
    (təmsiledici) qiymətlərdir — bax `approximation_notes()`.
    """
    nodes: Optional[np.ndarray] = None

    def __post_init__(self):
        super().__post_init__()                        # `dz`-ni (nz,) formasına salır
        if self.nodes is None:
            raise ValueError(
                "CornerPointGeometry `nodes` ((ncell,8,3)) TƏLƏB EDİR — "
                "Kartezian modeldən qurmaq üçün `from_cartesian()` işlədin.")
        nodes = np.asarray(self.nodes, float)
        if nodes.shape != (self.grid.ncell, 8, 3):
            raise ValueError(f"nodes ({self.grid.ncell},8,3) olmalıdır, alındı "
                             f"{nodes.shape}")
        object.__setattr__(self, "nodes", nodes)

        volumes, centroids = hex_metrics(nodes)
        object.__setattr__(self, "_volumes", volumes)
        object.__setattr__(self, "_centroids", centroids)
        # `Connections` obyekti başına keşlənmiş üz kəmiyyətləri — TPFA
        # hər zaman addımında `face_areas`/`face_half_distances` çağırır,
        # təkrar triangulyasiya israfdır (bax audit §17/§24).
        object.__setattr__(self, "_face_cache", {})

    # ── konstruktorlar ───────────────────────────────────────────────
    @classmethod
    def from_grdecl(cls, grid: CartesianGrid, coord: np.ndarray, zcorn: np.ndarray,
                    fix_orientation: bool = True
                    ) -> Tuple["CornerPointGeometry", Dict[str, object]]:
        """`COORD`/`ZCORN`-dan qurur; `(geometry, notes)` qaytarır.

        `notes` — `corner_point_nodes`-un diaqnostikası; çağıran
        (`GrdeclImporter`) onu `DiagnosticReport`-a yazır.
        """
        nodes, notes = corner_point_nodes(grid.nx, grid.ny, grid.nz, coord, zcorn,
                                          fix_orientation=fix_orientation)
        return cls.from_nodes(grid, nodes), notes

    @classmethod
    def from_nodes(cls, grid: CartesianGrid, nodes: np.ndarray) -> "CornerPointGeometry":
        """Hazır `(ncell,8,3)` təpələrdən — nominal `dx/dy/dz/top_depth`
        AVTOMATİK çıxarılır (bax `approximation_notes()`)."""
        nodes = np.asarray(nodes, float)
        if nodes.shape != (grid.ncell, 8, 3):
            raise ValueError(f"nodes ({grid.ncell},8,3) olmalıdır, alındı {nodes.shape}")
        cube = nodes.reshape(grid.nz, grid.ny, grid.nx, 8, 3)

        # Nominal ölçülər — YALNIZ köhnə skalyar API üçün (quyu
        # yerləşdirmə, dərinlik→lay). Həndəsi hesablamalarda İŞTİRAK
        # ETMİR.
        x_extent = cube[..., 0].max(axis=-1) - cube[..., 0].min(axis=-1)
        y_extent = cube[..., 1].max(axis=-1) - cube[..., 1].min(axis=-1)
        top_plane = cube[:, :, :, 0:4, 2].mean(axis=-1)         # (nz, ny, nx)
        bottom_plane = cube[:, :, :, 4:8, 2].mean(axis=-1)
        layer_thickness = (bottom_plane - top_plane).mean(axis=(1, 2))   # (nz,)

        column_top = top_plane[0].ravel()                        # (ny·nx,)
        return cls(grid=grid,
                   dx=float(max(x_extent.mean(), _EPS)),
                   dy=float(max(y_extent.mean(), _EPS)),
                   dz=np.maximum(layer_thickness, _EPS),
                   top_depth=float(column_top.min()),
                   top_depth_map=column_top.copy(),
                   nodes=nodes)

    @classmethod
    def from_cartesian(cls, geometry: CellGeometry) -> "CornerPointGeometry":
        """İstənilən Kartezian `CellGeometry`-ni EYNİ həndəsəni təsvir
        edən CPG-yə çevirir (tapşırıq §3 — Kartezian, CPG-nin xüsusi
        halı kimi). Həcm/sahə/mərkəz maşın dəqiqliyində üst-üstə düşür.
        """
        if isinstance(geometry, cls):
            return geometry
        nodes = cartesian_nodes(geometry.grid, geometry)
        return cls(grid=geometry.grid, dx=geometry.dx, dy=geometry.dy,
                   dz=np.asarray(geometry.dz, float).copy(),
                   top_depth=geometry.top_depth,
                   top_depth_map=(None if geometry.top_depth_map is None
                                  else np.asarray(geometry.top_depth_map, float).copy()),
                   nodes=nodes)

    # ── hüceyrə kəmiyyətləri (DƏQİQ) ─────────────────────────────────
    def volumes(self) -> np.ndarray:
        """`(ncell,)` — DƏQİQ çoxüzlü həcm, m³ (`dx·dy·dz` DEYİL)."""
        return self._volumes

    def cell_centroid(self) -> np.ndarray:
        """`(ncell, 3)` — HƏCM-ÇƏKİLİ mərkəz [X, Y, Z], m."""
        return self._centroids

    def cell_depths(self) -> np.ndarray:
        """`(ncell,)` — hüceyrə mərkəzinin dərinliyi = mərkəzin Z-i.

        Kartezian versiyada bu `top + Σdz + dz/2` idi; burada MAİLİ və
        ƏYRİ hüceyrələr üçün DƏQİQ həcm mərkəzidir."""
        return self._centroids[:, 2]

    def dz_per_cell(self) -> np.ndarray:
        """`(ncell,)` — hüceyrənin HƏQİQİ şaquli qalınlığı: daban üzünün
        mərkəzi ilə tavan üzünün mərkəzi arasındakı dərinlik fərqi
        (layın orta qalınlığı DEYİL)."""
        return (self.nodes[:, 4:8, 2].mean(axis=1)
                - self.nodes[:, 0:4, 2].mean(axis=1))

    def cell_nodes(self, cell: int) -> np.ndarray:
        """`(8, 3)` — bir hüceyrənin təpələri (görüntü/ixrac üçün)."""
        return self.nodes[int(cell)]

    # ── LOKAL ölçülər (Phase 5E) — HAMISI təpələrdən ─────────────────
    def cell_extents(self) -> np.ndarray:
        """`(ncell, 3)` — hər hüceyrənin X, Y, Z uzanması (təpələrin
        sərhəd qutusu), m.

        Bu, `dx`/`dy` NOMİNAL skalyarlarını ƏVƏZ EDİR: maili pillar,
        fay atımı və dəyişkən lay qalınlığı hər hüceyrəyə ÖZ ölçüsünü
        verir.

        DİQQƏT — bu, ox-boyu SƏRHƏD QUTUSUDUR (`max−min`), yəni
        fırlanmış/əyri hüceyrədə HƏQİQİ ölçüdən BÖYÜKDÜR (diaqonal
        uzanma). Ona görə Peaceman quyu indeksi ARTIQ bunu İŞLƏTMİR —
        o, `wellblock_geometry()`-dən keçir (bax həmin metod və
        `simulation/well_model.py`). Sərhəd qutusu qalır, çünki o,
        hüceyrənin qablaşdırma/axtarış (bounding) ölçüsü kimi hələ də
        doğru cavabdır."""
        return self.nodes.max(axis=1) - self.nodes.min(axis=1)

    def characteristic_length(self) -> np.ndarray:
        """`(ncell,)` — `V^(1/3)`, HƏQİQİ çoxüzlü həcmdən."""
        return np.cbrt(np.abs(self._volumes))

    def local_edge_vectors(self, cells=None) -> np.ndarray:
        """`(n, 3, 3)` — hər hüceyrənin YERLİ I/J/K ox vektorları, m.

        Hər ox, hüceyrənin həmin istiqamətdəki DÖRD kənarının
        ortalamasıdır; bu, eyni zamanda QARŞI ÜZLƏRİN MƏRKƏZLƏRİ
        arasındakı vektordur (üz mərkəzi 4 təpənin ortalamasıdır), yəni
        hüceyrənin həmin ox boyunca HƏQİQİ uzanması və istiqamətidir.
        Ox-boyu sərhəd qutusundan (`cell_extents`) fərqi budur: fırlanma
        və ya kəsilmə (shear) bu vektorun UZUNLUĞUNU şişirtmir, yalnız
        onu döndərir.
        """
        nodes = self.nodes if cells is None else self.nodes[
            np.atleast_1d(np.asarray(cells, dtype=int))]
        i_axis = (nodes[:, 1] - nodes[:, 0] + nodes[:, 2] - nodes[:, 3]
                  + nodes[:, 5] - nodes[:, 4] + nodes[:, 6] - nodes[:, 7]) / 4.0
        j_axis = (nodes[:, 3] - nodes[:, 0] + nodes[:, 2] - nodes[:, 1]
                  + nodes[:, 7] - nodes[:, 4] + nodes[:, 6] - nodes[:, 5]) / 4.0
        k_axis = (nodes[:, 4] - nodes[:, 0] + nodes[:, 5] - nodes[:, 1]
                  + nodes[:, 6] - nodes[:, 2] + nodes[:, 7] - nodes[:, 3]) / 4.0
        return np.stack([i_axis, j_axis, k_axis], axis=1)

    def wellblock_geometry(self, cells, direction="Z") -> WellblockGeometry:
        """HƏQİQİ corner-point wellblock həndəsəsi — Peaceman girişi.

        NİYƏ SƏRHƏD QUTUSU DEYİL
        ------------------------
        Əvvəl Peaceman `cell_extents()` (yəni `max(x)−min(x)`,
        `max(y)−min(y)`) oxuyurdu. Fırlanmış hüceyrədə bu ölçü hüceyrənin
        ÖZ eninə deyil, onun DİAQONALINA yaxınlaşır: 45° fırlanmış
        100×100 m blok üçün sərhəd qutusu 141×141 m verir, yəni `r_e`
        40 % şişir və WI süni olaraq ~8 % azalır — həndəsə isə heç
        dəyişməyib, sadəcə döndərilib. Kəsilmiş (skewed) hüceyrədə səhv
        daha böyükdür, çünki qutu HƏM şişir, HƏM də həqiqi sahəni itirir.

        HESABLAMA
        ---------
        1. Quyu oxu — hüceyrənin YERLİ oxu (`direction`): "Z" pillar
           istiqamətidir, "X"/"Y" isə I/J istiqaməti (üfüqi quyu).
           `length = ‖yerli ox‖` — perforasiyanın hüceyrə içindəki
           HƏQİQİ uzunluğu (dəyişkən qalınlıqda lay-lay fərqlidir).
        2. En kəsik sahəsi — `A = V / length`, burada `V` HƏQİQİ
           çoxüzlü həcmdir. Prizma üçün bu eynilikdir; maili/kəsilmiş
           hüceyrədə o, ayaq izinin quyu oxuna PROYEKSİYASINI verir —
           məhz Peaceman-ın tələb etdiyi kəmiyyət.
        3. Effektiv ölçülər — qalan iki yerli ox quyu oxuna
           perpendikulyar müstəviyə proyeksiya edilir; `d1/d2` onların
           uzunluqlarının NİSBƏTİNİ, `d1·d2` isə HƏQİQİ sahəni (`A`)
           saxlayır:
               d1 = √(A·‖e1⊥‖/‖e2⊥‖),   d2 = √(A·‖e2⊥‖/‖e1⊥‖)
           Beləliklə fırlanma WI-ni DƏYİŞMİR (ölçülər fırlanır, uzunluq
           qalır), kəsilmə isə həqiqi sahə qədər AZALDIR.
        4. `axis1`/`axis2` — həmin ölçülərin vahid istiqamətləri;
           anizotrop Peaceman onlara `Kx/Ky/Kz`-ni proyeksiya edir
           (bax `WellblockGeometry.directional_permeability`).

        Kartezian blokda addımların HAMISI eyniliklə keçir
        (`A = dx·dy`, `e1 = x̂·dx`, `e2 = ŷ·dy`) — nəticə maşın
        dəqiqliyində köhnə `dx`/`dy`/`dz_k` ilə üst-üstə düşür.
        """
        axis = completion_axis(direction)
        cells = np.atleast_1d(np.asarray(cells, dtype=int))
        edges = self.local_edge_vectors(cells)                  # (n, 3, 3)

        well = edges[:, axis, :]
        length = np.linalg.norm(well, axis=1)
        well_axis = well / np.maximum(length, _EPS)[:, None]

        # HƏQİQİ en kəsik sahəsi — nominal ayaq izi yox.
        area = np.abs(self._volumes[cells]) / np.maximum(length, _EPS)

        first, second = [a for a in (0, 1, 2) if a != axis]
        perpendicular = []
        for local in (first, second):
            edge = edges[:, local, :]
            along = np.einsum("ij,ij->i", edge, well_axis)
            perpendicular.append(edge - along[:, None] * well_axis)
        e1, e2 = perpendicular
        l1 = np.maximum(np.linalg.norm(e1, axis=1), _EPS)
        l2 = np.maximum(np.linalg.norm(e2, axis=1), _EPS)

        return WellblockGeometry(
            length=length,
            d1=np.sqrt(area * l1 / l2),
            d2=np.sqrt(area * l2 / l1),
            axis1=e1 / l1[:, None],
            axis2=e2 / l2[:, None],
            well_axis=well_axis)

    # ── quyu/dərinlik axtarışı — HƏQİQİ ayaq izi ─────────────────────
    def _column_footprints(self) -> np.ndarray:
        """`(ny, nx, 4, 2)` — hər sütunun ÜST layının tavan üzünün X/Y
        ayaq izi (dördbucaq, `HEX_FACE_VERTEX_INDICES["Z-"]` sırası ilə).

        Keşlənir: quyu axtarışı model başına bir neçə dəfə çağırılır,
        amma massiv `nodes`-dan asılıdır və o, DƏYİŞMİR."""
        cached = getattr(self, "_footprint_cache", None)
        if cached is None:
            top_layer = self.nodes[:self.grid.nx * self.grid.ny]
            quad = top_layer[:, list(HEX_FACE_VERTEX_INDICES["Z-"]), 0:2]
            cached = quad.reshape(self.grid.ny, self.grid.nx, 4, 2)
            object.__setattr__(self, "_footprint_cache", cached)
        return cached

    def locate_column(self, x: float, y: float) -> tuple:
        """`(x, y)` → `(i, j)` — HƏQİQİ ayaq izi üzərində axtarış.

        Nominal `int(x/dx)` bölməsi corner-point-də YANLIŞ hüceyrə verir:
        pillarlar maili ola bilər, fay atımı sütunları sürüşdürür, ayaq
        izi düzbucaqlı olmaya bilər. Ona görə nöqtə hər sütunun tavan
        dördbucağı ilə YOXLANILIR (sarğı işarəsi testi — dördbucaq
        qabarıq olmasa da işləyir, çünki iki üçbucağa bölünür).

        Nöqtə HEÇ BİR sütuna düşmürsə (grid-dən kənar, yaxud sütunlar
        arası mikroskopik boşluq) ƏN YAXIN sütun mərkəzi qaytarılır —
        `CellGeometry` ilə eyni müqavilə: bu metod HEÇ VAXT xəta atmır,
        "sərhəddən kənar" `validate_wells`-də AYRICA bildirilir.

        MƏHDUDİYYƏT (gizlədilmir): fay atımı areal olaraq ÜST-ÜSTƏ DÜŞƏN
        sütunlar yarada bilər; belə nöqtə üçün ƏN KİÇİK indeksli sütun
        qaytarılır. Bu, DETERMİNİSTİKDİR (təsadüfi deyil), amma "hansı
        fay bloku" sualını HƏLL ETMİR — onun üçün quyu trayektoriyası
        lazımdır, bu isə bu qatın məlumatı deyil.

        Axtarış ÜST layın tavan üzü üzərindədir: sütun boyu ayaq izi
        dərinliklə dəyişə bilər (maili pillar), ona görə çox dərin
        perforasiya üçün nəticə təxminidir.
        """
        point = np.array([float(x), float(y)])
        quads = self._column_footprints()                  # (ny, nx, 4, 2)
        inside = _point_in_quads(point, quads.reshape(-1, 4, 2))
        hits = np.flatnonzero(inside)
        if hits.size:
            flat = int(hits[0])
        else:
            centres = quads.mean(axis=2).reshape(-1, 2)    # (ny·nx, 2)
            flat = int(np.argmin(((centres - point) ** 2).sum(axis=1)))
        return int(flat % self.grid.nx), int(flat // self.grid.nx)

    def column_layer_edges(self, i: int, j: int) -> np.ndarray:
        """`(nz+1,)` — `(i, j)` sütununda HƏQİQİ lay sərhədləri.

        Hər sərhəd həmin layın tavan/daban üzünün MƏRKƏZ dərinliyidir,
        "tavan + kumulyativ orta `dz`" DEYİL — corner-point-də lay
        qalınlığı sütundan sütuna dəyişir və laylar maili ola bilər.
        """
        grid = self.grid
        cells = grid.index(int(i), int(j), 0) + np.arange(grid.nz) * grid.nx * grid.ny
        column = self.nodes[cells]                         # (nz, 8, 3)
        tops = column[:, 0:4, 2].mean(axis=1)
        bottoms = column[:, 4:8, 2].mean(axis=1)
        return np.concatenate([tops, bottoms[-1:]])

    # ── üz kəmiyyətləri (DƏQİQ) ──────────────────────────────────────
    def _faces(self, conn: Connections) -> Dict[str, np.ndarray]:
        """Bir `Connections` üçün üz sahəsi/mərkəzi/normalı — TƏK
        triangulyasiya keçidi, nəticə keşlənir.

        Paylaşılan fiziki üz HƏMİŞƏ `cell_a`-nın MÜSBƏT yerli üzüdür
        (bax `_AXIS_POSITIVE_FACE`), ona görə owner/neighbor arasında
        HEÇ BİR koordinat axtarışı LAZIM DEYİL və normal AVTOMATİK
        `cell_a → cell_b` istiqamətindədir.
        """
        # Keş açarı `id(conn)`-dur, ona görə keşlənən qeyd `conn`-un
        # ÖZÜNƏ də istinad saxlayır: əks halda `conn` zibil toplanandan
        # sonra BAŞQA obyekt eyni `id`-ni ala və KÖHNƏ üz massivlərini
        # oxuya bilərdi.
        key = id(conn)
        cached = self._face_cache.get(key)
        if cached is not None:
            return cached

        quads = np.empty((conn.count, 4, 3))
        for axis, name in _AXIS_POSITIVE_FACE.items():
            mask = conn.axis == axis
            if not np.any(mask):
                continue
            indices = list(HEX_FACE_VERTEX_INDICES[name])
            quads[mask] = self.nodes[np.asarray(conn.cell_a)[mask]][:, indices, :]

        area, centroid, area_vector = quad_metrics(quads)
        normal = unit_normals(area_vector)

        # Yarım-məsafələr: mərkəzdən üzə olan vektorun NORMAL üzərindəki
        # proyeksiyası — qeyri-ortoqonal hüceyrədə TPFA-nın tələb etdiyi
        # (və `dx/2` fərziyyəsindən fərqli) düzgün kəmiyyət budur.
        centroids = self._centroids
        half_a = np.abs(np.einsum("fj,fj->f", centroid - centroids[conn.cell_a], normal))
        half_b = np.abs(np.einsum("fj,fj->f", centroids[conn.cell_b] - centroid, normal))

        result = {"area": area, "centroid": centroid, "normal": normal,
                  "half_a": half_a, "half_b": half_b, "_conn": conn}
        self._face_cache[key] = result
        return result

    def face_areas(self, conn: Connections) -> np.ndarray:
        """`(nface,)` — ƏYRİ üzün HƏQİQİ 3D sahəsi, m² (qutu üzü DEYİL)."""
        return self._faces(conn)["area"]

    def face_normal(self, conn: Connections) -> np.ndarray:
        """`(nface, 3)` — HƏQİQİ vahid normal, `cell_a` → `cell_b`.

        SABİT `[1,0,0]`/`[0,1,0]`/`[0,0,1]` DEYİL: fay, maili lay və
        əyri hüceyrədə normal ox istiqamətindən sapır."""
        return self._faces(conn)["normal"]

    def face_centroid(self, conn: Connections) -> np.ndarray:
        """`(nface, 3)` — SAHƏ-ÇƏKİLİ üz mərkəzi."""
        return self._faces(conn)["centroid"]

    def face_half_distances(self, conn: Connections) -> Tuple[np.ndarray, np.ndarray]:
        """`(half_a, half_b)` — hər hüceyrə mərkəzindən üz müstəvisinə
        NORMAL boyunca məsafə (bax `_faces`)."""
        faces = self._faces(conn)
        return faces["half_a"], faces["half_b"]

    def areal_extent(self) -> tuple:
        """Grid-in HƏQİQİ X/Y əhatəsi (təpələrin sərhəd qutusundan)."""
        return (float(np.ptp(self.nodes[..., 0])), float(np.ptp(self.nodes[..., 1])))

    # ── MPFA-O körpüsü ───────────────────────────────────────────────
    def to_general_geometry(self, connections: Optional[Connections] = None):
        """`GeneralGridGeometry` — MPFA-O-nun (üz-əsaslı, owner/neighbor)
        gözlədiyi forma. Təpələr EYNİ massivdir, yenidən qurulmur."""
        from .general_grid_geometry import GeneralGridGeometry
        return GeneralGridGeometry(self.nodes, connections)

    # ── doğrulama / dürüstlük ────────────────────────────────────────
    def approximation_notes(self) -> List[str]:
        """İrs olunan SKALYAR sahələrin HARADA hələ də təxmini olduğunu
        AÇIQ sadalayır — sükutla gizlətməmək üçün (bax modul
        docstring-i, "GERİYƏ UYĞUNLUQ")."""
        return [
            f"`dx`≈{self.dx:.1f} m / `dy`≈{self.dy:.1f} m nominal ölçülərdir — "
            "YALNIZ `xy_to_ij()` (quyu → hüceyrə uyğunlaşdırması) işlədir. "
            "Həcm, üz sahəsi, normal və mərkəzlər BUNLARDAN ASILI DEYİL.",
            "`dz` lay üzrə ORTA qalınlıqdır — `layer_edges()`/`depth_to_k()` "
            "(dərinlik → K-lay) işlədir. Hüceyrənin öz qalınlığı üçün "
            "`dz_per_cell()` HƏQİQİ qiyməti verir.",
            "`top_depth_map` sütun tavanının ORTA dərinliyidir; hüceyrənin "
            "dörd tavan küncü fərqli dərinlikdə ola bilər (maili/faylı lay).",
        ]

    def validate(self) -> list:
        """Dejenerativ CPG həndəsəsini aşkarlayır — NaN/sonsuz təpə,
        sıfır/mənfi həcm, sıfır sahəli üz. HEÇ NƏ düzəldilmir."""
        issues = []
        grid_result = validate_grid_dimensions(self.grid.nx, self.grid.ny,
                                               self.grid.nz, self.dx, self.dy)
        issues.extend(grid_result.errors)
        if not np.all(np.isfinite(self.nodes)):
            bad = int(np.count_nonzero(~np.all(np.isfinite(self.nodes), axis=(1, 2))))
            issues.append(f"Corner-point təpələrində NaN/sonsuz var ({bad} hüceyrə).")
            return issues
        issues.extend(validate_cell_volumes(self.volumes(),
                                            label="corner-point hüceyrə həcmi").errors)
        return issues

    def quality_metrics(self, conn: Connections) -> Dict[str, float]:
        """Həndəsə keyfiyyəti — YALNIZ DİAQNOSTİKA, heç nəyi DƏYİŞMİR.

        `max_non_orthogonality_angle_deg` MPFA-O-nun nə qədər lazım
        olduğunun birbaşa ölçüsüdür: 0° = TPFA üçün ideal, böyüdükcə
        tək-nöqtəli axın xətası artır.
        """
        faces = self._faces(conn)
        d_ij = self._centroids[conn.cell_b] - self._centroids[conn.cell_a]
        distance = np.linalg.norm(d_ij, axis=1)
        cos_theta = np.abs(np.einsum("fj,fj->f", d_ij, faces["normal"]))
        cos_theta = np.divide(cos_theta, distance, out=np.zeros_like(cos_theta),
                              where=distance > _EPS)
        angle = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
        return {
            "min_cell_volume": float(self._volumes.min()),
            "max_cell_volume": float(self._volumes.max()),
            "min_face_area": float(faces["area"].min()) if conn.count else float("nan"),
            "max_face_area": float(faces["area"].max()) if conn.count else float("nan"),
            "max_non_orthogonality_angle_deg": (float(angle.max()) if conn.count
                                                else 0.0),
            "mean_non_orthogonality_angle_deg": (float(angle.mean()) if conn.count
                                                 else 0.0),
        }
