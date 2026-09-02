"""MPFA-O qarşılıqlı təsir bölgələri (interaction regions) — Phase 5A.

Bax `docs/mpfa_o_phase5a.md` §2/§3/§12 — bu modul HƏMİN spesifikasiyanın
implementasiyasıdır.

BU MODUL SAF HƏNDƏSƏ + TOPOLOGİYADIR: permeabilite, təzyiq, axın,
transmissivlik, Nyuton anlayışı YOXDUR (tapşırıq §21). Lokal riyazi
sistem `mpfa_o_local_system.py`-dədir.

Nə qurulur
----------
Hər grid TƏPƏSİ (node) üçün bir `MPFAOInteractionRegion`:

    node v
      ├── iştirakçı hüceyrələr  C(v)   (≤ 8)
      ├── sub-üzlər σ=(F,v)             (daxili node üçün 12)
      ├── hər hüceyrənin 3 sub-üzü      S(c,v)  (HƏMİŞƏ dəqiq 3)
      ├── sub-üz sahə vektorları  a_σ   (owner-dan KƏNARA)
      └── kəsilməzlik nöqtələri   x_σ = x_v + η(x_F − x_v)

`cell i` / `cell j` cütü BÖLGƏ DEYİL — bax tapşırıq §5.

Determinizm
-----------
Bölgə TAM TOPOLOJİ qurulur (`CartesianGrid` indeksləri + `HEX_FACE_
VERTEX_INDICES` konvensiyası) — koordinat müqayisəsi/floating-point
axtarış YOXDUR, eyni fəlsəfə: `GeneralGridGeometry._build_faces`. Ona
görə mürəkkəblik `O(N)`-dir (tapşırıq §31).

Struktursuz (corner-point) grid üçün bölgə qurucusu BU FAZADA YOXDUR
— bax `docs/mpfa_o_phase5a.md` §17 (Phase 5D).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from ..domain.general_grid_geometry import GeneralGridGeometry
from ..domain.grid import CartesianGrid
from ..domain.polyhedral_geometry import HEX_FACE_VERTEX_INDICES, Face

#: (di,dj,dk) oktantı -> hüceyrənin YERLİ təpə indeksi. `di=1` "hüceyrənin
#: x-böyük tərəfi" deməkdir. Bax `HEX_FACE_VERTEX_INDICES` docstring-indəki
#: ASCII diaqram: v0=(x0,y0,z0), v1=(x1,y0,z0), v2=(x1,y1,z0), v3=(x0,y1,z0),
#: v4..v7 — eyni sıra, üst üzdə.
_OCTANT_TO_LOCAL_VERTEX: Dict[Tuple[int, int, int], int] = {
    (0, 0, 0): 0, (1, 0, 0): 1, (1, 1, 0): 2, (0, 1, 0): 3,
    (0, 0, 1): 4, (1, 0, 1): 5, (1, 1, 1): 6, (0, 1, 1): 7,
}

#: YERLİ təpə indeksi -> həmin təpəni SAXLAYAN 3 yerli üz adı.
#: `HEX_FACE_VERTEX_INDICES`-dən DETERMİNİSTİK (açar sırası ilə)
#: törədilir — bax `docs/mpfa_o_phase5a.md` §12 (S(c,v) sırası).
_LOCAL_VERTEX_FACES: Dict[int, Tuple[str, ...]] = {
    vertex: tuple(name for name, idx in HEX_FACE_VERTEX_INDICES.items() if vertex in idx)
    for vertex in range(8)
}


@dataclass(frozen=True)
class MPFAOSubFace:
    """Bir sub-üz σ=(F,v) — `F` üzünün `v` təpəsinə bitişik dörddəbiri.

    Bax `docs/mpfa_o_phase5a.md` §2/§6.

    `area_vector` — `a_σ = A_σ n_σ`, HƏMİŞƏ `F`-in **owner**-indən
    KƏNARA (işarə konvensiyası §6). Qonşu tərəf üçün `-area_vector`
    işlədilir (`outward_area_vector`).
    """
    local_index: int                #: bölgə daxilində indeks
    face_index: int                 #: QLOBAL üz indeksi (GeneralGridGeometry)
    node_id: int                    #: bölgənin node identifikatoru
    owner: int                      #: QLOBAL hüceyrə indeksi
    neighbor: Optional[int]         #: `None` = sərhəd sub-üzü
    vertices: np.ndarray            #: (4,3) sub-üz poliqonu (owner-outward sırası)
    area: float
    area_vector: np.ndarray         #: (3,) = A_σ · n_σ
    centroid: np.ndarray            #: (3,) sub-üzün öz mərkəzi (diaqnostika)
    node_point: np.ndarray          #: (3,) təpə koordinatı x_v
    face_centroid: np.ndarray       #: (3,) ana üzün mərkəzi x_F
    continuity_point: np.ndarray    #: (3,) x_σ = x_v + η(x_F − x_v)

    @property
    def is_boundary(self) -> bool:
        return self.neighbor is None

    def outward_area_vector(self, cell: int) -> np.ndarray:
        """`a_σ^(c)` — `cell`-dən KƏNARA baxan sahə vektoru (§6)."""
        if cell == self.owner:
            return self.area_vector
        if cell == self.neighbor:
            return -self.area_vector
        raise ValueError(f"Hüceyrə {cell} sub-üz {self.face_index}-ə aid deyil "
                         f"(owner={self.owner}, neighbor={self.neighbor}).")

    def cells(self) -> Tuple[int, ...]:
        return (self.owner,) if self.neighbor is None else (self.owner, self.neighbor)


@dataclass
class MPFAOInteractionRegion:
    """Bir node ətrafındakı MPFA-O bölgəsi — bax `docs/mpfa_o_phase5a.md` §2.

    `cells` — iştirakçı QLOBAL hüceyrə indeksləri (artan sırada,
    deterministik). `cell_local` — qlobal -> bölgə-yerli indeks.
    `cell_sub_faces[a]` — `a` yerli hüceyrəsinin 3 sub-üzünün bölgə-yerli
    indeksləri (HƏMİŞƏ 3 ədəd, sərhəddə də — bax §2).
    """
    node_id: int
    node_ijk: Tuple[int, int, int]
    cells: List[int]
    cell_local: Dict[int, int]
    sub_faces: List[MPFAOSubFace]
    cell_sub_faces: List[List[int]]
    eta: float
    node_point: np.ndarray
    #: Bölgədə ƏN AZI bir sərhəd sub-üzü varmı (bax §10 — sərhəd
    #: bölgəsi AÇIQ fərqləndirilir, tapşırıq §20).
    is_boundary_region: bool = False
    #: İştirakçı hüceyrələrin təpə koordinatları arasındakı maksimum
    #: fərq — uyğun (conforming) grid üçün ≈0 (bax `validate`).
    node_point_mismatch: float = 0.0

    @property
    def n_cells(self) -> int:
        return len(self.cells)

    @property
    def n_sub_faces(self) -> int:
        return len(self.sub_faces)

    @property
    def interior_sub_faces(self) -> List[int]:
        return [s.local_index for s in self.sub_faces if not s.is_boundary]

    @property
    def boundary_sub_faces(self) -> List[int]:
        return [s.local_index for s in self.sub_faces if s.is_boundary]

    def closure_residual(self) -> np.ndarray:
        """`Σ_σ a_σ^(c)` HƏR hüceyrə üzrə DEYİL — bax `HexahedralCell.
        closure_residual`. Burada YALNIZ diaqnostika: bölgədəki bütün
        sub-üzlərin owner-outward sahə vektorlarının cəmi."""
        return np.sum([s.area_vector for s in self.sub_faces], axis=0)

    def describe(self) -> str:
        """İNSAN OXUYA BİLƏN dump — tapşırıq §8/§30 ("expose enough
        diagnostics to inspect one local interaction region")."""
        lines = [f"MPFAOInteractionRegion node={self.node_id} IJK={self.node_ijk} "
                 f"η={self.eta} sərhəd={self.is_boundary_region}",
                 f"  təpə x_v = {np.array2string(self.node_point, precision=6)}",
                 f"  hüceyrələr ({self.n_cells}): {self.cells}",
                 f"  sub-üzlər ({self.n_sub_faces}): "
                 f"daxili={len(self.interior_sub_faces)} sərhəd={len(self.boundary_sub_faces)}"]
        for s in self.sub_faces:
            lines.append(
                f"    σ{s.local_index}: üz={s.face_index} owner={s.owner} "
                f"neighbor={s.neighbor} A={s.area:.6g} "
                f"a_σ={np.array2string(s.area_vector, precision=6)} "
                f"x_σ={np.array2string(s.continuity_point, precision=6)}")
        for a, cell in enumerate(self.cells):
            lines.append(f"    S(c={cell}) = "
                         f"{[self.sub_faces[t].face_index for t in self.cell_sub_faces[a]]}")
        return "\n".join(lines)


def _sub_face_polygon(face_vertices: np.ndarray, corner: int,
                      face_centroid: np.ndarray) -> np.ndarray:
    """σ=(F,v) poliqonu — bax `docs/mpfa_o_phase5a.md` §2 düsturu.

    `face_vertices` (4,3) — ana üzün təpələri OWNER-OUTWARD sağ-əl
    sırasında; `corner` — `v`-nin bu sıradakı mövqeyi. Qaytarılan
    poliqon EYNİ fırlanma istiqamətindədir → normalı da owner-dan
    kənara baxır (Face-in sağ-əl qaydası).
    """
    n = face_vertices.shape[0]
    v_c = face_vertices[corner]
    v_next = face_vertices[(corner + 1) % n]
    v_prev = face_vertices[(corner - 1) % n]
    return np.array([v_c, 0.5 * (v_c + v_next), face_centroid, 0.5 * (v_prev + v_c)])


def build_interaction_regions(grid: CartesianGrid, geometry: GeneralGridGeometry,
                              eta: float = 1.0) -> List[MPFAOInteractionRegion]:
    """Bütün grid təpələri üçün MPFA-O bölgələri — `O(N)`.

    `eta` — kəsilməzlik nöqtəsi parametri η ∈ (0,1] (bax
    `docs/mpfa_o_phase5a.md` §4). DEFOLT `1.0`: x_σ = ana üzün mərkəzi,
    K-ortoqonal Kartezian gridd-də metodu DƏQİQ TPFA-ya reduksiya edən
    seçim (§16).

    `grid` TOPOLOGİYANI (hansı hüceyrə hansı node-u bölüşür), `geometry`
    isə HƏNDƏSƏNİ (təpə koordinatları, üz mərkəzləri) verir — ikisi
    `grid.ncell == geometry.ncell` ilə uzlaşmalıdır.
    """
    if not (0.0 < eta <= 1.0):
        raise ValueError(f"η ∈ (0, 1] olmalıdır, alındı {eta}")
    if grid.ncell != geometry.ncell:
        raise ValueError(f"grid.ncell ({grid.ncell}) != geometry.ncell "
                         f"({geometry.ncell}) — MPFA-O bölgələri qurula bilməz.")

    nx, ny, nz = grid.nx, grid.ny, grid.nz
    regions: List[MPFAOInteractionRegion] = []

    for kk in range(nz + 1):
        for jj in range(ny + 1):
            for ii in range(nx + 1):
                node_id = (kk * (ny + 1) + jj) * (nx + 1) + ii
                region = _build_region(grid, geometry, node_id, (ii, jj, kk), eta)
                if region is not None:
                    regions.append(region)
    return regions


def _build_region(grid: CartesianGrid, geometry: GeneralGridGeometry, node_id: int,
                  node_ijk: Tuple[int, int, int], eta: float
                  ) -> Optional[MPFAOInteractionRegion]:
    """TƏK bir node üçün bölgə (bax `build_interaction_regions`)."""
    ii, jj, kk = node_ijk

    # ── iştirakçı hüceyrələr + hər birinin YERLİ təpə indeksi ──────────
    participants: List[Tuple[int, int]] = []      # (qlobal hüceyrə, yerli təpə)
    for di in (1, 0):            # di=1 -> hüceyrə i=ii-1 (node onun x-böyük tərəfidir)
        i = ii - di
        if not (0 <= i < grid.nx):
            continue
        for dj in (1, 0):
            j = jj - dj
            if not (0 <= j < grid.ny):
                continue
            for dk in (1, 0):
                k = kk - dk
                if not (0 <= k < grid.nz):
                    continue
                cell = grid.index(i, j, k)
                participants.append((cell, _OCTANT_TO_LOCAL_VERTEX[(di, dj, dk)]))
    if not participants:
        return None

    participants.sort()
    cells = [cell for cell, _ in participants]
    cell_local = {cell: a for a, cell in enumerate(cells)}
    local_vertex_of = dict(participants)

    # ── sub-üzlər: (hüceyrə, yerli üz adı) -> QLOBAL üz -> TƏK sub-üz ──
    sub_faces: List[MPFAOSubFace] = []
    by_face: Dict[int, int] = {}                  # qlobal üz -> bölgə-yerli sub-üz indeksi
    cell_sub_faces: List[List[int]] = [[] for _ in cells]
    node_points: List[np.ndarray] = []

    for cell, local_vertex in participants:
        for local_name in _LOCAL_VERTEX_FACES[local_vertex]:
            face = geometry.face_index(cell, local_name)
            if face not in by_face:
                by_face[face] = len(sub_faces)
                sub_faces.append(_build_sub_face(geometry, face, node_id,
                                                 local_vertex_of, eta,
                                                 len(sub_faces)))
            cell_sub_faces[cell_local[cell]].append(by_face[face])
        node_points.append(geometry.cells[cell].vertices[local_vertex])

    for a, cell in enumerate(cells):
        if len(cell_sub_faces[a]) != 3:
            raise RuntimeError(
                f"Bölgə node={node_id}: hüceyrə {cell} üçün {len(cell_sub_faces[a])} "
                "sub-üz tapıldı, DƏQİQ 3 gözlənilirdi (hekzahedral fərziyyə pozulub "
                "— bax docs/mpfa_o_phase5a.md §2).")

    node_stack = np.array(node_points)
    mismatch = float(np.max(np.abs(node_stack - node_stack[0]))) if len(node_stack) > 1 else 0.0

    return MPFAOInteractionRegion(
        node_id=node_id, node_ijk=node_ijk, cells=cells, cell_local=cell_local,
        sub_faces=sub_faces, cell_sub_faces=cell_sub_faces, eta=eta,
        node_point=node_stack[0],
        is_boundary_region=any(s.is_boundary for s in sub_faces),
        node_point_mismatch=mismatch)


def _build_sub_face(geometry: GeneralGridGeometry, face: int, node_id: int,
                    local_vertex_of: Dict[int, int], eta: float,
                    local_index: int) -> MPFAOSubFace:
    grid_face = geometry.faces[face]
    owner = grid_face.owner
    face_vertices = grid_face.face.vertices              # OWNER-outward sırası
    face_centroid = grid_face.face.centroid()

    owner_vertex = local_vertex_of[owner]
    corner_indices = HEX_FACE_VERTEX_INDICES[grid_face.owner_local_name]
    corner = corner_indices.index(owner_vertex)

    polygon = _sub_face_polygon(face_vertices, corner, face_centroid)
    sub = Face(polygon)
    area = sub.area()
    node_point = face_vertices[corner]

    return MPFAOSubFace(
        local_index=local_index, face_index=face, node_id=node_id,
        owner=owner, neighbor=grid_face.neighbor,
        vertices=polygon, area=area, area_vector=area * sub.normal(),
        centroid=sub.centroid(), node_point=node_point, face_centroid=face_centroid,
        continuity_point=node_point + eta * (face_centroid - node_point))


def validate_interaction_regions(regions: List[MPFAOInteractionRegion],
                                 node_tolerance: float = 1e-9) -> List[str]:
    """Bölgələrin həndəsi uyğunluğu — HEÇ NƏ DÜZƏLDİLMİR, yalnız
    problemlərin siyahısı qaytarılır (eyni fəlsəfə: `GeneralGridGeometry.
    validate`).

    Yoxlanılanlar:
      1. `node_point_mismatch` — uyğun (conforming) grid-də iştirakçı
         hüceyrələrin həmin təpəsi EYNİ koordinatda olmalıdır. Bu
         pozulubsa MPFA-O kəsilməzlik nöqtələri fiziki cəhətdən
         mənasızdır (bax `docs/mpfa_o_phase5a.md` §2).
      2. Hər hüceyrənin DƏQİQ 3 sub-üzü (konstruksiyada da yoxlanılır).
      3. Degenerativ (sıfır sahəli) sub-üz.
    """
    issues: List[str] = []
    for region in regions:
        if region.node_point_mismatch > node_tolerance:
            issues.append(
                f"Bölgə node={region.node_id}: iştirakçı hüceyrələrin təpə "
                f"koordinatları uyğun gəlmir (fərq {region.node_point_mismatch:.3g} > "
                f"{node_tolerance:.3g}) — qeyri-uyğun (non-conforming) grid.")
        for a, cell in enumerate(region.cells):
            if len(region.cell_sub_faces[a]) != 3:
                issues.append(f"Bölgə node={region.node_id}: hüceyrə {cell} üçün "
                              f"{len(region.cell_sub_faces[a])} sub-üz (3 gözlənilir).")
        for s in region.sub_faces:
            if not np.isfinite(s.area) or s.area <= 1e-14:
                issues.append(f"Bölgə node={region.node_id}: sub-üz {s.face_index} "
                              f"degenerativdir (sahə {s.area:.3g}).")
    return issues
