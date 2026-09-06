"""GRID-səviyyəli ümumi (qeyri-ortoqonal ola bilən) həndəsə — MPFA-O üçün
son hazırlıq qatı (Phase 4: "Grid-Level General Geometry Integration").

Phase 3 `HexahedralCell`/`Face`-i (bax `polyhedral_geometry.py`) TƏK
hüceyrə səviyyəsində verirdi. Bu modul onları BİR ARAYA GƏTİRİR — çoxlu
hüceyrə, PAYLAŞILAN (deduplicated) daxili üzlər, sərhəd üzləri, mövcud
`Connections` (topologiya) ilə DETERMİNİSTİK xəritələmə.

Qat diaqramı (bax `ARCHITECTURE.md` §5.15):

    Geometry NÜVƏSİ (polyhedral_geometry.py — saf riyaziyyat, TƏK hüceyrə)
        ↓
    Grid Topology (`domain/grid.py::Connections` — HANSI hüceyrələr bağlıdır)
        ↓  cell/face ID-lər
    General Grid Geometry (BU FAYL — çoxlu hüceyrə, paylaşılan üzlər,
                            sərhəd üzləri, vektorlaşdırılmış massivlər)
        ↓
    Flux Discretization (TPFA indiki, MPFA-O gələcək)

BU MODUL DA (Phase 3 kimi) SAF HƏNDƏSƏDİR — transmissivlik/mobilite/
təzyiq/Nyuton/Jakobian YOXDUR (bax audit §21/Phase 4 hard rule).
`CellGeometry`/`TwoPointFluxDiscretization` DƏYİŞMİR, bu modul onları
İSTEHLAK ETMİR və ONLAR TƏRƏFİNDƏN İSTEHLAK EDİLMİR — TAM MÜSTƏQİL,
əlavə bir qatdır (bax audit §28: "must not change TPFA results unless
explicitly invoked").

**FAYL DEDUPLİKASİYASI (bax audit §5/§6)**: floating-point koordinat
uyğunluğuna GÜVƏNMİR. `Connections.axis` artıq HANSI iki hüceyrənin
bağlı olduğunu VƏ hansı OXDA olduğunu bilir (bax `CartesianGrid.
build_connections` — `cell_a` HƏMİŞƏ aşağı indeksdir) — ona görə hər
əlaqə DETERMİNİSTİK olaraq: `cell_a`-nın "+ox" yerli üzü === `cell_b`-nin
"-ox" yerli üzü (EYNİ FİZİKİ ÜZ). Bu, O(N) — heç bir geometrik axtarış
YOXDUR.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from .grid import CartesianGrid, Connections
from .polyhedral_geometry import (HEX_FACE_VERTEX_INDICES, Face, HexahedralCell,
                                  non_orthogonality_angle)
from .validation import ValidationResult

#: `Connections.axis` (0=X,1=Y,2=Z) -> (mənfi-ox yerli üz adı, müsbət-ox
#: yerli üz adı). `cell_a` HƏMİŞƏ aşağı indeksdir (bax `CartesianGrid.
#: build_connections`), ona görə `cell_a`-nın MÜSBƏT üzü `cell_b`-nin
#: MƏNFİ üzü ilə EYNİ fiziki üzdür.
_AXIS_LOCAL_NAMES = {0: ("X-", "X+"), 1: ("Y-", "Y+"), 2: ("Z-", "Z+")}
_OPPOSITE_LOCAL_NAME = {"X-": "X+", "X+": "X-", "Y-": "Y+", "Y+": "Y-", "Z-": "Z+", "Z+": "Z-"}


@dataclass
class GridFace:
    """Bir FİZİKİ üz — daxili (`neighbor` verilib) VƏ YA sərhəd
    (`neighbor is None`). YALNIZ TƏK `Face` saxlanılır (owner-a görə
    outward) — bax audit §5: "Do NOT create two independent geometric
    faces for every interior face"."""
    index: int
    owner: int
    owner_local_name: str
    neighbor: Optional[int]
    face: Face
    #: OWNER-OUTWARD düzəliş işarəsi (±1) — bax `GeneralGridGeometry.
    #: _orient_faces`. Təpə sırası (`HEX_FACE_VERTEX_INDICES`) düzgün
    #: olanda HƏMİŞƏ `+1`; sarğısı TƏRS olan (inverted) hüceyrədə `-1`.
    #: Həm `normal_from`, həm `area_vector_from` EYNİ işarəni işlədir —
    #: ikisi HEÇ VAXT bir-birindən ayrılmır.
    orientation_sign: float = 1.0

    @property
    def is_boundary(self) -> bool:
        return self.neighbor is None

    def normal(self) -> np.ndarray:
        """OWNER-dan KƏNARA baxan VAHİD normal."""
        return self.orientation_sign * self.face.normal()

    def area_vector(self) -> np.ndarray:
        """OWNER-dan KƏNARA baxan ORİYENTASİYALI sahə vektoru
        `A⃗ = Σ ½(b−a)×(c−a)` (bax `Face.oriented_area_vector`).

        DİQQƏT: bu, `face.area() * normal()` DEYİL — əyri (warped) üzdə
        ikisi FƏRQLİDİR və axın diskretizasiyası üçün DÜZGÜN olan
        BUDUR (bax FINDING-3)."""
        return self.orientation_sign * self.face.oriented_area_vector()

    def normal_from(self, cell: int) -> np.ndarray:
        """`cell`-ə görə OUTWARD normal — bax audit §9: "n(j,F) =
        -n(i,F)"."""
        return self._sign_for(cell) * self.normal()

    def area_vector_from(self, cell: int) -> np.ndarray:
        """`cell`-ə görə OUTWARD sahə vektoru. Daxili üz üçün
        `A⃗_owner = −A⃗_neighbor` KONSTRUKSİYA İLƏ təmin olunur — TƏK
        `Face` obyekti saxlanılır, qonşu yalnız İŞARƏNİ çevirir."""
        return self._sign_for(cell) * self.area_vector()

    def _sign_for(self, cell: int) -> float:
        if cell == self.owner:
            return 1.0
        if cell == self.neighbor:
            return -1.0
        raise ValueError(f"Hüceyrə {cell} bu üzə (owner={self.owner}, "
                         f"neighbor={self.neighbor}) aid deyil.")


class GeneralGridGeometry:
    """Çoxlu `HexahedralCell`-i idarə edən, `Connections`-la inteqrasiya
    olunan GRID-səviyyəli həndəsə (bax modul docstring-i).

    `vertices` — `(ncell, 8, 3)`, hər hüceyrənin öz `HEX_FACE_VERTEX_
    INDICES` konvensiyasına uyğun 8 təpəsi. Uniform dx/dy/dz VƏ ya
    oxa-uyğunlaşmış üz TƏLƏB OLUNMUR (bax audit §24) — istənilən
    hekzahedral təpə çoxluğu qəbul edilir.

    `connections` — İSTƏYƏ GÖRƏ (`None` = bütün üzlər sərhəddir, bax
    audit §25: "must not assume every cell... must necessarily
    participate in flow topology" — təcrid olunmuş/qeyri-aktiv hüceyrə
    dəstəyi üçün gələcək inteqrasiya nöqtəsi).
    """

    def __init__(self, vertices: np.ndarray, connections: Optional[Connections] = None):
        vertices = np.asarray(vertices, float)
        if vertices.ndim != 3 or vertices.shape[1:] != (8, 3):
            raise ValueError(
                f"vertices (ncell,8,3) olmalıdır, alındı {vertices.shape}")
        self.ncell = vertices.shape[0]
        self.connections = connections
        self.cells: List[HexahedralCell] = [HexahedralCell(vertices[c])
                                            for c in range(self.ncell)]

        self._cell_centroids = np.array([c.centroid() for c in self.cells])
        self._cell_volumes = np.array([c.volume() for c in self.cells])

        self.faces: List[GridFace] = []
        self._cell_face_indices: List[List[int]] = [[] for _ in range(self.ncell)]
        self._neighbor_lists: List[List[int]] = [[] for _ in range(self.ncell)]
        #: (hüceyrə, YERLİ üz adı) -> qlobal üz indeksi. Phase 5A ƏLAVƏSİ
        #: (bax `docs/mpfa_o_phase5a.md` §12) — MPFA-O bölgə qurucusu
        #: "hüceyrə c-nin X+ üzü hansı qlobal üzdür?" sualına O(1)
        #: cavab tələb edir. Qonşu üçün ƏKS ad qeyd olunur (EYNİ fiziki
        #: üz, bax `_OPPOSITE_LOCAL_NAME`). HEÇ BİR mövcud davranış
        #: dəyişmir — sadəcə əlavə axtarış cədvəli.
        self._cell_local_face: Dict[Tuple[int, str], int] = {}
        #: `Connections` sırası ilə hər əlaqənin qlobal üz indeksi.
        self._connection_faces: List[int] = []
        self._build_faces()
        #: Təpə sırasına KOR-KORANƏ güvənmirik (bax `_orient_faces`).
        self._reoriented_faces: List[int] = []
        self._orient_faces()

        self._face_areas = np.array([gf.face.area() for gf in self.faces])
        self._face_centroids = np.array([gf.face.centroid() for gf in self.faces])
        self._face_normals = np.array([gf.normal() for gf in self.faces])
        #: `(nface, 3)` — OWNER-outward ORİYENTASİYALI sahə vektorları.
        self._face_area_vectors = np.array([gf.area_vector() for gf in self.faces])
        self._face_owner = np.array([gf.owner for gf in self.faces], dtype=int)
        self._face_neighbor = np.array(
            [(-1 if gf.neighbor is None else gf.neighbor) for gf in self.faces], dtype=int)

    # ─────────────────────────────────────────── topologiya -> üz (O(N))
    def _build_faces(self) -> None:
        """Bax modul docstring-i, "FAYL DEDUPLİKASİYASI" — HEÇ bir
        geometrik (koordinat) axtarış YOXDUR, YALNIZ `Connections.axis`
        istifadə olunur."""
        covered: set = set()
        if self.connections is not None:
            conn = self.connections
            for k in range(conn.count):
                owner, neighbor, axis = (int(conn.cell_a[k]), int(conn.cell_b[k]),
                                         int(conn.axis[k]))
                # Sərhəd yoxlaması (bax audit §7 "invalid connectivity"):
                # MƏNFİ indeks Python-un dövri indeksləməsi ilə SƏSSİZ
                # olaraq SON hüceyrəyə bağlanardı, üstəlik `-1`
                # `face_neighbor` massivində "sərhəd" sentinelidir — yəni
                # EYNİ üz həm daxili, həm sərhəd görünərdi. Diapazondan
                # kənar müsbət indeks isə çılpaq `IndexError` verirdi.
                for cell in (owner, neighbor):
                    if not 0 <= cell < self.ncell:
                        raise ValueError(
                            f"Connections[{k}]: hüceyrə indeksi {cell} "
                            f"diapazondan kənardır (0..{self.ncell - 1}).")
                owner_name = _AXIS_LOCAL_NAMES[axis][1]        # owner-un + üzü
                neighbor_name = _AXIS_LOCAL_NAMES[axis][0]     # neighbor-un - üzü (EYNİ fiziki üz)
                self._connection_faces.append(len(self.faces))
                self._add_face(owner, owner_name, neighbor)
                covered.add((owner, owner_name))
                covered.add((neighbor, neighbor_name))

        for cell in range(self.ncell):
            for name in HEX_FACE_VERTEX_INDICES:
                if (cell, name) not in covered:
                    self._add_face(cell, name, None)

    def _orient_faces(self) -> None:
        """Hər üzün sahə vektorunun OWNER-dan KƏNARA baxdığını YOXLAYIR
        (FINDING-3 §3: "Normalın orientation-ını sadəcə vertex ordering-ə
        etibar etmə").

        Meyar YERLİ və yaxşı qoyulmuşdur: `A⃗ · (x_F − c_owner) > 0`.
        Bu, KONVEKS hüceyrənin hər üzü üçün doğrudur və YALNIZ sarğısı
        tərs (inverted, mənfi həcmli) hüceyrədə pozulur — belə halda
        həmin hüceyrənin BÜTÜN 6 üzü birdən çevrilir, yəni nəticə
        özlüyündə uzlaşmış qalır. İki-hüceyrəli `A⃗ · d_ij > 0` meyarı
        BURADA çevirmə üçün İŞLƏDİLMİR (güclü qeyri-ortoqonal griddə
        qanuni şəkildə sıfıra yaxınlaşa bilər) — o, `validate()`-də
        AYRICA bildirilir.

        Düzgün `HEX_FACE_VERTEX_INDICES` sırası ilə qurulmuş hər hüceyrədə
        bu metod HEÇ NƏYİ dəyişmir (`orientation_sign` `+1` qalır), ona
        görə mövcud Kartezian/corner-point davranışı DƏYİŞMİR — sadəcə
        səssiz yanlış istiqamət artıq mümkün deyil."""
        for gf in self.faces:
            outward = float(np.dot(gf.face.oriented_area_vector(),
                                   gf.face.centroid() - self._cell_centroids[gf.owner]))
            if outward < 0.0:
                gf.orientation_sign = -1.0
                self._reoriented_faces.append(gf.index)

    def _add_face(self, owner: int, owner_local_name: str, neighbor: Optional[int]) -> None:
        face_obj = self.cells[owner].faces()[owner_local_name]
        idx = len(self.faces)
        gf = GridFace(index=idx, owner=owner, owner_local_name=owner_local_name,
                     neighbor=neighbor, face=face_obj)
        self.faces.append(gf)
        self._cell_face_indices[owner].append(idx)
        self._cell_local_face[(owner, owner_local_name)] = idx
        if neighbor is not None:
            self._cell_face_indices[neighbor].append(idx)
            self._cell_local_face[(neighbor, _OPPOSITE_LOCAL_NAME[owner_local_name])] = idx
            self._neighbor_lists[owner].append(neighbor)
            self._neighbor_lists[neighbor].append(owner)

    # ─────────────────────────────────────── vektorlaşdırılmış girişlər
    @property
    def cell_centroids(self) -> np.ndarray:
        return self._cell_centroids

    @property
    def cell_volumes(self) -> np.ndarray:
        return self._cell_volumes

    @property
    def face_areas(self) -> np.ndarray:
        return self._face_areas

    @property
    def face_centroids(self) -> np.ndarray:
        return self._face_centroids

    @property
    def face_normals(self) -> np.ndarray:
        """Owner-a görə OUTWARD (bax `GridFace.normal_from`, neighbor
        üçün İŞARƏNİ ÖZÜNÜZ dəyişdirin: `-face_normals[i]`)."""
        return self._face_normals

    @property
    def face_area_vectors(self) -> np.ndarray:
        """`(nface, 3)` — OWNER-outward ORİYENTASİYALI sahə vektorları
        `A⃗_f` (qonşu üçün İŞARƏNİ ÖZÜNÜZ dəyişdirin: `-face_area_
        vectors[i]`, və ya `GridFace.area_vector_from`).

        Bu, `face_areas[:, None] * face_normals` DEYİL — əyri (warped)
        üzdə ikisi FƏRQLİDİR (bax `Face.oriented_area_vector`,
        FINDING-3). Divergensiya/qapanma və axın diskretizasiyası üçün
        DÜZGÜN kəmiyyət BUDUR."""
        return self._face_area_vectors

    @property
    def face_owner(self) -> np.ndarray:
        return self._face_owner

    @property
    def face_neighbor(self) -> np.ndarray:
        """`-1` = sərhəd üzü (qonşu yoxdur)."""
        return self._face_neighbor

    @property
    def is_boundary(self) -> np.ndarray:
        return self._face_neighbor < 0

    # ──────────────────────────────────────────────── sorğular (O(1)/O(deg))
    def neighbors(self, cell: int) -> List[int]:
        """Bax audit §11 — YALNIZ `Connections`-dan qurulmuş, ƏVVƏLCƏDƏN
        hesablanmış siyahı, heç bir geometrik axtarış YOXDUR."""
        return list(self._neighbor_lists[cell])

    def cell_faces(self, cell: int) -> List[int]:
        return list(self._cell_face_indices[cell])

    def is_boundary_face(self, face_index: int) -> bool:
        return self.faces[face_index].is_boundary

    def face_index(self, cell: int, local_name: str) -> int:
        """`cell` hüceyrəsinin YERLİ adı `local_name` ("X+", "Z-", ...)
        olan üzünün QLOBAL indeksi — O(1) (bax `_cell_local_face`).

        Phase 5A ƏLAVƏSİ (bax `docs/mpfa_o_phase5a.md` §12). Qonşu
        tərəfdən soruşulanda EYNİ fiziki üzü qaytarır (məs. `cell_a`-nın
        "X+" üzü ilə `cell_b`-nin "X-" üzü EYNİ indeksdir) — bu, MPFA-O
        bölgəsində sub-üzlərin DEDUPLİKASİYASINI koordinat müqayisəsi
        OLMADAN təmin edir."""
        try:
            return self._cell_local_face[(cell, local_name)]
        except KeyError:
            raise KeyError(f"Hüceyrə {cell} üçün '{local_name}' yerli üzü yoxdur "
                           f"(etibarlı adlar: {sorted(HEX_FACE_VERTEX_INDICES)}).") from None

    def connection_faces(self) -> np.ndarray:
        """`Connections` sırası ilə hər əlaqəyə uyğun QLOBAL üz indeksi
        — `(nconn,)`. `connections=None` verilibsə boş massiv.

        Bu, TPFA (`Connections` üzərində işləyir) ilə MPFA-O
        (`GeneralGridGeometry` üz indeksləri üzərində işləyir)
        nəticələrini MÜQAYİSƏ etmək üçün YEGANƏ düzgün körpüdür —
        indeks sırasının təsadüfən üst-üstə düşməsinə güvənmək YOX."""
        return np.array(self._connection_faces, dtype=int)

    # ────────────────────────────────────── MPFA-yönlü həndəsi vektorlar
    def d_ij(self, face_index: int) -> np.ndarray:
        """`c_j - c_i` — YALNIZ daxili üzlər üçün mənalıdır (bax audit §22)."""
        gf = self.faces[face_index]
        if gf.is_boundary:
            raise ValueError(f"Üz {face_index} sərhəd üzüdür — d_ij tərifsizdir.")
        return self._cell_centroids[gf.neighbor] - self._cell_centroids[gf.owner]

    def d_if(self, cell: int, face_index: int) -> np.ndarray:
        """`face_centroid - cell_centroid` — `cell` bu üzün owner-i VƏ YA
        qonşusu OLMALIDIR."""
        gf = self.faces[face_index]
        if cell not in (gf.owner, gf.neighbor):
            raise ValueError(f"Hüceyrə {cell} üz {face_index}-ə aid deyil.")
        return self._face_centroids[face_index] - self._cell_centroids[cell]

    # ──────────────────────────────────────────────────────── doğrulama
    def validate(self) -> ValidationResult:
        """Bax audit §15 — hüceyrələr, üzlər, DAXİLİ üzlərin owner/
        neighbor UYĞUNLUĞU (sahə/mərkəz/əks-normal). HEÇ NƏ düzəldilmir."""
        result = ValidationResult()
        self._validate_topology(result)
        for i, cell in enumerate(self.cells):
            cell_result = cell.validate(label=f"Hüceyrə {i}")
            result.errors.extend(cell_result.errors)
            result.warnings.extend(cell_result.warnings)

        for gf in self.faces:
            face_result = gf.face.validate(label=f"Üz {gf.index} (hüceyrə {gf.owner})")
            result.errors.extend(face_result.errors)
            result.warnings.extend(face_result.warnings)
            if gf.is_boundary:
                continue
            neighbor_name = _OPPOSITE_LOCAL_NAME[gf.owner_local_name]
            neighbor_face = self.cells[gf.neighbor].faces()[neighbor_name]
            if not np.isclose(gf.face.area(), neighbor_face.area(), rtol=1e-6):
                result.errors.append(
                    f"Üz {gf.index}: owner (hüceyrə {gf.owner}) və neighbor (hüceyrə "
                    f"{gf.neighbor}) sahələri uyğun gəlmir "
                    f"({gf.face.area():.6g} vs {neighbor_face.area():.6g}).")
            if not np.allclose(gf.face.centroid(), neighbor_face.centroid(), atol=1e-6):
                result.errors.append(
                    f"Üz {gf.index}: owner/neighbor mərkəzləri uyğun gəlmir.")
            if not np.allclose(gf.face.normal(), -neighbor_face.normal(), atol=1e-6):
                result.errors.append(
                    f"Üz {gf.index}: owner/neighbor normalları əks istiqamətdə DEYİL.")
            # FINDING-3 §3: sahə vektoru owner→neighbor istiqamətində
            # OLMALIDIR. Bu, `_orient_faces`-dəki YERLİ meyardan AYRI,
            # İKİ-hüceyrəli uzlaşma yoxlamasıdır — burada yalnız
            # BİLDİRİLİR, ÇEVRİLMİR (güclü qeyri-ortoqonal griddə qalıq
            # kiçik ola bilər, avtomatik çevirmə həndəsəni GİZLİCƏ
            # korlayardı).
            projection = float(np.dot(self._face_area_vectors[gf.index],
                                      self.d_ij(gf.index)))
            if projection <= 0.0:
                result.errors.append(
                    f"Üz {gf.index}: A⃗ · (c_neighbor − c_owner) = {projection:.3g} "
                    f"≤ 0 — sahə vektoru owner {gf.owner} → neighbor "
                    f"{gf.neighbor} istiqamətində DEYİL (tərs-yönümlü və ya "
                    "həddindən artıq deformasiya olunmuş hüceyrə).")
        return result

    def _validate_topology(self, result: ValidationResult) -> None:
        """Bax audit §2/§7 — `Connections` bu sinif üçün XARİCİ girişdir və
        `CartesianGrid.build_connections()`-dən BAŞQA mənbədən də gələ
        bilər (fault/NNC, importer, test). İki pozuntu HƏNDƏSƏ
        yoxlamalarından KEÇİB GEDƏ bilirdi, ona görə burada AÇIQ tutulur:

          · `owner == neighbor` (hüceyrənin özünə əlaqəsi);
          · EYNİ (hüceyrə, yerli üz) slotunun İKİ üzə düşməsi — məsələn
            eyni hüceyrə cütü `Connections`-da İKİ dəfə verildikdə.
            Nəticədə hüceyrə 6 yerinə 7 üz alır, `_cell_local_face`
            üzərinə yazılır və yığılmış `Σ A·n` SIFIR OLMUR (qeyri-
            konservativ həndəsə) — əvvəllər bu SƏSSİZ keçirdi.

        HEÇ NƏ DÜZƏLDİLMİR — yalnız bildirilir (bax `validate()` müqaviləsi).
        """
        for face_index in self._reoriented_faces:
            gf = self.faces[face_index]
            result.errors.append(
                f"Üz {face_index}: təpə sarğısı (winding) hüceyrə {gf.owner}-ə görə "
                "TƏRSDİR — sahə vektoru/normal owner-outward olsun deyə ÇEVRİLDİ "
                "(bax `_orient_faces`). Səbəb adətən tərs-yönümlü (mənfi həcmli) "
                "hüceyrədir.")
        slot_owner: Dict[Tuple[int, str], int] = {}
        for gf in self.faces:
            if gf.neighbor is not None and gf.neighbor == gf.owner:
                result.errors.append(
                    f"Üz {gf.index}: owner == neighbor ({gf.owner}) — hüceyrənin "
                    "ÖZÜNƏ əlaqəsi etibarsızdır.")
            slots = [(gf.owner, gf.owner_local_name)]
            if gf.neighbor is not None:
                slots.append((gf.neighbor, _OPPOSITE_LOCAL_NAME[gf.owner_local_name]))
            for slot in slots:
                first = slot_owner.setdefault(slot, gf.index)
                if first != gf.index:
                    result.errors.append(
                        f"Hüceyrə {slot[0]} üçün '{slot[1]}' yerli üzü İKİ DƏFƏ "
                        f"yaradılıb (üz {first} və üz {gf.index}) — təkrarlanan "
                        "və ya ziddiyyətli `Connections` girişi.")

    # ──────────────────────────────────────────────────────── diaqnostika
    def quality_metrics(self) -> Dict[str, float]:
        """Bax audit §23 — YALNIZ DİAQNOSTİKA, heç nəyi DƏYİŞMİR."""
        interior_distances: List[float] = []
        max_angle = 0.0
        min_normal_magnitude = float("inf")
        for gf in self.faces:
            min_normal_magnitude = min(min_normal_magnitude,
                                       float(np.linalg.norm(gf.face.normal())))
            if gf.is_boundary:
                continue
            d = self.d_ij(gf.index)
            interior_distances.append(float(np.linalg.norm(d)))
            angle = non_orthogonality_angle(d, gf.face.normal())
            if np.isfinite(angle):
                max_angle = max(max_angle, angle)
        return {
            "min_cell_volume": float(self._cell_volumes.min()),
            "max_cell_volume": float(self._cell_volumes.max()),
            "min_face_area": float(self._face_areas.min()),
            "max_face_area": float(self._face_areas.max()),
            "min_centroid_distance": (float(min(interior_distances))
                                      if interior_distances else float("nan")),
            "max_centroid_distance": (float(max(interior_distances))
                                      if interior_distances else float("nan")),
            "max_non_orthogonality_angle_deg": float(np.degrees(max_angle)),
            "min_face_normal_magnitude": min_normal_magnitude,
        }


def hexahedral_vertices(grid: CartesianGrid, geometry) -> np.ndarray:
    """`(ncell, 8, 3)` təpələr — həndəsənin HANSI növ olmasından ASILI
    OLMAYARAQ (Phase 5D: corner-point inteqrasiyası).

    Bu, MPFA-O-nun (və hər hansı üz-əsaslı istehlakçının) həndəsə
    mənbəyidir:

      · `CornerPointGeometry` → HAZIR `nodes` OLDUĞU KİMİ qaytarılır
        (COORD/ZCORN-dan gələn HƏQİQİ təpələr; yenidən qurulmur,
        Kartezianlaşdırılmır);
      · adi `CellGeometry`   → qutu təpələri qurulur (Kartezian model
        CPG-nin XÜSUSİ HALIDIR, bax `corner_point_geometry.py` §"GERİYƏ
        UYĞUNLUQ") — nəticə köhnə davranışla EYNİDİR.

    Yəni MPFA-O üçün TƏK kod yolu qalır və `MPFAODiscretization.build()`
    artıq həndəsəni Kartezian saymır.
    """
    nodes = getattr(geometry, "nodes", None)
    if nodes is not None:
        nodes = np.asarray(nodes, float)
        if nodes.shape != (grid.ncell, 8, 3):
            raise ValueError(
                f"Həndəsənin `nodes` sahəsi ({grid.ncell},8,3) olmalıdır, "
                f"alındı {nodes.shape}.")
        return nodes
    from .corner_point_geometry import cartesian_nodes
    return cartesian_nodes(grid, geometry)


def hexahedral_vertices_from_cartesian(grid: CartesianGrid, geometry) -> np.ndarray:
    """`CellGeometry`-dən (Kartezian) `GeneralGridGeometry`-yə TEST/KÖRPÜ
    funksiyası — hər hüceyrənin 8 təpəsini `geometry.dx/dy/dz/top_depth`-
    dən qurur, `HEX_FACE_VERTEX_INDICES` konvensiyasına uyğun.

    YALNIZ SINAQ/DEMONSTRASİYA üçündür (bax audit §28) — `CellGeometry`-ni
    ƏVƏZ ETMİR, TPFA bunu İSTİFADƏ ETMİR/ÇAĞIRMIR. Tam vektorlaşdırılıb
    (bax audit §29 — O(N), Python dövrü hüceyrə üzərində YOXDUR).
    """
    i, j, k = grid.ijk_array(np.arange(grid.ncell))
    dz_cum = np.concatenate(([0.0], np.cumsum(geometry.dz)))
    top = geometry.top_depth
    x0, x1 = i * geometry.dx, (i + 1) * geometry.dx
    y0, y1 = j * geometry.dy, (j + 1) * geometry.dy
    z0, z1 = top + dz_cum[k], top + dz_cum[k + 1]
    corners = (
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    )
    return np.stack([np.stack(corner, axis=-1) for corner in corners], axis=1)
