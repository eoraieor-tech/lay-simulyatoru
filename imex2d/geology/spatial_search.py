"""Anizotrop məkan axtarışı + PEŞƏKAR QONŞULUQ SEÇİMİ (A2).

Bu modulda İKİ qat var:

1. **Aşağı qat (Phase 4.1, DƏYİŞMƏYİB)** — `AnisotropicNeighborSearch` /
   `IncrementalAnisotropicSearch`: `scipy.spatial.cKDTree` üzərində sadə
   "radius + k-ən-yaxın" sorğusu. SGS/SIS (`sgs.py`, `facies.py`) bunları
   birbaşa işlədir; müqavilələri toxunulmazdır (bax
   `tests/test_spatial_search.py` — brute-force ilə PARİTET sübutu).

2. **Yuxarı qat (A2, YENİ)** — `NeighborhoodSelector`: məqsəd "ən yaxın N
   nöqtəni tap" DEYİL, STATİSTİK/HƏNDƏSİ MƏNALI qonşuluq qurmaqdır:

       anizotrop radius axtarışı
         → istiqamətli (sektor) balanslaşdırma
         → k-ən-yaxın ehtiyat yolu
         → radiusun genişləndirilməsi
         → (yalnız açıq icazə ilə) qlobal ehtiyat
       + hər hədəf üçün DƏSTƏK TƏSNİFATI (well/boundary/weak/extrapolated)

   Hər addımın nəticəsi `NeighborhoodResult.status`-da AÇIQ görünür —
   səssiz uğursuzluq YOXDUR.

KRİTİK QAYDA (A4.3): anizotrop axtarış Kriging-in ÖZÜNÜN işlətdiyi EYNİ
transformasiya fəzasında (`geology/anisotropy.AnisotropyParams.transform`)
aparılır — əks halda seçilən qonşuluq Kriging sisteminin hesabladığı ilə
UYĞUNSUZ olar. Ona görə bu modul XAM koordinatları DEYİL, transformasiya
edilmiş koordinatları indeksləyir (və ya heç bir anizotropluq
verilməyibsə, adi izotrop Evklid fəzasını).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
from scipy.spatial import cKDTree

from .anisotropy import AnisotropyParams, transform_points

# ── qonşuluq axtarışının STATUSLARI (A2.7 — hər ehtiyat yolu görünür) ──
STATUS_RADIUS = "radius"                    #: verilmiş radiusda tapıldı
STATUS_RADIUS_EXPANDED = "radius_expanded"  #: radius genişləndirildi
STATUS_KNN = "knn"                          #: k-ən-yaxın (radius yoxdur/boşdur)
STATUS_KNN_FALLBACK = "knn_fallback"        #: radius boş qaldı → k-ən-yaxın
STATUS_GLOBAL = "global"                    #: açıq icazə ilə bütün nöqtələr
STATUS_INSUFFICIENT = "insufficient"        #: min_neighbors təmin edilmədi
STATUS_EMPTY = "empty"                      #: heç bir nöqtə yoxdur

# ── dəstək təsnifatı (A2.8 — "qonşu sayı < N" DEYİL, HƏNDƏSƏ) ─────────
SUPPORT_WELL = "well_supported"
SUPPORT_BOUNDARY = "boundary"
SUPPORT_WEAK = "weak"
SUPPORT_EXTRAPOLATED = "extrapolated"

#: `SUPPORT_*` dəyərləri "yaxşıdan pisə" sıralı — hesabat/aqreqasiya üçün.
SUPPORT_ORDER: Tuple[str, ...] = (SUPPORT_WELL, SUPPORT_BOUNDARY, SUPPORT_WEAK,
                                  SUPPORT_EXTRAPOLATED)


class NeighborhoodError(ValueError):
    """Etibarsız qonşuluq konfiqurasiyası — səssiz düzəliş EDİLMİR."""


@dataclass(frozen=True)
class NeighborhoodConfig:
    """Qonşuluq seçiminin TAM konfiqurasiyası (A2.2/A2.3/A2.5/A2.6/A2.7).

    Bütün məsafələr TRANSFORMASİYA EDİLMİŞ (anizotrop) fəzadadır —
    `AnisotropyParams.transform()`-dan sonra, yəni `range_major` ilə
    eyni vahiddə (bax `geology/anisotropy.py`). YEGANƏ istisna
    `max_vertical_distance`-dir: o, XAM Z fərqidir (A2.6, aşağıya bax).

    * `min_neighbors` / `max_neighbors` — yerli sistemin ölçü hədləri.
    * `search_radius` — `None` olanda TƏMİZ k-ən-yaxın rejim.
    * `max_search_radius` — genişləndirmənin yuxarı həddi (`None` = hədsiz).
    * `radius_expansion_factor` / `max_radius_expansions` — sparse-data
      ehtiyat yolu (A2.7); `max_radius_expansions=0` genişləndirməni söndürür.
    * `sectors` — üfüqi (transformasiya edilmiş X,Y) müstəvidə bərabər
      bucaq sektorlarının sayı. `0` = balanslaşdırma YOX (defolt, əvvəlki
      davranış). `4` = kvadrantlar, `8` = 8 üfüqi sektor.
    * `vertical_sectors` — hər bucaq sektorunu Z-nin işarəsinə görə İKİYƏ
      bölür. `sectors=4, vertical_sectors=True` → əsl 3D OKTANTLAR.
    * `max_per_sector` — bir sektordan götürüləcək maksimum nöqtə
      (`None` = `ceil(max_neighbors / n_sector_cells)`).
    * `candidate_pool_factor` — sektor balanslaşdırması işlədiləndə
      k-ən-yaxın NAMİZƏD hovuzu `max_neighbors × bu əmsal` qədər
      genişləndirilir. Bu OLMASA balanslaşdırılacaq heç nə qalmır:
      `max_neighbors` ən yaxın nöqtə onsuz da bir klasterdən gələ bilər.
      `0` (defolt) = `n_sector_cells` (yəni hər sektor xanasına orta
      hesabla `max_neighbors` namizəd düşür). Radius rejimində hovuz
      onsuz da radiusla müəyyən olunur, bu əmsal TƏSİR ETMİR.
    * `max_vertical_distance` — XAM |ΔZ| həddi (A2.6): üfüqi cəhətdən
      yaxın, amma geoloji olaraq ÇOX uzaq (başqa lay) nöqtə anizotrop
      miqyaslanmadan ASILI OLMAYARAQ KƏSİLİR. `None` = hədd yox.
    * `allow_knn_fallback` — radius boş qalanda k-ən-yaxına keçilsinmi.
    * `allow_global_fallback` — hər şey uğursuz olanda BÜTÜN nöqtələr
      işlədilsinmi (defolt `False` — "qlobal ehtiyat" yalnız AÇIQ
      icazə ilə, çünki bu, uzaq/əlaqəsiz məlumatla proqnoz deməkdir).
    * `support_range` — dəstək təsnifatında "korrelyasiya radiusu" kimi
      işlədilən məsafə (adətən `range_major`). `None` olanda təsnifat
      məlumatın öz orta qonşu məsafəsindən çıxarılır (bax `_support_scale`).
    """

    min_neighbors: int = 1
    max_neighbors: Optional[int] = None
    search_radius: Optional[float] = None
    max_search_radius: Optional[float] = None
    radius_expansion_factor: float = 2.0
    max_radius_expansions: int = 0
    sectors: int = 0
    vertical_sectors: bool = False
    max_per_sector: Optional[int] = None
    candidate_pool_factor: int = 0
    max_vertical_distance: Optional[float] = None
    allow_knn_fallback: bool = False
    allow_global_fallback: bool = False
    support_range: Optional[float] = None

    def validate(self) -> None:
        if self.min_neighbors < 1:
            raise NeighborhoodError(
                f"min_neighbors ≥ 1 olmalıdır, alındı: {self.min_neighbors}")
        if self.max_neighbors is not None:
            if self.max_neighbors < 1:
                raise NeighborhoodError(
                    f"max_neighbors ≥ 1 olmalıdır, alındı: {self.max_neighbors}")
            if self.max_neighbors < self.min_neighbors:
                raise NeighborhoodError(
                    f"max_neighbors ({self.max_neighbors}) < min_neighbors "
                    f"({self.min_neighbors}) — qonşuluq heç vaxt qurula bilməz.")
        for name in ("search_radius", "max_search_radius", "max_vertical_distance",
                     "support_range"):
            value = getattr(self, name)
            if value is not None and (not np.isfinite(value) or value <= 0.0):
                raise NeighborhoodError(f"{name} müsbət və sonlu olmalıdır, alındı: {value!r}")
        if (self.search_radius is not None and self.max_search_radius is not None
                and self.max_search_radius < self.search_radius):
            raise NeighborhoodError(
                f"max_search_radius ({self.max_search_radius}) < search_radius "
                f"({self.search_radius}).")
        if self.max_radius_expansions < 0:
            raise NeighborhoodError(
                f"max_radius_expansions ≥ 0 olmalıdır, alındı: {self.max_radius_expansions}")
        if self.max_radius_expansions and self.radius_expansion_factor <= 1.0:
            raise NeighborhoodError(
                "radius_expansion_factor > 1 olmalıdır (əks halda genişləndirmə "
                f"heç vaxt bitməz), alındı: {self.radius_expansion_factor}")
        if self.sectors < 0 or self.sectors == 1:
            raise NeighborhoodError(
                f"sectors 0 (söndürülmüş) və ya ≥ 2 olmalıdır, alındı: {self.sectors}")
        if self.max_per_sector is not None and self.max_per_sector < 1:
            raise NeighborhoodError(
                f"max_per_sector ≥ 1 olmalıdır, alındı: {self.max_per_sector}")
        if self.candidate_pool_factor < 0:
            raise NeighborhoodError(
                f"candidate_pool_factor ≥ 0 olmalıdır, alındı: {self.candidate_pool_factor}")

    @property
    def n_sector_cells(self) -> int:
        """Sektor "xanalarının" ümumi sayı (şaquli bölmə daxil)."""
        if self.sectors <= 0:
            return 0
        return self.sectors * (2 if self.vertical_sectors else 1)

    @property
    def knn_pool_size(self) -> Optional[int]:
        """k-ən-yaxın rejimdə çəkiləcək NAMİZƏD sayı (`None` = hamısı)."""
        if self.max_neighbors is None:
            return None
        if not self.sectors:
            return self.max_neighbors
        factor = self.candidate_pool_factor or self.n_sector_cells
        return self.max_neighbors * max(int(factor), 1)


@dataclass
class NeighborhoodResult:
    """Bir hədəf üçün seçilmiş qonşuluq + NECƏ seçildiyi (A2.7/A2.8)."""

    indices: np.ndarray                 #: qlobal nöqtə indeksləri, məsafəyə görə artan
    distances: np.ndarray               #: anizotrop məsafələr (transformasiya fəzasında)
    status: str = STATUS_EMPTY          #: `STATUS_*`
    support: str = SUPPORT_EXTRAPOLATED  #: `SUPPORT_*`
    radius_used: Optional[float] = None
    n_sectors_total: int = 0            #: təsnifat üçün işlədilən kvadrant/oktant sayı
    n_sectors_occupied: int = 0
    n_candidates: int = 0               #: kəsilmədən ƏVVƏL tapılan namizəd sayı
    warnings: List[str] = field(default_factory=list)

    @property
    def count(self) -> int:
        return int(self.indices.size)

    @property
    def nearest_distance(self) -> float:
        """Ən yaxın seçilmiş qonşunun məsafəsi; qonşu yoxdursa `inf`."""
        return float(self.distances[0]) if self.distances.size else float("inf")

    @property
    def is_extrapolation(self) -> bool:
        return self.support == SUPPORT_EXTRAPOLATED

    @property
    def ok(self) -> bool:
        """Kriging sistemi qurula bilərmi (ən azı bir qonşu var)."""
        return self.indices.size > 0


def _empty_result(status: str = STATUS_EMPTY,
                  warnings: Optional[List[str]] = None) -> NeighborhoodResult:
    return NeighborhoodResult(np.array([], dtype=int), np.array([], dtype=float),
                              status=status, support=SUPPORT_EXTRAPOLATED,
                              warnings=list(warnings or []))


def _order_by_distance(indices: np.ndarray, distances: np.ndarray
                       ) -> Tuple[np.ndarray, np.ndarray]:
    """DETERMİNİSTİK sıralama: əvvəlcə məsafə, bərabərlikdə indeks.

    `cKDTree.query_ball_point` sırasız qaytarır; bərabər məsafəli
    nöqtələrdə "hansı qonşu seçilir" sualı təkrar icralarda EYNİ cavabı
    verməlidir (A2.3 determinizm tələbi)."""
    order = np.lexsort((indices, distances))
    return indices[order], distances[order]


def _sector_ids(offsets: np.ndarray, sectors: int, vertical: bool) -> np.ndarray:
    """Hədəfə nəzərən yerdəyişmə vektorlarından sektor nömrələri.

    Bucaq `atan2(Δy, Δx)` ilə [0, 2π) aralığına gətirilir və `sectors`
    bərabər hissəyə bölünür; `vertical=True` olanda `Δz ≥ 0` / `Δz < 0`
    ikiyə bölməsi əlavə olunur (nəticə: `sectors·2` xana)."""
    angle = np.arctan2(offsets[:, 1], offsets[:, 0]) % (2.0 * np.pi)
    ids = np.minimum((angle / (2.0 * np.pi / sectors)).astype(int), sectors - 1)
    if vertical:
        ids = ids * 2 + (offsets[:, 2] < 0.0).astype(int)
    return ids


def _balance_by_sector(indices: np.ndarray, offsets: np.ndarray,
                       config: NeighborhoodConfig, limit: int) -> np.ndarray:
    """Sektorlar arasında NÖVBƏLİ (round-robin) seçim (A2.5).

    Namizədlər artıq məsafəyə görə sıralıdır. Hər dövrədə HƏR sektordan
    ən yaxın hələ götürülməmiş nöqtə alınır — beləliklə 16 qonşu bir
    kiçik klasterdən DEYİL, mövcud sektorlara paylanmış şəkildə seçilir.
    BOŞ sektor üçün heç nə uydurulmur; sektorlar tükənəndə qalan yerlər
    yenə ən yaxın nöqtələrlə doldurulur.

    Qaytarır: `indices` massivindəki MÖVQE indekslərini (məsafə sırası
    qorunmuş şəkildə)."""
    sector_ids = _sector_ids(offsets, config.sectors, config.vertical_sectors)
    n_cells = config.n_sector_cells
    hard_cap = config.max_per_sector is not None
    per_sector = config.max_per_sector
    if per_sector is None:
        per_sector = int(np.ceil(limit / n_cells)) if n_cells else limit
    per_sector = max(int(per_sector), 1)

    buckets: dict = {}
    for position, sector in enumerate(sector_ids):
        buckets.setdefault(int(sector), []).append(position)

    picked: List[int] = []
    cursors = {sector: 0 for sector in buckets}
    taken = {sector: 0 for sector in buckets}
    while len(picked) < limit:
        progressed = False
        for sector in sorted(buckets):          # deterministik sıra
            if len(picked) >= limit:
                break
            if taken[sector] >= per_sector:
                continue
            bucket = buckets[sector]
            cursor = cursors[sector]
            if cursor >= len(bucket):
                continue
            picked.append(bucket[cursor])
            cursors[sector] = cursor + 1
            taken[sector] += 1
            progressed = True
        if not progressed:
            break

    if len(picked) < limit and not hard_cap:
        # Sektorlar tükəndi (boş sektor üçün heç nə UYDURULMUR) → qalan
        # yerlər ən yaxın nöqtələrlə doldurulur. Bu YALNIZ kvota AVTOMATİK
        # hesablananda edilir: `max_per_sector` AÇIQ veriləndə o, QƏTİ
        # məhdudiyyətdir — qonşuluq `limit`-dən kiçik qalsa da pozulmur.
        chosen = set(picked)
        for position in range(indices.size):
            if len(picked) >= limit:
                break
            if position not in chosen:
                picked.append(position)
    return np.sort(np.asarray(picked, dtype=int))


def classify_support(distances: np.ndarray, offsets: np.ndarray,
                     support_scale: float, is_3d: bool) -> Tuple[str, int, int]:
    """Hədəfin məlumat buludu ilə HƏQİQİ HƏNDƏSİ münasibətini təsnif edir.

    Bu, "qonşu sayı < N → ekstrapolyasiya" evristikası DEYİL (A2.8) —
    iki müstəqil ölçüdən istifadə edilir:

    1. `h = ən_yaxın_məsafə / support_scale` — hədəf ümumiyyətlə
       korrelyasiya radiusu daxilindədirmi.
    2. SEKTOR ƏHATƏSİ — qonşular hədəfin ƏTRAFINA paylanıb, yoxsa
       hamısı BİR TƏRƏFDƏDİR. 2D-də 4 kvadrant, 3D-də 8 oktant
       (təsnifat konfiqurasiyadakı `sectors`-dan ASILI DEYİL ki, nəticə
       balanslaşdırma parametrinə görə dəyişməsin).

    Qərar (deterministik, sıra ilə):
        qonşu yoxdur ....................... EXTRAPOLATED
        h > 1 (radiusdan kənar) ............ EXTRAPOLATED
        əhatə ≤ 1/2 (məlumat bir tərəfdə) .. BOUNDARY
        h > 1/2 və ya qonşu < 4 ............ WEAK
        qalan .............................. WELL_SUPPORTED

    Qaytarır: `(support, n_sectors_total, n_sectors_occupied)`.
    """
    n_cells = 8 if is_3d else 4
    if distances.size == 0:
        return SUPPORT_EXTRAPOLATED, n_cells, 0

    if is_3d:
        occupied_ids = ((offsets[:, 0] >= 0.0).astype(int) * 4
                        + (offsets[:, 1] >= 0.0).astype(int) * 2
                        + (offsets[:, 2] >= 0.0).astype(int))
    else:
        occupied_ids = ((offsets[:, 0] >= 0.0).astype(int) * 2
                        + (offsets[:, 1] >= 0.0).astype(int))
    occupied = int(np.unique(occupied_ids).size)

    scale = support_scale if (np.isfinite(support_scale) and support_scale > 0.0) else 1.0
    h = float(distances[0]) / scale
    coverage = occupied / n_cells

    if h > 1.0:
        support = SUPPORT_EXTRAPOLATED
    elif coverage <= 0.5:
        support = SUPPORT_BOUNDARY
    elif h > 0.5 or distances.size < 4:
        support = SUPPORT_WEAK
    else:
        support = SUPPORT_WELL
    return support, n_cells, occupied


@dataclass
class BatchNeighborhood:
    """ÇOXLU hədəf üçün qonşuluq nəticələri — massiv formasında (A7).

    Hər hədəf üçün ayrıca `NeighborhoodResult` obyekti yaratmaq 10⁴-10⁶
    hüceyrəli şəbəkədə Python obyekt xərcinə görə hesablamanın ÖZÜNDƏN
    baha başa gəlir. Bu quruluş eyni məlumatı sıx massivlərdə saxlayır:

        indices    (m, kmax) — boş yerlər -1
        distances  (m, kmax) — boş yerlər +inf
        counts     (m,)      — faktiki qonşu sayı

    `row(i)` lazım olanda tək sətri `NeighborhoodResult`-a çevirir, ona
    görə introspeksiya/hesabat yolu dəyişmir.
    """

    indices: np.ndarray
    distances: np.ndarray
    counts: np.ndarray
    support: np.ndarray
    status: np.ndarray
    warnings: List[str] = field(default_factory=list)

    def row(self, index: int) -> NeighborhoodResult:
        count = int(self.counts[index])
        return NeighborhoodResult(
            indices=self.indices[index, :count].copy(),
            distances=self.distances[index, :count].copy(),
            status=str(self.status[index]), support=str(self.support[index]),
            n_candidates=count)

    def neighbours(self, index: int) -> np.ndarray:
        return self.indices[index, :int(self.counts[index])]

    @property
    def n_targets(self) -> int:
        return int(self.counts.size)


def classify_support_batch(distances: np.ndarray, offsets: np.ndarray,
                           counts: np.ndarray, support_scale: float,
                           is_3d: bool) -> np.ndarray:
    """`classify_support`-un vektorlaşdırılmış (çoxlu hədəf) variantı.

    Qaydalar BİREBİR eynidir (bax `classify_support`) — burada yalnız
    hər hədəf üçün ayrıca Python çağırışı əvəzinə massiv əməliyyatları
    işlədilir. `offsets` (m, kmax, 3); `counts`-dan artıq sütunlar
    nəzərə alınmır.
    """
    m, kmax = distances.shape
    n_cells = 8 if is_3d else 4
    valid = np.arange(kmax)[None, :] < counts[:, None]

    if is_3d:
        ids = ((offsets[:, :, 0] >= 0.0).astype(np.int64) * 4
               + (offsets[:, :, 1] >= 0.0).astype(np.int64) * 2
               + (offsets[:, :, 2] >= 0.0).astype(np.int64))
    else:
        ids = ((offsets[:, :, 0] >= 0.0).astype(np.int64) * 2
               + (offsets[:, :, 1] >= 0.0).astype(np.int64))

    occupied_mask = np.zeros((m, n_cells), dtype=bool)
    rows = np.repeat(np.arange(m), kmax).reshape(m, kmax)
    occupied_mask[rows[valid], ids[valid]] = True
    coverage = occupied_mask.sum(axis=1) / n_cells

    scale = support_scale if (np.isfinite(support_scale) and support_scale > 0.0) else 1.0
    nearest = np.where(counts > 0, distances[:, 0], np.inf)
    h = nearest / scale

    support = np.full(m, SUPPORT_WELL, dtype=object)
    support[(h > 0.5) | (counts < 4)] = SUPPORT_WEAK
    support[coverage <= 0.5] = SUPPORT_BOUNDARY
    support[(h > 1.0) | (counts == 0)] = SUPPORT_EXTRAPOLATED
    return support


class NeighborhoodSelector:
    """Sabit nöqtə çoxluğu üzərində PEŞƏKAR qonşuluq seçicisi (A2).

    Nöqtələr bir dəfə `AnisotropyParams.transform()`-dan keçirilir və
    (defolt) `cKDTree`-yə verilir — hər hədəf sorğusu `O(log n)`-dir,
    HEÇ BİR hədəf üçün tam `O(n)` skan aparılmır (A2.1/A7).

    `anisotropy=None` — izotrop Evklid fəzası (nöqtələr (n,3)-ə padding
    edilir, məsafə dəyişmir). Kriging bu formada işlədir: nöqtələr
    ARTIQ transformasiya edilmiş halda verilir, ikiqat transformasiya
    riyazi səhv olardı.

    `index` — indeks strategiyası:
        ``"kdtree"`` (defolt) ağac qurur — ÇOX hədəf üçün;
        ``"brute"``  ağac QURMUR, hər sorğu `O(n)` skandır — AZ hədəf
                     üçün (məs. SGS/SIS hər hüceyrədə TƏK hədəflə
                     çağırır: orada `O(n log n)` ağac qurmaq bahadır).
    Strategiya nəticəyə TƏSİR ETMİR, yalnız xərcə (bax
    `tests/test_neighborhood_selection.py` — iki rejimin PARİTETİ).
    """

    def __init__(self, points: np.ndarray, anisotropy: Optional[AnisotropyParams] = None,
                 config: Optional[NeighborhoodConfig] = None, index: str = "kdtree"):
        raw = np.atleast_2d(np.asarray(points, float))
        if raw.size == 0:
            raw = raw.reshape(0, 3)
        self.anisotropy = anisotropy
        self.config = config or NeighborhoodConfig()
        self.config.validate()
        self.raw_points = raw
        self.transformed = (transform_points(raw, anisotropy) if raw.shape[0]
                            else np.zeros((0, 3)))
        self.n_points = int(self.transformed.shape[0])
        self._raw_z = (raw[:, 2].copy() if raw.shape[1] > 2
                       else np.zeros(raw.shape[0]))
        self._is_3d = bool(self.n_points and np.ptp(self.transformed[:, 2]) > 1e-12)
        self._tree: Optional[cKDTree] = None
        self._index_mode = "kdtree"
        self._support_scale: Optional[float] = None
        self.set_index_mode(index)

    # ── konfiqurasiya ─────────────────────────────────────────────────
    def set_index_mode(self, mode: str) -> None:
        """`"kdtree"` / `"brute"` — bax sinif docstring-i."""
        if mode not in ("kdtree", "brute"):
            raise NeighborhoodError(
                f"index rejimi 'kdtree' və ya 'brute' olmalıdır, alındı: {mode!r}")
        self._index_mode = mode
        if mode == "kdtree" and self.n_points and self._tree is None:
            self._tree = cKDTree(self.transformed)
        elif mode == "brute":
            self._tree = None

    def set_raw_vertical(self, raw_z: np.ndarray) -> None:
        """`max_vertical_distance` üçün XAM Z sütununu AÇIQ təyin edir.

        Kriging nöqtələri artıq transformasiya edilmiş halda verir, ona
        görə seçici xam Z-ni özü bilə bilməz — A2.6 kəsiyinin XAM
        vahiddə qalması üçün bu metod çağırılır."""
        raw_z = np.asarray(raw_z, float).ravel()
        if raw_z.size != self.n_points:
            raise NeighborhoodError(
                f"raw_z uzunluğu ({raw_z.size}) nöqtə sayına ({self.n_points}) "
                "bərabər olmalıdır.")
        self._raw_z = raw_z

    # ── daxili köməkçilər ─────────────────────────────────────────────
    def _support_scale_value(self) -> float:
        """Dəstək təsnifatının "korrelyasiya radiusu" (tənbəl hesablanır).

        `config.support_range` verilibsə odur (adətən `range_major`) —
        Kriging HƏMİŞƏ bunu ötürür, ona görə istehsal yolunda təsnifat
        FAKTİKİ variogram radiusuna görə aparılır.

        Verilməyəndə (seçici təkbaşına işlədiləndə) məlumatın ÖZ
        həndəsəsindən çıxarılır: əhatə qutusunun diaqonalı / 3 — məhz
        `OrdinaryKriging`-in radius verilmədikdə işlətdiyi `domen/3`
        evristikası (bax `interpolation._parameters`), ona görə iki modul
        arasında UYĞUNSUZLUQ yaranmır. Uydurulmuş sabit YOXDUR."""
        if self.config.support_range is not None:
            return float(self.config.support_range)
        if self._support_scale is not None:
            return self._support_scale
        if self.n_points < 2:
            self._support_scale = float("inf")   # tək nöqtə: radius təyin edilmir
            return self._support_scale
        extent = self.transformed.max(axis=0) - self.transformed.min(axis=0)
        diagonal = float(np.sqrt(np.sum(extent ** 2)))
        self._support_scale = max(diagonal / 3.0, 1e-12)
        return self._support_scale

    def _raw_radius(self, target_t: np.ndarray, radius: float
                    ) -> Tuple[np.ndarray, np.ndarray]:
        if self._tree is not None:
            found = np.asarray(self._tree.query_ball_point(target_t, r=float(radius)),
                               dtype=int)
            if found.size == 0:
                return found, np.array([], dtype=float)
            distances = np.linalg.norm(self.transformed[found] - target_t, axis=1)
        else:
            distances_all = np.linalg.norm(self.transformed - target_t, axis=1)
            found = np.where(distances_all <= float(radius))[0]
            distances = distances_all[found]
        return _order_by_distance(found, distances)

    def _raw_knn(self, target_t: np.ndarray, k: int) -> Tuple[np.ndarray, np.ndarray]:
        k = int(min(max(k, 1), self.n_points))
        if self._tree is not None:
            distances, found = self._tree.query(target_t, k=k)
            found = np.atleast_1d(found).astype(int)
            distances = np.atleast_1d(distances).astype(float)
        else:
            distances_all = np.linalg.norm(self.transformed - target_t, axis=1)
            found = np.argpartition(distances_all, k - 1)[:k] if k < self.n_points \
                else np.arange(self.n_points)
            distances = distances_all[found]
        return _order_by_distance(found, distances)

    def _knn_after_vertical_cut(self, target_t: np.ndarray, want: int,
                                target_raw_z: float) -> Tuple[np.ndarray, np.ndarray]:
        """Şaquli kəsikdən SONRA `want` qonşu qalana qədər k-ni böyüdür.

        Kəsik SÜZGƏCDİR: onu `max_neighbors` kəsimindən SONRA tətbiq
        etmək qonşuluğu süni şəkildə kiçildərdi (40 ən yaxın nöqtənin
        yalnız 4-ü öz layından ola bilər). Ona görə namizəd hovuzu
        süzgəcdən keçən sayı təmin edənə qədər (ən çoxu bütün məlumat
        çoxluğuna qədər) genişləndirilir — nöqtə UYDURULMUR, sadəcə
        mövcud olanlar tam nəzərdən keçirilir."""
        k = max(int(want), 1)
        while True:
            indices, distances = self._raw_knn(target_t, k)
            indices, distances = self._apply_vertical_cut(indices, distances,
                                                          target_raw_z)
            if indices.size >= want or k >= self.n_points:
                return indices, distances
            k = min(self.n_points, max(k * 4, k + 1))

    def _apply_vertical_cut(self, indices: np.ndarray, distances: np.ndarray,
                            target_raw_z: float) -> Tuple[np.ndarray, np.ndarray]:
        """A2.6 — XAM |ΔZ| həddi (anizotrop miqyaslanmadan ASILI OLMAYAN
        geoloji kəsik: "başqa layın nöqtəsi ümumiyyətlə iştirak etməsin")."""
        limit = self.config.max_vertical_distance
        if limit is None or indices.size == 0:
            return indices, distances
        keep = np.abs(self._raw_z[indices] - target_raw_z) <= float(limit)
        return indices[keep], distances[keep]

    # ── əsas API ──────────────────────────────────────────────────────
    def select(self, target, raw_vertical: Optional[float] = None) -> NeighborhoodResult:
        """Bir hədəf üçün tam qonşuluq qərarı (A2.7 ehtiyat zənciri).

        Zəncir: anizotrop radius → (genişləndirmə) → k-ən-yaxın ehtiyat →
        (açıq icazə ilə) qlobal → sektor balanslaşdırması → `max_neighbors`
        kəsimi → `min_neighbors` yoxlaması → dəstək təsnifatı.

        `raw_vertical` — hədəfin XAM Z-si (`max_vertical_distance` üçün);
        verilməyəndə hədəfin üçüncü sütunu işlədilir.
        """
        config = self.config
        if self.n_points == 0:
            return _empty_result(STATUS_EMPTY, ["Məlumat çoxluğu boşdur — qonşu yoxdur."])

        raw_target = np.atleast_2d(np.asarray(target, float))
        if not np.all(np.isfinite(raw_target[0])):
            return _empty_result(
                STATUS_EMPTY,
                ["Hədəf koordinatı NaN/sonsuzdur — qonşuluq təyin edilə bilməz."])
        target_t = transform_points(raw_target, self.anisotropy)[0]
        if raw_vertical is not None:
            target_raw_z = float(raw_vertical)
        else:
            target_raw_z = float(raw_target[0, 2]) if raw_target.shape[1] > 2 else 0.0

        warnings: List[str] = []
        radius_used: Optional[float] = None
        pool = config.knn_pool_size
        knn_k = pool if pool is not None else self.n_points

        if config.search_radius is None:
            status = STATUS_KNN
            if config.max_vertical_distance is None:
                indices, distances = self._raw_knn(target_t, knn_k)
            else:
                indices, distances = self._knn_after_vertical_cut(
                    target_t, knn_k, target_raw_z)
        else:
            radius_used = float(config.search_radius)
            indices, distances = self._raw_radius(target_t, radius_used)
            status = STATUS_RADIUS
            indices, distances = self._apply_vertical_cut(indices, distances, target_raw_z)

            expansions = 0
            while (indices.size < config.min_neighbors
                   and expansions < config.max_radius_expansions):
                next_radius = radius_used * config.radius_expansion_factor
                if config.max_search_radius is not None:
                    next_radius = min(next_radius, float(config.max_search_radius))
                    if next_radius <= radius_used:
                        break
                radius_used = next_radius
                indices, distances = self._raw_radius(target_t, radius_used)
                indices, distances = self._apply_vertical_cut(indices, distances,
                                                              target_raw_z)
                expansions += 1
                status = STATUS_RADIUS_EXPANDED
            if expansions:
                warnings.append(
                    f"Axtarış radiusu {expansions} dəfə genişləndirildi "
                    f"({config.search_radius:.4g} → {radius_used:.4g}) — bu hədəfdə "
                    "məlumat seyrəkdir.")

            if indices.size < config.min_neighbors and config.allow_knn_fallback:
                indices, distances = self._knn_after_vertical_cut(
                    target_t, knn_k, target_raw_z)
                status = STATUS_KNN_FALLBACK
                warnings.append(
                    "Radius axtarışı kifayət qədər qonşu tapmadı — k-ən-yaxın "
                    "ehtiyat yolu işlədildi (qonşular axtarış radiusundan KƏNAR ola bilər).")

        if indices.size < config.min_neighbors and config.allow_global_fallback:
            indices, distances = self._raw_knn(target_t, self.n_points)
            status = STATUS_GLOBAL
            warnings.append(
                "QLOBAL ehtiyat yolu işlədildi — proqnoz bütün məlumat nöqtələri ilə "
                "qurulub, o cümlədən korrelyasiya radiusundan çox uzaqları.")

        n_candidates = int(indices.size)

        limit = config.max_neighbors
        if limit is not None and indices.size > limit:
            if config.sectors:
                offsets = self.transformed[indices] - target_t
                positions = _balance_by_sector(indices, offsets, config, limit)
                indices, distances = indices[positions], distances[positions]
            else:
                indices, distances = indices[:limit], distances[:limit]

        if indices.size < config.min_neighbors:
            return _empty_result(
                STATUS_INSUFFICIENT,
                warnings + [
                    f"min_neighbors={config.min_neighbors} təmin edilmədi "
                    f"(tapıldı: {indices.size}) — dəyər UYDURULMUR."])

        offsets = self.transformed[indices] - target_t
        support, n_total, n_occupied = classify_support(
            distances, offsets, self._support_scale_value(), self._is_3d)
        return NeighborhoodResult(indices=indices, distances=distances, status=status,
                                  support=support, radius_used=radius_used,
                                  n_sectors_total=n_total, n_sectors_occupied=n_occupied,
                                  n_candidates=n_candidates, warnings=warnings)

    def _is_simple_knn(self) -> bool:
        """Konfiqurasiya SADƏ k-ən-yaxın rejimidirmi (vektorlaşdırıla bilər).

        Radius, sektor balanslaşdırması, şaquli kəsik və ehtiyat yolları
        hədəfdən-hədəfə FƏRQLİ qərar tələb edir, ona görə onlar üçün
        sətir-sətir `select()` işlədilir. Bu hal isə (defolt istehsal
        yolu) bütün hədəflər üçün eyni `k` deməkdir — TƏK toplu cKDTree
        sorğusu ilə həll olunur."""
        config = self.config
        return (config.search_radius is None and not config.sectors
                and config.max_vertical_distance is None
                and not config.allow_global_fallback
                and config.max_radius_expansions == 0)

    def select_batch(self, targets: np.ndarray) -> BatchNeighborhood:
        """Bütün hədəflər üçün qonşuluq — massiv formasında (A7).

        Sadə k-ən-yaxın konfiqurasiyasında TƏK toplu `cKDTree.query()`
        çağırılır (C səviyyəsində, hədəf başına Python xərci YOX);
        digər hallarda `select()` sətir-sətir işləyib nəticə eyni sıx
        massivlərə yığılır. NƏTİCƏ hər iki yolda EYNİDİR.
        """
        targets = np.atleast_2d(np.asarray(targets, float))
        m = targets.shape[0]
        finite = np.all(np.isfinite(targets), axis=1)
        if m and not np.all(finite):
            # cKDTree qeyri-sonlu koordinat qəbul etmir; belə hədəflər
            # üçün qonşuluq TƏYİN EDİLƏ BİLMƏZ — sətirlər BOŞ qalır
            # (dəyər uydurulmur, bax `interpolation._prepare`).
            partial = self.select_batch(targets[finite]) if np.any(finite) else None
            kmax = partial.indices.shape[1] if partial is not None else 1
            indices = np.full((m, kmax), -1, dtype=int)
            distances = np.full((m, kmax), np.inf)
            counts = np.zeros(m, dtype=int)
            support = np.full(m, SUPPORT_EXTRAPOLATED, dtype=object)
            status = np.full(m, STATUS_EMPTY, dtype=object)
            warnings: List[str] = []
            if partial is not None:
                indices[finite] = partial.indices
                distances[finite] = partial.distances
                counts[finite] = partial.counts
                support[finite] = partial.support
                status[finite] = partial.status
                warnings = list(partial.warnings)
            return BatchNeighborhood(indices, distances, counts, support, status,
                                     warnings)
        if self.n_points == 0 or m == 0:
            empty_i = np.full((m, 1), -1, dtype=int)
            empty_d = np.full((m, 1), np.inf)
            return BatchNeighborhood(
                empty_i, empty_d, np.zeros(m, dtype=int),
                np.full(m, SUPPORT_EXTRAPOLATED, dtype=object),
                np.full(m, STATUS_EMPTY, dtype=object),
                ["Məlumat çoxluğu boşdur — qonşu yoxdur."] if m else [])

        if not self._is_simple_knn():
            return self._select_batch_by_rows(targets)

        targets_t = transform_points(targets, self.anisotropy)
        pool = self.config.knn_pool_size
        k = int(min(pool if pool is not None else self.n_points, self.n_points))

        if self._tree is not None:
            distances, indices = self._tree.query(targets_t, k=k)
        else:
            diff = self.transformed[None, :, :] - targets_t[:, None, :]
            full = np.sqrt(np.sum(diff * diff, axis=-1))
            indices = np.argsort(full, axis=1, kind="stable")[:, :k]
            distances = np.take_along_axis(full, indices, axis=1)
        indices = np.atleast_2d(np.asarray(indices, dtype=int).reshape(m, k))
        distances = np.atleast_2d(np.asarray(distances, float).reshape(m, k))

        # DETERMİNİZM: bərabər məsafələrdə indeksə görə sırala (bax
        # `_order_by_distance` — burada eyni qayda, sətir-sətir vektorlu).
        order = np.lexsort((indices, distances), axis=-1)
        indices = np.take_along_axis(indices, order, axis=-1)
        distances = np.take_along_axis(distances, order, axis=-1)

        counts = np.full(m, k, dtype=int)
        status = np.full(m, STATUS_KNN, dtype=object)
        warnings: List[str] = []
        if k < self.config.min_neighbors:
            counts[:] = 0
            status[:] = STATUS_INSUFFICIENT
            warnings.append(
                f"min_neighbors={self.config.min_neighbors} təmin edilmədi "
                f"(mövcud nöqtə: {self.n_points}) — dəyər UYDURULMUR.")
            indices = np.full((m, max(k, 1)), -1, dtype=int)
            distances = np.full((m, max(k, 1)), np.inf)
            support = np.full(m, SUPPORT_EXTRAPOLATED, dtype=object)
            return BatchNeighborhood(indices, distances, counts, support, status,
                                     warnings)

        offsets = self.transformed[indices] - targets_t[:, None, :]
        support = classify_support_batch(distances, offsets, counts,
                                         self._support_scale_value(), self._is_3d)
        return BatchNeighborhood(indices, distances, counts, support, status, warnings)

    def _select_batch_by_rows_with_depth(self, targets: np.ndarray,
                                         raw_vertical: np.ndarray) -> BatchNeighborhood:
        """`_select_batch_by_rows`-un XAM Z-li variantı.

        `max_vertical_distance` XAM dərinlik fərqi ilə ölçülür; Kriging
        hədəfləri transformasiya edilmiş halda ötürdüyü üçün hər hədəfin
        xam Z-si AYRICA verilir (bax `select(raw_vertical=...)`)."""
        raw_vertical = np.asarray(raw_vertical, float).ravel()
        results = [self.select(targets[row:row + 1], raw_vertical=float(raw_vertical[row]))
                   for row in range(targets.shape[0])]
        return self._pack(results)

    def _select_batch_by_rows(self, targets: np.ndarray) -> BatchNeighborhood:
        """Mürəkkəb konfiqurasiya — `select()` sətir-sətir, nəticə sıx
        massivlərə yığılır (qonşu sayı hədəfdən-hədəfə dəyişə bilər)."""
        return self._pack([self.select(targets[row:row + 1])
                           for row in range(targets.shape[0])])

    @staticmethod
    def _pack(results: List[NeighborhoodResult]) -> BatchNeighborhood:
        """`NeighborhoodResult` siyahısını sıx massivlərə yığır."""
        m = len(results)
        kmax = max((r.count for r in results), default=0) or 1
        indices = np.full((m, kmax), -1, dtype=int)
        distances = np.full((m, kmax), np.inf)
        counts = np.zeros(m, dtype=int)
        support = np.empty(m, dtype=object)
        status = np.empty(m, dtype=object)
        warnings: List[str] = []
        for row, result in enumerate(results):
            count = result.count
            counts[row] = count
            if count:
                indices[row, :count] = result.indices
                distances[row, :count] = result.distances
            support[row] = result.support
            status[row] = result.status
            for message in result.warnings:
                if message not in warnings:
                    warnings.append(message)
        return BatchNeighborhood(indices, distances, counts, support, status, warnings)

    def select_many(self, targets: np.ndarray) -> List[NeighborhoodResult]:
        """Çoxlu hədəf — hər biri üçün `select()`. Ağac TƏKRAR qurulmur."""
        targets = np.atleast_2d(np.asarray(targets, float))
        return [self.select(targets[row:row + 1]) for row in range(targets.shape[0])]



class AnisotropicNeighborSearch:
    """SABİT (dəyişməyən) nöqtə çoxluğu üzərində qonşu axtarışı.

    `_distance_matrix`-in TAM brute-force analoqudur — YALNIZ sürətlə
    (cKDTree, O(log n) sorğu) — nəticə EYNİ olmalıdır (bax
    `tests/test_spatial_search.py`).

    Bu, aşağı qatdır: sadə radius/k-ən-yaxın. Sektor balanslaşdırması,
    radius genişləndirməsi və dəstək təsnifatı üçün `NeighborhoodSelector`
    işlədilir (SGS/SIS bu sinfin ARDICIL variantına bağlıdır, ona görə
    müqaviləsi dəyişdirilmir).
    """

    def __init__(self, points: np.ndarray, anisotropy: Optional[AnisotropyParams] = None):
        points = np.atleast_2d(np.asarray(points, float))
        self.anisotropy = anisotropy
        self.transformed = self._transform(points)
        self.n_points = self.transformed.shape[0]
        self._tree = cKDTree(self.transformed) if self.n_points > 0 else None

    def _transform(self, points: np.ndarray) -> np.ndarray:
        points = np.atleast_2d(np.asarray(points, float))
        return self.anisotropy.transform(points) if self.anisotropy is not None else points

    def query(self, target, search_radius: Optional[float] = None,
             max_neighbors: Optional[int] = None, min_neighbors: int = 1,
             sectors: int = 0, vertical_sectors: bool = False,
             max_per_sector: Optional[int] = None) -> np.ndarray:
        """Bir hədəf nöqtə üçün seçilmiş qonşu İNDEKSLƏRİ (məsafəyə görə
        ARTAN sıralı). `min_neighbors`-dan az qonşu tapılarsa BOŞ massiv
        qaytarır — NaN/dəyər UYDURULMUR, qərar çağırana aiddir (`facies.py`
        bunu qlobal nisbətə keçid kimi işlədir, bax `FACIES.md`).

        `sectors > 0` (defolt 0 = söndürülmüş, ƏVVƏLKİ davranış) —
        `max_neighbors` kəsimi ən yaxın N əvəzinə SEKTORLAR arasında
        balanslaşdırılmış şəkildə aparılır (A2.5)."""
        if self._tree is None or self.n_points == 0:
            return np.array([], dtype=int)
        target_t = self._transform(target)[0]

        if search_radius is not None:
            candidate = np.asarray(
                self._tree.query_ball_point(target_t, r=float(search_radius)), dtype=int)
            if candidate.size == 0:
                return candidate
            distances = np.linalg.norm(self.transformed[candidate] - target_t, axis=1)
            order = np.argsort(distances, kind="stable")
            candidate, distances = candidate[order], distances[order]
        else:
            k = self.n_points if max_neighbors is None else min(max_neighbors, self.n_points)
            distances, candidate = self._tree.query(target_t, k=k)
            candidate = np.atleast_1d(candidate).astype(int)
            distances = np.atleast_1d(distances).astype(float)

        if max_neighbors is not None and candidate.size > max_neighbors:
            if sectors:
                config = NeighborhoodConfig(sectors=sectors,
                                            vertical_sectors=vertical_sectors,
                                            max_per_sector=max_per_sector)
                offsets = self.transformed[candidate] - target_t
                if offsets.shape[1] < 3:
                    offsets = np.column_stack([offsets, np.zeros(offsets.shape[0])])
                positions = _balance_by_sector(candidate, offsets, config,
                                               max_neighbors)
                candidate = candidate[positions]
            else:
                candidate = candidate[:max_neighbors]
        if candidate.size < max(min_neighbors, 1):
            return np.array([], dtype=int)
        return candidate


class IncrementalAnisotropicSearch:
    """`AnisotropicNeighborSearch`-in ARDICIL (SIS kimi bir-bir nöqtə
    ƏLAVƏ OLUNAN) simulyasiyalar üçün genişlənməsi.

    cKDTree DƏYİŞMƏZDİR (immutable) — hər addımda YENİDƏN qurmaq
    `O(n log n)` olardı, YƏNİ TAM brute-force-dan (`O(n)`) belə PISDİR.
    Bunun əvəzinə: `rebuild_interval` nöqtədən bir ağac YENİDƏN qurulur;
    aralıqda (son qurulmadan bəri) əlavə olunan nöqtələr KİÇİK bufer
    kimi brute-force axtarılır və ağac nəticəsi ilə BİRLƏŞDİRİLİR.

    NƏTİCƏ TAM brute-force ilə EYNİDİR (approksimasiya DEYİL) — heç bir
    nöqtə axtarışdan kənarda qalmır, YALNIZ axtarış YÜKÜ ağac (sürətli)
    və bufer (yavaş, amma KİÇİK) arasında bölünür. Mürəkkəblik: axtarış
    hissəsi təxminən `O(n log n)` (ümumi, `rebuild_interval` addımda bir
    `O(m log m)` qurma, m ≈ o anki nöqtə sayı) — brute-force-un `O(n²)`-
    dən aşağı, bax `FACIES.md` "Mürəkkəblik (Phase 4.1)".
    """

    def __init__(self, initial_points: np.ndarray, anisotropy: Optional[AnisotropyParams] = None,
                rebuild_interval: int = 64):
        self.anisotropy = anisotropy
        self.rebuild_interval = max(int(rebuild_interval), 1)
        initial = np.atleast_2d(np.asarray(initial_points, float))
        self._all_transformed = self._transform(initial)
        self._tree_size = 0
        self._tree: Optional[cKDTree] = None
        self._rebuild_tree()

    def _transform(self, points: np.ndarray) -> np.ndarray:
        points = np.atleast_2d(np.asarray(points, float))
        return self.anisotropy.transform(points) if self.anisotropy is not None else points

    def _rebuild_tree(self) -> None:
        self._tree_size = self._all_transformed.shape[0]
        self._tree = cKDTree(self._all_transformed) if self._tree_size > 0 else None

    def add_point(self, point) -> None:
        transformed = self._transform(point)
        self._all_transformed = np.vstack([self._all_transformed, transformed])
        if self._all_transformed.shape[0] - self._tree_size >= self.rebuild_interval:
            self._rebuild_tree()

    @property
    def n_points(self) -> int:
        return self._all_transformed.shape[0]

    def query(self, target, search_radius: Optional[float] = None,
             max_neighbors: Optional[int] = None, min_neighbors: int = 1,
             sectors: int = 0, vertical_sectors: bool = False,
             max_per_sector: Optional[int] = None) -> np.ndarray:
        """`AnisotropicNeighborSearch.query`-lə EYNİ müqavilə — YALNIZ
        indeksləri iki mənbədən (ağac + gözləmə buferi) BİRLƏŞDİRİB
        sıralayır."""
        target_t = self._transform(target)[0]
        pending = self._all_transformed[self._tree_size:]

        tree_idx = np.array([], dtype=int)
        tree_dist = np.array([])
        if self._tree is not None and self._tree_size > 0:
            if search_radius is not None:
                tree_idx = np.asarray(
                    self._tree.query_ball_point(target_t, r=float(search_radius)), dtype=int)
                if tree_idx.size:
                    tree_dist = np.linalg.norm(
                        self._all_transformed[tree_idx] - target_t, axis=1)
            else:
                k = self._tree_size if max_neighbors is None else min(max_neighbors,
                                                                       self._tree_size)
                tree_dist, tree_idx = self._tree.query(target_t, k=k)
                tree_idx = np.atleast_1d(tree_idx).astype(int)
                tree_dist = np.atleast_1d(tree_dist)

        pending_idx = np.array([], dtype=int)
        pending_dist = np.array([])
        if pending.shape[0]:
            distances = np.linalg.norm(pending - target_t, axis=1)
            local_idx = np.arange(pending.shape[0])
            if search_radius is not None:
                keep = distances <= search_radius
                local_idx, distances = local_idx[keep], distances[keep]
            pending_idx = local_idx + self._tree_size
            pending_dist = distances

        combined_idx = np.concatenate([tree_idx, pending_idx])
        combined_dist = np.concatenate([tree_dist, pending_dist])
        if combined_idx.size == 0:
            return combined_idx
        order = np.argsort(combined_dist, kind="stable")
        combined_idx = combined_idx[order]
        combined_dist = combined_dist[order]
        if max_neighbors is not None and combined_idx.size > max_neighbors:
            if sectors:
                config = NeighborhoodConfig(sectors=sectors,
                                            vertical_sectors=vertical_sectors,
                                            max_per_sector=max_per_sector)
                offsets = self._all_transformed[combined_idx] - target_t
                if offsets.shape[1] < 3:
                    offsets = np.column_stack([offsets, np.zeros(offsets.shape[0])])
                positions = _balance_by_sector(combined_idx, offsets, config,
                                               max_neighbors)
                combined_idx = combined_idx[positions]
            else:
                combined_idx = combined_idx[:max_neighbors]
        if combined_idx.size < max(min_neighbors, 1):
            return np.array([], dtype=int)
        return combined_idx


def build_neighborhood_config(search_radius: Optional[float] = None,
                              max_neighbors: Optional[int] = None,
                              min_neighbors: int = 1,
                              **kwargs) -> NeighborhoodConfig:
    """`OrdinaryKriging`-in köhnə üç parametrindən `NeighborhoodConfig`.

    Geriyə uyğunluq körpüsü (A12): mövcud çağırışlar hələ də yalnız
    `search_radius`/`max_neighbors`/`min_neighbors` verir; yeni imkanlar
    (`sectors`, radius genişləndirməsi və s.) `kwargs` ilə əlavə olunur."""
    return NeighborhoodConfig(min_neighbors=min_neighbors, max_neighbors=max_neighbors,
                              search_radius=search_radius, **kwargs)
