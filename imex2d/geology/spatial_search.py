"""Anizotrop-uyğun məkan axtarışı — `scipy.spatial.cKDTree` əsaslı (Phase 4.1).

Niyə lazımdır: `geology/interpolation.OrdinaryKriging._solve_local`
(Phase 2-3) hər çağırışda TAM brute-force məsafə matrisi qurur — SIS-də
(Phase 4) hər hüceyrə üçün BÖYÜYƏN kondisioner çoxluğu üzərində bu
TƏKRAR-TƏKRAR çağırılır, nəticədə `O(K·n²)` mürəkkəblik yaranır (bax
`FACIES.md`). Bu modul YENİDƏN İSTİFADƏ OLUNA BİLƏN, sınanmış bir
alternativ təqdim edir — `OrdinaryKriging`-in ÖZÜ İSƏ TOXUNULMAYIB
(76+ Phase 2-3 testi qorunur, bax `tests/test_kriging_3d_anisotropy.py`
və s.) — bu axtarış YALNIZ `facies.py`-də, ƏLAVƏ (opt-in) sürətləndirmə
kimi işlədilir.

KRİTİK QAYDA (tapşırıq §5-6): anizotrop axtarış Kriging-in ÖZÜNÜN
işlətdiyi EYNİ transformasiya fəzasında (bax
`geology/variogram.AnisotropyParams.transform`) aparılmalıdır — əks
halda seçilən qonşuluq Kriging sisteminin hesabladığı ilə UYĞUNSUZ olar.
Ona görə bu sinif XAM koordinatları DEYİL, `AnisotropyParams.transform()`-
dan keçmiş koordinatları indeksləyir (və ya heç bir anizotropluq
verilməyibsə, adi izotrop Evklid fəzasını).

`tests/test_spatial_search.py` bunun `_solve_local`-ın brute-force
seçimi ilə EYNİ nəticəni verdiyini sübut edir — YALNIZ bundan sonra
`facies.py`-də əvəzləmə edilib (tapşırıq §6: "Only then replace").
"""

from __future__ import annotations

from typing import Optional

import numpy as np
from scipy.spatial import cKDTree

from .variogram import AnisotropyParams


class AnisotropicNeighborSearch:
    """SABİT (dəyişməyən) nöqtə çoxluğu üzərində qonşu axtarışı.

    `_distance_matrix`-in TAM brute-force analoqudur — YALNIZ sürətlə
    (cKDTree, O(log n) sorğu) — nəticə EYNİ olmalıdır (bax
    `tests/test_spatial_search.py`).
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
             max_neighbors: Optional[int] = None, min_neighbors: int = 1) -> np.ndarray:
        """Bir hədəf nöqtə üçün seçilmiş qonşu İNDEKSLƏRİ (məsafəyə görə
        ARTAN sıralı). `min_neighbors`-dan az qonşu tapılarsa BOŞ massiv
        qaytarır — NaN/dəyər UYDURULMUR, qərar çağırana aiddir (`facies.py`
        bunu qlobal nisbətə keçid kimi işlədir, bax `FACIES.md`)."""
        if self._tree is None or self.n_points == 0:
            return np.array([], dtype=int)
        target_t = self._transform(target)[0]

        if search_radius is not None:
            candidate = np.asarray(
                self._tree.query_ball_point(target_t, r=float(search_radius)), dtype=int)
            if candidate.size == 0:
                return candidate
            distances = np.linalg.norm(self.transformed[candidate] - target_t, axis=1)
            candidate = candidate[np.argsort(distances, kind="stable")]
        else:
            k = self.n_points if max_neighbors is None else min(max_neighbors, self.n_points)
            _, candidate = self._tree.query(target_t, k=k)
            candidate = np.atleast_1d(candidate).astype(int)

        if max_neighbors is not None:
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
             max_neighbors: Optional[int] = None, min_neighbors: int = 1) -> np.ndarray:
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
        if max_neighbors is not None:
            combined_idx = combined_idx[:max_neighbors]
        if combined_idx.size < max(min_neighbors, 1):
            return np.array([], dtype=int)
        return combined_idx
