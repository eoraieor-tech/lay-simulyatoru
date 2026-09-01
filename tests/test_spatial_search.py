"""Phase 4.1 — anizotrop məkan axtarışı: cKDTree vs brute-force PARİTETİ.

`OrdinaryKriging._solve_local` (Phase 2-3) hər hədəf üçün TAM brute-force
məsafə hesablayıb sıralayır/kəsir — bu fayl `AnisotropicNeighborSearch`/
`IncrementalAnisotropicSearch`-in EYNİ qonşuluğu seçdiyini sübut edir,
BUNDAN SONRA `facies.py` bu sürətli yolu işlədə bilər (tapşırıq §6:
"Only then replace the repeated brute-force search").

`OrdinaryKriging`-in ÖZÜ bu fayldan İSTİFADƏ OLUNMUR (import edilmir) —
`_solve_local`-ın brute-force MƏNTİQİ burada AYRICA (test məqsədilə)
təkrarlanır, çünki o, `interpolate()`-in daxili `_variogram`/`_solve_
global` çağırışları ilə QARIŞIQDIR (yalnız qonşu SEÇİMİNİ təcrid etmək
mümkün deyil metodun içindən). Bu TƏKRARLAMA DEYİL — `OrdinaryKriging`-in
kriging RİYAZİYYATI (variogram, xətti sistem) heç yerdə təkrarlanmır,
YALNIZ "hansı nöqtələr qonşudur" sual-cavabı (sadə həndəsə) təkrarlanır.
"""

from __future__ import annotations

import numpy as np

from imex2d.geology.spatial_search import (AnisotropicNeighborSearch,
                                            IncrementalAnisotropicSearch)
from imex2d.geology.variogram import AnisotropyParams


def _brute_force_candidates(points_t: np.ndarray, target_t: np.ndarray,
                            search_radius=None, max_neighbors=None, min_neighbors=1
                            ) -> np.ndarray:
    """`interpolation.py`-dəki `OrdinaryKriging._solve_local`-ın qonşu-
    seçim MƏNTİQİNİN eynisi (artıq transformasiya edilmiş fəzada)."""
    distances = np.linalg.norm(points_t - target_t, axis=1)
    candidate = np.arange(points_t.shape[0])
    if search_radius is not None:
        candidate = candidate[distances[candidate] <= search_radius]
    if candidate.size == 0:
        return candidate
    candidate = candidate[np.argsort(distances[candidate], kind="stable")]
    if max_neighbors is not None:
        candidate = candidate[:max_neighbors]
    if candidate.size < max(min_neighbors, 1):
        return np.array([], dtype=int)
    return candidate


def _random_points(n, seed, low=0.0, high=500.0, ndim=2):
    rng = np.random.default_rng(seed)
    return rng.uniform(low, high, size=(n, ndim))


# ── izotrop: sadə Evklid fəzası ────────────────────────────────────────
def test_isotropic_max_neighbors_matches_brute_force():
    points = _random_points(80, seed=1)
    target = np.array([[317.0, 183.0]])
    search = AnisotropicNeighborSearch(points)

    expected = _brute_force_candidates(points, target[0], max_neighbors=10)
    actual = search.query(target, max_neighbors=10)
    assert np.array_equal(actual, expected)


def test_isotropic_search_radius_matches_brute_force():
    points = _random_points(80, seed=2)
    target = np.array([[317.0, 183.0]])
    search = AnisotropicNeighborSearch(points)

    expected = _brute_force_candidates(points, target[0], search_radius=120.0)
    actual = search.query(target, search_radius=120.0)
    assert np.array_equal(actual, expected)


def test_isotropic_radius_and_max_neighbors_together_matches_brute_force():
    points = _random_points(80, seed=3)
    target = np.array([[317.0, 183.0]])
    search = AnisotropicNeighborSearch(points)

    expected = _brute_force_candidates(points, target[0], search_radius=150.0, max_neighbors=6)
    actual = search.query(target, search_radius=150.0, max_neighbors=6)
    assert np.array_equal(actual, expected)


def test_min_neighbors_rejects_when_too_few_found():
    points = _random_points(5, seed=4, high=50.0)
    target = np.array([[1000.0, 1000.0]])
    search = AnisotropicNeighborSearch(points)
    result = search.query(target, search_radius=10.0, min_neighbors=1)
    assert result.size == 0


# ── anizotrop: transformasiya edilmiş fəzada axtarış ────────────────────
def test_anisotropic_search_matches_brute_force_in_transformed_space():
    """Kriging-in ÖZÜ hansı fəzada axtarırsa (transformasiya edilmiş),
    bu axtarış da EYNİ fəzada işləməlidir — bax tapşırıq §5 "IMPORTANT"."""
    points = _random_points(100, seed=5)
    target = np.array([[317.0, 183.0]])
    aniso = AnisotropyParams(azimuth_deg=35.0, range_major=200.0, range_minor=40.0,
                             range_vertical=200.0)

    points3 = np.column_stack([points, np.zeros(points.shape[0])])
    target3 = np.column_stack([target, np.zeros(1)])
    points_t = aniso.transform(points3)
    target_t = aniso.transform(target3)[0]

    expected = _brute_force_candidates(points_t, target_t, max_neighbors=8)
    search = AnisotropicNeighborSearch(points3, anisotropy=aniso)
    actual = search.query(target3, max_neighbors=8)
    assert np.array_equal(actual, expected)


def test_anisotropic_search_with_radius_matches_brute_force():
    points = _random_points(100, seed=6)
    target = np.array([[317.0, 183.0]])
    aniso = AnisotropyParams(azimuth_deg=110.0, range_major=300.0, range_minor=25.0,
                             range_vertical=300.0)
    points3 = np.column_stack([points, np.zeros(points.shape[0])])
    target3 = np.column_stack([target, np.zeros(1)])
    points_t = aniso.transform(points3)
    target_t = aniso.transform(target3)[0]

    expected = _brute_force_candidates(points_t, target_t, search_radius=90.0, max_neighbors=15)
    search = AnisotropicNeighborSearch(points3, anisotropy=aniso)
    actual = search.query(target3, search_radius=90.0, max_neighbors=15)
    assert np.array_equal(actual, expected)


def test_isotropic_default_no_limits_returns_all_points_sorted():
    points = _random_points(30, seed=7)
    target = np.array([[317.0, 183.0]])
    search = AnisotropicNeighborSearch(points)
    expected = _brute_force_candidates(points, target[0])
    actual = search.query(target)
    assert np.array_equal(actual, expected)
    assert actual.size == 30


# ── inkremental (ARDICIL əlavə olunan) axtarış — SIS-in özü ─────────────
def test_incremental_search_matches_brute_force_as_points_are_added():
    """SIS-in simulyasiya etdiyi ssenari: nöqtələr BİR-BİR əlavə olunur,
    hər addımda sorğu aparılır — nəticə HƏR ADDIMDA tam brute-force ilə
    EYNİ olmalıdır (approksimasiya DEYİL)."""
    rng = np.random.default_rng(80)   # fərqli seed — `initial` ilə eyni axından QAÇINMAQ üçün
    initial = _random_points(10, seed=8)
    incremental = IncrementalAnisotropicSearch(initial, rebuild_interval=7)
    all_points = initial.copy()

    target = np.array([[317.0, 183.0]])
    for _ in range(40):
        new_point = rng.uniform(0.0, 500.0, size=(1, 2))
        incremental.add_point(new_point)
        all_points = np.vstack([all_points, new_point])

        expected = _brute_force_candidates(all_points, target[0], max_neighbors=9)
        actual = incremental.query(target, max_neighbors=9)
        assert np.array_equal(actual, expected), f"n_points={all_points.shape[0]}"


def test_incremental_search_matches_brute_force_with_radius_and_anisotropy():
    rng = np.random.default_rng(90)   # fərqli seed — `initial` ilə eyni axından QAÇINMAQ üçün
    aniso = AnisotropyParams(azimuth_deg=60.0, range_major=250.0, range_minor=50.0,
                             range_vertical=250.0)
    initial = _random_points(8, seed=9)
    initial3 = np.column_stack([initial, np.zeros(8)])
    incremental = IncrementalAnisotropicSearch(initial3, anisotropy=aniso, rebuild_interval=5)
    all_points3 = initial3.copy()

    target3 = np.array([[317.0, 183.0, 0.0]])
    target_t = aniso.transform(target3)[0]
    for _ in range(30):
        new_point = np.column_stack([rng.uniform(0.0, 500.0, size=(1, 2)), np.zeros((1, 1))])
        incremental.add_point(new_point)
        all_points3 = np.vstack([all_points3, new_point])

        points_t = aniso.transform(all_points3)
        expected = _brute_force_candidates(points_t, target_t, search_radius=140.0,
                                           max_neighbors=6)
        actual = incremental.query(target3, search_radius=140.0, max_neighbors=6)
        assert np.array_equal(actual, expected)


def test_incremental_search_n_points_tracks_additions():
    initial = _random_points(5, seed=10)
    incremental = IncrementalAnisotropicSearch(initial, rebuild_interval=3)
    assert incremental.n_points == 5
    incremental.add_point(np.array([[1.0, 1.0]]))
    assert incremental.n_points == 6
