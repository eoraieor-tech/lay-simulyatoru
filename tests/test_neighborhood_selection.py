"""A2 — peşəkar qonşuluq seçimi: k-ən-yaxın, radius, anizotropluq,
istiqamətli balanslaşdırma, şaquli ayırma, seyrək-data ehtiyatı və
ekstrapolyasiya təsnifatı.

`test_spatial_search.py` (Phase 4.1) aşağı qatın (cKDTree vs brute-force)
PARİTETİNİ sübut edir — bu fayl isə YUXARI qatın (`NeighborhoodSelector`)
QƏRAR məntiqini yoxlayır.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.geology.anisotropy import AnisotropyParams
from imex2d.geology.spatial_search import (STATUS_EMPTY, STATUS_GLOBAL,
                                           STATUS_INSUFFICIENT, STATUS_KNN,
                                           STATUS_KNN_FALLBACK, STATUS_RADIUS,
                                           STATUS_RADIUS_EXPANDED, SUPPORT_BOUNDARY,
                                           SUPPORT_EXTRAPOLATED, SUPPORT_WELL,
                                           NeighborhoodConfig, NeighborhoodError,
                                           NeighborhoodSelector)


def _random_points(n, seed, high=500.0, ndim=2):
    rng = np.random.default_rng(seed)
    return rng.uniform(0.0, high, size=(n, ndim))


def _brute_force_knn(points, target, k):
    distances = np.linalg.norm(points - np.asarray(target, float).ravel()[:points.shape[1]],
                               axis=1)
    return np.sort(np.argsort(distances, kind="stable")[:k])


# ── 13. k-ən-yaxın ─────────────────────────────────────────────────────
def test_knn_matches_brute_force():
    points = _random_points(200, seed=1)
    target = np.array([[213.0, 371.0]])
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(max_neighbors=12))
    result = selector.select(target)
    assert result.count == 12
    assert np.array_equal(np.sort(result.indices), _brute_force_knn(points, target, 12))
    assert result.status == STATUS_KNN
    assert np.all(np.diff(result.distances) >= -1e-12), "məsafəyə görə sıralı olmalıdır"


def test_knn_returns_all_points_when_fewer_than_requested():
    points = _random_points(5, seed=2)
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(max_neighbors=50))
    assert selector.select(np.array([[100.0, 100.0]])).count == 5


def test_kdtree_and_brute_index_modes_agree():
    """İndeks strategiyası PERFORMANS seçimidir — NƏTİCƏYƏ təsir etməməlidir."""
    points = _random_points(150, seed=3)
    targets = _random_points(20, seed=4)
    config = NeighborhoodConfig(max_neighbors=9, search_radius=180.0, min_neighbors=2)
    tree = NeighborhoodSelector(points, config=config, index="kdtree")
    brute = NeighborhoodSelector(points, config=config, index="brute")
    for row in range(targets.shape[0]):
        a = tree.select(targets[row:row + 1])
        b = brute.select(targets[row:row + 1])
        assert np.array_equal(a.indices, b.indices)
        assert np.allclose(a.distances, b.distances)
        assert a.status == b.status and a.support == b.support


def test_min_neighbors_rejects_and_reports_status():
    points = _random_points(3, seed=5)
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(min_neighbors=5))
    result = selector.select(np.array([[10.0, 10.0]]))
    assert result.count == 0
    assert result.status == STATUS_INSUFFICIENT
    assert any("min_neighbors" in w for w in result.warnings)


# ── 14. radius ─────────────────────────────────────────────────────────
def test_radius_search_returns_exactly_the_points_inside():
    points = _random_points(300, seed=6)
    target = np.array([[250.0, 250.0]])
    radius = 90.0
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(
        search_radius=radius, min_neighbors=1))
    result = selector.select(target)
    expected = np.where(np.linalg.norm(points - target[0], axis=1) <= radius)[0]
    assert np.array_equal(np.sort(result.indices), np.sort(expected))
    assert result.status == STATUS_RADIUS
    assert result.radius_used == pytest.approx(radius)


def test_radius_is_deterministic_across_repeated_queries():
    points = _random_points(120, seed=7)
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(
        search_radius=120.0, max_neighbors=10))
    target = np.array([[300.0, 200.0]])
    first = selector.select(target).indices
    second = selector.select(target).indices
    assert np.array_equal(first, second)


def test_tied_distances_are_broken_by_index_deterministically():
    """Tam bərabər məsafəli 4 nöqtə — seçim TƏKRAR icralarda eynidir."""
    points = np.array([[10., 0.], [-10., 0.], [0., 10.], [0., -10.]])
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(max_neighbors=2))
    target = np.array([[0.0, 0.0]])
    assert np.array_equal(selector.select(target).indices,
                          selector.select(target).indices)
    assert np.array_equal(selector.select(target).indices, np.array([0, 1]))


# ── 15. anizotrop qonşu sıralaması ─────────────────────────────────────
def test_anisotropy_changes_neighbour_ranking():
    """Major ox X (azimut 90°) boyunca 5 dəfə uzundur — X-də 200 m
    uzaqdakı nöqtə Y-də 60 m uzaqdakından DAHA yaxın sayılmalıdır."""
    points = np.array([[200., 0.], [0., 60.]])
    target = np.array([[0.0, 0.0]])

    isotropic = NeighborhoodSelector(points, config=NeighborhoodConfig(max_neighbors=1))
    assert isotropic.select(target).indices[0] == 1        # 60 < 200

    aniso = AnisotropyParams(azimuth_deg=90.0, range_major=500.0, range_minor=100.0,
                             range_vertical=500.0)
    anisotropic = NeighborhoodSelector(points, anisotropy=aniso,
                                       config=NeighborhoodConfig(max_neighbors=1))
    assert anisotropic.select(target).indices[0] == 0      # 200 vs 60·5 = 300


def test_anisotropic_radius_uses_transformed_distance():
    points = np.array([[400., 0.], [0., 120.]])
    aniso = AnisotropyParams(azimuth_deg=90.0, range_major=500.0, range_minor=100.0,
                             range_vertical=500.0)
    selector = NeighborhoodSelector(points, anisotropy=aniso,
                                    config=NeighborhoodConfig(search_radius=450.0))
    # transformasiya edilmiş məsafələr: 400 (major) və 120·5 = 600 (minor)
    result = selector.select(np.array([[0.0, 0.0]]))
    assert np.array_equal(result.indices, np.array([0]))


# ── 16-17. kvadrant / oktant balanslaşdırması ──────────────────────────
def _clustered_points():
    """Bir kvadrantda SIX klaster (20 nöqtə) + qalan üç kvadrantda 3-3."""
    rng = np.random.default_rng(99)
    cluster = rng.uniform(5.0, 25.0, size=(20, 2))                  # (+,+)
    left = np.column_stack([rng.uniform(-90.0, -40.0, 3), rng.uniform(20.0, 70.0, 3)])
    down = np.column_stack([rng.uniform(20.0, 70.0, 3), rng.uniform(-90.0, -40.0, 3)])
    both = np.column_stack([rng.uniform(-90.0, -40.0, 3), rng.uniform(-90.0, -40.0, 3)])
    return np.vstack([cluster, left, down, both])


def test_quadrant_balancing_spreads_neighbours_across_sectors():
    points = _clustered_points()
    target = np.array([[0.0, 0.0]])

    plain = NeighborhoodSelector(points, config=NeighborhoodConfig(max_neighbors=8))
    balanced = NeighborhoodSelector(points, config=NeighborhoodConfig(
        max_neighbors=8, sectors=4))

    def quadrants(indices):
        offsets = points[indices]
        return {(int(x >= 0), int(y >= 0)) for x, y in offsets}

    plain_quadrants = quadrants(plain.select(target).indices)
    balanced_quadrants = quadrants(balanced.select(target).indices)
    assert len(plain_quadrants) == 1, "balanslaşdırmasız hamısı klasterdən gəlir"
    assert len(balanced_quadrants) == 4, "balanslaşdırma dörd kvadrantı da əhatə etməlidir"
    assert balanced.select(target).count == 8


def test_sector_balancing_does_not_invent_data_for_empty_sectors():
    """İki kvadrant tamamilə boşdur — seçici oradan nöqtə UYDURMUR,
    qalan yerləri mövcud ən yaxınlarla doldurur."""
    rng = np.random.default_rng(5)
    points = np.vstack([rng.uniform(10.0, 90.0, size=(15, 2)),          # (+,+)
                        np.column_stack([rng.uniform(-90.0, -10.0, 5),
                                         rng.uniform(10.0, 90.0, 5)])])  # (-,+)
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(
        max_neighbors=10, sectors=4))
    result = selector.select(np.array([[0.0, 0.0]]))
    assert result.count == 10
    assert set(result.indices).issubset(set(range(points.shape[0])))


def test_octant_balancing_in_3d_spreads_above_and_below():
    """`sectors=4, vertical_sectors=True` → 8 oktant. Üstdə 20, altda 4
    nöqtə var; balanslaşdırma altdakıları da götürməlidir."""
    rng = np.random.default_rng(17)
    above = np.column_stack([rng.uniform(-60.0, 60.0, size=(20, 2)),
                             rng.uniform(5.0, 30.0, 20)])
    below = np.column_stack([rng.uniform(-60.0, 60.0, size=(4, 2)),
                             rng.uniform(-30.0, -5.0, 4)])
    points = np.vstack([above, below])
    target = np.array([[0.0, 0.0, 0.0]])

    plain = NeighborhoodSelector(points, config=NeighborhoodConfig(max_neighbors=8))
    balanced = NeighborhoodSelector(points, config=NeighborhoodConfig(
        max_neighbors=8, sectors=4, vertical_sectors=True))

    plain_below = int(np.sum(points[plain.select(target).indices, 2] < 0))
    balanced_below = int(np.sum(points[balanced.select(target).indices, 2] < 0))
    assert balanced_below > plain_below
    assert balanced_below >= 2


def test_max_per_sector_caps_each_sector():
    points = _clustered_points()
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(
        max_neighbors=12, sectors=4, max_per_sector=2))
    result = selector.select(np.array([[0.0, 0.0]]))
    chosen = points[result.indices]
    ids = [(int(x >= 0), int(y >= 0)) for x, y in chosen]
    for quadrant in set(ids):
        assert ids.count(quadrant) <= 2


# ── A2.6 şaquli ayırma ────────────────────────────────────────────────
def test_vertical_cut_removes_geologically_distant_points():
    """Üfüqi cəhətdən ÇOX yaxın (5 m), amma 300 m dərində olan nöqtə
    `max_vertical_distance=50` ilə TAMAMILƏ kəsilir."""
    points = np.array([[5., 0., 0.], [5., 0., 300.], [80., 0., 10.]])
    target = np.array([[0.0, 0.0, 0.0]])
    without = NeighborhoodSelector(points, config=NeighborhoodConfig(max_neighbors=3))
    with_cut = NeighborhoodSelector(points, config=NeighborhoodConfig(
        max_neighbors=3, max_vertical_distance=50.0))
    assert without.select(target).count == 3
    kept = with_cut.select(target)
    assert kept.count == 2
    assert 1 not in set(kept.indices)


def test_vertical_scaling_and_vertical_cut_are_independent():
    """Anizotrop şaquli miqyaslanma məsafəni DƏYİŞİR, `max_vertical_
    distance` isə XAM Z ilə KƏSİR — ikisi eyni şey deyil."""
    points = np.array([[10., 0., 0.], [0., 0., 40.]])
    aniso = AnisotropyParams(range_major=100.0, range_minor=100.0, range_vertical=10.0)
    scaled = NeighborhoodSelector(points, anisotropy=aniso,
                                  config=NeighborhoodConfig(max_neighbors=1))
    assert scaled.select(np.array([[0., 0., 0.]])).indices[0] == 0   # 10 < 40·10

    cut = NeighborhoodSelector(points, anisotropy=aniso, config=NeighborhoodConfig(
        max_neighbors=2, max_vertical_distance=20.0))
    assert cut.select(np.array([[0., 0., 0.]])).count == 1


# ── 18. seyrək-data ehtiyat zənciri ───────────────────────────────────
def test_radius_expansion_kicks_in_and_is_reported():
    points = np.array([[0., 0.], [300., 0.], [0., 300.], [300., 300.]])
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(
        search_radius=20.0, min_neighbors=2, max_radius_expansions=4,
        radius_expansion_factor=3.0))
    result = selector.select(np.array([[150.0, 150.0]]))
    assert result.count >= 2
    assert result.status == STATUS_RADIUS_EXPANDED
    assert result.radius_used > 20.0
    assert any("genişləndirildi" in w for w in result.warnings)


def test_max_search_radius_stops_the_expansion():
    points = np.array([[0., 0.], [900., 0.], [0., 900.], [900., 900.]])
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(
        search_radius=10.0, min_neighbors=1, max_radius_expansions=10,
        radius_expansion_factor=2.0, max_search_radius=40.0))
    result = selector.select(np.array([[450.0, 450.0]]))
    assert result.count == 0
    assert result.status == STATUS_INSUFFICIENT


def test_knn_fallback_is_opt_in_and_reported():
    points = np.array([[0., 0.], [500., 0.], [0., 500.]])
    target = np.array([[250.0, 250.0]])
    strict = NeighborhoodSelector(points, config=NeighborhoodConfig(
        search_radius=10.0, min_neighbors=1))
    lenient = NeighborhoodSelector(points, config=NeighborhoodConfig(
        search_radius=10.0, min_neighbors=1, max_neighbors=2, allow_knn_fallback=True))
    assert strict.select(target).count == 0
    fallback = lenient.select(target)
    assert fallback.count == 2
    assert fallback.status == STATUS_KNN_FALLBACK
    assert any("ehtiyat" in w for w in fallback.warnings)


def test_global_fallback_is_opt_in_and_loudly_reported():
    points = np.array([[0., 0.], [500., 0.], [0., 500.], [500., 500.]])
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(
        search_radius=5.0, min_neighbors=3, allow_global_fallback=True))
    result = selector.select(np.array([[250.0, 250.0]]))
    assert result.status == STATUS_GLOBAL
    assert result.count == 4
    assert any("QLOBAL" in w for w in result.warnings)


# ── 19. boş məlumat çoxluğu ───────────────────────────────────────────
def test_empty_dataset_returns_empty_result_not_exception():
    selector = NeighborhoodSelector(np.zeros((0, 3)))
    result = selector.select(np.array([[1.0, 2.0, 3.0]]))
    assert result.count == 0
    assert result.status == STATUS_EMPTY
    assert result.support == SUPPORT_EXTRAPOLATED
    assert np.isinf(result.nearest_distance)
    assert not result.ok


def test_single_point_dataset_works():
    selector = NeighborhoodSelector(np.array([[10.0, 20.0]]))
    result = selector.select(np.array([[0.0, 0.0]]))
    assert result.count == 1 and result.ok


# ── 20. ekstrapolyasiya aşkarlanması ──────────────────────────────────
def test_support_classification_uses_real_geometry_not_neighbour_count():
    """Eyni QONŞU SAYI, fərqli HƏNDƏSƏ → fərqli təsnifat."""
    rng = np.random.default_rng(23)
    points = rng.uniform(0.0, 400.0, size=(60, 2))
    config = NeighborhoodConfig(max_neighbors=8, support_range=200.0)
    selector = NeighborhoodSelector(points, config=config)

    inside = selector.select(np.array([[200.0, 200.0]]))
    edge = selector.select(np.array([[420.0, 200.0]]))
    far = selector.select(np.array([[3000.0, 3000.0]]))

    assert inside.count == edge.count == far.count == 8
    assert inside.support == SUPPORT_WELL
    assert edge.support == SUPPORT_BOUNDARY
    assert far.support == SUPPORT_EXTRAPOLATED
    assert far.is_extrapolation and not inside.is_extrapolation


def test_sector_occupancy_is_reported():
    points = _random_points(80, seed=31, high=400.0)
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(
        max_neighbors=12, support_range=200.0))
    inside = selector.select(np.array([[200.0, 200.0]]))
    assert inside.n_sectors_total == 4        # 2D → kvadrantlar
    assert inside.n_sectors_occupied == 4


def test_support_scale_falls_back_to_data_density_when_not_given():
    """`support_range` verilməyəndə təsnifat MƏLUMATIN öz sıxlığından
    çıxarılan miqyasla aparılır — sabit UYDURULMUR."""
    points = _random_points(60, seed=41, high=300.0)
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(max_neighbors=6))
    assert selector.select(np.array([[150.0, 150.0]])).support == SUPPORT_WELL
    assert selector.select(np.array([[9000.0, 9000.0]])).support == SUPPORT_EXTRAPOLATED


def test_select_many_matches_select():
    points = _random_points(70, seed=51)
    targets = _random_points(15, seed=52)
    selector = NeighborhoodSelector(points, config=NeighborhoodConfig(max_neighbors=6))
    many = selector.select_many(targets)
    for row, result in enumerate(many):
        assert np.array_equal(result.indices,
                              selector.select(targets[row:row + 1]).indices)


# ── konfiqurasiya doğrulaması ─────────────────────────────────────────
@pytest.mark.parametrize("kwargs", [
    {"min_neighbors": 0},
    {"max_neighbors": 0},
    {"min_neighbors": 5, "max_neighbors": 3},
    {"search_radius": -1.0},
    {"search_radius": 10.0, "max_search_radius": 5.0},
    {"sectors": 1},
    {"max_per_sector": 0},
    {"max_radius_expansions": 2, "radius_expansion_factor": 1.0},
    {"support_range": float("inf")},
])
def test_invalid_configuration_raises(kwargs):
    with pytest.raises(NeighborhoodError):
        NeighborhoodConfig(**kwargs).validate()


def test_invalid_index_mode_raises():
    with pytest.raises(NeighborhoodError):
        NeighborhoodSelector(_random_points(5, seed=1), index="octree")
