"""Phase 5 — `OrdinaryKriging.interpolate_with_variance()` (SGS üçün YENİ,
`interpolate()`-i DƏYİŞDİRMƏYƏN əlavə). Varians düsturu: σ²(x0) = Σ w_i·
γ(x_i,x0) + μ (Laqranj vuruğu) — ARTIQ həll edilmiş sistemdən oxunur."""

from __future__ import annotations

import numpy as np

from imex2d.geology.interpolation import OrdinaryKriging

POINTS = np.array([[0., 0.], [100., 0.], [0., 100.], [100., 100.], [50., 50.]])
VALUES = np.array([0.15, 0.25, 0.20, 0.30, 0.22])


def test_interpolate_unchanged_by_new_method_existing():
    """`interpolate()`-in ÖZÜ (Phase 2-3, 76+ test) dəyişməyib — eyni
    parametrlərlə eyni nəticəni verir, YENİ metoddan ASILI OLMADAN."""
    kriging = OrdinaryKriging(nugget=0.05)
    result = kriging.interpolate(POINTS, VALUES, np.array([[25., 25.], [75., 75.]]))
    estimate, _ = kriging.interpolate_with_variance(POINTS, VALUES,
                                                     np.array([[25., 25.], [75., 75.]]))
    assert np.allclose(result, estimate)


def test_variance_is_zero_at_exact_data_points_with_zero_nugget():
    kriging = OrdinaryKriging(nugget=0.0)
    estimate, variance = kriging.interpolate_with_variance(POINTS, VALUES, POINTS)
    assert np.allclose(estimate, VALUES, atol=1e-9)
    assert np.allclose(variance, 0.0, atol=1e-6)


def test_variance_increases_with_distance_from_data():
    """Verilənlərdən uzaqlaşdıqca kriging qeyri-müəyyənliyi (varians)
    ARTMALIDIR — standart kriging nəzəriyyəsi."""
    kriging = OrdinaryKriging(range_=200.0, sill=0.01, nugget=0.001)
    near = np.array([[10., 10.]])     # nöqtələrə yaxın
    far = np.array([[500., 500.]])    # nöqtələrdən uzaq
    _, var_near = kriging.interpolate_with_variance(POINTS, VALUES, near)
    _, var_far = kriging.interpolate_with_variance(POINTS, VALUES, far)
    assert var_far[0] > var_near[0]


def test_variance_is_non_negative():
    rng = np.random.default_rng(0)
    points = rng.uniform(0, 300, size=(15, 2))
    values = rng.uniform(0.1, 0.3, size=15)
    kriging = OrdinaryKriging(auto_fit=True, nugget=0.02)
    targets = rng.uniform(0, 300, size=(20, 2))
    _, variance = kriging.interpolate_with_variance(points, values, targets)
    assert np.all(variance >= -1e-9)


def test_variance_respects_local_search_neighbors():
    """Yerli axtarışla (`search_radius`/`max_neighbors`) da varians
    hesablanmalıdır — NaN yalnız kifayət qədər qonşu olmayanda."""
    kriging = OrdinaryKriging(range_=200.0, nugget=0.02, search_radius=60.0, min_neighbors=1)
    close_target = np.array([[10., 10.]])     # (0,0) nöqtəsinə yaxın, radiusda
    far_target = np.array([[900., 900.]])     # heç bir nöqtə radiusda yoxdur
    estimate, variance = kriging.interpolate_with_variance(POINTS, VALUES,
                                                           np.vstack([close_target, far_target]))
    assert np.isfinite(variance[0])
    assert np.isnan(estimate[1]) and np.isnan(variance[1])


def test_single_point_variance_uses_semivariogram_distance():
    kriging = OrdinaryKriging(range_=100.0, sill=0.02, nugget=0.0)
    single_point = np.array([[0., 0.]])
    single_value = np.array([0.2])
    estimate, variance = kriging.interpolate_with_variance(
        single_point, single_value, np.array([[0., 0.], [50., 0.], [1000., 0.]]))
    assert estimate[0] == 0.2 and estimate[1] == 0.2 and estimate[2] == 0.2
    assert variance[0] == 0.0                 # öz nöqtəsində
    assert variance[2] > variance[1] > 0.0    # uzaqlıqla artır
