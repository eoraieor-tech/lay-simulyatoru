"""Phase 5 §14/§15/§17 — SGS ansambl paylanması, mövzu-davamlılıq
(variogram) və sürətli/brute-force axtarış PARİTETİ yoxlanması.

Bu fayl SOFTWARE CORRECTNESS-i sınayır (bax `sgs.py` modul docstring-
indəki elmi çəkincə) — GEOLOJİ REALLIQ İDDİA EDİLMİR.
"""

from __future__ import annotations

import numpy as np

from imex2d.geology.distribution_analysis import summarize_distribution
from imex2d.geology.sgs import PropertyVariogramParams, run_realizations_sgs, simulate_sgs
from imex2d.geology.variogram import fit_variogram_from_data


def _grid_targets(nx, ny, dx=20.0, dy=20.0):
    xs = (np.arange(nx) + 0.5) * dx
    ys = (np.arange(ny) + 0.5) * dy
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


# ── §17: sürətli (cKDTree) vs brute-force PARİTETİ — SGS səviyyəsində ────
def test_fast_and_brute_force_search_agree_on_identical_realization():
    """QEYRİ-grid (təsadüfi) hədəflərlə — TAM (bit-bə-bit) uyğunluq.
    (Grid-uyğunlaşmış hədəflər üçün bax `test_fast_and_brute_force_
    differ_only_at_grid_tie_boundaries_not_wrong` — orada FƏRQLİ, amma
    SƏHV OLMAYAN bir hal sənədləşdirilib.)"""
    rng = np.random.default_rng(20)
    points = rng.uniform(0, 300, size=(12, 2))
    values = rng.normal(0.20, 0.03, size=12)
    targets = rng.uniform(0, 300, size=(100, 2))

    fast = simulate_sgs(points, values, targets, seed=5, max_neighbors=8, use_fast_search=True)
    slow = simulate_sgs(points, values, targets, seed=5, max_neighbors=8, use_fast_search=False)
    assert np.allclose(fast.values, slow.values, atol=1e-9)


def test_fast_and_brute_force_agree_with_explicit_anisotropy():
    vp = PropertyVariogramParams(range_=250.0, range_minor=40.0, azimuth_deg=30.0, nugget=0.02)
    rng = np.random.default_rng(21)
    points = rng.uniform(0, 300, size=(10, 2))
    values = rng.normal(150.0, 30.0, size=10)
    targets = rng.uniform(0, 300, size=(64, 2))

    fast = simulate_sgs(points, values, targets, variogram=vp, seed=3, search_radius=200.0,
                        max_neighbors=6, use_fast_search=True)
    slow = simulate_sgs(points, values, targets, variogram=vp, seed=3, search_radius=200.0,
                        max_neighbors=6, use_fast_search=False)
    assert np.allclose(fast.values, slow.values, atol=1e-9)


def test_fast_and_brute_force_differ_only_at_grid_tie_boundaries_not_wrong():
    """DÜRÜST TAPINTI: müntəzəm (grid) hədəflərdə fast/brute-force TAM
    üst-üstə düşməyə bilər — REGULYAR şəbəkədə `max_neighbors` sərhədində
    dəqiq/yaxın-bərabər məsafəli namizədlər TEZ-TEZ olur, cKDTree-nin
    daxili sıralaması ilə numpy-nin sabit `argsort`-u bu cür bərabərliyi
    FƏRQLİ (hər ikisi DÜZGÜN — eyni məsafədə) namizədlə həll edə bilər.
    Bu, `test_spatial_search.py`-də (qeyri-grid, təsadüfi nöqtələrlə)
    TAM uyğunluq sübut edildikdən SONRA aşkarlanan, sənədləşdirilmiş bir
    NÜANSDIR — YALNIZ qonşuluq SƏRHƏDİNDƏKİ bərabərliyə aiddir, alqoritm
    SƏHVİ deyil (bax FACIES.md/SGS.md üçün tam izah)."""
    rng = np.random.default_rng(20)
    points = rng.uniform(0, 300, size=(12, 2))
    values = rng.normal(0.20, 0.03, size=12)
    targets = _grid_targets(10, 10, dx=25.0, dy=25.0)

    fast = simulate_sgs(points, values, targets, seed=5, max_neighbors=8, use_fast_search=True)
    slow = simulate_sgs(points, values, targets, seed=5, max_neighbors=8, use_fast_search=False)
    # TAM bərabərlik GÖZLƏNİLMİR (yuxarı bax) — YALNIZ ki hər ikisi
    # MƏQBUL (kiçik) fərqlə eyni FİZİKİ diapazonda qalır
    assert np.allclose(fast.values, slow.values, atol=0.02)
    assert abs(np.mean(fast.values) - np.mean(slow.values)) < 0.01


# ── §14: ensambl paylanması — P10/P50/P90 müqayisəsi ────────────────────
def test_ensemble_p10_p50_p90_are_within_reasonable_range_of_target():
    rng = np.random.default_rng(22)
    points = rng.uniform(0, 300, size=(20, 2))
    values = rng.normal(0.22, 0.025, size=20)
    target_summary = summarize_distribution(values)
    targets = _grid_targets(12, 12, dx=20.0, dy=20.0)

    realizations = run_realizations_sgs(20, points, values, targets, seed=0)
    all_values = np.concatenate([r.values for r in realizations])
    ensemble_summary = summarize_distribution(all_values)

    print(f"\n[SGS ensemble] target P10/P50/P90 = {target_summary.p10:.4f}/"
         f"{target_summary.p50:.4f}/{target_summary.p90:.4f}  "
         f"ensemble = {ensemble_summary.p10:.4f}/{ensemble_summary.p50:.4f}/"
         f"{ensemble_summary.p90:.4f}")
    # geniş, DÜRÜST tolerantlıq — SGS kondisioner+kriging qeyri-
    # müəyyənliyini ƏLAVƏ edir, dəqiq bərabərlik GÖZLƏNİLMİR
    assert abs(ensemble_summary.p50 - target_summary.p50) < 0.03
    assert ensemble_summary.p10 < ensemble_summary.p50 < ensemble_summary.p90


# ── §15: mövzu-davamlılıq (variogram) reproduksiyası ────────────────────
def test_simulated_field_reproduces_approximate_variogram_range():
    """SGS-in nəticə SAHƏSİNDƏN hesablanmış deneysel variogram, HƏDƏF
    (giriş) variogram parametrlərinə YAXIN olmalıdır — DƏQİQ reproduksiya
    İDDİA EDİLMİR, YALNIZ kobud miqyas uyğunluğu yoxlanılır."""
    true_range = 80.0
    true_sill = 0.001
    rng = np.random.default_rng(23)
    # sıx bir sərt-data şəbəkəsi qurmaq üçün əvvəlcə GENİŞ bir SGS
    # icrası ilə "yer həqiqəti" sahəsi yaradılır, sonra bu sahədən
    # seyrək nöqtələr SEÇİLİR (dövri asılılıq YOXDUR: aşağıdakı test
    # YALNIZ bu sahədən yenidən fit edilən variogramın MƏNTİQİ miqyasda
    # olduğunu yoxlayır, DƏQİQ ƏDƏD gözləmir).
    dense_targets = _grid_targets(30, 30, dx=10.0, dy=10.0)
    seed_points = rng.uniform(0, 300, size=(10, 2))
    seed_values = rng.normal(0.20, 0.02, size=10)
    vp = PropertyVariogramParams(range_=true_range, sill=true_sill, nugget=0.0)
    ground_truth = simulate_sgs(seed_points, seed_values, dense_targets, variogram=vp, seed=1,
                                max_neighbors=20)

    exp_points_xy = dense_targets
    fit = fit_variogram_from_data(exp_points_xy, ground_truth.values, model="spherical")
    print(f"\n[SGS variogram check] true range={true_range}, fitted range={fit.range_:.1f}; "
         f"true sill={true_sill:.5f}, fitted sill={fit.sill:.6f}")
    # geniş tolerantlıq (sonlu nümunə + kriging-in özünün hamarlaşdırma
    # meyli) — YALNIZ eyni BÖYÜKLÜK dərəcəsində olduğunu yoxlayır
    assert 0.3 * true_range < fit.range_ < 3.0 * true_range
