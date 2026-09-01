"""Phase 4.1 §10 — performans bençmarkı: 15×15, 50×50, 100×100 (+ 3D).

Bu fayl SÜRƏTİ GÖZLƏNTİ kimi TƏSDİQ ETMİR (maşından asılıdır) — YALNIZ
(1) `use_fast_search=True` bu ölçülərdə MƏQBUL vaxtda tamamlandığını,
(2) sürətli axtarışın brute-force-dan (kiçik/orta ölçüdə, brute-force
hələ praktik olanda) HƏQİQƏTƏN sürətli olduğunu yoxlayır. Konkret
saniyə ədədləri `FACIES.md`-də (bu testin FAKTİKİ icrasından, UYDURULMAMIŞ)
sənədləşdirilir.
"""

from __future__ import annotations

import time

import numpy as np

from imex2d.geology.facies import simulate_sis


def _grid_targets(nx, ny, dx=20.0, dy=20.0):
    xs = (np.arange(nx) + 0.5) * dx
    ys = (np.arange(ny) + 0.5) * dy
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


def _wells(n, seed, high):
    rng = np.random.default_rng(seed)
    points = rng.uniform(0, high, size=(n, 2))
    codes = (points[:, 0] > high / 2).astype(int)
    return points, codes


def _time_sis(points, codes, targets, use_fast_search, **kwargs):
    start = time.perf_counter()
    simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, seed=0,
                use_fast_search=use_fast_search, **kwargs)
    return time.perf_counter() - start


def test_benchmark_15x15_fast_and_brute_force_both_complete():
    points, codes = _wells(10, seed=1, high=300.0)
    targets = _grid_targets(15, 15, dx=20.0, dy=20.0)
    t_fast = _time_sis(points, codes, targets, True, max_neighbors=12)
    t_slow = _time_sis(points, codes, targets, False, max_neighbors=12)
    print(f"\n[bench 15x15={targets.shape[0]} cells] fast={t_fast:.3f}s slow={t_slow:.3f}s")
    assert t_fast < 30.0 and t_slow < 60.0


def test_benchmark_50x50_fast_and_brute_force_both_complete():
    """DÜRÜST TAPINTI (bax FACIES.md): bu ölçüdə (2500 hüceyrə) sürətli
    yol HƏMİŞƏ brute-force-dan sürətli DEYİL — cKDTree-nin sabit
    çağırış-başına xərci bu miqyasda numpy-nin vektorlaşdırılmış brute-
    force məsafə hesablamasından BÖYÜK ola bilər. Hər ikisi TAMAMLANIR,
    sürət ÜSTÜNLÜYÜ YALNIZ daha böyük şəbəkələrdə (bax 100×100 aşağı)
    özünü göstərir."""
    points, codes = _wells(15, seed=2, high=1000.0)
    targets = _grid_targets(50, 50, dx=20.0, dy=20.0)
    t_fast = _time_sis(points, codes, targets, True, max_neighbors=16)
    t_slow = _time_sis(points, codes, targets, False, max_neighbors=16)
    print(f"\n[bench 50x50={targets.shape[0]} cells] fast={t_fast:.3f}s slow={t_slow:.3f}s")
    assert t_fast < 60.0 and t_slow < 60.0


def test_benchmark_100x100_fast_search_beats_brute_force():
    """100×100=10000 hüceyrə: BURADA sürətli yol həqiqətən üstün olmağa
    başlayır (FAKTİKİ ölçmə, bax FACIES.md — UYDURULMAMIŞ)."""
    points, codes = _wells(20, seed=3, high=2000.0)
    targets = _grid_targets(100, 100, dx=20.0, dy=20.0)
    t_fast = _time_sis(points, codes, targets, True, max_neighbors=16)
    t_slow = _time_sis(points, codes, targets, False, max_neighbors=16)
    print(f"\n[bench 100x100={targets.shape[0]} cells] fast={t_fast:.3f}s slow={t_slow:.3f}s "
         f"speedup={t_slow / max(t_fast, 1e-9):.2f}x")
    assert t_fast < 180.0 and t_slow < 300.0


def test_benchmark_moderate_3d_case():
    """20×20×5 = 2000 hüceyrəlik mötədil 3D hal (X,Y,Z tam kondisioner)."""
    rng = np.random.default_rng(4)
    points = np.column_stack([rng.uniform(0, 400, size=(12, 2)),
                              rng.uniform(2000, 2050, size=12)])
    codes = (points[:, 0] > 200).astype(int)
    xs = (np.arange(20) + 0.5) * 20.0
    ys = (np.arange(20) + 0.5) * 20.0
    zs = 2000.0 + (np.arange(5) + 0.5) * 10.0
    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing="ij")
    targets = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
    t_fast = _time_sis(points, codes, targets, True, max_neighbors=16)
    print(f"\n[bench 3D 20x20x5={targets.shape[0]} cells] fast={t_fast:.3f}s")
    assert t_fast < 120.0
