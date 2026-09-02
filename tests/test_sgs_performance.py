"""Phase 5 §20 — SGS performans bençmarkı: 15×15, 50×50, 100×100 (+ 3D).

`facies/test_facies_performance.py`-lə EYNİ dürüstlük qaydası: konkret
saniyə ədədləri BU İCRADAN (uydurulmamış) götürülür, "sürətli axtarış
HƏR ZAMAN daha sürətlidir" İDDİA EDİLMİR (bax FACIES_INTEGRATION.md-dəki
eyni tapıntı — kiçik/orta grid-də cKDTree-nin sabit xərci üstünlük təşkil
edə bilər).
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from imex2d.geology.sgs import simulate_sgs

pytestmark = pytest.mark.performance


def _grid_targets(nx, ny, dx=20.0, dy=20.0):
    xs = (np.arange(nx) + 0.5) * dx
    ys = (np.arange(ny) + 0.5) * dy
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


def _wells(n, seed, high):
    rng = np.random.default_rng(seed)
    points = rng.uniform(0, high, size=(n, 2))
    values = rng.normal(0.20, 0.03, size=n)
    return points, values


def _time(points, values, targets, use_fast_search, **kwargs):
    start = time.perf_counter()
    simulate_sgs(points, values, targets, seed=0, use_fast_search=use_fast_search, **kwargs)
    return time.perf_counter() - start


def test_benchmark_15x15():
    points, values = _wells(10, seed=1, high=300.0)
    targets = _grid_targets(15, 15, dx=20.0, dy=20.0)
    t_fast = _time(points, values, targets, True, max_neighbors=12)
    t_slow = _time(points, values, targets, False, max_neighbors=12)
    print(f"\n[SGS bench 15x15={targets.shape[0]} cells] fast={t_fast:.3f}s slow={t_slow:.3f}s")
    assert t_fast < 30.0 and t_slow < 60.0


def test_benchmark_50x50():
    points, values = _wells(15, seed=2, high=1000.0)
    targets = _grid_targets(50, 50, dx=20.0, dy=20.0)
    t_fast = _time(points, values, targets, True, max_neighbors=16)
    t_slow = _time(points, values, targets, False, max_neighbors=16)
    print(f"\n[SGS bench 50x50={targets.shape[0]} cells] fast={t_fast:.3f}s slow={t_slow:.3f}s")
    assert t_fast < 60.0 and t_slow < 60.0


def test_benchmark_100x100():
    points, values = _wells(20, seed=3, high=2000.0)
    targets = _grid_targets(100, 100, dx=20.0, dy=20.0)
    t_fast = _time(points, values, targets, True, max_neighbors=16)
    t_slow = _time(points, values, targets, False, max_neighbors=16)
    print(f"\n[SGS bench 100x100={targets.shape[0]} cells] fast={t_fast:.3f}s slow={t_slow:.3f}s "
         f"speedup={t_slow / max(t_fast, 1e-9):.2f}x")
    assert t_fast < 180.0 and t_slow < 300.0


def test_benchmark_moderate_3d_case():
    """20×20×5 = 2000 hüceyrəlik mötədil 3D hal."""
    rng = np.random.default_rng(4)
    points = np.column_stack([rng.uniform(0, 400, size=(12, 2)),
                              rng.uniform(2000, 2050, size=12)])
    values = rng.normal(0.20, 0.03, size=12)
    xs = (np.arange(20) + 0.5) * 20.0
    ys = (np.arange(20) + 0.5) * 20.0
    zs = 2000.0 + (np.arange(5) + 0.5) * 10.0
    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing="ij")
    targets = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
    t_fast = _time(points, values, targets, True, max_neighbors=16)
    print(f"\n[SGS bench 3D 20x20x5={targets.shape[0]} cells] fast={t_fast:.3f}s")
    assert t_fast < 120.0
