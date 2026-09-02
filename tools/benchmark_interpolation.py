"""İnterpolyasiya özəyinin performans ölçməsi (A7).

    python tools/benchmark_interpolation.py
    python tools/benchmark_interpolation.py --samples 100,500,1000,5000
    python tools/benchmark_interpolation.py --grids 41,101 --neighbors 16

Ölçülənlər (hər sərt-data ölçüsü × şəbəkə ölçüsü üçün):

    * QONŞULUQ vaxtı — `spatial_search.NeighborhoodSelector.select_batch`
      (cKDTree qurma + TOPLU sorğu/seçim — istehsal yolunun eynisi);
    * KRİGİNG vaxtı — yerli sistemlərin qurulması + toplu həlli;
    * TAM vaxt — `OrdinaryKriging.interpolate()` (boru xəttinin hamısı);
    * QLOBAL vaxt — müqayisə üçün, YALNIZ kiçik nöqtə sayında (böyük
      çoxluqda `O(n³)` və `O(n²)` yaddaş praktik deyil — məhz buna görə
      istehsal yolu YERLİDİR, bax `interpolation.OrdinaryKriging`
      `auto_local_threshold`);
    * PİK YADDAŞ — `tracemalloc` ilə (Python obyektləri + NumPy buferləri).

Ədədlər UYDURULMUR: skript nə ölçübsə, onu çap edir. Nəticələr
`PERFORMANCE.md`-də sənədləşdirilir.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import tracemalloc

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from imex2d.geology.anisotropy import AnisotropyParams
from imex2d.geology.interpolation import OrdinaryKriging
from imex2d.geology.spatial_search import NeighborhoodConfig, NeighborhoodSelector

DOMAIN = 5000.0
RANGE_MAJOR = 1200.0
RANGE_MINOR = 400.0
AZIMUTH = 35.0


def make_samples(n: int, seed: int = 0):
    """`n` sərt data nöqtəsi (X,Y,Z) + hamar sintetik dəyər."""
    rng = np.random.default_rng(seed)
    points = np.column_stack([rng.uniform(0.0, DOMAIN, size=(n, 2)),
                              rng.uniform(2000.0, 2100.0, n)])
    values = (0.20 + 0.05 * np.sin(points[:, 0] / 800.0)
              + 0.03 * np.cos(points[:, 1] / 600.0)
              + 0.0002 * (points[:, 2] - 2050.0))
    return points, values


def make_grid(nx: int, ny: int, nz: int = 1):
    xs = (np.arange(nx) + 0.5) * (DOMAIN / nx)
    ys = (np.arange(ny) + 0.5) * (DOMAIN / ny)
    zs = 2000.0 + (np.arange(nz) + 0.5) * (100.0 / nz)
    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])


def timed(callable_, *args, **kwargs):
    """`(nəticə, saniyə, pik_MB)` — pik yaddaş `tracemalloc` ilə."""
    tracemalloc.start()
    start = time.perf_counter()
    result = callable_(*args, **kwargs)
    elapsed = time.perf_counter() - start
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, elapsed, peak / (1024.0 * 1024.0)


def benchmark(n_samples: int, nx: int, ny: int, nz: int, max_neighbors: int,
              include_global: bool, global_limit: int):
    points, values = make_samples(n_samples, seed=n_samples)
    targets = make_grid(nx, ny, nz)
    anisotropy = AnisotropyParams(azimuth_deg=AZIMUTH, range_major=RANGE_MAJOR,
                                  range_minor=RANGE_MINOR, range_vertical=30.0)

    # ── 1. yalnız qonşuluq (indeks qurma + sorğular) ──────────────────
    transformed_points = anisotropy.transform(points)
    transformed_targets = anisotropy.transform(targets)
    config = NeighborhoodConfig(max_neighbors=max_neighbors, min_neighbors=1,
                                support_range=RANGE_MAJOR)

    def neighbourhood_pass():
        selector = NeighborhoodSelector(transformed_points, config=config)
        return selector.select_batch(transformed_targets)

    batch, t_neighbours, mem_neighbours = timed(neighbourhood_pass)
    mean_neighbours = float(np.mean(batch.counts))

    # ── 2. tam yerli boru xətti ───────────────────────────────────────
    local = OrdinaryKriging(range_=RANGE_MAJOR, range_minor=RANGE_MINOR,
                            range_v=30.0, azimuth_deg=AZIMUTH, sill=0.01,
                            nugget=0.0, max_neighbors=max_neighbors)
    _, t_local, mem_local = timed(local.interpolate, points, values, targets)

    # ── 3. qlobal sistem (yalnız kiçik n üçün) ────────────────────────
    t_global = mem_global = float("nan")
    if include_global and n_samples <= global_limit:
        glob = OrdinaryKriging(range_=RANGE_MAJOR, range_minor=RANGE_MINOR,
                               range_v=30.0, azimuth_deg=AZIMUTH, sill=0.01,
                               nugget=0.0, auto_local_threshold=10 ** 9)
        _, t_global, mem_global = timed(glob.interpolate, points, values, targets)

    return {
        "samples": n_samples,
        "cells": targets.shape[0],
        "grid": f"{nx}x{ny}x{nz}",
        "neighbours": mean_neighbours,
        "t_neighbourhood": t_neighbours,
        "t_kriging": max(t_local - t_neighbours, 0.0),
        "t_total": t_local,
        "t_global": t_global,
        "mem_local": max(mem_local, mem_neighbours),
        "mem_global": mem_global,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", default="100,500,1000,5000",
                        help="sərt data nöqtəsi sayları (vergüllə)")
    parser.add_argument("--grids", default="41,101",
                        help="kvadrat şəbəkə ölçüləri (vergüllə), nx=ny")
    parser.add_argument("--nz", type=int, default=1, help="lay sayı")
    parser.add_argument("--neighbors", type=int, default=24,
                        help="yerli sistemdə maksimum qonşu")
    parser.add_argument("--no-global", action="store_true",
                        help="qlobal sistemlə müqayisəni keç")
    parser.add_argument("--global-limit", type=int, default=1000,
                        help="qlobal sistem bu nöqtə sayına qədər sınanır")
    args = parser.parse_args()

    samples = [int(v) for v in args.samples.split(",") if v.strip()]
    grids = [int(v) for v in args.grids.split(",") if v.strip()]

    print(f"Domen {DOMAIN:.0f} m, anizotropluq: azimut {AZIMUTH}°, "
          f"major {RANGE_MAJOR:.0f} m, minor {RANGE_MINOR:.0f} m, "
          f"max_neighbors={args.neighbors}")
    header = (f"{'nöqtə':>7} {'şəbəkə':>12} {'hüceyrə':>9} {'qonşu':>6} "
              f"{'qonşuluq s':>11} {'kriging s':>10} {'CƏMİ s':>9} "
              f"{'µs/hüceyrə':>11} {'yaddaş MB':>10} {'qlobal s':>9}")
    print(header)
    print("-" * len(header))

    rows = []
    for n in samples:
        for size in grids:
            row = benchmark(n, size, size, args.nz, args.neighbors,
                            not args.no_global, args.global_limit)
            rows.append(row)
            per_cell = 1e6 * row["t_total"] / max(row["cells"], 1)
            global_text = ("—" if not np.isfinite(row["t_global"])
                           else f"{row['t_global']:.3f}")
            print(f"{row['samples']:>7} {row['grid']:>12} {row['cells']:>9} "
                  f"{row['neighbours']:>6.1f} {row['t_neighbourhood']:>11.3f} "
                  f"{row['t_kriging']:>10.3f} {row['t_total']:>9.3f} "
                  f"{per_cell:>11.1f} {row['mem_local']:>10.1f} {global_text:>9}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
