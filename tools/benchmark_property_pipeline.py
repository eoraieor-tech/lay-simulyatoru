"""Phase B boru xəttinin performans ölçməsi (B10).

    python tools/benchmark_property_pipeline.py
    python tools/benchmark_property_pipeline.py --samples 100,500,1000,5000
    python tools/benchmark_property_pipeline.py --grids 41,101 --skip-cv

Ölçülənlər (hər sərt-data ölçüsü üçün):

    QC          `data_quality.run_quality_control` (dublikat, hədd, kənar-dəyər)
    VARİOQRAM   deneysel + model fit (`variogram.fit_variogram_from_data`)
    İNTERP      TAM xassə boru xətti (`interpolate_property_field`)
    CV          leave-one-out çarpaz-doğrulama — QƏSDƏN BAHA, OFFLINE
    SGS         bir realizasiya (stoxastik budaq)
    YADDAŞ      `tracemalloc` pik

VACİB AYRIM (B10): LOOCV `n` dəfə tam boru xətti işlədir, yəni `O(n)`
dəfə bahadır. Bu, MODEL SEÇİMİ/kalibrləmə üçündür və İSTEHSAL şəbəkə
yolunda YOXDUR — `interpolate_property_field` şəbəkə hüceyrələri üçün
Phase A-nın YERLİ kriginqini bir dəfə çağırır. Cədvəldə ikisi AYRI
sütunlardır ki, qarışdırılmasın.

Ədədlər UYDURULMUR: skript nə ölçübsə, onu çap edir.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import tracemalloc

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from imex2d.geology.cross_validation import (ValidationDesign, ValidationKind,
                                             cross_validate_property)
from imex2d.geology.data_quality import run_quality_control
from imex2d.geology.property_config import resolve_strategy
from imex2d.geology.property_interpolation import interpolate_property_field
from imex2d.geology.sgs import PropertyVariogramParams
from imex2d.geology.sgs_ensemble import simulate_sgs_ensemble
from imex2d.geology.variogram import fit_variogram_from_data

DOMAIN = 5000.0
RANGE = 1200.0


def make_wells(n: int, seed: int = 0):
    """`n` quyu + loq-normal keçiricilik (ən bahalı yol: loq çevirməsi)."""
    rng = np.random.default_rng(seed)
    points = rng.uniform(0.0, DOMAIN, size=(n, 2))
    values = np.exp(4.0 + 1.2 * np.sin(points[:, 0] / 900.0)
                    + 0.4 * rng.standard_normal(n))
    return points, values


def make_grid(nx: int, ny: int):
    xs = (np.arange(nx) + 0.5) * (DOMAIN / nx)
    ys = (np.arange(ny) + 0.5) * (DOMAIN / ny)
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


def timed(callable_, *args, **kwargs):
    tracemalloc.start()
    start = time.perf_counter()
    result = callable_(*args, **kwargs)
    elapsed = time.perf_counter() - start
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, elapsed, peak / (1024.0 * 1024.0)


def benchmark(n_samples: int, nx: int, ny: int, cv_limit: int, skip_cv: bool,
              skip_sgs: bool, sgs_limit: int):
    points, values = make_wells(n_samples, seed=n_samples)
    targets = make_grid(nx, ny)
    strategy = resolve_strategy("PERMX")

    _, t_qc, mem_qc = timed(run_quality_control, points, values, strategy)
    _, t_variogram, _ = timed(fit_variogram_from_data, points, np.log(values))
    _, t_interp, mem_interp = timed(interpolate_property_field, points, values,
                                    targets, strategy=strategy)

    t_cv = float("nan")
    if not skip_cv and n_samples <= cv_limit:
        _, t_cv, _ = timed(
            cross_validate_property, points, values, strategy,
            ValidationDesign(kind=ValidationKind.LEAVE_ONE_OUT))

    t_sgs = float("nan")
    if not skip_sgs and n_samples <= sgs_limit:
        _, t_sgs, _ = timed(
            simulate_sgs_ensemble, 1, points, values, targets,
            variogram=PropertyVariogramParams(model="spherical", nugget=0.0,
                                              range_=RANGE),
            base_seed=0, max_neighbors=16)

    return {"samples": n_samples, "cells": targets.shape[0], "grid": f"{nx}x{ny}",
            "qc": t_qc, "variogram": t_variogram, "interp": t_interp,
            "cv": t_cv, "sgs": t_sgs, "memory": max(mem_qc, mem_interp)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", default="100,500,1000,5000")
    parser.add_argument("--grids", default="41,101")
    parser.add_argument("--cv-limit", type=int, default=1000,
                        help="LOOCV bu nöqtə sayına qədər ölçülür (O(n) baha)")
    parser.add_argument("--sgs-limit", type=int, default=500,
                        help="SGS bu nöqtə sayına qədər ölçülür")
    parser.add_argument("--skip-cv", action="store_true")
    parser.add_argument("--skip-sgs", action="store_true")
    args = parser.parse_args()

    samples = [int(v) for v in args.samples.split(",") if v.strip()]
    grids = [int(v) for v in args.grids.split(",") if v.strip()]

    print(f"Domen {DOMAIN:.0f} m · xassə PERMX (loq fəzası) · radius {RANGE:.0f} m")
    header = (f"{'nöqtə':>7} {'şəbəkə':>10} {'hüceyrə':>9} {'QC s':>8} "
              f"{'variogram s':>12} {'interp s':>9} {'µs/hüceyrə':>11} "
              f"{'LOOCV s':>9} {'SGS s':>8} {'yaddaş MB':>10}")
    print(header)
    print("-" * len(header))

    for n in samples:
        for size in grids:
            row = benchmark(n, size, size, args.cv_limit, args.skip_cv,
                            args.skip_sgs, args.sgs_limit)
            per_cell = 1e6 * row["interp"] / max(row["cells"], 1)
            cv = "—" if not np.isfinite(row["cv"]) else f"{row['cv']:.2f}"
            sgs = "—" if not np.isfinite(row["sgs"]) else f"{row['sgs']:.2f}"
            print(f"{row['samples']:>7} {row['grid']:>10} {row['cells']:>9} "
                  f"{row['qc']:>8.4f} {row['variogram']:>12.4f} "
                  f"{row['interp']:>9.3f} {per_cell:>11.1f} {cv:>9} {sgs:>8} "
                  f"{row['memory']:>10.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
