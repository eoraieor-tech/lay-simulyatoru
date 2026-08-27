"""Performans ölçmə aləti (C3).

    python tools/benchmark.py                 standart dəst
    python tools/benchmark.py --profile 41    bir ölçü üçün cProfile
    python tools/benchmark.py --sizes 21,41,61 --days 300

Məqsəd: A6 (fully implicit) mərhələsindən əvvəl vaxtın harada getdiyini
sənədləşdirmək və hər optimallaşdırmadan sonra müqayisə etmək.
"""

from __future__ import annotations

import argparse
import cProfile
import io
import os
import pstats
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from imex2d.application.config import (LinearSolverConfig, OutputConfig,
                                       SimulationConfig)
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.application.simulation_service import SimulationService
from imex2d.domain.scal import CoreyParameters
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter


def build_case(nx: int, ny: int = None, nz: int = 1):
    ny = ny or nx
    scal = CoreyParameters()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=20.0, dy=20.0, dz=10.0, porosity=0.22,
        permx_base=150.0, nz=nz)
    model = ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                          scal=scal)
    return model, scal


def make_service(scal, refresh: int = None) -> SimulationService:
    config = LinearSolverConfig()
    if refresh:
        config.preconditioner_refresh_steps = refresh
    return SimulationService(CoreyRelativePermeabilityAdapter(scal),
                             linear_solver=ScipyCgIluSolver(config))


def run_case(nx: int, nz: int = 1, days: float = 300.0, refresh: int = None,
             engine: str = "impes"):
    model, scal = build_case(nx, nz=nz)
    config = SimulationConfig(end_time=days,
                              output=OutputConfig(snapshot_count=3))
    if refresh:
        config.linear_solver.preconditioner_refresh_steps = refresh
    service = make_service(scal, refresh)
    if engine == "implicit":
        service = service.with_engine(FullyImplicitEngine)
    started = time.perf_counter()
    result = service.run(model, config)
    elapsed = time.perf_counter() - started
    return {
        "cells": model.ncell,
        "steps": result.steps,
        "seconds": elapsed,
        "ms_per_step": elapsed / max(result.steps, 1) * 1000.0,
        "us_per_cell_step": elapsed / max(result.steps * model.ncell, 1) * 1e6,
        "rf": result.final_recovery_factor,
    }


def scaling_table(sizes, days: float, nz: int = 1,
                  engine: str = "impes") -> None:
    print(f"\n{engine.upper()}")
    print(f"{'grid':>10} {'hüceyrə':>9} {'addım':>7} {'saniyə':>8} "
          f"{'ms/addım':>10} {'µs/hüc·addım':>13}")
    print("-" * 62)
    for size in sizes:
        data = run_case(size, nz=nz, days=days, engine=engine)
        label = f"{size}x{size}" + (f"x{nz}" if nz > 1 else "")
        print(f"{label:>10} {data['cells']:>9} {data['steps']:>7} "
              f"{data['seconds']:>8.2f} {data['ms_per_step']:>10.2f} "
              f"{data['us_per_cell_step']:>13.2f}")


def profile_case(nx: int, nz: int, days: float, top: int = 15) -> None:
    model, scal = build_case(nx, nz=nz)
    config = SimulationConfig(end_time=days,
                              output=OutputConfig(snapshot_count=3))
    engine = make_service(scal).create_engine(model, config)

    profiler = cProfile.Profile()
    profiler.enable()
    engine.run()
    profiler.disable()

    stream = io.StringIO()
    pstats.Stats(profiler, stream=stream).sort_stats("tottime").print_stats(top)
    print(f"\nProfil — {nx}x{nx}x{nz} ({model.ncell} hüceyrə)")
    print("\n".join(stream.getvalue().splitlines()[4:top + 8]))


def main() -> None:
    parser = argparse.ArgumentParser(description="IMEX-2D performans ölçməsi")
    parser.add_argument("--sizes", default="21,31,41",
                        help="vergüllə ayrılmış grid ölçüləri")
    parser.add_argument("--days", type=float, default=300.0)
    parser.add_argument("--nz", type=int, default=1)
    parser.add_argument("--profile", type=int, default=0,
                        help="verilən ölçü üçün cProfile hesabatı")
    parser.add_argument("--engine", default="impes",
                        choices=("impes", "implicit", "both"))
    arguments = parser.parse_args()

    if arguments.profile:
        profile_case(arguments.profile, arguments.nz, arguments.days)
        return

    sizes = [int(value) for value in arguments.sizes.split(",")]
    engines = (("impes", "implicit") if arguments.engine == "both"
               else (arguments.engine,))
    for engine in engines:
        scaling_table(sizes, arguments.days, arguments.nz, engine)


if __name__ == "__main__":
    main()
