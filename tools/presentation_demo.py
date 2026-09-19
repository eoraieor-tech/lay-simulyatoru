"""Təqdimat üçün nümunə qaçışı — slayd və səhifələrdəki BÜTÜN rəqəmlər (Seans 44).

    python tools/presentation_demo.py --out <qovluq>

Üç iş görür və nəticəni `<qovluq>`-a yazır:

1. `numune_quyular.csv` → adi Kriging → 41 × 41 × 3 geoloji model →
   five-spot (INJ-1 320 bar / PROD-1 150 bar) → `FullyImplicitEngine`,
   TPFA, 1500 gün. Hər Nyuton cəhdinin CNV tarixçəsi də yazılır.
2. 1D Buckley-Leverett (`tests/helpers.one_dimensional_model`, IMPES,
   250 gün) — metrikalar `MainWindow.run_validation` ilə EYNİDİR:
   RMS, orta nöqtə cəbhəsi, həcm balansı.
3. PNG qrafiklər (`--no-figs` ilə buraxılır).

Çıxış: `demo.json`, `bl_metrics.json`, `figs/*.png`.
Seans 44-də ölçülən dəyərlər `docs/teqdimat/README.md`-dədir.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

import numpy as np

from helpers import (bl_config, default_scal, make_service,
                     one_dimensional_model)
from imex2d.application.config import OutputConfig, SimulationConfig
from imex2d.application.geology_service import (GeologicalGridSpec,
                                                WellBasedGeologicalModelBuilder)
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import five_spot
from imex2d.geology.interpolation import OrdinaryKriging
from imex2d.geology.well_data_io import read_well_csv
from imex2d.simulation.analytical import buckley_leverett
from imex2d.simulation.implicit.engine import FullyImplicitEngine

BL_RATE, BL_TIME, BL_POROSITY, BL_DY, BL_DZ, BL_NX, BL_DX = 60.0, 250.0, 0.20, 100.0, 10.0, 120, 8.0


def run_example():
    """Nümunə yataq: geologiya + simulyasiya. `(data, massivlər)` qaytarır."""
    started = time.time()
    dataset = read_well_csv(os.path.join(ROOT, "numune_quyular.csv"))
    spec = GeologicalGridSpec(nx=41, ny=41, nz=3, dx=20.0, dy=20.0, dz=5.0,
                              top_depth=2000.0)
    geology, _ = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(dataset, spec)
    geology_seconds = time.time() - started
    shape = geology.grid.shape
    permx = geology.property_maps["PERMX"].values.reshape(shape)
    poro = geology.property_maps["PORO"].values.reshape(shape)

    model = ReservoirModelBuilder().build(
        geological_model=geology, wells=five_spot(geology.grid),
        scal=default_scal(), name="Nümunə five-spot")
    service = make_service().with_engine(FullyImplicitEngine)

    # Hər Nyuton CƏHDİNİ qeyd et (rədd olunanlar da daxil)
    attempts, holder = [], {}
    create = service.create_engine

    def create_logged(m, c):
        engine = create(m, c)
        solve = engine.newton.solve

        def logged(state, dt, *args, **kwargs):
            result = solve(state, dt, *args, **kwargs)
            attempts.append(dict(dt=float(dt), it=int(result.iterations),
                                 ok=bool(result.converged),
                                 cnv=[float(v) for v in result.cnv_history]))
            return result
        engine.newton.solve = logged
        holder["engine"] = engine
        return engine
    service.create_engine = create_logged

    config = SimulationConfig(end_time=1500.0, output=OutputConfig(snapshot_count=6))
    started = time.time()
    result = service.run(model, config)
    run_seconds = time.time() - started
    history = holder["engine"].time_stepper.history
    series = result.series
    data = dict(
        geology=dict(shape=list(shape), seconds=round(geology_seconds, 2),
                     permx_min=float(permx.min()), permx_max=float(permx.max()),
                     poro_min=float(poro.min()), poro_max=float(poro.max())),
        model=dict(ncell=int(model.ncell), issues=model.validate(),
                   ooip=float(result.ooip)),
        run=dict(seconds=round(run_seconds, 1), steps=int(result.steps),
                 converged=bool(result.converged), message=result.message,
                 rf_final=float(result.final_recovery_factor),
                 breakthrough=result.breakthrough_time,
                 newton_total=int(sum(h.iterations for h in history)),
                 repeats=int(sum(h.repeats for h in history))),
        series=dict(time=list(map(float, series.time)),
                    oil_rate=list(map(float, series.oil_rate)),
                    water_rate=list(map(float, series.water_rate)),
                    water_cut=list(map(float, series.water_cut)),
                    rf=list(map(float, series.recovery_factor))),
        steps=[dict(t=float(h.time), dt=float(h.dt), it=int(h.iterations),
                    rep=int(h.repeats), ds=float(h.max_saturation_change))
               for h in history],
        newton=attempts,
        snapshots=[float(s.time) for s in result.snapshots])
    arrays = dict(permx=permx,
                  sw=np.stack([s.water_saturation.reshape(shape) for s in result.snapshots]))
    return data, arrays


def run_buckley_leverett():
    """1D doğrulama — metrikalar `MainWindow.run_validation` ilə eyni."""
    scal = default_scal()
    model = one_dimensional_model()
    result = make_service().run(model, bl_config(BL_TIME))
    snapshot = result.snapshots[-1]
    sw = snapshot.water_saturation.ravel()
    x = (np.arange(BL_NX) + 0.5) * BL_DX
    solution = buckley_leverett(scal, model.fluids.water_viscosity,
                                model.fluids.oil_viscosity, BL_POROSITY,
                                BL_RATE, BL_DY * BL_DZ, snapshot.time,
                                length=BL_NX * BL_DX)
    exact = np.interp(x, solution.distance, solution.water_saturation)
    rms = float(np.sqrt(np.mean((sw - exact) ** 2)))
    level = 0.5 * (solution.shock_saturation + scal.swc)
    last = int(np.nonzero(sw >= level)[0][-1])
    weight = (sw[last] - level) / (sw[last] - sw[last + 1])
    front = float(x[last] + weight * (x[last + 1] - x[last]))
    swept = float(np.trapezoid(sw - scal.swc, x))
    injected = BL_RATE * snapshot.time / (BL_POROSITY * BL_DY * BL_DZ)
    return dict(time=float(snapshot.time), rms=rms, front_numeric=front,
                front_analytic=float(solution.front_position),
                front_error_pct=abs(front - solution.front_position)
                / solution.front_position * 100.0,
                balance_error_pct=abs(swept - injected) / injected * 100.0,
                shock_sw=float(solution.shock_saturation),
                x=x.tolist(), sw_numeric=sw.tolist(),
                x_analytic=list(map(float, solution.distance)),
                sw_analytic=list(map(float, solution.water_saturation)))


def ordinal_suffix(n: int) -> str:
    """Azərbaycan dilində sıra şəkilçisi: 257 → ci, 1500 → cü."""
    units = {1: "ci", 2: "ci", 3: "cü", 4: "cü", 5: "ci", 6: "cı", 7: "ci", 8: "ci", 9: "cu"}
    tens = {1: "cu", 2: "ci", 3: "cu", 4: "cı", 5: "ci", 6: "cı", 7: "ci", 8: "ci", 9: "cı"}
    if n % 10:
        return units[n % 10]
    if n % 100:
        return tens[(n // 10) % 10]
    return "ci" if n % 1000 == 0 else "cü"


def write_figures(data, arrays, bl, folder):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, LogNorm

    ink, muted, grid = "#13212B", "#4B5B66", "#C9D2D4"
    oil, water = "#C8742A", "#2A7FA0"
    plt.rcParams.update({
        "font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 15,
        "text.color": ink, "axes.labelcolor": ink, "axes.edgecolor": grid,
        "xtick.color": muted, "ytick.color": muted, "axes.titlesize": 17,
        "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False})
    os.makedirs(folder, exist_ok=True)

    def save(fig, name):
        fig.savefig(os.path.join(folder, name), dpi=150, transparent=True,
                    bbox_inches="tight", pad_inches=0.15)
        plt.close(fig)

    permx, sw = arrays["permx"], arrays["sw"]
    nz, ny, nx = permx.shape
    size = 20.0
    wells = {"W-1": (123, 123), "W-2": (656, 164), "W-3": (205, 615),
             "W-4": (697, 697), "W-5": (410, 410)}
    amber = LinearSegmentedColormap.from_list("amber", ["#FBF1E3", "#E9B77A", "#C8742A", "#6E3A12"])
    fig, axes = plt.subplots(1, nz, figsize=(15, 5.4), constrained_layout=True)
    norm = LogNorm(vmin=permx.min(), vmax=permx.max())
    for k, ax in enumerate(axes):
        image = ax.imshow(permx[k], origin="lower", extent=[0, nx * size, 0, ny * size], cmap=amber, norm=norm)
        for name, (wx, wy) in wells.items():
            ax.plot(wx, wy, "o", ms=9, mfc="white", mec=ink, mew=1.8)
            ax.annotate(name, (wx, wy), xytext=(7, 7), textcoords="offset points",
                        fontsize=12, weight="bold")
        ax.set_title(f"Qat {k + 1}")
        ax.set_xlabel("x, m")
    axes[0].set_ylabel("y, m")
    fig.colorbar(image, ax=axes, shrink=0.9, pad=0.01).set_label("PERMX, mD")
    save(fig, "permx.png")

    blue = LinearSegmentedColormap.from_list("water", ["#F6EBDD", "#D9C3A5", "#7FB5C9", "#2A7FA0", "#12455A"])
    times = data["snapshots"]
    fig, axes = plt.subplots(1, 4, figsize=(18, 5), constrained_layout=True)
    for ax, i in zip(axes, [1, 2, 4, 6]):
        image = ax.imshow(sw[i][0], origin="lower", extent=[0, nx * size, 0, ny * size],
                          cmap=blue, vmin=0.2, vmax=0.8)
        ax.plot(10, 10, "v", ms=13, mfc=water, mec="white", mew=1.5)
        ax.plot(nx * size - 10, ny * size - 10, "^", ms=13, mfc=oil, mec="white", mew=1.5)
        day = round(times[i])
        ax.set_title(f"{day}-{ordinal_suffix(day)} gün")
        ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(image, ax=axes, shrink=0.85, pad=0.01).set_label("Sw (su doyumluluğu)")
    save(fig, "sw_front.png")

    s = data["series"]
    t = np.array(s["time"])
    fig, (left, right) = plt.subplots(1, 2, figsize=(16, 5.6), constrained_layout=True)
    left.plot(t, s["oil_rate"], color=oil, lw=2.6, label="Neft debiti")
    left.plot(t, s["water_rate"], color=water, lw=2.6, label="Su debiti")
    left.set_xlabel("Zaman, gün"); left.set_ylabel("Debit, m³/gün"); left.legend(frameon=False)
    right.plot(t, s["rf"], color=ink, lw=2.8, label="Neftvermə (RF), %")
    right.plot(t, s["water_cut"], color=water, lw=2.2, label="Sulaşma, %")
    right.set_xlabel("Zaman, gün"); right.set_ylim(0, 100); right.legend(frameon=False)
    save(fig, "production.png")

    fig, ax = plt.subplots(figsize=(13, 5.6), constrained_layout=True)
    ax.plot(bl["x_analytic"], bl["sw_analytic"], color=ink, lw=2.6, label="Analitik həll (Buckley-Leverett)")
    ax.plot(bl["x"], bl["sw_numeric"], "o", ms=4.5, color=water, label="IMEX-2D (ədədi)")
    ax.set_xlim(0, 400); ax.set_xlabel("Vurucudan məsafə, m"); ax.set_ylabel("Sw")
    ax.legend(frameon=False)
    save(fig, "bl.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", required=True, help="nəticə qovluğu")
    parser.add_argument("--no-figs", action="store_true")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    data, arrays = run_example()
    bl = run_buckley_leverett()
    with open(os.path.join(args.out, "demo.json"), "w", encoding="utf-8") as handle:
        json.dump(data, handle, ensure_ascii=False, indent=1, default=str)
    with open(os.path.join(args.out, "bl_metrics.json"), "w", encoding="utf-8") as handle:
        json.dump({k: v for k, v in bl.items() if not isinstance(v, list)}, handle, indent=1)
    if not args.no_figs:
        write_figures(data, arrays, bl, os.path.join(args.out, "figs"))

    run = data["run"]
    print(f"Nümunə: {run['message']}")
    print(f"  RF = {run['rf_final']:.2f} %, su çatdı {run['breakthrough']:.0f}-ci gün, "
          f"{run['seconds']} s, OOIP {data['model']['ooip']:.0f} m³, "
          f"{len(data['newton'])} Nyuton cəhdi")
    print(f"BL: RMS {bl['rms']:.4f}, cəbhə {bl['front_numeric']:.1f} / "
          f"{bl['front_analytic']:.1f} m ({bl['front_error_pct']:.2f} %), "
          f"həcm balansı {bl['balance_error_pct']:.2f} %")


if __name__ == "__main__":
    main()
