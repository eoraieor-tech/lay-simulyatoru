"""SPE1CASE2: modeli qurur, 3650 gün işlədir, OPM Flow etalonu ilə müqayisə edir.

    python tools/spe1_compare.py <SPE1CASE2.DATA> <etalon_baza> [--max-dt 31]

`<etalon_baza>` — `SPE1CASE2.SMSPEC`/`.UNSMRY` fayllarının uzantısız yolu.
Faylların endirilməsi: `TEHVIL_TESLIM.md` §5 və `SPE1.md` §1.
"""

import argparse
import io
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from imex2d.application.simulation_service import ModelAwareSimulationService  # noqa: E402
from imex2d.benchmarks.spe1 import (build_spe1case2_model, compare_with_reference,  # noqa: E402
                                    first_time_below, gas_front_arrivals,
                                    simulated_series, spe1case2_config)
from imex2d.io.eclipse_summary import read_summary  # noqa: E402
from imex2d.simulation.implicit.engine import FullyImplicitEngine  # noqa: E402
from imex2d.simulation.linear_solver import ScipyCgIluSolver  # noqa: E402
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter  # noqa: E402


def main(argv=None):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")
    parser = argparse.ArgumentParser()
    parser.add_argument("deck")
    parser.add_argument("reference")
    parser.add_argument("--max-dt", type=float, default=31.0)
    args = parser.parse_args(argv)

    model = build_spe1case2_model(args.deck)
    service = ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(model.scal_parameters),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)
    started = time.time()
    result = service.run(model, spe1case2_config(max_dt=args.max_dt))
    print(f"{result.message} ({time.time() - started:.1f} san)")

    reference = read_summary(args.reference)
    print(f"{'t, gün':>7} {'kəmiyyət':<10} {'vahid':<9} {'bizdə':>10} "
          f"{'OPM Flow':>10} {'fərq, %':>8}")
    for row in compare_with_reference(model, result, reference):
        print(f"{row.time:7.0f} {row.quantity:<10} {row.unit:<9} "
              f"{row.simulated:10.4g} {row.reference:10.4g} "
              f"{100.0 * row.relative_error:8.1f}")

    time_, oil, _ = simulated_series(model, result)["FOPR"]
    # 0.5 % tolerans: səth debiti tənzimləyicisi ilk addımlarda hədəfdən
    # ~0.1 % sapır (ölçüldü: 19 981 STB/gün) — bu, rejim keçidi deyil.
    print("FOPR ilk dəfə 19 900-dən aşağı: bizdə",
          first_time_below(time_, oil, 19900.0), "gün, OPM",
          first_time_below(reference.series("TIME"), reference.series("FOPR"),
                           19900.0), "gün")

    print(f"\n{'blok':>5} {'(i,j,k)':<10} {'qaz çatıb, gün':>15} "
          f"{'OPM Flow':>10} {'fərq':>8}")
    for row in gas_front_arrivals(model, result, reference):
        i, j, k = row.ijk
        ours = "—" if row.simulated is None else f"{row.simulated:.0f}"
        opm = "—" if row.reference is None else f"{row.reference:.0f}"
        if row.simulated is None or row.reference is None:
            gap = "—"
        else:
            gap = f"{row.simulated - row.reference:+.0f}"
        print(f"{row.block:5d} {f'({i + 1},{j + 1},{k + 1})':<10} "
              f"{ours:>15} {opm:>10} {gap:>8}")


if __name__ == "__main__":
    main()
