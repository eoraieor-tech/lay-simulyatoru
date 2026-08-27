"""Qızıl etalon (golden files) — v68 konsolidasiya buraxılışı.

Məqsəd: qaz fazası və matplotlib 3D silinməzdən ƏVVƏL 2 fazalı nəticələri
qeyd etmək, silmədən SONRA isə eyni nəticələrin çıxdığını yoxlamaq.

Yalnız 2 FAZALI (neft-su) kəmiyyətlər qeyd olunur — qaza aid heç nə yoxdur,
ona görə də qaz silinəndən sonra da bu etalonlar keçərli qalır.

İstifadə:
    python tools/golden.py --write     # etalonu yaz (v67-də BİR DƏFƏ)
    python tools/golden.py             # cari kodu etalonla müqayisə et
    python tools/golden.py --case bl_1d --write   # tək keys

Nəticələr: tests/golden/<keys>.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

GOLDEN_DIR = ROOT / "tests" / "golden"

# Müqayisə dözümlülüyü. Təmiz refaktorinqdə nəticə praktiki olaraq
# eyni olmalıdır; bu hədd yalnız float toplama sırasının dəyişməsinə yer verir.
RTOL = 1e-9
ATOL = 1e-12


# ---------------------------------------------------------------- keyslər

def _cases():
    """Etalon keysləri. helpers.py-dakı qurucuları işlədir ki,
    mövcud testlərlə eyni yoldan getsin."""
    import helpers as h

    return {
        # 41x41 five-spot, qısa müddət — əsas sahə keysı
        "five_spot": lambda: (
            h.five_spot_model(),
            h.short_config(end_time=300.0, snapshots=10),
        ),
        # 1D Bakli-Leverett — analitik həllə yaxın, kəskin cəbhə
        "bl_1d": lambda: (
            h.one_dimensional_model(),
            h.bl_config(end_time=250.0),
        ),
        # Kiçik grid, uzun müddət — fərqli zaman addımı rejimi
        "five_spot_small": lambda: (
            h.five_spot_model(nx=21, ny=21, dx=40.0, dy=40.0),
            h.short_config(end_time=600.0, snapshots=6),
        ),
    }


# ---------------------------------------------------------------- çıxarış

def _digest(values, n_samples=7):
    """Bir seriyanın barmaq izi: ölçü, uc nöqtələr, statistika, nümunələr."""
    arr = np.asarray(values, dtype=float).ravel()
    if arr.size == 0:
        return {"n": 0}
    idx = np.unique(
        np.linspace(0, arr.size - 1, min(n_samples, arr.size)).astype(int)
    )
    return {
        "n": int(arr.size),
        "first": float(arr[0]),
        "last": float(arr[-1]),
        "min": float(arr.min()),
        "max": float(arr.max()),
        "mean": float(arr.mean()),
        "sum": float(arr.sum()),
        "samples": [float(arr[i]) for i in idx],
    }


# Yalnız 2 fazalı seriyalar — gas_rate, cumulative_gas, gas_oil_ratio YOXDUR
SERIES_FIELDS = (
    "time",
    "oil_rate",
    "water_rate",
    "water_injection_rate",
    "cumulative_oil",
    "cumulative_water",
    "water_cut",
    "average_pressure",
    "recovery_factor",
)


def extract(result):
    """SimulationResult-dan müqayisə oluna bilən quruluş çıxarır."""
    out = {
        "grid_shape": list(getattr(result, "grid_shape", ()) or ()),
        "steps": int(getattr(result, "steps", 0)),
        "converged": bool(getattr(result, "converged", True)),
        "ooip": float(getattr(result, "ooip", 0.0)),
    }

    series = getattr(result, "series", None)
    out["series"] = {
        name: _digest(getattr(series, name, []) or [])
        for name in SERIES_FIELDS
    }

    snaps = list(getattr(result, "snapshots", []) or [])
    out["snapshot_count"] = len(snaps)
    out["snapshots"] = [
        {
            "time": float(s.time),
            "pressure": _digest(s.pressure),
            "water_saturation": _digest(s.water_saturation),
        }
        for s in snaps
    ]

    # Quyu debitləri — yalnız neft və su (qaz yox)
    out["wells"] = {}
    for field in ("well_oil_rate", "well_water_rate"):
        table = getattr(result, field, {}) or {}
        out["wells"][field] = {
            str(k): _digest(v) for k, v in sorted(table.items())
        }

    return out


# ---------------------------------------------------------------- müqayisə

def compare(expected, actual, path=""):
    """Rekursiv müqayisə. Uyğunsuzluqların siyahısını qaytarır."""
    diffs = []

    if isinstance(expected, dict):
        if not isinstance(actual, dict):
            return [f"{path}: tip fərqi (dict gözlənilirdi)"]
        for key in expected:
            if key not in actual:
                diffs.append(f"{path}.{key}: yoxa çıxıb")
            else:
                diffs += compare(expected[key], actual[key], f"{path}.{key}")
        for key in actual:
            if key not in expected:
                diffs.append(f"{path}.{key}: yeni sahə")
        return diffs

    if isinstance(expected, list):
        if not isinstance(actual, list):
            return [f"{path}: tip fərqi (list gözlənilirdi)"]
        if len(expected) != len(actual):
            return [f"{path}: uzunluq {len(expected)} -> {len(actual)}"]
        for i, (e, a) in enumerate(zip(expected, actual)):
            diffs += compare(e, a, f"{path}[{i}]")
        return diffs

    if isinstance(expected, bool) or isinstance(actual, bool):
        if expected != actual:
            diffs.append(f"{path}: {expected} -> {actual}")
        return diffs

    if isinstance(expected, (int, float)) and isinstance(actual, (int, float)):
        if not np.isclose(expected, actual, rtol=RTOL, atol=ATOL,
                          equal_nan=True):
            delta = abs(actual - expected)
            rel = delta / max(abs(expected), 1e-30)
            diffs.append(
                f"{path}: {expected!r} -> {actual!r} "
                f"(mütləq {delta:.3e}, nisbi {rel:.3e})"
            )
        return diffs

    if expected != actual:
        diffs.append(f"{path}: {expected!r} -> {actual!r}")
    return diffs


# ---------------------------------------------------------------- işlətmə

def run_case(name, factory):
    import helpers as h

    model, config = factory()
    service = h.make_service()
    t0 = time.perf_counter()
    result = service.run(model, config)
    elapsed = time.perf_counter() - t0
    return extract(result), elapsed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true",
                        help="etalonu yaz (mövcudun üstünə yazır)")
    parser.add_argument("--case", default=None, help="yalnız bir keys")
    args = parser.parse_args()

    cases = _cases()
    if args.case:
        if args.case not in cases:
            print(f"Naməlum keys: {args.case}")
            print(f"Mövcud: {', '.join(cases)}")
            return 2
        cases = {args.case: cases[args.case]}

    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    failed = 0

    for name, factory in cases.items():
        print(f"\n=== {name} ===")
        try:
            data, elapsed = run_case(name, factory)
        except Exception as exc:
            print(f"  XƏTA: {type(exc).__name__}: {exc}")
            failed += 1
            continue

        print(f"  addım: {data['steps']}, yığıldı: {data['converged']}, "
              f"vaxt: {elapsed:.1f} s")
        rf = data["series"]["recovery_factor"]
        if rf.get("n"):
            print(f"  son RF: {rf['last']:.6f}")

        path = GOLDEN_DIR / f"{name}.json"

        if args.write:
            path.write_text(
                json.dumps(data, indent=2, sort_keys=True, ensure_ascii=False),
                encoding="utf-8",
            )
            print(f"  yazıldı -> {path.relative_to(ROOT)}")
            continue

        if not path.exists():
            print(f"  ETALON YOXDUR: {path.relative_to(ROOT)}")
            print("  əvvəlcə --write ilə yarat")
            failed += 1
            continue

        expected = json.loads(path.read_text(encoding="utf-8"))
        diffs = compare(expected, data)
        if diffs:
            print(f"  UYĞUNSUZLUQ ({len(diffs)} ədəd):")
            for d in diffs[:25]:
                print(f"    {d}")
            if len(diffs) > 25:
                print(f"    ... və daha {len(diffs) - 25}")
            failed += 1
        else:
            print("  UYĞUNDUR")

    print()
    if args.write:
        print("Etalonlar yazıldı. İndi commit et.")
        return 0
    if failed:
        print(f"{failed} keys uyğun gəlmədi.")
        return 1
    print("Bütün keyslər etalona uyğundur.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
