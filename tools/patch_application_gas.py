"""v69 addım 2 — application qatından qaz kodunu çıxarır.

Fayllar:
  imex2d/application/model_builder.py      gas_scal parametri
  imex2d/application/serialization.py      .imx-də gas_oil_contact açarı
  imex2d/application/simulation_service.py 3 fazalı mühərriyə keçid budağı

Bu addım `three_phase_engine.py` və `stone_relperm.py`-a olan YEGANƏ
işlək istinadı qoparır — ondan sonra o fayllar ölü koda çevrilir.

İstifadə:
    python tools/patch_application_gas.py --dry
    python tools/patch_application_gas.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
APP = ROOT / "imex2d" / "application"

# fayl -> (yoxlamalar, silinəcək aralıqlar)
PLAN = {
    "model_builder.py": {
        "checks": [
            (34, "scal: Optional[CoreyParameters] = None,"),
            (35, "gas_scal=None,"),
            (36, "capillary: Optional[CapillaryParameters] = None,"),
            (75, "scal_parameters=scal or CoreyParameters(),"),
            (76, "gas_scal_parameters=gas_scal,"),
            (77, "capillary_parameters=capillary or CapillaryParameters(),"),
        ],
        "deletions": [
            (76, 76, "gas_scal_parameters ötürməsi"),
            (35, 35, "gas_scal parametri"),
        ],
        "rewrites": [],
    },
    "serialization.py": {
        "checks": [
            (210, '"datum_depth", "datum_pressure", "water_saturation",'),
            (211, '"gas_oil_contact"'),
            (212, '"use_equilibration"]),'),
        ],
        "deletions": [],
        "rewrites": [
            (211, '"gas_oil_contact", ', "", "gas_oil_contact açarı"),
        ],
    },
    "simulation_service.py": {
        "checks": [
            (133, "def create_engine(self, model, config):"),
            (141, "if model.initial_conditions.use_equilibration else None)"),
            (143, "# A7"),
            (151, "self.pvt_provider.has_gas_phase()"),
            (153, "from ..simulation.implicit.three_phase_engine import ("),
            (155, 'gas_scal = getattr(model, "gas_scal_parameters", None)'),
            (159, "self.engine_factory = ThreePhaseSimulationEngine"),
            (160, "return super().create_engine(model, config)"),
        ],
        "deletions": [
            (143, 159, "A7 şərhləri + 3 fazalı mühərrik budağı"),
        ],
        "rewrites": [],
    },
}


def process(name, spec, dry):
    path = APP / name
    print(f"--- {name}")

    if not path.exists():
        print("  XƏTA: fayl tapılmadı")
        return None

    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    print(f"  {len(lines)} sətir")

    failed = 0
    for no, expected in spec["checks"]:
        if no > len(lines):
            print(f"  XƏTA sətir {no}: fayl qısadır")
            failed += 1
            continue
        actual = lines[no - 1].rstrip("\n")
        if expected not in actual:
            print(f"  XƏTA sətir {no}: gözlənilirdi {expected!r}")
            print(f"                  tapıldı      {actual!r}")
            failed += 1

    if failed:
        return None

    print(f"  {len(spec['checks'])} yoxlama keçdi")

    for no, old, new, label in spec["rewrites"]:
        cur = lines[no - 1]
        if old not in cur:
            print(f"  XƏTA sətir {no}: {old!r} tapılmadı")
            return None
        upd = cur.replace(old, new)
        print(f"  yenidən yazılır {no} — {label}")
        print(f"    {cur.rstrip()}")
        print(f" -> {upd.rstrip()}")
        lines[no - 1] = upd

    total = 0
    for start, end, label in sorted(spec["deletions"], reverse=True):
        count = end - start + 1
        total += count
        print(f"  silinir {start}-{end} ({count} sətir) — {label}")
        del lines[start - 1:end]

    print(f"  yeni ölçü: {len(lines)} sətir (-{total})")
    return path, "".join(lines)


def main():
    dry = "--dry" in sys.argv
    results = []

    for name, spec in PLAN.items():
        out = process(name, spec, dry)
        if out is None:
            print("\nBir fayl uğursuz oldu. HEÇ NƏ YAZILMADI.")
            return 1
        results.append(out)
        print()

    if dry:
        print("--dry rejimi: fayllar DƏYİŞDİRİLMƏDİ.")
        return 0

    for path, content in results:
        path.write_text(content, encoding="utf-8")
        print(f"Yazıldı: {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
