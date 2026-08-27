"""v69 addım 1b — imex2d/ui/panels.py faylından GOC (qaz-neft kontaktı) çıxarır.

Birinci patch-dən (patch_panels_gas.py) SONRA işlədilməlidir —
sətir nömrələri 904 sətirlik faylı nəzərdə tutur.

İstifadə:
    python tools/patch_panels_goc.py --dry
    python tools/patch_panels_goc.py

SAXLANILANLAR: owc (su-neft kontaktı), use_equilibration, datum_depth —
bunlar 2 fazalı neft-su modelində də lazımdır.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TARGET = ROOT / "imex2d" / "ui" / "panels.py"

EXPECTED_LEN = 904

CHECKS = [
    (857, "self.owc = _spin("),                    # qalır
    (858, "self.use_goc = QCheckBox("),
    (859, "(GOC)"),
    (860, "self.goc = _spin("),
    (861, "form.addRow(self.use_equilibration)"),  # qalır
    (872, "self.use_equilibration.stateChanged"),  # qalır
    (873, "form.addRow(self.use_goc)"),
    (874, "self.goc)"),
    (875, "self.use_goc.stateChanged"),
    (876, "self.goc.valueChanged"),
    (877, "note = QLabel("),                       # qalır
    (884, "equilibrate = self.use_equilibration.isChecked()"),
    (885, "use_goc = equilibrate and self.use_goc.isChecked()"),
    (886, "return InitialConditions("),
    (890, "oil_water_contact=self.owc.value()"),   # qalır
    (891, "gas_oil_contact=self.goc.value()"),
    (892, "use_equilibration=equilibrate)"),       # qalır
]

DELETIONS = [
    (858, 860, "use_goc checkbox + goc spin"),
    (873, 876, "GOC forma sətirləri + siqnallar"),
    (885, 885, "use_goc dəyişəni"),
    (891, 891, "gas_oil_contact arqumenti"),
]


def main():
    dry = "--dry" in sys.argv

    if not TARGET.exists():
        print(f"XƏTA: fayl tapılmadı: {TARGET}")
        return 2

    lines = TARGET.read_text(encoding="utf-8").splitlines(keepends=True)
    print(f"{TARGET.relative_to(ROOT)}: {len(lines)} sətir\n")

    if len(lines) != EXPECTED_LEN:
        print(f"XƏTA: {EXPECTED_LEN} sətir gözlənilirdi, {len(lines)} tapıldı.")
        print("Birinci patch (patch_panels_gas.py) tətbiq olunub?")
        print("HEÇ NƏ YAZILMADI.")
        return 1

    failed = 0
    for no, expected in CHECKS:
        actual = lines[no - 1].rstrip("\n")
        if expected not in actual:
            print(f"  XƏTA sətir {no}: gözlənilirdi {expected!r}")
            print(f"                  tapıldı      {actual!r}")
            failed += 1

    if failed:
        print(f"\n{failed} yoxlama uğursuz. HEÇ NƏ YAZILMADI.")
        return 1

    print(f"  {len(CHECKS)} yoxlamanın hamısı keçdi.\n")
    print("  Silinəcək aralıqlar:")

    total = 0
    for start, end, label in sorted(DELETIONS, reverse=True):
        count = end - start + 1
        total += count
        print(f"    {start}-{end} ({count} sətir) — {label}")
        del lines[start - 1:end]

    print(f"\n  Cəmi silinən: {total} sətir")
    print(f"  Yeni ölçü: {len(lines)} sətir")

    if dry:
        print("\n--dry rejimi: fayl DƏYİŞDİRİLMƏDİ.")
        return 0

    TARGET.write_text("".join(lines), encoding="utf-8")
    print(f"\nYazıldı: {TARGET.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
