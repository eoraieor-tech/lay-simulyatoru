"""v69 addım 2 — imex2d/ui/main_window.py faylından qaz kodunu çıxarır.

patch_panels_gas.py və patch_panels_goc.py-dan SONRA işlədilməlidir.

İstifadə:
    python tools/patch_mainwindow_gas.py --dry
    python tools/patch_mainwindow_gas.py

Dəyişikliklər:
  316   SCAL qrafiki 3 oxdan 2 oxa düşür (_figure(1,3) -> _figure(1,2))
  318   scal_gas_ax silinir
  1290-1298  rebuild_model: gas_active / gas_scal məntiqi silinir
  1304  model_builder.build(gas_scal=...) arqumenti silinir
  1565-1567  scal_renderer.draw_gas(...) çağırışı silinir
  1965-1966  Snapshot(gas_saturation=...) arqumenti silinir
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TARGET = ROOT / "imex2d" / "ui" / "main_window.py"

CHECKS = [
    (316, "_figure(1, 3)"),
    (317, "self.scal_axes = scal_axes[:2]"),
    (318, "self.scal_gas_ax = scal_axes[2]"),
    (319, "layout.addWidget(self.scal_canvas)"),
    (1289, "self.project.add_geological_model(geology)"),
    (1290, "gas_active = self.pvt_panel.gas_phase_active()"),
    (1291, "gas_scal = self.scal_panel.gas_values()"),
    (1297, "from ..domain.scal import GasCoreyParameters"),
    (1298, "gas_scal = GasCoreyParameters()"),
    (1299, "model = self.model_builder.build("),
    (1303, "scal=self.scal_panel.values(),"),
    (1304, "gas_scal=gas_scal,"),
    (1305, "capillary=self.scal_panel.capillary_values(),"),
    (1564, "capillary)"),
    (1565, "self.scal_renderer.draw_gas(self.scal_gas_ax,"),
    (1567, "self.scal_panel.gas_values())"),
    (1568, "self.scal_canvas.draw_idle()"),
    (1964, "snapshots=[Snapshot(time=s.time"),
    (1965, "water_saturation=s.water_saturation,"),
    (1966, "gas_saturation=s.gas_saturation)"),
    (1967, "for s in case.snapshots])"),
]

DELETIONS = [
    (1966, 1966, "Snapshot gas_saturation arqumenti"),
    (1565, 1567, "scal_renderer.draw_gas çağırışı"),
    (1304, 1304, "model_builder.build(gas_scal=...)"),
    (1290, 1298, "rebuild_model qaz məntiqi"),
    (318, 318, "scal_gas_ax"),
]


def main():
    dry = "--dry" in sys.argv

    if not TARGET.exists():
        print(f"XƏTA: fayl tapılmadı: {TARGET}")
        return 2

    lines = TARGET.read_text(encoding="utf-8").splitlines(keepends=True)
    print(f"{TARGET.relative_to(ROOT)}: {len(lines)} sətir\n")

    failed = 0
    for no, expected in CHECKS:
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
        print(f"\n{failed} yoxlama uğursuz. HEÇ NƏ YAZILMADI.")
        return 1

    print(f"  {len(CHECKS)} yoxlamanın hamısı keçdi.\n")

    # ---- yenidən yazılan sətirlər
    l316 = lines[315]
    new316 = l316.replace("_figure(1, 3)", "_figure(1, 2)")

    l1965 = lines[1964].rstrip("\n").rstrip()
    if not l1965.endswith(","):
        print(f"  XƏTA sətir 1965 vergüllə bitmir: {l1965!r}")
        return 1
    new1965 = l1965[:-1] + ")\n"

    print("  Dəyişəcək sətirlər:")
    print(f"    316:  {l316.rstrip()}")
    print(f"      -> {new316.rstrip()}")
    print(f"    1965: {l1965}")
    print(f"      -> {new1965.rstrip()}\n")

    lines[315] = new316
    lines[1964] = new1965

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
