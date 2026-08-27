"""v69 addım 1 — imex2d/ui/panels.py faylından qaz kodunu çıxarır.

TƏHLÜKƏSİZLİK: hər dəyişiklikdən əvvəl həmin sətrin gözlənilən mətnə
uyğunluğu yoxlanılır. Bir yoxlama belə uğursuz olsa, HEÇ NƏ yazılmır.

İstifadə:
    python tools/patch_panels_gas.py --dry     # yalnız yoxla, yazma
    python tools/patch_panels_gas.py           # tətbiq et

SAXLANILANLAR (bilərəkdən):
  - self.gas_gravity  -> qaz sıxlığı qara-neft korrelyasiyalarında neftin
                         öz xassələrinə (Rs, Bo, mu_o) daxildir, 2 fazalı
                         modeldə də lazımdır
  - self.bubble_point -> doymuş/doymamış ayrımı 2 fazalı black-oil-da da var
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TARGET = ROOT / "imex2d" / "ui" / "panels.py"


# (sətir_nömrəsi, gözlənilən_alt_sətir) — 1-dən başlayır
CHECKS = [
    (37, "GasCoreyParameters)"),
    (598, "Qaz-neft SCAL"),
    (599, "self.gas_enabled = QCheckBox("),
    (627, "self._on_gas_toggled()"),
    (629, "def _on_gas_toggled(self):"),
    (632, "widget.setEnabled(active)"),
    (634, "def values(self)"),
    (643, ""),                       # boş sətir
    (644, "def gas_values(self)"),
    (655, "self.nog.value())"),
    (688, "self.points"),
    (690, "self.gas_phase_enabled = QCheckBox("),
    (693, "self.gas_phase_enabled.stateChanged.connect"),
    (712, "form.addRow(note)"),
    (714, "form.addRow(self.gas_phase_enabled)"),
    (722, "form.addRow(gas_note)"),
    (724, "def is_enabled(self)"),
    (741, "bubble_point_bar=self.bubble_point.value(),"),
    (742, "include_gas=self.gas_phase_enabled.isChecked())"),
    (744, "def gas_phase_active(self)"),
    (745, "return self.enabled.isChecked()"),
    (748, "class WellPanel(QWidget):"),
]

# Silinəcək aralıqlar (1-dən başlayır, hər iki ucu daxil)
DELETIONS = [
    (36, 37, "idxal: GasCoreyParameters"),      # 36-nın vergülü ')' olur
    (598, 633, "ScalPanel qaz bloku + _on_gas_toggled"),
    (643, 655, "ScalPanel.gas_values"),
    (689, 693, "PvtPanel.gas_phase_enabled checkbox"),
    (713, 722, "PvtPanel qaz izahı + addRow"),
    (741, 742, "include_gas arqumenti"),        # 741 yenidən yazılır
    (743, 745, "PvtPanel.gas_phase_active"),
]


def main():
    dry = "--dry" in sys.argv

    if not TARGET.exists():
        print(f"XƏTA: fayl tapılmadı: {TARGET}")
        return 2

    lines = TARGET.read_text(encoding="utf-8").splitlines(keepends=True)
    print(f"{TARGET.relative_to(ROOT)}: {len(lines)} sətir\n")

    # ---- yoxlamalar
    failed = 0
    for no, expected in CHECKS:
        if no > len(lines):
            print(f"  XƏTA sətir {no}: fayl qısadır")
            failed += 1
            continue
        actual = lines[no - 1].rstrip("\n")
        if expected == "":
            ok = actual.strip() == ""
        else:
            ok = expected in actual
        if not ok:
            print(f"  XƏTA sətir {no}: gözlənilirdi {expected!r}")
            print(f"                  tapıldı      {actual!r}")
            failed += 1

    if failed:
        print(f"\n{failed} yoxlama uğursuz. HEÇ NƏ YAZILMADI.")
        print("Fayl gözlədiyimdən fərqlidir — mənə bu çıxışı göndər.")
        return 1

    print(f"  {len(CHECKS)} yoxlamanın hamısı keçdi.\n")

    # ---- iki xüsusi hal: vergülü bağlayan mötərizəyə çevir
    # 36: "..., " ilə bitir -> ")" ilə bitməlidir
    l36 = lines[35].rstrip("\n").rstrip()
    if not l36.endswith(","):
        print(f"  XƏTA sətir 36 vergüllə bitmir: {l36!r}")
        return 1
    new36 = l36[:-1] + ")\n"

    # 741: "...value()," -> "...value())"
    l741 = lines[740].rstrip("\n").rstrip()
    if not l741.endswith(","):
        print(f"  XƏTA sətir 741 vergüllə bitmir: {l741!r}")
        return 1
    new741 = l741[:-1] + ")\n"

    print("  Dəyişəcək sətirlər:")
    print(f"    36:  {l36}")
    print(f"      -> {new36.rstrip()}")
    print(f"    741: {l741}")
    print(f"      -> {new741.rstrip()}\n")

    lines[35] = new36
    lines[740] = new741

    # ---- silmələr: AŞAĞIDAN YUXARI (sətir nömrələri sürüşməsin)
    print("  Silinəcək aralıqlar:")
    total = 0
    for start, end, label in sorted(DELETIONS, reverse=True):
        # 36 və 741 yenidən yazıldı -> onları saxla, sonrakını sil
        s = start + 1 if start in (36, 741) else start
        count = end - s + 1
        total += count
        print(f"    {s}-{end} ({count} sətir) — {label}")
        del lines[s - 1:end]

    print(f"\n  Cəmi silinən: {total} sətir")
    print(f"  Yeni ölçü: {len(lines)} sətir")

    if dry:
        print("\n--dry rejimi: fayl DƏYİŞDİRİLMƏDİ.")
        return 0

    TARGET.write_text("".join(lines), encoding="utf-8")
    print(f"\nYazıldı: {TARGET.relative_to(ROOT)}")
    print("İndi yoxla:")
    print("  python -c \"import ast;ast.parse(open('imex2d/ui/panels.py',"
          "encoding='utf-8').read())\"")
    print("  findstr /n /i /c:\"gas\" imex2d\\ui\\panels.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
