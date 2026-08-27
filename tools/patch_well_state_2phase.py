"""v69 addım 4a — well_state.py-ı ThreePhaseState-dən ReservoirState-ə köçürür.

MƏQSƏD: quyu modelini (v61-66 işi) qaz kodundan AYIRMAQ, silmək yox.
Bundan sonra three_phase_state.py azad olur və bütün qaz faylları silinə bilər.

NƏ DƏYİŞİR:
  VARIABLES_PER_CELL  3 -> 2   (state.py-dan gəlir)
  reservoir sahəsi    ThreePhaseState -> ReservoirState
  from_vector()       is_saturated arqumenti çıxır (dəyişən keçidi
                      yalnız 3 fazalı black-oil anlayışıdır)

NƏ DƏYİŞMİR:
  updated() imzası eynidir -> coupled_newton.py toxunulmur
  WellUnknowns, BHP məntiqi, line search, max_stable_dt -> toxunulmur

İstifadə:
    python tools/patch_well_state_2phase.py --dry
    python tools/patch_well_state_2phase.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

WELL = ROOT / "imex2d" / "simulation" / "implicit" / "well_state.py"
TEST = ROOT / "tests" / "test_well_state.py"

# (sətir, köhnə_alt_sətir, yeni_alt_sətir, izah)
WELL_EDITS = [
    (48,
     "from .three_phase_state import VARIABLES_PER_CELL, ThreePhaseState",
     "from .state import VARIABLES_PER_CELL, ReservoirState",
     "idxal"),
    (119, "ThreePhaseState", "ReservoirState", "şərh"),
    (122, "reservoir: ThreePhaseState", "reservoir: ReservoirState",
     "sahə tipi"),
    (142, "ThreePhaseState", "ReservoirState", "şərh"),
    (158,
     "def from_vector(cls, vector: np.ndarray, is_saturated: np.ndarray,",
     "def from_vector(cls, vector: np.ndarray,",
     "from_vector imzası"),
    (166, "3-ə bölünmür", "2-ə bölünmür", "xəta mesajı"),
    (167,
     "ThreePhaseState.from_vector(vector[:offset], is_saturated)",
     "ReservoirState.from_vector(vector[:offset])",
     "rezervuar qurulması"),
    (177, "ThreePhaseState", "ReservoirState", "şərh"),
]

TEST_EDITS = [
    (14,
     "from imex2d.simulation.implicit.three_phase_state import ThreePhaseState",
     "from imex2d.simulation.implicit.state import ReservoirState",
     "idxal"),
    (31,
     "    reservoir = ThreePhaseState(np.full(n, 213.5), np.full(n, 0.35),",
     "    reservoir = ReservoirState(np.full(n, 213.5), np.full(n, 0.35))",
     "_setup rezervuar"),
    (97, "model.ncell * 3", "model.ncell * 2", "vektor ölçüsü"),
    (101, "3×3", "2×2", "şərh"),
    (105, "model.ncell * 3", "model.ncell * 2", "quyu ofseti"),
    (122,
     "    restored = CoupledState.from_vector(state.to_vector(),",
     "    restored = CoupledState.from_vector(state.to_vector(), wells.names)",
     "round-trip çağırışı"),
    (135,
     "CoupledState.from_vector(broken, reservoir.is_saturated, wells.names)",
     "CoupledState.from_vector(broken, wells.names)",
     "uzunluq testi"),
]

# rewrite-lardan SONRA silinəcək sətirlər (davamı olan sətirlər)
TEST_DELETIONS = [
    (123, "from_vector çağırışının davamı"),
    (32, "_setup rezervuar sətrinin davamı"),
]

# yoxlamalar: silinəcək sətirlərin doğru olduğunu təsdiqləyir
TEST_DELETE_CHECKS = [
    (32, "np.zeros(n), np.zeros(n, bool))"),
    (123, "reservoir.is_saturated, wells.names)"),
]


def apply_edits(path, edits, delete_checks=(), deletions=()):
    print(f"--- {path.relative_to(ROOT)}")
    if not path.exists():
        print("  XƏTA: fayl tapılmadı")
        return None

    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    print(f"  {len(lines)} sətir")

    failed = 0
    for no, old, _new, label in edits:
        if no > len(lines):
            print(f"  XƏTA sətir {no}: fayl qısadır")
            failed += 1
            continue
        if old not in lines[no - 1]:
            print(f"  XƏTA sətir {no} ({label}):")
            print(f"    gözlənilirdi: {old!r}")
            print(f"    tapıldı:      {lines[no - 1].rstrip()!r}")
            failed += 1

    for no, expected in delete_checks:
        if no > len(lines) or expected not in lines[no - 1]:
            got = lines[no - 1].rstrip() if no <= len(lines) else "(yoxdur)"
            print(f"  XƏTA silmə sətri {no}: gözlənilirdi {expected!r}")
            print(f"                        tapıldı      {got!r}")
            failed += 1

    if failed:
        return None

    print(f"  {len(edits) + len(delete_checks)} yoxlama keçdi")

    for no, old, new, label in edits:
        cur = lines[no - 1]
        # tam sətir əvəzlənməsi (girinti ilə başlayanlar)
        if old.startswith("    ") and old.strip() in cur.strip():
            upd = new + "\n"
        else:
            upd = cur.replace(old, new)
        print(f"    {no:>4} {label}")
        print(f"         - {cur.rstrip()}")
        print(f"         + {upd.rstrip()}")
        lines[no - 1] = upd

    for no, label in sorted(deletions, reverse=True):
        print(f"    sil {no} — {label}: {lines[no - 1].rstrip()}")
        del lines[no - 1]

    print(f"  yeni ölçü: {len(lines)} sətir\n")
    return path, "".join(lines)


def main():
    dry = "--dry" in sys.argv

    a = apply_edits(WELL, WELL_EDITS)
    if a is None:
        print("HEÇ NƏ YAZILMADI.")
        return 1

    b = apply_edits(TEST, TEST_EDITS, TEST_DELETE_CHECKS, TEST_DELETIONS)
    if b is None:
        print("HEÇ NƏ YAZILMADI.")
        return 1

    if dry:
        print("--dry rejimi: fayllar DƏYİŞDİRİLMƏDİ.")
        return 0

    for path, content in (a, b):
        path.write_text(content, encoding="utf-8")
        print(f"Yazıldı: {path.relative_to(ROOT)}")

    print("\nİndi yoxla:")
    print("  pytest -q tests/test_well_state.py")
    print("  python tools/golden.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
