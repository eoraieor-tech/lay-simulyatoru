"""Eclipse summary-ni konsola çap edir — `imex2d.io.eclipse_summary` üzərində.

    python tools/eclipse_summary.py <baza_yolu> "TIME,FOPR,FGOR,WBHP:PROD" [addım]

`<baza_yolu>` — uzantısız fayl adı (`.SMSPEC` və `.UNSMRY` yanında olmalıdır).
Oxuyucu Seans 38-də `imex2d/io/eclipse_summary.py`-yə köçürülüb.
"""

import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from imex2d.io.eclipse_summary import read_summary  # noqa: E402


def main(argv):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")
    summary = read_summary(argv[1])
    wanted = argv[2].split(",")
    print("mövcud:", ", ".join(summary.labels))
    missing = [w for w in wanted if w not in summary]
    if missing:
        print("TAPILMADI:", missing)
    wanted = [w for w in wanted if w in summary]
    print("\t".join(wanted))
    step = int(argv[3]) if len(argv) > 3 else 1
    columns = [summary.series(w) for w in wanted]
    count = summary.values.shape[0]
    for n in range(count):
        if n % step == 0 or n == count - 1:
            print("\t".join(f"{column[n]:.6g}" for column in columns))


if __name__ == "__main__":
    main(sys.argv)
