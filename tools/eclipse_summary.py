"""Minimal Eclipse summary oxuyucusu (SMSPEC + UNSMRY) — `resdata`-sız.

SPE1 etalonunu (OPM Flow nəticəsi) oxumaq üçün yazılıb (Seans 29, Q-23).
Format: big-endian Fortran "unformatted" qeydləri. Hər açar söz:

    [16] 8 simvol ad + int32 say + 4 simvol tip [16]
    sonra məlumat bloklarla (INTE/REAL/DOUB ≤ 1000, CHAR ≤ 105 element)

İstifadə:

    python tools/eclipse_summary.py <baza_yolu> "TIME,FOPR,FGOR,WBHP:PROD" [addım]

`<baza_yolu>` — uzantısız fayl adı (`.SMSPEC` və `.UNSMRY` yanında olmalıdır).
Quyu açarları `WBHP:PROD`, blok açarları `BPR:1` (hüceyrə nömrəsi) şəklindədir.

⏳ Hələ test olunmayıb və `imex2d`-yə qoşulmayıb — SPE1 müqayisə testi yazılanda
`imex2d/io/`-ya köçürülməli və testlə örtülməlidir.
"""

import io
import struct
import sys

SIZES = {b"INTE": (4, ">i"), b"REAL": (4, ">f"), b"DOUB": (8, ">d"),
         b"LOGI": (4, ">i"), b"CHAR": (8, None)}


def read_keywords(path):
    """`[(ad, dəyərlər), ...]` — faylda göründüyü sıra ilə."""
    data = open(path, "rb").read()
    pos, out = 0, []
    while pos < len(data):
        (n,) = struct.unpack(">i", data[pos:pos + 4])
        head = data[pos + 4:pos + 4 + n]
        pos += 8 + n
        name = head[:8].decode().strip()
        (count,) = struct.unpack(">i", head[8:12])
        kind = head[12:16]
        _, fmt = SIZES[kind]
        values, remaining = [], count
        block = 105 if kind == b"CHAR" else 1000
        while remaining > 0:
            (m,) = struct.unpack(">i", data[pos:pos + 4])
            chunk = data[pos + 4:pos + 4 + m]
            pos += 8 + m
            k = min(block, remaining)
            if kind == b"CHAR":
                values += [chunk[i * 8:(i + 1) * 8].decode().strip() for i in range(k)]
            else:
                values += list(struct.unpack(">" + fmt[1] * k, chunk))
            remaining -= k
        out.append((name, values))
    return out


def summary(base):
    """`(etiketlər, sətirlər, vahidlər)` — hər sətir bir hesabat addımıdır."""
    spec = {}
    for name, values in read_keywords(base + ".SMSPEC"):
        spec.setdefault(name, values)
    keys = spec["KEYWORDS"]
    wells = spec.get("WGNAMES", spec.get("NAMES"))
    nums = spec["NUMS"]
    labels = []
    for i, key in enumerate(keys):
        if key.startswith("W"):
            labels.append(f"{key}:{wells[i]}")
        elif key.startswith("B"):
            labels.append(f"{key}:{nums[i]}")
        else:
            labels.append(key)
    rows = [values for name, values in read_keywords(base + ".UNSMRY")
            if name == "PARAMS"]
    return labels, rows, spec.get("UNITS")


def main(argv):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")
    labels, rows, units = summary(argv[1])
    wanted = argv[2].split(",")
    index = {label: i for i, label in enumerate(labels)}
    print("mövcud:", ", ".join(labels))
    if units:
        print("vahidlər:", dict(zip(labels, units)))
    missing = [w for w in wanted if w not in index]
    if missing:
        print("TAPILMADI:", missing)
    wanted = [w for w in wanted if w in index]
    print("\t".join(wanted))
    step = int(argv[3]) if len(argv) > 3 else 1
    for n, row in enumerate(rows):
        if n % step == 0 or n == len(rows) - 1:
            print("\t".join(f"{row[index[w]]:.6g}" for w in wanted))


if __name__ == "__main__":
    main(sys.argv)
