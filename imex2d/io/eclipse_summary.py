"""Eclipse summary oxuyucusu (SMSPEC + UNSMRY) — `resdata`-sız.

SPE1 etalonunu (OPM Flow nəticəsi) oxumaq üçün yazılıb (Seans 29, Q-23),
Seans 38-də `tools/`-dan bura köçürülüb və testlə örtülüb.

Format: big-endian Fortran "unformatted" qeydləri. Hər açar söz:

    [16] 8 simvol ad + int32 say + 4 simvol tip [16]
    sonra məlumat bloklarla (INTE/REAL/DOUB/LOGI ≤ 1000, CHAR ≤ 105 element)

Etiketlər: sahə açarları `FOPR`, quyu açarları `WBHP:PROD`, blok açarları
`BPR:1` (1-dən başlayan qlobal hüceyrə nömrəsi).
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

_SIZES = {b"INTE": ">i", b"REAL": ">f", b"DOUB": ">d", b"LOGI": ">i"}


class EclipseSummaryError(Exception):
    """Summary faylı oxuna bilmədi — aydın, tutulan xəta."""


def read_keywords(path: str) -> List[Tuple[str, list]]:
    """`[(ad, dəyərlər), ...]` — faylda göründüyü sıra ilə."""
    with open(path, "rb") as handle:
        data = handle.read()
    pos, out = 0, []
    try:
        while pos < len(data):
            (length,) = struct.unpack(">i", data[pos:pos + 4])
            head = data[pos + 4:pos + 4 + length]
            pos += 8 + length
            name = head[:8].decode("ascii").strip()
            (count,) = struct.unpack(">i", head[8:12])
            kind = head[12:16]
            if kind != b"CHAR" and kind not in _SIZES:
                raise EclipseSummaryError(f"{path}: naməlum tip {kind!r} ({name})")
            values: list = []
            remaining = count
            block = 105 if kind == b"CHAR" else 1000
            while remaining > 0:
                (size,) = struct.unpack(">i", data[pos:pos + 4])
                chunk = data[pos + 4:pos + 4 + size]
                pos += 8 + size
                k = min(block, remaining)
                if kind == b"CHAR":
                    values += [chunk[i * 8:(i + 1) * 8].decode("ascii").strip()
                               for i in range(k)]
                else:
                    values += list(struct.unpack(">" + _SIZES[kind][1] * k, chunk))
                remaining -= k
            out.append((name, values))
    except struct.error as exc:
        raise EclipseSummaryError(f"{path}: fayl kəsilib və ya pozulub") from exc
    return out


@dataclass
class EclipseSummary:
    """Oxunmuş summary: hər etiket üçün zaman sırası."""

    labels: List[str]
    units: List[str]
    values: np.ndarray            # (addım sayı, etiket sayı)

    def __post_init__(self):
        self._index: Dict[str, int] = {}
        for position, label in enumerate(self.labels):
            self._index.setdefault(label, position)

    def __contains__(self, label: str) -> bool:
        return label in self._index

    def series(self, label: str) -> np.ndarray:
        if label not in self._index:
            raise KeyError(f"summary-də '{label}' yoxdur")
        return self.values[:, self._index[label]]

    def unit(self, label: str) -> str:
        return self.units[self._index[label]] if self.units else ""

    def at(self, label: str, time: float, time_label: str = "TIME") -> float:
        """`time` anındakı qiymət — zaman üzrə xətti interpolyasiya."""
        return float(np.interp(time, self.series(time_label), self.series(label)))


def read_summary(base: str) -> EclipseSummary:
    """`<base>.SMSPEC` + `<base>.UNSMRY` (uzantısız yol)."""
    spec: Dict[str, list] = {}
    for name, values in read_keywords(base + ".SMSPEC"):
        spec.setdefault(name, values)
    for required in ("KEYWORDS", "NUMS"):
        if required not in spec:
            raise EclipseSummaryError(f"{base}.SMSPEC: '{required}' açar sözü yoxdur")
    keys = spec["KEYWORDS"]
    wells = spec.get("WGNAMES", spec.get("NAMES"))
    nums = spec["NUMS"]
    labels = []
    for i, key in enumerate(keys):
        if key.startswith("W") and wells is not None:
            labels.append(f"{key}:{wells[i]}")
        elif key.startswith("B"):
            labels.append(f"{key}:{nums[i]}")
        else:
            labels.append(key)
    rows = [values for name, values in read_keywords(base + ".UNSMRY")
            if name == "PARAMS"]
    if not rows:
        raise EclipseSummaryError(f"{base}.UNSMRY: PARAMS qeydi yoxdur")
    width = len(labels)
    if any(len(row) != width for row in rows):
        raise EclipseSummaryError(f"{base}: PARAMS uzunluğu KEYWORDS ilə uyğun deyil")
    return EclipseSummary(labels, list(spec.get("UNITS") or []),
                          np.asarray(rows, dtype=float))
