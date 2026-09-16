"""Eclipse summary oxuyucusu (SMSPEC + UNSMRY) — Seans 38.

Test faylları burada yazılır (big-endian Fortran qeydləri), real OPM faylı
repoda deyil (lisenziya yoxlanılmayıb, bax `SPE1.md` §1).
"""

from __future__ import annotations

import os
import struct
import tempfile

import numpy as np
import pytest

from imex2d.io.eclipse_summary import (EclipseSummaryError, read_keywords,
                                       read_summary)


def _record(payload: bytes) -> bytes:
    return struct.pack(">i", len(payload)) + payload + struct.pack(">i", len(payload))


def _keyword(name: str, kind: str, values: list, block: int = None) -> bytes:
    head = name.ljust(8).encode() + struct.pack(">i", len(values)) + kind.encode()
    out = _record(head)
    block = block or (105 if kind == "CHAR" else 1000)
    for start in range(0, len(values), block):
        chunk = values[start:start + block]
        if kind == "CHAR":
            payload = b"".join(v.ljust(8).encode() for v in chunk)
        else:
            fmt = {"INTE": "i", "REAL": "f", "DOUB": "d"}[kind]
            payload = struct.pack(">" + fmt * len(chunk), *chunk)
        out += _record(payload)
    return out


def _write_case(rows, keywords=("TIME", "FOPR", "WBHP", "BPR"),
                wells=(":+:+:+:+", ":+:+:+:+", "PROD", ":+:+:+:+"),
                nums=(0, 0, 0, 300), units=("DAYS", "STB/DAY", "PSIA", "PSIA")):
    directory = tempfile.mkdtemp()
    base = os.path.join(directory, "CASE")
    with open(base + ".SMSPEC", "wb") as handle:
        handle.write(_keyword("KEYWORDS", "CHAR", list(keywords)))
        handle.write(_keyword("WGNAMES", "CHAR", list(wells)))
        handle.write(_keyword("NUMS", "INTE", list(nums)))
        handle.write(_keyword("UNITS", "CHAR", list(units)))
    with open(base + ".UNSMRY", "wb") as handle:
        for step, row in enumerate(rows):
            handle.write(_keyword("SEQHDR", "INTE", [step]))
            handle.write(_keyword("MINISTEP", "INTE", [step]))
            handle.write(_keyword("PARAMS", "REAL", list(row)))
    return base


ROWS = [(1.0, 20000.0, 2905.09, 5063.13),
        (1550.0, 19990.0, 1000.0, 4800.0),
        (3650.0, 5732.65, 1000.0, 4101.32)]


def test_labels_are_built_like_eclipse():
    summary = read_summary(_write_case(ROWS))
    assert summary.labels == ["TIME", "FOPR", "WBHP:PROD", "BPR:300"]
    assert summary.unit("WBHP:PROD") == "PSIA"


def test_series_values_are_read_in_order():
    summary = read_summary(_write_case(ROWS))
    assert summary.series("TIME") == pytest.approx([1.0, 1550.0, 3650.0])
    assert summary.series("FOPR")[-1] == pytest.approx(5732.65, rel=1e-6)
    assert "WBHP:PROD" in summary and "WBHP:INJ" not in summary


def test_value_at_time_is_interpolated():
    summary = read_summary(_write_case(ROWS))
    assert summary.at("FOPR", 3650.0) == pytest.approx(5732.65, rel=1e-6)
    midpoint = 0.5 * (1550.0 + 3650.0)
    assert summary.at("FOPR", midpoint) == pytest.approx(0.5 * (19990.0 + 5732.65),
                                                         rel=1e-6)


def test_long_arrays_are_split_into_blocks():
    """1000-dən uzun REAL massivi bir neçə qeyd blokunda yazılır."""
    directory = tempfile.mkdtemp()
    path = os.path.join(directory, "BIG")
    values = list(np.arange(2500, dtype=float))
    with open(path, "wb") as handle:
        handle.write(_keyword("PARAMS", "REAL", values))
        handle.write(_keyword("NAMES", "CHAR", [f"W{i}" for i in range(230)]))
    keywords = dict(read_keywords(path))
    assert keywords["PARAMS"] == pytest.approx(values)
    assert keywords["NAMES"][-1] == "W229"


def test_unknown_label_is_an_explicit_error():
    summary = read_summary(_write_case(ROWS))
    with pytest.raises(KeyError, match="WOPR:PROD"):
        summary.series("WOPR:PROD")


def test_truncated_file_is_an_explicit_error():
    base = _write_case(ROWS)
    with open(base + ".UNSMRY", "rb") as handle:
        data = handle.read()
    with open(base + ".UNSMRY", "wb") as handle:
        handle.write(data[:-10])
    with pytest.raises(EclipseSummaryError):
        read_summary(base)


def test_row_width_mismatch_is_an_explicit_error():
    base = _write_case([(1.0, 2.0, 3.0)])
    with pytest.raises(EclipseSummaryError, match="uyğun deyil"):
        read_summary(base)
