"""Fay tərifinin oxunması — B3.

İki format dəstəklənir:

CSV (əl ilə hazırlanmış siyahı)

    name,axis,plane_index,a_low,a_high,b_low,b_high,multiplier,sealing
    F1,I,10,0,40,0,4,0.1,0
    F2,J,20,,,,,0,1

`a_low..b_high` boş buraxılarsa, həmin ox üzrə bütün grid əhatə
olunur. `sealing=1` verilibsə çarpan nəzərə alınmır — fay tam bağlıdır.

Eclipse deck (`FAULTS` + `MULTFLT`)

    FAULTS
      'F1'  11 11  1 41  1 5  'I' /
    /
    MULTFLT
      'F1'  0.1 /
    /

Eclipse-də I/J/K 1-based və DİAPAZONDUR (I1 I2 J1 J2 K1 K2); bir
istiqamətdə fay tək bir sərhəddir, ona görə I1=I2 (və ya J1=J2/K1=K2)
gözlənilir. `MULTFLT` verilməyibsə çarpan 1.0 (şəffaf, yalnız
qeydiyyat) qəbul edilir.
"""

from __future__ import annotations

import csv
import re
from typing import Dict, List, Optional

from ..domain.structure import FaultReference

FACE_TO_AXIS = {"I": "I", "I+": "I", "I-": "I",
                "J": "J", "J+": "J", "J-": "J",
                "K": "K", "K+": "K", "K-": "K"}

NAME_KEYS = ("name", "ad")
AXIS_KEYS = ("axis", "ox")
PLANE_KEYS = ("plane_index", "plane", "muzavi")
A_LOW_KEYS = ("a_low", "alow")
A_HIGH_KEYS = ("a_high", "ahigh")
B_LOW_KEYS = ("b_low", "blow")
B_HIGH_KEYS = ("b_high", "bhigh")
MULT_KEYS = ("multiplier", "carpan", "çarpan", "mult")
SEAL_KEYS = ("sealing", "bagli", "bağlı")


class FaultFormatError(Exception):
    """Fayl strukturu gözləniləndən fərqlidir."""


def _pick(header: List[str], candidates) -> Optional[str]:
    lowered = {name.strip().lower(): name for name in header if name}
    for candidate in candidates:
        if candidate in lowered:
            return lowered[candidate]
    return None


def _to_int(text: Optional[str]) -> Optional[int]:
    text = (text or "").strip()
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _to_float(text: Optional[str], default: float = 1.0) -> float:
    text = (text or "").strip().replace(",", ".")
    if not text:
        return default
    try:
        return float(text)
    except ValueError:
        return default


def read_faults_csv(path: str) -> List[FaultReference]:
    with open(path, newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(handle, dialect=dialect)
        header = reader.fieldnames or []

        name_column = _pick(header, NAME_KEYS)
        axis_column = _pick(header, AXIS_KEYS)
        plane_column = _pick(header, PLANE_KEYS)
        if not (name_column and axis_column and plane_column):
            raise FaultFormatError(
                "CSV-də 'name', 'axis' və 'plane_index' sütunları tapılmadı. "
                f"Tapılan: {', '.join(header) or '—'}")

        a_low_c = _pick(header, A_LOW_KEYS)
        a_high_c = _pick(header, A_HIGH_KEYS)
        b_low_c = _pick(header, B_LOW_KEYS)
        b_high_c = _pick(header, B_HIGH_KEYS)
        mult_c = _pick(header, MULT_KEYS)
        seal_c = _pick(header, SEAL_KEYS)

        faults = []
        for row_number, row in enumerate(reader, start=2):
            name = (row.get(name_column) or "").strip()
            axis = (row.get(axis_column) or "").strip().upper()
            plane = _to_int(row.get(plane_column))
            if not name or axis not in ("I", "J", "K") or plane is None:
                continue

            def bounds(low_key, high_key):
                low = _to_int(row.get(low_key)) if low_key else None
                high = _to_int(row.get(high_key)) if high_key else None
                return (low, high) if low is not None and high is not None else None

            sealing_raw = (row.get(seal_c) or "").strip().lower() if seal_c else ""
            faults.append(FaultReference(
                name=name, source_id=name,
                transmissibility_multiplier=_to_float(
                    row.get(mult_c) if mult_c else None),
                sealing=sealing_raw in ("1", "true", "bəli", "evet", "yes"),
                axis=axis, plane_index=plane,
                range_a=bounds(a_low_c, a_high_c),
                range_b=bounds(b_low_c, b_high_c)))

    if not faults:
        raise FaultFormatError("Faylda oxunacaq fault sətri tapılmadı.")
    return faults


def write_faults_csv(path: str, faults: List[FaultReference]) -> str:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["name", "axis", "plane_index", "a_low", "a_high",
                         "b_low", "b_high", "multiplier", "sealing"])
        for fault in faults:
            a = fault.range_a or ("", "")
            b = fault.range_b or ("", "")
            writer.writerow([fault.name, fault.axis or "", fault.plane_index
                             if fault.plane_index is not None else "",
                             a[0], a[1], b[0], b[1],
                             f"{fault.transmissibility_multiplier:g}",
                             1 if fault.sealing else 0])
    return path


# ═══════════════════════════════════════════════════ Eclipse FAULTS/MULTFLT

def read_eclipse_faults(path: str) -> List[FaultReference]:
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        text = handle.read()

    faults_match = re.search(r"^\s*FAULTS\s*$", text, re.MULTILINE | re.IGNORECASE)
    if faults_match is None:
        raise FaultFormatError("Faylda FAULTS açar sözü tapılmadı.")

    geometry: Dict[str, dict] = {}
    body = text[faults_match.end():]
    stop = re.search(r"^\s*[A-Z][A-Z0-9_]{2,}\s*$", body, re.MULTILINE)
    section = body if stop is None else body[:stop.start()]

    for raw_line in section.splitlines():
        line = re.sub(r"--.*", "", raw_line).strip()
        if not line or line == "/":
            continue
        tokens = re.findall(r"'[^']*'|\S+", line.rstrip("/"))
        tokens = [token.strip("'") for token in tokens]
        if len(tokens) < 7:
            continue
        name = tokens[0]
        try:
            i1, i2, j1, j2, k1, k2 = (int(float(value)) for value in tokens[1:7])
        except ValueError:
            continue
        face = tokens[7].upper() if len(tokens) > 7 else "I"
        axis = FACE_TO_AXIS.get(face, "I")

        # Eclipse 1-based, daxildir. Bu, iki hüceyrə arasındakı sərhəddir:
        # tək sətirli oxda (i1==i2) plane_index = i1 - 1 (0-based, aşağı tərəf).
        if axis == "I":
            plane = i1 - 1
            range_a = (j1 - 1, j2 - 1)
            range_b = (k1 - 1, k2 - 1)
        elif axis == "J":
            plane = j1 - 1
            range_a = (i1 - 1, i2 - 1)
            range_b = (k1 - 1, k2 - 1)
        else:
            plane = k1 - 1
            range_a = (i1 - 1, i2 - 1)
            range_b = (j1 - 1, j2 - 1)

        geometry[name] = dict(axis=axis, plane_index=plane,
                              range_a=range_a, range_b=range_b)

    if not geometry:
        raise FaultFormatError("FAULTS bölməsində sətir tapılmadı.")

    multipliers: Dict[str, float] = {}
    sealing: Dict[str, bool] = {}
    mult_match = re.search(r"^\s*MULTFLT\s*$", text, re.MULTILINE | re.IGNORECASE)
    if mult_match is not None:
        mult_body = text[mult_match.end():]
        stop = re.search(r"^\s*[A-Z][A-Z0-9_]{2,}\s*$", mult_body, re.MULTILINE)
        mult_section = mult_body if stop is None else mult_body[:stop.start()]
        for raw_line in mult_section.splitlines():
            line = re.sub(r"--.*", "", raw_line).strip()
            if not line or line == "/":
                continue
            tokens = re.findall(r"'[^']*'|\S+", line.rstrip("/"))
            tokens = [token.strip("'") for token in tokens]
            if len(tokens) < 2:
                continue
            try:
                value = float(tokens[1])
            except ValueError:
                continue
            multipliers[tokens[0]] = value
            sealing[tokens[0]] = value <= 1e-12

    return [FaultReference(
        name=name, source_id=name,
        transmissibility_multiplier=multipliers.get(name, 1.0),
        sealing=sealing.get(name, False), **fields)
        for name, fields in geometry.items()]
