"""SCAL cədvəllərinin oxunması — B7.

İki format dəstəklənir:

CSV (laboratoriya hesabatlarının adi forması)

    region,sw,krw,kro,pc
    1,0.20,0.000,0.800,0.00
    1,0.30,0.015,0.520,0.00
    2,0.30,0.000,0.750,0.00

`region` sütunu olmasa, hamısı 1-ci regiona düşür.

ECLIPSE SWOF (deck-in `PROPS` bölməsindən)

    SWOF
    -- Sw     krw     kro     Pc
      0.20   0.000   0.800   0.0
      0.80   0.350   0.000   0.0
    /
      0.30   0.000   0.750   0.0
      ...
    /

Hər `/` yeni regionu bitirir — Eclipse-in öz qaydası.
"""

from __future__ import annotations

import csv
import re
from typing import Dict, List, Optional

import numpy as np

from ..domain.scal_tables import SaturationTable, SaturationTableSet

REGION_KEYS = ("region", "satnum", "zona", "zone")
SW_KEYS = ("sw", "s_w", "water_saturation", "su")
KRW_KEYS = ("krw", "kr_w", "krwater")
KRO_KEYS = ("kro", "kr_o", "kroil", "krow")
PC_KEYS = ("pc", "pcow", "capillary", "kapilyar")


class ScalFormatError(Exception):
    """Fayl strukturu gözləniləndən fərqlidir."""


def _pick(header: List[str], candidates) -> Optional[str]:
    lowered = {name.strip().lower(): name for name in header if name}
    for candidate in candidates:
        if candidate in lowered:
            return lowered[candidate]
    return None


def _to_float(text: Optional[str]) -> Optional[float]:
    text = (text or "").strip().replace(",", ".")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def read_scal_csv(path: str) -> SaturationTableSet:
    with open(path, newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(handle, dialect=dialect)
        header = reader.fieldnames or []

        sw_column = _pick(header, SW_KEYS)
        krw_column = _pick(header, KRW_KEYS)
        kro_column = _pick(header, KRO_KEYS)
        if not (sw_column and krw_column and kro_column):
            raise ScalFormatError(
                "CSV-də 'sw', 'krw' və 'kro' sütunları tapılmadı. "
                f"Tapılan: {', '.join(header) or '—'}")

        region_column = _pick(header, REGION_KEYS)
        pc_column = _pick(header, PC_KEYS)

        buckets: Dict[int, Dict[str, list]] = {}
        for row in reader:
            sw = _to_float(row.get(sw_column))
            krw = _to_float(row.get(krw_column))
            kro = _to_float(row.get(kro_column))
            if sw is None or krw is None or kro is None:
                continue
            region = 1
            if region_column:
                value = _to_float(row.get(region_column))
                region = int(value) if value is not None else 1
            bucket = buckets.setdefault(region,
                                        {"sw": [], "krw": [], "kro": [],
                                         "pc": []})
            bucket["sw"].append(sw)
            bucket["krw"].append(krw)
            bucket["kro"].append(kro)
            bucket["pc"].append(_to_float(row.get(pc_column))
                                if pc_column else None)

    if not buckets:
        raise ScalFormatError("Faylda oxunacaq sətir tapılmadı.")
    return _build(buckets, path)


def read_swof(path: str) -> SaturationTableSet:
    """Eclipse deck-dən `SWOF` açar sözünü oxuyur."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        text = handle.read()

    match = re.search(r"^\s*SWOF\s*$", text, re.MULTILINE | re.IGNORECASE)
    if match is None:
        raise ScalFormatError("Faylda SWOF açar sözü tapılmadı.")

    body = text[match.end():]
    stop = re.search(r"^\s*[A-Z][A-Z0-9_]{2,}\s*$", body, re.MULTILINE)
    if stop is not None:
        body = body[:stop.start()]

    buckets: Dict[int, Dict[str, list]] = {}
    region = 1
    for raw_line in body.splitlines():
        line = re.sub(r"--.*", "", raw_line).strip()
        if not line:
            continue
        if line.startswith("/"):
            if region in buckets:
                region += 1          # `/` regionu bitirir
            continue

        terminated = line.endswith("/")
        numbers = [float(token) for token in
                   re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?",
                              line.rstrip("/"))]
        if len(numbers) >= 3:
            bucket = buckets.setdefault(region,
                                        {"sw": [], "krw": [], "kro": [],
                                         "pc": []})
            bucket["sw"].append(numbers[0])
            bucket["krw"].append(numbers[1])
            bucket["kro"].append(numbers[2])
            bucket["pc"].append(numbers[3] if len(numbers) > 3 else None)
        if terminated and region in buckets:
            region += 1

    if not buckets:
        raise ScalFormatError("SWOF bölməsində rəqəm tapılmadı.")
    return _build(buckets, path)


def _build(buckets: Dict[int, dict], source: str) -> SaturationTableSet:
    table_set = SaturationTableSet()
    for region, bucket in sorted(buckets.items()):
        order = np.argsort(np.asarray(bucket["sw"], float))
        capillary = bucket["pc"]
        pc = (np.asarray([capillary[index] for index in order], float)
              if all(value is not None for value in capillary) else None)
        table_set.add(region, SaturationTable(
            sw=np.asarray(bucket["sw"], float)[order],
            krw=np.asarray(bucket["krw"], float)[order],
            kro=np.asarray(bucket["kro"], float)[order],
            pc=pc, name=f"Region {region}"))

    issues = table_set.validate()
    if issues:
        raise ScalFormatError("; ".join(issues))
    return table_set


def write_scal_csv(path: str, tables: SaturationTableSet) -> str:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["region", "sw", "krw", "kro", "pc"])
        for region, table in sorted(tables.tables.items()):
            pc = (table.pc if table.pc is not None
                  else np.zeros_like(table.sw))
            for index in range(table.sw.size):
                writer.writerow([region, f"{table.sw[index]:.5f}",
                                 f"{table.krw[index]:.6f}",
                                 f"{table.kro[index]:.6f}",
                                 f"{pc[index]:.5f}"])
    return path
