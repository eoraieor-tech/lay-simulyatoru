"""Müşahidə məlumatının CSV-dən oxunması.

Format — uzun (long) cədvəl:

    time,well,quantity,value[,uncertainty]
    30,PROD-1,OIL_RATE,142.5
    30,PROD-1,WATER_RATE,3.1
    30,,AVERAGE_PRESSURE,247.8

`quantity` sütununda qəbul edilən adlar `ALIASES` lüğətindədir —
həm ingilis (Eclipse SUMMARY üslubu: WOPR, WWPR, FPR), həm Azərbaycan
variantları tanınır.

Boş `well` sahəsi yataq səviyyəsi deməkdir.
"""

from __future__ import annotations

import csv
from typing import Dict, List, Optional

import numpy as np

from ..domain.observations import (ObservationSet, ObservedQuantity,
                                   ObservedSeries)

ALIASES: Dict[str, ObservedQuantity] = {}


def _register(quantity: ObservedQuantity, *names: str) -> None:
    for name in names:
        ALIASES[name.strip().upper()] = quantity


_register(ObservedQuantity.OIL_RATE,
          "OIL_RATE", "OILRATE", "QO", "WOPR", "FOPR", "NEFT", "NEFT_DEBITI")
_register(ObservedQuantity.WATER_RATE,
          "WATER_RATE", "WATERRATE", "QW", "WWPR", "FWPR", "SU", "SU_DEBITI")
_register(ObservedQuantity.WATER_INJECTION,
          "WATER_INJECTION", "QWI", "WWIR", "FWIR", "VURMA", "SU_VURMA")
_register(ObservedQuantity.WATER_CUT,
          "WATER_CUT", "WATERCUT", "WCT", "WWCT", "FWCT", "SULASMA", "SULAŞMA")
_register(ObservedQuantity.BOTTOM_HOLE_PRESSURE,
          "BHP", "WBHP", "BOTTOM_HOLE_PRESSURE", "QUYUDIBI_TEZYIQ")
_register(ObservedQuantity.AVERAGE_PRESSURE,
          "AVERAGE_PRESSURE", "FPR", "PRESSURE", "ORTA_TEZYIQ", "LAY_TEZYIQI")
_register(ObservedQuantity.CUMULATIVE_OIL,
          "CUMULATIVE_OIL", "CUM_OIL", "WOPT", "FOPT", "KUM_NEFT")

TIME_KEYS = ("time", "gun", "gün", "day", "days", "date", "t")
WELL_KEYS = ("well", "quyu", "name", "ad")
QUANTITY_KEYS = ("quantity", "kemiyyet", "kəmiyyət", "type", "tip", "parameter")
VALUE_KEYS = ("value", "qiymet", "qiymət", "val")
UNCERTAINTY_KEYS = ("uncertainty", "sigma", "error", "xeta", "xəta")

FIELD_MARKERS = {"", "FIELD", "FIELD-TOTAL", "YATAQ", "TOTAL", "-"}


class ObservationFormatError(Exception):
    """CSV strukturu gözləniləndən fərqlidir."""


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


def read_observations_csv(path: str) -> ObservationSet:
    with open(path, newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(handle, dialect=dialect)
        header = reader.fieldnames or []

        time_column = _pick(header, TIME_KEYS)
        quantity_column = _pick(header, QUANTITY_KEYS)
        value_column = _pick(header, VALUE_KEYS)
        if not (time_column and quantity_column and value_column):
            raise ObservationFormatError(
                "CSV-də 'time', 'quantity' və 'value' sütunları tapılmadı. "
                f"Tapılan: {', '.join(header) or '—'}")

        well_column = _pick(header, WELL_KEYS)
        uncertainty_column = _pick(header, UNCERTAINTY_KEYS)

        buckets: Dict[tuple, Dict[str, list]] = {}
        unknown: set = set()

        for row in reader:
            time = _to_float(row.get(time_column))
            value = _to_float(row.get(value_column))
            if time is None or value is None:
                continue

            raw = (row.get(quantity_column) or "").strip().upper()
            quantity = ALIASES.get(raw)
            if quantity is None:
                unknown.add(raw)
                continue

            well = (row.get(well_column) or "").strip() if well_column else ""
            if well.upper() in FIELD_MARKERS or quantity.is_field_level:
                well = ""

            bucket = buckets.setdefault((well, quantity),
                                        {"t": [], "v": [], "s": []})
            bucket["t"].append(time)
            bucket["v"].append(value)
            bucket["s"].append(_to_float(row.get(uncertainty_column))
                               if uncertainty_column else None)

    if not buckets:
        message = "Heç bir tanınan müşahidə tapılmadı."
        if unknown:
            message += (" Naməlum kəmiyyətlər: "
                        + ", ".join(sorted(unknown)[:8]))
        raise ObservationFormatError(message)

    series: List[ObservedSeries] = []
    for (well, quantity), bucket in buckets.items():
        order = np.argsort(np.asarray(bucket["t"], float))
        sigma = bucket["s"]
        uncertainty = (np.asarray([sigma[i] for i in order], float)
                       if all(value is not None for value in sigma) else None)
        series.append(ObservedSeries(
            well=well, quantity=quantity,
            time=np.asarray(bucket["t"], float)[order],
            values=np.asarray(bucket["v"], float)[order],
            uncertainty=uncertainty))

    dataset = ObservationSet(series=series, source=path)
    issues = dataset.validate()
    if issues:
        raise ObservationFormatError("; ".join(issues))
    return dataset


def write_observations_csv(path: str, observations: ObservationSet) -> str:
    """Nəticəni müşahidə formatında yazır — nümunə fayl yaratmaq üçün."""
    reverse = {}
    for name, quantity in ALIASES.items():
        reverse.setdefault(quantity, name)

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time", "well", "quantity", "value"])
        for item in observations.series:
            name = reverse[item.quantity]
            for time, value in zip(item.time, item.values):
                writer.writerow([f"{time:.3f}", item.well, name,
                                 f"{value:.5f}"])
    return path
