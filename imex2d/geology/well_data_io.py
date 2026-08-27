"""Quyu məlumatının CSV-dən oxunması və nümunə faylın yaradılması.

Gözlənilən format (sütun adları böyük-kiçik hərfə həssas deyil):

    well,x,y[,k][,depth],PORO,PERMX[,NTG,...]
    W-1,120,340,1,2005,0.21,180
    W-1,120,340,2,2011,0.19,140
    W-2,600,220,1,2008,0.24,320

Qaydalar:
  · `well`, `x`, `y` məcburidir
  · `k` (təbəqə, 1-dən başlayır) istəyə görədir — verilməsə xassə
    bütün təbəqələrə eyni tətbiq olunur
  · qalan bütün ədədi sütunlar xassə kimi qəbul edilir
"""

from __future__ import annotations

import csv
from typing import Dict, List, Optional

from ..domain.well_data import WellDataset, WellSample

WELL_KEYS = ("well", "quyu", "wellname", "name", "ad")
X_KEYS = ("x", "east", "easting", "x_m")
Y_KEYS = ("y", "north", "northing", "y_m")
LAYER_KEYS = ("k", "layer", "tebeqe", "təbəqə")
DEPTH_KEYS = ("depth", "md", "tvd", "derinlik", "dərinlik")

RESERVED = set(WELL_KEYS + X_KEYS + Y_KEYS + LAYER_KEYS + DEPTH_KEYS)


class WellDataFormatError(Exception):
    """CSV strukturu gözləniləndən fərqlidir."""


def _pick(header: List[str], candidates) -> Optional[str]:
    lowered = {name.strip().lower(): name for name in header}
    for candidate in candidates:
        if candidate in lowered:
            return lowered[candidate]
    return None


def _to_float(text: str) -> Optional[float]:
    text = (text or "").strip().replace(",", ".")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def read_well_csv(path: str) -> WellDataset:
    with open(path, newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(handle, dialect=dialect)
        header = reader.fieldnames or []

        well_column = _pick(header, WELL_KEYS)
        x_column = _pick(header, X_KEYS)
        y_column = _pick(header, Y_KEYS)
        if not (well_column and x_column and y_column):
            raise WellDataFormatError(
                "CSV-də 'well', 'x' və 'y' sütunları tapılmadı. "
                f"Tapılan sütunlar: {', '.join(header) or '—'}")

        layer_column = _pick(header, LAYER_KEYS)
        depth_column = _pick(header, DEPTH_KEYS)
        property_columns = [name for name in header
                            if name and name.strip().lower() not in RESERVED]

        samples: List[WellSample] = []
        for row_number, row in enumerate(reader, start=2):
            x = _to_float(row.get(x_column))
            y = _to_float(row.get(y_column))
            if x is None or y is None:
                continue
            values: Dict[str, float] = {}
            for column in property_columns:
                value = _to_float(row.get(column))
                if value is not None:
                    values[column.strip().upper()] = value
            if not values:
                continue
            layer = None
            if layer_column:
                raw = _to_float(row.get(layer_column))
                layer = int(raw) - 1 if raw is not None else None
            samples.append(WellSample(
                well=(row.get(well_column) or f"W-{row_number}").strip(),
                x=x, y=y, values=values, layer=layer,
                depth=_to_float(row.get(depth_column)) if depth_column else None))

    dataset = WellDataset(samples=samples, source=path)
    issues = dataset.validate()
    if issues:
        raise WellDataFormatError("; ".join(issues))
    return dataset


def write_example_csv(path: str, nx=41, ny=41, dx=20.0, dy=20.0, nz=1) -> str:
    """Nümunə quyu faylı — istifadəçinin formatı görməsi üçün."""
    import numpy as np

    rng = np.random.default_rng(3)
    length_x, length_y = nx * dx, ny * dy
    locations = [("W-1", 0.15, 0.15), ("W-2", 0.80, 0.20),
                 ("W-3", 0.25, 0.75), ("W-4", 0.85, 0.85),
                 ("W-5", 0.50, 0.50)]

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["well", "x", "y", "k", "PORO", "PERMX", "NTG"])
        for name, fx, fy in locations:
            x, y = fx * length_x, fy * length_y
            base_porosity = 0.16 + 0.10 * fx * fy + rng.normal(0, 0.01)
            base_permeability = 60.0 * np.exp(3.0 * base_porosity * 4.0)
            for k in range(1, nz + 1):
                factor = 1.0 - 0.12 * (k - 1)
                writer.writerow([
                    name, f"{x:.1f}", f"{y:.1f}", k,
                    f"{max(base_porosity * factor, 0.05):.4f}",
                    f"{max(base_permeability * factor ** 2, 1.0):.1f}",
                    f"{min(0.95, 0.7 + 0.05 * k):.2f}"])
    return path
