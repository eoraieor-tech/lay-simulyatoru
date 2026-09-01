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

VAHİD METADATASI (Phase 1, istəyə görə, GERİYƏ UYĞUN)
Sütun adının sonuna `[vahid]` əlavə edilə bilər:

    x[ft],y[ft],PERMX[D],PERMX[mD],PORO

Bu, sütunun BAZA adını (mötərizədən əvvəlki hissə, `x`/`PERMX`/...)
DƏYİŞMİR — yalnız həmin sütunun dəyərləri mühərrik vahidinə (m/mD)
ÇEVRİLİR. `[vahid]` verilməyən sütun ƏVVƏLKİ kimi mühərrik vahidində
qəbul edilir (dəyər DƏYİŞMİR) — YALNIZ keçiricilik (PERMX/PERMY/PERMZ)
üçün bu, açıq xəbərdarlıqla (`dataset.warnings`) bildirilir, çünki
mD/Darcy qarışıqlığı bu sahədə real və tez-tez rast gəlinən səhvdir.
X/Y/dərinlik üçün xəbərdarlıq YOXDUR — "metr" bu kod bazasında hər
yerdə sənədləşdirilmiş, mövcud DEFOLT konvensiyadır (bax UNITS.md),
hər faylda xatırlatmaq FAYDASIZ SƏS-KÜY olardı.
"""

from __future__ import annotations

import csv
import re
from typing import Dict, List, Optional, Tuple

from ..domain.unit_conversions import ENGINE_UNITS, to_engine_units
from ..domain.well_data import WellDataset, WellSample

WELL_KEYS = ("well", "quyu", "wellname", "name", "ad")
X_KEYS = ("x", "east", "easting", "x_m")
Y_KEYS = ("y", "north", "northing", "y_m")
LAYER_KEYS = ("k", "layer", "tebeqe", "təbəqə")
DEPTH_KEYS = ("depth", "md", "tvd", "derinlik", "dərinlik")

RESERVED = set(WELL_KEYS + X_KEYS + Y_KEYS + LAYER_KEYS + DEPTH_KEYS)

#: Sütun baza adı (böyük hərflə) -> (kəmiyyət növü, xəbərdarlıq lazımdırmı).
#: Yalnız keçiricilik üçün xəbərdarlıq — bax modul docstring-i.
_PROPERTY_QUANTITY: Dict[str, Tuple[str, bool]] = {
    "PERMX": ("permeability", True),
    "PERMY": ("permeability", True),
    "PERMZ": ("permeability", True),
}

_UNIT_SUFFIX = re.compile(r"^(.*?)\s*\[\s*([^\]]+?)\s*\]\s*$")


def _split_unit_suffix(column_name: str) -> Tuple[str, Optional[str]]:
    """`"PERMX[mD]"` -> `("PERMX", "mD")`; `"PERMX"` -> `("PERMX", None)`."""
    match = _UNIT_SUFFIX.match(column_name.strip())
    if match:
        return match.group(1).strip(), match.group(2).strip()
    return column_name.strip(), None


class WellDataFormatError(Exception):
    """CSV strukturu gözləniləndən fərqlidir."""


def _pick(parsed: Dict[str, Tuple[str, Optional[str]]], candidates
         ) -> Optional[Tuple[str, Optional[str]]]:
    for candidate in candidates:
        if candidate in parsed:
            return parsed[candidate]
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

        # sütun adı -> (baza_adı_lower, orijinal_baş_adı, vahid_ya_None)
        parsed: Dict[str, Tuple[str, Optional[str]]] = {}
        original_by_base: Dict[str, str] = {}
        for name in header:
            if not name:
                continue
            base, unit = _split_unit_suffix(name)
            parsed[base.lower()] = (name, unit)
            original_by_base[base.lower()] = name

        well_hit = _pick(parsed, WELL_KEYS)
        x_hit = _pick(parsed, X_KEYS)
        y_hit = _pick(parsed, Y_KEYS)
        if not (well_hit and x_hit and y_hit):
            raise WellDataFormatError(
                "CSV-də 'well', 'x' və 'y' sütunları tapılmadı. "
                f"Tapılan sütunlar: {', '.join(header) or '—'}")
        well_column, _ = well_hit
        (x_column, x_unit), (y_column, y_unit) = x_hit, y_hit

        layer_hit = _pick(parsed, LAYER_KEYS)
        depth_hit = _pick(parsed, DEPTH_KEYS)
        layer_column = layer_hit[0] if layer_hit else None
        depth_column, depth_unit = depth_hit if depth_hit else (None, None)

        reserved_bases = set(WELL_KEYS + X_KEYS + Y_KEYS + LAYER_KEYS + DEPTH_KEYS)
        # (baza_adı_upper, orijinal_sütun_adı, vahid_ya_None)
        property_columns = [
            (base_lower.upper(), original, unit)
            for base_lower, (original, unit) in parsed.items()
            if base_lower not in reserved_bases]

        warnings: List[str] = []
        warned_properties = set()
        samples: List[WellSample] = []
        for row_number, row in enumerate(reader, start=2):
            x = _to_float(row.get(x_column))
            y = _to_float(row.get(y_column))
            if x is None or y is None:
                continue
            if x_unit:
                x = to_engine_units(x, x_unit, "length")
            if y_unit:
                y = to_engine_units(y, y_unit, "length")
            values: Dict[str, float] = {}
            for prop_name, column, unit in property_columns:
                value = _to_float(row.get(column))
                if value is None:
                    continue
                quantity_info = _PROPERTY_QUANTITY.get(prop_name)
                if unit:
                    quantity = quantity_info[0] if quantity_info else None
                    if quantity is not None:
                        value = to_engine_units(value, unit, quantity)
                    elif prop_name not in warned_properties:
                        warnings.append(
                            f"'{prop_name}[{unit}]': bu xassə üçün vahid çevirməsi "
                            "dəstəklənmir — dəyər OLDUĞU KİMİ (çevrilmədən) qəbul edildi.")
                        warned_properties.add(prop_name)
                elif quantity_info and quantity_info[1] and prop_name not in warned_properties:
                    engine_unit = ENGINE_UNITS[quantity_info[0]]
                    warnings.append(
                        f"'{prop_name}': vahid göstərilməyib (məs. '{prop_name}[{engine_unit}]' "
                        f"və ya '{prop_name}[D]') — mövcud mühərrik vahidi ({engine_unit}) "
                        "qəbul edildi.")
                    warned_properties.add(prop_name)
                values[prop_name] = value
            if not values:
                continue
            layer = None
            if layer_column:
                raw = _to_float(row.get(layer_column))
                layer = int(raw) - 1 if raw is not None else None
            depth = _to_float(row.get(depth_column)) if depth_column else None
            if depth is not None and depth_unit:
                depth = to_engine_units(depth, depth_unit, "length")
            samples.append(WellSample(
                well=(row.get(well_column) or f"W-{row_number}").strip(),
                x=x, y=y, values=values, layer=layer, depth=depth))

    dataset = WellDataset(samples=samples, source=path, warnings=warnings)
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
