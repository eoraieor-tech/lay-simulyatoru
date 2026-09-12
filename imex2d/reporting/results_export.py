"""Simulyasiya nəticələrinin CSV / JSON ixracı — B5-a.

NİYƏ LAZIMDIR. Layihədə yalnız **PDF hesabat** vardı (`report.py`) — o,
insana baxmaq üçündür, maşına yox. B4-A `well_bhp` və `well_thp`
sıralarını yaratdıqdan sonra istifadəçinin əlində qrafikdən başqa heç
nə qalmırdı: rəqəmləri Excel-ə və ya başqa alətə çıxarmaq mümkün deyildi.

İKİ FORMAT, İKİ MƏQSƏD:

    CSV   — TƏMİZ cədvəl (başlıq + sətirlər), şərh sətri YOXDUR.
            Excel və `pandas.read_csv` onu heç bir parametr olmadan açır.
    JSON  — eyni məlumat + METADATA (model adı, OOIP/OGIP, addım sayı,
            yığılma vəziyyəti). Proqramla emal üçün.

Metadata CSV-yə QƏSDƏN salınmır: `#` ilə başlayan şərh sətirləri ciddi
CSV oxuyucularını sındırır və gedər-gələr testini mənasızlaşdırır.
"""

from __future__ import annotations

import csv
import json
import math
from typing import Dict, List, Optional, Sequence

from ..simulation.results import SimulationResult

#: Excel Windows-da UTF-8-i BOM olmadan tanımır və "ə/ş/ğ" hərfləri
#: korlanır. `utf-8-sig` BOM əlavə edir; `csv`/`pandas` oxuyanda BOM-u
#: özləri atır, yəni gedər-gələr pozulmur.
ENCODING = "utf-8-sig"

#: Sütun adları — vahid MÜTLƏQ başlıqdadır (bax `UNITS.md`).
#:
#: Vahid KVADRAT MÖTƏRİZƏDƏDİR, vergüllə deyil: vergül CSV-nin öz
#: ayırıcısıdır, başlığa salınsa hər başlıq dırnağa düşür və
#: `awk`/`cut` kimi sadə alətlər sınır.
SERIES_COLUMNS: Sequence[tuple] = (
    ("time", "t [gün]"),
    ("oil_rate", "q_neft [m³/gün]"),
    ("water_rate", "q_su [m³/gün]"),
    ("water_injection_rate", "q_vurulan_su [m³/gün]"),
    ("gas_rate", "q_qaz [m³/gün]"),
    ("cumulative_oil", "kum_neft [m³]"),
    ("cumulative_water", "kum_su [m³]"),
    ("cumulative_gas", "kum_qaz [m³]"),
    ("water_cut", "su_kəsri [%]"),
    ("average_pressure", "orta_P [bar]"),
    ("recovery_factor", "RF [%]"),
    ("gas_oil_ratio", "GOR [sm³/sm³]"),
)

#: Quyu üzrə sütunlar — `SimulationResult`-dakı lüğət adı → başlıq şəkli.
WELL_COLUMNS: Sequence[tuple] = (
    ("well_oil_rate", "{name}: q_neft [m³/gün]"),
    ("well_water_rate", "{name}: q_su [m³/gün]"),
    ("well_gas_rate", "{name}: q_qaz [m³/gün]"),
    ("well_bhp", "{name}: BHP [bar]"),
    ("well_thp", "{name}: THP [bar]"),
)


def _is_blank(value) -> bool:
    """`nan` və `None` BOŞ xana kimi yazılır — SIFIR kimi YOX.

    `well_thp`-də quyunun səthə axa bilmədiyi addımlar `nan`-dır
    (bax `simulation/wellbore/traverse.py` — bu, qəsdən belədir).
    Orada `0` yazmaq YANLIŞ olardı: 0 bar real ölçmə kimi oxunar və
    istifadəçini aldadar.
    """
    if value is None:
        return True
    try:
        return math.isnan(float(value))
    except (TypeError, ValueError):
        return True


def _cell(value) -> str:
    return "" if _is_blank(value) else repr(float(value))


def _well_names(result: SimulationResult) -> List[str]:
    """Nəticədə HƏR HANSI sırası olan quyular, əlifba sırası ilə."""
    names = set()
    for attribute, _ in WELL_COLUMNS:
        names.update(getattr(result, attribute, {}) or {})
    return sorted(names)


def _columns(result: SimulationResult) -> List[tuple]:
    """`(başlıq, dəyərlər)` cütləri — BOŞ sıralar sütun YARATMIR.

    İki fazalı nəticədə `gas_rate`/`cumulative_gas`/`gas_oil_ratio`
    boş siyahıdır (yalnız üç fazalı mühərrik onları doldurur). Boş
    sütun yazmaq əvəzinə onu BURAXIRIQ — beləcə faylın özü "burada
    qaz fazası yoxdur" deyir.

    Eyni məntiq quyulara da aiddir: `well_bhp`/`well_thp` yalnız BHP
    rejimli, lüləsi verilmiş istismarçılarda mövcuddur.
    """
    series = result.series
    columns: List[tuple] = []

    for attribute, header in SERIES_COLUMNS:
        values = getattr(series, attribute, None) or []
        if values:
            columns.append((header, list(values)))

    for name in _well_names(result):
        for attribute, template in WELL_COLUMNS:
            values = (getattr(result, attribute, {}) or {}).get(name, [])
            if values:
                columns.append((template.format(name=name), list(values)))
    return columns


def _row_count(columns: Sequence[tuple]) -> int:
    return max((len(values) for _, values in columns), default=0)


# ═══════════════════════════════ CSV ══════════════════════════════════

def write_csv(result: SimulationResult, path: str) -> str:
    """Nəticəni TƏMİZ CSV cədvəli kimi yazır. Qaytarır: `path`.

    Sətir sayı ən uzun sütuna görə götürülür; qısa sütunlarda qalan
    xanalar BOŞ qalır (uydurma dəyərlə doldurulmur).
    """
    columns = _columns(result)
    with open(path, "w", encoding=ENCODING, newline="") as handle:
        writer = csv.writer(handle)
        if not columns:
            writer.writerow(["(nəticə boşdur)"])
            return path
        writer.writerow([header for header, _ in columns])
        for index in range(_row_count(columns)):
            writer.writerow([
                _cell(values[index]) if index < len(values) else ""
                for _, values in columns])
    return path


def read_csv(path: str) -> Dict[str, List[float]]:
    """`write_csv`-nin TƏRSİ — sütun adı → dəyərlər.

    Boş xana `nan` kimi qayıdır, yəni gedər-gələr dövrü `nan`-ı
    saxlayır. Bu funksiya testlərin işini asanlaşdırmaq üçün deyil,
    formatın MÜQAVİLƏSİNİ açıq etmək üçündür: ixrac olunan fayl
    yenidən oxunanda eyni ədədləri verməlidir.
    """
    out: Dict[str, List[float]] = {}
    with open(path, "r", encoding=ENCODING, newline="") as handle:
        reader = csv.DictReader(handle)
        for field in reader.fieldnames or []:
            out[field] = []
        for row in reader:
            for field in reader.fieldnames or []:
                text = (row.get(field) or "").strip()
                out[field].append(float("nan") if text == "" else float(text))
    return out


# ═══════════════════════════════ JSON ═════════════════════════════════

def _metadata(result: SimulationResult) -> dict:
    return {
        "model_name": result.model_name,
        "grid_shape": list(result.grid_shape) if result.grid_shape else [],
        "steps": result.steps,
        "converged": bool(result.converged),
        "message": result.message,
        "ooip": float(result.ooip),
        "ogip": float(result.ogip),
    }


def _clean(values: Sequence) -> List[Optional[float]]:
    """`nan` → `None`, çünki JSON-da `NaN` standart deyil.

    `json.dump` defolt olaraq `NaN` yazır, lakin bu, JSON
    spesifikasiyasına ziddir və bir çox oxuyucu onu qəbul etmir.
    `null` isə "dəyər yoxdur" mənasını dəqiq verir.
    """
    return [None if _is_blank(v) else float(v) for v in values]


def write_json(result: SimulationResult, path: str) -> str:
    """Nəticəni metadata ilə birlikdə JSON kimi yazır."""
    series = result.series
    payload = {
        "metadata": _metadata(result),
        "series": {attribute: _clean(getattr(series, attribute, []) or [])
                   for attribute, _ in SERIES_COLUMNS
                   if getattr(series, attribute, None)},
        "wells": {},
    }
    for name in _well_names(result):
        well: Dict[str, List[Optional[float]]] = {}
        for attribute, _ in WELL_COLUMNS:
            values = (getattr(result, attribute, {}) or {}).get(name, [])
            if values:
                well[attribute] = _clean(values)
        if well:
            payload["wells"][name] = well

    with open(path, "w", encoding="utf-8", newline="") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2,
                  allow_nan=False)
    return path


def write(result: SimulationResult, path: str) -> str:
    """Uzantıya görə format seçir (`.json` → JSON, qalanı → CSV)."""
    return (write_json(result, path) if path.lower().endswith(".json")
            else write_csv(result, path))
