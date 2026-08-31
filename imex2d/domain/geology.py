"""Quyu cədvəli — geologiya bölməsinin yeganə giriş yolu.

Bu, istifadəçinin proqram daxilində birbaşa redaktə etdiyi quyu siyahısıdır
(əvvəlki CSV yükləməsinin əvəzidir). `None` və `0.0` fərqlidir: boş xana
`None`-dır, ölçülmüş sıfır deyil — Kriging/IDW `None` xassələri həmin quyuda
sadəcə nəzərə almır, `0.0`-ı isə real ölçü kimi qəbul edib nəticəni korlayar.

`top`/`bottom` PERFORASİYA DEYİL — layın bu nöqtədəki struktur dərinliyidir
(üst səthin interpolyasiyası üçün). Perforasiya rezervuar modelindədir
(bax `imex2d/domain/wells.py`).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional

from .geometry import CellGeometry, xy_to_ij


@dataclass
class GeologicalWell:
    name: str
    in_model: bool = True
    x: float = 0.0
    y: float = 0.0
    top: Optional[float] = None
    bottom: Optional[float] = None
    porosity: Optional[float] = None
    permeability: Optional[float] = None
    water_saturation: Optional[float] = None
    note: str = ""

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "in_model": self.in_model,
            "x": self.x,
            "y": self.y,
            "top": self.top,
            "bottom": self.bottom,
            "porosity": self.porosity,
            "permeability": self.permeability,
            "water_saturation": self.water_saturation,
            "note": self.note,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "GeologicalWell":
        return cls(
            name=data["name"],
            in_model=bool(data.get("in_model", True)),
            x=float(data.get("x", 0.0)),
            y=float(data.get("y", 0.0)),
            top=data.get("top"),
            bottom=data.get("bottom"),
            porosity=data.get("porosity"),
            permeability=data.get("permeability"),
            water_saturation=data.get("water_saturation"),
            note=data.get("note", ""),
        )


def method_minimum(method: str) -> int:
    """Seçilmiş interpolyasiya üsulu üçün minimum quyu sayı.

    Tapşırıqda "Kriging ≥ 3, IDW ≥ 1" yazılıb, amma paylaşılan
    `WellDataset.validate()` (CSV importunda da işlədilir) HƏR xassə üçün
    minimum 2 nöqtə tələb edir (bax `test_single_point_dataset_is_rejected`,
    `imex2d/domain/well_data.py`). Bunu dəyişmək CSV importunun mövcud
    sınağını sındırardı, ona görə faktiki minimum 2-yə qaldırılıb — bax
    `ISH_HESABATI.md`.
    """
    if "Kriging" in method:
        return 3
    return 2


@dataclass
class ValidationIssue:
    """Bir yoxlama nəticəsi. `error` interpolyasiyanı bloklayır."""
    level: str                    # "error" | "warning" | "info"
    message: str
    well: Optional[str] = None


_PROPERTY_LABELS = (("porosity", "φ"), ("permeability", "k"),
                    ("water_saturation", "Sw"))


def validate_wells(wells: List[GeologicalWell],
                   geometry: Optional[CellGeometry] = None,
                   method: str = "",
                   reservoir_well_names: Optional[Iterable[str]] = None
                   ) -> List[ValidationIssue]:
    """Geologiya cədvəlini yoxlayır. UI-dən asılı deyil, test yazıla bilər.

    `geometry` verilməyibsə (grid hələ qurulmayıb) sərhəd/hüceyrə
    yoxlamaları keçilir — yalnız cədvəlin öz daxili tutarlılığı yoxlanır.
    """
    issues: List[ValidationIssue] = []

    # ── ad: boş / təkrar ──────────────────────────────────────────────
    groups: dict = {}
    for well in wells:
        name = (well.name or "").strip()
        if not name:
            issues.append(ValidationIssue("error", "Quyu adı boşdur.", well.name))
            continue
        groups.setdefault(name, []).append(well)
    for name, group in groups.items():
        if len(group) > 1:
            issues.append(ValidationIssue(
                "error", f"Quyu adı təkrarlanır: '{name}' ({len(group)} dəfə).", name))

    # ── xassə diapazonları ───────────────────────────────────────────
    for well in wells:
        if well.porosity is not None and not (0.0 <= well.porosity <= 1.0):
            issues.append(ValidationIssue(
                "error", f"'{well.name}': φ = {well.porosity:g} [0, 1] aralığından kənardadır.",
                well.name))
        if well.water_saturation is not None and not (0.0 <= well.water_saturation <= 1.0):
            issues.append(ValidationIssue(
                "error", f"'{well.name}': Sw = {well.water_saturation:g} [0, 1] aralığından kənardadır.",
                well.name))
        if well.permeability is not None and well.permeability <= 0.0:
            issues.append(ValidationIssue(
                "error", f"'{well.name}': k = {well.permeability:g} ≤ 0.", well.name))
        if (well.top is not None and well.bottom is not None
                and well.bottom <= well.top):
            issues.append(ValidationIssue(
                "error", f"'{well.name}': lay altı ({well.bottom:g}) ≤ lay üstü ({well.top:g}).",
                well.name))
        missing = [label for attr, label in _PROPERTY_LABELS
                  if getattr(well, attr) is None]
        if well.top is None:
            missing.append("lay üstü")
        if well.bottom is None:
            missing.append("lay altı")
        if missing:
            issues.append(ValidationIssue(
                "info", f"'{well.name}': boş xanalar ({', '.join(missing)}) — "
                        "bu xassələr həmin quyunu nəzərə almadan hesablanacaq.",
                well.name))

    # ── üsul üçün quyu sayı (xassə-üzrə, BLOKLAMIR — bax ISH_HESABATI.md) ──
    if method:
        required = method_minimum(method)
        any_sufficient = False
        for attr, label in _PROPERTY_LABELS + (("top", "lay üstü"), ("bottom", "lay altı")):
            count = sum(1 for w in wells if getattr(w, attr) is not None)
            if count >= required:
                any_sufficient = True
            elif 0 < count < required:
                issues.append(ValidationIssue(
                    "warning",
                    f"'{label}' üçün quyu sayı azdır ({count}/{required}, üsul: {method}) "
                    "— bu xassə interpolyasiya olunmayacaq.", None))
        if wells and not any_sufficient:
            issues.append(ValidationIssue(
                "error", "Heç bir xassə üçün seçilmiş üsula kifayət qədər quyu yoxdur "
                        f"(tələb olunan: {required}).", None))

    # ── grid-lə bağlı yoxlamalar (yalnız grid qurulubsa) ────────────────
    if geometry is not None:
        x_max, y_max = geometry.areal_extent()
        cells: dict = {}
        for well in wells:
            in_bounds = 0.0 <= well.x <= x_max and 0.0 <= well.y <= y_max
            if not in_bounds:
                issues.append(ValidationIssue(
                    "error",
                    f"'{well.name}': (X={well.x:g}, Y={well.y:g}) grid sərhədlərindən "
                    f"kənardadır (0–{x_max:g}, 0–{y_max:g}).", well.name))
                continue
            i, j = xy_to_ij(well.x, well.y, geometry)
            cells.setdefault((i, j), []).append(well.name)

            near_x = min(well.x % geometry.dx, geometry.dx - (well.x % geometry.dx))
            near_y = min(well.y % geometry.dy, geometry.dy - (well.y % geometry.dy))
            if near_x < 1e-6 or near_y < 1e-6:
                issues.append(ValidationIssue(
                    "warning",
                    f"'{well.name}': hüceyrə kənarına düşür — Peaceman quyu indeksi "
                    "etibarsızlaşa bilər.", well.name))
        for (i, j), names in cells.items():
            if len(names) > 1:
                issues.append(ValidationIssue(
                    "warning",
                    f"İki və ya daha çox quyu eyni hüceyrədədir ({i}, {j}): "
                    f"{', '.join(names)}.", None))

    # ── 7-ci bölmə ilə əlaqə ─────────────────────────────────────────
    if reservoir_well_names is not None:
        reservoir_names = set(reservoir_well_names)
        geology_names = {w.name for w in wells}
        for well in wells:
            if well.in_model and well.name not in reservoir_names:
                issues.append(ValidationIssue(
                    "warning",
                    f"'{well.name}': modeldə işarələnib, amma 7-ci bölmədə rejim "
                    "təyin olunmayıb.", well.name))
        for name in reservoir_names - geology_names:
            issues.append(ValidationIssue(
                "warning",
                f"'{name}': 7-ci bölmədə var, amma geologiya cədvəlində yoxdur "
                "(sahibsiz quyu).", name))

    return issues
