"""Quyu tərifi — YALNIZ məlumat. Peaceman hesablaması burada deyil."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional

from .tubing import TubingGeometry
from .validation import validate_pressure, validate_well_rate


class WellType(Enum):
    PRODUCER = "PROD"
    INJECTOR = "INJ"


class ControlMode(Enum):
    BHP = "BHP"
    RATE = "RATE"
    #: Quyu başı təzyiqi, bar (B4-B). Mühərrik onu HƏR ADDIMDA tərs
    #: traversdən BHP-yə çevirir — bax `simulation/wellbore/thp_control.py`.
    #: Lülə həndəsəsi (`Well.tubing`) tələb edir; yalnız istismarçılar.
    THP = "THP"


class Phase(Enum):
    WATER = "WATER"
    OIL = "OIL"
    GAS = "GAS"


@dataclass
class Perforation:
    """Bir hüceyrədə açılmış interval."""
    i: int
    j: int
    k: int = 0
    open: bool = True
    skin: float = 0.0
    #: Perforasiyanın hüceyrəni HANSI ox boyunca deldiyi (Eclipse
    #: `COMPDAT` "DIR" sütunu ilə eyni hərflər): "Z"/"K" şaquli (defolt,
    #: mövcud bütün modellər), "X"/"I" və "Y"/"J" üfüqi tamamlama.
    #: Peaceman hesablaması bunu hüceyrənin YERLİ oxuna çevirir — bax
    #: `domain/geometry.py::completion_axis` və
    #: `CornerPointGeometry.wellblock_geometry`.
    direction: str = "Z"


@dataclass
class WellControl:
    mode: ControlMode = ControlMode.BHP
    target: float = 150.0
    injected_phase: Phase = Phase.WATER

    def validate(self) -> List[str]:
        """`target` mənası `mode`-dan asılıdır: BHP üçün mütləq təzyiq
        (bar), RATE üçün HƏMİŞƏ müsbət debit böyüklüyü (m³/gün) — bax
        `simulation/implicit/standard_well.py:_signed_rate_target`."""
        if self.mode is ControlMode.BHP:
            return validate_pressure([self.target], label="BHP hədəfi").errors
        if self.mode is ControlMode.THP:
            return validate_pressure([self.target], label="THP hədəfi").errors
        return validate_well_rate(self.target, label="debit hədəfi").errors

    def validate_warnings(self) -> List[str]:
        if self.mode is ControlMode.BHP:
            return validate_pressure([self.target], label="BHP hədəfi").warnings
        if self.mode is ControlMode.THP:
            return validate_pressure([self.target], label="THP hədəfi").warnings
        return validate_well_rate(self.target, label="debit hədəfi").warnings


@dataclass
class Well:
    name: str
    well_type: WellType = WellType.PRODUCER
    control: WellControl = field(default_factory=WellControl)
    perforations: List[Perforation] = field(default_factory=list)
    radius: float = 0.1
    active: bool = True
    # Metrlə perforasiya intervalı — YALNIZ geologiya cədvəli ilə əlaqəli
    # quyular üçün doldurulur. `perforations` (i,j,k) bunlardan HESABLANIR
    # və mühərrikin gördüyü yeganə şeydir; bu iki sahə yalnız metr girişinin
    # grid ölçüsü dəyişəndə itməməsi üçün saxlanılır (bax `depth_to_k`).
    perf_top: Optional[float] = None
    perf_bottom: Optional[float] = None

    #: Lülə borusu — THP hesabatı üçün (B4). `None` olanda quyu üçün
    #: THP hesablanmır; mühərrik onsuz da yalnız BHP ilə işləyir, ona
    #: görə bu sahə HEÇ NƏYİ məcbur etmir (geriyə tam uyğundur).
    tubing: Optional[TubingGeometry] = None

    @classmethod
    def vertical(cls, name, i, j, well_type=WellType.PRODUCER,
                 mode=ControlMode.BHP, target=150.0, radius=0.1,
                 skin=0.0, k=0) -> "Well":
        return cls(name=name, well_type=well_type,
                   control=WellControl(mode, target),
                   perforations=[Perforation(i, j, k, True, skin)],
                   radius=radius)

    @property
    def is_injector(self) -> bool:
        return self.well_type is WellType.INJECTOR

    def open_perforations(self) -> List[Perforation]:
        return [p for p in self.perforations if p.open]
