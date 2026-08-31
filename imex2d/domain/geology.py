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
from typing import Optional


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
