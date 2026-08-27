"""Vahid sistemi. Sabitlər koda səpələnmək əvəzinə bir yerdə saxlanılır."""

from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class UnitSystem:
    """Darsi sabiti vahid sistemindən asılıdır — hardcode edilməməlidir."""
    name: str
    darcy_constant: float
    length: str
    permeability: str
    viscosity: str
    pressure: str
    time: str
    rate: str


METRIC = UnitSystem("METRIC", 0.008527, "m", "mD", "cP", "bar", "day", "m3/day")
FIELD = UnitSystem("FIELD", 1.127e-3, "ft", "mD", "cP", "psi", "day", "bbl/day")
DEFAULT_UNITS = METRIC
