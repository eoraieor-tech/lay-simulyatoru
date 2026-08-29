"""İlkin şərtlər — PLACEHOLDER.

Hazırda yalnız bərabər paylanmış təzyiq və doyumluluq saxlanılır.
A3-dən sonra: `use_equilibration` və `oil_water_contact` verildikdə
EquilibriumInitializationProvider bu sinfi GİRİŞ məlumatı kimi oxuyur
və dərinlikdən asılı təzyiq/doyumluluq sahələri qurur. Bayraq
söndürülübsə, bərabər paylanmış ilkin şərtlər işlədilir (köhnə davranış).
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional


@dataclass
class InitialConditions:
    datum_depth: float = 0.0
    datum_pressure: float = 250.0
    water_saturation: float = 0.20
    oil_water_contact: Optional[float] = None
    equilibration_region: int = 1
    use_equilibration: bool = False
