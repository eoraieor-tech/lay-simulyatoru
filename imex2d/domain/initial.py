"""İlkin şərtlər.

İki müstəqil bayraq ilkin təzyiq və doyumluluğun HARADAN gəldiyini
müəyyən edir. Bayraqlar bir-birini əvəz etmir — biri TƏZYİQƏ, digəri
DOYUMLULUĞA aiddir:

    use_equilibration    — təzyiq: hidrostatik profil (bax
                           `simulation/initialization/equilibrium.py`).
    use_saturation_map   — doyumluluq: `property_maps["SW"]` xəritəsi
                           (geologiya cədvəlindən interpolyasiya olunub,
                           bax `simulation/initialization/saturation_map.py`).

PRİORİTET CƏDVƏLİ

    use_equilibration  use_saturation_map  Təzyiq        Doyumluluq
    ─────────────────  ──────────────────  ────────────  ─────────────────────
    False              False               datum (bərabər)  `water_saturation` skalyar
    False              True                datum (bərabər)  `SW` xəritəsi
    True               False               hidrostatik      kontakt/keçid zonası
    True               True                hidrostatik      `SW` xəritəsi

Son sətir: təzyiq fiziki MODELDƏN (hidrostatika) gəlir, doyumluluq isə
ÖLÇÜLMÜŞ məlumatdır — ölçmə modelin üstündədir, ona görə xəritə
kontakt/keçid zonası hesabını ƏVƏZ EDİR (jurnala açıq INFO yazılır).

Hər iki bayraq defolt `False`-dur — köhnə `.imx` faylları və mövcud
etalonlar əvvəlki yolu gedir.
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
    use_saturation_map: bool = False
