"""Xassə növü reyestri — KATEQORİK vs KƏSİLMƏZ (Phase 4.1).

Niyə lazımdır: audit (bax `FACIES_INTEGRATION.md`) təsdiqlədi ki,
`geology_service.py`-nin mövcud iş axını HƏR sütunu (adı nə olursa olsun)
eyni kəsilməz Kriging/IDW yolundan keçirir — FACIES kimi kateqorik sütun
sükutla ədədi kimi interpolyasiya olunub 1.4 kimi mənasız aralıq dəyər
verə bilər. Bu reyestr həmin qərarı BİR yerdə, AÇIQ və GENİŞLƏNƏ BİLƏN
şəkildə tutur — sütun adı ilə səpələnmiş `if name == "FACIES"` yoxlaması
YOX.

`classify_property()` HEÇ VAXT sükutla "kateqorikdir" qərarı VERMİR —
naməlum ad defolt olaraq KƏSİLMƏZ sayılır (bu, mövcud davranışdır,
GERİYƏ UYĞUNLUĞU qoruyur) — YALNIZ AÇIQ reyestrdə olan və ya çağıranın
`overrides`-də bildirdiyi adlar KATEQORİK sayılır.
"""

from __future__ import annotations

from enum import Enum
from typing import Dict, Optional


class PropertyType(Enum):
    CONTINUOUS = "continuous"
    CATEGORICAL = "categorical"


#: Defolt reyestr — YALNIZ bu adlar (böyük hərflə) müəyyən bir tipə
#: bağlanır. Yeni kateqorik xassə əlavə etmək üçün BURAYA (və ya
#: `classify_property(overrides=...)`-ə) bir sətir əlavə etmək kifayətdir
#: — `geology_service.py`-nin özünə TOXUNMADAN.
DEFAULT_PROPERTY_TYPES: Dict[str, PropertyType] = {
    "PORO": PropertyType.CONTINUOUS,
    "PERMX": PropertyType.CONTINUOUS,
    "PERMY": PropertyType.CONTINUOUS,
    "PERMZ": PropertyType.CONTINUOUS,
    "NTG": PropertyType.CONTINUOUS,
    "SW": PropertyType.CONTINUOUS,
    "VSH": PropertyType.CONTINUOUS,
    "PRESSURE": PropertyType.CONTINUOUS,
    "FACIES": PropertyType.CATEGORICAL,
    "LITHOLOGY": PropertyType.CATEGORICAL,
    "ROCKTYPE": PropertyType.CATEGORICAL,
}


def classify_property(name: str, overrides: Optional[Dict[str, PropertyType]] = None
                      ) -> PropertyType:
    """`name` (böyük-kiçik hərfə həssas deyil) üçün xassə növünü qaytarır.

    `overrides` çağıranın öz layihəsinə məxsus əlavə/dəyişdirilmiş
    təyinatlarıdır (dəyişdirmədən DEFAULT reyestri toxunulmaz qalır).
    Naməlum ad → `CONTINUOUS` (mövcud, geriyə-uyğun davranış).
    """
    registry = dict(DEFAULT_PROPERTY_TYPES)
    if overrides:
        registry.update({key.upper(): value for key, value in overrides.items()})
    return registry.get(name.upper(), PropertyType.CONTINUOUS)


def is_categorical(name: str, overrides: Optional[Dict[str, PropertyType]] = None) -> bool:
    return classify_property(name, overrides) is PropertyType.CATEGORICAL
