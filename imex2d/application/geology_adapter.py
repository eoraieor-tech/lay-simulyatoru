"""`GeologicalWell` cədvəlini `geology_service.py`-nin gözlədiyi
`WellDataset`-ə çevirir.

`geology_service.py` TOXUNULMUR (bax `ISH_HESABATI.md`, bölmə 9) — bu
modul yalnız girişi CSV oxuyucunun yerinə keçirir. Xassə-üzrə davranış
("bir quyuda dəyər yoxdursa yalnız o xassə düşür, sətir yox") elə burada
əldə olunur: `WellSample.values` yalnız MÖVCUD olan xassələri daşıyır,
`WellDataset.points()` isə həmin xassəni istəməyən sətirləri özü süzür.

Seçilmiş üsul üçün quyu sayı kifayət etmirsə (bax
`imex2d.domain.geology.method_minimum`), həmin xassə dataset-dən
TAMAMİLƏ ÇIXARILIR ki, `geology_service.py` onu heç görməsin və
digər xassələr toxunulmadan hesablansın (`WellBasedGeologicalModelBuilder`
`dataset.property_names()` üzərində dövr edir).
"""

from __future__ import annotations

from typing import Dict, List, Tuple

from ..domain.geology import GeologicalWell, method_minimum
from ..domain.well_data import WellDataset, WellSample

# GeologicalWell atributu -> WellDataset xassə açarı (DEFAULT_RULES ilə uyğun)
PROPERTY_MAP = {
    "porosity": "PORO",
    "permeability": "PERMX",
    "water_saturation": "SW",
    "top": "TOP",
    "bottom": "BOTTOM",
}


def wells_to_dataset(wells: List[GeologicalWell],
                     method: str = "") -> Tuple[WellDataset, Dict[str, str]]:
    """`(dataset, skipped)` qaytarır.

    `skipped`: xassə adı → səbəb mesajı. Bu xassələr `dataset`-ə
    daxil edilməyib (quyu sayı kifayət etmədiyi üçün).
    """
    required = method_minimum(method)
    counts: Dict[str, int] = {}
    for well in wells:
        for attr, target in PROPERTY_MAP.items():
            if getattr(well, attr) is not None:
                counts[target] = counts.get(target, 0) + 1

    skipped = {
        target: (f"'{target}' üçün {count}/{required} quyu var "
                 f"(üsul: {method or 'seçilməyib'}) — bu xassə interpolyasiya olunmadı.")
        for target, count in counts.items() if count < required
    }

    samples: List[WellSample] = []
    for well in wells:
        values = {}
        for attr, target in PROPERTY_MAP.items():
            value = getattr(well, attr)
            if value is not None and target not in skipped:
                values[target] = value
        if values:
            samples.append(WellSample(well=well.name, x=well.x, y=well.y,
                                      values=values))
    return WellDataset(samples=samples, source="wells"), skipped
