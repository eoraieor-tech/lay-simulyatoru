"""QUYU İNTERVALI → GRID K-TƏBƏQƏSİ uyğunlaşdırması VƏ lay-üzrə
məlumat mövcudluğunun HESABLANMASI.

Bu modul `domain/data_availability.py`-nin (SAF data modeli) ÜSTÜNDƏ
duran HESABLAMA qatıdır: `WellDataset` + `CellGeometry` → `Model
DataAvailability`.

İKİ ANLAYIŞ BURADA DA AYRIDIR (tapşırıq §7):

    `well_interval_layers()`  — quyunun FİZİKİ intervalının KƏSDİYİ
                                laylar ("effective layers"). Bu, YALNIZ
                                HƏNDƏSƏDİR.
    `compute_availability()`  — həmin xassənin HANSI layda HƏQİQƏTƏN
                                ölçüldüyü ("data availability").

Birincisi ikincisini AVTOMATİK VERMİR. `top=2000, bottom=2210` beş layı
kəssə də, PERMX yalnız L1–L3-də ölçülmüş ola bilər — bu fərq
`LayerDataPolicy` ilə AÇIQ idarə olunur, SÜKUTLA fərz edilmir.
"""

from __future__ import annotations

from enum import Enum
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np

from ..domain.data_availability import (DataStatus, ModelDataAvailability,
                                        PropertyAvailability)
from ..domain.geometry import CellGeometry, depth_to_k, interval_layers, xy_to_ij
from ..domain.well_data import WellDataset, WellSample


class LayerDataPolicy(str, Enum):
    """LAY ETİKETİ OLMAYAN (`layer is None`, `depth is None`) nümunə
    hansı laylara aid sayılsın.

    Bu, TAPŞIRIĞIN ƏSAS QƏRAR NÖQTƏSİDİR — hər üç variant AÇIQ seçimdir,
    heç biri "sükutla fərz edilən" deyil:

    `BROADCAST` — köhnə (geriyə-uyğun) davranış: etiketsiz nümunə BÜTÜN
        laylara aiddir. UI cədvəlindən gələn tək-dəyərli quyular indiyə
        qədər məhz belə işlənib (bax `tests/test_layer_aware_kriging_
        leak.py::test_ui_table_wells_broadcast_single_value_to_every_
        layer`). Lay-məlumatlı rejim SÖNDÜRÜLÜ olanda DEFOLTdur.

    `STRICT` — etiketsiz nümunə HEÇ BİR laya aid deyil: mövcudluq YALNIZ
        açıq `layer`/`depth` (və ya UI-dakı "Data layları" sütunu) ilə
        müəyyən olunur. Lay-məlumatlı rejimin DEFOLTudur (tapşırıq §10:
        "sistem səssiz fərziyyə yaratmamalıdır").

    `INTERVAL` — etiketsiz nümunə quyunun `top/bottom` intervalının
        KƏSDİYİ laylara aid sayılır. Bu, AÇIQ İSTİFADƏÇİ SEÇİMİDİR (və
        hesabatda xəbərdarlıq kimi görünür), çünki "interval = data"
        fərziyyəsi ELMİ CƏHƏTDƏN ZƏMANƏTLİ DEYİL.
    """

    BROADCAST = "broadcast"
    STRICT = "strict"
    INTERVAL = "interval"


def well_interval_layers(x: float, y: float, top: Optional[float],
                         bottom: Optional[float], geometry: CellGeometry) -> List[int]:
    """Quyu intervalının KƏSDİYİ K-təbəqələri (0-əsaslı, artan sıra).

    HƏNDƏSƏ `domain/geometry.interval_layers()`-dədir (SAF domain, kənar
    halların hamısı orada sənədləşib) — bu funksiya yalnız "verilməyib"
    (`None`) halını əlavə edir: `top`/`bottom` boşdursa FƏRZİYYƏ
    QURULMUR, boş siyahı qayıdır.
    """
    if top is None or bottom is None:
        return []
    return interval_layers(x, y, float(top), float(bottom), geometry)


def sample_layers(sample: WellSample, geometry: CellGeometry,
                  policy: LayerDataPolicy = LayerDataPolicy.STRICT,
                  interval: Optional[Sequence[int]] = None) -> List[int]:
    """Bir `WellSample` HANSI laylara aid sayılır.

    Sıra (ilk uyğun gələn qazanır):
      0. `sample.areal`                           → BÜTÜN laylar
         (struktur səth — TƏBİƏTƏN laydan asılı deyil, bax `WellSample.areal`)
      1. AÇIQ `sample.layer`                      → yalnız o lay
      2. ÖLÇÜLMÜŞ `sample.depth`                  → `depth_to_k`
      3. `policy` (bax `LayerDataPolicy`)

    Grid-dən kənar dərinlik `[]` qaytarır — SƏSSİZCƏ ən yaxın laya
    "sancılmır" (§26: silent corruption yoxdur); çağıran bunu
    xəbərdarlıq kimi göstərir.
    """
    nz = geometry.grid.nz
    if sample.areal and sample.layer is None:
        return list(range(nz))
    if sample.layer is not None:
        return [int(sample.layer)] if 0 <= int(sample.layer) < nz else []
    if sample.depth is not None:
        k = depth_to_k(sample.x, sample.y, float(sample.depth), geometry)
        return [] if k is None else [int(k)]
    if policy is LayerDataPolicy.BROADCAST:
        return list(range(nz))
    if policy is LayerDataPolicy.INTERVAL:
        return list(interval or [])
    return []


def compute_availability(dataset: WellDataset, geometry: CellGeometry,
                         policy: LayerDataPolicy = LayerDataPolicy.STRICT,
                         properties: Optional[Iterable[str]] = None,
                         intervals: Optional[Dict[str, Sequence[int]]] = None,
                         ) -> ModelDataAvailability:
    """`WellDataset` → XASSƏ-SPESİFİK lay mövcudluğu.

    `intervals` — quyu adı → həmin quyunun interval laylarının siyahısı
    (yalnız `policy=INTERVAL` üçün; `wells_to_dataset`/çağıran hazırlayır).

    Nəticədə HƏR xassənin ÖZ `PropertyAvailability`-si olur: hansı layda
    neçə sərt nöqtə var (`n_data`) və status `MEASURED`/`MISSING`. Bu,
    İNTERPOLYASİYADAN ƏVVƏLKİ (giriş) mənzərədir — interpolyasiya/
    completion sonra bu obyekti YENİLƏYİR (bax `geology_service`).
    """
    nz = geometry.grid.nz
    names = sorted(properties) if properties is not None else dataset.property_names()
    availability = ModelDataAvailability(nz=nz)
    counts: Dict[str, np.ndarray] = {name: np.zeros(nz, dtype=int) for name in names}

    for sample in dataset.samples:
        layers = sample_layers(sample, geometry, policy,
                               (intervals or {}).get(sample.well))
        if not layers:
            continue
        for name in names:
            if name in sample.values and np.isfinite(sample.values[name]):
                counts[name][layers] += 1

    for name in names:
        entry = PropertyAvailability(name=name, nz=nz)
        for k in range(nz):
            n = int(counts[name][k])
            entry.set(k, n_data=n,
                      status=DataStatus.MEASURED if n else DataStatus.MISSING)
        availability.properties[name] = entry
    return availability


def hard_data_cells(dataset: WellDataset, geometry: CellGeometry, source: str,
                    policy: LayerDataPolicy = LayerDataPolicy.STRICT,
                    intervals: Optional[Dict[str, Sequence[int]]] = None) -> np.ndarray:
    """`(ncell,)` bool — hansı hüceyrədə HƏQİQİ ölçmə var.

    `MEASURED` statusunu HÜCEYRƏ səviyyəsində vermək üçün lazımdır:
    interpolyasiya olunmuş layda quyu hüceyrəsinin özü ÖLÇÜLMÜŞDÜR
    (sərt-data honoring sayəsində dəyər ölçmə ilə üst-üstə düşür),
    qonşu hüceyrələr isə İNTERPOLYASİYADIR.
    """
    grid = geometry.grid
    mask = np.zeros(grid.ncell, dtype=bool)
    for sample in dataset.samples:
        if source not in sample.values or not np.isfinite(sample.values[source]):
            continue
        i, j = xy_to_ij(sample.x, sample.y, geometry)
        for k in sample_layers(sample, geometry, policy,
                               (intervals or {}).get(sample.well)):
            mask[np.ravel_multi_index((k, j, i), grid.shape)] = True
    return mask


def unassigned_samples(dataset: WellDataset, geometry: CellGeometry,
                       policy: LayerDataPolicy = LayerDataPolicy.STRICT,
                       intervals: Optional[Dict[str, Sequence[int]]] = None
                       ) -> List[WellSample]:
    """Heç bir laya aid edilə bilməyən nümunələr — çağıran bunları AÇIQ
    xəbərdarlıq kimi göstərməlidir (SƏSSİZCƏ atılmır)."""
    return [s for s in dataset.samples
            if not sample_layers(s, geometry, policy, (intervals or {}).get(s.well))]
