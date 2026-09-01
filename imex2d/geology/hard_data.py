"""Sərt (hard) fasiya datasının hüceyrə xəritələnməsi + ziddiyyət aşkarlanması.

Niyə lazımdır (Phase 4.1 §4): iki quyu nümunəsi eyni grid hüceyrəsinə
düşüb FƏRQLİ fasiya kodu bildirirsə, bu, SƏSSİZCƏ (məs. "sonuncu qazanır"
massiv təyinatı ilə) HƏLL EDİLMƏMƏLİDİR — istifadəçi bunu BİLMƏLİDİR və
ya AÇIQ strategiya seçməlidir. Bu modul `domain/geometry.xy_to_ij`/
`depth_to_k`-nı (mövcud, TƏKRARLANMAYAN) istifadə edərək hər nümunənin
hansı `(i,j,k)` hüceyrəsinə düşdüyünü tapır və ziddiyyətləri toplayır.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..domain.geometry import CellGeometry, depth_to_k, xy_to_ij
from ..domain.well_data import WellSample

Cell = Tuple[int, int, int]

#: Qəbul edilən ziddiyyət-həll strategiyaları. "average" YALNIZ kəsilməz
#: rejimdə (`tolerance` verilib) dəstəklənir (bax `resolve_hard_data`).
CONFLICT_STRATEGIES = ("raise", "majority", "keep_first", "keep_last", "average")


def find_exact_matches(points: np.ndarray, targets: np.ndarray, tolerance: float) -> np.ndarray:
    """Hər hədəf üçün TAM (tolerantlıq daxilində) üst-üstə düşən sərt
    data nöqtəsinin indeksini qaytarır, tapılmasa -1.

    Phase 4 (SIS) və Phase 5 (SGS) ARASINDA PAYLAŞILAN, xassə növündən
    ASILI OLMAYAN sadə həndəsə — hər ikisi bunu İDXAL edir, TƏKRAR
    YAZMIR (əvvəllər `geology/facies.py`-də idi, bura köçürülüb).
    """
    points = np.asarray(points, float)
    targets = np.asarray(targets, float)
    result = np.full(targets.shape[0], -1, dtype=int)
    if points.shape[0] == 0:
        return result
    ndim = min(points.shape[1], targets.shape[1])
    for row in range(targets.shape[0]):
        diff = points[:, :ndim] - targets[row, :ndim]
        distances = np.sqrt(np.sum(diff * diff, axis=1))
        closest = int(np.argmin(distances))
        if distances[closest] <= tolerance:
            result[row] = closest
    return result


class HardDataConflictError(ValueError):
    """`on_conflict="raise"` (defolt) və ya həll mümkün olmayanda (məs.
    bərabər səs) atılır — heç bir ziddiyyət SƏSSİZ qalmır."""


@dataclass(frozen=True)
class HardDataConflict:
    cell: Cell
    wells: Tuple[str, ...]
    codes: Tuple[float, ...]   # kateqorik rejimdə tam ədəd, kəsilməzdə float

    def __str__(self) -> str:
        pairs = ", ".join(f"{w}={c}" for w, c in zip(self.wells, self.codes))
        return f"hüceyrə {self.cell}: {pairs}"


def map_samples_to_cells(samples: Sequence[WellSample], source: str, grid,
                         geometry: CellGeometry
                         ) -> Tuple[Dict[Cell, List[int]], List[WellSample]]:
    """`source` xassəsini daşıyan nümunələri `(i,j,k)` hüceyrələrinə xəritələyir.

    Qaytarır: `(cell -> [samples_for indeksləri], samples_for)`.
    `sample.layer` verilibsə birbaşa K kimi işlədilir (artıq 0-based,
    bax `well_data_io.py`); yoxdursa `sample.depth` + `depth_to_k` ilə
    tapılır; heç biri yoxdursa (nümunə "bütün laylara aiddir") həmin
    nümunə xəritələnmir — ziddiyyət yoxlamasından KƏNARDA qalır (onun
    K-sı qəsdən qeyri-müəyyəndir, süni K UYDURULMUR).
    """
    samples_for = [s for s in samples if source in s.values]
    mapping: Dict[Cell, List[int]] = {}
    for idx, sample in enumerate(samples_for):
        i, j = xy_to_ij(sample.x, sample.y, geometry)
        if sample.layer is not None:
            k = sample.layer
        elif sample.depth is not None:
            k = depth_to_k(sample.x, sample.y, sample.depth, geometry)
        else:
            k = None
        if k is None:
            continue
        mapping.setdefault((i, j, int(k)), []).append(idx)
    return mapping, samples_for


def detect_hard_data_conflicts(samples: Sequence[WellSample], source: str, grid,
                               geometry: CellGeometry,
                               tolerance: Optional[float] = None) -> List[HardDataConflict]:
    """Eyni hüceyrəyə düşən, AMMA UYUŞMAYAN dəyər bildirən nümunə
    qruplarını tapır. Eyni dəyəri bildirən təkrarlar (duplikat, ziddiyyət
    DEYİL) siyahıya DAXİL EDİLMİR.

    `tolerance=None` (defolt) — KATEQORİK rejim: dəyərlər tam ədədə
    çevrilir, İSTƏNİLƏN fərq ziddiyyətdir (Phase 4.1, FACIES üçün).
    `tolerance=<float>` — KƏSİLMƏZ rejim (Phase 5, PORO/PERMX üçün):
    dəyərlər RAW float saxlanılır, YALNIZ `max-min > tolerance` olanda
    ziddiyyət sayılır (kiçik ölçmə fərqi ziddiyyət DEYİL).
    """
    mapping, samples_for = map_samples_to_cells(samples, source, grid, geometry)
    conflicts = []
    for cell, indices in mapping.items():
        if len(indices) < 2:
            continue
        raw = [float(samples_for[i].values[source]) for i in indices]
        if tolerance is None:
            codes = tuple(int(v) for v in raw)
            is_conflict = len(set(codes)) > 1
        else:
            codes = tuple(raw)
            is_conflict = (max(raw) - min(raw)) > tolerance
        if is_conflict:
            wells = tuple(samples_for[i].well for i in indices)
            conflicts.append(HardDataConflict(cell=cell, wells=wells, codes=codes))
    return conflicts


def resolve_hard_data(samples: Sequence[WellSample], source: str, grid,
                      geometry: CellGeometry, on_conflict: str = "raise",
                      tolerance: Optional[float] = None) -> List[WellSample]:
    """Ziddiyyəti HƏLL EDİLMİŞ nümunə alt-çoxluğunu qaytarır (X/Y/dərinlik
    HƏLƏ orijinal `WellSample`-dədir — çağıran (`geology_service.py`)
    3D (X,Y,Z) nöqtə massivini ÖZ mövcud `_sample_depth` məntiqi ilə
    qurur, bu funksiya YALNIZ "hansı nümunələr saxlanılır" qərarını verir,
    dərinlik-qurma MƏNTİQİNİ TƏKRARLAMIR).

    `tolerance` — bax `detect_hard_data_conflicts` (None=kateqorik,
    float=kəsilməz). `on_conflict`:
        "raise"      (defolt) — HƏR ziddiyyətdə `HardDataConflictError`.
        "majority"   — (YALNIZ kateqorik) ən çox təkrarlanan kodu
                       saxlayır; SƏS BƏRABƏRDİRSƏ yenə atılır.
        "keep_first" — nümunələr siyahısındakı İLK görünüşü saxlayır.
        "keep_last"  — SONUNCU görünüşü saxlayır.
        "average"    — (YALNIZ kəsilməz) ziddiyyətli dəyərlərin ORTASI
                       yeni (sintetik) nümunə kimi işlədilir.
    Heç biri SƏSSİZCƏ "sonuncu qazanır" DEYİL — seçim HƏMİŞƏ AÇIQ və
    sənədləşdirilmiş bir strategiyanın NƏTİCƏSİDİR.
    """
    if on_conflict not in CONFLICT_STRATEGIES:
        raise ValueError(f"Naməlum on_conflict: {on_conflict!r}. Dəstəklənən: {CONFLICT_STRATEGIES}")
    if on_conflict == "majority" and tolerance is not None:
        raise ValueError("'majority' yalnız kateqorik rejimdə (tolerance=None) mənalıdır.")
    if on_conflict == "average" and tolerance is None:
        raise ValueError("'average' yalnız kəsilməz rejimdə (tolerance verilməlidir) mənalıdır.")

    mapping, samples_for = map_samples_to_cells(samples, source, grid, geometry)
    conflicts = detect_hard_data_conflicts(samples, source, grid, geometry, tolerance)
    if conflicts and on_conflict == "raise":
        detail = "; ".join(str(c) for c in conflicts)
        raise HardDataConflictError(
            f"'{source}': {len(conflicts)} hüceyrədə ziddiyyətli sərt data (fərqli quyular "
            f"uyuşmayan dəyər bildirir): {detail}. Açıq strategiya seçin "
            f"(`on_conflict`={CONFLICT_STRATEGIES[1:]}) və ya məlumatı düzəldin.")

    result: List[WellSample] = []
    conflicted_cells = {c.cell for c in conflicts}
    for cell, indices in mapping.items():
        if cell not in conflicted_cells:
            result.extend(samples_for[i] for i in indices)
            continue
        if on_conflict == "keep_first":
            result.append(samples_for[indices[0]])
        elif on_conflict == "keep_last":
            result.append(samples_for[indices[-1]])
        elif on_conflict == "average":
            raw = [float(samples_for[i].values[source]) for i in indices]
            base = samples_for[indices[0]]
            averaged = dataclasses.replace(
                base, values={**base.values, source: float(np.mean(raw))})
            result.append(averaged)
        elif on_conflict == "majority":
            codes = [int(samples_for[i].values[source]) for i in indices]
            unique, counts = np.unique(codes, return_counts=True)
            top = counts.max()
            winners = unique[counts == top]
            if winners.size > 1:
                raise HardDataConflictError(
                    f"'{source}': hüceyrə {cell} üçün səs BƏRABƏRDİR ({dict(zip(unique, counts))}) "
                    "— 'majority' bunu həll edə bilmir, `on_conflict='keep_first'/'keep_last'` "
                    "seçin və ya məlumatı düzəldin.")
            winner_code = int(winners[0])
            result.append(next(samples_for[i] for i in indices
                               if int(samples_for[i].values[source]) == winner_code))

    # xəritələnə bilməyən (K qeyri-müəyyən) nümunələri də DAXİL et —
    # onlar heç bir hüceyrə ilə ziddiyyətə düşə bilməzlər (fərqli laylar üçün keçərlidirlər)
    mapped_indices = {i for indices in mapping.values() for i in indices}
    result.extend(samples_for[i] for i in range(len(samples_for)) if i not in mapped_indices)
    return result
