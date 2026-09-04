"""`GeologicalWell` cədvəlini `geology_service.py`-nin gözlədiyi
`WellDataset`-ə çevirir.

Qeyd: `geology_service.py`-a "toxunulmazdır" qərarı (bax `ISH_HESABATI.md`,
bölmə 9) həmin fazanın (CSV → cədvəl keçidi) öz əhatəsi üçün idi. Layer
sızması düzəlişi (bax `ISH_HESABATI.md`, M1) üçün `_interpolate_volume`
dəyişdirildi — Kriging riyaziyyatına toxunmadan, yalnız boş laya hansı
nöqtələrin daxil ediləcəyi məntiqinə.

LAY-MƏLUMATLI REJİM (yeni). Bu modul indi İKİ rejimdə işləyir:

  · `LayerDataPolicy.BROADCAST` (DEFOLT, geriyə-uyğun) — `WellSample.layer`
    doldurulmur, tək dəyər bütün K-lara yayılır. Bu, kod bazasının
    indiyə qədərki DAVRANIŞIDIR (bax `tests/test_layer_aware_kriging_
    leak.py::test_ui_table_wells_broadcast_single_value_to_every_layer`)
    və `geometry` verilmədikdə YEGANƏ mümkün rejimdir.

  · `LayerDataPolicy.STRICT` / `INTERVAL` — hər quyu üçün `Geological
    Well.data_layers_text` ("Data layları" sütunu) oxunur və HƏR LAY
    ÜÇÜN AYRICA `WellSample` yaradılır (`layer=k`, `depth` = həmin
    hüceyrənin mərkəz dərinliyi). Beləliklə `dataset.is_layered()` DOĞRU
    olur və `geology_service._interpolate_volume` hər layı YALNIZ ÖZ
    məlumatı ilə hesablayır.

KRİTİK FƏRQ (tapşırıq §1): `top`/`bottom` (quyunun FİZİKİ intervalı) HEÇ
BİR REJİMDƏ lay indeksi kimi İŞLƏDİLMİR. `INTERVAL` siyasəti bunu YALNIZ
istifadəçi AÇIQ seçəndə edir və nəticə hesabatda xəbərdarlıq kimi görünür.
`TOP`/`BOTTOM` isə STRUKTUR SƏTHLƏRİDİR (petrofizika deyil) — `WellSample.
areal=True` ilə işarələnir, ona görə hər layda keçərlidir və lay-üzrə
mövcudluq yoxlamasını bloklamır.

Xassə-üzrə davranış ("bir quyuda dəyər yoxdursa yalnız o xassə düşür, sətir
yox") elə burada əldə olunur: `WellSample.values` yalnız MÖVCUD olan
xassələri daşıyır, `WellDataset.points()` isə həmin xassəni istəməyən
sətirləri özü süzür.

Seçilmiş üsul üçün quyu sayı kifayət etmirsə (bax
`imex2d.domain.geology.method_minimum`), həmin xassə dataset-dən
TAMAMİLƏ ÇIXARILIR ki, `geology_service.py` onu heç görməsin və
digər xassələr toxunulmadan hesablansın (`WellBasedGeologicalModelBuilder`
`dataset.property_names()` üzərində dövr edir).
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple



from ..domain.geology import GeologicalWell, method_minimum, well_effective_layers
from ..domain.geometry import CellGeometry, xy_to_ij
from ..domain.well_data import WellDataset, WellSample
from ..geology.layer_availability import LayerDataPolicy

# GeologicalWell atributu -> WellDataset xassə açarı (DEFAULT_RULES ilə uyğun)
PROPERTY_MAP = {
    "porosity": "PORO",
    "permeability": "PERMX",
    "water_saturation": "SW",
    "top": "TOP",
    "bottom": "BOTTOM",
}

#: STRUKTUR SƏTHLƏRİ — petrofiziki ölçmə DEYİL, ona görə lay-üzrə
#: mövcudluq anlayışına TABE DEYİL (bax `WellSample.areal`). Bu, "PORO
#: üçün xüsusi kod" tipli hack DEYİL: burada AYRILAN şey xassənin ADI
#: yox, ÖLÇMƏNİN TƏBİƏTİDİR (səth vs. həcm xassəsi).
AREAL_TARGETS = ("TOP", "BOTTOM")


def wells_to_dataset(wells: List[GeologicalWell],
                     method: str = "",
                     geometry: Optional[CellGeometry] = None,
                     policy: LayerDataPolicy = LayerDataPolicy.BROADCAST,
                     ) -> Tuple[WellDataset, Dict[str, str]]:
    """`(dataset, skipped)` qaytarır.

    `skipped`: xassə adı → səbəb mesajı. Bu xassələr `dataset`-ə
    daxil edilməyib (quyu sayı kifayət etmədiyi üçün).

    `geometry`/`policy` — lay-məlumatlı rejim üçün. `geometry=None` və ya
    `policy=BROADCAST` olanda davranış ƏVVƏLKİ İLƏ EYNİDİR (nə `layer`,
    nə `depth` doldurulur). Lay bəyanı ilə bağlı qeydlər `dataset.
    warnings`-ə yazılır — SƏSSİZ ATILMA YOXDUR.
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

    layer_aware = geometry is not None and policy is not LayerDataPolicy.BROADCAST
    if not layer_aware:
        return _broadcast_dataset(wells, skipped), skipped
    return _layered_dataset(wells, skipped, geometry, policy), skipped


# ─────────────────────────────────────────────── köhnə (laysız) yol
def _broadcast_dataset(wells: List[GeologicalWell],
                       skipped: Dict[str, str]) -> WellDataset:
    """ƏVVƏLKİ davranış — hər quyu üçün BİR nümunə, lay etiketi YOXDUR."""
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
    return WellDataset(samples=samples, source="wells")


# ──────────────────────────────────────────── lay-məlumatlı yol
def _layered_dataset(wells: List[GeologicalWell], skipped: Dict[str, str],
                     geometry: CellGeometry, policy: LayerDataPolicy) -> WellDataset:
    nz = geometry.grid.nz
    depths = geometry.cell_depths().reshape(geometry.grid.shape)   # (nz, ny, nx)
    samples: List[WellSample] = []
    warnings: List[str] = []

    for well in wells:
        i, j = xy_to_ij(well.x, well.y, geometry)
        try:
            default_layers, per_property = well.data_layer_sets(nz)
        except ValueError as error:
            raise ValueError(f"'{well.name}': data layları oxunmadı — {error}") from error

        effective = _effective_layers(well, geometry, warnings)
        per_layer: Dict[int, Dict[str, float]] = {}
        areal_values: Dict[str, float] = {}

        for attr, target in PROPERTY_MAP.items():
            value = getattr(well, attr)
            if value is None or target in skipped:
                continue
            if target in AREAL_TARGETS:
                areal_values[target] = value
                continue
            layers = _resolve_layers(well, target, default_layers, per_property,
                                     effective, policy, warnings)
            for k in layers:
                per_layer.setdefault(k, {})[target] = value

        for k, values in sorted(per_layer.items()):
            samples.append(WellSample(
                well=well.name, x=well.x, y=well.y, values=values, layer=k,
                depth=float(depths[k, j, i])))
        if areal_values:
            samples.append(WellSample(well=well.name, x=well.x, y=well.y,
                                      values=areal_values, areal=True))

    return WellDataset(samples=samples, source="wells", warnings=warnings)


def _effective_layers(well: GeologicalWell, geometry: CellGeometry,
                      warnings: List[str]) -> List[int]:
    try:
        return well_effective_layers(well, geometry)
    except ValueError as error:
        warnings.append(f"'{well.name}': interval oxunmadı — {error}")
        return []


def _resolve_layers(well: GeologicalWell, target: str,
                    default_layers: Optional[Sequence[int]],
                    per_property: Dict[str, Sequence[int]],
                    effective: Sequence[int], policy: LayerDataPolicy,
                    warnings: List[str]) -> List[int]:
    """Bu quyunun bu xassəsi HANSI laylara MƏLUMAT verir.

    Üstünlük sırası — hamısı AÇIQ bəyandır:
      1. `PORO:1-3` kimi XASSƏ-SPESİFİK bəyan (tapşırıq §4/§15)
      2. quyunun ümumi `1-3` bəyanı
      3. `policy=INTERVAL` isə quyunun FİZİKİ intervalı (AÇIQ seçim,
         xəbərdarlıqla)
      4. heç nə — bu quyu bu xassə üçün heç bir laya məlumat VERMİR
    """
    layers = per_property.get(target)
    if layers is not None:
        return list(layers)
    if default_layers is not None:
        return list(default_layers)
    if policy is LayerDataPolicy.INTERVAL:
        if not effective:
            warnings.append(
                f"'{well.name}' ({target}): interval siyasəti seçilib, amma quyunun "
                "lay üstü/altı verilməyib (və ya grid-dən kənardadır) — bu quyu bu "
                "xassə üçün heç bir laya məlumat vermir.")
            return []
        warnings.append(
            f"'{well.name}' ({target}): data layı bəyan edilməyib — İSTİFADƏÇİNİN "
            f"seçdiyi 'interval' siyasəti ilə quyunun fiziki intervalının kəsdiyi "
            f"laylar məlumat sayıldı. Bu, FƏRZİYYƏDİR: interval ölçmə ilə eyni "
            f"şey deyil.")
        return list(effective)
    warnings.append(
        f"'{well.name}' ({target}): data layı bəyan edilməyib — lay-məlumatlı "
        "rejimdə bu quyu bu xassə üçün HEÇ BİR laya məlumat vermir (səssiz "
        "yayılma tətbiq edilmir). 'Data layları' sütununu doldurun.")
    return []


def well_layer_summary(wells: List[GeologicalWell], geometry: CellGeometry,
                       targets: Sequence[str] = ("PORO", "PERMX", "SW"),
                       policy: LayerDataPolicy = LayerDataPolicy.STRICT
                       ) -> Dict[str, Dict[str, List[int]]]:
    """UI/hesabat üçün quyu-üzrə xülasə (tapşırıq §14).

    `{quyu adı: {"effective": [...], "data": [...], "PORO": [...], ...}}`
    — 0-əsaslı K indeksləri. HEÇ BİR interpolyasiya aparmır, yalnız
    BƏYANLARI oxuyur; ona görə UI-də hər cədvəl redaktəsindən sonra
    ucuz şəkildə çağırıla bilər.
    """
    nz = geometry.grid.nz
    summary: Dict[str, Dict[str, List[int]]] = {}
    for well in wells:
        entry: Dict[str, List[int]] = {}
        try:
            entry["effective"] = well_effective_layers(well, geometry)
        except ValueError:
            entry["effective"] = []
        try:
            default_layers, per_property = well.data_layer_sets(nz)
        except ValueError:
            summary[well.name] = {**entry, "data": [], "error": []}
            continue
        union: set = set()
        for attr, target in PROPERTY_MAP.items():
            if target in AREAL_TARGETS or target not in targets:
                continue
            if getattr(well, attr) is None:
                entry[target] = []
                continue
            layers = _resolve_layers(well, target, default_layers, per_property,
                                     entry["effective"], policy, [])
            entry[target] = layers
            union.update(layers)
        entry["data"] = sorted(union)
        summary[well.name] = entry
    return summary
