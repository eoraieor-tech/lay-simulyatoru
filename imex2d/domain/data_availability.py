"""LAY-ÜZRƏ MƏLUMAT MÖVCUDLUĞU — "hansı xassə hansı K-təbəqəsində
HƏQİQƏTƏN ölçülüb" sualının YEGANƏ cavab yeri.

BU MODULUN ƏSAS QAYDASI (tapşırıq §2): aşağıdakı ALTI anlayış BİR-BİRİNƏ
QARIŞDIRILMIR —

    1. quyunun FİZİKİ intervalı            (top/bottom, metr)
    2. grid təbəqə həndəsəsi               (`CellGeometry`, K sərhədləri)
    3. MƏLUMAT MÖVCUDLUĞU                  (bu modul — `PropertyAvailability`)
    4. İNTERPOLYASİYA HƏDƏFİ               (istifadəçinin seçdiyi laylar)
    5. TAMAMLAMA (completion) ÜSULU        (məlumatsız lay necə doldurulur)
    6. QEYRİ-MÜƏYYƏNLİK / ETİBARLILIQ      (`confidence`)

Xüsusilə: `top=2000, bottom=2210` HEÇ VAXT "L1–L5 üçün petrofiziki
məlumat var" DEMƏK DEYİL. Birincisi HƏNDƏSƏDİR (quyu hansı layları
KƏSİR), ikincisi MƏLUMATDIR (hansı layda ÖLÇMƏ var). Bu modul yalnız
(3)-ü təmsil edir; (1)→(2) uyğunlaşdırması `geology/layer_availability.py`
-dədir, (4)/(5) isə `application/geology_service.py`-dəki
`LayerInterpolationConfig`-dədir.

Bu modul SAF DOMAIN-dir: numpy-dan başqa asılılığı yoxdur, heç bir
interpolyasiya/geostatistika alqoritmini TANIMIR.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


class DataStatus(str, Enum):
    """Bir hüceyrənin/layın dəyəri HARADAN gəlir.

    QƏTİ QAYDA (tapşırıq §5/§26): `ESTIMATED`/`SIMULATED`/`EXTRAPOLATED`
    nəticə HEÇ VAXT `MEASURED` kimi qeyd edilmir. Status DƏYƏRLƏ
    BİRLİKDƏ saxlanılır (bax `properties.PropertyProvenance`) ki, model
    simulyatora ötürüləndə mənşəyi İTMƏSİN.
    """

    MEASURED = "measured"           #: quyuda ölçülüb (sərt data hüceyrəsi)
    INTERPOLATED = "interpolated"   #: həmin layın ÖZ sərt datası ilə hesablanıb
    ESTIMATED = "estimated"         #: məlumatsız lay, AÇIQ completion üsulu ilə
    EXTRAPOLATED = "extrapolated"   #: məlumat zərfindən KƏNARDA qiymətləndirilib
    SIMULATED = "simulated"         #: stoxastik realizasiya (SGS/SIS)
    PRESERVED = "preserved"         #: ORİJİNAL sahədən toxunulmadan gətirilib
    MISSING = "missing"             #: məlumat YOXDUR və tamamlanmayıb


#: "ən çox məlumatla dəstəklənəndən" "ən azına" — hesabat/aqreqasiya sırası.
STATUS_ORDER: Tuple[DataStatus, ...] = (
    DataStatus.MEASURED, DataStatus.INTERPOLATED, DataStatus.ESTIMATED,
    DataStatus.SIMULATED, DataStatus.EXTRAPOLATED, DataStatus.PRESERVED,
    DataStatus.MISSING)

#: 3D görüntüləmə üçün ƏDƏDİ kod (rəng şkalası). SIRALI DEYİL, KATEQORİKDİR
#: — "2 > 1 deməli daha yaxşıdır" MƏNASI YOXDUR, yalnız rəng ayırd etmək üçün.
STATUS_CODE: Dict[str, int] = {status.value: index
                               for index, status in enumerate(STATUS_ORDER)}

#: İstifadəçiyə göstərilən qısa etiketlər.
STATUS_LABEL: Dict[str, str] = {
    DataStatus.MEASURED.value: "ÖLÇÜLÜB",
    DataStatus.INTERPOLATED.value: "İNTERPOLYASİYA",
    DataStatus.ESTIMATED.value: "QİYMƏTLƏNDİRİLİB",
    DataStatus.EXTRAPOLATED.value: "EKSTRAPOLYASİYA",
    DataStatus.SIMULATED.value: "SİMULYASİYA",
    DataStatus.PRESERVED.value: "ORİJİNAL (toxunulmayıb)",
    DataStatus.MISSING.value: "MƏLUMAT YOXDUR",
}


@dataclass(frozen=True)
class LayerStatus:
    """Bir xassənin bir K-təbəqəsindəki vəziyyəti.

    `confidence` — `[0, 1]` ORDİNAL dəstək balı və ya `None`. `None`
    "hesablanmadı" DEMƏKDİR, "sıfır etibar" DEYİL: tapşırıq §18-ə görə
    ƏSASLANDIRILA BİLMƏYƏN halda SAXTA rəqəm YARADILMIR (bax
    `geology_service._layer_confidence`).
    """

    k: int
    status: DataStatus = DataStatus.MISSING
    method: str = ""
    confidence: Optional[float] = None
    n_data: int = 0
    note: str = ""

    @property
    def has_data(self) -> bool:
        """Bu layda HƏQİQİ (sərt) ölçmə varmı — `status`-dan ASILI DEYİL."""
        return self.n_data > 0


@dataclass
class PropertyAvailability:
    """BİR xassənin BÜTÜN K-təbəqələri üzrə vəziyyəti.

    Xassə-spesifikdir (tapşırıq §4): PORO L4-də ola bilər, PERMX olmaya
    bilər — ona görə hər xassənin ÖZ obyekti var, ortaq "lay siyahısı"
    YOXDUR.
    """

    name: str
    nz: int
    layers: Dict[int, LayerStatus] = field(default_factory=dict)

    def __post_init__(self):
        for k in range(self.nz):
            self.layers.setdefault(k, LayerStatus(k=k))

    # ---------------------------------------------------------- sorğular
    def status(self, k: int) -> DataStatus:
        return self.layers[k].status

    def set(self, k: int, **changes) -> None:
        """Bir layın vəziyyətini əvəz edir (frozen `LayerStatus` → yeni obyekt)."""
        if not 0 <= k < self.nz:
            raise ValueError(f"{self.name}: K={k} grid diapazonundan kənardadır "
                             f"(0..{self.nz - 1}).")
        current = self.layers[k]
        self.layers[k] = LayerStatus(
            k=k,
            status=changes.get("status", current.status),
            method=changes.get("method", current.method),
            confidence=changes.get("confidence", current.confidence),
            n_data=changes.get("n_data", current.n_data),
            note=changes.get("note", current.note))

    def data_layers(self) -> List[int]:
        """SƏRT DATASI OLAN laylar — `status`-dan ASILI DEYİL (giriş faktı)."""
        return [k for k in range(self.nz) if self.layers[k].has_data]

    def layers_with(self, *statuses: DataStatus) -> List[int]:
        wanted = {s.value for s in statuses}
        return [k for k in range(self.nz) if self.layers[k].status.value in wanted]

    def missing_layers(self) -> List[int]:
        return self.layers_with(DataStatus.MISSING)

    def status_codes(self) -> np.ndarray:
        """`(nz,)` ədədi kod massivi — 3D görüntü/filtr üçün."""
        return np.asarray([STATUS_CODE[self.layers[k].status.value]
                           for k in range(self.nz)], dtype=float)

    def confidences(self) -> np.ndarray:
        """`(nz,)` — hesablanmayan laylar `NaN` (SAXTA rəqəm YOXDUR)."""
        return np.asarray([np.nan if self.layers[k].confidence is None
                           else float(self.layers[k].confidence)
                           for k in range(self.nz)], dtype=float)

    def as_text(self) -> str:
        lines = [f"{self.name}:"]
        for k in range(self.nz):
            entry = self.layers[k]
            confidence = ("—" if entry.confidence is None
                          else f"{entry.confidence:.2f}")
            method = f" [{entry.method}]" if entry.method else ""
            lines.append(
                f"  L{k + 1} (K={k}): {STATUS_LABEL[entry.status.value]}{method}"
                f"  data={entry.n_data}  etibar={confidence}"
                + (f"  — {entry.note}" if entry.note else ""))
        return "\n".join(lines)


@dataclass
class ModelDataAvailability:
    """Xassə adı → `PropertyAvailability`. Model ilə BİRLİKDƏ saxlanılır."""

    nz: int
    properties: Dict[str, PropertyAvailability] = field(default_factory=dict)

    def __contains__(self, name: str) -> bool:
        return name in self.properties

    def __getitem__(self, name: str) -> PropertyAvailability:
        return self.properties[name]

    def get(self, name: str) -> Optional[PropertyAvailability]:
        return self.properties.get(name)

    def require(self, name: str) -> PropertyAvailability:
        if name not in self.properties:
            self.properties[name] = PropertyAvailability(name=name, nz=self.nz)
        return self.properties[name]

    def names(self) -> List[str]:
        return sorted(self.properties)

    def as_text(self) -> str:
        return "\n".join(self.properties[name].as_text() for name in self.names())


# ═════════════════════════════════════════════ lay seçimi mətn formatı
def format_layers(layers: Iterable[int], one_based: bool = True) -> str:
    """`[0, 1, 2, 4]` → `"1-3,5"` (defolt 1-əsaslı, İSTİFADƏÇİ görünüşü).

    Boş siyahı → `"—"` (BOŞ SƏTİR DEYİL) ki, UI-da "yoxdur" halı
    təsadüfən "hamısı" kimi oxunmasın.
    """
    values = sorted({int(k) for k in layers})
    if not values:
        return "—"
    offset = 1 if one_based else 0
    parts: List[str] = []
    start = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        parts.append(_range_text(start + offset, previous + offset))
        start = previous = value
    parts.append(_range_text(start + offset, previous + offset))
    return ",".join(parts)


def _range_text(low: int, high: int) -> str:
    return str(low) if low == high else f"{low}-{high}"


def parse_layers(text: str, nz: int, one_based: bool = True) -> List[int]:
    """`"1-3,5"` → `[0, 1, 2, 4]` (0-əsaslı MÜHƏRRIK indeksləri).

    QƏBUL EDİLƏNLƏR: `"1-3"`, `"1,2,3"`, `"1-3, 5"`, `"*"`/`"all"`/
    `"hamısı"` (bütün laylar), boş mətn (`[]`).

    XƏTA (SƏSSİZ DÜZƏLİŞ YOXDUR — tapşırıq §23.2): diapazondan kənar
    indeks, tərsinə diapazon (`3-1`), rəqəm olmayan simvol.
    """
    text = (text or "").strip()
    if not text:
        return []
    if text.lower() in {"*", "all", "hamisi", "hamısı"}:
        return list(range(nz))
    offset = 1 if one_based else 0
    result: set = set()
    for token in text.replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue
        low_text, separator, high_text = token.partition("-")
        try:
            low = int(low_text.strip())
            high = int(high_text.strip()) if separator else low
        except ValueError:
            raise ValueError(
                f"Lay seçimi oxunmadı: {token!r} (gözlənilən format: '1-3' və ya '1,2,5').")
        if high < low:
            raise ValueError(
                f"Lay diapazonu tərsinədir: {token!r} (aşağı hədd yuxarıdan böyükdür).")
        for value in range(low - offset, high - offset + 1):
            if not 0 <= value < nz:
                raise ValueError(
                    f"Lay {value + offset} grid diapazonundan kənardadır "
                    f"(1..{nz}) — NZ = {nz}.")
            result.add(value)
    return sorted(result)


def parse_property_layers(text: str, nz: int, one_based: bool = True
                          ) -> Tuple[Optional[List[int]], Dict[str, List[int]]]:
    """Quyu cədvəlinin "Data layları" sütununu oxuyur.

    İKİ format DƏSTƏKLƏNİR:

        "1-3"                     → BÜTÜN xassələr üçün L1–L3
        "PORO:1-5; PERMX:1-3"     → XASSƏ-SPESİFİK (tapşırıq §4/§15)

    Qaytarır `(default_layers, per_property)`; `default_layers` `None`
    olanda "bu quyu üçün ümumi lay siyahısı verilməyib" deməkdir (yalnız
    adı çəkilən xassələr üçün məlumat var).
    """
    text = (text or "").strip()
    if not text:
        return None, {}
    if ":" not in text:
        return parse_layers(text, nz, one_based), {}

    default: Optional[List[int]] = None
    per_property: Dict[str, List[int]] = {}
    for chunk in _split_property_chunks(text):
        name, separator, layer_text = chunk.partition(":")
        if not separator:
            default = parse_layers(chunk, nz, one_based)
            continue
        key = name.strip().upper()
        if not key:
            raise ValueError(f"Xassə adı boşdur: {chunk!r}.")
        per_property[key] = parse_layers(layer_text, nz, one_based)
    return default, per_property


def _split_property_chunks(text: str) -> List[str]:
    """`"PORO:1-3,5; PERMX:1-2"` → `["PORO:1-3,5", "PERMX:1-2"]`.

    Vergül HƏM xassələri, HƏM DƏ lay siyahısını ayıra bildiyi üçün sadə
    `split(",")` KİFAYƏT ETMİR: yeni parça YALNIZ növbəti `":"` görünəndə
    başlayır.
    """
    chunks: List[str] = []
    current: List[str] = []
    for token in text.replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue
        if ":" in token and current:
            chunks.append(",".join(current))
            current = [token]
        else:
            current.append(token)
    if current:
        chunks.append(",".join(current))
    return chunks


def describe_layer_sets(effective: Sequence[int], data: Sequence[int],
                        interpolated: Sequence[int], missing: Sequence[int]) -> str:
    """UI/hesabat üçün BİR sətir — dörd anlayış AÇIQ AYRILIR (tapşırıq §14)."""
    return (f"kəsir: {format_layers(effective)} | "
            f"data: {format_layers(data)} | "
            f"interp: {format_layers(interpolated)} | "
            f"yoxdur: {format_layers(missing)}")
