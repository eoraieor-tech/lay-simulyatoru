"""Struktur elementlərə İSTİNADLAR və regionlar.

Əsas qərar: geoloq və mühəndis fərqli məlumatla işləyir. Geoloq
`GeologicalModel.Fault`-u yaradır — hansı I/J/K müstəvisində, hansı
diapazonda. Mühəndis isə simulyasiya üçün YALNIZ transmissivlik
çarpanını dəyişir, geoloji detalları bilməli deyil.

`FaultReference` bu ikisinin BİŞMİŞ (resolved) nəticəsidir: modelin
qurulma anında `GeologicalModel.fault_references()` vasitəsilə
həndəsə köçürülür, mühəndis sonra yalnız `transmissibility_multiplier`
və `sealing`-i dəyişir. Diskretizasiya birbaşa bunu oxuyur — geoloji
modelə geri qayıtmır, çünki `ReservoirModel` onu saxlamır (yalnız adı,
`source_geological_model`).
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from .properties import PropertyMap

FAULT_AXES = ("I", "J", "K")


@dataclass(frozen=True)
class HorizonReference:
    """Geoloji modeldəki horizonta istinad."""
    name: str
    source_id: str
    role: str = "top"


@dataclass(frozen=True)
class FaultReference:
    """Geoloji modeldəki faya istinad + simulyasiya parametri.

    `axis`, `plane_index`, `range_a`, `range_b` — geoloji modeldən
    bişmiş həndəsə (bax modul sənədləşməsi). Köhnə layihə fayllarında
    bu sahələr yox idi, ona görə hamısı `None` defoltludur və belə
    fay heç bir bağlantıya təsir etmir (`matches()` həmişə False).
    """
    name: str
    source_id: str
    transmissibility_multiplier: float = 1.0
    sealing: bool = False
    axis: Optional[str] = None
    plane_index: Optional[int] = None
    range_a: Optional[Tuple[int, int]] = None
    range_b: Optional[Tuple[int, int]] = None

    def __post_init__(self):
        if self.axis is not None and self.axis.upper() not in FAULT_AXES:
            raise ValueError(f"'{self.name}': ox 'I', 'J' və ya 'K' olmalıdır "
                             f"(alındı: '{self.axis}').")
        if self.transmissibility_multiplier < 0.0:
            raise ValueError(f"'{self.name}': çarpan mənfi ola bilməz.")

    @property
    def has_geometry(self) -> bool:
        return self.axis is not None and self.plane_index is not None

    @property
    def effective_multiplier(self) -> float:
        """`sealing=True` çarpandan asılı olmayaraq axını sıfırlayır."""
        return 0.0 if self.sealing else self.transmissibility_multiplier

    @property
    def is_sealing(self) -> bool:
        return self.sealing or self.transmissibility_multiplier <= 1e-12

    def validate(self, grid) -> List[str]:
        if not self.has_geometry:
            return []
        issues = []
        limits = {"I": grid.nx, "J": grid.ny, "K": grid.nz}
        axis = self.axis.upper()
        limit = limits[axis]
        if self.plane_index < 0 or self.plane_index >= limit - 1:
            issues.append(
                f"'{self.name}': müstəvi indeksi ({self.plane_index}) "
                f"grid-dən kənardadır ({axis} oxunda maks {max(limit - 2, 0)}).")
        other_axes = [a for a in FAULT_AXES if a != axis]
        for label, bounds in zip(other_axes, (self.range_a, self.range_b)):
            if bounds is None:
                continue
            low, high = bounds
            axis_limit = limits[label]
            if low > high:
                issues.append(f"'{self.name}': {label} diapazonu tərsdir.")
            elif low < 0 or high >= axis_limit:
                issues.append(
                    f"'{self.name}': {label} diapazonu ({low}..{high}) "
                    f"grid-dən kənardadır (0..{axis_limit - 1}).")
        return issues

    def matches(self, connection_axis_code: int, coordinate_a,
               coordinate_b) -> np.ndarray:
        """Bu fayın üzərində olan bağlantıların maskası.

        `connection_axis_code` — `Connections.axis` kodu (0=I, 1=J, 2=K).
        `coordinate_a`/`coordinate_b` digər iki oxun (bağlantının `cell_a`
        tərəfindən çıxarılmış) koordinat massivləridir.
        """
        if not self.has_geometry:
            return np.zeros(np.shape(coordinate_a), dtype=bool)
        if FAULT_AXES.index(self.axis.upper()) != connection_axis_code:
            return np.zeros(np.shape(coordinate_a), dtype=bool)
        mask = np.ones(np.shape(coordinate_a), dtype=bool)
        for bounds, coordinate in ((self.range_a, coordinate_a),
                                   (self.range_b, coordinate_b)):
            if bounds is not None:
                low, high = bounds
                mask &= (coordinate >= low) & (coordinate <= high)
        return mask

    def summary(self) -> str:
        if not self.has_geometry:
            return f"{self.name}: həndəsəsiz (fault, çarpan={self.transmissibility_multiplier:g})"
        span = lambda bounds: "hamısı" if bounds is None else f"{bounds[0]}-{bounds[1]}"
        other_axes = [a for a in FAULT_AXES if a != self.axis.upper()]
        return (f"{self.name}: {self.axis}={self.plane_index}|"
               f"{self.plane_index + 1}  {other_axes[0]}={span(self.range_a)}  "
               f"{other_axes[1]}={span(self.range_b)}  "
               f"çarpan={self.effective_multiplier:g}"
               f"{'  (BAĞLI)' if self.is_sealing else ''}")


@dataclass
class RegionSet:
    """Hüceyrələrin regionlara bölgüsü.

    Region ID-lər gələcəkdə SCAL və PVT cədvəllərinin seçilməsi üçün
    provider-lərə ötürüləcək (rock region / PVT region).
    """
    region_id: PropertyMap
    names: Dict[int, str] = field(default_factory=dict)

    @classmethod
    def single(cls, ncell: int, name: str = "REGION-1") -> "RegionSet":
        return cls(PropertyMap.uniform("REGION_ID", 1, ncell), {1: name})

    @property
    def ids(self) -> np.ndarray:
        return np.unique(self.region_id.values.astype(int))
