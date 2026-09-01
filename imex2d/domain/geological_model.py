"""GEOLOJİ MODEL — yerin statik təsviri.

Nə var: grid topologiyası, həndəsə, xassə xəritələri, regionlar,
horizontlar, faylar.

Nə YOXDUR (qəsdən): quyular, flüidlər, ilkin şərtlər, simulyasiya
parametrləri. Geoloq bu obyekti quyu və istismar məlumatı olmadan
yaradır və dəyişir.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .facies_field import FaciesField
from .geometry import CellGeometry
from .grid import CartesianGrid
from .properties import PropertyMap
from .structure import FaultReference, HorizonReference, RegionSet
from .validation import validate_permeability, validate_porosity


@dataclass
class Horizon:
    """Horizontun özü (geoloji modelə məxsusdur)."""
    name: str
    depth_map: Optional[PropertyMap] = None


@dataclass
class Fault:
    """Fayın özü (geoloji modelə məxsusdur).

    `polyline`/`throw`/`dip` — TƏSVİRİ məlumatdır, xəritədə göstərmək
    üçün; hesablamaya təsir etmir. Axın hesablaması `axis`/`plane_index`/
    `range_a`/`range_b` — grid oxlarına düz bucaqlı müstəvi — ilə
    aparılır. Bu, Eclipse-in `FAULTS` açar sözü ilə eyni fikirdir və
    hazırkı bərabər-blok həndəsəsi (`CellGeometry`) ilə uyğundur;
    corner-point-in real əyri fay səthini yalnız approksimasiya edir
    (bax `ECLIPSE_IO.md`).

    `axis`/`plane_index` verilməyibsə (`None`), fay heç bir bağlantıya
    təsir etmir — mövcud olması, lakin axına toxunmaması qanunidir
    (məsələn yalnız xəritədə göstərmək üçün).
    """
    name: str
    polyline: Optional[list] = None
    throw: float = 0.0
    dip: float = 90.0
    axis: Optional[str] = None
    plane_index: Optional[int] = None
    range_a: Optional[tuple] = None
    range_b: Optional[tuple] = None
    transmissibility_multiplier: float = 1.0
    sealing: bool = False


@dataclass
class GeologicalModel:
    name: str
    grid: CartesianGrid
    geometry: CellGeometry
    property_maps: Dict[str, PropertyMap] = field(default_factory=dict)
    regions: Optional[RegionSet] = None
    horizons: List[Horizon] = field(default_factory=list)
    faults: List[Fault] = field(default_factory=list)
    coordinate_system: str = "LOCAL"
    #: Kateqorik (SIS ilə yaradılan) fasiya sahələri — `regions`-dan
    #: (SATNUM/SCAL region, stoxastik DEYİL) AYRICA saxlanılır, bax
    #: `domain/facies_field.py` modul docstring-i.
    facies_fields: Dict[str, FaciesField] = field(default_factory=dict)

    def __post_init__(self):
        if self.regions is None:
            self.regions = RegionSet.single(self.grid.ncell)

    def add_facies_field(self, facies: FaciesField) -> None:
        if facies.ncell != self.grid.ncell:
            raise ValueError(
                f"{facies.name}: hüceyrə sayı uyğun gəlmir ({facies.ncell} != {self.grid.ncell})")
        self.facies_fields[facies.name] = facies

    def add_property(self, prop: PropertyMap) -> None:
        if prop.values.size != self.grid.ncell:
            raise ValueError(f"{prop.name}: hüceyrə sayı uyğun gəlmir")
        self.property_maps[prop.name] = prop

    def require(self, name: str) -> PropertyMap:
        if name not in self.property_maps:
            raise KeyError(f"Geoloji modeldə '{name}' xassəsi yoxdur")
        return self.property_maps[name]

    def horizon_references(self) -> List[HorizonReference]:
        return [HorizonReference(h.name, h.name) for h in self.horizons]

    def fault_references(self) -> List[FaultReference]:
        """Geoloji fayları simulyasiya üçün BİŞİRİR.

        Bu, `ReservoirModel` qurularkən bir dəfə çağırılır (bax
        `ReservoirModelBuilder`) — sonrakı dəyişikliklər (çarpan,
        sealing) birbaşa `FaultReference` üzərində edilir, geoloji
        modelə qayıtmır.
        """
        return [FaultReference(
            name=f.name, source_id=f.name,
            transmissibility_multiplier=f.transmissibility_multiplier,
            sealing=f.sealing, axis=f.axis, plane_index=f.plane_index,
            range_a=f.range_a, range_b=f.range_b) for f in self.faults]

    def validate(self) -> list:
        """Sərt xətalar: məcburi xəritələr, dejenerativ həndəsə, fiziki
        cəhətdən qeyri-mümkün PORO/PERM dəyəri (bax `domain/validation.py`,
        Phase 1). Qeyri-adi-amma-mümkün diapazonlar üçün `validate_
        warnings()`-ə bax — burada RƏDD EDİLMİR."""
        issues = []
        for req in ("PORO", "PERMX"):
            if req not in self.property_maps:
                issues.append(f"Geoloji modeldə məcburi '{req}' xəritəsi yoxdur.")
        issues.extend(self.geometry.validate())
        if "PORO" in self.property_maps:
            issues.extend(validate_porosity(self.property_maps["PORO"].values, "PORO").errors)
        for key in ("PERMX", "PERMY", "PERMZ"):
            if key in self.property_maps:
                issues.extend(validate_permeability(self.property_maps[key].values, key).errors)
        return issues

    def validate_warnings(self) -> list:
        """Rədd edilməyən, qeyri-adi diapazon xəbərdarlıqları."""
        warnings = []
        if "PORO" in self.property_maps:
            warnings += validate_porosity(self.property_maps["PORO"].values, "PORO").warnings
        for key in ("PERMX", "PERMY", "PERMZ"):
            if key in self.property_maps:
                warnings += validate_permeability(self.property_maps[key].values, key).warnings
        return warnings
