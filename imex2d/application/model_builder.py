"""İŞ AXINININ ADDIMI: Geoloji model → Rezervuar modeli.

Bu keçid əvvəllər ümumiyyətlə yox idi — simulyator birbaşa UI-dən
gələn skalyar dəyərlərdən özü üçün massivlər qururdu. İndi keçid
açıq şəkildə burada baş verir və izlənə bilir.
"""

from __future__ import annotations
from typing import Dict, List, Optional

from ..domain.geological_model import GeologicalModel
from ..domain.initial import InitialConditions
from ..domain.properties import FluidProperties, PropertyMap, RockProperties
from ..domain.pvt import PVTTable
from ..domain.reservoir_model import ReservoirModel
from ..domain.scal import CapillaryParameters, CoreyParameters
from ..domain.unit_conversions import to_engine_units
from ..domain.units import DEFAULT_UNITS, UnitSystem
from ..domain.wells import Well


class ReservoirModelBuilder:
    """Geoloji modeli mühəndis məlumatı ilə zənginləşdirib simulyasiyaya hazır model qurur."""

    PORO_KEY = "PORO"
    PERMX_KEY = "PERMX"
    PERMY_KEY = "PERMY"
    PERMZ_KEY = "PERMZ"
    NTG_KEY = "NTG"

    def build(self,
              geological_model: GeologicalModel,
              wells: List[Well],
              fluids: Optional[FluidProperties] = None,
              scal: Optional[CoreyParameters] = None,
              capillary: Optional[CapillaryParameters] = None,
              initial: Optional[InitialConditions] = None,
              pvt_table: Optional[PVTTable] = None,
              scal_tables=None,
              fault_references: Optional[List] = None,
              rock_compressibility: float = 4.5e-5,
              rock_compressibility_unit: str = "bar",
              name: Optional[str] = None,
              units: UnitSystem = DEFAULT_UNITS) -> ReservoirModel:
        """`rock_compressibility_unit` — `rock_compressibility`-nin HANSI
        vahiddə verildiyi (təzyiq vahidi, sıxılma = 1/[bu vahid]). Defolt
        `"bar"` — mühərrikin öz vahidi, ona görə `rock_compressibility_unit`
        DƏYİŞMƏYƏNDƏ dəyər ƏVVƏLKİ kimi HEÇ ÇEVRİLMİR (bax
        `to_engine_units`, `unit==target` olanda `convert()` no-op-dur)."""
        # BU, SİMULYATOR ÖNCƏSİ VALİDASİYA QAPISIDIR (tapşırıq §20):
        # `GeologicalModel.validate()` indi lay-üzrə MISSING xassələri də
        # AÇIQ, lay-adlı mesajla bildirir — yəni "L4-L5 üçün PORO yoxdur"
        # modeli simulyasiyaya BURAXILMIR. Tamamlanmış (estimated/
        # simulated) sahələr isə keçir, amma mənşəyi (`provenance`)
        # rezervuar modelinə DAŞINIR, itmir.
        issues = geological_model.validate()
        if issues:
            raise ValueError("Geoloji model natamamdır: " + "; ".join(issues))

        rock_compressibility = to_engine_units(
            rock_compressibility, rock_compressibility_unit, "compressibility")

        maps: Dict[str, PropertyMap] = dict(geological_model.property_maps)
        permx = geological_model.require(self.PERMX_KEY)
        permy = maps.get(self.PERMY_KEY, permx)

        rock = RockProperties(
            porosity=geological_model.require(self.PORO_KEY),
            permx=permx,
            permy=permy,
            permz=maps.get(self.PERMZ_KEY),
            net_to_gross=maps.get(self.NTG_KEY),
            compressibility=rock_compressibility,
        )

        return ReservoirModel(
            name=name or f"{geological_model.name} — rezervuar modeli",
            grid=geological_model.grid,
            geometry=geological_model.geometry,
            rock=rock,
            fluids=fluids or FluidProperties(),
            property_maps=maps,
            regions=geological_model.regions,
            fault_references=(fault_references if fault_references is not None
                              else geological_model.fault_references()),
            horizon_references=geological_model.horizon_references(),
            wells=list(wells),
            initial_conditions=initial or InitialConditions(),
            scal_parameters=scal or CoreyParameters(),
            capillary_parameters=capillary or CapillaryParameters(),
            pvt_table=pvt_table,
            scal_tables=scal_tables,
            units=units,
            source_geological_model=geological_model.name,
            provenance=dict(geological_model.provenance),
        )
