"""SPE1CASE2 (Odeh 1981, OPM `opm-tests/spe1`) — rezervuar modeli və etalon müqayisəsi.

Cədvəllər (SWOF, SGOF, PVTW, PVDG, PVTO) deck faylından OXUNUR. Qalan
parametrlər deck-dən sətir-sətir yoxlanılıb (`SPE1.md` §2) və burada FIELD
vahidlərində, deck açar sözünün adı ilə yazılıb — çevirmə yalnız
`domain/unit_conversions.py` ilə aparılır.

Deck faylı repoda DEYİL (lisenziya yoxlanılmayıb) — bax `SPE1.md` §1.

⏳ MƏLUM FƏRQLƏR (açıq yazılır):
  * Üç fazalı kro — Stone II (Eclipse-in defoltu başqa modeldir). SPE1-də
    su hərəkətsiz qaldığı üçün (Sw ≈ Swc) hər ikisi kro = krog verir.
  * Doymuş qolun 5014.7 psia-dan yuxarı davranışı (G7) — plato.
  * Quyu lüləsində hidrostatik hədd yoxdur; SPE1-də hər quyunun bir
    perforasiyası var və BHP istinad dərinliyi perforasiyanın mərkəzidir.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np

from ..application.config import SimulationConfig, TimeSteppingConfig, OutputConfig
from ..application.model_builder import ReservoirModelBuilder
from ..domain.geological_model import GeologicalModel
from ..domain.geometry import CellGeometry
from ..domain.grid import CartesianGrid
from ..domain.initial import InitialConditions
from ..domain.properties import FluidProperties, PropertyMap
from ..domain.reservoir_model import ReservoirModel
from ..domain.scal import CoreyParameters, GasCoreyParameters
from ..domain.structure import RegionSet
from ..domain.unit_conversions import convert
from ..domain.wells import (ControlMode, Perforation, Phase, RateBasis, Well,
                            WellControl, WellType)
from ..io.pvt_io import read_deck_pvt
from ..io.scal_io import read_sgof, read_swof

#: `DIMENS 10 10 3`
DIMENS = (10, 10, 3)
#: `DX = DY = 1000` ft, `DZ 20 30 50` ft, `TOPS 8325` ft
DX_FT, DY_FT, DZ_FT, TOPS_FT = 1000.0, 1000.0, (20.0, 30.0, 50.0), 8325.0
#: `PERMX = PERMY = PERMZ` qatlar üzrə, mD; `PORO 0.3`
PERM_MD = (500.0, 50.0, 200.0)
PORO = 0.3
#: `ROCK 14.7 3E-6` — istinad psia, sıxılma 1/psi
ROCK_REFERENCE_PSIA, ROCK_COMPRESSIBILITY_PSI = 14.7, 3.0e-6
#: `DENSITY 53.66 64.49 0.0533` lb/ft³ (neft, su, qaz)
DENSITY_LB_FT3 = (53.66, 64.49, 0.0533)
#: `EQUIL 8400 4800 8450 0 8300 0`
DATUM_FT, DATUM_PSIA, OWC_FT, GOC_FT = 8400.0, 4800.0, 8450.0, 8300.0
#: `RSVD` 8300–8450 ft: 1.270 Mscf/STB
INITIAL_RS_MSCF_STB = 1.270
#: `COMPDAT` diametr 0.5 ft
WELL_RADIUS_FT = 0.25
#: `WCONPROD 'PROD' 'OPEN' 'ORAT' 20000 4* 1000`, `COMPDAT 'PROD' 10 10 3 3`
PRODUCER_ORAT_STB, PRODUCER_MIN_BHP_PSIA, PRODUCER_IJK = 20000.0, 1000.0, (9, 9, 2)
#: `WCONINJE 'INJ' 'GAS' 'OPEN' 'RATE' 100000 1* 9014`, `COMPDAT 'INJ' 1 1 1 1`
INJECTOR_RATE_MSCF, INJECTOR_MAX_BHP_PSIA, INJECTOR_IJK = 100000.0, 9014.0, (0, 0, 0)
#: `TSTEP` — 120 ay, cəmi 3650 gün
END_TIME_DAYS = 3650.0


def _ft(value):
    return convert(value, "ft", "m", "length")


def _psi(value):
    return convert(value, "psi", "bar", "pressure")


def build_spe1case2_model(deck_path: str) -> ReservoirModel:
    """SPE1CASE2 rezervuar modeli — mühərrik (METRIC) vahidlərində."""
    nx, ny, nz = DIMENS
    grid = CartesianGrid(nx, ny, nz)
    geometry = CellGeometry(grid, _ft(DX_FT), _ft(DY_FT),
                            [_ft(dz) for dz in DZ_FT], top_depth=_ft(TOPS_FT))
    n = grid.ncell
    layer_perm = np.repeat(np.asarray(PERM_MD, float), nx * ny)

    geology = GeologicalModel(name="SPE1CASE2", grid=grid, geometry=geometry,
                              regions=RegionSet.single(n))
    geology.add_property(PropertyMap.from_array("PORO", PORO, n))
    for key in ("PERMX", "PERMY", "PERMZ"):
        geology.add_property(PropertyMap.from_array(key, layer_perm, n, "mD"))

    oil_density, water_density, gas_density = (
        convert(value, "lb/ft3", "kg/m3", "density") for value in DENSITY_LB_FT3)
    deck = read_deck_pvt(deck_path)
    water_oil = read_swof(deck_path)
    gas_oil = read_sgof(deck_path)
    table = water_oil.tables[min(water_oil.tables)]
    swc = float(table.sw[0])
    # Corey parametrləri YALNIZ doyma hədləri üçündür (diaqnostika və ilkin
    # tarazlıq onları oxuyur) — əyrilər cədvəldən gəlir.
    immobile_oil = table.sw[table.kro <= 0.0]
    sor = 1.0 - float(immobile_oil[0]) if immobile_oil.size else 0.0

    wells = [
        Well("PROD", WellType.PRODUCER,
             WellControl(ControlMode.RATE,
                         convert(PRODUCER_ORAT_STB, "stb/day", "m3/day", "rate"),
                         bhp_limit=_psi(PRODUCER_MIN_BHP_PSIA),
                         rate_basis=RateBasis.SURFACE),
             [Perforation(*PRODUCER_IJK)], radius=_ft(WELL_RADIUS_FT)),
        Well("INJ", WellType.INJECTOR,
             WellControl(ControlMode.RATE,
                         convert(INJECTOR_RATE_MSCF, "Mscf/day", "m3/day", "rate"),
                         Phase.GAS, bhp_limit=_psi(INJECTOR_MAX_BHP_PSIA),
                         rate_basis=RateBasis.SURFACE),
             [Perforation(*INJECTOR_IJK)], radius=_ft(WELL_RADIUS_FT)),
    ]

    initial = InitialConditions(
        datum_depth=_ft(DATUM_FT), datum_pressure=_psi(DATUM_PSIA),
        water_saturation=swc, oil_water_contact=_ft(OWC_FT),
        gas_oil_contact=_ft(GOC_FT), use_equilibration=True,
        solution_gor=convert(INITIAL_RS_MSCF_STB, "Mscf/stb", "sm3/sm3",
                             "solution_gor"))

    return ReservoirModelBuilder().build(
        geology, wells,
        fluids=FluidProperties(oil_density=oil_density,
                               water_density=water_density,
                               gas_density=gas_density),
        scal=CoreyParameters(swc=swc, sor=sor),
        gas_scal=GasCoreyParameters(),
        initial=initial,
        pvt_table=deck.to_pvt_table(source="SPE1CASE2.DATA"),
        pvt_oil_branches=deck.oil.branches,
        scal_tables=water_oil, gas_scal_tables=gas_oil,
        rock_compressibility=ROCK_COMPRESSIBILITY_PSI,
        rock_compressibility_unit="psi",
        rock_compressibility_reference=ROCK_REFERENCE_PSIA,
        rock_compressibility_reference_unit="psi",
        name="SPE1CASE2")


def spe1case2_config(max_dt: float = 31.0, snapshots: int = 120) -> SimulationConfig:
    return SimulationConfig(
        end_time=END_TIME_DAYS,
        time_stepping=TimeSteppingConfig(max_dt=max_dt),
        output=OutputConfig(snapshot_count=snapshots))


# ═══════════════════════════════ müqayisə ════════════════════════════

@dataclass
class ComparisonRow:
    time: float
    quantity: str
    unit: str
    simulated: float
    reference: float

    @property
    def relative_error(self) -> float:
        scale = max(abs(self.reference), 1e-12)
        return (self.simulated - self.reference) / scale


def simulated_series(model: ReservoirModel, result) -> Dict[str, tuple]:
    """Nəticəni etalon etiketlərinə və FIELD vahidlərinə çevirir.

    Qaytarır: `{etiket: (zaman, qiymətlər, vahid)}`.
    """
    series = result.series
    time = np.asarray(series.time, float)
    out: Dict[str, tuple] = {}
    if series.oil_rate:
        out["FOPR"] = (time, convert(np.asarray(series.oil_rate, float),
                                     "m3/day", "stb/day", "rate"), "STB/DAY")
    if series.gas_oil_ratio:
        out["FGOR"] = (time, convert(np.asarray(series.gas_oil_ratio, float),
                                     "sm3/sm3", "Mscf/stb", "solution_gor"),
                       "MSCF/STB")
    for well, values in result.well_bhp.items():
        if len(values) == time.size:
            out[f"WBHP:{well}"] = (time, convert(np.asarray(values, float),
                                                 "bar", "psi", "pressure"), "PSIA")
    if result.snapshots:
        snap_time = np.asarray([s.time for s in result.snapshots], float)
        for cell in (0, model.ncell - 1):
            values = np.asarray([s.pressure[cell] for s in result.snapshots], float)
            out[f"BPR:{cell + 1}"] = (snap_time, convert(values, "bar", "psi",
                                                         "pressure"), "PSIA")
    return out


def compare_with_reference(model: ReservoirModel, result, reference,
                           times=(1.0, 304.0, 1034.0, 1399.0, 1580.0, 2129.0, 3650.0),
                           labels=("FOPR", "FGOR", "WBHP:PROD", "WBHP:INJ",
                                   "BPR:1", "BPR:300")) -> List[ComparisonRow]:
    """`SPE1.md` §3-dəki nöqtələrdə simulyasiya ↔ OPM Flow etalonu."""
    simulated = simulated_series(model, result)
    rows: List[ComparisonRow] = []
    for label in labels:
        if label not in simulated or label not in reference:
            continue
        sim_time, sim_values, unit = simulated[label]
        for time in times:
            if time < sim_time[0] or time > sim_time[-1]:
                continue
            rows.append(ComparisonRow(
                time, label, unit,
                float(np.interp(time, sim_time, sim_values)),
                reference.at(label, time)))
    return rows


def first_time_below(time: np.ndarray, values: np.ndarray,
                     threshold: float) -> Optional[float]:
    """`values < threshold` olan ilk an (məs. FOPR 20 000-dən aşağı)."""
    below = np.nonzero(np.asarray(values) < threshold)[0]
    return float(np.asarray(time)[below[0]]) if below.size else None
