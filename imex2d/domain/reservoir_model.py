"""REZERVUAR MODELİ — hər simulyasiya üçün YEGANƏ HƏQİQƏT MƏNBƏYİ.

Hesablama mühərriki bundan başqa heç bir yerdən məlumat götürmür və
öz daxili modelini QURMUR. UI, fayl oxuyucusu, skript — hamısı eyni
obyekti hazırlayıb mühərrikə ötürür.

Tərkibi (tapşırıqda tələb olunan siyahı):
    grid, geometry, rock, fluids (placeholder), property_maps, regions,
    fault_references, horizon_references, wells, initial_conditions
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from .diagnostics import DiagnosticReport, Severity
from .geometry import CellGeometry
from .grid import CartesianGrid, Connections
from .initial import InitialConditions
from .properties import FluidProperties, PropertyMap, RockProperties
from .pvt import PVTTable
from .scal import CapillaryParameters, CoreyParameters
from .structure import FaultReference, HorizonReference, RegionSet
from .units import DEFAULT_UNITS, UnitSystem
from .validation import validate_query_range
from .wells import ControlMode, Well


@dataclass
class ReservoirModel:
    name: str
    grid: CartesianGrid
    geometry: CellGeometry
    rock: RockProperties
    fluids: FluidProperties = field(default_factory=FluidProperties)
    property_maps: Dict[str, PropertyMap] = field(default_factory=dict)
    regions: Optional[RegionSet] = None
    fault_references: List[FaultReference] = field(default_factory=list)
    horizon_references: List[HorizonReference] = field(default_factory=list)
    wells: List[Well] = field(default_factory=list)
    initial_conditions: InitialConditions = field(default_factory=InitialConditions)
    scal_parameters: CoreyParameters = field(default_factory=CoreyParameters)
    capillary_parameters: CapillaryParameters = field(default_factory=CapillaryParameters)
    pvt_table: Optional[PVTTable] = None
    scal_tables: Optional[object] = None
    """Laboratoriya SCAL cədvəlləri (`SaturationTableSet`).

    Verilibsə, `scal_parameters`-dəki Corey düsturunun yerinə işlədilir.
    Tip annotasiyası `object`-dir ki, domain qatı `scal_tables` modulundan
    məcburi asılı olmasın — köhnə layihə faylları onsuz da yüklənməlidir.
    """
    units: UnitSystem = DEFAULT_UNITS
    source_geological_model: str = ""
    #: Xassə adı -> `PropertyProvenance` — geoloji modeldən GƏTİRİLİR
    #: (`ReservoirModelBuilder`), simulyasiya HESABLAMASINA HEÇ CÜR təsir
    #: ETMİR. Məqsəd: "bu PORO ölçülüb, yoxsa qiymətləndirilib" sualının
    #: cavabı simulyasiya modelinə çatanda İTMƏSİN (tapşırıq §20 —
    #: provenance saxlanılmalıdır) və 3D görüntü onu göstərə bilsin.
    #: Defolt boş lüğət — köhnə modellər üçün HEÇ NƏ dəyişmir.
    provenance: Dict[str, object] = field(default_factory=dict)

    _connections: Optional[Connections] = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        if self.regions is None:
            self.regions = RegionSet.single(self.grid.ncell)

    @property
    def ncell(self) -> int:
        return self.grid.ncell

    def connections(self) -> Connections:
        """Qonşuluq qrafı bir dəfə qurulur və saxlanılır."""
        if self._connections is None:
            self._connections = self.grid.build_connections()
        return self._connections

    def active_wells(self) -> List[Well]:
        return [w for w in self.wells if w.active and w.open_perforations()]

    def pore_volume(self) -> np.ndarray:
        pv = self.rock.porosity.values * self.geometry.volumes()
        if self.rock.net_to_gross is not None:
            pv = pv * self.rock.net_to_gross.values
        return pv

    # təxmini çatlama qradiyenti (bar/m) — vurucu BHP-nin üst həddi üçün
    FRACTURE_GRADIENT = 0.160
    # bundan dayaz modellərdə çatlama təxmini mənasızdır (sintetik/test modelləri)
    MINIMUM_DEPTH_FOR_FRACTURE_CHECK = 100.0

    def validate(self) -> list:
        """Yalnız BLOKLAYAN problemlərin siyahısı (geriyə uyğunluq üçün)."""
        return self.diagnose().messages(Severity.ERROR)

    def diagnose(self) -> DiagnosticReport:
        """Tam diaqnostika: xətalar + xəbərdarlıqlar.

        Bütün yoxlamalar BURADA — UI-də deyil. Əvvəllər bir hissəsi
        QMessageBox çağırışları ilə interfeysdə idi, yəni skriptdən
        istifadə edəndə heç işə düşmürdü.
        """
        report = DiagnosticReport()

        for message in self.rock.validate():
            report.error(message, "süxur")
        for message in self.scal_parameters.validate():
            report.error(message, "SCAL")
        for message in self.capillary_parameters.validate():
            report.error(message, "kapilyar")
        if self.pvt_table is not None:
            for message in self.pvt_table.validate():
                report.error(message, "PVT")

        self._check_initial_saturation(report)
        self._check_wells(report)
        self._check_faults(report)
        self._check_pvt_scal_ranges(report)
        return report

    # ------------------------------------------------------ yoxlamalar
    def _check_initial_saturation(self, report: DiagnosticReport) -> None:
        sw = self.initial_conditions.water_saturation
        scal = self.scal_parameters
        if not (scal.swc - 1e-9 <= sw <= 1.0 - scal.sor + 1e-9):
            report.error(
                "İlkin Sw hərəkətli doyumluluq intervalından kənardadır.",
                "ilkin şərtlər",
                f"Swc = {scal.swc:.3f} … 1−Sor = {1.0 - scal.sor:.3f}")

    def _check_faults(self, report: DiagnosticReport) -> None:
        seen = {}
        for fault in self.fault_references:
            for message in fault.validate(self.grid):
                report.error(message, fault.name)
            seen[fault.name] = seen.get(fault.name, 0) + 1
        for name, count in seen.items():
            if count > 1:
                report.warning(f"'{name}' adlı fault {count} dəfə təkrarlanır — "
                               f"çarpanları vurulacaq.", name)

    def _check_wells(self, report: DiagnosticReport) -> None:
        wells = self.active_wells()
        if not wells:
            report.error("Modeldə aktiv quyu yoxdur.", "quyular")
            return
        if not any(w.control.mode is ControlMode.BHP for w in wells):
            report.error(
                "Ən azı bir quyu BHP ilə idarə olunmalıdır — "
                "əks halda təzyiq səviyyəsi qeyri-müəyyəndir.", "quyular")

        seen = {}
        for well in wells:
            seen[well.name] = seen.get(well.name, 0) + 1
        for name, count in seen.items():
            if count > 1:
                report.warning(
                    f"'{name}' adı {count} quyuda təkrarlanır — debitlər "
                    f"hesabatda birləşdiriləcək.", name,
                    "Adları fərqləndir (INJ-1, INJ-2 …)")

        reference = self._reference_pressure()
        fracture = self._fracture_pressure()

        for well in wells:
            for perforation in well.open_perforations():
                if not (0 <= perforation.i < self.grid.nx
                        and 0 <= perforation.j < self.grid.ny
                        and 0 <= perforation.k < self.grid.nz):
                    report.error(
                        f"{well.name}: perforasiya grid-dən kənardadır "
                        f"(i={perforation.i}, j={perforation.j}, "
                        f"k={perforation.k + 1}).", well.name,
                        f"Grid: {self.grid.nx}×{self.grid.ny}×{self.grid.nz}")
            self._check_well_control(well, reference, fracture, report)

    def _check_well_control(self, well, reference, fracture,
                            report: DiagnosticReport) -> None:
        if well.control.mode is not ControlMode.BHP:
            return
        target = well.control.target

        if well.is_injector:
            if target <= reference:
                report.warning(
                    f"{well.name}: vurucu quyunun BHP-si ({target:.0f} bar) "
                    f"lay təzyiqindən ({reference:.0f} bar) yüksək deyil — "
                    f"quyu heç nə vurmayacaq.", well.name,
                    f"BHP-ni {reference:.0f} bardan yuxarı qaldır")
            elif fracture and target > fracture:
                report.warning(
                    f"{well.name}: BHP ({target:.0f} bar) təxmini çatlama "
                    f"təzyiqindən ({fracture:.0f} bar) yüksəkdir — real layda "
                    f"hidravlik yarılma riski var.", well.name,
                    "Çatlama modeli hesablanmır, nəticələr nikbin ola bilər")
        else:
            if target >= reference:
                report.warning(
                    f"{well.name}: hasilat quyusunun BHP-si ({target:.0f} bar) "
                    f"lay təzyiqindən ({reference:.0f} bar) aşağı deyil — "
                    f"quyu hasilat verməyəcək.", well.name,
                    f"BHP-ni {reference:.0f} bardan aşağı sal")
            elif (self.pvt_table is not None
                  and self.pvt_table.bubble_point > 0):
                # v69: qaz modelləşdirilmir, ona görə doyma təzyiqindən
                # aşağı BHP həmişə xəbərdarlıq doğurur.
                bubble = self.pvt_table.bubble_point
                if target < bubble:
                    report.warning(
                        f"{well.name}: BHP ({target:.0f} bar) doyma "
                        f"təzyiqindən ({bubble:.0f} bar) aşağıdır — quyudibində "
                        f"qaz ayrılacaq.", well.name,
                        "Qaz fazası modelləşdirilmir, nəticələr nikbin ola bilər")

    def _check_pvt_scal_ranges(self, report: DiagnosticReport) -> None:
        """Simulyasiya BAŞLAMAZDAN ƏVVƏL PVT/SCAL cədvəllərinin gözlənilən
        işləmə diapazonunu (ilkin təzyiq, quyu BHP hədəfləri, ilkin Sw)
        əhatə etdiyini yoxlayır — Phase 1 (giriş boru xətti).

        Newton həlqəsi bunu HƏR İTERASİYADA YOXLAMIR (performans/davranış
        riski, bax UNITS.md) — YALNIZ burada, model qurularkən bir dəfə.
        Yüngül kənarlaşma xəbərdarlıqdır (cədvəl sərhədə kəsir, nəticə
        dəyişmir); HƏDDİNDƏN ARTIQ kənarlaşma (adətən vahid qarışıqlığı
        əlaməti) sərt xətadır (bax `validate_query_range`).
        """
        if self.pvt_table is not None and self.pvt_table.size >= 2:
            pressures = [float(self.initial_conditions.datum_pressure)]
            pressures += [float(w.control.target) for w in self.active_wells()
                         if w.control.mode is ControlMode.BHP]
            result = validate_query_range(
                pressures, float(self.pvt_table.pressure[0]),
                float(self.pvt_table.pressure[-1]), "PVT təzyiq sorğusu")
            for message in result.errors:
                report.error(message, "PVT")
            for message in result.warnings:
                report.warning(message, "PVT")

        tables = getattr(self.scal_tables, "tables", None)
        if tables:
            sw = float(self.initial_conditions.water_saturation)
            for region, table in sorted(tables.items()):
                if table.sw.size < 2:
                    continue
                result = validate_query_range(
                    [sw], float(table.sw[0]), float(table.sw[-1]),
                    f"SCAL (region {region}) Sw sorğusu")
                for message in result.errors:
                    report.error(message, "SCAL")
                for message in result.warnings:
                    report.warning(message, "SCAL")

    def _reference_pressure(self) -> float:
        return float(self.initial_conditions.datum_pressure)

    def _fracture_pressure(self) -> Optional[float]:
        depths = self.geometry.cell_depths()
        mean_depth = float(np.mean(depths))
        if mean_depth < self.MINIMUM_DEPTH_FOR_FRACTURE_CHECK:
            return None
        return mean_depth * self.FRACTURE_GRADIENT

    def summary(self) -> dict:
        return {
            "name": self.name,
            "cells": self.ncell,
            "wells": len(self.active_wells()),
            "regions": int(self.regions.ids.size),
            "faults": len(self.fault_references),
            "horizons": len(self.horizon_references),
            "units": self.units.name,
            "pvt": self.pvt_table.source if self.pvt_table else "statik (PVT yoxdur)",
            "from": self.source_geological_model,
        }
