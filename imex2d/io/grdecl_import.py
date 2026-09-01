"""GRDECL -> GeologicalModel çevirməsi.

Oxuma (`grdecl.py`) və çevirmə (bu fayl) qəsdən ayrıdır: birincisi
faylın sintaksisini bilir, ikincisi domain modelini. Beləliklə yeni
format əlavə edəndə yalnız oxuyucu yazılır.
"""

from __future__ import annotations
from typing import Optional, Tuple

import numpy as np

from ..domain.diagnostics import DiagnosticReport
from ..domain.geological_model import GeologicalModel
from ..domain.geometry import CellGeometry
from ..domain.grid import CartesianGrid
from ..domain.properties import PropertyMap
from ..domain.structure import RegionSet
from ..domain.unit_conversions import to_engine_units
from ..logging_setup import get_logger
from .grdecl import GrdeclDeck, GrdeclError

LOG = get_logger(__name__)

DEFAULT_CELL_SIZE = 100.0
DEFAULT_THICKNESS = 10.0
DEFAULT_TOP_DEPTH = 2000.0

#: Eclipse-in hər vahid sistemində UZUNLUĞUN native vahidi (METRIC/FIELD/
#: LAB-ın RUNSPEC-də bəyan olunan üç dəsti). Keçiricilik (mD) VƏ lözlük
#: (cP) HƏR ÜÇÜNDƏ EYNİDİR (Eclipse-in özünəməxsus konvensiyası) — ona
#: görə bu cədvəldə YOXDUR, çevrilməyə EHTİYAC yoxdur.
LENGTH_UNIT_BY_DECK_SYSTEM = {"METRIC": "m", "FIELD": "ft", "LAB": "cm"}


class GrdeclImporter:
    """Oxunmuş deck-dən geoloji model qurur."""

    def build(self, deck: GrdeclDeck,
              report: Optional[DiagnosticReport] = None,
              name: str = "GRDECL modeli") -> GeologicalModel:
        report = report if report is not None else DiagnosticReport()
        unit_system = self._declared_unit_system(deck)
        length_unit = LENGTH_UNIT_BY_DECK_SYSTEM[unit_system]
        self._check_declared_units(unit_system, report)
        nx, ny, nz = deck.dimensions
        grid = CartesianGrid(nx, ny, nz)

        geometry = self._geometry(deck, grid, report, length_unit)
        model = GeologicalModel(name=name, grid=grid, geometry=geometry,
                                regions=self._regions(deck, grid, report),
                                coordinate_system="ECLIPSE")

        self._add_properties(deck, model, grid, report)
        self._check_inactive_cells(deck, report)

        issues = model.validate()
        if issues:
            raise GrdeclError("Qurulan model natamamdır: " + "; ".join(issues))
        return model

    @staticmethod
    def _declared_unit_system(deck: GrdeclDeck) -> str:
        """RUNSPEC-də bəyan edilmiş vahid sistemi. Heç biri yoxdursa,
        Eclipse-in öz defoltu (METRIC) qəbul edilir — bu, bu oxucunun
        DEFOLT davranışı ilə üst-üstə düşdüyü üçün YEGANƏ sükutla
        qəbul edilən haldır."""
        if "FIELD" in deck.keywords_seen:
            return "FIELD"
        if "LAB" in deck.keywords_seen:
            return "LAB"
        return "METRIC"

    @staticmethod
    def _check_declared_units(unit_system: str, report: DiagnosticReport) -> None:
        """`FIELD`/`LAB` bəyan edilibsə NƏYİN çevrildiyini, NƏYİN YOX
        olduğunu AÇIQ bildirir — heç vaxt sükutla METRIC kimi oxumur
        (bax Phase 1 audit: əvvəllər BÜTÜN dəyərlər sükutla METRIC
        sayılırdı, bu, ft dəyərinin metr kimi oxunması demək idi).

        Uzunluq (DX/DY/DZ/TOPS/COORD/ZCORN): FIELD→ft, LAB→cm-dən "m"-ə
        FAKTİKİ ÇEVRİLİR (bax `_uniform`/`_top_surface`/`_from_corner_point`).
        Keçiricilik (PERMX/PERMY/PERMZ): ÇEVRİLMİR — bu, çevirmənin
        YAZILMAMASI deyil, Eclipse-in ÜÇ vahid sistemində də (METRIC/
        FIELD/LAB) keçiriciliyin HƏMİŞƏ mD olması faktına görə EHTİYAC
        YOXLUĞUDUR. Təzyiq/PVT açar sözləri bu oxucuda HƏLƏ ÜMUMİYYƏTLƏ
        YOXDUR — gələcəkdə əlavə ediləndə FIELD/LAB üçün AYRICA çevirmə
        YAZILMALIDIR, bu funksiya onu HƏLƏ ETMİR.
        """
        if unit_system == "METRIC":
            return
        label = {"FIELD": "FIELD (ft/psi/mD)", "LAB": "LAB (sm/atm/mD)"}[unit_system]
        report.warning(
            f"Fayl {label} vahid sistemini bəyan edir. Uzunluq (DX/DY/DZ/TOPS/"
            "COORD/ZCORN) mühərrik vahidinə (m) ÇEVRİLDİ. Keçiricilik (PERMX/"
            "PERMY/PERMZ) çevrilmədi — Eclipse-də bu üç sistemdə də mD-dir, "
            "çevrilməyə EHTİYAC YOXDUR. Bu oxucu təzyiq/PVT açar sözlərini "
            "HƏLƏ oxumur.", "GRDECL",
            "Gələcəkdə uzunluqdan başqa yeni açar söz əlavə edilsə, onun da "
            "vahidi FIELD/LAB üçün yoxlanılmalıdır — bax UNITS.md.")

    # ═══════════════════════════════════════════════════ həndəsə
    def _geometry(self, deck: GrdeclDeck, grid: CartesianGrid,
                  report: DiagnosticReport, length_unit: str = "m") -> CellGeometry:
        if deck.has("COORD") and deck.has("ZCORN"):
            return self._from_corner_point(deck, grid, report, length_unit)
        return self._from_block_centred(deck, grid, report, length_unit)

    def _from_block_centred(self, deck: GrdeclDeck, grid: CartesianGrid,
                            report: DiagnosticReport, length_unit: str = "m"
                            ) -> CellGeometry:
        """DX/DY/DZ/TOPS — birbaşa uyğunluq (deck vahidindən "m"-ə çevrilir)."""
        dx = self._uniform(deck, "DX", DEFAULT_CELL_SIZE, report, length_unit)
        dy = self._uniform(deck, "DY", dx, report, length_unit)
        dz = self._uniform(deck, "DZ", DEFAULT_THICKNESS, report, length_unit)

        surface, top_depth = self._top_surface(deck, grid, dz, report, length_unit)
        return CellGeometry(grid=grid, dx=dx, dy=dy, dz=dz,
                            top_depth=top_depth, top_depth_map=surface)

    def _from_corner_point(self, deck: GrdeclDeck, grid: CartesianGrid,
                           report: DiagnosticReport, length_unit: str = "m"
                           ) -> CellGeometry:
        """COORD/ZCORN — bərabər bloka APPROKSİMASİYA.

        `CellGeometry` hazırda dəyişkən hüceyrə ölçüsü saxlamır, ona görə
        orta ölçü hesablanır. Bu, dürüst şəkildə xəbərdarlıq tələb edir —
        nəticələr orijinal həndəsə ilə tam üst-üstə düşməyəcək.
        """
        coord = deck.get("COORD")
        zcorn = deck.get("ZCORN")
        if length_unit != "m":
            coord = to_engine_units(coord, length_unit, "length")
            zcorn = to_engine_units(zcorn, length_unit, "length")
        nx, ny, nz = grid.nx, grid.ny, grid.nz

        expected_coord = 6 * (nx + 1) * (ny + 1)
        expected_zcorn = 8 * nx * ny * nz
        if coord.size != expected_coord or zcorn.size != expected_zcorn:
            raise GrdeclError(
                f"COORD/ZCORN ölçüsü uyğun gəlmir: {coord.size}/{zcorn.size}, "
                f"gözlənilən {expected_coord}/{expected_zcorn}.")

        pillars = coord.reshape((ny + 1), (nx + 1), 6)
        x_pillars = pillars[:, :, 0]
        y_pillars = pillars[:, :, 1]
        dx = float(np.mean(np.abs(np.diff(x_pillars, axis=1))))
        dy = float(np.mean(np.abs(np.diff(y_pillars, axis=0))))

        corners = zcorn.reshape(2 * nz, 2 * ny, 2 * nx)
        cell_tops = corners[0::2].reshape(nz, 2 * ny, 2 * nx)
        cell_bottoms = corners[1::2].reshape(nz, 2 * ny, 2 * nx)
        thickness = float(np.mean(np.abs(cell_bottoms - cell_tops)))

        top_layer = cell_tops[0].reshape(ny, 2, nx, 2).mean(axis=(1, 3))

        report.warning(
            "Fayl corner-point həndəsəsindədir. Model bərabər ölçülü "
            f"bloklara approksimasiya olundu (DX≈{dx:.1f}, DY≈{dy:.1f}, "
            f"DZ≈{thickness:.1f} m).", "GRDECL",
            "Nəticələr orijinal həndəsə ilə tam üst-üstə düşməyəcək")

        return CellGeometry(grid=grid, dx=max(dx, 1e-3), dy=max(dy, 1e-3),
                            dz=max(thickness, 1e-3),
                            top_depth=float(np.mean(top_layer)),
                            top_depth_map=top_layer.ravel())

    def _top_surface(self, deck: GrdeclDeck, grid: CartesianGrid, dz: float,
                     report: DiagnosticReport, length_unit: str = "m"
                     ) -> Tuple[Optional[np.ndarray], float]:
        tops = deck.get("TOPS")
        if tops is None:
            report.info(f"TOPS verilməyib — tavan {DEFAULT_TOP_DEPTH:.0f} m "
                        f"qəbul edildi.", "GRDECL")
            return None, DEFAULT_TOP_DEPTH
        if length_unit != "m":
            tops = to_engine_units(tops, length_unit, "length")

        areal = grid.nx * grid.ny
        surface = tops[:areal] if tops.size >= areal else None
        if surface is None:
            return None, float(np.nanmean(tops))
        if np.allclose(surface, surface[0], equal_nan=True):
            return None, float(surface[0])
        return surface.copy(), float(np.nanmin(surface))

    @staticmethod
    def _uniform(deck: GrdeclDeck, keyword: str, fallback: float,
                 report: DiagnosticReport, length_unit: str = "m") -> float:
        """Hüceyrə ölçüsü massivdirsə, orta qiymət götürülür (deck-in öz
        uzunluq vahidindən — `length_unit` — mühərrik vahidinə, "m",
        ÇEVRİLİR). `fallback` artıq "m"-dədir (verilməmiş DX/DY/DZ üçün
        bu modulun özünün defoltu), ona görə HEÇ vaxt çevrilmir."""
        values = deck.get(keyword)
        if values is None:
            report.info(f"{keyword} verilməyib — {fallback:g} m "
                        f"qəbul edildi.", "GRDECL")
            return float(fallback)
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            return float(fallback)
        if not np.allclose(finite, finite[0]):
            report.warning(
                f"{keyword} hüceyrədən hüceyrəyə dəyişir "
                f"({finite.min():.1f}–{finite.max():.1f} {length_unit}); "
                f"orta qiymət {finite.mean():.1f} {length_unit} işlədilir.", "GRDECL",
                "Dəyişkən hüceyrə ölçüsü hələ dəstəklənmir")
        mean_value = float(finite.mean())
        if length_unit != "m":
            mean_value = to_engine_units(mean_value, length_unit, "length")
        return mean_value

    # ═══════════════════════════════════════════════════ xassələr
    def _add_properties(self, deck: GrdeclDeck, model: GeologicalModel,
                        grid: CartesianGrid,
                        report: DiagnosticReport) -> None:
        porosity = deck.get("PORO")
        if porosity is None:
            raise GrdeclError("Faylda PORO massivi yoxdur.")
        permx = deck.get("PERMX")
        if permx is None:
            raise GrdeclError("Faylda PERMX massivi yoxdur.")

        model.add_property(PropertyMap.from_array(
            "PORO", self._clean(porosity, 0.01, 0.6), grid.ncell))
        model.add_property(PropertyMap.from_array(
            "PERMX", self._clean(permx, 1e-4, 1e6), grid.ncell, "mD"))

        for keyword, fallback_key, factor in (("PERMY", "PERMX", 1.0),
                                              ("PERMZ", "PERMX", 0.1)):
            values = deck.get(keyword)
            if values is None:
                values = deck.get(fallback_key) * factor
                report.info(f"{keyword} verilməyib — {fallback_key}×{factor:g} "
                            f"işlədildi.", "GRDECL")
            model.add_property(PropertyMap.from_array(
                keyword, self._clean(values, 1e-6, 1e6), grid.ncell, "mD"))

        net_to_gross = deck.get("NTG")
        if net_to_gross is not None:
            model.add_property(PropertyMap.from_array(
                "NTG", self._clean(net_to_gross, 0.0, 1.0), grid.ncell))

    @staticmethod
    def _clean(values: np.ndarray, minimum: float, maximum: float
               ) -> np.ndarray:
        """NaN (Eclipse-in `n*` defoltu) və hədləri təmizləyir."""
        cleaned = np.asarray(values, float).copy()
        finite = cleaned[np.isfinite(cleaned)]
        replacement = float(np.median(finite)) if finite.size else minimum
        cleaned[~np.isfinite(cleaned)] = replacement
        return np.clip(cleaned, minimum, maximum)

    # ═══════════════════════════════════════════════════ regionlar
    @staticmethod
    def _regions(deck: GrdeclDeck, grid: CartesianGrid,
                 report: DiagnosticReport) -> RegionSet:
        satnum = deck.get("SATNUM")
        if satnum is None:
            return RegionSet.single(grid.ncell)
        identifiers = np.nan_to_num(satnum, nan=1.0).astype(int)
        names = {int(value): f"SATNUM-{int(value)}"
                 for value in np.unique(identifiers)}
        report.info(f"SATNUM oxundu: {len(names)} region.", "GRDECL")
        return RegionSet(
            PropertyMap.from_array("REGION_ID", identifiers, grid.ncell),
            names)

    @staticmethod
    def _check_inactive_cells(deck: GrdeclDeck,
                              report: DiagnosticReport) -> None:
        """ACTNUM hələ dəstəklənmir — susmaq təhlükəlidir."""
        actnum = deck.get("ACTNUM")
        if actnum is None:
            return
        inactive = int(np.sum(np.nan_to_num(actnum, nan=1.0) < 0.5))
        if inactive:
            report.warning(
                f"Faylda {inactive} qeyri-aktiv hüceyrə var (ACTNUM = 0). "
                f"Hazırkı model bütün hüceyrələri aktiv sayır.", "GRDECL",
                "Həcm və ehtiyat hesabları böyük çıxacaq")
