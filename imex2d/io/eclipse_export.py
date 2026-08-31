"""Eclipse `.DATA` deck yazıcısı.

Məqsəd: modelin CMG/Eclipse/tNavigator-a ötürülməsi. Yazılan deck
tam işlək olmaya bilər (hər simulyatorun öz xüsusiyyətləri var), lakin
strukturu düzgündür və əl ilə tamamlanmağa hazırdır.

Bölmələr: RUNSPEC, GRID, PROPS, SOLUTION, SUMMARY, SCHEDULE.

VAHİDLƏR: METRIC (m, bar, m3, gün) — modelin öz sistemi ilə eynidir,
ona görə heç bir çevirmə aparılmır və çevirmə səhvi riski yoxdur.
"""

from __future__ import annotations

import os
from datetime import datetime
from typing import List, Optional

import numpy as np

from ..domain.reservoir_model import ReservoirModel
from ..domain.wells import ControlMode, WellType
from ..logging_setup import get_logger
from ..version import VERSION

LOG = get_logger(__name__)

VALUES_PER_LINE = 8


def _compress(values: np.ndarray, decimals: int = 5) -> List[str]:
    """Ardıcıl eyni dəyərləri `n*value` şəklində sıxır.

    200 000 hüceyrəli modeldə bu, faylı onlarla dəfə kiçildir.
    """
    rounded = np.round(np.asarray(values, float), decimals)
    tokens: List[str] = []
    count = 1
    for index in range(1, rounded.size + 1):
        if index < rounded.size and rounded[index] == rounded[index - 1]:
            count += 1
            continue
        value = rounded[index - 1]
        text = f"{value:g}"
        tokens.append(f"{count}*{text}" if count > 1 else text)
        count = 1
    return tokens


def _array_block(keyword: str, values: np.ndarray, decimals: int = 5) -> str:
    tokens = _compress(values, decimals)
    lines = [keyword]
    for start in range(0, len(tokens), VALUES_PER_LINE):
        lines.append("  " + " ".join(tokens[start:start + VALUES_PER_LINE]))
    lines.append("/")
    return "\n".join(lines)


class EclipseDeckWriter:
    """ReservoirModel -> Eclipse `.DATA` mətn faylı."""

    def __init__(self, end_time: float = 1500.0, report_steps: int = 20):
        self.end_time = float(end_time)
        self.report_steps = int(report_steps)

    def write(self, model: ReservoirModel, path: str) -> str:
        text = self.render(model)
        directory = os.path.dirname(os.path.abspath(path))
        if directory:
            os.makedirs(directory, exist_ok=True)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
        LOG.info("Eclipse deck yazıldı: %s (%d hüceyrə)", path, model.ncell)
        return path

    def render(self, model: ReservoirModel) -> str:
        return "\n\n".join([
            self._header(model),
            self._runspec(model),
            self._grid(model),
            self._props(model),
            self._solution(model),
            self._summary(),
            self._schedule(model),
        ]) + "\n"

    # ══════════════════════════════════════════════════════ bölmələr
    @staticmethod
    def _header(model: ReservoirModel) -> str:
        stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
        return "\n".join([
            "-- " + "=" * 70,
            f"-- IMEX-2D v{VERSION} tərəfindən yaradıldı  ·  {stamp}",
            f"-- Model: {model.name}",
            f"-- Grid:  {model.grid.nx} x {model.grid.ny} x {model.grid.nz}"
            f"  =  {model.ncell} hüceyrə",
            "--",
            "-- QEYD: bu deck avtomatik yaradılıb. Hər simulyatorun öz",
            "-- xüsusiyyətləri var — işə salmazdan əvvəl yoxlayın.",
            "-- " + "=" * 70,
        ])

    def _runspec(self, model: ReservoirModel) -> str:
        grid = model.grid
        wells = model.active_wells()
        max_perforations = max((len(w.open_perforations()) for w in wells),
                               default=1)
        return "\n".join([
            "RUNSPEC",
            "",
            "TITLE",
            f"  {model.name}",
            "/",
            "",
            "DIMENS",
            f"  {grid.nx} {grid.ny} {grid.nz} /",
            "",
            "-- iki fazalı: su və neft",
            "OIL",
            "WATER",
            "",
            "METRIC",
            "",
            "START",
            "  1 'JAN' 2020 /",
            "",
            "WELLDIMS",
            f"  {max(len(wells), 1)} {max_perforations} "
            f"{max(len(wells), 1)} {max(len(wells), 1)} /",
        ])

    def _grid(self, model: ReservoirModel) -> str:
        grid, geometry = model.grid, model.geometry
        rock = model.rock
        depths = geometry.cell_depths().reshape(grid.shape)
        dz_cell = geometry.dz_per_cell()
        tops = depths[0] - dz_cell.reshape(grid.shape)[0] * 0.5

        parts = [
            "GRID",
            "",
            "-- bərabər ölçülü bloklar",
            _array_block("DX", np.full(model.ncell, geometry.dx), 3),
            "",
            _array_block("DY", np.full(model.ncell, geometry.dy), 3),
            "",
            _array_block("DZ", dz_cell, 3),
            "",
            _array_block("TOPS", tops.ravel(), 3),
            "",
            _array_block("PORO", rock.porosity.values, 5),
            "",
            _array_block("PERMX", rock.permx.values, 4),
            "",
            _array_block("PERMY", rock.permy.values, 4),
        ]
        if rock.permz is not None:
            parts += ["", _array_block("PERMZ", rock.permz.values, 4)]
        if rock.net_to_gross is not None:
            parts += ["", _array_block("NTG", rock.net_to_gross.values, 4)]

        regions = model.regions
        if regions is not None and regions.ids.size > 1:
            parts += ["", _array_block("SATNUM",
                                       regions.region_id.values, 0)]
        parts += ["", "INIT"]
        return "\n".join(parts)

    def _props(self, model: ReservoirModel) -> str:
        scal = model.scal_parameters
        fluids = model.fluids
        saturations = np.linspace(scal.swc, 1.0 - scal.sor, 12)

        capillary = model.capillary_parameters
        pc_values = np.zeros_like(saturations)
        if capillary.enabled:
            from ..simulation.capillary import BrooksCoreyCapillaryProvider
            pc_values = BrooksCoreyCapillaryProvider(capillary, scal).pcow(
                saturations)

        rows = ["SWOF", "-- Sw        krw        kro        Pc"]
        for sw, pc in zip(saturations, pc_values):
            rows.append(f"  {sw:8.5f} {float(scal.krw(sw)):10.6f} "
                        f"{float(scal.kro(sw)):10.6f} {pc:10.5f}")
        rows.append("/")

        parts = ["PROPS", "", "\n".join(rows), ""]
        parts.append(self._pvt_tables(model, fluids))
        parts += ["", "ROCK", f"  {model.initial_conditions.datum_pressure:.1f} "
                  f"{model.rock.compressibility:.3e} /"]
        parts += ["", "DENSITY",
                  f"  {fluids.oil_density:.1f} {fluids.water_density:.1f} "
                  f"1.0 /"]
        return "\n".join(parts)

    @staticmethod
    def _pvt_tables(model: ReservoirModel, fluids) -> str:
        reference = model.initial_conditions.datum_pressure
        table = model.pvt_table

        if table is None:
            return "\n".join([
                "PVCDO",
                f"  {reference:.1f} {fluids.oil_fvf:.5f} "
                f"{fluids.oil_compressibility:.4e} {fluids.oil_viscosity:.4f} "
                f"0.0 /",
                "",
                "PVTW",
                f"  {reference:.1f} {fluids.water_fvf:.5f} "
                f"{fluids.water_compressibility:.4e} "
                f"{fluids.water_viscosity:.4f} 0.0 /",
            ])

        rows = ["PVDO", "--  P, bar      Bo        muo"]
        for pressure, fvf, viscosity in zip(table.pressure, table.oil_fvf,
                                            table.oil_viscosity):
            rows.append(f"  {pressure:9.3f} {fvf:9.5f} {viscosity:9.5f}")
        rows.append("/")
        rows += ["", "PVTW",
                 f"  {reference:.1f} {float(np.interp(reference, table.pressure, table.water_fvf)):.5f} "
                 f"{fluids.water_compressibility:.4e} "
                 f"{float(np.interp(reference, table.pressure, table.water_viscosity)):.4f} 0.0 /"]
        return "\n".join(rows)

    def _solution(self, model: ReservoirModel) -> str:
        initial = model.initial_conditions
        depths = model.geometry.cell_depths()
        datum = initial.datum_depth or float(np.min(depths))
        contact = (initial.oil_water_contact
                   if initial.oil_water_contact is not None
                   else float(np.max(depths)) + 100.0)
        return "\n".join([
            "SOLUTION",
            "",
            "EQUIL",
            "--  datum    P       OWC    Pc(OWC)  GOC  Pc(GOC)",
            f"  {datum:8.2f} {initial.datum_pressure:8.2f} "
            f"{contact:8.2f}  0.0  0.0  0.0  1 /",
        ])

    @staticmethod
    def _summary() -> str:
        return "\n".join([
            "SUMMARY",
            "",
            "FOPR", "FWPR", "FWIR", "FOPT", "FWPT", "FPR", "FWCT",
            "WBHP", "  / ",
        ])

    def _schedule(self, model: ReservoirModel) -> str:
        parts = ["SCHEDULE", "", "WELSPECS"]
        for well in model.active_wells():
            perforations = well.open_perforations()
            if not perforations:
                continue
            phase = "WATER" if well.well_type is WellType.INJECTOR else "OIL"
            parts.append(f"  '{well.name}' 'G1' {perforations[0].i + 1} "
                         f"{perforations[0].j + 1} 1* '{phase}' /")
        parts.append("/")

        parts += ["", "COMPDAT"]
        for well in model.active_wells():
            for perforation in well.open_perforations():
                parts.append(
                    f"  '{well.name}' {perforation.i + 1} "
                    f"{perforation.j + 1} {perforation.k + 1} "
                    f"{perforation.k + 1} 'OPEN' 2* {well.radius * 2:.4f} "
                    f"3* {perforation.skin:.2f} /")
        parts.append("/")

        producers = [w for w in model.active_wells()
                     if w.well_type is WellType.PRODUCER]
        injectors = [w for w in model.active_wells()
                     if w.well_type is WellType.INJECTOR]

        if producers:
            parts += ["", "WCONPROD"]
            for well in producers:
                if well.control.mode is ControlMode.BHP:
                    parts.append(f"  '{well.name}' 'OPEN' 'BHP' 5* "
                                 f"{well.control.target:.2f} /")
                else:
                    parts.append(f"  '{well.name}' 'OPEN' 'LRAT' 2* "
                                 f"{abs(well.control.target):.2f} 2* 1.0 /")
            parts.append("/")

        if injectors:
            parts += ["", "WCONINJE"]
            for well in injectors:
                if well.control.mode is ControlMode.BHP:
                    parts.append(f"  '{well.name}' 'WATER' 'OPEN' 'BHP' 2* "
                                 f"{well.control.target:.2f} /")
                else:
                    parts.append(f"  '{well.name}' 'WATER' 'OPEN' 'RATE' "
                                 f"{abs(well.control.target):.2f} 1* "
                                 f"1000.0 /")
            parts.append("/")

        step = self.end_time / max(self.report_steps, 1)
        parts += ["", "TSTEP",
                  f"  {self.report_steps}*{step:.4f} /", "", "END"]
        return "\n".join(parts)
