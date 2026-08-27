"""İŞ AXINININ YENİ ADDIMI: quyu məlumatı → geoloji model.

    Karotaj interpretasiyası (xaricdə)
        ↓  CSV
    WellDataset
        ↓  interpolyasiya (IDW / Kriging / ən yaxın qonşu)
    GeologicalModel  →  ReservoirModel  →  Simulyasiya

Bu qat interpolyasiya alqoritmini TANIMIR — yalnız IPropertyInterpolator
interfeysini bilir. Alqoritm konstruktora inject edilir.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Optional

import numpy as np

from ..domain.geological_model import GeologicalModel
from ..domain.geometry import CellGeometry
from ..domain.grid import CartesianGrid
from ..domain.properties import PropertyMap
from ..domain.structure import RegionSet
from ..domain.well_data import WellDataset
from ..geology.interpolation import interpolate_property
from ..interfaces.interpolation import IPropertyInterpolator


@dataclass
class PropertyRule:
    """Bir xassənin necə interpolyasiya olunacağı."""
    target: str                      # grid açarı: PORO, PERMX, NTG…
    log_transform: bool = False
    minimum: Optional[float] = None
    maximum: Optional[float] = None


DEFAULT_RULES: Dict[str, PropertyRule] = {
    "PORO": PropertyRule("PORO", False, 0.01, 0.45),
    "PERMX": PropertyRule("PERMX", True, 0.01, 1e5),
    "PERMY": PropertyRule("PERMY", True, 0.01, 1e5),
    "PERMZ": PropertyRule("PERMZ", True, 0.001, 1e5),
    "NTG": PropertyRule("NTG", False, 0.01, 1.0),
    "SW": PropertyRule("SW", False, 0.0, 1.0),
    "VSH": PropertyRule("VSH", False, 0.0, 1.0),
}


@dataclass
class GeologicalGridSpec:
    """Interpolyasiyanın aparılacağı grid."""
    nx: int = 41
    ny: int = 41
    nz: int = 1
    dx: float = 20.0
    dy: float = 20.0
    dz: float = 10.0
    top_depth: float = 2000.0
    dip_x: float = 0.0
    dip_y: float = 0.0


@dataclass
class InterpolationReport:
    """Nəyin necə hesablandığı — istifadəçiyə göstərilir."""
    method: str = ""
    entries: list = field(default_factory=list)

    def add(self, target: str, source: str, log_transform: bool, values: np.ndarray):
        self.entries.append({
            "target": target, "source": source, "log": log_transform,
            "min": float(values.min()), "max": float(values.max()),
            "mean": float(values.mean())})

    def as_text(self) -> str:
        lines = [f"Üsul: {self.method}"]
        for entry in self.entries:
            lines.append(
                f"  {entry['target']:<6} ← {entry['source']:<6} "
                f"{'(log)' if entry['log'] else '     '}  "
                f"min {entry['min']:.4g}  orta {entry['mean']:.4g}  "
                f"maks {entry['max']:.4g}")
        return "\n".join(lines)


class WellBasedGeologicalModelBuilder:
    """Quyu nöqtələrindən grid xassələri qurur."""

    def __init__(self, interpolator: IPropertyInterpolator,
                 rules: Optional[Dict[str, PropertyRule]] = None):
        self.interpolator = interpolator
        self.rules = dict(DEFAULT_RULES)
        if rules:
            self.rules.update(rules)

    # ---------------------------------------------------------- public
    def build(self, dataset: WellDataset, spec: GeologicalGridSpec,
              ky_over_kx: float = 1.0, kv_over_kh: float = 0.1,
              name: str = "Quyu məlumatından geoloji model"):
        issues = dataset.validate()
        if issues:
            raise ValueError("Quyu məlumatı yararsızdır: " + "; ".join(issues))

        grid = CartesianGrid(spec.nx, spec.ny, spec.nz)
        geometry = CellGeometry(grid, spec.dx, spec.dy, spec.dz,
                                top_depth=spec.top_depth,
                                top_depth_map=self._surface(grid, spec))
        model = GeologicalModel(name=name, grid=grid, geometry=geometry,
                                regions=RegionSet.single(grid.ncell))
        report = InterpolationReport(method=self.interpolator.describe())

        targets = self._cell_centres(grid, spec)
        available = dataset.property_names()

        for source in available:
            rule = self.rules.get(source, PropertyRule(source))
            values = self._interpolate_volume(dataset, source, rule, targets, grid)
            model.add_property(PropertyMap.from_array(rule.target, values,
                                                      grid.ncell))
            report.add(rule.target, source, rule.log_transform, values)

        self._fill_missing_permeability(model, grid, ky_over_kx, kv_over_kh, report)
        issues = model.validate()
        if issues:
            raise ValueError("Qurulan geoloji model natamamdır: " + "; ".join(issues))
        return model, report

    # -------------------------------------------------------- internal
    @staticmethod
    def _cell_centres(grid: CartesianGrid, spec: GeologicalGridSpec) -> np.ndarray:
        x = (np.arange(grid.nx) + 0.5) * spec.dx
        y = (np.arange(grid.ny) + 0.5) * spec.dy
        yy, xx = np.meshgrid(y, x, indexing="ij")
        return np.column_stack([xx.ravel(), yy.ravel()])

    @staticmethod
    def _surface(grid: CartesianGrid, spec: GeologicalGridSpec):
        if abs(spec.dip_x) < 1e-12 and abs(spec.dip_y) < 1e-12:
            return None
        i = np.arange(grid.nx)
        j = np.arange(grid.ny)
        jj, ii = np.meshgrid(j, i, indexing="ij")
        return spec.top_depth + ii * spec.dip_x + jj * spec.dip_y

    def _interpolate_volume(self, dataset, source, rule, targets, grid) -> np.ndarray:
        """Hər təbəqə üçün areal interpolyasiya, sonra həcmə yığılır."""
        layers = []
        for k in range(grid.nz):
            layer = k if dataset.is_layered() else None
            points, values = dataset.points(source, layer)
            if values.size == 0:                      # bu təbəqədə məlumat yoxdur
                points, values = dataset.points(source, None)
            if values.size == 0:
                raise ValueError(f"'{source}' üçün nöqtə tapılmadı.")
            layers.append(interpolate_property(
                self.interpolator, points, values, targets,
                log_transform=rule.log_transform,
                minimum=rule.minimum, maximum=rule.maximum))
        return np.concatenate(layers)

    @staticmethod
    def _fill_missing_permeability(model, grid, ky_over_kx, kv_over_kh, report):
        """PERMY/PERMZ verilməyibsə anizotropluq əmsalları ilə qurulur."""
        if "PERMX" not in model.property_maps:
            return
        permx = model.property_maps["PERMX"].values
        for key, factor in (("PERMY", ky_over_kx), ("PERMZ", kv_over_kh)):
            if key in model.property_maps:
                continue
            values = permx * factor
            model.add_property(PropertyMap.from_array(key, values, grid.ncell, "mD"))
            report.add(key, "PERMX", False, values)
