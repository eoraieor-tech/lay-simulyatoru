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
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from ..domain.geological_model import GeologicalModel
from ..domain.geometry import CellGeometry
from ..domain.grid import CartesianGrid
from ..domain.properties import PropertyMap
from ..domain.structure import RegionSet
from ..domain.well_data import WellDataset
from ..geology.cross_validation import CrossValidationResult, k_fold, leave_one_out
from ..geology.interpolation import interpolate_property
from ..interfaces.interpolation import IPropertyInterpolator

#: `cross_validate_all`-un baxdığı xassələr — PORO və PERM* istiqamətləri.
_CROSS_VALIDATED_PROPERTIES = ("PORO", "PERMX", "PERMY", "PERMZ")


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
    dz: Union[float, Sequence[float]] = 10.0
    top_depth: float = 2000.0
    dip_x: float = 0.0
    dip_y: float = 0.0


@dataclass
class InterpolationReport:
    """Nəyin necə hesablandığı — istifadəçiyə göstərilir."""
    method: str = ""
    entries: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    def add(self, target: str, source: str, log_transform: bool, values: np.ndarray):
        self.entries.append({
            "target": target, "source": source, "log": log_transform,
            "min": float(values.min()), "max": float(values.max()),
            "mean": float(values.mean())})

    def warn(self, message: str) -> None:
        self.warnings.append(message)

    def as_text(self) -> str:
        lines = [f"Üsul: {self.method}"]
        for entry in self.entries:
            lines.append(
                f"  {entry['target']:<6} ← {entry['source']:<6} "
                f"{'(log)' if entry['log'] else '     '}  "
                f"min {entry['min']:.4g}  orta {entry['mean']:.4g}  "
                f"maks {entry['max']:.4g}")
        for message in self.warnings:
            lines.append(f"  ⚠ {message}")
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
              name: str = "Quyu məlumatından geoloji model",
              allow_cross_layer_fallback: bool = False):
        """`allow_cross_layer_fallback` — bir K-təbəqəsində (dataset laylı
        olanda) heç bir quyu nöqtəsi yoxdursa nə edilsin.

        Defolt (`False`): AÇIQ XƏTA atılır — sükutla başqa laylardan
        (məs. lay 1-3) məlumat "sızdırılmır". `True` verilsə köhnə davranış
        (bütün laylardan hovuzlanmış nöqtələrlə kriging) işə düşür, amma
        `report.warnings`-ə açıq xəbərdarlıq yazılır ki, istifadəçi bunun
        bilərəkdən ekstrapolyasiya olduğunu bilsin (bax M0 sınaqları,
        `tests/test_layer_aware_kriging_leak.py`).
        """
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
            values = self._interpolate_volume(dataset, source, rule, targets, grid,
                                              allow_cross_layer_fallback, report, geometry)
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

    @staticmethod
    def _sample_depth(sample, current_k: int, layer_mean_depth: np.ndarray) -> float:
        """Bir nümunənin Z-si (3D/anizotrop kriging üçün).

        Ölçülmüş `sample.depth` varsa (CSV `depth`/`md`/`tvd` sütunu)
        HƏMİŞƏ ona üstünlük verilir — bu, əsl ölçmədir. Yoxdursa, nümunə
        öz layına aid olan orta dərinliklə (`layer_mean_depth[layer]`)
        əvəz olunur; laysız (hər yerə aid) nümunə isə SORĞU olunan K-nın
        öz dərinliyini alır ki, "hər yerə aiddir" mənası qorunsun (heç bir
        süni şaquli məsafə yaranmır). Bu, fiziki dəyər UYDURMAQ deyil —
        yalnız grid həndəsəsindən artıq məlum olan dərinliyi işlədir.
        """
        if sample.depth is not None:
            return float(sample.depth)
        source_layer = sample.layer if sample.layer is not None else current_k
        return float(layer_mean_depth[source_layer])

    def _interpolate_volume(self, dataset, source, rule, targets, grid,
                            allow_cross_layer_fallback, report, geometry) -> np.ndarray:
        """Hər təbəqə üçün areal (və Kriging 3D dəstəkləyirsə — X,Y,Z)
        interpolyasiya, sonra həcmə yığılır.

        `dataset.samples_for(source, layer=k)` düzgün süzür: laya bağlı
        olmayan (`sample.layer is None`) nümunələr HƏR K üçün daxil
        edilir, yalnız k-ya bağlı nümunələr öz K-sına məhdudlaşır. Əgər
        bir K-də NƏ laysız, NƏ DƏ ona bağlı nümunə yoxdursa, defolt
        olaraq bu "digər layların nöqtələrini sükutla hovuzla" demək
        DEYİL (əvvəlki nöqsan, bax M1) — açıq xəta atılır, yalnız
        `allow_cross_layer_fallback=True` ilə bilərəkdən hovuzlanır.

        İnterpolyator 3D dəstəkləyirsə (`supports_z`, yalnız
        `OrdinaryKriging`) hər nümunənin Z-si (bax `_sample_depth`) və
        hədəfin öz K-sının həqiqi hüceyrə-mərkəzi dərinliyi ötürülür.
        Bunun sayəsində `allow_cross_layer_fallback` işə düşəndə fərqli
        laylardan gələn nöqtələr artıq BƏRABƏR yox, öz dərinlik
        fərqlərinə görə (range_v vasitəsilə) ÇƏKİLİ qatılır — yaxın lay
        uzaq laydan çox təsir edir (M2: geoloji cəhətdən əsaslandırılmış
        borclanma, kor-koranə hovuzlama yox).
        """
        use_z = self.interpolator.supports_z
        layer_mean_depth = target_depths = None
        if use_z:
            depths_grid = geometry.cell_depths().reshape(grid.shape)   # (nz, ny, nx)
            layer_mean_depth = depths_grid.mean(axis=(1, 2))
            target_depths = depths_grid

        layers = []
        for k in range(grid.nz):
            layer = k if dataset.is_layered() else None
            samples = dataset.samples_for(source, layer)
            if not samples and layer is not None:
                if not allow_cross_layer_fallback:
                    raise ValueError(
                        f"'{source}' üçün {k + 1}-ci laydan (K={k}) heç bir quyu "
                        "nöqtəsi yoxdur. Bu laya sükutla başqa layların "
                        "dəyərləri işlədilmir (bax ISH_HESABATI.md, M1). "
                        "Ya bu laya məlumat əlavə et, ya da bilərəkdən "
                        "allow_cross_layer_fallback=True keçir.")
                samples = dataset.samples_for(source, None)
                report.warn(
                    f"'{source}': {k + 1}-ci layda (K={k}) məlumat yoxdur — "
                    "allow_cross_layer_fallback=True ilə bütün laylardan "
                    "hovuzlanmış nöqtələr işlədildi (ekstrapolyasiya).")
            if not samples:
                raise ValueError(f"'{source}' üçün nöqtə tapılmadı.")

            values = np.asarray([s.values[source] for s in samples], float)
            if use_z:
                depths = np.asarray(
                    [self._sample_depth(s, k, layer_mean_depth) for s in samples], float)
                points = np.column_stack(
                    [[s.x for s in samples], [s.y for s in samples], depths])
                target_points = np.column_stack([targets, target_depths[k].ravel()])
            else:
                points = np.asarray([(s.x, s.y) for s in samples], float).reshape(-1, 2)
                target_points = targets

            layers.append(interpolate_property(
                self.interpolator, points, values, target_points,
                log_transform=rule.log_transform,
                minimum=rule.minimum, maximum=rule.maximum))
        return np.concatenate(layers)

    # ------------------------------------------------------- M4: cross-validation
    def cross_validate(self, dataset: WellDataset, source: str,
                       method: str = "loo", k: int = 5, seed: int = 42
                       ) -> Tuple[Dict[Optional[int], CrossValidationResult], Dict[Optional[int], str]]:
        """`source` üçün REAL dəqiqliyi ölçür — "100% dəqiq" vəd ETMİR.

        Laylı dataset-də hər lay ÖZ nöqtələri ilə (`samples_for` — məhz
        `_interpolate_volume`-un istifadə etdiyi eyni süzgəc) ayrıca
        doğrulanır, çünki M1/M2-dən sonra hər lay öz kriging səthinə
        malikdir — bir "qlobal model" yoxdur ki, tək ədədlə doğrulansın.
        Laysız dataset-də (məs. UI cədvəli) bir dəfə, bütün nöqtələrlə.

        Qaytarır: `(nəticələr, buraxılanlar)` — `buraxılanlar` 3 nöqtədən
        az olan laylar üçün səbəb mesajıdır (CV üçün bu laylar keçilir,
        səbəbsiz gizlədilmir).
        """
        rule = self.rules.get(source, PropertyRule(source))
        layers: Sequence[Optional[int]] = (
            [l for l in dataset.layers if l is not None] if dataset.is_layered() else [None])

        results: Dict[Optional[int], CrossValidationResult] = {}
        skipped: Dict[Optional[int], str] = {}
        for layer in layers:
            samples = dataset.samples_for(source, layer)
            if len(samples) < 3:
                label = f"K={layer}" if layer is not None else "bütün model"
                skipped[layer] = (
                    f"{label}: {len(samples)} nöqtə var, cross-validation üçün "
                    "ən azı 3 lazımdır — bu lay/model üçün doğrulama aparılmadı.")
                continue
            points = np.asarray([(s.x, s.y) for s in samples], float)
            values = np.asarray([s.values[source] for s in samples], float)
            runner = k_fold if method == "k-fold" else leave_one_out
            kwargs = {"k": k, "seed": seed} if method == "k-fold" else {}
            results[layer] = runner(
                self.interpolator, points, values, log_transform=rule.log_transform,
                compute_log_metrics=rule.log_transform, **kwargs)
        return results, skipped

    def cross_validate_all(self, dataset: WellDataset, method: str = "loo",
                           k: int = 5, seed: int = 42
                           ) -> Dict[str, Tuple[Dict[Optional[int], CrossValidationResult],
                                               Dict[Optional[int], str]]]:
        """PORO və mövcud PERMX/PERMY/PERMZ üçün `cross_validate` icra edir."""
        available = set(dataset.property_names())
        return {
            source: self.cross_validate(dataset, source, method=method, k=k, seed=seed)
            for source in _CROSS_VALIDATED_PROPERTIES if source in available
        }

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


_DISPLAY_LABEL = {"PORO": "POROSITY", "PERMX": "PERMEABILITY (PERMX)",
                  "PERMY": "PERMEABILITY (PERMY)", "PERMZ": "PERMEABILITY (PERMZ)"}


def format_cross_validation_report(
        all_results: Dict[str, Tuple[Dict[Optional[int], CrossValidationResult],
                                    Dict[Optional[int], str]]]) -> str:
    """`cross_validate_all()` nəticəsini istifadəçiyə göstərilən mətnə çevirir.

    Format tapşırıqdakı kimi (bax M4):

        POROSITY:
          RMSE = ...
          MAE  = ...
          R²   = ...

        PERMEABILITY (PERMX):
          ...

    Heç bir xassə üçün "100% dəqiq" iddiası YOXDUR — mövcud R²/RMSE
    olduğu kimi göstərilir, pis nəticə gizlədilmir.
    """
    if not all_results:
        return "Cross-validation üçün heç bir xassə tapılmadı (ən azı 3 nöqtə lazımdır)."
    blocks = []
    for source, (results, skipped) in all_results.items():
        label = _DISPLAY_LABEL.get(source, source)
        lines = [f"{label}:"]
        if not results:
            lines.append("  (heç bir lay üçün kifayət qədər nöqtə yoxdur)")
        for layer, result in sorted(results.items(), key=lambda kv: (kv[0] is None, kv[0])):
            layer_label = f"lay K={layer}" if layer is not None else ""
            lines.append(result.as_text(layer_label))
        for message in skipped.values():
            lines.append(f"  ⚠ {message}")
        blocks.append("\n".join(lines))
    return "\n\n".join(blocks)
