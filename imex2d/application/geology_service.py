"""İŞ AXINININ YENİ ADDIMI: quyu məlumatı → geoloji model.

    Karotaj interpretasiyası (xaricdə)
        ↓  CSV
    WellDataset
        ↓  xassə növü: KƏSİLMƏZ → interpolyasiya (IDW/Kriging)
        ↓             KATEQORİK → Sequential Indicator Simulation (Phase 4.1)
    GeologicalModel  →  ReservoirModel  →  Simulyasiya

Bu qat interpolyasiya alqoritmini TANIMIR — yalnız IPropertyInterpolator
interfeysini bilir. Alqoritm konstruktora inject edilir.

KRİTİK ELMİ QAYDA (Phase 4.1): heç bir KATEQORİK xassə (bax
`geology/property_types.py`) `interpolate_property()`-dən (kəsilməz
Kriging/IDW) KEÇMİR. Bu, `_interpolate_volume` daxilində AÇIQ yoxlanılır
(bax `build()`) — kateqorik sütun aşkarlansa, KƏSİLMƏZ yol heç
ÇAĞIRILMIR, bunun əvəzinə `_simulate_categorical_field()` (SIS,
`geology/facies.py`) işə düşür.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Optional, Sequence, Tuple, Union

import numpy as np

from ..domain.facies_field import FaciesField
from ..domain.geological_model import GeologicalModel
from ..domain.geometry import CellGeometry, depth_to_k, xy_to_ij
from ..domain.grid import CartesianGrid
from ..domain.properties import PropertyMap
from ..domain.structure import RegionSet
from ..domain.well_data import WellDataset
from ..geology.cross_validation import CrossValidationResult, k_fold, leave_one_out
from ..geology.distribution_analysis import log_transform_is_justified
from ..geology.facies import (FaciesVariogramParams, observed_proportions, simulate_sis)
from ..geology.hard_data import resolve_hard_data
from ..geology.interpolation import interpolate_property
from ..geology.property_types import PropertyType, classify_property
from ..geology.sgs import (DEFAULT_MIN_HARD_DATA_FOR_OWN_MODEL, FaciesPropertyConfig,
                           PropertyVariogramParams, simulate_sgs, simulate_sgs_facies_conditioned)
from ..interfaces.interpolation import IPropertyInterpolator

#: `cross_validate_all`-un baxdığı xassələr — PORO və PERM* istiqamətləri.
_CROSS_VALIDATED_PROPERTIES = ("PORO", "PERMX", "PERMY", "PERMZ")


@dataclass
class FaciesBuildConfig:
    """Bir kateqorik sütunun (məs. FACIES) necə simulyasiya olunacağı.

    `proportions` verilməyibsə (`None`) sərt datadan MÜŞAHİDƏ OLUNAN
    nisbətlər avtomatik hesablanır — bu, İSTİFADƏÇİNİN AÇIQ seçimi
    DEYİL, ona görə `report.warnings`-ə AÇIQ qeyd düşülür (bax
    `WellBasedGeologicalModelBuilder._simulate_categorical_field`).
    """
    proportions: Optional[Dict[int, float]] = None
    category_names: Optional[Dict[int, str]] = None
    variograms: Optional[Dict[int, FaciesVariogramParams]] = None
    seed: int = 0
    realization_id: int = 0
    search_radius: Optional[float] = None
    max_neighbors: Optional[int] = 24
    min_neighbors: int = 1
    on_conflict: str = "raise"


@dataclass
class ContinuousSGSConfig:
    """Bir KƏSİLMƏZ sütunun (PORO/PERMX/...) Sequential Gaussian
    Simulation (Phase 5) ilə (deterministik Kriging ƏVƏZİNƏ) necə
    simulyasiya olunacağı. Konfiqurasiya edilməyən sütunlar ƏVVƏLKİ
    (deterministik Kriging/IDW) yolu ilə davam edir — bu, TAM opt-in-dir.

    `facies_field_name` verilibsə (məs. `"FACIES"`), simulyasiya
    `model.facies_fields[facies_field_name]`-ə ŞƏRTLƏNİR (bax
    `geology/sgs.simulate_sgs_facies_conditioned`) — hər fasiya üçün
    (istəyə görə) AYRI `facies_configs` konfiqurasiyası.

    `log_space`/`bounds` verilməyəndə (`None`) MÜVAFİQ olaraq
    `distribution_analysis.log_transform_is_justified()` (data-əsaslı)
    və `DEFAULT_RULES[source]`-in hədləri İSTİFADƏ OLUNUR (mövcud
    kəsilməz yolla EYNİ hədlər, TƏKRARLANMIR).
    """
    variogram: Optional[PropertyVariogramParams] = None
    log_space: Optional[bool] = None
    bounds: Optional[Tuple[Optional[float], Optional[float]]] = None
    seed: int = 0
    realization_id: int = 0
    search_radius: Optional[float] = None
    max_neighbors: Optional[int] = 24
    min_neighbors: int = 1
    on_conflict: str = "raise"
    conflict_tolerance: float = 0.0
    facies_field_name: Optional[str] = None
    facies_configs: Optional[Dict[int, FaciesPropertyConfig]] = None
    min_hard_data_for_own_model: int = DEFAULT_MIN_HARD_DATA_FOR_OWN_MODEL


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
              allow_cross_layer_fallback: bool = False,
              facies_config: Optional[Dict[str, FaciesBuildConfig]] = None,
              property_type_overrides: Optional[Dict[str, PropertyType]] = None,
              sgs_config: Optional[Dict[str, ContinuousSGSConfig]] = None):
        """`allow_cross_layer_fallback` — bir K-təbəqəsində (dataset laylı
        olanda) heç bir quyu nöqtəsi yoxdursa nə edilsin.

        Defolt (`False`): AÇIQ XƏTA atılır — sükutla başqa laylardan
        (məs. lay 1-3) məlumat "sızdırılmır". `True` verilsə köhnə davranış
        (bütün laylardan hovuzlanmış nöqtələrlə kriging) işə düşür, amma
        `report.warnings`-ə açıq xəbərdarlıq yazılır ki, istifadəçi bunun
        bilərəkdən ekstrapolyasiya olduğunu bilsin (bax M0 sınaqları,
        `tests/test_layer_aware_kriging_leak.py`).

        `facies_config` — hər KATEQORİK sütun (bax `geology/property_
        types.py`) üçün `FaciesBuildConfig` (Phase 4.1). Kateqorik sütun
        HEÇ VAXT `_interpolate_volume`-dan (kəsilməz Kriging/IDW) KEÇMİR
        — bunun əvəzinə `_simulate_categorical_field` (SIS) çağırılır və
        nəticə `model.facies_fields`-ə (PropertyMap-DAN AYRICA) yazılır.

        `sgs_config` — hər KƏSİLMƏZ sütun üçün (istəyə görə)
        `ContinuousSGSConfig` (Phase 5). Konfiqurasiya edilməyən kəsilməz
        sütun ƏVVƏLKİ deterministik Kriging/IDW yolu ilə davam edir —
        SGS TAM opt-in-dir, mövcud davranışı DƏYİŞMİR. Kateqorik sütunlar
        HƏMİŞƏ kəsilməzlərdən ƏVVƏL emal olunur ki, `sgs_config`-in
        `facies_field_name` istinadı artıq qurulmuş olsun.
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
        categorical_sources = [s for s in available
                               if classify_property(s, property_type_overrides)
                               is PropertyType.CATEGORICAL]
        continuous_sources = [s for s in available if s not in categorical_sources]

        for source in categorical_sources:
            config = (facies_config or {}).get(source)
            facies_field = self._simulate_categorical_field(
                dataset, source, targets, grid, geometry, config, report)
            model.add_facies_field(facies_field)

        for source in continuous_sources:
            rule = self.rules.get(source, PropertyRule(source))
            sgs = (sgs_config or {}).get(source)
            if sgs is not None:
                values = self._simulate_continuous_sgs_field(
                    dataset, source, targets, grid, geometry, sgs, report, model)
            else:
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

    def _simulate_categorical_field(self, dataset: WellDataset, source: str,
                                    targets: np.ndarray, grid: CartesianGrid,
                                    geometry: CellGeometry,
                                    config: Optional[FaciesBuildConfig],
                                    report: "InterpolationReport") -> FaciesField:
        """KATEQORİK sütunu Sequential Indicator Simulation ilə (Phase
        4.1) simulyasiya edir — `_interpolate_volume`-un (kəsilməz
        Kriging/IDW) ƏVƏZİNƏ, ONUN YERİNƏ DEYİL (bu sütun `build()`-də
        heç vaxt `_interpolate_volume`-a ötürülmür).

        Bütün laylar BİR simulyasiyada (tam 3D X,Y,Z kondisioner +
        hədəf) işlənir — hər lay üçün AYRICA (əvvəlki `allow_cross_layer_
        fallback` kimi) DEYİL, çünki şaquli variogram (`range_v`) artıq
        laylar arası davamlılığı DOĞRU modelləşdirir (bax FACIES.md).
        """
        config = config or FaciesBuildConfig()
        raw_samples = [s for s in dataset.samples if source in s.values]
        resolved_samples = resolve_hard_data(raw_samples, source, grid, geometry,
                                             on_conflict=config.on_conflict)
        if not resolved_samples:
            raise ValueError(f"'{source}' üçün istifadə edilə bilən sərt data tapılmadı.")

        depths_grid = geometry.cell_depths().reshape(grid.shape)   # (nz, ny, nx)

        # Sərt datanı öz EV HÜCEYRƏSİNİN mərkəzinə "sancırıq" (snap).
        # SƏBƏB: real quyu demək olar HEÇ VAXT dəqiq hüceyrə mərkəzində
        # deyil — `simulate_sis`-in sərt-data hörməti dəqiq KOORDİNAT
        # üst-üstə düşməsinə əsaslanır (bax `facies._find_hard_data_
        # matches`), ona görə kondisioner nöqtəni HƏDƏF massivindəki
        # (`full_targets`) EYNİ hüceyrənin mərkəzi ilə eyniləşdiririk —
        # əks halda "sərt data honored" QORUNMUR (yalnız TƏSADÜFƏN quyu
        # mərkəzdə olanda işləyərdi). `resolve_hard_data` artıq eyni
        # hüceyrəyə düşən ziddiyyətli nümunələri BLOKLAYIB, ona görə bu
        # sancma YENİ ziddiyyət yaratmır.
        xs, ys, zs, codes = [], [], [], []
        skipped = 0
        for sample in resolved_samples:
            i, j = xy_to_ij(sample.x, sample.y, geometry)
            if sample.layer is not None:
                k = sample.layer
            elif sample.depth is not None:
                k = depth_to_k(sample.x, sample.y, sample.depth, geometry)
            else:
                k = None
            if k is None:
                skipped += 1
                continue
            xs.append((i + 0.5) * geometry.dx)
            ys.append((j + 0.5) * geometry.dy)
            zs.append(float(depths_grid[k, j, i]))
            codes.append(int(sample.values[source]))
        if skipped:
            report.warn(
                f"'{source}': {skipped} nümunə nə lay, nə dərinlik daşıyır — 3D mövqeyi "
                "müəyyən edilə bilmədi, simulyasiyaya daxil edilmədi.")
        if not codes:
            raise ValueError(f"'{source}' üçün 3D mövqeyi müəyyən edilə bilən sərt data yoxdur.")

        points = np.column_stack([xs, ys, zs])
        codes_array = np.asarray(codes, int)

        proportions = config.proportions
        if proportions is None:
            proportions = observed_proportions(codes_array)
            report.warn(
                f"'{source}': fasiya nisbətləri verilməyib — quyu datasından MÜŞAHİDƏ "
                "OLUNAN nisbətlər işlədildi (istifadəçi AÇIQ seçim ETMƏYİB): "
                + ", ".join(f"{k}={v:.3f}" for k, v in sorted(proportions.items())))
            proportion_source = "observed"
        else:
            proportion_source = "user"

        full_targets = np.concatenate(
            [np.column_stack([targets, depths_grid[k].ravel()]) for k in range(grid.nz)], axis=0)

        realization = simulate_sis(
            points, codes_array, full_targets, proportions, variograms=config.variograms,
            seed=config.seed, realization_id=config.realization_id,
            search_radius=config.search_radius, max_neighbors=config.max_neighbors,
            min_neighbors=config.min_neighbors)
        for message in realization.warnings:
            report.warn(f"'{source}' (SIS): {message}")

        variogram_metadata = {k: {"model": v.model, "nugget": v.nugget, "sill": v.sill,
                                  "range_": v.range_, "range_v": v.range_v,
                                  "azimuth_deg": v.azimuth_deg, "range_minor": v.range_minor}
                              for k, v in (config.variograms or {}).items()}

        return FaciesField(
            name=source,
            codes=realization.codes,
            category_names=dict(config.category_names or {}),
            realization_id=realization.realization_id,
            seed=realization.seed,
            requested_proportions=realization.requested_proportions,
            realized_proportions=realization.realized_proportions,
            variogram_metadata=variogram_metadata,
            conditioning_data_stats={
                "n_hard_points": int(codes_array.size),
                "n_wells": len({s.well for s in resolved_samples}),
                "n_skipped_unplaceable": skipped,
                "proportion_source": proportion_source,
                "negative_probability_events": realization.diagnostics.negative_probability_events,
                "excess_probability_events": realization.diagnostics.excess_probability_events,
                "nan_fallback_cells": realization.diagnostics.nan_fallback_cells,
                "zero_sum_fallback_cells": realization.diagnostics.zero_sum_fallback_cells,
                "n_cells_simulated": realization.diagnostics.n_cells_simulated,
            },
            warnings=list(realization.warnings))

    def _simulate_continuous_sgs_field(self, dataset: WellDataset, source: str,
                                       targets: np.ndarray, grid: CartesianGrid,
                                       geometry: CellGeometry, config: ContinuousSGSConfig,
                                       report: "InterpolationReport",
                                       model: GeologicalModel) -> np.ndarray:
        """KƏSİLMƏZ sütunu Sequential Gaussian Simulation ilə (Phase 5)
        simulyasiya edir — `_interpolate_volume`-un (deterministik
        Kriging) YERİNƏ, YALNIZ `sgs_config`-də AÇIQ tələb olunanda.

        `_simulate_categorical_field`-lə EYNİ hüceyrə-mərkəzinə-sancma
        (snap) qaydası (bax orada tapılan HƏQİQİ səhv) VƏ EYNİ 3D
        (X,Y,Z) tərtib məntiqi işlədilir — TƏKRARLANMIR (kod paylaşımı
        praktik olaraq mümkün olduğu qədər, iki metodun fərqli sərt-data
        həll strategiyası (`resolve_hard_data(tolerance=...)`) səbəbilə
        tam ortaq funksiyaya çıxarılmayıb).
        """
        raw_samples = [s for s in dataset.samples if source in s.values]
        resolved_samples = resolve_hard_data(raw_samples, source, grid, geometry,
                                             on_conflict=config.on_conflict,
                                             tolerance=config.conflict_tolerance)
        if not resolved_samples:
            raise ValueError(f"'{source}' üçün istifadə edilə bilən sərt data tapılmadı.")

        depths_grid = geometry.cell_depths().reshape(grid.shape)

        xs, ys, zs, values, cell_indices, used_samples = [], [], [], [], [], []
        skipped = 0
        for sample in resolved_samples:
            i, j = xy_to_ij(sample.x, sample.y, geometry)
            if sample.layer is not None:
                k = sample.layer
            elif sample.depth is not None:
                k = depth_to_k(sample.x, sample.y, sample.depth, geometry)
            else:
                k = None
            if k is None:
                skipped += 1
                continue
            used_samples.append(sample)
            xs.append((i + 0.5) * geometry.dx)
            ys.append((j + 0.5) * geometry.dy)
            zs.append(float(depths_grid[k, j, i]))
            values.append(float(sample.values[source]))
            cell_indices.append((i, j, k))
        if skipped:
            report.warn(f"'{source}': {skipped} nümunə nə lay, nə dərinlik daşıyır — "
                       "3D mövqeyi müəyyən edilə bilmədi, simulyasiyaya daxil edilmədi.")
        if not values:
            raise ValueError(f"'{source}' üçün 3D mövqeyi müəyyən edilə bilən sərt data yoxdur.")

        points = np.column_stack([xs, ys, zs])
        values_array = np.asarray(values, float)
        log_space = (config.log_space if config.log_space is not None
                    else log_transform_is_justified(values_array))
        rule = self.rules.get(source, PropertyRule(source))
        bounds = config.bounds if config.bounds is not None else (rule.minimum, rule.maximum)

        full_targets = np.concatenate(
            [np.column_stack([targets, depths_grid[k].ravel()]) for k in range(grid.nz)], axis=0)
        common_kwargs = dict(seed=config.seed, realization_id=config.realization_id,
                             search_radius=config.search_radius,
                             max_neighbors=config.max_neighbors,
                             min_neighbors=config.min_neighbors)

        if config.facies_field_name:
            facies_field = model.facies_fields.get(config.facies_field_name)
            if facies_field is None:
                raise ValueError(
                    f"'{source}': facies_field_name={config.facies_field_name!r} "
                    "model.facies_fields-də yoxdur (kateqorik sütun əvvəlcə emal olunmalıdır).")
            facies_at_points = np.array([
                int(sample.values[config.facies_field_name])
                if config.facies_field_name in sample.values
                else int(facies_field.codes[np.ravel_multi_index(cell, grid.shape)])
                for sample, cell in zip(used_samples, cell_indices)], int)
            realization = simulate_sgs_facies_conditioned(
                points, values_array, facies_at_points, full_targets, facies_field.codes,
                facies_configs=config.facies_configs, facies_reference=config.facies_field_name,
                min_hard_data_for_own_model=config.min_hard_data_for_own_model, **common_kwargs)
        else:
            realization = simulate_sgs(points, values_array, full_targets, log_space=log_space,
                                       bounds=bounds, variogram=config.variogram, **common_kwargs)

        for message in realization.warnings:
            report.warn(f"'{source}' (SGS): {message}")
        diag = realization.diagnostics
        if diag.n_cells_simulated:
            report.warn(
                f"'{source}' (SGS) diaqnostika: {diag.bound_corrections} hədd kəsilməsi "
                f"({diag.rate(diag.bound_corrections) * 100:.1f}%), "
                f"{diag.nan_fallback_cells} NaN-geri-dönüş "
                f"({diag.rate(diag.nan_fallback_cells) * 100:.1f}%).")
        return realization.values

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
