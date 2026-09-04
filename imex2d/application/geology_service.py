"""İŞ AXINININ YENİ ADDIMI: quyu məlumatı → geoloji model.

    Karotaj interpretasiyası (xaricdə)
        ↓  CSV
    WellDataset
        ↓  xassə növü: KƏSİLMƏZ → Phase B xassə-strategiyalı interpolyasiya
        ↓             KATEQORİK → Sequential Indicator Simulation (Phase 4.1,
        ↓                        defolt) / Phase B indikator kriginq (opt-in)
    GeologicalModel  →  ReservoirModel  →  Simulyasiya

Bu qat interpolyasiya alqoritmini TANIMIR — yalnız IPropertyInterpolator
interfeysini bilir. Alqoritm konstruktora inject edilir; bu, HƏLƏ DƏ
sərt-data honoring/qonşuluq/variogram ÖZƏYİ üçün doğrudur (`_kriging_
overrides` vasitəsilə inject edilmiş `OrdinaryKriging`-in AÇIQ dəyişdirdiyi
sahələr Phase B mühərrikinə ötürülür). Amma XASSƏNİN STATİSTİK TƏBİƏTİ
(çevirmə/hədlər/interpolyasiya NÖVÜ) artıq `geology/property_config.
PropertyStrategy` reyestrindən HƏLL OLUNUR (B-INTEGRATION-FIX) — bax
`_resolve_property_strategy`/`_interpolate_volume`.

KRİTİK ELMİ QAYDA (Phase 4.1, B-INTEGRATION-FIX-də DƏYİŞMƏDƏN qalıb): heç
bir KATEQORİK xassə (bax `geology/property_types.py`) kəsilməz Kriging/IDW-
dən KEÇMİR. Bu, `_interpolate_volume` daxilində AÇIQ yoxlanılır (bax
`build()`) — kateqorik sütun aşkarlansa, KƏSİLMƏZ yol heç ÇAĞIRILMIR,
bunun əvəzinə `_simulate_categorical_field()` (SIS) və ya (bax
`FaciesBuildConfig.deterministic`) Phase B `interpolate_categorical_field()`
işə düşür — hər ikisi indikator-əsaslıdır, HEÇ BİRİ fasiya kodunu ədədi
kəsilməz dəyər kimi kriging etmir (GATE B4).

KÖHNƏ Phase A `geology.interpolation.interpolate_property()` ARTIQ bu
modulda ÇAĞIRILMIR (bax B-INTEGRATION-FIX hesabatı) — funksiya özü
`geology/interpolation.py`-də SİLİNMƏYİB (başqa çağıranlar üçün geriyə-
uyğun qalır, bax `tests/test_geology_import.py`-in bir hissəsi onu
BİRBAŞA, bu moduldan keçmədən sınayır), sadəcə bu istehsalat boru
xəttinin İCRA YOLUNDAN ÇIXARILIB.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from ..domain.data_availability import (DataStatus, ModelDataAvailability,
                                        PropertyAvailability, format_layers)
from ..domain.diagnostics import DiagnosticReport
from ..domain.facies_field import FaciesField
from ..domain.geological_model import GeologicalModel
from ..domain.geometry import CellGeometry, depth_to_k, xy_to_ij
from ..domain.grid import CartesianGrid
from ..domain.properties import (CategoricalUncertainty, PropertyMap, PropertyProvenance,
                                 PropertyUncertainty)
from ..domain.structure import RegionSet
from ..domain.well_data import WellDataset
from ..geology.layer_availability import (LayerDataPolicy, compute_availability,
                                          hard_data_cells, unassigned_samples)
from ..geology.cross_validation import CrossValidationResult, k_fold, leave_one_out
from ..geology.distribution_analysis import log_transform_is_justified
from ..geology.facies import (FaciesVariogramParams, observed_proportions, simulate_sis)
from ..geology.hard_data import resolve_hard_data
from ..geology.interpolation import OrdinaryKriging
from ..geology.property_config import (PropertyStrategy, VariableType, resolve_strategy)
from ..geology.property_interpolation import (CategoricalEstimate, PropertyEstimate,
                                              interpolate_categorical_field,
                                              interpolate_property_field)
from ..geology.property_types import PropertyType, classify_property
from ..geology.sgs import (DEFAULT_MIN_HARD_DATA_FOR_OWN_MODEL, FaciesPropertyConfig,
                           PropertyVariogramParams, simulate_sgs, simulate_sgs_facies_conditioned)
from ..geology.spatial_search import (SUPPORT_BOUNDARY, SUPPORT_EXTRAPOLATED, SUPPORT_WEAK,
                                      SUPPORT_WELL)
from ..geology.transforms import IDENTITY_TRANSFORM, BackTransform, LogTransform
from ..interfaces.interpolation import IPropertyInterpolator

#: `cross_validate_all`-un baxdığı xassələr — PORO və PERM* istiqamətləri.
_CROSS_VALIDATED_PROPERTIES = ("PORO", "PERMX", "PERMY", "PERMZ")

#: `Confidence` (ORDİNAL interpretasiya kateqoriyası, bax `geology/
#: property_interpolation.Confidence`) → `[0,1]` ədədi bal. Bu, EHTİMAL
#: DEYİL: mövcud, ƏSASLANDIRILMIŞ kateqoriyanın (qonşu sayı + məsafə +
#: nisbi kriginq variansı) MONOTON ədədi əksidir, ona görə "saxta rəqəm"
#: deyil — amma `PropertyProvenance.confidence_kind` ilə AÇIQ şəkildə
#: "ordinal_support_score" kimi etiketlənir.
_CONFIDENCE_SCORE = {"high": 0.90, "medium": 0.60, "low": 0.30,
                     "extrapolated": 0.10}

#: Tamamlama (completion) HEÇ VAXT ölçmə/interpolyasiya səviyyəsində
#: etibarlı sayılmır — zərfin İÇİNDƏ də tavanı var, KƏNARINDA daha aşağı.
_TREND_INSIDE_CEILING = 0.50
_TREND_OUTSIDE_CEILING = 0.30


def _confidence_scores(confidence) -> np.ndarray:
    """`Confidence` massivini `[0,1]` bal massivinə çevirir; tanınmayan
    dəyər `NaN` (uydurulmur)."""
    return np.asarray([_CONFIDENCE_SCORE.get(str(value), np.nan)
                       for value in np.asarray(confidence, dtype=object)], float)


def _mean_or_none(scores) -> Optional[float]:
    """Lay üzrə orta bal — hamısı NaN olanda `None` ("hesablanmadı")."""
    values = np.asarray(scores, float)
    finite = np.isfinite(values)
    return float(values[finite].mean()) if finite.any() else None


@dataclass
class FaciesBuildConfig:
    """Bir kateqorik sütunun (məs. FACIES) necə simulyasiya olunacağı.

    `proportions` verilməyibsə (`None`) sərt datadan MÜŞAHİDƏ OLUNAN
    nisbətlər avtomatik hesablanır — bu, İSTİFADƏÇİNİN AÇIQ seçimi
    DEYİL, ona görə `report.warnings`-ə AÇIQ qeyd düşülür (bax
    `WellBasedGeologicalModelBuilder._simulate_categorical_field`).

    `deterministic` (B-INTEGRATION-FIX, TAM opt-in, DEFOLT `False`):
    `True` olanda bu sütun SIS (`simulate_sis`) ƏVƏZİNƏ Phase B-nin
    deterministik indikator-kriginq mühərriki
    (`property_interpolation.interpolate_categorical_field`) ilə
    hesablanır — stoxastik realizasiya YOX, ən-ehtimallı kateqoriya +
    tam ehtimal vektoru (bax `WellBasedGeologicalModelBuilder.
    _estimate_categorical_field_phase_b`). `seed`/`realization_id`/
    `variograms`/`proportions` bu rejimdə İSTİFADƏ OLUNMUR (SIS-ə
    məxsusdur); əvəzinə `WellBasedGeologicalModelBuilder`-in konstruktor
    zamanı aldığı Kriging parametrləri (bax `_kriging_overrides`)
    işlədilir. DEFOLT (`False`) ilə davranış TAM ƏVVƏLKİ kimi qalır
    (SIS) — mövcud testlər DƏYİŞMİR."""
    proportions: Optional[Dict[int, float]] = None
    category_names: Optional[Dict[int, str]] = None
    variograms: Optional[Dict[int, FaciesVariogramParams]] = None
    seed: int = 0
    realization_id: int = 0
    search_radius: Optional[float] = None
    max_neighbors: Optional[int] = 24
    min_neighbors: int = 1
    on_conflict: str = "raise"
    deterministic: bool = False


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


class CompletionMethod(str, Enum):
    """MƏLUMATSIZ layın (data yoxdur) necə tamamlanacağı — tapşırıq §9.

    HEÇ BİRİ nəticəni `MEASURED` kimi qeyd ETMİR; hər biri provenance-də
    öz adı ilə görünür.

    `NONE` — DEFOLT. Lay TAMAMLANMIR: orijinal sahə varsa `PRESERVED`,
        yoxdursa `MISSING` qalır və simulyator ÖNCƏSİ validasiya
        BLOKLAYIR (§8/§10 — "sistem səssiz fərziyyə yaratmamalıdır").
    `PRESERVE_ORIGINAL` — mövcud geoloji prior (`original_fields`) olduğu
        kimi saxlanılır → `PRESERVED`.
    `VERTICAL_TREND` — məlumatlı layların lay-ortaları üzrə dərinliyə
        görə XƏTTİ trend qurulur, məlumatsız laya YALNIZ LAY ORTASI
        verilir. LATERAL struktur UYDURULMUR (qonşu layın xəritəsi
        KOPYALANMIR — §26). → `ESTIMATED` / zərfdən kənarda `EXTRAPOLATED`.
    `GEOSTATISTICAL_3D` — mövcud Kriging mühərriki BÜTÜN layların sərt
        datası ilə, HƏQİQİ 3D (X,Y,Z) məsafə və şaquli range ilə həmin
        laya qiymət verir (yaxın lay uzaq laydan çox təsir edir).
        → `ESTIMATED` / `EXTRAPOLATED`.
    `SGS` — Sequential Gaussian Simulation realizasiyası → `SIMULATED`
        (§17: HEÇ VAXT `MEASURED` deyil).
    `CONSTANT` — istifadəçinin AÇIQ verdiyi lay dəyəri → `ESTIMATED`.
    """

    NONE = "none"
    PRESERVE_ORIGINAL = "preserve_original"
    VERTICAL_TREND = "vertical_trend"
    GEOSTATISTICAL_3D = "geostatistical_3d"
    SGS = "sgs"
    CONSTANT = "constant"


@dataclass(frozen=True)
class CompletionSpec:
    """Bir xassə üçün tamamlama qaydası."""

    method: CompletionMethod = CompletionMethod.NONE
    #: `CONSTANT` üçün MƏCBURİ dəyər (fiziki vahiddə).
    value: Optional[float] = None
    #: İstifadəçinin AÇIQ bəyan etdiyi etibarlılıq `[0,1]`. Verilməyəndə
    #: (`None`) HESABLANA BİLƏNDƏ hesablanır, bilinməyəndə `NaN` qalır —
    #: SAXTA rəqəm YARADILMIR (§18).
    confidence: Optional[float] = None
    #: `SGS` üçün konfiqurasiya (verilməyibsə defolt `ContinuousSGSConfig`).
    sgs: Optional["ContinuousSGSConfig"] = None
    #: YALNIZ bu laylar tamamlansın (0-əsaslı). `None` — bütün məlumatsız
    #: laylar.
    layers: Optional[Sequence[int]] = None


@dataclass
class LayerInterpolationConfig:
    """LAY-MƏLUMATLI rejimin BÜTÜN parametrləri — `build(layer_config=...)`.

    Bu obyekt VERİLMƏYƏNDƏ (`None`, defolt) `build()` TAM ƏVVƏLKİ kimi
    işləyir (§25 geriyə-uyğunluq). Verildikdə isə §2-nin altı anlayışı
    ayrı-ayrı idarə olunur:

        `policy`                  → məlumat mövcudluğu necə oxunur
        `target_layers`           → İNTERPOLYASİYA HƏDƏFİ (0-əsaslı K!)
        `property_target_layers`  → xassə-üzrə hədəf (üstünlüyü var)
        `completion`              → məlumatsız layların tamamlanması
        `property_completion`     → xassə-üzrə tamamlama (üstünlüyü var)
        `original_fields`         → ORİJİNAL sahə (varsa) — `final`-ın bazası

    `target_layers=None` → "məlumatı OLAN bütün laylar" (səssiz genişlənmə
    YOXDUR: məlumatsız laya heç vaxt interpolyasiya edilmir).

    QEYD (UI 1-əsaslı, mühərrik 0-əsaslı): burada indekslər HƏMİŞƏ
    0-əsaslıdır. Çevirmə UI-də (`panels.py`) aparılır — `domain/
    data_availability.parse_layers()` bunu bir yerdə edir, hər çağıranda
    təkrarlanmır.
    """

    policy: LayerDataPolicy = LayerDataPolicy.STRICT
    target_layers: Optional[Sequence[int]] = None
    property_target_layers: Dict[str, Sequence[int]] = field(default_factory=dict)
    completion: CompletionSpec = field(default_factory=CompletionSpec)
    property_completion: Dict[str, CompletionSpec] = field(default_factory=dict)
    original_fields: Dict[str, np.ndarray] = field(default_factory=dict)
    #: Tamamlanmayan laylarda orijinal sahə (varsa) SAXLANILSIN mı.
    #: `False` — orijinal olsa belə `MISSING` qalır (ən sərt rejim).
    preserve_original_when_missing: bool = True

    def targets_for(self, source: str, nz: int) -> Optional[List[int]]:
        """`source` üçün istənilən hədəf laylar, yoxlanmış (0-əsaslı).

        `None` → "hədəf açıq verilməyib" (çağıran məlumatlı layları
        işlədəcək). Diapazondan kənar indeks AÇIQ `ValueError`-dur
        (§23.2 — səssiz kəsmə YOXDUR).
        """
        requested = self.property_target_layers.get(source, self.target_layers)
        if requested is None:
            return None
        layers = sorted({int(k) for k in requested})
        outside = [k for k in layers if not 0 <= k < nz]
        if outside:
            raise ValueError(
                f"'{source}': interpolyasiya üçün seçilmiş lay indeksi grid-dən "
                f"kənardadır: K={outside} (icazə verilən: 0..{nz - 1}, NZ={nz}).")
        if not layers:
            raise ValueError(
                f"'{source}': interpolyasiya üçün seçilmiş lay siyahısı BOŞDUR — "
                "ən azı bir lay seçin, ya da seçimi tamamilə boş buraxın "
                "(o halda məlumatı olan bütün laylar işlədilir).")
        return layers

    def completion_for(self, source: str) -> CompletionSpec:
        return self.property_completion.get(source, self.completion)


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
    #: BLOKLAYAN problemlər (lay-məlumatlı rejim) — model QAYTARILIR ki,
    #: istifadəçi onu görüb completion strategiyası seçə bilsin, amma
    #: `ReservoirModelBuilder`/simulyator onu QƏBUL ETMİR (bax
    #: `GeologicalModel.completeness_issues()`).
    blocking: list = field(default_factory=list)
    #: Lay-üzrə mövcudluq/status mənzərəsi (varsa) — hesabat mətnində
    #: göstərilir və UI-də cədvələ çevrilir.
    availability: Optional[ModelDataAvailability] = None

    def add(self, target: str, source: str, log_transform: bool, values: np.ndarray):
        """Statistika NaN-ə DÖZÜMLÜDÜR: lay-məlumatlı rejimdə MISSING
        hüceyrələr `NaN` daşıyır — onları statistikaya qatmaq bütün
        sətri "nan" edərdi. MISSING sayı AYRICA göstərilir, gizlədilmir."""
        values = np.asarray(values, float)
        finite = np.isfinite(values)
        entry = {"target": target, "source": source, "log": log_transform,
                 "missing": int(np.sum(~finite)), "n": int(values.size)}
        if finite.any():
            entry.update({"min": float(values[finite].min()),
                          "max": float(values[finite].max()),
                          "mean": float(values[finite].mean())})
        else:
            entry.update({"min": float("nan"), "max": float("nan"),
                          "mean": float("nan")})
        self.entries.append(entry)

    def warn(self, message: str) -> None:
        self.warnings.append(message)

    def block(self, message: str) -> None:
        self.blocking.append(message)

    @property
    def has_blocking(self) -> bool:
        return bool(self.blocking)

    def as_text(self) -> str:
        lines = [f"Üsul: {self.method}"]
        for entry in self.entries:
            missing = (f"  MISSING {entry['missing']}/{entry['n']}"
                       if entry.get("missing") else "")
            lines.append(
                f"  {entry['target']:<6} ← {entry['source']:<6} "
                f"{'(log)' if entry['log'] else '     '}  "
                f"min {entry['min']:.4g}  orta {entry['mean']:.4g}  "
                f"maks {entry['max']:.4g}{missing}")
        if self.availability is not None:
            lines.append("")
            lines.append("Lay-üzrə vəziyyət:")
            lines.append(self.availability.as_text())
        for message in self.warnings:
            lines.append(f"  ⚠ {message}")
        for message in self.blocking:
            lines.append(f"  ⛔ {message}")
        return "\n".join(lines)


class WellBasedGeologicalModelBuilder:
    """Quyu nöqtələrindən grid xassələri qurur."""

    def __init__(self, interpolator: IPropertyInterpolator,
                 rules: Optional[Dict[str, PropertyRule]] = None):
        self.interpolator = interpolator
        self.rules = dict(DEFAULT_RULES)
        #: İSTİFADƏÇİNİN konstruktora AÇIQ verdiyi (modul `DEFAULT_RULES`-
        #: dən FƏRQLİ) qaydalar — `self.rules`-dan AYRICA saxlanılır ki,
        #: `_resolve_property_strategy` "istifadəçi bilərəkdən bunu
        #: dəyişdi" ilə "bu, sadəcə modul defoltudur" fərqini bilsin (bax
        #: orada, B-INTEGRATION-FIX §5: köhnə `PropertyRule` geriyə-uyğun
        #: qalır, amma Phase B `PropertyStrategy` reyestri ARTIQ hədlərin/
        #: çevirmənin ƏSAS mənbəyidir).
        self._explicit_rules: Dict[str, PropertyRule] = dict(rules) if rules else {}
        if rules:
            self.rules.update(rules)

    # ---------------------------------------------------------- public
    def build(self, dataset: WellDataset, spec: GeologicalGridSpec,
              ky_over_kx: float = 1.0, kv_over_kh: float = 0.1,
              name: str = "Quyu məlumatından geoloji model",
              allow_cross_layer_fallback: bool = False,
              facies_config: Optional[Dict[str, FaciesBuildConfig]] = None,
              property_type_overrides: Optional[Dict[str, PropertyType]] = None,
              sgs_config: Optional[Dict[str, ContinuousSGSConfig]] = None,
              calibrated_strategies: Optional[Dict[str, PropertyStrategy]] = None,
              layer_config: Optional[LayerInterpolationConfig] = None):
        """`allow_cross_layer_fallback` — bir K-təbəqəsində (dataset laylı
        olanda) heç bir quyu nöqtəsi yoxdursa nə edilsin.

        `calibrated_strategies` (PHASE C, TAM opt-in) — `calibrate_
        property()`-nin (spatial-block CV-əsaslı model seçimi) qaytardığı
        `ModelSelectionReport.selected.candidate.strategy`-ni birbaşa
        production interpolyasiyasında İSTİFADƏ ETMƏK üçün: `{source:
        PropertyStrategy}`. Verilməyən xassələr ƏVVƏLKİ kimi `_resolve_
        property_strategy()`-dən (reyestr + `rules`) keçir — CV çəkisi
        HEÇ VAXT production yolunu MƏCBURİ dəyişmir, yalnız çağıran AÇIQ
        şəkildə "bu xassə üçün kalibrlənmiş modeli işlət" desə.

        Defolt (`False`): AÇIQ XƏTA atılır — sükutla başqa laylardan
        (məs. lay 1-3) məlumat "sızdırılmır". `True` verilsə köhnə davranış
        (bütün laylardan hovuzlanmış nöqtələrlə kriging) işə düşür, amma
        `report.warnings`-ə açıq xəbərdarlıq yazılır ki, istifadəçi bunun
        bilərəkdən ekstrapolyasiya olduğunu bilsin (bax M0 sınaqları,
        `tests/test_layer_aware_kriging_leak.py`).

        `facies_config` — hər KATEQORİK sütun (bax `geology/property_
        types.py`) üçün `FaciesBuildConfig` (Phase 4.1). Kateqorik sütun
        HEÇ VAXT `_interpolate_volume`-dan (kəsilməz Kriging/IDW) KEÇMİR
        — bunun əvəzinə DEFOLT olaraq `_simulate_categorical_field` (SIS)
        çağırılır, `FaciesBuildConfig.deterministic=True` verilibsə isə
        Phase B-nin indikator-kriginq mühərriki (`_estimate_categorical_
        field_phase_b`, B-INTEGRATION-FIX) — hər ikisi ehtimal-əsaslı
        indikator yoludur, HEÇ BİRİ fasiya kodunu ədədi kəsilməz dəyər
        kimi kriging etmir. Nəticə `model.facies_fields`-ə (PropertyMap-
        DAN AYRICA) yazılır.

        `layer_config` (LAY-MƏLUMATLI REJİM, TAM opt-in, bax
        `LAYER_AWARE_MODELING.md`) — verilməyəndə (`None`, DEFOLT) bu
        metod TAM ƏVVƏLKİ kimi işləyir: nə `model.provenance`, nə
        `model.availability` doldurulur, `report.blocking` boş qalır.
        Verildikdə isə:

          · xassə mövcudluğu HƏR LAY ÜÇÜN AYRICA hesablanır
            (`compute_availability`) — PORO L4-də ola bilər, PERMX olmaya;
          · YALNIZ seçilmiş VƏ məlumatlı laylar interpolyasiya olunur
            (`_build_layer_aware_field`) — `allow_cross_layer_fallback`
            bu rejimdə İSTİFADƏ OLUNMUR (onun yerini AÇIQ `CompletionSpec`
            tutur, iki fərqli məntiq eyni layı doldurmasın deyə);
          · məlumatsız lay YALNIZ AÇIQ `CompletionSpec` ilə doldurulur,
            əks halda `MISSING` qalır və `report.blocking`-ə düşür;
          · nəticə `PropertyProvenance` (status/üsul/etibarlılıq) ilə
            birlikdə modeldə saxlanılır.

        DİQQƏT (davranış fərqi): bu rejimdə MISSING lay `build()`-dan
        XƏTA ATMIR — model QAYTARILIR ki, istifadəçi onu 3D-də görüb
        completion seçə bilsin. Simulyasiya qapısı `ReservoirModelBuilder`
        -dədir (`GeologicalModel.completeness_issues()`), yəni natamam
        model HEÇ VAXT simulyatora keçmir.

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
        report = InterpolationReport(
            method=f"Phase B (xassə-strategiyası reyestri) — konfiqurasiya edilmiş "
                   f"interpolyator: {self.interpolator.describe()}")
        if not isinstance(self.interpolator, OrdinaryKriging):
            report.warn(
                f"Seçilmiş üsul ('{self.interpolator.describe()}') qeydiyyatlı (PORO, "
                "PERMX/Y/Z, SW, NTG, VSH, FACIES, ...) kəsilməz xassələr üçün İSTİFADƏ "
                "OLUNMUR — bu xassələrin interpolyasiya NÖVÜ (Kriging/log-kriging/"
                "logit-kriging) `PropertyStrategy` reyestrindən gəlir (B-INTEGRATION-FIX, "
                "GATE: ad → strategiya → çevirmə → üsul, `if property==...` YOXDUR). "
                "Yalnız Kriging seçilibsə onun parametrləri (range/sill/nugget/axtarış "
                "radiusu) Phase B mühərrikinə ötürülür.")

        targets = self._cell_centres(grid, spec)
        available = dataset.property_names()
        categorical_sources = [s for s in available
                               if classify_property(s, property_type_overrides)
                               is PropertyType.CATEGORICAL]
        continuous_sources = [s for s in available if s not in categorical_sources]

        for source in categorical_sources:
            config = (facies_config or {}).get(source)
            if config is not None and config.deterministic:
                facies_field = self._estimate_categorical_field_phase_b(
                    dataset, source, targets, grid, geometry, config, report, model)
            else:
                facies_field = self._simulate_categorical_field(
                    dataset, source, targets, grid, geometry, config, report)
            model.add_facies_field(facies_field)

        availability = None
        if layer_config is not None:
            availability = self._prepare_availability(
                dataset, geometry, continuous_sources, layer_config, report)
            model.availability = availability
            report.availability = availability
            self._record_categorical_availability(
                dataset, geometry, categorical_sources, layer_config, availability,
                report)

        for source in continuous_sources:
            rule = self.rules.get(source, PropertyRule(source))
            sgs = (sgs_config or {}).get(source)
            if layer_config is not None:
                values = self._build_layer_aware_field(
                    dataset, source, rule, targets, grid, geometry, model, report,
                    availability, layer_config, sgs,
                    (calibrated_strategies or {}).get(source))
            elif sgs is not None:
                values = self._simulate_continuous_sgs_field(
                    dataset, source, targets, grid, geometry, sgs, report, model)
            else:
                values = self._interpolate_volume(dataset, source, rule, targets, grid,
                                                  allow_cross_layer_fallback, report, geometry,
                                                  model,
                                                  (calibrated_strategies or {}).get(source))
            model.add_property(PropertyMap.from_array(rule.target, values,
                                                      grid.ncell))
            report.add(rule.target, source, rule.log_transform, values)

        self._fill_missing_permeability(model, grid, ky_over_kx, kv_over_kh, report)
        issues = model.validate()
        if layer_config is None:
            if issues:
                raise ValueError("Qurulan geoloji model natamamdır: " + "; ".join(issues))
            return model, report

        # LAY-MƏLUMATLI rejim: MISSING lay AÇIQ, GÖZLƏNİLƏN nəticədir
        # (istifadəçi hələ completion strategiyası seçməyib) — model
        # QAYTARILIR ki, 3D-də/hesabatda GÖRÜNSÜN, amma `report.blocking`
        # doludur və `ReservoirModelBuilder` onu QƏBUL ETMİR (§20).
        # QALAN hər cür problem (fiziki cəhətdən qeyri-mümkün dəyər,
        # dejenerativ həndəsə) ƏVVƏLKİ kimi DƏRHAL xətadır.
        completeness = set(model.completeness_issues())
        other = [issue for issue in issues if issue not in completeness]
        if other:
            raise ValueError("Qurulan geoloji model natamamdır: " + "; ".join(other))
        for message in sorted(completeness):
            report.block(message)
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
    def _gather_categorical_hard_data(dataset: WellDataset, source: str, grid: CartesianGrid,
                                      geometry: CellGeometry, on_conflict: str,
                                      report: "InterpolationReport"):
        """Kateqorik sütun üçün sərt datanı 3D (X,Y,Z) nöqtələrə "sancır"
        (bax aşağıdakı qeyd) — həm SIS (`_simulate_categorical_field`),
        həm Phase B (`_estimate_categorical_field_phase_b`) EYNİ
        gathering məntiqini işlədir (TƏKRARLANMIR).

        Sərt datanı öz EV HÜCEYRƏSİNİN mərkəzinə "sancırıq" (snap).
        SƏBƏB: real quyu demək olar HEÇ VAXT dəqiq hüceyrə mərkəzində
        deyil — indikator-əsaslı mühərriklərin (SIS, Phase B kriginq)
        sərt-data hörməti dəqiq KOORDİNAT üst-üstə düşməsinə əsaslanır,
        ona görə kondisioner nöqtəni HƏDƏF massivindəki EYNİ hüceyrənin
        mərkəzi ilə eyniləşdiririk — əks halda "sərt data honored"
        QORUNMUR (yalnız TƏSADÜFƏN quyu mərkəzdə olanda işləyərdi).
        `resolve_hard_data` artıq eyni hüceyrəyə düşən ziddiyyətli
        nümunələri BLOKLAYIB, ona görə bu sancma YENİ ziddiyyət yaratmır.
        """
        raw_samples = [s for s in dataset.samples if source in s.values]
        resolved_samples = resolve_hard_data(raw_samples, source, grid, geometry,
                                             on_conflict=on_conflict)
        if not resolved_samples:
            raise ValueError(f"'{source}' üçün istifadə edilə bilən sərt data tapılmadı.")

        depths_grid = geometry.cell_depths().reshape(grid.shape)   # (nz, ny, nx)
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
        return resolved_samples, points, codes_array, skipped

    def _estimate_categorical_field_phase_b(self, dataset: WellDataset, source: str,
                                            targets: np.ndarray, grid: CartesianGrid,
                                            geometry: CellGeometry,
                                            config: FaciesBuildConfig,
                                            report: "InterpolationReport",
                                            model: GeologicalModel) -> FaciesField:
        """KATEQORİK sütunu Phase B-nin DETERMİNİSTİK indikator-kriginq
        mühərriki (`property_interpolation.interpolate_categorical_field`)
        ilə hesablayır — `FaciesBuildConfig.deterministic=True` olanda,
        `_simulate_categorical_field` (SIS) ƏVƏZİNƏ (bax `build()`).

        SIS kimi bu da ƏSLA kəsilməz kriginqə keçmir (GATE B4) — hər
        kateqoriya üçün İNDİKATOR kriginq → [0,1]-ə kəsilib normallaşan
        ehtimal vektoru → ən-ehtimallı kod (bax `interpolate_categorical_
        field` docstring-i). Nəticə STOXASTİK REALİZASİYA DEYİL (seed/
        realization_id mənasız — `0` qoyulur), ona görə `FaciesField`-ə
        çevriləndə `requested_proportions` BOŞ qalır (heç bir hədəf
        nisbət YOXDUR, yalnız DATA-dan gələn ən-ehtimallı təsnifat).
        """
        resolved_samples, points, codes_array, _skipped = self._gather_categorical_hard_data(
            dataset, source, grid, geometry, config.on_conflict, report)
        depths_grid = geometry.cell_depths().reshape(grid.shape)
        full_targets = np.concatenate(
            [np.column_stack([targets, depths_grid[k].ravel()]) for k in range(grid.nz)], axis=0)

        estimate = interpolate_categorical_field(
            points, codes_array, full_targets, property_name=source,
            kriging_overrides=self._kriging_overrides())
        for message in estimate.warnings:
            report.warn(f"'{source}' (Phase B kateqorik): {message}")

        realized_proportions = {
            int(code): float(np.mean(estimate.most_probable == code))
            for code in estimate.categories}

        model.add_uncertainty(source, CategoricalUncertainty(
            name=source, categories=estimate.categories,
            probabilities=estimate.probabilities, entropy=estimate.entropy,
            normalized_entropy=estimate.normalized_entropy,
            max_probability=estimate.max_probability, confidence=estimate.confidence,
            support=np.asarray(estimate.support, dtype=object),
            neighbor_count=estimate.neighbor_count, nearest_distance=estimate.nearest_distance,
            extrapolated=estimate.extrapolated,
            n_probability_corrections=estimate.n_probability_corrections,
            warnings=list(estimate.warnings)))

        return FaciesField(
            name=source, codes=estimate.most_probable,
            category_names=dict(config.category_names or {}),
            realization_id=0, seed=0,
            requested_proportions={}, realized_proportions=realized_proportions,
            conditioning_data_stats={
                "n_hard_points": int(codes_array.size),
                "n_wells": len({s.well for s in resolved_samples}),
                "method": "phase_b_indicator_kriging",
                "mean_normalized_entropy": float(np.mean(estimate.normalized_entropy)),
                "n_probability_corrections": estimate.n_probability_corrections,
            },
            warnings=list(estimate.warnings))

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
        resolved_samples, points, codes_array, skipped = self._gather_categorical_hard_data(
            dataset, source, grid, geometry, config.on_conflict, report)
        depths_grid = geometry.cell_depths().reshape(grid.shape)   # (nz, ny, nx)

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
                                       model: GeologicalModel,
                                       layers: Optional[Sequence[int]] = None) -> np.ndarray:
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

        # `layers` verilibsə YALNIZ həmin K-ların hüceyrələri simulyasiya
        # olunur (§24 — lazımsız hesablama yoxdur). ŞƏRTLƏNDİRMƏ (sərt
        # data) BÜTÜN laylardan gəlməkdə DAVAM EDİR: SGS-in şaquli
        # kəsilməzliyi məhz buna əsaslanır (§17), yalnız HƏDƏF dəyişir.
        selected = list(range(grid.nz)) if layers is None else [int(k) for k in layers]
        full_targets = np.concatenate(
            [np.column_stack([targets, depths_grid[k].ravel()]) for k in selected], axis=0)
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
            # hədəf məhdudlaşdırılıbsa fasiya kodları DA eyni hüceyrələrə
            # kəsilməlidir — əks halda kod/hədəf uyğunluğu sürüşərdi.
            target_codes = np.concatenate(
                [facies_field.codes.reshape(grid.shape)[k].ravel() for k in selected])
            realization = simulate_sgs_facies_conditioned(
                points, values_array, facies_at_points, full_targets, target_codes,
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

    # ---------------------------------------------- B-INTEGRATION-FIX: Phase B körpüsü
    def _resolve_property_strategy(self, source: str) -> PropertyStrategy:
        """`source` üçün `PropertyStrategy` — reyestr (`property_config.
        resolve_strategy`) ƏSAS mənbədir (B-INTEGRATION-FIX §4: ad →
        strategiya → çevirmə → üsul → hədlər → qeyri-müəyyənlik, HEÇ bir
        `if property_name == ...` YOXDUR).

        `self._explicit_rules[source]` (konstruktora İSTİFADƏÇİNİN AÇIQ
        verdiyi `PropertyRule`, modul `DEFAULT_RULES`-dən FƏRQLİ olaraq)
        varsa, YALNIZ o, strategiyanın hədlərini/loq-çevirməsini
        `derive()` ilə ƏVƏZ edir — geriyə-uyğunluq (§5): köhnə
        `rules={...}` konstruktor parametri İNDİ DƏ işləyir (bax
        `tests/test_geology_import.py::test_custom_rule_overrides_
        default`), amma İNDİ Phase B strategiyasının ÜSTÜNDƏN keçir,
        onu ƏVƏZ ETMİR (QC/uncertainty/anizotropluq strategiyadan
        qalır)."""
        strategy = resolve_strategy(source)
        rule = self._explicit_rules.get(source)
        if rule is None:
            return strategy
        overrides: Dict[str, object] = {}
        if rule.minimum is not None or rule.maximum is not None:
            # DİQQƏT: köhnə `PropertyRule.minimum/maximum` Phase A-da YALNIZ
            # NƏTİCƏNİ kəsirdi (`interpolate_property(minimum=..., maximum=...)`
            # — girişə TOXUNMURDU). Ona görə YALNIZ `output_bounds` əvəz
            # olunur; `physical_bounds` (Phase B-nin QC girişi RƏDD ETMƏ
            # həddi) strategiyanın ÖZ dəyərində qalır — əks halda dar bir
            # `maximum` HƏQİQİ, fiziki etibarlı sərt datanı QC-də səssizcə
            # ATAR (tutulub, bax `tests/test_geology_import.py::
            # test_custom_rule_overrides_default` — geniş fiziki hədd, dar
            # NƏTİCƏ həddi gözləyir).
            overrides["output_bounds"] = (rule.minimum, rule.maximum)
        if rule.log_transform and strategy.transform.is_identity:
            overrides["variable_type"] = VariableType.LOGNORMAL
            overrides["transform"] = LogTransform()
            overrides["legacy_log_transform"] = True
        elif (not rule.log_transform and not strategy.transform.is_identity
              and strategy.legacy_log_transform):
            overrides["variable_type"] = VariableType.CONTINUOUS
            overrides["transform"] = IDENTITY_TRANSFORM
            overrides["back_transform"] = BackTransform.MEDIAN
            overrides["legacy_log_transform"] = False
        return strategy.derive(**overrides) if overrides else strategy

    def _kriging_overrides(self) -> Optional[Dict[str, object]]:
        """`self.interpolator` (istifadəçinin/UI-nin konstruktora ötürdüyü
        `IPropertyInterpolator`) `OrdinaryKriging` olub AÇIQ (defoltdan
        FƏRQLİ) dəyişdirdiyi sahələri Phase B-nin `kriging_overrides`
        mexanizminə körpüləyir (bax `property_interpolation.
        interpolate_property_field`).

        YALNIZ defoltdan FƏRQLİ sahələr ötürülür — əks halda, məs.
        istifadəçi HEÇ TOXUNMAYIB deyə `auto_fit=False`/`model=
        "spherical"` kimi `OrdinaryKriging()`-in xam defoltları
        strategiyanın öz ("auto" model fit, 24-qonşulu yerli axtarış)
        AĞILLI defoltlarını SƏSSİZCƏ əzərdi (B-INTEGRATION-FIX §12/§13).
        Əksinə, `range_`/`range_v`/`nugget`/`search_radius`/... kimi
        AÇIQ istifadəçi seçimi (məs. UI panelindəki Kriging sahələri)
        HƏMİŞƏ ötürülür — bu, variogram parametrinin dəyişməsinin
        nəticəyə TƏSİR ETDİYİNİ sübut edən test üçün vacibdir (§12)."""
        interp = self.interpolator
        if not isinstance(interp, OrdinaryKriging):
            return None
        defaults = OrdinaryKriging()
        fields = ("range_", "range_v", "sill", "nugget", "model", "auto_fit",
                 "auto_fit_nugget", "azimuth_deg", "range_minor", "dip_deg",
                 "search_radius", "max_neighbors", "min_neighbors", "sectors",
                 "honor_hard_data",
                 # PHASE C: istifadəçi (məs. `OrdinaryKriging(auto_detect_
                 # anisotropy=True)`) data-əsaslı anizotropluq aşkarlanmasına
                 # AÇIQ opt-in edə bilsin — defolt (`False`) DƏYİŞMİR, izotrop
                 # qalır (bax C§6 — "mümkün olduqda data-driven", MƏCBURİ deyil).
                 "auto_detect_anisotropy")
        overrides = {f: getattr(interp, f) for f in fields
                    if getattr(interp, f) != getattr(defaults, f)}
        return overrides or None

    @staticmethod
    def _attach_continuous_uncertainty(model: GeologicalModel, target: str,
                                       estimates: "list[PropertyEstimate]") -> None:
        """Hər laydan gələn `PropertyEstimate`-ləri BİR grid halına yığıb
        `model.uncertainty[target]`-ə yazır (B-INTEGRATION-FIX §9) —
        qeyri-müəyyənlik `report.warnings`-ə YAZILAN mətnlə YANAŞI, ƏDƏDİ
        massiv kimi DƏ istehsalat modelindən əlçatan olur, itmir."""
        if not estimates:
            return
        model.add_uncertainty(target, PropertyUncertainty(
            name=target,
            variance=np.concatenate([e.variance for e in estimates]),
            std=np.concatenate([e.std for e in estimates]),
            confidence=np.concatenate(
                [np.asarray(e.confidence, dtype=object) for e in estimates]),
            support=np.concatenate(
                [np.asarray(e.support, dtype=object) for e in estimates]),
            neighbor_count=np.concatenate([e.neighbor_count for e in estimates]),
            nearest_distance=np.concatenate([e.nearest_distance for e in estimates]),
            data_density=np.concatenate([e.data_density for e in estimates]),
            extrapolated=np.concatenate([e.extrapolated for e in estimates]),
            variance_kind=estimates[0].variance_kind.value,
            warnings=[w for e in estimates for w in e.warnings]))

    def _interpolate_volume(self, dataset, source, rule, targets, grid,
                            allow_cross_layer_fallback, report, geometry,
                            model: GeologicalModel,
                            calibrated_strategy: Optional[PropertyStrategy] = None
                            ) -> np.ndarray:
        """Hər təbəqə üçün Phase B xassə-strategiyalı 3D (X,Y,Z) kriginq,
        sonra həcmə yığılır (B-INTEGRATION-FIX — KÖHNƏ Phase A
        `geology.interpolation.interpolate_property()` ARTIQ burada
        ÇAĞIRILMIR, bax modul docstring-i).

        `dataset.samples_for(source, layer=k)` düzgün süzür: laya bağlı
        olmayan (`sample.layer is None`) nümunələr HƏR K üçün daxil
        edilir, yalnız k-ya bağlı nümunələr öz K-sına məhdudlaşır. Əgər
        bir K-də NƏ laysız, NƏ DƏ ona bağlı nümunə yoxdursa, defolt
        olaraq bu "digər layların nöqtələrini sükutla hovuzla" demək
        DEYİL (əvvəlki nöqsan, bax M1) — açıq xəta atılır, yalnız
        `allow_cross_layer_fallback=True` ilə bilərəkdən hovuzlanır.

        Hər nümunənin Z-si (bax `_sample_depth`) və hədəfin öz K-sının
        həqiqi hüceyrə-mərkəzi dərinliyi HƏR ZAMAN ötürülür (Phase B
        mühərriki HƏMİŞƏ `OrdinaryKriging` — bax `property_config.
        PropertyStrategy.interpolation` defoltu — ona görə 3D/anizotrop
        dəstək artıq `self.interpolator`-un növündən ASILI DEYİL, M2
        davranışı BÜTÜN qeydiyyatlı xassələrə YAYILIB). Bunun sayəsində
        `allow_cross_layer_fallback` işə düşəndə fərqli laylardan gələn
        nöqtələr artıq BƏRABƏR yox, öz dərinlik fərqlərinə görə
        (range_v vasitəsilə) ÇƏKİLİ qatılır — yaxın lay uzaq laydan çox
        təsir edir (M2: geoloji cəhətdən əsaslandırılmış borclanma).
        """
        strategy = (calibrated_strategy if calibrated_strategy is not None
                   else self._resolve_property_strategy(source))
        overrides = self._kriging_overrides()

        depths_grid = geometry.cell_depths().reshape(grid.shape)   # (nz, ny, nx)
        layer_mean_depth = depths_grid.mean(axis=(1, 2))
        target_depths = depths_grid

        layers = []
        estimates: "list[PropertyEstimate]" = []
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

            estimate = self._estimate_layer(samples, source, k, strategy, overrides,
                                            targets, target_depths, layer_mean_depth)
            layer_label = f" (K={k})" if grid.nz > 1 else ""
            for message in estimate.warnings:
                report.warn(f"'{source}'{layer_label}: {message}")
            layers.append(estimate.estimate)
            estimates.append(estimate)

        self._attach_continuous_uncertainty(model, rule.target, estimates)
        return np.concatenate(layers)

    def _estimate_layer(self, samples, source: str, k: int,
                        strategy: PropertyStrategy, overrides, targets: np.ndarray,
                        depths_grid: np.ndarray,
                        layer_mean_depth: np.ndarray) -> "PropertyEstimate":
        """BİR K-təbəqəsi üçün Phase B interpolyasiyası — HƏM köhnə
        (`_interpolate_volume`), HƏM DƏ lay-məlumatlı yol EYNİ bu metodu
        çağırır (Kriging riyaziyyatı BİR YERDƏDİR, təkrarlanmır)."""
        values = np.asarray([s.values[source] for s in samples], float)
        depths = np.asarray(
            [self._sample_depth(s, k, layer_mean_depth) for s in samples], float)
        points = np.column_stack(
            [[s.x for s in samples], [s.y for s in samples], depths])
        target_points = np.column_stack([targets, depths_grid[k].ravel()])
        return interpolate_property_field(points, values, target_points,
                                          strategy=strategy, kriging_overrides=overrides)

    # ═══════════════════════════════════════ LAY-MƏLUMATLI (layer-aware) YOL
    def _prepare_availability(self, dataset: WellDataset, geometry: CellGeometry,
                              sources: Sequence[str],
                              config: LayerInterpolationConfig,
                              report: "InterpolationReport") -> ModelDataAvailability:
        """Giriş mənzərəsi: hansı xassə hansı layda HƏQİQƏTƏN ölçülüb."""
        availability = compute_availability(dataset, geometry, config.policy, sources)
        for message in dataset.warnings:
            report.warn(message)
        stray = unassigned_samples(dataset, geometry, config.policy)
        if stray:
            names = ", ".join(sorted({s.well for s in stray}))
            report.warn(
                f"{len(stray)} quyu nöqtəsi ({names}) heç bir K-laya aid edilə bilmədi "
                "(nə lay indeksi, nə dərinlik, nə də 'Data layları' bəyanı var) — bu "
                "nöqtələr HEÇ BİR laya məlumat vermir. Səssiz yayılma tətbiq edilmir "
                "(bax LayerDataPolicy).")
        return availability

    def _record_categorical_availability(self, dataset: WellDataset,
                                         geometry: CellGeometry,
                                         sources: Sequence[str],
                                         config: LayerInterpolationConfig,
                                         availability: ModelDataAvailability,
                                         report: "InterpolationReport") -> None:
        """KATEQORİK sütunlar üçün lay-üzrə mövcudluğu QEYDƏ ALIR.

        Kateqorik yola (SIS / indikator kriginq) LAY MASKASI TƏTBİQ
        EDİLMİR — bu, bilərəkdəndir (§15: kəsilməz interpolyasiya
        qaydalarını fasiyaya KOR-KORANƏ tətbiq etmə; SIS onsuz da 3D
        şərtlənmiş STOXASTİK üsuldur və şaquli kəsilməzlik onun öz
        modelindən gəlir).

        Amma "sərt datası olmayan layda SIS nəticəsi var" faktı GİZLİ
        QALMAMALIDIR: həmin laylar `SIMULATED` kimi qeyd olunur (ölçülmüş
        KİMİ YOX) və hesabatda AÇIQ xəbərdarlıq verilir.
        """
        if not sources:
            return
        categorical = compute_availability(dataset, geometry, config.policy, sources)
        for source in sources:
            entry = categorical[source]
            simulated = [k for k in range(geometry.grid.nz) if not entry.layers[k].has_data]
            for k in simulated:
                entry.set(k, status=DataStatus.SIMULATED, method="sis",
                          note="bu layda sərt fasiya datası yoxdur — stoxastik realizasiya")
            for k in entry.data_layers():
                entry.set(k, status=DataStatus.SIMULATED, method="sis",
                          note="sərt data ilə şərtlənmiş realizasiya")
            availability.properties[source] = entry
            if simulated:
                report.warn(
                    f"'{source}' (kateqorik): L{format_layers(simulated)} laylarında "
                    "sərt fasiya datası YOXDUR — nəticə SIS realizasiyasıdır "
                    "(status SIMULATED), ölçmə DEYİL. Kateqorik sütuna lay maskası "
                    "tətbiq edilmir (bax §15).")

    @staticmethod
    def _layer_sample_index(dataset: WellDataset, source: str,
                            geometry: CellGeometry,
                            policy: LayerDataPolicy) -> Dict[int, list]:
        """`{K: [həmin laya məlumat VERƏN nümunələr]}` — BİR keçidlə.

        `WellDataset.samples_for()`-dan FƏRQİ: o, `layer is None` olan HƏR
        nümunəni HƏR laya daxil edir (köhnə "hər yerə aiddir" semantikası).
        Lay-məlumatlı rejimdə bu, məhz QADAĞAN olunan səssiz yayılmadır —
        ona görə burada `layer_availability.sample_layers()` işlədilir.

        İndeks BİR DƏFƏ qurulur (O(n·nz)) və bütün laylar/tamamlama
        addımları onu paylaşır — hər lay üçün dataset-i yenidən gəzmək
        (O(n·nz²)) böyük NZ-də lazımsız yükdür (§24).
        """
        from ..geology.layer_availability import sample_layers
        index: Dict[int, list] = {k: [] for k in range(geometry.grid.nz)}
        for sample in dataset.samples:
            value = sample.values.get(source)
            if value is None or not np.isfinite(value):
                continue
            for k in sample_layers(sample, geometry, policy):
                index[k].append(sample)
        return index

    def _build_layer_aware_field(self, dataset: WellDataset, source: str,
                                 rule: PropertyRule, targets: np.ndarray,
                                 grid: CartesianGrid, geometry: CellGeometry,
                                 model: GeologicalModel, report: "InterpolationReport",
                                 availability: ModelDataAvailability,
                                 config: LayerInterpolationConfig,
                                 sgs: Optional[ContinuousSGSConfig],
                                 calibrated_strategy: Optional[PropertyStrategy]
                                 ) -> np.ndarray:
        """`final` sahəni AÇIQ, izlənə bilən addımlarla qurur (§11):

            final ← original (varsa)              → status PRESERVED
            final[seçilmiş VƏ məlumatlı laylar]   ← interpolyasiya  → INTERPOLATED
            final[sərt data hüceyrələri]          ← ölçmə            → MEASURED
            final[məlumatsız laylar]              ← completion       → ESTIMATED/
                                                    (seçilibsə)        EXTRAPOLATED/
                                                                       SIMULATED
            qalan                                 ← toxunulmur       → MISSING

        HEÇ BİR ADDIM digərinin nəticəsini SÜKUTLA əzmir; hər hüceyrənin
        son statusu `PropertyProvenance`-da saxlanılır.
        """
        nz, ncell = grid.nz, grid.ncell
        areal = grid.nx * grid.ny
        entry = availability.require(source)
        data_layers = entry.data_layers()

        requested = config.targets_for(source, nz)
        if requested is None:
            requested = list(data_layers)
        interp_layers = [k for k in requested if k in data_layers]
        requested_without_data = [k for k in requested if k not in data_layers]
        if requested_without_data:
            report.warn(
                f"'{source}': seçilmiş laylardan L{format_layers(requested_without_data)} "
                "üçün sərt data YOXDUR — bu laylar İNTERPOLYASİYA EDİLMİR. Onlara yalnız "
                "AÇIQ tamamlama (completion) strategiyası tətbiq oluna bilər; "
                "seçilməyibsə MISSING qalırlar.")
        excluded_with_data = [k for k in data_layers if k not in interp_layers]
        if excluded_with_data:
            report.warn(
                f"'{source}': L{format_layers(excluded_with_data)} laylarında data VAR, "
                "amma istifadəçi onları interpolyasiya hədəfinə salmayıb — həmin laylar "
                "DƏYİŞDİRİLMİR (orijinal/MISSING olduğu kimi qalır).")

        original = config.original_fields.get(source)
        if original is not None:
            original = np.asarray(original, float).ravel()
            if original.size != ncell:
                raise ValueError(
                    f"'{source}': original_fields ölçüsü grid ilə uyğun gəlmir "
                    f"({original.size} != {ncell}).")

        final = np.full(ncell, np.nan)
        status = np.full(ncell, DataStatus.MISSING.value, dtype=object)
        method = np.full(ncell, "", dtype=object)
        confidence = np.full(ncell, np.nan)
        interpolated_field = np.full(ncell, np.nan)
        estimated_field = np.full(ncell, np.nan)
        layer_methods: Dict[int, str] = {}

        if original is not None and config.preserve_original_when_missing:
            final = original.copy()
            status[:] = DataStatus.PRESERVED.value
            method[:] = "original"

        # ── 1. seçilmiş VƏ məlumatlı laylar: mövcud Kriging boru xətti ──
        strategy = (calibrated_strategy if calibrated_strategy is not None
                    else self._resolve_property_strategy(source))
        overrides = self._kriging_overrides()
        depths_grid = geometry.cell_depths().reshape(grid.shape)
        layer_mean_depth = depths_grid.mean(axis=(1, 2))
        measured = hard_data_cells(dataset, geometry, source,
                                   config.policy).reshape(grid.shape)
        # lay → nümunə indeksi BİR DƏFƏ qurulur və bütün addımlar
        # (interpolyasiya + tamamlama) onu paylaşır (§24).
        layer_index = self._layer_sample_index(dataset, source, geometry, config.policy)
        estimates: Dict[int, PropertyEstimate] = {}

        for k in interp_layers:
            samples = layer_index[k]
            if not samples:                      # mövcudluq hesabı ilə ziddiyyət
                raise ValueError(
                    f"'{source}': K={k} üçün mövcudluq cədvəli data göstərir, amma "
                    "nümunə tapılmadı — daxili uyğunsuzluq.")
            if sgs is not None:
                continue                          # SGS aşağıda TOPLU işlənir
            estimate = self._estimate_layer(samples, source, k, strategy, overrides,
                                            targets, depths_grid, layer_mean_depth)
            for message in estimate.warnings:
                report.warn(f"'{source}' (K={k}): {message}")
            index = np.arange(k * areal, (k + 1) * areal)
            final[index] = estimate.estimate
            interpolated_field[index] = estimate.estimate
            status[index] = DataStatus.INTERPOLATED.value
            method[index] = "kriging"
            confidence[index] = _confidence_scores(estimate.confidence)
            hard = index[measured[k].ravel()]
            status[hard] = DataStatus.MEASURED.value
            method[hard] = "measured"
            confidence[hard] = 1.0
            estimates[k] = estimate
            layer_methods[k] = "kriging"
            entry.set(k, status=DataStatus.INTERPOLATED, method="kriging",
                      confidence=_mean_or_none(confidence[index]))

        if sgs is not None and interp_layers:
            values = self._simulate_continuous_sgs_field(
                dataset, source, targets, grid, geometry, sgs, report, model,
                layers=interp_layers)
            for position, k in enumerate(interp_layers):
                index = np.arange(k * areal, (k + 1) * areal)
                final[index] = values[position * areal:(position + 1) * areal]
                interpolated_field[index] = final[index]
                status[index] = DataStatus.SIMULATED.value
                method[index] = "sgs"
                hard = index[measured[k].ravel()]
                status[hard] = DataStatus.MEASURED.value
                method[hard] = "measured"
                confidence[hard] = 1.0
                layer_methods[k] = "sgs"
                entry.set(k, status=DataStatus.SIMULATED, method="sgs",
                          note=f"realization={sgs.realization_id}, seed={sgs.seed}")

        if estimates:
            self._attach_partial_uncertainty(model, rule.target, estimates, grid)

        # ── 2. məlumatsız laylar: YALNIZ AÇIQ completion ────────────────
        completion_layers = [k for k in range(nz)
                             if k not in interp_layers and k not in data_layers]
        spec = config.completion_for(source)
        if spec.layers is not None:
            allowed = {int(k) for k in spec.layers}
            skipped = [k for k in completion_layers if k not in allowed]
            completion_layers = [k for k in completion_layers if k in allowed]
            if skipped:
                report.warn(
                    f"'{source}': L{format_layers(skipped)} layları tamamlama "
                    "siyahısına (CompletionSpec.layers) daxil deyil — MISSING qalır.")

        if completion_layers:
            self._complete_missing_layers(
                dataset, source, spec, completion_layers, data_layers, layer_index,
                grid, geometry, targets, strategy, overrides, depths_grid,
                layer_mean_depth, config, original, model, report, entry,
                final, status, method, confidence, estimated_field, layer_methods)

        # ── 3. provenance + mövcudluq cədvəlinin yekunlaşdırılması ──────
        for k in range(nz):
            if k in interp_layers or k in completion_layers:
                continue
            index = np.arange(k * areal, (k + 1) * areal)
            # İKİ FƏRQLİ SƏBƏB, İKİ FƏRQLİ MESAJ (səbəb gizlədilmir):
            #   · layda data VAR, amma istifadəçi onu hədəfə salmayıb;
            #   · layda data YOXDUR və tamamlama seçilməyib.
            reason = ("məlumat var, amma interpolyasiya hədəfinə daxil deyil"
                      if k in data_layers else
                      "məlumat yoxdur, tamamlama strategiyası seçilməyib")
            if original is not None and config.preserve_original_when_missing:
                entry.set(k, status=DataStatus.PRESERVED, method="original",
                          note=reason)
            else:
                status[index] = DataStatus.MISSING.value
                entry.set(k, status=DataStatus.MISSING, method="", note=reason)

        model.add_provenance(PropertyProvenance(
            name=rule.target, status=status, method=method, confidence=confidence,
            final=final, original=original, interpolated=interpolated_field,
            estimated=estimated_field, layer_methods=layer_methods))
        return final

    def _complete_missing_layers(self, dataset, source, spec: CompletionSpec,
                                 completion_layers, data_layers, layer_index,
                                 grid, geometry, targets, strategy, overrides,
                                 depths_grid, layer_mean_depth, config, original,
                                 model, report, entry: PropertyAvailability,
                                 final, status, method, confidence, estimated_field,
                                 layer_methods) -> None:
        """Tapşırıq §9-un tamamlama strategiyaları. HEÇ BİRİ `MEASURED`
        vermir və heç biri DEFOLT deyil — `CompletionMethod.NONE` ilə
        laylar toxunulmadan MISSING/PRESERVED qalır."""
        areal = grid.nx * grid.ny
        label = format_layers(completion_layers)

        if spec.method is CompletionMethod.NONE:
            preserved = original is not None and config.preserve_original_when_missing
            for k in completion_layers:
                # Hüceyrə-səviyyəli status ARTIQ doğrudur (PRESERVED və ya
                # MISSING); burada YALNIZ mövcudluq cədvəli eyni həqiqəti
                # təkrarlayır ki, iki mənbə bir-birinə ZİDD OLMASIN.
                entry.set(k,
                          status=DataStatus.PRESERVED if preserved else DataStatus.MISSING,
                          method="original" if preserved else "",
                          note=("orijinal sahə saxlanıldı — tamamlama seçilməyib"
                                if preserved else
                                "məlumat yoxdur, tamamlama strategiyası seçilməyib"))
            report.warn(
                f"'{source}': L{label} üçün sərt data YOXDUR və tamamlama strategiyası "
                + ("SEÇİLMƏYİB — orijinal sahə OLDUĞU KİMİ saxlanıldı (§10 defolt "
                   "davranışı: preserve)."
                   if preserved else
                   "SEÇİLMƏYİB — bu laylar MISSING olaraq qalır (§10 defolt davranışı). "
                   "Simulyasiyadan əvvəl ya data əlavə edin, ya da completion seçin."))
            return

        if spec.method is CompletionMethod.PRESERVE_ORIGINAL:
            if original is None:
                raise ValueError(
                    f"'{source}': completion='preserve_original' seçilib, amma "
                    "`original_fields` verilməyib — saxlanacaq orijinal sahə yoxdur.")
            for k in completion_layers:
                index = np.arange(k * areal, (k + 1) * areal)
                final[index] = original[index]
                status[index] = DataStatus.PRESERVED.value
                method[index] = "preserve_original"
                confidence[index] = (np.nan if spec.confidence is None
                                     else float(spec.confidence))
                layer_methods[k] = "preserve_original"
                entry.set(k, status=DataStatus.PRESERVED, method="preserve_original",
                          confidence=spec.confidence,
                          note="mövcud geoloji prior saxlanıldı")
            return

        if spec.method is CompletionMethod.CONSTANT:
            if spec.value is None or not np.isfinite(spec.value):
                raise ValueError(
                    f"'{source}': completion='constant' üçün `value` verilməlidir "
                    "(sonlu ədəd).")
            value = float(spec.value)
            for k in completion_layers:
                index = np.arange(k * areal, (k + 1) * areal)
                final[index] = value
                estimated_field[index] = value
                status[index] = DataStatus.ESTIMATED.value
                method[index] = f"constant={value:g}"
                confidence[index] = (np.nan if spec.confidence is None
                                     else float(spec.confidence))
                layer_methods[k] = "constant"
                entry.set(k, status=DataStatus.ESTIMATED, method=f"constant={value:g}",
                          confidence=spec.confidence,
                          note="istifadəçinin AÇIQ verdiyi lay dəyəri")
            report.warn(
                f"'{source}': L{label} istifadəçinin verdiyi sabit dəyərlə ({value:g}) "
                "tamamlandı — status ESTIMATED, ölçmə DEYİL.")
            return

        if spec.method is CompletionMethod.VERTICAL_TREND:
            self._complete_by_vertical_trend(
                source, spec, completion_layers, data_layers, layer_index, grid,
                strategy, layer_mean_depth, report, entry,
                final, status, method, confidence, estimated_field, layer_methods)
            return

        if spec.method is CompletionMethod.GEOSTATISTICAL_3D:
            self._complete_by_3d_kriging(
                source, spec, completion_layers, data_layers, layer_index, grid,
                targets, strategy, overrides, depths_grid, layer_mean_depth,
                report, entry, final, status, method, confidence, estimated_field,
                layer_methods)
            return

        if spec.method is CompletionMethod.SGS:
            self._complete_by_sgs(
                dataset, source, spec, completion_layers, grid, geometry, targets,
                model, report, entry, final, status, method, confidence,
                estimated_field, layer_methods)
            return

        raise ValueError(f"'{source}': naməlum completion üsulu {spec.method!r}.")

    # ---------------------------------------------------- şaquli trend
    def _complete_by_vertical_trend(self, source, spec, completion_layers,
                                    data_layers, layer_index, grid, strategy,
                                    layer_mean_depth, report, entry,
                                    final, status, method, confidence,
                                    estimated_field, layer_methods) -> None:
        """Məlumatlı layların LAY-ORTALARINDAN dərinliyə görə xətti trend.

        NƏ EDİLİR: `m(z) = a + b·z` (çevrilmiş fəzada — PERMX üçün loq)
        məlumatlı layların ortaları ilə ən-kiçik-kvadratlarla qurulur,
        məlumatsız laya YALNIZ ONUN ORTASI verilir.

        NƏ EDİLMİR (§26): qonşu layın LATERAL xəritəsi KOPYALANMIR. Trend
        şaquli məlumat verir, üfüqi struktur haqqında MƏLUMAT VERMİR —
        onu uydurmaq elmi cəhətdən müdafiə oluna bilməz. Bu, nəticənin
        lay daxilində SABİT olması deməkdir və hesabatda AÇIQ deyilir.
        """
        points = []
        for k in data_layers:
            samples = layer_index[k]
            values = np.asarray([s.values[source] for s in samples], float)
            if values.size == 0:
                continue
            try:
                transformed = strategy.transform.forward(values)
            except Exception as error:                       # TransformError və s.
                raise ValueError(
                    f"'{source}': şaquli trend üçün çevirmə uğursuz oldu (K={k}): {error}"
                ) from error
            points.append((float(layer_mean_depth[k]), float(np.mean(transformed))))

        if len(points) < 2:
            raise ValueError(
                f"'{source}': şaquli trend üçün ən azı İKİ məlumatlı lay lazımdır "
                f"(tapıldı: {len(points)}). Tək layla trend qurmaq onu SABİTƏ çevirir "
                "— bu, trend deyil, gizli ekstrapolyasiya olardı.")

        depths = np.asarray([p[0] for p in points], float)
        means = np.asarray([p[1] for p in points], float)
        slope, intercept = np.polyfit(depths, means, 1)
        low, high = float(depths.min()), float(depths.max())
        lo_bound, hi_bound = strategy.output_bounds
        range_v = self._vertical_range()
        areal = grid.nx * grid.ny

        for k in completion_layers:
            z = float(layer_mean_depth[k])
            value = float(strategy.transform.inverse(
                np.asarray([intercept + slope * z], float))[0])
            if lo_bound is not None:
                value = max(value, lo_bound)
            if hi_bound is not None:
                value = min(value, hi_bound)
            inside = low - 1e-9 <= z <= high + 1e-9
            gap = 0.0 if inside else min(abs(z - low), abs(z - high))
            score = self._extrapolation_confidence(gap, range_v, inside)
            if spec.confidence is not None:
                score = float(spec.confidence)
            index = np.arange(k * areal, (k + 1) * areal)
            final[index] = value
            estimated_field[index] = value
            state = DataStatus.ESTIMATED if inside else DataStatus.EXTRAPOLATED
            status[index] = state.value
            method[index] = "vertical_trend"
            confidence[index] = np.nan if score is None else score
            layer_methods[k] = "vertical_trend"
            entry.set(k, status=state, method="vertical_trend", confidence=score,
                      note=(f"trend: {strategy.transform.describe()} fəzasında "
                            f"b={slope:.4g}/m; şaquli məsafə {gap:.1f} m"))
        report.warn(
            f"'{source}': L{format_layers(completion_layers)} şaquli trendlə tamamlandı "
            "— lay daxilində SABİT dəyər (lateral struktur UYDURULMUR), status "
            "ESTIMATED/EXTRAPOLATED, ölçmə DEYİL."
            + ("" if range_v else " Şaquli korrelyasiya radiusu (range_v) verilmədiyi "
                                 "üçün etibarlılıq balı HESABLANMADI (NaN) — saxta "
                                 "rəqəm yaradılmır."))

    # ------------------------------------------------- 3D geostatistika
    def _complete_by_3d_kriging(self, source, spec, completion_layers,
                                data_layers, layer_index, grid, targets, strategy,
                                overrides, depths_grid, layer_mean_depth,
                                report, entry, final, status, method, confidence,
                                estimated_field, layer_methods) -> None:
        """BÜTÜN məlumatlı layların sərt datası ilə, HƏQİQİ 3D (X,Y,Z)
        məsafə üzərindən həmin laya qiymət verir — mövcud Kriging
        mühərriki (`interpolate_property_field`) DƏYİŞMƏDƏN işlədilir.

        Bu, köhnə `allow_cross_layer_fallback` davranışının ELMİ CƏHƏTDƏN
        DÜZGÜN, ETİKETLƏNMİŞ variantıdır: nəticə `INTERPOLATED` deyil,
        `ESTIMATED`/`EXTRAPOLATED` kimi qeyd olunur və provenance-da
        görünür."""
        samples = []
        for k in data_layers:
            samples.extend(layer_index[k])
        if len(samples) < 2:
            raise ValueError(
                f"'{source}': 3D geostatistik tamamlama üçün ən azı iki sərt nöqtə "
                f"lazımdır (tapıldı: {len(samples)}).")
        depths = [float(layer_mean_depth[k]) for k in data_layers]
        low, high = min(depths), max(depths)
        range_v = self._vertical_range()
        areal = grid.nx * grid.ny

        for k in completion_layers:
            estimate = self._estimate_layer(samples, source, k, strategy, overrides,
                                            targets, depths_grid, layer_mean_depth)
            for message in estimate.warnings:
                report.warn(f"'{source}' (tamamlama K={k}): {message}")
            z = float(layer_mean_depth[k])
            inside = low - 1e-9 <= z <= high + 1e-9
            state = DataStatus.ESTIMATED if inside else DataStatus.EXTRAPOLATED
            index = np.arange(k * areal, (k + 1) * areal)
            final[index] = estimate.estimate
            estimated_field[index] = estimate.estimate
            status[index] = state.value
            method[index] = "geostatistical_3d"
            scores = _confidence_scores(estimate.confidence)
            if not inside:
                gap = min(abs(z - low), abs(z - high))
                penalty = self._extrapolation_confidence(gap, range_v, False)
                scores = (np.full(scores.size, np.nan) if penalty is None
                          else np.minimum(scores, penalty))
            if spec.confidence is not None:
                scores = np.full(scores.size, float(spec.confidence))
            confidence[index] = scores
            layer_methods[k] = "geostatistical_3d"
            entry.set(k, status=state, method="geostatistical_3d",
                      confidence=_mean_or_none(scores),
                      note="bütün məlumatlı layların 3D sərt datası ilə")
        report.warn(
            f"'{source}': L{format_layers(completion_layers)} 3D geostatistik "
            "qiymətləndirmə ilə tamamlandı — status ESTIMATED/EXTRAPOLATED, "
            "ÖLÇMƏ DEYİL.")

    # ---------------------------------------------------------- SGS
    def _complete_by_sgs(self, dataset, source, spec, completion_layers, grid,
                         geometry, targets, model, report, entry, final, status,
                         method, confidence, estimated_field, layer_methods) -> None:
        """Mövcud SGS mühərriki (`geology/sgs.py`) DƏYİŞMƏDƏN çağırılır;
        yalnız HƏDƏF HÜCEYRƏLƏR məlumatsız laylarla məhdudlaşdırılır
        (§24: lazımsız hesablama yoxdur). Nəticə HƏMİŞƏ `SIMULATED`."""
        config = spec.sgs if spec.sgs is not None else ContinuousSGSConfig()
        values = self._simulate_continuous_sgs_field(
            dataset, source, targets, grid, geometry, config, report, model,
            layers=completion_layers)
        areal = grid.nx * grid.ny
        for position, k in enumerate(completion_layers):
            index = np.arange(k * areal, (k + 1) * areal)
            block = values[position * areal:(position + 1) * areal]
            final[index] = block
            estimated_field[index] = block
            status[index] = DataStatus.SIMULATED.value
            method[index] = "sgs"
            confidence[index] = (np.nan if spec.confidence is None
                                 else float(spec.confidence))
            layer_methods[k] = "sgs"
            entry.set(k, status=DataStatus.SIMULATED, method="sgs",
                      confidence=spec.confidence,
                      note=(f"realization={config.realization_id}, seed={config.seed}"
                            " — TƏK realizasiya; qeyri-müəyyənlik üçün ansambl lazımdır"))
        report.warn(
            f"'{source}': L{format_layers(completion_layers)} SGS realizasiyası ilə "
            f"tamamlandı (seed={config.seed}, realization={config.realization_id}) — "
            "status SIMULATED, ÖLÇMƏ DEYİL. Tək realizasiyadan etibarlılıq balı "
            "hesablanmır (NaN); ansambl üçün `sgs_ensemble` işlədin.")

    # --------------------------------------------------------- köməkçilər
    def _vertical_range(self) -> Optional[float]:
        """Şaquli korrelyasiya radiusu (m) — YALNIZ istifadəçi AÇIQ verəndə.

        `OrdinaryKriging(range_v=...)` verilməyibsə `None` qaytarır və
        çağıran etibarlılıq balını HESABLAMIR (§18: əsassız rəqəm
        yaradılmır)."""
        interpolator = self.interpolator
        if isinstance(interpolator, OrdinaryKriging) and interpolator.range_v:
            return float(interpolator.range_v)
        return None

    @staticmethod
    def _extrapolation_confidence(gap: float, range_v: Optional[float],
                                  inside: bool) -> Optional[float]:
        """Şaquli ekstrapolyasiya məsafəsinə görə ORDİNAL etibarlılıq balı.

        `exp(-gap / range_v)` — variogram dəstəyinə ƏSASLANIR: şaquli
        korrelyasiya radiusu qədər uzaqlaşanda bal `1/e`-yə düşür. Zərfin
        İÇİNDƏ (interpolyasiya) `gap=0` → 1.0-a yaxın, amma sabit
        `_TREND_INSIDE_CEILING` ilə məhdudlaşdırılır ki, HEÇ VAXT ölçmə
        (1.0) ilə eyni səviyyəyə çıxmasın.

        `range_v` bilinmirsə `None` — SAXTA rəqəm YOXDUR (§18).
        """
        if range_v is None or range_v <= 0.0:
            return None
        score = float(np.exp(-max(gap, 0.0) / range_v))
        ceiling = _TREND_INSIDE_CEILING if inside else _TREND_OUTSIDE_CEILING
        return min(score, ceiling)

    @staticmethod
    def _attach_partial_uncertainty(model: GeologicalModel, target: str,
                                    estimates: Dict[int, "PropertyEstimate"],
                                    grid: CartesianGrid) -> None:
        """`_attach_continuous_uncertainty`-nin lay-məlumatlı variantı:
        YALNIZ hesablanmış laylar üçün nəticə var, qalan laylar NaN/boş
        qalır — uydurma qiymət YAZILMIR."""
        if not estimates:
            return
        ncell = grid.ncell
        areal = grid.nx * grid.ny
        sample = next(iter(estimates.values()))
        variance = np.full(ncell, np.nan)
        std = np.full(ncell, np.nan)
        confidence = np.full(ncell, "", dtype=object)
        support = np.full(ncell, "", dtype=object)
        neighbor_count = np.zeros(ncell, dtype=float)
        nearest = np.full(ncell, np.nan)
        density = np.zeros(ncell, dtype=float)
        extrapolated = np.zeros(ncell, dtype=bool)
        warnings: list = []
        for k, estimate in estimates.items():
            index = np.arange(k * areal, (k + 1) * areal)
            variance[index] = estimate.variance
            std[index] = estimate.std
            confidence[index] = np.asarray(estimate.confidence, dtype=object)
            support[index] = np.asarray(estimate.support, dtype=object)
            neighbor_count[index] = estimate.neighbor_count
            nearest[index] = estimate.nearest_distance
            density[index] = estimate.data_density
            extrapolated[index] = estimate.extrapolated
            warnings.extend(estimate.warnings)
        model.add_uncertainty(target, PropertyUncertainty(
            name=target, variance=variance, std=std, confidence=confidence,
            support=support, neighbor_count=neighbor_count, nearest_distance=nearest,
            data_density=density, extrapolated=extrapolated,
            variance_kind=sample.variance_kind.value, warnings=warnings))

    # ------------------------------------------------- PHASE C: model calibration
    def calibrate_property(self, dataset: WellDataset, source: str,
                           candidates=None, design=None, weights=None, qc=None):
        """`source` üçün variogram-model seçimini SIZMASIZ, MƏKAN-BLOKLU
        (spatial-block) çarpaz-doğrulama ilə aparır və `ModelSelectionReport`
        qaytarır (PHASE C §2/§3) — CV/model-selection infrastrukturu
        (`geology.cross_validation.select_property_model`) artıq mövcud
        idi, bura YALNIZ onu production builder-dən ÇAĞIRILA BİLƏN,
        DEFOLT olaraq spatial-block dizaynlı bir addım kimi körpüləyir.

        NƏTİCƏ AVTOMATİK production interpolyasiyasına TƏTBİQ OLUNMUR —
        `build(calibrated_strategies={source: report.selected.candidate.
        strategy})` ilə İSTİFADƏÇİ AÇIQ tətbiq edir (bax `build()`
        docstring-i). Bu, §19-un tələbini tərcümə edir: "CV/LOOCV production
        interpolyasiya yoluna DAXİL OLMAMALIDIR" — kalibrasiya BİR DƏFƏLİK,
        AYRICA addımdır, hər grid hüceyrəsi üçün TƏKRARLANMIR.

        Laylı dataset-də bütün laylardan HOVUZLANMIŞ nöqtələrlə (layer
        SÜZGƏCSİZ) kalibrasiya aparılır — model SEÇİMİ (hansı variogram
        AİLƏSİ ən yaxşı fit olur) layların cüzi fərqinə HƏSSAS statistik
        sual deyil, bu, `_interpolate_volume`-un öz laylı kriginqindən
        (hər K ÖZ sistemi ilə) FƏRQLİDİR və qəsdən sadələşdirilib."""
        from ..geology.cross_validation import ValidationDesign, ValidationKind, default_candidates, select_property_model

        samples = [s for s in dataset.samples if source in s.values]
        if not samples:
            raise ValueError(f"'{source}' üçün kalibrasiya ediləcək sərt data yoxdur.")
        points = np.asarray([(s.x, s.y) for s in samples], float)
        values = np.asarray([s.values[source] for s in samples], float)
        candidates = candidates if candidates is not None else default_candidates(source)
        design = design or ValidationDesign(kind=ValidationKind.SPATIAL_BLOCK)
        return select_property_model(points, values, candidates, property_name=source,
                                     design=design, weights=weights, qc=qc)

    # ------------------------------------------------------- M4: cross-validation
    def cross_validate(self, dataset: WellDataset, source: str,
                       method: str = "loo", k: int = 5, seed: int = 42,
                       nz: Optional[int] = None
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

        `nz` (LAY-MƏLUMATLI rejim, tapşırıq §19) — verilibsə BÜTÜN K
        təbəqələri gəzilir, təkcə dataset-də NÖQTƏSİ OLANLAR yox. Belə
        olanda məlumatı OLMAYAN lay (məs. L4/L5) hesabatda AÇIQ "doğrulama
        məlumatı yoxdur" kimi görünür — sükutla siyahıdan DÜŞMÜR və HEÇ
        VAXT "RMSE = 0" kimi SAXTA uğur nəticəsi yaratmır.
        """
        rule = self.rules.get(source, PropertyRule(source))
        if dataset.is_layered():
            layers: Sequence[Optional[int]] = (
                list(range(nz)) if nz is not None
                else [l for l in dataset.layers if l is not None])
        else:
            layers = [None]

        results: Dict[Optional[int], CrossValidationResult] = {}
        skipped: Dict[Optional[int], str] = {}
        for layer in layers:
            samples = dataset.samples_for(source, layer)
            if len(samples) < 3:
                label = f"K={layer}" if layer is not None else "bütün model"
                reason = ("bu layda HEÇ BİR doğrulama məlumatı yoxdur"
                          if not samples else
                          f"{len(samples)} nöqtə var, cross-validation üçün ən azı 3 lazımdır")
                skipped[layer] = (
                    f"{label}: {reason} — bu lay/model üçün doğrulama APARILMADI "
                    "(nəticə 'mükəmməl' kimi göstərilmir).")
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
                           k: int = 5, seed: int = 42, nz: Optional[int] = None
                           ) -> Dict[str, Tuple[Dict[Optional[int], CrossValidationResult],
                                               Dict[Optional[int], str]]]:
        """PORO və mövcud PERMX/PERMY/PERMZ üçün `cross_validate` icra edir.

        `nz` — bax `cross_validate` (məlumatsız laylar da hesabatda
        görünsün deyə)."""
        available = set(dataset.property_names())
        return {
            source: self.cross_validate(dataset, source, method=method, k=k, seed=seed,
                                        nz=nz)
            for source in _CROSS_VALIDATED_PROPERTIES if source in available
        }

    @staticmethod
    def _fill_missing_permeability(model, grid, ky_over_kx, kv_over_kh, report):
        """PERMY/PERMZ verilməyibsə anizotropluq əmsalları ilə qurulur.

        LAY-MƏLUMATLI rejimdə PERMX-in MƏNŞƏYİ (provenance) də TÖRƏMƏ
        sahələrə KEÇİRİLİR: PERMX-in MISSING olduğu hüceyrədə PERMY/PERMZ
        də MISSING-dir (NaN × əmsal = NaN) — status bunu AÇIQ göstərir,
        "əmsalla doldurduq, deməli məlumatlıdır" TƏƏSSÜRATI YARANMIR.
        """
        if "PERMX" not in model.property_maps:
            return
        permx = model.property_maps["PERMX"].values
        origin = model.provenance.get("PERMX")
        for key, factor in (("PERMY", ky_over_kx), ("PERMZ", kv_over_kh)):
            if key in model.property_maps:
                continue
            values = permx * factor
            model.add_property(PropertyMap.from_array(key, values, grid.ncell, "mD"))
            report.add(key, "PERMX", False, values)
            if origin is None:
                continue
            model.add_provenance(PropertyProvenance(
                name=key,
                status=np.array(origin.status, dtype=object),
                method=np.asarray([f"{m}→{key}(×{factor:g})" if m else ""
                                   for m in origin.method], dtype=object),
                confidence=np.array(origin.confidence, float),
                final=values,
                original=None if origin.original is None else origin.original * factor,
                interpolated=(None if origin.interpolated is None
                              else origin.interpolated * factor),
                estimated=(None if origin.estimated is None
                           else origin.estimated * factor),
                layer_methods=dict(origin.layer_methods)))
            if model.availability is not None and "PERMX" in model.availability:
                source_entry = model.availability["PERMX"]
                derived = model.availability.require(key)
                for k in range(grid.nz):
                    state = source_entry.layers[k]
                    derived.set(k, status=state.status, n_data=state.n_data,
                                confidence=state.confidence,
                                method=(f"{state.method}→{key}" if state.method else ""),
                                note=f"PERMX × {factor:g} (anizotropluq əmsalı)")


# ═══════════════════════════════════════════════════ TƏSİR (impact) ANALİZİ
@dataclass
class ImpactResult:
    """"Fərz edək ki..." ssenarisinin ƏSAS modelə TƏSİRİ (tapşırıq §12).

    QƏTİ QAYDA: bu obyekt YALNIZ HESABLAMA NƏTİCƏSİDİR — nə `original`,
    nə də `hypothetical` model DƏYİŞDİRİLİR (massivlər KOPYALANIR).
    Təsir HEÇ VAXT `final_field`-ə YAZILMIR; UI onu AYRICA təbəqə kimi
    göstərir.

    `delta = hypothetical − original` (fiziki vahiddə),
    `relative = delta / |original|` (orijinal sıfır/NaN olan hüceyrədə NaN).
    """

    name: str
    original: np.ndarray
    hypothetical: np.ndarray
    delta: np.ndarray
    relative: np.ndarray
    shape: Tuple[int, int, int]

    @property
    def changed_cells(self) -> int:
        return int(np.sum(np.abs(np.nan_to_num(self.delta, nan=0.0)) > 0.0))

    def layer_mean_delta(self) -> np.ndarray:
        """`(nz,)` — hər layın orta təsiri (NaN-lar nəzərə alınmır)."""
        grid = self.delta.reshape(self.shape)
        with np.errstate(invalid="ignore"):
            return np.asarray([
                float(np.nanmean(grid[k])) if np.isfinite(grid[k]).any() else np.nan
                for k in range(self.shape[0])], float)

    def as_text(self) -> str:
        lines = [f"TƏSİR — {self.name}:",
                 f"  dəyişən hüceyrə: {self.changed_cells}/{self.delta.size}"]
        for k, value in enumerate(self.layer_mean_delta()):
            lines.append(f"  L{k + 1} (K={k}): Δ orta = "
                         + ("—" if not np.isfinite(value) else f"{value:+.5g}"))
        return "\n".join(lines)


def compute_property_impact(original: GeologicalModel, hypothetical: GeologicalModel,
                            name: str) -> ImpactResult:
    """İki modelin EYNİ xassəsi arasındakı fərq.

    HEÇ BİR modeli DƏYİŞMİR (§12) — massivlər `copy()` ilə götürülür,
    ona görə nəticə üzərində aparılan hər hansı əməliyyat mənbəyə
    QAYITMIR (bax `tests/test_layer_data_availability.py`, TEST H).
    """
    if original.grid.shape != hypothetical.grid.shape:
        raise ValueError(
            f"Təsir analizi üçün grid ölçüləri eyni olmalıdır: "
            f"{original.grid.shape} != {hypothetical.grid.shape}")
    if name not in original.property_maps or name not in hypothetical.property_maps:
        raise KeyError(f"'{name}' hər iki modeldə olmalıdır (təsir analizi).")
    base = np.array(original.property_maps[name].values, float, copy=True)
    other = np.array(hypothetical.property_maps[name].values, float, copy=True)
    delta = other - base
    with np.errstate(divide="ignore", invalid="ignore"):
        relative = np.where(np.abs(base) > 0.0, delta / np.abs(base), np.nan)
    return ImpactResult(name=name, original=base, hypothetical=other, delta=delta,
                        relative=relative, shape=original.grid.shape)


def compute_impact(original: GeologicalModel, hypothetical: GeologicalModel,
                   names: Optional[Sequence[str]] = None) -> Dict[str, ImpactResult]:
    """Hər ortaq (və ya `names`-də sadalanan) xassə üçün `ImpactResult`."""
    shared = (list(names) if names is not None
              else sorted(set(original.property_maps) & set(hypothetical.property_maps)))
    return {name: compute_property_impact(original, hypothetical, name)
            for name in shared}


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


# ═══════════════════════════════════════════════ PHASE C: diagnostic report + gate
_SUPPORT_LABELS = {
    SUPPORT_WELL: "GOOD_SUPPORT", SUPPORT_BOUNDARY: "LIMITED_SUPPORT",
    SUPPORT_WEAK: "POOR_SUPPORT", SUPPORT_EXTRAPOLATED: "EXTRAPOLATION",
}


def _label_distribution(values: np.ndarray,
                        mapping: Optional[Dict[str, str]] = None) -> Dict[str, int]:
    if values is None or np.asarray(values).size == 0:
        return {}
    labels = [mapping.get(str(v), str(v)) if mapping else str(v) for v in values]
    unique, counts = np.unique(labels, return_counts=True)
    return {str(u): int(c) for u, c in zip(unique, counts)}


@dataclass
class PropertyQualityReport:
    """Bir xassənin PRODUCTION interpolyasiya diaqnostikası (PHASE C §24).

    `build_quality_report()` bunu `GeologicalModel.uncertainty`-dən
    (B-INTEGRATION-FIX-də doldurulan) qurur — heç bir yeni interpolyasiya
    aparmır, sadəcə artıq hesablanmış nəticəni İCMALLAŞDIRIR."""
    property_name: str
    kind: str                                    #: "continuous" | "categorical"
    strategy_summary: str
    sample_count: int
    neighbor_count_min: Optional[int]
    neighbor_count_mean: Optional[float]
    neighbor_count_max: Optional[int]
    mean_uncertainty: Optional[float]             #: kəsilməz: orta std; kateqorik: orta norm. entropiya
    max_uncertainty: Optional[float]
    extrapolated_fraction: float
    support_distribution: Dict[str, int] = field(default_factory=dict)
    confidence_distribution: Dict[str, int] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    def as_text(self) -> str:
        lines = [f"{self.property_name} ({self.kind}): {self.strategy_summary}",
                 f"  nümunə: {self.sample_count}"]
        if self.neighbor_count_mean is not None:
            lines.append(f"  qonşu sayı: min {self.neighbor_count_min} "
                         f"orta {self.neighbor_count_mean:.1f} maks {self.neighbor_count_max}")
        if self.mean_uncertainty is not None:
            lines.append(f"  qeyri-müəyyənlik: orta {self.mean_uncertainty:.4g} "
                         f"maks {self.max_uncertainty:.4g}")
        lines.append(f"  ekstrapolyasiya: {self.extrapolated_fraction * 100:.1f}%")
        if self.support_distribution:
            lines.append("  dəstək: " + ", ".join(
                f"{k}={v}" for k, v in sorted(self.support_distribution.items())))
        if self.confidence_distribution:
            lines.append("  etimad: " + ", ".join(
                f"{k}={v}" for k, v in sorted(self.confidence_distribution.items())))
        lines.extend(f"  ⚠ {w}" for w in self.warnings)
        return "\n".join(lines)


def build_quality_report(builder: "WellBasedGeologicalModelBuilder", model: GeologicalModel,
                         dataset: WellDataset) -> Dict[str, PropertyQualityReport]:
    """`model.uncertainty`-dən avtomatik diaqnostik hesabat (PHASE C §24).

    Yalnız Phase B mühərrikindən HƏQİQƏTƏN keçmiş xassələr (`model.
    uncertainty`-də olanlar) daxildir. Anizotropluq əmsalı ilə DOLDURULAN
    PERMY/PERMZ (bax `_fill_missing_permeability`) buraya DAXİL DEYİL —
    onlar interpolyasiya EDİLMƏYİB, bu ÖZÜ diaqnostik həqiqətdir,
    gizlədilmir (sadəcə "hesabatda yoxdur" deyil, `describe_missing_
    properties`-dən görünə bilər, bax aşağıda)."""
    reports: Dict[str, PropertyQualityReport] = {}
    for name, unc in model.uncertainty.items():
        sample_count = len(dataset.samples_for(name, None))
        neighbor = np.asarray(unc.neighbor_count, dtype=float)
        neighbor_min = int(neighbor.min()) if neighbor.size else None
        neighbor_mean = float(neighbor.mean()) if neighbor.size else None
        neighbor_max = int(neighbor.max()) if neighbor.size else None
        extrapolated_fraction = (float(np.mean(unc.extrapolated))
                                 if np.asarray(unc.extrapolated).size else 0.0)
        support_dist = _label_distribution(unc.support, _SUPPORT_LABELS)
        confidence_dist = _label_distribution(unc.confidence)

        if isinstance(unc, PropertyUncertainty):
            strategy = builder._resolve_property_strategy(name)
            finite_std = unc.std[np.isfinite(unc.std)]
            mean_unc = float(finite_std.mean()) if finite_std.size else None
            max_unc = float(finite_std.max()) if finite_std.size else None
            kind = "continuous"
        elif isinstance(unc, CategoricalUncertainty):
            strategy = resolve_strategy(name)
            entropy = unc.normalized_entropy
            mean_unc = float(np.mean(entropy)) if np.asarray(entropy).size else None
            max_unc = float(np.max(entropy)) if np.asarray(entropy).size else None
            kind = "categorical"
        else:
            continue

        reports[name] = PropertyQualityReport(
            property_name=name, kind=kind, strategy_summary=strategy.describe(),
            sample_count=sample_count, neighbor_count_min=neighbor_min,
            neighbor_count_mean=neighbor_mean, neighbor_count_max=neighbor_max,
            mean_uncertainty=mean_unc, max_uncertainty=max_unc,
            extrapolated_fraction=extrapolated_fraction,
            support_distribution=support_dist, confidence_distribution=confidence_dist,
            warnings=list(unc.warnings))
    return reports


def quality_report_as_text(reports: Dict[str, PropertyQualityReport]) -> str:
    if not reports:
        return "Diaqnostik hesabat üçün heç bir Phase B nəticəsi yoxdur."
    return "\n\n".join(r.as_text() for r in reports.values())


def run_validation_gate(model: GeologicalModel, report: InterpolationReport,
                        quality: Optional[Dict[str, PropertyQualityReport]] = None,
                        max_extrapolated_fraction: float = 0.5) -> DiagnosticReport:
    """Production model-in "validated" statusu almazdan ƏVVƏLKİ son qapı
    (PHASE C §25/26) — mövcud `domain.diagnostics.DiagnosticReport`/
    `Severity` semantikasını (XƏTA=bloklayır, XƏBƏRDARLIQ=bloklamır)
    TƏKRAR İCAD ETMİR, birbaşa işlədir.

    XƏTA (bloklayır): NaN/Inf property/fasiya dəyəri. `build()` artıq
    `model.validate()` ilə əsas fiziki hədd pozmalarını (mənfi PORO/PERM,
    PORO>1) bloklayıb — bura MÜSTƏQİL, İKİNCİ təsdiq üçündür (gate
    `model.validate()`-in nəticəsinə kor-koranə etibar etmir).
    XƏBƏRDARLIQ (bloklamır): interpolyasiya zamanı toplanmış `report.
    warnings` (seyrək data, QC düzəlişləri, fallback...) VƏ yüksək
    ekstrapolyasiya nisbəti (`max_extrapolated_fraction`-dan çox).

    Qayıdan `DiagnosticReport.has_errors is False` ⇒ "validated" statusu
    ALINA BİLƏR; `True` ⇒ ALINA BİLMƏZ."""
    gate = DiagnosticReport()
    for name, prop in model.property_maps.items():
        values = prop.values
        bad = ~np.isfinite(values)
        if np.any(bad):
            gate.error(f"{int(np.sum(bad))}/{values.size} hüceyrədə NaN/Inf dəyər.",
                      source=name, hint="QC/interpolyasiya konfiqurasiyasını yenidən yoxla.")
    for name, facies in model.facies_fields.items():
        if not np.all(np.isfinite(facies.codes)):
            gate.error("NaN/Inf fasiya kodu.", source=name)
    for message in report.warnings:
        gate.warning(message, source="interpolation")
    for name, q in (quality or {}).items():
        if q.extrapolated_fraction > max_extrapolated_fraction:
            gate.warning(
                f"hüceyrələrin {q.extrapolated_fraction * 100:.1f}%-i ekstrapolyasiyadır "
                f"(həddi {max_extrapolated_fraction * 100:.0f}%) — sərt data ilə zəif örtülüb.",
                source=name, hint="əlavə quyu/nöqtə əlavə et və ya `max_extrapolated_fraction` dəyiş.")
        if q.sample_count < 3:
            gate.warning(f"cəmi {q.sample_count} sərt data nöqtəsi — statistik cəhətdən zəif.",
                        source=name)
    return gate
