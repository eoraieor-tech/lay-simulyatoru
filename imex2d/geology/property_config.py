"""XASSƏ STRATEGİYASI — hər rezervuar xassəsi üçün TƏK, AVTORİTAR konfiqurasiya (B1.7).

Bu modulun BÜTÜN məqsədi bir cümlədə: **POROSITY ≠ PERMEABILITY ≠
SATURATION ≠ NTG ≠ FACIES**, və bu fərq kodun HƏR YERİNƏ səpələnmiş
`if name == "PERMX"` yoxlamaları ilə DEYİL, BİR reyestrdə saxlanılır.

Hər `PropertyStrategy` bir xassə üçün bunları təyin edir:

    dəyişən növü          kəsilməz / loq-normal / hədli / kateqorik
    çevirmə               `transforms.ValueTransform` (loq, logit, eynilik…)
    geri-çevirmə mənası   median / orta / OK-düzəlişli orta (B1.3)
    fiziki hədlər         nə etibarlıdır (QC üçün) və nə tətbiq olunur
    hədd siyasəti         pozanı rədd et / kəs / toxunma
    interpolyasiya        kriging / IDW / indikator (kateqorik)
    variogram strategiyası model adı və ya "auto"
    anizotropluq          `AnisotropyParams` və ya `None` (avtomatik)
    sərt data siyasəti    `honor_hard_data` (bax `interpolation.py`)
    dublikat siyasəti     orta / median / ilk / son / ayrı saxla / xəta
    kənar-dəyər analizi   metod + hədd (SİLMİR — işarələyir, bax `data_quality.py`)
    kateqoriyalar         yalnız kateqorik xassələr üçün

Mövcud `property_types.classify_property()` (Phase 4.1) TOXUNULMUR —
bu modul onun üzərində qurulur: kateqorik/kəsilməz ilkin təsnifatı
oradan gəlir, statistik təfərrüat isə buradadır.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from .anisotropy import AnisotropyParams
from .property_types import PropertyType, classify_property
from .transforms import (IDENTITY_TRANSFORM, BackTransform, LogitTransform,
                         LogTransform, NormalScoreValueTransform, ValueTransform)
from .variogram import KNOWN_MODELS


class VariableType(str, Enum):
    """Xassənin STATİSTİK təbiəti — çevirmə/uncertainty seçimini idarə edir."""

    CONTINUOUS = "continuous"    #: kəsilməz, çevirməsiz (PORO, təzyiq)
    LOGNORMAL = "lognormal"      #: ciddi müsbət, çarpıq (PERMX/Y/Z)
    BOUNDED = "bounded"          #: [lo, hi] ilə məhdud (SW, NTG, VSH)
    CATEGORICAL = "categorical"  #: kod çoxluğu (FACIES, LITHOLOGY)


class BoundPolicy(str, Enum):
    """Fiziki hədləri pozan NƏTİCƏ ilə nə edilsin.

    DİQQƏT: bu, GİRİŞ datasına deyil, İNTERPOLYASİYA NƏTİCƏSİNƏ aiddir.
    Giriş datasının hədd pozması ayrıca məsələdir (bax `data_quality.py`)."""

    NONE = "none"      #: heç nə — nəticə olduğu kimi
    CLIP = "clip"      #: hədlərə kəsilir, KƏSİLƏN hüceyrələr SAYILIR və bildirilir
    FLAG = "flag"      #: kəsilmir, yalnız işarələnir (istifadəçi qərar verir)


class InterpolationKind(str, Enum):
    KRIGING = "kriging"
    IDW = "idw"
    NEAREST = "nearest"
    INDICATOR = "indicator"     #: kateqorik — indikator kriginq + ehtimallar


class DuplicatePolicy(str, Enum):
    """Eyni koordinatda birdən çox müşahidə (B4.2)."""

    MEAN = "mean"
    MEDIAN = "median"
    KEEP_FIRST = "keep_first"
    KEEP_LAST = "keep_last"
    KEEP_SEPARATE = "keep_separate"   #: birləşdirmə — solver jitter ilə həll edir
    RAISE = "raise"
    MAJORITY = "majority"             #: yalnız kateqorik


class OutlierMethod(str, Enum):
    """Kənar-dəyər DİAQNOSTİKASI (B4.5) — heç biri avtomatik SİLMİR."""

    NONE = "none"
    MAD = "mad"          #: robust z = 0.6745·|z − median| / MAD
    IQR = "iqr"          #: Tukey çəpərləri
    SPATIAL = "spatial"  #: yerli qonşuluğa nəzərən robust z (məkan kənarı)


class UncertaintyKind(str, Enum):
    """Qaytarılan variansın MƏNASI (B2.1) — qarışdırılmasın deyə."""

    KRIGING_VARIANCE = "kriging_variance"        #: orijinal fəzada, çevirmə yox
    BACK_TRANSFORMED = "back_transformed"        #: çevrilmiş fəzadan köçürülüb
    CATEGORICAL_ENTROPY = "categorical_entropy"  #: ehtimal paylanmasının entropiyası


class PropertyConfigError(ValueError):
    """Ziddiyyətli/etibarsız xassə strategiyası — səssiz düzəliş YOX."""


@dataclass(frozen=True)
class PropertyStrategy:
    """Bir xassənin TAM interpolyasiya reseptidir (B1.7).

    `frozen=True` — strategiya DƏYİŞMƏZDİR; variant qurmaq üçün
    `derive(...)` işlədilir. Bu, çarpaz-doğrulamada namizəd modelləri
    təhlükəsiz saxlamağa imkan verir (bir namizədin dəyişməsi digərinə
    təsir etmir).
    """

    name: str
    variable_type: VariableType = VariableType.CONTINUOUS
    #: interpolyasiyanın aparılacağı fəza
    transform: ValueTransform = IDENTITY_TRANSFORM
    back_transform: BackTransform = BackTransform.MEDIAN
    #: FİZİKİ etibarlılıq aralığı (QC üçün); `None` = məhdudiyyət yox
    physical_bounds: Tuple[Optional[float], Optional[float]] = (None, None)
    #: NƏTİCƏYƏ tətbiq olunan hədlər (adətən `physical_bounds` ilə eyni)
    output_bounds: Tuple[Optional[float], Optional[float]] = (None, None)
    bound_policy: BoundPolicy = BoundPolicy.CLIP
    interpolation: InterpolationKind = InterpolationKind.KRIGING
    variogram_model: str = "auto"
    anisotropy: Optional[AnisotropyParams] = None
    honor_hard_data: str = "auto"
    duplicate_policy: DuplicatePolicy = DuplicatePolicy.MEAN
    duplicate_tolerance: float = 1e-9
    outlier_method: OutlierMethod = OutlierMethod.MAD
    outlier_threshold: float = 3.5
    remove_outliers: bool = False          #: DEFOLT: SİLMİR, yalnız işarələyir
    max_neighbors: Optional[int] = 24
    min_neighbors: int = 1
    search_radius: Optional[float] = None
    sectors: int = 0
    #: kateqorik xassələr üçün icazə verilən kodlar (`None` = datadan)
    categories: Optional[Tuple[int, ...]] = None
    #: `interpolate_property`-nin köhnə `log_transform` bayrağı ilə körpü
    legacy_log_transform: bool = False
    notes: str = ""

    # ── doğrulama ─────────────────────────────────────────────────────
    def __post_init__(self):
        if self.variogram_model != "auto" and self.variogram_model not in KNOWN_MODELS:
            raise PropertyConfigError(
                f"{self.name}: naməlum variogram modeli {self.variogram_model!r}. "
                f"Dəstəklənən: {KNOWN_MODELS + ('auto',)}")
        if self.min_neighbors < 1:
            raise PropertyConfigError(f"{self.name}: min_neighbors ≥ 1 olmalıdır.")
        if self.max_neighbors is not None and self.max_neighbors < self.min_neighbors:
            raise PropertyConfigError(
                f"{self.name}: max_neighbors ({self.max_neighbors}) < min_neighbors "
                f"({self.min_neighbors}).")
        if self.honor_hard_data not in ("auto", "always", "never"):
            raise PropertyConfigError(
                f"{self.name}: honor_hard_data 'auto'/'always'/'never' olmalıdır.")
        lo, hi = self.physical_bounds
        if lo is not None and hi is not None and hi <= lo:
            raise PropertyConfigError(
                f"{self.name}: physical_bounds üst hədd alt həddən böyük olmalıdır.")
        if self.is_categorical and self.interpolation is not InterpolationKind.INDICATOR:
            raise PropertyConfigError(
                f"{self.name}: KATEQORİK xassə kəsilməz interpolyasiyadan KEÇƏ BİLMƏZ "
                f"({self.interpolation.value}) — 'FACIES=1.73' mənasızdır. "
                "InterpolationKind.INDICATOR işlədin (B1.6/GATE B4).")
        if (self.variable_type is VariableType.LOGNORMAL
                and self.transform.is_identity and not self.legacy_log_transform):
            raise PropertyConfigError(
                f"{self.name}: LOGNORMAL xassə üçün eynilik çevirməsi seçilib — "
                "xam keçiricilik üzərində xətti kriginq DEFOLT ola bilməz (GATE B2). "
                "Bunu bilərəkdən istəyirsinizsə `variable_type=CONTINUOUS` verin.")
        if (self.back_transform is not BackTransform.MEDIAN
                and self.transform.is_identity):
            raise PropertyConfigError(
                f"{self.name}: çevirmə yoxdursa geri-çevirmə modu MEDIAN olmalıdır "
                f"(alındı: {self.back_transform.value}).")

    # ── törəmə xassələr ───────────────────────────────────────────────
    @property
    def is_categorical(self) -> bool:
        return self.variable_type is VariableType.CATEGORICAL

    @property
    def is_bounded(self) -> bool:
        return self.variable_type is VariableType.BOUNDED

    @property
    def uncertainty_kind(self) -> UncertaintyKind:
        if self.is_categorical:
            return UncertaintyKind.CATEGORICAL_ENTROPY
        if self.transform.is_identity:
            return UncertaintyKind.KRIGING_VARIANCE
        return UncertaintyKind.BACK_TRANSFORMED

    def derive(self, **changes) -> "PropertyStrategy":
        """Dəyişdirilmiş NÜSXƏ (orijinal toxunulmaz qalır)."""
        return replace(self, **changes)

    # ── giriş datasının fiziki etibarlılığı ───────────────────────────
    def invalid_value_mask(self, values: np.ndarray) -> np.ndarray:
        """FİZİKİ olaraq etibarsız GİRİŞ dəyərlərinin maskası (B4.4).

        DİQQƏT — bu, kənar-dəyər (outlier) DEYİL: "fiziki cəhətdən
        mümkünsüz" (məs. mənfi məsaməlik) ilə "nadir amma mümkün"
        (məs. 5000 mD keçiricilik) FƏRQLİ problemlərdir və ayrı-ayrı
        işlənir (bax `data_quality.py`)."""
        values = np.asarray(values, float)
        invalid = ~np.isfinite(values)
        lo, hi = self.physical_bounds
        if lo is not None:
            invalid |= values < lo
        if hi is not None:
            invalid |= values > hi
        if self.variable_type is VariableType.LOGNORMAL:
            offset = getattr(self.transform, "offset", 0.0)
            with np.errstate(invalid="ignore"):
                invalid |= (values + offset) <= 0.0
        if self.is_categorical and self.categories is not None:
            allowed = np.asarray(self.categories, float)
            invalid |= ~np.isin(values, allowed)
        return invalid

    def fit_transform(self, values: np.ndarray) -> ValueTransform:
        """Datadan asılı çevirmələri (normal-score) FİT edir, digərlərini
        olduğu kimi qaytarır. Çarpaz-doğrulama bunu HƏR qat üçün ayrıca
        çağırır — sızmanın qarşısı (B3.1)."""
        return self.transform.fit(np.asarray(values, float))

    def apply_output_bounds(self, values: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """`(dəyər, düzəldilmiş_maska)` — `bound_policy`-yə görə.

        `CLIP` olanda kəsilir, `FLAG`/`NONE` olanda dəyər DƏYİŞMİR;
        hər üç halda hansı hüceyrələrin hədd pozduğu QAYTARILIR ki,
        səssiz kəsmə mümkün olmasın (tapşırıq: "do not silently clip")."""
        values = np.asarray(values, float)
        lo, hi = self.output_bounds
        violated = np.zeros(values.shape, dtype=bool)
        if lo is not None:
            violated |= values < lo
        if hi is not None:
            violated |= values > hi
        if self.bound_policy is BoundPolicy.CLIP and np.any(violated):
            values = values.copy()
            if lo is not None:
                values = np.maximum(values, lo)
            if hi is not None:
                values = np.minimum(values, hi)
        return values, violated

    def describe(self) -> str:
        lo, hi = self.physical_bounds
        bounds = ("—" if lo is None and hi is None
                  else f"[{'-∞' if lo is None else f'{lo:g}'}, "
                       f"{'+∞' if hi is None else f'{hi:g}'}]")
        return (f"{self.name}: {self.variable_type.value} | çevirmə "
                f"{self.transform.describe()} | geri {self.back_transform.value} | "
                f"hədlər {bounds} | {self.interpolation.value} | "
                f"variogram {self.variogram_model}")


# ── DEFOLT REYESTR ────────────────────────────────────────────────────
def _porosity(name: str = "PORO") -> PropertyStrategy:
    """Məsaməlik: kəsilməz, praktikada [0, 1], təxminən simmetrik.

    Çevirmə YOXDUR — məsaməlik loq-normal deyil və `[0,1]` hədlərinə
    real datada nadir hallarda yaxınlaşır; logit tətbiqi burada faydadan
    çox təhrif gətirər. Hədlər NƏTİCƏYƏ kəsilir və kəsmə SAYILIR."""
    return PropertyStrategy(
        name=name, variable_type=VariableType.CONTINUOUS,
        transform=IDENTITY_TRANSFORM, back_transform=BackTransform.MEDIAN,
        physical_bounds=(0.0, 1.0), output_bounds=(0.0, 1.0),
        bound_policy=BoundPolicy.CLIP, interpolation=InterpolationKind.KRIGING,
        notes="Kəsilməz; çevirməsiz kriginq; nəticə [0,1]-ə kəsilir və sayılır.")


def _permeability(name: str, offset: float = 0.0) -> PropertyStrategy:
    """Keçiricilik: LOQ-NORMAL. Xam kriginq DEFOLT DEYİL (GATE B2).

    Geri-çevirmə DEFOLT `MEDIAN`-dır: `exp(ŷ)` şərti mediandır və
    monoton çevirmə altında dəyişməzdir. Şərti ORTA istəyən çağıran
    `back_transform=BackTransform.MEAN` (və ya adi kriginq üçün riyazi
    olaraq daha düzgün `MEAN_OK`) seçir — bu, AÇIQ qərardır, sükutla
    tətbiq edilmir (tapşırıq B1.3)."""
    return PropertyStrategy(
        name=name, variable_type=VariableType.LOGNORMAL,
        transform=LogTransform(offset=offset), back_transform=BackTransform.MEDIAN,
        physical_bounds=(0.0, None), output_bounds=(1e-6, None),
        bound_policy=BoundPolicy.CLIP, interpolation=InterpolationKind.KRIGING,
        legacy_log_transform=True,
        notes="ln(K) fəzasında kriginq; geri-çevirmə defolt MEDİAN, "
              "orta üçün back_transform=MEAN/MEAN_OK.")


def _bounded_unit(name: str, notes: str) -> PropertyStrategy:
    """[0,1] hədli xassə (SW, NTG, VSH): logit fəzasında kriginq.

    Geri çevirmə RİYAZİ OLARAQ hədləri poza bilmir (logistik funksiya),
    ona görə "Gauss kriginq 1.07 doyma verdi" halı KÖKDƏN mümkün deyil."""
    return PropertyStrategy(
        name=name, variable_type=VariableType.BOUNDED,
        transform=LogitTransform(lower=0.0, upper=1.0, eps=1e-4),
        back_transform=BackTransform.MEDIAN,
        physical_bounds=(0.0, 1.0), output_bounds=(0.0, 1.0),
        bound_policy=BoundPolicy.CLIP, interpolation=InterpolationKind.KRIGING,
        notes=notes)


def _facies(name: str = "FACIES") -> PropertyStrategy:
    """Kateqorik: indikator kriginq → kateqoriya EHTİMALLARI (B1.6)."""
    return PropertyStrategy(
        name=name, variable_type=VariableType.CATEGORICAL,
        transform=IDENTITY_TRANSFORM, back_transform=BackTransform.MEDIAN,
        physical_bounds=(None, None), output_bounds=(None, None),
        bound_policy=BoundPolicy.NONE, interpolation=InterpolationKind.INDICATOR,
        duplicate_policy=DuplicatePolicy.MAJORITY,
        outlier_method=OutlierMethod.NONE,
        notes="İndikator kriginq; nəticə kateqoriya ehtimalları + ən ehtimallı kod.")


DEFAULT_STRATEGIES: Dict[str, PropertyStrategy] = {
    "PORO": _porosity("PORO"),
    "PERMX": _permeability("PERMX"),
    "PERMY": _permeability("PERMY"),
    "PERMZ": _permeability("PERMZ"),
    "SW": _bounded_unit("SW", "Su doyumluluğu: logit fəzasında kriginq, [0,1] təmin edilir."),
    "NTG": _bounded_unit("NTG", "Net-to-gross: logit fəzasında kriginq, [0,1] təmin edilir."),
    "VSH": _bounded_unit("VSH", "Gil həcmi: logit fəzasında kriginq, [0,1] təmin edilir."),
    "FACIES": _facies("FACIES"),
    "LITHOLOGY": _facies("LITHOLOGY"),
    "ROCKTYPE": _facies("ROCKTYPE"),
    "PRESSURE": PropertyStrategy(
        name="PRESSURE", variable_type=VariableType.CONTINUOUS,
        physical_bounds=(0.0, None), output_bounds=(0.0, None),
        notes="Təzyiq: kəsilməz, çevirməsiz; yalnız müsbətlik yoxlanılır."),
    "TOP": PropertyStrategy(
        name="TOP", variable_type=VariableType.CONTINUOUS,
        bound_policy=BoundPolicy.NONE, outlier_method=OutlierMethod.NONE,
        notes="Struktur səthi: hədsiz, kənar-dəyər analizi söndürülüb."),
    "BOTTOM": PropertyStrategy(
        name="BOTTOM", variable_type=VariableType.CONTINUOUS,
        bound_policy=BoundPolicy.NONE, outlier_method=OutlierMethod.NONE,
        notes="Struktur səthi: hədsiz, kənar-dəyər analizi söndürülüb."),
}


def resolve_strategy(name: str,
                     overrides: Optional[Dict[str, PropertyStrategy]] = None
                     ) -> PropertyStrategy:
    """`name` üçün strategiya — açıq `overrides`, sonra defolt reyestr.

    Naməlum ad üçün strategiya UYDURULMUR: `property_types.
    classify_property()` (Phase 4.1 reyestri) ilə kateqorik/kəsilməz
    təsnifatı alınır və NEYTRAL (çevirməsiz, hədsiz) kəsilməz strategiya
    qaytarılır — mövcud, geriyə-uyğun davranış. Kateqorik olduğu bilinən
    naməlum ad isə indikator yoluna göndərilir (kəsilməz kriginqə YOX).
    """
    key = name.upper()
    if overrides and key in {k.upper() for k in overrides}:
        for candidate_name, strategy in overrides.items():
            if candidate_name.upper() == key:
                return strategy
    if key in DEFAULT_STRATEGIES:
        return DEFAULT_STRATEGIES[key]
    if classify_property(name) is PropertyType.CATEGORICAL:
        return _facies(key)
    return PropertyStrategy(
        name=key, variable_type=VariableType.CONTINUOUS,
        bound_policy=BoundPolicy.NONE,
        notes="Reyestrdə olmayan ad — neytral kəsilməz strategiya (geriyə-uyğun).")


def normal_score_strategy(base: PropertyStrategy) -> PropertyStrategy:
    """`base`-in normal-score fəzasında işləyən variantı (SGS üçün).

    SGS Gauss fəzası TƏLƏB EDİR; bu funksiya eyni hədləri/siyasətləri
    saxlayaraq yalnız çevirməni əvəz edir, beləliklə SGS və deterministik
    kriginq EYNİ strategiya obyektindən törəyir (B7 ardıcıllıq qaydası)."""
    if base.is_categorical:
        raise PropertyConfigError(
            f"{base.name}: kateqorik xassə normal-score fəzasına çevrilə bilməz — "
            "kateqorik stoxastik modelləşdirmə üçün SIS (`facies.py`) işlədilir.")
    return base.derive(transform=NormalScoreValueTransform(),
                       back_transform=BackTransform.MEDIAN,
                       variable_type=VariableType.CONTINUOUS,
                       legacy_log_transform=False)


def strategy_table(names: Sequence[str] = tuple(DEFAULT_STRATEGIES)) -> str:
    """Reyestrin mətn cədvəli — hesabat/sənəd üçün."""
    lines = [f"{'Xassə':<10} {'Növ':<12} {'Çevirmə':<26} {'Geri':<8} "
             f"{'Hədlər':<12} {'Üsul':<10}"]
    lines.append("-" * len(lines[0]))
    for name in names:
        s = resolve_strategy(name)
        lo, hi = s.physical_bounds
        bounds = ("—" if lo is None and hi is None
                  else f"[{'-∞' if lo is None else f'{lo:g}'},"
                       f"{'∞' if hi is None else f'{hi:g}'}]")
        lines.append(f"{s.name:<10} {s.variable_type.value:<12} "
                     f"{s.transform.describe():<26} {s.back_transform.value:<8} "
                     f"{bounds:<12} {s.interpolation.value:<10}")
    return "\n".join(lines)


#: Hesabatda göstərilən standart xassə sırası.
REPORTED_PROPERTIES: Tuple[str, ...] = ("PORO", "PERMX", "PERMY", "PERMZ",
                                        "SW", "NTG", "FACIES")


__all__ = ["VariableType", "BoundPolicy", "InterpolationKind", "DuplicatePolicy",
           "OutlierMethod", "UncertaintyKind", "PropertyConfigError",
           "PropertyStrategy", "DEFAULT_STRATEGIES", "resolve_strategy",
           "normal_score_strategy", "strategy_table", "REPORTED_PROPERTIES"]
