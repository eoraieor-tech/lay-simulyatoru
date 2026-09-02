"""Mərkəzləşdirilmiş fiziki yoxlama qatı (Phase 1).

Hər funksiya `ValidationResult(errors, warnings)` qaytarır:

    errors   — FİZİKİ CƏHƏTDƏN QEYRİ-MÜMKÜN dəyərlər (mənfi keçiricilik,
               1-dən böyük doyumluluq, sıfır/mənfi hüceyrə həcmi...).
               Çağıran bunları ADƏTƏN `ValueError` kimi yuxarı ötürməlidir.
    warnings — "QEYRİ-ADİ, AMMA FİZİKİ CƏHƏTDƏN MÜMKÜN" dəyərlər (məs.
               20 D keçiricilik, 50000 cP lözlük). RƏDD EDİLMİR — real
               yataq məlumatı tez-tez "qeyri-adi" olur, ona görə burada
               yalnız İSTİFADƏÇİYƏ BİLDİRİLİR, hesablama DAYANDIRILMIR.

Bu modul MÖVCUD `validate()` metodlarını (bax `PVTTable.validate`,
`SaturationTable.validate`, `WellDataset.validate`, `RockProperties.
validate`) TƏKRARLAMIR və DƏYİŞDİRMİR (geriyə uyğunluq) — yalnız onların
İNDİYƏ QƏDƏR VERMƏDİYİ iki şeyi əlavə edir: (1) aydın NaN/inf mesajı,
(2) "qeyri-adi diapazon" xəbərdarlığı. Hansı obyektlərin bu funksiyaları
faktiki çağırdığı üçün bax hər faylın "Wired" qeydini (`pvt.py`,
`geometry.py`, `properties.py`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np


@dataclass
class ValidationResult:
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def extend(self, other: "ValidationResult") -> "ValidationResult":
        self.errors.extend(other.errors)
        self.warnings.extend(other.warnings)
        return self


def _finite_issue(values: np.ndarray, label: str) -> Optional[str]:
    values = np.asarray(values, float)
    n_nan = int(np.sum(np.isnan(values)))
    n_inf = int(np.sum(np.isinf(values)))
    if n_nan or n_inf:
        parts = []
        if n_nan:
            parts.append(f"{n_nan} NaN")
        if n_inf:
            parts.append(f"{n_inf} sonsuz")
        return f"{label}: etibarsız dəyər tapıldı ({', '.join(parts)})."
    return None


# ── petrofiziki xassələr ─────────────────────────────────────────────────
def validate_porosity(values, label: str = "PORO") -> ValidationResult:
    """0 <= phi <= 1 sərt tələbdir. phi > 0.40 (çox nadir hətta boşluqlu
    qumda) və ya phi < 0.02 (praktik olaraq qeyri-kollektor) XƏBƏRDARLIQ —
    RƏDD EDİLMİR, bəzi karbonat/fraktura modellərində bu diapazondan
    kənar dəyər qanuni ola bilər."""
    result = ValidationResult()
    values = np.asarray(values, float)
    finite = _finite_issue(values, label)
    if finite:
        result.errors.append(finite)
        return result
    if np.any(values < 0.0):
        result.errors.append(f"{label}: mənfi məsaməlilik fiziki cəhətdən qeyri-mümkündür.")
    if np.any(values >= 1.0):
        result.errors.append(f"{label}: məsaməlilik >= 1.0 fiziki cəhətdən qeyri-mümkündür.")
    if np.any((values > 0.40) & (values < 1.0)):
        result.warnings.append(
            f"{label}: {int(np.sum(values > 0.40))} hüceyrədə məsaməlilik 0.40-dan "
            "yüksəkdir — qeyri-adidir, boşluqlu qum/karst istisna olmaqla.")
    if np.any((values > 0.0) & (values < 0.02)):
        result.warnings.append(
            f"{label}: {int(np.sum((values > 0.0) & (values < 0.02)))} hüceyrədə "
            "məsaməlilik 0.02-dən aşağıdır — praktik olaraq qeyri-kollektor ola bilər.")
    return result


def validate_saturation(values, label: str = "Sw") -> ValidationResult:
    result = ValidationResult()
    values = np.asarray(values, float)
    finite = _finite_issue(values, label)
    if finite:
        result.errors.append(finite)
        return result
    if np.any(values < -1e-9) or np.any(values > 1.0 + 1e-9):
        result.errors.append(f"{label}: doyumluluq [0, 1] intervalından kənardadır.")
    return result


def validate_permeability(values, label: str = "PERMX", unit: str = "mD") -> ValidationResult:
    """k > 0 sərt tələbdir (log-keçiricilik interpolyasiyası da bunu
    tələb edir — bax `geology/interpolation.py`). k > 20000 mD (vuqlu
    karbonat) və ya k < 0.01 mD (sıx/qeyri-ənənəvi) XƏBƏRDARLIQ."""
    result = ValidationResult()
    values = np.asarray(values, float)
    finite = _finite_issue(values, label)
    if finite:
        result.errors.append(finite)
        return result
    if np.any(values <= 0.0):
        result.errors.append(
            f"{label}: sıfır və ya mənfi keçiricilik fiziki cəhətdən qeyri-mümkündür.")
    positive = values[values > 0.0]
    if unit == "mD" and positive.size:
        if np.any(positive > 20000.0):
            result.warnings.append(
                f"{label}: {int(np.sum(positive > 20000.0))} hüceyrədə keçiricilik "
                "20000 mD-dən yüksəkdir — qeyri-adidir (vuqlu karbonat istisna olmaqla).")
        if np.any(positive < 0.01):
            result.warnings.append(
                f"{label}: {int(np.sum(positive < 0.01))} hüceyrədə keçiricilik "
                "0.01 mD-dən aşağıdır — sıx/qeyri-ənənəvi kollektor ola bilər.")
    return result


def validate_viscosity(value, label: str = "lözlük", unit: str = "cP") -> ValidationResult:
    result = ValidationResult()
    values = np.atleast_1d(np.asarray(value, float))
    finite = _finite_issue(values, label)
    if finite:
        result.errors.append(finite)
        return result
    if np.any(values <= 0.0):
        result.errors.append(f"{label}: sıfır və ya mənfi lözlük fiziki cəhətdən qeyri-mümkündür.")
    positive = values[values > 0.0]
    if unit == "cP" and positive.size:
        if np.any(positive > 50000.0):
            result.warnings.append(
                f"{label}: 50000 cP-dən yüksək — ağır/bitum nefti istisna olmaqla qeyri-adidir.")
        if np.any(positive < 0.05):
            result.warnings.append(f"{label}: 0.05 cP-dən aşağı — qeyri-adidir (qazdan da azdır).")
    return result


def validate_density(value, label: str = "sıxlıq", unit: str = "kg/m3") -> ValidationResult:
    result = ValidationResult()
    values = np.atleast_1d(np.asarray(value, float))
    finite = _finite_issue(values, label)
    if finite:
        result.errors.append(finite)
        return result
    if np.any(values <= 0.0):
        result.errors.append(f"{label}: sıfır və ya mənfi sıxlıq fiziki cəhətdən qeyri-mümkündür.")
    positive = values[values > 0.0]
    if unit == "kg/m3" and positive.size and (np.any(positive < 200.0) or np.any(positive > 1500.0)):
        result.warnings.append(f"{label}: [200, 1500] kg/m³ diapazonundan kənardadır — qeyri-adidir.")
    return result


def validate_compressibility(value, label: str = "sıxılma", unit: str = "1/bar") -> ValidationResult:
    result = ValidationResult()
    values = np.atleast_1d(np.asarray(value, float))
    finite = _finite_issue(values, label)
    if finite:
        result.errors.append(finite)
        return result
    if np.any(values <= 0.0):
        result.errors.append(
            f"{label}: sıfır və ya mənfi sıxılma fiziki cəhətdən qeyri-mümkündür "
            "(stabil sistemdə sıxılma müsbətdir).")
    positive = values[values > 0.0]
    if unit == "1/bar" and positive.size and (np.any(positive < 1e-6) or np.any(positive > 1e-2)):
        result.warnings.append(
            f"{label}: [1e-6, 1e-2] 1/bar diapazonundan kənardadır — qeyri-adidir.")
    return result


def validate_pressure(values, label: str = "təzyiq", unit: str = "bar",
                      depth_m: Optional[float] = None,
                      fracture_gradient_bar_per_m: float = 0.160) -> ValidationResult:
    """p > 0 (mütləq təzyiq) sərt tələbdir. `depth_m` verilibsə,
    hidrostatik+çatlama qradiyenti təxminini keçən təzyiq XƏBƏRDARLIQ
    doğurur (`fracture_gradient_bar_per_m` defolt dəyəri
    `reservoir_model.FRACTURE_GRADIENT` ilə eynidir — anormal-yüksək-
    təzyiqli yataqlarda bu qanuni şəkildə aşıla bilər, ona görə XƏTA yox,
    XƏBƏRDARLIQdır)."""
    result = ValidationResult()
    values = np.atleast_1d(np.asarray(values, float))
    finite = _finite_issue(values, label)
    if finite:
        result.errors.append(finite)
        return result
    if np.any(values <= 0.0):
        result.errors.append(f"{label}: sıfır və ya mənfi mütləq təzyiq fiziki cəhətdən qeyri-mümkündür.")
    if depth_m is not None and unit == "bar" and depth_m > 0:
        limit = fracture_gradient_bar_per_m * depth_m
        if np.any(values > limit):
            result.warnings.append(
                f"{label}: {depth_m:.0f} m dərinlikdə təxmini çatlama həddini "
                f"({limit:.0f} bar) aşır — anormal-yüksək-təzyiqli yataq deyilsə yoxlanılmalıdır.")
    return result


# ── həndəsə ────────────────────────────────────────────────────────────
def validate_grid_dimensions(nx: int, ny: int, nz: int,
                             dx: float, dy: float) -> ValidationResult:
    result = ValidationResult()
    for name, value in (("NX", nx), ("NY", ny), ("NZ", nz)):
        if value < 1:
            result.errors.append(f"{name}: müsbət tam ədəd olmalıdır (alındı: {value}).")
    for name, value in (("DX", dx), ("DY", dy)):
        if not np.isfinite(value) or value <= 0.0:
            result.errors.append(f"{name}: müsbət olmalıdır (alındı: {value}).")
    return result


def validate_thickness(dz, label: str = "DZ") -> ValidationResult:
    result = ValidationResult()
    values = np.atleast_1d(np.asarray(dz, float))
    finite = _finite_issue(values, label)
    if finite:
        result.errors.append(finite)
        return result
    if np.any(values <= 0.0):
        result.errors.append(f"{label}: hər təbəqənin qalınlığı müsbət olmalıdır.")
    return result


def validate_cell_volumes(volumes, label: str = "hüceyrə həcmi") -> ValidationResult:
    """Sıfır/mənfi həcm — dejenerativ hüceyrə, sərt xəta. Çox kiçik/çox
    böyük həcm bəzən VAHİD QARIŞIQLIĞININ (məs. ft yerinə m) əlaməti ola
    bilər — buna görə XƏBƏRDARLIQ (rədd yox, yoxlama tövsiyəsi)."""
    result = ValidationResult()
    values = np.atleast_1d(np.asarray(volumes, float))
    finite = _finite_issue(values, label)
    if finite:
        result.errors.append(finite)
        return result
    if np.any(values <= 0.0):
        n_bad = int(np.sum(values <= 0.0))
        result.errors.append(
            f"{label}: {n_bad} hüceyrədə sıfır/mənfi həcm (dejenerativ həndəsə).")
    positive = values[values > 0.0]
    if positive.size:
        if np.any(positive < 1e-3):
            result.warnings.append(
                f"{label}: {int(np.sum(positive < 1e-3))} hüceyrə çox kiçikdir "
                "(< 0.001 m³) — vahid qarışıqlığı (məs. ft əvəzinə m) yoxlanılmalıdır.")
        if np.any(positive > 1e9):
            result.warnings.append(
                f"{label}: {int(np.sum(positive > 1e9))} hüceyrə çox böyükdür "
                "(> 1e9 m³) — vahid qarışıqlığı yoxlanılmalıdır.")
    return result


# ── quyu ──────────────────────────────────────────────────────────────
def validate_well_rate(value: float, label: str = "debit") -> ValidationResult:
    """Mühərrik konvensiyası: istifadəçi RATE hədəfini HƏMİŞƏ müsbət
    böyüklük kimi verir (bax `standard_well.py:_signed_rate_target`) —
    mənfi dəyər burada XƏTADIR. Sıfır debit (bağlı quyu niyyəti) mümkün
    ola bilər, ona görə YALNIZ XƏBƏRDARLIQdır."""
    result = ValidationResult()
    if not np.isfinite(value):
        result.errors.append(f"{label}: etibarsız dəyər (NaN/sonsuz).")
        return result
    if value < 0.0:
        result.errors.append(
            f"{label}: mənfi debit qəbul edilmir — istiqamət `WellType`/`Phase` "
            "ilə müəyyənləşir, böyüklük həmişə müsbətdir.")
    elif value == 0.0:
        result.warnings.append(f"{label}: debit sıfırdır — bağlı quyu nəzərdə tutulursa qanunidir.")
    elif value > 100000.0:
        result.warnings.append(f"{label}: {value:g} m³/gün — tək quyu üçün qeyri-adi yüksəkdir.")
    return result


def validate_query_range(values, table_min: float, table_max: float, label: str = "sorğu",
                         severe_factor: float = 0.5) -> ValidationResult:
    """Sorğu dəyərini cədvəl diapazonu ilə müqayisə edir — DƏRƏCƏLİ:

    * diapazon daxilində — problem yoxdur.
    * diapazondan bir az kənar (< `severe_factor`·diapazon eni qədər) —
      XƏBƏRDARLIQ: "açıq icazə verilən" yüngül ekstrapolyasiya (cədvəl
      sərhədə kəsir, nəticə DƏYİŞMİR, amma bilinməlidir).
    * diapazondan ÇOX kənar (>= `severe_factor`·diapazon eni) — SƏRT
      XƏTA: bu, adətən VAHİD QARIŞIQLIĞININ əlamətidir (məs. psi dəyəri
      bar cədvəlinə sorğulanıb) — "qeyri-müəyyən ekstrapolyasiya" kimi
      İCAZƏSİZ sayılır.
    """
    result = ValidationResult()
    values = np.atleast_1d(np.asarray(values, float))
    span = max(table_max - table_min, 1e-9)
    severe_lo = table_min - severe_factor * span
    severe_hi = table_max + severe_factor * span

    severe = int(np.sum((values < severe_lo) | (values > severe_hi)))
    mild = int(np.sum(((values < table_min) & (values >= severe_lo))
                      | ((values > table_max) & (values <= severe_hi))))
    if severe:
        result.errors.append(
            f"{label}: {severe} dəyər cədvəl diapazonundan ([{table_min:g}, "
            f"{table_max:g}]) HƏDDİNDƏN ARTIQ kənardadır — bu, VAHİD "
            "QARIŞIQLIĞI əlaməti ola bilər (məs. psi/bar). Qeyri-müəyyən "
            "ekstrapolyasiya İCAZƏ VERİLMİR.")
    if mild:
        result.warnings.append(
            f"{label}: {mild} dəyər cədvəl diapazonundan ([{table_min:g}, "
            f"{table_max:g}]) bir qədər kənardadır — sərhədə kəsilib "
            "(açıq icazə verilən yüngül ekstrapolyasiya).")
    return result


# ── fasiya nisbətləri (Phase 4) ──────────────────────────────────────────
def validate_facies_proportions(proportions: Dict[int, float],
                                label: str = "fasiya nisbətləri",
                                tol: float = 1e-6) -> ValidationResult:
    """`{fasiya_kodu: nisbət}` lüğəti — cəmi 1.0 olmalıdır (tolerantlıq
    daxilində), hər nisbət [0, 1]-dədir. Bunlar İSTİFADƏÇİNİN AÇIQ
    seçdiyi hədəf nisbətlərdir (məlumatdan çıxarılmır) — SƏSSİZCƏ
    dəyişdirilmir, yalnız fiziki/riyazi cəhətdən qeyri-mümkün olanda
    rədd edilir."""
    result = ValidationResult()
    if not proportions:
        result.errors.append(f"{label}: boşdur.")
        return result
    values = np.array(list(proportions.values()), float)
    finite = _finite_issue(values, label)
    if finite:
        result.errors.append(finite)
        return result
    if np.any(values < -1e-9):
        result.errors.append(f"{label}: mənfi nisbət qəbul edilmir.")
    if np.any(values > 1.0 + 1e-9):
        result.errors.append(f"{label}: 1.0-dan böyük tək nisbət qəbul edilmir.")
    total = float(values.sum())
    if abs(total - 1.0) > tol:
        result.errors.append(
            f"{label}: cəm {total:.6f}, 1.0 olmalıdır (tolerantlıq {tol:g}). "
            "Nisbətlər avtomatik normallaşdırılmır — özünüz düzəldin.")
    return result


def compare_observed_vs_requested_proportions(
        observed: Dict[int, float], requested: Dict[int, float],
        label: str = "fasiya nisbətləri", warn_threshold: float = 0.15) -> List[str]:
    """Müşahidə olunan (sərt data-dan hesablanmış) və tələb olunan
    nisbətləri müqayisə edir — FƏRQ OLANDA YALNIZ XƏBƏRDARLIQ verir,
    heç bir dəyəri DƏYİŞDİRMİR (istifadəçinin açıq seçimi toxunulmazdır)."""
    warnings: List[str] = []
    for code in sorted(set(observed) | set(requested)):
        o = observed.get(code, 0.0)
        r = requested.get(code, 0.0)
        if abs(o - r) > warn_threshold:
            warnings.append(
                f"{label}: fasiya {code} üçün quyu məlumatında müşahidə olunan "
                f"nisbət ({o:.3f}) tələb olunan nisbətdən ({r:.3f}) çox fərqlənir "
                f"(fərq > {warn_threshold:g}).")
    return warnings


# ── PVT/SCAL — mövcud `validate()`-in ÜSTÜNƏ əlavə diaqnostika ─────────
def check_extrapolation_range(query_values, table_min: float, table_max: float,
                              label: str = "sorğu") -> List[str]:
    """Cədvəlin laboratoriya diapazonundan KƏNAR sorğunu bildirir.

    Mövcud `BlackOilPVTProvider`/`SaturationTable.interpolate_*` sərhədə
    KƏSİR (clamp) — bu, Nyutonu dayandırmaq üçün DOĞRU mühəndislik
    seçimidir (bax `A6_PLAN.md`), amma İSTİFADƏÇİYƏ "nəticə ekstrapolyasiya
    olunub" bildirilmir. Bu funksiya HƏLƏ heç bir avtomatik axına
    bağlanmayıb (bax GEOSTATISTICS.md-yə bənzər qeyd, hesabat sonunda) —
    çağıran (məs. `ReservoirModelBuilder`) İSTƏYƏ görə çağıra bilər.
    """
    values = np.atleast_1d(np.asarray(query_values, float))
    below = int(np.sum(values < table_min))
    above = int(np.sum(values > table_max))
    warnings = []
    if below:
        warnings.append(
            f"{label}: {below} dəyər cədvəlin aşağı həddindən ({table_min:g}) kiçikdir "
            "— sərhədə kəsilib (ekstrapolyasiya EDİLMİR, sabit dəyər saxlanılır).")
    if above:
        warnings.append(
            f"{label}: {above} dəyər cədvəlin yuxarı həddindən ({table_max:g}) böyükdür "
            "— sərhədə kəsilib (ekstrapolyasiya EDİLMİR, sabit dəyər saxlanılır).")
    return warnings
