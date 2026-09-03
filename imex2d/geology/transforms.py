"""Xassə-fəzası ÇEVİRMƏLƏRİ — interpolyasiya HANSI fəzada aparılsın (B1).

Phase A-nın məkan/Kriging özəyi (`interpolation.py`, `variogram.py`,
`anisotropy.py`, `spatial_search.py`) dəyərin FİZİKİ MƏNASINI bilmir və
bilməməlidir. Rezervuar xassələrinin statistik təbiəti isə eyni DEYİL:

    PORO    kəsilməz, praktikada [0, 1] arasında, təxminən simmetrik
    PERMX   LOQ-NORMAL, ciddi müsbət, onlarla dəfə diapazon
    SW/NTG  HƏDLİ (bounded) [0, 1] — Gauss kriging hədləri poza bilər
    FACIES  KATEQORİK — "1.73 fasiya" mənasızdır

Bu modul həmin fərqi BİR yerdə, riyazi olaraq təyin edir. Hər çevirmə:

    forward(z)            → çevrilmiş (Gauss-a daha yaxın) fəza
    inverse(y)            → geri
    inverse_variance(...) → çevrilmiş fəzadakı VARİANSI orijinal fəzaya
                            köçürür (dəqiq düsturla, mümkün olmayanda
                            delta metodu ilə — HANSI olduğu bildirilir)

GERİ-ÇEVİRMƏ SEMANTİKASI (B1.3) — loq-normal üçün BİR cavab yoxdur::

    MEDIAN : K = exp(ŷ)                        şərti MEDİAN
    MEAN   : K = exp(ŷ + σ²/2)                 şərti ORTA (loq-normal)
    MEAN_OK: K = exp(ŷ + σ²/2 − μ)             adi kriging üçün DÜZƏLİŞLİ
                                               şərti orta (μ — Laqranj vuruğu)

Bu üçü FƏRQLİ kəmiyyətlərdir və heç biri "hər yerdə doğru" deyil, ona
görə seçim `PropertyStrategy`-də AÇIQ saxlanılır (bax `property_config.py`),
sükutla `exp(ŷ + σ²/2)` HƏR YERDƏ tətbiq EDİLMİR.

Niyə `MEAN_OK`: adi kriging (OK) Laqranj vuruğu `μ` ilə yansızlıq şərtini
məcbur edir; loq-fəzada yansız olan qiymət geri çevriləndə ORTA üçün
düzəliş `σ²_OK/2 − μ`-dir (Journel & Huijbregts, "Mining Geostatistics",
§7.4 — sadə kriginq üçün `μ = 0`, ona görə `MEAN`-ə çevrilir). Phase A-nın
`KrigingResult.lagrange` sahəsi məhz bunun üçün ötürülür.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional, Tuple

import numpy as np
from scipy.special import expit

from .gaussian_transform import NormalScoreTransform


class BackTransform(str, Enum):
    """Geri-çevirmənin STATİSTİK mənası (B1.3)."""

    MEDIAN = "median"      #: şərti median — `exp(ŷ)`
    MEAN = "mean"          #: şərti orta — `exp(ŷ + σ²/2)`
    MEAN_OK = "mean_ok"    #: adi kriginq üçün düzəlişli orta — `exp(ŷ + σ²/2 − μ)`


class VarianceKind(str, Enum):
    """`inverse_variance()`-in nə qaytardığı — səhv oxunmasın deyə AÇIQ."""

    EXACT = "exact"        #: qapalı düsturla dəqiq (loq-normal)
    DELTA = "delta"        #: birinci tərtib (delta metodu) YAXINLAŞMASI
    IDENTITY = "identity"  #: çevirmə yoxdur — varians olduğu kimi
    UNDEFINED = "undefined"  #: riyazi olaraq təyin edilməyib → NaN


class TransformError(ValueError):
    """Çevirmə tətbiq edilə bilməyən data — SƏSSİZ düzəliş EDİLMİR."""


class ValueTransform:
    """Baza sinif = EYNİLİK (identity) çevirməsi.

    Konkret sinif deyil, abstrakt deyil — qəsdən: `ValueTransform()` özü
    işlək eynilik çevirməsidir (Phase A-dakı davranış birəbir qorunur),
    alt-siniflər isə `forward`/`inverse`-i əvəz edir.

    `fit()` DATADAN asılı çevirmələr üçündür (normal-score, avtomatik
    hədlər). Eynilik/loq/logit üçün o, `self`-i qaytarır — yəni çevirmə
    məlumatdan ASILI DEYİL, deməli çarpaz-doğrulamada SIZMA (leakage)
    mənbəyi ola bilməz (bax `cross_validation.py`).
    """

    name = "identity"
    is_identity = True
    #: `fit()` təlim datasından statistika öyrənirmi (B3.1 sızma auditi)
    data_dependent = False

    def fit(self, values: np.ndarray) -> "ValueTransform":
        return self

    def validate(self, values: np.ndarray) -> None:
        """Çevirmənin tətbiq edilə bilməsini yoxlayır; ola bilməzsə atır."""
        values = np.asarray(values, float)
        if np.any(~np.isfinite(values)):
            raise TransformError(
                f"'{self.name}' çevirməsi NaN/sonsuz dəyər qəbul etmir.")

    def forward(self, values) -> np.ndarray:
        return np.asarray(values, float)

    def inverse(self, values, variance=None, lagrange=None,
                mode: BackTransform = BackTransform.MEDIAN) -> np.ndarray:
        return np.asarray(values, float)

    def inverse_variance(self, values, variance) -> Tuple[np.ndarray, VarianceKind]:
        return np.asarray(variance, float), VarianceKind.IDENTITY

    def describe(self) -> str:
        return self.name


# ── loq fəzası (PERMX/PERMY/PERMZ) ────────────────────────────────────
@dataclass
class LogTransform(ValueTransform):
    """`Y = ln(Z + offset)` — loq-normal xassələr (keçiricilik) üçün.

    `offset` (defolt 0) sıfıra-bərabər dəyərlərə icazə vermək üçün AÇIQ
    sürüşmədir; sıfırdan fərqlidirsə geri çevirmə də onu çıxır, yəni
    `inverse(forward(z)) == z`. `offset=0` olanda BÜTÜN dəyərlər ciddi
    müsbət olmalıdır — əks halda `TransformError` (səssizcə kiçik müsbət
    ədədə çevirmə YOXDUR, tapşırıq B1.2).
    """

    offset: float = 0.0
    name: str = "log"
    is_identity: bool = False
    data_dependent: bool = False

    def validate(self, values: np.ndarray) -> None:
        values = np.asarray(values, float)
        if np.any(~np.isfinite(values)):
            raise TransformError("Loq çevirməsi NaN/sonsuz dəyər qəbul etmir.")
        shifted = values + self.offset
        if np.any(shifted <= 0.0):
            bad = int(np.sum(shifted <= 0.0))
            raise TransformError(
                f"Loq çevirməsi üçün bütün dəyərlər müsbət olmalıdır "
                f"(offset={self.offset:g} ilə {bad} dəyər ≤ 0). Sıfır/mənfi "
                "keçiricilik SƏSSİZCƏ müsbət ədədə çevrilmir — ya `offset` verin, "
                "ya da məlumatı düzəldin (bax `data_quality.py` siyasətləri).")

    def forward(self, values) -> np.ndarray:
        values = np.asarray(values, float)
        self.validate(values)
        return np.log(values + self.offset)

    def inverse(self, values, variance=None, lagrange=None,
                mode: BackTransform = BackTransform.MEDIAN) -> np.ndarray:
        """`mode`-a görə median/orta/OK-düzəlişli orta (bax modul docstring-i)."""
        y = np.asarray(values, float)
        if mode == BackTransform.MEDIAN or variance is None:
            corrected = y
        else:
            sigma2 = np.clip(np.asarray(variance, float), 0.0, None)
            corrected = y + 0.5 * sigma2
            if mode == BackTransform.MEAN_OK:
                if lagrange is None:
                    raise TransformError(
                        "BackTransform.MEAN_OK üçün Laqranj vuruğu (`lagrange`) "
                        "lazımdır — `KrigingResult.lagrange` ötürün.")
                corrected = corrected - np.asarray(lagrange, float)
        return np.exp(corrected) - self.offset

    def inverse_variance(self, values, variance) -> Tuple[np.ndarray, VarianceKind]:
        """Loq-normalın DƏQİQ variansı::

            Var[K] = exp(2ŷ + σ²)·(exp(σ²) − 1)

        Yaxınlaşma DEYİL — `Y ~ N(ŷ, σ²)` fərziyyəsi altında qapalı düstur.
        """
        y = np.asarray(values, float)
        sigma2 = np.clip(np.asarray(variance, float), 0.0, None)
        with np.errstate(over="ignore"):
            result = np.exp(2.0 * y + sigma2) * np.expm1(sigma2)
        return result, VarianceKind.EXACT

    def describe(self) -> str:
        return f"ln(z + {self.offset:g})" if self.offset else "ln(z)"


# ── logit fəzası (SW/NTG və digər hədli xassələr) ─────────────────────
@dataclass
class LogitTransform(ValueTransform):
    """Hədli `[lower, upper]` xassələr üçün logit çevirməsi (B1.4/B1.5)::

        p = (z − lower) / (upper − lower)                    ∈ [0, 1]
        p̃ = eps + p·(1 − 2·eps)                              ∈ (0, 1)
        y = ln( p̃ / (1 − p̃) )

    Geri::

        p̃ = 1 / (1 + exp(−y))
        p  = (p̃ − eps) / (1 − 2·eps)   → [0,1]-ə kəsilir
        z  = lower + p·(upper − lower)

    NİYƏ `eps` VAR VƏ TƏSİRİ NƏDİR (tapşırıq B1.4: "do not introduce
    arbitrary epsilon without documenting its effect"):

    * `p = 0` və ya `1` olanda logit ±sonsuzdur — Kriging matrisi
      hesablana bilməz. `eps` bu iki nöqtəni SONLU dəyərə gətirir:
      `eps = 1e-4` üçün `y(0) = ln(1e-4/0.9999) ≈ −9.21`.
    * SIXMA GERİ AÇILIR (`inverse` `eps`-i çıxır), ona görə
      `inverse(forward(z)) == z` MAŞIN DƏQİQLİYİ ilə — DAXİL OLMAQLA
      dəqiq `lower`/`upper` dəyərləri. Yəni sərt datanın DƏQİQ honor
      edilməsi hədlərdə də POZULMUR.
    * Praktiki nəticə: `eps` yalnız hədlərdəki nöqtələrin logit fəzasında
      NƏ QƏDƏR uzağa düşdüyünü təyin edir. Kiçik `eps` → daha ekstremal
      `y` → variogram sillini şişirdir. Böyük `eps` → hədlər bir-birinə
      yaxınlaşır. `1e-4` bu ikisi arasında sənədləşdirilmiş kompromisdir.
    * Geri çevirmə RİYAZİ OLARAQ hədləri POZA BİLMİR (logistik funksiya
      (0,1)-dədir, sıxma isə [0,1]-ə kəsilir) — bu, "Gauss kriging
      fiziki cəhətdən mümkünsüz doyma verdi" probleminin KÖKDƏN həllidir.
    """

    lower: float = 0.0
    upper: float = 1.0
    eps: float = 1e-4
    name: str = "logit"
    is_identity: bool = False
    data_dependent: bool = False

    def __post_init__(self):
        if not np.isfinite(self.lower) or not np.isfinite(self.upper):
            raise TransformError("Logit hədləri sonlu olmalıdır.")
        if self.upper <= self.lower:
            raise TransformError(
                f"Logit üçün upper > lower olmalıdır ({self.upper} <= {self.lower}).")
        if not 0.0 < self.eps < 0.25:
            raise TransformError(f"eps (0, 0.25) aralığında olmalıdır, alındı: {self.eps}")

    @property
    def span(self) -> float:
        return float(self.upper - self.lower)

    def validate(self, values: np.ndarray) -> None:
        values = np.asarray(values, float)
        if np.any(~np.isfinite(values)):
            raise TransformError("Logit çevirməsi NaN/sonsuz dəyər qəbul etmir.")
        outside = (values < self.lower - 1e-12) | (values > self.upper + 1e-12)
        if np.any(outside):
            raise TransformError(
                f"Logit çevirməsi üçün bütün dəyərlər [{self.lower:g}, {self.upper:g}] "
                f"aralığında olmalıdır ({int(np.sum(outside))} dəyər kənardadır). "
                "Hədd pozan dəyər SƏSSİZCƏ kəsilmir — bax `data_quality.py`.")

    def forward(self, values) -> np.ndarray:
        values = np.asarray(values, float)
        self.validate(values)
        p = np.clip((values - self.lower) / self.span, 0.0, 1.0)
        squeezed = self.eps + p * (1.0 - 2.0 * self.eps)
        return np.log(squeezed / (1.0 - squeezed))

    def inverse(self, values, variance=None, lagrange=None,
                mode: BackTransform = BackTransform.MEDIAN) -> np.ndarray:
        """Logistik geri çevirmə.

        `mode` BURADA nəticəni DƏYİŞMİR: logit-normal paylanmanın ortası
        qapalı formada YOXDUR (elementar funksiyalarla ifadə edilmir), ona
        görə "orta" adı ilə YANLIŞ düstur tətbiq etmək əvəzinə MEDİAN
        qaytarılır. Bu, monoton çevirmə altında medianın DƏYİŞMƏZ qalması
        xassəsinə görə RİYAZİ OLARAQ DÜZGÜNDÜR (`median(g(Y)) =
        g(median(Y))`), sadəcə şərti orta deyil — və bu fərq burada AÇIQ
        yazılır (tapşırıq B1.3: "do not automatically use this formula
        everywhere").
        """
        y = np.asarray(values, float)
        # `expit` = 1/(1+e⁻ʸ), ƏDƏDİ DAYANIQLI (|y| ≫ 1-də daşma yoxdur)
        squeezed = expit(y)
        p = (squeezed - self.eps) / (1.0 - 2.0 * self.eps)
        return self.lower + np.clip(p, 0.0, 1.0) * self.span

    def inverse_variance(self, values, variance) -> Tuple[np.ndarray, VarianceKind]:
        """DELTA metodu (birinci tərtib yaxınlaşma)::

            dz/dy = span · (1 − 2ε) · p̃ · (1 − p̃)
            Var[z] ≈ (dz/dy)² · σ²_y

        Dəqiq deyil və belə də bildirilir (`VarianceKind.DELTA`): logit-
        normalın variansı qapalı formada yoxdur. Yaxınlaşma `σ_y` kiçik
        olanda yaxşı, hədlərə çox yaxın nöqtələrdə isə zəifdir.
        """
        y = np.asarray(values, float)
        sigma2 = np.clip(np.asarray(variance, float), 0.0, None)
        # `expit` = 1/(1+e⁻ʸ), ƏDƏDİ DAYANIQLI (|y| ≫ 1-də daşma yoxdur)
        squeezed = expit(y)
        derivative = self.span * (1.0 - 2.0 * self.eps) * squeezed * (1.0 - squeezed)
        return derivative ** 2 * sigma2, VarianceKind.DELTA

    def describe(self) -> str:
        return f"logit(z; [{self.lower:g}, {self.upper:g}], eps={self.eps:g})"


# ── normal-score (SGS-in tələb etdiyi fəza) ───────────────────────────
@dataclass
class NormalScoreValueTransform(ValueTransform):
    """`gaussian_transform.NormalScoreTransform`-un `ValueTransform`
    interfeysinə uyğunlaşdırılması.

    DATADAN ASILIDIR (`data_dependent=True`): cədvəl məhz fit edildiyi
    nümunədən qurulur. Buna görə çarpaz-doğrulamada HƏR qat üçün YENİDƏN
    fit edilməlidir — əks halda gizlədilmiş nöqtə çevirmə statistikasına
    sızardı (B3.1). `cross_validation.py` bunu MƏCBUR edir.
    """

    table: Optional[NormalScoreTransform] = None
    name: str = "normal_score"
    is_identity: bool = False
    data_dependent: bool = True

    def fit(self, values: np.ndarray) -> "NormalScoreValueTransform":
        return NormalScoreValueTransform(table=NormalScoreTransform.fit(values))

    def _table(self) -> NormalScoreTransform:
        if self.table is None:
            raise TransformError(
                "Normal-score çevirməsi işlədilməzdən ƏVVƏL `fit(values)` "
                "çağırılmalıdır (cədvəl datadan qurulur).")
        return self.table

    def validate(self, values: np.ndarray) -> None:
        values = np.asarray(values, float)
        if np.any(~np.isfinite(values)):
            raise TransformError("Normal-score çevirməsi NaN/sonsuz qəbul etmir.")

    def forward(self, values) -> np.ndarray:
        return self._table().forward(values)

    def inverse(self, values, variance=None, lagrange=None,
                mode: BackTransform = BackTransform.MEDIAN) -> np.ndarray:
        return self._table().inverse(values)

    def inverse_variance(self, values, variance) -> Tuple[np.ndarray, VarianceKind]:
        """Delta metodu — tərs cədvəlin ƏDƏDİ törəməsi ilə.

        Normal-score tərs çevirməsi PARÇALI XƏTTİ empirik funksiyadır,
        qapalı varians düsturu yoxdur. Törəmə mərkəzi fərqlə (Gauss
        fəzasında ±0.01 addım) hesablanır."""
        y = np.asarray(values, float)
        sigma2 = np.clip(np.asarray(variance, float), 0.0, None)
        table = self._table()
        if table.is_constant:
            return np.zeros_like(sigma2), VarianceKind.EXACT
        step = 0.01
        derivative = (table.inverse(y + step) - table.inverse(y - step)) / (2.0 * step)
        return derivative ** 2 * sigma2, VarianceKind.DELTA

    def describe(self) -> str:
        return "normal-score (empirik CDF → Gauss kvantili)"


#: Paylaşılan eynilik nüsxəsi (obyekt yaratmadan yoxlamaq üçün).
IDENTITY_TRANSFORM = ValueTransform()
LOG_TRANSFORM = LogTransform()


def apply_back_transform(transform: ValueTransform, estimate: np.ndarray,
                         variance: Optional[np.ndarray],
                         lagrange: Optional[np.ndarray],
                         mode: BackTransform
                         ) -> Tuple[np.ndarray, np.ndarray, VarianceKind]:
    """`(dəyər, varians, varians_növü)` — çevrilmiş fəzadan orijinala.

    Bir yerdə cəmlənib ki, "hansı mod hansı düsturla geri çevrildi"
    sualı bütün çağıranlar üçün EYNİ cavabı versin (B7 ardıcıllıq qaydası).
    """
    values = transform.inverse(estimate, variance=variance, lagrange=lagrange, mode=mode)
    if variance is None:
        return values, np.full(np.shape(estimate), np.nan), VarianceKind.UNDEFINED
    back_variance, kind = transform.inverse_variance(estimate, variance)
    return values, back_variance, kind
