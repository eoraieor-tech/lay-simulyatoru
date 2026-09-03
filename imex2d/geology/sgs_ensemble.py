"""SGS ANSAMBLI + DOĞRULAMA (B5.4 / B5.6 / B5.7).

`sgs.py` BİR realizasiya istehsal edir (normal-score → Gauss kriginq →
şərti nümunə → tərs çevirmə). Bu modul onun üzərində iki şey qurur:

    1. ANSAMBL       — `n` müstəqil realizasiya + ansambl statistikası
                       (orta, dispersiya, kvantillər).
    2. DOĞRULAMA     — realizasiyanın SGS-dən GÖZLƏNİLƏNİ verib-vermədiyi:
                       sərt datanı hörmət edirmi, marjinal paylanmanı
                       təkrarlayırmı, variogramı təkrarlayırmı.

TERMİNOLOGİYA (B5.7 — qəsdən dəqiq): `P10/P50/P90` ANSAMBL KVANTİLLƏRİDİR,
"etibar intervalı" DEYİL. Onlar modelin ÖZ fərziyyələri (variogram,
marjinal paylanma, kondisiyalaşdırma) daxilindəki dəyişkənliyi göstərir;
model səhvdirsə kvantillər də səhv olacaq. `SGSEnsemble` heç bir yerdə
"confidence interval" ifadəsini işlətmir.

TƏKRARLAMA YOXDUR: simulyasiyanın özü `sgs.simulate_sgs`-dədir, variogram
`variogram.py`-dədir, kriginq `interpolation.py`-dədir (B5.3 — İKİNCİ,
uyğunsuz variogram implementasiyası YARADILMIR).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from .distribution_analysis import DistributionSummary, summarize_distribution
from .sgs import PropertyRealization, PropertyVariogramParams, simulate_sgs
from .variogram import (MODEL_FUNCS, VariogramParameters, experimental_variogram,
                        fit_variogram)

#: `run_realizations_sgs`-in Phase 5 konvensiyası — dəyişdirilmir ki,
#: mövcud realizasiya identifikatorları/seed-ləri EYNİ qalsın.
SEED_STRIDE = 1000


@dataclass
class SGSEnsemble:
    """`n` realizasiyanın ansamblı + statistikası (B5.7)."""

    realizations: List[PropertyRealization]
    values: np.ndarray                  #: (n_real, m)
    hard_data_mask: np.ndarray          #: (m,) — bütün realizasiyalarda eyni
    base_seed: int
    warnings: List[str] = field(default_factory=list)

    @property
    def n_realizations(self) -> int:
        return int(self.values.shape[0])

    @property
    def n_cells(self) -> int:
        return int(self.values.shape[1])

    @property
    def mean(self) -> np.ndarray:
        """Ansambl ORTASI — hər hüceyrədə realizasiyalar üzrə orta.

        DİQQƏT: bu, KRİGİNQ qiyməti DEYİL (ona yaxınlaşır, amma
        realizasiya sayı sonlu olduğu üçün bərabər deyil) və HAMARDIR —
        heç bir realizasiyanın məkan davamlılığını daşımır. Simulyasiyanın
        məqsədi məhz hamarlıqdan qaçmaqdır, ona görə ansambl ortası
        MODEL kimi işlədilməməlidir."""
        return np.nanmean(self.values, axis=0)

    @property
    def variance(self) -> np.ndarray:
        """Ansambl DİSPERSİYASI — realizasiyalar arası dəyişkənlik.

        Bu, kriginq variansından FƏRQLİ kəmiyyətdir (B2.1 №7): kriginq
        variansı XƏTTİ qiymətin nəzəri xətasıdır, ansambl dispersiyası isə
        modelin verdiyi SİMULYASİYA qeyri-müəyyənliyidir."""
        return np.nanvar(self.values, axis=0, ddof=0)

    @property
    def std(self) -> np.ndarray:
        return np.sqrt(self.variance)

    def percentile(self, q: float) -> np.ndarray:
        """Ansambl kvantili (`q` faizlə). Realizasiya sayı azdırsa
        kvantil kobuddur — bu, xəbərdarlıqla bildirilir."""
        return np.nanpercentile(self.values, q, axis=0)

    @property
    def p10(self) -> np.ndarray:
        return self.percentile(10.0)

    @property
    def p50(self) -> np.ndarray:
        return self.percentile(50.0)

    @property
    def p90(self) -> np.ndarray:
        return self.percentile(90.0)

    def as_grids(self) -> Dict[str, np.ndarray]:
        return {"mean": self.mean, "variance": self.variance, "std": self.std,
                "p10": self.p10, "p50": self.p50, "p90": self.p90,
                "hard_data": self.hard_data_mask.astype(float)}

    def summary(self) -> str:
        lines = [f"SGS ansamblı: {self.n_realizations} realizasiya × "
                 f"{self.n_cells} hüceyrə (base_seed={self.base_seed})",
                 f"  ansambl ortası: [{np.nanmin(self.mean):.5g}, "
                 f"{np.nanmax(self.mean):.5g}]",
                 f"  ansambl std: orta {np.nanmean(self.std):.5g}",
                 f"  sərt data hüceyrəsi: {int(np.sum(self.hard_data_mask))}"]
        lines.extend(f"  ⚠ {w}" for w in self.warnings)
        return "\n".join(lines)


def simulate_sgs_ensemble(n_realizations: int, points, values, targets,
                          variogram: Optional[PropertyVariogramParams] = None,
                          base_seed: int = 0, **kwargs) -> SGSEnsemble:
    """`n_realizations` MÜSTƏQİL realizasiya (B5.4).

    Seed konvensiyası `sgs.run_realizations_sgs`-in EYNİSİDİR
    (`base_seed + i·1000`), ona görə mövcud çağıranlarla eyni nəticələr
    alınır. Zəmanətlər:

    * eyni `base_seed` + eyni konfiqurasiya → EYNİ ansambl (bit-bit);
    * fərqli seed → fərqli realizasiyalar;
    * hər realizasiya sərt datanı ÖZ konfiqurasiyasına görə hörmət edir.
    """
    if n_realizations < 1:
        raise ValueError(f"n_realizations ≥ 1 olmalıdır, alındı: {n_realizations}")
    realizations = [
        simulate_sgs(points, values, targets, variogram=variogram,
                     seed=base_seed + i * SEED_STRIDE, realization_id=i, **kwargs)
        for i in range(n_realizations)]

    stacked = np.vstack([r.values for r in realizations])
    warnings: List[str] = []
    if n_realizations < 10:
        warnings.append(
            f"Yalnız {n_realizations} realizasiya — P10/P90 kvantilləri KOBUDDUR "
            "(statistik etibarlı kvantil üçün adətən ≥ 50 realizasiya lazımdır).")
    for realization in realizations:
        for message in realization.warnings:
            if message not in warnings:
                warnings.append(f"[realizasiya {realization.realization_id}] {message}")

    return SGSEnsemble(realizations=realizations, values=stacked,
                       hard_data_mask=realizations[0].hard_data_mask.copy(),
                       base_seed=base_seed, warnings=warnings)


# ── doğrulama (B5.6) ──────────────────────────────────────────────────
@dataclass
class SGSValidationReport:
    """Bir realizasiyanın (və ya ansamblın) SGS keyfiyyət hesabatı."""

    n_cells: int
    n_hard_data: int
    #: sərt data hüceyrələrində maksimum MÜTLƏQ fərq (ideal 0)
    hard_data_max_error: float = float("nan")
    hard_data_honored: bool = False
    #: marjinal paylanma müqayisəsi
    data_summary: Optional[DistributionSummary] = None
    simulated_summary: Optional[DistributionSummary] = None
    mean_shift: float = float("nan")        #: (sim.orta − data.orta) / data.std
    std_ratio: float = float("nan")         #: sim.std / data.std (ideal 1)
    ks_statistic: float = float("nan")      #: Kolmoqorov-Smirnov D (ideal 0)
    #: variogram müqayisəsi (hədəf modelə nəzərən)
    target_variogram: Optional[VariogramParameters] = None
    realized_range: float = float("nan")
    range_ratio: float = float("nan")       #: realizə/hədəf (ideal 1)
    variogram_rmse: float = float("nan")    #: γ_sim vs γ_model, sillə normallanmış
    warnings: List[str] = field(default_factory=list)

    def as_dict(self) -> Dict[str, object]:
        return {"n_cells": self.n_cells, "n_hard_data": self.n_hard_data,
                "hard_data_max_error": self.hard_data_max_error,
                "hard_data_honored": self.hard_data_honored,
                "mean_shift": self.mean_shift, "std_ratio": self.std_ratio,
                "ks_statistic": self.ks_statistic,
                "realized_range": self.realized_range,
                "range_ratio": self.range_ratio,
                "variogram_rmse": self.variogram_rmse,
                "warnings": list(self.warnings)}

    def as_text(self) -> str:
        lines = [f"SGS doğrulaması ({self.n_cells} hüceyrə, "
                 f"{self.n_hard_data} sərt data)"]
        lines.append(f"  sərt data: maks xəta {self.hard_data_max_error:.3g} → "
                     f"{'HÖRMƏT EDİLİB' if self.hard_data_honored else 'POZULUB'}")
        if np.isfinite(self.ks_statistic):
            lines.append(f"  paylanma: KS D={self.ks_statistic:.4f} · "
                         f"orta sürüşməsi={self.mean_shift:+.3f}σ · "
                         f"std nisbəti={self.std_ratio:.3f}")
        if np.isfinite(self.range_ratio):
            lines.append(f"  variogram: realizə radius={self.realized_range:.4g} · "
                         f"nisbət={self.range_ratio:.3f} · "
                         f"γ-RMSE={self.variogram_rmse:.4f}")
        lines.extend(f"  ⚠ {w}" for w in self.warnings)
        return "\n".join(lines)


def _ks_statistic(a: np.ndarray, b: np.ndarray) -> float:
    """İki nümunə üçün Kolmoqorov-Smirnov `D` (maksimum CDF fərqi).

    `scipy.stats.ks_2samp` p-dəyər də hesablayır; burada YALNIZ `D`
    lazımdır (hipotez testi DEYİL — SGS-in marjinal paylanmanı nə qədər
    təkrarladığını göstərən ölçü)."""
    a = np.sort(np.asarray(a, float))
    b = np.sort(np.asarray(b, float))
    if a.size == 0 or b.size == 0:
        return float("nan")
    grid = np.concatenate([a, b])
    cdf_a = np.searchsorted(a, grid, side="right") / a.size
    cdf_b = np.searchsorted(b, grid, side="right") / b.size
    return float(np.max(np.abs(cdf_a - cdf_b)))


def _comparison_max_lag(points_xy: np.ndarray, target_range: float) -> float:
    """Realizə variogramını hədəf modellə müqayisə etmək üçün lag həddi.

    İKİ standart qayda kəsişir və KİÇİYİ götürülür:

    1. **Domenin yarısı** — deneysel variogram yalnız domen ölçüsünün
       təxminən yarısına qədər etibarlıdır; ondan uzaqda cüt sayı azalır
       və kənar effektləri üstünlük təşkil edir.
    2. **Hədəf radiusun 2.5 misli** — model onsuz da sillə çatıb;
       daha uzaq laglar müqayisəyə heç nə əlavə etmir.

    NİYƏ VACİBDİR (ölçülüb): 800 m domendə, hədəf radius 220 m olan
    realizasiya üçün laglar 591 m-ə (defolt 75%-lik kvantil) qədər
    aparılanda fit `735 m` radius verirdi — çünki deneysel `γ` uzaq
    laglarda sillin ÜSTÜNƏ qalxır (zonal/kənar effekti) və fitter bunu
    uzun radiusla izah edir. Etibarlı pəncərədə (≤ 550 m) eyni realizasiya
    `229 m` verir — yəni SGS məkan davamlılığını DÜZGÜN təkrarlayır, qüsur
    ÖLÇMƏ pəncərəsində idi.
    """
    extent = points_xy.max(axis=0) - points_xy.min(axis=0)
    half_domain = 0.5 * float(np.sqrt(np.sum(extent ** 2)))
    return float(max(min(2.5 * float(target_range), half_domain), 1e-9))


def validate_realization(realization: PropertyRealization, points, values, targets,
                         target_variogram: Optional[VariogramParameters] = None,
                         hard_data_tolerance: float = 1e-9,
                         distribution_ks_threshold: float = 0.25,
                         range_ratio_tolerance: float = 0.6) -> SGSValidationReport:
    """Bir realizasiyanı ÜÇ meyarla yoxlayır (B5.6).

    1. **SƏRT DATA** — `hard_data_mask` işarəli hüceyrələrdə simulyasiya
       dəyəri müşahidə ilə ÜST-ÜSTƏ düşməlidir (yalnız qlobal histoqram
       uyğunluğu KİFAYƏT DEYİL — B5.2).
    2. **MARJİNAL PAYLANMA** — simulyasiya edilmiş dəyərlərin paylanması
       sərt datanınkına yaxın olmalıdır (KS statistikası + moment
       müqayisəsi). Tam bərabərlik GÖZLƏNİLMİR: kondisiyalaşdırma və
       sonlu şəbəkə paylanmanı bir qədər dəyişir.
    3. **VARİOGRAM** — realizasiyadan fit edilən radius hədəf modelin
       radiusuna yaxın olmalıdır. SGS-in məkan davamlılığı təkrarlaması
       məhz budur.

    Hədlər (`*_threshold`) SƏNƏDLƏŞDİRİLMİŞ, konfiqurasiya edilə bilən
    tolerantlıqlardır — "keçdi/keçmədi" qərarı çağırana aiddir, bu
    funksiya ÖLÇÜR və xəbərdarlıq yazır.
    """
    values = np.asarray(values, float).ravel()
    targets = np.atleast_2d(np.asarray(targets, float))
    simulated = np.asarray(realization.values, float)
    hard_mask = np.asarray(realization.hard_data_mask, bool)
    report = SGSValidationReport(n_cells=int(simulated.size),
                                 n_hard_data=int(np.sum(hard_mask)))

    # ── 1. sərt data ──────────────────────────────────────────────────
    if np.any(hard_mask):
        from .hard_data import find_exact_matches
        points_array = np.atleast_2d(np.asarray(points, float))
        matches = find_exact_matches(points_array, targets, 1e-6)
        matched = matches >= 0
        error = np.abs(simulated[matched] - values[matches[matched]])
        report.hard_data_max_error = float(np.max(error)) if error.size else 0.0
        report.hard_data_honored = bool(report.hard_data_max_error <= hard_data_tolerance)
        if not report.hard_data_honored:
            report.warnings.append(
                f"Sərt data POZULUB: maksimum fərq {report.hard_data_max_error:.4g} > "
                f"tolerantlıq {hard_data_tolerance:.4g}.")
    else:
        report.hard_data_honored = True
        report.warnings.append(
            "Heç bir hədəf hüceyrəsi sərt data ilə üst-üstə düşmür — "
            "kondisiyalaşdırma yoxlanıla bilmədi.")

    # ── 2. marjinal paylanma ──────────────────────────────────────────
    finite = np.isfinite(simulated)
    if np.sum(finite) >= 3 and values.size >= 3:
        report.data_summary = summarize_distribution(values)
        report.simulated_summary = summarize_distribution(simulated[finite])
        data_std = max(report.data_summary.std, 1e-12)
        report.mean_shift = float(
            (report.simulated_summary.mean - report.data_summary.mean) / data_std)
        report.std_ratio = float(report.simulated_summary.std / data_std)
        report.ks_statistic = _ks_statistic(simulated[finite], values)
        if report.ks_statistic > distribution_ks_threshold:
            report.warnings.append(
                f"Marjinal paylanma fərqi böyükdür (KS D={report.ks_statistic:.3f} > "
                f"{distribution_ks_threshold:.2f}) — normal-score çevirməsi və ya "
                "kondisiyalaşdırma paylanmanı gözləniləndən çox dəyişib.")

    # ── 3. variogram ──────────────────────────────────────────────────
    # YALNIZ SİMULYASİYA EDİLMİŞ hüceyrələr: sərt-data hüceyrələri SGS-in
    # MƏHSULU DEYİL (onlar müşahidədir və ölçüsü/paylanması fərqlidir).
    # Onları da fitə qatmaq realizə radiusunu süni şəkildə uzadırdı —
    # ölçülüb: hədəf 220 m üçün qarışıq fit 780 m, yalnız simulyasiya
    # edilmiş hüceyrələrlə 211-229 m (tutulmuş qüsur).
    simulated_only = finite & ~hard_mask
    if target_variogram is not None and np.sum(simulated_only) >= 20:
        try:
            comparison_lag = _comparison_max_lag(
                targets[simulated_only][:, :2], target_variogram.range_)
            experimental = experimental_variogram(
                targets[simulated_only][:, :2], simulated[simulated_only],
                n_lags=12, max_lag=comparison_lag)
            # SİLL SABİT, YALNIZ RADİUS sərbəst — bax `_comparison_max_lag`
            # və `variogram.detect_anisotropy`-dəki EYNİ səbəb: nugget/sill/
            # radiusu birlikdə fit etmək zəif müəyyəndir. Burada realizə
            # sahəsinin dispersiyası ONSUZ DA məlumdur, ona görə onu sabit
            # saxlamaq həm düzgün, həm də sabitdir (ölçülüb: sərbəst fit
            # hədəf 400 m üçün 3542 m verirdi, sabit sill 534 m).
            realized = fit_variogram(
                experimental, model=target_variogram.model, fix_nugget=0.0,
                fix_sill=float(np.var(simulated[simulated_only])))
            report.realized_range = float(realized.range_)
            report.range_ratio = float(realized.range_ / max(target_variogram.range_, 1e-12))
            model_func = MODEL_FUNCS[target_variogram.model]
            expected = model_func(experimental.lags, target_variogram.nugget,
                                  target_variogram.sill, target_variogram.range_)
            filled = experimental.counts > 0
            total_sill = max(target_variogram.nugget + target_variogram.sill, 1e-12)
            report.variogram_rmse = float(np.sqrt(np.mean(
                ((experimental.gamma[filled] - expected[filled]) / total_sill) ** 2)))
            low, high = range_ratio_tolerance, 1.0 / range_ratio_tolerance
            if not low <= report.range_ratio <= high:
                report.warnings.append(
                    f"Realizə edilmiş variogram radiusu hədəfdən kənardır "
                    f"(nisbət {report.range_ratio:.2f}, gözlənilən "
                    f"[{low:.2f}, {high:.2f}]) — məkan davamlılığı tam təkrarlanmayıb. "
                    "Bu, sonlu şəbəkə/qonşuluq ölçüsünün məlum təsiri ola bilər.")
        except ValueError as exc:
            report.warnings.append(f"Realizasiyadan variogram fit alınmadı: {exc}")

    return report


def validate_ensemble(ensemble: SGSEnsemble, points, values, targets,
                      target_variogram: Optional[VariogramParameters] = None,
                      **kwargs) -> List[SGSValidationReport]:
    """Hər realizasiya üçün `validate_realization` — ansambl səviyyəsində."""
    return [validate_realization(realization, points, values, targets,
                                 target_variogram=target_variogram, **kwargs)
            for realization in ensemble.realizations]


def ensemble_statistics(ensemble: SGSEnsemble,
                        quantiles: Sequence[float] = (10.0, 50.0, 90.0)
                        ) -> Dict[str, np.ndarray]:
    """`mean`/`variance`/`std` + istənilən kvantillər (B5.7).

    Açarlar `p10`, `p50`, `p90` formasındadır. Bunlar ANSAMBL
    KVANTİLLƏRİDİR — "etibar intervalı" DEYİL (bax modul docstring-i)."""
    stats: Dict[str, np.ndarray] = {"mean": ensemble.mean,
                                    "variance": ensemble.variance,
                                    "std": ensemble.std}
    for q in quantiles:
        stats[f"p{int(round(q))}"] = ensemble.percentile(float(q))
    return stats
