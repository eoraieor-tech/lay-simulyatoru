"""Sequential Gaussian Simulation (SGS) — kəsilməz xassələr üçün (Phase 5).

`geology/facies.py`-dəki Sequential Indicator Simulation-un (Phase 4)
KƏSİLMƏZ analoqudur — EYNİ memarlıq (təsadüfi yol, yerli kriging, nəticəni
kondisioner çoxluğuna əlavə et), YALNIZ "indikator + kateqorik nümunələmə"
əvəzinə "normal-score çevirmə + Gauss kriging + normal nümunələmə"
işlədilir:

    sərt data
        → (istəyə görə) log-fəza (PERMX kimi çarpıq xassələr üçün)
        → normal-score çevirməsi (`gaussian_transform.py`) → Gauss fəzası
        → HƏR hədəf üçün: yerli kriging (`OrdinaryKriging.
          interpolate_with_variance`, Phase 5-də ƏLAVƏ EDİLİB) →
          N(mean_kriging, variance_kriging)-dən NÜMUNƏ (DEYİL ki, sadəcə
          kriging qiymətini yaz — bu, HAMARLANMIŞ DETERMİNİSTİK sahə
          verərdi, SGS-in ƏSAS QADAĞASI)
        → tərs normal-score çevirməsi → (istəyə görə) exp() → fiziki hədd

TƏKRARLAMA YOXDUR: variogram (`variogram.py`), kriging + anizotropluq
(`interpolation.OrdinaryKriging`), məkan axtarışı (`spatial_search.py`),
sərt-data həndəsəsi (`hard_data.py`) BİRBAŞA İSTİFADƏ OLUNUR — bu modul
YALNIZ Gauss çevirməsini, ardıcıl nümunələməni və hədd/diaqnostika
idarəetməsini əlavə edir.

ELMİ ÇƏKİNCƏ: bu modul PROQRAM CƏHƏTDƏN düzgün SGS-i tətbiq edir (Gauss
çevirməsi doğrudur, kriging variansı düzgün hesablanır, nümunə ŞƏRTİ
paylanmadan ÇƏKİLİR, sərt data hörmət olunur). Bu, YERİN gerçək xassə
paylanmasını əks etdirdiyi demək DEYİL (bax `facies.py`-dəki eyni
çəkincə) — YALNIZ kifayət qədər, təmsiledici quyu sıxlığı və DOĞRU
seçilmiş variogram/paylanma modeli ilə mümkündür.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from .distribution_analysis import log_transform_is_justified
from .gaussian_transform import NormalScoreTransform
from .hard_data import find_exact_matches
from .interpolation import OrdinaryKriging
from .spatial_search import IncrementalAnisotropicSearch
from .variogram import MODEL_SPHERICAL, AnisotropyParams, fit_variogram_from_data

DEFAULT_CORRECTION_WARN_THRESHOLD = 0.30

#: Fasiyaya-xas variogram/paylanma qurmaq üçün minimum sərt nöqtə sayı
#: (bax `simulate_sgs_facies_conditioned`, tapşırıq §5). Bundan az olan
#: fasiya ÜÇÜN AYRICA model UYDURULMUR — sənədləşdirilmiş FALLBACK-a keçilir.
DEFAULT_MIN_HARD_DATA_FOR_OWN_MODEL = 8


@dataclass
class PropertyVariogramParams:
    """Bir (xassə[, fasiya]) cütü üçün variogram + anizotropluq
    parametrləri — `facies.FaciesVariogramParams`-ın kəsilməz analoqu.
    Heç biri verilməyəndə (`range_=None`) sərt datadan (Gauss fəzasında)
    AVTOMATİK fit edilir."""
    model: str = MODEL_SPHERICAL
    nugget: float = 0.0
    sill: Optional[float] = None
    range_: Optional[float] = None
    range_v: Optional[float] = None
    azimuth_deg: Optional[float] = None
    range_minor: Optional[float] = None


@dataclass
class PropertyDiagnostics:
    """Ehtimal/hədd düzəlişlərinin AYRICA sayğacları (tapşırıq §6/§7:
    "do not silently clip... report number of corrections, rate")."""
    n_cells_simulated: int = 0
    bound_corrections: int = 0     # fiziki hədlərə görə kəsilən (clip) hüceyrə sayı
    nan_fallback_cells: int = 0    # yerli kriging NaN verdikdə qlobal N(0,1)-ə keçid

    def rate(self, count: int) -> float:
        return count / self.n_cells_simulated if self.n_cells_simulated else 0.0

    def summary_warnings(self, warn_threshold: float = DEFAULT_CORRECTION_WARN_THRESHOLD
                         ) -> List[str]:
        if self.n_cells_simulated == 0:
            return []
        messages = []
        for label, count in (("fiziki hədd kəsilməsi (clip)", self.bound_corrections),
                             ("NaN → qlobal N(0,1) geri-dönüşü", self.nan_fallback_cells)):
            if count == 0:
                continue
            rate = self.rate(count)
            messages.append(f"{label}: {count} hadisə ({rate * 100:.1f}% hüceyrə).")
            if rate > warn_threshold:
                messages.append(
                    f"GÜCLÜ XƏBƏRDARLIQ: \"{label}\" nisbəti ({rate * 100:.1f}%) "
                    f"həddi ({warn_threshold * 100:.0f}%) aşır — variogram/paylanma modeli "
                    "ƏDƏDİ CƏHƏTDƏN QEYRİ-SABİT ola bilər. Nəticəyə ehtiyatla yanaşın.")
        return messages


@dataclass
class PropertyRealization:
    """Bir SGS icrasının nəticəsi. `facies_reference` — hansı
    `FaciesField.name`-ə (varsa) əsaslanıb (tapşırıq §13: "reuse Phase 4
    realization concept where appropriate")."""
    realization_id: int
    seed: int
    values: np.ndarray
    hard_data_mask: np.ndarray
    facies_reference: Optional[str] = None
    metadata: Dict[str, object] = field(default_factory=dict)
    diagnostics: PropertyDiagnostics = field(default_factory=PropertyDiagnostics)
    warnings: List[str] = field(default_factory=list)


def _resolve_property_variogram(gaussian_points_xy: np.ndarray, gaussian_values: np.ndarray,
                                vp: Optional[PropertyVariogramParams], span: float
                                ) -> Tuple[dict, Optional[str]]:
    """`facies._resolve_facies_variogram`-ın kəsilməz analoqu — EYNİ
    fallback fəlsəfəsi (AÇIQ verilməyibsə fit et, mümkün olmasa AÇIQ
    xəbərdarlıqla domen/3 evristikasına keç).

    KRİTİK: `sill` HƏMİŞƏ konkret ədədə həll edilir, `None` saxlanılmır
    — bax `facies._resolve_facies_variogram`-dakı EYNİ qeyd (tutulmuş
    HƏQİQİ səhv: `sill=None` yerli axtarış alt-çoxluğundan asılı olaraq
    ÇAĞIRIŞDAN ÇAĞIRIŞA dəyişən, qeyri-stabil sill yaradırdı).
    """
    full_variance = max(float(np.var(gaussian_values)), 1e-12)

    if vp is not None and vp.range_ is not None:
        sill = vp.sill if vp.sill is not None else full_variance
        return dict(model=vp.model, nugget=vp.nugget, sill=sill, range_=vp.range_,
                    range_v=vp.range_v, azimuth_deg=vp.azimuth_deg,
                    range_minor=vp.range_minor), None
    model = vp.model if vp is not None else "auto"
    try:
        fit = fit_variogram_from_data(gaussian_points_xy, gaussian_values, model=model)
        return dict(model=fit.model, nugget=fit.nugget, sill=fit.sill, range_=fit.range_,
                    range_v=(vp.range_v if vp else None),
                    azimuth_deg=(vp.azimuth_deg if vp else None),
                    range_minor=(vp.range_minor if vp else None)), None
    except ValueError as exc:
        fallback_range = max(span / 3.0, 1e-6)
        warning = (f"Gauss-fəza variogramı fit alınmadı ({exc}) — ehtiyat evristika "
                  f"(domen/3 = {fallback_range:.3g}) işlədildi.")
        return dict(model=MODEL_SPHERICAL, nugget=0.0, sill=full_variance, range_=fallback_range,
                    range_v=(vp.range_v if vp else None),
                    azimuth_deg=(vp.azimuth_deg if vp else None),
                    range_minor=(vp.range_minor if vp else None)), warning


def _domain_span(points_xy: np.ndarray) -> float:
    if points_xy.shape[0] < 2:
        return 1.0
    lo, hi = points_xy.min(axis=0), points_xy.max(axis=0)
    return float(np.sqrt(np.sum((hi - lo) ** 2)))


def simulate_sgs(points, values, targets, variogram: Optional[PropertyVariogramParams] = None,
                 seed: int = 0, realization_id: int = 0,
                 search_radius: Optional[float] = None,
                 max_neighbors: Optional[int] = 24, min_neighbors: int = 1,
                 hard_data_tolerance: float = 1e-6,
                 bounds: Optional[Tuple[Optional[float], Optional[float]]] = None,
                 log_space: bool = False,
                 use_fast_search: bool = True, rebuild_interval: int = 64,
                 correction_warn_threshold: float = DEFAULT_CORRECTION_WARN_THRESHOLD,
                 facies_reference: Optional[str] = None) -> PropertyRealization:
    """Bir SGS realizasiyası (tək fasiya/paylanma, bax modul docstring-i
    üçün tam alqoritm).

    `bounds=(lo, hi)` — tərs çevirmədən sonra tətbiq olunan fiziki hədd
    (məs. PORO üçün `(0.0, 1.0)`-a yaxın, YAXUD layihənin öz
    `DEFAULT_RULES` hədləri) — kəsilən (clip) hər dəyər `diagnostics.
    bound_corrections`-da SAYILIR, SƏSSİZ DEYİL (tapşırıq §6).

    `log_space=True` — dəyərlər ƏVVƏLCƏ `log()`-a keçirilir (bütün
    dəyərlər müsbət OLMALIDIR, əks halda `ValueError`), SGS log-fəzada
    aparılır, sonda `exp()` ilə geri qaytarılır — nəticə HƏMİŞƏ müsbətdir
    (tapşırıq §7: "k > 0 must always be preserved").
    """
    points = np.atleast_2d(np.asarray(points, float))
    raw_values = np.asarray(values, float).ravel()
    targets = np.atleast_2d(np.asarray(targets, float))
    if points.shape[0] != raw_values.shape[0]:
        raise ValueError(
            f"points ({points.shape[0]}) və values ({raw_values.shape[0]}) uzunluğu uyğun gəlmir.")
    if points.shape[0] < 1:
        raise ValueError("Ən azı 1 sərt data nöqtəsi lazımdır.")
    if np.any(~np.isfinite(raw_values)):
        raise ValueError("Sərt data dəyərləri NaN/sonsuz ola bilməz.")
    if log_space and np.any(raw_values <= 0):
        raise ValueError("log_space=True yalnız BÜTÜN dəyərlər müsbət olanda işlədilə bilər.")

    warnings: List[str] = []
    n_targets = targets.shape[0]
    hard_index = find_exact_matches(points, targets, hard_data_tolerance)
    hard_mask = hard_index >= 0
    simulated = np.full(n_targets, np.nan)
    simulated[hard_mask] = raw_values[hard_index[hard_mask]]
    to_simulate = np.where(~hard_mask)[0]

    diag = PropertyDiagnostics()
    metadata: Dict[str, object] = {"log_space": log_space, "facies_reference": facies_reference}

    if to_simulate.size == 0:
        return PropertyRealization(realization_id, seed, simulated, hard_mask,
                                   facies_reference, metadata, diag, warnings)

    working_values = np.log(raw_values) if log_space else raw_values.copy()
    if np.ptp(working_values) < 1e-12:
        # SABİT XASSƏ: Gauss çevirməsi degenerativdir, SIMULYASİYAYA
        # EHTİYAC YOXDUR — sabit dəyər UYDURULMUŞ dəyişkənlik OLMADAN yayılır.
        constant = float(raw_values[0])
        simulated[to_simulate] = constant
        warnings.append(
            "Sərt data SABİTDİR (dəyişkənlik yoxdur) — SGS keçildi, sabit dəyər tətbiq edildi.")
        return PropertyRealization(realization_id, seed, simulated, hard_mask,
                                   facies_reference, metadata, diag, warnings)

    transform = NormalScoreTransform.fit(working_values)
    gaussian_values = transform.forward(working_values)

    points3 = points if points.shape[1] >= 3 else np.column_stack(
        [points, np.zeros(points.shape[0])])
    targets3 = targets if targets.shape[1] >= 3 else np.column_stack(
        [targets, np.zeros(targets.shape[0])])
    points_xy = points[:, :2] if points.shape[1] >= 2 else points

    span = _domain_span(points_xy)
    kwargs, warn = _resolve_property_variogram(points_xy, gaussian_values, variogram, span)
    if warn:
        warnings.append(warn)
    metadata["variogram"] = dict(kwargs)

    range_h = kwargs["range_"]
    range_v = kwargs["range_v"] if kwargs.get("range_v") is not None else range_h
    range_minor = kwargs["range_minor"] if kwargs.get("range_minor") is not None else range_h
    azimuth = kwargs["azimuth_deg"] if kwargs.get("azimuth_deg") is not None else 0.0
    aniso = AnisotropyParams(azimuth_deg=azimuth, range_major=range_h,
                             range_minor=range_minor, range_vertical=range_v)

    if use_fast_search:
        search = IncrementalAnisotropicSearch(points3, anisotropy=aniso,
                                              rebuild_interval=rebuild_interval)
        estimator = OrdinaryKriging(search_radius=None, max_neighbors=None,
                                    min_neighbors=1, **kwargs)
    else:
        estimator = OrdinaryKriging(search_radius=search_radius, max_neighbors=max_neighbors,
                                    min_neighbors=min_neighbors, **kwargs)

    rng = np.random.default_rng(seed)
    path = rng.permutation(to_simulate)

    # ƏVVƏLCƏDƏN AYRILMIŞ bufer: kondisioner çoxluğu hər addımda bir
    # sətir/dəyər böyüyür. `np.vstack`/`np.append` hər dəfə BÜTÜN massivi
    # YENİDƏN köçürürdü (O(n²) məcmu köçürmə — profillə təsdiqlənib: 150×150
    # şəbəkədə ~%8 ümumi vaxt və ölçü ilə superxətti artır). Son ölçü
    # ƏVVƏLCƏDƏN MƏLUMDUR (sərt data + hədəf sayı), ona görə tək dəfə
    # ayrılıb doldurulur. Nəticə EYNİDİR: eyni sətir/dəyər sırası, eyni RNG
    # çağırış ardıcıllığı — YALNIZ yaddaş idarəetməsi dəyişib.
    max_conditioning = points3.shape[0] + to_simulate.size
    sim_points = np.empty((max_conditioning, points3.shape[1]), dtype=points3.dtype)
    sim_points[:points3.shape[0]] = points3
    sim_gaussian = np.empty(max_conditioning, dtype=gaussian_values.dtype)
    sim_gaussian[:gaussian_values.size] = gaussian_values
    n_conditioning = points3.shape[0]

    for idx in path:
        target_point = targets3[idx:idx + 1]
        diag.n_cells_simulated += 1
        active_points = sim_points[:n_conditioning]
        active_gaussian = sim_gaussian[:n_conditioning]

        if use_fast_search:
            neighbor_idx = search.query(target_point, search_radius=search_radius,
                                        max_neighbors=max_neighbors, min_neighbors=min_neighbors)
            if neighbor_idx.size == 0:
                mean_g, var_g = 0.0, 1.0   # Gauss fəzasında marjinal N(0,1) — bax modul docstring-i
                diag.nan_fallback_cells += 1
            else:
                est, var = estimator.interpolate_with_variance(
                    active_points[neighbor_idx], active_gaussian[neighbor_idx], target_point)
                mean_g, var_g = float(est[0]), float(var[0])
        else:
            est, var = estimator.interpolate_with_variance(
                active_points, active_gaussian, target_point)
            mean_g, var_g = float(est[0]), float(var[0])
            if not np.isfinite(mean_g) or not np.isfinite(var_g):
                mean_g, var_g = 0.0, 1.0
                diag.nan_fallback_cells += 1

        sample_g = rng.normal(mean_g, np.sqrt(max(var_g, 0.0)))
        simulated_g_value = sample_g
        sim_points[n_conditioning] = target_point[0]
        sim_gaussian[n_conditioning] = simulated_g_value
        n_conditioning += 1
        if use_fast_search:
            search.add_point(target_point)
        simulated[idx] = simulated_g_value   # müvəqqəti: tərs çevirmə AŞAĞIDA toplu aparılır

    # tərs çevirmə: Gauss -> (log-fəza) -> orijinal vahid
    simulated_mask = ~hard_mask
    inv = transform.inverse(simulated[simulated_mask])
    if log_space:
        inv = np.exp(inv)
    simulated[simulated_mask] = inv

    if bounds is not None:
        lo, hi = bounds
        below = (simulated_mask & (simulated < lo) if lo is not None
                else np.zeros_like(hard_mask))
        above = (simulated_mask & (simulated > hi) if hi is not None
                else np.zeros_like(hard_mask))
        n_corrected = int(np.sum(below)) + int(np.sum(above))
        if lo is not None:
            simulated[below] = lo
        if hi is not None:
            simulated[above] = hi
        diag.bound_corrections = n_corrected

    warnings.extend(diag.summary_warnings(correction_warn_threshold))
    return PropertyRealization(realization_id, seed, simulated, hard_mask, facies_reference,
                               metadata, diag, warnings)


def run_realizations_sgs(n_realizations: int, points, values, targets,
                         variogram: Optional[PropertyVariogramParams] = None,
                         seed: int = 0, **kwargs) -> List[PropertyRealization]:
    """`n_realizations` MÜSTƏQİL SGS realizasiyası — `facies.run_
    realizations`-la EYNİ seed konvensiyası (`seed + i*1000`)."""
    return [simulate_sgs(points, values, targets, variogram=variogram, seed=seed + i * 1000,
                         realization_id=i, **kwargs)
            for i in range(n_realizations)]


@dataclass
class FaciesPropertyConfig:
    """Bir (fasiya) üçün xassə konfiqurasiyası — `facies × property`
    matrisinin bir xanası (tapşırıq §5)."""
    variogram: Optional[PropertyVariogramParams] = None
    log_space: Optional[bool] = None   # None -> `distribution_analysis`-dan AVTOMATİK qərar
    bounds: Optional[Tuple[Optional[float], Optional[float]]] = None


def _global_fallback_structural_model(points: np.ndarray, values: np.ndarray
                                      ) -> Tuple[Optional["PropertyVariogramParams"],
                                                 Optional[NormalScoreTransform]]:
    """Fasiyaya-xas model qurmaq mümkün olmayanda (çox az/sıfır ÖZ sərt
    data) istifadə olunan PAYLAŞILAN (bütün fasiyalar üzrə pooled) STRUKTUR
    modeli — variogram forması/aralığı VƏ marjinal paylanma forması.

    KRİTİK AYRIM: bu, YALNIZ struktur/forma məlumatıdır. Heç bir fasiyanın
    FAKTİKİ sərt-data DƏYƏRİ bu yolla başqa fasiyaya "hard-conditioning"
    kimi keçmir — çağıran (`simulate_sgs_facies_conditioned`) hər zaman
    yalnız fasiyanın ÖZ nöqtələrini `simulate_sgs`-ə sərt-data kimi verir,
    bu funksiya isə yalnız (əgər fasiyaya-xas fit mümkün deyilsə) variogram
    aralığı/nisbətini və ya (heç bir öz nöqtə olmayanda) marjinal paylanma
    formasını borc verir — standart geostatistik təcrübə (bax: seyrək
    alt-populyasiya üçün qlobal variogramdan struktur borc almaq).
    """
    if points.shape[0] < 2 or np.ptp(values) < 1e-12:
        return None, None
    points_xy = points[:, :2] if points.shape[1] >= 2 else points
    transform = NormalScoreTransform.fit(values)
    gaussian_values = transform.forward(values)
    span = _domain_span(points_xy)
    kwargs, _ = _resolve_property_variogram(points_xy, gaussian_values, None, span)
    vp = PropertyVariogramParams(model=kwargs["model"], nugget=kwargs["nugget"],
                                 sill=kwargs["sill"], range_=kwargs["range_"],
                                 range_v=kwargs.get("range_v"), azimuth_deg=kwargs.get("azimuth_deg"),
                                 range_minor=kwargs.get("range_minor"))
    return vp, transform


def simulate_sgs_facies_conditioned(
        points, values, facies_at_points, targets, facies_at_targets,
        facies_configs: Optional[Dict[int, FaciesPropertyConfig]] = None,
        seed: int = 0, realization_id: int = 0,
        min_hard_data_for_own_model: int = DEFAULT_MIN_HARD_DATA_FOR_OWN_MODEL,
        facies_reference: Optional[str] = None, **kwargs) -> PropertyRealization:
    """Fasiya-şərtli SGS: hər fasiya ÖZ paylanması/variogramı ilə, YALNIZ
    öz sərt datası VƏ öz hədəf hüceyrələri üzərində simulyasiya olunur.

    QAT-İ QAYDA (elmi düzgünlük): bir fasiyanın hard-conditioning-i HEÇ
    VAXT başqa fasiyanın sərt-data DƏYƏRLƏRİNDƏN istifadə etmir — bu,
    fasiyalar arası "sızma" (cross-facies contamination) yaradardı (məs.
    Facies A-nın hədəf hüceyrəsi Facies B-nin PORO=100-130 kimi dəyərini
    "dəqiq" sərt-data kimi ala bilərdi, halbuki bu nöqtə fiziki olaraq
    Facies A-ya aid deyil).

    Kifayət qədər sərt data (< `min_hard_data_for_own_model`) OLMAYAN
    fasiya üçün fasiyaya-xas variogram/paylanma AYRICA qurula bilmir.
    Bu halda:

    - `0 < n_own < min_hard_data_for_own_model` — hard-conditioning YENƏ
      DƏ yalnız bu fasiyanın öz nöqtələri ilə aparılır; YALNIZ variogram
      STRUKTURU (aralıq/nisbət) bütün fasiyalar üzrə PAYLAŞILAN (qlobal)
      modeldən borc alınır (`_global_fallback_structural_model`) —
      standart geostatistik təcrübə (seyrək alt-populyasiya üçün qlobal
      struktur borcu), FAKTİKİ DƏYƏR borcu DEYİL.
    - `n_own == 0` — heç bir kondisiyalaşdırma mümkün deyil (kondisiyalaşdırmaq
      üçün başqa fasiyanın dəyərini götürmək QADAĞANDIR); ƏVƏZİNƏ bütün
      fasiyalar üzrə paylaşılan qlobal marjinal paylanmadan ŞƏRTSİZ
      (məkan kondisiyası olmadan) nümunə çəkilir.

    Hər iki fallback AÇIQ xəbərdarlıqla bildirilir (tapşırıq §5: "do not
    fabricate a variogram... report that fallback").
    """
    points = np.atleast_2d(np.asarray(points, float))
    values = np.asarray(values, float).ravel()
    facies_at_points = np.asarray(facies_at_points, int).ravel()
    targets = np.atleast_2d(np.asarray(targets, float))
    facies_at_targets = np.asarray(facies_at_targets, int).ravel()
    if points.shape[0] != values.shape[0] or points.shape[0] != facies_at_points.shape[0]:
        raise ValueError("points/values/facies_at_points uzunluqları uyğun gəlmir.")
    if targets.shape[0] != facies_at_targets.shape[0]:
        raise ValueError("targets/facies_at_targets uzunluqları uyğun gəlmir.")

    facies_configs = facies_configs or {}
    result_values = np.full(targets.shape[0], np.nan)
    result_hard_mask = np.zeros(targets.shape[0], dtype=bool)
    all_warnings: List[str] = []
    combined_diag = PropertyDiagnostics()
    per_facies_metadata: Dict[int, dict] = {}

    # PAYLAŞILAN struktur modeli (variogram forması + marjinal paylanma) —
    # bax `_global_fallback_structural_model` docstring-i: YALNIZ struktur
    # borc alınır, heç bir fasiyanın DƏYƏRİ başqasına hard-data kimi keçmir.
    global_variogram, global_transform = _global_fallback_structural_model(points, values)

    facies_list = sorted(set(np.unique(facies_at_targets).tolist()))
    for facies in facies_list:
        target_mask = facies_at_targets == facies
        if not np.any(target_mask):
            continue
        point_mask = facies_at_points == facies
        facies_points, facies_values = points[point_mask], values[point_mask]
        config = facies_configs.get(facies, FaciesPropertyConfig())
        n_own = facies_points.shape[0]

        if n_own == 0:
            # QAYDA: başqa fasiyanın sərt-data DƏYƏRLƏRİ bu fasiya üçün
            # HEÇ VAXT hard-conditioning kimi istifadə OLUNMUR — heç bir öz
            # nöqtə olmadıqda kondisiyalı simulyasiya sadəcə MÜMKÜN DEYİL.
            # Yeganə elmi cəhətdən düzgün seçim: bütün fasiyalar üzrə
            # PAYLAŞILAN qlobal marjinal paylanmadan ŞƏRTSİZ (heç bir məkan
            # kondisiyası olmadan) nümunə çəkmək.
            n_cells = int(np.sum(target_mask))
            all_warnings.append(
                f"Fasiya {facies}: HEÇ bir öz sərt nöqtə yoxdur — fasiyaya-xas VƏ YA "
                "kondisiyalı simulyasiya mümkün deyil. Başqa fasiyanın sərt-data "
                "dəyərləri hard-conditioning kimi İSTİFADƏ OLUNMUR (elmi qayda). "
                "ƏVƏZİNƏ bütün fasiyalar üzrə paylaşılan qlobal marjinal "
                "paylanmadan ŞƏRTSİZ nümunə çəkildi — nəticə heç bir konkret "
                "quyu ilə MƏKANCA əlaqələndirilməyib.")
            rng = np.random.default_rng(seed + int(facies) * 7919)
            if global_transform is None:
                fallback_value = float(values[0]) if values.size else 0.0
                facies_values_out = np.full(n_cells, fallback_value)
            else:
                facies_values_out = global_transform.inverse(rng.normal(0.0, 1.0, size=n_cells))
            n_corrected = 0
            if config.bounds is not None:
                lo, hi = config.bounds
                if lo is not None:
                    below = facies_values_out < lo
                    n_corrected += int(np.sum(below))
                    facies_values_out[below] = lo
                if hi is not None:
                    above = facies_values_out > hi
                    n_corrected += int(np.sum(above))
                    facies_values_out[above] = hi
            result_values[target_mask] = facies_values_out
            combined_diag.n_cells_simulated += n_cells
            combined_diag.bound_corrections += n_corrected
            per_facies_metadata[int(facies)] = {
                "used_cross_facies_fallback": False,
                "used_unconditional_global_fallback": True,
                "n_hard_points": 0,
            }
            continue

        used_fallback = False
        variogram_override = config.variogram
        if n_own < min_hard_data_for_own_model:
            used_fallback = True
            all_warnings.append(
                f"Fasiya {facies}: yalnız {n_own} sərt nöqtə var "
                f"(minimum {min_hard_data_for_own_model} tələb olunur) — fasiyaya-xas "
                "variogram etibarlı fit edilə bilmədi. ƏVƏZİNƏ bütün fasiyalar üzrə "
                "PAYLAŞILAN (qlobal) variogram STRUKTURU (aralıq/nisbət) istifadə "
                f"olundu — AMMA hard-conditioning YALNIZ bu fasiyanın öz {n_own} "
                "nöqtəsi ilə aparılır, başqa fasiyanın sərt-data DƏYƏRLƏRİ heç vaxt "
                "bu fasiyaya keçmir.")
            if variogram_override is None:
                variogram_override = global_variogram

        log_space = config.log_space
        if log_space is None:
            log_space = log_transform_is_justified(facies_values)

        sub = simulate_sgs(facies_points, facies_values, targets[target_mask],
                           variogram=variogram_override, seed=seed, realization_id=realization_id,
                           bounds=config.bounds, log_space=log_space,
                           facies_reference=facies_reference, **kwargs)
        result_values[target_mask] = sub.values
        result_hard_mask[target_mask] = sub.hard_data_mask
        combined_diag.n_cells_simulated += sub.diagnostics.n_cells_simulated
        combined_diag.bound_corrections += sub.diagnostics.bound_corrections
        combined_diag.nan_fallback_cells += sub.diagnostics.nan_fallback_cells
        per_facies_metadata[int(facies)] = {
            **sub.metadata, "used_cross_facies_fallback": used_fallback,
            "n_hard_points": int(n_own),
            "log_space": log_space,
        }
        all_warnings.extend(f"[fasiya {facies}] {w}" for w in sub.warnings)

    return PropertyRealization(
        realization_id=realization_id, seed=seed, values=result_values,
        hard_data_mask=result_hard_mask, facies_reference=facies_reference,
        metadata={"per_facies": per_facies_metadata}, diagnostics=combined_diag,
        warnings=all_warnings)
