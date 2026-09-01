"""Fasiya modelləşdirməsi — Sequential Indicator Simulation (SIS), Phase 4.

Niyə lazımdır: `geology_service.py`-dəki mövcud iş axını (Kriging/IDW)
KATEQORİK dəyəri (fasiya kodu 0, 1, 2...) ADİ RƏQƏMSƏL kimi
interpolyasiya edir — nəticədə iki fasiya arasında "1.4" kimi mənasız
aralıq dəyər çıxa bilər. Bu modul əvəzinə HƏR fasiya üçün indikator
funksiyası

    I_k(x) = 1  əgər fasiya(x) = k
    I_k(x) = 0  əks halda

qurub, HƏR birini AYRICA Kriging edir (indikator kriging), sonra
ehtimallardan TƏSADÜFİ NÜMUNƏ götürür — Deutsch & Journel (GSLIB)
"Sequential Indicator Simulation" alqoritminin standart formasıdır.

TƏKRARLAMA YOXDUR: variogram riyaziyyatı (`geology/variogram.py`) və
kriging həlli (`geology/interpolation.OrdinaryKriging`) BİRBAŞA İSTİFADƏ
OLUNUR — bu modul YALNIZ (1) hər fasiya üçün indikator çevirməsini,
(2) ardıcıl simulyasiya yolunu/nümunələməni, (3) ehtimal normallaşdırma
diaqnostikasını əlavə edir.

ELMİ ÇƏKİNCƏ (tapşırığın "Phase 26/15" qaydası): bu modul PROQRAM
CƏHƏTDƏN düzgün SIS-i tətbiq edir (indikator çevirməsi doğrudur, kriging
ehtimalları normallaşdırılır, sərt data hörmət olunur, təkrarlana
bilər). Bu, YERİN ÖZÜNÜN gerçək fasiya paylanmasını əks etdirdiyi demək
DEYİL — real geoloji "doğruluq" YALNIZ kifayət qədər, təmsiledici quyu
sıxlığı və DOĞRU seçilmiş variogram modeli ilə mümkündür. Testlər
proqram düzgünlüyünü yoxlayır (bax `tests/test_facies.py`) — "geoloji
reallıq"ı YOX, bunu YOXLAMAQ MÜMKÜN DEYİL sintetik testlə.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..domain.validation import (ValidationResult, compare_observed_vs_requested_proportions,
                                 validate_facies_proportions)
from .hard_data import find_exact_matches
from .interpolation import OrdinaryKriging
from .spatial_search import IncrementalAnisotropicSearch
from .variogram import MODEL_SPHERICAL, AnisotropyParams, fit_variogram_from_data

#: Ehtimal cəmi bundan çox kənarlaşarsa "nəzərəçarpan" sayılır və
#: xəbərdarlıqda qeyd olunur. Kriging çəkiləri mənfi/1-dən çox ola
#: bildiyi üçün kiçik kənarlaşma NORMALDIR — yalnız BÖYÜK kənarlaşma
#: diaqnostik dəyərə malikdir.
LARGE_RENORMALIZATION_THRESHOLD = 0.05

#: `FaciesDiagnostics.summary_warnings()`-in defolt həddi (§7) — bu
#: nisbətdən çox hüceyrə düzəliş/geri-dönüş tələb edərsə GÜCLÜ
#: xəbərdarlıq verilir (variogram/qonşuluq modeli ədədi cəhətdən
#: qeyri-sabit ola bilər).
DEFAULT_CORRECTION_WARN_THRESHOLD = 0.30


@dataclass
class FaciesDiagnostics:
    """Ehtimal düzəlişlərinin AYRICA sayğacları (tapşırıq §7) — heç biri
    digərinə YIĞILMIR, çünki hərəsi FƏRQLİ ədədi problemi göstərir."""
    n_cells_simulated: int = 0
    negative_probability_events: int = 0     # neçə HÜCEYRƏDƏ ən azı bir fasiya mənfi ehtimal verdi
    excess_probability_events: int = 0       # neçə hüceyrədə kliplənmədən ƏVVƏL cəm 1-dən çox idi
    nan_fallback_cells: int = 0              # neçə hüceyrə NaN kriging nəticəsi ilə qlobal nisbətə keçdi
    zero_sum_fallback_cells: int = 0         # neçə hüceyrə sıfır-cəm ilə qlobal nisbətə keçdi

    @property
    def global_proportion_fallback_cells(self) -> int:
        """§7-də ayrıca tələb olunan "qlobal nisbət geri-dönüşü" ÜMUMİ sayı."""
        return self.nan_fallback_cells + self.zero_sum_fallback_cells

    def rate(self, count: int) -> float:
        return count / self.n_cells_simulated if self.n_cells_simulated else 0.0

    def summary_warnings(self, warn_threshold: float = DEFAULT_CORRECTION_WARN_THRESHOLD
                         ) -> List[str]:
        """Hər kateqoriya üçün say+faiz, HƏDDİ AŞANDA əlavə GÜCLÜ mesaj.
        Heç bir düzəliş baş verməyibsə boş siyahı (səssiz-uğurlu)."""
        if self.n_cells_simulated == 0:
            return []
        messages = []
        checks = [
            ("mənfi kriging ehtimalı", self.negative_probability_events),
            ("kliplənmədən əvvəl ehtimal cəmi > 1", self.excess_probability_events),
            ("NaN ehtimal → qlobal/regional nisbətə keçid", self.nan_fallback_cells),
            ("sıfır-cəm ehtimal → qlobal/regional nisbətə keçid", self.zero_sum_fallback_cells),
        ]
        for label, count in checks:
            if count == 0:
                continue
            rate = self.rate(count)
            messages.append(f"{label}: {count} hadisə ({rate * 100:.1f}% hüceyrə).")
            if rate > warn_threshold:
                messages.append(
                    f"GÜCLÜ XƏBƏRDARLIQ: \"{label}\" nisbəti ({rate * 100:.1f}%) "
                    f"həddi ({warn_threshold * 100:.0f}%) aşır — variogram/qonşuluq/ehtimal "
                    "modeli ƏDƏDİ CƏHƏTDƏN QEYRİ-SABİT ola bilər. Nəticəyə ehtiyatla yanaşın.")
        return messages


@dataclass
class FaciesVariogramParams:
    """Bir fasiya üçün indikator variogram + anizotropluq parametrləri.

    Heç biri verilməyəndə (hamısı `None`) `simulate_sis` sərt datadan
    AVTOMATİK fit edir (bax `_resolve_facies_variogram`) — məlumat
    kifayət deyilsə açıq xəbərdarlıqla domen/3 evristikasına keçir.
    """
    model: str = MODEL_SPHERICAL
    nugget: float = 0.0
    sill: Optional[float] = None
    range_: Optional[float] = None
    range_v: Optional[float] = None
    azimuth_deg: Optional[float] = None
    range_minor: Optional[float] = None


@dataclass
class FaciesRealization:
    """Bir SIS icrasının nəticəsi — `realization_id`/`seed` təkrarlana
    bilənliyi TƏMİN EDİR (eyni giriş + eyni seed = eyni nəticə)."""
    realization_id: int
    seed: int
    codes: np.ndarray                          # (n_targets,) simulyasiya edilmiş fasiya kodları
    requested_proportions: Dict[int, float]
    realized_proportions: Dict[int, float]
    hard_data_mask: np.ndarray                 # (n_targets,) bool — sərt data ilə üst-üstə düşən hədəflər
    warnings: List[str] = field(default_factory=list)
    diagnostics: FaciesDiagnostics = field(default_factory=FaciesDiagnostics)


@dataclass
class FaciesProportions:
    """Qlobal + (istəyə görə) region/lay-asaslı fasiya nisbətləri.

    Prioritet: `layer_proportions` (varsa) > `region_proportions` (varsa)
    > `global_proportions`. Hər üçü (verilmişsə) AYRICA yoxlanılır (cəm=1).
    """
    global_proportions: Dict[int, float]
    region_proportions: Optional[Dict[int, Dict[int, float]]] = None
    layer_proportions: Optional[Dict[int, Dict[int, float]]] = None

    def validate(self) -> ValidationResult:
        result = ValidationResult()
        result.extend(validate_facies_proportions(self.global_proportions, "qlobal nisbətlər"))
        for region, props in (self.region_proportions or {}).items():
            result.extend(validate_facies_proportions(props, f"region {region} nisbətləri"))
        for layer, props in (self.layer_proportions or {}).items():
            result.extend(validate_facies_proportions(props, f"lay {layer} nisbətləri"))
        return result

    def for_cell(self, region: Optional[int] = None, layer: Optional[int] = None
                ) -> Dict[int, float]:
        if layer is not None and self.layer_proportions and layer in self.layer_proportions:
            return self.layer_proportions[layer]
        if region is not None and self.region_proportions and region in self.region_proportions:
            return self.region_proportions[region]
        return self.global_proportions


def indicator(codes: np.ndarray, facies: int) -> np.ndarray:
    """I_k(x) — `codes == facies` olan yerdə 1.0, əks halda 0.0."""
    return (np.asarray(codes) == facies).astype(float)


def observed_proportions(codes: np.ndarray) -> Dict[int, float]:
    """Sərt datadan müşahidə olunan fasiya nisbətləri (yalnız məlumat
    üçün/müqayisə üçün — modeli İDARƏ ETMİR)."""
    codes = np.asarray(codes)
    n = codes.size
    if n == 0:
        return {}
    unique, counts = np.unique(codes, return_counts=True)
    return {int(code): float(count) / n for code, count in zip(unique, counts)}


def _domain_span(points_xy: np.ndarray) -> float:
    if points_xy.shape[0] < 2:
        return 1.0
    lo = points_xy.min(axis=0)
    hi = points_xy.max(axis=0)
    return float(np.sqrt(np.sum((hi - lo) ** 2)))


def _resolve_facies_variogram(points_xy: np.ndarray, codes: np.ndarray, facies: int,
                              vp: Optional[FaciesVariogramParams], span: float
                              ) -> Tuple[dict, Optional[str]]:
    """Bir fasiyanın Kriging arqumentlərini müəyyənləşdirir.

    İstifadəçi `range_` veribsə (AÇIQ seçim) — birbaşa işlədilir.
    Verilməyibsə sərt datadan fit edilir; kifayət qədər nöqtə/dəyişkənlik
    yoxdursa (bax `fit_variogram_from_data`) AÇIQ xəbərdarlıqla domen/3
    evristikasına keçilir — SƏSSİZ heç nə UYDURULMUR.

    KRİTİK: `sill` BURADA HƏMİŞƏ KONKRET (sabit) ədədə HƏLL EDİLİR —
    heç vaxt `None` saxlanılmır. Səbəb: `OrdinaryKriging._parameters()`
    `sill=None` olanda onu HƏR ÇAĞIRIŞDA ötürülən nöqtələrin öz
    variansından (`np.var(values)`) YENİDƏN hesablayır — SIS-in yerli
    axtarışında bu, ÇAĞIRIŞDAN ÇAĞIRIŞA (qonşuluq alt-çoxluğundan asılı
    olaraq) DƏYİŞƏN, YƏNİ QEYRİ-STABİL sill deməkdir (stasionar variogram
    fərziyyəsini pozur) — məhz bu, sürətli (cKDTree) və brute-force
    axtarışın FƏRQLİ nəticə verməsinə səbəb olan HƏQİQİ SƏHV idi (tutulub
    və düzəldilib, bax `tests/test_facies_integration.py`).
    """
    ind = indicator(codes, facies)
    full_variance = max(float(np.var(ind)), 1e-12)

    if vp is not None and vp.range_ is not None:
        sill = vp.sill if vp.sill is not None else full_variance
        return dict(model=vp.model, nugget=vp.nugget, sill=sill, range_=vp.range_,
                    range_v=vp.range_v, azimuth_deg=vp.azimuth_deg,
                    range_minor=vp.range_minor), None
    model = vp.model if vp is not None else "auto"
    try:
        fit = fit_variogram_from_data(points_xy, ind, model=model)
        return dict(model=fit.model, nugget=fit.nugget, sill=fit.sill, range_=fit.range_,
                    range_v=(vp.range_v if vp else None),
                    azimuth_deg=(vp.azimuth_deg if vp else None),
                    range_minor=(vp.range_minor if vp else None)), None
    except ValueError as exc:
        fallback_range = max(span / 3.0, 1e-6)
        warning = (f"Fasiya {facies}: indikator variogram fit alınmadı ({exc}) — "
                  f"ehtiyat evristika (domen/3 = {fallback_range:.3g}) işlədildi.")
        return dict(model=MODEL_SPHERICAL, nugget=0.0, sill=full_variance, range_=fallback_range,
                    range_v=(vp.range_v if vp else None),
                    azimuth_deg=(vp.azimuth_deg if vp else None),
                    range_minor=(vp.range_minor if vp else None)), warning


#: Phase 5-də (SGS) da paylaşılan, xassə növündən asılı olmayan həndəsə
#: — bax `hard_data.find_exact_matches` modul docstring-i.
_find_hard_data_matches = find_exact_matches


def simulate_sis(points, codes, targets, proportions: Dict[int, float],
                 variograms: Optional[Dict[int, FaciesVariogramParams]] = None,
                 seed: int = 0, realization_id: int = 0,
                 search_radius: Optional[float] = None,
                 max_neighbors: Optional[int] = 24,
                 min_neighbors: int = 1,
                 hard_data_tolerance: float = 1e-6,
                 proportion_warn_threshold: float = 0.15,
                 correction_warn_threshold: float = DEFAULT_CORRECTION_WARN_THRESHOLD,
                 use_fast_search: bool = True,
                 rebuild_interval: int = 64) -> FaciesRealization:
    """Bir Sequential Indicator Simulation realizasiyası.

    Alqoritm (tapşırıqda tələb olunan addımlar):
        1. `seed`-dən qurulan RNG ilə TƏSADÜFİ simulyasiya yolu.
        2-8. Hər hədəf üçün: yerli kondisioner axtarışı, HƏR fasiya üçün
             indikator kriging sistemi, şərti ehtimal, [0,1]-ə kəsmə +
             normallaşdırma, kateqorik nümunələmə, nəticənin kondisioner
             çoxluğuna ƏLAVƏ edilməsi.
        9. Bütün hədəflər simulyasiya olunana qədər davam.

    Sərt data (targets-in points-a TAM üst-üstə düşən sətirləri)
    HEÇ VAXT simulyasiya YOLUNA daxil edilmir — birbaşa müşahidə
    dəyərini alır (bax `_find_hard_data_matches`).

    `use_fast_search=True` (defolt) — qonşu axtarışı `geology/
    spatial_search.IncrementalAnisotropicSearch` (cKDTree, Phase 4.1)
    ilə aparılır, YALNIZ seçilmiş yerli kondisioner alt-çoxluğu
    `OrdinaryKriging`-ə (qlobal həll rejimində) ötürülür — kriging
    RİYAZİYYATININ ÖZÜ dəyişmir, YALNIZ qonşu SEÇİMİ sürətlənir (bax
    `tests/test_spatial_search.py` — brute-force ilə EYNİ nəticə sübut
    edilib). `use_fast_search=False` — Phase 4-ün ƏSL yolu:
    `OrdinaryKriging`-in öz daxili (brute-force) yerli axtarışı. Hər
    ikisi EYNİ nəticəni verməlidir (bax
    `tests/test_facies_integration.py::test_fast_and_brute_force_search_agree`).

    Mürəkkəblik: `use_fast_search=False` təxminən `O(K·n²)` (bax Phase 4
    hesabatı); `use_fast_search=True` axtarış hissəsini təxminən
    `O(K·n log n)`-ə endirir (bax `FACIES.md` "Mürəkkəblik (Phase 4.1)").
    """
    facies_list = sorted(proportions)
    global_check = validate_facies_proportions(proportions, "qlobal nisbətlər")
    if not global_check.ok:
        raise ValueError("Fasiya nisbətləri etibarsızdır: " + "; ".join(global_check.errors))

    points = np.atleast_2d(np.asarray(points, float))
    codes = np.asarray(codes, int).ravel()
    targets = np.atleast_2d(np.asarray(targets, float))
    if points.shape[0] != codes.shape[0]:
        raise ValueError(
            f"points ({points.shape[0]}) və codes ({codes.shape[0]}) uzunluğu uyğun gəlmir.")
    unknown = sorted(set(np.unique(codes)) - set(facies_list))
    if unknown:
        raise ValueError(
            f"Sərt datada `proportions`-da olmayan fasiya kodu var: {unknown}.")

    warnings: List[str] = list(global_check.warnings)
    if codes.size:
        warnings.extend(compare_observed_vs_requested_proportions(
            observed_proportions(codes), proportions,
            warn_threshold=proportion_warn_threshold))

    n_targets = targets.shape[0]
    hard_index = _find_hard_data_matches(points, targets, hard_data_tolerance)
    hard_mask = hard_index >= 0
    simulated = np.full(n_targets, -1, dtype=int)
    simulated[hard_mask] = codes[hard_index[hard_mask]]

    to_simulate = np.where(~hard_mask)[0]
    if to_simulate.size == 0:
        realized = observed_proportions(simulated)
        return FaciesRealization(realization_id, seed, simulated, dict(proportions),
                                 realized, hard_mask, warnings)

    if len(facies_list) == 1:
        # tək fasiya: qeyri-müəyyənlik yoxdur, kriging/nümunələmə lazımsızdır
        simulated[to_simulate] = facies_list[0]
        realized = observed_proportions(simulated)
        return FaciesRealization(realization_id, seed, simulated, dict(proportions),
                                 realized, hard_mask, warnings)

    # `OrdinaryKriging.supports_z=True` — (n,2) verilsə Z=0 ilə (n,3)-ə
    # PADDING edir (`_as_points`), amma bizim öz axtarış strukturumuz
    # (`spatial_search.py`) bunu ETMİR — həmişə eyni sütun sayını
    # qoruyaq deyə burda ÖZÜMÜZ eyni şəkildə paddinq edirik.
    points3 = points if points.shape[1] >= 3 else np.column_stack(
        [points, np.zeros(points.shape[0])])
    targets3 = targets if targets.shape[1] >= 3 else np.column_stack(
        [targets, np.zeros(targets.shape[0])])
    points_xy = points[:, :2] if points.shape[1] >= 2 else points

    span = _domain_span(points_xy)
    krigers: Dict[int, OrdinaryKriging] = {}
    global_krigers: Dict[int, OrdinaryKriging] = {}
    searches: Dict[int, IncrementalAnisotropicSearch] = {}
    for k in facies_list:
        vp = (variograms or {}).get(k)
        kwargs, warn = _resolve_facies_variogram(points_xy, codes, k, vp, span)
        if warn:
            warnings.append(warn)
        if use_fast_search:
            range_h = kwargs["range_"]
            range_v = kwargs["range_v"] if kwargs.get("range_v") is not None else range_h
            range_minor = (kwargs["range_minor"] if kwargs.get("range_minor") is not None
                          else range_h)
            azimuth = kwargs["azimuth_deg"] if kwargs.get("azimuth_deg") is not None else 0.0
            aniso = AnisotropyParams(azimuth_deg=azimuth, range_major=range_h,
                                     range_minor=range_minor, range_vertical=range_v)
            searches[k] = IncrementalAnisotropicSearch(points3, anisotropy=aniso,
                                                       rebuild_interval=rebuild_interval)
            global_krigers[k] = OrdinaryKriging(search_radius=None, max_neighbors=None,
                                               min_neighbors=1, **kwargs)
        else:
            krigers[k] = OrdinaryKriging(search_radius=search_radius, max_neighbors=max_neighbors,
                                         min_neighbors=min_neighbors, **kwargs)

    rng = np.random.default_rng(seed)
    path = rng.permutation(to_simulate)

    sim_points = points3.copy()
    sim_codes = codes.copy()
    diag = FaciesDiagnostics()

    for idx in path:
        target_point = targets3[idx:idx + 1]
        diag.n_cells_simulated += 1
        raw_probs = np.full(len(facies_list), np.nan)
        insufficient_neighbors = False

        for i, k in enumerate(facies_list):
            ind_values = (sim_codes == k).astype(float)
            if use_fast_search:
                neighbor_idx = searches[k].query(target_point, search_radius=search_radius,
                                                 max_neighbors=max_neighbors,
                                                 min_neighbors=min_neighbors)
                if neighbor_idx.size == 0:
                    insufficient_neighbors = True
                    continue
                raw_probs[i] = global_krigers[k].interpolate(
                    sim_points[neighbor_idx], ind_values[neighbor_idx], target_point)[0]
            else:
                raw_probs[i] = krigers[k].interpolate(sim_points, ind_values, target_point)[0]

        if insufficient_neighbors or np.any(~np.isfinite(raw_probs)):
            probs = np.array([proportions[k] for k in facies_list])
            diag.nan_fallback_cells += 1
        else:
            if np.any(raw_probs < -1e-9):
                diag.negative_probability_events += 1
            clipped = np.clip(raw_probs, 0.0, None)
            total = float(clipped.sum())
            if total > 1.0 + 1e-9:
                diag.excess_probability_events += 1
            if total <= 1e-12:
                probs = np.array([proportions[k] for k in facies_list])
                diag.zero_sum_fallback_cells += 1
            else:
                probs = clipped / total

        # `rng.choice` numpy-nin daxili tolerantlığı (~1.49e-8) bizim
        # `validate_facies_proportions`-ın tolerantlığından (1e-6) SƏRTDİR
        # — son anda YENİDƏN (ehtiyat) normallaşdırma ilə "cəm 1 deyil"
        # xətasının qarşısı alınır (dəyər artıq [0,1]-dədir, YALNIZ son-bit
        # dəqiqləşdirmədir).
        probs = probs / probs.sum()
        chosen = facies_list[int(rng.choice(len(facies_list), p=probs))]
        simulated[idx] = chosen
        sim_points = np.vstack([sim_points, target_point])
        sim_codes = np.append(sim_codes, chosen)
        if use_fast_search:
            for k in facies_list:
                searches[k].add_point(target_point)

    warnings.extend(diag.summary_warnings(correction_warn_threshold))

    realized = observed_proportions(simulated)
    return FaciesRealization(realization_id, seed, simulated, dict(proportions),
                             realized, hard_mask, warnings, diag)


def run_realizations(n_realizations: int, points, codes, targets, proportions: Dict[int, float],
                     seed: int = 0, **kwargs) -> List[FaciesRealization]:
    """`n_realizations` MÜSTƏQİL SIS realizasiyası — hər biri
    `seed + i*1000` ilə (bax `application/scenarios.py`-dəki eyni
    konvensiya). Memarlıq N=1-dən N=100+-a qədər DƏYİŞİKLİKSİZ işləyir
    — burada PARALLELLƏŞDİRMƏ YOXDUR (bilərəkdən, bax FACIES.md)."""
    return [simulate_sis(points, codes, targets, proportions, seed=seed + i * 1000,
                         realization_id=i, **kwargs)
            for i in range(n_realizations)]


def summarize_realized_proportions(realizations: Sequence[FaciesRealization]
                                   ) -> Dict[int, Dict[str, float]]:
    """Hər fasiya üçün: tələb olunan nisbət + reallaşan nisbətlərin
    mean/std/min/max (bax tapşırıq §12 — TƏK realizasiya tələb olunan
    nisbəti DƏQİQ vermir, ORTALAMA yaxınlaşmalıdır)."""
    if not realizations:
        raise ValueError("Ən azı bir realizasiya lazımdır.")
    facies_list = sorted(realizations[0].requested_proportions)
    summary = {}
    for k in facies_list:
        values = np.array([r.realized_proportions.get(k, 0.0) for r in realizations])
        summary[k] = {
            "requested": realizations[0].requested_proportions[k],
            "mean": float(values.mean()),
            "std": float(values.std()),
            "min": float(values.min()),
            "max": float(values.max()),
        }
    return summary
