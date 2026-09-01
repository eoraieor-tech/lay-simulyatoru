"""Məkan interpolyasiyası — IPropertyInterpolator implementasiyaları.

Üç üsul, artan mürəkkəblik sırası ilə:

    NearestNeighbour  — ən yaxın quyunun dəyəri. Sürətli, pilləli.
    InverseDistance   — 1/d^p çəkiləri. Hamar, amma quyular arasında
                        həmişə orta dəyərə meyl edir ("öküz gözü" effekti).
    OrdinaryKriging   — variogram əsaslı, statistik optimal xətti qiymət.
                        Məkan korrelyasiyasını nəzərə alır.

Keçiricilik üçün `log_transform=True` tövsiyə olunur: keçiricilik
log-normal paylanır, ona görə interpolyasiya ln(k) fəzasında aparılır.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from ..interfaces.interpolation import IPropertyInterpolator
from .variogram import (MODEL_FUNCS, MODEL_SPHERICAL, AnisotropyDetectionResult,
                        AnisotropyParams, VariogramParameters, detect_anisotropy,
                        fit_variogram_from_data)


def _distance_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """(n,2) və (m,2) arasında Evklid məsafələri -> (n,m)."""
    diff = a[:, None, :] - b[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=-1))


class NearestNeighbour(IPropertyInterpolator):
    name = "Ən yaxın qonşu"

    def interpolate(self, points, values, targets) -> np.ndarray:
        points = np.asarray(points, float).reshape(-1, 2)
        values = np.asarray(values, float)
        targets = np.asarray(targets, float).reshape(-1, 2)
        distances = _distance_matrix(targets, points)
        return values[np.argmin(distances, axis=1)]


@dataclass
class InverseDistance(IPropertyInterpolator):
    """IDW: w_i = 1 / d_i^p."""
    power: float = 2.0
    search_radius: Optional[float] = None
    name: str = "Əks məsafə (IDW)"

    def interpolate(self, points, values, targets) -> np.ndarray:
        points = np.asarray(points, float).reshape(-1, 2)
        values = np.asarray(values, float)
        targets = np.asarray(targets, float).reshape(-1, 2)

        distances = _distance_matrix(targets, points)
        exact = distances < 1e-9
        distances = np.maximum(distances, 1e-9)

        weights = 1.0 / distances ** self.power
        if self.search_radius:
            weights = np.where(distances <= self.search_radius, weights, 0.0)
            empty = weights.sum(axis=1) <= 0.0
            if np.any(empty):
                # radius daxilində quyu yoxdursa, ən yaxınına qayıdırıq
                nearest = np.argmin(distances[empty], axis=1)
                weights[empty] = 0.0
                weights[np.where(empty)[0], nearest] = 1.0

        result = (weights * values[None, :]).sum(axis=1) / weights.sum(axis=1)

        # quyu nöqtəsində dəqiq dəyər qaytarılır
        hit_rows, hit_columns = np.where(exact)
        result[hit_rows] = values[hit_columns]
        return result


@dataclass
class OrdinaryKriging(IPropertyInterpolator):
    """Adi kriging — çoxlu variogram modeli, avtomatik fit, tam 3D
    anizotropluq (azimut + major/minor + şaquli), yerli axtarış dəstəyi ilə.

        γ(h) = c0 + c·[1.5·(h/a) − 0.5·(h/a)³],  h < a      (sferik, defolt)
        γ(h) = c0 + c,                            h ≥ a

    Eksponensial/qauss modelləri də var (bax `geology/variogram.py`).
    c0 — nugget, c — sill, a — major üfüqi radius (`range_`).

    **Parametr mənbəyi** — `range_`/`sill` birbaşa verilə bilər; verilməyib
    `auto_fit=True` olsa `geology/variogram.py`-dəki deneysel variogram +
    çəkili ən kiçik kvadrat fit işə düşür (bax `fit_variogram_from_data`);
    fit mümkün olmasa (nöqtə azdır) AÇIQ xəbərdarlıqla köhnə `domen/3`
    evristikasına geri qayıdılır. `auto_fit=False` (defolt) olanda dəyişən
    YOXDUR — köhnə `domen/3`/`var(dəyərlər)` evristikası birbaşa işlədilir,
    tamamilə əvvəlki davranış.

    **3D/tam anizotropluq** — `points`/`targets` (n,2) [X,Y] və ya (n,3)
    [X,Y,Z] ola bilər. `range_v` (şaquli radius) verilməyibsə `range_`-ə
    bərabər qəbul edilir. `azimuth_deg`/`range_minor` verilməyəndə (defolt)
    üfüqi müstəvi İZOTROPDUR — yalnız Z miqyaslanır (əvvəlki M2 davranışı,
    bax `geology/variogram.AnisotropyParams.transform` docstring-i: bu hal
    ədəd-ədəd əvvəlki nəticəni verir). `azimuth_deg` verilsə (və ya
    `auto_detect_anisotropy=True` etibarlı nəticə tapsa) üfüqi müstəvidə də
    major/minor radiuslar fərqli işlədilir — tam geometrik anizotropluq
    transformu (`AnisotropyParams.transform`).

    **Yerli axtarış (moving neighbourhood)** — `search_radius` və/və ya
    `max_neighbors` veriləndə hər hədəf üçün YALNIZ ən yaxın uyğun
    nöqtələrlə AYRICA yerli kriging sistemi qurulur (qlobal sistemin
    əvəzinə). `min_neighbors`-dan az uyğun nöqtə olan hədəf üçün NaN
    qaytarılır — dəyər UYDURULMUR, sadəcə "bu nöqtədə etibarlı proqnoz
    yoxdur" bildirilir. Hamısı `None`/defolt olanda nəticə köhnə QLOBAL
    kriging ilə BİRƏBİR eynidir (bax `_solve_global`) — bu yol reqressiya
    testləri ilə qorunur.
    """
    range_: Optional[float] = None
    range_v: Optional[float] = None
    sill: Optional[float] = None
    nugget: float = 0.0
    model: str = MODEL_SPHERICAL
    auto_fit: bool = False
    auto_fit_nugget: bool = False
    azimuth_deg: Optional[float] = None
    range_minor: Optional[float] = None
    auto_detect_anisotropy: bool = False
    search_radius: Optional[float] = None
    min_neighbors: int = 1
    max_neighbors: Optional[int] = None
    name: str = "Kriging (adi)"

    supports_z = True

    #: sonuncu `interpolate()` çağırışının introspeksiya/hesabat nəticələri
    #: (istəyə görə UI/PDF hesabatda göstərilə bilər) — bax `_parameters`.
    last_fit_: Optional[VariogramParameters] = field(default=None, repr=False, compare=False)
    last_anisotropy_: Optional[AnisotropyDetectionResult] = field(
        default=None, repr=False, compare=False)
    last_warnings_: list = field(default_factory=list, repr=False, compare=False)

    def __post_init__(self):
        if self.model == "auto" and not self.auto_fit:
            raise ValueError("model='auto' üçün auto_fit=True lazımdır.")
        if self.model != "auto" and self.model not in MODEL_FUNCS:
            raise ValueError(
                f"Naməlum variogram modeli: {self.model!r}. "
                f"Dəstəklənən: {tuple(MODEL_FUNCS) + ('auto',)}")

    def _variogram(self, h: np.ndarray, range_: float, sill: float, nugget: float,
                   model: str, zero_at_origin: bool = True) -> np.ndarray:
        """`zero_at_origin` yalnız məlumat-məlumat matrisi üçün doğrudur:
        nöqtənin özü ilə arasındakı variogram sıfırdır. Sağ tərəfdə
        (məlumat-hədəf) isə nugget saxlanılır — məhz buna görə nugget
        sıfırdan böyük olanda kriging DƏQİQ interpolyator olmaqdan
        çıxır və ölçmə səhvini süzür.
        """
        func = MODEL_FUNCS[model]
        gamma = func(h, nugget, sill, range_)
        if zero_at_origin:
            return np.where(h <= 1e-12, 0.0, gamma)
        return gamma

    def _dedupe_conflicting_points(self, points3: np.ndarray, values: np.ndarray):
        """Tam üst-üstə düşən (məsafə < 1e-9) giriş nöqtələrini araşdırır.

        Əvvəlki sükut davranış: `_solve_global`-da tam-üst-üstə-düşən
        hədəflər üçün dəyər `result[rows] = values[columns]` fancy-
        indexing ilə bərpa edilirdi — bir neçə eyni-koordinatlı GİRİŞ
        nöqtəsi ziddiyyətli dəyərlə verildikdə, hansının "qazanacağı"
        NumPy-ın sənədləşdirilməmiş sırasına görə həll olunurdu.

        İndi: ziddiyyətli dublikatlar DETERMİNİSTİK olaraq ORTALANIR və
        `last_warnings_`-ə yazılır; dəyərləri praktik eyni olan
        dublikatlar səssizcə (xəbərdarlıqsız) birləşdirilir.

        Bu funksiya HƏR `interpolate()`/`interpolate_with_variance()`
        çağırışında (o cümlədən SGS/SIS-in artan-nöqtə simulyasiya
        dövründə hər addımda, bax `sgs.py`/`facies.py`) işə düşür, ona
        görə İKİ PİLLƏLİDİR: əvvəlcə ucuz `O(n log n)` lexsort-əsaslı
        "heç bir dublikat yoxdurmu?" yoxlaması (adi, dublikatsız hal —
        demək olar bütün çağırışlar) sürətlə keçir; yalnız faktiki
        dublikat tapılan NADIR halda daha ətraflı (yenə `O(n log n)`,
        `np.unique` hash-əsaslı) qruplaşdırma işə düşür. Tam cüt-cüt
        (n,n) məsafə matrisi HEÇ VAXT hesablanmır — bu, əvvəlki versiyada
        `test_facies_performance.py`-ni O(n³)-ə qədər yavaşladan real
        performans reqressiyası idi (performans auditində tapılıb və
        canlı təsdiqlənib), buna görə iki pilləli struktur seçildi.
        """
        n = points3.shape[0]
        if n < 2:
            return points3, values, []

        order = np.lexsort((points3[:, 2], points3[:, 1], points3[:, 0]))
        adjacent_close = np.all(
            np.abs(np.diff(points3[order], axis=0)) < 1e-9, axis=1)
        if not np.any(adjacent_close):
            return points3, values, []

        rounded = np.round(np.ascontiguousarray(points3), 9)
        view = rounded.view([("", rounded.dtype)] * rounded.shape[1])
        _, first_idx, inverse, counts = np.unique(
            view, return_index=True, return_inverse=True, return_counts=True)
        inverse = inverse.reshape(-1)

        keep = np.zeros(n, dtype=bool)
        keep[first_idx] = True
        resolved = values.copy()
        warnings: list = []
        for group_id in np.where(counts > 1)[0]:
            idx = np.where(inverse == group_id)[0]
            group_values = values[idx]
            if not np.allclose(group_values, group_values[0], atol=1e-9, rtol=0.0):
                xy = points3[idx[0], :2]
                warnings.append(
                    f"Eyni koordinatda ({xy[0]:.4g}, {xy[1]:.4g}) {idx.size} "
                    f"ziddiyyətli sərt-data nöqtəsi tapıldı (dəyərlər: "
                    f"{group_values.tolist()}) — deterministik ORTA dəyər "
                    "istifadə olundu (əvvəlki sükut 'son yazı qazanır' "
                    "davranışı əvəzinə).")
            resolved[first_idx[group_id]] = float(np.mean(group_values))

        return points3[keep], resolved[keep], warnings

    def _parameters(self, points_xy: np.ndarray, values: np.ndarray):
        """`points_xy` — YALNIZ X,Y (Z auto-range təxminini korlamasın).

        Qaytarır: `(range_h, range_v, range_minor, azimuth_deg, sill,
        nugget, model)` — hamısı bu çağırış üçün İSTİFADƏ OLUNACAQ faktiki
        dəyərlər (fit/aşkarlanma nəticələri daxil).
        """
        warnings: list = []
        fit: Optional[VariogramParameters] = None
        model = self.model if self.model != "auto" else MODEL_SPHERICAL

        if self.range_ is not None:
            range_h = float(self.range_)
        elif self.auto_fit:
            try:
                fit = fit_variogram_from_data(points_xy, values, model=self.model)
                range_h = fit.range_
                model = fit.model
            except ValueError as exc:
                warnings.append(
                    f"Avtomatik variogram fit alınmadı ({exc}) — ehtiyat evristika "
                    "(domen/3) işlədildi.")
                span = _distance_matrix(points_xy, points_xy).max()
                range_h = max(span / 3.0, 1e-6)
        else:
            span = _distance_matrix(points_xy, points_xy).max()
            range_h = max(span / 3.0, 1e-6)

        range_v = float(self.range_v) if self.range_v is not None else range_h
        sill = (max(fit.sill, 1e-12) if (fit is not None and self.sill is None)
               else max(float(self.sill) if self.sill is not None else float(np.var(values)),
                        1e-12))
        nugget = fit.nugget if (fit is not None and self.auto_fit_nugget) else self.nugget
        range_minor = float(self.range_minor) if self.range_minor is not None else range_h
        azimuth_deg = float(self.azimuth_deg) if self.azimuth_deg is not None else 0.0

        anisotropy: Optional[AnisotropyDetectionResult] = None
        if self.auto_detect_anisotropy and self.azimuth_deg is None:
            try:
                anisotropy = detect_anisotropy(points_xy, values)
                if anisotropy.reliable:
                    azimuth_deg = anisotropy.azimuth_deg
                    if self.range_ is None:
                        range_h = anisotropy.range_major
                    if self.range_minor is None:
                        range_minor = anisotropy.range_minor
                else:
                    warnings.extend(anisotropy.warnings)
            except ValueError as exc:
                warnings.append(f"Anizotropluq aşkarlanması alınmadı: {exc}")

        self.last_fit_ = fit
        self.last_anisotropy_ = anisotropy
        self.last_warnings_ = warnings
        return range_h, range_v, range_minor, azimuth_deg, sill, nugget, model

    @staticmethod
    def _as_points(arr) -> np.ndarray:
        """(n,2) -> Z=0 sütunu ilə (n,3); (n,3) olduğu kimi qalır."""
        arr = np.asarray(arr, float)
        if arr.ndim != 2:
            raise ValueError("Nöqtələr 2D massiv olmalıdır: (n,2) və ya (n,3)")
        if arr.shape[1] == 2:
            return np.column_stack([arr, np.zeros(arr.shape[0])])
        if arr.shape[1] == 3:
            return arr
        raise ValueError(f"Nöqtələr (n,2) və ya (n,3) olmalıdır, alındı: {arr.shape}")

    def interpolate(self, points, values, targets) -> np.ndarray:
        points3 = self._as_points(points)
        values = np.asarray(values, float)
        targets3 = self._as_points(targets)
        points3, values, dup_warnings = self._dedupe_conflicting_points(points3, values)
        n = points3.shape[0]
        if n == 1:
            self.last_warnings_ = dup_warnings
            return np.full(targets3.shape[0], values[0])

        (range_h, range_v, range_minor, azimuth_deg, sill, nugget,
         model) = self._parameters(points3[:, :2], values)
        self.last_warnings_ = dup_warnings + self.last_warnings_
        anisotropy = AnisotropyParams(azimuth_deg=azimuth_deg, range_major=range_h,
                                      range_minor=range_minor, range_vertical=range_v)
        scaled_points = anisotropy.transform(points3)
        scaled_targets = anisotropy.transform(targets3)

        if self.search_radius is None and self.max_neighbors is None:
            return self._solve_global(scaled_points, values, scaled_targets, range_h,
                                      sill, nugget, model)
        return self._solve_local(scaled_points, values, scaled_targets, range_h, sill,
                                 nugget, model)

    def _solve_global(self, points, values, targets, range_h, sill, nugget, model,
                      return_variance: bool = False):
        """Bütün nöqtələrlə TƏK kriging sistemi — parametrsiz köhnə yol.

        `return_variance=True` (defolt `False`, mövcud çağırışlar
        DƏYİŞMİR) — Phase 5 (SGS) üçün: kriging variansı `σ²(x0) =
        Σ w_i·γ(x_i,x0) + μ` ARTIQ HƏLL EDİLMİŞ sistemdən (`weights`,
        `right`, Laqranj vuruğu `solution[n,:]`) hesablanır — YENİDƏN
        HƏLL YOXDUR, sadəcə eyni nəticədən ƏLAVƏ oxuma.
        """
        n = points.shape[0]
        # sol tərəf: variogram matrisi + Laqranj sətri/sütunu
        left = np.ones((n + 1, n + 1))
        left[:n, :n] = self._variogram(_distance_matrix(points, points), range_h,
                                       sill, nugget, model, zero_at_origin=True)
        left[n, n] = 0.0

        right = np.ones((n + 1, targets.shape[0]))
        right[:n, :] = self._variogram(_distance_matrix(points, targets), range_h,
                                       sill, nugget, model, zero_at_origin=nugget <= 0.0)

        try:
            solution = np.linalg.solve(left, right)
        except np.linalg.LinAlgError:
            solution = np.linalg.lstsq(left, right, rcond=None)[0]

        weights = solution[:n, :]
        result = (weights * values[:, None]).sum(axis=0)

        # nugget sıfırdırsa kriging dəqiq interpolyatordur — nöqtələri bərpa edirik
        exact_rows = None
        if nugget <= 0.0:
            hits = _distance_matrix(targets, points) < 1e-9
            rows, columns = np.where(hits)
            result[rows] = values[columns]
            exact_rows = rows

        if not return_variance:
            return result

        variance = np.sum(weights * right[:n, :], axis=0) + solution[n, :]
        variance = np.clip(variance, 0.0, None)
        if exact_rows is not None:
            variance[exact_rows] = 0.0
        return result, variance

    def _solve_local(self, points, values, targets, range_h, sill, nugget, model,
                     return_variance: bool = False):
        """Hər hədəf üçün ayrıca yerli (moving neighbourhood) sistem.

        `points`/`targets` artıq anizotrop-miqyaslanmışdır, ona görə
        məsafə/radius bu miqyaslanmış fəzada ölçülür — variogramın özü
        işlətdiyi fəza ilə eynidir.
        """
        distances = _distance_matrix(targets, points)   # (m, n)
        result = np.full(targets.shape[0], np.nan)
        variance = np.full(targets.shape[0], np.nan) if return_variance else None
        for row in range(targets.shape[0]):
            dist_row = distances[row]
            candidate = np.arange(points.shape[0])
            if self.search_radius is not None:
                candidate = candidate[dist_row[candidate] <= self.search_radius]
            if candidate.size == 0:
                continue
            candidate = candidate[np.argsort(dist_row[candidate])]
            if self.max_neighbors is not None:
                candidate = candidate[: self.max_neighbors]
            if candidate.size < max(self.min_neighbors, 1):
                continue
            if return_variance:
                local, local_var = self._solve_global(
                    points[candidate], values[candidate], targets[row:row + 1],
                    range_h, sill, nugget, model, return_variance=True)
                variance[row] = local_var[0]
            else:
                local = self._solve_global(points[candidate], values[candidate],
                                           targets[row:row + 1], range_h, sill, nugget, model)
            result[row] = local[0]
        if return_variance:
            return result, variance
        return result

    def interpolate_with_variance(self, points, values, targets):
        """`interpolate()`-in EYNİSİ, ƏLAVƏ olaraq kriging VARİANSINI da
        qaytarır (Phase 5/SGS: şərti Gauss paylanması `N(estimate,
        variance)` üçün tələb olunur). `interpolate()`-in ÖZÜ bu metodla
        DƏYİŞMƏYİB — tamamilə yeni, əlavə giriş nöqtəsidir.

        Qaytarır: `(estimate, variance)`, hər ikisi `(m,)` massiv.
        """
        points3 = self._as_points(points)
        values = np.asarray(values, float)
        targets3 = self._as_points(targets)
        points3, values, dup_warnings = self._dedupe_conflicting_points(points3, values)
        n = points3.shape[0]

        (range_h, range_v, range_minor, azimuth_deg, sill, nugget,
         model) = self._parameters(points3[:, :2], values)
        self.last_warnings_ = dup_warnings + self.last_warnings_
        anisotropy = AnisotropyParams(azimuth_deg=azimuth_deg, range_major=range_h,
                                      range_minor=range_minor, range_vertical=range_v)
        scaled_points = anisotropy.transform(points3)
        scaled_targets = anisotropy.transform(targets3)

        if n == 1:
            # tək nöqtə: `interpolate()` bu halda sistemi keçib birbaşa
            # dəyəri qaytarır (riyazi cəhətdən `_solve_global`-ın n=1
            # üçün verdiyi ilə EYNİDİR — çəki=1, Laqranj vuruğu=γ). AMMA
            # varians üçün bu QISA YOLU İZLƏMİRİK: `_solve_global`-ın
            # ÖZÜNÜ çağırırıq ki, `_solve_local`-ın (yerli axtarışda
            # dəqiq 1 qonşu tapılanda) verdiyi NƏTİCƏ İLƏ TAM UYĞUN olsun
            # — ilk versiya bunları FƏRQLİ (yanlış: sadə-kriging-bənzər
            # `γ`, düzgünü: `σ²=w·γ+μ=2γ`, n=1-də) hesablayırdı, bu,
            # sürətli/brute-force axtarışın FƏRQLİ nəticə verməsinə səbəb
            # olan HƏQİQİ SƏHV idi (tutulub, bax `tests/test_sgs_
            # validation.py`).
            return self._solve_global(scaled_points, values, scaled_targets, range_h,
                                      sill, nugget, model, return_variance=True)

        if self.search_radius is None and self.max_neighbors is None:
            return self._solve_global(scaled_points, values, scaled_targets, range_h,
                                      sill, nugget, model, return_variance=True)
        return self._solve_local(scaled_points, values, scaled_targets, range_h, sill,
                                 nugget, model, return_variance=True)


INTERPOLATORS = {
    "Ən yaxın qonşu": NearestNeighbour,
    "Əks məsafə (IDW)": InverseDistance,
    "Kriging (adi)": OrdinaryKriging,
}


def interpolate_property(interpolator: IPropertyInterpolator,
                         points: np.ndarray, values: np.ndarray,
                         targets: np.ndarray,
                         log_transform: bool = False,
                         minimum: Optional[float] = None,
                         maximum: Optional[float] = None) -> np.ndarray:
    """İnterpolyasiya + log çevirmə + hədlərin tətbiqi.

    Log çevirmə keçiricilik üçündür: ln(k) fəzasında interpolyasiya
    həm mənfi dəyərin qarşısını alır, həm də fiziki cəhətdən daha
    doğrudur, çünki keçiricilik log-normal paylanır.
    """
    values = np.asarray(values, float)
    if log_transform:
        if np.any(values <= 0):
            raise ValueError("Log interpolyasiya üçün bütün dəyərlər müsbət olmalıdır.")
        result = np.exp(interpolator.interpolate(points, np.log(values), targets))
    else:
        result = interpolator.interpolate(points, values, targets)

    if minimum is not None:
        result = np.maximum(result, minimum)
    if maximum is not None:
        result = np.minimum(result, maximum)
    return result
