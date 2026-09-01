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
from dataclasses import dataclass
from typing import Optional

import numpy as np

from ..interfaces.interpolation import IPropertyInterpolator


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
    """Sferik variogramlı adi kriging — 3D, anizotrop, yerli axtarış dəstəyi ilə.

        γ(h) = c0 + c·[1.5·(h/a) − 0.5·(h/a)³],  h < a
        γ(h) = c0 + c,                            h ≥ a

    c0 — nugget (ölçmə səhvi / kiçik miqyaslı dəyişkənlik),
    c  — sill töhfəsi, a — üfüqi korrelyasiya radiusu (`range_`).

    `range_` verilməyibsə məlumatın (yalnız X,Y) ölçüsündən təxmin edilir.

    **3D/anizotropluq** — `points`/`targets` (n,2) [X,Y] və ya (n,3)
    [X,Y,Z] ola bilər. (n,2) veriləndə Z=0 qəbul edilir — köhnə 2D
    davranış BİRƏBİR qorunur (bax `_as_points`). `range_v` (şaquli
    radius) verilməyibsə `range_`-ə bərabər qəbul edilir (izotrop
    defolt) — bu, fiziki dəyər UYDURMAQ deyil, yalnız "əks halda nə
    fərz edək" sualına ən neytral cavabdır. `range_v != range_`
    olanda Z oxu `range_/range_v` əmsalı ilə miqyaslanır, sonra TƏK
    izotrop variogram miqyaslanmış məsafəyə tətbiq olunur — bu,
    geostatistikada standart "geometric anisotropy" transformudur.

    **Yerli axtarış (moving neighbourhood)** — `search_radius` və/və ya
    `max_neighbors` veriləndə hər hədəf üçün YALNIZ ən yaxın uyğun
    nöqtələrlə AYRICA yerli kriging sistemi qurulur (qlobal sistemin
    əvəzinə). `min_neighbors`-dan az uyğun nöqtə olan hədəf üçün NaN
    qaytarılır — dəyər UYDURULMUR, sadəcə "bu nöqtədə etibarlı proqnoz
    yoxdur" bildirilir. Hər ikisi `None`/defolt olanda (heç bir axtarış
    məhdudiyyəti yoxdursa) nəticə köhnə QLOBAL kriging ilə BİRƏBİR
    eynidir (bax `_solve_global`) — bu yol reqressiya testləri ilə
    qorunur.
    """
    range_: Optional[float] = None
    range_v: Optional[float] = None
    sill: Optional[float] = None
    nugget: float = 0.0
    search_radius: Optional[float] = None
    min_neighbors: int = 1
    max_neighbors: Optional[int] = None
    name: str = "Kriging (adi)"

    supports_z = True

    def _variogram(self, h: np.ndarray, range_: float, sill: float,
                   zero_at_origin: bool = True) -> np.ndarray:
        """Sferik variogram.

        `zero_at_origin` yalnız məlumat-məlumat matrisi üçün doğrudur:
        nöqtənin özü ilə arasındakı variogram sıfırdır. Sağ tərəfdə
        (məlumat-hədəf) isə nugget saxlanılır — məhz buna görə nugget
        sıfırdan böyük olanda kriging DƏQİQ interpolyator olmaqdan
        çıxır və ölçmə səhvini süzür.
        """
        ratio = np.clip(h / max(range_, 1e-9), 0.0, 1.0)
        gamma = self.nugget + sill * (1.5 * ratio - 0.5 * ratio ** 3)
        if zero_at_origin:
            return np.where(h <= 1e-12, 0.0, gamma)
        return gamma

    def _parameters(self, points_xy: np.ndarray, values: np.ndarray):
        """`points_xy` — YALNIZ X,Y (Z auto-range təxminini korlamasın)."""
        if self.range_ is not None:
            range_h = float(self.range_)
        else:
            span = _distance_matrix(points_xy, points_xy).max()
            range_h = max(span / 3.0, 1e-6)
        range_v = float(self.range_v) if self.range_v is not None else range_h
        sill = float(self.sill) if self.sill is not None else float(np.var(values))
        return range_h, range_v, max(sill, 1e-12)

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

    @staticmethod
    def _apply_anisotropy(points: np.ndarray, range_h: float, range_v: float) -> np.ndarray:
        if abs(range_v - range_h) < 1e-12:
            return points
        scaled = points.copy()
        scaled[:, 2] *= range_h / max(range_v, 1e-9)
        return scaled

    def interpolate(self, points, values, targets) -> np.ndarray:
        points3 = self._as_points(points)
        values = np.asarray(values, float)
        targets3 = self._as_points(targets)
        n = points3.shape[0]
        if n == 1:
            return np.full(targets3.shape[0], values[0])

        range_h, range_v, sill = self._parameters(points3[:, :2], values)
        scaled_points = self._apply_anisotropy(points3, range_h, range_v)
        scaled_targets = self._apply_anisotropy(targets3, range_h, range_v)

        if self.search_radius is None and self.max_neighbors is None:
            return self._solve_global(scaled_points, values, scaled_targets, range_h, sill)
        return self._solve_local(scaled_points, values, scaled_targets, range_h, sill)

    def _solve_global(self, points, values, targets, range_h, sill) -> np.ndarray:
        """Bütün nöqtələrlə TƏK kriging sistemi — parametrsiz köhnə yol."""
        n = points.shape[0]
        # sol tərəf: variogram matrisi + Laqranj sətri/sütunu
        left = np.ones((n + 1, n + 1))
        left[:n, :n] = self._variogram(_distance_matrix(points, points), range_h,
                                       sill, zero_at_origin=True)
        left[n, n] = 0.0

        right = np.ones((n + 1, targets.shape[0]))
        right[:n, :] = self._variogram(_distance_matrix(points, targets), range_h,
                                       sill, zero_at_origin=self.nugget <= 0.0)

        try:
            solution = np.linalg.solve(left, right)
        except np.linalg.LinAlgError:
            solution = np.linalg.lstsq(left, right, rcond=None)[0]

        weights = solution[:n, :]
        result = (weights * values[:, None]).sum(axis=0)

        # nugget sıfırdırsa kriging dəqiq interpolyatordur — nöqtələri bərpa edirik
        if self.nugget <= 0.0:
            hits = _distance_matrix(targets, points) < 1e-9
            rows, columns = np.where(hits)
            result[rows] = values[columns]
        return result

    def _solve_local(self, points, values, targets, range_h, sill) -> np.ndarray:
        """Hər hədəf üçün ayrıca yerli (moving neighbourhood) sistem.

        `points`/`targets` artıq anizotrop-miqyaslanmışdır, ona görə
        məsafə/radius bu miqyaslanmış fəzada ölçülür — variogramın özü
        işlətdiyi fəza ilə eynidir.
        """
        distances = _distance_matrix(targets, points)   # (m, n)
        result = np.full(targets.shape[0], np.nan)
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
            local = self._solve_global(points[candidate], values[candidate],
                                       targets[row:row + 1], range_h, sill)
            result[row] = local[0]
        return result


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
