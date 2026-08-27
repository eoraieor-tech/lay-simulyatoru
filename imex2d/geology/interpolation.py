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
    """Sferik variogramlı adi kriging.

        γ(h) = c0 + c·[1.5·(h/a) − 0.5·(h/a)³],  h < a
        γ(h) = c0 + c,                            h ≥ a

    c0 — nugget (ölçmə səhvi / kiçik miqyaslı dəyişkənlik),
    c  — sill töhfəsi, a — korrelyasiya radiusu (range).

    `range_` verilməyibsə məlumatın ölçüsündən təxmin edilir.
    """
    range_: Optional[float] = None
    sill: Optional[float] = None
    nugget: float = 0.0
    name: str = "Kriging (adi)"

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

    def _parameters(self, points: np.ndarray, values: np.ndarray):
        if self.range_ is not None:
            range_ = float(self.range_)
        else:
            span = _distance_matrix(points, points).max()
            range_ = max(span / 3.0, 1e-6)
        sill = float(self.sill) if self.sill is not None else float(np.var(values))
        return range_, max(sill, 1e-12)

    def interpolate(self, points, values, targets) -> np.ndarray:
        points = np.asarray(points, float).reshape(-1, 2)
        values = np.asarray(values, float)
        targets = np.asarray(targets, float).reshape(-1, 2)
        n = points.shape[0]
        if n == 1:
            return np.full(targets.shape[0], values[0])

        range_, sill = self._parameters(points, values)

        # sol tərəf: variogram matrisi + Laqranj sətri/sütunu
        left = np.ones((n + 1, n + 1))
        left[:n, :n] = self._variogram(_distance_matrix(points, points), range_,
                                       sill, zero_at_origin=True)
        left[n, n] = 0.0

        right = np.ones((n + 1, targets.shape[0]))
        right[:n, :] = self._variogram(_distance_matrix(points, targets), range_,
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
