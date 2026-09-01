"""İnterpolyasiya interfeysi.

Geoloji modelləşdirmə qatı konkret alqoritmdən deyil, bu müqavilədən
asılıdır. Gələcəkdə maşın öyrənmə əsaslı və ya çoxdəyişənli
(co-kriging, seysmik trendlə) interpolyator əlavə etmək üçün yalnız
yeni implementasiya yazmaq kifayətdir.
"""

from __future__ import annotations
from abc import ABC, abstractmethod

import numpy as np


class IPropertyInterpolator(ABC):

    name: str = "interpolator"

    #: `True` olan implementasiyalar (n,3) [X,Y,Z] nöqtə/hədəf qəbul edir
    #: (bax `OrdinaryKriging`). `False` qalanlar (NearestNeighbour, IDW)
    #: yalnız (n,2)-ni dəstəkləyir — geology_service.py bunu yoxlayıb
    #: hansı fəzada nöqtə quracağına qərar verir.
    supports_z: bool = False

    @abstractmethod
    def interpolate(self, points: np.ndarray, values: np.ndarray,
                    targets: np.ndarray) -> np.ndarray:
        """points (n,2) və ya (n,3) [supports_z], values (n,),
        targets (m,2)/(m,3) -> (m,) dəyərlər."""

    def describe(self) -> str:
        return self.name
