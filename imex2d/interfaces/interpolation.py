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

    @abstractmethod
    def interpolate(self, points: np.ndarray, values: np.ndarray,
                    targets: np.ndarray) -> np.ndarray:
        """points (n,2), values (n,), targets (m,2) -> (m,) dəyərlər."""

    def describe(self) -> str:
        return self.name
