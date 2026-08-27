"""MÜVƏQQƏTİ adapter: mövcud Corey kodunu provider interfeysinə bağlayır.

Bu YENİ modul deyil — domain.scal.CoreyParameters-dəki (əvvəlki nüvədən
köçürülmüş) düsturları IRelativePermeabilityProvider müqaviləsinə
uyğunlaşdırır ki, mühərrik konkret sinifdən asılı olmasın.

Real SCAL modulu (cədvəl oxuma, region üzrə fərqli əyrilər, histerezis)
yazılanda bu fayl silinir və yerinə həmin modul inject edilir.
"""

from __future__ import annotations
from typing import Optional

import numpy as np

from ..domain.scal import CoreyParameters
from ..interfaces.providers import IRelativePermeabilityProvider


class CoreyRelativePermeabilityAdapter(IRelativePermeabilityProvider):

    def __init__(self, parameters: CoreyParameters):
        self._p = parameters

    def krw(self, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._p.krw(sw)

    def kro(self, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._p.kro(sw)

    def krw_derivative(self, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        """dkrw/dSw — analitik.

            krw = krw_end · Sn^nw,   Sn = (Sw − Swc) / (1 − Swc − Sor)
            dkrw/dSw = krw_end · nw · Sn^(nw−1) / (1 − Swc − Sor)

        Hədlərdən kənarda (Sn = 0 və ya 1-dən kənar) törəmə sıfırdır,
        çünki krw orada sabitdir.
        """
        p = self._p
        span = max(1.0 - p.swc - p.sor, 1e-8)
        sw = np.asarray(sw, float)
        sn = (sw - p.swc) / span
        inside = (sn > 0.0) & (sn < 1.0)
        result = np.zeros_like(sn)
        clipped = np.clip(sn, 1e-12, 1.0)
        result[inside] = (p.krw_end * p.nw * clipped[inside] ** (p.nw - 1.0)
                          / span)
        return result

    def kro_derivative(self, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        """dkro/dSw — analitik (mənfi işarəli)."""
        p = self._p
        span = max(1.0 - p.swc - p.sor, 1e-8)
        sw = np.asarray(sw, float)
        sn = (sw - p.swc) / span
        inside = (sn > 0.0) & (sn < 1.0)
        result = np.zeros_like(sn)
        clipped = np.clip(1.0 - sn, 1e-12, 1.0)
        result[inside] = (-p.kro_end * p.no * clipped[inside] ** (p.no - 1.0)
                          / span)
        return result

    def saturation_limits(self, region: Optional[int] = None) -> tuple:
        return self._p.swc, 1.0 - self._p.sor

    def endpoint_water_mobility(self, water_viscosity: float,
                                region: Optional[int] = None) -> float:
        return self._p.krw_end / water_viscosity

    def max_fractional_flow_derivative(self, water_viscosity: float,
                                       oil_viscosity: float,
                                       region: Optional[int] = None) -> float:
        lo, hi = self.saturation_limits()
        s = np.linspace(lo, hi, 400)
        lw = self._p.krw(s) / water_viscosity
        loil = self._p.kro(s) / oil_viscosity
        fw = lw / np.maximum(lw + loil, 1e-30)
        return float(np.nanmax(np.abs(np.gradient(fw, s))))
