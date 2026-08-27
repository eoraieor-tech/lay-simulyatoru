"""BrooksCoreyCapillaryProvider — A4.

    Pc(Sw) = Pe · Sw_norm^(-1/λ),    Sw_norm = (Sw − Swc) / (1 − Swc − Sor)

Sw → Swc olduqda Pc sonsuzluğa gedir, ona görə yuxarı hədd qoyulur.
Törəmə analitikdir — gələcək fully implicit sxem (A6) üçün Jakobianda
lazım olacaq.
"""

from __future__ import annotations
from typing import Optional

import numpy as np

from ..domain.scal import CapillaryParameters, CoreyParameters
from ..interfaces.providers import ICapillaryPressureProvider


class BrooksCoreyCapillaryProvider(ICapillaryPressureProvider):

    def __init__(self, capillary: CapillaryParameters, scal: CoreyParameters):
        issues = capillary.validate()
        if issues:
            raise ValueError("Kapilyar model yararsızdır: " + "; ".join(issues))
        self.params = capillary
        self.scal = scal
        self._span = max(1.0 - scal.swc - scal.sor, 1e-8)
        # Pc yuxarı həddinə uyğun minimal normallaşdırılmış doyumluluq
        self._sn_min = (capillary.entry_pressure /
                        capillary.max_pressure) ** capillary.lambda_exponent \
            if capillary.enabled else 1e-6

    def _normalized(self, sw) -> np.ndarray:
        sn = (np.asarray(sw, float) - self.scal.swc) / self._span
        return np.clip(sn, self._sn_min, 1.0)

    def pcow(self, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        p = self.params
        sn = self._normalized(sw)
        return np.minimum(p.entry_pressure * sn ** (-1.0 / p.lambda_exponent),
                          p.max_pressure)

    def dpcow_dsw(self, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        p = self.params
        sn = self._normalized(sw)
        raw = p.entry_pressure * sn ** (-1.0 / p.lambda_exponent)
        derivative = (-raw / (p.lambda_exponent * np.maximum(sn, 1e-12))
                      / self._span)
        return np.where(raw >= p.max_pressure, 0.0, derivative)

    def saturation_from_pc(self, pc) -> np.ndarray:
        """Pc → Sw tərs çevirmə (equilibration keçid zonası üçün)."""
        p = self.params
        pc = np.clip(np.asarray(pc, float), 1e-12, p.max_pressure)
        sn = np.where(pc <= p.entry_pressure, 1.0,
                      (p.entry_pressure / pc) ** p.lambda_exponent)
        return self.scal.swc + np.clip(sn, 0.0, 1.0) * self._span
