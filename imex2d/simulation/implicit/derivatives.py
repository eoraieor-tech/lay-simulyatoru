"""Törəmə təchizatçısı — Jakobian üçün lazım olan bütün törəmələr.

Prinsip: provider analitik törəmə verirsə, o işlədilir; verməyibsə
mərkəzi sonlu fərqə keçilir. Bu sayədə provider interfeysləri
dəyişmir (Open/Closed), lakin analitik törəmə mövcud olduqda dəqiqlik
və sürət qazanılır.

Hazırda analitik olanlar:
    dkrw/dSw, dkro/dSw   — Corey düsturundan
    dPc/dSw              — Brooks-Corey (A4-də yazılıb)

Ədədi olanlar:
    dBw/dp, dBo/dp, dμ/dp — PVT cədvəli parçalı xəttidir, ona görə
                            mərkəzi fərq kifayətdir
"""

from __future__ import annotations

import numpy as np

PRESSURE_STEP = 1e-4      # bar
SATURATION_STEP = 1e-6


class DerivativeProvider:
    """Provider-lərin törəmələrini vahid interfeys altında toplayır."""

    def __init__(self, relperm, pvt=None, capillary=None,
                 fluids=None, reference_pressure: float = 0.0,
                 pressure_step: float = PRESSURE_STEP,
                 saturation_step: float = SATURATION_STEP):
        self.relperm = relperm
        self.pvt = pvt
        self.capillary = capillary
        self.fluids = fluids
        self.reference_pressure = float(reference_pressure)
        self.dp = pressure_step
        self.ds = saturation_step

    # ══════════════════════════ statik sıxılma modeli (PVT olmayanda)
    def _static_fvf_derivative(self, pressure, reference_fvf: float,
                               compressibility: float) -> np.ndarray:
        """B(p) = B_ref / (1 + c·(p − p_ref))  ->  dB/dp = −c·B² / B_ref."""
        pressure = np.asarray(pressure, float)
        if self.fluids is None or compressibility <= 0.0:
            return np.zeros_like(pressure)
        factor = np.maximum(
            1.0 + compressibility * (pressure - self.reference_pressure), 1e-6)
        fvf = reference_fvf / factor
        return -compressibility * fvf ** 2 / reference_fvf

    # ═══════════════════════════════════════════ nisbi keçiricilik
    def dkrw_dsw(self, sw: np.ndarray) -> np.ndarray:
        analytic = getattr(self.relperm, "krw_derivative", None)
        if analytic is not None:
            return np.asarray(analytic(sw), float)
        return self._central(lambda s: self.relperm.krw(s), sw, self.ds)

    def dkro_dsw(self, sw: np.ndarray) -> np.ndarray:
        analytic = getattr(self.relperm, "kro_derivative", None)
        if analytic is not None:
            return np.asarray(analytic(sw), float)
        return self._central(lambda s: self.relperm.kro(s), sw, self.ds)

    # ═══════════════════════════════════════════════ kapilyar təzyiq
    def dpc_dsw(self, sw: np.ndarray) -> np.ndarray:
        if self.capillary is None:
            return np.zeros_like(np.asarray(sw, float))
        return np.asarray(self.capillary.dpcow_dsw(sw), float)

    # ═══════════════════════════════════════════════════════════ PVT
    def dbw_dp(self, pressure: np.ndarray) -> np.ndarray:
        if self.pvt is None:
            if self.fluids is None:
                return np.zeros_like(np.asarray(pressure, float))
            return self._static_fvf_derivative(
                pressure, self.fluids.water_fvf,
                self.fluids.water_compressibility)
        analytic = getattr(self.pvt, "water_fvf_derivative", None)
        if analytic is not None:
            return np.asarray(analytic(pressure), float)
        return self._central(self.pvt.water_fvf, pressure, self.dp)

    def dbo_dp(self, pressure: np.ndarray) -> np.ndarray:
        if self.pvt is None:
            if self.fluids is None:
                return np.zeros_like(np.asarray(pressure, float))
            return self._static_fvf_derivative(
                pressure, self.fluids.oil_fvf,
                self.fluids.oil_compressibility)
        analytic = getattr(self.pvt, "oil_fvf_derivative", None)
        if analytic is not None:
            return np.asarray(analytic(pressure), float)
        return self._central(self.pvt.oil_fvf, pressure, self.dp)

    def dmuw_dp(self, pressure: np.ndarray) -> np.ndarray:
        if self.pvt is None:
            return np.zeros_like(np.asarray(pressure, float))
        analytic = getattr(self.pvt, "water_viscosity_derivative", None)
        if analytic is not None:
            return np.asarray(analytic(pressure), float)
        return self._central(self.pvt.water_viscosity, pressure, self.dp)

    def dmuo_dp(self, pressure: np.ndarray) -> np.ndarray:
        if self.pvt is None:
            return np.zeros_like(np.asarray(pressure, float))
        analytic = getattr(self.pvt, "oil_viscosity_derivative", None)
        if analytic is not None:
            return np.asarray(analytic(pressure), float)
        return self._central(self.pvt.oil_viscosity, pressure, self.dp)

    # ═══════════════════════════════════ birləşmiş hərəkətlilik/B
    def water_transport_derivatives(self, pressure, sw, fluid):
        """M_w = krw / (μw·Bw) üçün (∂M/∂p, ∂M/∂Sw)."""
        krw = np.asarray(self.relperm.krw(sw), float)
        denominator = fluid.mu_w * fluid.bw
        d_denominator = (self.dmuw_dp(pressure) * fluid.bw
                         + fluid.mu_w * self.dbw_dp(pressure))
        d_dp = -krw * d_denominator / denominator ** 2
        d_dsw = self.dkrw_dsw(sw) / denominator
        return d_dp, d_dsw

    def oil_transport_derivatives(self, pressure, sw, fluid):
        """M_o = kro / (μo·Bo) üçün (∂M/∂p, ∂M/∂Sw)."""
        kro = np.asarray(self.relperm.kro(sw), float)
        denominator = fluid.mu_o * fluid.bo
        d_denominator = (self.dmuo_dp(pressure) * fluid.bo
                         + fluid.mu_o * self.dbo_dp(pressure))
        d_dp = -kro * d_denominator / denominator ** 2
        d_dsw = self.dkro_dsw(sw) / denominator
        return d_dp, d_dsw

    # ═══════════════════════════════════════════════════ köməkçi
    @staticmethod
    def _central(function, x: np.ndarray, step: float) -> np.ndarray:
        x = np.asarray(x, float)
        forward = np.asarray(function(x + step), float)
        backward = np.asarray(function(x - step), float)
        return (forward - backward) / (2.0 * step)
