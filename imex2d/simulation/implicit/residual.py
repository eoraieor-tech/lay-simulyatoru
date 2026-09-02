"""Qalıq (residual) vektoru — A6, mərhələ 1.

Hər hüceyrə üçün iki kütlə balansı tənliyi yazılır (su və neft).
Səth həcmi vahidlərində, çünki B_p ilə bölünmə orada kütləni saxlayır:

    R_p,c = (PV_c / Δt) · [ (S_p / B_p)ⁿ⁺¹ − (S_p / B_p)ⁿ ]
            − Σ_üzlər  T · (λ_p / B_p)_upstream · ΔΦ_p
            − q_p,c

    Φ_o = p            − ρ_o·g·D
    Φ_w = p − Pc(Sw)   − ρ_w·g·D

Tənliklərin işarəsi belədir ki, R = 0 tarazlıq deməkdir. Bütün
axın və quyu üzvləri IMPES mühərriki ilə EYNİ formullardan istifadə
edir — fərq yalnız ondadır ki, burada onlar YENİ zaman qatında
qiymətləndirilir (implicit), IMPES-də isə köhnə qatda idi.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

import numpy as np

from ...domain.reservoir_model import ReservoirModel
from ...domain.unit_conversions import PRESSURE_TO_PA, STANDARD_GRAVITY_M_S2
from ...domain.wells import ControlMode
from ...interfaces.providers import (ICapillaryPressureProvider, IPVTProvider,
                                     IRelativePermeabilityProvider)
from ..discretization import DiscretizedGrid
from .state import ReservoirState, VARIABLES_PER_CELL

#: Mərkəzləşdirilmiş sabitlər (bax `domain/unit_conversions.py`) — dəyərlər
#: ƏVVƏLKİ lokal literallarla BİT-BƏ-BİT eynidir, ədədi nəticə DƏYİŞMİR.
GRAVITY = STANDARD_GRAVITY_M_S2
PA_TO_BAR = 1.0 / PRESSURE_TO_PA["bar"]

WATER = 0
OIL = 1


@dataclass
class FluidState:
    """Verilmiş təzyiq/doyumluluqda flüid xassələri (hüceyrə üzrə)."""
    mu_w: np.ndarray
    mu_o: np.ndarray
    bw: np.ndarray
    bo: np.ndarray
    lam_w: np.ndarray
    lam_o: np.ndarray
    pc: Optional[np.ndarray] = None


@dataclass
class WellRates:
    """Qalıq hesablamasının yan məhsulu — hesabat üçün lazımdır."""
    water: np.ndarray          # hüceyrə üzrə, səth m3/gün (+ vurma, − hasilat)
    oil: np.ndarray
    per_well_water: dict
    per_well_oil: dict


class ResidualAssembler:
    """Verilmiş vəziyyət üçün qalıq vektorunu qurur."""

    def __init__(self, model: ReservoirModel, grid: DiscretizedGrid,
                 wells: list, relperm: IRelativePermeabilityProvider,
                 pvt: Optional[IPVTProvider] = None,
                 capillary: Optional[ICapillaryPressureProvider] = None):
        self.model = model
        self.reference_pressure = float(model.initial_conditions.datum_pressure)
        self.grid = grid
        self.wells = wells
        self.relperm = relperm
        self.pvt = pvt
        self.capillary = capillary

        self.connections = grid.connections
        self.transmissibility = grid.transmissibility
        self.pore_volume = grid.pore_volume
        self.ncell = model.ncell

        depths = model.geometry.cell_depths()
        self._depth_difference = (depths[self.connections.cell_a]
                                  - depths[self.connections.cell_b])
        self._has_gravity = bool(np.any(np.abs(self._depth_difference) > 1e-12))

        self._well_index = {c.well_name: c for c in wells}
        self._producer_names = sorted({c.well_name for c in wells
                                       if not c.is_injector})
        self._injector_names = sorted({c.well_name for c in wells
                                       if c.is_injector})

    # ═════════════════════════════════════════════════ flüid vəziyyəti
    def fluid_state(self, state: ReservoirState) -> FluidState:
        fluids = self.model.fluids
        pressure, sw = state.pressure, state.water_saturation

        if self.pvt is None:
            # PVT provider yoxdursa, sadə sıxılma modeli:
            #     B(p) = B_ref / (1 + c·(p − p_ref))
            # Bu, IMPES-in `ct` yanaşması ilə ekvivalentdir. Onsuz
            # implicit sxem sıxılmayan sistem həll edərdi və hasilat
            # yalnız vurulan həcmlə məhdudlaşardı — IMPES ilə fərq
            # yaradan məhz bu idi.
            shape = np.shape(pressure)
            delta = pressure - self.reference_pressure
            mu_w = np.full(shape, fluids.water_viscosity)
            mu_o = np.full(shape, fluids.oil_viscosity)
            bw = fluids.water_fvf / np.maximum(
                1.0 + fluids.water_compressibility * delta, 1e-6)
            bo = fluids.oil_fvf / np.maximum(
                1.0 + fluids.oil_compressibility * delta, 1e-6)
        else:
            mu_w = np.asarray(self.pvt.water_viscosity(pressure), float)
            mu_o = np.asarray(self.pvt.oil_viscosity(pressure), float)
            bw = np.asarray(self.pvt.water_fvf(pressure), float)
            bo = np.asarray(self.pvt.oil_fvf(pressure), float)

        pc = None if self.capillary is None else np.asarray(
            self.capillary.pcow(sw), float)

        return FluidState(
            mu_w=mu_w, mu_o=mu_o, bw=bw, bo=bo,
            lam_w=np.asarray(self.relperm.krw(sw), float) / mu_w,
            lam_o=np.asarray(self.relperm.kro(sw), float) / mu_o,
            pc=pc)

    # ═══════════════════════════════════════════════ faza potensialları
    def potentials(self, state: ReservoirState, fluid: FluidState):
        """Üzlər üzrə (ΔΦ_w, ΔΦ_o) — bar."""
        conn = self.connections
        dp = state.pressure[conn.cell_a] - state.pressure[conn.cell_b]
        d_phi_w = dp.copy()
        d_phi_o = dp.copy()

        if self._has_gravity:
            rho_w = self.model.fluids.water_density / np.maximum(fluid.bw, 1e-9)
            rho_o = self.model.fluids.oil_density / np.maximum(fluid.bo, 1e-9)
            head = GRAVITY * self._depth_difference * PA_TO_BAR
            d_phi_w -= 0.5 * (rho_w[conn.cell_a] + rho_w[conn.cell_b]) * head
            d_phi_o -= 0.5 * (rho_o[conn.cell_a] + rho_o[conn.cell_b]) * head

        if fluid.pc is not None:
            d_phi_w -= fluid.pc[conn.cell_a] - fluid.pc[conn.cell_b]

        return d_phi_w, d_phi_o

    # ═══════════════════════════════════════════════════ üzlər üzrə axın
    def face_fluxes(self, state: ReservoirState, fluid: FluidState):
        """A → B istiqamətində səth həcmi axını (m3/gün)."""
        conn = self.connections
        d_phi_w, d_phi_o = self.potentials(state, fluid)

        up_w = np.where(d_phi_w >= 0, conn.cell_a, conn.cell_b)
        up_o = np.where(d_phi_o >= 0, conn.cell_a, conn.cell_b)

        water = self.transmissibility * (fluid.lam_w[up_w] / fluid.bw[up_w]) * d_phi_w
        oil = self.transmissibility * (fluid.lam_o[up_o] / fluid.bo[up_o]) * d_phi_o
        return water, oil

    def net_influx(self, state: ReservoirState, fluid: FluidState):
        """Hüceyrəyə daxil olan xalis axın (m3/gün, səth)."""
        conn = self.connections
        water_flux, oil_flux = self.face_fluxes(state, fluid)

        water = np.zeros(self.ncell)
        oil = np.zeros(self.ncell)
        np.add.at(water, conn.cell_a, -water_flux)
        np.add.at(water, conn.cell_b, +water_flux)
        np.add.at(oil, conn.cell_a, -oil_flux)
        np.add.at(oil, conn.cell_b, +oil_flux)
        return water, oil

    # ══════════════════════════════════════════════════════════ quyular
    def well_rates(self, state: ReservoirState, fluid: FluidState) -> WellRates:
        """Quyu debitləri, səth həcmi. Müsbət = laya daxil olur."""
        water = np.zeros(self.ncell)
        oil = np.zeros(self.ncell)
        per_well_water = {name: 0.0 for name in
                          self._producer_names + self._injector_names}
        per_well_oil = dict(per_well_water)

        endpoint_mobility = self.relperm.endpoint_water_mobility(1.0)

        for connection in self.wells:
            cell = connection.cell
            if connection.is_injector:
                mobility = endpoint_mobility / fluid.mu_w[cell]
                if connection.mode is ControlMode.BHP:
                    rate = (connection.well_index * mobility
                            * (connection.target - state.pressure[cell]))
                else:
                    rate = abs(connection.target)
                rate = max(rate, 0.0) / fluid.bw[cell]
                water[cell] += rate
                per_well_water[connection.well_name] += rate
            else:
                lam_w = fluid.lam_w[cell]
                lam_o = fluid.lam_o[cell]
                if connection.mode is ControlMode.BHP:
                    drawdown = connection.target - state.pressure[cell]
                    qw = connection.well_index * lam_w * drawdown
                    qo = connection.well_index * lam_o * drawdown
                else:
                    total = -abs(connection.target)
                    fraction = lam_w / max(lam_w + lam_o, 1e-30)
                    qw, qo = total * fraction, total * (1.0 - fraction)
                qw = min(qw, 0.0) / fluid.bw[cell]
                qo = min(qo, 0.0) / fluid.bo[cell]
                water[cell] += qw
                oil[cell] += qo
                per_well_water[connection.well_name] += qw
                per_well_oil[connection.well_name] += qo

        return WellRates(water, oil, per_well_water, per_well_oil)

    # ═══════════════════════════════════════════════════════ akkumulyasiya
    def pore_volume_at(self, pressure: np.ndarray) -> np.ndarray:
        """Süxur sıxılması: PV(p) = PV_ref · [1 + c_r·(p − p_ref)]."""
        compressibility = self.model.rock.compressibility
        if compressibility <= 0.0:
            return self.pore_volume
        factor = 1.0 + compressibility * (pressure - self.reference_pressure)
        return self.pore_volume * np.maximum(factor, 1e-6)

    def accumulation(self, state: ReservoirState, fluid: FluidState):
        """PV(p) · S_p / B_p — səth həcmi kimi hüceyrədəki flüid miqdarı."""
        sw = state.water_saturation
        pore_volume = self.pore_volume_at(state.pressure)
        water = pore_volume * sw / fluid.bw
        oil = pore_volume * (1.0 - sw) / fluid.bo
        return water, oil

    # ═══════════════════════════════════════════════════════════ qalıq
    def residual(self, state: ReservoirState, previous: ReservoirState,
                 dt: float, previous_fluid: Optional[FluidState] = None):
        """R(x) vektoru — ölçüsü 2N.

        Qaytarır: (residual, fluid_state, well_rates)
        """
        fluid = self.fluid_state(state)
        old_fluid = previous_fluid or self.fluid_state(previous)

        water_new, oil_new = self.accumulation(state, fluid)
        water_old, oil_old = self.accumulation(previous, old_fluid)

        influx_water, influx_oil = self.net_influx(state, fluid)
        rates = self.well_rates(state, fluid)

        residual_water = (water_new - water_old) / dt - influx_water - rates.water
        residual_oil = (oil_new - oil_old) / dt - influx_oil - rates.oil

        vector = np.empty(self.ncell * VARIABLES_PER_CELL)
        vector[WATER::VARIABLES_PER_CELL] = residual_water
        vector[OIL::VARIABLES_PER_CELL] = residual_oil
        return vector, fluid, rates

    # ═══════════════════════════════════════════════ konvergensiya ölçüsü
    def scaled_residual_norm(self, residual: np.ndarray, dt: float) -> float:
        """Normallaşdırılmış qalıq — CNV tipli ölçü.

        Xam qalıq m3/gün vahidindədir və hüceyrə həcmindən asılıdır.
        Məsamə həcminə bölmək onu ölçüsüz edir və müxtəlif grid
        ölçülərində eyni konvergensiya həddini işlətməyə imkan verir.
        """
        scale = np.maximum(self.pore_volume / dt, 1e-30)
        water = np.abs(residual[WATER::VARIABLES_PER_CELL]) / scale
        oil = np.abs(residual[OIL::VARIABLES_PER_CELL]) / scale
        return float(max(water.max(), oil.max()))

    def material_balance_error(self, residual: np.ndarray, dt: float) -> tuple:
        """Qalıqların cəmi = ümumi kütlə balansı səhvi (m3, səth).

        Fiziki mənası: bütün hüceyrələr üzrə toplayanda daxili axınlar
        bir-birini yeyir və yalnız [dəyişiklik − quyu debiti] qalır.
        """
        water = float(np.sum(residual[WATER::VARIABLES_PER_CELL]) * dt)
        oil = float(np.sum(residual[OIL::VARIABLES_PER_CELL]) * dt)
        return water, oil
