"""EquilibriumInitializationProvider — A3.

Hidrostatik tarazlıqdan ilkin təzyiq və doyumluluq sahələrini qurur:

    p(z) = p_datum + ρ·g·(z − z_datum) / 1e5      [bar]

Kontaktdan yuxarı neft zonası (ρ = ρo), aşağı su zonası (ρ = ρw).

A4-dən sonra: kapilyar provider verilibsə, kontaktdan yuxarıda
kapilyar-cazibə tarazlığından **keçid zonası** qurulur —

    Pc(D) = (ρw − ρo)·g·(D_owc − D) / 1e5,   Sw = Pc⁻¹(Pc(D))

Provider verilmədikdə kontakt kəskin qalır (köhnə davranış).

A7-dən sonra: `gas_oil_contact` verilib PVT provider `has_gas_phase()`
qaytarırsa, GOC-dan yuxarı QAZ PAPAĞI (gas cap) qurulur —
üç zonalı tarazlıq:

    GOC-dan yuxarı   qaz zonası     Sw = Swc,  Sg = 1 − Swc
    GOC — OWC arası  neft zonası    Sg = 0,    Sw əvvəlki kimi
    OWC-dən aşağı    su zonası      Sw = 1 − Sor

Sadələşdirmə: qaz-neft kontaktında kapilyar keçid zonası hələ
qurulmur (kəskin sərhəd) — su-neft kontaktındakı A4 keçid zonası
kimi hamar keçid gələcək təkmilləşdirmədir. Təzyiq profilinə qazın
sütun çəkisi daxildir (adətən çox kiçikdir, çünki ρ_qaz ≪ ρ_neft).

Sıxlıqlar lay şəraitinə çevrilir: ρ_lay = ρ_səth / B. PVT provider
verilibsə B(p) ondan, verilməyibsə modelin sabit dəyərindən alınır.
"""

from __future__ import annotations
from typing import Optional

import numpy as np

from ...domain.reservoir_model import ReservoirModel
from ...interfaces.providers import (ICapillaryPressureProvider,
                                     IInitializationProvider, InitialState,
                                     IPVTProvider)

GRAVITY = 9.80665          # m/s2
PA_TO_BAR = 1.0e-5


class EquilibriumInitializationProvider(IInitializationProvider):

    def __init__(self, pvt: Optional[IPVTProvider] = None,
                 capillary: Optional[ICapillaryPressureProvider] = None,
                 max_iterations: int = 20, tolerance: float = 1e-6):
        self.pvt = pvt
        self.capillary = capillary
        self.max_iterations = int(max_iterations)
        self.tolerance = float(tolerance)

    # ------------------------------------------------------------ public
    def initialize(self, model: ReservoirModel) -> InitialState:
        ic = model.initial_conditions
        depths = model.geometry.cell_depths()
        scal = model.scal_parameters

        contact = ic.oil_water_contact
        if contact is None:
            contact = float(np.max(depths)) + 1.0   # bütün lay neft zonasında

        has_gas_cap = (ic.gas_oil_contact is not None and self.pvt is not None
                       and self.pvt.has_gas_phase())
        goc = ic.gas_oil_contact if has_gas_cap else None

        below = depths >= contact
        pressure = self._pressure_profile(model, depths, contact, goc)

        if self.capillary is None:
            sw = np.where(below, 1.0 - scal.sor, scal.swc)
        else:
            sw = self._transition_zone(model, depths, contact, pressure, below)

        if not has_gas_cap:
            sw = np.clip(sw, scal.swc, 1.0 - scal.sor)
            return InitialState(pressure=pressure, water_saturation=sw)

        # qaz papağı: GOC-dan yuxarı Sw -> Swc, Sg -> 1 − Swc
        in_gas_cap = depths < goc
        sw = np.where(in_gas_cap, scal.swc, sw)
        sg = np.where(in_gas_cap, 1.0 - scal.swc, 0.0)
        sw = np.clip(sw, scal.swc, 1.0 - scal.sor)
        sg = np.clip(sg, 0.0, 1.0 - scal.swc)
        return InitialState(pressure=pressure, water_saturation=sw,
                            gas_saturation=sg)

    def _transition_zone(self, model, depths, contact, pressure, below):
        """Kapilyar-cazibə tarazlığı: kontaktdan yuxarı hamar keçid."""
        scal = model.scal_parameters
        rho_oil = self._phase_density(model, "oil", pressure)
        rho_water = self._phase_density(model, "water", pressure)
        height = np.clip(contact - depths, 0.0, None)          # kontaktdan yuxarı, m
        pc = (rho_water - rho_oil) * GRAVITY * height * PA_TO_BAR

        if hasattr(self.capillary, "saturation_from_pc"):
            sw = self.capillary.saturation_from_pc(pc)
        else:
            sw = self._invert_numerically(pc, scal)
        return np.where(below, 1.0 - scal.sor, sw)

    def _invert_numerically(self, pc, scal):
        """Provider tərs funksiya vermirsə, cədvəl üzrə axtarış."""
        grid = np.linspace(scal.swc, 1.0 - scal.sor, 512)
        curve = np.asarray(self.capillary.pcow(grid), float)
        order = np.argsort(curve)
        return np.interp(np.asarray(pc, float), curve[order], grid[order])

    # ----------------------------------------------------------- internal
    def _phase_density(self, model: ReservoirModel, phase: str,
                       pressure) -> np.ndarray:
        """Lay şəraitində sıxlıq, kg/m3."""
        fluids = model.fluids
        if phase == "oil":
            surface, fvf = fluids.oil_density, (
                self.pvt.oil_fvf(pressure) if self.pvt else fluids.oil_fvf)
        elif phase == "gas":
            surface = fluids.gas_density
            fvf = (self.pvt.gas_fvf(pressure) if self.pvt
                   and self.pvt.has_gas_phase() else 1.0)
        else:
            surface, fvf = fluids.water_density, (
                self.pvt.water_fvf(pressure) if self.pvt else fluids.water_fvf)
        fvf = np.full(np.shape(pressure), fvf, dtype=float) \
            if np.ndim(fvf) == 0 else np.asarray(fvf, float)
        return surface / np.maximum(fvf, 1e-9)

    def _pressure_profile(self, model: ReservoirModel, depths: np.ndarray,
                          contact: float,
                          goc: Optional[float] = None) -> np.ndarray:
        """Hidrostatik təzyiq. Sıxlıq təzyiqdən asılı olduğu üçün iterasiya.

        `goc` verilibsə, datum ilə GOC arasında qaz sütununun çəkisi
        də əlavə olunur (adətən çox kiçikdir, çünki ρ_qaz ≪ ρ_neft,
        lakin dərin qaz papaqlarında hiss oluna bilər).
        """
        ic = model.initial_conditions
        datum_depth = ic.datum_depth if ic.datum_depth else float(np.min(depths))
        pressure = np.full(depths.size, float(ic.datum_pressure))
        gas_top = min(goc, datum_depth) if goc is not None else datum_depth

        for _ in range(self.max_iterations):
            rho_oil = self._phase_density(model, "oil", pressure)
            rho_water = self._phase_density(model, "water", pressure)

            oil_column = (np.clip(depths, gas_top, contact)
                         - max(gas_top, min(datum_depth, contact)))
            water_column = np.clip(depths - contact, 0.0, None)

            new_pressure = (float(ic.datum_pressure)
                            + rho_oil * GRAVITY * oil_column * PA_TO_BAR
                            + rho_water * GRAVITY * water_column * PA_TO_BAR)

            if goc is not None:
                rho_gas = self._phase_density(model, "gas", pressure)
                gas_column = (np.clip(depths, gas_top, goc) - gas_top)
                new_pressure = new_pressure + (rho_gas * GRAVITY * gas_column
                                               * PA_TO_BAR)

            if np.max(np.abs(new_pressure - pressure)) < self.tolerance:
                pressure = new_pressure
                break
            pressure = new_pressure
        return pressure
