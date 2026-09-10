"""BlackOilPVTProvider — IPVTProvider-in ilk implementasiyası.

Yeganə hesablama üsulu: cədvəl üzrə xətti interpolyasiya (np.interp).
Cədvəldən kənarda sərhəd dəyəri saxlanılır (np.interp-in defolt davranışı) —
bu, ekstrapolyasiyadan daha təhlükəsizdir.

Sıxılma cədvəldən ədədi törəmə ilə alınır və bir dəfə əvvəlcədən
hesablanır ki, hər zaman addımında yenidən hesablanmasın.
"""

from __future__ import annotations
from typing import Optional

import numpy as np

from ...domain.pvt import PVTTable
from ...interfaces.providers import IPVTProvider


class BlackOilPVTProvider(IPVTProvider):

    def __init__(self, table: PVTTable):
        issues = table.validate()
        if issues:
            raise ValueError("PVT cədvəli yararsızdır: " + "; ".join(issues))
        self.table = table
        self._co = table.compressibility("oil_fvf")
        self._cw = table.compressibility("water_fvf")
        self._cr = table.rock_compressibility
        # Törəmələr bir dəfə hesablanır. Cədvəl PARÇALI XƏTTİ olduğuna
        # görə törəmə hər intervalda sabitdir — `np.gradient`-in verdiyi
        # hamar qiymət deyil, məhz interval meyli. Bu seçim ölçmə ilə
        # təsdiqləndi: hamar törəmə Nyutonun iterasiya sayını
        # dəyişmir, lakin Jakobianı sonlu fərqdən 10⁶ dəfə uzaqlaşdırır.
        pressure = table.pressure
        columns = [("oil_fvf", table.oil_fvf), ("water_fvf", table.water_fvf),
                  ("oil_viscosity", table.oil_viscosity),
                  ("water_viscosity", table.water_viscosity),
                  ("solution_gor", table.solution_gor)]
        if table.has_gas_phase:
            columns += [("gas_fvf", table.gas_fvf),
                       ("gas_viscosity", table.gas_viscosity)]
        self._slopes = {name: np.diff(values) / np.diff(pressure)
                        for name, values in columns}

    # ------------------------------------------------------ interpolyasiya
    def _interp(self, values: np.ndarray, pressure) -> np.ndarray:
        return np.interp(np.asarray(pressure, float), self.table.pressure, values)

    def oil_fvf(self, pressure, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._interp(self.table.oil_fvf, pressure)

    def oil_viscosity(self, pressure, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._interp(self.table.oil_viscosity, pressure)

    def water_fvf(self, pressure, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._interp(self.table.water_fvf, pressure)

    def water_viscosity(self, pressure, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._interp(self.table.water_viscosity, pressure)

    def solution_gor(self, pressure, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._interp(self.table.solution_gor, pressure)

    def total_compressibility(self, pressure, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        sw = np.asarray(sw, float)
        co = self._interp(self._co, pressure)
        cw = self._interp(self._cw, pressure)
        return cw * sw + co * (1.0 - sw) + self._cr

    def bubble_point(self, region: Optional[np.ndarray] = None) -> float:
        return self.table.bubble_point

    def has_gas_phase(self, region: Optional[np.ndarray] = None) -> bool:
        return self.table.has_gas_phase

    def gas_fvf(self, pressure, region: Optional[np.ndarray] = None) -> np.ndarray:
        if not self.table.has_gas_phase:
            raise NotImplementedError(
                "Bu PVT cədvəlində qaz xassələri yoxdur "
                "(build_pvt_table(..., include_gas=True) işlədin).")
        return self._interp(self.table.gas_fvf, pressure)

    def gas_viscosity(self, pressure, region: Optional[np.ndarray] = None) -> np.ndarray:
        if not self.table.has_gas_phase:
            raise NotImplementedError(
                "Bu PVT cədvəlində qaz xassələri yoxdur "
                "(build_pvt_table(..., include_gas=True) işlədin).")
        return self._interp(self.table.gas_viscosity, pressure)

    # ─────────────────────────────────────────── analitik törəmələr
    def _slope(self, name: str, pressure) -> np.ndarray:
        """Parçalı xətti cədvəlin dəqiq törəməsi.

        Cədvəldən kənarda sıfırdır, çünki `np.interp` orada sərhəd
        dəyərini saxlayır (funksiya sabitdir).
        """
        pressure = np.atleast_1d(np.asarray(pressure, float))
        nodes = self.table.pressure
        index = np.clip(np.searchsorted(nodes, pressure, side="right") - 1,
                        0, nodes.size - 2)
        slopes = self._slopes[name][index]
        outside = (pressure < nodes[0]) | (pressure > nodes[-1])
        return np.where(outside, 0.0, slopes)

    def oil_fvf_derivative(self, pressure, region=None) -> np.ndarray:
        return self._slope("oil_fvf", pressure)

    def water_fvf_derivative(self, pressure, region=None) -> np.ndarray:
        return self._slope("water_fvf", pressure)

    def oil_viscosity_derivative(self, pressure, region=None) -> np.ndarray:
        return self._slope("oil_viscosity", pressure)

    def water_viscosity_derivative(self, pressure, region=None) -> np.ndarray:
        return self._slope("water_viscosity", pressure)

    def gas_fvf_derivative(self, pressure, region=None) -> np.ndarray:
        return self._slope("gas_fvf", pressure)

    def gas_viscosity_derivative(self, pressure, region=None) -> np.ndarray:
        return self._slope("gas_viscosity", pressure)

    def solution_gor_derivative(self, pressure, region=None) -> np.ndarray:
        """dRs_sat/dp — Jakobianda doymuş hüceyrələr üçün lazımdır (A7/6c).

        Doyma təzyiqindən yuxarıda Rs sabitdir, ona görə bu, sıfıra
        düşür — cədvəlin özündə bu sabitlik artıq mövcuddur, əlavə
        şərtə ehtiyac yoxdur.
        """
        return self._slope("solution_gor", pressure)
