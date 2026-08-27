"""Analitik həllər — validasiya üçün. Mövcud koddan köçürülüb.

Simulyasiya mühərrikindən ayrı saxlanılır: bu, verifikasiya alətidir,
hesablama zəncirinin hissəsi deyil.
"""

from __future__ import annotations
from dataclasses import dataclass

import numpy as np

from ..domain.scal import CoreyParameters


@dataclass
class BuckleyLeverettSolution:
    distance: np.ndarray
    water_saturation: np.ndarray
    shock_saturation: float
    front_position: float


def buckley_leverett(scal: CoreyParameters, mu_w: float, mu_o: float,
                     porosity: float, total_rate: float, area: float,
                     time: float, sw_initial: float = None) -> BuckleyLeverettSolution:
    swi = scal.swc if sw_initial is None else sw_initial
    s = np.linspace(swi + 1e-6, 1.0 - scal.sor - 1e-6, 2000)

    lw = scal.krw(s) / mu_w
    lo = scal.kro(s) / mu_o
    f = lw / np.maximum(lw + lo, 1e-30)
    lw_i = scal.krw(np.array([swi]))[0] / mu_w
    lo_i = scal.kro(np.array([swi]))[0] / mu_o
    f_i = lw_i / max(lw_i + lo_i, 1e-30)

    chord = (f - f_i) / (s - swi)
    k = int(np.argmax(chord))
    sw_shock, slope_shock = s[k], chord[k]

    dfds = np.gradient(f, s)
    velocity = total_rate / (area * porosity)
    x_front = velocity * slope_shock * time

    s_prof = s[s >= sw_shock]
    x_prof = velocity * dfds[s >= sw_shock] * time
    x_prof = np.maximum.accumulate(x_prof[::-1])[::-1]

    x = np.concatenate([[0.0], x_prof, [x_front, x_front * 1.6]])
    sw = np.concatenate([[1.0 - scal.sor], s_prof, [sw_shock, swi]])
    order = np.argsort(x)
    return BuckleyLeverettSolution(x[order], sw[order], float(sw_shock), float(x_front))
