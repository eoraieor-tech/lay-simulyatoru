"""Black-oil korrelyasiyaları — PVT CƏDVƏLİ GENERATORU.

Bu modul mühərrikə heç vaxt birbaşa qoşulmur. Onun yeganə işi
laboratoriya məlumatı olmayanda `PVTTable` istehsal etməkdir.
Bu ayrılıq sayəsində mühərrik korrelyasiyaların mövcudluğundan
xəbərsizdir və həmişə cədvəllə işləyir.

İstifadə olunan korrelyasiyalar (klassik, dəyişdirilməyib):
    Standing (1947)        — Pb və Rs
    Vazquez-Beggs (1980)   — Bo, doyma altı neftin sıxılması
    Beggs-Robinson (1975)  — ölü və doymuş neftin lözlüyü
    Meehan / McCain        — su üçün Bw və μw
    Sutton (1985)          — psevdo-kritik Tpc/Ppc (qaz sıxlığından)
    Beggs-Brill (1973)     — Z-faktoru (Standing-Katz əyrisinin
                             açıq (iterasiyasız) approksimasiyası)
    Lee-Gonzalez-Eakin (1966) — qaz lözlüyü
"""

from __future__ import annotations

import numpy as np

from ...domain.pvt import PVTTable

BAR_TO_PSI = 14.5037744
SM3M3_TO_SCFSTB = 5.61458
C_TO_F = lambda c: c * 9.0 / 5.0 + 32.0
C_TO_RANKINE = lambda c: (c * 9.0 / 5.0 + 32.0) + 459.67

# standart şərait: 1 atm, 60°F (Eclipse METRIC-in özü ilə uyğun)
STANDARD_PRESSURE_BAR = 1.01325
STANDARD_TEMPERATURE_RANKINE = 519.67


def api_to_specific_gravity(api: float) -> float:
    return 141.5 / (api + 131.5)


def standing_solution_gor(pressure_bar, api: float, gas_gravity: float,
                          temperature_c: float) -> np.ndarray:
    """Rs(p), sm3/sm3. Standing (1947)."""
    p_psi = np.asarray(pressure_bar, float) * BAR_TO_PSI
    t_f = C_TO_F(temperature_c)
    x = 0.0125 * api - 0.00091 * t_f
    rs_scf_stb = gas_gravity * ((p_psi / 18.2 + 1.4) * 10.0 ** x) ** 1.2048
    return rs_scf_stb / SM3M3_TO_SCFSTB


def standing_bubble_point(rs_sm3_sm3: float, api: float, gas_gravity: float,
                          temperature_c: float) -> float:
    """Pb, bar. Standing korrelyasiyasının tərs formu."""
    rs = rs_sm3_sm3 * SM3M3_TO_SCFSTB
    t_f = C_TO_F(temperature_c)
    x = 0.0125 * api - 0.00091 * t_f
    p_psi = 18.2 * ((rs / gas_gravity) ** 0.83 * 10.0 ** (-x) - 1.4)
    return max(p_psi / BAR_TO_PSI, 1.0)


def vazquez_beggs_oil_fvf(rs_sm3_sm3, api: float, gas_gravity: float,
                          temperature_c: float) -> np.ndarray:
    """Bo(Rs) doyma nöqtəsinə qədər. Vazquez-Beggs (1980)."""
    rs = np.asarray(rs_sm3_sm3, float) * SM3M3_TO_SCFSTB
    t_f = C_TO_F(temperature_c)
    if api <= 30.0:
        c1, c2, c3 = 4.677e-4, 1.751e-5, -1.811e-8
    else:
        c1, c2, c3 = 4.670e-4, 1.100e-5, 1.337e-9
    return 1.0 + c1 * rs + (t_f - 60.0) * (api / gas_gravity) * (c2 + c3 * rs)


def vazquez_beggs_undersaturated_compressibility(rs_sm3_sm3: float, api: float,
                                                 gas_gravity: float,
                                                 temperature_c: float,
                                                 pressure_bar: float) -> float:
    """co, 1/bar — doyma təzyiqindən yuxarı."""
    rs = rs_sm3_sm3 * SM3M3_TO_SCFSTB
    t_f = C_TO_F(temperature_c)
    p_psi = max(pressure_bar * BAR_TO_PSI, 1.0)
    co_psi = (-1433.0 + 5.0 * rs + 17.2 * t_f - 1180.0 * gas_gravity
              + 12.61 * api) / (1e5 * p_psi)
    return max(co_psi, 1e-9) * BAR_TO_PSI


def beggs_robinson_dead_oil_viscosity(api: float, temperature_c: float) -> float:
    """μod, cP."""
    t_f = C_TO_F(temperature_c)
    z = 3.0324 - 0.02023 * api
    y = 10.0 ** z
    x = y * t_f ** -1.163
    return 10.0 ** x - 1.0


def beggs_robinson_saturated_viscosity(mu_dead: float, rs_sm3_sm3) -> np.ndarray:
    """μo(Rs), cP — doymuş neft."""
    rs = np.asarray(rs_sm3_sm3, float) * SM3M3_TO_SCFSTB
    a = 10.715 * (rs + 100.0) ** -0.515
    b = 5.44 * (rs + 150.0) ** -0.338
    return a * mu_dead ** b


def water_fvf(pressure_bar, temperature_c: float) -> np.ndarray:
    """Bw(p, T) — McCain tipli sadə korrelyasiya."""
    p_psi = np.asarray(pressure_bar, float) * BAR_TO_PSI
    t_f = C_TO_F(temperature_c)
    d_vwt = -1.0001e-2 + 1.33391e-4 * t_f + 5.50654e-7 * t_f ** 2
    d_vwp = (-1.95301e-9 * p_psi * t_f - 1.72834e-13 * p_psi ** 2 * t_f
             - 3.58922e-7 * p_psi - 2.25341e-10 * p_psi ** 2)
    return (1.0 + d_vwt) * (1.0 + d_vwp)


def water_viscosity(temperature_c: float, salinity_ppm: float = 0.0) -> float:
    """μw(T), cP — Meehan korrelyasiyası."""
    t_f = C_TO_F(temperature_c)
    s = salinity_ppm / 1e4
    a = 109.574 - 8.40564 * s + 0.313314 * s ** 2 + 8.72213e-3 * s ** 3
    b = -1.12166 + 2.63951e-2 * s - 6.79461e-4 * s ** 2 - 5.47119e-5 * s ** 3
    return max(a * t_f ** b, 0.05)


def sutton_pseudo_critical(gas_gravity: float) -> tuple:
    """Tpc, Ppc (Rankine, psia) — qaz nisbi sıxlığından.

    Sutton (1985) korrelyasiyası: laboratoriya kompozisiya analizi
    olmadan, yalnız qazın havaya nisbətən sıxlığından (γg) psevdo-
    kritik temperatur və təzyiqi təxmin edir. Quru qaz üçün sənayə
    standartıdır.
    """
    t_pc = 169.2 + 349.5 * gas_gravity - 74.0 * gas_gravity ** 2
    p_pc = 756.8 - 131.0 * gas_gravity - 3.6 * gas_gravity ** 2
    return t_pc, p_pc


def beggs_brill_z_factor(pressure_bar, gas_gravity: float,
                         temperature_c: float) -> np.ndarray:
    """Z(p) — Standing-Katz əyrisinin açıq approksimasiyası.

    Dranchuk-Abou-Kassem kimi iterativ tənliklərdən fərqli olaraq
    Beggs-Brill (1973) BİRBAŞA (qapalı formada) hesablanır — vektor
    əməliyyatlarına və Nyuton-daxili çağırışlara ideal uyğundur, çünki
    hər addımda min hüceyrə üçün iterasiyasız qiymətləndirilir.
    """
    p_psi = np.atleast_1d(np.asarray(pressure_bar, float)) * BAR_TO_PSI
    t_pc, p_pc = sutton_pseudo_critical(gas_gravity)
    t_pr = C_TO_RANKINE(temperature_c) / t_pc
    p_pr = p_psi / p_pc

    a = 1.39 * np.sqrt(np.maximum(t_pr - 0.92, 1e-6)) - 0.36 * t_pr - 0.101
    b = ((0.62 - 0.23 * t_pr) * p_pr
         + (0.066 / np.maximum(t_pr - 0.86, 1e-6) - 0.037) * p_pr ** 2
         + 0.32 * p_pr ** 6 / 10.0 ** (9.0 * (t_pr - 1.0)))
    c = 0.132 - 0.32 * np.log10(np.maximum(t_pr, 1e-6))
    d = 10.0 ** (0.3106 - 0.49 * t_pr + 0.1824 * t_pr ** 2)

    z = a + (1.0 - a) / np.exp(np.clip(b, -50, 50)) + c * p_pr ** d
    return np.maximum(z, 0.2)          # fiziki alt hədd (çox aşağı P-də)


def gas_fvf(pressure_bar, gas_gravity: float,
           temperature_c: float) -> np.ndarray:
    """Bg(p), rm³/sm³ — həqiqi qaz qanunu.

        Bg = (Z·T/P) / (Zsc·Tsc/Psc),   Zsc ≈ 1

    Neft/su FVF-dən fərqli olaraq Bg təzyiqlə KƏSKİN azalır (qaz
    sıxıla bilən) — aşağı təzyiqdə 100-dən çox ola bilər.
    """
    pressure_bar = np.maximum(np.atleast_1d(np.asarray(pressure_bar, float)),
                              1e-3)
    z = beggs_brill_z_factor(pressure_bar, gas_gravity, temperature_c)
    t_r = C_TO_RANKINE(temperature_c)
    p_psi = pressure_bar * BAR_TO_PSI
    return 0.02827 * z * t_r / p_psi * (
        STANDARD_PRESSURE_BAR / BAR_TO_PSI * BAR_TO_PSI
        / STANDARD_TEMPERATURE_RANKINE) * STANDARD_TEMPERATURE_RANKINE         / (STANDARD_PRESSURE_BAR * BAR_TO_PSI) * (STANDARD_PRESSURE_BAR * BAR_TO_PSI)


def gas_viscosity(pressure_bar, gas_gravity: float,
                  temperature_c: float) -> np.ndarray:
    """μg(p), cP — Lee, Gonzalez, Eakin (1966).

    Neft/su lözlüyündən fərqli olaraq qaz lözlüyü təzyiqlə ARTIR
    (sıxlıq artdıqca molekullararası toqquşma çoxalır) — neft/su-nun
    əks istiqamətdəki davranışı ilə müqayisə üçün faydalı yoxlamadır.
    """
    pressure_bar = np.atleast_1d(np.asarray(pressure_bar, float))
    z = beggs_brill_z_factor(pressure_bar, gas_gravity, temperature_c)
    t_r = C_TO_RANKINE(temperature_c)
    p_psi = pressure_bar * BAR_TO_PSI
    molar_mass = 28.97 * gas_gravity                 # lb/lb-mol

    density = 1.4935e-3 * p_psi * molar_mass / np.maximum(z * t_r, 1e-6)  # g/cm3

    k = (9.4 + 0.02 * molar_mass) * t_r ** 1.5 / (209.0 + 19.0 * molar_mass + t_r)
    x = 3.5 + 986.0 / t_r + 0.01 * molar_mass
    y = 2.4 - 0.2 * x
    return k * 1e-4 * np.exp(x * np.maximum(density, 1e-6) ** y)


def build_pvt_table(api: float = 32.0,
                    gas_gravity: float = 0.75,
                    temperature_c: float = 70.0,
                    salinity_ppm: float = 30000.0,
                    pressure_min: float = 1.0,
                    pressure_max: float = 400.0,
                    n_points: int = 40,
                    bubble_point_bar: float = None,
                    rock_compressibility: float = 4.5e-5,
                    include_gas: bool = False) -> PVTTable:
    """Korrelyasiyalardan tam PVT cədvəli qurur.

    `include_gas=True` — Bg və μg də hesablanır (A7, üç fazalı model
    üçün). Defolt `False`-dur ki, mövcud iki fazalı testlər və
    ssenarilər dəyişməsin.
    """
    pressure = np.linspace(pressure_min, pressure_max, int(n_points))

    rs_saturated = standing_solution_gor(pressure, api, gas_gravity, temperature_c)
    if bubble_point_bar is None:
        bubble_point_bar = float(pressure_max * 0.6)
    rs_at_pb = float(standing_solution_gor(np.array([bubble_point_bar]), api,
                                           gas_gravity, temperature_c)[0])

    # Pb-dən yuxarı Rs sabit qalır (qaz artıq həll olub)
    rs = np.where(pressure < bubble_point_bar, rs_saturated, rs_at_pb)

    bo = vazquez_beggs_oil_fvf(rs, api, gas_gravity, temperature_c)
    bo_at_pb = float(vazquez_beggs_oil_fvf(np.array([rs_at_pb]), api,
                                           gas_gravity, temperature_c)[0])
    co = vazquez_beggs_undersaturated_compressibility(
        rs_at_pb, api, gas_gravity, temperature_c, bubble_point_bar)
    above = pressure >= bubble_point_bar
    bo[above] = bo_at_pb * np.exp(-co * (pressure[above] - bubble_point_bar))

    mu_dead = beggs_robinson_dead_oil_viscosity(api, temperature_c)
    mu_oil = beggs_robinson_saturated_viscosity(mu_dead, rs)
    mu_at_pb = float(beggs_robinson_saturated_viscosity(
        mu_dead, np.array([rs_at_pb]))[0])
    # Pb-dən yuxarı lözlük təzyiqlə artır (Vazquez-Beggs)
    ratio = np.maximum(pressure[above] / max(bubble_point_bar, 1e-9), 1.0)
    mu_oil[above] = mu_at_pb * ratio ** 0.278

    return PVTTable(
        pressure=pressure,
        oil_fvf=bo,
        oil_viscosity=mu_oil,
        solution_gor=rs,
        water_fvf=water_fvf(pressure, temperature_c),
        water_viscosity=np.full_like(pressure,
                                     water_viscosity(temperature_c, salinity_ppm)),
        bubble_point=float(bubble_point_bar),
        rock_compressibility=rock_compressibility,
        source=f"correlation(API={api:g}, γg={gas_gravity:g}, T={temperature_c:g}°C)",
        gas_fvf=(gas_fvf(pressure, gas_gravity, temperature_c)
                 if include_gas else None),
        gas_viscosity=(gas_viscosity(pressure, gas_gravity, temperature_c)
                      if include_gas else None),
    )
