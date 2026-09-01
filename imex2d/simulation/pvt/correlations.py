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
from ...domain.unit_conversions import convert_temperature

BAR_TO_PSI = 14.5037744
SM3M3_TO_SCFSTB = 5.61458
C_TO_F = lambda c: c * 9.0 / 5.0 + 32.0


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


def build_pvt_table(api: float = 32.0,
                    gas_gravity: float = 0.75,
                    temperature_c: float = 70.0,
                    salinity_ppm: float = 30000.0,
                    pressure_min: float = 1.0,
                    pressure_max: float = 400.0,
                    n_points: int = 40,
                    bubble_point_bar: float = None,
                    rock_compressibility: float = 4.5e-5,
                    temperature_unit: str = "C") -> PVTTable:
    """Korrelyasiyalardan tam PVT cədvəli qurur (iki fazalı: neft-su).

    `gas_gravity` qaz fazası üçün DEYİL — Standing (Rs) və
    Vazquez-Beggs (Bo) korrelyasiyaları həll olmuş qazın nisbi
    sıxlığını arqument kimi tələb edir, ona görə iki fazalı modeldə də
    lazımdır.

    `temperature_unit` — `temperature_c` parametrinin FAKTİKİ vahidi
    ("C"/"F"/"K"; defolt "C" — ad `temperature_c` olsa da, korrelyasiyalar
    daxildə onsuz da Fahrenheit-ə çevirir, `C_TO_F` ilə). Defolt DƏYİŞMİR
    — yalnız "C"-dən fərqli vahid verildikdə çevirmə tətbiq olunur.
    Korrelyasiyaların ÖZÜ (aşağı) TOXUNULMAYIB.
    """
    if temperature_unit != "C":
        temperature_c = convert_temperature(temperature_c, temperature_unit, "C")
    pressure = np.linspace(pressure_min, pressure_max, int(n_points))

    rs_saturated = standing_solution_gor(pressure, api, gas_gravity, temperature_c)
    if bubble_point_bar is None:
        bubble_point_bar = float(pressure_max * 0.6)
    rs_at_pb = float(standing_solution_gor(np.array([bubble_point_bar]), api,
                                           gas_gravity, temperature_c)[0])

    # Rs YALNIZ hesabat/diaqnostika üçündür (bax `ReservoirModel.diagnose`,
    # `rendering/renderers.py`) — mühərrikin qalıq/Jakobian hesablamasında
    # HEÇ YERDƏ istifadə OLUNMUR (bu, iki fazalı modeldə sərbəst qaz
    # fazasının izlənmədiyinin birbaşa nəticəsidir). Ona görə Pb-də
    # doymuş qaz-neft nisbətinin "qırılması" saxlanılır — bu, YALNIZ
    # görüntüləmə üçündür və Nyutona TƏSİR ETMİR.
    rs = np.where(pressure < bubble_point_bar, rs_saturated, rs_at_pb)

    # Bo(p) və μo(p) İSƏ Nyutonun HƏLL ETDİYİ tənliklərə birbaşa girir
    # (bax `BlackOilPVTProvider`/`ResidualAssembler`). TAPILAN SƏHV: bu
    # ikisi əvvəllər Pb-də FƏRQLİ düsturlara keçirdi (doymuş qaz-neft
    # qarışığı ↔ doymamış maye) — DƏYƏR kəsilməzdir, lakin TÖRƏMƏ Pb-də
    # sıçrayır (∂Bo/∂p işarə dəyişir). Bu, Nyutonu Pb ətrafında sonsuz
    # OSSİLYASİYAYA sala bilir (ölçüldü və sənədləşdirilib —
    # `test_line_search_prevents_infinite_oscillation_near_a_well`).
    #
    # Qərar: qaz fazası onsuz da modelləşdirilmədiyi üçün (istifadəçiyə
    # UI-də artıq bildirilir: "nəticələr nikbin ola bilər"), Pb-dən
    # AŞAĞIDA da EYNİ (doymamış maye) düsturunu davam etdiririk —
    # doymuş qarışığa KEÇMİRİK. Bu, Bo/μo-nu BÜTÜN təzyiq oblastında
    # HAMAR (C¹) edir və qırılmanı kökündən aradan qaldırır — "yumşaq
    # uğursuzluq" kimi bir SONRAKI TƏDBİR deyil, məhz SƏBƏBİN özünün
    # düzəldilməsidir.
    bo_at_pb = float(vazquez_beggs_oil_fvf(np.array([rs_at_pb]), api,
                                           gas_gravity, temperature_c)[0])
    co = vazquez_beggs_undersaturated_compressibility(
        rs_at_pb, api, gas_gravity, temperature_c, bubble_point_bar)
    bo = bo_at_pb * np.exp(-co * (pressure - bubble_point_bar))

    mu_dead = beggs_robinson_dead_oil_viscosity(api, temperature_c)
    mu_at_pb = float(beggs_robinson_saturated_viscosity(
        mu_dead, np.array([rs_at_pb]))[0])
    ratio = np.maximum(pressure / max(bubble_point_bar, 1e-9), 1e-6)
    mu_oil = mu_at_pb * ratio ** 0.278

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
    )
