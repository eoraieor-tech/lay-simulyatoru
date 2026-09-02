"""FINAL CONSOLIDATION audit §19 — sadə 1D homogen rezervuar üçün əlavə
ədədi keyfiyyət benchmark-ları, analitik gözlənti ilə müqayisə.

Bunlar `test_physics.py`-dəki Bukley-Leverett (cəbhə mövqeyi/shock
doyumluluğu — DOYUMLULUQ NƏQLİYYATI) və `test_discretization.py`-dəki
tək-üz transmissivlik əl-hesabı (BİR ÜZ düsturu) testlərini TƏKRARLAMIR —
onların ARASINDAKI boşluğu doldurur:

    1) YIĞILMIŞ (bütün şəbəkə, N hüceyrə) TPFA sisteminin özünün
       kontinuum Darcy həllini DƏQİQ (maşın dəqiqliyi ilə) reproduksiya
       etdiyini sübut edir — bu, tək-üz düsturunun düzgünlüyündən
       FƏRQLİ bir sual (əlaqələndirmə/indeksləşdirmə düzgünlüyü).
    2) Zaman addımının INCƏLƏŞDİRİLMƏSİNİN nəticəni "partlatmadığını"
       (instability), əksinə kiçik, sərhədli (converging) şəkildə
       dəyişdiyini göstərir.
"""

from __future__ import annotations

import numpy as np

from helpers import default_scal, make_service, one_dimensional_model
from imex2d.application.config import OutputConfig, SimulationConfig, TimeSteppingConfig
from imex2d.simulation.discretization import TwoPointFluxDiscretization


def _steady_single_phase_pressure(model, rate: float, mu: float, bhp: float) -> np.ndarray:
    """TPFA transmissivlik ŞƏBƏKƏSİNDƏN (relperm/doyumluluq OLMADAN,
    təmiz tək-fazalı) sabit-hal təzyiq sahəsini BİRBAŞA həll edir.

    Bu, `TwoPointFluxDiscretization`-ın çıxardığı `transmissibility`
    massivinin ÖZÜNÜ (artıq `test_discretization.py`-də tək üz üçün əl
    hesabı ilə doğrulanıb) N-hüceyrəlik bir xətti sistemə YIĞARAQ
    işlədir — yəni ƏLAQƏLƏNDİRMƏ (`Connections.cell_a`/`cell_b`
    indeksləşdirməsi) düzgünlüyünü sınayır, tək üzün düsturunu YOX.
    """
    grid = TwoPointFluxDiscretization().build(model)
    conn = grid.connections
    n = model.grid.ncell
    system = np.zeros((n, n))
    for face in range(conn.count):
        i, j = int(conn.cell_a[face]), int(conn.cell_b[face])
        t = grid.transmissibility[face] / mu
        system[i, i] += t
        system[j, j] += t
        system[i, j] -= t
        system[j, i] -= t
    rhs = np.zeros(n)
    rhs[0] += rate                      # inyeksiya: hüceyrə 0-a daxil olan sabit debit
    system[n - 1, :] = 0.0
    system[n - 1, n - 1] = 1.0
    rhs[n - 1] = bhp                    # hasilat: son hüceyrədə sabit BHP (Dirichlet)
    return np.linalg.solve(system, rhs)


def test_tpfa_network_reproduces_analytical_linear_darcy_pressure_profile():
    """1D homogen xətti axın: yığılmış TPFA şəbəkəsinin sabit-hal həlli
    kontinuum Darcy həlli ilə maşın-dəqiqliyi səviyyəsində üst-üstə
    düşməlidir: P(x) = P_L + q·μ/(k·A)·(L−x).
    """
    scal = default_scal()
    nx, dx, dy, dz, porosity, permeability, rate, bhp = 60, 8.0, 100.0, 10.0, 0.20, 200.0, 60.0, 200.0
    model = one_dimensional_model(nx=nx, dx=dx, dy=dy, dz=dz, porosity=porosity,
                                  permeability=permeability, injection_rate=rate, scal=scal)
    mu = model.fluids.oil_viscosity

    pressure = _steady_single_phase_pressure(model, rate, mu, bhp)
    x_cells = (np.arange(nx) + 0.5) * dx
    slope_numeric, _intercept = np.polyfit(x_cells, pressure, 1)

    darcy_constant = model.units.darcy_constant
    face_area = dy * dz
    slope_analytical = -rate * mu / (darcy_constant * permeability * face_area)

    # profil dəqiq XƏTTİ olmalıdır (residual demək olar sıfır)
    residual = pressure - np.polyval([slope_numeric, _intercept], x_cells)
    assert np.max(np.abs(residual)) < 1e-6, "Təzyiq profili XƏTTİ deyil (gözlənilməz)"

    relative_error = abs(slope_numeric - slope_analytical) / abs(slope_analytical)
    assert relative_error < 1e-8, (
        f"Yığılmış TPFA sistemi analitik Darcy meyllindən {relative_error:.3e} "
        f"fərqlənir (analitik={slope_analytical:.6f}, ədədi={slope_numeric:.6f})")


def test_timestep_refinement_converges_not_diverges():
    """Adaptiv (CFL-əsaslı) zaman addımı artıq həll üçün YETƏRLİ addım
    ölçüsünü seçirsə, `max_dt`-i daha da azaltmaq nəticəni DƏYİŞMƏMƏLİDİR
    (CFL nəzarətçisi artıq məhdudlaşdırıcı amildir). `max_dt`-i CFL-in
    seçəcəyindən DAHA KİÇİK məcbur etdikdə isə nəticə YALNIZ KİÇİK,
    SƏRHƏDLİ şəkildə dəyişməlidir (yığılma), heç bir halda "partlamamalı"
    (kəskin sıçrayış/NaN) — bu, ƏDƏDİ SABİTLİYİN birbaşa sübutudur.
    """
    scal = default_scal()
    rate, end_time = 60.0, 250.0
    nx, dx, dy, dz, porosity = 120, 8.0, 100.0, 10.0, 0.20

    def _run(max_dt):
        model = one_dimensional_model(nx=nx, dx=dx, dy=dy, dz=dz, porosity=porosity,
                                      injection_rate=rate, scal=scal)
        config = SimulationConfig(
            end_time=end_time,
            time_stepping=TimeSteppingConfig(max_dt=max_dt, cfl_factor=0.4),
            output=OutputConfig(snapshot_count=2))
        result = make_service(scal).run(model, config)
        return result.series.cumulative_oil[-1], result.steps

    # max_dt CFL-in seçdiyindən BÖYÜK olanda (16→2) CFL nəzarətçisi artıq
    # məhdudlaşdırıcıdır — nəticə TAM sabit qalmalıdır.
    coarse_values = [_run(dt)[0] for dt in (16.0, 8.0, 4.0, 2.0)]
    assert np.allclose(coarse_values, coarse_values[0], rtol=0, atol=1e-9), (
        "CFL-məhdud rejimdə fərqli max_dt fərqli nəticə verdi — adaptiv "
        "nəzarətçi gözlənildiyi kimi işləmir")

    # max_dt-i CFL-in seçəcəyindən DAHA KİÇİK məcbur etsək (250 addım),
    # nəticə YALNIZ kiçik miqdarda dəyişməlidir — ƏDƏDİ PARTLAMA yoxdur.
    fine_value, fine_steps = _run(1.0)
    assert fine_steps > 198, "max_dt=1.0 gözlənildiyi kimi daha çox addım tələb etmədi"
    relative_change = abs(fine_value - coarse_values[0]) / abs(coarse_values[0])
    assert relative_change < 0.01, (
        f"Zaman addımını incələşdirmək kumulyativ hasilatı {relative_change * 100:.3f}% "
        "dəyişdi — gözlənilən kiçik (<1%) yığılma sərhədini aşır")
    assert np.isfinite(fine_value), "İncə zaman addımında nəticə sonsuz/NaN oldu"
