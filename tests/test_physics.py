"""Fiziki doğruluq testləri: material balansı və Bukley-Leverett.

Bunlar simulyatorun "düzgün işlədiyini" sübut edən yeganə testlərdir —
qalanları yalnız kodun quruluşunu yoxlayır.
"""

import numpy as np

from helpers import (SKIP_SLOW, bl_config, default_scal, five_spot_model,
                     make_service, one_dimensional_model, short_config)
from imex2d.simulation.analytical import buckley_leverett


def _trapezoid(y, x):
    return float(np.trapezoid(np.asarray(y), np.asarray(x)))


def test_water_material_balance_is_conserved():
    """Layda toplanan su = vurulan − çıxarılan. Xəta < 0.5 % olmalıdır."""
    scal = default_scal()
    model = five_spot_model(scal=scal)
    service = make_service(scal)
    engine = service.create_engine(model, short_config(end_time=400.0))
    result = engine.run()

    pore_volume = model.pore_volume()
    initial_sw = model.initial_conditions.water_saturation
    accumulated = float(np.sum(pore_volume * (engine.sw - initial_sw)))

    series = result.series
    injected = _trapezoid(series.water_injection_rate, series.time)
    produced = _trapezoid(series.water_rate, series.time) * model.fluids.water_fvf
    net = injected - produced

    error = abs(accumulated - net) / max(injected, 1e-9) * 100.0
    assert error < 0.5, f"Material balans xətası {error:.3f} % (limit 0.5 %)"


def test_saturation_stays_within_physical_limits():
    scal = default_scal()
    result = make_service(scal).run(five_spot_model(scal=scal),
                                    short_config(end_time=400.0))
    for snapshot in result.snapshots:
        sw = snapshot.water_saturation
        assert sw.min() >= scal.swc - 1e-9, "Sw bağlı sudan aşağı düşdü"
        assert sw.max() <= 1.0 - scal.sor + 1e-9, "Sw qalıq neft həddini keçdi"


def test_cumulative_production_is_monotonic():
    scal = default_scal()
    result = make_service(scal).run(five_spot_model(scal=scal),
                                    short_config(end_time=400.0))
    cumulative = np.array(result.series.cumulative_oil)
    assert np.all(np.diff(cumulative) >= -1e-9), "Kumulyativ hasilat azaldı"
    assert np.all(np.array(result.series.recovery_factor) <= 100.0)


def test_buckley_leverett_front_position():
    """Ədədi cəbhə analitik cəbhə ilə 12 % daxilində üst-üstə düşməlidir.

    Tam üst-üstə düşmə gözlənilmir: upstream çəkilənmə ədədi diffuziya
    yaradır və cəbhəni yayır.
    """
    scal = default_scal()
    rate, end_time = 60.0, 250.0
    nx, dx, dy, dz, porosity = 120, 8.0, 100.0, 10.0, 0.20

    model = one_dimensional_model(nx=nx, dx=dx, dy=dy, dz=dz,
                                  porosity=porosity, injection_rate=rate,
                                  scal=scal)
    result = make_service(scal).run(model, bl_config(end_time))

    analytical = buckley_leverett(scal, model.fluids.water_viscosity,
                                  model.fluids.oil_viscosity, porosity,
                                  rate, dy * dz, end_time)
    x_cells = (np.arange(nx) + 0.5) * dx
    sw = result.snapshots[-1].water_saturation.ravel()
    swept = sw > scal.swc + 0.01
    assert swept.any(), "Su heç irəliləmədi"
    numerical_front = x_cells[swept][-1]

    error = abs(numerical_front - analytical.front_position) / analytical.front_position
    assert error < 0.12, (f"Cəbhə mövqeyi fərqi {error * 100:.1f} % "
                          f"(analitik {analytical.front_position:.1f} m, "
                          f"ədədi {numerical_front:.1f} m)")


def test_buckley_leverett_shock_saturation_is_reached():
    """Ədədi həllin maksimal doyumluluğu shock qiymətindən aşağı olmamalıdır."""
    scal = default_scal()
    model = one_dimensional_model(scal=scal)
    result = make_service(scal).run(model, bl_config())
    analytical = buckley_leverett(scal, model.fluids.water_viscosity,
                                  model.fluids.oil_viscosity, 0.20, 60.0,
                                  100.0 * 10.0, 250.0)
    sw_max = result.snapshots[-1].water_saturation.max()
    assert sw_max >= analytical.shock_saturation - 0.02
    assert scal.swc < analytical.shock_saturation < 1.0 - scal.sor


def test_higher_oil_viscosity_lowers_recovery():
    """Fiziki gözlənti: lözlük nisbəti pisləşdikcə RF düşür."""
    if SKIP_SLOW:
        return
    scal = default_scal()
    results = []
    for viscosity in (3.0, 40.0):
        model = five_spot_model(nx=21, ny=21, scal=scal)
        model.fluids.oil_viscosity = viscosity
        results.append(make_service(scal).run(
            model, short_config(end_time=600.0)).final_recovery_factor)
    assert results[1] < results[0], "Lözlü neftdə RF azalmadı"
