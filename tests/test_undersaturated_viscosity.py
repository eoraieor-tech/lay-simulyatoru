"""G2 (1-ci hissə) — doymamış neft özlülüyü μo(p, Rs), PROVIDER səviyyəsi.

    μo(p, Rs) = μo_sat(Pb(Rs)) · (p / Pb(Rs))^n

Bu, `oil_fvf_undersaturated` (B3-B) qolunun güzgüsüdür. Doymamış hüceyrədə
neftin tərkibi sabitdir, ona görə özlülük həmin Rs-in doyma təzyiqindən
başlayan qola aiddir — cədvəlin doymuş qoluna YOX.

ÜSTƏL HARADAN GƏLİR (ölçüldü, Seans 34): korrelyasiya sabit 0.278 işlədir,
SPE1 deck-inin qolları isə 0.4602 və 0.5802 verir (Seans 37-də düzəldilib). Ona görə üstəl CƏDVƏLİN ÖZ
doymamış sətirlərindən fit olunur; fit mümkün olmayanda korrelyasiya qiymətinə
düşülür VƏ xəbərdarlıq verilir.

BU ADDIMDA MÜHƏRRİK TOXUNULMUR — qoşulma ayrı commit-dir (G2b). Sonuncu test
məhz bunu kilidləyir.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from imex2d.domain.pvt import PVTTable
from imex2d.simulation.pvt.black_oil import (CORRELATION_VISCOSITY_EXPONENT,
                                             BlackOilPVTProvider)
from imex2d.simulation.pvt.correlations import build_pvt_table

BUBBLE_POINT = 200.0


def _synthetic_table(exponent: float, above_points: int = 6) -> PVTTable:
    """Doymamış qolu MƏLUM üstəllə qurulmuş cədvəl.

    Doymuş hissə (p < Pb): Rs və Bo artır, μo azalır — real davranış.
    Doymamış hissə (p ≥ Pb): μo = μo_b·(p/Pb)^n, Bo = Bo_b·exp(c_o·(Pb−p)).
    """
    below = np.linspace(20.0, BUBBLE_POINT, 9)
    above = np.linspace(BUBBLE_POINT, 400.0, above_points + 1)[1:]
    pressure = np.concatenate([below, above])

    rs_b, bo_b, mu_b, co = 120.0, 1.45, 0.62, 1.8e-4
    rs = np.concatenate([rs_b * below / BUBBLE_POINT,
                         np.full(above.size, rs_b)])
    bo = np.concatenate([1.05 + (bo_b - 1.05) * below / BUBBLE_POINT,
                         bo_b * np.exp(co * (BUBBLE_POINT - above))])
    mu = np.concatenate([mu_b * (BUBBLE_POINT / below) ** 0.25,
                         mu_b * (above / BUBBLE_POINT) ** exponent])
    return PVTTable(
        pressure=pressure, oil_fvf=bo, oil_viscosity=mu, solution_gor=rs,
        water_fvf=np.full(pressure.size, 1.02),
        water_viscosity=np.full(pressure.size, 0.5),
        bubble_point=BUBBLE_POINT,
        gas_fvf=0.9 * BUBBLE_POINT / pressure,
        gas_viscosity=np.full(pressure.size, 0.02))


def _correlation_provider() -> BlackOilPVTProvider:
    return BlackOilPVTProvider(build_pvt_table(
        pressure_min=1.0, pressure_max=400.0, n_points=40,
        bubble_point_bar=240.0, include_gas=True))


# ═══════════════════════ üstəlin fit olunması ════════════════════════

def test_exponent_is_fitted_from_the_table_not_hardcoded():
    """SPE1-in qolları 0.46-0.58 verir — sabit 0.278 YARAMIR."""
    provider = BlackOilPVTProvider(_synthetic_table(0.46))
    assert provider.viscosity_exponent_fitted
    assert provider.undersaturated_viscosity_exponent == pytest.approx(0.46, rel=1e-6)


def test_a_different_exponent_is_recovered_too():
    provider = BlackOilPVTProvider(_synthetic_table(0.5802))
    assert provider.undersaturated_viscosity_exponent == pytest.approx(0.5802,
                                                                       rel=1e-6)


def test_correlation_table_reproduces_its_own_exponent():
    """Korrelyasiya cədvəli 0.278 ilə qurulur — fit onu geri tapmalıdır."""
    provider = _correlation_provider()
    assert provider.viscosity_exponent_fitted
    assert provider.undersaturated_viscosity_exponent == pytest.approx(
        CORRELATION_VISCOSITY_EXPONENT, rel=1e-6)


def test_fallback_warns_when_the_table_has_no_undersaturated_rows():
    """Səssiz ehtiyat qiymət YOXDUR — xəbərdarlıq yazılır (Q-27 qaydası)."""
    table = _synthetic_table(0.46)
    thin = PVTTable(
        pressure=table.pressure, oil_fvf=table.oil_fvf,
        oil_viscosity=table.oil_viscosity, solution_gor=table.solution_gor,
        water_fvf=table.water_fvf, water_viscosity=table.water_viscosity,
        bubble_point=float(table.pressure[-1]),      # Pb = maksimum təzyiq
        gas_fvf=table.gas_fvf, gas_viscosity=table.gas_viscosity)

    logger = logging.getLogger("imex2d.simulation.pvt.black_oil")
    messages = []

    class Collector(logging.Handler):
        def emit(self, record):
            messages.append(record.getMessage())

    handler = Collector(level=logging.WARNING)
    logger.addHandler(handler)
    try:
        provider = BlackOilPVTProvider(thin)
    finally:
        logger.removeHandler(handler)

    assert not provider.viscosity_exponent_fitted
    assert provider.undersaturated_viscosity_exponent == pytest.approx(
        CORRELATION_VISCOSITY_EXPONENT)
    assert any("özlülük üstəli" in message for message in messages), messages


# ═══════════════════════ fiziki davranış ═════════════════════════════

def test_branch_is_continuous_at_the_bubble_point():
    """Rs = Rs_sat(p) olanda Pb(Rs) = p, yəni nəticə doymuş qolun özüdür."""
    provider = _correlation_provider()

    # ÖLÇÜLDÜ (Seans 36): cədvəl DÜYÜNLƏRİNDƏ uyğunluq TAM DƏQİQDİR.
    # Düyünlər arasında, xüsusən son düyünlə Pb arasında, xam sütun sınığı
    # xətti kəsir və fərq 7.6e-3 nisbi səviyyəyə çatır — bu, qolun deyil,
    # cədvəlin öz interpolyasiyasının xətasıdır.
    nodes = provider.table.pressure
    nodes = nodes[(nodes >= 100.0) & (nodes <= 239.0)]
    assert np.allclose(provider.oil_viscosity_undersaturated(
        nodes, provider.solution_gor(nodes)),
        provider.oil_viscosity(nodes), rtol=1e-12)

    pressure = np.array([120.0, 200.0, 239.0])
    rs = provider.solution_gor(pressure)
    assert np.allclose(provider.oil_viscosity_undersaturated(pressure, rs),
                       provider.oil_viscosity(pressure), rtol=1e-2)


def test_undersaturated_oil_thickens_with_pressure():
    """SPE1-də 0.51 → 0.74 cP: doymamış neft təzyiqlə QATILAŞIR."""
    provider = BlackOilPVTProvider(_synthetic_table(0.46))
    rs = np.full(4, 120.0)
    pressure = np.array([200.0, 260.0, 330.0, 400.0])
    values = provider.oil_viscosity_undersaturated(pressure, rs)
    assert np.all(np.diff(values) > 0.0)


def test_saturated_branch_would_miss_the_thickening():
    """Doymuş qoldan oxumaq artımı İTİRİR — boşluğun fiziki mənası budur."""
    provider = BlackOilPVTProvider(_synthetic_table(0.46))
    pressure = np.array([400.0])
    rs = np.array([60.0])                     # Pb ≈ 100 bar, yəni çox doymamış
    branch = provider.oil_viscosity_undersaturated(pressure, rs)
    assert branch[0] > provider.oil_viscosity(pressure)[0] * 1.2


# ═══════════════════════ törəmələr — sonlu fərq ══════════════════════

def test_pressure_derivative_matches_finite_difference():
    provider = _correlation_provider()
    pressure = np.array([250.0, 300.0, 380.0])
    rs = np.array([150.0, 150.0, 150.0])
    step = 1e-4
    analytic, _ = provider.oil_viscosity_undersaturated_derivatives(pressure, rs)
    numeric = (provider.oil_viscosity_undersaturated(pressure + step, rs)
               - provider.oil_viscosity_undersaturated(pressure - step, rs)) / (2 * step)
    assert np.allclose(analytic, numeric, rtol=1e-6)


def test_solution_gor_derivative_matches_finite_difference():
    """∂μo/∂Rs ≠ 0 — doymamış hüceyrədə 3-cü dəyişən məhz Rs-dir."""
    provider = _correlation_provider()
    pressure = np.array([300.0, 350.0])
    rs = np.array([100.0, 140.0])      # cədvəlin Rs diapazonu: 0.70…153.87
    step = 1e-4
    _, analytic = provider.oil_viscosity_undersaturated_derivatives(pressure, rs)
    numeric = (provider.oil_viscosity_undersaturated(pressure, rs + step)
               - provider.oil_viscosity_undersaturated(pressure, rs - step)) / (2 * step)
    assert np.allclose(analytic, numeric, rtol=1e-5)
    assert np.all(np.abs(analytic) > 0.0)


def test_derivative_is_zero_outside_the_tables_solution_gor_range():
    """QƏSDİ DAVRANIŞ (Bo qolu ilə EYNİ qayda): cədvəlin Rs diapazonundan
    kənarda `dPb/dRs = 0`, ona görə ∂μo/∂Rs də sıfırdır — ekstrapolyasiya
    edilmir.

    ÖLÇÜLDÜ: korrelyasiya cədvəlinin doymuş Rs diapazonu 0.702…153.869
    sm³/sm³-dir; Rs = 170 ondan kənardadır. Sonlu fərq də sıfır verir,
    yəni analitik törəmə cədvəlin öz davranışı ilə UYĞUNDUR.
    """
    provider = _correlation_provider()
    pressure, rs, step = np.array([300.0]), np.array([170.0]), 1e-4
    _, analytic = provider.oil_viscosity_undersaturated_derivatives(pressure, rs)
    numeric = (provider.oil_viscosity_undersaturated(pressure, rs + step)
               - provider.oil_viscosity_undersaturated(pressure, rs - step)
               ) / (2 * step)
    assert analytic[0] == pytest.approx(0.0, abs=1e-15)
    assert numeric[0] == pytest.approx(0.0, abs=1e-15)


def test_derivatives_are_consistent_with_the_exponent():
    """∂μ/∂p = n·μ/p — analitik eyniliyin özü."""
    provider = BlackOilPVTProvider(_synthetic_table(0.46))
    pressure = np.array([260.0, 340.0])
    rs = np.full(2, 120.0)
    mu = provider.oil_viscosity_undersaturated(pressure, rs)
    d_dp, _ = provider.oil_viscosity_undersaturated_derivatives(pressure, rs)
    assert np.allclose(d_dp, 0.46 * mu / pressure, rtol=1e-10)


# ═══════════════════════ mühərrik HƏLƏ dəyişmir ══════════════════════

def test_engine_wiring_happened_in_the_next_step():
    """G2a-da bu test QOŞULMAMAĞI yoxlayırdı — G2b-də qoşulma baş verdi.

    Test silinmir, TƏRSİNƏ çevrilir: staging-in tamamlandığını sənədləşdirir
    (provider ayrı commit-də yazıldı, mühərriyə qoşulma ayrı commit-də).
    """
    import inspect
    from imex2d.simulation.implicit import three_phase_newton
    source = inspect.getsource(three_phase_newton.ThreePhaseNewtonSolver)
    assert "_oil_viscosity" in source
    assert "oil_viscosity_undersaturated" in source
