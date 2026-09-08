"""Bakli-Leverett ANALİTİK ETALONUNUN doğruluq testləri.

`tests/test_physics.py` mühərriki etalonla müqayisə edir; bu fayl isə
ETALONUN ÖZÜNÜ yoxlayır. Sınıq etalon fizika səhvindən təhlükəlidir:
o, mühərrikin gələcək səhvlərini maskalayır.

Testlərin əsas dayağı — Bakli-Leverett həllində DƏQİQ ödənən kütlə
eyniliyidir:

    ∫₀^∞ (Sw(x) − Swi) dx  =  q·t / (φ·A)  =  velocity · t

Bu, modulun öz daxili hesabına yox, sıxılmaz iki-fazalı sıxışdırmanın
həcm balansına söykənir — yəni müstəqil ölçüdür.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import brentq

from helpers import (SKIP_SLOW, bl_config, default_scal, make_service,
                     one_dimensional_model)
from imex2d.domain.scal import CoreyParameters
from imex2d.simulation.analytical import buckley_leverett

# ── Etalon məsələ: A3 tapşırığında ölçülmüş konfiqurasiya ──────────────
MU_W, MU_O = 0.5, 3.0
POROSITY, RATE, AREA, TIME = 0.2, 60.0, 1000.0, 250.0
VELOCITY = RATE / (AREA * POROSITY)          # m/gün, frontal sürət
EXPECTED_INTEGRAL = VELOCITY * TIME          # = 75.0 m


def _solution(scal=None, **kwargs):
    scal = scal or default_scal()
    return buckley_leverett(scal, MU_W, MU_O, POROSITY, RATE, AREA, TIME,
                            **kwargs)


def _swept_volume(distance, saturation, swi) -> float:
    """∫ (Sw − Swi) dx — profil altındakı süpürülmüş su həcmi (vahid en)."""
    return float(np.trapezoid(np.asarray(saturation) - swi,
                              np.asarray(distance)))


# ══════════════════════════════════════════════════ müstəqil Welge həlli
# Modulun nəticəsini modulun özü ilə müqayisə edən test heç nə yoxlamır.
# Ona görə fraksion axın və onun törəməsi burada ANALİTİK (np.gradient-siz)
# yenidən qurulur, kök isə brentq ilə tapılır.

def _fractional_flow(scal: CoreyParameters, s, mu_w, mu_o):
    lw = scal.krw(s) / mu_w
    lo = scal.kro(s) / mu_o
    return lw / (lw + lo)


def _fractional_flow_derivative(scal: CoreyParameters, s, mu_w, mu_o):
    """df/ds — zəncir qaydası ilə, ədədi diferensiallaşdırma olmadan."""
    den = 1.0 - scal.swc - scal.sor
    sn = (s - scal.swc) / den
    lw = scal.krw_end * sn ** scal.nw / mu_w
    lo = scal.kro_end * (1.0 - sn) ** scal.no / mu_o
    dlw = scal.krw_end * scal.nw * sn ** (scal.nw - 1.0) / (den * mu_w)
    dlo = -scal.kro_end * scal.no * (1.0 - sn) ** (scal.no - 1.0) / (den * mu_o)
    return (dlw * lo - lw * dlo) / (lw + lo) ** 2


def _welge_shock(scal: CoreyParameters, swi, mu_w, mu_o):
    """Welge toxunanı: df/ds(s) = (f(s) − f_i)/(s − swi). (s_shock, chord)."""
    f_i = _fractional_flow(scal, swi, mu_w, mu_o)

    def residual(s):
        chord = (_fractional_flow(scal, s, mu_w, mu_o) - f_i) / (s - swi)
        return _fractional_flow_derivative(scal, s, mu_w, mu_o) - chord

    # Kobud tarama ilə işarə dəyişməsini tap, sonra dəqiq kök.
    grid = np.linspace(swi + 1e-4, 1.0 - scal.sor - 1e-4, 400)
    values = np.array([residual(s) for s in grid])
    sign_change = np.nonzero(np.diff(np.sign(values)) != 0)[0]
    assert sign_change.size, "Welge toxunanı üçün işarə dəyişməsi tapılmadı"
    index = sign_change[-1]
    root = brentq(residual, grid[index], grid[index + 1], xtol=1e-14)
    chord = (_fractional_flow(scal, root, mu_w, mu_o) - f_i) / (root - swi)
    return float(root), float(chord)


# ══════════════════════════════════════════════════════════ §5.1 ƏSAS TEST
def test_swept_volume_matches_injected_volume():
    """KÜTLƏ EYNİLİYİ — bu modulun əsas reqressiya testi.

    Profil altındakı sahə vurulan həcmə bərabər olmalıdır. Şok mailli
    xətlə "hamarlanarsa" sahə şişir: sınıq kod 24.75 % xəta verirdi.
    """
    scal = default_scal()
    solution = _solution(scal)
    integral = _swept_volume(solution.distance, solution.water_saturation,
                             scal.swc)
    error = abs(integral - EXPECTED_INTEGRAL) / EXPECTED_INTEGRAL * 100.0
    assert error < 0.1, (
        f"Kütlə eyniliyi pozuldu: ∫(Sw−Swi)dx = {integral:.4f}, "
        f"gözlənilən q·t/(φ·A) = {EXPECTED_INTEGRAL:.4f} → xəta "
        f"{error:.3f} % (limit 0.1 %)")


# ═══════════════════════════════════════════════════════════ §5.2 quyruq
def test_saturation_behind_front_is_exactly_initial():
    """Cəbhədən sonra toxunulmamış zonadır: Sw DƏQİQ Swi (tolerantlıq yox)."""
    scal = default_scal()
    solution = _solution(scal)
    tail = solution.distance > solution.front_position
    assert tail.any(), "Cəbhədən sonra heç bir nöqtə yoxdur"
    excess = float(solution.water_saturation[tail].max()) - scal.swc
    assert np.all(solution.water_saturation[tail] == scal.swc), (
        "Cəbhədən sonra Sw ≠ Swi — şok mailli keçidlə əvəz olunub: "
        f"maksimum artıqlıq {excess:.6f}")


# ═════════════════════════════════════════════ §5.3 ciddi artan məsafə
def test_distance_is_strictly_increasing():
    """`np.interp` və matplotlib ciddi artan x tələb edir."""
    solution = _solution()
    diff = np.diff(solution.distance)
    assert np.all(diff > 0), (
        f"{int(np.count_nonzero(diff <= 0))} nöqtədə x artmır — şok üçün "
        "`np.nextafter` işlədilməlidir")


# ══════════════════════════════════════════════════════ §5.4 monotonluq
def test_saturation_is_non_increasing_with_distance():
    """Sıxışdırma profili məsafə boyu azalmayan ola bilməz."""
    solution = _solution()
    assert np.all(np.diff(solution.water_saturation) <= 0), (
        "Sw(x) hansısa nöqtədə artdı — profil sıralaması pozulub")


# ═══════════════════════════════════════════ §5.5 müstəqil təsdiq (Welge)
def test_shock_matches_independent_welge_root():
    """Şok doyumluluğu və cəbhə mövqeyi müstəqil brentq kökü ilə üst-üstə."""
    scal = default_scal()
    solution = _solution(scal)
    root, chord = _welge_shock(scal, scal.swc, MU_W, MU_O)
    expected_front = VELOCITY * chord * TIME

    assert solution.shock_saturation == pytest.approx(root, rel=1e-4), (
        f"shock Sw: modul {solution.shock_saturation:.6f}, "
        f"müstəqil kök {root:.6f}")
    assert solution.front_position == pytest.approx(expected_front, rel=1e-4), (
        f"x_front: modul {solution.front_position:.4f}, "
        f"müstəqil {expected_front:.4f}")


# ═════════════════════════════════════════════════════ §5.6 zaman miqyası
def test_front_position_scales_linearly_with_time():
    """BL xəttidir: t iki dəfə artanda x_front da dəqiq iki dəfə artır."""
    scal = default_scal()
    single = buckley_leverett(scal, MU_W, MU_O, POROSITY, RATE, AREA, TIME)
    double = buckley_leverett(scal, MU_W, MU_O, POROSITY, RATE, AREA, 2 * TIME)
    assert double.front_position == pytest.approx(2 * single.front_position,
                                                  rel=1e-12)
    assert double.shock_saturation == pytest.approx(single.shock_saturation,
                                                    rel=1e-12)


# ═══════════════════════════════════════════════════════ §5.7 breakthrough
def test_breakthrough_flag_is_reported():
    """Cəbhə modeldən çıxıbsa, bu AÇIQ bildirilməlidir — səssiz faiz yox."""
    scal = default_scal()
    inside = _solution(scal, length=400.0)
    assert inside.breakthrough is False
    assert inside.front_position < 400.0

    outside = _solution(scal, length=100.0)
    assert outside.breakthrough is True, (
        f"x_front = {outside.front_position:.1f} m > length = 100 m, "
        "amma breakthrough bayrağı qalxmadı")


def test_profile_ends_at_requested_length():
    """`length` verildikdə profil ora qədər uzanır — sehrli 1.6 sabiti yox."""
    solution = _solution(length=960.0)
    assert solution.distance[-1] == pytest.approx(960.0, rel=1e-12)


# ══════════════════════════════════════════════════ §5.8 giriş yoxlaması
@pytest.mark.parametrize("kwargs", [
    {"time": 0.0},
    {"time": -10.0},
    {"total_rate": 0.0},
    {"total_rate": -60.0},
    {"area": 0.0},
    {"porosity": 0.0},
    {"porosity": -0.2},
])
def test_invalid_inputs_raise_value_error(kwargs):
    scal = default_scal()
    arguments = dict(mu_w=MU_W, mu_o=MU_O, porosity=POROSITY,
                     total_rate=RATE, area=AREA, time=TIME)
    arguments.update(kwargs)
    with pytest.raises(ValueError):
        buckley_leverett(scal, **arguments)


@pytest.mark.parametrize("sw_initial", [0.1, 0.75, 0.8])
def test_invalid_initial_saturation_raises(sw_initial):
    """swc ≤ Swi < 1 − sor olmalıdır (defoltlarda 0.20 ≤ Swi < 0.75)."""
    scal = default_scal()
    with pytest.raises(ValueError):
        buckley_leverett(scal, MU_W, MU_O, POROSITY, RATE, AREA, TIME,
                         sw_initial=sw_initial)


# ═══════════════════════════════════ diaqnostika: səssiz düzəliş qadağandır
def test_monotonicity_correction_count_is_exposed():
    """`np.maximum.accumulate` heç nə düzəltməməlidir — sayı görünsün."""
    solution = _solution()
    assert solution.monotonicity_corrections == 0, (
        f"Rarefaksiyada {solution.monotonicity_corrections} nöqtə zorla "
        "monotonlaşdırıldı — SCAL parametrləri qeyri-adidir")


# ═════════════════════════════════════════ §5.9 GRID YAXINSAMASI (ƏSL SÜBUT)
LENGTH = 960.0


def _numerical_profile(nx, length=LENGTH):
    """Eyni fiziki məsələ, fərqli şəbəkə addımı ilə."""
    scal = default_scal()
    dx, dy, dz = length / nx, 100.0, 10.0
    model = one_dimensional_model(nx=nx, dx=dx, dy=dy, dz=dz,
                                  porosity=POROSITY, injection_rate=RATE,
                                  scal=scal)
    result = make_service(scal).run(model, bl_config(TIME))
    x_cells = (np.arange(nx) + 0.5) * dx
    return x_cells, result.snapshots[-1].water_saturation.ravel()


def _rms_against_analytical(x_cells, sw_numeric, length=LENGTH):
    scal = default_scal()
    solution = buckley_leverett(scal, MU_W, MU_O, POROSITY, RATE,
                                100.0 * 10.0, TIME, length=length)
    sw_exact = np.interp(x_cells, solution.distance, solution.water_saturation)
    return float(np.sqrt(np.mean((sw_numeric - sw_exact) ** 2)))


def test_refining_the_grid_reduces_error_against_analytical():
    """Şəbəkə xırdalandıqca ədədi həll analitikə YAXINLAŞMALIDIR.

    Sınıq etalonla xəta əksinə ARTIRDI (0.0498 → 0.0564 → 0.0606) — çünki
    ölçü cihazının özü səhv idi. Bu test düzəlişin ən güclü sübutudur.
    """
    if SKIP_SLOW:
        pytest.skip("IMEX_SKIP_SLOW")
    errors = []
    for nx in (60, 120, 240):
        x_cells, sw_numeric = _numerical_profile(nx)
        errors.append(_rms_against_analytical(x_cells, sw_numeric))
    trace = " → ".join(f"{e:.4f}" for e in errors)
    assert errors[1] < errors[0] and errors[2] < errors[1], (
        f"RMS xəta şəbəkə xırdalandıqca azalmadı: {trace}")


def test_numerical_profile_conserves_injected_volume():
    """Mühərrikdən asılı olmayan yoxlama: ədədi profilin altındakı sahə.

    Bu, BL nəzəriyyəsindən asılı deyil — xalis həcm balansıdır və
    mühərrikin özünü sınayır.
    """
    if SKIP_SLOW:
        pytest.skip("IMEX_SKIP_SLOW")
    scal = default_scal()
    x_cells, sw_numeric = _numerical_profile(120)
    integral = _swept_volume(x_cells, sw_numeric, scal.swc)
    error = abs(integral - EXPECTED_INTEGRAL) / EXPECTED_INTEGRAL * 100.0
    assert error < 3.0, (
        f"Ədədi profil {integral:.3f} m, vurulan {EXPECTED_INTEGRAL:.3f} m "
        f"→ {error:.2f} % (hüceyrə mərkəzi trapesiya xətası daxil)")


# ══════════════════════════════════ §4 ölçmə metrikası: orta nöqtə keçidi
def _midpoint_front(x_cells, saturation, level):
    """`main_window` metrikası — Qt pəncərəsi qurmadan (statik metoddur)."""
    QtWidgets = pytest.importorskip("PyQt5.QtWidgets")  # noqa: F841
    from imex2d.ui.main_window import MainWindow
    return MainWindow._front_by_midpoint(x_cells, saturation, level)


def test_midpoint_metric_recovers_a_symmetrically_smeared_front():
    """Simmetrik yayılmış cəbhədə orta nöqtə DƏQİQ mərkəzi qaytarır."""
    scal = default_scal()
    shock, swi = 0.54, scal.swc
    level = 0.5 * (shock + swi)
    centre = 200.0
    x_cells = np.linspace(0.0, 400.0, 401)
    # tanh — mərkəzə görə tam simmetrik yayılma
    saturation = swi + (shock - swi) * 0.5 * (1.0 - np.tanh((x_cells - centre) / 12.0))
    assert _midpoint_front(x_cells, saturation, level) == pytest.approx(centre,
                                                                       abs=0.5)


def test_midpoint_metric_is_less_biased_than_leading_edge():
    """Ön kənar metrikası ədədi diffuziyanın quyruğuna görə ŞİŞİRDİR."""
    if SKIP_SLOW:
        pytest.skip("IMEX_SKIP_SLOW")
    scal = default_scal()
    x_cells, sw_numeric = _numerical_profile(120)
    solution = buckley_leverett(scal, MU_W, MU_O, POROSITY, RATE,
                                100.0 * 10.0, TIME, length=LENGTH)

    leading_edge = float(x_cells[sw_numeric > scal.swc + 0.01][-1])
    level = 0.5 * (solution.shock_saturation + scal.swc)
    midpoint = _midpoint_front(x_cells, sw_numeric, level)

    front = solution.front_position
    assert leading_edge > front, "Ön kənar cəbhəni qabaqlamalıdır (qərəz)"
    assert abs(midpoint - front) < abs(leading_edge - front), (
        f"orta nöqtə {midpoint:.1f} m, ön kənar {leading_edge:.1f} m, "
        f"analitik {front:.1f} m — orta nöqtə daha yaxın olmalıdır")


# ═══════════════════════════ §5.10 kapilyar/cazibə izolyasiyası
def test_validation_model_has_no_gravity_and_no_capillary():
    """Validasiya modeli 1 təbəqəli və üfüqidir; BL bunları nəzərə almır."""
    scal = default_scal()
    service = make_service(scal)
    model = one_dimensional_model(scal=scal)
    engine = service.create_engine(model, bl_config())
    assert engine._has_gravity is False, (
        "Validasiya modelində cazibə aktivləşdi — BL üfüqi nəzəriyyədir")
    assert engine.capillary is None, (
        "Validasiya modelinə kapilyar provider sızdı — BL Pc = 0 qəbul edir")
