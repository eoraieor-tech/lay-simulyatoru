"""Phase 1 — vahid çevirmə qatı: dəqiq dəyərlər, dəyirmi-səyahət, çarpaz yoxlama."""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.domain import unit_conversions as uc


# ── konkret, əl ilə hesablanmış dəyərlər (tapşırıqda tələb olunan cütlər) ─
def test_psi_to_pa_known_value():
    assert uc.psi_to_pa(1.0) == pytest.approx(6894.757293168361, rel=1e-12)


def test_pa_to_psi_known_value():
    assert uc.pa_to_psi(6894.757293168361) == pytest.approx(1.0, rel=1e-12)


def test_bar_to_pa_known_value():
    assert uc.bar_to_pa(1.0) == pytest.approx(100000.0, rel=1e-12)


def test_pa_to_bar_known_value():
    assert uc.pa_to_bar(100000.0) == pytest.approx(1.0, rel=1e-12)


def test_bar_to_psi_matches_standard_constant():
    # standart neft-mühəndisliyi sabiti: 1 bar = 14.5037744 psi
    assert uc.bar_to_psi(1.0) == pytest.approx(14.5037744, rel=1e-6)


def test_ft_to_m_known_value():
    assert uc.ft_to_m(1.0) == pytest.approx(0.3048, rel=1e-12)


def test_m_to_ft_known_value():
    assert uc.m_to_ft(1.0) == pytest.approx(1.0 / 0.3048, rel=1e-12)


def test_md_to_m2_known_value():
    assert uc.md_to_m2(1.0) == pytest.approx(9.869232667160128e-16, rel=1e-12)


def test_darcy_to_m2_known_value():
    assert uc.darcy_to_m2(1.0) == pytest.approx(9.869232667160128e-13, rel=1e-12)


def test_cp_to_pas_known_value():
    assert uc.cp_to_pas(1.0) == pytest.approx(1.0e-3, rel=1e-12)


def test_pas_to_cp_known_value():
    assert uc.pas_to_cp(1.0e-3) == pytest.approx(1.0, rel=1e-12)


def test_bbl_to_m3_known_value():
    assert uc.bbl_to_m3(1.0) == pytest.approx(0.158987294928, rel=1e-12)


def test_m3_to_bbl_known_value():
    assert uc.m3_to_bbl(1.0) == pytest.approx(1.0 / 0.158987294928, rel=1e-12)


def test_stb_per_day_to_m3_per_day_matches_bbl_volume():
    # stb fiziki olaraq bbl ilə eyni həcmdir (42 qallon) — bax modul docstring-i
    assert uc.stb_per_day_to_m3_per_day(1.0) == pytest.approx(0.158987294928, rel=1e-12)


# ── dəyirmi-səyahət: dəyər -> vahid A -> SI -> vahid A -> orijinal dəyər ──
@pytest.mark.parametrize("quantity,unit,value", [
    ("pressure", "psi", 2500.0),
    ("pressure", "bar", 200.0),
    ("pressure", "kPa", 15000.0),
    ("pressure", "MPa", 12.5),
    ("length", "ft", 3280.8),
    ("length", "m", 1000.0),
    ("area", "ft2", 435600.0),
    ("area", "acre", 10.0),
    ("permeability", "mD", 150.0),
    ("permeability", "D", 2.5),
    ("viscosity", "cP", 3.2),
    ("volume", "bbl", 1000.0),
    ("volume", "stb", 5000.0),
    ("rate", "bbl/day", 500.0),
    ("rate", "stb/day", 1200.0),
    ("density", "lb/ft3", 55.0),
])
def test_round_trip_through_si(quantity, unit, value):
    si = uc.to_si(value, unit, quantity)
    back = uc.from_si(si, unit, quantity)
    assert back == pytest.approx(value, rel=1e-9)


@pytest.mark.parametrize("unit,value", [("C", 85.0), ("F", 185.0), ("K", 350.0)])
def test_round_trip_temperature(unit, value):
    kelvin = uc.convert_temperature(value, unit, "K")
    back = uc.convert_temperature(kelvin, "K", unit)
    assert back == pytest.approx(value, rel=1e-9)


@pytest.mark.parametrize("unit,value", [("bar", 4.5e-5), ("psi", 3.1e-6), ("Pa", 4.5e-10)])
def test_round_trip_compressibility(unit, value):
    si = uc.convert_compressibility(value, unit, "Pa")
    back = uc.convert_compressibility(si, "Pa", unit)
    assert back == pytest.approx(value, rel=1e-9)


def test_round_trip_through_engine_units():
    """`to_engine_units`/`from_engine_units` cütü də dəyirmi-səyahətdə itki verməməlidir."""
    original_psi = 3200.0
    engine_value = uc.to_engine_units(original_psi, "psi", "pressure")   # -> bar
    assert engine_value == pytest.approx(uc.psi_to_bar(original_psi), rel=1e-12)
    back_psi = uc.from_engine_units(engine_value, "psi", "pressure")
    assert back_psi == pytest.approx(original_psi, rel=1e-9)


# ── temperatur konkret dəyərlər (ofset yoxlaması) ───────────────────────
def test_celsius_to_fahrenheit_known_points():
    assert uc.convert_temperature(0.0, "C", "F") == pytest.approx(32.0)
    assert uc.convert_temperature(100.0, "C", "F") == pytest.approx(212.0)


def test_celsius_to_kelvin_known_point():
    assert uc.convert_temperature(0.0, "C", "K") == pytest.approx(273.15)


# ── sıxılma: TƏRS miqyaslama, mövcud kod ilə çarpaz yoxlama ─────────────
def test_compressibility_matches_existing_correlations_conversion():
    """`simulation/pvt/correlations.py`-dəki
    `vazquez_beggs_undersaturated_compressibility` psi⁻¹-i `* BAR_TO_PSI`
    ilə bar⁻¹-ə çevirir. Bu, sıxılmanın TƏRS miqyaslı olması səbəbindən
    DOĞRUDUR (audit bunu şübhəli sayıb, əl hesabı ilə YOXLANILIB): `co
    [1/bar] = co [1/psi] * 14.5037744`. Bu test elə bu faktı qoruyur."""
    co_psi = 3.0e-6
    BAR_TO_PSI = 14.5037744
    co_bar_existing_formula = co_psi * BAR_TO_PSI
    co_bar_via_layer = uc.convert_compressibility(co_psi, "psi", "bar")
    assert co_bar_via_layer == pytest.approx(co_bar_existing_formula, rel=1e-6)
    # sağlamlıq yoxlaması: nəticə tipik süxur/flüid sıxılması diapazonundadır
    assert 1e-6 < co_bar_via_layer < 1e-3


def test_compressibility_larger_pressure_unit_gives_larger_numeric_value():
    """1 bar = 14.5 psi (bar daha 'böyük' vahiddir) -> eyni fiziki
    sıxılma bar⁻¹-də PSİ⁻¹-DƏN böyük ədədlə ifadə olunur."""
    co_psi = 3.0e-6
    co_bar = uc.convert_compressibility(co_psi, "psi", "bar")
    assert co_bar > co_psi


# ── xəta halları ─────────────────────────────────────────────────────────
def test_convert_rejects_unknown_quantity():
    with pytest.raises(ValueError, match="kəmiyyət"):
        uc.convert(1.0, "m", "ft", "speed")


def test_convert_rejects_unknown_unit():
    with pytest.raises(ValueError, match="vahidi"):
        uc.convert(1.0, "furlong", "m", "length")


def test_convert_temperature_rejects_unknown_unit():
    with pytest.raises(ValueError, match="temperatur"):
        uc.convert_temperature(1.0, "R", "K")


def test_convert_compressibility_rejects_non_pressure_unit():
    with pytest.raises(ValueError, match="təzyiq"):
        uc.convert_compressibility(1.0, "m", "bar")


def test_quantity_wrapper_converts_to_engine_units():
    assert uc.Quantity(3000.0, "psi", "pressure").to_engine() == pytest.approx(
        uc.psi_to_bar(3000.0), rel=1e-12)


def test_known_units_lists_expected_members():
    assert "psi" in uc.known_units("pressure")
    assert "mD" in uc.known_units("permeability")
    assert set(uc.known_units("temperature")) == {"K", "C", "F"}


# ── vahiddən asılı olmayan fiziki nəticə: eyni ssenari iki vahid sistemində ─
def test_reservoir_model_builder_converts_rock_compressibility_unit():
    """Phase 1 (giriş boru xətti): `ReservoirModelBuilder.build()` indi
    `rock_compressibility_unit` qəbul edir — defolt ('bar') dəyişməzliyi
    qoruyur, fərqli vahid (psi) verilsə mühərrik vahidinə (1/bar) çevrilir."""
    from imex2d.application.model_builder import ReservoirModelBuilder
    from imex2d.application.scenarios import SyntheticGeologicalModelBuilder, five_spot

    geology = SyntheticGeologicalModelBuilder().build(
        nx=3, ny=3, dx=20.0, dy=20.0, dz=10.0, porosity=0.2, permx_base=100.0)
    wells = five_spot(geology.grid)

    default_model = ReservoirModelBuilder().build(geology, wells)
    assert default_model.rock.compressibility == pytest.approx(4.5e-5)

    co_psi = 3.1e-6
    model_from_psi = ReservoirModelBuilder().build(
        geology, wells, rock_compressibility=co_psi, rock_compressibility_unit="psi")
    assert model_from_psi.rock.compressibility == pytest.approx(
        uc.convert_compressibility(co_psi, "psi", "bar"))
    assert model_from_psi.rock.compressibility == pytest.approx(4.5e-5, rel=1e-3)


def test_darcy_flow_estimate_is_unit_invariant():
    """q = k·A·Δp / (μ·L) — eyni fiziki ssenarini METRIC (mD/m²/bar/m/cP)
    və FIELD (mD/ft²/psi/ft/cP) girişindən eyni SI axın sürətinə gətirməli.

    Bu, `discretization.py`-nin özünü DƏYİŞDİRMİR — yalnız YENİ çevirmə
    qatının Darsi düsturunun hər tərəfini DOĞRU çevirdiyini yoxlayır."""
    k_md, area_m2, dp_bar, length_m, mu_cp = 150.0, 200.0, 20.0, 50.0, 2.0

    def darcy_flux_si(k_m2, area_m2_, dp_pa, length_m_, mu_pas):
        return k_m2 * area_m2_ * dp_pa / (mu_pas * length_m_)

    flux_metric = darcy_flux_si(
        uc.md_to_m2(k_md), area_m2, uc.bar_to_pa(dp_bar), length_m, uc.cp_to_pas(mu_cp))

    # eyni fiziki ssenari FIELD girişi kimi (ft²/psi/ft) ifadə olunub
    area_ft2 = uc.convert(area_m2, "m2", "ft2", "area")
    dp_psi = uc.bar_to_psi(dp_bar)
    length_ft = uc.m_to_ft(length_m)

    flux_field = darcy_flux_si(
        uc.md_to_m2(k_md),
        uc.convert(area_ft2, "ft2", "m2", "area"),
        uc.psi_to_pa(dp_psi),
        uc.convert(length_ft, "ft", "m", "length"),
        uc.cp_to_pas(mu_cp))

    assert flux_field == pytest.approx(flux_metric, rel=1e-9)
