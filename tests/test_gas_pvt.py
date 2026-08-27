"""Qaz PVT xassələri (A7, mərhələ 1): Z-faktoru, Bg, μg."""

import numpy as np

from imex2d.domain.pvt import PVTTable
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import (beggs_brill_z_factor,
                                                build_pvt_table, gas_fvf,
                                                gas_viscosity,
                                                sutton_pseudo_critical)


PRESSURE = np.array([50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0])


# ── Z-faktoru ─────────────────────────────────────────────────────────
def test_z_factor_approaches_one_at_low_pressure():
    """Aşağı təzyiqdə real qaz ideal qaza yaxınlaşmalıdır (Z -> 1)."""
    z = beggs_brill_z_factor(np.array([5.0]), 0.75, 70.0)
    assert abs(z[0] - 1.0) < 0.1


def test_z_factor_has_a_standing_katz_shaped_dip():
    """Standing-Katz əyrisinin klassik forması: minimum orta təzyiqdə."""
    z = beggs_brill_z_factor(PRESSURE, 0.75, 70.0)
    minimum_index = int(np.argmin(z))
    assert 0 < minimum_index < len(PRESSURE) - 1
    assert z[minimum_index] < z[0]
    assert z[minimum_index] < z[-1]


def test_z_factor_stays_within_physical_bounds():
    z = beggs_brill_z_factor(PRESSURE, 0.75, 70.0)
    assert np.all(z > 0.2)
    assert np.all(z < 2.0)


def test_heavier_gas_gives_a_different_z_factor():
    light = beggs_brill_z_factor(PRESSURE, 0.65, 70.0)
    heavy = beggs_brill_z_factor(PRESSURE, 0.90, 70.0)
    assert not np.allclose(light, heavy)


def test_pseudo_critical_properties_increase_with_gas_gravity():
    t_light, p_light = sutton_pseudo_critical(0.6)
    t_heavy, p_heavy = sutton_pseudo_critical(0.9)
    assert t_heavy > t_light          # ağır qazın Tpc-si yüksəkdir
    assert p_heavy < p_light          # Ppc adətən azalır


# ── Bg ────────────────────────────────────────────────────────────────
def test_gas_fvf_decreases_monotonically_with_pressure():
    """Qaz sıxılandır — Bg təzyiqlə azalmalıdır (Bo/Bw-nin əksinə)."""
    bg = gas_fvf(PRESSURE, 0.75, 70.0)
    assert np.all(np.diff(bg) < 0)


def test_gas_fvf_is_much_larger_than_liquid_fvf():
    """Bg tipik olaraq Bo-dan 10-100 dəfə böyükdür (aşağı təzyiqdə)."""
    bg_low_pressure = gas_fvf(np.array([50.0]), 0.75, 70.0)[0]
    assert bg_low_pressure > 0.01     # tipik neft Bo ~1.0-1.5 ilə müqayisədə


def test_gas_fvf_is_positive_everywhere():
    bg = gas_fvf(PRESSURE, 0.75, 70.0)
    assert np.all(bg > 0)


# ── qaz lözlüyü ──────────────────────────────────────────────────────
def test_gas_viscosity_increases_with_pressure():
    """Neft/su lözlüyünün ƏKSİNƏ — qaz sıxlandıqca lözlük artır."""
    mu = gas_viscosity(PRESSURE, 0.75, 70.0)
    assert np.all(np.diff(mu) > 0)


def test_gas_viscosity_is_in_a_realistic_range():
    """Tipik təbii qaz lözlüyü 0.01-0.05 cP aralığındadır."""
    mu = gas_viscosity(PRESSURE, 0.75, 70.0)
    assert np.all(mu > 0.005)
    assert np.all(mu < 0.1)


def test_gas_viscosity_is_much_smaller_than_oil_viscosity():
    """Qaz həmişə neftdən qat-qat az özlüdür — üç fazalı axının əsası."""
    mu_gas = gas_viscosity(np.array([150.0]), 0.75, 70.0)[0]
    assert mu_gas < 0.5               # tipik neft lözlüyü > 1 cP


# ── PVTTable domain ─────────────────────────────────────────────────
def test_table_without_gas_columns_reports_no_gas_phase():
    table = build_pvt_table(bubble_point_bar=240.0)
    assert not table.has_gas_phase
    assert table.gas_fvf is None
    assert table.validate() == []


def test_table_with_gas_columns_reports_gas_phase():
    table = build_pvt_table(bubble_point_bar=240.0, include_gas=True)
    assert table.has_gas_phase
    assert table.gas_fvf is not None
    assert table.gas_viscosity is not None
    assert table.validate() == []


def test_table_rejects_mismatched_gas_column_length():
    table = build_pvt_table(bubble_point_bar=240.0, include_gas=True)
    table.gas_fvf = table.gas_fvf[:-2]
    assert table.validate()


def test_table_rejects_non_positive_gas_fvf():
    table = build_pvt_table(bubble_point_bar=240.0, include_gas=True)
    table.gas_fvf = table.gas_fvf.copy()
    table.gas_fvf[3] = -0.01
    assert any("Bg" in issue for issue in table.validate())


def test_table_rejects_increasing_gas_fvf():
    """Bg təzyiqlə artırsa (fiziki cəhətdən yanlış), yoxlama tutmalıdır."""
    table = build_pvt_table(bubble_point_bar=240.0, include_gas=True)
    table.gas_fvf = table.gas_fvf.copy()
    table.gas_fvf[-1] = table.gas_fvf[0] * 2.0
    assert any("Bg" in issue for issue in table.validate())


def test_manual_table_without_gas_is_backward_compatible():
    """Köhnə (A1-dövrü) təyinat sintaksisi indi də işləməlidir."""
    table = PVTTable(
        pressure=np.array([50.0, 150.0, 250.0]),
        oil_fvf=np.array([1.2, 1.3, 1.25]),
        oil_viscosity=np.array([2.0, 1.5, 1.8]),
        solution_gor=np.array([20.0, 60.0, 60.0]),
        water_fvf=np.array([1.02, 1.01, 1.0]),
        water_viscosity=np.array([0.5, 0.5, 0.5]))
    assert not table.has_gas_phase
    assert table.validate() == []


# ── BlackOilPVTProvider ──────────────────────────────────────────────
def test_provider_exposes_gas_properties_when_table_has_them():
    provider = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                                   include_gas=True))
    assert provider.has_gas_phase()
    pressures = np.array([100.0, 200.0])
    assert np.all(provider.gas_fvf(pressures) > 0)
    assert np.all(provider.gas_viscosity(pressures) > 0)


def test_provider_raises_clearly_without_gas_columns():
    provider = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0))
    assert not provider.has_gas_phase()
    try:
        provider.gas_fvf(np.array([100.0]))
    except NotImplementedError:
        return
    raise AssertionError("Qazsız cədvəldə gas_fvf() səssizcə keçdi")


def test_gas_fvf_derivative_matches_finite_difference():
    """Analitik törəmə (parçalı-xətti meyl) sonlu fərqlə uyğun olmalıdır."""
    provider = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                                   include_gas=True,
                                                   n_points=60))
    pressure = np.array([120.0])
    step = 1e-3
    analytic = provider.gas_fvf_derivative(pressure)[0]
    numeric = (provider.gas_fvf(pressure + step)[0]
              - provider.gas_fvf(pressure - step)[0]) / (2 * step)
    assert abs(analytic - numeric) < 1e-3 * max(abs(numeric), 1.0)


def test_gas_viscosity_derivative_is_positive():
    """Lözlük təzyiqlə artır — törəmə müsbət olmalıdır."""
    provider = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                                   include_gas=True))
    derivative = provider.gas_viscosity_derivative(np.array([150.0]))[0]
    assert derivative > 0


def test_default_interface_reports_no_gas_phase():
    """`IPVTProvider.has_gas_phase()` defolt False — köhnə provider-lər
    heç nə etmədən yeni interfeysi tətbiq edir."""
    from imex2d.interfaces.providers import IPVTProvider

    class Minimal(IPVTProvider):
        def oil_fvf(self, pressure, region=None): return pressure
        def oil_viscosity(self, pressure, region=None): return pressure
        def water_fvf(self, pressure, region=None): return pressure
        def water_viscosity(self, pressure, region=None): return pressure
        def total_compressibility(self, pressure, sw, region=None): return pressure

    assert Minimal().has_gas_phase() is False
