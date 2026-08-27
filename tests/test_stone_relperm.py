"""Stone II üç fazalı nisbi keçiricilik (A7, mərhələ 3)."""

import numpy as np

from imex2d.domain.scal import CoreyParameters, GasCoreyParameters
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.stone_relperm import StoneRelativePermeabilityProvider


def _stone(**gas_kwargs):
    wo = CoreyParameters()
    go = GasCoreyParameters(**gas_kwargs)
    return StoneRelativePermeabilityProvider.from_corey(wo, go), wo, go


# ── GasCoreyParameters ─────────────────────────────────────────────────
def test_krg_increases_with_gas_saturation():
    go = GasCoreyParameters()
    sg = np.linspace(0.0, 0.6, 10)
    krg = go.krg(sg, swc=0.2)
    assert np.all(np.diff(krg) >= -1e-12)


def test_krog_decreases_with_gas_saturation():
    """Qaz-neft sistemində Sg artdıqca neftin keçiriciliyi azalmalıdır."""
    go = GasCoreyParameters()
    sg = np.linspace(0.0, 0.6, 10)
    krog = go.krog(sg, swc=0.2, kro_end=0.9)
    assert np.all(np.diff(krog) <= 1e-12)


def test_krg_is_zero_below_critical_saturation():
    go = GasCoreyParameters(sgc=0.1)
    assert go.krg(np.array([0.0, 0.05]), swc=0.2)[0] == 0.0
    assert go.krg(np.array([0.0, 0.05]), swc=0.2)[1] == 0.0


def test_gas_parameters_validate_the_saturation_span():
    go = GasCoreyParameters(sgc=0.5, sorg=0.5)
    assert go.validate(swc=0.2)          # 0.2+0.5+0.5 >= 0.95 -> xəta


def test_gas_parameters_reject_bad_endpoints():
    assert GasCoreyParameters(krg_end=1.5).validate(swc=0.2)
    assert GasCoreyParameters(sgc=-0.1).validate(swc=0.2)


# ── Stone II — sərhəd hallara reduksiya (ƏSAS FİZİKİ SINAQ) ────────────
def test_stone_reduces_exactly_to_two_phase_water_oil_at_zero_gas():
    """Sg=0-da Stone iki fazalı su-neft əyrisinə DƏQİQ uyğun gəlməlidir."""
    stone, wo, _ = _stone()
    sw = np.linspace(wo.swc, 1.0 - wo.sor, 20)
    stone_kro = stone.kro_three_phase(sw, np.zeros_like(sw))
    assert np.allclose(stone_kro, wo.kro(sw), atol=1e-10)


def test_stone_reduces_exactly_to_gas_oil_curve_at_connate_water():
    """Sw=Swc-də Stone qaz-neft əyrisinə DƏQİQ uyğun gəlməlidir."""
    stone, wo, go = _stone()
    sg = np.linspace(0.0, 1.0 - wo.swc - go.sorg, 20)
    stone_kro = stone.kro_three_phase(np.full_like(sg, wo.swc), sg)
    assert np.allclose(stone_kro, go.krog(sg, wo.swc, wo.kro_end), atol=1e-10)


def test_two_phase_kro_method_matches_the_wrapped_provider():
    """`kro(sw)` (iki fazalı metod) sarğılanan provider ilə eyni olmalıdır."""
    stone, wo, _ = _stone()
    sw = np.linspace(wo.swc, 1.0 - wo.sor, 10)
    assert np.allclose(stone.kro(sw), wo.kro(sw))


# ── monotonluq və fiziki hədlər ─────────────────────────────────────────
def test_kro_decreases_as_gas_saturation_increases_at_fixed_water():
    stone, wo, go = _stone()
    sw = np.full(8, 0.5)
    sg = np.linspace(0.0, 0.4, 8)
    kro = stone.kro_three_phase(sw, sg)
    assert np.all(np.diff(kro) <= 1e-9)


def test_kro_never_negative_across_the_full_saturation_plane():
    """Stone I-in məlum qüsuru (mənfi kro) Stone II-də olmamalıdır."""
    stone, wo, go = _stone()
    sw_grid, sg_grid = np.meshgrid(np.linspace(wo.swc, 1 - wo.sor, 20),
                                   np.linspace(0.0, 0.5, 20))
    kro = stone.kro_three_phase(sw_grid.ravel(), sg_grid.ravel())
    assert np.all(kro >= -1e-12)


def test_kro_never_exceeds_the_endpoint():
    stone, wo, go = _stone()
    sw_grid, sg_grid = np.meshgrid(np.linspace(wo.swc, 1 - wo.sor, 20),
                                   np.linspace(0.0, 0.5, 20))
    kro = stone.kro_three_phase(sw_grid.ravel(), sg_grid.ravel())
    assert np.all(kro <= wo.kro_end + 1e-9)


def test_kro_decreases_as_water_saturation_increases_at_fixed_gas():
    """İki fazalı su-neft davranışı üç fazalı kəsikdə də qorunmalıdır."""
    stone, wo, go = _stone()
    sg = np.full(8, 0.1)
    sw = np.linspace(wo.swc, 1.0 - wo.sor - 0.1, 8)
    kro = stone.kro_three_phase(sw, sg)
    assert np.all(np.diff(kro) <= 1e-9)


# ── iki fazalı metodların ötürülməsi ────────────────────────────────────
def test_krw_is_forwarded_unchanged():
    stone, wo, _ = _stone()
    sw = np.linspace(wo.swc, 1.0 - wo.sor, 10)
    assert np.allclose(stone.krw(sw), wo.krw(sw))


def test_saturation_limits_are_forwarded():
    stone, wo, _ = _stone()
    assert stone.saturation_limits() == (wo.swc, 1.0 - wo.sor)


def test_cfl_derivative_is_forwarded():
    stone, wo, _ = _stone()
    adapter = CoreyRelativePermeabilityAdapter(wo)
    assert abs(stone.max_fractional_flow_derivative(0.5, 3.0)
              - adapter.max_fractional_flow_derivative(0.5, 3.0)) < 1e-9


def test_gas_saturation_limits():
    stone, wo, go = _stone(sgc=0.05, sorg=0.1)
    low, high = stone.gas_saturation_limits()
    assert abs(low - 0.05) < 1e-9
    assert abs(high - (1.0 - wo.swc - 0.1)) < 1e-9


# ── has_gas_phase / geriyə uyğunluq ─────────────────────────────────────
def test_stone_reports_gas_phase():
    stone, _, _ = _stone()
    assert stone.has_gas_phase()


def test_two_phase_adapter_reports_no_gas_phase():
    adapter = CoreyRelativePermeabilityAdapter(CoreyParameters())
    assert not adapter.has_gas_phase()


def test_two_phase_adapter_raises_clearly_for_gas_methods():
    adapter = CoreyRelativePermeabilityAdapter(CoreyParameters())
    try:
        adapter.krg(np.array([0.1]))
    except NotImplementedError:
        return
    raise AssertionError("İki fazalı adapter krg() çağırışını səssizcə keçdi")


def test_kro_end_is_taken_from_water_oil_parameters_automatically():
    """`from_corey` iki əyrini eyni son nöqtəyə bağlamalıdır."""
    wo = CoreyParameters(kro_end=0.75)
    go = GasCoreyParameters()
    stone = StoneRelativePermeabilityProvider.from_corey(wo, go)
    assert abs(stone.kro_end - 0.75) < 1e-12


def test_invalid_gas_parameters_are_rejected_at_construction():
    wo = CoreyParameters()
    go = GasCoreyParameters(sgc=0.6, sorg=0.6)     # 0.2+0.6+0.6 > 0.95
    try:
        StoneRelativePermeabilityProvider.from_corey(wo, go)
    except ValueError:
        return
    raise AssertionError("Yararsız qaz parametrləri qəbul edildi")


# ── mühərriklə inteqrasiya (qaz hələ bağlanmayıb, iki fazalı yol) ──────
def test_wrapped_provider_still_runs_a_normal_two_phase_simulation():
    """Stone provider mühərrikə İKİ FAZALI kimi verilsə problemsiz işləməlidir.

    `kro`, `krw`, CFL — hamısı ötürülür, mühərrik `kro_three_phase`-i
    hələ çağırmır (mərhələ 4-6-nın işidir).
    """
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
    from helpers import default_scal, five_spot_model, make_service, short_config

    scal = default_scal()
    stone, _, _ = _stone()
    model = five_spot_model(nx=9, ny=9, scal=scal)
    result = make_service(scal).run(model, short_config(end_time=100.0))
    # baza (Corey adapteri) ilə müqayisə — Stone-un iki fazalı yolu
    # eyni nəticəni verməlidir, çünki krw/kro birbaşa ötürülür
    from imex2d.application.simulation_service import SimulationService
    stone_result = SimulationService(relperm_provider=stone).run(
        five_spot_model(nx=9, ny=9, scal=scal), short_config(end_time=100.0))
    assert abs(result.final_recovery_factor
              - stone_result.final_recovery_factor) < 1e-6
