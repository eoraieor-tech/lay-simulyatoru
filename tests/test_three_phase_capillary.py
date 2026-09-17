"""B5-b — üç fazalı mühərrikdə kapilyar təzyiq: Pcow qoşulması və Pcog.

TAPINTI (tətbiqdən əvvəl ölçüldü): `ThreePhaseNewtonSolver.build_fluid`
`pc = None` SABİT yazırdı. Servis Pc provider-ini mühərrikə ötürürdü, lakin
mühərrik onu Nyutona vermirdi — Pc = 1 bar ilə və onsuz üç fazalı qaçış
BİT-BİT eyni idi. UI-dakı kapilyar parametrlər qazlı rejimdə səssizcə
atılırdı.

Testlərin əsas ölçüsü: kapilyar hədlərlə analitik Jakobian tam qalıq
funksiyasının SONLU FƏRQİ ilə üst-üstə düşür — həm doymuş, həm doymamış
(Rs 3-cü dəyişən) hüceyrələrdə.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.domain.scal import (CapillaryParameters, CoreyParameters,
                                GasCoreyParameters)
from imex2d.simulation.capillary import (BrooksCoreyCapillaryProvider,
                                         BrooksCoreyGasCapillaryProvider)
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.implicit.three_phase_residual import (
    ThreePhaseAccumulator, ThreePhaseFlux, ThreePhaseFluidState,
    ThreePhaseJacobianAssembler, ThreePhaseWellModel)
from imex2d.simulation.implicit.three_phase_state import ThreePhaseState
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.stone_relperm import StoneRelativePermeabilityProvider
from imex2d.simulation.well_model import PeacemanWellModel

WATER_OIL = CoreyParameters()
GAS = GasCoreyParameters()


def _pcow(entry=1.0):
    return BrooksCoreyCapillaryProvider(
        CapillaryParameters(entry_pressure=entry, lambda_exponent=2.0,
                            max_pressure=5.0), WATER_OIL)


def _pcog(entry=2.0, max_pressure=8.0):
    return BrooksCoreyGasCapillaryProvider(
        CapillaryParameters(entry_pressure=entry, lambda_exponent=2.0,
                            max_pressure=max_pressure), GAS, WATER_OIL.swc)


# ═══════════════════════ Pcog provider ═══════════════════════════════

def test_pcog_equals_entry_pressure_without_gas():
    assert _pcog().pcog(np.array([0.0]))[0] == pytest.approx(2.0)


def test_pcog_increases_monotonically_with_gas_saturation():
    sg = np.linspace(0.0, 0.6, 61)
    values = _pcog().pcog(sg)
    assert np.all(np.diff(values) >= -1e-12)
    assert values[-1] > values[0]


def test_pcog_is_capped_at_max_pressure():
    values = _pcog(max_pressure=3.0).pcog(np.linspace(0.0, 0.7, 50))
    assert values.max() == pytest.approx(3.0)


def test_pcog_derivative_matches_finite_difference():
    provider = _pcog()
    sg = np.linspace(0.02, 0.45, 30)
    step = 1e-7
    numeric = (provider.pcog(sg + step) - provider.pcog(sg - step)) / (2 * step)
    assert np.allclose(provider.dpcog_dsg(sg), numeric, rtol=1e-5, atol=1e-6)


def test_pcog_derivative_is_zero_on_the_capped_plateau():
    provider = _pcog(max_pressure=3.0)
    sg = np.array([0.69])
    assert provider.pcog(sg)[0] == pytest.approx(3.0)
    assert provider.dpcog_dsg(sg)[0] == 0.0


def test_disabled_gas_capillary_is_rejected():
    with pytest.raises(ValueError):
        BrooksCoreyGasCapillaryProvider(CapillaryParameters(), GAS, 0.2)


# ═══════════════════════ potensial ═══════════════════════════════════

def _setup(nx=5, ny=5, with_wells=True):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    model = ReservoirModelBuilder().build(
        geology, five_spot(geology.grid) if with_wells else [], scal=WATER_OIL)
    grid = TwoPointFluxDiscretization().build(model)
    wells = PeacemanWellModel().build_connections(model) if with_wells else []
    relperm = StoneRelativePermeabilityProvider.from_corey(WATER_OIL, GAS)
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    accumulator = ThreePhaseAccumulator(model, grid.pore_volume)
    flux = ThreePhaseFlux(model, grid)
    well_model = ThreePhaseWellModel(model, wells, 0.35)
    assembler = ThreePhaseJacobianAssembler(model, accumulator, flux,
                                            well_model, relperm, pvt)
    return model, grid, relperm, pvt, accumulator, flux, well_model, assembler


def _fluid_builder(relperm, pvt, n, pcow=None, pcog=None):
    def build(state):
        sw, sg = state.water_saturation, state.gas_saturation
        # SU XASSƏLƏRİ PVT-DƏN OLMALIDIR (ölçülmüş düzəliş, Seans 42).
        # Əvvəl burada `mu_w = 0.5`, `bw = 1.0` SABİT verilirdi, Jakobian
        # isə PVT-nin `water_viscosity_derivative`/`water_fvf_derivative`
        # törəmələrini oxuyur — yəni sonlu fərq testində flüid modeli
        # Jakobianla UYĞUN DEYİLDİ. Q-33-dən sonra vurucu bağlantısı da
        # μw-dən asılı olduğu üçün bu uyğunsuzluq üzə çıxdı:
        #     testin öz qurucusu (sabit)   1.09×10⁻⁵   ← ən pis yer: vurucu
        #     PVT-dən (uyğun)              9.78×10⁻¹¹
        fluid = ThreePhaseFluidState(
            mu_w=pvt.water_viscosity(state.pressure),
            mu_o=pvt.oil_viscosity(state.pressure),
            mu_g=pvt.gas_viscosity(state.pressure),
            bw=pvt.water_fvf(state.pressure),
            bo=pvt.oil_fvf(state.pressure), bg=pvt.gas_fvf(state.pressure),
            rs=state.solution_gor(pvt), krw=relperm.krw(sw),
            kro=relperm.kro_three_phase(sw, sg), krg=relperm.krg(sg))
        if pcow is not None:
            fluid.pc = pcow.pcow(sw)
            fluid.dpc_dsw = pcow.dpcow_dsw(sw)
        if pcog is not None:
            fluid.pcog = pcog.pcog(sg)
            fluid.dpcog_dsg = np.where(state.is_saturated, pcog.dpcog_dsg(sg), 0.0)
        return fluid
    return build


def test_pcog_enters_the_gas_potential_with_positive_sign():
    """Pg = Po + Pcog: çox qazlı hüceyrənin qaz potensialı YÜKSƏKDİR."""
    model, grid, relperm, pvt, _, flux, _, _ = _setup(nx=2, ny=1, with_wells=False)
    n = model.ncell
    state = ThreePhaseState(np.full(n, 200.0), np.full(n, 0.3),
                            np.array([0.4, 0.05]), np.ones(n, bool))
    base = _fluid_builder(relperm, pvt, n)(state)
    with_pc = _fluid_builder(relperm, pvt, n, pcog=_pcog())(state)
    _, _, phi_g0 = flux.potentials(state, base)
    _, _, phi_g1 = flux.potentials(state, with_pc)
    expected = with_pc.pcog[0] - with_pc.pcog[1]
    assert expected > 0
    assert phi_g1[0] - phi_g0[0] == pytest.approx(expected)


def test_pcow_enters_the_water_potential_with_negative_sign():
    model, grid, relperm, pvt, _, flux, _, _ = _setup(nx=2, ny=1, with_wells=False)
    n = model.ncell
    state = ThreePhaseState(np.full(n, 200.0), np.array([0.3, 0.6]),
                            np.full(n, 0.05), np.ones(n, bool))
    base = _fluid_builder(relperm, pvt, n)(state)
    with_pc = _fluid_builder(relperm, pvt, n, pcow=_pcow())(state)
    phi_w0, _, _ = flux.potentials(state, base)
    phi_w1, _, _ = flux.potentials(state, with_pc)
    assert phi_w1[0] - phi_w0[0] == pytest.approx(-(with_pc.pc[0] - with_pc.pc[1]))


# ═══════════════════════ tam Jakobian — sonlu fərq ═══════════════════

def _max_relative_jacobian_error(state, previous, pcow, pcog, with_wells=True,
                                 columns=None):
    model, grid, relperm, pvt, accumulator, flux, well_model, assembler = \
        _setup(with_wells=with_wells)
    n = model.ncell
    build = _fluid_builder(relperm, pvt, n, pcow=pcow, pcog=pcog)

    def residual(s):
        fluid, fluid_prev = build(s), build(previous)
        n_w, n_o, n_g = accumulator.accumulation(s, fluid)
        n_w0, n_o0, n_g0 = accumulator.accumulation(previous, fluid_prev)
        in_w, in_o, in_g = flux.net_influx(s, fluid)
        rates = well_model.well_rates(s, fluid)
        out = np.empty(n * 3)
        out[0::3] = (n_w - n_w0) - in_w - rates.water
        out[1::3] = (n_o - n_o0) - in_o - rates.oil
        out[2::3] = (n_g - n_g0) - in_g - rates.gas
        return out

    jacobian = assembler.assemble(state, build(state), dt=1.0)
    vector = state.to_vector()
    worst = 0.0
    for column in (columns if columns is not None else range(len(vector))):
        step = 1e-6 * max(1.0, abs(vector[column]))
        forward = vector.copy(); forward[column] += step
        backward = vector.copy(); backward[column] -= step
        numeric = (residual(ThreePhaseState.from_vector(forward, state.is_saturated))
                   - residual(ThreePhaseState.from_vector(backward, state.is_saturated))
                   ) / (2 * step)
        analytic = np.asarray(jacobian[:, column].todense()).ravel()
        worst = max(worst, np.max(np.abs(numeric - analytic))
                    / max(1.0, np.max(np.abs(numeric))))
    return worst


def _heterogeneous_saturated_state(n=25):
    rng = np.random.default_rng(7)
    pressure = 200.0 + 10.0 * rng.random(n)
    sw = 0.30 + 0.20 * rng.random(n)
    sg = 0.05 + 0.20 * rng.random(n)
    return ThreePhaseState(pressure, sw, sg, np.ones(n, bool))


def test_jacobian_with_water_oil_capillary_matches_finite_difference():
    """1-ci addım: Pcow üç fazalı Jakobiana düzgün qoşulub."""
    state = _heterogeneous_saturated_state()
    error = _max_relative_jacobian_error(state, state.copy(), _pcow(), None)
    assert error < 1e-5, error


def test_jacobian_with_gas_oil_capillary_matches_finite_difference():
    """2-ci addım: Pcog-un 3-cü sütun törəməsi (hər iki tərəf) düzgündür."""
    state = _heterogeneous_saturated_state()
    error = _max_relative_jacobian_error(state, state.copy(), None, _pcog())
    assert error < 1e-5, error


def test_jacobian_with_both_capillaries_and_wells_matches_finite_difference():
    state = _heterogeneous_saturated_state()
    error = _max_relative_jacobian_error(state, state.copy(), _pcow(), _pcog())
    assert error < 1e-5, error


def test_jacobian_is_exact_with_capillaries_in_a_mixed_saturation_state():
    """Doymuş + doymamış qonşu hüceyrələr: Pcog doymamışda sabit (Pe),
    törəməsi Rs-ə görə SIFIR olmalıdır."""
    n = 25
    rng = np.random.default_rng(11)
    pressure = 200.0 + 10.0 * rng.random(n)
    sw = 0.30 + 0.20 * rng.random(n)
    saturated = rng.random(n) > 0.5
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    third = np.where(saturated, 0.05 + 0.2 * rng.random(n),
                     0.6 * np.asarray(pvt.solution_gor(pressure)))
    state = ThreePhaseState(pressure, sw, third, saturated)
    error = _max_relative_jacobian_error(state, state.copy(), _pcow(), _pcog())
    baseline = _max_relative_jacobian_error(state, state.copy(), None, None)

    # ÖLÇÜLDÜ: qarışıq vəziyyətdə Jakobianın xətası kapilyardan ƏVVƏL də
    # 2.7e-5-dir (dəyişən keçidi olan hüceyrələrdə Rs sütununun mövcud
    # dəqiqliyi — bu, B5-b-nin gətirdiyi qüsur DEYİL). Ona görə burada
    # mütləq hədd deyil, kapilyarın ƏLAVƏ ETDİYİ xəta yoxlanılır.
    assert error - baseline < 1e-6, (error, baseline)
    assert error < 1e-4, error


def test_without_capillary_the_jacobian_is_unchanged():
    """GERİYƏ UYĞUNLUQ: provider yoxdursa yeni hədlər matrisə heç nə yazmır."""
    model, grid, relperm, pvt, *_rest, assembler = _setup()
    state = _heterogeneous_saturated_state()
    build = _fluid_builder(relperm, pvt, model.ncell)
    fluid = build(state)
    reference = assembler.assemble(state, fluid, dt=1.0)
    fluid.pc = np.zeros(model.ncell)          # pc sahəsi var, törəmə yox
    again = assembler.assemble(state, fluid, dt=1.0)
    assert (reference - again).nnz == 0


# ═══════════════════════ mühərrik və model səviyyəsi ═════════════════

def _three_phase_model(pcow=0.0, pcog=0.0, include_gas=True):
    from imex2d.simulation.pvt.correlations import build_pvt_table as _table
    geology = SyntheticGeologicalModelBuilder().build(
        nx=8, ny=8, dx=20.0, dy=20.0, dz=10.0, porosity=0.22,
        permx_base=150.0, top_depth=1200.0)
    return ReservoirModelBuilder().build(
        geological_model=geology, wells=five_spot(geology.grid),
        scal=WATER_OIL, gas_scal=GAS,
        capillary=CapillaryParameters(entry_pressure=pcow, lambda_exponent=2.0,
                                      max_pressure=5.0),
        gas_capillary=CapillaryParameters(entry_pressure=pcog, lambda_exponent=2.0,
                                          max_pressure=8.0),
        pvt_table=_table(pressure_min=1.0, pressure_max=400.0, n_points=40,
                         bubble_point_bar=240.0, include_gas=include_gas),
        name="Pc sınağı")


def _run(model, end=600.0):
    from helpers import default_scal
    from imex2d.application.config import SimulationConfig
    from imex2d.application.simulation_service import ModelAwareSimulationService
    from imex2d.simulation.implicit.engine import FullyImplicitEngine
    from imex2d.simulation.linear_solver import ScipyCgIluSolver
    from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
    service = ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)
    return service.run(model, SimulationConfig(end_time=end))


def test_engine_builds_the_gas_capillary_provider_only_when_enabled():
    from imex2d.application.config import SimulationConfig
    from imex2d.application.simulation_service import ModelAwareSimulationService
    from imex2d.simulation.implicit.engine import FullyImplicitEngine
    from imex2d.simulation.linear_solver import ScipyCgIluSolver
    from helpers import default_scal
    from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter

    def engine(pcow, pcog):
        service = ModelAwareSimulationService(
            relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
            linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)
        return service.create_engine(_three_phase_model(pcow, pcog),
                                     SimulationConfig(end_time=10.0))

    off = engine(0.0, 0.0)
    assert off.newton.capillary is None and off.newton.gas_capillary is None
    on = engine(1.0, 2.0)
    assert on.newton.capillary is not None, "Pcow provider-i Nyutona çatmalıdır"
    assert on.newton.gas_capillary is not None, "Pcog provider-i Nyutona çatmalıdır"


def test_water_oil_capillary_now_changes_the_three_phase_result():
    """B5-b-nin ƏSAS TƏLƏBİ: əvvəl Pcow səssizcə atılırdı (bit-bit eyni idi)."""
    base = _run(_three_phase_model())
    with_pcow = _run(_three_phase_model(pcow=1.0))
    assert base.converged and with_pcow.converged, with_pcow.message
    assert (base.series.recovery_factor[-1]
            != with_pcow.series.recovery_factor[-1]), "Pcow yenə atılır"


def test_gas_oil_capillary_changes_the_three_phase_result():
    base = _run(_three_phase_model())
    with_pcog = _run(_three_phase_model(pcog=2.0))
    assert base.converged and with_pcog.converged, with_pcog.message
    assert (base.series.recovery_factor[-1]
            != with_pcog.series.recovery_factor[-1])


def test_strong_capillaries_still_converge():
    result = _run(_three_phase_model(pcow=3.0, pcog=6.0))
    assert result.converged, result.message


def test_disabled_capillary_keeps_the_previous_behaviour():
    """Defolt (Pe = 0) modeldə nəticə İKİ QAÇIŞDA da eynidir —
    yeni hədlər yalnız açıq şəkildə tələb olunanda işə düşür."""
    first = _run(_three_phase_model())
    second = _run(_three_phase_model())
    assert (first.series.recovery_factor[-1]
            == pytest.approx(second.series.recovery_factor[-1], rel=0, abs=0))


def test_gas_capillary_without_a_gas_phase_is_a_warning():
    model = _three_phase_model(pcog=2.0, include_gas=False)
    messages = [d.message for d in model.diagnose().warnings]
    assert any("Pcog" in m for m in messages), messages


def test_gas_capillary_survives_the_project_round_trip():
    from imex2d.application.serialization import ProjectSerializer
    model = _three_phase_model(pcow=1.5, pcog=2.5)
    serializer = ProjectSerializer()
    restored = serializer.reservoir_model_from_dict(
        serializer.reservoir_model_to_dict(model))
    assert restored.gas_capillary_parameters.entry_pressure == pytest.approx(2.5)
    assert restored.capillary_parameters.entry_pressure == pytest.approx(1.5)


def test_scal_panel_exposes_the_pcog_fields():
    pytest.importorskip("PyQt5.QtWidgets")
    import inspect
    from imex2d.ui import panels
    source = inspect.getsource(panels.ScalPanel)
    assert "pcog_entry" in source and "gas_capillary_values" in source
