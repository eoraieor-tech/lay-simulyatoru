"""Analitik Jakobian (A6, mərhələ 2).

Əsas yoxlama üsulu: analitik Jakobianı mərkəzi sonlu fərq ilə
müqayisə etmək. Jakobian səhv olsa, Nyuton yığılmaz — ona görə bu
testlər A6-nın qalan mərhələlərinin təməlidir.
"""

import dataclasses

import numpy as np

from helpers import default_scal, five_spot_model
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.domain.scal import CapillaryParameters, CoreyParameters
from imex2d.domain.wells import ControlMode, WellControl
from imex2d.simulation.capillary import BrooksCoreyCapillaryProvider
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.implicit.derivatives import DerivativeProvider
from imex2d.simulation.implicit.jacobian import JacobianAssembler
from imex2d.simulation.implicit.residual import ResidualAssembler
from imex2d.simulation.implicit.state import (PRESSURE, VARIABLES_PER_CELL,
                                              WATER_SATURATION, ReservoirState)
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.well_model import PeacemanWellModel

# Cazibə üzvündə ∂ρ/∂p buraxılıb (standart sadələşdirmə). Sıxlıq
# ρ = ρ_səth / B(p) olduğuna görə bu, dərinlik fərqi olan HƏR modeldə
# özünü göstərir — həm PVT ilə, həm sıxılma modeli ilə.
EXACT = 1e-8
WITH_GRAVITY = 1e-3


def _pvt():
    return BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0))


def _capillary(scal):
    return BrooksCoreyCapillaryProvider(
        CapillaryParameters(entry_pressure=0.4), scal)


def _assemblers(model, scal, pvt=None, capillary=None):
    grid = TwoPointFluxDiscretization().build(model)
    wells = PeacemanWellModel().build_connections(model)
    residual = ResidualAssembler(model, grid, wells,
                                 CoreyRelativePermeabilityAdapter(scal),
                                 pvt=pvt, capillary=capillary)
    return residual, JacobianAssembler(residual)


def _random_state(model, scal, seed=1):
    rng = np.random.default_rng(seed)
    previous = ReservoirState(np.full(model.ncell, 250.0),
                              np.full(model.ncell, 0.40))
    state = ReservoirState(
        250.0 + rng.normal(0.0, 6.0, model.ncell),
        np.clip(0.45 + rng.normal(0.0, 0.06, model.ncell),
                scal.swc + 0.02, 1.0 - scal.sor - 0.02))
    return state, previous


def _max_relative_error(model, scal, pvt=None, capillary=None, dt=5.0, seed=1):
    residual_assembler, jacobian = _assemblers(model, scal, pvt, capillary)
    state, previous = _random_state(model, scal, seed)
    _, fluid, _ = residual_assembler.residual(state, previous, dt)
    analytic = jacobian.assemble(state, fluid, dt).toarray()
    numeric = jacobian.numerical(state, previous, dt)
    scale = np.maximum(np.abs(numeric).max(axis=0), 1e-30)
    return float(np.max(np.abs(analytic - numeric) / scale))


def _dipping(model, dip=6.0):
    surface = np.add.outer(np.zeros(model.grid.ny),
                           np.arange(model.grid.nx) * dip) + 2000.0
    model.geometry = dataclasses.replace(model.geometry,
                                         top_depth_map=surface)
    return model


# ── analitik törəmələr ────────────────────────────────────────────────
def test_corey_derivatives_match_finite_differences():
    scal = CoreyParameters()
    adapter = CoreyRelativePermeabilityAdapter(scal)
    sw = np.linspace(scal.swc + 0.02, 1.0 - scal.sor - 0.02, 15)
    step = 1e-6
    numeric_w = (adapter.krw(sw + step) - adapter.krw(sw - step)) / (2 * step)
    numeric_o = (adapter.kro(sw + step) - adapter.kro(sw - step)) / (2 * step)
    assert np.max(np.abs(adapter.krw_derivative(sw) - numeric_w)) < 1e-6
    assert np.max(np.abs(adapter.kro_derivative(sw) - numeric_o)) < 1e-6


def test_relative_permeability_derivatives_have_expected_signs():
    scal = CoreyParameters()
    adapter = CoreyRelativePermeabilityAdapter(scal)
    sw = np.linspace(scal.swc + 0.02, 1.0 - scal.sor - 0.02, 10)
    assert np.all(adapter.krw_derivative(sw) > 0)     # krw artır
    assert np.all(adapter.kro_derivative(sw) < 0)     # kro azalır


def test_derivatives_vanish_outside_the_mobile_range():
    scal = CoreyParameters()
    adapter = CoreyRelativePermeabilityAdapter(scal)
    outside = np.array([scal.swc - 0.05, 1.0 - scal.sor + 0.05])
    assert np.allclose(adapter.krw_derivative(outside), 0.0)
    assert np.allclose(adapter.kro_derivative(outside), 0.0)


def test_derivative_provider_prefers_analytic_over_numeric():
    scal = CoreyParameters()
    adapter = CoreyRelativePermeabilityAdapter(scal)
    provider = DerivativeProvider(adapter)
    sw = np.linspace(scal.swc + 0.05, 1.0 - scal.sor - 0.05, 8)
    assert np.allclose(provider.dkrw_dsw(sw), adapter.krw_derivative(sw))


def test_derivative_provider_falls_back_to_finite_differences():
    """Provider analitik törəmə verməsə də işləməlidir."""
    class _NoDerivatives:
        def __init__(self, inner):
            self._inner = inner

        def krw(self, sw, region=None):
            return self._inner.krw(sw)

        def kro(self, sw, region=None):
            return self._inner.kro(sw)

    scal = CoreyParameters()
    adapter = CoreyRelativePermeabilityAdapter(scal)
    provider = DerivativeProvider(_NoDerivatives(adapter))
    sw = np.linspace(scal.swc + 0.05, 1.0 - scal.sor - 0.05, 8)
    assert np.allclose(provider.dkrw_dsw(sw), adapter.krw_derivative(sw),
                       atol=1e-5)


def test_pvt_derivatives_are_zero_without_a_provider():
    provider = DerivativeProvider(
        CoreyRelativePermeabilityAdapter(CoreyParameters()))
    pressure = np.array([200.0, 250.0])
    assert np.allclose(provider.dbo_dp(pressure), 0.0)
    assert np.allclose(provider.dmuo_dp(pressure), 0.0)


def test_pvt_derivatives_have_expected_signs_below_bubble_point():
    provider = DerivativeProvider(
        CoreyRelativePermeabilityAdapter(CoreyParameters()), pvt=_pvt())
    pressure = np.array([150.0, 200.0])
    assert np.all(provider.dbo_dp(pressure) > 0)      # Bo Pb-yə qədər artır
    assert np.all(provider.dmuo_dp(pressure) < 0)     # μo azalır


# ── Jakobianın strukturu ──────────────────────────────────────────────
def test_jacobian_shape_and_sparsity():
    scal = default_scal()
    model = five_spot_model(nx=5, ny=5, scal=scal)
    residual_assembler, jacobian = _assemblers(model, scal)
    state, previous = _random_state(model, scal)
    _, fluid, _ = residual_assembler.residual(state, previous, 5.0)
    matrix = jacobian.assemble(state, fluid, 5.0)

    size = model.ncell * VARIABLES_PER_CELL
    assert matrix.shape == (size, size)
    assert matrix.nnz < size * size          # seyrək olmalıdır
    assert matrix.nnz > 0


def test_jacobian_pattern_is_reused_between_assemblies():
    scal = default_scal()
    model = five_spot_model(nx=5, ny=5, scal=scal)
    residual_assembler, jacobian = _assemblers(model, scal)
    state, previous = _random_state(model, scal)
    _, fluid, _ = residual_assembler.residual(state, previous, 5.0)
    first = jacobian.assemble(state, fluid, 5.0)
    second = jacobian.assemble(state, fluid, 5.0)
    assert first is second, "Matris hər yığımda yenidən yaradılır"


def test_accumulation_diagonal_has_expected_sign():
    """∂R_su/∂Sw = PV/(Δt·Bw) > 0,  ∂R_neft/∂Sw = −PV/(Δt·Bo) < 0."""
    scal = default_scal()
    model = five_spot_model(nx=4, ny=4, scal=scal)
    residual_assembler, jacobian = _assemblers(model, scal)
    state, previous = _random_state(model, scal)
    dt = 10.0
    _, fluid, _ = residual_assembler.residual(state, previous, dt)
    matrix = jacobian.assemble(state, fluid, dt).toarray()

    cell = model.ncell // 2
    water_row = cell * VARIABLES_PER_CELL + 0
    oil_row = cell * VARIABLES_PER_CELL + 1
    column = cell * VARIABLES_PER_CELL + WATER_SATURATION
    assert matrix[water_row, column] > 0
    assert matrix[oil_row, column] < 0


# ── sonlu fərqlə müqayisə (əsas yoxlama) ──────────────────────────────
def test_matches_finite_differences_base_case():
    scal = default_scal()
    assert _max_relative_error(five_spot_model(nx=5, ny=5, scal=scal),
                               scal) < EXACT


def test_matches_finite_differences_with_capillary_pressure():
    scal = default_scal()
    error = _max_relative_error(five_spot_model(nx=5, ny=5, scal=scal), scal,
                                capillary=_capillary(scal))
    assert error < EXACT


def test_matches_finite_differences_with_pvt():
    scal = default_scal()
    error = _max_relative_error(five_spot_model(nx=5, ny=5, scal=scal), scal,
                                pvt=_pvt())
    assert error < EXACT


def test_matches_finite_differences_with_pvt_and_capillary():
    scal = default_scal()
    error = _max_relative_error(five_spot_model(nx=5, ny=5, scal=scal), scal,
                                pvt=_pvt(), capillary=_capillary(scal))
    assert error < EXACT


def test_gravity_stays_within_documented_tolerance():
    """Cazibə üzvündə ∂ρ/∂p buraxılıb — sənaye praktikasında standartdır.

    Xəta 1e-3-dən kiçikdir, yəni Nyuton hələ də sürətlə yığılır:
    konvergensiya QALIĞA görə yoxlanılır, Jakobiana görə yox.
    """
    scal = default_scal()
    model = _dipping(five_spot_model(nx=5, ny=5, scal=scal))
    assert _max_relative_error(model, scal) < WITH_GRAVITY
    assert _max_relative_error(model, scal, pvt=_pvt()) < WITH_GRAVITY


def test_horizontal_model_jacobian_is_exact():
    """Dərinlik fərqi olmayanda cazibə üzvü yoxdur — Jakobian dəqiqdir."""
    scal = default_scal()
    assert _max_relative_error(five_spot_model(nx=5, ny=5, scal=scal),
                               scal, pvt=_pvt()) < EXACT


def test_matches_finite_differences_in_three_dimensions():
    """3D-də təbəqələr arasında dərinlik fərqi var -> cazibə toleransı."""
    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=4, ny=4, dx=25.0, dy=25.0, dz=5.0, porosity=0.2,
        permx_base=150.0, nz=3, top_depth=2000.0)
    model = ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                          scal=scal)
    assert _max_relative_error(model, scal,
                               capillary=_capillary(scal)) < WITH_GRAVITY


def test_matches_finite_differences_with_rate_controlled_wells():
    scal = default_scal()
    model = five_spot_model(nx=5, ny=5, scal=scal)
    for well in model.wells:
        well.control = WellControl(ControlMode.RATE,
                                   40.0 if well.is_injector else 30.0)
    assert _max_relative_error(model, scal, pvt=_pvt()) < EXACT


def test_matches_finite_differences_for_several_random_states():
    scal = default_scal()
    model = five_spot_model(nx=4, ny=4, scal=scal)
    for seed in (2, 3, 4):
        error = _max_relative_error(model, scal, pvt=_pvt(),
                                    capillary=_capillary(scal), seed=seed)
        assert error < EXACT, f"seed={seed}: {error:.2e}"


def test_matches_finite_differences_for_several_time_steps():
    scal = default_scal()
    model = five_spot_model(nx=4, ny=4, scal=scal)
    for dt in (0.1, 30.0):
        assert _max_relative_error(model, scal, dt=dt) < EXACT
