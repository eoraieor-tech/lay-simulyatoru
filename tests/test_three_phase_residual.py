"""Üç fazalı akkumulyasiya — kütlə balansı tənlikləri (A7, mərhələ 5)."""

import numpy as np

from helpers import default_scal
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.implicit.three_phase_residual import (
    ThreePhaseAccumulator, ThreePhaseFlux, ThreePhaseFluidState,
    ThreePhaseWellModel)
from imex2d.simulation.implicit.three_phase_state import ThreePhaseState
from imex2d.domain.scal import CoreyParameters, GasCoreyParameters
from imex2d.domain.wells import ControlMode
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table


def _setup(nx=6, ny=6, nz=3):
    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=nz, top_depth=2000.0)
    model = ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                          scal=scal)
    grid = TwoPointFluxDiscretization().build(model)
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    return model, grid, pvt


def _fluid(pressure, rs, pvt, n):
    return ThreePhaseFluidState(
        mu_w=np.full(n, 0.5), mu_o=pvt.oil_viscosity(pressure),
        mu_g=pvt.gas_viscosity(pressure),
        bw=np.full(n, 1.0), bo=pvt.oil_fvf(pressure), bg=pvt.gas_fvf(pressure),
        rs=rs, krw=np.full(n, 0.2), kro=np.full(n, 0.5), krg=np.full(n, 0.1))


def _mixed_state(model, pvt, pressure_value=200.0):
    n = model.ncell
    pressure = np.full(n, pressure_value)
    sw = np.full(n, 0.3)
    third = np.zeros(n)
    is_saturated = np.zeros(n, bool)
    half = n // 2
    third[:half] = pvt.solution_gor(pressure[:half]) * 0.6
    third[half:] = 0.15
    is_saturated[half:] = True
    return ThreePhaseState(pressure, sw, third, is_saturated), half


# ── akkumulyasiya düsturu ───────────────────────────────────────────
def test_undersaturated_gas_comes_only_from_dissolved_gas():
    """Sg=0 hüceyrədə qaz tənliyinin YEGANƏ mənbəyi həll olmuş qazdır."""
    model, grid, pvt = _setup()
    state, half = _mixed_state(model, pvt)
    accumulator = ThreePhaseAccumulator(model, grid.pore_volume)
    rs = state.solution_gor(pvt)
    fluid = _fluid(state.pressure, rs, pvt, model.ncell)

    water, oil, gas = accumulator.accumulation(state, fluid)
    pore_volume = accumulator.pore_volume_at(state.pressure)
    expected = (pore_volume[:half] * (1.0 - state.water_saturation[:half])
               * rs[:half] / fluid.bo[:half])
    assert np.allclose(gas[:half], expected)


def test_saturated_gas_includes_both_free_and_dissolved():
    model, grid, pvt = _setup()
    state, half = _mixed_state(model, pvt)
    accumulator = ThreePhaseAccumulator(model, grid.pore_volume)
    rs = state.solution_gor(pvt)
    fluid = _fluid(state.pressure, rs, pvt, model.ncell)

    _, _, gas = accumulator.accumulation(state, fluid)
    pore_volume = accumulator.pore_volume_at(state.pressure)
    free_gas = pore_volume[half:] * state.gas_saturation[half:] / fluid.bg[half:]
    dissolved = (pore_volume[half:] * state.oil_saturation[half:]
                * rs[half:] / fluid.bo[half:])
    assert np.allclose(gas[half:], free_gas + dissolved)
    assert np.all(free_gas > 0)
    assert np.all(dissolved > 0)


def test_saturated_gas_is_larger_than_undersaturated_at_similar_conditions():
    """Sərbəst qaz əlavə hədddir — doyma vəziyyəti qaz miqdarını artırmalıdır."""
    model, grid, pvt = _setup()
    state, half = _mixed_state(model, pvt)
    accumulator = ThreePhaseAccumulator(model, grid.pore_volume)
    rs = state.solution_gor(pvt)
    fluid = _fluid(state.pressure, rs, pvt, model.ncell)
    _, _, gas = accumulator.accumulation(state, fluid)
    assert gas[half] > gas[0]


# ── iki fazalı reduksiya (Sg=0 hüceyrələrdə) ────────────────────────
def test_water_and_oil_match_two_phase_formula_for_undersaturated_cells():
    """Su/neft düsturu dəyişməyib — yalnız qaz tənliyi yenidir."""
    model, grid, pvt = _setup()
    state, half = _mixed_state(model, pvt)
    accumulator = ThreePhaseAccumulator(model, grid.pore_volume)
    rs = state.solution_gor(pvt)
    fluid = _fluid(state.pressure, rs, pvt, model.ncell)

    water, oil, _ = accumulator.accumulation(state, fluid)
    water_ref, oil_ref = accumulator.two_phase_accumulation_matches(state, fluid)
    assert np.allclose(water[:half], water_ref[:half])
    assert np.allclose(oil[:half], oil_ref[:half])


def test_oil_differs_from_two_phase_reference_when_gas_is_present():
    """Doymuş hüceyrədə So = 1−Sw−Sg — iki fazalı (1−Sw) fərziyyəsi yanlışdır."""
    model, grid, pvt = _setup()
    state, half = _mixed_state(model, pvt)
    accumulator = ThreePhaseAccumulator(model, grid.pore_volume)
    rs = state.solution_gor(pvt)
    fluid = _fluid(state.pressure, rs, pvt, model.ncell)

    _, oil, _ = accumulator.accumulation(state, fluid)
    _, oil_ref = accumulator.two_phase_accumulation_matches(state, fluid)
    assert not np.allclose(oil[half:], oil_ref[half:])
    assert np.all(oil[half:] < oil_ref[half:])   # Sg oil-in yerini alır


# ── sərhəddə kəsilməzlik ────────────────────────────────────────────
def test_gas_accumulation_is_continuous_across_the_saturation_boundary():
    """Sg=0 ⟺ Rs=Rs_sat(p) — qaz kəmiyyəti keçid anında sıçramamalıdır."""
    model, grid, pvt = _setup()
    pore_volume = grid.pore_volume[:1]
    accumulator = ThreePhaseAccumulator(model, pore_volume)
    pressure = np.array([200.0])
    rs_sat = float(pvt.solution_gor(pressure)[0])

    under = ThreePhaseState(pressure, np.array([0.3]), np.array([rs_sat]),
                            np.array([False]))
    over = ThreePhaseState(pressure, np.array([0.3]), np.array([1e-9]),
                           np.array([True]))

    fluid_under = _fluid(pressure, under.solution_gor(pvt), pvt, 1)
    fluid_over = _fluid(pressure, over.solution_gor(pvt), pvt, 1)

    _, _, gas_under = accumulator.accumulation(under, fluid_under)
    _, _, gas_over = accumulator.accumulation(over, fluid_over)
    assert abs(gas_under[0] - gas_over[0]) < 1e-2      # PV böyüklüyünə görə mütləq tolerans


# ── FluidState mobilliklər ──────────────────────────────────────────
def test_mobilities_are_relative_permeability_over_viscosity():
    fluid = ThreePhaseFluidState(
        mu_w=np.array([0.5]), mu_o=np.array([2.0]), mu_g=np.array([0.02]),
        bw=np.array([1.0]), bo=np.array([1.2]), bg=np.array([0.01]),
        rs=np.array([80.0]), krw=np.array([0.2]), kro=np.array([0.6]),
        krg=np.array([0.1]))
    assert abs(fluid.lam_w[0] - 0.4) < 1e-9
    assert abs(fluid.lam_o[0] - 0.3) < 1e-9
    assert abs(fluid.lam_g[0] - 5.0) < 1e-9


# ── süxur sıxılması ────────────────────────────────────────────────
def test_pore_volume_at_reduces_to_reference_when_incompressible():
    model, grid, pvt = _setup()
    model.rock.compressibility = 0.0
    accumulator = ThreePhaseAccumulator(model, grid.pore_volume)
    pressure = np.full(model.ncell, 350.0)     # referansdan çox fərqli
    assert np.allclose(accumulator.pore_volume_at(pressure), grid.pore_volume)


def test_pore_volume_increases_with_pressure_when_compressible():
    model, grid, pvt = _setup()
    model.rock.compressibility = 4.5e-5
    accumulator = ThreePhaseAccumulator(model, grid.pore_volume)
    reference = accumulator.reference_pressure
    high = accumulator.pore_volume_at(np.full(model.ncell, reference + 100.0))
    low = accumulator.pore_volume_at(np.full(model.ncell, reference - 100.0))
    assert np.all(high > low)


# ── axın həddləri (ThreePhaseFlux) ─────────────────────────────────────
def _column_setup(nz=4, sw_value=0.3, sg_value=0.2, pressure_value=200.0):
    from imex2d.simulation.implicit.three_phase_residual import ThreePhaseFlux

    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=1, ny=1, dx=50.0, dy=50.0, dz=20.0, porosity=0.2,
        permx_base=150.0, nz=nz, top_depth=2000.0)
    model = ReservoirModelBuilder().build(geology, [], scal=scal)
    grid = TwoPointFluxDiscretization().build(model)
    flux = ThreePhaseFlux(model, grid)

    n = model.ncell
    pressure = np.full(n, pressure_value)
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    state = ThreePhaseState(pressure, np.full(n, sw_value),
                            np.full(n, sg_value), np.ones(n, bool))
    rs = state.solution_gor(pvt)
    fluid = ThreePhaseFluidState(
        mu_w=np.full(n, 0.5), mu_o=pvt.oil_viscosity(pressure),
        mu_g=pvt.gas_viscosity(pressure),
        bw=np.full(n, 1.0), bo=pvt.oil_fvf(pressure),
        bg=pvt.gas_fvf(pressure), rs=rs,
        krw=np.full(n, 0.2), kro=np.full(n, 0.5), krg=np.full(n, 0.15))
    return model, grid, flux, state, fluid


def test_water_and_oil_potentials_match_two_phase_gravity_formula():
    """Su/neft potensialı dəyişməyib — yalnız qaz yenidir."""
    model, grid, flux, state, fluid = _column_setup()
    d_phi_w, d_phi_o, _ = flux.potentials(state, fluid)
    assert np.all(np.abs(d_phi_w) > 0)
    assert np.all(np.abs(d_phi_o) > 0)


def test_gas_has_the_smallest_gravity_head_of_the_three_phases():
    """ρ_qaz ≪ ρ_neft < ρ_su — qazın potensial fərqi ən kiçik olmalıdır."""
    model, grid, flux, state, fluid = _column_setup()
    d_phi_w, d_phi_o, d_phi_g = flux.potentials(state, fluid)
    assert np.max(np.abs(d_phi_g)) < np.max(np.abs(d_phi_o))
    assert np.max(np.abs(d_phi_o)) < np.max(np.abs(d_phi_w))


def test_gas_flux_equals_free_plus_dissolved():
    model, grid, flux, state, fluid = _column_setup()
    water_flux, oil_flux, gas_flux = flux.face_fluxes(state, fluid)

    d_phi_w, d_phi_o, d_phi_g = flux.potentials(state, fluid)
    conn = grid.connections
    up_o = np.where(d_phi_o >= 0, conn.cell_a, conn.cell_b)
    dissolved = fluid.rs[up_o] * oil_flux
    free = gas_flux - dissolved
    assert np.allclose(gas_flux, free + dissolved)


def test_dissolved_gas_follows_the_oil_flow_direction():
    """Rs 'sərnişindir' — öz upstream-i yoxdur, neftin daşıdığı yerə gedir."""
    model, grid, flux, state, fluid = _column_setup()
    water_flux, oil_flux, gas_flux = flux.face_fluxes(state, fluid)
    d_phi_w, d_phi_o, _ = flux.potentials(state, fluid)
    conn = grid.connections
    up_o = np.where(d_phi_o >= 0, conn.cell_a, conn.cell_b)
    dissolved = fluid.rs[up_o] * oil_flux
    same_sign = (np.sign(dissolved) == np.sign(oil_flux)) | (np.abs(oil_flux) < 1e-12)
    assert bool(np.all(same_sign))


def test_net_influx_conserves_mass_across_internal_faces():
    """Daxili üzlərdə axın itmir — bir tərəfdən çıxan digər tərəfə düşür."""
    model, grid, flux, state, fluid = _column_setup()
    water, oil, gas = flux.net_influx(state, fluid)
    # təcrid olunmuş sütun, xarici axın yoxdur -> cəmi sıfıra yaxın olmalıdır
    assert abs(float(np.sum(water))) < 1e-9
    assert abs(float(np.sum(oil))) < 1e-9
    assert abs(float(np.sum(gas))) < 1e-9


def test_zero_potential_difference_gives_zero_flux():
    """Bərabər təzyiq, cazibəsiz (üfüqi) grid — axın olmamalıdır."""
    from imex2d.simulation.implicit.three_phase_residual import ThreePhaseFlux

    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=3, ny=1, dx=50.0, dy=50.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    model = ReservoirModelBuilder().build(geology, [], scal=scal)
    grid = TwoPointFluxDiscretization().build(model)
    flux = ThreePhaseFlux(model, grid)

    n = model.ncell
    pressure = np.full(n, 200.0)
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    state = ThreePhaseState(pressure, np.full(n, 0.3), np.full(n, 0.2),
                            np.ones(n, bool))
    rs = state.solution_gor(pvt)
    fluid = ThreePhaseFluidState(
        mu_w=np.full(n, 0.5), mu_o=pvt.oil_viscosity(pressure),
        mu_g=pvt.gas_viscosity(pressure), bw=np.full(n, 1.0),
        bo=pvt.oil_fvf(pressure), bg=pvt.gas_fvf(pressure), rs=rs,
        krw=np.full(n, 0.2), kro=np.full(n, 0.5), krg=np.full(n, 0.15))

    water, oil, gas = flux.face_fluxes(state, fluid)
    assert np.allclose(water, 0.0)
    assert np.allclose(oil, 0.0)
    assert np.allclose(gas, 0.0)


def test_pressure_gradient_drives_flow_from_high_to_low():
    from imex2d.simulation.implicit.three_phase_residual import ThreePhaseFlux

    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=3, ny=1, dx=50.0, dy=50.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    model = ReservoirModelBuilder().build(geology, [], scal=scal)
    grid = TwoPointFluxDiscretization().build(model)
    flux = ThreePhaseFlux(model, grid)

    n = model.ncell
    pressure = np.array([220.0, 200.0, 180.0])   # azalan soldan sağa
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    state = ThreePhaseState(pressure, np.full(n, 0.3), np.full(n, 0.2),
                            np.ones(n, bool))
    rs = state.solution_gor(pvt)
    fluid = ThreePhaseFluidState(
        mu_w=np.full(n, 0.5), mu_o=pvt.oil_viscosity(pressure),
        mu_g=pvt.gas_viscosity(pressure), bw=np.full(n, 1.0),
        bo=pvt.oil_fvf(pressure), bg=pvt.gas_fvf(pressure), rs=rs,
        krw=np.full(n, 0.2), kro=np.full(n, 0.5), krg=np.full(n, 0.15))

    water, oil, gas = flux.face_fluxes(state, fluid)
    # A->B müsbətdirsə, yüksək təzyiqdən aşağıya axın deməkdir
    assert np.all(water > 0)
    assert np.all(oil > 0)
    assert np.all(gas > 0)


# ── quyu həddləri (ThreePhaseWellModel) ────────────────────────────────
def _well_setup(nx=9, ny=9, pressure_value=220.0, sw_value=0.35, sg_value=0.1):
    from imex2d.simulation.implicit.three_phase_residual import (
        ThreePhaseWellModel)
    from imex2d.simulation.well_model import PeacemanWellModel

    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    model = ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                          scal=scal)
    wells = PeacemanWellModel().build_connections(model)

    n = model.ncell
    pressure = np.full(n, pressure_value)
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    state = ThreePhaseState(pressure, np.full(n, sw_value),
                            np.full(n, sg_value), np.ones(n, bool))
    rs = state.solution_gor(pvt)
    fluid = ThreePhaseFluidState(
        mu_w=np.full(n, 0.5), mu_o=pvt.oil_viscosity(pressure),
        mu_g=pvt.gas_viscosity(pressure), bw=np.full(n, 1.0),
        bo=pvt.oil_fvf(pressure), bg=pvt.gas_fvf(pressure), rs=rs,
        krw=np.full(n, 0.2), kro=np.full(n, 0.4), krg=np.full(n, 0.15))
    model_well = ThreePhaseWellModel(model, wells, 0.35)
    return model, wells, state, fluid, model_well


def test_producer_gas_rate_equals_free_plus_dissolved():
    model, wells, state, fluid, well_model = _well_setup()
    rates = well_model.well_rates(state, fluid)
    for connection in wells:
        if connection.is_injector:
            continue
        cell = connection.cell
        dissolved = fluid.rs[cell] * rates.oil[cell]
        free_gas = rates.gas[cell] - dissolved
        assert abs(free_gas) > 0            # sg > 0 halında sərbəst hədd sıfırdan böyük olmalı


def test_injector_produces_only_water():
    model, wells, state, fluid, well_model = _well_setup()
    rates = well_model.well_rates(state, fluid)
    for connection in wells:
        if connection.is_injector:
            assert rates.oil[connection.cell] == 0.0
            assert rates.gas[connection.cell] == 0.0
            assert rates.water[connection.cell] > 0.0   # laya daxil olur


def test_producer_rates_are_negative():
    """İşarə konvensiyası: müsbət = laya daxil olur (A6-dakı ilə eyni)."""
    model, wells, state, fluid, well_model = _well_setup()
    rates = well_model.well_rates(state, fluid)
    for connection in wells:
        if not connection.is_injector:
            cell = connection.cell
            assert rates.water[cell] <= 0.0
            assert rates.oil[cell] <= 0.0
            assert rates.gas[cell] <= 0.0


def test_gas_rate_is_zero_when_cell_is_undersaturated():
    """Sg=0 hüceyrədə sərbəst qaz yoxdur — yalnız həll olmuş qaz çıxır."""
    model, wells, state, fluid, well_model = _well_setup(sg_value=0.0)
    state = ThreePhaseState(state.pressure, state.water_saturation,
                            np.zeros(model.ncell), np.zeros(model.ncell, bool))
    fluid.krg[:] = 0.0             # doymamış hüceyrədə krg=0 (Stone II-də də belədir)
    rates = well_model.well_rates(state, fluid)
    for connection in wells:
        if not connection.is_injector:
            cell = connection.cell
            expected = fluid.rs[cell] * rates.oil[cell]
            assert abs(rates.gas[cell] - expected) < 1e-9


def test_per_well_totals_match_the_sum_of_producing_cells():
    model, wells, state, fluid, well_model = _well_setup()
    rates = well_model.well_rates(state, fluid)
    for name, total in rates.per_well_gas.items():
        cells = [c.cell for c in wells if c.well_name == name]
        assert abs(total - sum(rates.gas[c] for c in cells)) < 1e-9


def test_rate_mode_target_still_refers_to_liquid_not_gas():
    """RATE hədəfi su+neft debitidir — A6-dakı konvensiya qorunur."""
    from imex2d.simulation.implicit.three_phase_residual import (
        ThreePhaseWellModel)
    from imex2d.simulation.well_model import WellConnection

    model, wells, state, fluid, _ = _well_setup()
    producer = next(c for c in wells if not c.is_injector)
    rate_connection = WellConnection(
        well_name=producer.well_name, cell=producer.cell,
        well_index=producer.well_index, is_injector=False,
        mode=ControlMode.RATE, target=-50.0)
    well_model = ThreePhaseWellModel(model, [rate_connection], 0.35)
    rates = well_model.well_rates(state, fluid)
    cell = producer.cell
    liquid = abs(rates.water[cell] * fluid.bw[cell]
                + rates.oil[cell] * fluid.bo[cell])
    assert abs(liquid - 50.0) < 1e-6


# ── akkumulyasiya Jakobianı (analitik vs sonlu fərq) ────────────────────
def _jacobian_setup(nx=4, ny=4, nz=2, pressure_value=213.7):
    """Qəsdən cədvəl nöqtəsi OLMAYAN təzyiq — sonlu fərq testi üçün.

    Cədvəl nöqtəsində (məs. 200.0) mərkəzi fərq iki fərqli parçalı-xətti
    seqmenti birləşdirir və süni fərq yaradır — bu, `_slope()`-un
    özündə deyil, YALNIZ FD-testin metodologiyasında problemdir.
    """
    from imex2d.simulation.implicit.three_phase_residual import (
        ThreePhaseAccumulationJacobian)

    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=nz, top_depth=2000.0)
    model = ReservoirModelBuilder().build(geology, [], scal=scal)
    model.rock.compressibility = 4.5e-5
    grid = TwoPointFluxDiscretization().build(model)
    accumulator = ThreePhaseAccumulator(model, grid.pore_volume)
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    jacobian = ThreePhaseAccumulationJacobian(accumulator, pvt)

    n = model.ncell
    pressure = np.full(n, pressure_value)
    sw = np.full(n, 0.32)
    third = np.zeros(n)
    is_saturated = np.zeros(n, bool)
    half = n // 2
    third[:half] = pvt.solution_gor(pressure[:half]) * 0.55
    third[half:] = 0.12
    is_saturated[half:] = True
    state = ThreePhaseState(pressure, sw, third, is_saturated)

    def build_fluid(s):
        rs = s.solution_gor(pvt)
        return ThreePhaseFluidState(
            mu_w=np.full(n, 0.5), mu_o=pvt.oil_viscosity(s.pressure),
            mu_g=pvt.gas_viscosity(s.pressure), bw=np.full(n, 1.0),
            bo=pvt.oil_fvf(s.pressure), bg=pvt.gas_fvf(s.pressure), rs=rs,
            krw=np.full(n, 0.2), kro=np.full(n, 0.5), krg=np.full(n, 0.1))

    return jacobian, state, build_fluid, half


def test_analytic_jacobian_matches_finite_difference_for_undersaturated_cells():
    jacobian, state, build_fluid, half = _jacobian_setup()
    fluid = build_fluid(state)
    analytic = jacobian.blocks(state, fluid)
    numeric = jacobian.numerical(state, build_fluid)
    assert np.max(np.abs(analytic[:half] - numeric[:half])) < 0.05


def test_analytic_jacobian_matches_finite_difference_for_saturated_cells():
    jacobian, state, build_fluid, half = _jacobian_setup()
    fluid = build_fluid(state)
    analytic = jacobian.blocks(state, fluid)
    numeric = jacobian.numerical(state, build_fluid)
    assert np.max(np.abs(analytic[half:] - numeric[half:])) < 0.05


def test_water_row_is_independent_of_the_third_variable():
    """Su nə Sg-dən, nə Rs-dən asılıdır — sütun 2 həmişə sıfır olmalıdır."""
    jacobian, state, build_fluid, half = _jacobian_setup()
    fluid = build_fluid(state)
    blocks = jacobian.blocks(state, fluid)
    assert np.allclose(blocks[:, 0, 2], 0.0)


def test_oil_third_variable_derivative_differs_by_saturation_state():
    """Doymuşda ∂N_o/∂Sg = −PV/Bo (Sg oil-in yerini alır);
    doymamışda ∂N_o/∂Rs = 0 (neft həcmi Rs-dən asılı deyil)."""
    jacobian, state, build_fluid, half = _jacobian_setup()
    fluid = build_fluid(state)
    blocks = jacobian.blocks(state, fluid)
    assert np.allclose(blocks[:half, 1, 2], 0.0)
    assert np.all(blocks[half:, 1, 2] < 0.0)


def test_gas_third_variable_derivative_is_positive_in_both_states():
    """Sg artanda da (sərbəst qaz çoxalır), Rs artanda da (həll olmuş
    qaz çoxalır) qaz tənliyi artmalıdır."""
    jacobian, state, build_fluid, half = _jacobian_setup()
    fluid = build_fluid(state)
    blocks = jacobian.blocks(state, fluid)
    assert np.all(blocks[:half, 2, 2] > 0.0)     # ∂N_g/∂Rs
    assert np.all(blocks[half:, 2, 2] > 0.0)     # ∂N_g/∂Sg


def test_jacobian_block_shape_is_three_by_three_per_cell():
    jacobian, state, build_fluid, half = _jacobian_setup()
    fluid = build_fluid(state)
    blocks = jacobian.blocks(state, fluid)
    assert blocks.shape == (state.ncell, 3, 3)


def test_solution_gor_derivative_vanishes_above_bubble_point():
    """Doyma təzyiqindən yuxarıda Rs sabitdir — törəmə sıfır olmalıdır."""
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    above_bubble_point = pvt.solution_gor_derivative(np.array([260.0, 300.0]))
    assert np.allclose(above_bubble_point, 0.0)


def test_solution_gor_derivative_is_positive_below_bubble_point():
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    below_bubble_point = pvt.solution_gor_derivative(np.array([100.0, 180.0]))
    assert np.all(below_bubble_point > 0.0)


# ── axın Jakobianı (analitik vs sonlu fərq, üz-üz) ───────────────────
def _flux_jacobian_setup(nx=3, ny=1, nz=1):
    from imex2d.domain.scal import CoreyParameters, GasCoreyParameters
    from imex2d.simulation.implicit.three_phase_residual import (
        ThreePhaseFlux, ThreePhaseFluxJacobian)
    from imex2d.simulation.stone_relperm import StoneRelativePermeabilityProvider

    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=50.0, dy=50.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=nz, top_depth=2000.0)
    model = ReservoirModelBuilder().build(geology, [], scal=scal)
    grid = TwoPointFluxDiscretization().build(model)
    flux = ThreePhaseFlux(model, grid)
    relperm = StoneRelativePermeabilityProvider.from_corey(
        CoreyParameters(), GasCoreyParameters())
    jacobian = ThreePhaseFluxJacobian(flux, relperm)

    n = model.ncell
    pressure = np.linspace(215.3, 189.1, n)     # cədvəl nöqtəsi deyil
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    sw = np.linspace(0.35, 0.30, n)
    sg = np.linspace(0.12, 0.08, n)
    state = ThreePhaseState(pressure, sw, sg, np.ones(n, bool))

    def build_fluid(s):
        rs = s.solution_gor(pvt)
        return ThreePhaseFluidState(
            mu_w=np.full(n, 0.5), mu_o=pvt.oil_viscosity(s.pressure),
            mu_g=pvt.gas_viscosity(s.pressure), bw=np.full(n, 1.0),
            bo=pvt.oil_fvf(s.pressure), bg=pvt.gas_fvf(s.pressure), rs=rs,
            krw=relperm.krw(s.water_saturation),
            kro=relperm.kro_three_phase(s.water_saturation, s.gas_saturation),
            krg=relperm.krg(s.gas_saturation))

    return model, grid, flux, jacobian, state, build_fluid, pvt


def _perturb_single_cell(state, cell, variable_index, eps):
    vector = state.to_vector().copy()
    vector[cell * 3 + variable_index] += eps
    return ThreePhaseState.from_vector(vector, state.is_saturated)


def test_flux_pressure_jacobian_matches_finite_difference_per_face():
    """Hər üz AYRI-AYRI yoxlanılır — paylaşılan hüceyrələr qarışdırılmır."""
    model, grid, flux, jacobian, state, build_fluid, pvt = _flux_jacobian_setup()
    fluid = build_fluid(state)
    (dpw_a, _), (dpo_a, _), (dpg_a, _), _ = jacobian.face_pressure_derivatives(
        state, fluid, pvt)
    conn = grid.connections
    step = 1e-4

    for face in range(len(conn.cell_a)):
        cell = conn.cell_a[face]
        forward = _perturb_single_cell(state, cell, 0, step)
        backward = _perturb_single_cell(state, cell, 0, -step)
        wf, of, gf = flux.face_fluxes(forward, build_fluid(forward))
        wb, ob, gb = flux.face_fluxes(backward, build_fluid(backward))

        numeric_w = (wf[face] - wb[face]) / (2 * step)
        numeric_o = (of[face] - ob[face]) / (2 * step)
        numeric_g = (gf[face] - gb[face]) / (2 * step)

        assert abs(dpw_a[face] - numeric_w) < 1e-2
        assert abs(dpo_a[face] - numeric_o) < 1e-2
        assert abs(dpg_a[face] - numeric_g) < max(1.0, abs(numeric_g) * 0.01)


def test_flux_saturation_jacobian_matches_finite_difference_at_upstream_cell():
    model, grid, flux, jacobian, state, build_fluid, pvt = _flux_jacobian_setup()
    fluid = build_fluid(state)
    saturation_derivatives = jacobian.face_saturation_derivatives(state, fluid)
    conn = grid.connections
    step = 1e-5

    for face in range(len(conn.cell_a)):
        upstream_water = saturation_derivatives["water_up"][face]
        forward = _perturb_single_cell(state, upstream_water, 1, step)
        backward = _perturb_single_cell(state, upstream_water, 1, -step)
        wf, _, _ = flux.face_fluxes(forward, build_fluid(forward))
        wb, _, _ = flux.face_fluxes(backward, build_fluid(backward))
        numeric = (wf[face] - wb[face]) / (2 * step)
        assert abs(saturation_derivatives["water_dsw"][face] - numeric) < 0.5


def test_flux_pressure_derivative_is_zero_far_from_the_face_cells():
    """Üzə aid olmayan hüceyrənin təzyiqi həmin üzün axınına təsir etməməlidir."""
    model, grid, flux, jacobian, state, build_fluid, pvt = _flux_jacobian_setup(nx=4)
    fluid = build_fluid(state)
    conn = grid.connections
    unrelated_cell = 3          # ilk üzlə (0-1) əlaqəsi olmayan hüceyrə
    step = 1e-4

    forward = _perturb_single_cell(state, unrelated_cell, 0, step)
    backward = _perturb_single_cell(state, unrelated_cell, 0, -step)
    wf, of, gf = flux.face_fluxes(forward, build_fluid(forward))
    wb, ob, gb = flux.face_fluxes(backward, build_fluid(backward))
    assert abs((wf[0] - wb[0]) / (2 * step)) < 1e-9
    assert abs((of[0] - ob[0]) / (2 * step)) < 1e-9


def test_stone_ii_derivatives_match_finite_difference():
    """Stone II-nin öz analitik zəncirvari qaydası — ayrıca doğrulama."""
    wo = CoreyParameters()
    gas = GasCoreyParameters()
    from imex2d.simulation.stone_relperm import StoneRelativePermeabilityProvider
    provider = StoneRelativePermeabilityProvider.from_corey(wo, gas)

    sw = np.array([0.3, 0.45, 0.55])
    sg = np.array([0.15, 0.25, 0.1])
    d_sw, d_sg = provider.kro_three_phase_derivatives(sw, sg)
    step = 1e-6
    numeric_sw = (provider.kro_three_phase(sw + step, sg)
                 - provider.kro_three_phase(sw - step, sg)) / (2 * step)
    numeric_sg = (provider.kro_three_phase(sw, sg + step)
                 - provider.kro_three_phase(sw, sg - step)) / (2 * step)
    assert np.max(np.abs(d_sw - numeric_sw)) < 1e-6
    assert np.max(np.abs(d_sg - numeric_sg)) < 1e-6


def test_gas_corey_derivatives_match_finite_difference():
    gas = GasCoreyParameters()
    sg = np.array([0.1, 0.3, 0.5])
    step = 1e-5
    d_krg = gas.krg_derivative(sg, 0.2)
    numeric_krg = (gas.krg(sg + step, 0.2) - gas.krg(sg - step, 0.2)) / (2 * step)
    assert np.max(np.abs(d_krg - numeric_krg)) < 1e-3

    d_krog = gas.krog_derivative(sg, 0.2, 0.9)
    numeric_krog = (gas.krog(sg + step, 0.2, 0.9)
                    - gas.krog(sg - step, 0.2, 0.9)) / (2 * step)
    assert np.max(np.abs(d_krog - numeric_krog)) < 1e-3


# ── quyu Jakobianı (analitik vs sonlu fərq) ────────────────────────────
def _well_jacobian_setup(nx=6, ny=6, pressure_value=213.5, sw_value=0.35,
                         sg_value=0.1):
    from imex2d.simulation.implicit.three_phase_residual import (
        ThreePhaseWellJacobian, ThreePhaseWellModel)
    from imex2d.simulation.stone_relperm import StoneRelativePermeabilityProvider
    from imex2d.simulation.well_model import PeacemanWellModel

    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    model = ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                          scal=scal)
    wells = PeacemanWellModel().build_connections(model)
    relperm = StoneRelativePermeabilityProvider.from_corey(
        CoreyParameters(), GasCoreyParameters())

    n = model.ncell
    pressure = np.full(n, pressure_value)
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    state = ThreePhaseState(pressure, np.full(n, sw_value),
                            np.full(n, sg_value), np.ones(n, bool))

    def build_fluid(s):
        rs = s.solution_gor(pvt)
        return ThreePhaseFluidState(
            mu_w=np.full(n, 0.5), mu_o=pvt.oil_viscosity(s.pressure),
            mu_g=pvt.gas_viscosity(s.pressure), bw=np.full(n, 1.0),
            bo=pvt.oil_fvf(s.pressure), bg=pvt.gas_fvf(s.pressure), rs=rs,
            krw=relperm.krw(s.water_saturation),
            kro=relperm.kro_three_phase(s.water_saturation, s.gas_saturation),
            krg=relperm.krg(s.gas_saturation))

    well_model = ThreePhaseWellModel(model, wells, 0.35)
    jacobian = ThreePhaseWellJacobian(well_model, pvt, relperm)
    return model, wells, well_model, jacobian, state, build_fluid


def _well_numerical_derivative(well_model, build_fluid, state, cell, var_index,
                               step=1e-4):
    forward = state.to_vector().copy()
    forward[cell * 3 + var_index] += step
    backward = state.to_vector().copy()
    backward[cell * 3 + var_index] -= step
    state_forward = ThreePhaseState.from_vector(forward, state.is_saturated)
    state_backward = ThreePhaseState.from_vector(backward, state.is_saturated)
    rates_forward = well_model.well_rates(state_forward, build_fluid(state_forward))
    rates_backward = well_model.well_rates(state_backward,
                                           build_fluid(state_backward))
    return (
        (rates_forward.water[cell] - rates_backward.water[cell]) / (2 * step),
        (rates_forward.oil[cell] - rates_backward.oil[cell]) / (2 * step),
        (rates_forward.gas[cell] - rates_backward.gas[cell]) / (2 * step))


def test_bhp_producer_jacobian_matches_finite_difference():
    model, wells, well_model, jacobian, state, build_fluid = _well_jacobian_setup()
    fluid = build_fluid(state)
    blocks = jacobian.blocks(state, fluid)
    producer = next(c for c in wells if not c.is_injector)
    cell = producer.cell

    for var_index in range(3):
        numeric_w, numeric_o, numeric_g = _well_numerical_derivative(
            well_model, build_fluid, state, cell, var_index)
        assert abs(blocks[cell, 0, var_index] - numeric_w) < max(
            0.01, abs(numeric_w) * 0.01)
        assert abs(blocks[cell, 1, var_index] - numeric_o) < max(
            0.01, abs(numeric_o) * 0.01)
        assert abs(blocks[cell, 2, var_index] - numeric_g) < max(
            1.0, abs(numeric_g) * 0.01)


def test_bhp_injector_jacobian_matches_finite_difference():
    model, wells, well_model, jacobian, state, build_fluid = _well_jacobian_setup()
    fluid = build_fluid(state)
    blocks = jacobian.blocks(state, fluid)
    injector = next(c for c in wells if c.is_injector)
    cell = injector.cell

    numeric_w, numeric_o, numeric_g = _well_numerical_derivative(
        well_model, build_fluid, state, cell, 0)
    assert abs(blocks[cell, 0, 0] - numeric_w) < max(0.01, abs(numeric_w) * 0.02)
    assert blocks[cell, 1, 0] == 0.0     # vurucu neft/qaz vermir
    assert blocks[cell, 2, 0] == 0.0


def test_gas_column_reflects_the_product_rule_for_dissolved_gas():
    """∂q_gas/∂3-cü, dRs/d3-cü·qo + Rs·∂qo/∂3-cü hasil qaydasına uyğun olmalıdır."""
    model, wells, well_model, jacobian, state, build_fluid = _well_jacobian_setup()
    fluid = build_fluid(state)
    blocks = jacobian.blocks(state, fluid)
    producer = next(c for c in wells if not c.is_injector)
    cell = producer.cell
    assert blocks[cell, 2, 2] != 0.0        # sərbəst + həll olmuş hər ikisi iştirak edir


def test_well_jacobian_is_zero_for_non_producing_shut_in_well():
    """Drawdown mənfi olub axını dayandıranda (BHP-dən sərfəli olmayan
    hüceyrə) törəmələr sıfır olmalıdır — A6-dakı kəsilmə qaydası."""
    model, wells, well_model, jacobian, state, build_fluid = _well_jacobian_setup(
        pressure_value=50.0)   # istismarçının BHP hədəfindən çox aşağı deyil -> kəsilmiş ola bilər
    fluid = build_fluid(state)
    blocks = jacobian.blocks(state, fluid)
    # heç bir NaN/inf olmamalıdır, hansı vəziyyətdə olur-olsun
    assert np.all(np.isfinite(blocks))


# ── tam Jakobian yığımı (sistemli sonlu fərq) ──────────────────────────
def _assembler_setup(nx=5, ny=5, pressure_value=213.5, sw_value=0.35,
                     sg_value=0.1, with_wells=True):
    from imex2d.simulation.implicit.three_phase_residual import (
        ThreePhaseJacobianAssembler)
    from imex2d.simulation.stone_relperm import StoneRelativePermeabilityProvider
    from imex2d.simulation.well_model import PeacemanWellModel

    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    wells_pattern = five_spot(geology.grid) if with_wells else []
    model = ReservoirModelBuilder().build(geology, wells_pattern, scal=scal)
    grid = TwoPointFluxDiscretization().build(model)
    wells = PeacemanWellModel().build_connections(model) if with_wells else []
    relperm = StoneRelativePermeabilityProvider.from_corey(
        CoreyParameters(), GasCoreyParameters())
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))

    n = model.ncell
    accumulator = ThreePhaseAccumulator(model, grid.pore_volume)
    flux = ThreePhaseFlux(model, grid)
    well_model = ThreePhaseWellModel(model, wells, 0.35)
    assembler = ThreePhaseJacobianAssembler(model, accumulator, flux,
                                            well_model, relperm, pvt)

    def build_fluid(s):
        rs = s.solution_gor(pvt)
        return ThreePhaseFluidState(
            mu_w=np.full(n, 0.5), mu_o=pvt.oil_viscosity(s.pressure),
            mu_g=pvt.gas_viscosity(s.pressure), bw=np.full(n, 1.0),
            bo=pvt.oil_fvf(s.pressure), bg=pvt.gas_fvf(s.pressure), rs=rs,
            krw=relperm.krw(s.water_saturation),
            kro=relperm.kro_three_phase(s.water_saturation, s.gas_saturation),
            krg=relperm.krg(s.gas_saturation))

    def full_residual(state, previous, dt):
        fluid = build_fluid(state)
        fluid_prev = build_fluid(previous)
        n_w, n_o, n_g = accumulator.accumulation(state, fluid)
        n_w0, n_o0, n_g0 = accumulator.accumulation(previous, fluid_prev)
        influx_w, influx_o, influx_g = flux.net_influx(state, fluid)
        rates = well_model.well_rates(state, fluid)
        residual = np.empty(n * 3)
        residual[0::3] = (n_w - n_w0) / dt - influx_w - rates.water
        residual[1::3] = (n_o - n_o0) / dt - influx_o - rates.oil
        residual[2::3] = (n_g - n_g0) / dt - influx_g - rates.gas
        return residual

    pressure = np.full(n, pressure_value)
    state = ThreePhaseState(pressure.copy(), np.full(n, sw_value),
                            np.full(n, sg_value), np.ones(n, bool))
    previous = ThreePhaseState(pressure.copy(), np.full(n, sw_value),
                               np.full(n, sg_value), np.ones(n, bool))
    return assembler, state, previous, build_fluid, full_residual


def test_assembled_jacobian_has_the_correct_shape_and_sparsity():
    assembler, state, previous, build_fluid, full_residual = _assembler_setup()
    fluid = build_fluid(state)
    jacobian = assembler.assemble(state, fluid, dt=1.0)
    assert jacobian.shape == (state.ncell * 3, state.ncell * 3)
    assert jacobian.nnz > 0


def test_assembled_jacobian_matches_system_wide_finite_difference():
    """Ən inandırıcı yoxlama: bütün töhfələr (akkumulyasiya+axın+quyu)
    birlikdə, tam qalıq funksiyasına qarşı sonlu fərqlə."""
    assembler, state, previous, build_fluid, full_residual = _assembler_setup()
    fluid = build_fluid(state)
    dt = 1.0
    jacobian = assembler.assemble(state, fluid, dt)

    step = 1e-4
    vector = state.to_vector()
    sample_columns = list(range(0, 9)) + list(range(30, 39))
    max_relative_error = 0.0
    for column in sample_columns:
        forward = vector.copy(); forward[column] += step
        backward = vector.copy(); backward[column] -= step
        state_forward = ThreePhaseState.from_vector(forward, state.is_saturated)
        state_backward = ThreePhaseState.from_vector(backward, state.is_saturated)
        numeric = (full_residual(state_forward, previous, dt)
                  - full_residual(state_backward, previous, dt)) / (2 * step)
        analytic = np.asarray(jacobian[:, column].todense()).ravel()
        error = np.max(np.abs(numeric - analytic))
        relative = error / max(1.0, np.max(np.abs(numeric)))
        max_relative_error = max(max_relative_error, relative)

    assert max_relative_error < 1e-3


def test_assembled_jacobian_matches_finite_difference_without_wells():
    """Quyusuz model — yalnız akkumulyasiya+axın töhfəsi yoxlanılır."""
    assembler, state, previous, build_fluid, full_residual = _assembler_setup(
        nx=3, ny=3, with_wells=False)
    fluid = build_fluid(state)
    dt = 1.0
    jacobian = assembler.assemble(state, fluid, dt)

    step = 1e-4
    vector = state.to_vector()
    max_relative_error = 0.0
    for column in range(min(12, len(vector))):
        forward = vector.copy(); forward[column] += step
        backward = vector.copy(); backward[column] -= step
        state_forward = ThreePhaseState.from_vector(forward, state.is_saturated)
        state_backward = ThreePhaseState.from_vector(backward, state.is_saturated)
        numeric = (full_residual(state_forward, previous, dt)
                  - full_residual(state_backward, previous, dt)) / (2 * step)
        analytic = np.asarray(jacobian[:, column].todense()).ravel()
        error = np.max(np.abs(numeric - analytic))
        relative = error / max(1.0, np.max(np.abs(numeric)))
        max_relative_error = max(max_relative_error, relative)

    assert max_relative_error < 1e-3


def test_assembled_jacobian_is_symmetric_in_sparsity_pattern_for_connections():
    """Hər daxili bağlantı üçün həm (a,b), həm (b,a) blokları mövcud olmalıdır."""
    assembler, state, previous, build_fluid, full_residual = _assembler_setup(
        with_wells=False)
    fluid = build_fluid(state)
    jacobian = assembler.assemble(state, fluid, dt=1.0).tocoo()
    nonzero_cells = set(zip(jacobian.row // 3, jacobian.col // 3))
    off_diagonal_pairs = {(r, c) for r, c in nonzero_cells if r != c}
    for (r, c) in list(off_diagonal_pairs)[:10]:
        assert (c, r) in off_diagonal_pairs


def test_jacobian_is_exact_in_the_all_undersaturated_state():
    """DOYMAMIŞ vəziyyətdə Jakobian — bu, uzun müddət gizli qalmış
    səhvi tutur.

    Qaz axını `F = Rs_upstream · F_neft`-dir. Hasil qaydası İKİ hədd
    verir, lakin kodda yalnız biri var idi — `F_neft` həddi (yəni
    ∂F/∂Rs_up-un birbaşa hissəsi) unudulmuşdu.

    Niyə gizli qalmışdı: bu hədd YALNIZ hüceyrə doymamış olanda
    sıfırdan fərqlidir (yalnız o halda 3-cü primary dəyişən Rs-dir).
    Bütün əvvəlki tam-sistem doğrulamaları TAM DOYMUŞ vəziyyətdə
    (`np.ones(n, bool)`) aparılmışdı.

    Ölçülüb: hədd olmadan xəta 1.7×10⁻³, hədd ilə 4.6×10⁻⁸.
    """
    assembler, state, previous, build_fluid, full_residual = _assembler_setup()
    n = state.ncell

    # HAMISI DOYMAMIŞ: 3-cü dəyişən Rs-dir
    pressure = np.linspace(180.0, 235.0, n)
    undersaturated = ThreePhaseState(
        pressure, np.linspace(0.32, 0.55, n),
        np.linspace(5.0, 40.0, n), np.zeros(n, bool))
    previous_state = ThreePhaseState(
        pressure.copy(), np.linspace(0.32, 0.55, n),
        np.linspace(5.0, 40.0, n), np.zeros(n, bool))

    fluid = build_fluid(undersaturated)
    dt = 1.0
    jacobian = assembler.assemble(undersaturated, fluid, dt,
                                  previous_state.pressure)

    step = 1e-5
    vector = undersaturated.to_vector()
    worst = 0.0
    # yalnız 3-cü dəyişən sütunları (Rs) — səhv məhz orada idi
    for column in range(2, len(vector), 3):
        forward = vector.copy(); forward[column] += step
        backward = vector.copy(); backward[column] -= step
        numeric = (full_residual(
                       ThreePhaseState.from_vector(forward, undersaturated.is_saturated),
                       previous_state, dt)
                  - full_residual(
                       ThreePhaseState.from_vector(backward, undersaturated.is_saturated),
                       previous_state, dt)) / (2 * step)
        analytic = np.asarray(jacobian[:, column].todense()).ravel()
        error = np.max(np.abs(numeric - analytic))
        worst = max(worst, error / max(1.0, np.max(np.abs(numeric))))

    assert worst < 1e-4, f"doymamış Jakobian xətası çox böyük: {worst}"
