"""PVT cədvəli, korrelyasiyalar və provider testləri (A1)."""

import numpy as np

from helpers import (SKIP_SLOW, default_scal, five_spot_model, make_service,
                     short_config)
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder, five_spot
from imex2d.application.simulation_service import SimulationService
from imex2d.domain.pvt import PVTTable
from imex2d.simulation.pvt import correlations as C
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter


def _table(**kwargs):
    return C.build_pvt_table(**kwargs)


# ── cədvəl yoxlaması ──────────────────────────────────────────────────
def test_generated_table_is_valid():
    assert _table().validate() == []


def test_table_rejects_non_monotonic_pressure():
    table = _table()
    table.pressure[5], table.pressure[6] = table.pressure[6], table.pressure[5]
    assert any("artan" in issue for issue in table.validate())


def test_table_rejects_mismatched_column_length():
    table = _table()
    table.oil_fvf = table.oil_fvf[:-3]
    assert table.validate()


def test_table_rejects_negative_viscosity():
    table = _table()
    table.oil_viscosity[0] = -1.0
    assert table.validate()


# ── korrelyasiyaların fiziki davranışı ────────────────────────────────
def test_solution_gor_increases_then_flattens_at_bubble_point():
    pb = 240.0
    table = _table(bubble_point_bar=pb, pressure_max=400.0)
    below = table.pressure < pb
    above = table.pressure >= pb
    assert np.all(np.diff(table.solution_gor[below]) > 0), "Rs doyma altında artmır"
    assert np.ptp(table.solution_gor[above]) < 1e-9, "Rs doyma üstündə sabit deyil"


def test_oil_fvf_peaks_at_bubble_point():
    pb = 240.0
    table = _table(bubble_point_bar=pb, pressure_max=400.0)
    peak_pressure = table.pressure[int(np.argmax(table.oil_fvf))]
    assert abs(peak_pressure - pb) < (table.pressure[1] - table.pressure[0]) * 2


def test_oil_viscosity_is_minimum_at_bubble_point():
    pb = 240.0
    table = _table(bubble_point_bar=pb, pressure_max=400.0)
    min_pressure = table.pressure[int(np.argmin(table.oil_viscosity))]
    assert abs(min_pressure - pb) < (table.pressure[1] - table.pressure[0]) * 2


def test_heavier_oil_has_higher_viscosity():
    light = C.beggs_robinson_dead_oil_viscosity(api=40.0, temperature_c=70.0)
    heavy = C.beggs_robinson_dead_oil_viscosity(api=18.0, temperature_c=70.0)
    assert heavy > light


def test_bubble_point_correlation_round_trip():
    """Rs(Pb) hesabla, sonra həmin Rs-dən Pb-ni geri al."""
    api, gas_gravity, temperature = 32.0, 0.75, 70.0
    pb_expected = 220.0
    rs = float(C.standing_solution_gor(np.array([pb_expected]), api,
                                       gas_gravity, temperature)[0])
    pb_back = C.standing_bubble_point(rs, api, gas_gravity, temperature)
    assert abs(pb_back - pb_expected) / pb_expected < 0.02


# ── provider ──────────────────────────────────────────────────────────
def test_provider_reproduces_table_values_at_nodes():
    table = _table()
    provider = BlackOilPVTProvider(table)
    assert np.allclose(provider.oil_fvf(table.pressure), table.oil_fvf)
    assert np.allclose(provider.oil_viscosity(table.pressure), table.oil_viscosity)
    assert np.allclose(provider.water_fvf(table.pressure), table.water_fvf)


def test_provider_clamps_outside_table_range():
    """Ekstrapolyasiya yox, sərhəd dəyəri saxlanılır."""
    table = _table(pressure_min=50.0, pressure_max=300.0)
    provider = BlackOilPVTProvider(table)
    assert abs(provider.oil_fvf(1.0) - table.oil_fvf[0]) < 1e-12
    assert abs(provider.oil_fvf(9999.0) - table.oil_fvf[-1]) < 1e-12


def test_provider_rejects_invalid_table():
    table = _table()
    table.pressure = table.pressure[::-1]
    try:
        BlackOilPVTProvider(table)
    except ValueError:
        return
    raise AssertionError("Yararsız cədvəl qəbul edildi")


def test_total_compressibility_is_positive_and_saturation_weighted():
    provider = BlackOilPVTProvider(_table())
    ct_water = provider.total_compressibility(200.0, 1.0)
    ct_oil = provider.total_compressibility(200.0, 0.0)
    assert ct_water > 0 and ct_oil > 0
    assert abs(provider.total_compressibility(200.0, 0.5)
               - 0.5 * (ct_water + ct_oil)) < 1e-12


# ── mühərriklə inteqrasiya ────────────────────────────────────────────
def _pvt_service(scal, table):
    return SimulationService(
        relempl := CoreyRelativePermeabilityAdapter(scal),
        pvt_provider=BlackOilPVTProvider(table))


def test_engine_runs_with_pvt_provider():
    scal = default_scal()
    service = _pvt_service(scal, _table())
    result = service.run(five_spot_model(nx=15, ny=15, scal=scal),
                         short_config(end_time=200.0))
    assert result.converged
    assert result.final_recovery_factor > 0.0


def test_pvt_changes_ooip_versus_static_model():
    """Bo təzyiqdən asılı olduğu üçün OOIP statik haldan fərqlənməlidir."""
    scal = default_scal()
    model = five_spot_model(nx=11, ny=11, scal=scal)
    static_engine = make_service(scal).create_engine(model, short_config(end_time=10.0))
    pvt_engine = _pvt_service(scal, _table()).create_engine(model, short_config(end_time=10.0))
    assert abs(pvt_engine.original_oil_in_place()
               - static_engine.original_oil_in_place()) > 1.0


def test_pvt_table_attached_to_model_is_validated():
    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=7, ny=7, dx=20, dy=20, dz=10, porosity=0.2, permx_base=150.0)
    broken = _table()
    broken.oil_fvf[3] = -1.0
    model = ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=scal, pvt_table=broken)
    assert any("həcm əmsalı" in issue for issue in model.validate())


def test_static_path_unchanged_when_no_pvt_provider():
    """A1-in əsas zəmanəti: PVT verilmədikdə köhnə davranış qorunur."""
    if SKIP_SLOW:
        return
    scal = default_scal()
    result = make_service(scal).run(five_spot_model(nx=21, ny=21, scal=scal),
                                    SimulationConfig(end_time=500.0))
    engine = make_service(scal).create_engine(five_spot_model(nx=21, ny=21, scal=scal),
                                              SimulationConfig(end_time=500.0))
    assert engine.pvt is None
    assert engine._static is not None, "Statik flüid yolu işləmir"
    assert result.converged
