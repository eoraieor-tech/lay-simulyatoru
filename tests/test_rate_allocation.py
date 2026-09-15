"""RATE hədəfinin perforasiyalara bölünməsi (Seans 27, B7 addım 2-dən əvvəl).

ÖLÇÜLMÜŞ SƏHV: qalıq, Jakobian və IMPES RATE hədəfini HƏR perforasiyaya
TAM yazırdı — 3 perforasiyalı quyu 50 m³/gün hədəfdə 150 m³/gün hasil
edirdi. İndi hədəf `WI·λ` nisbətində bölünür (bax
`simulation/well_constraints.py`); pay addımın əvvəlində verilir və Nyuton
daxilində sabitdir, ona görə Jakobianın strukturu dəyişmir.
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers import default_scal
from test_gas_injection import _jacobian_error as _three_phase_jacobian_error
from test_implicit_jacobian import EXACT, _assemblers, _pvt, _random_state
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.scal import CoreyParameters, GasCoreyParameters
from imex2d.domain.wells import ControlMode, Perforation, Well, WellControl, WellType
from imex2d.simulation.impes_engine import ImpesEngine
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.well_constraints import (assign_rate_shares,
                                                 needs_rate_allocation)
from imex2d.simulation.well_model import PeacemanWellModel, WellConnection

TARGET = 50.0


def _model(producer_cells, injector_cells=((0, 0, 0),), nz=1,
           producer_mode=ControlMode.RATE, injector_mode=ControlMode.BHP,
           injector_target=320.0, include_gas=False, scal=None):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=5, ny=5, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=nz, top_depth=2000.0)
    wells = [
        Well("INJ", WellType.INJECTOR, WellControl(injector_mode, injector_target),
             [Perforation(*cell) for cell in injector_cells]),
        Well("PROD", WellType.PRODUCER, WellControl(producer_mode, TARGET),
             [Perforation(*cell) for cell in producer_cells]),
    ]
    return ReservoirModelBuilder().build(
        geology, wells, scal=scal or default_scal(),
        gas_scal=GasCoreyParameters() if include_gas else None,
        pvt_table=build_pvt_table(pressure_min=1.0, pressure_max=400.0,
                                  n_points=40, bubble_point_bar=240.0,
                                  include_gas=include_gas),
        name="RATE payı sınağı")


def _layers(i, j, nz):
    return [(i, j, k) for k in range(nz)]


def _service(engine=FullyImplicitEngine):
    return ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=engine)


def _connection(cell, well_index, name="P"):
    return WellConnection(name, cell, well_index, False, ControlMode.RATE, TARGET)


# ═══════════════════════ pay qaydası ═════════════════════════════════

def test_single_perforation_share_is_exactly_one():
    """GERİYƏ UYĞUNLUQ: `target * 1.0` — tək perforasiyada bit-bit eyni."""
    connections = PeacemanWellModel().build_connections(_model([(4, 4, 0)]))
    assert all(c.rate_share == 1.0 for c in connections)


def test_initial_shares_follow_the_well_index_and_sum_to_one():
    connections = [_connection(0, 2.0), _connection(1, 6.0)]
    assign_rate_shares(connections)
    assert [c.rate_share for c in connections] == pytest.approx([0.25, 0.75])


def test_shares_follow_well_index_times_mobility():
    connections = [_connection(0, 2.0), _connection(1, 6.0), _connection(2, 1.0, "Q")]
    assign_rate_shares(connections, [3.0, 1.0, 0.7])
    assert connections[0].rate_share == pytest.approx(0.5)
    assert connections[1].rate_share == pytest.approx(0.5)
    assert connections[2].rate_share == 1.0, "başqa quyunun payı ayrıca hesablanır"


def test_zero_mobility_everywhere_falls_back_to_the_well_index():
    connections = [_connection(0, 1.0), _connection(1, 3.0)]
    assign_rate_shares(connections, [0.0, 0.0])
    assert [c.rate_share for c in connections] == pytest.approx([0.25, 0.75])
    assert sum(c.rate_share for c in connections) == pytest.approx(1.0)


def test_allocation_is_needed_only_for_multi_perforation_rate_wells():
    single = PeacemanWellModel().build_connections(_model([(4, 4, 0)]))
    layered_bhp = PeacemanWellModel().build_connections(
        _model(_layers(4, 4, 3), nz=3, producer_mode=ControlMode.BHP,
               injector_cells=[(0, 0, 0)]))
    layered_rate = PeacemanWellModel().build_connections(
        _model(_layers(4, 4, 3), nz=3))
    assert not needs_rate_allocation(single)
    assert not needs_rate_allocation(layered_bhp)
    assert needs_rate_allocation(layered_rate)


# ═══════════════════════ quyu hədəfi BİR DƏFƏ çatdırılır ══════════════

def test_two_phase_rate_wells_deliver_the_target_once():
    """ƏSAS DÜZƏLİŞ: 3 perforasiya → hədəfin 3 qatı DEYİL, özü."""
    model = _model(_layers(4, 4, 3), _layers(0, 0, 3), nz=3,
                   injector_mode=ControlMode.RATE, injector_target=TARGET)
    # Yalnız RATE quyulu model servisin yoxlamasından keçmir (təzyiq
    # səviyyəsi qeyri-müəyyəndir) — debit hesabını yoxlamaq üçün mühərrik
    # birbaşa qurulur, addım atılmır.
    engine = FullyImplicitEngine(model, SimulationConfig(end_time=5.0),
                                 CoreyRelativePermeabilityAdapter(default_scal()),
                                 pvt=_pvt())
    engine._update_rate_shares()
    assembler = engine.residual_assembler
    fluid = assembler.fluid_state(engine.state)
    rates = assembler.well_rates(engine.state, fluid)

    produced = -sum(rates.water[c.cell] * fluid.bw[c.cell]
                    + rates.oil[c.cell] * fluid.bo[c.cell]
                    for c in assembler.wells if c.well_name == "PROD")
    injected = sum(rates.water[c.cell] * fluid.bw[c.cell]
                   for c in assembler.wells if c.well_name == "INJ")
    assert produced == pytest.approx(TARGET, rel=1e-12)
    assert injected == pytest.approx(TARGET, rel=1e-12)


def test_three_phase_rate_well_delivers_the_target_once():
    model = _model(_layers(4, 4, 3), nz=3, include_gas=True)
    engine = _service().create_engine(model, SimulationConfig(end_time=5.0))
    engine._update_rate_shares()
    newton = engine.newton
    fluid = newton.build_fluid(engine.state)
    rates = newton.well_model.well_rates(engine.state, fluid)
    produced = -sum(rates.water[c.cell] * fluid.bw[c.cell]
                    + rates.oil[c.cell] * fluid.bo[c.cell]
                    for c in newton.well_model.wells if c.well_name == "PROD")
    assert produced == pytest.approx(TARGET, rel=1e-12)


def test_engine_shares_follow_the_mobility_of_the_converged_state():
    model = _model(_layers(4, 4, 3), nz=3)
    engine = _service().create_engine(model, SimulationConfig(end_time=5.0))
    engine.state.water_saturation[model.grid.index(4, 4, 0)] = 0.6   # sulanmış təbəqə
    engine._update_rate_shares()
    assembler = engine.residual_assembler
    mobility = assembler.connection_mobilities(assembler.fluid_state(engine.state))
    producer = [(c, m) for c, m in zip(assembler.wells, mobility)
                if c.well_name == "PROD"]
    weights = np.array([c.well_index * m for c, m in producer])
    assert [c.rate_share for c, _ in producer] == pytest.approx(
        list(weights / weights.sum()), rel=1e-12)
    assert len({round(c.rate_share, 9) for c, _ in producer}) > 1


@pytest.mark.parametrize("engine, include_gas", [
    (FullyImplicitEngine, False), (ImpesEngine, False), (FullyImplicitEngine, True)],
    ids=["tam-implicit", "impes", "uc-fazali"])
def test_layered_rate_well_produces_like_a_single_perforation(engine, include_gas):
    """UC-UCA: təbəqə sayı quyunun CƏMİ debitini dəyişməməlidir.

    Düzəlişdən əvvəl 3 perforasiyalı quyu ~3 dəfə çox hasil edirdi.
    """
    def cumulative_liquid(nz):
        model = _model(_layers(4, 4, nz), _layers(0, 0, nz), nz=nz,
                       include_gas=include_gas)
        result = _service(engine).run(model, SimulationConfig(end_time=10.0))
        assert result.converged, result.message
        return result.series.cumulative_oil[-1] + result.series.cumulative_water[-1]

    ratio = cumulative_liquid(3) / cumulative_liquid(1)
    assert ratio == pytest.approx(1.0, abs=0.05), ratio


# ═══════════════════════ Jakobian — sonlu fərq ═══════════════════════

def test_two_phase_jacobian_with_split_rate_targets_is_exact():
    """Pay Nyuton daxilində SABİTDİR — törəmələr sadəcə paya vurulur."""
    scal = default_scal()
    model = _model([(4, 4, 0), (3, 4, 0)], [(0, 0, 0), (1, 0, 0)],
                   injector_mode=ControlMode.RATE, injector_target=40.0, scal=scal)
    residual, jacobian = _assemblers(model, scal, pvt=_pvt())
    assign_rate_shares(residual.wells, [1.0, 3.0, 2.0, 5.0])
    assert {round(c.rate_share, 6) for c in residual.wells} == {0.25, 0.75,
                                                               round(2 / 7, 6),
                                                               round(5 / 7, 6)}
    state, previous = _random_state(model, scal)
    _, fluid, _ = residual.residual(state, previous, 5.0)
    analytic = jacobian.assemble(state, fluid, 5.0).toarray()
    numeric = jacobian.numerical(state, previous, 5.0)
    scale = np.maximum(np.abs(numeric).max(axis=0), 1e-30)
    assert float(np.max(np.abs(analytic - numeric) / scale)) < EXACT


def test_three_phase_jacobian_with_split_rate_target_is_not_worse():
    """Üç fazalı RATE istismarçısının Jakobianı ƏVVƏLDƏN təqribidir (sərbəst
    qazın törəməsi buraxılır — köhnə sadələşdirmə). Pay onu PİSLƏŞDİRMƏMƏLİDİR:
    iki perforasiyalı quyunun xətası tək perforasiyalınınkından böyük deyil."""
    single = _three_phase_jacobian_error(
        _model([(4, 4, 0)], include_gas=True, scal=CoreyParameters()),
        columns=range(75))
    split = _three_phase_jacobian_error(
        _model([(4, 4, 0), (3, 4, 0)], include_gas=True, scal=CoreyParameters()),
        columns=range(75))
    assert split <= single * (1.0 + 1e-9), (split, single)
