"""Qaz fazasının UI qoşulması — A7, mərhələ 6d (UI hissəsi)."""

import numpy as np

from helpers import default_scal
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.scal import GasCoreyParameters
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.implicit.three_phase_engine import (
    ThreePhaseSimulationEngine)
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter


def _model(with_gas_pvt=True, with_gas_scal=True):
    scal = default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=4, ny=4, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    pvt_table = (build_pvt_table(bubble_point_bar=240.0, include_gas=True)
                if with_gas_pvt else None)
    gas_scal = GasCoreyParameters() if with_gas_scal else None
    return ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=scal, pvt_table=pvt_table,
        gas_scal=gas_scal)


def _service():
    return ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()))


# ── model sahəsi ────────────────────────────────────────────────────
def test_model_carries_gas_scal_parameters():
    model = _model()
    assert model.gas_scal_parameters is not None
    assert isinstance(model.gas_scal_parameters, GasCoreyParameters)


def test_model_without_gas_scal_defaults_to_none():
    model = _model(with_gas_scal=False)
    assert model.gas_scal_parameters is None


# ── mühərrik seçimi ─────────────────────────────────────────────────
def test_service_selects_three_phase_engine_when_gas_pvt_present():
    from imex2d.application.config import SimulationConfig

    model = _model(with_gas_pvt=True, with_gas_scal=True)
    service = _service()
    engine = service.create_engine(model, SimulationConfig(end_time=10.0))
    assert isinstance(engine, ThreePhaseSimulationEngine)


def test_service_selects_two_phase_engine_without_gas_pvt():
    from imex2d.application.config import SimulationConfig

    model = _model(with_gas_pvt=False, with_gas_scal=False)
    service = _service().with_engine(FullyImplicitEngine)
    engine = service.create_engine(model, SimulationConfig(end_time=10.0))
    assert isinstance(engine, FullyImplicitEngine)


def test_three_phase_engine_selection_ignores_impes_choice():
    """Qaz aktivdirsə IMPES seçimi olsa belə üç fazalı mühərrik işlədilir —
    IMPES qazı dəstəkləmir."""
    from imex2d.application.config import SimulationConfig
    from imex2d.simulation.impes_engine import ImpesEngine

    model = _model(with_gas_pvt=True, with_gas_scal=True)
    service = _service().with_engine(ImpesEngine)
    engine = service.create_engine(model, SimulationConfig(end_time=10.0))
    assert isinstance(engine, ThreePhaseSimulationEngine)


# ── mühərrikin özü — kiçik griddə tam iş axını ─────────────────────
def test_three_phase_engine_runs_without_crashing_on_a_small_grid():
    """Ən vacib UI-qoşulma testi: nəticə HƏMİŞƏ SimulationResult
    olmalıdır (converged True ya da False) — heç vaxt istisna."""
    from imex2d.application.config import SimulationConfig

    model = _model()
    model.initial_conditions.datum_pressure = 213.5
    model.initial_conditions.water_saturation = 0.35
    service = _service()
    config = SimulationConfig(end_time=5.0)
    engine = service.create_engine(model, config)
    result = engine.run()
    assert result is not None
    assert hasattr(result, "converged")
    assert hasattr(result, "message")
    assert result.message != ""


def test_three_phase_engine_uses_its_own_compatible_linear_solver():
    """`ScipyCgIluSolver` (A6-nın 2×2 CPR-i üçün) ötürülsə belə,
    üç fazalı mühərrik onu İŞLƏTMİR — öz uyğun həlledicisini saxlayır."""
    from imex2d.application.config import SimulationConfig
    from imex2d.simulation.linear_solver import ScipyCgIluSolver

    model = _model()
    incompatible_solver = ScipyCgIluSolver()
    from imex2d.simulation.stone_relperm import StoneRelativePermeabilityProvider
    from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
    relperm = StoneRelativePermeabilityProvider.from_corey(
        model.scal_parameters, model.gas_scal_parameters)
    pvt = BlackOilPVTProvider(model.pvt_table)
    engine = ThreePhaseSimulationEngine(
        model, SimulationConfig(end_time=5.0), relperm, pvt,
        linear_solver=incompatible_solver)
    assert not isinstance(engine.newton.linear_solver, ScipyCgIluSolver)


# ── nəticə obyektinin qaz sahələri ──────────────────────────────────
def test_partial_result_still_has_series_up_to_failure_point():
    """Yığılma tam bitməsə belə, o nöqtəyə qədər olan seriyalar mövcud
    olmalıdır — istifadəçi bunlara baxa bilməlidir."""
    from imex2d.application.config import SimulationConfig

    model = _model()
    model.initial_conditions.datum_pressure = 213.5
    model.initial_conditions.water_saturation = 0.35
    service = _service()
    config = SimulationConfig(end_time=5.0)
    engine = service.create_engine(model, config)
    result = engine.run()
    # nəticə converged olsun-olmasın, series strukturu mövcuddur
    assert hasattr(result.series, "gas_rate")
    assert hasattr(result.series, "gas_oil_ratio")
    assert hasattr(result.series, "cumulative_gas")
