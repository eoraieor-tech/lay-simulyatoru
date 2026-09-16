"""G3 (2-ci hissə) — çox qollu PVT provider-i MÜHƏRRİKDƏ.

G3a provider-ə `dc_o/dRs` və `dn/dRs` hədlərini əlavə etdi. Mühərrik bu
törəmələri `bo_rs`/`mu_o_rs` kimi flüid vəziyyətindən alır, yəni qalıq və
Jakobian kodu DƏYİŞMƏMƏLİDİR. Bu testlər bunu sonlu fərqlə yoxlayır:
doymamış hüceyrələrin Rs-i iki deck qolunun ARASINDADIR (hədlər sıfırdan
fərqli olan yeganə yer).
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers import default_scal
from test_deck_oil_branches import _deck_pvt, _undersaturated_branches
from test_undersaturated_viscosity_engine import _column_kind_errors
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder, five_spot
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.scal import CoreyParameters, GasCoreyParameters
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.implicit.three_phase_state import ThreePhaseState
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter

#: Hər iki Rs qolunun doyma təzyiqi arasında (276.80 … 345.75 bar)
PRESSURE = 330.0


class _DeckBranchService(ModelAwareSimulationService):
    """Servisin PVT seçimini çox qollu provider ilə əvəz edir."""

    deck = None

    def _build_pvt_provider(self, model):
        return BlackOilPVTProvider(model.pvt_table,
                                   oil_branches=self.deck.oil.branches)


def _engine(end_time: float = 60.0):
    deck = _deck_pvt()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=4, ny=4, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    model = ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=CoreyParameters(),
        gas_scal=GasCoreyParameters(), pvt_table=deck.to_pvt_table(),
        name="G3b sınağı")
    service = _DeckBranchService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)
    service.deck = deck
    return service.create_engine(model, SimulationConfig(end_time=end_time))


def _between_branches_state(engine):
    """Yarısı doymuş (Sg), yarısı doymamış — Rs İKİ QOLUN ARASINDA."""
    low, high = _undersaturated_branches(_deck_pvt())
    n = engine.model.ncell
    saturated = np.zeros(n, bool)
    saturated[: n // 2] = True
    rs = np.linspace(low.solution_gor + 5.0, 250.0, n)
    third = np.where(saturated, 0.08, rs)
    return ThreePhaseState(np.full(n, PRESSURE), np.full(n, 0.3), third, saturated)


def _slopes_zeroed(engine):
    pvt = engine.newton.pvt
    pvt._branch_co_slope = np.zeros_like(pvt._branch_co_slope)
    pvt._branch_n_slope = np.zeros_like(pvt._branch_n_slope)


def test_the_engine_uses_the_branch_provider():
    engine = _engine()
    assert engine.newton.pvt._branch_rs is not None


def test_undersaturated_cells_sit_between_the_branches():
    """Sınaq vəziyyətinin özü: Pb(Rs) < p və Rs qollar arasındadır."""
    engine = _engine()
    state = _between_branches_state(engine)
    under = ~state.is_saturated
    rs = state.third_variable[under]
    pvt = engine.newton.pvt
    assert np.all(pvt.saturation_pressure(rs) < PRESSURE)
    assert np.all((rs > pvt._branch_rs[0]) & (rs < pvt._branch_rs[-1]))


def test_third_column_matches_finite_difference_between_branches():
    engine = _engine()
    errors = _column_kind_errors(engine, _between_branches_state(engine))
    assert errors["Sw"] < 1e-8, errors
    assert errors["third"] < 1e-8, errors


def test_the_branch_slope_terms_are_necessary_in_the_engine():
    engine = _engine()
    state = _between_branches_state(engine)
    with_terms = _column_kind_errors(engine, state)["third"]
    _slopes_zeroed(engine)
    without_terms = _column_kind_errors(engine, state)["third"]
    assert without_terms > with_terms * 1000.0, (with_terms, without_terms)


def test_run_converges_with_the_branch_provider():
    engine = _engine(end_time=120.0)
    result = engine.run()
    assert result.converged, result.message


def test_service_reads_the_branches_from_the_model():
    """Seans 38: qollar `ReservoirModel.pvt_oil_branches` ilə — alt-sinif lazım deyil."""
    deck = _deck_pvt()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=3, ny=3, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    model = ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=CoreyParameters(),
        gas_scal=GasCoreyParameters(), pvt_table=deck.to_pvt_table(),
        pvt_oil_branches=deck.oil.branches, name="qollar modeldə")
    service = ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)
    engine = service.create_engine(model, SimulationConfig(end_time=10.0))
    assert engine.newton.pvt._branch_rs is not None
    assert engine.newton.pvt._branch_co == pytest.approx([2.0564e-4, 1.8317e-4],
                                                         rel=1e-3)
