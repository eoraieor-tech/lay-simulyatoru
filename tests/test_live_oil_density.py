"""G9 — canlı neftin sıxlığı və cazibə həddinin Jakobianı.

G9a (bu faylın ilk hissəsi): cazibə həddi `ΔΦ = Δp − ½(ρ_a+ρ_b)·g·ΔD`-də
sıxlığın törəmələri əvvəl ATILIRDI ("cazibə üzvündə sıxlığın təzyiqdən
asılılığı nəzərə alınmır"). Ölçüldü (3 laylı model):

    doymuş vəziyyət, təzyiq sütunu     8.2e-5  →  3.2e-7
    qarışıq vəziyyət, 3-cü sütun       7.7e-3  →  1.6e-10   (ρo = ρo_s/Bo(Rs))

Qarışıq vəziyyətdəki təzyiq sütunu (~0.03) DƏYİŞMİR — TB-2-nin mənbəyi
cazibə deyil (tək laylı modeldə də var).
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers import default_scal
from test_deck_oil_branches import _deck_pvt
from test_deck_oil_branches_engine import _between_branches_state, _DeckBranchService
from test_undersaturated_viscosity_engine import _column_kind_errors
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder, five_spot
from imex2d.domain.scal import CoreyParameters, GasCoreyParameters
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.implicit.three_phase_residual import ThreePhaseFluxJacobian
from imex2d.simulation.implicit.three_phase_state import ThreePhaseState
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter

PRESSURE = 330.0


def _engine(nz: int = 3, end_time: float = 10.0):
    deck = _deck_pvt()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=3, ny=3, dx=25.0, dy=25.0, dz=10.0, porosity=0.2, permx_base=150.0,
        nz=nz, top_depth=2000.0, kv_over_kh=0.5)
    model = ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=CoreyParameters(),
        gas_scal=GasCoreyParameters(), pvt_table=deck.to_pvt_table(),
        name="G9 sınağı")
    service = _DeckBranchService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)
    service.deck = deck
    return service.create_engine(model, SimulationConfig(end_time=end_time))


def _saturated_state(engine):
    n = engine.model.ncell
    return ThreePhaseState(np.full(n, PRESSURE), np.full(n, 0.3),
                           np.full(n, 0.08), np.ones(n, bool))


@pytest.fixture
def without_density_derivatives(monkeypatch):
    """G9a-dan ƏVVƏLKİ davranış: sıxlıq törəmələri sıfır."""
    def zero_pressure(self, state, fluid, pvt, bw_p, bo_p, bg_p):
        zeros = np.zeros(state.ncell)
        return zeros, zeros, zeros

    def zero_third(self, state, fluid):
        return np.zeros(state.ncell)

    monkeypatch.setattr(ThreePhaseFluxJacobian, "density_pressure_derivatives",
                        zero_pressure)
    monkeypatch.setattr(ThreePhaseFluxJacobian, "density_third_derivatives",
                        zero_third)


# ═══════════════════════ G9a — cazibə həddinin Jakobianı ══════════════

def test_layered_model_has_gravity_and_flat_model_does_not():
    assert _engine(nz=3).newton.flux.has_gravity
    assert not _engine(nz=1).newton.flux.has_gravity


def test_pressure_column_matches_finite_difference_with_gravity():
    engine = _engine()
    errors = _column_kind_errors(engine, _saturated_state(engine))
    assert errors["p"] < 5e-6, errors
    assert errors["Sw"] < 1e-8 and errors["third"] < 1e-8, errors


def test_third_column_matches_finite_difference_with_gravity():
    engine = _engine()
    errors = _column_kind_errors(engine, _between_branches_state(engine))
    assert errors["third"] < 1e-8, errors
    assert errors["Sw"] < 1e-8, errors


def test_density_derivatives_are_necessary(without_density_derivatives):
    """ZƏRURİLİK: törəmələr sıfırlananda hər iki sütun ölçülmüş səviyyəyə qayıdır."""
    engine = _engine()
    saturated = _column_kind_errors(engine, _saturated_state(engine))
    mixed = _column_kind_errors(engine, _between_branches_state(engine))
    assert saturated["p"] > 5e-5, saturated
    assert mixed["third"] > 1e-3, mixed


def test_mixed_pressure_column_debt_is_not_from_gravity(without_density_derivatives):
    """TB-2: qarışıq vəziyyətin təzyiq sütunu cazibə törəmələrindən asılı deyil."""
    engine = _engine()
    without = _column_kind_errors(engine, _between_branches_state(engine))["p"]
    flat = _engine(nz=1)
    flat_error = _column_kind_errors(flat, _between_branches_state(flat))["p"]
    assert without > 1e-2 and flat_error > 1e-2, (without, flat_error)
