"""G6 — süxur sıxılmasının İSTİNAD TƏZYİQİ (SPE1: `ROCK 14.7 3E-6`).

    PV(p) = PV_ref · [1 + c_r · (p − p_istinad)]

Əvvəl `p_istinad` HƏMİŞƏ `datum_pressure` idi. Eclipse `ROCK` açar sözü isə
istinadı ayrıca verir — SPE1-də 14.7 psia, datum isə 4800 psia. Fərq məsamə
həcmində ~1.4 %-dir.

ƏN VACİB TEST (sahibkarın xüsusi tələbi): kodda İKİ AYRI istinad təzyiqi var
və G6 yalnız BİRİNCİSİNƏ toxunur:

  * məsamə həcmi  — `PV(p)`      → G6 bunu ayırır;
  * FLÜİDİN statik sıxılma modeli — `B(p) = B_ref/(1+c·Δp)`, PVT provider
    olmayanda — **TOXUNULMUR**, istinadı datum olaraq qalır.

`test_fluid_static_model_keeps_the_datum_reference` bunu kilidləyir.
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder, five_spot
from imex2d.application.serialization import ProjectSerializer
from imex2d.domain.properties import RockProperties
from imex2d.domain.scal import CoreyParameters, GasCoreyParameters
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.implicit.residual import (ResidualAssembler,
                                                 rock_reference_pressure)
from imex2d.simulation.implicit.state import ReservoirState
from imex2d.simulation.implicit.three_phase_residual import (
    ThreePhaseAccumulationJacobian, ThreePhaseAccumulator, ThreePhaseFluidState)
from imex2d.simulation.implicit.three_phase_state import ThreePhaseState
from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.well_model import PeacemanWellModel

#: SPE1 `ROCK 14.7 3E-6` — FIELD vahidlərindən metrikə
SPE1_REFERENCE_BAR = 14.7 / 14.5037738
SPE1_COMPRESSIBILITY_PER_BAR = 3.0e-6 * 14.5037738

DATUM = 250.0
COMPRESSIBILITY = 4.5e-5


def _model(reference=None, datum=DATUM, compressibility=COMPRESSIBILITY,
           reference_unit="bar", include_gas=False):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=4, ny=4, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=1, top_depth=2000.0)
    from imex2d.domain.initial import InitialConditions
    return ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=default_scal(),
        gas_scal=GasCoreyParameters() if include_gas else None,
        initial=InitialConditions(datum_pressure=datum),
        pvt_table=build_pvt_table(pressure_min=1.0, pressure_max=400.0,
                                  n_points=40, bubble_point_bar=240.0,
                                  include_gas=include_gas),
        rock_compressibility=compressibility,
        rock_compressibility_reference=reference,
        rock_compressibility_reference_unit=reference_unit,
        name="G6 sınağı")


def _assembler(model):
    grid = TwoPointFluxDiscretization().build(model)
    wells = PeacemanWellModel().build_connections(model)
    return ResidualAssembler(model, grid, wells,
                             CoreyRelativePermeabilityAdapter(default_scal()))


def _accumulator(model):
    grid = TwoPointFluxDiscretization().build(model)
    return ThreePhaseAccumulator(model, grid.pore_volume)


# ═══════════════════════ domen və ötürücülər ═════════════════════════

def test_default_is_none_so_nothing_changes():
    assert RockProperties.__dataclass_fields__[
        "compressibility_reference_pressure"].default is None
    assert _model().rock.compressibility_reference_pressure is None


def test_helper_falls_back_to_the_datum_pressure():
    assert rock_reference_pressure(_model()) == pytest.approx(DATUM)


def test_helper_uses_the_explicit_reference():
    assert rock_reference_pressure(_model(reference=1.0)) == pytest.approx(1.0)


def test_builder_converts_the_unit():
    """SPE1 dəyəri psi-dədir — mühərrik bar işlədir."""
    model = _model(reference=14.7, reference_unit="psi")
    assert (model.rock.compressibility_reference_pressure
            == pytest.approx(SPE1_REFERENCE_BAR, rel=1e-6))


def test_reference_survives_the_project_round_trip():
    serializer = ProjectSerializer()
    for reference in (None, 1.0135):
        model = _model(reference=reference)
        restored = serializer.reservoir_model_from_dict(
            serializer.reservoir_model_to_dict(model))
        if reference is None:
            assert restored.rock.compressibility_reference_pressure is None
        else:
            assert (restored.rock.compressibility_reference_pressure
                    == pytest.approx(reference))


# ═══════════════════════ məsamə həcmi — iki fazalı ═══════════════════

def test_two_phase_pore_volume_is_unchanged_without_a_reference():
    """GERİYƏ UYĞUNLUQ: `None` → köhnə düstur (datum), BİT-BİT."""
    assembler = _assembler(_model())
    pressure = np.linspace(180.0, 320.0, assembler.ncell)
    expected = assembler.pore_volume * (
        1.0 + COMPRESSIBILITY * (pressure - DATUM))
    assert np.array_equal(assembler.pore_volume_at(pressure), expected)


def test_two_phase_pore_volume_uses_the_explicit_reference():
    reference = 1.0135
    assembler = _assembler(_model(reference=reference))
    pressure = np.linspace(180.0, 320.0, assembler.ncell)
    expected = assembler.pore_volume * (
        1.0 + COMPRESSIBILITY * (pressure - reference))
    assert np.allclose(assembler.pore_volume_at(pressure), expected, rtol=0, atol=0)


def test_reference_equal_to_datum_reproduces_the_old_result():
    """Açıq şəkildə datum verilsə, nəticə `None` halı ilə EYNİ olmalıdır."""
    pressure = np.linspace(180.0, 320.0, 16)
    without = _assembler(_model()).pore_volume_at(pressure)
    explicit = _assembler(_model(reference=DATUM)).pore_volume_at(pressure)
    assert np.array_equal(without, explicit)


# ═══════════════════════ məsamə həcmi — üç fazalı ════════════════════

def test_three_phase_pore_volume_uses_the_explicit_reference():
    reference = 1.0135
    accumulator = _accumulator(_model(reference=reference, include_gas=True))
    pressure = np.linspace(180.0, 320.0, accumulator.model.ncell)
    expected = accumulator.pore_volume * (
        1.0 + COMPRESSIBILITY * (pressure - reference))
    assert np.allclose(accumulator.pore_volume_at(pressure), expected, rtol=0, atol=0)


def test_three_phase_pore_volume_is_unchanged_without_a_reference():
    accumulator = _accumulator(_model(include_gas=True))
    pressure = np.linspace(180.0, 320.0, accumulator.model.ncell)
    expected = accumulator.pore_volume * (
        1.0 + COMPRESSIBILITY * (pressure - DATUM))
    assert np.array_equal(accumulator.pore_volume_at(pressure), expected)


# ═══════════════════════ İKİ İSTİNADIN AYRILMASI ═════════════════════

def test_fluid_static_model_keeps_the_datum_reference():
    """ƏSAS TƏLƏB: G6 YALNIZ məsamə həcminə toxunur.

    PVT provider olmayanda flüidin `B(p)` modeli də istinad təzyiqi işlədir
    (`B_ref/(1+c·Δp)`). O, DATUM olaraq qalmalıdır — əks halda süxur üçün
    verilən istinad səssizcə FLÜİDİN sıxılmasını da dəyişərdi.
    """
    plain = _assembler(_model())
    shifted = _assembler(_model(reference=1.0135))

    assert plain.reference_pressure == pytest.approx(DATUM)
    assert shifted.reference_pressure == pytest.approx(DATUM), \
        "flüidin istinadı süxurunkundan ASILI OLMAMALIDIR"
    assert shifted.rock_reference_pressure == pytest.approx(1.0135)

    state = ReservoirState(np.full(plain.ncell, 300.0),
                           np.full(plain.ncell, 0.35))
    first, second = plain.fluid_state(state), shifted.fluid_state(state)
    assert np.array_equal(first.bw, second.bw)
    assert np.array_equal(first.bo, second.bo)


# ═══════════════════════ Jakobian dəyişmir ═══════════════════════════

def test_accumulation_jacobian_matches_finite_difference_with_a_reference():
    """`d(PV)/dp = PV_ref·c_r` — istinaddan ASILI DEYİL, ona görə Jakobian
    kodu toxunulmadı. Bu test həmin iddianı kilidləyir."""
    model = _model(reference=1.0135, include_gas=True)
    accumulator = _accumulator(model)
    pvt = BlackOilPVTProvider(build_pvt_table(bubble_point_bar=240.0,
                                              include_gas=True))
    n = model.ncell

    def build_fluid(state):
        return ThreePhaseFluidState(
            mu_w=np.full(n, 0.5), mu_o=pvt.oil_viscosity(state.pressure),
            mu_g=pvt.gas_viscosity(state.pressure),
            # Bw PVT-DƏN götürülməlidir: analitik blok `water_fvf_derivative`
            # işlədir, sabit Bw ilə sonlu fərq uyğunsuz olardı (test qüsuru).
            bw=pvt.water_fvf(state.pressure),
            bo=pvt.oil_fvf(state.pressure), bg=pvt.gas_fvf(state.pressure),
            rs=state.solution_gor(pvt), krw=np.full(n, 0.2),
            kro=np.full(n, 0.5), krg=np.full(n, 0.1))

    state = ThreePhaseState(np.full(n, 230.0), np.full(n, 0.3),
                            np.full(n, 0.1), np.ones(n, bool))
    blocks = ThreePhaseAccumulationJacobian(accumulator, pvt).blocks(
        state, build_fluid(state))

    step = 1e-4
    vector = state.to_vector()
    for cell in range(min(n, 4)):
        forward, backward = vector.copy(), vector.copy()
        forward[cell * 3] += step
        backward[cell * 3] -= step
        state_f = ThreePhaseState.from_vector(forward, state.is_saturated)
        state_b = ThreePhaseState.from_vector(backward, state.is_saturated)
        numeric = [(a[cell] - b[cell]) / (2 * step) for a, b in zip(
            accumulator.accumulation(state_f, build_fluid(state_f)),
            accumulator.accumulation(state_b, build_fluid(state_b)))]
        for equation in range(3):
            assert blocks[cell, equation, 0] == pytest.approx(
                numeric[equation], rel=1e-4), (cell, equation)


# ═══════════════════════ SPE1 böyüklüyü ══════════════════════════════

def test_spe1_reference_changes_the_pore_volume_by_about_one_and_a_half_percent():
    """SPE1-in öz rəqəmləri: `ROCK 14.7 3E-6`, datum 4800 psia."""
    datum_bar = 4800.0 / 14.5037738
    with_reference = _model(reference=SPE1_REFERENCE_BAR, datum=datum_bar,
                            compressibility=SPE1_COMPRESSIBILITY_PER_BAR)
    without = _model(datum=datum_bar, compressibility=SPE1_COMPRESSIBILITY_PER_BAR)

    pressure = np.full(with_reference.ncell, datum_bar)
    ratio = (_assembler(with_reference).pore_volume_at(pressure)
             / _assembler(without).pore_volume_at(pressure))
    assert np.all(ratio > 1.013) and np.all(ratio < 1.016), ratio[0]


def test_rock_fluid_panel_exposes_the_reference_field():
    pytest.importorskip("PyQt5.QtWidgets")
    import inspect
    from imex2d.ui import panels
    source = inspect.getsource(panels.RockFluidPanel)
    assert "rock_reference_enabled" in source
    assert "rock_compressibility_reference_value" in source
