# -*- coding: utf-8 -*-
"""ACTNUM — qeyri-aktiv hüceyrələrin simulyasiyadan ÇIXARILMASI.

Bu fayl dörd audit məsələsini bağlayır:

    #11  ACTNUM oxunurdu, amma qeyri-aktiv hüceyrə simulyasiyadan
         çıxarılmırdı
    #12  qeyri-aktiv hüceyrə PV/OOIP/material balansına daxil idi
    #13  aktiv↔qeyri-aktiv transmissivlik bağlantısı qurulurdu
    #19  OPM idxalı `active_count`-u yalnız xəbərdarlıq üçün işlədirdi

Ən güclü test `test_padding_cells_marked_inactive_reproduce_the_small_model`-dir:
qeyri-aktiv "doldurucu" hüceyrələr əlavə etmək simulyasiyanın NƏTİCƏSİNİ
DƏYİŞMƏMƏLİDİR — nə OOIP, nə hasilat, nə də təzyiq. Fizikanın qeyri-aktiv
hüceyrələrdən TAM təcrid olunduğunu yalnız belə bir bərabərlik sübut edir.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from helpers import default_scal, make_service, short_config
from imex2d.application.config import (OutputConfig, SimulationConfig,
                                       TimeSteppingConfig)
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder
from imex2d.domain.diagnostics import Severity
from imex2d.domain.grid import ActiveMap, CartesianGrid, Connections
from imex2d.domain.initial import InitialConditions
from imex2d.domain.wells import (ControlMode, Perforation, Well, WellControl,
                                 WellType)
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.impes_engine import ImpesEngine
from imex2d.simulation.well_model import PeacemanWellModel


# ═══════════════════════════════════════════════════════ model qurucular
def _model(nx, inactive=(), scal=None, ny=1, dx=10.0, permeability=200.0):
    """1D model; `inactive` — qeyri-aktiv qlobal hüceyrə indeksləri.

    Quyular HƏMİŞƏ 0 və `nx_active_end` hüceyrələrindədir — aşağıdakı
    bərabərlik testində iki modelin quyuları eyni fiziki yerdə olsun.
    """
    scal = scal or default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=dx, dy=100.0, dz=10.0,
        porosity=0.20, permx_base=permeability)

    if inactive:
        actnum = np.ones(geology.grid.ncell, dtype=np.int8)
        actnum[list(inactive)] = 0
        grid = geology.grid.with_actnum(actnum)
        geology.grid = grid
        geology.geometry = dataclasses.replace(geology.geometry, grid=grid)

    return geology


def _built(geology, producer_i, scal=None, injection_rate=40.0):
    scal = scal or default_scal()
    wells = [
        Well("INJ", WellType.INJECTOR,
             WellControl(ControlMode.RATE, injection_rate),
             [Perforation(0, 0, 0)]),
        Well("PROD", WellType.PRODUCER, WellControl(ControlMode.BHP, 180.0),
             [Perforation(producer_i, 0, 0)]),
    ]
    return ReservoirModelBuilder().build(
        geological_model=geology, wells=wells, scal=scal,
        initial=InitialConditions(datum_pressure=220.0,
                                  water_saturation=scal.swc),
        name="ACTNUM testi")


# ═══════════════════════════════════════════════════════ ActiveMap (§1)
def test_active_map_is_a_two_way_bijection():
    actnum = np.array([1, 0, 1, 1, 0, 1])
    active = ActiveMap.from_actnum(actnum, 6)

    assert active.n_global == 6
    assert active.n_active == 4
    assert active.n_inactive == 2
    assert list(active.global_to_active) == [0, -1, 1, 2, -1, 3]
    assert list(active.active_to_global) == [0, 2, 3, 5]

    # tərs çevirmələr bir-birini bərpa edir
    for a in range(active.n_active):
        assert active.global_to_active[active.active_to_global[a]] == a
    for g in active.active_to_global:
        assert active.active_to_global[active.global_to_active[g]] == g


def test_inactive_cells_are_marked_minus_one_not_zero():
    """`global_to_active` -də qeyri-aktiv hüceyrə −1-dir; 0 OLSAYDI o,
    birinci aktiv hüceyrə ilə səssizcə qarışardı."""
    active = ActiveMap.from_actnum([0, 1, 1], 3)
    assert active.global_to_active[0] == -1
    assert active.global_to_active[1] == 0


def test_all_active_map_is_the_identity():
    active = ActiveMap.all_active(5)
    assert active.all_cells_active
    assert list(active.global_to_active) == [0, 1, 2, 3, 4]
    assert list(active.active_to_global) == [0, 1, 2, 3, 4]


def test_nan_actnum_counts_as_active():
    """Deck-də `n*` defoltu ilə yazılmış hüceyrəni səssizcə YOX ETMƏK
    təhlükəlidir — NaN aktiv sayılır."""
    active = ActiveMap.from_actnum([np.nan, 1.0, 0.0], 3)
    assert active.n_active == 2
    assert active.actnum[0] == 1


def test_active_map_rejects_a_size_mismatch():
    with pytest.raises(ValueError):
        ActiveMap.from_actnum([1, 1, 1], 5)


def test_grid_without_actnum_behaves_exactly_as_before():
    grid = CartesianGrid(4, 3, 2)
    assert grid.n_active == grid.ncell == 24
    assert not grid.has_inactive_cells
    assert grid.active.all_cells_active
    # bərabərlik/hash əvvəlki kimi yalnız ölçülərə görədir
    assert grid == CartesianGrid(4, 3, 2)
    assert len({grid, CartesianGrid(4, 3, 2)}) == 1


def test_grid_equality_ignores_actnum_but_keeps_the_mapping():
    grid = CartesianGrid(2, 2, 1, [1, 1, 0, 1])
    assert grid.n_active == 3
    assert grid.has_inactive_cells
    assert grid == CartesianGrid(2, 2, 1)          # __eq__ yalnız ölçülər
    assert hash(grid) == hash(CartesianGrid(2, 2, 1))


# ═════════════════════════════════════════════ bağlantılar / TPFA (§3)
def test_no_connection_touches_an_inactive_cell():
    grid = CartesianGrid(3, 3, 1, [1, 1, 1, 1, 0, 1, 1, 1, 1])
    conn = grid.build_connections()

    assert 4 not in conn.cell_a
    assert 4 not in conn.cell_b
    # 3x3-də 12 üz var; mərkəzi hüceyrə 4 üz aparır
    assert conn.count == 8
    assert np.all(grid.active.active_face_mask(conn.cell_a, conn.cell_b))


def test_inactive_cells_never_enter_the_sparse_matrix_topology():
    """T_AB = 0 yazmaq KİFAYƏT DEYİL — üz ÜMUMİYYƏTLƏ olmamalıdır,
    yoxsa qeyri-aktiv hüceyrə matrisin strukturunda naməlum kimi qalır."""
    scal = default_scal()
    model = _built(_model(6, inactive=(3,), scal=scal), producer_i=5, scal=scal)
    grid = TwoPointFluxDiscretization().build(model)

    assert grid.connections.count == 3          # 5 üzdən 2-si atılıb
    assert 3 not in grid.connections.cell_a
    assert 3 not in grid.connections.cell_b
    assert grid.transmissibility.size == grid.connections.count


def test_discretization_rejects_hand_built_active_inactive_connections():
    """Müqavilə qoruyucusu: kimsə `Connections`-ı əl ilə qurub ötürsə,
    səssiz yanlış nəticə əvəzinə AÇIQ xəta alsın."""
    scal = default_scal()
    model = _built(_model(4, inactive=(2,), scal=scal), producer_i=3, scal=scal)
    model._connections = Connections(np.array([1]), np.array([2]),
                                     np.array([0], dtype=np.int8))
    with pytest.raises(ValueError, match="qeyri-aktiv"):
        TwoPointFluxDiscretization().build(model)


def test_a_barrier_of_inactive_cells_splits_the_grid_in_two():
    """Qeyri-aktiv sütun MANEƏDİR — axın onu KEÇƏ BİLMƏZ."""
    nx, ny = 5, 3
    barrier = [i for i in range(nx * ny) if i % nx == 2]
    actnum = np.ones(nx * ny, dtype=np.int8)
    actnum[barrier] = 0
    conn = CartesianGrid(nx, ny, 1, actnum).build_connections()

    left = {c for c in range(nx * ny) if c % nx < 2}
    right = {c for c in range(nx * ny) if c % nx > 2}
    for a, b in zip(conn.cell_a, conn.cell_b):
        assert not ({int(a), int(b)} & left and {int(a), int(b)} & right), \
            "maneədən axın sızır"


# ═════════════════════════════════════════════ məsamə həcmi / OOIP (§2)
def test_pore_volume_is_exactly_zero_in_inactive_cells():
    scal = default_scal()
    model = _built(_model(5, inactive=(1, 3), scal=scal), producer_i=4,
                   scal=scal)
    pv = model.pore_volume()

    assert pv[1] == 0.0
    assert pv[3] == 0.0
    assert np.all(pv[[0, 2, 4]] > 0.0)
    assert model.bulk_volume()[1] == 0.0


def test_ooip_ignores_inactive_cells():
    """OOIP qeyri-aktiv hüceyrələrin payını SAYMAMALIDIR (#12)."""
    scal = default_scal()
    full = _built(_model(8, scal=scal), producer_i=7, scal=scal)
    holed = _built(_model(8, inactive=(4,), scal=scal), producer_i=7, scal=scal)

    engine_full = ImpesEngine(full, short_config(),
                              *_providers(scal))
    engine_holed = ImpesEngine(holed, short_config(), *_providers(scal))

    ooip_full = engine_full.original_oil_in_place()
    ooip_holed = engine_holed.original_oil_in_place()

    assert ooip_holed < ooip_full
    assert ooip_holed == pytest.approx(ooip_full * 7.0 / 8.0, rel=1e-9)


def _providers(scal):
    from imex2d.simulation.linear_solver import ScipyCgIluSolver
    from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
    return CoreyRelativePermeabilityAdapter(scal), ScipyCgIluSolver()


# ═════════════════════════════════════════════════════ xətti sistem (§1)
def test_impes_pressure_matrix_has_exactly_n_active_unknowns():
    scal = default_scal()
    model = _built(_model(10, inactive=(4, 5), scal=scal), producer_i=9,
                   scal=scal)
    engine = ImpesEngine(model, short_config(), *_providers(scal))

    assert model.ncell == 10
    assert model.n_active == 8
    assert engine._matrix.shape == (8, 8), \
        "matris QLOBAL ölçüdə qalıb — qeyri-aktiv sətirlər sinqulyardır"


def test_newton_system_is_reduced_to_active_degrees_of_freedom():
    from imex2d.simulation.implicit.active_reduction import ActiveDofReduction

    active = ActiveMap.from_actnum([1, 0, 1], 3)
    reduction = ActiveDofReduction(active)
    assert reduction.size == 6                  # 3 hüceyrə × 2 dəyişən
    assert reduction.reduced_size == 4           # 2 aktiv × 2 dəyişən
    assert list(reduction.kept_dofs) == [0, 1, 4, 5]

    delta = reduction.expand(np.array([1.0, 2.0, 3.0, 4.0]))
    assert list(delta) == [1.0, 2.0, 0.0, 0.0, 3.0, 4.0]


def test_reduction_keeps_trailing_well_unknowns():
    from imex2d.simulation.implicit.active_reduction import ActiveDofReduction

    active = ActiveMap.from_actnum([1, 0], 2)
    reduction = ActiveDofReduction(active, extra_dofs=2)
    assert reduction.size == 6                   # 2×2 hüceyrə + 2 quyu
    assert list(reduction.kept_dofs) == [0, 1, 4, 5]


def test_reduction_is_a_no_op_when_every_cell_is_active():
    from imex2d.simulation.implicit.active_reduction import ActiveDofReduction
    import scipy.sparse as sp

    reduction = ActiveDofReduction(ActiveMap.all_active(3))
    assert reduction.is_identity
    matrix = sp.eye(6, format="csr")
    reduced, rhs = reduction.restrict(matrix, np.arange(6.0))
    assert reduced is matrix                    # KOPYALANMIR
    assert list(rhs) == list(np.arange(6.0))


# ═══════════════════════════════════════════════════════ quyular (§4)
def test_a_perforation_in_an_inactive_cell_is_disabled():
    scal = default_scal()
    model = _built(_model(6, inactive=(5,), scal=scal), producer_i=5,
                   scal=scal)
    connections = PeacemanWellModel().build_connections(model)

    assert not any(c.cell == 5 for c in connections), "WI sıfırlanmayıb"
    assert [c.well_name for c in connections] == ["INJ"]


def test_a_partly_inactive_well_warns_but_keeps_its_open_perforations():
    scal = default_scal()
    geology = _model(6, inactive=(5,), scal=scal)
    wells = [
        Well("INJ", WellType.INJECTOR, WellControl(ControlMode.BHP, 320.0),
             [Perforation(0, 0, 0)]),
        Well("PROD", WellType.PRODUCER, WellControl(ControlMode.BHP, 150.0),
             [Perforation(4, 0, 0), Perforation(5, 0, 0)]),
    ]
    model = ReservoirModelBuilder().build(
        geological_model=geology, wells=wells, scal=scal,
        initial=InitialConditions(datum_pressure=220.0,
                                  water_saturation=scal.swc))

    report = model.diagnose()
    assert any("qeyri-aktiv" in message for message in
               report.messages(Severity.WARNING))
    assert not report.messages(Severity.ERROR)

    cells = [c.cell for c in PeacemanWellModel().build_connections(model)]
    assert 4 in cells and 5 not in cells


def test_a_fully_inactive_well_is_an_error_not_a_warning():
    scal = default_scal()
    geology = _model(6, inactive=(5,), scal=scal)
    wells = [
        Well("INJ", WellType.INJECTOR, WellControl(ControlMode.BHP, 320.0),
             [Perforation(0, 0, 0)]),
        Well("PROD", WellType.PRODUCER, WellControl(ControlMode.BHP, 150.0),
             [Perforation(5, 0, 0)]),
    ]
    model = ReservoirModelBuilder().build(
        geological_model=geology, wells=wells, scal=scal,
        initial=InitialConditions(datum_pressure=220.0,
                                  water_saturation=scal.swc))

    errors = model.diagnose().messages(Severity.ERROR)
    assert any("BÜTÜN perforasiyalar" in message for message in errors)


# ══════════════════════════════════════════ UÇTAN-UCA BƏRABƏRLİK (§1-§4)
def test_padding_cells_marked_inactive_reproduce_the_small_model():
    """ƏSAS TEST: qeyri-aktiv "doldurucu" hüceyrələr nəticəyə TƏSİR ETMİR.

    İki model:
        A — 10 hüceyrəli, hamısı aktiv
        B — 14 hüceyrəli, son 4-ü qeyri-aktiv

    B-nin aktiv hissəsi A ilə eynidir (eyni həndəsə, eyni xassələr, eyni
    quyular), ona görə OOIP, hasilat və təzyiq BİRƏBİR üst-üstə
    düşməlidir. Fərq çıxsa, qeyri-aktiv hüceyrələr fizikaya sızır.
    """
    scal = default_scal()
    config = SimulationConfig(
        end_time=400.0,
        time_stepping=TimeSteppingConfig(max_dt=5.0),
        output=OutputConfig(snapshot_count=4))

    small = _built(_model(10, scal=scal), producer_i=9, scal=scal)
    padded = _built(_model(14, inactive=(10, 11, 12, 13), scal=scal),
                    producer_i=9, scal=scal)

    assert padded.ncell == 14 and padded.n_active == 10

    result_small = make_service(scal).run(small, config)
    result_padded = make_service(scal).run(padded, config)

    assert result_padded.ooip == pytest.approx(result_small.ooip, rel=1e-12)
    assert result_padded.steps == result_small.steps
    assert result_padded.final_recovery_factor == pytest.approx(
        result_small.final_recovery_factor, rel=1e-9)
    assert result_padded.series.average_pressure[-1] == pytest.approx(
        result_small.series.average_pressure[-1], rel=1e-9)
    assert result_padded.series.cumulative_oil[-1] == pytest.approx(
        result_small.series.cumulative_oil[-1], rel=1e-9)


def test_the_simulation_runs_and_stays_finite_with_inactive_cells():
    """Reduksiya olmasaydı matris SİNQULYAR olardı və həll NaN verərdi."""
    scal = default_scal()
    model = _built(_model(12, inactive=(3, 7), scal=scal), producer_i=11,
                   scal=scal)
    result = make_service(scal).run(model, short_config(end_time=200.0))

    assert result.converged
    assert np.all(np.isfinite(result.series.average_pressure))
    assert result.ooip > 0.0
    assert result.steps > 0


def test_inactive_cells_do_not_collapse_the_cfl_time_step():
    """PV = 0 olan hüceyrə CFL nisbətini 0 edərdi və simulyasiya
    dayanardı — dt yalnız AKTİV hüceyrələr üzrə hesablanır."""
    scal = default_scal()
    model = _built(_model(12, inactive=(6,), scal=scal), producer_i=11,
                   scal=scal)
    result = make_service(scal).run(model, short_config(end_time=200.0))

    assert result.converged
    assert result.steps < 10_000, "dt sıfıra yaxın qalıb"


def test_average_pressure_is_taken_over_active_cells_only():
    scal = default_scal()
    model = _built(_model(8, inactive=(4,), scal=scal), producer_i=7,
                   scal=scal)
    engine = ImpesEngine(model, short_config(), *_providers(scal))
    engine.pressure[4] = 1.0e6            # qeyri-aktiv "zibil" dəyər

    assert float(np.mean(model.active_values(engine.pressure))) < 1000.0


# ════════════════════════════════ FULLY IMPLICIT / Nyuton yolu (§1, §3)
def _newton(model, scal):
    from imex2d.simulation.implicit.newton import NewtonSolver
    from imex2d.simulation.implicit.residual import ResidualAssembler
    from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter

    grid = TwoPointFluxDiscretization().build(model)
    connections = PeacemanWellModel().build_connections(model)
    residual = ResidualAssembler(model, grid, connections,
                                 CoreyRelativePermeabilityAdapter(scal))
    return NewtonSolver(residual), residual


def _state(model):
    from imex2d.simulation.implicit.state import ReservoirState
    ic = model.initial_conditions
    return ReservoirState(np.full(model.ncell, ic.datum_pressure),
                          np.full(model.ncell, ic.water_saturation))


def test_newton_converges_with_inactive_cells():
    """Daraltma olmasaydı Jakobianın qeyri-aktiv sətirləri SIFIR olardı
    və xətti həll NaN qaytarardı."""
    from imex2d.simulation.implicit.newton import NewtonStatus

    scal = default_scal()
    model = _built(_model(12, inactive=(4, 8), scal=scal), producer_i=11,
                   scal=scal)
    solver, _ = _newton(model, scal)

    assert solver.reduction.reduced_size == 2 * model.n_active == 20
    result = solver.solve(_state(model), dt=5.0)

    assert result.status is NewtonStatus.CONVERGED
    assert np.all(np.isfinite(result.state.pressure))
    assert np.all(np.isfinite(result.state.water_saturation))


def test_inactive_cell_state_never_changes_during_newton():
    """Qeyri-aktiv hüceyrənin naməlumu YOXDUR — dəyəri toxunulmaz qalır."""
    scal = default_scal()
    model = _built(_model(10, inactive=(5,), scal=scal), producer_i=9,
                   scal=scal)
    solver, _ = _newton(model, scal)

    state = _state(model)
    before = float(state.pressure[5])
    result = solver.solve(state, dt=5.0)

    assert float(result.state.pressure[5]) == pytest.approx(before, abs=1e-12)


def test_fully_implicit_engine_runs_with_inactive_cells():
    from imex2d.simulation.implicit.engine import FullyImplicitEngine
    from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter

    scal = default_scal()
    model = _built(_model(12, inactive=(3, 7), scal=scal), producer_i=11,
                   scal=scal)
    config = SimulationConfig(end_time=120.0,
                              output=OutputConfig(snapshot_count=3))
    result = FullyImplicitEngine(
        model, config, CoreyRelativePermeabilityAdapter(scal)).run()

    assert np.all(np.isfinite(result.series.average_pressure))
    assert result.ooip > 0.0


def test_material_balance_excludes_inactive_cells():
    """Akkumulyasiya PV-yə vurulur, PV isə qeyri-aktiv hüceyrədə 0-dır —
    material balansının məxrəci qeyri-aktiv payı SAYMAMALIDIR (#12)."""
    scal = default_scal()
    full = _built(_model(8, scal=scal), producer_i=7, scal=scal)
    holed = _built(_model(8, inactive=(4,), scal=scal), producer_i=7, scal=scal)

    _, residual_full = _newton(full, scal)
    _, residual_holed = _newton(holed, scal)

    water_full, oil_full = residual_full.accumulation(
        _state(full), residual_full.fluid_state(_state(full)))
    water_holed, oil_holed = residual_holed.accumulation(
        _state(holed), residual_holed.fluid_state(_state(holed)))

    assert water_holed[4] == 0.0
    assert oil_holed[4] == 0.0
    assert oil_holed.sum() == pytest.approx(oil_full.sum() * 7.0 / 8.0,
                                            rel=1e-9)


# ═══════════════════════════════════════════════════════════ MPFA-O (§3)
def test_mpfa_o_refuses_inactive_cells_instead_of_ignoring_them():
    """MPFA-O interaction-region stensilindən hüceyrə çıxara bilmir —
    SƏSSİZCƏ davam etmək #13-ü MPFA yolunda gizlədərdi."""
    from imex2d.discretization.mpfa_o import MPFAODiscretization
    from imex2d.simulation.implicit.residual import ResidualAssembler
    from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter

    scal = default_scal()
    model = _built(_model(6, inactive=(3,), scal=scal, ny=3), producer_i=5,
                   scal=scal)
    grid = MPFAODiscretization().build(model)

    assert any("ACTNUM" in feature for feature in grid.unsupported_features)
    with pytest.raises(NotImplementedError, match="ACTNUM"):
        ResidualAssembler(model, grid, [],
                          CoreyRelativePermeabilityAdapter(scal))


def test_tpfa_reports_inactive_cells_in_its_warnings():
    scal = default_scal()
    model = _built(_model(6, inactive=(3,), scal=scal), producer_i=5,
                   scal=scal)
    warnings = TwoPointFluxDiscretization().build(model).warnings
    assert any("ACTNUM" in w for w in warnings)


def test_tpfa_emits_no_actnum_warning_when_all_cells_are_active():
    scal = default_scal()
    model = _built(_model(6, scal=scal), producer_i=5, scal=scal)
    warnings = TwoPointFluxDiscretization().build(model).warnings
    assert not any("ACTNUM" in w for w in warnings)
