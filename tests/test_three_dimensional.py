"""3D modelləşdirmə (A5) testləri."""

import matplotlib
matplotlib.use("Agg")

import numpy as np
from matplotlib.figure import Figure

from helpers import default_scal, make_service, short_config
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot, line_drive)
from imex2d.application.simulation_service import SimulationService
from imex2d.domain.grid import CartesianGrid
from imex2d.domain.initial import InitialConditions
from imex2d.rendering import renderers as R
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.initialization.equilibrium import EquilibriumInitializationProvider
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.simulation.well_model import PeacemanWellModel


def _model(nx=9, ny=9, nz=4, dz=5.0, top=2000.0, owc=None, kv=0.1,
           equilibrate=False, scal=None):
    scal = scal or default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=dz, porosity=0.22,
        permx_base=200.0, nz=nz, kv_over_kh=kv, top_depth=top)
    initial = InitialConditions(datum_depth=top, datum_pressure=250.0,
                                oil_water_contact=owc,
                                use_equilibration=equilibrate)
    return ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                         scal=scal, initial=initial)


# ── grid və həndəsə ───────────────────────────────────────────────────
def test_three_dimensional_cell_count():
    model = _model(nx=9, ny=9, nz=4)
    assert model.ncell == 9 * 9 * 4


def test_vertical_connections_are_created():
    model = _model(nx=5, ny=5, nz=3)
    conn = model.connections()
    vertical = int((conn.axis == 2).sum())
    assert vertical == 5 * 5 * 2


def test_layer_depths_increase_downward():
    model = _model(nx=4, ny=4, nz=5, dz=6.0, top=2000.0)
    depths = model.geometry.cell_depths().reshape(model.grid.shape)
    per_layer = depths.mean(axis=(1, 2))
    assert np.all(np.diff(per_layer) > 0)
    assert abs(per_layer[0] - 2003.0) < 1e-9


def test_permz_controls_vertical_transmissibility():
    """Kv/Kh azaldıqda şaquli üzlərin transmissivliyi mütənasib azalır."""
    high = _model(nz=3, kv=0.5)
    low = _model(nz=3, kv=0.05)
    conn = high.connections()
    vertical = conn.axis == 2
    t_high = TwoPointFluxDiscretization().build(high).transmissibility[vertical]
    t_low = TwoPointFluxDiscretization().build(low).transmissibility[vertical]
    ratio = float(np.mean(t_high / t_low))
    assert abs(ratio - 10.0) / 10.0 < 1e-6


def test_horizontal_transmissibility_unaffected_by_permz():
    high = _model(nz=3, kv=0.5)
    low = _model(nz=3, kv=0.05)
    conn = high.connections()
    horizontal = conn.axis < 2
    t_high = TwoPointFluxDiscretization().build(high).transmissibility[horizontal]
    t_low = TwoPointFluxDiscretization().build(low).transmissibility[horizontal]
    assert np.allclose(t_high, t_low)


# ── quyular ───────────────────────────────────────────────────────────
def test_wells_perforate_every_layer():
    grid = CartesianGrid(7, 7, 5)
    for pattern in (five_spot, line_drive):
        for well in pattern(grid):
            assert len(well.perforations) == 5
            assert sorted(p.k for p in well.perforations) == list(range(5))


def test_well_connections_created_per_perforation():
    model = _model(nx=5, ny=5, nz=4)
    connections = PeacemanWellModel().build_connections(model)
    assert len(connections) == 2 * 4
    assert len({c.cell for c in connections}) == 8


def test_two_dimensional_model_still_has_single_perforation():
    """Reqressiya qoruyucusu: nz = 1 halda davranış dəyişmir."""
    grid = CartesianGrid(9, 9, 1)
    for well in five_spot(grid):
        assert len(well.perforations) == 1
        assert well.perforations[0].k == 0


# ── cazibə və equilibration ───────────────────────────────────────────
def test_gravity_activates_automatically_in_three_dimensions():
    """Şaquli ölçü olan kimi dərinlik fərqi yaranır — cazibə üzvü açılır."""
    engine_3d = make_service(default_scal()).create_engine(
        _model(nz=4), short_config(end_time=10.0))
    engine_2d = make_service(default_scal()).create_engine(
        _model(nz=1), short_config(end_time=10.0))
    assert engine_3d._has_gravity
    assert not engine_2d._has_gravity


def test_equilibration_places_water_in_lower_layers():
    scal = default_scal()
    model = _model(nx=7, ny=7, nz=6, dz=5.0, top=2000.0, owc=2018.0,
                   equilibrate=True, scal=scal)
    state = EquilibriumInitializationProvider().initialize(model)
    sw = state.water_saturation.reshape(model.grid.shape)
    per_layer = sw.mean(axis=(1, 2))
    assert np.all(np.diff(per_layer) >= -1e-12), "Sw aşağı təbəqələrdə artmır"
    assert per_layer[0] < per_layer[-1]


def test_pressure_increases_with_layer_depth():
    model = _model(nx=5, ny=5, nz=5, owc=None, equilibrate=True)
    state = EquilibriumInitializationProvider().initialize(model)
    per_layer = state.pressure.reshape(model.grid.shape).mean(axis=(1, 2))
    assert np.all(np.diff(per_layer) > 0)


# ── simulyasiya ───────────────────────────────────────────────────────
def test_three_dimensional_simulation_runs():
    scal = default_scal()
    model = _model(nx=9, ny=9, nz=4, owc=2018.0, equilibrate=True, scal=scal)
    service = SimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(scal),
        initialization_provider=EquilibriumInitializationProvider())
    result = service.run(model, short_config(end_time=200.0))
    assert result.converged
    assert result.snapshots[-1].water_saturation.shape == model.grid.shape


def test_material_balance_holds_in_three_dimensions():
    scal = default_scal()
    model = _model(nx=9, ny=9, nz=4, scal=scal)
    engine = make_service(scal).create_engine(model, short_config(end_time=200.0))
    initial_water = float(np.sum(model.pore_volume() * engine.sw))
    result = engine.run()
    final_water = float(np.sum(model.pore_volume() * engine.sw))
    series = result.series
    injected = float(np.trapezoid(series.water_injection_rate, series.time))
    produced = float(np.trapezoid(series.water_rate, series.time)) * model.fluids.water_fvf
    error = abs((final_water - initial_water) - (injected - produced)) / max(injected, 1e-9)
    assert error < 0.005, f"Material balans xətası {error * 100:.3f} %"


def test_vertical_permeability_controls_gravity_segregation():
    """Homogen modeldə şaquli fərqi yaradan yeganə mexanizm cazibədir.

    Kv böyük olanda su aşağı təbəqələrə çökür və təbəqələr arasında
    doyumluluq fərqi artır. Kv kiçik olanda təbəqələr bir-birindən
    təcrid olunur; hər biri eyni xassələrə və eyni perforasiyaya malik
    olduğu üçün eyni cür davranır və fərq itir.
    """
    scal = default_scal()
    spreads, profiles = [], []
    for kv in (0.5, 0.01):
        model = _model(nx=9, ny=9, nz=4, kv=kv, scal=scal)
        result = make_service(scal).run(model, short_config(end_time=300.0))
        sw = result.snapshots[-1].water_saturation.reshape(model.grid.shape)
        per_layer = sw.mean(axis=(1, 2))
        spreads.append(float(np.ptp(per_layer)))
        profiles.append(per_layer)

    assert spreads[0] > spreads[1] * 5, \
        "Yüksək Kv-də cazibə seqreqasiyası güclənmədi"
    assert profiles[0][-1] > profiles[0][0], \
        "Su alt təbəqələrdə toplanmadı"
    assert spreads[1] < 0.01, "Şaquli təcrid təbəqələri bərabərləşdirmədi"


# ── vizuallaşdırma ────────────────────────────────────────────────────
def test_map_renderer_draws_each_layer():
    model = _model(nx=7, ny=7, nz=3)
    figure = Figure()
    axes = figure.subplots()
    for layer in range(3):
        assert R.MapRenderer().draw(figure=figure, ax=axes, model=model,
                                    property_key=R.PERMEABILITY,
                                    layer=layer) is not None


def test_cross_section_renderer_draws_both_orientations():
    model = _model(nx=7, ny=7, nz=3)
    figure = Figure()
    axes = figure.subplots()
    for axis in ("J", "I"):
        assert R.CrossSectionRenderer().draw(axes, figure, model, R.SATURATION,
                                             None, None, axis=axis,
                                             index=2) is not None
