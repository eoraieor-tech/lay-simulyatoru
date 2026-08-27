"""Rendering qatı Qt olmadan işləməlidir — bu, ayrılığın sübutudur."""

import matplotlib
matplotlib.use("Agg")

import numpy as np
from matplotlib.figure import Figure

from helpers import default_scal, five_spot_model, make_service, short_config
from imex2d.rendering import renderers as R
from imex2d.simulation.analytical import buckley_leverett


def _result():
    scal = default_scal()
    return scal, make_service(scal).run(five_spot_model(nx=15, ny=15, scal=scal),
                                        short_config(end_time=150.0))


def test_map_renderer_draws_without_qt():
    scal, result = _result()
    model = five_spot_model(nx=15, ny=15, scal=scal)
    figure = Figure()
    axes = figure.subplots()
    for key in (R.SATURATION, R.PRESSURE, R.PERMEABILITY, R.POROSITY):
        bar = R.MapRenderer().draw(axes, figure, model, key, result.snapshots[-1])
        assert bar is not None


def test_map_renderer_works_without_simulation_result():
    model = five_spot_model(nx=9, ny=9)
    figure = Figure()
    R.MapRenderer().draw(figure.subplots(), figure, model, R.SATURATION, None)


def test_production_curve_renderer_draws_without_qt():
    _, result = _result()
    figure = Figure()
    R.ProductionCurveRenderer().draw(figure.subplots(2, 2), result)


def test_scal_renderer_draws_without_qt():
    scal = default_scal()
    figure = Figure()
    R.ScalRenderer().draw(figure.subplots(1, 2), scal, 0.5, 3.0)


def test_repeated_draws_do_not_shrink_the_axes():
    """Reqressiya: hər çəkilişdə yeni colorbar oxu daraldırdı.

    `figure.colorbar(..., ax=ax)` əsas oxdan yer alır və şkala silinəndə
    həmin yeri qaytarmır. 60 çəkilişdən sonra ox eni 97 % kiçilirdi.
    Sabit `cax` ilə ox toxunulmaz qalmalıdır.
    """
    model = five_spot_model(nx=11, ny=11)
    figure = Figure()
    spec = figure.add_gridspec(1, 2, width_ratios=[40, 1])
    axes = figure.add_subplot(spec[0, 0])
    cax = figure.add_subplot(spec[0, 1])

    initial_width = axes.get_position().width
    for _ in range(60):
        R.MapRenderer().draw(axes, figure, model, R.PERMEABILITY, None, cax=cax)
    assert abs(axes.get_position().width - initial_width) < 1e-9

    for _ in range(30):
        R.CrossSectionRenderer().draw(axes, figure, model, R.SATURATION, None,
                                      axis="J", index=0, cax=cax)
    assert abs(axes.get_position().width - initial_width) < 1e-9


def test_comparison_renderer_handles_multiple_and_zero_runs():
    scal, result = _result()
    figure = Figure()
    axes = figure.subplots(2, 2)
    R.RunComparisonRenderer().draw(axes, [("RUN-001", result),
                                          ("RUN-002", result)])
    R.RunComparisonRenderer().draw(axes, [])   # boş siyahı çökməməlidir


def test_pvt_renderer_draws_without_qt():
    from imex2d.simulation.pvt.correlations import build_pvt_table
    figure = Figure()
    R.PvtRenderer().draw(figure.subplots(2, 2), build_pvt_table())


def test_validation_renderer_draws_without_qt():
    scal = default_scal()
    analytical = buckley_leverett(scal, 0.5, 3.0, 0.2, 60.0, 1000.0, 250.0)
    x = np.linspace(0, 900, 120)
    sw = np.full(120, scal.swc)
    figure = Figure()
    R.ValidationRenderer().draw(figure.subplots(), analytical, x, sw, 250.0, 120)
