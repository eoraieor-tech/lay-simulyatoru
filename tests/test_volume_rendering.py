"""3D həcm görüntüsü."""

import matplotlib
matplotlib.use("Agg")

import numpy as np
from matplotlib.figure import Figure

from helpers import default_scal
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.rendering.volume import VolumeFilter, VolumeRenderer


def _model(nx=8, ny=8, nz=4):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=5.0, porosity=0.2,
        permx_base=200.0, nz=nz, top_depth=2000.0)
    return ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                         scal=default_scal())


def _axes():
    figure = Figure()
    spec = figure.add_gridspec(1, 2, width_ratios=[40, 1])
    return figure, figure.add_subplot(spec[0, 0], projection="3d"), \
        figure.add_subplot(spec[0, 1])


def _faces(model, values, volume_filter=None):
    """(çoxbucaqlılar, dəyərlər) — parlaqlıq ayrıca testlərdə yoxlanılır."""
    mask = (volume_filter or VolumeFilter()).mask(values, model.grid.shape)
    polygons, face_values, _ = VolumeRenderer()._visible_faces(
        model, values, mask, 1.0)
    return polygons, face_values


# ── görünən üz çıxarışı ───────────────────────────────────────────────
def test_only_boundary_faces_are_drawn_without_filtering():
    """Filtrsiz halda yalnız modelin xarici səthi çəkilməlidir.

    nx×ny×nz bloku üçün xarici üz sayı analitikdir:
        2(nx·ny + ny·nz + nx·nz)
    """
    model = _model(nx=6, ny=5, nz=4)
    values = model.rock.permx.values
    polygons, _ = _faces(model, values)

    nz, ny, nx = model.grid.shape
    expected = 2 * (nx * ny + ny * nz + nx * nz)
    assert len(polygons) == expected


def test_face_extraction_scales_far_below_the_naive_count():
    model = _model(nx=20, ny=20, nz=10)
    polygons, _ = _faces(model, model.rock.permx.values)
    naive = model.ncell * 6
    assert len(polygons) < naive / 10, f"{len(polygons)} vs {naive}"


def test_every_face_has_four_corners_in_three_dimensions():
    model = _model()
    polygons, values = _faces(model, model.rock.permx.values)
    assert len(values) == len(polygons)
    for polygon in polygons[:20]:
        assert np.shape(polygon) == (4, 3)


def test_face_values_come_from_the_owning_cell():
    """Üzün rəngi öz hüceyrəsinin dəyəri olmalıdır."""
    model = _model(nx=4, ny=4, nz=2)
    values = np.arange(model.ncell, dtype=float)
    _, face_values = _faces(model, values)
    assert set(np.unique(face_values)).issubset(set(values))


# ── filtr və kəsim ────────────────────────────────────────────────────
def test_threshold_filter_hides_low_values():
    model = _model(nx=8, ny=8, nz=4)
    values = np.linspace(0.0, 1.0, model.ncell)
    mask = VolumeFilter(value_min=0.5).mask(values, model.grid.shape)
    assert mask.sum() < model.ncell
    assert np.all(values.reshape(model.grid.shape)[mask] >= 0.5)


def test_filtering_exposes_interior_faces():
    """Hüceyrə gizlədiləndə qonşusunun üzü açılmalıdır — daxili struktur."""
    model = _model(nx=8, ny=8, nz=4)
    values = np.random.default_rng(0).random(model.ncell)
    outside = _faces(model, values)[0]
    filtered = _faces(model, values, VolumeFilter(value_min=0.5))[0]
    assert len(filtered) > len(outside) * 0.5
    assert len(filtered) != len(outside)


def test_layer_range_limits_the_visible_volume():
    model = _model(nx=6, ny=6, nz=5)
    values = model.rock.permx.values
    mask = VolumeFilter(k_range=(1, 2)).mask(values, model.grid.shape)
    per_layer = mask.sum(axis=(1, 2))
    assert per_layer[0] == 0 and per_layer[3] == 0
    assert per_layer[1] > 0 and per_layer[2] > 0


def test_single_layer_slice_produces_a_flat_sheet():
    model = _model(nx=5, ny=5, nz=4)
    values = model.rock.permx.values
    polygons, _ = _faces(model, values, VolumeFilter(k_range=(2, 2)))
    nz, ny, nx = model.grid.shape
    expected = 2 * (nx * ny + ny * 1 + nx * 1)
    assert len(polygons) == expected


def test_empty_filter_produces_no_faces():
    model = _model()
    values = model.rock.permx.values
    polygons, _ = _faces(model, values,
                         VolumeFilter(value_min=values.max() * 10))
    assert polygons == []


# ── çəkmə ─────────────────────────────────────────────────────────────
def test_renderer_draws_without_qt():
    model = _model()
    figure, axes, cax = _axes()
    bar = VolumeRenderer().draw(figure=figure, ax=axes, model=model,
                                values=model.rock.permx.values,
                                label="Kx", cax=cax)
    assert bar is not None
    assert axes.collections


def test_renderer_handles_every_property():
    from imex2d.rendering import renderers as R

    model = _model()
    figure, axes, cax = _axes()
    renderer = VolumeRenderer()
    for key in (R.PERMEABILITY, R.POROSITY, R.DEPTH):
        values, colormap, low, high = R.MapRenderer()._select_volume(
            model, key, None)
        renderer.draw(axes, figure, model, np.asarray(values).ravel(),
                      label=key, colormap=colormap, cax=cax)
    assert axes.collections


def test_repeated_draws_do_not_shrink_the_axes():
    """Colorbar reqressiyası — 2D-də olduğu kimi 3D-də də yoxlanılır."""
    model = _model()
    figure, axes, cax = _axes()
    renderer = VolumeRenderer()
    width = axes.get_position().width
    for _ in range(20):
        renderer.draw(axes, figure, model, model.rock.permx.values,
                      label="Kx", cax=cax)
    assert abs(axes.get_position().width - width) < 1e-9


def test_coordinates_are_always_true_depths():
    """Mübaliğə KOORDİNATLARA toxunmamalıdır.

    Əvvəl hüceyrələr `dz × exaggeration` hündürlüyündə çəkilirdi, ox
    etiketləri isə həqiqi dərinliyi göstərirdi — ikisi bir-birinə zidd
    idi. 1599 m tavanlı, 100 m qalınlıqlı lay `Z×3` ilə 1349–1949 m
    kimi görünürdü.
    """
    model = _model(nx=5, ny=5, nz=3)
    values = model.rock.permx.values
    mask = VolumeFilter().mask(values, model.grid.shape)
    renderer = VolumeRenderer()

    flat = np.array(renderer._visible_faces(model, values, mask, 1.0)[0])
    stretched = np.array(renderer._visible_faces(model, values, mask, 8.0)[0])
    assert np.allclose(flat, stretched)


def test_depth_axis_matches_the_real_layer_interval():
    """Ox modelin həqiqi dərinlik intervalını göstərməlidir."""
    model = _model(nx=4, ny=4, nz=1)
    figure, axes, cax = _axes()
    VolumeRenderer().draw(axes, figure, model, model.rock.permx.values,
                          vertical_exaggeration=5.0, cax=cax)

    depths = model.geometry.cell_depths()
    half = model.geometry.dz * 0.5
    top, base = depths.min() - half, depths.max() + half

    low, high = sorted(axes.get_zlim())
    span = base - top
    assert abs(low - top) < span * 0.2
    assert abs(high - base) < span * 0.2


def test_exaggeration_changes_only_the_box_aspect():
    model = _model(nx=6, ny=6, nz=2)
    figure, axes, cax = _axes()
    renderer = VolumeRenderer()

    aspects, limits = [], []
    for exaggeration in (1.0, 6.0):
        renderer.draw(axes, figure, model, model.rock.permx.values,
                      vertical_exaggeration=exaggeration, cax=cax)
        aspects.append(axes.get_box_aspect()[2])
        limits.append(sorted(axes.get_zlim()))

    assert aspects[1] > aspects[0] * 3          # görüntü uzanır
    assert np.allclose(limits[0], limits[1])    # ox dəyişmir


def test_wells_are_drawn_at_their_grid_positions():
    model = _model(nx=6, ny=6, nz=3)
    figure, axes, cax = _axes()
    VolumeRenderer().draw(axes, figure, model, model.rock.permx.values,
                          show_wells=True, cax=cax)
    assert axes.lines, "Quyu lüləsi çəkilmədi"


def test_wells_can_be_hidden():
    model = _model(nx=6, ny=6, nz=3)
    figure, axes, cax = _axes()
    VolumeRenderer().draw(axes, figure, model, model.rock.permx.values,
                          show_wells=False, cax=cax)
    assert not axes.lines


def test_two_dimensional_model_still_renders():
    """nz = 1 halda 3D görüntü tək təbəqəli lövhə olmalıdır."""
    model = _model(nx=6, ny=6, nz=1)
    figure, axes, cax = _axes()
    VolumeRenderer().draw(axes, figure, model, model.rock.permx.values,
                          cax=cax)
    assert axes.collections


def test_simulation_snapshot_can_be_rendered():
    from helpers import make_service, short_config
    from imex2d.domain.wells import ControlMode, WellControl

    scal = default_scal()
    model = _model(nx=7, ny=7, nz=3)
    rate = model.pore_volume()[0] * 0.05
    for well in model.wells:
        well.control = (WellControl(ControlMode.RATE, rate) if well.is_injector
                        else WellControl(ControlMode.BHP, 150.0))
    result = make_service(scal).run(model, short_config(end_time=100.0))

    figure, axes, cax = _axes()
    VolumeRenderer().draw(axes, figure, model,
                          result.snapshots[-1].water_saturation.ravel(),
                          label="Sw", value_limits=(scal.swc, 1 - scal.sor),
                          cax=cax)
    assert axes.collections


# ── işıqlandırma və görünüş ───────────────────────────────────────────
def test_face_shading_depends_on_orientation():
    """Altı üz istiqamətinin altı fərqli parlaqlığı olmalıdır.

    Bütün üzlər eyni rəngdə olanda 3D model yastı görünür — dərinlik
    hissi məhz kölgədən yaranır.
    """
    model = _model(nx=5, ny=5, nz=3)
    values = model.rock.permx.values
    mask = VolumeFilter().mask(values, model.grid.shape)
    _, _, shades = VolumeRenderer()._visible_faces(model, values, mask, 1.0)
    assert len(np.unique(np.round(shades, 6))) == 6


def test_shading_does_not_change_the_number_of_faces():
    model = _model()
    figure, axes, cax = _axes()
    counts = []
    for shading in (0.0, 1.0):
        axes.clear()
        VolumeRenderer().draw(axes, figure, model, model.rock.permx.values,
                              shading=shading, cax=cax)
        counts.append(len(axes.collections[0].get_paths())
                      if hasattr(axes.collections[0], "get_paths")
                      else len(axes.collections))
    assert counts[0] == counts[1]


def test_stronger_shading_darkens_the_faces():
    model = _model(nx=5, ny=5, nz=3)
    values = model.rock.permx.values
    mask = VolumeFilter().mask(values, model.grid.shape)
    renderer = VolumeRenderer()
    _, face_values, shades = renderer._visible_faces(model, values, mask, 1.0)

    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    norm = Normalize(float(np.min(face_values)), float(np.max(face_values)))
    mappable = ScalarMappable(norm=norm, cmap="viridis")

    flat = renderer._shaded_colours(mappable, face_values, norm, shades,
                                    0.0, 1.0)
    shaded = renderer._shaded_colours(mappable, face_values, norm, shades,
                                      0.9, 1.0)
    assert shaded[:, :3].mean() < flat[:, :3].mean()
    assert np.all(shaded[:, :3] >= 0.0)


def test_opacity_is_applied_to_the_alpha_channel():
    model = _model(nx=4, ny=4, nz=2)
    values = model.rock.permx.values
    mask = VolumeFilter().mask(values, model.grid.shape)
    renderer = VolumeRenderer()
    _, face_values, shades = renderer._visible_faces(model, values, mask, 1.0)

    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    norm = Normalize(float(np.min(face_values)), float(np.max(face_values)))
    mappable = ScalarMappable(norm=norm, cmap="viridis")
    colours = renderer._shaded_colours(mappable, face_values, norm, shades,
                                       0.45, 0.6)
    assert np.allclose(colours[:, 3], 0.6)


def test_view_presets_set_the_camera_angle():
    from imex2d.rendering.volume import VIEW_ANGLES

    model = _model()
    figure, axes, cax = _axes()
    renderer = VolumeRenderer()
    for name, (elevation, azimuth) in VIEW_ANGLES.items():
        renderer.draw(axes, figure, model, model.rock.permx.values,
                      view=name, cax=cax)
        assert abs(axes.elev - elevation) < 1e-6, name
        assert abs(axes.azim - azimuth) < 1e-6, name


def test_free_view_does_not_reset_the_camera():
    """`view=None` halda kamera toxunulmaz qalmalıdır."""
    model = _model()
    figure, axes, cax = _axes()
    renderer = VolumeRenderer()
    axes.view_init(elev=17.0, azim=-33.0)
    renderer.draw(axes, figure, model, model.rock.permx.values, view=None,
                  cax=cax)
    assert abs(axes.elev - 17.0) < 1e-6
    assert abs(axes.azim + 33.0) < 1e-6


# ── dərinlik oxu və şaquli miqyas ─────────────────────────────────────
def test_depth_ticks_sit_on_layer_boundaries():
    """Ox işarələri təbəqə sərhədlərində olmalıdır.

    İxtiyari addımlar (2050, 2150 …) təbəqələrlə üst-üstə düşmür və
    hansı hüceyrənin harada bitdiyini görmək olmur. 5 təbəqə × 80 m,
    tavan 2000 -> 2000, 2080, 2160, 2240, 2320, 2400.
    """
    import dataclasses

    model = _model(nx=6, ny=6, nz=5)
    model.geometry = dataclasses.replace(model.geometry, dz=80.0,
                                         top_depth=2000.0)
    figure, axes, cax = _axes()
    VolumeRenderer().draw(axes, figure, model, model.rock.permx.values,
                          cax=cax)
    ticks = [float(value) for value in axes.get_zticks()]
    assert ticks == [2000.0, 2080.0, 2160.0, 2240.0, 2320.0, 2400.0], ticks


def test_depth_ticks_are_thinned_for_many_layers():
    """20 təbəqədə hər sərhədi işarələmək oxu oxunmaz edir."""
    model = _model(nx=4, ny=4, nz=20)
    figure, axes, cax = _axes()
    VolumeRenderer().draw(axes, figure, model, model.rock.permx.values,
                          cax=cax)
    assert len(axes.get_zticks()) <= 12


def test_zoom_magnifies_without_changing_the_geometry():
    """Yaxınlaşdırma yalnız görüntü böyütməsidir.

    Ox hədləri, dərinlik işarələri və çəkilən çoxbucaqlılar
    dəyişməməlidir — yalnız model çərçivədə daha böyük görünür.
    """
    model = _model(nx=6, ny=6, nz=3)
    figure, axes, cax = _axes()
    renderer = VolumeRenderer()

    limits, ticks = [], []
    for zoom in (1.0, 2.5):
        renderer.draw(axes, figure, model, model.rock.permx.values,
                      zoom=zoom, cax=cax)
        limits.append(sorted(axes.get_zlim()))
        ticks.append([float(value) for value in axes.get_zticks()])

    assert np.allclose(limits[0], limits[1])
    assert ticks[0] == ticks[1]


def test_default_zoom_fills_the_frame():
    """100 % modelin çərçivəni doldurması deməkdir.

    matplotlib-in 3D oxu defolt olaraq geniş boş kənar buraxır — ona
    görə baza əmsalı tətbiq olunur.
    """
    from imex2d.rendering.volume import BASE_FIT

    assert BASE_FIT > 1.0


def test_unit_scale_gives_the_true_geometric_ratio():
    """Z× 1 həqiqi nisbət olmalıdır — heç bir gizli mübaliğə yoxdur."""
    model = _model(nx=10, ny=10, nz=4)
    figure, axes, cax = _axes()
    VolumeRenderer().draw(axes, figure, model, model.rock.permx.values,
                          vertical_exaggeration=1.0, cax=cax)
    aspect = axes.get_box_aspect()

    plan = max(model.grid.nx * model.geometry.dx,
               model.grid.ny * model.geometry.dy)
    thickness = float(model.geometry.dz.sum())
    assert abs(aspect[2] / aspect[0] - thickness / plan) < 1e-6




def test_zoom_scales_the_view_without_touching_the_data():
    """Yaxınlaşdırma yalnız kameraya təsir edir.

    Koordinatlar, ox hədləri və çəkilən üzlər dəyişməməlidir —
    yalnız model kadrda böyük görünür.
    """
    model = _model(nx=8, ny=8, nz=3)
    figure, axes, cax = _axes()
    renderer = VolumeRenderer()

    states = []
    for zoom in (1.0, 2.5):
        renderer.draw(axes, figure, model, model.rock.permx.values,
                      zoom=zoom, cax=cax)
        states.append((sorted(axes.get_zlim()), axes.get_xlim(),
                       len(axes.collections)))

    assert np.allclose(states[0][0], states[1][0])
    assert np.allclose(states[0][1], states[1][1])
    assert states[0][2] == states[1][2]


def test_default_zoom_is_neutral():
    """Defolt 100 % — heç bir yaxınlaşdırma tətbiq olunmur."""
    model = _model(nx=6, ny=6, nz=2)
    figure, axes, cax = _axes()
    renderer = VolumeRenderer()

    renderer.draw(axes, figure, model, model.rock.permx.values, cax=cax)
    implicit = axes.get_box_aspect()
    renderer.draw(axes, figure, model, model.rock.permx.values, zoom=1.0,
                  cax=cax)
    explicit = axes.get_box_aspect()
    assert np.allclose(implicit, explicit)


def test_zoom_is_clamped_to_a_positive_value():
    model = _model(nx=5, ny=5, nz=2)
    figure, axes, cax = _axes()
    VolumeRenderer().draw(axes, figure, model, model.rock.permx.values,
                          zoom=0.0, cax=cax)
    assert axes.collections


# ── faultlar (3D görüntüdə) ─────────────────────────────────────────────
def _model_with_fault(nx=10, ny=8, nz=3, **fault_kwargs):
    from imex2d.domain.structure import FaultReference

    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=15.0, porosity=0.2,
        permx_base=200.0, nz=nz, top_depth=2000.0)
    defaults = dict(name="F1", source_id="F1", axis="I", plane_index=5,
                    transmissibility_multiplier=0.1)
    defaults.update(fault_kwargs)
    return ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=default_scal(),
        fault_references=[FaultReference(**defaults)])


def test_fault_without_geometry_produces_no_polygons():
    from imex2d.domain.structure import FaultReference

    geology = SyntheticGeologicalModelBuilder().build(
        nx=8, ny=8, dx=25.0, dy=25.0, dz=15.0, porosity=0.2,
        permx_base=200.0, nz=3, top_depth=2000.0)
    model = ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=default_scal(),
        fault_references=[FaultReference(name="F1", source_id="F1")])
    polygons, multipliers, labels = VolumeRenderer()._fault_polygons(model, 1.0)
    assert polygons == []


def test_i_axis_fault_produces_a_vertical_plane_at_the_correct_x():
    model = _model_with_fault(axis="I", plane_index=5, range_a=(0, 7),
                              range_b=(0, 2))
    polygons, multipliers, labels = VolumeRenderer()._fault_polygons(model, 1.0)
    assert len(polygons) == 1
    quad = polygons[0]
    expected_x = 6 * 25.0             # (plane_index + 1) * dx
    assert np.allclose(quad[:, 0], expected_x)
    # Y oxu tam diapazonu (0..7) əhatə etməlidir
    assert abs(quad[:, 1].min() - 0.0) < 1e-6
    assert abs(quad[:, 1].max() - 8 * 25.0) < 1e-6


def test_j_axis_fault_produces_a_vertical_plane_at_the_correct_y():
    model = _model_with_fault(axis="J", plane_index=3, range_a=(0, 9),
                              range_b=(0, 2))
    polygons, _, _ = VolumeRenderer()._fault_polygons(model, 1.0)
    quad = polygons[0]
    expected_y = 4 * 25.0
    assert np.allclose(quad[:, 1], expected_y)


def test_k_axis_fault_produces_a_horizontal_plane():
    model = _model_with_fault(axis="K", plane_index=1, range_a=(0, 9),
                              range_b=(0, 7))
    polygons, _, _ = VolumeRenderer()._fault_polygons(model, 1.0)
    assert len(polygons) == 1
    quad = polygons[0]
    assert np.ptp(quad[:, 2]) < 1e-6      # bütün küncləri eyni dərinlikdə


def test_fault_out_of_grid_range_is_skipped():
    model = _model_with_fault(axis="I", plane_index=999)
    polygons, _, _ = VolumeRenderer()._fault_polygons(model, 1.0)
    assert polygons == []


def test_fault_label_sits_at_the_plane_centre():
    model = _model_with_fault(axis="I", plane_index=5, range_a=(0, 7),
                              range_b=(0, 2))
    _, _, labels = VolumeRenderer()._fault_polygons(model, 1.0)
    assert len(labels) == 1
    name, x, y, z = labels[0]
    assert name == "F1"
    assert abs(x - 6 * 25.0) < 1e-6
    assert abs(y - 4 * 25.0) < 1e-6       # (0+7+1)/2 * dy


def test_sealing_fault_reports_zero_effective_multiplier():
    model = _model_with_fault(axis="I", plane_index=5, sealing=True,
                              transmissibility_multiplier=0.8)
    _, multipliers, _ = VolumeRenderer()._fault_polygons(model, 1.0)
    assert multipliers[0] == 0.0


def test_multiple_faults_all_produce_polygons():
    from imex2d.domain.structure import FaultReference

    geology = SyntheticGeologicalModelBuilder().build(
        nx=10, ny=8, dx=25.0, dy=25.0, dz=15.0, porosity=0.2,
        permx_base=200.0, nz=3, top_depth=2000.0)
    model = ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=default_scal(),
        fault_references=[
            FaultReference(name="F1", source_id="F1", axis="I",
                          plane_index=3, transmissibility_multiplier=0.1),
            FaultReference(name="F2", source_id="F2", axis="J",
                          plane_index=2, sealing=True)])
    polygons, multipliers, labels = VolumeRenderer()._fault_polygons(model, 1.0)
    assert len(polygons) == 2
    assert {label[0] for label in labels} == {"F1", "F2"}


def test_draw_with_faults_does_not_raise():
    model = _model_with_fault()
    figure, axes, cax = _axes()
    VolumeRenderer().draw(axes, figure, model, model.rock.permx.values,
                          show_faults=True, cax=cax)
    assert axes.lines, "Fault konturu (Line3D) çəkilmədi"


def test_faults_can_be_hidden():
    model = _model_with_fault()
    figure, axes, cax = _axes()
    VolumeRenderer().draw(axes, figure, model, model.rock.permx.values,
                          show_wells=False, show_faults=False, cax=cax)
    assert not axes.lines


def test_model_without_faults_draws_no_extra_lines():
    model = _model()
    figure, axes, cax = _axes()
    VolumeRenderer().draw(axes, figure, model, model.rock.permx.values,
                          show_wells=False, show_faults=True, cax=cax)
    assert not axes.lines
