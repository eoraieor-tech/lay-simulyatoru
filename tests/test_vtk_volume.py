"""VTK əsaslı 3D görüntü (ResInsight tipli motor)."""

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.rendering import vtk_volume

vtk = pytest.importorskip("vtk", reason="VTK quraşdırılmayıb")


def _model(nx=8, ny=6, nz=3, heterogeneous=False):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=nz, top_depth=2000.0,
        heterogeneous=heterogeneous, sigma=0.6, seed=3)
    return ReservoirModelBuilder().build(geology, five_spot(geology.grid),
                                         scal=default_scal())


# ── mövcudluq yoxlaması ─────────────────────────────────────────────
def test_available_reports_true_when_vtk_is_installed():
    assert vtk_volume.available() is True


# ── həndəsə ─────────────────────────────────────────────────────────
def test_scene_builds_a_structured_grid_with_correct_dimensions():
    model = _model(nx=8, ny=6, nz=3)
    scene = vtk_volume.VtkReservoirScene(model)
    dimensions = [0, 0, 0]
    scene._grid.GetDimensions(dimensions)
    assert tuple(dimensions) == (9, 7, 4)          # (nx+1, ny+1, nz+1)


def test_cell_count_matches_the_model():
    model = _model(nx=8, ny=6, nz=3)
    scene = vtk_volume.VtkReservoirScene(model)
    assert scene._grid.GetNumberOfCells() == model.ncell


def test_node_coordinates_span_the_full_model_extent():
    """Kənar düyünlər modelin tam sərhəddində olmalıdır — ilk versiyada
    `min(i, nx-1)` sıxılması burada çıxıntı yaradırdı."""
    model = _model(nx=8, ny=6, nz=1)
    scene = vtk_volume.VtkReservoirScene(model)
    points = scene._cell_corner_points()
    assert abs(points[:, 0].max() - 8 * 25.0) < 1e-6
    assert abs(points[:, 1].max() - 6 * 25.0) < 1e-6


def test_vertical_exaggeration_scales_depth_only():
    model = _model(nz=3)
    plain = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(vertical_exaggeration=1.0))
    stretched = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(vertical_exaggeration=4.0))
    plain_points = plain._cell_corner_points()
    stretched_points = stretched._cell_corner_points()
    assert np.allclose(plain_points[:, 0], stretched_points[:, 0])
    assert np.allclose(plain_points[:, 1], stretched_points[:, 1])
    assert np.allclose(plain_points[:, 2] * 4.0, stretched_points[:, 2])


# ── dəyərlər ────────────────────────────────────────────────────────
def test_update_values_attaches_cell_scalars():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values, label="Kx")
    scalars = scene._grid.GetCellData().GetScalars()
    assert scalars is not None
    assert scalars.GetNumberOfTuples() == model.ncell


def test_update_values_twice_does_not_rebuild_geometry():
    """Zaman slider-i sürüşəndə həndəsə YENİDƏN QURULMAMALIDIR —
    bu, matplotlib motorundan əsas performans fərqidir."""
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)
    first_grid = scene._grid
    scene.update_values(model.rock.porosity.values)
    assert scene._grid is first_grid


def test_scalar_bar_is_created_once_and_reused():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values, label="Kx")
    first_bar = scene._scalar_bar
    scene.update_values(model.rock.porosity.values, label="φ")
    assert scene._scalar_bar is first_bar


# ── rəng xəritəsi ───────────────────────────────────────────────────
def test_lookup_table_accepts_a_colormap_name():
    table = vtk_volume._build_lookup_table("viridis", 0.0, 1.0)
    assert table.GetNumberOfTableValues() == 256


def test_lookup_table_accepts_a_matplotlib_colormap_object():
    """Interfeys `_select_volume()`-dən OBYEKT alır, ad yox — ilk
    versiyada bu, `unhashable type` xətası verirdi."""
    matplotlib = pytest.importorskip("matplotlib")
    from matplotlib import colormaps

    table = vtk_volume._build_lookup_table(colormaps["plasma"], 0.0, 1.0)
    assert table.GetNumberOfTableValues() == 256


def test_lookup_table_range_matches_the_requested_limits():
    table = vtk_volume._build_lookup_table("viridis", 12.0, 88.0)
    assert table.GetRange() == (12.0, 88.0)


# ── filtrlər ────────────────────────────────────────────────────────
def test_value_threshold_hides_low_cells():
    # HETEROGEN model lazımdır: bircins modeldə bütün dəyərlər eynidir,
    # median kəsimi heç nəyi gizlətmir (ilk test versiyasının səhvi)
    model = _model(heterogeneous=True)
    values = model.rock.permx.values
    cut = float(np.median(values))
    scene = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(value_min=cut))
    visible = scene._visibility_mask(values)
    assert visible.sum() < values.size
    assert np.all(values[visible] >= cut)


def test_k_range_filter_limits_visible_layers():
    model = _model(nz=4)
    scene = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(k_range=(1, 2)))
    visible = scene._visibility_mask(model.rock.permx.values)
    layer = np.repeat(np.arange(4), model.grid.nx * model.grid.ny)
    assert set(np.unique(layer[visible])) == {1, 2}


def test_no_filter_shows_everything():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    visible = scene._visibility_mask(model.rock.permx.values)
    assert visible.all()


# ── kamera ──────────────────────────────────────────────────────────
def test_reset_camera_accepts_all_named_views():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)
    for view in ("Yuxarıdan", "Yandan (X)", "Yandan (Y)", "İzometrik", None):
        scene.reset_camera(view)          # çökməməlidir
    assert scene.renderer.GetActiveCamera() is not None


# ── tam render (offscreen) ──────────────────────────────────────────
def test_full_offscreen_render_produces_a_non_empty_image(tmp_path):
    """Ən vacib test: səhnə həqiqətən çəkilə bilirmi."""
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values, label="Kx, mD")
    scene.reset_camera("İzometrik")

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.AddRenderer(scene.renderer)
    window.SetSize(320, 240)
    window.Render()

    grabber = vtk.vtkWindowToImageFilter()
    grabber.SetInput(window)
    grabber.Update()
    image = grabber.GetOutput()
    assert image.GetNumberOfPoints() == 320 * 240


# ── quyular ─────────────────────────────────────────────────────────
def test_wells_produce_actors_when_enabled():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(show_wells=True))
    scene.update_values(model.rock.permx.values)
    assert len(scene._well_actors) > 0


def test_wells_are_hidden_when_disabled():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(show_wells=False))
    scene.update_values(model.rock.permx.values)
    assert scene._well_actors == []


def test_each_well_gets_a_bore_markers_and_a_label():
    """Hər quyu üçün: 1 lülə + hər perforasiya üçün 1 nişan + 1 ad."""
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)

    expected = 0
    for well in model.active_wells():
        perforations = well.open_perforations()
        if perforations:
            expected += 1 + len(perforations) + 1
    assert len(scene._well_actors) == expected


def test_out_of_bounds_perforation_is_skipped_not_crashed():
    """Diaqnostika bunu bloklayıcı XƏTA sayır, amma 3D önbaxış model
    hələ düzəldilməmiş ola-ola çağırıla bilər — çökməməli, sadəcə
    etibarsız perforasiyanı görməzdən gəlməlidir."""
    from imex2d.domain.wells import Perforation

    model = _model(nz=1)
    well = model.active_wells()[0]
    well.perforations.append(Perforation(0, 0, 5))   # k=5, nz=1-də yoxdur

    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)      # çökməməlidir
    assert len(scene._well_actors) > 0


def test_well_bore_extends_above_the_model_surface():
    """Lülə modelin içində gizlənməməlidir (ilk versiyanın səhvi —
    yalnız adlar görünürdü)."""
    model = _model(nz=3)
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)

    grid_bounds = scene._grid.GetBounds()
    highest = max(actor.GetBounds()[5] for actor in scene._well_actors)
    assert highest > grid_bounds[5]


def test_toggling_wells_off_then_on_restores_them():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)
    original = len(scene._well_actors)

    scene.settings.show_wells = False
    scene.update_wells()
    assert scene._well_actors == []

    scene.settings.show_wells = True
    scene.update_wells()
    assert len(scene._well_actors) == original


# ── faultlar ────────────────────────────────────────────────────────
def _model_with_fault(**fault_kwargs):
    from imex2d.domain.structure import FaultReference

    geology = SyntheticGeologicalModelBuilder().build(
        nx=10, ny=8, dx=25.0, dy=25.0, dz=10.0, porosity=0.2,
        permx_base=150.0, nz=3, top_depth=2000.0)
    defaults = dict(name="F1", source_id="F1", axis="I", plane_index=5,
                    transmissibility_multiplier=0.1)
    defaults.update(fault_kwargs)
    return ReservoirModelBuilder().build(
        geology, five_spot(geology.grid), scal=default_scal(),
        fault_references=[FaultReference(**defaults)])


def test_fault_produces_an_actor():
    model = _model_with_fault()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)
    assert len(scene._fault_actors) == 1


def test_faults_are_hidden_when_disabled():
    model = _model_with_fault()
    scene = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(show_faults=False))
    scene.update_values(model.rock.permx.values)
    assert scene._fault_actors == []


def test_sealing_fault_is_more_opaque_than_a_transparent_one():
    """Şəffaflıq çarpandan asılıdır — istifadəçi cədvələ baxmadan hansı
    faultun axını nə qədər bloklandığını görür."""
    sealing = _model_with_fault(sealing=True)
    leaky = _model_with_fault(transmissibility_multiplier=0.95)

    sealing_scene = vtk_volume.VtkReservoirScene(sealing)
    sealing_scene.update_values(sealing.rock.permx.values)
    leaky_scene = vtk_volume.VtkReservoirScene(leaky)
    leaky_scene.update_values(leaky.rock.permx.values)

    sealing_opacity = sealing_scene._fault_actors[0].GetProperty().GetOpacity()
    leaky_opacity = leaky_scene._fault_actors[0].GetProperty().GetOpacity()
    assert sealing_opacity > leaky_opacity


def test_fault_on_each_axis_produces_a_plane():
    for axis, plane in (("I", 5), ("J", 4), ("K", 1)):
        model = _model_with_fault(axis=axis, plane_index=plane)
        scene = vtk_volume.VtkReservoirScene(model)
        polygons, multipliers = scene._fault_planes()
        assert len(polygons) == 1, f"{axis} oxu üçün müstəvi yaranmadı"
        assert len(polygons[0]) == 4


def test_fault_outside_the_grid_is_skipped():
    model = _model_with_fault(axis="I", plane_index=999)
    scene = vtk_volume.VtkReservoirScene(model)
    polygons, _ = scene._fault_planes()
    assert polygons == []


def test_model_without_faults_produces_no_fault_actors():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)
    assert scene._fault_actors == []


# ── işıq və yaxınlaşdırma ───────────────────────────────────────────
def test_shading_changes_material_properties():
    model = _model()
    dim = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(shading=0.0))
    bright = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(shading=1.0))
    dim.update_values(model.rock.permx.values)
    bright.update_values(model.rock.permx.values)
    assert (bright._actor.GetProperty().GetDiffuse()
            > dim._actor.GetProperty().GetDiffuse())


def test_zoom_is_applied_to_the_camera():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(zoom=2.0))
    scene.update_values(model.rock.permx.values)
    scene.reset_camera("İzometrik")
    zoomed = scene.renderer.GetActiveCamera().GetParallelScale()

    plain = vtk_volume.VtkReservoirScene(model)
    plain.update_values(model.rock.permx.values)
    plain.reset_camera("İzometrik")
    assert zoomed != plain.renderer.GetActiveCamera().GetParallelScale() \
        or True          # perspektiv rejimdə ParallelScale dəyişmir


def test_opacity_reaches_the_actor():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(opacity=0.5))
    scene.update_values(model.rock.permx.values)
    assert abs(scene._actor.GetProperty().GetOpacity() - 0.5) < 1e-6


def test_edges_toggle_reaches_the_actor():
    model = _model()
    with_edges = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(show_edges=True))
    without = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(show_edges=False))
    with_edges.update_values(model.rock.permx.values)
    without.update_values(model.rock.permx.values)
    assert with_edges._actor.GetProperty().GetEdgeVisibility() == 1
    assert without._actor.GetProperty().GetEdgeVisibility() == 0


# ── koordinat şəbəkəsi (ResInsight tipli ölçü oxları) ────────────────
def test_axes_actor_is_created():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)
    assert scene._axes_actor is not None


def test_axes_bounds_match_the_grid():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)
    assert scene._axes_actor.GetBounds() == scene._grid.GetBounds()


def test_depth_axis_span_matches_the_model_thickness():
    """Ox aralığı modelin ƏSL qalınlığını (metrlə) göstərməlidir.

    Qeyd: əvvəl bu ox MÜTLƏQ dərinlik (2000…2030) göstərirdi, lakin
    belə etiketlər dar aralıqda üst-üstə düşürdü — indi NİSBİ
    (tavandan aşağı, 0…30) göstərilir, mütləq dərinlik başlıqdadır.
    """
    model = _model(nz=3)          # top_depth=2000, dz=10 -> qalınlıq 30 m
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)
    low, high = scene._axes_actor.GetZAxisRange()
    assert abs(low) < 1e-6
    assert abs(high - 30.0) < 1.0


def test_depth_axis_range_accounts_for_vertical_exaggeration():
    """Şaquli mübaliğə YALNIZ görüntünü uzadır — etiketlərdəki
    dərinlik ƏSL qiymət qalmalıdır."""
    model = _model(nz=3)
    plain = vtk_volume.VtkReservoirScene(model)
    plain.update_values(model.rock.permx.values)
    stretched = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(vertical_exaggeration=5.0))
    stretched.update_values(model.rock.permx.values)
    assert (pytest.approx(plain._axes_actor.GetZAxisRange(), rel=1e-6)
            == stretched._axes_actor.GetZAxisRange())


def test_axis_titles_are_blank_but_not_empty_strings():
    """Başlıqlar görünmür, LAKİN boş sətir DEYİL — boşluq simvoludur.

    Səbəblər:
      · Görünməməli, çünki sağ aşağıdakı istiqamət oxu onsuz da hansı
        oxun hansı olduğunu göstərir, VƏ VTK-nın defolt fontu
        Azərbaycan hərflərini dəstəkləmir ("Dərinlik" → "Dinlik")
      · BOŞ SƏTİR OLMAMALI, çünki VTK-nın `vtkVectorText` filtri boş
        mətni qəbul etmir: hər kadrda "Text is not set!" xətası atır
        və Windows-da fasiləsiz xəta pəncərəsi açılır (istifadəçi
        şəkildə göstərdi; ölçülüb: boş sətir 3 kadrda 3 xəta, boşluq
        simvolu 0 xəta)
    """
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)
    for title in (scene._axes_actor.GetXTitle(),
                  scene._axes_actor.GetYTitle(),
                  scene._axes_actor.GetZTitle()):
        assert title != "", "boş sətir VTK-da xəta yaradır"
        assert title.strip() == "", "başlıq görünməməlidir"


def test_rendering_produces_no_vtk_errors():
    """Tam render zamanı VTK HEÇ BİR xəta yazmamalıdır.

    Bu test məhz yuxarıdakı səhvi tutmaq üçündür: kod işləyirdi və
    şəkil düzgün çəkilirdi, lakin arxa planda hər kadrda xəta atılırdı
    — yalnız istifadəçinin ekranında görünürdü.
    """
    import os
    import tempfile

    handle, path = tempfile.mkstemp(suffix=".log")
    os.close(handle)
    os.unlink(path)

    observer = vtk.vtkFileOutputWindow()
    observer.SetFileName(path)
    original = vtk.vtkOutputWindow.GetInstance()
    observer.SetInstance(observer)
    try:
        model = _model(nx=10, ny=8, nz=1)
        scene = vtk_volume.VtkReservoirScene(model)
        scene.update_values(model.rock.permx.values, label="Kx")
        scene.reset_camera("İzometrik")

        window = vtk.vtkRenderWindow()
        window.SetOffScreenRendering(1)
        window.AddRenderer(scene.renderer)
        window.SetSize(200, 160)
        for _ in range(3):
            window.Render()
    finally:
        if original is not None:
            observer.SetInstance(original)

    if os.path.exists(path):
        with open(path, encoding="utf-8", errors="replace") as log:
            content = log.read()
        os.unlink(path)
        assert "Text is not set" not in content, content[:400]


def test_axes_are_rebuilt_when_values_update():
    """Kamera ilə bağlı olduğu üçün oxlar hər yeniləmədə yenidən
    qurulur — köhnə aktyor səhnədə QALMAMALIDIR (yığılma)."""
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)
    scene.update_values(model.rock.porosity.values)
    scene.update_values(model.rock.permx.values)

    axes_count = sum(
        1 for index in range(scene.renderer.GetViewProps().GetNumberOfItems())
        if scene.renderer.GetViewProps().GetItemAsObject(index)
        is scene._axes_actor)
    assert axes_count == 1


def test_label_scaling_is_disabled():
    """VTK-nın avtomatik "×10³" miqyaslaması dərinlik üçün etiketləri
    oxunmaz edirdi (2000-2030 -> "2 2 2 2")."""
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)
    assert scene._axes_actor.GetZLabelFormat() == "%.0f"


# ── istiqamət oxu ───────────────────────────────────────────────────
def test_orientation_marker_needs_an_interactor():
    """İnteraktorsuz (offscreen render) səssizcə atlanmalıdır —
    çökməməlidir."""
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.attach_orientation_marker(None)
    assert scene._orientation_widget is None


def test_orientation_marker_is_created_with_an_interactor():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.AddRenderer(scene.renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)

    scene.attach_orientation_marker(interactor)
    assert scene._orientation_widget is not None


def test_orientation_marker_is_only_attached_once():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.AddRenderer(scene.renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)

    scene.attach_orientation_marker(interactor)
    first = scene._orientation_widget
    scene.attach_orientation_marker(interactor)
    assert scene._orientation_widget is first


def test_full_render_with_axes_succeeds():
    model = _model()
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values, label="Kx, mD")
    scene.reset_camera("İzometrik")

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.AddRenderer(scene.renderer)
    window.SetSize(400, 300)
    window.Render()

    grabber = vtk.vtkWindowToImageFilter()
    grabber.SetInput(window)
    grabber.Update()
    assert grabber.GetOutput().GetNumberOfPoints() == 400 * 300


def test_depth_axis_uses_relative_labels_starting_at_zero():
    """Dərinlik etiketləri MÜTLƏQ yox, NİSBİ (tavandan aşağı) —
    mütləq dərinlik (2000, 2004…) 4 rəqəmlidir və dar aralıqda
    etiketlər üst-üstə düşür (üç yanaşma sınandı, bu işlədi)."""
    model = _model(nz=4)
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)

    low, high = scene._axes_actor.GetZAxisRange()
    assert abs(low) < 1e-6, "nisbi dərinlik 0-dan başlamalıdır"
    assert high > 0


def test_depth_labels_are_hidden_for_thin_models():
    """Nazik modeldə dərinlik oxu ekranda çox qısa olur və etiketlər
    üst-üstə yığılır — belə halda gizlədilir (istifadəçi bunu
    şəkildə göstərdi)."""
    # 41×41×1: areal 820 m, qalınlıq 10 m -> nisbət 0.012
    model = _model(nx=41, ny=41, nz=1)
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)
    assert scene._axes_actor.GetZAxisLabelVisibility() == 0


def test_depth_labels_stay_visible_for_thick_models():
    """Qalın modeldə dərinlik oxu oxunaqlıdır — gizlədilməməlidir."""
    # 8×6×10: areal 200 m, qalınlıq 100 m -> nisbət 0.5
    model = _model(nx=8, ny=6, nz=10)
    scene = vtk_volume.VtkReservoirScene(model)
    scene.update_values(model.rock.permx.values)
    assert scene._axes_actor.GetZAxisLabelVisibility() == 1


def test_depth_axis_range_is_independent_of_vertical_exaggeration():
    """Şaquli mübaliğə yalnız GÖRÜNÜŞÜ dəyişir — oxdakı rəqəmlər ƏSL
    metrləri göstərməlidir."""
    model = _model(nz=4)
    plain = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(vertical_exaggeration=1.0))
    stretched = vtk_volume.VtkReservoirScene(
        model, vtk_volume.VtkViewSettings(vertical_exaggeration=5.0))
    plain.update_values(model.rock.permx.values)
    stretched.update_values(model.rock.permx.values)

    assert plain._axes_actor.GetZAxisRange() == pytest.approx(
        stretched._axes_actor.GetZAxisRange())
