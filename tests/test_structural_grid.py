"""A2 — quyulardan gələn lay üstü/altı səthlərinin grid həndəsəsinə
qoşulması.

Testlər üç qatı əhatə edir:

  §1  `domain/structural_grid.py` — saf təpə qurma və yoxlama;
  §2  `application/geology_service.py` — `structure_from_wells` yolu,
      siyasətlər, hesabat və İKİ MƏRHƏLƏLİ qurma sırası;
  §3  FİZİKAYA TƏSİRİN SÜBUTU — məsamə həcmi, cazibə, perforasiya→K,
      equilibration, transmissibillik.

Ən vacib test qrupu §3-dür: dəyişikliyin HƏQİQƏTƏN işlədiyinin sübutu
"kod çağırıldı"-da yox, "rəqəm dəyişdi"-dədir.
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers import default_scal, make_service, short_config

from imex2d.application.geology_service import (GeologicalGridSpec,
                                                InterpolationReport,
                                                WellBasedGeologicalModelBuilder)
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.project import Project
from imex2d.application.scenarios import five_spot
from imex2d.application.serialization import ProjectSerializer
from imex2d.domain.corner_point_geometry import CornerPointGeometry, cartesian_nodes
from imex2d.domain.geological_model import GeologicalModel
from imex2d.domain.geometry import CellGeometry, depth_to_k, interval_layers
from imex2d.domain.grid import CartesianGrid
from imex2d.domain.initial import InitialConditions
from imex2d.domain.properties import PropertyMap
from imex2d.domain.structural_grid import (structural_nodes, structure_statistics,
                                           validate_structure)
from imex2d.domain.structure import RegionSet
from imex2d.domain.well_data import WellDataset, WellSample
from imex2d.domain.wells import (ControlMode, Perforation, Well, WellControl,
                                 WellType)
from imex2d.geology.interpolation import OrdinaryKriging
from imex2d.simulation.discretization import TwoPointFluxDiscretization
from imex2d.simulation.initialization.equilibrium import EquilibriumInitializationProvider


# ═══════════════════════════════════════════════ §1 saf təpə həndəsəsi
def test_flat_surfaces_reproduce_cartesian_nodes_exactly():
    """Sabit TOP/BOTTOM → `cartesian_nodes` ilə BİRƏBİR eyni massiv.

    Bu, təpə SIRASININ (`HEX_FACE_VERTEX_INDICES` konvensiyası) düzgün
    olmasının sübutudur — sıra səhv olsaydı massivlər fərqlənərdi."""
    grid = CartesianGrid(4, 3, 5)
    geometry = CellGeometry(grid, dx=20.0, dy=25.0, dz=8.0, top_depth=2000.0)
    reference = cartesian_nodes(grid, geometry)

    areal = grid.nx * grid.ny
    top = np.full(areal, 2000.0)
    bottom = np.full(areal, 2000.0 + 5 * 8.0)
    nodes = structural_nodes(grid, top, bottom, dx=20.0, dy=25.0)

    assert nodes.shape == (grid.ncell, 8, 3)
    assert np.array_equal(nodes, reference)


def test_flat_structure_matches_cartesian_geometry_to_machine_precision():
    """Degenerativ hal: qurulan CPG `CellGeometry` ilə eyni həcm,
    dərinlik və üz sahələri verir (rtol=1e-12)."""
    grid = CartesianGrid(5, 4, 3)
    geometry = CellGeometry(grid, dx=30.0, dy=30.0, dz=12.0, top_depth=1800.0)
    areal = grid.nx * grid.ny
    nodes = structural_nodes(grid, np.full(areal, 1800.0),
                             np.full(areal, 1800.0 + 3 * 12.0), 30.0, 30.0)
    cpg = CornerPointGeometry.from_nodes(grid, nodes)
    conn = grid.build_connections()

    assert np.allclose(cpg.volumes(), geometry.volumes(), rtol=1e-12)
    assert np.allclose(cpg.cell_depths(), geometry.cell_depths(), rtol=1e-12)
    assert np.allclose(cpg.face_areas(conn), geometry.face_areas(conn), rtol=1e-12)
    half_a, half_b = cpg.face_half_distances(conn)
    ref_a, ref_b = geometry.face_half_distances(conn)
    assert np.allclose(half_a, ref_a, rtol=1e-12)
    assert np.allclose(half_b, ref_b, rtol=1e-12)


def test_variable_thickness_gives_proportional_cell_volumes():
    """İki sütun (10 m və 30 m) → həcm nisbəti DƏQİQ 1:3."""
    grid = CartesianGrid(2, 1, 2)
    top = np.array([2000.0, 2000.0])
    bottom = np.array([2010.0, 2030.0])
    cpg = CornerPointGeometry.from_nodes(
        grid, structural_nodes(grid, top, bottom, 20.0, 20.0))

    volumes = cpg.volumes().reshape(grid.nz, -1)
    # hər layda nazik/qalın sütun nisbəti 1:3
    assert np.allclose(volumes[:, 1] / volumes[:, 0], 3.0, rtol=1e-12)
    # `proportional` bölgü: hər lay sütunun yarısıdır
    assert np.allclose(volumes[0], volumes[1], rtol=1e-12)
    assert np.isclose(cpg.volumes().sum(), 20.0 * 20.0 * (10.0 + 30.0), rtol=1e-12)


def test_column_layer_edges_follow_the_interpolated_surfaces():
    """Sütun-üzrə lay sərhədləri məhz TOP/BOTTOM səthlərindən gəlir."""
    grid = CartesianGrid(2, 1, 4)
    top = np.array([1500.0, 1520.0])
    bottom = np.array([1540.0, 1620.0])
    cpg = CornerPointGeometry.from_nodes(
        grid, structural_nodes(grid, top, bottom, 10.0, 10.0))

    assert np.allclose(cpg.column_layer_edges(0, 0),
                       np.linspace(1500.0, 1540.0, 5), rtol=1e-12)
    assert np.allclose(cpg.column_layer_edges(1, 0),
                       np.linspace(1520.0, 1620.0, 5), rtol=1e-12)


def test_structural_nodes_rejects_unknown_layering():
    grid = CartesianGrid(2, 2, 1)
    with pytest.raises(ValueError, match="layering"):
        structural_nodes(grid, np.zeros(4), np.ones(4), 10.0, 10.0,
                         layering="equal")


def test_structural_nodes_accepts_full_volume_arrays():
    """`geology_service` xassələri `(ncell,)` saxlayır — hər iki forma
    qəbul edilir, üçüncü forma AÇIQ xəta verir."""
    grid = CartesianGrid(2, 2, 3)
    top = np.tile(np.full(4, 1000.0), grid.nz)
    bottom = np.tile(np.full(4, 1030.0), grid.nz)
    nodes = structural_nodes(grid, top, bottom, 10.0, 10.0)
    assert nodes.shape == (grid.ncell, 8, 3)
    with pytest.raises(ValueError, match="gözlənilən"):
        structural_nodes(grid, np.zeros(7), np.ones(7), 10.0, 10.0)


# ────────────────────────────────────────────────── validate_structure
def test_validate_structure_reports_crossing_surfaces_with_indices():
    grid = CartesianGrid(3, 2, 1)
    top = np.full(6, 2000.0)
    bottom = np.full(6, 2020.0)
    bottom[4] = 1990.0                    # daban tavandan YUXARI
    issues = validate_structure(top, bottom, min_thickness=0.1, grid=grid)
    assert len(issues) == 1
    assert "KƏSİŞİR" in issues[0]
    assert "1 sütunda" in issues[0]
    assert "(1, 1)" in issues[0]          # flat 4 → i=1, j=1


def test_validate_structure_reports_nan_and_thin_columns():
    grid = CartesianGrid(4, 1, 1)
    top = np.array([2000.0, 2000.0, np.nan, 2000.0])
    bottom = np.array([2020.0, 2000.05, 2020.0, 2020.0])
    issues = validate_structure(top, bottom, min_thickness=0.1, grid=grid)
    text = "\n".join(issues)
    assert "NaN" in text
    assert "qalınlıq minimumdan" in text
    assert "1 sütunda qalınlıq" in text


def test_validate_structure_is_silent_for_healthy_surfaces():
    grid = CartesianGrid(3, 3, 2)
    top = np.full(9, 1000.0)
    bottom = top + np.linspace(5.0, 50.0, 9)
    assert validate_structure(top, bottom, min_thickness=0.1, grid=grid) == []


def test_structure_statistics_reports_thickness_and_relief():
    grid = CartesianGrid(2, 2, 1)
    top = np.array([1000.0, 1010.0, 1020.0, 1030.0])
    bottom = top + np.array([10.0, 20.0, 30.0, 40.0])
    stats = structure_statistics(top, bottom, grid)
    assert stats["min_thickness"] == pytest.approx(10.0)
    assert stats["max_thickness"] == pytest.approx(40.0)
    assert stats["mean_thickness"] == pytest.approx(25.0)
    assert stats["relief"] == pytest.approx(30.0)
    assert stats["nan_columns"] == 0


# ═════════════════════════════════ §1b sütun-həssas dərinlik→lay uyğunluğu
def test_depth_to_k_is_column_sensitive_under_variable_thickness():
    """EYNİ metr dərinliyi qalın və nazik sütunda FƏRQLİ K verir.

    Perforasiyanın metrdən K-təbəqəsinə çevrilməsi məhz budur —
    düz müstəvi həndəsədə hər iki sütun eyni K verərdi (səhv)."""
    grid = CartesianGrid(2, 1, 4)
    top = np.array([2000.0, 2000.0])
    bottom = np.array([2040.0, 2008.0])          # 40 m və 8 m sütun
    cpg = CornerPointGeometry.from_nodes(
        grid, structural_nodes(grid, top, bottom, 100.0, 100.0))

    depth = 2007.0
    k_thick = depth_to_k(50.0, 50.0, depth, cpg)     # sütun (0,0), 10 m/lay
    k_thin = depth_to_k(150.0, 50.0, depth, cpg)     # sütun (1,0), 2 m/lay
    assert k_thick == 0
    assert k_thin == 3
    assert k_thick != k_thin


# ═══════════════════════════════ §2 servis yolu (`structure_from_wells`)
_WELL_XY = ((60.0, 60.0), (390.0, 60.0), (60.0, 390.0), (390.0, 390.0),
            (225.0, 225.0))


def _structural_dataset(tops, bottoms) -> WellDataset:
    """Beş quyu — TOP/BOTTOM + minimal petrofizika.

    `geology_adapter._broadcast_dataset` ilə eyni forma (hər quyu üçün
    BİR nümunə, lay etiketi yoxdur) — yəni UI-nin defolt yolu ilə eyni.
    """
    poro = (0.20, 0.24, 0.22, 0.18, 0.21)
    permx = (150.0, 400.0, 220.0, 90.0, 260.0)
    samples = []
    for index, ((x, y), top, bottom) in enumerate(zip(_WELL_XY, tops, bottoms)):
        samples.append(WellSample(
            well=f"W-{index + 1}", x=x, y=y,
            values={"TOP": float(top), "BOTTOM": float(bottom),
                    "PORO": poro[index], "PERMX": permx[index]}))
    return WellDataset(samples=samples, source="test")


def _spec(**kwargs) -> GeologicalGridSpec:
    base = dict(nx=9, ny=9, nz=3, dx=50.0, dy=50.0, dz=10.0,
                top_depth=2000.0, dip_x=0.0, dip_y=0.0)
    base.update(kwargs)
    return GeologicalGridSpec(**base)


def _build(spec, dataset):
    return WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(dataset, spec)


def test_default_path_is_untouched_when_structure_is_off():
    """`structure_from_wells=False` → HƏNDƏSƏ TİPİ və HƏCMLƏR əvvəlki kimi."""
    dataset = _structural_dataset((2000, 2010, 2020, 2030, 2015),
                                  (2040, 2055, 2070, 2085, 2062))
    model, report = _build(_spec(), dataset)

    assert type(model.geometry) is CellGeometry
    assert np.allclose(model.geometry.volumes(), 50.0 * 50.0 * 10.0)
    assert model.geometry.top_depth == 2000.0
    assert model.structural_issues == []
    assert not report.has_blocking
    # TOP/BOTTOM hələ də XASSƏ kimi interpolyasiya olunur (köhnə davranış)
    assert "TOP" in model.property_maps and "BOTTOM" in model.property_maps


def test_structure_from_wells_builds_corner_point_geometry_from_surfaces():
    dataset = _structural_dataset((2000, 2010, 2020, 2030, 2015),
                                  (2040, 2055, 2070, 2085, 2062))
    model, report = _build(_spec(structure_from_wells=True), dataset)

    assert isinstance(model.geometry, CornerPointGeometry)
    assert not report.has_blocking
    # Həndəsə MƏHZ interpolyasiya edilmiş səthdən gəlir: hər sütunun
    # tavan/daban sərhədi TOP/BOTTOM xəritəsinin həmin sütun dəyəridir.
    areal = model.grid.nx * model.grid.ny
    top = model.property_maps["TOP"].values[:areal]
    bottom = model.property_maps["BOTTOM"].values[:areal]
    for i, j in ((0, 0), (4, 4), (8, 8), (0, 8)):
        edges = model.geometry.column_layer_edges(i, j)
        flat = j * model.grid.nx + i
        assert edges[0] == pytest.approx(top[flat], rel=1e-12)
        assert edges[-1] == pytest.approx(bottom[flat], rel=1e-12)
    # Struktur HƏQİQƏTƏN dəyişkəndir (test özü degenerativ deyil)
    assert (bottom - top).max() - (bottom - top).min() > 1.0


def test_structural_mode_says_out_loud_which_parameters_are_unused():
    """§3.5 — `top_depth`/`dip_x`/`dip_y`/`dz` artıq iştirak etmir və bu,
    hesabatda AÇIQ yazılır (səssiz "parametr işləmir" halı QADAĞANDIR)."""
    dataset = _structural_dataset((2000, 2010, 2020, 2030, 2015),
                                  (2040, 2055, 2070, 2085, 2062))
    _model, report = _build(_spec(structure_from_wells=True), dataset)
    text = report.as_text()
    assert "İSTİFADƏ OLUNMUR" in text
    assert "tavan dərinliyi" in text and "maillik X" in text and "DZ" in text
    assert "Struktur statistikası" in text
    assert "Həndəsə keyfiyyəti" in text


def test_structural_geometry_ignores_top_depth_and_dip_completely():
    """`top_depth`/`dip_x`/`dip_y`/`dz` DƏYİŞSƏ DƏ həndəsə eyni qalır —
    "iştirak etmir" iddiasının SÜBUTU."""
    dataset = _structural_dataset((2000, 2010, 2020, 2030, 2015),
                                  (2040, 2055, 2070, 2085, 2062))
    first, _ = _build(_spec(structure_from_wells=True), dataset)
    second, _ = _build(
        _spec(structure_from_wells=True, top_depth=500.0, dip_x=7.0, dip_y=-3.0,
              dz=99.0), dataset)
    assert np.allclose(first.geometry.nodes, second.geometry.nodes, rtol=1e-12)
    assert np.allclose(first.geometry.volumes(), second.geometry.volumes(), rtol=1e-12)


def test_structural_surfaces_are_interpolated_once_not_per_layer():
    """§3.2 — TOP/BOTTOM BİR DƏFƏ (areal) hesablanır; hər K eyni səthi
    daşıyır, yəni həndəsədən sonra TƏKRAR interpolyasiya YOXDUR."""
    dataset = _structural_dataset((2000, 2010, 2020, 2030, 2015),
                                  (2040, 2055, 2070, 2085, 2062))
    model, _ = _build(_spec(structure_from_wells=True), dataset)
    grid = model.grid
    for key in ("TOP", "BOTTOM"):
        volume = model.property_maps[key].values.reshape(grid.shape)
        for k in range(1, grid.nz):
            assert np.array_equal(volume[0], volume[k])


def test_other_properties_are_interpolated_after_the_real_geometry():
    """§3.2 sırasının SÜBUTU: PORO/PERMX 3D kriginqi hüceyrə mərkəzinin
    Z-ini işlədir, ona görə struktur rejimində nəticə düz-müstəvi
    rejimindən FƏRQLƏNMƏLİDİR (eyni quyu məlumatı ilə)."""
    dataset = _structural_dataset((2000, 2040, 2080, 2120, 2060),
                                  (2030, 2075, 2120, 2165, 2095))
    flat, _ = _build(_spec(), dataset)
    structural, _ = _build(_spec(structure_from_wells=True), dataset)
    assert not np.allclose(structural.property_maps["PORO"].values,
                           flat.property_maps["PORO"].values)


def test_thickness_source_constant_takes_top_from_wells_and_dz_for_thickness():
    dataset = _structural_dataset((2000, 2010, 2020, 2030, 2015),
                                  (2040, 2055, 2070, 2085, 2062))
    model, report = _build(
        _spec(structure_from_wells=True, thickness_source="constant", nz=3, dz=7.0),
        dataset)

    assert np.allclose(model.geometry.dz_per_cell(), 7.0, rtol=1e-12)
    # tavan hələ də QUYULARDAN (əyri səth) gəlir
    tops = model.geometry.nodes[:model.grid.nx * model.grid.ny, 0:4, 2]
    assert tops.max() - tops.min() > 5.0
    assert "thickness_source='constant'" in report.as_text()


def test_thickness_source_constant_accepts_per_layer_dz():
    dataset = _structural_dataset((2000, 2010, 2020, 2030, 2015),
                                  (2040, 2055, 2070, 2085, 2062))
    model, _ = _build(
        _spec(structure_from_wells=True, thickness_source="constant", nz=3,
              dz=[4.0, 6.0, 8.0]), dataset)
    # `proportional` bölgü Σdz = 18 m-i BƏRABƏR üç hissəyə bölür
    assert np.allclose(model.geometry.dz_per_cell(), 6.0, rtol=1e-12)


def test_missing_bottom_column_is_an_explicit_error_not_a_silent_fallback():
    samples = [WellSample(well=f"W-{n + 1}", x=x, y=y,
                          values={"TOP": 2000.0 + 10.0 * n, "PORO": 0.2,
                                  "PERMX": 150.0})
               for n, (x, y) in enumerate(_WELL_XY)]
    dataset = WellDataset(samples=samples, source="test")
    with pytest.raises(ValueError, match="lay altı"):
        _build(_spec(structure_from_wells=True), dataset)
    # `constant` rejimində eyni dataset İŞLƏYİR — tavan kifayətdir
    model, _ = _build(
        _spec(structure_from_wells=True, thickness_source="constant"), dataset)
    assert isinstance(model.geometry, CornerPointGeometry)


def test_deactivate_policy_refuses_explicitly():
    dataset = _structural_dataset((2000, 2010, 2020, 2030, 2015),
                                  (2040, 2055, 2070, 2085, 2062))
    with pytest.raises(NotImplementedError) as error:
        _build(_spec(structure_from_wells=True, on_zero_thickness="deactivate"),
               dataset)
    assert "ACTNUM" in str(error.value)


def test_unknown_policies_are_rejected():
    dataset = _structural_dataset((2000, 2010, 2020, 2030, 2015),
                                  (2040, 2055, 2070, 2085, 2062))
    for kwargs, needle in ((dict(thickness_source="guess"), "thickness_source"),
                           (dict(layering="equal"), "layering"),
                           (dict(on_zero_thickness="ignore"), "on_zero_thickness")):
        with pytest.raises(ValueError, match=needle):
            _build(_spec(structure_from_wells=True, **kwargs), dataset)


# ── xəta siyasətləri ─────────────────────────────────────────────────
def _crossing_dataset() -> WellDataset:
    """Daban tavandan YUXARI olan quyular — səthlər kəsişir."""
    return _structural_dataset((2000, 2010, 2020, 2030, 2015),
                               (1990, 1995, 2005, 2010, 1998))


def test_error_policy_blocks_and_reservoir_builder_refuses_the_model():
    model, report = _build(_spec(structure_from_wells=True), _crossing_dataset())
    assert report.has_blocking
    assert any("KƏSİŞİR" in message for message in report.blocking)
    assert model.structural_issues                     # modeldə də qalır
    assert model.validate()                            # simulyasiya qapısı bağlıdır
    with pytest.raises(ValueError, match="Geoloji model natamamdır"):
        ReservoirModelBuilder().build(
            geological_model=model, wells=five_spot(model.grid),
            scal=default_scal(), name="rədd edilməli")


def test_thin_columns_block_even_when_volume_is_positive():
    """Qalınlıq `min_thickness`-dən azdır, amma MÜSBƏTDİR — köhnə həcm
    yoxlaması bunu tutmazdı, ona görə AYRICA bloklanır."""
    dataset = _structural_dataset((2000, 2000, 2000, 2000, 2000),
                                  (2000.02, 2000.03, 2000.02, 2000.03, 2000.02))
    model, report = _build(
        _spec(structure_from_wells=True, min_thickness=0.5), dataset)
    assert report.has_blocking
    assert any("qalınlıq minimumdan" in message for message in report.blocking)
    assert np.all(model.geometry.volumes() > 0.0)      # həcm MÜSBƏTDİR
    with pytest.raises(ValueError, match="Geoloji model natamamdır"):
        ReservoirModelBuilder().build(
            geological_model=model, wells=five_spot(model.grid),
            scal=default_scal(), name="rədd edilməli")


def test_clamp_policy_reports_the_number_of_repaired_columns():
    dataset = _structural_dataset((2000, 2010, 2020, 2030, 2015),
                                  (2000.01, 2010.01, 2020.01, 2030.01, 2015.01))
    model, report = _build(
        _spec(structure_from_wells=True, on_zero_thickness="clamp",
              min_thickness=0.5), dataset)

    assert not report.has_blocking
    clamp_lines = [w for w in report.warnings if "on_zero_thickness='clamp'" in w]
    assert len(clamp_lines) == 1
    # RƏQƏM görünür: neçə sütun düzəldilib / neçədən
    columns = model.grid.nx * model.grid.ny
    assert f"{columns} sütunda lay altı" in clamp_lines[0]
    assert f"ümumi {columns} sütundan" in clamp_lines[0]
    assert np.allclose(model.geometry.dz_per_cell(), 0.5 / model.grid.nz, rtol=1e-9)


def test_clamp_does_not_hide_nan_surfaces():
    """`clamp` NaN-ı DÜZƏLTMİR — NaN yox olan məlumatdır, nazik lay deyil."""
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    spec = _spec(structure_from_wells=True, on_zero_thickness="clamp")
    grid = CartesianGrid(spec.nx, spec.ny, spec.nz)
    report = InterpolationReport()
    surfaces = {"TOP": np.full(grid.nx * grid.ny, 2000.0),
                "BOTTOM": np.full(grid.nx * grid.ny, 2030.0)}
    surfaces["BOTTOM"][5] = np.nan
    _geometry, blocking = builder._build_structural_geometry(
        grid, spec, surfaces, report)
    assert blocking
    assert any("NaN" in message for message in blocking)


def test_single_well_surface_is_reported_as_constant():
    """§4.2 — 3 quyudan az olanda səthin mənasızlığı AÇIQ yazılır."""
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    samples = [WellSample(well="W-1", x=60.0, y=60.0,
                          values={"TOP": 2000.0, "BOTTOM": 2030.0}),
               WellSample(well="W-1", x=60.0, y=60.0, layer=1,
                          values={"TOP": 2000.0, "BOTTOM": 2030.0})]
    dataset = WellDataset(samples=samples, source="test")
    grid = CartesianGrid(4, 4, 2)
    geometry = CellGeometry(grid, 50.0, 50.0, 10.0, top_depth=2000.0)
    model = GeologicalModel(name="tək quyu", grid=grid, geometry=geometry,
                            regions=RegionSet.single(grid.ncell))
    report = InterpolationReport()
    builder._interpolate_areal(dataset, ["TOP"],
                               builder._cell_centres(grid, geometry), grid,
                               geometry, model, report)
    assert any("yalnız 1 quyuda" in message for message in report.warnings)
    assert any("SABİTDİR" in message for message in report.warnings)


# ═══════════════════════════════════ §3 FİZİKAYA TƏSİRİN SÜBUTU
def _structural_geometry(top, bottom, nx, ny, nz, dx=50.0, dy=50.0):
    grid = CartesianGrid(nx, ny, nz)
    return CornerPointGeometry.from_nodes(
        grid, structural_nodes(grid, top, bottom, dx, dy))


def _geological_model(geometry, porosity=0.22, permx=150.0) -> GeologicalModel:
    grid = geometry.grid
    model = GeologicalModel(name="struktur", grid=grid, geometry=geometry,
                            regions=RegionSet.single(grid.ncell))
    model.add_property(PropertyMap.from_array(
        "PORO", np.full(grid.ncell, porosity), grid.ncell))
    model.add_property(PropertyMap.from_array(
        "PERMX", np.full(grid.ncell, permx), grid.ncell))
    return model


def _reservoir_from_geometry(geometry, initial=None, wells=None):
    """Verilmiş həndəsə ilə minimal rezervuar modeli."""
    grid = geometry.grid
    scal = default_scal()
    if wells is None:
        wells = [
            Well("INJ", WellType.INJECTOR, WellControl(ControlMode.RATE, 40.0),
                 [Perforation(0, 0, k) for k in range(grid.nz)]),
            Well("PROD", WellType.PRODUCER, WellControl(ControlMode.BHP, 200.0),
                 [Perforation(grid.nx - 1, grid.ny - 1, k) for k in range(grid.nz)]),
        ]
    return ReservoirModelBuilder().build(
        geological_model=_geological_model(geometry), wells=wells, scal=scal,
        initial=initial or InitialConditions(datum_pressure=250.0,
                                             water_saturation=scal.swc),
        name="struktur modeli")


def test_pore_volume_ratio_follows_the_thickness_ratio():
    """§5.1 — dəyişən qalınlıqda sütunlar arası PV nisbəti qalınlıq
    nisbətinə BƏRABƏRDİR (məsamə həcmi struktura tabe oldu)."""
    nx, ny, nz = 4, 1, 2
    top = np.full(nx * ny, 2000.0)
    bottom = top + np.array([10.0, 20.0, 30.0, 40.0])
    reservoir = _reservoir_from_geometry(_structural_geometry(top, bottom, nx, ny, nz))

    pv = reservoir.pore_volume().reshape(nz, ny, nx)
    column_pv = pv.sum(axis=0).ravel()
    assert np.allclose(column_pv / column_pv[0], [1.0, 2.0, 3.0, 4.0], rtol=1e-12)


def test_flat_structure_reproduces_the_cartesian_pore_volume():
    """Düz struktur → köhnə Kartezian PV ilə maşın dəqiqliyində eyni."""
    nx, ny, nz = 5, 4, 2
    top = np.full(nx * ny, 1800.0)
    bottom = np.full(nx * ny, 1800.0 + nz * 10.0)
    structural = _reservoir_from_geometry(_structural_geometry(top, bottom, nx, ny, nz))
    grid = CartesianGrid(nx, ny, nz)
    cartesian = _reservoir_from_geometry(
        CellGeometry(grid, 50.0, 50.0, 10.0, top_depth=1800.0))
    assert np.allclose(structural.pore_volume(), cartesian.pore_volume(), rtol=1e-12)


def test_transmissibility_matches_cartesian_for_a_flat_structure():
    """§5.5 — düz struktur halında CPG transmissibillikləri Kartezian
    ilə maşın dəqiqliyində üst-üstə düşür."""
    nx, ny, nz = 5, 4, 3
    top = np.full(nx * ny, 1800.0)
    bottom = np.full(nx * ny, 1800.0 + nz * 10.0)
    structural = _reservoir_from_geometry(_structural_geometry(top, bottom, nx, ny, nz))
    grid = CartesianGrid(nx, ny, nz)
    cartesian = _reservoir_from_geometry(
        CellGeometry(grid, 50.0, 50.0, 10.0, top_depth=1800.0))

    scheme = TwoPointFluxDiscretization()
    assert np.allclose(scheme.build(structural).transmissibility,
                       scheme.build(cartesian).transmissibility, rtol=1e-12)


def test_dipping_structure_creates_lateral_gravity_and_runs():
    """§5.2 — cazibə `cell_depths()` fərqindən qurulur (`_setup_gravity`).

    ŞAQULİ əlaqədə dərinlik fərqi struktursuz modeldə də var (laylar
    üst-üstədir), ona görə `_has_gravity` tək başına heç nə sübut etmir.
    Struktur məhz YANAL (X üzrə) əlaqələrdə cazibə yaradır — düz layda
    o fərq SIFIRDIR, maili layda deyil. Test bunu ölçür."""
    nx, ny, nz = 6, 1, 2
    top = 2000.0 + 8.0 * np.arange(nx * ny, dtype=float)     # maili tavan
    dipping = _reservoir_from_geometry(
        _structural_geometry(top, top + 20.0, nx, ny, nz))
    flat_top = np.full(nx * ny, 2000.0)
    flat = _reservoir_from_geometry(
        _structural_geometry(flat_top, flat_top + 20.0, nx, ny, nz))

    service = make_service()
    config = short_config(end_time=20.0, snapshots=2)
    engine_dipping = service.create_engine(dipping, config)
    engine_flat = service.create_engine(flat, config)
    assert engine_dipping._has_gravity is True

    lateral = engine_dipping._connections.axis == 0
    assert np.all(np.abs(engine_dipping._depth_difference[lateral]) > 1.0)
    lateral_flat = engine_flat._connections.axis == 0
    assert np.allclose(engine_flat._depth_difference[lateral_flat], 0.0, atol=1e-12)

    result = service.run(dipping, config)
    assert result.converged
    assert result.steps > 0


def test_equilibration_oil_column_follows_the_structural_relief():
    """§5.4 — kontakt SABİT dərinlikdədirsə, struktur qalxıqlığında neft
    zonası QALIN, çökəkliyində NAZİK olur."""
    nx, ny, nz = 2, 1, 6
    top = np.array([2000.0, 2060.0])          # sol sütun QALXIQ, sağ ÇÖKƏK
    geometry = _structural_geometry(top, top + 60.0, nx, ny, nz)
    scal = default_scal()
    reservoir = _reservoir_from_geometry(
        geometry,
        initial=InitialConditions(datum_pressure=250.0, datum_depth=2000.0,
                                  use_equilibration=True,
                                  oil_water_contact=2075.0),
        wells=[Well("PROD", WellType.PRODUCER, WellControl(ControlMode.BHP, 200.0),
                    [Perforation(0, 0, 0)])])

    state = EquilibriumInitializationProvider().initialize(reservoir)
    sw = np.asarray(state.water_saturation, float).reshape(nz, ny, nx)
    oil_threshold = 0.5 * (scal.swc + (1.0 - scal.sor))
    oil_cells_high = int(np.sum(sw[:, 0, 0] < oil_threshold))
    oil_cells_low = int(np.sum(sw[:, 0, 1] < oil_threshold))
    assert oil_cells_high > oil_cells_low


def test_perforation_interval_maps_to_different_layers_per_column():
    """§5.3 — eyni metr intervalı qalın və nazik sütunda FƏRQLİ
    K-təbəqələri verir."""
    nx, ny, nz = 2, 1, 4
    top = np.full(nx * ny, 2000.0)
    bottom = top + np.array([40.0, 8.0])
    geometry = _structural_geometry(top, bottom, nx, ny, nz)

    assert interval_layers(25.0, 25.0, 2000.5, 2007.5, geometry) == [0]
    assert interval_layers(75.0, 25.0, 2000.5, 2007.5, geometry) == [0, 1, 2, 3]


# ═══════════════════════════════════════════ §4 serializasiya round-trip
def test_structural_geometry_survives_the_imx_round_trip(tmp_path):
    """`nodes` fayla yazılır və geri oxunur — həcmlər MAŞIN
    DƏQİQLİYİNDƏ eynidir. Fayl ölçüsü də qeyd olunur ki, `nodes`
    (hüceyrə başına 24 float) gələcəkdə böyüsə görünsün."""
    dataset = _structural_dataset((2000, 2010, 2020, 2030, 2015),
                                  (2040, 2055, 2070, 2085, 2062))
    model, _ = _build(_spec(structure_from_wells=True), dataset)
    assert isinstance(model.geometry, CornerPointGeometry)

    project = Project(name="struktur layihəsi")
    project.add_geological_model(model)
    path = tmp_path / "struktur.imx"
    serializer = ProjectSerializer()
    serializer.save(project, str(path))
    restored = serializer.load(str(path)).geological_models[model.name]

    assert isinstance(restored.geometry, CornerPointGeometry)
    assert np.array_equal(restored.geometry.nodes, model.geometry.nodes)
    assert np.array_equal(restored.geometry.volumes(), model.geometry.volumes())
    assert np.array_equal(restored.geometry.cell_depths(),
                          model.geometry.cell_depths())

    # Ölçü qeydi: 9×9×3 = 243 hüceyrə → 243·24 = 5832 ədəd `nodes`.
    # Fayl gzip-lidir (bax `serialization.py`), ona görə praktikada
    # kiçikdir; hədd yalnız "gözlənilmədən böyüdü" halını tutmaq üçündür.
    size = path.stat().st_size
    assert size < 400_000, f"fayl ölçüsü gözləniləndən böyükdür: {size} bayt"


def test_blocking_structural_issues_are_not_lost_on_save_and_load(tmp_path):
    """Yararsız həndəsəli layihəni saxlayıb açmaq simulyasiya qapısını
    SÜKUTLA AÇMAMALIDIR."""
    model, report = _build(_spec(structure_from_wells=True), _crossing_dataset())
    assert report.has_blocking

    project = Project(name="yararsız")
    project.add_geological_model(model)
    path = tmp_path / "yararsiz.imx"
    serializer = ProjectSerializer()
    serializer.save(project, str(path))
    restored = serializer.load(str(path)).geological_models[model.name]
    assert restored.structural_issues == model.structural_issues
    assert restored.validate()


def test_old_projects_without_the_key_still_load(tmp_path):
    """Köhnə `.imx` faylında `structural_issues` açarı YOXDUR — geriyə
    uyğunluq: boş siyahı kimi oxunur, xəta atılmır."""
    dataset = _structural_dataset((2000, 2010, 2020, 2030, 2015),
                                  (2040, 2055, 2070, 2085, 2062))
    model, _ = _build(_spec(), dataset)
    serializer = ProjectSerializer()
    data = serializer.geological_model_to_dict(model)
    data.pop("structural_issues")
    restored = serializer.geological_model_from_dict(data)
    assert restored.structural_issues == []
    assert type(restored.geometry) is CellGeometry


# ═══════════════════════════════════════════════════════════ §5 UI qatı
# `test_ui_units.py` ilə eyni üslub: real (ekransız) Qt widget-ləri —
# burada YALNIZ struktur deyil, faktiki söndürmə/ötürmə davranışı
# yoxlanılır. `MainWindow` QURULMUR (VTK render pəncərəsi ekransız
# mühitdə işləmir), ona görə pəncərə qatından yalnız SAF köməkçilər
# sınanır.
QtWidgets = pytest.importorskip("PyQt5.QtWidgets")


@pytest.fixture(scope="module")
def qapp():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield app


def test_grid_panel_default_keeps_the_old_behaviour(qapp):
    from imex2d.ui.panels import GridGeometryPanel
    panel = GridGeometryPanel()
    assert panel.structure_values() == {"structure_from_wells": False,
                                        "thickness_source": "wells"}
    assert panel.top_depth.isEnabled() and panel.dip_x.isEnabled()
    assert panel.dz.isEnabled()
    assert panel.top_depth.toolTip() == ""
    # `values()` sintetik qurucuya BİRBAŞA açıldığı üçün yeni açar
    # ORAYA əlavə OLUNMAMALIDIR.
    assert set(panel.values()) == {"nx", "ny", "nz", "dx", "dy", "dz",
                                   "top_depth", "dip_x", "dip_y"}


def test_grid_panel_disables_the_parameters_that_stop_working(qapp):
    """§4.3 — struktur rejimində iştirak etməyən sahələr SÖNDÜRÜLÜR və
    tooltip səbəbi yazır."""
    from imex2d.ui.panels import GridGeometryPanel
    panel = GridGeometryPanel()
    panel.structure_from_wells.setChecked(True)

    for widget in (panel.top_depth, panel.dip_x, panel.dip_y, panel.dz,
                   panel.thickness_mode, panel.base_depth, panel.per_layer):
        assert not widget.isEnabled()
        assert "istifadə olunmur" in widget.toolTip()
    assert panel.thickness_source.isEnabled()
    assert panel.structure_values()["structure_from_wells"] is True


def test_grid_panel_keeps_dz_alive_for_the_constant_thickness_source(qapp):
    """`thickness_source="constant"` halında qalınlıq MƏHZ DZ-dən gəlir,
    ona görə DZ söndürülmür — yalnız tavan/maillik söndürülür."""
    from imex2d.ui.panels import GridGeometryPanel
    panel = GridGeometryPanel()
    panel.structure_from_wells.setChecked(True)
    panel.thickness_source.setCurrentIndex(
        panel.thickness_source.findData("constant"))

    assert panel.dz.isEnabled()
    assert panel.thickness_mode.isEnabled()
    assert panel.dz.toolTip() == ""
    assert not panel.top_depth.isEnabled()
    assert panel.structure_values()["thickness_source"] == "constant"


def test_grid_panel_restores_everything_when_structure_is_switched_off(qapp):
    from imex2d.ui.panels import GridGeometryPanel
    panel = GridGeometryPanel()
    panel.structure_from_wells.setChecked(True)
    panel.structure_from_wells.setChecked(False)

    assert panel.top_depth.isEnabled() and panel.dip_x.isEnabled()
    assert panel.dz.isEnabled() and panel.thickness_mode.isEnabled()
    assert not panel.thickness_source.isEnabled()
    for widget in (panel.top_depth, panel.dip_x, panel.dz):
        assert widget.toolTip() == ""


def test_grid_panel_values_feed_the_spec_without_manual_mapping(qapp):
    """`grid_panel.values()` + `structure_values()` birlikdə DƏQİQ
    `GeologicalGridSpec` konstruktoruna uyğun gəlir (main_window bunu
    `**` ilə açır — açar adı fərqlənsəydi TypeError olardı)."""
    from imex2d.ui.panels import GridGeometryPanel
    panel = GridGeometryPanel()
    panel.structure_from_wells.setChecked(True)
    values = panel.values()
    spec = GeologicalGridSpec(
        nx=values["nx"], ny=values["ny"], nz=values["nz"],
        dx=values["dx"], dy=values["dy"], dz=values["dz"],
        top_depth=values["top_depth"], dip_x=values["dip_x"],
        dip_y=values["dip_y"], **panel.structure_values())
    assert spec.structure_from_wells is True
    assert spec.thickness_source == "wells"


def test_main_window_reports_only_the_wells_that_actually_moved():
    """§4.4 tələsi — həndəsə dəyişəndə quyu i/j/K yenidən hesablanır və
    DƏYİŞƏNLƏR istifadəçiyə göstərilir (dəyişməyənlər səs-küy yaratmır)."""
    from imex2d.ui.main_window import MainWindow
    before = {"W-1": (2, 3, 0, 1), "W-2": (5, 5, 0, 0), "W-3": (1, 1, 2, 2)}
    after = {"W-1": (2, 3, 0, 3), "W-2": (5, 5, 0, 0), "W-3": (4, 1, 2, 2)}
    messages = MainWindow._describe_moved_wells(before, after)

    assert len(messages) == 2
    assert any(m.startswith("• W-1:") and "K 1–2" in m and "K 1–4" in m
               for m in messages)
    assert any(m.startswith("• W-3:") and "(1, 1)" in m and "(4, 1)" in m
               for m in messages)
    assert not any("W-2" in m for m in messages)
    assert MainWindow._describe_moved_wells(before, before) == []


def test_main_window_passes_the_structure_options_into_the_spec():
    """`_interpolate_geology` `structure_values()`-i spec-ə ötürür və
    `NotImplementedError`-i də tutur (AST yoxlaması — `MainWindow`
    ekransız mühitdə qurula bilmir, bax modul şərhi)."""
    import ast
    import inspect
    from imex2d.ui import main_window as module

    source = inspect.getsource(module)
    tree = ast.parse(source)
    function = next(node for node in ast.walk(tree)
                    if isinstance(node, ast.FunctionDef)
                    and node.name == "_interpolate_geology")
    body = ast.dump(function)
    assert "structure_values" in body
    assert "NotImplementedError" in body
    assert "set_geology_context" in body
    assert "_describe_moved_wells" in body


# ═════════════════════════════ §6 UI adapteri ilə uçdan-uca (§4.6, §7)
def _ui_rows():
    """GUI-də əl ilə doldurulan cədvəlin eynisi — fərqli lay üstü/altı
    dərinlikləri ilə 4 quyu (bax qəbul meyarları, §7)."""
    from imex2d.domain.geology import GeologicalWell
    return [
        GeologicalWell(name="W-1", in_model=True, x=60.0, y=60.0, porosity=0.20,
                       permeability=150.0, water_saturation=0.25,
                       top=2000.0, bottom=2030.0),
        GeologicalWell(name="W-2", in_model=True, x=380.0, y=60.0, porosity=0.24,
                       permeability=400.0, water_saturation=0.22,
                       top=2020.0, bottom=2075.0),
        GeologicalWell(name="W-3", in_model=True, x=60.0, y=380.0, porosity=0.22,
                       permeability=220.0, water_saturation=0.27,
                       top=2040.0, bottom=2080.0),
        GeologicalWell(name="W-4", in_model=True, x=380.0, y=380.0, porosity=0.18,
                       permeability=90.0, water_saturation=0.30,
                       top=2060.0, bottom=2130.0),
    ]


def test_ui_adapter_carries_top_and_bottom_through_to_the_geometry():
    """§4.6 — `wells_to_dataset` TOP/BOTTOM-u İTİRMİR və struktur yolu
    onları həndəsəyə çevirir (adapterdə dəyişiklik LAZIM DEYİL)."""
    from imex2d.application.geology_adapter import AREAL_TARGETS, wells_to_dataset

    dataset, skipped = wells_to_dataset(_ui_rows(), "Kriging")
    assert not skipped
    for target in AREAL_TARGETS:
        assert target in dataset.property_names()

    spec = GeologicalGridSpec(nx=11, ny=11, nz=3, dx=40.0, dy=40.0, dz=10.0,
                              top_depth=2000.0, structure_from_wells=True)
    model, report = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(
        dataset, spec)
    assert isinstance(model.geometry, CornerPointGeometry)
    assert not report.has_blocking


def test_structural_mode_changes_ooip_and_switching_it_off_restores_it():
    """§7 qəbul meyarı: struktur rejimi OOIP-i DƏYİŞİR, söndürüləndə isə
    KÖHNƏ rəqəmlərə qayıdır (çünki köhnə yol toxunulmayıb)."""
    from imex2d.application.geology_adapter import wells_to_dataset

    dataset, _ = wells_to_dataset(_ui_rows(), "Kriging")

    def _reservoir(structural: bool):
        spec = GeologicalGridSpec(nx=11, ny=11, nz=3, dx=40.0, dy=40.0, dz=10.0,
                                  top_depth=2000.0, structure_from_wells=structural)
        geology, _report = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(
            dataset, spec)
        return geology, ReservoirModelBuilder().build(
            geological_model=geology, wells=five_spot(geology.grid),
            scal=default_scal(), name="qəbul meyarı")

    flat_geology, flat = _reservoir(False)
    structural_geology, structural = _reservoir(True)

    assert type(flat_geology.geometry) is CellGeometry
    assert isinstance(structural_geology.geometry, CornerPointGeometry)
    # Düz modeldə qalınlıq HƏR YERDƏ eynidir, strukturda deyil
    assert np.allclose(flat_geology.geometry.dz_per_cell(), 10.0)
    thickness = structural_geology.geometry.dz_per_cell()
    assert thickness.max() - thickness.min() > 5.0
    # Məsamə həcmi (deməli OOIP və RF) HƏQİQƏTƏN dəyişir
    assert structural.pore_volume().sum() > flat.pore_volume().sum() * 1.2

    # Söndürüləndə köhnə rəqəm BİT-BƏ-BİT qayıdır
    repeat_geology, repeat = _reservoir(False)
    assert np.array_equal(repeat.pore_volume(), flat.pore_volume())
    assert repeat_geology.geometry.top_depth == 2000.0
