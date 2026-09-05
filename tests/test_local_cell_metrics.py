"""Phase 5E — "Variable Grid / Local Cell Metrics".

Tapşırıq §1–§9: hər hüceyrə ÖZ real corner-point həndəsəsindən istifadə
edir; heç bir istehsal hesablamasında gizli global/nominal `DX/DY/DZ`
fallback qalmır.

    Corner-Point Geometry → Local Cell Metrics → Real Cell Volume
        → Real Geometric Distances → Downstream Calculations

Test modelləri (tapşırıq §8): uniform, dəyişkən qalınlıq, maili (dipping),
əyri (warped), pinch-out, ACTNUM + dəyişkən həndəsə, qeyri-ortoqonal.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.domain.corner_point_geometry import CornerPointGeometry
from imex2d.domain.geometry import (CellGeometry, depth_to_k, interval_layers,
                                    layer_edges, xy_to_ij)
from imex2d.domain.grid import CartesianGrid

NX, NY, NZ = 4, 3, 2
DX, DY, DZ, TOP = 100.0, 80.0, 10.0, 2000.0


# ══════════════════════════════════ model qurucuları (§8) ═══════════════
def _coord(nx=NX, ny=NY, dx=DX, dy=DY, top=TOP, base=3200.0,
           tilt_x=0.0, tilt_wave=0.0) -> np.ndarray:
    values = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            x, y = i * dx, j * dy
            shift_x = tilt_x + tilt_wave * np.sin(1.3 * i + 0.7 * j)
            shift_y = tilt_wave * np.cos(0.9 * i - 1.1 * j)
            values += [x, y, top, x + shift_x, y + shift_y, base]
    return np.array(values, float)


def _zcorn(nx=NX, ny=NY, nz=NZ, top=TOP, dip_x=0.0, dip_y=0.0,
           layer_dz=None, warp=0.0, pinch_layer=None) -> np.ndarray:
    """`layer_dz` — lay başına qalınlıq (dəyişkən qalınlıq modeli);
    `warp` — künc-indeksindən asılı hamar səth (konformal əyrilik);
    `pinch_layer` — həmin layın qalınlığı SIFIR (pinch-out)."""
    layer_dz = list(layer_dz) if layer_dz is not None else [DZ] * nz
    z = np.zeros((nz, 2, ny, 2, nx, 2))
    for j in range(ny):
        for cj in range(2):
            for i in range(nx):
                for ci in range(2):
                    gi, gj = i + ci, j + cj
                    surface = (top + dip_x * gi + dip_y * gj
                               + warp * np.sin(1.7 * gi) * np.cos(1.1 * gj))
                    depth = surface
                    for k in range(nz):
                        thickness = 0.0 if k == pinch_layer else layer_dz[k]
                        z[k, 0, j, cj, i, ci] = depth
                        z[k, 1, j, cj, i, ci] = depth + thickness
                        depth += thickness
    return z.ravel()


def _grid(nx=NX, ny=NY, nz=NZ, actnum=None) -> CartesianGrid:
    return CartesianGrid(nx, ny, nz, actnum)


def _cpg(grid=None, **kwargs) -> CornerPointGeometry:
    """Ölçülər HƏMİŞƏ `grid`-dən götürülür (`nx/ny/nz` kwargs kimi
    ötürülmür) — qalan açarlar COORD və ZCORN arasında bölüşdürülür."""
    grid = grid if grid is not None else _grid()
    coord_keys = ("tilt_x", "tilt_wave")
    coord = _coord(nx=grid.nx, ny=grid.ny,
                   **{k: v for k, v in kwargs.items() if k in coord_keys})
    zcorn = _zcorn(nx=grid.nx, ny=grid.ny, nz=grid.nz,
                   **{k: v for k, v in kwargs.items()
                      if k not in coord_keys and k not in ("nx", "ny", "nz")})
    geometry, _ = CornerPointGeometry.from_grdecl(grid, coord, zcorn)
    return geometry


def _uniform() -> CornerPointGeometry:
    return _cpg()


def _variable_thickness() -> CornerPointGeometry:
    return _cpg(layer_dz=[6.0, 21.0])


def _dipping() -> CornerPointGeometry:
    return _cpg(dip_x=25.0, dip_y=9.0)


def _warped() -> CornerPointGeometry:
    return _cpg(tilt_x=300.0, tilt_wave=120.0, warp=7.0)


def _non_orthogonal() -> CornerPointGeometry:
    return _cpg(tilt_x=450.0, dip_x=20.0)


def _pinched() -> CornerPointGeometry:
    return _cpg(nz=NZ, layer_dz=[12.0, 12.0], pinch_layer=0)


ALL_MODELS = {
    "uniform": _uniform,
    "variable_thickness": _variable_thickness,
    "dipping": _dipping,
    "warped": _warped,
    "non_orthogonal": _non_orthogonal,
}


# ══════════════════════════════════ §2 — lokal cell metrics ═════════════
@pytest.mark.parametrize("name", sorted(ALL_MODELS))
def test_every_cell_metric_is_per_cell_and_finite(name):
    geometry = ALL_MODELS[name]()
    ncell = geometry.grid.ncell
    extents = geometry.cell_extents()
    assert extents.shape == (ncell, 3)
    assert geometry.cell_thickness().shape == (ncell,)
    assert geometry.characteristic_length().shape == (ncell,)
    assert np.all(np.isfinite(extents)) and np.all(extents > 0.0)
    assert np.all(geometry.characteristic_length() > 0.0)


def test_variable_thickness_gives_each_layer_its_own_thickness():
    """§1/§2 — `mean(DZ)` DEYİL: hər layın öz qalınlığı."""
    geometry = _variable_thickness()
    thickness = geometry.cell_thickness().reshape(NZ, NY * NX)
    assert np.allclose(thickness[0], 6.0)
    assert np.allclose(thickness[1], 21.0)
    # Nominal orta (13.5) HEÇ BİR hüceyrənin həqiqi qalınlığı deyil.
    assert not np.any(np.isclose(geometry.cell_thickness(), 13.5))


def test_warped_model_has_genuinely_different_metrics_per_cell():
    """§9 — əyri modeldə hüceyrələr EYNİ ölçüdə OLMAMALIDIR."""
    geometry = _warped()
    extents = geometry.cell_extents()
    for axis in range(3):
        assert np.ptp(extents[:, axis]) > 1e-6, f"ox {axis} sabit qalıb"
    assert np.ptp(geometry.characteristic_length()) > 1e-6


def test_characteristic_length_is_the_cube_root_of_the_real_volume():
    geometry = _dipping()
    assert np.allclose(geometry.characteristic_length(),
                       np.cbrt(np.abs(geometry.volumes())))


def test_cell_extents_are_bounding_box_of_that_cells_own_nodes():
    geometry = _warped()
    for cell in (0, 5, geometry.grid.ncell - 1):
        nodes = geometry.cell_nodes(cell)
        expected = nodes.max(axis=0) - nodes.min(axis=0)
        assert np.allclose(geometry.cell_extents()[cell], expected)


def test_cartesian_metrics_match_the_scalar_values_exactly():
    """§ geriyə uyğunluq — bərabər bloklarda lokal ölçü = köhnə skalyar."""
    grid = _grid()
    box = CellGeometry(grid=grid, dx=DX, dy=DY, dz=[6.0, 21.0], top_depth=TOP)
    extents = box.cell_extents()
    assert np.allclose(extents[:, 0], DX)
    assert np.allclose(extents[:, 1], DY)
    assert np.allclose(extents[:, 2], box.dz_per_cell())
    assert np.allclose(box.cell_thickness(), box.dz_per_cell())

    cpg = CornerPointGeometry.from_cartesian(box)
    assert np.allclose(cpg.cell_extents(), extents)
    assert np.allclose(cpg.cell_thickness(), box.cell_thickness())
    assert np.allclose(cpg.characteristic_length(), box.characteristic_length())


# ══════════════════════════════════ §3 — pore volume ════════════════════
def _rock(ncell, poro=0.24, ntg=None):
    from imex2d.domain.properties import PropertyMap, RockProperties
    return RockProperties(
        porosity=PropertyMap.uniform("PORO", poro, ncell),
        permx=PropertyMap.uniform("PERMX", 150.0, ncell),
        permy=PropertyMap.uniform("PERMY", 150.0, ncell),
        permz=PropertyMap.uniform("PERMZ", 15.0, ncell),
        net_to_gross=(None if ntg is None
                      else PropertyMap.uniform("NTG", ntg, ncell)))


def _reservoir(geometry, ntg=None):
    from imex2d.domain.reservoir_model import ReservoirModel
    grid = geometry.grid
    return ReservoirModel(name="phase5e", grid=grid, geometry=geometry,
                          rock=_rock(grid.ncell, ntg=ntg), wells=[])


@pytest.mark.parametrize("name", sorted(ALL_MODELS))
def test_pore_volume_is_real_volume_times_poro_times_ntg(name):
    """§3 — PV = REAL corner-point həcm × PORO × NTG."""
    geometry = ALL_MODELS[name]()
    model = _reservoir(geometry, ntg=0.8)
    assert np.allclose(model.pore_volume(), geometry.volumes() * 0.24 * 0.8)
    assert np.all(model.pore_volume() > 0.0)


@pytest.mark.parametrize("name", sorted(ALL_MODELS))
def test_volume_never_uses_the_global_nominal_box(name):
    """§1 — həcm `mean(DX)·mean(DY)·mean(DZ)` QLOBAL qutusundan gəlmir.

    QEYD: xələnmə (maili lay, maili pillar) həcmi Kavalyeri prinsipinə
    görə DƏYİŞMİR — ona görə "həcm qutudan fərqlidir" UNİVERSAL bir
    iddia DEYİL və burada iddia edilmir. Universal olan budur: həcm
    QLOBAL orta ölçülərdən qurulmuş tək bir qutu dəyəri deyil, hər
    hüceyrənin ÖZ təpələrindən gəlir (bunu `warped` modeli, yəni ayaq
    izi dərinliklə dəyişən yeganə model, birbaşa göstərir)."""
    geometry = ALL_MODELS[name]()
    global_box = float(geometry.dx) * float(geometry.dy) * float(np.mean(geometry.dz))
    if name in ("variable_thickness", "warped"):
        assert not np.allclose(geometry.volumes(), global_box)
    # Hər halda: həcm hüceyrənin ÖZ təpələrindən yenidən hesablana bilir.
    from imex2d.domain.polyhedral_geometry import HexahedralCell
    for cell in (0, geometry.grid.ncell - 1):
        assert np.isclose(geometry.volumes()[cell],
                          HexahedralCell(geometry.cell_nodes(cell)).volume())


def test_warped_footprint_makes_volume_depart_from_the_box_formula():
    """Ayaq izi dərinliklə dəyişəndə həcm `sahə × qalınlıq` DEYİL."""
    geometry = _warped()
    box = DX * DY * geometry.cell_thickness()
    assert not np.allclose(geometry.volumes(), box)


def test_pore_volume_without_ntg_omits_the_factor():
    geometry = _dipping()
    model = _reservoir(geometry, ntg=None)
    assert np.allclose(model.pore_volume(), geometry.volumes() * 0.24)


def test_inactive_cells_have_zero_pore_volume_and_bulk_volume():
    """§3/§7 — ACTNUM = 0 → PV = 0, amma HƏNDƏSƏ (volumes()) dəyişmir."""
    actnum = np.ones(NX * NY * NZ, dtype=int)
    actnum[[0, 5, 17]] = 0
    grid = _grid(actnum=actnum)
    geometry = _cpg(grid=grid, dip_x=20.0, layer_dz=[8.0, 15.0])
    model = _reservoir(geometry, ntg=0.9)

    pv = model.pore_volume()
    assert np.allclose(pv[[0, 5, 17]], 0.0)
    assert np.all(pv[actnum > 0] > 0.0)
    assert np.allclose(model.bulk_volume()[[0, 5, 17]], 0.0)
    # Xalis həndəsə ACTNUM-dan XƏBƏRSİZDİR (qəsdən).
    assert np.all(geometry.volumes() > 0.0)


def test_pinched_out_cells_carry_zero_pore_volume():
    """§8 pinch-out — sıfır qalınlıqlı lay sıfır PV verir, NaN yox."""
    geometry = _pinched()
    model = _reservoir(geometry)
    pv = model.pore_volume()
    areal = NX * NY
    assert np.allclose(pv[:areal], 0.0)          # pinch-out layı
    assert np.all(pv[areal:] > 0.0)
    assert np.all(np.isfinite(pv))


# ══════════════════════════════════ §4/§5 — transmissibility ════════════
def _discretize(geometry, ntg=None):
    from imex2d.simulation.discretization import TwoPointFluxDiscretization
    return TwoPointFluxDiscretization().build(_reservoir(geometry, ntg))


@pytest.mark.parametrize("name", sorted(ALL_MODELS))
def test_transmissibility_uses_real_areas_and_real_half_distances(name):
    """§4/§5 — TPFA girişləri HƏNDƏSƏDƏN gəlir, nominal `dx/2` DEYİL."""
    geometry = ALL_MODELS[name]()
    conn = geometry.grid.build_connections()
    half_a, half_b = geometry.face_half_distances(conn)

    # Yarım-məsafə = |(üz mərkəzi − hüceyrə mərkəzi) · n̂| — birbaşa yoxla.
    centroids = geometry.cell_centroid()
    face_c = geometry.face_centroid(conn)
    normal = geometry.face_normal(conn)
    assert np.allclose(half_a, np.abs(np.einsum(
        "fj,fj->f", face_c - centroids[conn.cell_a], normal)))
    assert np.allclose(half_b, np.abs(np.einsum(
        "fj,fj->f", centroids[conn.cell_b] - face_c, normal)))
    assert np.all(half_a > 0.0) and np.all(half_b > 0.0)

    grid = _discretize(geometry)
    assert np.all(np.isfinite(grid.connections.cell_a))
    assert np.all(grid.transmissibility > 0.0)


def test_half_distances_are_not_the_nominal_half_cell_size_when_skewed():
    """§5 — `dx/2`/`dy/2`/`dz/2` fərziyyəsi POZULMALIDIR."""
    geometry = _non_orthogonal()
    conn = geometry.grid.build_connections()
    half_a, _ = geometry.face_half_distances(conn)
    nominal = np.where(conn.axis == 0, DX / 2.0,
                       np.where(conn.axis == 1, DY / 2.0, DZ / 2.0))
    assert not np.allclose(half_a, nominal)


def test_variable_thickness_changes_vertical_transmissibility():
    """§4 — qalın/nazik lay arasındakı Z üzü nominal `dz/2` DEYİL."""
    geometry = _variable_thickness()
    conn = geometry.grid.build_connections()
    z_faces = conn.axis == 2
    half_a, half_b = geometry.face_half_distances(conn)
    assert np.allclose(half_a[z_faces], 3.0)     # 6 m layın yarısı
    assert np.allclose(half_b[z_faces], 10.5)    # 21 m layın yarısı


def test_transmissibility_differs_between_uniform_and_dipping_models():
    """Maili modeldə üz sahəsi/məsafə dəyişir → T dəyişir. Əks halda
    həndəsə heç yerə ÖTÜRÜLMÜRDÜ demək olardı (§6)."""
    uniform = _discretize(_uniform()).transmissibility
    dipping = _discretize(_dipping()).transmissibility
    assert not np.allclose(uniform, dipping)


def test_cartesian_transmissibility_is_unchanged_by_the_cpg_path():
    """Bərabər blok modelində CPG yolu köhnə qutu yolu ilə EYNİ T verir."""
    grid = _grid()
    box = CellGeometry(grid=grid, dx=DX, dy=DY, dz=[6.0, 21.0], top_depth=TOP)
    cpg = CornerPointGeometry.from_cartesian(box)
    assert np.allclose(_discretize(box).transmissibility,
                       _discretize(cpg).transmissibility)


# ══════════════════════════════════ §6 — quyu həndəsəsi/indeksi ═════════
def _well_connections(geometry, i=1, j=1):
    from imex2d.domain.wells import (ControlMode, Perforation, Well, WellControl,
                                     WellType)
    from imex2d.simulation.well_model import PeacemanWellModel
    grid = geometry.grid
    well = Well(name="P1", well_type=WellType.PRODUCER,
                control=WellControl(ControlMode.BHP, 150.0),
                perforations=[Perforation(i, j, k, True, 0.0)
                              for k in range(grid.nz)], radius=0.1)
    model = _reservoir(geometry)
    object.__setattr__(model, "wells", [well])
    return PeacemanWellModel().build_connections(model)


def test_well_index_uses_each_cells_own_extents_not_a_global_dx():
    """§6 — dəyişkən qalınlıqda WI lay-lay FƏRQLİ olmalıdır (WI ∝ dz)."""
    geometry = _variable_thickness()
    connections = _well_connections(geometry)
    assert len(connections) == NZ
    wi = [c.well_index for c in connections]
    assert not np.isclose(wi[0], wi[1])
    assert np.isclose(wi[1] / wi[0], 21.0 / 6.0, rtol=1e-9)


def test_well_index_differs_between_uniform_and_warped_geometry():
    """§1 — global `dx`/`dy` işlədilsəydi bu iki model EYNİ WI verərdi."""
    uniform = [c.well_index for c in _well_connections(_uniform())]
    warped = [c.well_index for c in _well_connections(_warped())]
    assert not np.allclose(uniform, warped)


def test_well_index_on_cartesian_model_is_unchanged():
    """Bərabər bloklarda lokal ölçü = köhnə skalyar → WI DƏYİŞMİR."""
    grid = _grid()
    box = CellGeometry(grid=grid, dx=DX, dy=DY, dz=[6.0, 21.0], top_depth=TOP)
    cpg = CornerPointGeometry.from_cartesian(box)
    assert np.allclose([c.well_index for c in _well_connections(box)],
                       [c.well_index for c in _well_connections(cpg)])


def test_well_in_inactive_cell_is_dropped_with_variable_geometry():
    """§7 — ACTNUM + dəyişkən həndəsə: qeyri-aktiv laya düşən
    perforasiya WI siyahısına DÜŞMÜR (indeks sürüşməsi olmadan)."""
    actnum = np.ones(NX * NY * NZ, dtype=int)
    grid_probe = _grid()
    dead = grid_probe.index(1, 1, 0)
    actnum[dead] = 0
    grid = _grid(actnum=actnum)
    geometry = _cpg(grid=grid, layer_dz=[6.0, 21.0], dip_x=15.0)

    connections = _well_connections(geometry, i=1, j=1)
    assert len(connections) == NZ - 1
    assert all(c.cell != dead for c in connections)
    # Qalan bağlantı 2-ci layın HƏQİQİ qalınlığını görür.
    assert connections[0].cell == grid.index(1, 1, 1)


# ══════════════════════════════════ §6 — dərinlik/sütun axtarışı ════════
def test_layer_edges_come_from_the_real_column_not_average_dz():
    """§1 — `layer_edges` sütunun ÖZ künc dərinliklərindən qurulur."""
    geometry = _dipping()
    x, y = 3.5 * DX, 0.5 * DY                    # i = 3 sütunu
    edges = layer_edges(x, y, geometry)
    assert edges.size == NZ + 1
    # i = 3 sütunu dip_x·3.5 ≈ 87.5 m daha dərindədir (künc ortalaması).
    shallow = layer_edges(0.5 * DX, y, geometry)
    assert edges[0] > shallow[0] + 50.0


def test_layer_edges_track_variable_thickness_per_layer():
    geometry = _variable_thickness()
    edges = layer_edges(1.5 * DX, 1.5 * DY, geometry)
    assert np.allclose(np.diff(edges), [6.0, 21.0])


def test_depth_to_k_uses_the_real_column_edges():
    geometry = _variable_thickness()
    x, y = 1.5 * DX, 1.5 * DY
    assert depth_to_k(x, y, TOP + 3.0, geometry) == 0        # 6 m layın içi
    assert depth_to_k(x, y, TOP + 15.0, geometry) == 1       # 21 m layın içi
    assert depth_to_k(x, y, TOP + 500.0, geometry) is None   # grid-dən kənar


def test_interval_layers_tolerance_follows_the_columns_own_thickness():
    geometry = _variable_thickness()
    x, y = 1.5 * DX, 1.5 * DY
    assert interval_layers(x, y, TOP + 1.0, TOP + 4.0, geometry) == [0]
    assert interval_layers(x, y, TOP + 1.0, TOP + 20.0, geometry) == [0, 1]
    # Lay SƏRHƏDİ üzərindəki sıfır-örtmə qonşu layı DAXİL ETMİR.
    assert interval_layers(x, y, TOP + 6.0, TOP + 20.0, geometry) == [1]


def test_locate_column_follows_the_real_footprint_on_a_shifted_grid():
    """§1 — nominal `int(x/dx)` sürüşmüş ayaq izində yanlış sütun verir."""
    geometry = _cpg(tilt_x=0.0)
    # Şaquli pillarlarda ayaq izi nominal şəbəkə ilə üst-üstə düşür.
    assert xy_to_ij(2.5 * DX, 1.5 * DY, geometry) == (2, 1)
    assert xy_to_ij(0.0, 0.0, geometry) == (0, 0)
    # Grid-dən kənar nöqtə ən yaxın sütuna düşür, XƏTA ATILMIR.
    assert xy_to_ij(-500.0, -500.0, geometry) == (0, 0)
    assert xy_to_ij(1e6, 1e6, geometry) == (NX - 1, NY - 1)


def test_locate_column_matches_cartesian_for_an_unshifted_grid():
    grid = _grid()
    box = CellGeometry(grid=grid, dx=DX, dy=DY, dz=DZ, top_depth=TOP)
    cpg = CornerPointGeometry.from_cartesian(box)
    for x, y in ((10.0, 10.0), (250.0, 100.0), (399.0, 239.0), (0.0, 0.0)):
        assert cpg.locate_column(x, y) == box.locate_column(x, y)


def test_locate_column_finds_the_cell_containing_its_own_centroid():
    """Ən sərt yoxlama: hər sütunun mərkəzi MƏHZ o sütuna düşməlidir."""
    geometry = _warped()
    centroids = geometry.cell_centroid()
    grid = geometry.grid
    for j in range(grid.ny):
        for i in range(grid.nx):
            centre = centroids[grid.index(i, j, 0)]
            assert geometry.locate_column(centre[0], centre[1]) == (i, j)


# ══════════════════════════════════ §7 — indeks ardıcıllığı ═════════════
def test_geometry_property_and_active_indexing_stay_aligned():
    """§7 — global cell → active cell → geometry cell → property cell.

    Hər hüceyrəyə UNİKAL qalınlıq verilir, sonra aktiv altmassivin
    həndəsəsi ilə xassəsi eyni hüceyrəni göstərməlidir."""
    from imex2d.domain.properties import PropertyMap

    actnum = np.ones(NX * NY * NZ, dtype=int)
    actnum[[2, 9, 20]] = 0
    grid = _grid(actnum=actnum)
    geometry = _cpg(grid=grid, layer_dz=[7.0, 19.0], dip_x=12.0)

    # Xassə = hüceyrənin QLOBAL indeksi (izlənə bilən marker).
    marker = PropertyMap.from_array("MARK", np.arange(grid.ncell, dtype=float),
                                    grid.ncell)
    active = grid.active
    assert np.array_equal(active.to_active(marker.values),
                          active.active_to_global.astype(float))

    # Aktiv altmassivdə həndəsə də EYNİ hüceyrələri verir.
    thickness_global = geometry.cell_thickness()
    assert np.allclose(active.to_active(thickness_global),
                       thickness_global[active.active_to_global])
    volumes_global = geometry.volumes()
    assert np.allclose(active.to_active(volumes_global),
                       volumes_global[active.active_to_global])

    # Geri yayılma: aktiv → qlobal, qeyri-aktiv 0.
    restored = active.to_global(active.to_active(marker.values))
    assert np.allclose(restored[active.mask], marker.values[active.mask])
    assert np.allclose(restored[~active.mask], 0.0)


def test_cell_metrics_are_indexed_by_global_cell_id():
    """Bütün lokal ölçülər QLOBAL (ncell) indekslidir — xassələrlə eyni."""
    actnum = np.ones(NX * NY * NZ, dtype=int)
    actnum[0] = 0
    grid = _grid(actnum=actnum)
    geometry = _cpg(grid=grid, layer_dz=[7.0, 19.0])
    for values in (geometry.cell_thickness(), geometry.characteristic_length(),
                   geometry.volumes(), geometry.cell_depths()):
        assert values.shape == (grid.ncell,)
    assert geometry.cell_extents().shape == (grid.ncell, 3)
    assert geometry.cell_centroid().shape == (grid.ncell, 3)
    # Qeyri-aktiv hüceyrənin HƏNDƏSƏSİ hələ də mövcuddur (silinmir).
    assert geometry.volumes()[0] > 0.0


# ══════════════════════════════════ §6 — CFL / zaman addımı ═════════════
def test_cfl_ratio_is_driven_by_the_real_pore_volume():
    """§6 — CFL addımı `PV / throughput`-dur; PV real həcmdən gəldiyi
    üçün dəyişkən həndəsə CFL-ə AVTOMATİK ötürülür."""
    uniform = _reservoir(_uniform()).pore_volume()
    variable = _reservoir(_variable_thickness()).pore_volume()
    assert not np.allclose(uniform, variable)
    # Nazik lay daha kiçik PV → daha sərt CFL həddi verir.
    areal = NX * NY
    assert variable[:areal].max() < uniform[:areal].min()


# ══════════════════════════════════ §9 — bütöv zəncir ═══════════════════
def test_full_chain_carries_real_geometry_end_to_end():
    """§9 — Corner-Point → Local Metrics → Real Volume → Real Distances
    → Downstream. Hər addımda qutu düsturundan FƏRQ olmalıdır."""
    geometry = _non_orthogonal()
    grid = geometry.grid
    conn = grid.build_connections()

    # Lokal ölçü pillar meylini TUTUR: hüceyrənin X uzanması deck-in
    # pillar addımından (DX) böyükdür, çünki hüceyrə xələnib.
    extents = geometry.cell_extents()
    assert np.all(extents[:, 0] > DX + 1e-6)
    # Üz həndəsəsi nominal qutu üzündən FƏRQLİDİR (xələnmə həcmi
    # saxlayır, amma sahə/məsafə/normalı SAXLAMIR — bax
    # `test_volume_never_uses_the_global_nominal_box`).
    box_area = np.where(conn.axis == 0, DY * DZ,
                        np.where(conn.axis == 1, DX * DZ, DX * DY))
    assert not np.allclose(geometry.face_areas(conn), box_area)

    model = _reservoir(geometry, ntg=0.85)
    assert np.allclose(model.pore_volume(), geometry.volumes() * 0.24 * 0.85)

    discretized = _discretize(geometry, ntg=0.85)
    assert np.allclose(discretized.cell_volume, geometry.volumes())
    assert np.all(discretized.transmissibility > 0.0)

    # Üz normalları ox vektorları DEYİL → həqiqi qeyri-ortoqonallıq.
    normals = geometry.face_normal(conn)
    assert np.any(np.abs(normals).max(axis=1) < 0.999)
    assert geometry.quality_metrics(conn)["max_non_orthogonality_angle_deg"] > 1.0
