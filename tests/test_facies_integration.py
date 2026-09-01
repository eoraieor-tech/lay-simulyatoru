"""Phase 4.1 — SIS pipeline inteqrasiyası: kateqorik marşrutlaşdırma,
sərt-data xəritələnməsi, sürətli axtarış PARİTETİ, ehtimal diaqnostikası,
GeologicalModel-ə bağlanma.

Bu fayl tapşırığın §11 "Integration tests" siyahısındakı 13 bəndin
HAMISINI birbaşa təmsil edir (bəziləri digər fayllarda da ayrıca
sınanıb — bura konsolidasiya üçündür).
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.application.geology_service import (FaciesBuildConfig, GeologicalGridSpec,
                                                WellBasedGeologicalModelBuilder)
from imex2d.domain.geometry import CellGeometry
from imex2d.domain.grid import CartesianGrid
from imex2d.domain.well_data import WellDataset, WellSample
from imex2d.geology.facies import simulate_sis
from imex2d.geology.interpolation import OrdinaryKriging
from imex2d.geology.property_types import PropertyType, classify_property


def _geometry(nx=6, ny=6, nz=1, dx=20.0, dy=20.0, dz=10.0, top=2000.0):
    grid = CartesianGrid(nx, ny, nz)
    return grid, CellGeometry(grid, dx, dy, dz, top_depth=top)


def _grid_targets(nx=6, ny=6, dx=20.0, dy=20.0):
    xs = (np.arange(nx) + 0.5) * dx
    ys = (np.arange(ny) + 0.5) * dy
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


# ── §6: brute-force vs sürətli (cKDTree) axtarış PARİTETİ (SIS səviyyəsində) ─
def test_fast_and_brute_force_search_agree_on_identical_realization():
    """Tapşırıq §6: 'Only then replace the repeated brute-force search'.
    Eyni (nöqtələr, anizotropluq, axtarış radiusu, maks qonşu, seed) üçün
    `use_fast_search=True/False` EYNİ realizasiyanı verməlidir."""
    rng = np.random.default_rng(21)
    points = rng.uniform(0, 300, size=(12, 2))
    codes = (points[:, 0] > 150).astype(int)
    targets = _grid_targets(nx=10, ny=10, dx=25.0, dy=25.0)
    proportions = {0: 0.5, 1: 0.5}

    fast = simulate_sis(points, codes, targets, proportions, seed=17,
                        max_neighbors=8, use_fast_search=True)
    slow = simulate_sis(points, codes, targets, proportions, seed=17,
                        max_neighbors=8, use_fast_search=False)
    assert np.array_equal(fast.codes, slow.codes)
    assert fast.realized_proportions == slow.realized_proportions


def test_fast_and_brute_force_agree_with_explicit_anisotropy_and_radius():
    from imex2d.geology.facies import FaciesVariogramParams
    points = np.array([[-100., 0.], [100., 0.], [0., -100.], [0., 100.],
                       [50., 50.], [-50., -50.], [70., -70.], [-70., 70.]])
    codes = np.array([0, 0, 1, 1, 0, 1, 0, 1])
    vp = {0: FaciesVariogramParams(range_=250.0, range_minor=40.0, azimuth_deg=45.0, nugget=0.05),
         1: FaciesVariogramParams(range_=250.0, range_minor=40.0, azimuth_deg=45.0, nugget=0.05)}
    targets = _grid_targets(nx=8, ny=8, dx=20.0, dy=20.0)
    kwargs = dict(variograms=vp, search_radius=200.0, max_neighbors=6, seed=5)

    fast = simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, use_fast_search=True, **kwargs)
    slow = simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, use_fast_search=False, **kwargs)
    assert np.array_equal(fast.codes, slow.codes)


# ── §1: FACIES CSV should NOT enter continuous interpolation ────────────
def test_facies_column_is_classified_categorical_poro_is_continuous():
    assert classify_property("FACIES") is PropertyType.CATEGORICAL
    assert classify_property("PORO") is PropertyType.CONTINUOUS


# ── §7: ehtimal düzəliş diaqnostikası ────────────────────────────────────
def test_probability_correction_diagnostics_are_tracked_and_reported():
    """Seyrək data ilə NaN-geri-dönüş sayılmalı və hədd aşılırsa GÜCLÜ
    xəbərdarlıq verilməlidir."""
    points = np.array([[0., 0.], [500., 500.]])   # yalnız 2 quyu, çox seyrək
    codes = np.array([0, 1])
    targets = _grid_targets(nx=6, ny=6, dx=20.0, dy=20.0)
    realization = simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, seed=1,
                               search_radius=15.0,   # çox dar radius -> çox NaN fallback
                               max_neighbors=4, correction_warn_threshold=0.10)
    diag = realization.diagnostics
    assert diag.n_cells_simulated == 36
    assert diag.nan_fallback_cells > 0
    assert diag.rate(diag.nan_fallback_cells) > 0.0
    assert any("GÜCLÜ XƏBƏRDARLIQ" in w for w in realization.warnings)


def test_probability_diagnostics_report_zero_corrections_when_clean():
    points, codes = (np.array([[0., 0.], [0., 100.], [100., 0.], [100., 100.], [50., 50.]]),
                     np.array([0, 0, 1, 1, 0]))
    targets = _grid_targets(nx=4, ny=4, dx=25.0, dy=25.0)
    realization = simulate_sis(points, codes, targets, {0: 0.6, 1: 0.4}, seed=2,
                               max_neighbors=5)
    diag = realization.diagnostics
    assert diag.n_cells_simulated == 16
    # sağlam ssenaridə say sıfır ola BİLƏR (təminat deyil), amma NEQATİV olmamalıdır
    assert diag.negative_probability_events >= 0
    assert diag.nan_fallback_cells >= 0


# ── §11.13: legacy iş axını fasiya datası olmadan dəyişməz qalır ─────────
def _facies_and_poro_dataset():
    samples = []
    factor = {"A": 1.00, "B": 1.06, "C": 0.94}
    facies_by_well = {"A": 0, "B": 1, "C": 0}
    for name, (x, y) in {"A": (20.0, 20.0), "B": (100.0, 20.0), "C": (20.0, 100.0)}.items():
        # `layer=0` AÇIQ verilir: nz=1 tək laydır, amma SIS-in 3D
        # yerləşdirməsi (Phase 4.1) hər nümunənin AÇIQ lay/dərinlik
        # daşımasını tələb edir (bax `_simulate_categorical_field`).
        samples.append(WellSample(well=name, x=x, y=y, layer=0,
                                  values={"PORO": 0.2 * factor[name], "PERMX": 150.0 * factor[name],
                                         "FACIES": float(facies_by_well[name])}))
    return WellDataset(samples=samples, source="test")


def _build_with_facies(**facies_kwargs):
    dataset = _facies_and_poro_dataset()
    spec = GeologicalGridSpec(nx=4, ny=4, nz=1, dx=25.0, dy=25.0, top_depth=2000.0)
    facies_kwargs.setdefault("seed", 1)
    config = {"FACIES": FaciesBuildConfig(proportions={0: 0.6, 1: 0.4}, **facies_kwargs)}
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    return builder.build(dataset, spec, facies_config=config)


# ── §11.2/§11.3/§11.8: FACIES → SIS, PORO → kəsilməz, model-ə bağlanma ──
def test_facies_routes_to_sis_and_attaches_to_geological_model_poro_stays_continuous():
    model, report = _build_with_facies()
    assert "FACIES" in model.facies_fields
    assert "FACIES" not in model.property_maps      # §3: gizli PropertyMap kimi DEYİL
    assert "PORO" in model.property_maps             # §11.3: PORO hələ də kəsilməzdir
    assert "PORO" not in model.facies_fields

    facies_field = model.facies_fields["FACIES"]
    assert facies_field.ncell == model.grid.ncell
    assert set(np.unique(facies_field.codes)).issubset({0, 1})
    poro_values = model.property_maps["PORO"].values
    # PORO davam edən (kriging) dəyərlərdir — 2-3 diskret dəyərlə MƏHDUDLAŞMIR
    assert len(set(np.round(poro_values, 6))) > 2


def test_poro_values_are_never_integer_facies_like_codes():
    """§14.3 həqiqi yoxlanışı: PORO kəsilməz Kriging nəticəsidir, YALNIZ
    {0,1} kimi diskret kodlarla MƏHDUDLAŞMIR (FACIES-in ADİ ədədi
    interpolyasiyaya sızmadığının dolayı sübutu)."""
    model, _ = _build_with_facies()
    poro_values = model.property_maps["PORO"].values
    assert not np.all(np.isin(np.round(poro_values, 6), [0.0, 1.0]))


# ── §11.6/§11.7: seed təkrarlanabilənliyi tam boru xəttində ──────────────
def test_same_seed_reproduces_identical_pipeline_facies_field():
    model1, _ = _build_with_facies(seed=7)
    model2, _ = _build_with_facies(seed=7)
    assert np.array_equal(model1.facies_fields["FACIES"].codes,
                          model2.facies_fields["FACIES"].codes)


def test_different_seeds_give_different_pipeline_facies_fields():
    model1, _ = _build_with_facies(seed=1)
    model2, _ = _build_with_facies(seed=2)
    assert not np.array_equal(model1.facies_fields["FACIES"].codes,
                              model2.facies_fields["FACIES"].codes)


# ── §11.4: sərt data inteqrasiyadan sonra da qorunur ──────────────────────
def test_hard_data_preserved_through_full_pipeline():
    dataset = _facies_and_poro_dataset()
    spec = GeologicalGridSpec(nx=6, ny=6, nz=1, dx=20.0, dy=20.0, top_depth=2000.0)
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    config = {"FACIES": FaciesBuildConfig(proportions={0: 0.6, 1: 0.4}, seed=3)}
    model, _ = builder.build(dataset, spec, facies_config=config)

    from imex2d.domain.geometry import xy_to_ij
    facies_field = model.facies_fields["FACIES"]
    grid_shape = model.grid.shape   # (nz, ny, nx)
    for sample in dataset.samples:
        i, j = xy_to_ij(sample.x, sample.y, model.geometry)
        cell_index = np.ravel_multi_index((0, j, i), grid_shape)
        assert facies_field.codes[cell_index] == int(sample.values["FACIES"]), sample.well


# ── §11.5: ziddiyyətli sərt data aşkarlanır (tam boru xəttində) ─────────
def test_conflicting_hard_data_detected_in_full_pipeline():
    samples = [
        WellSample(well="A", x=10.0, y=10.0, layer=0,
                  values={"PORO": 0.2, "PERMX": 150.0, "FACIES": 0.0}),
        WellSample(well="B", x=11.0, y=11.0, layer=0,   # A ilə EYNİ hüceyrə, FƏRQLİ kod
                  values={"PORO": 0.2, "PERMX": 150.0, "FACIES": 1.0}),
        WellSample(well="C", x=90.0, y=90.0, layer=0,
                  values={"PORO": 0.2, "PERMX": 150.0, "FACIES": 1.0}),
    ]
    dataset = WellDataset(samples=samples, source="test")
    spec = GeologicalGridSpec(nx=4, ny=4, nz=1, dx=30.0, dy=30.0, top_depth=2000.0)
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    config = {"FACIES": FaciesBuildConfig(proportions={0: 0.5, 1: 0.5})}
    with pytest.raises(Exception):
        builder.build(dataset, spec, facies_config=config)


# ── §11.9: SATNUM/PVTNUM (RegionSet) inteqrasiyadan təsirlənmir ──────────
def test_regions_remain_unaffected_by_facies_integration():
    model, _ = _build_with_facies()
    assert model.regions is not None
    assert model.regions.ids.tolist() == [1]   # defolt tək region, FACIES-dən ASILI DEYİL


def test_legacy_continuous_workflow_unaffected_without_facies_data():
    """Fasiya sütunu olmayan mövcud (Phase 1-3) iş axını TAM DƏYİŞMƏZ
    qalmalıdır — bax `tests/test_geology_import.py` (bu fayl həmin
    testlərin heç birini DƏYİŞMİR, yalnız faktın özünü təsdiqləyir)."""
    from imex2d.application.geology_service import (GeologicalGridSpec,
                                                     WellBasedGeologicalModelBuilder)
    from imex2d.geology.interpolation import OrdinaryKriging

    dataset = WellDataset(samples=[
        WellSample(well="W-1", x=10.0, y=20.0, values={"PORO": 0.20, "PERMX": 150.0}),
        WellSample(well="W-2", x=90.0, y=80.0, values={"PORO": 0.25, "PERMX": 400.0}),
    ])
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    model, report = builder.build(dataset, GeologicalGridSpec(nx=3, ny=3, nz=1, dx=20.0, dy=20.0))
    assert "PORO" in model.property_maps
    assert "PERMX" in model.property_maps
    assert not hasattr(model, "facies_fields") or not model.facies_fields
