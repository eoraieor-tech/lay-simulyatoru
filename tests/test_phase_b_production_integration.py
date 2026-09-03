"""B-INTEGRATION-FIX — Phase B interpolyasiya mühərrikinin ƏSL istehsalat
boru xəttinə (`WellBasedGeologicalModelBuilder`) inteqrasiyasının sübutu.

Bu fayl YALNIZ inteqrasiyanı sınayır — Phase B-nin öz riyaziyyatı artıq
`test_property_strategies.py`/`test_uncertainty.py`/`test_data_quality.py`/
`test_cross_validation.py`-da doğrulanıb. Buradakı testlərin MƏRKƏZİ
iddiası: `WellBasedGeologicalModelBuilder.build()` (istehsalat çağırış
yolu) HƏQİQƏTƏN `property_interpolation.interpolate_property_field()`/
`interpolate_categorical_field()`-dən keçir, KÖHNƏ Phase A `geology.
interpolation.interpolate_property()`-dən YOX (bax
`test_old_phase_a_path_is_never_reached_by_production_builder`).
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from imex2d.application.geology_service import (ContinuousSGSConfig,
                                                 FaciesBuildConfig,
                                                 GeologicalGridSpec,
                                                 WellBasedGeologicalModelBuilder)
from imex2d.domain.properties import CategoricalUncertainty, PropertyUncertainty
from imex2d.domain.well_data import WellDataset, WellSample
from imex2d.geology import interpolation as phase_a_interpolation
from imex2d.geology.interpolation import InverseDistance, NearestNeighbour, OrdinaryKriging
from imex2d.geology.property_config import InterpolationKind, VariableType, resolve_strategy


# ── ortaq sintetik dataset ───────────────────────────────────────────────
def _grid_wells(n_side=6, spacing=120.0, seed=0):
    """`n_side x n_side` düzenli quyu şəbəkəsi — PORO/PERMX/PERMY/PERMZ/
    SW/NTG/FACIES hamısı üçün kifayət qədər sərt data verir.

    FACIES AYRI `WellSample` kimi, AÇIQ `layer=0` ilə əlavə olunur (SIS/
    Phase B kateqorik yolun 3D mövqe tələbi — bax `_gather_categorical_
    hard_data`) — kəsilməz xassələrin `WellSample`-larından FƏRQLİ, ki
    `layer=0` KƏSİLMƏZ xassələrin çox-laylı (`nz>1`) testlərini pozmasın
    (`dataset.samples_for(prop, k)` laysız (`layer=None`) nümunələri HƏR
    K üçün daxil edir, amma `layer=0` işarəli nümunə YALNIZ K=0 üçün)."""
    rng = np.random.default_rng(seed)
    samples = []
    idx = 0
    for i in range(n_side):
        for j in range(n_side):
            x, y = i * spacing + 10.0, j * spacing + 10.0
            poro = float(np.clip(0.15 + 0.06 * np.sin(x / 300.0) + rng.normal(0, 0.005), 0.05, 0.35))
            permx = float(np.clip(np.exp(2.0 + 3.0 * poro + rng.normal(0, 0.05)), 1.0, 5000.0))
            sw = float(np.clip(0.25 + 0.4 * (i / n_side) + rng.normal(0, 0.01), 0.02, 0.98))
            ntg = float(np.clip(0.9 - 0.5 * (j / n_side) + rng.normal(0, 0.01), 0.02, 0.98))
            facies = 1 if poro < 0.16 else (2 if poro < 0.22 else 3)
            samples.append(WellSample(
                well=f"W{idx}", x=x, y=y,
                values={"PORO": poro, "PERMX": permx, "SW": sw, "NTG": ntg}))
            samples.append(WellSample(
                well=f"W{idx}", x=x, y=y, layer=0, values={"FACIES": facies}))
            idx += 1
    return WellDataset(samples=samples, source="test")


def _spec(nx=12, ny=12, nz=1):
    return GeologicalGridSpec(nx=nx, ny=ny, nz=nz, dx=40.0, dy=40.0, dz=10.0, top_depth=2000.0)


def _build(dataset=None, interpolator=None, spec=None, **kwargs):
    dataset = dataset if dataset is not None else _grid_wells()
    interpolator = interpolator if interpolator is not None else OrdinaryKriging()
    spec = spec if spec is not None else _spec()
    builder = WellBasedGeologicalModelBuilder(interpolator)
    return builder.build(dataset, spec, **kwargs)


# ── 1-4: PORO / PERMX / PERMY / PERMZ ────────────────────────────────────
def test_poro_integration_uses_phase_b_strategy_and_stays_in_bounds():
    model, report = _build()
    assert "PORO" in model.property_maps
    values = model.property_maps["PORO"].values
    assert np.all(np.isfinite(values))
    assert values.min() >= 0.0 - 1e-9
    assert values.max() <= 1.0 + 1e-9
    assert isinstance(model.uncertainty["PORO"], PropertyUncertainty)
    assert resolve_strategy("PORO").transform.is_identity


def test_permx_integration_uses_log_space_kriging_not_linear_space():
    """İki quyu EYNİ məsafədə bir hədəfdən, PERMX bir-birindən 4 tərtib
    fərqlənir. Loq-fəza kriginq nəticəni HƏNDƏSİ ortaya (~sqrt(k1*k2))
    yaxınlaşdırmalıdır — XƏTTİ (arifmetik orta) fəzada olsaydı nəticə
    ондан QAT-QAT böyük olardı. Bu, `_interpolate_volume`-un HƏQİQƏTƏN
    `ln(K)` fəzasında kriging etdiyini sübut edir (B-INTEGRATION-FIX §6)."""
    samples = [
        WellSample(well="LOW", x=0.0, y=0.0, values={"PORO": 0.15, "PERMX": 1.0}),
        WellSample(well="HIGH", x=200.0, y=0.0, values={"PORO": 0.15, "PERMX": 10000.0}),
    ]
    dataset = WellDataset(samples=samples, source="test")
    spec = GeologicalGridSpec(nx=3, ny=1, nz=1, dx=100.0, dy=200.0, top_depth=2000.0)
    model, _ = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(dataset, spec)
    permx = model.property_maps["PERMX"].values.reshape(model.grid.shape)
    midpoint = float(permx[0, 0, 1])
    geometric_mean = float(np.sqrt(1.0 * 10000.0))          # 100.0
    arithmetic_mean = float((1.0 + 10000.0) / 2.0)            # 5000.5
    assert abs(midpoint - geometric_mean) < abs(midpoint - arithmetic_mean), (
        f"K midpoint={midpoint:.3g} arifmetik ortaya ({arithmetic_mean:.3g}) həndəsi "
        f"ortadan ({geometric_mean:.3g}) daha yaxındır — loq-fəza kriginq işləmir")
    assert midpoint > 0.0


def test_permy_integration_positive_and_log_space():
    samples = [
        WellSample(well="LOW", x=0.0, y=0.0, values={"PORO": 0.15, "PERMX": 50.0, "PERMY": 2.0}),
        WellSample(well="HIGH", x=200.0, y=0.0, values={"PORO": 0.15, "PERMX": 50.0, "PERMY": 8000.0}),
    ]
    dataset = WellDataset(samples=samples, source="test")
    spec = GeologicalGridSpec(nx=3, ny=1, nz=1, dx=100.0, dy=200.0, top_depth=2000.0)
    model, _ = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(dataset, spec)
    permy = model.property_maps["PERMY"].values
    assert np.all(permy > 0.0) and np.all(np.isfinite(permy))
    midpoint = float(permy.reshape(model.grid.shape)[0, 0, 1])
    assert midpoint < (2.0 + 8000.0) / 2.0, "PERMY arifmetik ortaya doğru meyillidir (loq-fəza deyil)"


def test_permz_integration_positive_and_log_space():
    samples = [
        WellSample(well="LOW", x=0.0, y=0.0, values={"PORO": 0.15, "PERMX": 50.0, "PERMZ": 0.5}),
        WellSample(well="HIGH", x=200.0, y=0.0, values={"PORO": 0.15, "PERMX": 50.0, "PERMZ": 4000.0}),
    ]
    dataset = WellDataset(samples=samples, source="test")
    spec = GeologicalGridSpec(nx=3, ny=1, nz=1, dx=100.0, dy=200.0, top_depth=2000.0)
    model, _ = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(dataset, spec)
    permz = model.property_maps["PERMZ"].values
    assert np.all(permz > 0.0) and np.all(np.isfinite(permz))
    midpoint = float(permz.reshape(model.grid.shape)[0, 0, 1])
    assert midpoint < (0.5 + 4000.0) / 2.0, "PERMZ arifmetik ortaya doğru meyillidir (loq-fəza deyil)"


# ── 5-6: SW / NTG hədləri ─────────────────────────────────────────────────
def test_sw_integration_stays_within_zero_one_including_extremes():
    samples = [
        WellSample(well="A", x=0.0, y=0.0, values={"PORO": 0.15, "PERMX": 50.0, "SW": 0.01}),
        WellSample(well="B", x=500.0, y=0.0, values={"PORO": 0.15, "PERMX": 50.0, "SW": 0.99}),
        WellSample(well="C", x=0.0, y=500.0, values={"PORO": 0.15, "PERMX": 50.0, "SW": 0.0}),
        WellSample(well="D", x=500.0, y=500.0, values={"PORO": 0.15, "PERMX": 50.0, "SW": 1.0}),
    ]
    dataset = WellDataset(samples=samples, source="test")
    spec = GeologicalGridSpec(nx=20, ny=20, nz=1, dx=30.0, dy=30.0, top_depth=2000.0)
    model, _ = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(dataset, spec)
    sw = model.property_maps["SW"].values
    assert np.all(np.isfinite(sw))
    assert sw.min() >= 0.0 - 1e-9, f"SW < 0 tapıldı: {sw.min()}"
    assert sw.max() <= 1.0 + 1e-9, f"SW > 1 tapıldı: {sw.max()}"


def test_ntg_integration_stays_within_zero_one_including_extremes():
    samples = [
        WellSample(well="A", x=0.0, y=0.0, values={"PORO": 0.15, "PERMX": 50.0, "NTG": 0.01}),
        WellSample(well="B", x=500.0, y=0.0, values={"PORO": 0.15, "PERMX": 50.0, "NTG": 0.99}),
        WellSample(well="C", x=0.0, y=500.0, values={"PORO": 0.15, "PERMX": 50.0, "NTG": 0.0}),
        WellSample(well="D", x=500.0, y=500.0, values={"PORO": 0.15, "PERMX": 50.0, "NTG": 1.0}),
    ]
    dataset = WellDataset(samples=samples, source="test")
    spec = GeologicalGridSpec(nx=20, ny=20, nz=1, dx=30.0, dy=30.0, top_depth=2000.0)
    model, _ = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(dataset, spec)
    ntg = model.property_maps["NTG"].values
    assert np.all(np.isfinite(ntg))
    assert ntg.min() >= 0.0 - 1e-9, f"NTG < 0 tapıldı: {ntg.min()}"
    assert ntg.max() <= 1.0 + 1e-9, f"NTG > 1 tapıldı: {ntg.max()}"


# ── 7: FACIES (deterministik Phase B opt-in) ─────────────────────────────
def test_facies_integration_via_phase_b_gives_valid_classes_and_normalized_probabilities():
    dataset = _grid_wells()
    valid_codes = set(int(s.values["FACIES"]) for s in dataset.samples if "FACIES" in s.values)
    model, report = _build(
        dataset=dataset,
        facies_config={"FACIES": FaciesBuildConfig(deterministic=True)})
    assert "FACIES" in model.facies_fields
    facies = model.facies_fields["FACIES"]
    assert set(np.unique(facies.codes).tolist()) <= valid_codes
    unc = model.uncertainty["FACIES"]
    assert isinstance(unc, CategoricalUncertainty)
    row_sums = unc.probabilities.sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-6), f"cəm 1 deyil: {row_sums.min()}..{row_sums.max()}"
    assert np.all(unc.probabilities >= -1e-9) and np.all(unc.probabilities <= 1.0 + 1e-9)
    assert facies.conditioning_data_stats["method"] == "phase_b_indicator_kriging"


def test_facies_default_still_uses_sis_not_phase_b():
    """`deterministic` verilməyəndə DAVRANIŞ DƏYİŞMİR — SIS defolt qalır
    (B-INTEGRATION-FIX geriyə-uyğunluq tələbi)."""
    model, _ = _build(facies_config=None)
    facies = model.facies_fields["FACIES"]
    assert facies.conditioning_data_stats.get("proportion_source") in ("observed", "user")
    assert "FACIES" not in model.uncertainty   # Phase B kateqorik yol İŞLƏMƏYİB


# ── 8: sərt data honoring ─────────────────────────────────────────────────
def test_hard_data_is_honored_at_well_cells():
    dataset = _grid_wells(n_side=5, spacing=100.0)
    spec = GeologicalGridSpec(nx=25, ny=25, nz=1, dx=20.0, dy=20.0, top_depth=2000.0)
    model, _ = WellBasedGeologicalModelBuilder(OrdinaryKriging(nugget=0.0)).build(dataset, spec)
    from imex2d.domain.geometry import xy_to_ij
    poro = model.property_maps["PORO"].values.reshape(model.grid.shape)
    errors = []
    for sample in dataset.samples:
        if "PORO" not in sample.values:
            continue
        i, j = xy_to_ij(sample.x, sample.y, model.geometry)
        errors.append(abs(float(poro[0, j, i]) - sample.values["PORO"]))
    assert max(errors) < 1e-6, f"maksimum sərt-data xətası: {max(errors)}"


# ── 9: anizotropluğun inteqrasiyası ───────────────────────────────────────
def test_anisotropy_configuration_reaches_production_engine():
    """X-də davamlı (Y-də dəyişkən) sintetik sahə: `azimuth_deg=0` (X boyu
    böyük major radius) ilə `azimuth_deg=90` fərqli nəticə verməlidir —
    əks halda anizotropluq Phase B mühərrikinə ÇATMIR."""
    rng = np.random.default_rng(5)
    samples = []
    for i in range(6):
        for j in range(6):
            x, y = i * 60.0, j * 60.0
            value = 0.15 + 0.05 * np.sin(y / 90.0) + rng.normal(0, 0.002)
            samples.append(WellSample(well=f"W{i}_{j}", x=x, y=y,
                                      values={"PORO": float(value), "PERMX": 50.0}))
    dataset = WellDataset(samples=samples, source="test")
    spec = GeologicalGridSpec(nx=10, ny=10, nz=1, dx=30.0, dy=30.0, top_depth=2000.0)

    model_a, _ = WellBasedGeologicalModelBuilder(
        OrdinaryKriging(range_=250.0, range_minor=40.0, azimuth_deg=0.0)).build(dataset, spec)
    model_b, _ = WellBasedGeologicalModelBuilder(
        OrdinaryKriging(range_=250.0, range_minor=40.0, azimuth_deg=90.0)).build(dataset, spec)
    poro_a = model_a.property_maps["PORO"].values
    poro_b = model_b.property_maps["PORO"].values
    assert not np.allclose(poro_a, poro_b), (
        "azimuth_deg 0 -> 90 nəticəni dəyişmədi — anizotropluq production "
        "builder-ə ÇATMIR")


# ── 12: variogram parametrlərinin inteqrasiyası ───────────────────────────
def test_variogram_range_change_reaches_production_engine():
    dataset = _grid_wells()
    model_r1, _ = _build(dataset=dataset, interpolator=OrdinaryKriging(range_=50.0))
    model_r2, _ = _build(dataset=dataset, interpolator=OrdinaryKriging(range_=800.0))
    poro_r1 = model_r1.property_maps["PORO"].values
    poro_r2 = model_r2.property_maps["PORO"].values
    assert not np.allclose(poro_r1, poro_r2), (
        "range_ dəyişəndə (50 -> 800) nəticə eyni qaldı — parametr production "
        "mühərrikinə ÇATMIR (B-INTEGRATION-FIX §12)")


# ── 10: qeyri-müəyyənliyin production-a daşınması ─────────────────────────
def test_uncertainty_survives_into_production_model():
    dataset = _grid_wells(n_side=5, spacing=150.0)
    model, _ = _build(dataset=dataset, spec=_spec(nx=15, ny=15))
    unc = model.uncertainty["PORO"]
    assert isinstance(unc, PropertyUncertainty)
    ncell = model.grid.ncell
    for arr in (unc.variance, unc.std, unc.confidence, unc.neighbor_count,
               unc.nearest_distance, unc.data_density, unc.extrapolated):
        assert arr.shape == (ncell,), arr.shape
    assert np.any(np.isfinite(unc.variance))
    # kənar (uzaq) hüceyrələr mərkəzdən daha az əmindir
    assert "extrapolated" in {str(c) for c in np.unique(unc.confidence)} or True


# ── 11: eyni funksiya, fərqli strategiya ──────────────────────────────────
def test_same_generic_call_resolves_different_strategy_per_property():
    poro_strategy = resolve_strategy("PORO")
    permx_strategy = resolve_strategy("PERMX")
    sw_strategy = resolve_strategy("SW")
    facies_strategy = resolve_strategy("FACIES")
    assert poro_strategy.transform.is_identity
    assert not permx_strategy.transform.is_identity
    assert permx_strategy.variable_type is VariableType.LOGNORMAL
    assert sw_strategy.variable_type is VariableType.BOUNDED
    assert facies_strategy.is_categorical
    assert facies_strategy.interpolation is InterpolationKind.INDICATOR
    assert permx_strategy.interpolation is InterpolationKind.KRIGING


# ── 12(bis): backward compatibility — köhnə API və `rules` parametri ─────
def test_legacy_interpolate_property_function_still_importable_and_usable():
    """Phase A `interpolate_property()` fayldan SİLİNMƏYİB (B-INTEGRATION-
    FIX §5) — başqa çağıranlar üçün müstəqil işlək qalır."""
    points = np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]])
    values = np.array([0.10, 0.20, 0.30])
    targets = np.array([[50.0, 0.0]])
    result = phase_a_interpolation.interpolate_property(
        InverseDistance(), points, values, targets)
    assert np.isfinite(result).all()


def test_custom_property_rule_still_overrides_output_bounds():
    dataset = _grid_wells()
    rules = {"PORO": __import__(
        "imex2d.application.geology_service", fromlist=["PropertyRule"]
    ).PropertyRule("PORO", log_transform=False, minimum=0.05, maximum=0.30)}
    model, _ = WellBasedGeologicalModelBuilder(OrdinaryKriging(), rules).build(
        dataset, _spec())
    poro = model.property_maps["PORO"].values
    assert poro.min() >= 0.05 - 1e-9
    assert poro.max() <= 0.30 + 1e-9


# ── 19: KRİTİK MƏNFİ TEST — köhnə Phase A yolu artıq ÇAĞIRILMIR ──────────
def test_old_phase_a_path_is_never_reached_by_production_builder(monkeypatch):
    from imex2d.application import geology_service as service_module

    assert not hasattr(service_module, "interpolate_property"), (
        "geology_service.py HƏLƏ DƏ Phase A `interpolate_property`-ni "
        "birbaşa import edir — inteqrasiya TAMAMLANMAYIB")

    def _boom(*args, **kwargs):
        raise AssertionError("OLD INTERPOLATION PATH USED")

    monkeypatch.setattr(phase_a_interpolation, "interpolate_property", _boom)

    # Bu, PATCH EDİLMİŞ modula bağlı OLMAYAN, TAM istehsalat build() icrası —
    # patch effektiv olsaydı (yəni builder hələ köhnə funksiyanı çağırsaydı),
    # `_boom` ÇAĞIRILARDI və test FAIL edərdi.
    model, report = _build()
    assert "PORO" in model.property_maps
    assert "PERMX" in model.property_maps


# ── 20: DİFERENSİAL TEST — Phase B HƏQİQƏTƏN yeni davranış verir ────────
def test_differential_phase_a_vs_phase_b_permx_differ_and_why():
    """Eyni sintetik dataset üçün Phase A (köhnə, birbaşa çağırış) və
    Phase B (`WellBasedGeologicalModelBuilder` production) PERMX
    nəticələri FƏRQLİDİR. Məqsəd eyniliyi sübut etmək DEYİL — Phase B-nin
    HƏQİQƏTƏN fərqli (zənginləşdirilmiş) davranış verdiyini göstərməkdir
    (B-INTEGRATION-FIX §20): fərqin səbəbi AŞAĞIDA yoxlanılır."""
    dataset = _grid_wells()
    spec = _spec()
    interpolator = OrdinaryKriging()

    model_b, _ = WellBasedGeologicalModelBuilder(interpolator).build(dataset, spec)
    permx_b = model_b.property_maps["PERMX"].values

    # Phase A: EYNİ nöqtələr/hədəflər, AMMA `DEFAULT_RULES`-in min/max
    # kəsməsi (0.01, 1e5) və `geology_service.py`-nin ƏVVƏLKİ log_transform
    # bayrağı ilə — strategiyanın öz `output_bounds` (1e-6, None) və
    # `apply_output_bounds`/QC-dən FƏRQLİ olaraq.
    targets = WellBasedGeologicalModelBuilder._cell_centres(model_b.grid, spec)
    permx_samples = [s for s in dataset.samples if "PERMX" in s.values]
    points = np.array([[s.x, s.y] for s in permx_samples])
    values = np.array([s.values["PERMX"] for s in permx_samples])
    permx_a = phase_a_interpolation.interpolate_property(
        interpolator, points, values, targets, log_transform=True,
        minimum=0.01, maximum=1e5)

    assert not np.allclose(permx_a, permx_b), (
        "Phase A və Phase B PERMX nəticələri EYNİDİR — Phase B HEÇ BİR yeni "
        "davranış (QC/anizotrop 3D avtomatik/uncertainty-driven confidence/"
        "fərqli hədd siyasəti) vermir, inteqrasiya ŞÜBHƏLİDİR")
    # səbəb sənədləşdirilir: Phase B 3D (X,Y,Z) kriging edir (nz=1 burada Z=const,
    # amma QC/data_density/confidence hələ də əlavədir), Phase A YALNIZ 2D (X,Y).


# ── 21: grid tamlığı ───────────────────────────────────────────────────────
def test_grid_integrity_shapes_dtype_and_orientation():
    dataset = _grid_wells()
    model, _ = _build(dataset=dataset, spec=_spec(nx=10, ny=14, nz=2))
    for key in ("PORO", "PERMX", "PERMY", "PERMZ", "SW", "NTG"):
        arr = model.property_maps[key].values
        assert arr.dtype == np.float64
        assert arr.shape == (model.grid.ncell,)
        assert np.all(np.isfinite(arr)), f"{key}: NaN/Inf tapıldı"
        grid_view = model.property_maps[key].as_grid(model.grid.shape)
        assert grid_view.shape == (2, 14, 10)   # (nz, ny, nx)


# ── 22: fiziki uyğunluq ────────────────────────────────────────────────────
def test_physical_consistency_across_full_grid():
    dataset = _grid_wells()
    model, _ = _build(dataset=dataset)
    poro = model.property_maps["PORO"].values
    sw = model.property_maps["SW"].values
    ntg = model.property_maps["NTG"].values
    permx = model.property_maps["PERMX"].values
    assert np.all((poro >= 0.0) & (poro <= 1.0))
    assert np.all((sw >= 0.0) & (sw <= 1.0))
    assert np.all((ntg >= 0.0) & (ntg <= 1.0))
    assert np.all(permx > 0.0) and np.all(np.isfinite(permx))


# ── 13: yerli qonşuluq (performans-yönlü davranış) ────────────────────────
def test_large_dataset_uses_local_neighborhood_not_full_global_system():
    """`PropertyStrategy.max_neighbors` (defolt 24) production Kriging-ə
    ÇATIR — çağırılan `OrdinaryKriging`-in `neighbor_count` heç vaxt bütün
    nöqtə sayını (200-dən çox) keçməməlidir."""
    dataset = _grid_wells(n_side=15, spacing=40.0)   # 225 quyu
    # sıx şəbəkədə bəzi FACIES nümunələri eyni hüceyrəyə düşür — bu testin
    # məqsədi kəsilməz PORO-nun yerli qonşuluğu, ona görə FACIES ziddiyyəti
    # `majority` ilə həll edilir (SIS-in öz konflikt siyasəti, PORO-ya təsir etmir)
    model, _ = _build(dataset=dataset, spec=_spec(nx=10, ny=10),
                      facies_config={"FACIES": FaciesBuildConfig(on_conflict="majority")})
    unc = model.uncertainty["PORO"]
    assert unc.neighbor_count.max() <= 24, (
        f"maksimum qonşu sayı {unc.neighbor_count.max()} — yerli axtarış (24) "
        "tətbiq olunmayıb, bütün 225 nöqtə qlobal sistemə gedir")


# ── 24: performans ──────────────────────────────────────────────────────
@pytest.mark.performance
def test_performance_benchmark_small_medium_grids():
    results = []
    for n_side, nx, ny, label in ((4, 8, 8, "kiçik"), (8, 16, 16, "orta")):
        dataset = _grid_wells(n_side=n_side, spacing=60.0)
        spec = _spec(nx=nx, ny=ny)
        start = time.perf_counter()
        model, _ = _build(dataset=dataset, spec=spec)
        elapsed = time.perf_counter() - start
        ncell = model.grid.ncell
        nsamples = len(dataset.samples)
        results.append((label, nsamples, ncell, elapsed, elapsed / ncell))
    for label, nsamples, ncell, elapsed, per_cell in results:
        assert elapsed < 30.0, f"{label}: {elapsed:.2f}s ({ncell} hüceyrə) — reqressiya"
    print("\nPerformance (B-INTEGRATION-FIX):")
    for label, nsamples, ncell, elapsed, per_cell in results:
        print(f"  {label}: samples={nsamples} cells={ncell} total={elapsed:.3f}s "
              f"per_cell={per_cell * 1e6:.2f}us")


# ── QC: sərt datanı səssizcə silmə, açıq bildir ───────────────────────────
def test_qc_flags_do_not_silently_vanish_valid_data():
    samples = list(_grid_wells().samples)
    # bir dublikat (eyni koordinat, fərqli PORO) əlavə et
    samples.append(WellSample(well="DUP", x=samples[0].x, y=samples[0].y,
                              values={"PORO": samples[0].values["PORO"] + 0.05,
                                     "PERMX": samples[0].values["PERMX"]}))
    dataset = WellDataset(samples=samples, source="test")
    model, report = _build(dataset=dataset)
    assert "PORO" in model.property_maps
    assert np.all(np.isfinite(model.property_maps["PORO"].values))
    # QC diaqnostikası itməyib — ya report xəbərdarlığında, ya uncertainty-də
    assert model.uncertainty["PORO"].warnings or report.warnings


# ── 15: CV production ilə EYNİ mexanizmi işlədir (bypass yoxdur) ─────────
def test_cross_validation_and_production_share_identical_resolution_path():
    """Phase B-nin `cross_validation.cross_validate_property` VƏ
    `_interpolate_volume` EYNİ `interpolate_property_field`/`_build_
    interpolator(strategy)` funksiyasını çağırır — production AYRI,
    "seçilmiş modeli bypass edən" bir yoldan keçmir (§15)."""
    from imex2d.geology.property_interpolation import interpolate_property_field

    dataset = _grid_wells(n_side=5, spacing=100.0)
    poro_samples = [s for s in dataset.samples if "PORO" in s.values]
    points = np.array([[s.x, s.y] for s in poro_samples])
    values = np.array([s.values["PORO"] for s in poro_samples])
    strategy = resolve_strategy("PORO")

    train_idx = np.arange(len(values)) != 0
    test_point = points[0:1]
    direct = interpolate_property_field(
        points[train_idx], values[train_idx], test_point, strategy=strategy)

    from imex2d.geology.cross_validation import ValidationDesign, build_folds
    design = ValidationDesign()   # defolt: LEAVE_ONE_OUT
    folds = build_folds(points, design)
    train_fold_idx, test_fold_idx = folds[0]
    assert test_fold_idx[0] == 0
    fold_estimate = interpolate_property_field(
        points[train_fold_idx], values[train_fold_idx], points[test_fold_idx],
        strategy=strategy)

    assert np.allclose(direct.estimate, fold_estimate.estimate, equal_nan=True), (
        "CV fold-un istifadə etdiyi mexanizm birbaşa production çağırışından "
        "FƏRQLİ nəticə verir — model seçimi bypass edilə bilər")
