"""LAY-ÜZRƏ MƏLUMAT MÖVCUDLUĞU + İNTERPOLYASİYA + TAMAMLAMA +
QEYRİ-MÜƏYYƏNLİK — tapşırıq §22-nin TEST A–O siyahısı və §23 kənar halları.

FİZİKİ SSENARİ (tapşırıq §1): 5 laylı model, quyu məlumatı YALNIZ
L1–L3 üçün var. Quyunun fiziki intervalı (`top`/`bottom`) beş layı KƏSSƏ
DƏ, bu, "beş layda petrofiziki ölçmə var" DEMƏK DEYİL.

Bu faylın hər testi bir ELMİ QAYDANI qoruyur; heç biri "UI-da gizlət"
tipli həlli qəbul etmir (bax §26).
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.application.geology_adapter import (well_layer_summary, wells_to_dataset)
from imex2d.application.geology_service import (CompletionMethod, CompletionSpec,
                                                ContinuousSGSConfig,
                                                GeologicalGridSpec,
                                                LayerInterpolationConfig,
                                                WellBasedGeologicalModelBuilder,
                                                compute_property_impact,
                                                format_cross_validation_report)
from imex2d.domain.data_availability import (DataStatus, PropertyAvailability,
                                             format_layers, parse_layers,
                                             parse_property_layers)
from imex2d.domain.geology import GeologicalWell, well_effective_layers
from imex2d.domain.geometry import CellGeometry, interval_layers
from imex2d.domain.grid import CartesianGrid
from imex2d.domain.well_data import WellDataset, WellSample
from imex2d.geology.interpolation import OrdinaryKriging
from imex2d.geology.layer_availability import (LayerDataPolicy, compute_availability,
                                               hard_data_cells, well_interval_layers)

# ── ortaq quraşdırma ─────────────────────────────────────────────────
_WELL_XY = {"A": (25.0, 25.0), "B": (125.0, 25.0), "C": (25.0, 125.0)}
_FACTOR = {"A": 1.00, "B": 1.06, "C": 0.94}
_PERMX_BY_LAYER = {0: 500.0, 1: 100.0, 2: 10.0}
_PORO_BY_LAYER = {0: 0.25, 1: 0.20, 2: 0.15}
_TOP_DEPTH = 2000.0
_DZ = 10.0


def _geometry(nz: int = 5, nx: int = 3, ny: int = 3) -> CellGeometry:
    return CellGeometry(CartesianGrid(nx, ny, nz), 50.0, 50.0, _DZ,
                        top_depth=_TOP_DEPTH)


def _spec(nz: int = 5, nx: int = 3, ny: int = 3) -> GeologicalGridSpec:
    return GeologicalGridSpec(nx=nx, ny=ny, nz=nz, dx=50.0, dy=50.0, dz=_DZ,
                              top_depth=_TOP_DEPTH)


def _dataset(data_layers=(0, 1, 2), properties=("PERMX", "PORO")) -> WellDataset:
    """3 quyu × verilmiş laylar — QALAN laylarda HEÇ BİR məlumat yoxdur."""
    samples = []
    for name, (x, y) in _WELL_XY.items():
        factor = _FACTOR[name]
        for k in data_layers:
            values = {}
            if "PERMX" in properties:
                values["PERMX"] = _PERMX_BY_LAYER.get(k, 50.0) * factor
            if "PORO" in properties:
                values["PORO"] = _PORO_BY_LAYER.get(k, 0.18) * factor
            samples.append(WellSample(well=name, x=x, y=y, layer=k, values=values))
    return WellDataset(samples=samples, source="test")


def _build(dataset: WellDataset, nz: int = 5, config=None, interpolator=None,
           **kwargs):
    builder = WellBasedGeologicalModelBuilder(interpolator or OrdinaryKriging())
    return builder.build(dataset, _spec(nz),
                         layer_config=config if config is not None
                         else LayerInterpolationConfig(), **kwargs)


def _grid_of(model, name: str) -> np.ndarray:
    return model.property_maps[name].as_grid(model.grid.shape)


def _status_grid(model, name: str) -> np.ndarray:
    return np.asarray(model.provenance[name].status,
                      dtype=object).reshape(model.grid.shape)


def _original(nz: int = 5, nx: int = 3, ny: int = 3, value: float = 0.111):
    """Süni "mövcud geoloji prior" — hər layda FƏRQLİ sabit."""
    field = np.zeros((nz, ny, nx))
    for k in range(nz):
        field[k] = value + 0.01 * k
    return field.ravel()


# ═══════════════════════════════════════════════════════════ TEST A
def test_a_five_layers_data_only_in_first_three():
    """5 lay, data yalnız L1–L3: L1–L3 interpolyasiya olunur, L4–L5
    TOXUNULMUR (MISSING qalır) — sükutla doldurulmur."""
    model, report = _build(_dataset(), nz=5)
    poro = _grid_of(model, "PORO")
    status = _status_grid(model, "PORO")

    for k in (0, 1, 2):
        assert np.all(np.isfinite(poro[k])), f"L{k + 1} hesablanmalı idi"
        assert set(status[k].ravel()) <= {DataStatus.INTERPOLATED.value,
                                          DataStatus.MEASURED.value}
    for k in (3, 4):
        assert np.all(np.isnan(poro[k])), f"L{k + 1} SƏSSİZCƏ doldurulub"
        assert set(status[k].ravel()) == {DataStatus.MISSING.value}
    assert model.availability["PORO"].data_layers() == [0, 1, 2]
    assert model.availability["PORO"].missing_layers() == [3, 4]
    assert report.has_blocking, "MISSING lay üçün bloklayan mesaj olmalıdır"


def test_a_missing_layers_are_not_a_copy_of_the_last_data_layer():
    """§26: L4/L5 nə L3-ün surəti, nə sıfır, nə də ixtiyari sabitdir —
    onlar MISSING-dir və bu, statusda AÇIQ görünür."""
    model, _ = _build(_dataset(), nz=5)
    poro = _grid_of(model, "PORO")
    status = _status_grid(model, "PORO")
    for k in (3, 4):
        assert np.isnan(poro[k]).all(), "MISSING lay ədədi dəyər daşımamalıdır"
        assert not np.allclose(np.nan_to_num(poro[k]), poro[2]), "L3-ün surəti"
        assert not np.allclose(np.nan_to_num(poro[k], nan=np.nan), 0.0), "sıfırla dolub"
        assert set(status[k].ravel()) == {DataStatus.MISSING.value}, (
            "problem NaN-ın arxasında GİZLƏDİLMƏMƏLİDİR — status MISSING olmalıdır")
    assert model.completeness_issues(), "MISSING vəziyyət validasiyada görünməlidir"


# ═══════════════════════════════════════════════════════════ TEST B
def test_b_selected_layers_one_to_three_do_not_overwrite_four_and_five():
    """İstifadəçi AÇIQ L1–L3 seçir: L4–L5 orijinal sahədən DƏYİŞMİR."""
    original = _original()
    config = LayerInterpolationConfig(target_layers=[0, 1, 2],
                                      original_fields={"PORO": original,
                                                       "PERMX": original * 1000.0})
    model, _ = _build(_dataset(), nz=5, config=config)
    poro = _grid_of(model, "PORO")
    base = original.reshape(model.grid.shape)
    status = _status_grid(model, "PORO")

    assert np.allclose(poro[3], base[3]), "L4 orijinaldan dəyişdirilib"
    assert np.allclose(poro[4], base[4]), "L5 orijinaldan dəyişdirilib"
    assert set(status[3].ravel()) == {DataStatus.PRESERVED.value}
    assert not np.allclose(poro[0], base[0]), "L1 interpolyasiya olunmalı idi"
    # mövcudluq cədvəli hüceyrə-səviyyəli provenance ilə ZİDD OLMAMALIDIR
    assert model.availability["PORO"].layers[3].status is DataStatus.PRESERVED
    assert model.missing_layers("PORO") == [], (
        "orijinal sahə saxlanılıbsa lay MISSING sayılmamalıdır")


# ═══════════════════════════════════════════════════════════ TEST C
def test_c_selected_layers_two_to_four_change_only_those_layers():
    """Seçim L2–L4, data L1–L3 → YALNIZ L2–L3 dəyişir (L4-də data yoxdur,
    ona görə interpolyasiya EDİLMİR); L1 və L5 TOXUNULMUR."""
    original = _original()
    config = LayerInterpolationConfig(target_layers=[1, 2, 3],
                                      original_fields={"PORO": original,
                                                       "PERMX": original * 1000.0})
    model, report = _build(_dataset(), nz=5, config=config)
    poro = _grid_of(model, "PORO")
    base = original.reshape(model.grid.shape)

    assert np.allclose(poro[0], base[0]), "L1 seçilməyib, dəyişməməli idi"
    assert np.allclose(poro[4], base[4]), "L5 seçilməyib, dəyişməməli idi"
    assert np.allclose(poro[3], base[3]), "L4-də data yoxdur — dəyişməməli idi"
    assert not np.allclose(poro[1], base[1])
    assert not np.allclose(poro[2], base[2])
    assert any("L4" in message or "data YOXDUR" in message
               for message in report.warnings)


def test_c_layer_selection_is_not_hard_coded_to_three_layers():
    """§6: hər NZ üçün işləyir — seçim 7 laylı gridin 5-ci layındadır."""
    dataset = _dataset(data_layers=(4,))
    config = LayerInterpolationConfig(target_layers=[4])
    model, _ = _build(dataset, nz=7, config=config)
    poro = _grid_of(model, "PORO")
    assert np.all(np.isfinite(poro[4]))
    for k in (0, 1, 2, 3, 5, 6):
        assert np.all(np.isnan(poro[k])), f"L{k + 1} toxunulmamalı idi"


# ═══════════════════════════════════════════════════════════ TEST D
def test_d_effective_interval_is_not_data_availability():
    """Quyu intervalı L1–L5-i kəsir, amma data yalnız L1–L3-dədir."""
    geometry = _geometry(nz=5)
    well = GeologicalWell(name="W-1", x=25.0, y=25.0, top=2000.0, bottom=2050.0,
                          porosity=0.25, permeability=500.0,
                          data_layers_text="1-3")
    effective = well_effective_layers(well, geometry)
    declared, _per_property = well.data_layer_sets(5)

    assert effective == [0, 1, 2, 3, 4], "interval beş layı kəsir"
    assert declared == [0, 1, 2], "məlumat yalnız üç laydadır"
    assert effective != declared, "iki anlayış EYNİ sayılmamalıdır"

    summary = well_layer_summary([well], geometry)
    assert summary["W-1"]["effective"] == [0, 1, 2, 3, 4]
    assert summary["W-1"]["data"] == [0, 1, 2]


def test_d_interval_does_not_leak_into_the_dataset_as_data():
    """`top/bottom` HEÇ VAXT lay indeksi kimi işlədilmir (§26)."""
    geometry = _geometry(nz=5)
    wells = [GeologicalWell(name=name, x=x, y=y, top=2000.0, bottom=2050.0,
                            porosity=0.25, permeability=500.0,
                            data_layers_text="1-3")
             for name, (x, y) in _WELL_XY.items()]
    dataset, _skipped = wells_to_dataset(wells, "Kriging (adi)", geometry,
                                         LayerDataPolicy.STRICT)
    availability = compute_availability(dataset, geometry, LayerDataPolicy.STRICT,
                                        ["PORO", "PERMX"])
    assert availability["PORO"].data_layers() == [0, 1, 2]
    assert availability["PERMX"].missing_layers() == [3, 4]


# ═══════════════════════════════════════════════════════════ TEST E
def test_e_without_completion_strategy_layers_stay_missing():
    """Tamamlama seçilməyibsə L4–L5 MISSING qalır və AÇIQ bloklayır."""
    model, report = _build(_dataset(), nz=5)
    assert model.missing_layers("PORO") == [3, 4]
    issues = model.completeness_issues()
    assert issues and all("MISSING" in issue for issue in issues)
    assert report.has_blocking
    # heç bir səssiz doldurma: nə 0, nə də sabit
    poro = _grid_of(model, "PORO")
    assert np.isnan(poro[3]).all()


def test_e_missing_model_is_rejected_before_the_simulator():
    """§20: natamam model rezervuar modelinə ÇEVRİLMİR."""
    from imex2d.application.model_builder import ReservoirModelBuilder
    model, _ = _build(_dataset(), nz=5)
    with pytest.raises(ValueError) as error:
        ReservoirModelBuilder().build(model, wells=[])
    assert "məlumat YOXDUR" in str(error.value)


# ═══════════════════════════════════════════════════════════ TEST F
def test_f_vertical_trend_marks_layers_estimated_with_method_metadata():
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.VERTICAL_TREND))
    model, report = _build(_dataset(), nz=5, config=config,
                           interpolator=OrdinaryKriging(range_v=15.0))
    poro = _grid_of(model, "PORO")
    status = _status_grid(model, "PORO")
    availability = model.availability["PORO"]

    assert np.all(np.isfinite(poro[3])) and np.all(np.isfinite(poro[4]))
    for k in (3, 4):
        assert set(status[k].ravel()) <= {DataStatus.ESTIMATED.value,
                                          DataStatus.EXTRAPOLATED.value}
        assert availability.layers[k].method == "vertical_trend"
        assert availability.layers[k].status is not DataStatus.MEASURED
    assert model.provenance["PORO"].layer_methods[3] == "vertical_trend"
    assert not report.has_blocking, "tamamlanmış model bloklanmamalıdır"


def test_f_vertical_trend_needs_at_least_two_data_layers():
    """Tək laydan "trend" qurmaq gizli ekstrapolyasiya olardı — rədd edilir."""
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.VERTICAL_TREND))
    with pytest.raises(ValueError) as error:
        _build(_dataset(data_layers=(0,)), nz=5, config=config)
    assert "ən azı İKİ məlumatlı lay" in str(error.value)


def test_f_vertical_trend_does_not_copy_the_neighbouring_layer_map():
    """§26: tamamlanan lay L3-ün LATERAL xəritəsinin surəti DEYİL."""
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.VERTICAL_TREND))
    model, _ = _build(_dataset(), nz=5, config=config)
    poro = _grid_of(model, "PORO")
    assert not np.allclose(poro[3], poro[2])
    assert np.allclose(poro[3], poro[3].ravel()[0]), (
        "şaquli trend lay daxilində SABİT dəyər verməlidir (lateral struktur "
        "uydurulmur)")


# ═══════════════════════════════════════════════════════════ TEST G
def test_g_sgs_completion_is_marked_simulated_with_realization_metadata():
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.SGS,
                                  sgs=ContinuousSGSConfig(seed=7, realization_id=2)))
    model, _ = _build(_dataset(), nz=5, config=config)
    status = _status_grid(model, "PORO")
    availability = model.availability["PORO"]

    for k in (3, 4):
        assert set(status[k].ravel()) == {DataStatus.SIMULATED.value}
        assert availability.layers[k].method == "sgs"
        assert "realization=2" in availability.layers[k].note
        assert "seed=7" in availability.layers[k].note
    assert DataStatus.MEASURED.value not in set(status[3].ravel())


def test_g_sgs_completion_leaves_data_layers_untouched():
    """SGS tamamlaması YALNIZ məlumatsız layları toxundurur."""
    plain, _ = _build(_dataset(), nz=5)
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.SGS))
    simulated, _ = _build(_dataset(), nz=5, config=config)
    for k in (0, 1, 2):
        assert np.allclose(_grid_of(plain, "PORO")[k], _grid_of(simulated, "PORO")[k])


# ═══════════════════════════════════════════════════════════ TEST H
def test_h_impact_analysis_does_not_modify_either_model():
    base, _ = _build(_dataset(), nz=3)
    hypothetical, _ = _build(_dataset(), nz=3,
                             config=LayerInterpolationConfig(
                                 completion=CompletionSpec(method=CompletionMethod.NONE)))
    # süni "fərz edək ki" ssenarisi: yalnız kopyada dəyişiklik
    hypothetical.property_maps["PORO"].values[:] *= 1.20

    before = np.array(base.property_maps["PORO"].values, copy=True)
    impact = compute_property_impact(base, hypothetical, "PORO")
    impact.delta[:] = 999.0                      # nəticəni korlamağa cəhd

    assert np.allclose(base.property_maps["PORO"].values, before), (
        "təsir analizi ƏSAS modeli dəyişdirdi")
    fresh = compute_property_impact(base, hypothetical, "PORO")
    assert np.allclose(fresh.delta, before * 0.20, rtol=1e-9)
    assert fresh.changed_cells == before.size
    assert "TƏSİR" in fresh.as_text()


# ═══════════════════════════════════════════════════════════ TEST I
def test_i_availability_is_property_specific():
    """PORO L1–L5, PERMX L1–L3 — hər xassə ÖZ mövcudluq xəritəsini işlədir."""
    geometry = _geometry(nz=5)
    wells = [GeologicalWell(name=name, x=x, y=y, top=2000.0, bottom=2050.0,
                            porosity=0.20 + 0.01 * index,
                            permeability=400.0 + 20.0 * index,
                            data_layers_text="PORO:1-5; PERMX:1-3")
             for index, (name, (x, y)) in enumerate(_WELL_XY.items())]
    dataset, _skipped = wells_to_dataset(wells, "Kriging (adi)", geometry,
                                         LayerDataPolicy.STRICT)
    model, _report = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(
        dataset, _spec(5), layer_config=LayerInterpolationConfig())

    assert model.availability["PORO"].data_layers() == [0, 1, 2, 3, 4]
    assert model.availability["PERMX"].data_layers() == [0, 1, 2]
    assert model.missing_layers("PORO") == []
    assert model.missing_layers("PERMX") == [3, 4]
    poro = _grid_of(model, "PORO")
    permx = _grid_of(model, "PERMX")
    assert np.all(np.isfinite(poro))
    assert np.isnan(permx[3]).all() and np.isnan(permx[4]).all()


def test_i_derived_permeability_inherits_permx_provenance():
    """PERMY/PERMZ PERMX-dən törəyir — MISSING statusu da miras qalır."""
    model, _ = _build(_dataset(), nz=5)
    for key, factor in (("PERMY", 1.0), ("PERMZ", 0.1)):
        assert key in model.provenance
        assert model.missing_layers(key) == [3, 4]
        methods = np.asarray(model.provenance[key].method, dtype=object)
        assert all(text.endswith(f"→{key}(×{factor:g})") for text in methods[:9]), (
            f"{key} mənşəyi PERMX-dən törədiyini AÇIQ göstərməlidir")
        assert model.availability[key].layers[3].status is DataStatus.MISSING


# ═══════════════════════════════════════════════════════════ TEST J
def test_j_hard_data_is_honoured_inside_the_layer_aware_path():
    """Quyu hüceyrəsində nəticə ÖLÇMƏ ilə üst-üstə düşür və MEASURED
    kimi işarələnir (interpolyasiya kimi YOX)."""
    geometry = _geometry(nz=3)
    dataset = _dataset(data_layers=(0, 1, 2))
    model, _ = _build(dataset, nz=3)
    poro = _grid_of(model, "PORO")
    status = _status_grid(model, "PORO")

    assert poro[0, 0, 0] == pytest.approx(_PORO_BY_LAYER[0] * _FACTOR["A"])
    assert poro[0, 0, 2] == pytest.approx(_PORO_BY_LAYER[0] * _FACTOR["B"])
    assert status[0, 0, 0] == DataStatus.MEASURED.value
    assert status[0, 1, 1] == DataStatus.INTERPOLATED.value

    mask = hard_data_cells(dataset, geometry, "PORO",
                           LayerDataPolicy.STRICT).reshape(model.grid.shape)
    assert mask[0, 0, 0] and not mask[0, 1, 1]


# ═══════════════════════════════════════════════════════════ TEST K
def test_k_samples_carry_true_cell_centre_depth():
    """3D Z koordinatı grid həndəsəsindən gəlir, uydurulmur."""
    geometry = _geometry(nz=5)
    depths = geometry.cell_depths().reshape(geometry.grid.shape)
    wells = [GeologicalWell(name=name, x=x, y=y, top=2000.0, bottom=2050.0,
                            porosity=0.25, permeability=500.0,
                            data_layers_text="1-3")
             for name, (x, y) in _WELL_XY.items()]
    dataset, _ = wells_to_dataset(wells, "Kriging (adi)", geometry,
                                  LayerDataPolicy.STRICT)
    layered = [s for s in dataset.samples if s.layer is not None]
    assert layered, "lay-məlumatlı rejimdə lay etiketli nümunə olmalıdır"
    for sample in layered:
        i = int(sample.x // geometry.dx)
        j = int(sample.y // geometry.dy)
        assert sample.depth == pytest.approx(float(depths[sample.layer, j, i]))


def test_k_variable_layer_thickness_changes_the_mapped_layer():
    """Dəyişən DZ ilə dərinlik→lay uyğunlaşdırması sürüşür."""
    geometry = CellGeometry(CartesianGrid(3, 3, 3), 50.0, 50.0, [5.0, 20.0, 5.0],
                            top_depth=2000.0)
    assert interval_layers(25.0, 25.0, 2000.0, 2004.0, geometry) == [0]
    assert interval_layers(25.0, 25.0, 2006.0, 2020.0, geometry) == [1]
    assert interval_layers(25.0, 25.0, 2004.0, 2026.0, geometry) == [0, 1, 2]


# ═══════════════════════════════════════════════════════════ TEST L
def test_l_anisotropic_kriging_still_affects_the_layer_aware_result():
    """Anizotropluq parametrləri lay maskası ilə BİRLİKDƏ işləyir."""
    dataset = _dataset(data_layers=(0, 1, 2))
    isotropic, _ = _build(dataset, nz=3, interpolator=OrdinaryKriging(range_=200.0))
    anisotropic, _ = _build(dataset, nz=3,
                            interpolator=OrdinaryKriging(range_=200.0,
                                                         range_minor=20.0,
                                                         azimuth_deg=45.0))
    assert not np.allclose(_grid_of(isotropic, "PORO"), _grid_of(anisotropic, "PORO")), (
        "anizotropluq lay-məlumatlı yolda təsirsiz qalıb"
    )


# ═══════════════════════════════════════════════════════════ TEST M
def test_m_categorical_property_still_goes_through_the_indicator_path():
    """Fasiya sütunu lay-məlumatlı rejimdə də kəsilməz Kriging-dən KEÇMİR."""
    samples = []
    for name, (x, y) in _WELL_XY.items():
        for k in (0, 1, 2):
            samples.append(WellSample(
                well=name, x=x, y=y, layer=k,
                values={"PORO": _PORO_BY_LAYER[k] * _FACTOR[name],
                        "PERMX": _PERMX_BY_LAYER[k] * _FACTOR[name],
                        "FACIES": float(1 if name == "A" else 2)}))
    dataset = WellDataset(samples=samples, source="test")
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.VERTICAL_TREND))
    model, _ = _build(dataset, nz=5, config=config)

    assert "FACIES" in model.facies_fields
    assert "FACIES" not in model.property_maps, (
        "kateqorik sütun kəsilməz xassə xəritəsinə çevrilməməlidir")
    codes = np.unique(model.facies_fields["FACIES"].codes)
    assert set(codes.tolist()) <= {1, 2}, "fasiya kodu ədədi olaraq interpolyasiya olunub"


# ═══════════════════════════════════════════════════════════ TEST N
def test_n_sgs_as_the_interpolation_engine_keeps_layer_masking():
    """SGS interpolyasiya mühərriki kimi seçiləndə də məlumatsız lay
    AVTOMATİK doldurulmur."""
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    model, _ = builder.build(
        _dataset(), _spec(5), layer_config=LayerInterpolationConfig(),
        sgs_config={"PORO": ContinuousSGSConfig(seed=3)})
    status = _status_grid(model, "PORO")
    assert set(status[0].ravel()) <= {DataStatus.SIMULATED.value,
                                      DataStatus.MEASURED.value}
    assert set(status[3].ravel()) == {DataStatus.MISSING.value}


# ═══════════════════════════════════════════════ ÇARPAZ-DOĞRULAMA (§19)
def test_cross_validation_reports_layers_without_validation_data():
    """L4–L5 üçün "RMSE = 0" kimi saxta nəticə YARADILMIR."""
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    results, skipped = builder.cross_validate(_dataset(), "PORO", nz=5)
    assert sorted(results) == [0, 1, 2]
    assert sorted(skipped) == [3, 4]
    for message in skipped.values():
        assert "doğrulama məlumatı yoxdur" in message
    text = format_cross_validation_report({"PORO": (results, skipped)})
    assert "K=3" in text and "K=4" in text


# ═══════════════════════════════════════════ ETİBARLILIQ / PROVENANCE
def test_confidence_is_absent_when_it_cannot_be_justified():
    """§18: şaquli korrelyasiya radiusu bilinmirsə SAXTA rəqəm yaradılmır."""
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.VERTICAL_TREND))
    model, _ = _build(_dataset(), nz=5, config=config,
                      interpolator=OrdinaryKriging())      # range_v verilməyib
    availability = model.availability["PORO"]
    assert availability.layers[3].confidence is None
    confidence = model.provenance["PORO"].confidence.reshape(model.grid.shape)
    assert np.isnan(confidence[3]).all()


def test_confidence_decreases_with_vertical_extrapolation_distance():
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.VERTICAL_TREND))
    model, _ = _build(_dataset(), nz=5, config=config,
                      interpolator=OrdinaryKriging(range_v=15.0))
    availability = model.availability["PORO"]
    assert availability.layers[3].confidence > availability.layers[4].confidence
    assert availability.layers[0].confidence > availability.layers[4].confidence


def test_provenance_keeps_original_interpolated_and_final_apart():
    """§11: sahələr AYRI saxlanılır, bir-birini əvəz etmir."""
    original = _original()
    config = LayerInterpolationConfig(
        target_layers=[0, 1, 2],
        original_fields={"PORO": original, "PERMX": original * 1000.0},
        completion=CompletionSpec(method=CompletionMethod.CONSTANT, value=0.17))
    model, _ = _build(_dataset(), nz=5, config=config)
    provenance = model.provenance["PORO"]
    shape = model.grid.shape

    assert np.allclose(provenance.original, original)
    assert np.isnan(provenance.interpolated.reshape(shape)[3]).all(), (
        "interpolyasiya sahəsi tamamlanmış layda dolu olmamalıdır")
    assert np.allclose(provenance.estimated.reshape(shape)[3], 0.17)
    assert np.isnan(provenance.estimated.reshape(shape)[0]).all()
    assert np.allclose(provenance.final.reshape(shape)[3], 0.17)
    assert provenance.confidence_kind == "ordinal_support_score"


def test_provenance_reaches_the_reservoir_model():
    """§20: mənşə simulyasiya modelinə çatanda İTMİR."""
    from imex2d.application.model_builder import ReservoirModelBuilder
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.CONSTANT, value=0.17))
    model, _ = _build(_dataset(), nz=5, config=config)
    reservoir = ReservoirModelBuilder().build(model, wells=[])
    assert "PORO" in reservoir.provenance
    status = np.asarray(reservoir.provenance["PORO"].status,
                        dtype=object).reshape(reservoir.grid.shape)
    assert set(status[3].ravel()) == {DataStatus.ESTIMATED.value}


# ══════════════════════════════════════════════ GERİYƏ-UYĞUNLUQ (§25)
def test_default_mode_is_completely_unchanged():
    """`layer_config` verilmədikdə köhnə davranış (bütün laylara yayılma)."""
    wells = [GeologicalWell(name=name, x=x, y=y, porosity=0.25, permeability=500.0)
             for name, (x, y) in _WELL_XY.items()]
    dataset, _skipped = wells_to_dataset(wells, "Kriging (adi)")
    assert not dataset.is_layered()
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    model, report = builder.build(dataset, _spec(3))
    poro = _grid_of(model, "PORO")
    assert np.allclose(poro[0], poro[1]) and np.allclose(poro[1], poro[2])
    assert not model.provenance and model.availability is None
    assert not report.has_blocking


def test_broadcast_policy_keeps_the_old_dataset_shape():
    geometry = _geometry(nz=5)
    wells = [GeologicalWell(name=name, x=x, y=y, porosity=0.25, permeability=500.0,
                            data_layers_text="1-3")
             for name, (x, y) in _WELL_XY.items()]
    dataset, _ = wells_to_dataset(wells, "Kriging (adi)", geometry,
                                  LayerDataPolicy.BROADCAST)
    assert not dataset.is_layered(), (
        "BROADCAST siyasəti lay etiketi yaratmamalıdır (geriyə-uyğunluq)")


# ═════════════════════════════════════════════════ KƏNAR HALLAR (§23)
def test_edge_empty_layer_selection_is_rejected():
    config = LayerInterpolationConfig(target_layers=[])
    with pytest.raises(ValueError) as error:
        config.targets_for("PORO", 5)
    assert "BOŞDUR" in str(error.value)


def test_edge_layer_index_outside_the_grid_is_rejected():
    config = LayerInterpolationConfig(target_layers=[0, 7])
    with pytest.raises(ValueError) as error:
        config.targets_for("PORO", 5)
    assert "kənardadır" in str(error.value)
    with pytest.raises(ValueError):
        parse_layers("1-9", 5)


def test_edge_top_greater_than_bottom_is_rejected():
    geometry = _geometry(nz=5)
    with pytest.raises(ValueError):
        interval_layers(25.0, 25.0, 2050.0, 2000.0, geometry)
    well = GeologicalWell(name="W", x=25.0, y=25.0, top=2050.0, bottom=2000.0)
    assert well_effective_layers(well, geometry) == []


def test_edge_interval_completely_outside_the_grid():
    geometry = _geometry(nz=5)
    assert well_interval_layers(25.0, 25.0, 1900.0, 1950.0, geometry) == []
    assert well_interval_layers(25.0, 25.0, 2100.0, 2150.0, geometry) == []


def test_edge_interval_touching_a_single_layer_and_a_boundary():
    geometry = _geometry(nz=5)
    assert well_interval_layers(25.0, 25.0, 2010.0, 2020.0, geometry) == [1]
    # sərhəd ÜZƏRİNDƏ bitən interval qonşu layı KƏSMİR
    assert well_interval_layers(25.0, 25.0, 2000.0, 2010.0, geometry) == [0]


@pytest.mark.parametrize("nz", [1, 2, 5, 8])
def test_edge_various_nz_values(nz):
    """NZ = 1, 2, 5, >5 — heç birində sabit "3 lay" fərziyyəsi yoxdur."""
    data_layers = tuple(range(min(nz, 3)))
    model, _ = _build(_dataset(data_layers=data_layers), nz=nz)
    poro = _grid_of(model, "PORO")
    for k in data_layers:
        assert np.all(np.isfinite(poro[k]))
    for k in range(len(data_layers), nz):
        assert np.isnan(poro[k]).all()


def test_edge_data_in_a_single_layer_only():
    model, _ = _build(_dataset(data_layers=(2,)), nz=5)
    assert model.availability["PORO"].data_layers() == [2]
    assert model.missing_layers("PORO") == [0, 1, 3, 4]


def test_edge_data_in_every_layer_leaves_nothing_missing():
    model, _ = _build(_dataset(data_layers=(0, 1, 2, 3, 4)), nz=5)
    assert model.missing_layers("PORO") == []
    assert not model.completeness_issues()


def test_edge_no_layer_declaration_in_strict_mode_yields_no_data():
    """Bəyan yoxdursa STRICT rejimdə heç bir laya məlumat getmir və bu,
    AÇIQ xəbərdarlıqla bildirilir (səssiz yayılma YOXDUR)."""
    geometry = _geometry(nz=5)
    wells = [GeologicalWell(name=name, x=x, y=y, porosity=0.25, permeability=500.0)
             for name, (x, y) in _WELL_XY.items()]
    dataset, _ = wells_to_dataset(wells, "Kriging (adi)", geometry,
                                  LayerDataPolicy.STRICT)
    assert [s for s in dataset.samples if s.layer is not None] == []
    assert any("HEÇ BİR laya məlumat vermir" in w for w in dataset.warnings)


def test_edge_interval_policy_is_an_explicit_opt_in_with_a_warning():
    geometry = _geometry(nz=5)
    wells = [GeologicalWell(name=name, x=x, y=y, top=2000.0, bottom=2020.0,
                            porosity=0.25, permeability=500.0)
             for name, (x, y) in _WELL_XY.items()]
    dataset, _ = wells_to_dataset(wells, "Kriging (adi)", geometry,
                                  LayerDataPolicy.INTERVAL)
    availability = compute_availability(dataset, geometry, LayerDataPolicy.INTERVAL,
                                        ["PORO"])
    assert availability["PORO"].data_layers() == [0, 1]
    assert any("FƏRZİYYƏDİR" in w for w in dataset.warnings)


def test_edge_duplicate_samples_do_not_break_the_layer_mask():
    dataset = _dataset()
    dataset.samples.extend([WellSample(well="A", x=25.0, y=25.0, layer=0,
                                       values={"PORO": 0.25, "PERMX": 500.0})])
    model, _ = _build(dataset, nz=5)
    assert model.availability["PORO"].layers[0].n_data == 4
    assert np.all(np.isfinite(_grid_of(model, "PORO")[0]))


def test_edge_nan_values_are_not_counted_as_data():
    samples = []
    for name, (x, y) in _WELL_XY.items():
        samples.append(WellSample(well=name, x=x, y=y, layer=0,
                                  values={"PORO": _PORO_BY_LAYER[0] * _FACTOR[name],
                                          "PERMX": _PERMX_BY_LAYER[0] * _FACTOR[name]}))
        samples.append(WellSample(well=name, x=x, y=y, layer=1,
                                  values={"PORO": float("nan"), "PERMX": float("nan")}))
    geometry = _geometry(nz=3)
    dataset = WellDataset(samples=samples, source="test")
    availability = compute_availability(dataset, geometry, LayerDataPolicy.STRICT,
                                        ["PORO", "PERMX"])
    assert availability["PORO"].data_layers() == [0], (
        "NaN dəyər 'məlumat var' kimi sayılmamalıdır")


def test_edge_constant_completion_requires_a_value():
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.CONSTANT))
    with pytest.raises(ValueError) as error:
        _build(_dataset(), nz=5, config=config)
    assert "value" in str(error.value)


def test_edge_preserve_original_requires_an_original_field():
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.PRESERVE_ORIGINAL))
    with pytest.raises(ValueError) as error:
        _build(_dataset(), nz=5, config=config)
    assert "original_fields" in str(error.value)


def test_edge_original_field_size_must_match_the_grid():
    config = LayerInterpolationConfig(original_fields={"PORO": np.zeros(7)})
    with pytest.raises(ValueError) as error:
        _build(_dataset(), nz=5, config=config)
    assert "uyğun gəlmir" in str(error.value)


# ══════════════════════════════════════════════════ MƏTN FORMATLARI
def test_layer_text_round_trip():
    assert format_layers([0, 1, 2, 4]) == "1-3,5"
    assert format_layers([]) == "—"
    assert parse_layers("1-3,5", 5) == [0, 1, 2, 4]
    assert parse_layers("*", 3) == [0, 1, 2]
    assert parse_property_layers("PORO:1-5; PERMX:1-3", 5) == (
        None, {"PORO": [0, 1, 2, 3, 4], "PERMX": [0, 1, 2]})
    assert parse_property_layers("1-3", 5) == ([0, 1, 2], {})


def test_availability_rejects_layer_index_outside_the_grid():
    availability = PropertyAvailability(name="PORO", nz=3)
    with pytest.raises(ValueError):
        availability.set(5, status=DataStatus.MEASURED)


# ═════════════════════════════════════════════════ PERFORMANS (§24)
def test_only_selected_layers_are_actually_kriged(monkeypatch):
    """20 laylı gridin yalnız 3 layında data var → Kriging mühərriki
    CƏMİ 3 dəfə (xassə başına) çağırılır, 20 dəfə YOX."""
    from imex2d.application import geology_service

    calls = []
    original = geology_service.interpolate_property_field

    def spy(points, values, target_points, **kwargs):
        calls.append(int(np.asarray(target_points).shape[0]))
        return original(points, values, target_points, **kwargs)

    monkeypatch.setattr(geology_service, "interpolate_property_field", spy)
    model, _ = _build(_dataset(data_layers=(0, 1, 2)), nz=20)
    areal = model.grid.nx * model.grid.ny

    assert len(calls) == 6, (
        f"iki kəsilməz xassə × üç lay = 6 çağırış gözlənilirdi, oldu {len(calls)}")
    assert all(count == areal for count in calls), (
        "hər çağırış YALNIZ bir layın hüceyrələrini hədəfləməlidir")


def test_sgs_completion_only_targets_missing_layer_cells(monkeypatch):
    """SGS tamamlaması bütün həcmi deyil, YALNIZ məlumatsız layları
    simulyasiya edir (§24)."""
    from imex2d.application import geology_service

    sizes = []
    original = geology_service.simulate_sgs

    def spy(points, values, targets, **kwargs):
        sizes.append(int(np.asarray(targets).shape[0]))
        return original(points, values, targets, **kwargs)

    monkeypatch.setattr(geology_service, "simulate_sgs", spy)
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.SGS))
    model, _ = _build(_dataset(), nz=5, config=config)
    areal = model.grid.nx * model.grid.ny

    assert sizes, "SGS çağırılmalı idi"
    assert all(size == 2 * areal for size in sizes), (
        f"yalnız iki məlumatsız layın hüceyrələri gözlənilirdi, oldu {sizes}")


# ═══════════════════════════════════════════════════ SERİALİZASİYA
def test_data_layer_declaration_survives_save_and_load(tmp_path):
    from imex2d.application.project import Project
    from imex2d.application.serialization import ProjectSerializer

    project = Project(name="lay-məlumatlı")
    project.geology_wells = [
        GeologicalWell(name="W-1", x=25.0, y=25.0, top=2000.0, bottom=2050.0,
                       porosity=0.25, permeability=500.0,
                       data_layers_text="PORO:1-5; PERMX:1-3")]
    path = tmp_path / "layers.imx"
    serializer = ProjectSerializer()
    serializer.save(project, str(path))
    loaded = serializer.load(str(path))

    assert loaded.geology_wells[0].data_layers_text == "PORO:1-5; PERMX:1-3"
    assert loaded.geology_wells[0].data_layer_sets(5) == (
        None, {"PORO": [0, 1, 2, 3, 4], "PERMX": [0, 1, 2]})


def test_old_project_without_layer_declaration_still_loads():
    """v1/v2 faylında bu açar YOXDUR — boş mətn = "bəyan edilməyib"."""
    well = GeologicalWell.from_dict({"name": "W", "x": 1.0, "y": 2.0})
    assert well.data_layers_text == ""
    assert well.data_layer_sets(5) == (None, {})


# ══════════════════ SİMULYATOR / MPFA-TPFA İNTEQRASİYASI (§20/§21)
def _completed_reservoir_model(nz: int = 5):
    from helpers import default_scal
    from imex2d.application.model_builder import ReservoirModelBuilder
    from imex2d.domain.initial import InitialConditions
    from imex2d.domain.wells import (ControlMode, Perforation, Well, WellControl,
                                     WellType)
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.VERTICAL_TREND))
    model, _ = _build(_dataset(), nz=nz, config=config,
                      interpolator=OrdinaryKriging(range_v=15.0))
    scal = default_scal()
    wells = [
        Well("INJ", WellType.INJECTOR, WellControl(ControlMode.RATE, 20.0),
             [Perforation(0, 0, k) for k in range(nz)]),
        Well("PROD", WellType.PRODUCER, WellControl(ControlMode.BHP, 190.0),
             [Perforation(2, 2, k) for k in range(nz)])]
    reservoir = ReservoirModelBuilder().build(
        model, wells, scal=scal,
        initial=InitialConditions(datum_pressure=200.0, water_saturation=scal.swc),
        name="lay-məlumatlı tamamlanmış model")
    return model, reservoir


def test_completed_model_reaches_the_tpfa_and_mpfa_discretisations():
    """Tamamlanmış (ESTIMATED) sahə simulyatora ÇATIR və hər iki axın
    diskretizasiyası onu qəbul edir — MPFA/TPFA nüvəsinə TOXUNULMAYIB."""
    from imex2d.discretization import MPFAOBoundaryClosure, MPFAODiscretization
    from imex2d.simulation.discretization import TwoPointFluxDiscretization
    _model, reservoir = _completed_reservoir_model(nz=5)

    assert np.all(np.isfinite(reservoir.rock.porosity.values))
    assert np.all(np.isfinite(reservoir.rock.permx.values))

    tpfa = TwoPointFluxDiscretization().build(reservoir)
    assert np.all(np.isfinite(tpfa.transmissibility))
    assert np.all(tpfa.transmissibility >= 0.0)

    # MPFA-O `transmissibility` skalyar sahəsi DEYİL, üz-əsaslı əmsallar
    # qurur — nüvəsinə TOXUNULMAYIB, burada yalnız tamamlanmış sahənin
    # ONA DA problemsiz çatdığı yoxlanılır.
    mpfa = MPFAODiscretization(closure=MPFAOBoundaryClosure.NEUMANN_ZERO).build(reservoir)
    assert np.all(np.isfinite(mpfa.pore_volume))
    assert mpfa.coefficients is not None
    assert not mpfa.unsupported_features, mpfa.unsupported_features


def test_completed_model_runs_an_end_to_end_simulation():
    """Ucdan-uca: quyu cədvəli -> lay maskası -> tamamlama -> simulyasiya."""
    from helpers import default_scal, make_service, short_config
    _model, reservoir = _completed_reservoir_model(nz=3)
    engine = make_service(default_scal()).create_engine(reservoir, short_config())
    result = engine.run()
    assert result.snapshots, "simulyasiya heç bir nəticə vermədi"
    assert np.all(np.isfinite(result.snapshots[-1].pressure))


def test_incomplete_model_never_reaches_the_flow_solver():
    """MISSING lay qalıbsa boru xətti SİMULYASİYADAN ƏVVƏL dayanır."""
    from imex2d.application.model_builder import ReservoirModelBuilder
    model, _ = _build(_dataset(), nz=5)          # tamamlama YOXDUR
    with pytest.raises(ValueError):
        ReservoirModelBuilder().build(model, wells=[])


# ═══════════════════════════════════════════════════════ UI (§6/§14)
@pytest.fixture(scope="module")
def qapp():
    QtWidgets = pytest.importorskip("PyQt5.QtWidgets")
    yield QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def test_ui_panel_reads_layer_selection_and_builds_a_config(qapp):
    """UI 1-əsaslıdır, mühərrik 0-əsaslı — çevirmə paneldə olur."""
    from imex2d.ui.panels import GeologyPanel
    panel = GeologyPanel()
    panel.set_geometry(_geometry(nz=5))
    assert panel.layer_config(5) is None, "rejim söndürülü ikən config olmamalıdır"

    panel.layer_aware.setChecked(True)
    panel.target_layers.setText("1-3")
    config = panel.layer_config(5)
    assert config.target_layers == [0, 1, 2], "1-əsaslı giriş 0-əsaslı olmalıdır"
    assert config.policy is LayerDataPolicy.STRICT
    assert config.completion.method is CompletionMethod.NONE
    assert panel.cross_layer_fallback_allowed() is False

    panel.completion_method.setCurrentIndex(
        panel.completion_method.findData(CompletionMethod.CONSTANT.value))
    panel.completion_value.setValue(0.17)
    spec = panel.layer_config(5).completion
    assert spec.method is CompletionMethod.CONSTANT
    assert spec.value == pytest.approx(0.17)


def test_ui_panel_rejects_an_out_of_range_layer_selection(qapp):
    from imex2d.ui.panels import GeologyPanel
    panel = GeologyPanel()
    panel.layer_aware.setChecked(True)
    panel.target_layers.setText("1-9")
    with pytest.raises(ValueError):
        panel.layer_config(5)


def test_ui_table_round_trips_the_data_layer_column(qapp):
    """"Data layları" sütunu cədvəldən oxunur və "Kəsdiyi laylar"
    sütunu AYRICA, top/bottom-dan hesablanır (§14)."""
    from imex2d.ui.panels import GeologyPanel
    panel = GeologyPanel()
    panel.set_geometry(_geometry(nz=5))
    panel.add_row(GeologicalWell(name="W-1", x=25.0, y=25.0, top=2000.0,
                                 bottom=2050.0, porosity=0.25, permeability=500.0,
                                 data_layers_text="1-3"))
    wells = panel.wells()
    assert wells[0].data_layers_text == "1-3"
    effective_item = panel.table.item(0, GeologyPanel.COL_EFFECTIVE)
    assert effective_item.text() == "1-5", (
        "kəsdiyi laylar data laylarından FƏRQLİ göstərilməlidir")
    assert "kəsir 1-5" in panel.layer_summary.toPlainText()
    assert "data 1-3" in panel.layer_summary.toPlainText()


# ══════════════════════════════════ KATEQORİK MÖVCUDLUQ QEYDİ (§15/§17)
def test_categorical_layers_without_hard_data_are_recorded_as_simulated():
    samples = []
    for name, (x, y) in _WELL_XY.items():
        for k in (0, 1, 2):
            samples.append(WellSample(
                well=name, x=x, y=y, layer=k,
                values={"PORO": _PORO_BY_LAYER[k] * _FACTOR[name],
                        "PERMX": _PERMX_BY_LAYER[k] * _FACTOR[name],
                        "FACIES": float(1 if name == "A" else 2)}))
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.VERTICAL_TREND))
    model, report = _build(WellDataset(samples=samples, source="t"), nz=5,
                           config=config)
    facies = model.availability["FACIES"]
    assert facies.layers[3].status is DataStatus.SIMULATED
    assert facies.layers[3].method == "sis"
    assert facies.layers[0].status is DataStatus.SIMULATED, (
        "SIS nəticəsi şərtlənmiş olsa belə 'ölçülmüş' deyil")
    assert any("kateqorik" in message for message in report.warnings)


# ═══════════════════════════════════════ 3D GÖRÜNTÜ (§13) — data yolu
def test_provenance_is_available_to_the_renderer_as_numeric_codes():
    from imex2d.application.model_builder import ReservoirModelBuilder
    from imex2d.rendering import renderers as R
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.CONSTANT, value=0.17))
    model, _ = _build(_dataset(), nz=5, config=config)
    reservoir = ReservoirModelBuilder().build(model, wells=[])
    volume, _cmap, low, high = R.MapRenderer()._select_volume(
        reservoir, R.provenance_key("PORO"), None)
    assert volume.shape == reservoir.grid.shape
    assert low == 0.0 and high > 0.0
    assert len(set(volume[0].ravel().tolist())) >= 1
    assert volume[0, 0, 0] != volume[3, 0, 0], (
        "ölçülmüş və qiymətləndirilmiş lay eyni rəng koduna düşməməlidir")
    confidence, _cmap, _low, _high = R.MapRenderer()._select_volume(
        reservoir, R.confidence_key("PORO"), None)
    assert confidence.shape == reservoir.grid.shape


def test_original_and_impact_views_are_available_without_touching_final():
    """§12: UI Original / Final / Impact-i AYRI göstərə bilir və TƏSİR
    baxışı əsas sahəni DƏYİŞMİR."""
    from imex2d.application.model_builder import ReservoirModelBuilder
    from imex2d.rendering import renderers as R
    original = _original()
    config = LayerInterpolationConfig(
        original_fields={"PORO": original, "PERMX": original * 1000.0},
        completion=CompletionSpec(method=CompletionMethod.CONSTANT, value=0.17))
    model, _ = _build(_dataset(), nz=5, config=config)
    reservoir = ReservoirModelBuilder().build(model, wells=[])
    final_before = np.array(reservoir.rock.porosity.values, copy=True)

    base, _cmap, _low, _high = R.MapRenderer()._select_volume(
        reservoir, R.original_key("PORO"), None)
    impact, _cmap, low, high = R.MapRenderer()._select_volume(
        reservoir, R.impact_key("PORO"), None)

    assert np.allclose(base.ravel(), original)
    assert np.allclose(impact.ravel(),
                       reservoir.rock.porosity.values - original, equal_nan=True)
    assert low == -high and high > 0.0, "təsir şkalası simmetrik olmalıdır"
    assert np.allclose(reservoir.rock.porosity.values, final_before), (
        "təsir baxışı FİNAL sahəni dəyişdirdi")


def test_status_filter_hides_cells_with_a_different_provenance():
    """§13: "yalnız qiymətləndirilmiş hüceyrələri göstər" filtri."""
    from imex2d.application.model_builder import ReservoirModelBuilder
    from imex2d.rendering.volume import VolumeFilter
    config = LayerInterpolationConfig(
        completion=CompletionSpec(method=CompletionMethod.CONSTANT, value=0.17))
    model, _ = _build(_dataset(), nz=5, config=config)
    reservoir = ReservoirModelBuilder().build(model, wells=[])

    mask = reservoir.provenance["PORO"].mask(DataStatus.ESTIMATED.value)
    visible = VolumeFilter(cell_mask=mask).mask(
        reservoir.rock.porosity.values, reservoir.grid.shape)
    assert visible[3].all() and visible[4].all(), "tamamlanmış laylar görünməlidir"
    assert not visible[0].any(), "interpolyasiya olunmuş lay gizlədilməli idi"
    # filtrsiz hər şey görünür — filtr ƏSAS sahəni DƏYİŞMİR
    assert VolumeFilter().mask(reservoir.rock.porosity.values,
                               reservoir.grid.shape).all()
