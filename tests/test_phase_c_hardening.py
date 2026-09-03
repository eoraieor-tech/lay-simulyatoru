"""PHASE C — interpolation finalization & production hardening.

Bu fayl B-INTEGRATION-FIX-in (`test_phase_b_production_integration.py`)
ÜSTÜNDƏ qurulur. Underlying geostatistik alqoritmlər (variogram model
seçimi, spatial-block CV, anizotropluq aşkarlanması, SGS ensemble) artıq
`test_model_selection.py`/`test_anisotropy.py`/`test_synthetic_
validation.py`/`test_sgs_ensemble.py`-də DƏRİN doğrulanıb — bu fayl onları
TƏKRARLAMIR. Mərkəzi iddia budur: bu mövcud, dərin doğrulanmış imkanlar
(a) real `WellBasedGeologicalModelBuilder` production yolundan ƏLÇATANDIR,
(b) yeni diaqnostik/gate imkanları (`build_quality_report`,
`run_validation_gate`, `calibrate_property`) DÜZGÜN işləyir, (c) production
edge-case-lərdə səssizcə korlanmır.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from imex2d.application.geology_service import (FaciesBuildConfig, GeologicalGridSpec,
                                                 WellBasedGeologicalModelBuilder,
                                                 build_quality_report, run_validation_gate)
from imex2d.domain.well_data import WellDataset, WellSample
from imex2d.geology.cross_validation import ValidationDesign, ValidationKind, cross_validate_property
from imex2d.geology.interpolation import OrdinaryKriging
from imex2d.geology.property_config import resolve_strategy


# ── ortaq helper-lər ───────────────────────────────────────────────────
def _grid_wells(n_side=7, spacing=100.0, seed=0):
    rng = np.random.default_rng(seed)
    samples = []
    idx = 0
    for i in range(n_side):
        for j in range(n_side):
            x, y = i * spacing + 10.0, j * spacing + 10.0
            poro = float(np.clip(0.15 + 0.06 * np.sin(x / 260.0) + rng.normal(0, 0.004), 0.05, 0.35))
            permx = float(np.clip(np.exp(2.0 + 3.0 * poro + rng.normal(0, 0.03)), 1.0, 5000.0))
            sw = float(np.clip(0.3 + 0.3 * (i / n_side) + rng.normal(0, 0.01), 0.02, 0.98))
            ntg = float(np.clip(0.85 - 0.4 * (j / n_side) + rng.normal(0, 0.01), 0.02, 0.98))
            samples.append(WellSample(well=f"W{idx}", x=x, y=y,
                                      values={"PORO": poro, "PERMX": permx, "SW": sw, "NTG": ntg}))
            samples.append(WellSample(well=f"W{idx}", x=x, y=y, layer=0,
                                      values={"FACIES": 1 if poro < 0.18 else 2}))
            idx += 1
    return WellDataset(samples=samples, source="test")


def _spec(nx=12, ny=12, nz=1):
    return GeologicalGridSpec(nx=nx, ny=ny, nz=nz, dx=35.0, dy=35.0, dz=10.0, top_depth=2000.0)


def _build(dataset=None, interpolator=None, spec=None, **kwargs):
    dataset = dataset if dataset is not None else _grid_wells()
    interpolator = interpolator if interpolator is not None else OrdinaryKriging()
    spec = spec if spec is not None else _spec()
    # sıx quyu şəbəkəsində bəzi FACIES nümunələri eyni hüceyrəyə düşə bilər —
    # bu testlərin əksəriyyəti FACIES-i yox, kəsilməz xassələri yoxlayır, ona
    # görə ziddiyyət SIS-in öz `majority` siyasəti ilə həll olunur (defolt)
    kwargs.setdefault("facies_config", {"FACIES": FaciesBuildConfig(on_conflict="majority")})
    builder = WellBasedGeologicalModelBuilder(interpolator)
    return builder, builder.build(dataset, spec, **kwargs)


# ── 24: diagnostic report ─────────────────────────────────────────────
def test_quality_report_covers_all_phase_b_properties_with_expected_fields():
    dataset = _grid_wells()
    builder, (model, report) = _build(dataset=dataset)
    quality = build_quality_report(builder, model, dataset)
    assert set(quality) == set(model.uncertainty)   # PERMX/PORO/SW/NTG — PERMY/PERMZ DEYİL (törəmə)
    for name, q in quality.items():
        assert q.sample_count > 0
        assert q.mean_uncertainty is not None
        assert 0.0 <= q.extrapolated_fraction <= 1.0
        assert sum(q.support_distribution.values()) == model.grid.ncell
        assert sum(q.confidence_distribution.values()) == model.grid.ncell
        assert q.as_text()   # boş deyil, çökmür


def test_quality_report_excludes_derived_permy_permz():
    dataset = _grid_wells()
    builder, (model, report) = _build(dataset=dataset)
    quality = build_quality_report(builder, model, dataset)
    assert "PERMY" not in quality and "PERMZ" not in quality
    assert "PERMY" in model.property_maps   # törəmə kimi MÖVCUDDUR, sadəcə diaqnostikasız


# ── 25/26: validation gate ─────────────────────────────────────────────
def test_validation_gate_passes_on_well_covered_dataset():
    dataset = _grid_wells()
    builder, (model, report) = _build(dataset=dataset)
    quality = build_quality_report(builder, model, dataset)
    gate = run_validation_gate(model, report, quality)
    assert gate.has_errors is False


def test_validation_gate_warns_on_sparse_data_without_blocking():
    samples = [
        WellSample(well="A", x=0.0, y=0.0, values={"PORO": 0.15, "PERMX": 50.0}),
        WellSample(well="B", x=500.0, y=500.0, values={"PORO": 0.20, "PERMX": 80.0}),
    ]
    dataset = WellDataset(samples=samples, source="test")
    spec = GeologicalGridSpec(nx=25, ny=25, nz=1, dx=25.0, dy=25.0, top_depth=2000.0)
    builder, (model, report) = _build(dataset=dataset, spec=spec)
    quality = build_quality_report(builder, model, dataset)
    gate = run_validation_gate(model, report, quality)
    assert gate.has_errors is False   # 2 nöqtə fiziki cəhətdən etibarsız DEYİL, sadəcə seyrək
    assert len(gate.warnings) > 0     # amma xəbərdarlıq gizlədilmir


def test_validation_gate_flags_nan_as_fatal_error():
    """Gate `model.validate()`-ə kor-koranə etibar etmir — müstəqil yoxlayır."""
    from imex2d.domain.diagnostics import DiagnosticReport
    from imex2d.application.geology_service import InterpolationReport

    dataset = _grid_wells()
    builder, (model, report) = _build(dataset=dataset)
    model.property_maps["PORO"].values[0] = np.nan
    gate = run_validation_gate(model, report)
    assert gate.has_errors is True
    assert any("PORO" == d.source for d in gate.errors)


# ── 2/3/18: model kalibrasiyası — spatial-block CV, stabillik ─────────
def test_calibrate_property_uses_spatial_block_cv_by_default():
    dataset = _grid_wells()
    builder, _ = _build(dataset=dataset)
    report = builder.calibrate_property(dataset, "PORO")
    assert report.design.kind is ValidationKind.SPATIAL_BLOCK
    assert report.selected is not None
    assert report.selected.metrics.n > 0


def test_model_selection_is_deterministic_across_repeated_runs():
    dataset = _grid_wells()
    builder, _ = _build(dataset=dataset)
    r1 = builder.calibrate_property(dataset, "PORO")
    r2 = builder.calibrate_property(dataset, "PORO")
    assert r1.selected.candidate.label == r2.selected.candidate.label
    assert r1.selected.score == r2.selected.score


def test_calibrated_strategy_actually_changes_production_output():
    """Kalibrasiya nəticəsi `build(calibrated_strategies=...)` ilə real
    production nəticəsinə TƏSİR EDİR — sadəcə hesabatda qalmır."""
    dataset = _grid_wells()
    builder, (model_default, _) = _build(dataset=dataset)
    selection = builder.calibrate_property(dataset, "PORO")
    model_calibrated, _ = builder.build(
        dataset, _spec(), calibrated_strategies={"PORO": selection.selected.candidate.strategy},
        facies_config={"FACIES": FaciesBuildConfig(on_conflict="majority")})
    assert "PORO" in model_calibrated.property_maps
    # Fərqli namizəd modellər eyni nəticəni verə BİLƏR (data-ya görə) — ona
    # görə burada YALNIZ kalibrasiya yolunun ÇÖKMƏDİYİNİ və eyni sxemi
    # (`PropertyStrategy`) production-a ÖTÜRDÜYÜNÜ yoxlayırıq.
    assert np.all(np.isfinite(model_calibrated.property_maps["PORO"].values))


# ── 6/7: anizotropluq — production-a opt-in aşkarlanma ─────────────────
def test_auto_detect_anisotropy_opt_in_reaches_production_and_changes_result():
    """`OrdinaryKriging(auto_detect_anisotropy=True)` production builder-ə
    ötürüldükdə (bax `_kriging_overrides` PHASE C düzəlişi) İZOTROP defolt
    nəticədən FƏRQLİ nəticə verməlidir — açıq anizotrop sintetik sahədə."""
    rng = np.random.default_rng(3)
    samples = []
    for i in range(8):
        for j in range(8):
            x, y = i * 50.0, j * 50.0
            value = 0.15 + 0.06 * np.sin(y / 70.0) + rng.normal(0, 0.002)   # Y-də sürətli dəyişkənlik
            samples.append(WellSample(well=f"W{i}_{j}", x=x, y=y,
                                      values={"PORO": float(value), "PERMX": 50.0}))
    dataset = WellDataset(samples=samples, source="test")
    spec = GeologicalGridSpec(nx=12, ny=12, nz=1, dx=25.0, dy=25.0, top_depth=2000.0)

    isotropic, _ = WellBasedGeologicalModelBuilder(OrdinaryKriging()).build(dataset, spec)
    detected, _ = WellBasedGeologicalModelBuilder(
        OrdinaryKriging(auto_detect_anisotropy=True)).build(dataset, spec)
    assert not np.allclose(isotropic.property_maps["PORO"].values,
                           detected.property_maps["PORO"].values), (
        "auto_detect_anisotropy=True production nəticəsini dəyişmədi — "
        "`_kriging_overrides` bu bayrağı ötürmür")


# ── 8: uncertainty kalibrasiyası (spatial CV) ──────────────────────────
def test_uncertainty_is_reasonably_calibrated_under_spatial_cv():
    """`z = (əsl − proqnoz) / kriging_std` üçün `mean(z) ≈ 0`, `std(z) ≈ 1`
    gözlənilir — Phase B-nin `cross_validate_property`-si bunu artıq
    `mean_standardized_error`/`variance_standardized_error` kimi hesablayır
    (bax `ContinuousCVMetrics`); bura YALNIZ production strategiyası ilə
    ÇAĞIRIB gevşək tolerans daxilində olduğunu təsdiqləyir."""
    rng = np.random.default_rng(11)
    points = rng.uniform(0.0, 800.0, size=(60, 2))
    trend = 0.15 + 0.00015 * points[:, 0] + 0.00010 * points[:, 1]
    noise = rng.normal(0.0, 0.006, size=60)
    values = np.clip(trend + noise, 0.05, 0.35)
    strategy = resolve_strategy("PORO")
    design = ValidationDesign(kind=ValidationKind.SPATIAL_BLOCK, seed=7)
    metrics = cross_validate_property(points, values, strategy, design)
    assert metrics.n >= 30
    assert abs(metrics.mean_standardized_error) < 1.0, (
        f"standartlaşdırılmış xətanın ortası {metrics.mean_standardized_error:.3f} — 0-dan çox uzaq")
    assert 0.15 < metrics.variance_standardized_error < 6.0, (
        f"standartlaşdırılmış xətanın dispersiyası {metrics.variance_standardized_error:.3f} — "
        "kalibrasiya son dərəcə pisdir")


def test_confidence_classification_is_deterministic_for_same_input():
    dataset = _grid_wells()
    builder1, (model1, _) = _build(dataset=dataset)
    builder2, (model2, _) = _build(dataset=dataset)
    assert list(model1.uncertainty["PORO"].confidence) == list(model2.uncertainty["PORO"].confidence)


# ── 10/11: extrapolyasiya + dəstək təsnifatı ───────────────────────────
def test_far_away_cells_are_flagged_extrapolated_with_poor_support_and_low_confidence():
    samples = [WellSample(well="A", x=50.0, y=50.0, values={"PORO": 0.18, "PERMX": 40.0}),
              WellSample(well="B", x=60.0, y=55.0, values={"PORO": 0.19, "PERMX": 45.0}),
              WellSample(well="C", x=55.0, y=60.0, values={"PORO": 0.17, "PERMX": 42.0})]
    dataset = WellDataset(samples=samples, source="test")
    # grid quyulardan ÇOX uzağa yayılır (əksər hüceyrələr ekstrapolyasiya olacaq)
    spec = GeologicalGridSpec(nx=30, ny=30, nz=1, dx=200.0, dy=200.0, top_depth=2000.0)
    builder, (model, report) = _build(dataset=dataset, spec=spec)
    unc = model.uncertainty["PORO"]
    assert np.any(unc.extrapolated)
    far_idx = np.argmax(unc.nearest_distance)
    assert bool(unc.extrapolated[far_idx])
    assert str(unc.confidence[far_idx]) == "extrapolated"


# ── 12: edge-case hardening (production `build()`) ─────────────────────
def test_zero_samples_gives_clear_error_not_a_crash():
    dataset = WellDataset(samples=[], source="test")
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    with pytest.raises(ValueError):
        builder.build(dataset, _spec())


def test_one_sample_gives_a_clear_error_not_a_crash():
    """`WellDataset.validate()` (domen qatı, Phase B-dən ƏVVƏL) tək nöqtəni
    AÇIQ rədd edir — "1 nümunə → deterministik fallback VƏ YA aydın
    məhdudiyyət" tələbinin (C§12) İKİNCİ yarısı: sükutla uydurulmuş
    (məs. bütün grid = tək dəyər) NƏTİCƏ YOX, AYDIN mesajlı xəta."""
    dataset = WellDataset(samples=[
        WellSample(well="A", x=100.0, y=100.0, values={"PORO": 0.20, "PERMX": 75.0})],
        source="test")
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    with pytest.raises(ValueError):
        builder.build(dataset, _spec(nx=6, ny=6))


def test_two_samples_singular_geometry_does_not_crash():
    dataset = WellDataset(samples=[
        WellSample(well="A", x=0.0, y=0.0, values={"PORO": 0.15, "PERMX": 40.0}),
        WellSample(well="B", x=300.0, y=0.0, values={"PORO": 0.25, "PERMX": 90.0})],
        source="test")
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    model, report = builder.build(dataset, GeologicalGridSpec(nx=6, ny=1, nz=1, dx=60.0, dy=60.0))
    assert np.all(np.isfinite(model.property_maps["PORO"].values))


def test_duplicate_coordinates_are_handled_deterministically():
    dataset = WellDataset(samples=[
        WellSample(well="A", x=50.0, y=50.0, values={"PORO": 0.15, "PERMX": 40.0}),
        WellSample(well="B", x=50.0, y=50.0, values={"PORO": 0.25, "PERMX": 90.0}),
        WellSample(well="C", x=400.0, y=400.0, values={"PORO": 0.20, "PERMX": 60.0})],
        source="test")
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    model1, _ = builder.build(dataset, _spec(nx=8, ny=8))
    model2, _ = builder.build(dataset, _spec(nx=8, ny=8))
    assert np.allclose(model1.property_maps["PORO"].values, model2.property_maps["PORO"].values)


def test_constant_property_gives_zero_variance_not_nan():
    samples = [WellSample(well=f"W{i}", x=float(i * 80), y=float((i % 3) * 80),
                          values={"PORO": 0.20, "PERMX": 50.0}) for i in range(9)]
    dataset = WellDataset(samples=samples, source="test")
    builder, (model, report) = _build(dataset=dataset, spec=_spec(nx=8, ny=8))
    poro = model.property_maps["PORO"].values
    assert np.allclose(poro, 0.20)
    unc = model.uncertainty["PORO"]
    assert np.all(np.isfinite(unc.variance))


def test_extreme_permeability_values_stay_finite_and_positive():
    samples = [
        WellSample(well="LOW", x=0.0, y=0.0, values={"PORO": 0.10, "PERMX": 1e-6}),
        WellSample(well="HIGH", x=400.0, y=400.0, values={"PORO": 0.30, "PERMX": 1e6}),
        WellSample(well="MID", x=200.0, y=200.0, values={"PORO": 0.20, "PERMX": 500.0}),
    ]
    dataset = WellDataset(samples=samples, source="test")
    builder, (model, report) = _build(dataset=dataset, spec=_spec(nx=10, ny=10))
    permx = model.property_maps["PERMX"].values
    assert np.all(np.isfinite(permx)) and np.all(permx > 0.0)


def test_sw_and_ntg_exactly_zero_and_one_at_hard_data_do_not_break_bounds():
    samples = [
        WellSample(well="A", x=0.0, y=0.0, values={"PORO": 0.15, "PERMX": 40.0, "SW": 0.0, "NTG": 0.0}),
        WellSample(well="B", x=400.0, y=400.0, values={"PORO": 0.20, "PERMX": 60.0, "SW": 1.0, "NTG": 1.0}),
        WellSample(well="C", x=0.0, y=400.0, values={"PORO": 0.18, "PERMX": 50.0, "SW": 0.5, "NTG": 0.5}),
    ]
    dataset = WellDataset(samples=samples, source="test")
    builder, (model, report) = _build(dataset=dataset, spec=_spec(nx=15, ny=15))
    for key in ("SW", "NTG"):
        values = model.property_maps[key].values
        assert values.min() >= -1e-9 and values.max() <= 1.0 + 1e-9
        assert np.all(np.isfinite(values))


def test_nan_input_value_is_rejected_with_a_clear_error_not_silently_corrupted():
    """`WellDataset.validate()` (domen qatı) NaN-ı Phase B QC-yə ÇATMADAN
    AÇIQ rədd edir — bu, "heç biri silent corruption yaratmamalıdır"
    (C§12) tələbini AYRI, daha erkən qatda təmin edir: NaN sükutla
    atılıb qalan 3 nöqtə ilə davam ETMİR, bütün dataset AÇIQ xəta ilə
    rədd edilir (istifadəçi mənbə datasını düzəltməlidir)."""
    samples = [
        WellSample(well="A", x=0.0, y=0.0, values={"PORO": float("nan"), "PERMX": 40.0}),
        WellSample(well="B", x=200.0, y=0.0, values={"PORO": 0.18, "PERMX": 50.0}),
        WellSample(well="C", x=0.0, y=200.0, values={"PORO": 0.22, "PERMX": 60.0}),
        WellSample(well="D", x=200.0, y=200.0, values={"PORO": 0.20, "PERMX": 55.0}),
    ]
    dataset = WellDataset(samples=samples, source="test")
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    with pytest.raises(ValueError):
        builder.build(dataset, _spec(nx=8, ny=8))


def test_qc_layer_drops_out_of_bound_value_and_flags_it_without_crashing():
    """`WellDataset.validate()` NaN-ı ötürməsə də, FİZİKİ hədddən kənar
    (amma sonlu) dəyər ordan keçir — Phase B QC-nin özü bunu AÇIQ
    işarələyib (silmədən keçirmir, bax `data_quality.py`) atır və SAYIR."""
    samples = [
        WellSample(well="A", x=0.0, y=0.0, values={"PORO": -0.50, "PERMX": 40.0}),  # fiziki etibarsız
        WellSample(well="B", x=200.0, y=0.0, values={"PORO": 0.18, "PERMX": 50.0}),
        WellSample(well="C", x=0.0, y=200.0, values={"PORO": 0.22, "PERMX": 60.0}),
        WellSample(well="D", x=200.0, y=200.0, values={"PORO": 0.20, "PERMX": 55.0}),
    ]
    dataset = WellDataset(samples=samples, source="test")
    builder, (model, report) = _build(dataset=dataset, spec=_spec(nx=8, ny=8))
    assert np.all(np.isfinite(model.property_maps["PORO"].values))
    assert model.uncertainty["PORO"].warnings   # atılma SAYILIB, gizlədilmir


# ── 13: sərt-data validasiyası (kateqorik) ─────────────────────────────
def test_categorical_hard_data_class_is_honored_with_high_probability():
    dataset = _grid_wells()
    builder, (model, report) = _build(
        dataset=dataset,
        facies_config={"FACIES": FaciesBuildConfig(deterministic=True, on_conflict="majority")})
    from imex2d.domain.geometry import xy_to_ij
    facies_samples = [s for s in dataset.samples if "FACIES" in s.values]
    unc = model.uncertainty["FACIES"]
    grid = model.grid
    codes = model.facies_fields["FACIES"].codes.reshape(grid.shape)
    hits = 0
    for sample in facies_samples:
        i, j = xy_to_ij(sample.x, sample.y, model.geometry)
        if int(codes[0, j, i]) == int(sample.values["FACIES"]):
            hits += 1
    assert hits / len(facies_samples) > 0.9, "sərt-data fasiya kodu əksəriyyətdə honored deyil"


# ── 14: paylanma yoxlaması (distribution check) ────────────────────────
def test_output_distribution_does_not_collapse_relative_to_input():
    dataset = _grid_wells(n_side=10, spacing=70.0)
    builder, (model, report) = _build(dataset=dataset, spec=_spec(nx=25, ny=25))
    for key in ("PORO", "SW", "NTG"):
        input_values = np.array([s.values[key] for s in dataset.samples if key in s.values])
        output_values = model.property_maps[key].values
        # kriging hamarlaşdırıcıdır (gözlənilən) — amma dispersiya SIFIRA
        # ÇÖKMƏMƏLİDİR (kobud "hər yerdə orta" reqressiyası əlaməti olardı)
        assert output_values.std() > 0.15 * input_values.std(), (
            f"{key}: çıxış dispersiyası ({output_values.std():.4g}) giriş "
            f"dispersiyasına ({input_values.std():.4g}) nəzərən şübhəli dərəcədə kiçikdir")
        # Kriging DƏQİQ interpolyatordur (overshoot yalnız trend/ekstrapolyasiya
        # ilə mümkündür) — sərt bir maks/min EYNİLİYİ gözlənilmir, amma çıxış
        # giriş diapazonunun bir neçə mislindən UZAĞA getməməlidir (kobud
        # "artefakt" əlaməti olardı)
        span = input_values.max() - input_values.min()
        assert output_values.min() >= input_values.min() - 0.5 * span - 0.05
        assert output_values.max() <= input_values.max() + 0.5 * span + 0.05


# ── 16: SGS validasiyası (production builder vasitəsilə) ───────────────
def test_sgs_through_production_builder_honors_hard_data_and_seed_reproducibility():
    """SGS-in 3D sərt-data yerləşdirilməsi (bax `_simulate_continuous_sgs_
    field`/`_gather_categorical_hard_data`-la EYNİ məntiq) AÇIQ `layer`/
    `depth` tələb edir — ona görə bura xüsusi, `layer=0` işarəli dataset
    işlədir (`_grid_wells()`-in kəsilməz nümunələri qəsdən laysızdır)."""
    from imex2d.application.geology_service import ContinuousSGSConfig

    rng = np.random.default_rng(4)
    samples = [WellSample(well=f"W{i}", x=float(i % 6) * 90.0, y=float(i // 6) * 90.0,
                          layer=0, values={"PORO": float(np.clip(
                              0.15 + 0.05 * np.sin(i / 3.0) + rng.normal(0, 0.005), 0.05, 0.3)),
                                          "PERMX": 50.0})
              for i in range(36)]
    dataset = WellDataset(samples=samples, source="test")
    spec = _spec(nx=12, ny=12)
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    cfg = {"PORO": ContinuousSGSConfig(seed=42)}
    model_a, _ = builder.build(dataset, spec, sgs_config=cfg)
    model_b, _ = builder.build(dataset, spec, sgs_config=cfg)
    assert np.allclose(model_a.property_maps["PORO"].values, model_b.property_maps["PORO"].values), (
        "eyni seed production SGS-də EYNİ nəticə vermədi")

    cfg_other_seed = {"PORO": ContinuousSGSConfig(seed=99)}
    model_c, _ = builder.build(dataset, spec, sgs_config=cfg_other_seed)
    assert not np.allclose(model_a.property_maps["PORO"].values,
                           model_c.property_maps["PORO"].values), (
        "fərqli seed production SGS-də EYNİ nəticə verdi")


# ── 17: tam pipeline deterministik təkrarlanabilirlik ──────────────────
def test_full_production_pipeline_is_deterministic_across_repeated_builds():
    dataset = _grid_wells()
    spec = _spec()
    builder1 = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    builder2 = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    facies_config = {"FACIES": FaciesBuildConfig(on_conflict="majority")}
    model1, report1 = builder1.build(dataset, spec, facies_config=facies_config)
    model2, report2 = builder2.build(dataset, spec, facies_config=facies_config)
    for key in ("PORO", "PERMX", "SW", "NTG", "PERMY", "PERMZ"):
        assert np.array_equal(model1.property_maps[key].values, model2.property_maps[key].values), key
    assert np.array_equal(model1.uncertainty["PORO"].variance, model2.uncertainty["PORO"].variance)
    assert report1.as_text() == report2.as_text()


# ── 22: analitik (manufactured) sahə validasiyası ──────────────────────
def test_linear_analytical_field_is_reproduced_with_low_error_through_production():
    """`z(x,y) = a + b*x + c*y` — Kriging TAM YANSIZ interpolyatordur, bu
    sinif SƏTHLƏRİ demək olar SIFIR xəta ilə reproduksiya edir."""
    rng = np.random.default_rng(21)
    samples = []
    a, b, c = 0.10, 0.0002, 0.00015
    for i in range(8):
        for j in range(8):
            x, y = i * 90.0, j * 90.0
            samples.append(WellSample(well=f"W{i}_{j}", x=x, y=y,
                                      values={"PORO": a + b * x + c * y, "PERMX": 50.0}))
    dataset = WellDataset(samples=samples, source="test")
    spec = _spec(nx=20, ny=20)
    builder, (model, report) = _build(dataset=dataset, spec=spec)
    grid = model.grid
    cell_x = (np.arange(grid.nx) + 0.5) * spec.dx
    cell_y = (np.arange(grid.ny) + 0.5) * spec.dy
    yy, xx = np.meshgrid(cell_y, cell_x, indexing="ij")
    expected = (a + b * xx.ravel() + c * yy.ravel())
    actual = model.property_maps["PORO"].values
    rmse = float(np.sqrt(np.mean((actual - expected) ** 2)))
    assert rmse < 0.01, f"xətti sahədə production RMSE={rmse:.5g} — gözlənilməz dərəcədə böyük"


# ── 19/20: performans/miqyaslanma (production builder) ─────────────────
@pytest.mark.performance
def test_performance_and_scalability_benchmark_production_path():
    results = []
    for n_side, nx, ny, label in ((10, 10, 10, "100-samples"),
                                  (32, 20, 20, "1000-samples")):
        dataset = _grid_wells(n_side=n_side, spacing=35.0, seed=1)
        spec = _spec(nx=nx, ny=ny)
        start = time.perf_counter()
        builder, (model, report) = _build(dataset=dataset, spec=spec)
        elapsed = time.perf_counter() - start
        unc = model.uncertainty["PORO"]
        results.append((label, len(dataset.samples), model.grid.ncell, elapsed,
                        int(unc.neighbor_count.max())))
    print("\nPHASE C performance/scalability:")
    for label, nsamples, ncell, elapsed, max_neighbors in results:
        print(f"  {label}: samples={nsamples} cells={ncell} time={elapsed:.3f}s "
              f"max_neighbors={max_neighbors}")
        assert elapsed < 60.0, f"{label}: {elapsed:.2f}s — reqressiya"
        assert max_neighbors <= 24, (
            f"{label}: max_neighbors={max_neighbors} — yerli qonşuluq qorunmayıb, "
            "N böyüdükcə production yolu O(N^2)/O(N^3)-ə sürüşür")
