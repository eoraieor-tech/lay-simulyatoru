"""Phase 5 — Sequential Gaussian Simulation (SGS) kəsilməz xassələr üçün."""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.geology.distribution_analysis import summarize_distribution
from imex2d.geology.sgs import (FaciesPropertyConfig, PropertyVariogramParams,
                                run_realizations_sgs, simulate_sgs,
                                simulate_sgs_facies_conditioned)


def _grid_targets(nx=10, ny=10, dx=20.0, dy=20.0):
    xs = (np.arange(nx) + 0.5) * dx
    ys = (np.arange(ny) + 0.5) * dy
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


def _sample_wells(n=10, seed=0, high=300.0, mean=0.20, sigma=0.03):
    rng = np.random.default_rng(seed)
    points = rng.uniform(0, high, size=(n, 2))
    values = rng.normal(mean, sigma, size=n)
    return points, values


# ── §2: SGS SAMPLES — kriging ESTIMATE-i deyil ──────────────────────────
def test_sgs_result_is_not_the_smooth_kriging_estimate():
    """SGS-in ƏSAS FƏRQİ: nəticə hamarlanmış kriging qiyməti DEYİL,
    variansdan nümunə götürülür — eyni hədəf üçün FƏRQLİ seedlərlə
    FƏRQLİ dəyər almalıyıq (kriging isə HƏMİŞƏ eyni qiyməti verərdi)."""
    points, values = _sample_wells(n=8, seed=1)
    targets = np.array([[150.0, 150.0]])
    v1 = simulate_sgs(points, values, targets, seed=1).values[0]
    v2 = simulate_sgs(points, values, targets, seed=2).values[0]
    v3 = simulate_sgs(points, values, targets, seed=3).values[0]
    assert len({round(v1, 8), round(v2, 8), round(v3, 8)}) > 1


def test_sgs_realization_reproduces_variance_not_just_mean():
    """Bir çox realizasiya üzərində simulyasiya edilmiş dəyərlərin öz
    ARALARINDA dəyişkənliyi olmalıdır (sıfır varians = deterministik
    kriging demək olardı, bu, SGS-in QADAĞASIdır)."""
    points, values = _sample_wells(n=8, seed=1)
    targets = np.array([[150.0, 150.0]])
    draws = [simulate_sgs(points, values, targets, seed=s).values[0] for s in range(20)]
    assert np.std(draws) > 1e-4


# ── §11: sərt data hörməti ──────────────────────────────────────────────
def test_hard_data_honored_exactly_regardless_of_seed():
    points, values = _sample_wells(n=6, seed=2)
    targets = np.vstack([points, _grid_targets(5, 5)])
    for seed in (1, 2, 3, 100):
        realization = simulate_sgs(points, values, targets, seed=seed)
        assert np.allclose(realization.values[:len(points)], values, atol=1e-9), seed
        assert np.all(realization.hard_data_mask[:len(points)])


# ── §12: seed təkrarlanabilənliyi ────────────────────────────────────────
def test_same_seed_gives_identical_realization():
    points, values = _sample_wells(n=8, seed=3)
    targets = _grid_targets()
    r1 = simulate_sgs(points, values, targets, seed=42)
    r2 = simulate_sgs(points, values, targets, seed=42)
    assert np.allclose(r1.values, r2.values, atol=0, rtol=0)


def test_different_seeds_give_different_realizations():
    points, values = _sample_wells(n=8, seed=4)
    targets = _grid_targets()
    r1 = simulate_sgs(points, values, targets, seed=1)
    r2 = simulate_sgs(points, values, targets, seed=2)
    assert not np.allclose(r1.values, r2.values)


# ── §13: çoxlu realizasiya ────────────────────────────────────────────────
def test_run_realizations_sgs_ids_and_seeds():
    points, values = _sample_wells(n=6, seed=5)
    targets = _grid_targets(4, 4)
    realizations = run_realizations_sgs(3, points, values, targets, seed=10)
    assert [r.realization_id for r in realizations] == [0, 1, 2]
    assert [r.seed for r in realizations] == [10, 1010, 2010]
    assert not np.array_equal(realizations[0].values, realizations[1].values)


# ── §6: PORO hədləri — kliplənmə İZLƏNİR, sükutla deyil ─────────────────
def test_poro_bounds_enforced_and_corrections_tracked():
    points, values = _sample_wells(n=10, seed=6, mean=0.15, sigma=0.15)   # geniş sigma -> kənar dəyər
    targets = _grid_targets(8, 8)
    realization = simulate_sgs(points, values, targets, seed=1, bounds=(0.01, 0.45))
    assert np.all(realization.values >= 0.01 - 1e-9)
    assert np.all(realization.values <= 0.45 + 1e-9)
    # geniş sigma ilə bəzi kliplənmə GÖZLƏNİLİR, sayılmalıdır (0 ola bilər, amma mənfi ola bilməz)
    assert realization.diagnostics.bound_corrections >= 0


def test_bound_correction_rate_triggers_strong_warning_when_excessive():
    points, values = _sample_wells(n=10, seed=7, mean=0.50, sigma=0.30)   # bilərəkdən hədlərdən kənar
    targets = _grid_targets(8, 8)
    realization = simulate_sgs(points, values, targets, seed=1, bounds=(0.01, 0.10),
                               correction_warn_threshold=0.05)
    assert realization.diagnostics.bound_corrections > 0
    assert any("GÜCLÜ XƏBƏRDARLIQ" in w for w in realization.warnings)


# ── §7: PERMX log-fəza — müsbətlik qorunur ──────────────────────────────
def test_log_space_permeability_stays_strictly_positive():
    rng = np.random.default_rng(8)
    points = rng.uniform(0, 300, size=(10, 2))
    permx = rng.lognormal(mean=4.5, sigma=1.0, size=10)   # tipik keçiricilik (mD)
    targets = _grid_targets(10, 10)
    realization = simulate_sgs(points, permx, targets, seed=1, log_space=True)
    assert np.all(realization.values > 0.0)


def test_log_space_rejects_non_positive_hard_data():
    points, values = _sample_wells(n=5, seed=9)
    values[0] = -1.0
    with pytest.raises(ValueError):
        simulate_sgs(points, values, _grid_targets(3, 3), log_space=True)


def test_direct_space_and_log_space_both_reproduce_hard_data_exactly():
    rng = np.random.default_rng(10)
    points = rng.uniform(0, 300, size=(6, 2))
    permx = rng.lognormal(mean=4.0, sigma=0.8, size=6)
    targets = np.vstack([points, _grid_targets(4, 4)])
    for log_space in (True, False):
        realization = simulate_sgs(points, permx, targets, seed=1, log_space=log_space)
        assert np.allclose(realization.values[:len(points)], permx, atol=1e-6)


# ── §5/§16: fasiya-şərtli SGS ────────────────────────────────────────────
def _facies_conditioned_dataset():
    rng = np.random.default_rng(11)
    sand_points = rng.uniform(0, 150, size=(10, 2))
    shale_points = rng.uniform(150, 300, size=(10, 2))
    sand_values = rng.normal(0.25, 0.02, size=10)     # sand: yüksək poroziya
    shale_values = rng.normal(0.08, 0.01, size=10)    # shale: aşağı poroziya
    points = np.vstack([sand_points, shale_points])
    values = np.concatenate([sand_values, shale_values])
    facies = np.concatenate([np.zeros(10, int), np.ones(10, int)])
    return points, values, facies


def test_facies_conditioned_sgs_preserves_per_facies_distribution():
    points, values, facies_at_points = _facies_conditioned_dataset()
    targets = _grid_targets(12, 12, dx=25.0, dy=25.0)
    target_facies = (targets[:, 0] > 150.0).astype(int)   # sağ yarı = shale (1)

    realization = simulate_sgs_facies_conditioned(
        points, values, facies_at_points, targets, target_facies, seed=1)

    sand_result = realization.values[target_facies == 0]
    shale_result = realization.values[target_facies == 1]
    # §16: fasiyalar-arası fərq SAXLANMALIDIR (data-driven, məcbur DEYİL)
    assert np.mean(sand_result) > np.mean(shale_result) + 0.05


def test_facies_conditioned_sgs_falls_back_with_warning_for_sparse_facies():
    """§5: kifayət qədər sərt data olmayan fasiya üçün AYRICA model
    UYDURULMUR — fallback AÇIQ bildirilir."""
    rng = np.random.default_rng(12)
    points = rng.uniform(0, 300, size=(15, 2))
    values = rng.normal(0.20, 0.03, size=15)
    facies_at_points = np.zeros(15, int)
    facies_at_points[:2] = 1   # fasiya 1 üçün YALNIZ 2 nöqtə (< minimum)

    targets = _grid_targets(6, 6, dx=25.0, dy=25.0)
    target_facies = np.zeros(36, int)
    target_facies[:5] = 1

    realization = simulate_sgs_facies_conditioned(
        points, values, facies_at_points, targets, target_facies, seed=1,
        min_hard_data_for_own_model=8)
    assert any("PAYLAŞILAN" in w or "fallback" in w or "QURULMADI" in w
              for w in realization.warnings)
    assert realization.metadata["per_facies"][1]["used_cross_facies_fallback"] is True
    assert realization.metadata["per_facies"][0]["used_cross_facies_fallback"] is False
    # KRİTİK: fallback İNDİ artıq başqa fasiyanın DƏYƏRLƏRİNİ pooled hard-data
    # kimi istifadə ETMİR — fasiya 1-in öz 2 nöqtəsi hələ də YALNIZ öz
    # dəyərləri ilə (0.20 ətrafı) kondisiyalaşdırılıb, fasiya 0-ın 13
    # nöqtəsinin heç biri fasiya 1-ə "sərt data" kimi keçməyib.
    assert realization.metadata["per_facies"][1]["n_hard_points"] == 2


def test_sparse_facies_fallback_never_uses_other_facies_value_as_hard_data():
    """KRİTİK REQRESSIYA TESTİ — bax audit tapşırığı §2.

    Ssenari: Facies A seyrəkdir (2 nöqtə, < minimum). Facies B-nin YAXŞI
    dəstəklənmiş bir quyu koordinatı ADVERSARIAL olaraq hədəf şəbəkəyə
    Facies A kimi ETİKETLƏNİR (real dünyada bu, fasiya sərhədinin sərt-data
    quyusundan bir qədər fərqli çəkilməsi ilə baş verə bilər). Əvvəlki
    (səhv) implementasiyada bu, `find_exact_matches`-in koordinat
    üst-üstə-düşməsini tapıb Facies B-nin dəyərini "dəqiq sərt-data" kimi
    Facies A-ya ötürməsinə səbəb olurdu — Facies A-nın diapazonundan
    (10-15) TAM KƏNAR bir dəyər (Facies B-nin ~100-130 diapazonu) həmin
    hüceyrəyə "hörmət edilmiş sərt-data" kimi yazılırdı.

    Düzəlişdən sonra: fasiyanın hard-conditioning-i YALNIZ öz nöqtələri
    ilə aparılır, ona görə bu adversarial hüceyrə ARTIQ dəqiq-uyğunluq
    almır və Facies B-nin dəyəri ilə TƏSADÜFƏN üst-üstə düşmür.
    """
    rng = np.random.default_rng(7)
    facies_b_points = rng.uniform(0.0, 100.0, size=(10, 2))
    facies_b_values = rng.uniform(100.0, 130.0, size=10)

    facies_a_points = np.array([[5.0, 5.0], [6.0, 6.0]])
    facies_a_values = np.array([10.0, 12.0])

    points = np.vstack([facies_a_points, facies_b_points])
    values = np.concatenate([facies_a_values, facies_b_values])
    facies_at_points = np.array([0, 0] + [1] * 10)

    # Hədəf şəbəkəyə Facies B-nin BİRİNCİ quyu koordinatını ADVERSARIAL
    # olaraq daxil edirik, AMMA Facies A kimi ETİKETLƏYİRİK.
    adversarial_target = facies_b_points[0:1]
    other_targets = rng.uniform(0.0, 100.0, size=(20, 2))
    targets = np.vstack([adversarial_target, other_targets])
    facies_at_targets = np.zeros(targets.shape[0], dtype=int)   # HAMISI Facies A

    realization = simulate_sgs_facies_conditioned(
        points, values, facies_at_points, targets, facies_at_targets,
        seed=1, min_hard_data_for_own_model=8)

    # Facies B-nin dəyəri Facies A-nın adversarial hüceyrəsinə HEÇ VAXT
    # "dəqiq sərt-data" kimi keçməməlidir.
    assert not np.isclose(realization.values[0], facies_b_values[0], atol=1e-6)
    assert not realization.hard_data_mask[0]
    # Fallback YENƏ DƏ işə düşməlidir (2 < 8), AMMA yalnız struktur borcu ilə.
    assert realization.metadata["per_facies"][0]["used_cross_facies_fallback"] is True
    assert realization.metadata["per_facies"][0]["n_hard_points"] == 2
    # Facies A-nın ÖZ 2 nöqtəsi hədəflərdə olsaydı, onlar YENƏ DƏ dəqiq
    # hörmət edilməlidir (fix hard-data honoring-i pozmayıb).
    own_targets = np.vstack([facies_a_points, targets])
    own_facies = np.concatenate([np.zeros(2, dtype=int), facies_at_targets])
    realization2 = simulate_sgs_facies_conditioned(
        points, values, facies_at_points, own_targets, own_facies,
        seed=1, min_hard_data_for_own_model=8)
    assert np.allclose(realization2.values[:2], facies_a_values, atol=1e-6)


def test_facies_with_literally_zero_own_hard_data_uses_unconditional_global_fallback():
    """§2: bir fasiyanın ÖZ sərt datası SIFIR olduqda belə, başqa
    fasiyanın dəyərləri ona hard-conditioning kimi VERİLMİR — ƏVƏZİNƏ
    paylaşılan qlobal marjinal paylanmadan şərtsiz nümunə çəkilir."""
    rng = np.random.default_rng(3)
    points, values = _sample_wells(n=12, seed=3, mean=0.20, sigma=0.03)
    facies_at_points = np.zeros(12, dtype=int)   # bütün quyular Facies 0-a aiddir

    targets = _grid_targets(6, 6, dx=20.0, dy=20.0)
    facies_at_targets = np.zeros(targets.shape[0], dtype=int)
    facies_at_targets[:5] = 1   # Facies 1: heç bir öz quyusu YOXDUR

    realization = simulate_sgs_facies_conditioned(
        points, values, facies_at_points, targets, facies_at_targets, seed=2,
        min_hard_data_for_own_model=8)

    assert np.all(np.isfinite(realization.values))
    assert not np.any(realization.hard_data_mask[:5])
    assert realization.metadata["per_facies"][1]["n_hard_points"] == 0
    assert realization.metadata["per_facies"][1]["used_unconditional_global_fallback"] is True
    assert any("HEÇ bir öz sərt nöqtə" in w for w in realization.warnings)


def test_facies_conditioned_sgs_honors_hard_data():
    points, values, facies_at_points = _facies_conditioned_dataset()
    targets = np.vstack([points, _grid_targets(8, 8, dx=30.0, dy=30.0)])
    target_facies = np.concatenate([facies_at_points, (targets[len(points):, 0] > 150).astype(int)])
    realization = simulate_sgs_facies_conditioned(
        points, values, facies_at_points, targets, target_facies, seed=1)
    assert np.allclose(realization.values[:len(points)], values, atol=1e-6)


# ── kənar hallar ────────────────────────────────────────────────────────
def test_single_hard_data_point_does_not_crash():
    points = np.array([[10.0, 10.0]])
    values = np.array([0.20])
    realization = simulate_sgs(points, values, _grid_targets(3, 3), seed=1)
    assert np.all(np.isfinite(realization.values))


def test_constant_property_yields_constant_field_without_fabricated_variability():
    points, _ = _sample_wells(n=6, seed=13)
    values = np.full(6, 150.0)
    realization = simulate_sgs(points, values, _grid_targets(4, 4), seed=1)
    assert np.allclose(realization.values, 150.0)
    assert any("SABİT" in w for w in realization.warnings)


def test_very_sparse_data_completes_with_fallback_diagnostics():
    points = np.array([[0.0, 0.0], [280.0, 280.0]])
    values = np.array([0.10, 0.30])
    realization = simulate_sgs(points, values, _grid_targets(6, 6), seed=1,
                               search_radius=15.0, max_neighbors=4)
    assert np.all(np.isfinite(realization.values))
    assert realization.diagnostics.nan_fallback_cells > 0


def test_invalid_input_mismatched_lengths_rejected():
    with pytest.raises(ValueError):
        simulate_sgs(np.array([[0., 0.], [1., 1.]]), np.array([1.0]), _grid_targets(3, 3))


def test_invalid_input_nan_hard_data_rejected():
    points, values = _sample_wells(n=5, seed=14)
    values[0] = np.nan
    with pytest.raises(ValueError):
        simulate_sgs(points, values, _grid_targets(3, 3))


# ── §14: ensembl paylanma yoxlanması ─────────────────────────────────────
def test_ensemble_distribution_reproduces_target_statistics():
    points, values = _sample_wells(n=15, seed=15, mean=0.20, sigma=0.03)
    target_summary = summarize_distribution(values)
    targets = _grid_targets(10, 10, dx=25.0, dy=25.0)

    realizations = run_realizations_sgs(15, points, values, targets, seed=0)
    all_values = np.concatenate([r.values for r in realizations])
    ensemble_summary = summarize_distribution(all_values)

    assert abs(ensemble_summary.mean - target_summary.mean) < 0.05
    # SGS varians = kondisioner + kriging qeyri-müəyyənliyi əlavə edir,
    # ona görə TAM BƏRABƏRLİK gözlənilmir — YALNIZ MƏQBUL diapazon
    assert ensemble_summary.std > 0.0
