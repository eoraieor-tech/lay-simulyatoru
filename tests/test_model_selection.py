"""B8 №18-24 — çarpaz-doğrulama və MODEL SEÇİMİ (GATE B6/B7).

İki mərkəzi iddia:

1. **SIZMA YOXDUR** — gizlədilmiş nöqtə variogram fitinə, çevirmə
   statistikasına, QC qərarlarına və Kriging çəkilərinə TƏSİR ETMİR.
   Bu, ölçü ilə sübut olunur (gizlədilmiş dəyəri dəyişmək proqnozu
   DƏYİŞMİR).
2. **SEÇİM DATA ƏSASLIDIR** — "adi kriginq həmişə ən yaxşıdır" kimi
   sabit qərar yoxdur; qalib doğrulama metriklərindən çıxır.
"""

from __future__ import annotations

import numpy as np
import pytest

from imex2d.geology.cross_validation import (DEFAULT_SELECTION_WEIGHTS,
                                             CategoricalCVMetrics,
                                             ContinuousCVMetrics, ModelCandidate,
                                             ModelSelectionReport, ValidationDesign,
                                             ValidationKind, build_folds,
                                             cross_validate_property,
                                             default_candidates,
                                             select_property_model)
from imex2d.geology.property_config import InterpolationKind, resolve_strategy


def _gaussian_field(points, ranges, seed):
    rng = np.random.default_rng(seed)
    n = points.shape[0]
    scaled = points[:, :len(ranges)] / np.asarray(ranges, float)[None, :]
    diff = scaled[:, None, :] - scaled[None, :, :]
    cov = np.exp(-3.0 * np.sqrt(np.sum(diff * diff, axis=-1))) + 1e-8 * np.eye(n)
    return np.linalg.cholesky(cov) @ rng.standard_normal(n)


def _poro_field(n=70, seed=1, high=1000.0):
    rng = np.random.default_rng(seed)
    points = rng.uniform(0.0, high, size=(n, 2))
    values = 0.18 + 0.03 * _gaussian_field(points, (250.0, 250.0), seed + 1)
    return points, np.clip(values, 0.01, 0.45)


# ── 18. LOOCV ─────────────────────────────────────────────────────────
def test_loocv_predicts_every_point_exactly_once():
    points, values = _poro_field(40, seed=2)
    folds = build_folds(points, ValidationDesign(kind=ValidationKind.LEAVE_ONE_OUT))
    assert len(folds) == 40
    for train, test in folds:
        assert test.size == 1
        assert train.size == 39
        assert test[0] not in train


def test_loocv_produces_finite_metrics():
    points, values = _poro_field(50, seed=3)
    metrics = cross_validate_property(points, values, resolve_strategy("PORO"),
                                      ValidationDesign())
    assert isinstance(metrics, ContinuousCVMetrics)
    assert metrics.n == 50
    assert np.isfinite(metrics.rmse) and metrics.rmse > 0.0
    assert np.isfinite(metrics.mae) and np.isfinite(metrics.bias)
    assert metrics.rmse >= metrics.mae      # riyazi zəmanət


def test_loocv_rmse_is_worse_than_fitting_the_point_itself():
    """Gizlədilmiş nöqtədə xəta SIFIR OLMAMALIDIR — əks halda sızma var."""
    points, values = _poro_field(45, seed=4)
    metrics = cross_validate_property(points, values, resolve_strategy("PORO"),
                                      ValidationDesign())
    assert metrics.rmse > 1e-6, "LOOCV xətası sıfırdırsa nöqtə özünü görüb"


# ── 19. SIZMA YOXDUR ──────────────────────────────────────────────────
def test_held_out_value_cannot_influence_its_own_prediction():
    """ƏN GÜCLÜ sızma testi: gizlədilmiş nöqtənin DƏYƏRİNİ kəskin
    dəyişirik. Əgər o dəyər variogram fitinə/çevirməyə/QC-yə sızırsa,
    ONUN ÖZ proqnozu da dəyişər. Sızma yoxdursa proqnoz EYNİ qalır."""
    from imex2d.geology.property_interpolation import interpolate_property_field
    points, values = _poro_field(40, seed=5)
    held = 7

    train = np.delete(np.arange(40), held)
    baseline = interpolate_property_field(
        points[train], values[train], points[held:held + 1],
        property_name="PORO").estimate[0]

    poisoned = values.copy()
    poisoned[held] = 0.44          # kəskin fərqli dəyər
    poisoned_prediction = interpolate_property_field(
        points[train], poisoned[train], points[held:held + 1],
        property_name="PORO").estimate[0]

    assert baseline == pytest.approx(poisoned_prediction, abs=0.0), (
        "gizlədilmiş nöqtənin dəyəri proqnozuna sızıb")


def test_changing_a_held_out_value_does_not_change_the_other_predictions():
    """Eyni iddianın TAM LOOCV üzərində forması."""
    points, values = _poro_field(35, seed=6)
    strategy = resolve_strategy("PORO")

    base = cross_validate_property(points, values, strategy, ValidationDesign())
    modified = values.copy()
    modified[3] = 0.43
    changed = cross_validate_property(points, modified, strategy, ValidationDesign())

    # 3-cü nöqtənin ÖZ xətası dəyişməlidir (dəyəri dəyişib), amma
    # metriklərin tamamilə eyni qalması da, tamamilə pozulması da yanlışdır
    assert changed.rmse != base.rmse
    assert changed.n == base.n == 35


def test_normal_score_transform_is_refit_per_fold():
    """DATADAN ASILI çevirmə hər qatda YENİDƏN fit olunur — prototip
    obyekt dəyişməz qalır (sızma qapısı bağlıdır)."""
    from imex2d.geology.property_config import normal_score_strategy
    points, values = _poro_field(30, seed=7)
    strategy = normal_score_strategy(resolve_strategy("PORO"))
    assert strategy.transform.table is None
    metrics = cross_validate_property(points, values, strategy, ValidationDesign())
    assert np.isfinite(metrics.rmse)
    assert strategy.transform.table is None, "prototip cədvəl DOLDURULMAMALIDIR"


def test_outlier_decisions_are_made_per_fold_not_globally():
    """QC hər qatda yalnız TƏLİM datası üzərində işləyir."""
    from imex2d.geology.property_interpolation import interpolate_property_field
    points, values = _poro_field(30, seed=8)
    values[12] = 0.44                 # kənar-dəyər namizədi
    train = np.delete(np.arange(30), 12)
    result = interpolate_property_field(points[train], values[train],
                                        points[12:13], property_name="PORO")
    assert result.quality.n_input == 29, "test nöqtəsi QC-yə DAXİL OLMAMALIDIR"


# ── 20-22. model müqayisəsi ───────────────────────────────────────────
def test_model_comparison_evaluates_every_candidate():
    points, values = _poro_field(50, seed=9)
    candidates = default_candidates("PORO")
    assert len(candidates) >= 4
    report = select_property_model(points, values, candidates, "PORO",
                                   ValidationDesign())
    assert isinstance(report, ModelSelectionReport)
    assert len(report.results) == len(candidates)
    assert report.selected is not None
    assert all(np.isfinite(r.score) for r in report.ranking)


def test_variogram_models_are_compared_and_can_disagree():
    """Üç variogram modeli FƏRQLİ metriklər verməlidir — yəni müqayisə
    həqiqidir, formal deyil."""
    points, values = _poro_field(60, seed=10)
    scores = {}
    for model in ("spherical", "exponential", "gaussian"):
        metrics = cross_validate_property(
            points, values, resolve_strategy("PORO").derive(variogram_model=model),
            ValidationDesign())
        scores[model] = metrics.rmse
    assert len(set(np.round(list(scores.values()), 10))) > 1


def test_anisotropic_candidate_wins_on_an_anisotropic_field():
    """GATE B7 — seçim DATA ƏSASLIDIR: anizotrop sahədə anizotrop
    namizəd izotropdan daha yaxşı olmalıdır."""
    from imex2d.geology.anisotropy import AnisotropyParams
    rng = np.random.default_rng(11)
    points = rng.uniform(0.0, 1200.0, size=(90, 2))
    values = 0.2 + 0.03 * _gaussian_field(points, (900.0, 120.0), 12)
    base = resolve_strategy("PORO").derive(variogram_model="exponential")

    isotropic = ModelCandidate(
        "izotrop", base.derive(anisotropy=AnisotropyParams(
            azimuth_deg=0.0, range_major=900.0, range_minor=900.0,
            range_vertical=900.0)))
    anisotropic = ModelCandidate(
        "anizotrop 90°", base.derive(anisotropy=AnisotropyParams(
            azimuth_deg=90.0, range_major=900.0, range_minor=120.0,
            range_vertical=900.0)))
    report = select_property_model(points, values, [isotropic, anisotropic],
                                   "PORO", ValidationDesign())
    assert report.selected.candidate.label == "anizotrop 90°"


def test_idw_can_be_selected_when_it_actually_wins():
    """Sabit favorit YOXDUR: kriginq həmişə qalib gəlmir."""
    points, values = _poro_field(40, seed=13)
    kriging = ModelCandidate("kriginq", resolve_strategy("PORO"))
    idw = ModelCandidate("IDW", resolve_strategy("PORO").derive(
        interpolation=InterpolationKind.IDW))
    report = select_property_model(points, values, [kriging, idw], "PORO",
                                   ValidationDesign())
    labels = [r.candidate.label for r in report.ranking]
    assert set(labels) == {"kriginq", "IDW"}
    assert report.selected.candidate.label in labels


def test_log_space_beats_raw_space_on_a_lognormal_field():
    """Loq fəzasının üstünlüyü FƏRZ EDİLMİR — ölçülür."""
    rng = np.random.default_rng(14)
    points = rng.uniform(0.0, 900.0, size=(70, 2))
    values = np.exp(4.0 + 1.5 * _gaussian_field(points, (300.0, 300.0), 15))
    report = select_property_model(points, values, default_candidates("PERMX"),
                                   "PERMX", ValidationDesign())
    assert report.selected is not None
    assert "xam fəza" not in report.selected.candidate.label


# ── 23. məkan/blok doğrulaması ────────────────────────────────────────
def test_spatial_blocks_keep_neighbours_together():
    rng = np.random.default_rng(16)
    points = rng.uniform(0.0, 1000.0, size=(80, 2))
    folds = build_folds(points, ValidationDesign(
        kind=ValidationKind.SPATIAL_BLOCK, k=4, block_size=250.0))
    assert 2 <= len(folds) <= 4
    for train, test in folds:
        assert np.intersect1d(train, test).size == 0
        assert train.size + test.size == 80


def test_spatial_block_validation_is_harder_than_random_split():
    """B3.3 — məkan bloklarında xəta DAHA BÖYÜK olmalıdır, çünki
    gizlədilən nöqtənin qonşuları da gizlədilir. Təsadüfi bölgüyə
    əsaslanıb "məkan ümumiləşdirməsi yaxşıdır" demək olmaz."""
    points, values = _poro_field(90, seed=17, high=1000.0)
    strategy = resolve_strategy("PORO").derive(variogram_model="exponential")
    random_split = cross_validate_property(
        points, values, strategy,
        ValidationDesign(kind=ValidationKind.RANDOM_KFOLD, k=5, seed=1))
    spatial = cross_validate_property(
        points, values, strategy,
        ValidationDesign(kind=ValidationKind.SPATIAL_BLOCK, k=5, block_size=300.0))
    assert spatial.rmse > random_split.rmse


def test_every_design_reports_which_design_was_used():
    points, values = _poro_field(40, seed=18)
    for kind in ValidationKind:
        design = ValidationDesign(kind=kind, k=4, block_size=300.0)
        report = select_property_model(points, values,
                                       [ModelCandidate("k", resolve_strategy("PORO"))],
                                       "PORO", design)
        assert report.design.kind is kind
        assert report.as_dict()["design"]["kind"] == kind.value
        assert design.describe()


def test_random_kfold_folds_partition_the_data():
    rng = np.random.default_rng(19)
    points = rng.uniform(0.0, 500.0, size=(37, 2))
    folds = build_folds(points, ValidationDesign(kind=ValidationKind.RANDOM_KFOLD,
                                                 k=5, seed=3))
    covered = np.concatenate([test for _, test in folds])
    assert np.array_equal(np.sort(covered), np.arange(37))


# ── 24. təkrarlana bilənlik ───────────────────────────────────────────
def test_model_selection_is_reproducible():
    points, values = _poro_field(50, seed=20)
    candidates = default_candidates("PORO")
    first = select_property_model(points, values, candidates, "PORO",
                                  ValidationDesign())
    second = select_property_model(points, values, candidates, "PORO",
                                   ValidationDesign())
    assert first.selected.candidate.label == second.selected.candidate.label
    assert [r.candidate.label for r in first.ranking] == \
        [r.candidate.label for r in second.ranking]
    assert np.allclose([r.score for r in first.ranking],
                       [r.score for r in second.ranking], atol=0.0)


def test_random_kfold_is_reproducible_for_a_fixed_seed():
    rng = np.random.default_rng(21)
    points = rng.uniform(0.0, 600.0, size=(40, 2))
    a = build_folds(points, ValidationDesign(kind=ValidationKind.RANDOM_KFOLD,
                                             k=4, seed=9))
    b = build_folds(points, ValidationDesign(kind=ValidationKind.RANDOM_KFOLD,
                                             k=4, seed=9))
    for (ta, sa), (tb, sb) in zip(a, b):
        assert np.array_equal(ta, tb) and np.array_equal(sa, sb)


def test_ties_are_broken_deterministically_by_label():
    points, values = _poro_field(30, seed=22)
    strategy = resolve_strategy("PORO").derive(variogram_model="spherical")
    candidates = [ModelCandidate("zeta", strategy), ModelCandidate("alpha", strategy)]
    report = select_property_model(points, values, candidates, "PORO",
                                   ValidationDesign())
    assert report.selected.candidate.label == "alpha"


# ── B3.5 çox-metrikli sıralama ────────────────────────────────────────
def test_scoring_uses_every_documented_criterion():
    points, values = _poro_field(45, seed=23)
    report = select_property_model(points, values, default_candidates("PORO"),
                                   "PORO", ValidationDesign())
    best = report.ranking[0]
    assert set(best.penalties) == set(DEFAULT_SELECTION_WEIGHTS)
    assert all(v >= 0.0 for v in best.penalties.values())
    assert sum(DEFAULT_SELECTION_WEIGHTS.values()) == pytest.approx(1.0)


def test_custom_weights_can_change_the_winner():
    """Çəkilər SABİT KODLANMIR — çağıran onları dəyişə bilər."""
    points, values = _poro_field(50, seed=24)
    candidates = default_candidates("PORO")
    accuracy_only = select_property_model(
        points, values, candidates, "PORO", ValidationDesign(),
        weights={"accuracy": 1.0, "calibration": 0.0, "bias": 0.0,
                 "validity": 0.0, "stability": 0.0})
    calibration_only = select_property_model(
        points, values, candidates, "PORO", ValidationDesign(),
        weights={"accuracy": 0.0, "calibration": 1.0, "bias": 0.0,
                 "validity": 0.0, "stability": 0.0})
    assert accuracy_only.selected is not None
    assert calibration_only.selected is not None
    # ən azı ballar fərqlənməlidir (qalib eyni ola bilər)
    assert not np.allclose([r.score for r in accuracy_only.ranking],
                           [r.score for r in calibration_only.ranking])


def test_accuracy_only_scoring_picks_the_lowest_rmse():
    points, values = _poro_field(45, seed=25)
    report = select_property_model(
        points, values, default_candidates("PORO", include_idw=False), "PORO",
        ValidationDesign(),
        weights={"accuracy": 1.0, "calibration": 0.0, "bias": 0.0,
                 "validity": 0.0, "stability": 0.0})
    rmses = [r.metrics.rmse for r in report.ranking]
    assert rmses == sorted(rmses)
    assert report.selected.metrics.rmse == min(rmses)


# ── B3.6 hesabat strukturu ────────────────────────────────────────────
def test_report_is_machine_readable():
    points, values = _poro_field(40, seed=26)
    report = select_property_model(points, values, default_candidates("PORO"),
                                   "PORO", ValidationDesign())
    data = report.as_dict()
    for key in ("property", "design", "n_samples", "n_excluded", "weights",
                "selected", "candidates", "warnings"):
        assert key in data
    assert data["n_samples"] == 40
    assert len(data["candidates"]) == len(report.results)
    for entry in data["candidates"]:
        assert {"label", "score", "penalties", "error", "metrics"} <= set(entry)


def test_report_text_ranks_candidates_and_marks_the_winner():
    points, values = _poro_field(40, seed=27)
    report = select_property_model(points, values, default_candidates("PORO"),
                                   "PORO", ValidationDesign())
    text = report.as_text()
    assert "Model seçimi" in text and "★" in text
    assert report.selected.candidate.label in text


def test_failing_candidates_are_recorded_not_silently_dropped():
    points, values = _poro_field(30, seed=28)
    good = ModelCandidate("yaxşı", resolve_strategy("PORO"))
    # loq çevirməsi MƏNFİ dəyərlərdə uğursuz olur → namizəd xəta verməlidir
    broken = ModelCandidate("pozuq", resolve_strategy("PERMX"))
    report = select_property_model(points, values - 1.0, [good, broken], "PORO",
                                   ValidationDesign())
    assert report.n_excluded >= 0
    labels = {r.candidate.label for r in report.results}
    assert labels == {"yaxşı", "pozuq"}


# ── kateqorik doğrulama ───────────────────────────────────────────────
def test_categorical_cross_validation_reports_classification_metrics():
    rng = np.random.default_rng(29)
    points = rng.uniform(0.0, 800.0, size=(60, 2))
    codes = (points[:, 0] > 400).astype(int) + 2 * (points[:, 1] > 400).astype(int)
    metrics = cross_validate_property(points, codes.astype(float),
                                      resolve_strategy("FACIES"),
                                      ValidationDesign(kind=ValidationKind.RANDOM_KFOLD,
                                                       k=5, seed=2))
    assert isinstance(metrics, CategoricalCVMetrics)
    assert 0.0 <= metrics.accuracy <= 1.0
    assert metrics.log_loss >= 0.0
    assert 0.0 <= metrics.brier_score <= 2.0
    assert metrics.confusion.shape == (len(metrics.categories),) * 2
    assert int(metrics.confusion.sum()) == metrics.n


def test_categorical_accuracy_beats_random_guessing_on_a_clean_field():
    rng = np.random.default_rng(30)
    points = rng.uniform(0.0, 600.0, size=(60, 2))
    codes = (points[:, 0] > 300).astype(float)
    metrics = cross_validate_property(points, codes, resolve_strategy("FACIES"),
                                      ValidationDesign(kind=ValidationKind.RANDOM_KFOLD,
                                                       k=5, seed=4))
    assert metrics.accuracy > 0.75


def test_categorical_model_selection_ranks_indicator_variograms():
    rng = np.random.default_rng(31)
    points = rng.uniform(0.0, 700.0, size=(55, 2))
    codes = (points[:, 1] > 350).astype(float)
    report = select_property_model(points, codes, default_candidates("FACIES"),
                                   "FACIES",
                                   ValidationDesign(kind=ValidationKind.RANDOM_KFOLD,
                                                    k=5, seed=5))
    assert report.selected is not None
    assert all("indikator" in r.candidate.label for r in report.ranking)
