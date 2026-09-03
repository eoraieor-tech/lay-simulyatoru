"""B9 — SİNTETİK REZERVUAR DOĞRULAMASI (7 hal).

Hər halda HƏQİQƏT MƏLUMDUR (sahə nəzarət altında qurulur), ona görə
"işlədi" iddiası ölçülə bilir:

    Hal 1  məsaməlik         məlum kəsilməz sahə
    Hal 2  loq-keçiricilik   məlum loq-normal sahə (median/orta fərqi)
    Hal 3  doyma             [0,1] hədli sahə
    Hal 4  NTG               [0,1] hədli sahə + CV
    Hal 5  fasiya            kateqorik sahə (ehtimallar)
    Hal 6  anizotrop rezervuar (güclü üfüqi, zəif şaquli davamlılıq)
    Hal 7  seyrək quyular    qeyri-müəyyənlik məsafə ilə ARTIR

Sahələr kovariasiya matrisinin Xolesskiy parçalanması ilə qurulur —
"korrelyasiyalı görünən" səs-küy DEYİL.
"""

from __future__ import annotations

import numpy as np

from imex2d.geology.cross_validation import (ValidationDesign, ValidationKind,
                                             cross_validate_property,
                                             default_candidates,
                                             select_property_model)
from imex2d.geology.property_config import (BackTransform, VariableType,
                                            resolve_strategy)
from imex2d.geology.property_interpolation import (Confidence,
                                                   interpolate_categorical_field,
                                                   interpolate_property_field)
from imex2d.geology.sgs import PropertyVariogramParams
from imex2d.geology.sgs_ensemble import simulate_sgs_ensemble, validate_realization
from imex2d.geology.variogram import VariogramParameters, experimental_variogram


def _field(points, ranges, seed):
    """`exp(−3·d_ani)` kovariasiyalı standart qauss sahəsi."""
    rng = np.random.default_rng(seed)
    n = points.shape[0]
    scaled = points[:, :len(ranges)] / np.asarray(ranges, float)[None, :]
    diff = scaled[:, None, :] - scaled[None, :, :]
    cov = np.exp(-3.0 * np.sqrt(np.sum(diff * diff, axis=-1))) + 1e-8 * np.eye(n)
    return np.linalg.cholesky(cov) @ rng.standard_normal(n)


def _sample_and_grid(n_wells, seed, high=1200.0, n_grid=22):
    rng = np.random.default_rng(seed)
    wells = rng.uniform(0.0, high, size=(n_wells, 2))
    axis = np.linspace(0.0, high, n_grid)
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    grid = np.column_stack([xx.ravel(), yy.ravel()])
    return wells, grid


# ══ Hal 1 — MƏSAMƏLİK ═════════════════════════════════════════════════
def test_case1_porosity_is_interpolated_and_validated():
    wells, grid = _sample_and_grid(70, seed=101)
    truth = 0.18 + 0.035 * _field(wells, (300.0, 300.0), 102)
    values = np.clip(truth, 0.02, 0.40)

    result = interpolate_property_field(wells, values, grid, property_name="PORO")
    assert np.all(np.isfinite(result.estimate))
    assert np.all(result.estimate >= 0.0) and np.all(result.estimate <= 1.0)
    # qiymət sərt datanın diapazonundan çox uzağa getməməlidir
    span = values.max() - values.min()
    assert result.estimate.min() > values.min() - 0.5 * span
    assert result.estimate.max() < values.max() + 0.5 * span


def test_case1_porosity_variance_is_smallest_near_the_wells():
    wells, grid = _sample_and_grid(50, seed=103)
    values = np.clip(0.18 + 0.03 * _field(wells, (300.0, 300.0), 104), 0.02, 0.40)
    result = interpolate_property_field(wells, values, grid, property_name="PORO")
    near = result.nearest_distance <= np.percentile(result.nearest_distance, 25)
    far = result.nearest_distance >= np.percentile(result.nearest_distance, 75)
    assert float(np.mean(result.variance[near])) < float(np.mean(result.variance[far]))


def test_case1_porosity_cross_validation_beats_the_global_mean():
    """Modelin heç bir dəyəri yoxdursa CV RMSE-si sadə qlobal ortadan
    yaxşı olmalıdır — bu, minimum məqbulluq həddidir."""
    wells, _ = _sample_and_grid(70, seed=105)
    values = np.clip(0.18 + 0.035 * _field(wells, (300.0, 300.0), 106), 0.02, 0.40)
    metrics = cross_validate_property(wells, values, resolve_strategy("PORO"),
                                      ValidationDesign())
    baseline_rmse = float(np.std(values))
    assert metrics.rmse < baseline_rmse
    assert metrics.r2 > 0.0


# ══ Hal 2 — LOQ-KEÇİRİCİLİK ═══════════════════════════════════════════
def test_case2_lognormal_permeability_is_recovered_in_log_space():
    wells, grid = _sample_and_grid(80, seed=201)
    log_truth = 4.5 + 1.3 * _field(wells, (350.0, 350.0), 202)
    values = np.exp(log_truth)

    result = interpolate_property_field(wells, values, grid, property_name="PERMX")
    assert np.all(result.estimate > 0.0)
    assert result.strategy.transform.name == "log"
    # çevrilmiş fəza həqiqətən ln(K)-dır və data diapazonundadır
    assert result.transformed_estimate.min() > log_truth.min() - 3.0
    assert result.transformed_estimate.max() < log_truth.max() + 3.0


def test_case2_median_and_mean_back_transforms_differ_as_theory_predicts():
    """`E[K]/median(K) = exp(σ²/2)` — geri çevirmə fərqi NƏZƏRİ düsturla
    uyğun gəlməlidir (B1.3)."""
    wells, grid = _sample_and_grid(60, seed=203)
    values = np.exp(4.0 + 1.2 * _field(wells, (300.0, 300.0), 204))

    median = interpolate_property_field(wells, values, grid, property_name="PERMX")
    mean = interpolate_property_field(
        wells, values, grid,
        strategy=resolve_strategy("PERMX").derive(back_transform=BackTransform.MEAN))

    ratio = mean.raw_estimate / median.raw_estimate
    expected = np.exp(0.5 * median.transformed_variance)
    assert np.allclose(ratio, expected, rtol=1e-9)
    assert float(np.mean(ratio)) > 1.0


def test_case2_log_space_cross_validation_beats_raw_space():
    """Loq fəzasının üstünlüyü DOĞRULAMA ilə təsdiqlənir."""
    wells, _ = _sample_and_grid(80, seed=205)
    values = np.exp(4.0 + 1.5 * _field(wells, (320.0, 320.0), 206))
    report = select_property_model(wells, values, default_candidates("PERMX"),
                                   "PERMX", ValidationDesign())
    assert report.selected is not None
    raw = next(r for r in report.results if "xam fəza" in r.candidate.label)
    assert report.selected.score <= raw.score


def test_case2_permeability_honours_wells_exactly():
    wells, _ = _sample_and_grid(40, seed=207)
    values = np.exp(4.0 + _field(wells, (300.0, 300.0), 208))
    result = interpolate_property_field(wells, values, wells, property_name="PERMX")
    assert np.allclose(result.estimate, values, rtol=1e-8)


# ══ Hal 3 — DOYMA ═════════════════════════════════════════════════════
def test_case3_saturation_stays_physically_bounded_everywhere():
    wells, grid = _sample_and_grid(60, seed=301)
    raw = 0.5 + 0.28 * _field(wells, (280.0, 280.0), 302)
    values = np.clip(raw, 0.0, 1.0)

    result = interpolate_property_field(wells, values, grid, property_name="SW")
    assert np.all(result.raw_estimate >= 0.0)
    assert np.all(result.raw_estimate <= 1.0)
    assert not np.any(result.bound_adjusted), "logit yolunda kəsməyə ehtiyac yoxdur"


def test_case3_saturation_extrapolation_remains_bounded():
    """Məlumat buludundan ÇOX UZAQDA da hədlər pozulmur."""
    wells, _ = _sample_and_grid(40, seed=303, high=600.0)
    values = np.clip(0.5 + 0.3 * _field(wells, (200.0, 200.0), 304), 0.0, 1.0)
    far = np.array([[-5000.0, -5000.0], [9000.0, 9000.0], [0.0, 9000.0]])
    result = interpolate_property_field(wells, values, far, property_name="SW")
    finite = np.isfinite(result.raw_estimate)
    assert np.all(result.raw_estimate[finite] >= 0.0)
    assert np.all(result.raw_estimate[finite] <= 1.0)


def test_case3_saturation_cross_validation_has_no_bound_violations():
    wells, _ = _sample_and_grid(60, seed=305)
    values = np.clip(0.45 + 0.3 * _field(wells, (280.0, 280.0), 306), 0.0, 1.0)
    metrics = cross_validate_property(wells, values, resolve_strategy("SW"),
                                      ValidationDesign())
    assert metrics.bound_violations == 0
    assert np.isfinite(metrics.rmse)


# ══ Hal 4 — NTG ═══════════════════════════════════════════════════════
def test_case4_ntg_is_bounded_and_cross_validates():
    """Real (hədlərdə YIĞILMAYAN) NTG sahəsi: logit yolu həm hədləri
    saxlayır, həm də qlobal ortadan yaxşı proqnoz verir."""
    wells, grid = _sample_and_grid(80, seed=401)
    values = np.clip(0.55 + 0.12 * _field(wells, (400.0, 400.0), 402), 0.0, 1.0)
    assert float(np.mean((values <= 1e-12) | (values >= 1.0 - 1e-12))) == 0.0

    result = interpolate_property_field(wells, values, grid, property_name="NTG")
    assert np.all(result.estimate >= 0.0) and np.all(result.estimate <= 1.0)

    metrics = cross_validate_property(wells, values, resolve_strategy("NTG"),
                                      ValidationDesign())
    assert metrics.rmse < float(np.std(values))
    assert metrics.r2 > 0.0
    assert metrics.bound_violations == 0


def test_case4_censored_ntg_is_a_known_limitation_that_model_selection_detects():
    """ÖLÇÜLMÜŞ MƏHDUDİYYƏT (gizlədilmir): dəyərlərin bir hissəsi DƏQİQ
    hədddə "yığılıbsa" (censored), logit onları eyni ekstremal nöqtəyə
    (`±ln((1−ε)/ε) ≈ ±9.2`) göndərir. Çevrilmiş sahə iki modallı və ağır
    quyruqlu olur, kriginq isə bunu pis idarə edir.

    Bu, çevirmənin SƏHV olduğunu DEYİL, HƏR DATAYA UYĞUN OLMADIĞINI
    göstərir — məhz buna görə `default_candidates()` hədli xassələr üçün
    də ÇEVİRMƏSİZ namizədi daxil edir və qərar DOĞRULAMAYA buraxılır
    (B3.4). Test bu davranışı QORUYUR.

    Ölçülmüş (60 quyu, 8% dəyər `NTG = 1.0`-da yığılıb, R²):
        logit + exponential  −0.120
        xam   + exponential  +0.128     ← çevirməsiz DAHA YAXŞI
    """
    wells, _ = _sample_and_grid(60, seed=401)
    values = np.clip(0.72 + 0.20 * _field(wells, (330.0, 330.0), 402), 0.0, 1.0)
    censored = float(np.mean(values >= 1.0 - 1e-12))
    assert censored > 0.05, "bu ssenari qəsdən 'yığılmış' datadır"

    from imex2d.geology.transforms import IDENTITY_TRANSFORM
    logit = cross_validate_property(
        wells, values, resolve_strategy("NTG").derive(variogram_model="exponential"),
        ValidationDesign())
    raw = cross_validate_property(
        wells, values,
        resolve_strategy("NTG").derive(variogram_model="exponential",
                                       transform=IDENTITY_TRANSFORM,
                                       variable_type=VariableType.CONTINUOUS),
        ValidationDesign())
    assert raw.rmse < logit.rmse, (
        "yığılmış datada çevirməsiz variant daha yaxşı olmalıdır")
    # hər iki halda fiziki hədlər POZULMUR (xam variantda kəsmə işə düşür)
    assert logit.bound_violations == 0


def test_case4_ntg_uncertainty_is_reported_as_a_delta_approximation():
    wells, grid = _sample_and_grid(45, seed=403)
    values = np.clip(0.7 + 0.2 * _field(wells, (300.0, 300.0), 404), 0.0, 1.0)
    result = interpolate_property_field(wells, values, grid, property_name="NTG")
    assert result.variance_kind.value == "delta"
    assert np.all(result.variance[np.isfinite(result.variance)] >= 0.0)


# ══ Hal 5 — FASİYA ════════════════════════════════════════════════════
def _facies_field(n_wells, seed, high=1200.0, ranges=(700.0, 400.0)):
    """İki fasiyalı sahə: qauss sahəsinin işarəsinə görə (məkanca
    davamlı zolaqlar, təsadüfi səpələnmə DEYİL)."""
    wells, grid = _sample_and_grid(n_wells, seed, high=high)
    everything = np.vstack([wells, grid])
    # Radiuslar quyu aralığından BÖYÜK seçilir: 150 m-lik zolaqları 80
    # quyu ilə bərpa etmək prinsipcə mümkün deyil (nümunələmə həddi,
    # kodun qüsuru DEYİL) — bax `test_case5_facies_accuracy_degrades_...`.
    latent = _field(everything, ranges, seed + 1)
    codes = (latent > 0.0).astype(int)
    return wells, grid, codes[:n_wells], codes[n_wells:]


def test_case5_facies_probabilities_are_valid_and_predictive():
    wells, grid, well_codes, grid_truth = _facies_field(90, seed=501)
    result = interpolate_categorical_field(wells, well_codes.astype(float), grid,
                                           property_name="FACIES")
    assert np.allclose(result.probabilities.sum(axis=1), 1.0)
    assert np.all(result.probabilities >= 0.0) and np.all(result.probabilities <= 1.0)
    accuracy = float(np.mean(result.most_probable == grid_truth))
    assert accuracy > 0.7, f"şəbəkə dəqiqliyi {accuracy:.3f}"
    assert accuracy > float(np.mean(grid_truth == np.bincount(grid_truth).argmax())), \
        "ən çox rast gələn kodu təxmin etməkdən yaxşı olmalıdır"


def test_case5_facies_confidence_tracks_prediction_accuracy():
    """Yüksək etimadlı hüceyrələr HƏQİQƏTƏN daha dəqiq olmalıdır —
    entropiya mənalı ölçüdür, bəzək deyil."""
    wells, grid, well_codes, grid_truth = _facies_field(90, seed=503)
    result = interpolate_categorical_field(wells, well_codes.astype(float), grid,
                                           property_name="FACIES")
    correct = result.most_probable == grid_truth
    confident = result.normalized_entropy <= np.percentile(result.normalized_entropy, 30)
    uncertain = result.normalized_entropy >= np.percentile(result.normalized_entropy, 70)
    assert float(np.mean(correct[confident])) > float(np.mean(correct[uncertain]))


def test_case5_facies_cross_validation_reports_classification_metrics():
    wells, _, well_codes, _ = _facies_field(90, seed=505)
    metrics = cross_validate_property(
        wells, well_codes.astype(float), resolve_strategy("FACIES"),
        ValidationDesign(kind=ValidationKind.RANDOM_KFOLD, k=5, seed=7))
    assert metrics.accuracy > 0.6
    assert metrics.brier_score < 0.5
    # Log-loss AĞIR QUYRUQLUDUR: indikator kriginq bəzi hüceyrələrdə
    # SƏHV kateqoriyaya demək olar 1 ehtimal verir və `−ln(p)` partlayır.
    # Ona görə müqayisə üçün BRIER (məhdud) işlədilir, log-loss isə yalnız
    # hesabatda saxlanılır — bu, metriklərin dürüst oxunmasıdır.
    assert np.isfinite(metrics.log_loss) and metrics.log_loss > 0.0
    assert metrics.confusion.sum() == metrics.n


def test_case5_facies_accuracy_degrades_when_structure_is_finer_than_well_spacing():
    """NÜMUNƏLƏMƏ HƏDDİ (kod qüsuru deyil): fasiya zolaqları quyu
    aralığından KİÇİK olanda heç bir interpolyator onları bərpa edə
    bilməz. Test bunu ÖLÇÜR ki, gələcəkdə "dəqiqlik aşağıdır" şikayəti
    səhv yerdə axtarılmasın."""
    _, _, coarse_codes, coarse_truth = _facies_field(90, seed=507,
                                                     ranges=(700.0, 400.0))
    wells_c, grid_c, codes_c, truth_c = _facies_field(90, seed=507,
                                                      ranges=(700.0, 400.0))
    wells_f, grid_f, codes_f, truth_f = _facies_field(90, seed=507,
                                                      ranges=(250.0, 90.0))
    coarse = interpolate_categorical_field(wells_c, codes_c.astype(float), grid_c,
                                           property_name="FACIES")
    fine = interpolate_categorical_field(wells_f, codes_f.astype(float), grid_f,
                                         property_name="FACIES")
    coarse_accuracy = float(np.mean(coarse.most_probable == truth_c))
    fine_accuracy = float(np.mean(fine.most_probable == truth_f))
    assert coarse_accuracy > fine_accuracy
    # incə sahədə qeyri-müəyyənlik də DAHA BÖYÜK olmalıdır (dürüst siqnal)
    assert float(np.mean(fine.normalized_entropy)) > \
        float(np.mean(coarse.normalized_entropy))


# ══ Hal 6 — ANİZOTROP REZERVUAR ═══════════════════════════════════════
def test_case6_anisotropic_kriging_beats_isotropic_on_a_layered_reservoir():
    """Güclü üfüqi, zəif şaquli davamlılıq — deterministik kriginq."""
    from imex2d.geology.anisotropy import AnisotropyParams
    rng = np.random.default_rng(601)
    points = np.column_stack([rng.uniform(0.0, 1500.0, size=(120, 2)),
                              rng.uniform(0.0, 40.0, 120)])
    values = 0.2 + 0.03 * _field(points, (1000.0, 1000.0, 15.0), 602)

    def loo_rmse(strategy):
        return cross_validate_property(points, values, strategy,
                                       ValidationDesign()).rmse

    base = resolve_strategy("PORO").derive(variogram_model="exponential")
    isotropic = base.derive(anisotropy=AnisotropyParams(
        range_major=1000.0, range_minor=1000.0, range_vertical=1000.0))
    layered = base.derive(anisotropy=AnisotropyParams(
        range_major=1000.0, range_minor=1000.0, range_vertical=15.0))
    assert loo_rmse(layered) < loo_rmse(isotropic)


def test_case6_sgs_reproduces_directional_continuity():
    """Eyni anizotropluq STOXASTİK yolda da işləyir (SGS)."""
    rng = np.random.default_rng(603)
    wells = rng.uniform(0.0, 1500.0, size=(45, 2))
    values = 0.2 + 0.03 * _field(wells, (800.0, 120.0), 604)
    axis = np.linspace(0.0, 1500.0, 26)
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    grid = np.column_stack([xx.ravel(), yy.ravel()])
    targets = np.vstack([wells, grid])

    realization = simulate_sgs_ensemble(
        1, wells, values, targets, base_seed=605, max_neighbors=20,
        variogram=PropertyVariogramParams(model="spherical", nugget=0.0,
                                          range_=800.0, range_minor=120.0,
                                          azimuth_deg=90.0)).realizations[0]
    simulated = realization.values[len(values):]
    along = experimental_variogram(grid, simulated, n_lags=8, max_lag=600.0,
                                   azimuth_deg=90.0, azimuth_tolerance_deg=20.0)
    across = experimental_variogram(grid, simulated, n_lags=8, max_lag=600.0,
                                    azimuth_deg=0.0, azimuth_tolerance_deg=20.0)
    assert along.gamma[4] < across.gamma[4]


def test_case6_sgs_honours_wells_in_the_anisotropic_reservoir():
    rng = np.random.default_rng(606)
    wells = rng.uniform(0.0, 1200.0, size=(35, 2))
    values = 0.2 + 0.03 * _field(wells, (700.0, 150.0), 607)
    axis = np.linspace(0.0, 1200.0, 20)
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    targets = np.vstack([wells, np.column_stack([xx.ravel(), yy.ravel()])])
    variogram = PropertyVariogramParams(model="spherical", nugget=0.0, range_=700.0,
                                        range_minor=150.0, azimuth_deg=90.0)
    realization = simulate_sgs_ensemble(1, wells, values, targets,
                                        variogram=variogram, base_seed=608,
                                        max_neighbors=20).realizations[0]
    target_model = VariogramParameters("spherical", 0.0, float(np.var(values)),
                                       700.0, 0.0, 0)
    report = validate_realization(realization, wells, values, targets,
                                  target_variogram=target_model)
    assert report.hard_data_honored
    assert report.ks_statistic < 0.3


# ══ Hal 7 — SEYRƏK QUYULAR ════════════════════════════════════════════
def test_case7_uncertainty_grows_with_distance_from_the_wells():
    wells = np.array([[600.0, 600.0], [650.0, 600.0], [600.0, 650.0],
                      [650.0, 650.0], [625.0, 625.0]])
    values = np.array([0.20, 0.21, 0.19, 0.22, 0.205])
    distances = np.array([0.0, 100.0, 400.0, 1500.0, 6000.0])
    targets = np.column_stack([625.0 + distances, np.full(5, 625.0)])

    result = interpolate_property_field(wells, values, targets,
                                        property_name="PORO",
                                        strategy=resolve_strategy("PORO").derive(
                                            variogram_model="spherical"))
    finite = np.isfinite(result.variance)
    assert np.all(np.diff(result.variance[finite]) >= -1e-12), (
        "varians məsafə ilə MONOTON artmalıdır")
    assert result.confidence[-1] == Confidence.EXTRAPOLATED.value
    assert result.data_density[0] > result.data_density[-1]


def test_case7_sparse_wells_produce_wider_ensemble_spread_far_from_data():
    wells = np.array([[300.0, 300.0], [900.0, 300.0], [300.0, 900.0],
                      [900.0, 900.0], [600.0, 600.0]])
    values = np.array([0.18, 0.22, 0.19, 0.23, 0.205])
    axis = np.linspace(0.0, 1200.0, 18)
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    grid = np.column_stack([xx.ravel(), yy.ravel()])
    targets = np.vstack([wells, grid])

    ensemble = simulate_sgs_ensemble(
        8, wells, values, targets, base_seed=701, max_neighbors=8,
        variogram=PropertyVariogramParams(model="spherical", nugget=0.0,
                                          range_=400.0))
    spread = ensemble.std[len(values):]
    from scipy.spatial import cKDTree
    distance, _ = cKDTree(wells).query(grid, k=1)
    near = distance <= np.percentile(distance, 25)
    far = distance >= np.percentile(distance, 75)
    assert float(np.mean(spread[far])) > float(np.mean(spread[near]))


def test_case7_extrapolated_cells_are_labelled_and_not_hidden():
    wells = np.array([[500.0, 500.0], [520.0, 500.0], [500.0, 520.0],
                      [520.0, 520.0], [510.0, 510.0]])
    values = np.array([0.20, 0.21, 0.19, 0.22, 0.205])
    axis = np.linspace(0.0, 4000.0, 12)
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    grid = np.column_stack([xx.ravel(), yy.ravel()])
    result = interpolate_property_field(wells, values, grid, property_name="PORO")
    assert np.any(result.extrapolated)
    labels = np.asarray(result.confidence).astype(str)
    assert np.sum(labels == Confidence.EXTRAPOLATED.value) == int(
        np.sum(result.extrapolated))
