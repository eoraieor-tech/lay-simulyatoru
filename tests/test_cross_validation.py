"""M4 — cross-validation: real dəqiqliyi ölçür, 100% vəd etmir."""

from __future__ import annotations

import numpy as np

from imex2d.application.geology_service import (GeologicalGridSpec,
                                                 WellBasedGeologicalModelBuilder,
                                                 format_cross_validation_report)
from imex2d.domain.well_data import WellDataset, WellSample
from imex2d.geology.cross_validation import k_fold, leave_one_out
from imex2d.geology.interpolation import InverseDistance, NearestNeighbour, OrdinaryKriging


def test_loo_rejects_fewer_than_three_points():
    try:
        leave_one_out(InverseDistance(), np.array([[0., 0.], [1., 1.]]),
                      np.array([1.0, 2.0]))
    except ValueError as error:
        assert "3 nöqtə" in str(error)
        return
    raise AssertionError("2 nöqtə ilə ValueError gözlənilirdi")


def test_loo_recovers_a_smooth_linear_trend_well():
    """Səthin özü hamar (xətti) olanda LOO xətası kiçik, R² yüksək olmalıdır."""
    rng = np.random.default_rng(0)
    xy = rng.uniform(0, 100, size=(12, 2))
    values = 0.15 + 0.001 * xy[:, 0] + 0.0005 * xy[:, 1]   # aydın xətti trend
    result = leave_one_out(InverseDistance(power=2.0), xy, values)
    assert result.rmse < 0.02
    assert result.r2 > 0.5


def test_poor_model_reports_low_or_negative_r2_without_hiding_it():
    """Az və nizamsız nöqtə + ən yaxın qonşu ilə pis proqnoz — R² gizlədilmir,
    mənfi ola bilər (bax tələb: '100% dəqiq' vəd edilmir)."""
    rng = np.random.default_rng(1)
    xy = rng.uniform(0, 100, size=(6, 2))
    values = rng.uniform(0, 1000, size=6)   # heç bir məkan korrelyasiyası yoxdur
    result = leave_one_out(NearestNeighbour(), xy, values)
    assert result.r2 < 0.8, "tamamilə təsadüfi datada süni yüksək R² şübhəlidir"
    # nəticə nə olursa olsun as_text()-də görünməlidir, gizlədilmir
    assert "R²" in result.as_text()


def test_mape_safely_excludes_near_zero_actual_values():
    """Ölçülmüş dəyərlərdən biri ~0-dırsa MAPE sıfıra bölünmə ilə çökməməli,
    o nöqtəni çıxarıb sayını bildirməlidir."""
    xy = np.array([[0., 0.], [10., 0.], [0., 10.], [10., 10.], [5., 5.], [3., 7.]])
    values = np.array([0.0, 1.0, 2.0, 3.0, 2.5, 1.5])
    result = leave_one_out(InverseDistance(), xy, values)
    assert result.mape_excluded >= 1
    assert result.mape is not None and np.isfinite(result.mape)


def test_r2_is_nan_not_crashing_when_all_values_identical():
    xy = np.array([[0., 0.], [10., 0.], [0., 10.], [10., 10.], [5., 5.]])
    values = np.full(5, 0.2)
    result = leave_one_out(InverseDistance(), xy, values)
    assert np.isnan(result.r2)
    assert any("R²" in w for w in result.warnings)


def test_k_fold_falls_back_to_loo_when_too_few_points_for_requested_k():
    xy = np.array([[0., 0.], [10., 0.], [0., 10.], [10., 10.]])
    values = np.array([0.1, 0.2, 0.15, 0.25])
    result = k_fold(InverseDistance(), xy, values, k=5)
    assert result.method == "leave-one-out"
    assert any("leave-one-out" in w for w in result.warnings)


def test_k_fold_runs_with_enough_points():
    rng = np.random.default_rng(2)
    xy = rng.uniform(0, 100, size=(20, 2))
    values = 100.0 + 0.5 * xy[:, 0]
    result = k_fold(OrdinaryKriging(nugget=0.01), xy, values, k=4)
    assert result.method == "4-fold"
    assert result.n_points == 20


def test_log_space_metrics_populated_for_permeability_like_data():
    """PERMX kimi log-normal keçiricilik: `compute_log_metrics=True` +
    `log_transform=True` ilə HƏM mD, HƏM log(mD) fəzasında metrik çıxmalıdır."""
    rng = np.random.default_rng(3)
    xy = rng.uniform(0, 100, size=(10, 2))
    permx = 50.0 * np.exp(0.02 * xy[:, 0])   # log-normal xarakterli
    result = leave_one_out(OrdinaryKriging(nugget=0.0), xy, permx,
                           log_transform=True, compute_log_metrics=True)
    assert result.rmse_log is not None and np.isfinite(result.rmse_log)
    assert result.mae_log is not None
    assert "log fəzasında" in result.as_text()


def test_search_radius_nan_predictions_are_dropped_not_fabricated():
    """`OrdinaryKriging(search_radius=...)` bəzi hədəflər üçün NaN
    qaytara bilər (yerli axtarışda uyğun nöqtə tapılmayanda) — CV bunları
    fabrikasiya etmədən metrikdən çıxarmalı və xəbərdarlıq verməlidir."""
    xy = np.array([[0., 0.], [1., 0.], [0., 1.], [1., 1.], [500., 500.]])
    values = np.array([0.1, 0.12, 0.11, 0.13, 0.9])   # son nöqtə təcrid olunub
    interpolator = OrdinaryKriging(nugget=0.0, search_radius=5.0, min_neighbors=1)
    result = leave_one_out(interpolator, xy, values)
    assert result.n_points < 5, "təcrid olunmuş nöqtənin NaN proqnozu metrikdən düşməli idi"
    assert any("NaN" in w for w in result.warnings)


# ── inteqrasiya: geology_service.py boru xəttindən keçərək ─────────────
def _layered_dataset_for_cv():
    samples = []
    xy = {"A": (25., 25.), "B": (125., 25.), "C": (25., 125.), "D": (125., 125.)}
    permx_by_layer = {0: 500.0, 1: 100.0, 2: 10.0}
    for name, (x, y) in xy.items():
        for k, permx in permx_by_layer.items():
            samples.append(WellSample(well=name, x=x, y=y, layer=k,
                                      values={"PERMX": permx * (1.0 + 0.05 * hash(name) % 3),
                                             "PORO": 0.1 + 0.05 * k}))
    return WellDataset(samples=samples, source="test")


def test_builder_cross_validate_runs_per_layer_independently():
    dataset = _layered_dataset_for_cv()
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    results, skipped = builder.cross_validate(dataset, "PERMX")
    assert set(results.keys()) == {0, 1, 2}
    assert skipped == {}
    for result in results.values():
        assert result.rmse_log is not None, "PERMX log_transform=True -> log metrikləri gəlməlidir"


def test_builder_cross_validate_skips_layers_with_too_few_points_with_a_reason():
    samples = [
        WellSample(well="A", x=25., y=25., layer=0, values={"PERMX": 500.0}),
        WellSample(well="B", x=125., y=25., layer=0, values={"PERMX": 480.0}),
        # K=1: cəmi 1 nöqtə - CV üçün kifayət etmir
        WellSample(well="A", x=25., y=25., layer=1, values={"PERMX": 100.0}),
    ]
    dataset = WellDataset(samples=samples, source="test")
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    results, skipped = builder.cross_validate(dataset, "PERMX")
    assert 0 not in results and 0 in skipped   # 2 nöqtə də 3-dən azdır
    assert 1 not in results and 1 in skipped


def test_cross_validate_all_and_report_format_cover_poro_and_permx():
    dataset = _layered_dataset_for_cv()
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    all_results = builder.cross_validate_all(dataset)
    assert set(all_results.keys()) == {"PORO", "PERMX"}
    text = format_cross_validation_report(all_results)
    assert "POROSITY:" in text
    assert "PERMEABILITY (PERMX):" in text
    assert "RMSE" in text and "MAE" in text and "R²" in text


def test_empty_cross_validation_report_does_not_crash():
    assert "tap" in format_cross_validation_report({}).lower() or format_cross_validation_report({})
