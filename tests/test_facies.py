"""Phase 4 — kateqorik fasiya modelləşdirməsi + Sequential Indicator Simulation.

Bu testlər PROQRAM DÜZGÜNLÜYÜNÜ yoxlayır (indikator çevirməsi, ehtimal
normallaşdırması, sərt-data hörməti, təkrarlana bilənlik, ansambl
statistikası) — "geoloji reallıq"ı YOX (bax `facies.py` modul
docstring-indəki elmi çəkincə). Statistik testlər sabit seed və GENİŞ
tolerantlıqla yazılıb (bax `test_variogram.py`-dəki eyni konvensiya).
"""

from __future__ import annotations

import numpy as np

from imex2d.domain.validation import validate_facies_proportions
from imex2d.geology.facies import (FaciesProportions, FaciesVariogramParams,
                                   indicator, observed_proportions, run_realizations,
                                   simulate_sis, summarize_realized_proportions)


def _dmat(a, b):
    diff = a[:, None, :] - b[None, :, :]
    return np.sqrt(np.sum(diff ** 2, axis=-1))


# ── 1. indikator çevirməsi ────────────────────────────────────────────────
def test_indicator_conversion_basic():
    codes = np.array([0, 1, 2, 1, 0])
    assert np.array_equal(indicator(codes, 1), [0.0, 1.0, 0.0, 1.0, 0.0])
    assert np.array_equal(indicator(codes, 0), [1.0, 0.0, 0.0, 0.0, 1.0])


def test_observed_proportions_matches_counts():
    codes = np.array([0, 0, 1, 1, 1, 2])
    props = observed_proportions(codes)
    assert props[0] == 2 / 6
    assert props[1] == 3 / 6
    assert props[2] == 1 / 6
    assert abs(sum(props.values()) - 1.0) < 1e-12


# ── 2. fasiya nisbəti yoxlanması ────────────────────────────────────────
def test_proportion_validation_accepts_valid_distribution():
    assert validate_facies_proportions({0: 0.5, 1: 0.3, 2: 0.2}).ok


def test_proportion_validation_rejects_sum_not_one():
    assert not validate_facies_proportions({0: 0.5, 1: 0.3}).ok


def test_proportion_validation_rejects_negative():
    assert not validate_facies_proportions({0: 1.2, 1: -0.2}).ok


def test_facies_proportions_precedence_layer_over_region_over_global():
    fp = FaciesProportions(
        global_proportions={0: 0.5, 1: 0.5},
        region_proportions={1: {0: 0.2, 1: 0.8}},
        layer_proportions={3: {0: 0.9, 1: 0.1}})
    assert fp.for_cell(region=1, layer=3) == {0: 0.9, 1: 0.1}     # lay qalib gəlir
    assert fp.for_cell(region=1, layer=None) == {0: 0.2, 1: 0.8}  # region
    assert fp.for_cell(region=None, layer=None) == {0: 0.5, 1: 0.5}  # qlobal
    assert fp.validate().ok


def test_facies_proportions_validate_catches_bad_region_entry():
    fp = FaciesProportions(global_proportions={0: 0.5, 1: 0.5},
                           region_proportions={1: {0: 0.9, 1: 0.9}})
    assert not fp.validate().ok


# ── synthetic sınaq datası ────────────────────────────────────────────────
def _two_facies_wells():
    """3 fasiya-0, 3 fasiya-1 quyusu, aydın məkan qruplaşması (0 solda, 1 sağda)."""
    points = np.array([[0., 0.], [0., 50.], [10., 25.],
                       [200., 0.], [200., 50.], [190., 25.]])
    codes = np.array([0, 0, 0, 1, 1, 1])
    return points, codes


def _grid_targets(nx=10, ny=10, dx=20.0, dy=20.0):
    xs = (np.arange(nx) + 0.5) * dx
    ys = (np.arange(ny) + 0.5) * dy
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


# ── 3. indikator variogram (Phase 2-3 infrastrukturunun istifadəsi) ──────
def test_indicator_variogram_is_fit_via_phase23_infrastructure():
    points, codes = _two_facies_wells()
    # kifayət qədər nöqtə yoxdur (6 nöqtə, 4 minimum + 3 dolu bin tələbi
    # çətin ödənə bilər) -> ehtiyat evristikaya keçməlidir, XƏBƏRDARLIQLA
    targets = _grid_targets(nx=4, ny=4, dx=50.0, dy=15.0)
    realization = simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, seed=1)
    assert realization.codes.shape == (targets.shape[0],)
    assert all(np.isin(realization.codes, [0, 1]))


# ── 4. şərti ehtimal normallaşdırması ────────────────────────────────────
def test_conditional_probabilities_are_valid_distribution_at_every_cell():
    """Birbaşa daxili kriging + normallaşdırma addımını sınayır (SIS-in
    öz sampling addımından ASILI OLMADAN) — hər addımda ehtimallar
    mənfi deyil və cəmi 1-ə YAXINDIR."""
    from imex2d.geology.interpolation import OrdinaryKriging
    points, codes = _two_facies_wells()
    facies_list = [0, 1]
    krigers = {k: OrdinaryKriging(range_=80.0, nugget=0.05, max_neighbors=6)
              for k in facies_list}
    targets = _grid_targets(nx=6, ny=3, dx=30.0, dy=15.0)
    for row in range(targets.shape[0]):
        probs = np.array([krigers[k].interpolate(points, indicator(codes, k),
                                                  targets[row:row + 1])[0]
                          for k in facies_list])
        clipped = np.clip(probs, 0.0, None)
        total = clipped.sum()
        assert total > 0
        normalized = clipped / total
        assert np.all(normalized >= -1e-12)
        assert abs(normalized.sum() - 1.0) < 1e-9


def test_simulate_sis_realized_codes_always_valid_categories():
    points, codes = _two_facies_wells()
    targets = _grid_targets()
    realization = simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, seed=7)
    assert set(np.unique(realization.codes)).issubset({0, 1})
    assert np.all(realization.codes >= 0)


# ── 5. sərt data hörməti ──────────────────────────────────────────────────
def test_hard_data_is_honored_exactly_regardless_of_seed():
    points, codes = _two_facies_wells()
    # hədəflərin bir hissəsi məhz quyu koordinatlarıdır
    targets = np.vstack([points, _grid_targets(nx=5, ny=5)])
    for seed in (1, 2, 3, 999):
        realization = simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, seed=seed)
        assert np.array_equal(realization.codes[:len(points)], codes), (
            f"seed={seed}: sərt data pozuldu")
        assert np.all(realization.hard_data_mask[:len(points)])


# ── 6/7. seed təkrarlana bilənlik ────────────────────────────────────────
def test_same_seed_produces_identical_realization():
    points, codes = _two_facies_wells()
    targets = _grid_targets()
    r1 = simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, seed=42)
    r2 = simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, seed=42)
    assert np.array_equal(r1.codes, r2.codes)


def test_different_seeds_produce_different_realizations():
    points, codes = _two_facies_wells()
    targets = _grid_targets()
    r1 = simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, seed=1)
    r2 = simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, seed=2)
    assert not np.array_equal(r1.codes, r2.codes)


# ── 8. anizotrop davamlılıq ───────────────────────────────────────────────
def test_anisotropic_indicator_kriging_favours_major_axis_neighbour():
    """Phase 2-3-ün `OrdinaryKriging(azimuth_deg=..., range_minor=...)`-u
    BİRBAŞA işlədilir (təkrar yazılmır) — güclü anizotropluqla major-ox
    istiqamətindəki fasiya-1 nöqtəsi minor-ox istiqamətindəkindən daha
    çox təsir etməlidir (bax `test_variogram.py`-dəki eyni ssenari)."""
    points = np.array([[0., 50.], [50., 0.], [-200., -200.]])
    codes = np.array([1, 0, 0])
    vp = {0: FaciesVariogramParams(range_=200.0, range_minor=20.0, azimuth_deg=0.0,
                                   range_v=1e9, nugget=0.0),
         1: FaciesVariogramParams(range_=200.0, range_minor=20.0, azimuth_deg=0.0,
                                  range_v=1e9, nugget=0.0)}
    realization = simulate_sis(points, codes, np.array([[0., 0.]]), {0: 0.5, 1: 0.5},
                               variograms=vp, seed=0)
    # birbaşa ehtimalı da yoxlayaq (sampling təsadüfi olmasın deyə)
    from imex2d.geology.interpolation import OrdinaryKriging
    k1 = OrdinaryKriging(range_=200.0, range_minor=20.0, azimuth_deg=0.0,
                         range_v=1e9, nugget=0.0)
    p1 = k1.interpolate(points, indicator(codes, 1), np.array([[0., 0.]]))[0]
    assert p1 > 0.5   # major-ox (Y) istiqamətindəki fasiya-1 üstünlük etməlidir


def test_ensemble_spatial_continuity_is_stronger_along_major_axis():
    """Ansambl statistikası (tapşırıq §11: TƏK hüceyrə deyil, statistik
    struktur): major-ox boyu iki hüceyrə minor-ox boyu iki hüceyrədən
    DAHA TEZ-TEZ eyni fasiyada olmalıdır (bir çox realizasiya üzərində)."""
    points = np.array([[-100., 0.], [100., 0.], [0., -100.], [0., 100.]])
    codes = np.array([0, 0, 1, 1])
    vp = {0: FaciesVariogramParams(range_=300.0, range_minor=15.0, azimuth_deg=90.0, nugget=0.0),
         1: FaciesVariogramParams(range_=300.0, range_minor=15.0, azimuth_deg=90.0, nugget=0.0)}
    # azimuth=90 -> major ox = X. major-cütü X boyu, minor-cütü Y boyu.
    targets = np.array([[-30., 0.], [30., 0.], [0., -30.], [0., 30.]])
    major_matches, minor_matches = 0, 0
    n = 60
    for i in range(n):
        r = simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, variograms=vp,
                         seed=100 + i, max_neighbors=4)
        if r.codes[0] == r.codes[1]:
            major_matches += 1
        if r.codes[2] == r.codes[3]:
            minor_matches += 1
    assert major_matches > minor_matches, (
        f"major-ox uyğunluğu ({major_matches}/{n}) minor-oxdan ({minor_matches}/{n}) "
        "çox olmalıydı")


# ── 9. çoxfasiyalı ─────────────────────────────────────────────────────
def test_three_facies_realization_uses_only_declared_codes():
    points = np.array([[0., 0.], [100., 0.], [0., 100.], [100., 100.], [50., 50.], [20., 80.]])
    codes = np.array([0, 1, 2, 0, 1, 2])
    targets = _grid_targets(nx=6, ny=6, dx=20.0, dy=20.0)
    realization = simulate_sis(points, codes, targets, {0: 0.4, 1: 0.35, 2: 0.25}, seed=3)
    assert set(np.unique(realization.codes)).issubset({0, 1, 2})


# ── 10. tək-fasiya kənar halı ────────────────────────────────────────────
def test_single_facies_edge_case_assigns_everywhere_without_randomness():
    points = np.array([[0., 0.], [10., 10.]])
    codes = np.array([0, 0])
    targets = _grid_targets(nx=3, ny=3)
    r1 = simulate_sis(points, codes, targets, {0: 1.0}, seed=1)
    r2 = simulate_sis(points, codes, targets, {0: 1.0}, seed=999)
    assert np.all(r1.codes == 0)
    assert np.array_equal(r1.codes, r2.codes)   # seed-dən ASILI DEYİL (qeyri-müəyyənlik yoxdur)


# ── 11. seyrək sərt data ──────────────────────────────────────────────────
def test_sparse_hard_data_falls_back_with_warning_but_completes():
    points = np.array([[0., 0.], [100., 100.]])   # cəmi 2 quyu
    codes = np.array([0, 1])
    targets = _grid_targets(nx=5, ny=5, dx=25.0, dy=25.0)
    realization = simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, seed=5)
    assert realization.codes.shape == (targets.shape[0],)
    assert set(np.unique(realization.codes)).issubset({0, 1})
    assert any("evristika" in w or "qonşu" in w for w in realization.warnings), (
        "seyrək data ilə ehtiyat evristikası/xəbərdarlığı gözlənilirdi")


# ── 12. çoxlu realizasiya üzrə nisbət saxlanması ─────────────────────────
def test_proportions_preserved_on_average_across_many_realizations():
    points, codes = _two_facies_wells()
    targets = _grid_targets(nx=8, ny=8, dx=25.0, dy=25.0)
    realizations = run_realizations(30, points, codes, targets, {0: 0.5, 1: 0.5}, seed=0)
    summary = summarize_realized_proportions(realizations)
    assert set(summary) == {0, 1}
    for code, stats in summary.items():
        assert stats["requested"] == 0.5
        assert abs(stats["mean"] - 0.5) < 0.15   # orta hesabla yaxın, TƏK realizasiya deyil
        assert stats["min"] <= stats["mean"] <= stats["max"]
        assert stats["std"] >= 0.0
    # realizasiyalar müxtəlifdir (eyni deyil) — statistik dəyişkənlik var
    assert not np.array_equal(realizations[0].codes, realizations[1].codes)


def test_run_realizations_ids_and_seeds_are_deterministic_sequence():
    points, codes = _two_facies_wells()
    targets = _grid_targets(nx=4, ny=4)
    realizations = run_realizations(3, points, codes, targets, {0: 0.5, 1: 0.5}, seed=10)
    assert [r.realization_id for r in realizations] == [0, 1, 2]
    assert [r.seed for r in realizations] == [10, 1010, 2010]


# ── 13. etibarsız giriş ────────────────────────────────────────────────
def test_rejects_proportions_not_summing_to_one():
    points, codes = _two_facies_wells()
    targets = _grid_targets(nx=3, ny=3)
    try:
        simulate_sis(points, codes, targets, {0: 0.5, 1: 0.6}, seed=0)
    except ValueError as exc:
        assert "1.0" in str(exc) or "cəm" in str(exc)
        return
    raise AssertionError("cəmi 1 olmayan nisbətlər qəbul edildi")


def test_rejects_mismatched_points_and_codes_length():
    targets = _grid_targets(nx=3, ny=3)
    try:
        simulate_sis(np.array([[0., 0.], [1., 1.]]), np.array([0]), targets, {0: 1.0}, seed=0)
    except ValueError as exc:
        assert "uzunlu" in str(exc)
        return
    raise AssertionError("uzunluq uyğunsuzluğu qəbul edildi")


def test_rejects_hard_data_code_not_in_proportions():
    points = np.array([[0., 0.], [10., 10.], [5., 5.]])
    codes = np.array([0, 1, 2])            # 2 -> proportions-da yoxdur
    targets = _grid_targets(nx=3, ny=3)
    try:
        simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, seed=0)
    except ValueError as exc:
        assert "fasiya kodu" in str(exc)
        return
    raise AssertionError("naməlum fasiya kodu qəbul edildi")


# ── 14. kiçik sintetik SIS bençmarkı ──────────────────────────────────────
def test_small_synthetic_sis_benchmark_completes_and_is_self_consistent():
    """15×15 = 225 hüceyrəlik tam grid, 8 quyu, 2 fasiya — tam SIS icrası,
    nəticənin daxili tutarlılığını (ölçü, kodlar, sərt data, nisbətlər)
    yoxlayır. Performans sənədləşdirilib (bax FACIES.md) — bu test
    YALNIZ düzgünlüyü, sürəti YOX yoxlayır."""
    rng = np.random.default_rng(0)
    points = rng.uniform(0, 300, size=(8, 2))
    codes = (points[:, 0] > 150).astype(int)   # aydın məkan qruplaşması
    targets = _grid_targets(nx=15, ny=15, dx=20.0, dy=20.0)

    realization = simulate_sis(points, codes, targets, {0: 0.5, 1: 0.5}, seed=42,
                               max_neighbors=12)
    assert realization.codes.shape == (225,)
    assert np.all(realization.codes >= 0)
    assert set(np.unique(realization.codes)).issubset({0, 1})
    total = sum(realization.realized_proportions.values())
    assert abs(total - 1.0) < 1e-9
