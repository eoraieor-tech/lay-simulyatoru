"""M2 — Kriging-in 3D/anizotrop, axtarış-radiuslu genişlənməsi.

Fon: M0/M1 (bax `test_layer_aware_kriging_leak.py`) laylar arasındakı
sızmanı bağladı, amma bunu edərkən boş lay üçün "açıq razılıqla"
(`allow_cross_layer_fallback=True`) digər laylardan borc alanda HAMISINI
BƏRABƏR çəki ilə hovuzlayırdı — yaxın lay ilə uzaq lay arasında fərq
qoymurdu. M2 bunu düzəldir: `OrdinaryKriging` indi əsl 3D (X,Y,Z) və
anizotrop (üfüqi/şaquli range ayrı) işləyə bilir, `geology_service.py`
hər nümunənin Z-sini (ölçülmüş `depth` və ya öz layının orta dərinliyi)
hədəfin öz K-sının həqiqi dərinliyi ilə müqayisə edir. Nəticə: boş layın
qonşu laydan borc alması artıq GEOLOJİ CƏHƏTDƏN ƏSASLANDIRILMIŞDIR (yaxın
lay çox, uzaq lay az təsir edir) — kor-koranə bərabər hovuzlama YOX.

Kriging riyaziyyatının özü (sferik variogram, Laqranj sistemi) DƏYİŞMƏYİB
— yalnız neçə ölçülü nöqtə qəbul etdiyi və (istəyə görə) hər hədəf üçün
YEREL sistemi qurub-qurmadığı dəyişib. `range_v`/`search_radius`/
`max_neighbors` heç biri verilməyəndə nəticə köhnə 2D qlobal kriging ilə
BİRƏBİR eynidir (bax `test_default_parameters_reproduce_pre_m2_2d_behaviour`).
"""

from __future__ import annotations

import numpy as np

from imex2d.application.geology_service import (GeologicalGridSpec,
                                                 WellBasedGeologicalModelBuilder)
from imex2d.domain.well_data import WellDataset, WellSample
from imex2d.geology.interpolation import OrdinaryKriging

POINTS = np.array([[0., 0.], [100., 0.], [0., 100.], [100., 100.], [50., 50.]])
VALUES = np.array([0.15, 0.25, 0.20, 0.30, 0.22])


# ── vahid səviyyəli: birbaşa OrdinaryKriging üzərində ───────────────────
def test_default_parameters_reproduce_pre_m2_2d_behaviour():
    """`range_v`/`search_radius`/`max_neighbors` verilməyəndə (defolt)
    nəticə M2-dən ƏVVƏLKİ 2D qlobal kriging ilə ədəd-ədəd eynidir."""
    result = OrdinaryKriging(nugget=0.05).interpolate(
        POINTS, VALUES, np.array([[25., 25.], [75., 75.]]))
    # 2D nöqtələr veriləndə daxili Z=0 ilə eyni nəticəni verməlidir
    points3d = np.column_stack([POINTS, np.zeros(5)])
    result3d = OrdinaryKriging(nugget=0.05).interpolate(
        points3d, VALUES, np.array([[25., 25., 0.], [75., 75., 0.]]))
    assert np.allclose(result, result3d)


def test_exact_interpolator_property_preserved_in_3d():
    points3d = np.column_stack([POINTS, np.array([10., 20., 15., 25., 18.])])
    result = OrdinaryKriging(nugget=0.0).interpolate(points3d, VALUES, points3d)
    assert np.allclose(result, VALUES), "nugget=0 ilə öz nöqtələrində dəqiq olmalıdır"


def test_large_vertical_range_makes_far_depth_nearly_isotropic_with_surface():
    """`range_v` çox böyükdürsə (demək olar limitsiz), Z fərqi nəticəyə
    demək olar təsir etməməlidir — anizotropluq SÖNÜR."""
    points3d = np.column_stack([POINTS, np.zeros(5)])
    target_shallow = np.array([[50., 50., 0.0]])
    target_deep = np.array([[50., 50., 5000.0]])
    near = OrdinaryKriging(nugget=0.0, range_v=1e9).interpolate(points3d, VALUES, target_shallow)
    far = OrdinaryKriging(nugget=0.0, range_v=1e9).interpolate(points3d, VALUES, target_deep)
    assert abs(near[0] - far[0]) < 1e-3


def test_small_vertical_range_makes_far_depth_drift_towards_global_mean():
    """`range_v` kiçikdirsə, uzaq Z-dəki hədəf lokal nöqtələrdən deyil,
    (kriging-in qərəzsizlik şərtinə görə) qlobal ortaya yaxınlaşır."""
    points3d = np.column_stack([POINTS, np.zeros(5)])
    target_deep = np.array([[50., 50., 5000.0]])
    far = OrdinaryKriging(nugget=0.0, range_v=10.0).interpolate(points3d, VALUES, target_deep)
    assert abs(far[0] - VALUES.mean()) < 1e-2


def test_search_radius_returns_nan_when_nothing_within_reach():
    result = OrdinaryKriging(nugget=0.0, search_radius=5.0, min_neighbors=1).interpolate(
        POINTS, VALUES, np.array([[500., 500.]]))
    assert np.isnan(result[0]), "radiusdan kənar hədəf üçün dəyər UYDURULMAMALIDIR"


def test_max_neighbors_limits_local_system_size():
    """`max_neighbors=1` təcrübədə ən yaxın qonşuya bərabər olmalıdır."""
    nn_like = OrdinaryKriging(nugget=0.0, max_neighbors=1).interpolate(
        POINTS, VALUES, np.array([[1., 1.]]))
    assert abs(nn_like[0] - 0.15) < 1e-9   # (0,0) nöqtəsinin dəyəri


# ── inteqrasiya: geology_service.py boru xəttindən keçərək ─────────────
_WELL_XY = {"A": (25.0, 25.0), "B": (125.0, 25.0), "C": (25.0, 125.0)}


def _spec(nz):
    return GeologicalGridSpec(nx=3, ny=3, nz=nz, dx=50.0, dy=50.0, dz=10.0,
                              top_depth=2000.0)


def _layered_dataset_with_explicit_depth(depth_offset_a_layer0: float = 0.0):
    """3 quyu x 3 lay (K=0,1,2), K=3 laysız — hər nümunənin ÖLÇÜLMÜŞ
    `depth`-i var (layer mərkəzindən `depth_offset_a_layer0` qədər sürüşən
    A quyusunun K=0 nöqtəsi istisna olmaqla)."""
    samples = []
    factor = {"A": 1.00, "B": 1.06, "C": 0.94}
    layer_centre_depth = {0: 2005.0, 1: 2015.0, 2: 2025.0}   # dz=10, top=2000
    permx_by_layer = {0: 500.0, 1: 100.0, 2: 10.0}
    for name, (x, y) in _WELL_XY.items():
        for k, permx in permx_by_layer.items():
            depth = layer_centre_depth[k]
            if name == "A" and k == 0:
                depth += depth_offset_a_layer0
            samples.append(WellSample(
                well=name, x=x, y=y, layer=k, depth=depth,
                values={"PERMX": permx * factor[name], "PORO": 0.2 * factor[name]}))
    return WellDataset(samples=samples, source="test")


def _build(dataset, nz, interpolator=None, allow_cross_layer_fallback=False):
    builder = WellBasedGeologicalModelBuilder(interpolator or OrdinaryKriging())
    return builder.build(dataset, _spec(nz), kv_over_kh=0.2,
                         allow_cross_layer_fallback=allow_cross_layer_fallback)


def test_explicit_depth_column_changes_result_within_its_own_layer():
    """A quyusunun K=0 nöqtəsinin `depth`-i lay mərkəzindən çox uzaqlaşanda
    (məs. lay 1-ə daha yaxın) — kiçik `range_v` ilə bu, K=0-ın nəticəsinə
    təsir etməlidir. `range_v` böyükdürsə (izotrop-dan uzaq, demək olar
    şaquli fərqi görməzdən gəlir) təsir kiçik olmalıdır."""
    baseline, _ = _build(_layered_dataset_with_explicit_depth(0.0), nz=3,
                         interpolator=OrdinaryKriging(range_=200.0, range_v=20.0))
    shifted, _ = _build(_layered_dataset_with_explicit_depth(9.0), nz=3,
                        interpolator=OrdinaryKriging(range_=200.0, range_v=20.0))
    permx_baseline = baseline.property_maps["PERMX"].as_grid(baseline.grid.shape)
    permx_shifted = shifted.property_maps["PERMX"].as_grid(shifted.grid.shape)
    assert not np.allclose(permx_baseline[0], permx_shifted[0]), (
        "A quyusunun ölçülmüş depth-i dəyişəndə (kiçik range_v ilə) K=0 "
        "nəticəsi dəyişmədi — 3D Z-awareness işləmir")


def test_cross_layer_fallback_weights_nearest_layer_more_than_distant_layer():
    """M2-nin əsas nəticəsi: K=3 (laysız) M1-in açıq-razılıq ehtiyatı ilə
    K=0,1,2-dən borc alanda, ən yaxın lay (K=2, dərinlik 2025) ən uzaq
    laydan (K=0, dərinlik 2005) DAHA ÇOX çəki almalıdır — bərabər orta
    DEYİL. Bunu iki ssenari ilə fərqləndiririk: yalnız K=0-ın PERMX-i
    dəyişəndə K=3 az dəyişməli, yalnız K=2-nin PERMX-i (eyni nisbətdə)
    dəyişəndə isə DAHA ÇOX dəyişməlidir."""
    def _dataset_with_layer_scaled(k, factor):
        samples = []
        base_factor = {"A": 1.00, "B": 1.06, "C": 0.94}
        layer_centre_depth = {0: 2005.0, 1: 2015.0, 2: 2025.0}
        permx_by_layer = {0: 500.0, 1: 100.0, 2: 10.0}
        for name, (x, y) in _WELL_XY.items():
            for kk, permx in permx_by_layer.items():
                value = permx * base_factor[name] * (factor if kk == k else 1.0)
                samples.append(WellSample(
                    well=name, x=x, y=y, layer=kk, depth=layer_centre_depth[kk],
                    values={"PERMX": value, "PORO": 0.2 * base_factor[name]}))
        return WellDataset(samples=samples, source="test")

    interpolator = OrdinaryKriging(range_=200.0, range_v=20.0, nugget=0.0)

    base, _ = _build(_dataset_with_layer_scaled(0, 1.0), nz=4,
                     interpolator=interpolator, allow_cross_layer_fallback=True)
    layer0_changed, _ = _build(_dataset_with_layer_scaled(0, 3.0), nz=4,
                               interpolator=interpolator, allow_cross_layer_fallback=True)
    layer2_changed, _ = _build(_dataset_with_layer_scaled(2, 3.0), nz=4,
                               interpolator=interpolator, allow_cross_layer_fallback=True)

    permx_base = base.property_maps["PERMX"].as_grid(base.grid.shape)[3].mean()
    permx_layer0_changed = layer0_changed.property_maps["PERMX"].as_grid(
        layer0_changed.grid.shape)[3].mean()
    permx_layer2_changed = layer2_changed.property_maps["PERMX"].as_grid(
        layer2_changed.grid.shape)[3].mean()

    shift_from_layer0 = abs(permx_layer0_changed - permx_base)
    shift_from_layer2 = abs(permx_layer2_changed - permx_base)

    assert shift_from_layer2 > shift_from_layer0, (
        f"K=3 (dərinlik ~2035) ən yaxın laydan (K=2, 2025) ən uzaq laydan "
        f"(K=0, 2005) DAHA ÇOX təsirlənməlidir; alınan sürüşmələr: "
        f"K=0 dəyişəndə {shift_from_layer0:.3g}, K=2 dəyişəndə {shift_from_layer2:.3g}")
