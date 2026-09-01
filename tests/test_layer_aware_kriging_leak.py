"""M0/M1 — quyu-lay sızmasının sınan testi VƏ düzəlişin reqressiya testi.

M0 mərhələsində bu fayl YALNIZ nöqsanı sənədləşdirirdi (kod dəyişməmişdi).
M1-də `geology_service.py`-a düzəliş edildikdən sonra iki test yeniləndi
(`test_layer0_change_must_not_leak_into_unsampled_layer` və
`test_unsampled_layer_result_equals_pooled_fallback_mechanism` əvəzinə
aşağıdakı iki test) ki, DÜZƏLİŞİ də sübut etsin, yalnız nöqsanı yox. M2
(3D/anizotrop Kriging) `tests/test_kriging_3d_anisotropy.py`-dadır.

FƏRZ EDİLƏN NÖQSAN (istifadəçi hesabatı): quyuda lay 1-3 üçün qeyd edilmiş
PORO/PERMX dəyərləri, məlumatı olmayan lay 5-ə "sızır".

KÖK SƏBƏB (Mərhələ 0 analizi ilə tapılıb, bax söhbətdəki hesabat):

  imex2d/application/geology_service.py:151-152 (`_interpolate_volume`)
      if values.size == 0:                      # bu təbəqədə məlumat yoxdur
          points, values = dataset.points(source, None)   # <-- BURADA
  imex2d/domain/well_data.py:62 (`WellDataset.points`)
      if layer is not None and sample.layer is not None and sample.layer != layer:
          continue

  `layer=None` sorğusu iki fərqli məna daşıyır və bu modul onları
  qarışdırır: (1) "bu nümunə heç bir laya aid deyil, hər yerdə keçərlidir"
  (sample.layer is None-un həqiqi mənası) və (2) "bu K-də məlumat yoxdur,
  ən azı bir şey tap" (boş-lay ehtiyat mexanizmi, geology_service.py:152).
  (2) münasibətində filter TAM SÖNDÜRÜLÜR, ona görə açıq lay etiketli
  (layer=0,1,2) nümunələr, aid olmadıqları laya (məs. K=4) sızır.

  İKİNCİ, MÜSTƏQİL NÖQSAN (A tipi): `imex2d/application/geology_adapter.py`
  (UI-dən redaktə olunan quyu cədvəli — CSV-nin əvəzi, bax commit e70fa49)
  `WellSample`-ı HEÇ VAXT `layer`-lə doldurmur (`domain/geology.py`-də
  `GeologicalWell`-in per-lay porosity/permeability sahəsi ümumiyyətlə
  yoxdur — bir quyuda cəmi BİR φ, BİR k dəyəri var). Ona görə
  `dataset.is_layered()` bu yolda həmişə `False`-dur və
  `geology_service.py:149` hər zaman `layer=None` göndərir — tək dəyər
  bütün K-lara SİNKRON YAYILIR (bug deyil, sxemin özündə lay anlayışı
  yoxdur). Bu, aşağıda `test_ui_table_wells_broadcast_single_value_to_every_layer`
  ilə sənədləşdirilib.

"""

from __future__ import annotations

import numpy as np

from helpers import default_scal, make_service, short_config
from imex2d.application.config import (OutputConfig, SimulationConfig,
                                       TimeSteppingConfig)
from imex2d.application.geology_adapter import wells_to_dataset
from imex2d.application.geology_service import (GeologicalGridSpec,
                                                 WellBasedGeologicalModelBuilder)
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.domain.geology import GeologicalWell
from imex2d.domain.geometry import xy_to_ij
from imex2d.domain.initial import InitialConditions
from imex2d.domain.structure import FaultReference
from imex2d.domain.well_data import WellDataset, WellSample
from imex2d.domain.wells import ControlMode, Perforation, Well, WellControl, WellType
from imex2d.geology.interpolation import OrdinaryKriging
from imex2d.simulation.discretization import TwoPointFluxDiscretization

# 3 quyu, X/Y-də fərqli yerlərdə — 50x3 = 150 m tutan 3x3 arealda
_WELL_XY = {"A": (25.0, 25.0), "B": (125.0, 25.0), "C": (25.0, 125.0)}

# lay üzrə "əsl" dəyərlər (tapşırıqdakı ədədlərlə üst-üstə düşür)
_PERMX_BY_LAYER = {0: 500.0, 1: 100.0, 2: 10.0}   # lay 1, 2, 3 (1-based)
_PORO_BY_LAYER = {0: 0.25, 1: 0.20, 2: 0.15}


def _layered_dataset(permx_layer0_a: float = 500.0) -> WellDataset:
    """3 quyu x 3 lay (K=0,1,2) — K=3,4 üçün HEÇ BİR YERDƏ məlumat yoxdur.

    `permx_layer0_a` — yalnız A quyusunun 0-cı lay PERMX dəyəri dəyişdirilə
    bilsin deyə parametrləşdirilib (sızma üçün diferensial test).
    """
    samples = []
    for name, (x, y) in _WELL_XY.items():
        # kiçik quyu-üzrə vurma əmsalı (əlavə YOX) ki, kriging səthi sabit
        # olmasın, amma lay 3-ün (10 mD) balası mənfi olmasın
        factor = {"A": 1.00, "B": 1.06, "C": 0.94}[name]
        for k, permx in _PERMX_BY_LAYER.items():
            value = permx_layer0_a if (name == "A" and k == 0) else permx * factor
            samples.append(WellSample(
                well=name, x=x, y=y, layer=k,
                values={"PERMX": value, "PORO": _PORO_BY_LAYER[k] * factor}))
    return WellDataset(samples=samples, source="test")


def _build(dataset: WellDataset, nz: int = 5, allow_cross_layer_fallback: bool = False):
    spec = GeologicalGridSpec(nx=3, ny=3, nz=nz, dx=50.0, dy=50.0, dz=10.0,
                              top_depth=2000.0)
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    model, report = builder.build(dataset, spec, kv_over_kh=0.2,
                                  allow_cross_layer_fallback=allow_cross_layer_fallback)
    return model, report


# ── M1 DÜZƏLİŞİNDƏN SONRA: boş lay artıq sükutla sızdırmır ─────────────
# (əvvəlki `test_layer0_change_must_not_leak_into_unsampled_layer` və
# `test_unsampled_layer_result_equals_pooled_fallback_mechanism` bu iki
# testlə əvəzləndi — bax ISH_HESABATI.md, M1)
def test_unsampled_layer_raises_instead_of_silently_pooling_other_layers():
    """K=3 və K=4-də (heç bir quyuda məlumat yoxdur) `build()` artıq
    sükutla başqa layların (K=0,1,2) nöqtələrini hovuzlamır — açıq
    `ValueError` atır. Bu, M1 düzəlişinin əsas reqressiya qoruyucusudur:
    əvvəllər bu ssenari sükutla keçib lay 1-in dəyərini lay 5-ə
    sızdırırdı (bax M0: `git log`-da bu commitdən əvvəlki versiya)."""
    dataset = _layered_dataset()
    try:
        _build(dataset, nz=5)   # K=3, K=4 üçün heç bir quyuda məlumat yoxdur
    except ValueError as error:
        assert "heç bir quyu nöqtəsi yoxdur" in str(error)
        return
    raise AssertionError("Boş lay (K=3/4) üçün ValueError gözlənilirdi, tapılmadı")


def test_layer0_change_no_longer_affects_unsampled_layer_when_fallback_is_explicit():
    """`allow_cross_layer_fallback=True` ilə (bilərəkdən açılan köhnə
    davranış) K=4 nəticəsi HƏLƏ DƏ K=0-dan asılıdır — bu gözlənilir,
    çünki istifadəçi ekstrapolyasiyaya açıq razılıq verib. Fərq: bu artıq
    DEFOLT deyil və `report.warnings`-də açıq qeyd olunur."""
    baseline, report_a = _build(_layered_dataset(permx_layer0_a=500.0), nz=5,
                                allow_cross_layer_fallback=True)
    changed, report_b = _build(_layered_dataset(permx_layer0_a=5000.0), nz=5,
                               allow_cross_layer_fallback=True)

    permx_baseline = baseline.property_maps["PERMX"].as_grid(baseline.grid.shape)
    permx_changed = changed.property_maps["PERMX"].as_grid(changed.grid.shape)
    assert not np.allclose(permx_baseline[4], permx_changed[4]), (
        "allow_cross_layer_fallback=True ilə K=4 hələ də K=0-dan asılı "
        "olmalıdır (bu, bilərəkdən seçilmiş ekstrapolyasiyadır)")
    assert any("PERMX" in w and "K=3" in w for w in report_a.warnings), (
        "cross-layer fallback işə düşəndə report.warnings-də açıq qeyd olmalıdır")


# ── müsbət nəzarət: hər layın öz məlumatı olanda sızma YOXDUR ──────────
def test_each_fully_sampled_layer_keeps_its_own_value_when_no_fallback_triggers():
    """Hər 3 lay üçün də məlumat olanda (K=0,1,2, nz=3) per-K dövrü heç vaxt
    boş-lay ehtiyatına düşmür — kriging riyaziyyatı düzgündür, nöqsan yalnız
    ehtiyat budağındadır. Lay 1-in dəyəri lay 2/3-ə kopyalanmır."""
    dataset = _layered_dataset()
    model, _report = _build(dataset, nz=3)     # yalnız K=0,1,2 - hamısında məlumat var
    permx = model.property_maps["PERMX"].as_grid(model.grid.shape)

    layer0_mean = float(permx[0].mean())
    layer1_mean = float(permx[1].mean())
    layer2_mean = float(permx[2].mean())

    assert layer0_mean > layer1_mean > layer2_mean, (
        f"laylar üzrə orta PERMX ({layer0_mean:.1f}, {layer1_mean:.1f}, "
        f"{layer2_mean:.1f}) gözlənilən 500>100>10 sırasını saxlamır")
    # lay 1 lay 2-yə ədəd-ədəd bərabər olmamalıdır (kopyalanma yoxdur)
    assert not np.allclose(permx[0], permx[1])
    assert not np.allclose(permx[1], permx[2])
    # amma hər laydakı dəyər öz təbəqəsinin əsl dəyərlərinə yaxın olmalıdır
    assert abs(layer0_mean - 500.0) / 500.0 < 0.1
    assert abs(layer2_mean - 10.0) / 10.0 < 0.1


# ── A-tipi: UI cədvəli hər quyu üçün cəmi BİR dəyər daşıyır ─────────────
def test_ui_table_wells_broadcast_single_value_to_every_layer():
    """`GeologicalWell` (proqramda indi işlədilən quyu cədvəli, bax
    domain/geology.py) per-lay φ/k sahəsi daşımır — bir quyu üçün cəmi BİR
    porosity, BİR permeability var. `geology_adapter.wells_to_dataset` bunu
    heç vaxt `WellSample.layer`-ə yazmır, ona görə `is_layered()` HƏMİŞƏ
    False-dur və `geology_service.py:149` hər K üçün eyni 2D nəticəni
    işlədir — tək dəyər bütün laylara SİNKRON yayılır.

    Bu, "sızma" deyil, İKİNCİ, MÜSTƏQİL bir nöqsan/çatışmazlıqdır: UI
    cədvəlində ümumiyyətlə per-lay giriş yolu yoxdur. Bu test bugünkü
    kodla KEÇİR — mexanizmi sənədləşdirir."""
    wells = [
        GeologicalWell(name="A", x=25.0, y=25.0, porosity=0.25, permeability=500.0),
        GeologicalWell(name="B", x=125.0, y=25.0, porosity=0.22, permeability=480.0),
        GeologicalWell(name="C", x=25.0, y=125.0, porosity=0.20, permeability=520.0),
    ]
    dataset, skipped = wells_to_dataset(wells, method="Kriging (adi)")
    assert skipped == {}
    assert not dataset.is_layered()

    model = _build_from_dataset_no_layer(dataset, nz=3)
    permx = model.property_maps["PERMX"].as_grid(model.grid.shape)
    poro = model.property_maps["PORO"].as_grid(model.grid.shape)

    assert np.allclose(permx[0], permx[1]) and np.allclose(permx[1], permx[2]), (
        "UI cədvəlindən qurulan modeldə PERMX laylar üzrə fərqlənir — "
        "gözlənilməz, çünki GeologicalWell-də per-lay sahə yoxdur")
    assert np.allclose(poro[0], poro[1]) and np.allclose(poro[1], poro[2])


def _build_from_dataset_no_layer(dataset: WellDataset, nz: int):
    spec = GeologicalGridSpec(nx=3, ny=3, nz=nz, dx=50.0, dy=50.0, dz=10.0,
                              top_depth=2000.0)
    builder = WellBasedGeologicalModelBuilder(OrdinaryKriging())
    model, _ = builder.build(dataset, spec, kv_over_kh=0.2)
    return model


# ── inteqrasiya: kriging-dən transmissivliyə/PV-yə/kütlə balansına ─────
def _reservoir_model(nz: int = 3, faults=None):
    dataset = _layered_dataset()
    geology, _report = _build(dataset, nz=nz)
    wells = [
        Well("INJ", WellType.INJECTOR, WellControl(ControlMode.RATE, 20.0),
             [Perforation(0, 0, k) for k in range(nz)]),
        Well("PROD", WellType.PRODUCER, WellControl(ControlMode.BHP, 190.0),
             [Perforation(2, 2, k) for k in range(nz)]),
    ]
    scal = default_scal()
    return ReservoirModelBuilder().build(
        geology, wells, scal=scal, fault_references=faults,
        initial=InitialConditions(datum_pressure=200.0, water_saturation=scal.swc),
        name="M0 kriging pipeline testi")


def test_well_xy_maps_to_expected_grid_cell():
    """Quyu koordinatı düzgün (i, j) hüceyrəsinə uyğunlaşır (dx=dy=50 m,
    3x3 grid, quyu A (25,25) -> (0,0) hüceyrəsi, quyu B (125,25) -> (2,0))."""
    model = _reservoir_model(nz=3)
    i, j = xy_to_ij(25.0, 25.0, model.geometry)
    assert (i, j) == (0, 0)
    i, j = xy_to_ij(125.0, 25.0, model.geometry)
    assert (i, j) == (2, 0)


def test_pore_volume_uses_per_cell_kriged_porosity():
    """PV = phi * V_bulk — laylar üzrə fərqli porosity fərqli PV verir."""
    model = _reservoir_model(nz=3)
    pv = model.pore_volume().reshape(model.grid.shape)
    bulk = model.geometry.volumes().reshape(model.grid.shape)
    poro = model.rock.porosity.values.reshape(model.grid.shape)
    assert np.allclose(pv, poro * bulk)
    # lay 1 (yüksək porosity) PV-si lay 3-dən (aşağı porosity) böyükdür
    assert pv[0].sum() > pv[2].sum()


def test_transmissibility_uses_permx_permy_horizontally_and_permz_vertically():
    model = _reservoir_model(nz=3)
    grid = TwoPointFluxDiscretization().build(model)
    conn = grid.connections
    assert np.all(grid.transmissibility[conn.axis < 2] > 0)
    assert np.all(grid.transmissibility[conn.axis == 2] > 0)
    # PERMZ kv_over_kh=0.2 ilə PERMX-dən kiçikdir -> şaquli T horizontaldan
    # sistematik fərqli olmalıdır (eyni geometriyada eyni olsaydı kv/kh
    # tətbiq olunmurdu deməkdir)
    assert not np.allclose(grid.transmissibility[conn.axis == 2].mean(),
                           grid.transmissibility[conn.axis == 0].mean(), rtol=0.2)


def test_sealing_fault_zeroes_transmissibility_even_with_high_permeability():
    """Lay 1-in (K=0) PERMX-i 500 mD-dir (yüksək), amma sealing fay
    olduğu üzdə T=0 olmalıdır — keçiricilik yüksək olsa belə."""
    model = _reservoir_model(nz=3, faults=[FaultReference(
        name="F1", source_id="F1", axis="I", plane_index=1, sealing=True)])
    grid = TwoPointFluxDiscretization().build(model)
    i_a, j_a, k_a = model.grid.ijk_array(grid.connections.cell_a)
    on_fault = (grid.connections.axis == 0) & (i_a == 1) & (k_a == 0)
    assert on_fault.any()
    assert np.allclose(grid.transmissibility[on_fault], 0.0)


def test_mass_balance_holds_for_kriging_derived_reservoir_model():
    """Well data -> Kriging -> ReservoirModel boru xəttinin ucdan-uca kütlə
    balansı. Bu, mövcud `test_material_balance_holds_in_three_dimensions`
    (SyntheticGeologicalModelBuilder ilə) ilə eyni yoxlamadır, sadəcə
    modeli kriging boru xəttindən qurur."""
    scal = default_scal()
    model = _reservoir_model(nz=3)
    # kiçik 3x3x3 grid — default max_dt=20 gün bu ölçüdə həddindən artıq
    # iri addımdır (CFL-i pozur), ona görə 1D BL testindəki kimi kiçik
    # addım seçilir (bax tests/helpers.py bl_config). end_time qısa
    # olanda (150) producer hələ heç nə vermir (su cəbhəsi çatmayıb) —
    # bu halda trapezoid inteqrasiyasının son-addım artefaktı (~1 günlük
    # sərf) "istehsal ~0"-a nisbətdə süni şəkildə böyük görünür; bu, kriging
    # sızması ilə əlaqəsizdir, ona görə axın həqiqətən başlayana qədər
    # işlədilir.
    config = SimulationConfig(
        end_time=800.0,
        time_stepping=TimeSteppingConfig(max_dt=2.0, cfl_factor=0.4),
        output=OutputConfig(snapshot_count=20))
    engine = make_service(scal).create_engine(model, config)
    initial_water = float(np.sum(model.pore_volume() * engine.sw))
    result = engine.run()
    final_water = float(np.sum(model.pore_volume() * engine.sw))
    series = result.series
    injected = float(np.trapezoid(series.water_injection_rate, series.time))
    produced = float(np.trapezoid(series.water_rate, series.time)) * model.fluids.water_fvf
    error = abs((final_water - initial_water) - (injected - produced)) / max(injected, 1e-9)
    assert error < 0.005, f"Kütlə balansı xətası {error * 100:.3f} %"
