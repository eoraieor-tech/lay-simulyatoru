# Testlər

## İşə salmaq

```bash
pip install pytest        # bir dəfə
pytest                    # bütün testlər
pytest -q                 # qısa
IMEX_SKIP_SLOW=1 pytest   # yavaş testləri keç
pytest tests/test_physics.py -v      # yalnız fizika
```

`pytest` yoxdursa, ehtiyat runner işləyir (əlavə paket tələb etmir):

```bash
python run_tests.py
IMEX_SKIP_SLOW=1 python run_tests.py -q
```

## Fayllar

| Fayl | Test | Nəyi qoruyur |
|---|---|---|
| `test_domain.py` | 12 | Grid indeksləşməsi, xassələr, SCAL son nöqtələri, model yoxlaması |
| `test_discretization.py` | 5 | Harmonik orta, Peaceman quyu indeksi, skin təsiri |
| `test_config_and_service.py` | 9 | Konfiqurasiya, DI davranışı, layihə obyekti |
| `test_physics.py` | 6 | **Material balansı, Bukley-Leverett, fiziki hədlər** |
| `test_regression.py` | 2 | **Etalon nəticə: RF 16.840 %, 4314 addım** |
| `test_rendering.py` | 8 | Rendering-in Qt-dən asılı olmaması, ox ölçüsünün sabitliyi |
| `test_pvt.py` | 17 | **A1: PVT cədvəli, korrelyasiyalar, provider** |
| `test_gas_pvt.py` | 22 | **A7/1: Z-faktoru, Bg, μg, geriyə uyğunluq** |
| `test_three_phase.py` | 17 | **A7/2: üç fazalı doyumluluq, qaz papağı equilibration** |
| `test_stone_relperm.py` | 22 | **A7/3: Stone II, iki fazalı reduksiya, mənfi kro yoxdur** |
| `test_variable_switching.py` | 20 | **A7/4: dəyişən keçid, kütlə balansı, kəsilməzlik** |
| `test_three_phase_residual.py` | 43 | **A7/5-6c(4): tam Jakobian yığımı, sistemli sonlu fərqlə doğrulama** |
| `test_standard_well.py` | 26 | **OPM quyu modeli, mərhələ 2: BHP-dən debitlər, idarəetmə tənlikləri, HAMARLIQ** |
| `test_well_state.py` | 16 | **OPM tipli quyu modeli, mərhələ 1: BHP naməlum dəyişən kimi** |
| `test_three_phase_newton.py` | 18 | **A7/6d: üç fazalı Nyuton döngəsi + çökməyə qarşı mühafizə (heç vaxt istisna atmır)** |
| `test_gas_ui_wiring.py` | 8 | **A7/6d(UI): mühərrik seçimi, uyğun xətti həlledici, təhlükəsiz uğursuzluq** |
| `test_diagnostics.py` | 21 | +1: köhnə "qaz modelləşdirilmir" xəbərdarlığı qaz aktivdirsə göstərilmir |
| `test_vtk_volume.py` | 47 | **VTK 3D motoru — pytest, `vtk` tələb edir** (həndəsə, filtrlər, rəng xəritəsi, tam offscreen render) |
| `test_opm_import.py` | 12 | **OPM Flow idxalı — pytest, `resdata` tələb edir** (sintetik round-trip Eclipse halı, öz VolumeRenderer-imizlə çəkilmə) |
| `test_initialization.py` | 14 | **A3: equilibration, hidrostatik qradiyent, OWC** |
| `test_capillary_gravity.py` | 16 | **A4: Brooks-Corey Pc, keçid zonası, cazibə potensialı** |
| `test_three_dimensional.py` | 16 | **A5: 3D grid, şaquli axın, Kv/Kh, kəsik vizuallaşdırma** |
| `test_perforation.py` | 6 | Perforasiya intervalı, qismən açılmış quyular |
| `test_serialization.py` | 13 | **B1: layihə faylı (.imx), eyni nəticənin bərpası** |
| `test_geology_import.py` | 23 | **B2: quyu CSV, IDW/Kriging, geoloji model qurulması** |
| `test_diagnostics.py` | 20 | Xəta/xəbərdarlıq ayrımı, quyu rejimi, loglama, matris keşi |
| `test_ui_static.py` | 9 | UI qatının AST yoxlanışı (import, atribut, qurucu, tab indeksi, tab siyahısı) |
| `test_ui_wiring.py` | 7 | UI panel–model bağlantısı |
| `test_implicit_residual.py` | 16 | **A6/1: qalıq vektoru, kütlə balansı, IMPES ilə ardıcıllıq** |
| `test_implicit_jacobian.py` | 20 | **A6/2: analitik Jakobian vs sonlu fərq (8 konfiqurasiya)** |
| `test_implicit_newton.py` | 16 | **A6/3: Nyuton döngəsi, konvergensiya, IMPES ilə uyğunluq** |
| `test_implicit_engine.py` | 17 | **A6/4: adaptiv Δt, FullyImplicitEngine** |
| `test_cpr.py` | 16 | **A6/5: CPR dekuplinqi, blok-Jakobi, yaddaş üstünlüyü** |
| `test_volume_rendering.py` | 32 | **3D görüntü: üz çıxarışı, filtr, kəsim, işıqlandırma, baxış** |

Cəmi **722 test**.py`).

Performans ölçmələri: `PERFORMANCE.md` və `tools/benchmark.py`.
Fully implicit sxemin planı və nəticələri: `A6_PLAN.md`.

## Etalon dəyərlər

`helpers.py` içindəki `REFERENCE_FIVE_SPOT` refaktorinqdən ƏVVƏLKİ
`core.py`-dən götürülüb:

| Göstərici | Dəyər | Tolerans |
|---|---|---|
| Recovery Factor | 16.840 % | ±0.01 |
| Addım sayı | 4314 | dəqiq |
| OOIP | 1 029 064.3 m³ | ±1 m³ |

Bu rəqəmlərdən biri dəyişirsə, hesablama davranışı dəyişib.
Dəyişiklik qəsdəndirsə, etalon yenilənir və səbəbi qeyd edilir.

## Diaqnostika: xəta vs xəbərdarlıq

`model.diagnose()` iki səviyyə qaytarır:

| Səviyyə | Davranış |
|---|---|
| `ERROR` | Model işə salınmır (grid-dən kənar perforasiya, quyu yoxdur) |
| `WARNING` | İşə salınır, amma istifadəçidən təsdiq istənilir |

`model.validate()` geriyə uyğunluq üçün saxlanılıb — yalnız `ERROR`
mesajlarını qaytarır.

## UI qatı niyə statik yoxlanılır

UI kodu Qt tələb etdiyi üçün avtomatik testlərlə əhatə olunmurdu və
iki dəfə eyni tipli səhv buraxıldı: yaradılmadan çağırılan atribut
(`pvt_axes`) və import edilməyən ad (`QSpinBox`). `test_ui_static.py`
AST təhlili ilə bu sinif səhvləri Qt işə salmadan tutur.

Yoxlanıldı: `QSpinBox` importu qəsdən silindikdə test dərhal xəta verir.

3D tab əlavə olunanda iki yeni sinif səhv üzə çıxdı və hər ikisi üçün
test yazıldı:

| Səhv | Test |
|---|---|
| `_build_volume_tab()` yazılıb, amma çağırılmayıb | `test_builder_methods_are_actually_called` |
| `__init__`-də metod atributundan əvvəl çağırılıb | `test_init_does_not_call_methods_before_their_attributes_exist` |
| Sərt tab indeksləri yeni tab əlavə olunanda sürüşür | `test_tabs_are_selected_by_name_not_index` |

Birincisi xüsusilə maraqlıdır: atribut testi onu TUTMADI, çünki
`self.volume_time` metod daxilində təyin olunmuşdu — sadəcə həmin
metod heç vaxt icra olunmurdu.

## Mutasiya yoxlaması

Testlərin həqiqətən xəta tutduğu qəsdən pozuntularla yoxlanılıb:

| Pozuntu | Tutuldu |
|---|---|
| Harmonik orta → arifmetik orta | `test_transmissibility_uses_harmonic_mean` |
| Peaceman sabiti 0.28 → 0.20 | `test_peaceman_well_index_matches_hand_calculation` |
| Upstream → downstream çəkilənmə | `test_water_material_balance_is_conserved` (xəta 99 %) |
| `QSpinBox` importunun silinməsi | `test_every_qt_name_used_in_ui_is_imported` |

## Modul əlavə edəndə

Qayda: mövcud reqressiya testləri **dəyişmədən keçməlidir**.
Keçmirsə, yeni funksiya köhnə yolu sındırıb.

Hər yeni provider (PVT, kapilyar, initialization) verilmədikdə
davranış əvvəlki kimi qalmalıdır — bu, ayrıca testlə qorunur.
