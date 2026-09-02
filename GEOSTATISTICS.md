# Variogram model fitting və tam 3D anizotropluq (M5)

## Niyə lazımdır

`OrdinaryKriging` əvvəllər YALNIZ sabit sferik variogram işlədirdi.
`range_`/`sill` verilməyəndə evristika idi:

    range = (nöqtələrin əhatə etdiyi sahənin) diaqonalı / 3
    sill  = var(dəyərlər)

Bu, "əks halda nə fərz edək" sualına neytral cavabdır, amma MƏLUMATIN
ÖZÜNDƏN heç nə öyrənmir — real korrelyasiya radiusu 3 dəfə fərqli ola
bilər, evristika bunu görməz. Anizotropluq da tək bir şaquli
nisbətlə (`range_v`) məhdud idi — üfüqi müstəvidə istiqamətlilik
(məs. çay kanalı, fay xətti boyu uzanan keçiricilik) tamamilə
nəzərə alınmırdı.

## Yeni iş axını (`imex2d/geology/variogram.py`)

    quyu nöqtələri
        → deneysel (empirical) variogram — lag-binlənmiş yarım-dəyişkənlik
        → model fit (sferik/eksponensial/qauss, çəkili ən kiçik kvadrat)
        → (istəyə görə) 6 istiqamətdə fit → ən böyük/kiçik radiuslu ox
        → OrdinaryKriging (auto_fit=True / auto_detect_anisotropy=True)

Bütün model funksiyaları "praktiki radius" konvensiyasında yazılıb:
γ(range_) ≈ 0.95·sill hər üç modeldə — ona görə `range_` modellər arasında
birbaşa müqayisə edilə bilir.

## Model seçimi — vizual DEYİL, doğrulama əsaslı

`select_best_variogram_model()` hər namizəd modeli (sferik/eksponensial/
qauss) fit edir, sonra HƏR BİRİ ilə əsl Kriging qurub `cross_validation.py`
(leave-one-out/k-fold) ilə REAL proqnoz RMSE-sini ölçür. Seçim ən kiçik
CV RMSE-yə görədir — deneysel variogrammla ən "hamar" uyğunlaşan model
YOX (bu, tapşırığın "vizual ən hamar model seçmə" qadağasına uyğundur).

## Tam 3D anizotropluq

`AnisotropyParams(azimuth_deg, range_major, range_minor, range_vertical)`
nöqtələri fırladıb miqyaslayaraq (X,Y,Z)-i vahid `range_major` radiuslu
izotrop fəzaya çevirir — standart geometrik anizotropluq transformu.
`azimuth_deg` — major oxun istiqaməti, +Y (Şimal) oxundan SAAT ƏQRƏBİ
istiqamətində (0°=Şimal, 90°=Şərq).

`azimuth_deg`/`range_minor` verilməyəndə (defolt) transform ƏVVƏLKİ
M2 davranışını ƏDƏD-ƏDƏD verir (yalnız Z miqyaslanır) — bax
`test_anisotropy_transform_reduces_to_z_only_scaling_when_isotropic_horizontal`.

`detect_anisotropy()` 6 üfüqi istiqamətdə sferik radius fit edib ən böyük
radiuslu istiqaməti major ox seçir. **Tipik reservoir modelində 5-20
quyu ilə bu ƏKSƏR HALLARDA etibarsız olacaq** — istiqamətli binlərdə
kifayət qədər cüt olmayacaq. Bu GİZLƏDİLMİR: `reliable=False` və izahlı
`warnings` qaytarılır, izotrop defolta qayıdılır. Sınanıb:
`test_detect_anisotropy_reports_unreliable_with_sparse_wells` (5 quyu).

Sintetik 180 nöqtəlik anizotrop sahədə (məlum azimut=30°, major=150,
minor=40) aşkarlanma: azimut ±35° daxilində, nisbət < 0.85 — bax
`test_detect_anisotropy_recovers_known_azimuth_and_ratio`.

## `OrdinaryKriging` üzərində yeni parametrlər

| Parametr | Defolt | Təsir |
|---|---|---|
| `model` | `"spherical"` | `"exponential"`, `"gaussian"`, ya da `"auto"` (yalnız `auto_fit=True` ilə) |
| `auto_fit` | `False` | `True` olanda `range_`/`sill` verilməyibsə deneysel variogramdan fit edilir |
| `auto_fit_nugget` | `False` | `True` olanda nugget də fit edilir (defoltda `nugget=0.0` saxlanılır) |
| `azimuth_deg` | `None` | verilsə üfüqi müstəvidə tam anizotrop transform işə düşür |
| `range_minor` | `None` | üfüqi minor radius (`range_` = major) |
| `auto_detect_anisotropy` | `False` | `True` və `azimuth_deg=None` olanda `detect_anisotropy()` çağırılır (etibarsızdırsa xəbərdarlıqla izotropa qayıdır) |

Hamısının defoltu KÖHNƏ davranışı verir — 76 mövcud kriging/geologiya
testi dəyişməz qaldı (`test_kriging_3d_anisotropy.py`,
`test_layer_aware_kriging_leak.py`, `test_cross_validation.py`, s.).

Fit/aşkarlanma nəticələri introspeksiya üçün saxlanılır:
`kriging.last_fit_`, `kriging.last_anisotropy_`, `kriging.last_warnings_`
(hesabat/UI-ya bağlamaq üçün hazırdır, hələ bağlanmayıb — bax aşağı).

## Sınanmadı / qalan iş

- **UI-da (`ui/panels.py`) yeni parametrlər hələ ekspoz edilməyib** —
  yalnız kod səviyyəsində istifadə oluna bilər. UI cəhdi bilərəkdən bu
  mərhələdə edilmədi (tapşırığın prioritet sırası: fiziki düzgünlük >
  ədədi stabillik > doğrulama > qeyri-müəyyənlik > performans > UI).
- Universal Kriging (trend-li), Simple Kriging, Matérn modeli —
  Phase 2-də tələb olunub, bu mərhələdə YAZILMADI.
- Şaquli istiqamətli variogram (`experimental_variogram(vertical=True)`)
  yazılıb və sınanıb, amma `OrdinaryKriging`-ə (auto-fit üçün) hələ
  bağlanmayıb — hazırkı `range_v` hələ də əl ilə verilir və ya `range_`-ə
  bərabər qəbul edilir.
- **Tam tenzor permeabilite (`PermeabilityTensor`, bax `ARCHITECTURE.md`
  §5.13) bu boru xəttindən HEÇ VAXT DOLDURULMUR** — co-kriging və ya
  digər çoxdəyişənli/tenzor-əsaslı interpolyasiya bu fazada YAZILMADI
  (bilərəkdən, MPFA-O-dan əvvəlki hazırlıq mərhələsidir). Kxy/Kxz/Kyz
  komponentləri hazırda YALNIZ kodda ƏLLƏ (birbaşa `PropertyMap`
  qurularaq) verilə bilər — heç bir geoloji məlumat İCAD EDİLMİR.
