# Sequential Gaussian Simulation (SGS) — kəsilməz xassələr (Phase 5)

## Audit: mövcud kəsilməz xassə iş axını

`geology_service.py::_interpolate_volume` → `interpolate_property()` →
`OrdinaryKriging`/`InverseDistance`. Bu, HƏMİŞƏ TƏK, DETERMİNİSTİK
(hamarlanmış) qiymət qaytarır — eyni giriş həmişə eyni nəticəni verir,
heç bir təsadüfilik/qeyri-müəyyənlik REALİZASİYASI yoxdur. `PropertyRule`
(log-transform + min/max klip) VƏ `OrdinaryKriging` (variogram,
anizotropluq, yerli axtarış — Phase 2-3) DƏYİŞMƏDİ, TƏKRARLANMADI —
SGS bunları BİRBAŞA istifadə edir (aşağı bax).

## Fayllar

**Yeni:** `imex2d/geology/gaussian_transform.py` (normal-score çevirmə),
`imex2d/geology/distribution_analysis.py` (paylanma statistikası, log-
fəza qərarı), `imex2d/geology/sgs.py` (əsas SGS mühərriki), + 6 test faylı.
**Dəyişdirilmiş:** `imex2d/geology/interpolation.py` (`interpolate_
with_variance()` YENİ, `interpolate()` DƏYİŞMƏDİ), `imex2d/geology/
hard_data.py` (kəsilməz ziddiyyət rejimi — `tolerance` parametri, YENİ,
kateqorik rejim DƏYİŞMƏDİ), `imex2d/geology/facies.py` (sill sabitliyi
düzəlişi, bax aşağı), `imex2d/application/geology_service.py`
(`ContinuousSGSConfig`, opt-in SGS marşrutlaşdırması).

## SGS alqoritmi

```
sərt data → (istəyə görə log) → normal-score çevirmə (Gauss fəzası)
    → hər hədəf üçün: yerli kriging (mean, variance) Gauss fəzasında
    → N(mean, variance)-dən NÜMUNƏ (Kriging ESTİMATE-i DEYİL!)
    → nəticə kondisioner çoxluğuna əlavə olunur
    → (bütün hədəflər bitəndən sonra) tərs çevirmə + (log-fəza isə) exp()
    → fiziki hədd (bounds) — KLİPLƏNMƏ SAYILIR, sükutla deyil
```

Sınanıb: eyni hədəf üçün fərqli seed-lər FƏRQLİ dəyər verir (kriging
ESTİMATE-i sabit olardı) — `test_sgs_result_is_not_the_smooth_kriging_
estimate`.

## Kriging variansı — YENİ, `interpolate()`-i DƏYİŞMƏYƏN əlavə

`OrdinaryKriging.interpolate_with_variance()` — σ²(x0) = Σw_i·γ(x_i,x0) + μ
(Laqranj vuruğu) ARTIQ HƏLL EDİLMİŞ xətti sistemdən oxunur, YENİDƏN HƏLL
YOXDUR. `interpolate()`-in ÖZÜ TOXUNULMAYIB (76+ Phase 2-3 testi
qorunur, YENİDƏN işlədilib təsdiqləndi).

## Normal-score çevirmə

Sıra-əsaslı (rank-based, Deutsch & Journel): `scipy.stats.rankdata
(method="average")` (TIES üçün), Hazen kvantili `(rank-0.5)/n`,
`scipy.stats.norm.ppf`. Tərs/irəli — DƏQİQ CƏDVƏLİN xətti
interpolyasiyası, ekstrapolyasiya YOXDUR (sərhədə kəsilir). Sabit
xassə (bütün dəyərlər eyni) AYRICA yoxlanılır — SGS keçilir, sabit
dəyər UYDURULMUŞ dəyişkənlik OLMADAN yayılır.

## Paylanma təhlili — log-fəza DATA-ƏSASLI qərarlaşdırılır

`log_transform_is_justified(values)`: `log(dəyər)`-in çarpıqlığı
ORİJİNAL dəyərin çarpıqlığından KİÇİKSƏ (mütləq qiymətcə) log fəzası
üstünlük təşkil edir — SABİT ehtimal (yalnız PERMX-ə görə) İLƏ DEYİL.

## Fasiya-şərtli SGS

`simulate_sgs_facies_conditioned()`: hər fasiya ÖZ sərt datası VƏ öz
hədəf hüceyrələri üzərində AYRI SGS. Kifayət qədər sərt data (<8,
`DEFAULT_MIN_HARD_DATA_FOR_OWN_MODEL`) OLMAYAN fasiya üçün AYRICA
model UYDURULMUR — BÜTÜN fasiyalar-arası data ilə FALLBACK, AÇIQ
xəbərdarlıqla (`used_cross_facies_fallback` metadata sahəsi).

## PORO/PERMX/PERMY/PERMZ

PORO: `bounds=(rule.minimum, rule.maximum)` defolt olaraq MÖVCUD
`DEFAULT_RULES["PORO"]` (0.01-0.45) hədlərindən GÖTÜRÜLÜR (TƏKRARLANMIR).
PERMX/Y/Z: `log_space=True/False/None` (None=avtomatik). Log-fəzada
`exp()` tərs çevirmədən sonra tətbiq olunur, nəticə HƏMİŞƏ müsbətdir
(sınanıb, `test_log_space_permeability_stays_strictly_positive`).
Geoloji modelin ÖZÜ verdiyi anizotropluq nisbətləri (`ky_over_kx`/
`kv_over_kh`, `_fill_missing_permeability`) SGS TƏRƏFİNDƏN YENİDƏN
İCAD EDİLMİR — SGS YALNIZ AÇIQ konfiqurasiya edilən sütunlar üçün işə
düşür, `_fill_missing_permeability` DƏYİŞMƏDƏN qalır.

## Sərt data + ziddiyyət

`hard_data.py` GENİŞLƏNDİRİLDİ (kateqorik rejim DƏYİŞMƏDƏN): `tolerance`
parametri ilə KƏSİLMƏZ rejim (`max-min > tolerance` = ziddiyyət).
`on_conflict="average"` YENİ (yalnız kəsilməz) — ziddiyyətli dəyərlərin
ORTASI sintetik nümunə kimi işlədilir. Defolt `conflict_tolerance=0.0`
(ƏN SƏRT) — istifadəçi AÇIQ boşluq vermədikcə HƏR fərq ziddiyyət sayılır.

**Sərt data → hüceyrə mərkəzinə sancma**: Phase 4.1-də tapılan EYNİ
prinsip (real quyu nadir hallarda dəqiq hüceyrə mərkəzindədir) burada
da tətbiq edilib — `_simulate_continuous_sgs_field` sərt nöqtəni öz
ev hüceyrəsinin mərkəzinə köçürür ki, hədəf massivi ilə DƏQİQ üst-üstə
düşsün.

## TAPILAN VƏ DÜZƏLDİLƏN HƏQİQİ SƏHVLƏR (bu fazada)

1. **Qeyri-sabit `sill`**: `_resolve_facies_variogram`/`_resolve_
   property_variogram`-ın FALLBACK budağı (variogram fit alınmayanda)
   `sill=None` saxlayırdı — `OrdinaryKriging._parameters()` `sill=None`
   olanda onu HƏR ÇAĞIRIŞDA ötürülən nöqtələrin öz variansından YENİDƏN
   hesablayır, bu, YERLİ axtarışın (qonşuluq alt-çoxluğu) ÖLÇÜSÜNDƏN
   ASILI OLARAQ DƏYİŞƏN, qeyri-stabil sill deməkdir — stasionar variogram
   fərziyyəsini POZUR. Fast/brute-force müqayisə testi (`tests/test_sgs_
   validation.py`) bunu TUTDU. Düzəliş: `sill` HƏMİŞƏ tam datadan BİR
   DƏFƏ hesablanan konkret ədədə HƏLL EDİLİR (`facies.py`-də DƏ eyni
   düzəliş edildi, kateqorik SIS-ə də aid idi).
2. **Yanlış n=1 kriging variansı**: `interpolate_with_variance()`-in
   ilk versiyası tək-nöqtə halında `variance=γ(məsafə)` (sadə-kriging-
   bənzər) hesablayırdı — DÜZGÜN ordinar kriging düsturu isə (Laqranj
   məhdudiyyətinin "naməlum ortalama" cəzası ilə) `variance=2γ(məsafə)`
   verir (əl hesabı ilə TƏSDİQLƏNİB, bax `_solve_global`-ın ÖZÜ). Bu,
   yerli axtarışda DƏQİQ 1 qonşu tapılanda sürətli/brute-force yolların
   FƏRQLİ nəticə verməsinə səbəb olurdu. Düzəliş: `interpolate_with_
   variance()` n=1 halında da BİRBAŞA `_solve_global`-ı çağırır (xüsusi
   düstur YOXDUR).
3. **DÜRÜST QALAN NÜANS (səhv DEYİL)**: müntəzəm (grid) hədəflərdə
   `max_neighbors` sərhədində dəqiq/yaxın-bərabər məsafəli namizədlər
   TEZ-TEZ olur — cKDTree və brute-force bu cür bərabərliyi FƏRQLİ (hər
   ikisi DÜZGÜN) namizədlə həll edə bilər. QEYRİ-grid (təsadüfi) hədəflərlə
   TAM (bit-bə-bit) uyğunluq sübut edilib (`test_fast_and_brute_force_
   search_agree_on_identical_realization`); grid hədəflərdə YALNIZ KİÇİK,
   sənədləşdirilmiş fərq (bax `test_fast_and_brute_force_differ_only_
   at_grid_tie_boundaries_not_wrong`).

## Performans (FAKTİKİ ölçmə, UYDURULMAMIŞ)

| Grid | Hüceyrə | sürətli | brute-force |
|---|---|---|---|
| 15×15 | 225 | 0.062 s | 0.044 s |
| 50×50 | 2500 | 0.677 s | 0.528 s |
| 100×100 | 10000 | 2.391 s | 3.538 s (1.48× sürətli) |
| 20×20×5 (3D) | 2000 | 0.442 s | (ölçülmədi) |

Facies/SIS-lə EYNİ tapıntı: kiçik/orta şəbəkədə brute-force rəqabətlidir,
sürət üstünlüyü ~10⁴ hüceyrədən başlayır.

## Elmi çəkincə

Bu modul PROQRAM CƏHƏTDƏN düzgün SGS-i tətbiq edir (Gauss çevirməsi
doğrudur, kriging variansı DÜZGÜN düsturla hesablanır — YUXARIDA
tapılan səhv daxil olmaqla YOXLANILIB, nümunə ŞƏRTİ paylanmadan çəkilir,
sərt data hörmət olunur). Bu, YERİN gerçək xassə paylanmasını əks
etdirdiyi demək DEYİL. Ensembl statistikası (P10/P50/P90, mean/std)
HƏDƏF paylanmaya YAXINLAŞIR (dəqiq bərabərlik gözlənilmir, sonlu
nümunə + kriging qeyri-müəyyənliyi). Deneysel variogram reproduksiyası
YALNIZ BÖYÜKLÜK dərəcəsində yoxlanılıb (`test_simulated_field_
reproduces_approximate_variogram_range`) — DƏQİQ ədədi uyğunluq İDDİA
EDİLMİR.

## Qalan iş / NOT DONE / NOT VALIDATED YET

* `simulate_sgs_facies_conditioned`-in `geology_service.py` inteqrasiyası
  YALNIZ `model.facies_fields`-dən oxuyur — çoxlu FACIES REALİZASİYASI
  üzərində SGS ansamblının necə birləşəcəyi (hər fasiya realizasiyası
  üçün ayrı PORO ansamblı) HƏLƏ AVTOMATLAŞDIRILMAYIB.
* UI-da SGS parametrləri (Phase 4.1-in `FaciesPanel`-inə bənzər) əlavə
  EDİLMƏDİ bu fazada — **NOT DONE**.
* Ayrıca "cross-facies müqayisə hesabatı" (mean/variance/quantile
  cədvəli) reusable funksiya kimi İXRAC EDİLMƏDİ — YALNIZ testlərdə
  (`test_sgs.py`, `test_sgs_integration.py`) ASSERT edilib.
* 100,000+ hüceyrəli tam sahə grid performansı ÖLÇÜLMƏDİ (yalnız 10⁴-ə
  qədər) — **NOT VALIDATED YET**.
* `OrdinaryKriging`-in ÖZÜ (kontinual Kriging/IDW yolu) performans
  baxımından TOXUNULMADI.
