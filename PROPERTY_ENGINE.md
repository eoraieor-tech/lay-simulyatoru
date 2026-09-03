# Phase B — rezervuar xassəsi interpolyasiyası + qeyri-müəyyənlik (B1–B5)

Phase A (`INTERPOLATION_CORE.md`) MƏKAN özəyini qurdu: anizotrop məsafə,
qonşuluq, variogram, adi kriginq, kriginq variansı. Bu sənəd onun
ÜZƏRİNDƏKİ XASSƏ qatını yazır.

Phase B-nin bir cümləlik məqsədi:
**POROSITY ≠ PERMEABILITY ≠ SATURATION ≠ NTG ≠ FACIES.**

## 0. Modul xəritəsi

| Modul | Məsuliyyət |
|---|---|
| `geology/transforms.py` | dəyər-fəzası çevirmələri: loq, logit, normal-score + geri-çevirmə semantikası və varians köçürməsi (B1.2-B1.5) |
| `geology/property_config.py` | hər xassə üçün TƏK avtoritar `PropertyStrategy` reyestri (B1.7) |
| `geology/data_quality.py` | QC boru xətti: koordinat, sonluluq, dublikat, fiziki hədd, kənar-dəyər (B4) |
| `geology/property_interpolation.py` | xassə-yönlü interpolyasiya + qeyri-müəyyənlik sistemi (B1/B2) |
| `geology/cross_validation.py` | (genişləndirilib) sızmasız CV + model seçimi (B3) |
| `geology/sgs_ensemble.py` | SGS ansamblı + doğrulama (B5) |

Phase A-nın heç bir riyazi mühərriki TƏKRARLANMIR — məsafə
`anisotropy.py`-dən, qonşuluq `spatial_search.py`-dən, variogram
`variogram.py`-dən, kriginq `interpolation.py`-dən gəlir.

---

## 1. Xassə strategiyası cədvəli

| Xassə | Dəyişən növü | Çevirmə | İnterpolyasiya | Hədlər | Qeyri-müəyyənlik |
|---|---|---|---|---|---|
| PORO | kəsilməz | eynilik | adi kriginq | [0, 1] kəsilir + **sayılır** | kriginq variansı (`identity`) |
| PERMX | loq-normal | `ln(K)` | adi kriginq loq fəzasında | K > 0 | `Var[K] = e^{2ŷ+σ²}(e^{σ²}−1)` — **DƏQİQ** |
| PERMY | loq-normal | `ln(K)` | eyni | K > 0 | eyni (dəqiq) |
| PERMZ | loq-normal | `ln(K)` | eyni | K > 0 | eyni (dəqiq) |
| SW | hədli | `logit` | kriginq logit fəzasında | [0, 1] **riyazi zəmanət** | delta metodu (yaxınlaşma, belə bildirilir) |
| NTG | hədli | `logit` | eyni | [0, 1] **riyazi zəmanət** | delta metodu |
| FACIES | kateqorik | yox | **indikator kriginq** | kateqoriya çoxluğu | ehtimal paylanması + entropiya |

Reyestr: `property_config.DEFAULT_STRATEGIES`; mətn cədvəli
`property_config.strategy_table()`.

---

## 2. Çevirmə riyaziyyatı

### Loq (keçiricilik) — B1.2/B1.3

    Y = ln(K)

Geri-çevirmənin ÜÇ FƏRQLİ mənası var və heç biri "hər yerdə doğru" deyil:

| Mod | Düstur | Nədir |
|---|---|---|
| `MEDIAN` (defolt) | `exp(ŷ)` | şərti **median** — monoton çevirmə altında dəyişməzdir |
| `MEAN` | `exp(ŷ + σ²/2)` | loq-normal şərti **orta** |
| `MEAN_OK` | `exp(ŷ + σ²/2 − μ)` | **adi kriginq** üçün düzəlişli orta (`μ` — Laqranj vuruğu) |

`MEAN_OK` üçün Phase A-nın `KrigingResult.lagrange` sahəsi əlavə edildi
(Journel & Huijbregts §7.4; sadə kriginqdə `μ = 0`, onda `MEAN`-ə çevrilir).

Defolt `MEDIAN`-dır — `exp(ŷ + σ²/2)` HƏR YERDƏ sükutla TƏTBİQ EDİLMİR.

### Logit (hədli xassələr) — B1.4/B1.5

    p  = (z − lo)/(hi − lo)
    p̃  = ε + p·(1 − 2ε)
    y  = ln(p̃/(1 − p̃))

Geri çevirmə sıxmanı **AÇIR**, ona görə `inverse(forward(z)) == z` maşın
dəqiqliyi ilə — DAXİL OLMAQLA dəqiq `0` və `1`. Yəni sərt datanın dəqiq
honor edilməsi hədlərdə də pozulmur.

`ε = 1e-4`-in TƏSİRİ (sənədləşdirilib): yalnız hədlərin logit fəzasında
nə qədər uzağa düşdüyünü təyin edir (`y(0) ≈ −9.21`). Geri çevirmə
logistik funksiya olduğu üçün `[lo, hi]`-dən **riyazi olaraq çıxa bilmir**
— "Gauss kriginq SW = 1.07 verdi" problemi kökündən yoxdur.

**ÖLÇÜLMÜŞ MƏHDUDİYYƏT:** dəyərlərin bir hissəsi DƏQİQ hədddə yığılıbsa
(censored), logit onları eyni ekstremal nöqtəyə göndərir və nəticə
pisləşir. 60 quyu, 8% dəyər `NTG = 1.0`-da: logit R² = −0.12, çevirməsiz
R² = +0.13. Məhz buna görə `default_candidates()` hədli xassələr üçün də
ÇEVİRMƏSİZ namizədi daxil edir və qərarı DOĞRULAMAYA buraxır.

### Normal-score (SGS)

`gaussian_transform.NormalScoreTransform` `ValueTransform` interfeysinə
uyğunlaşdırılıb. **DATADAN ASILIDIR** (`data_dependent=True`) — ona görə
çarpaz-doğrulamada hər qat üçün YENİDƏN fit olunur.

---

## 3. Qeyri-müəyyənlik semantikası (B2.1)

Yeddi fərqli kəmiyyət QARIŞDIRILMIR:

| Sahə | Nədir |
|---|---|
| `transformed_variance` | ÇEVRİLMİŞ fəzada ƏSL kriginq variansı `σ²_y` |
| `variance` | orijinal fəzaya köçürülmüş varians |
| `variance_kind` | köçürmənin növü: `exact` (loq-normal), `delta` (logit/normal-score), `identity`, `undefined` |
| `nearest_distance` | anizotrop fəzada ən yaxın sərt data |
| `data_density` | korrelyasiya radiusundakı nöqtə SAYI — **diaqnostika**, ehtimal deyil |
| `support` | Phase A həndəsi dəstək təsnifatı |
| `confidence` | HIGH/MEDIUM/LOW/EXTRAPOLATED — **İNTERPRETASİYA**, kalibrlənmiş ehtimal DEYİL |
| ansambl `variance` (SGS) | realizasiyalar arası dəyişkənlik — kriginq variansından FƏRQLİ |

`confidence` qaydası (deterministik, `classify_confidence`):

    dəstək EXTRAPOLATED / qonşu yoxdur ..... EXTRAPOLATED
    dəstək BOUNDARY / WEAK ................. LOW
    σ²/(nugget+sill) > 0.70 ................ LOW
    σ²/(nugget+sill) > 0.30 ................ MEDIUM
    qalan .................................. HIGH

Variansı olmayan üsul (IDW/NearestNeighbour) üçün `variance` **NaN**
qalır və bu, xəbərdarlıqla bildirilir — ədəd UYDURULMUR.

Standartlaşdırılmış xəta (B2.2) ÇEVRİLMİŞ fəzada hesablanır (varians məhz
orada əsl kriginq variansıdır): `e = (z − ẑ)/σ`; ideal `orta ≈ 0`,
`dispersiya ≈ 1`. `coverage_68`/`coverage_95` ideal `0.683`/`0.954`.

---

## 4. Məlumat keyfiyyəti (B4)

Sıra (kənar-dəyər ƏN SONDA — qalan təmiz paylanmaya nəzərən ölçülməlidir):

    koordinat → sonluluq → dublikat → fiziki hədd → kənar-dəyər

**ÜÇ FƏRQLİ PROBLEM AYRILIR:**

1. **Fiziki cəhətdən etibarsız** (mənfi məsaməlik, `SW = 1.4`) → defolt
   ÇIXARILIR, səbəb yazılır.
2. **Statistik kənar-dəyər** (5000 mD) → defolt YALNIZ İŞARƏLƏNİR.
   `remove_outliers=True` AÇIQ seçim tələb edir.
3. **Nadir amma etibarlı** (çat zonası) → heç bir avtomatik qərar yoxdur.

Dublikat siyasətləri: `mean`/`median`/`keep_first`/`keep_last`/
`keep_separate`/`majority`/`raise`. Ziddiyyətli dublikatlar ayrıca sayılır.

Kənar-dəyər metodları: `MAD`, `IQR`, `SPATIAL` (yerli qonşuluğa nəzərən —
Phase A-nın `NeighborhoodSelector`-u ilə, EYNİ anizotrop həndəsədə).

---

## 5. Çarpaz-doğrulama və model seçimi (B3)

**SIZMA MEMARLIQLA qarşısı alınır:** hər qatda
`interpolate_property_field()` YALNIZ təlim alt-çoxluğu ilə çağırılır və
o, QC-ni, çevirməni, variogram fitini ÖZÜ həmin nöqtələrdən qurur. Test
nöqtəsi heç bir mərhələyə fiziki olaraq daxil olmur.

Doğrulama dizaynları: `loo`, `random_kfold`, `spatial_block`.
Məkan blokları BÜTÖV bölgələri gizlədir — təsadüfi bölgüyə əsaslanıb
"məkan ümumiləşdirməsi yaxşıdır" demək olmaz (ölçülüb: eyni datada
məkan-blok RMSE-si təsadüfi bölgüdən BÖYÜKDÜR).

Çox-metrikli sıralama (`DEFAULT_SELECTION_WEIGHTS`, cəmi 1.0):

| Meyar | Çəki | Səbəb |
|---|---|---|
| accuracy | 0.50 | proqnoz xətası əsas məqsəddir, amma tək meyar deyil |
| calibration | 0.25 | qeyri-müəyyənlik yanlışdırsa aşağı RMSE də aldadıcıdır |
| bias | 0.15 | sistematik sürüşmə həcm hesablamasını təhrif edir |
| validity | 0.07 | fiziki hədd pozması uyğunsuzluq əlamətidir |
| stability | 0.03 | alınmayan proqnozlar üçün kiçik, amma sıfır olmayan cəza |

Bütün cərimələr ƏN YAXŞI namizədə NİSBİDİR, ona görə metriklərin miqyası
(mD vs kəsr) nəticəyə təsir etmir. Çəkilər sabit kodlanmır.

**Ölçülə bilməyən meyar → ƏN PİS müşahidə olunan cərimə** (neytral 0
DEYİL). Səbəb ölçülüb: neytral 0 ilə IDW (kriginq variansı yoxdur) balın
25%-ini pulsuz qazanırdı və PORO ssenarisində RMSE-si DAHA PİS olmasına
baxmayaraq (0.02684 vs 0.02594) QALİB gəlirdi — yəni qeyri-müəyyənliyi
ümumiyyətlə verə bilməyən üsul, verə bilən üsuldan üstün tutulurdu.
"Yoxlaya bilmirik" → "kredit vermirik".

Faktiki nəticələr (LOOCV, sintetik sahələr):

    PORO (70 quyu)                      bal      RMSE      var(e)
      1. kriginq + exponential        0.2309   0.02594     1.912
      2. kriginq + spherical          0.6414   0.02643     3.522
      3. IDW                         11.0664   0.02684     yox
      4. kriginq + gaussian          11.2600   0.03710    45.180

    PERMX (80 quyu, loq-normal sahə)    bal      RMSE      var(e)
      1. kriginq + exponential (log)  0.0652   401.79      1.164
      2. kriginq + spherical  (log)   0.2927   400.91      2.095
      3. kriginq + spherical (xam)    2.3156   437.83     10.058
      4. kriginq + gaussian   (log)   5.7689   447.27     23.809
      5. IDW                          5.7758   439.78      yox

Qalib hər iki halda ƏN YAXŞI KALİBRLƏNMİŞ modeldir (`var(e)` 1-ə ən
yaxın) — sırf ən kiçik RMSE deyil. Loq fəzası xam fəzadan üstün gəlir,
qauss modeli isə pis kalibrləməyə görə cəzalanır.

Kateqorik metriklər: dəqiqlik, log-loss, Brier, qarışıqlıq matrisi.
**Qeyd:** indikator kriginq bəzən SƏHV kateqoriyaya ~1 ehtimal verir və
`−ln(p)` partlayır — ona görə müqayisədə BRIER (məhdud) üstünlük təşkil
edir, log-loss isə hesabatda saxlanılır.

---

## 6. SGS finalizasiyası (B5)

`sgs_ensemble.simulate_sgs_ensemble(n, …, base_seed)`:

* seed konvensiyası `run_realizations_sgs`-in EYNİSİDİR (`base_seed + i·1000`);
* eyni seed → **bit-bit eyni** ansambl; fərqli seed → fərqli realizasiyalar;
* `mean`/`variance`/`std`/`p10`/`p50`/`p90`.

**Terminologiya (B5.7):** `P10/P50/P90` ANSAMBL KVANTİLLƏRİDİR, "etibar
intervalı" DEYİL.

`validate_realization()` üç meyarı ÖLÇÜR: sərt data (dəqiq), marjinal
paylanma (KS + momentlər), variogram (realizə/hədəf radius nisbəti).

---

## 7. Aşkarlanmış və düzəldilmiş REAL qüsurlar

Phase B-nin testləri Phase A-da və öz kodunda dörd həqiqi səhv tapdı:

1. **Qauss variogramı + sıfıra yaxın nugget = yararsız kriginq.** Qauss
   modeli başlanğıcda parabolikdir; sıfır nugget-lə matris demək olar
   təkil olur, çəkilər ossilyasiya edir. Sistem "uğurla" həll olunur
   (`Σw = 1` ödənir, solver `direct` bildirir), amma qiymət məlumat
   diapazonundan çıxır. **Ölçülmüş:** məlumat `[0.151, 0.210]`, nəticə
   `[−0.300, +1.126]`. Düzəliş: `variogram.stabilizing_nugget()` — qauss
   modelində nugget sillin ən azı 1%-inə qaldırılır (GSLIB tövsiyəsi),
   AÇIQ xəbərdarlıqla. Eksponensial/sferik modellər hər nugget üçün
   sabitdir (ölçülüb).
2. **Ədədi sabitləşdirici nugget sərt datanın honor edilməsini
   söndürürdü.** `honor_hard_data="auto"` "nugget varsa honor etmə"
   deməkdir — amma SIRF ƏDƏDİ requlyarlaşdırma "ölçmədə səhv var"
   demək DEYİL. Düzəliş: `model_nugget` (sabitləşdirmədən ƏVVƏLKİ) ayrıca
   saxlanılır və honor qərarı ONA baxır.
3. **MAD = 0 kənar-dəyər aşkarlanmasını məhz TƏMİZ datada söndürürdü.**
   29 nöqtədə `K = 50`, birində `5000` → MAD = 0 → "kənar-dəyər yoxdur".
   Düzəliş: `_robust_scale()` iki pilləlidir (MAD, sonra orta mütləq
   kənarlanma); yalnız BÜTÜN dəyərlər eyni olanda sıfır qaytarır.
4. **SGS variogram doğrulaması yanlış pəncərədə ölçürdü.** Realizə
   variogramı (a) sərt-data hüceyrələrini də daxil edirdi, (b) domenin
   74%-inə qədər lag işlədirdi, (c) nugget/sill/radiusu birlikdə fit
   edirdi. Nəticə: hədəf 220 m üçün "735 m", hədəf 400 m üçün "3542 m".
   Düzəliş: yalnız simulyasiya edilmiş hüceyrələr + etibarlı lag pəncərəsi
   (`min(2.5·radius, domenin yarısı)`) + SABİT sill. Nəticə: 220 → 243,
   120 → 110, 400 → 534.

---

## 8. Performans (B10)

`python tools/benchmark_property_pipeline.py`

PERMX (loq fəzası), domen 5000 m, radius 1200 m — ölçülmüş:

| nöqtə | şəbəkə | hüceyrə | QC s | variogram s | interp s | µs/hüceyrə | LOOCV s | SGS s | yaddaş MB |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 41×41 | 1 681 | 0.0015 | 0.035 | 0.115 | 68 | 4.19 | 1.45 | 56 |
| 100 | 101×101 | 10 201 | 0.0006 | 0.035 | 0.493 | 48 | 3.89 | 8.97 | 122 |
| 500 | 41×41 | 1 681 | 0.0008 | 0.052 | 0.133 | 79 | 29.2 | 1.46 | 56 |
| 500 | 101×101 | 10 201 | 0.0008 | 0.052 | 0.512 | 50 | 27.8 | 9.03 | 122 |
| 1 000 | 41×41 | 1 681 | 0.0011 | 0.092 | 0.179 | 107 | 99.2 | — | 56 |
| 1 000 | 101×101 | 10 201 | 0.0013 | 0.094 | 0.564 | 55 | 96.0 | — | 122 |
| 5 000 | 41×41 | 1 681 | 0.0032 | 0.456 | 0.543 | 323 | — | — | 302 |
| 5 000 | 101×101 | 10 201 | 0.0028 | 0.451 | 0.960 | 94 | — | — | 302 |

Oxunuş:

* **QC praktiki olaraq pulsuzdur** (≤ 3 ms hətta 5000 nöqtədə).
* **İstehsal yolu YERLİ qalır**: 10 201 hüceyrədə 48-94 µs/hüceyrə,
  nöqtə sayından zəif asılı (Phase A-nın xassəsi qorunub).
* **LOOCV `O(n)` dəfə bahadır** (1000 nöqtə → 96 s) — bu, MODEL SEÇİMİ
  üçündür və istehsal şəbəkə yolunda YOXDUR. Cədvəldə qəsdən ayrı sütun.
* **Cüt alt-nümunəsi** (`variogram.MAX_VARIOGRAM_PAIRS = 2·10⁶`) 5000
  nöqtədə variogramı 1.57 s → 0.46 s və yaddaşı 858 MB → 302 MB endirdi
  (deterministik, `n ≲ 2000` üçün heç nə dəyişmir — bütün mövcud testlər
  o həddin altındadır).
* Phase A-nın öz bençmarkı DƏYİŞMƏYİB: 5000 nöqtə / 10 201 hüceyrə üçün
  0.566 s, 55.5 µs/hüceyrə, 150 MB.

---

## 9. Geriyə uyğunluq (B11)

Bütün Phase A/Phase 4-5 ictimai adları proqramla yoxlanılıb və yerindədir:

* `interpolation`: `NearestNeighbour`, `InverseDistance`, `OrdinaryKriging`,
  `INTERPOLATORS`, `interpolate_property`, `KrigingResult`,
  `ValueTransform`, `LogTransform` (sonuncu ikisi `transforms.py`-dən
  RE-EKSPORTDUR — **iki implementasiya saxlanmır**);
* `variogram`: 14 ad, `AnisotropyParams` daxil;
* `spatial_search`: `AnisotropicNeighborSearch`, `IncrementalAnisotropicSearch`;
* `sgs`: `simulate_sgs`, `run_realizations_sgs`,
  `simulate_sgs_facies_conditioned`, `PropertyVariogramParams`;
* `facies`: `simulate_sis`, `observed_proportions`, `indicator`;
* `cross_validation`: `leave_one_out`, `k_fold`, `CrossValidationResult`.

`application/geology_service.py` və `application/geology_adapter.py`
**dəyişdirilməyib**. `ui/panels.py` dəyişdirilməyib.

---

## 10. Qalan iş (Part C/D)

* Tam keçiricilik tenzoru interpolyasiyası + MPFA-O inteqrasiyası (Part D).
  `anisotropy.metric_tensor()` genişlənmə nöqtəsi hazırdır.
* Ko-kriginq / universal kriginq (trendli) — Part C/D.
* Sensor/təzyiq datası ilə çoxdəyişənli şərtləndirmə.
* Yeni xassə mühərrikinin `geology_service.py`/UI-a bağlanması: hazırda
  `PropertyEstimate`/`CategoricalEstimate` backend-də tam hazırdır, amma
  `WellBasedGeologicalModelBuilder` hələ Phase A-nın `interpolate_property()`
  yolundan istifadə edir (geriyə uyğunluq üçün qəsdən toxunulmayıb).
* Yığılmış (censored) hədli xassələr üçün truncated/qarışıq model.
