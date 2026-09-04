# Phase A — peşəkar geostatistik interpolyasiya özəyi (A1–A4)

Bu sənəd `imex2d/geology/` altındakı interpolyasiya özəyinin RİYAZİ
müqaviləsini və qərarların səbəblərini yazır. `GEOSTATISTICS.md` (M5)
variogram fitinin tarixçəsidir; bura Phase A-nın nəticəsidir.

## 0. Modul xəritəsi

| Modul | Məsuliyyət |
|---|---|
| `geology/anisotropy.py` | **VAHİD** anizotrop həndəsə: dönmə + miqyaslanma, məsafə, metrik tenzor (A4) |
| `geology/variogram.py` | deneysel γ(h), istiqamətli/şaquli variogram, model fit + doğrulama (A3) |
| `geology/spatial_search.py` | cKDTree indeks, qonşuluq siyasəti, sektor balanslaşdırması, dəstək təsnifatı (A2) |
| `geology/interpolation.py` | Adi Kriging sistemi, dayanıqlı solver, boru xətti, `KrigingResult` (A1/A5) |
| `geology/layer_availability.py` | quyu intervalı → K-lay uyğunlaşdırması, XASSƏ-ÜZRƏ lay mövcudluğu (bax `LAYER_AWARE_MODELING.md`) |

**Lay maskası bu özəyə TOXUNMUR:** lay-məlumatlı rejim (bax
`LAYER_AWARE_MODELING.md`) yalnız HANSI layın HANSI nöqtələrlə
hesablanacağını müəyyən edir; Kriging riyaziyyatı, variogram, anizotropluq
və sərt-data honoring DƏYİŞMƏDƏN qalır.

**Tək həndəsə qaydası (Gate 5):** anizotrop məsafə YALNIZ
`anisotropy.AnisotropyParams.transform()`-dan gəlir. Kriging matrisi,
qonşuluq axtarışı və (istiqamətli) variogram qiymətləndirməsi eyni
obyektin verdiyi fəzada işləyir; ikinci, gizli bir məsafə düsturu YOXDUR.

---

## 1. Kriging formulyasiyası (A1.1)

N qonşu üçün adi kriging sistemi (yarım-dəyişkənlik formasında):

    [Γ  1] [w]   [γ₀]
    [1ᵀ 0] [μ] = [ 1]

* `Γ[i,j] = γ(hᵢⱼ)`, diaqonal **0** (nöqtənin özü ilə fərqi sıfırdır);
* `γ₀[i] = γ(hᵢ₀)` — nugget SAXLANILIR (nugget > 0 olanda kriging dəqiq
  interpolyator olmaqdan çıxır və ölçmə səhvini süzür — bu, qəsdəndir);
* `Σ wᵢ = 1` Laqranj vuruğu `μ` ilə MƏCBUR edilir;
* qiymət `ẑ = Σ wᵢ zᵢ`;
* varians `σ² = Σ wᵢ γ₀ᵢ + μ` — **eyni həll edilmiş sistemdən** oxunur,
  ikinci həll YOXDUR.

Bütün `h` məsafələri anizotrop transformasiyadan sonrakı fəzadadır.

## 2. Variogram konvensiyası (A1.6/A3.4)

    γ(h) = nugget + sill·g(h/a),   g(0)=0, g(∞)=1
    C(h) = (nugget + sill) − γ(h)  (h > 0);  C(0) = nugget + sill

`sill` parametri **quruluşlu** hissədir (nugget DAXİL DEYİL).
`range_` hər üç modeldə **praktiki radiusdur**: γ(range_) ≈ 0.95·(nugget+sill)
(sferikdə dəqiq sill). Modelin öz e-qatlama radiusu `effective_range()`
ilə alınır (eksponensial: `range_/3`, qauss: `range_/√3`) — bu ikisi
QARIŞDIRILMAMALIDIR.

Kriging və variogram eyni `MODEL_FUNCS` funksiyalarını çağırır, ona görə
konvensiya uyğunsuzluğu struktur olaraq mümkün deyil.

## 3. Anizotrop məsafə (A4.1)

    u  = R₀ x        azimut dönməsi  (X,Y → major, minor)
    u' = R_dip u     dip dönməsi     (major, vertical müstəvisində)
    x''= S u'        S = diag(1, a_maj/a_min, a_maj/a_v)
    d_ani = ‖x''(p) − x''(q)‖₂

`azimuth_deg` +Y (Şimal) oxundan saat əqrəbi ilə (0°=Şimal, 90°=Şərq).
`dip_deg` major oxu azimut istiqamətində Z-nin artan tərəfinə əyir.
`dip_deg == 0` (defolt) olanda transform ƏVVƏLKİ (M2) davranışı **bit-bit**
təkrarlayır — dönmə qısa qapanma ilə heç tətbiq edilmir.

`matrix()` (`M = S·R`) və `metric_tensor()` (`G = MᵀM`, `d² = ΔxᵀGΔx`)
eyni həndəsənin cəbri formasıdır — Part D-nin tenzor işi üçün genişlənmə
nöqtəsi. Burada tam tenzor interpolyasiyası **implementasiya edilməyib**
(qəsdən — yanlış tenzor modeli yazmaqdansa təmiz interfeys saxlanılıb).

## 4. Qonşuluq siyasəti (A2)

Ehtiyat zənciri, hər addım `NeighborhoodResult.status`-da görünür:

    anizotrop radius axtarışı            STATUS_RADIUS
      → radius genişləndirməsi           STATUS_RADIUS_EXPANDED
      → k-ən-yaxın ehtiyatı (opt-in)     STATUS_KNN_FALLBACK
      → qlobal ehtiyat (AÇIQ icazə ilə)  STATUS_GLOBAL
      → alınmadı                         STATUS_INSUFFICIENT  → NaN

Sonra: **sektor balanslaşdırması** (kvadrant/oktant), `max_neighbors`
kəsimi, `min_neighbors` yoxlaması.

Süzgəclərin sırası vacibdir: `max_vertical_distance` (XAM |ΔZ| kəsiyi)
namizəd hovuzuna `max_neighbors` kəsimindən ƏVVƏL tətbiq olunur və hovuz
kifayət qədər nöqtə qalana qədər genişləndirilir — əks halda "40 ən
yaxının yalnız 4-ü öz layındandır" vəziyyəti qonşuluğu süni kiçildərdi.
Eyni səbəbdən sektor balanslaşdırması namizəd hovuzunu
`max_neighbors × sektor sayı`-na genişləndirir (`candidate_pool_factor`):
əks halda balanslaşdırılacaq artıq namizəd qalmır.

**Dəstək təsnifatı (A2.8)** "qonşu sayı < N" evristikası DEYİL:

    h = ən_yaxın_məsafə / korrelyasiya_radiusu
    əhatə = dolu kvadrant (2D) / oktant (3D) sayı ÷ ümumi

    qonşu yoxdur ......................... EXTRAPOLATED
    h > 1 ................................ EXTRAPOLATED
    əhatə ≤ ½ (məlumat bir tərəfdə) ...... BOUNDARY
    h > ½ və ya qonşu < 4 ................ WEAK
    qalan ................................ WELL_SUPPORTED

## 5. Sərt datanın honor edilməsi (A1.4)

Siyasət `honor_hard_data` ilə AÇIQDIR:

* `"auto"` (defolt) — yalnız `nugget == 0` olanda dəqiq honor. Nugget > 0
  modelin ölçmə səhvini SÜZMƏSİ deməkdir; ölçülmüş dəyərə geri qaytarmaq
  modelin öz fərziyyəsini pozardı.
* `"always"` — nugget-dən asılı olmayaraq dəqiq honor.
* `"never"` — yalnız sistemin öz nəticəsi.

**Dublikat siyasəti:** eyni koordinatda ziddiyyətli dəyərlər DETERMİNİSTİK
olaraq ORTALANIR və xəbərdarlıq yazılır (əvvəlki sükut "son yazı qazanır"
davranışı əvəzinə). Honor edilən dəyər məhz bu ortalanmış təkil dəyərdir.

Yerli yolda üst-üstə düşmə yoxlaması ƏLAVƏ hesablama tələb etmir —
qonşuluq seçimi ən yaxın qonşunu və məsafəsini onsuz da qaytarır.

## 6. Ədədi ehtiyat zənciri (A1.3)

`_solve_single_robust` — hər pillə statusla görünür:

1. `np.linalg.solve` (LAPACK LU, qismən pivotlama) — `inv()` İSTİFADƏ
   EDİLMİR. → `direct`
2. **Jitter**: məlumat blokunun diaqonalına izə nisbi `ε` (Tixonov).
   Riyazi olaraq sonsuz kiçik nugget əlavə etməyə bərabərdir; dublikat/
   demək olar üst-üstə düşən nöqtələrin təkliyini aradan qaldırır.
   `1e-12 … 1e-4` sırası ilə, İLK uğurlu dayandırır. → `jitter`
3. `lstsq` (minimal-norma/psevdo-tərs); yansızlıq pozularsa çəkilər AÇIQ
   yenidən normallanır. → `lstsq` / `renormalized`
4. **IDW ehtiyatı** `1/d²` — bu, Kriging DEYİL: varians a-priori sill
   kimi qaytarılır, çünki həqiqi kriging variansı mövcud deyil. → `idw_fallback`

Hər həll `Σwᵢ = 1` (`UNBIASED_TOLERANCE = 1e-6`) yoxlamasından keçir.
Toplu (batched) həlldə YALNIZ yoxlamadan keçməyən sistemlər tək-tək
ehtiyat zəncirinə düşür — normal halda heç bir yavaşlama yoxdur.

Qeyri-sonlu (NaN/±inf) sərt data AÇIQ xəbərdarlıqla çıxarılır
(`drop_non_finite=False` isə xəta atır); qeyri-sonlu hədəf üçün NaN
qaytarılır.

## 7. Boru xətti (A5)

    XAM SƏRT DATA
      → DOĞRULAMA (NaN/±inf, uzunluq, dublikat siyasəti)
      → KOORDİNAT NORMALLAŞDIRMASI ((n,2) → (n,3))
      → DENEYSEL VARİOGRAM + MODEL FİT        (auto_fit)
      → PARAMETR DOĞRULAMASI                  (etibarsız model solver-ə ÇATMIR)
      → ANİZOTROP TRANSFORMASİYA
      → MƏKAN İNDEKSİ (cKDTree) + QONŞULUQ SEÇİMİ
      → YERLİ KRİGİNG (toplu həll)
      → QİYMƏT + VARİANS + DİAQNOSTİKA (`KrigingResult`)

## 8. Qlobal → yerli avtomatik keçid (Gate 3)

Heç bir yerli parametr verilməyəndə də nöqtə sayı `auto_local_threshold`
(defolt 100) həddini keçərsə sistem AVTOMATİK yerli rejimə keçir
(`auto_local_max_neighbors`, defolt 40) və bunu `warnings`-ə yazır.

Bu həddən AŞAĞIDA (tipik quyu çoxluğu: 5-50 quyu) nəticə əvvəlki QLOBAL
kriging ilə **birəbir** eynidir — mövcud reqressiya testləri bunu qoruyur.
Qlobal yol da hər hüceyrə üçün YENİDƏN həll etmir: sistem BİR DƏFƏ qurulub
`m` sağ tərəflə həll olunur.

## 9. Aşkarlanmış və düzəldilmiş REAL qüsurlar

Bu iş zamanı testlərlə üzə çıxan (əvvəlki koddakı) səhvlər:

1. **`max_lag`-dan uzaq cütlər sonuncu binə yığılırdı** (`min(bin, n_lags-1)`).
   Sonuncu bin öz mərkəzindən qat-qat böyük məsafələri təmsil edir və cüt
   sayı — yəni fit ÇƏKİSİ — nəhəng olurdu; fit süni uzun radiusa çəkilirdi.
   İndi belə cütlər ATILIR. Bu, istiqamətli radius müqayisəsini tamamilə
   etibarsız edən əsas səbəb idi.
2. **`detect_anisotropy` hər istiqaməti sərbəst (nugget+sill+range) fit
   edirdi.** Sill/nugget radiusu kompensasiya etdiyi üçün fit zəif müəyyən
   idi: izotrop sahədə 100 dəfə fərqli "radiuslar", dönmüş sahədə 45-75°
   azimut xətası (ölçülüb). İndi ümumi sill omnidireksional fitdən SABİT
   götürülür, nugget istiqamətli fitdə 0 saxlanılır — yeganə sərbəst
   parametr radiusdur. 45°/60°/135° dönmələri DƏQİQ bərpa olunur.
3. **`_parameters` `domen/3` evristikası üçün tam (n,n) məsafə matrisi
   qururdu** — 5000 nöqtədə 200 MB. İndi qabarıq örtük diametri (dəqiq,
   `O(n log n)`) və ya parçalı hesablama işlədilir; `n ≤ 1500` üçün əvvəlki
   yol bit-bit saxlanılır.
4. **"Üfüqi" variogramda şaquli sızma.** Laylı məlumatda 500 m aralıqda
   ±5° dip pəncərəsi ±44 m şaquli fərqə, yəni bir neçə LAYA icazə verir.
   `vertical_tolerance` (şaquli bant eni) əlavə edildi.

## 10. Performans (A7)

`python tools/benchmark_interpolation.py --grids 41,101,201`

Domen 5000 m, anizotrop (azimut 35°, major 1200 m, minor 400 m),
`max_neighbors=24`, ölçülmüş (uydurulmamış) nəticələr:

| nöqtə | şəbəkə | hüceyrə | qonşuluq s | kriging s | CƏMİ s | µs/hüceyrə | pik MB | qlobal s |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 41×41 | 1 681 | 0.009 | 0.069 | 0.078 | 46.6 | 69.7 | 0.025 |
| 100 | 101×101 | 10 201 | 0.049 | 0.434 | 0.483 | 47.4 | 149.0 | 0.146 |
| 100 | 201×201 | 40 401 | 0.198 | 1.597 | 1.796 | 44.4 | 171.8 | 0.528 |
| 500 | 41×41 | 1 681 | 0.010 | 0.072 | 0.082 | 48.8 | 69.7 | 0.129 |
| 500 | 101×101 | 10 201 | 0.061 | 0.425 | 0.485 | 47.6 | 149.3 | 0.666 |
| 500 | 201×201 | 40 401 | 0.227 | 1.669 | 1.896 | 46.9 | 172.9 | 2.525 |
| 1 000 | 41×41 | 1 681 | 0.012 | 0.073 | 0.085 | 50.4 | 69.8 | 0.289 |
| 1 000 | 101×101 | 10 201 | 0.065 | 0.437 | 0.503 | 49.3 | 149.4 | 1.370 |
| 1 000 | 201×201 | 40 401 | 0.247 | 1.743 | 1.990 | 49.3 | 173.3 | 5.544 |
| 5 000 | 41×41 | 1 681 | 0.013 | 0.074 | 0.088 | 52.2 | 69.9 | — |
| 5 000 | 101×101 | 10 201 | 0.070 | 0.460 | 0.530 | 51.9 | 149.6 | — |
| 5 000 | 201×201 | 40 401 | 0.262 | 1.811 | 2.073 | 51.3 | 173.6 | — |

Oxunuş:

* **Yerli yolun vaxtı nöqtə sayından demək olar ASILI DEYİL** (100 → 5000
  nöqtə: 1.80 s → 2.07 s, 40 401 hüceyrədə). Mürəkkəblik `O(m·k³ + m log n)`,
  `O(n²)` per-hüceyrə darboğazı YOXDUR (Gate 10).
* Qlobal yol nöqtə sayı ilə kəskin artır (40 401 hüceyrə: 100 → 1000 nöqtə
  = 0.53 s → 5.54 s) və 5000 nöqtədə praktik deyil — məhz buna görə
  istehsal defoltu YERLİDİR.
* Qonşuluq axtarışı ümumi vaxtın ~13%-idir; qalanı yerli sistemlərin
  qurulması/həllidir (gözlənilən nisbət).
* Toplu (`select_batch`) qonşuluq seçimi hədəf başına Python xərcini
  aradan qaldırır: eyni iş üçün 190 µs/hüceyrə → 45 µs/hüceyrə (4.2×,
  ölçülüb).

SGS/SIS bençmarkları (`pytest -m performance`) dəyişməyib — bu axınlar
`spatial_search`-in aşağı qatını (`IncrementalAnisotropicSearch`) işlədir
və müqavilələri toxunulmayıb.

## 11. Qalan iş (Phase A-ya AİD DEYİL)

* Universal/Simple Kriging, co-kriging, Matérn modeli.
* Şaquli variogramın `OrdinaryKriging`-ə avtomatik bağlanması (`range_v`
  hazırda əl ilə verilir və ya `range_`-ə bərabər qəbul edilir).
* Tam keçiricilik tenzoru interpolyasiyası (Part D) — `metric_tensor()`
  genişlənmə nöqtəsi hazırdır, tenzor iş axını YAZILMAYIB.
* Xassə-yönlü çevirmələr (log-keçiricilik, hədli SW/NTG, kateqorik
  fasiya) — `ValueTransform` interfeysi hazırdır, Part B doldurur.
* `ui/panels.py`-də yeni parametrlərin (sektor, dip, ehtiyat siyasəti)
  ekspoz edilməsi.
