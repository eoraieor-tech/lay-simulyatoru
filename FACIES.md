# Fasiya modelləşdirməsi — Sequential Indicator Simulation (Phase 4)

## Audit nəticəsi: əvvəllər nə var idi?

Heç nə. Kod bazasında `facies`/`SIS`/indikator simulyasiyası sözünün
BİR DƏFƏ də keçmədiyi təsdiqləndi (bu sessiyanın əvvəlki audit-ləri +
təzə axtarış). Kateqorik xassə üçün ən yaxın analoq `domain/structure.py`-
dəki `RegionSet` idi (SATNUM-bənzər tam ədəd kod + ad lüğəti) — amma bu,
SCAL/PVT region seçimi üçündür, fasiya üçün İSTİFADƏ EDİLMİR (fərqli
məqsəd, ayrıca saxlanıldı). `geology_service.py`-nin mövcud iş axını
(`WellBasedGeologicalModelBuilder`) fərz edilən bir `"FACIES"` sütununu
DİGƏR bütün rəqəmsal sütunlar kimi (`DEFAULT_RULES.get(source, Property
Rule(source))` — heç bir xüsusi qayda yoxdur) birbaşa Kriging/IDW-yə
göndərəcəkdi, yəni fasiya kodu 0/1/2 arasında "1.4" kimi mənasız aralıq
dəyər çıxara bilərdi. **Bu, DƏYİŞDİRİLMƏDİ** — `geology_service.py`-a
TOXUNULMADI (tapşırıq flow solver/PVT/SCAL-a toxunmamağı tələb etdi, və
facies-in tam iş axınına inteqrasiyası Phase 5-in işi kimi görünür). Bu
o deməkdir ki, əgər kimsə bu gün `WellBasedGeologicalModelBuilder`
vasitəsilə "FACIES" sütunlu CSV yükləsə, HƏLƏ DƏ köhnə (səhv) yol
işləyəcək — YENİ SIS modulu bundan TAMAMILƏ AYRI, paralel bir imkandır.
**NOT DONE**: SIS-in `geology_service.py`/`GeologicalModel.property_maps`
iş axınına drop-in əvəz kimi bağlanması.

## Memarlıq: nə YENİDƏN YAZILMADI

`imex2d/geology/facies.py` YALNIZ üç şey əlavə edir:

1. İndikator çevirməsi (`indicator()`), nisbət yoxlaması/müqayisəsi.
2. Ardıcıl simulyasiya yolu + kateqorik nümunələmə (SIS-in özəl hissəsi).
3. Ehtimal normallaşdırma diaqnostikası (kliplənmə/renormallaşma sayğacı).

Aşağıdakılar **BİRBAŞA, DƏYİŞİKLİKSİZ** Phase 2-3-dən istifadə olunur:

* **Variogram riyaziyyatı və fit** — `geology/variogram.py`-dəki
  `fit_variogram_from_data()` (deneysel variogram + sferik/eksponensial/
  qauss model fit, çəkili ən kiçik kvadrat). Hər fasiyanın indikator
  variogramı bu funksiya ilə fit edilir — YENİ variogram düsturu
  YAZILMADI.
* **Kriging həlli və anizotropluq** — `geology/interpolation.
  OrdinaryKriging` BİRBAŞA çağırılır (hər fasiya üçün ayrıca instans).
  `azimuth_deg`/`range_minor`/`range_v` (tam 3D anizotropluq),
  `search_radius`/`max_neighbors`/`min_neighbors` (yerli axtarış) — HAMISI
  Phase 2-3-dən miras qalır, TƏKRAR YAZILMADI.

## SIS alqoritmi (tapşırıqda tələb olunan addımlar)

```
1. seed-dən np.random.default_rng(seed) qurulur
2. simulyasiya ediləcək hüceyrələr (sərt-data ilə üst-üstə düşməyənlər)
   TƏSADÜFİ sıra ilə (rng.permutation) seçilir
3. hər hüceyrə üçün:
   a. hər fasiya k üçün: I_k(mövcud kondisioner nöqtələr) hesablanır
   b. OrdinaryKriging(range_, azimuth_deg, range_minor, search_radius,
      max_neighbors, ...) bu hüceyrədə şərti ehtimalı EVALƏNDİRİR
      (bu, "yerli axtarış" + "indikator kriging sistemi" addımlarıdır)
   c. ehtimallar [0, ∞)-dan [0, 1]-ə KƏSİLİR (mənfi kriging çəkiləri
      mümkündür — klassik indikator-kriging artefaktı), sonra CƏMƏ
      bölünərək normallaşdırılır
   d. rng.choice(fasiyalar, p=normallaşmış_ehtimallar) ilə NÜMUNƏ
      GÖTÜRÜLÜR (ən yüksək ehtimallı fasiya AVTOMATİK seçilmir!)
   e. seçilmiş fasiya kondisioner çoxluğuna ƏLAVƏ olunur (növbəti
      hüceyrələr bunu "məlumat" kimi görür)
4. bütün hüceyrələr simulyasiya olunana qədər davam
```

`P(nümunə = ən-yüksək-ehtimallı-fasiya) ≠ 1` — sınanıb
(`test_run_realizations_ids_and_seeds_are_deterministic_sequence` və
`test_proportions_preserved_on_average_across_many_realizations`
realizasiyaların FƏRQLİ olduğunu göstərir, "həmişə ən yüksək ehtimallı"
olsaydı hamısı EYNİ olardı).

## Ehtimal etibarlılığı — "səssizcə kliplənmə" YOXDUR

`domain/validation.py`-dəki "heç vaxt gizli problemi göstərmə" qaydasına
uyğun olaraq:

* Kriging NaN qaytarsa (yerli axtarışda kifayət qədər qonşu yoxdursa) —
  o hüceyrə üçün QLOBAL/REGIONAL nisbətlərə keçilir, sayğac artırılır,
  sonda ÜMUMİ say XƏBƏRDARLIQ kimi bildirilir (`fallback_events`).
* Ehtimalların cəmi kliplənmədən ƏVVƏL 1-dən **böyük** kənarlaşırsa
  (>`LARGE_RENORMALIZATION_THRESHOLD=0.05`) — bu da sayılır və
  bildirilir (`renormalization_events`) — KİÇİK kənarlaşma (məs. 1.02)
  normaldır (kriging-in qərəzsizlik şərti dəqiq deyil, tolerantlıqla),
  YALNIZ BÖYÜK kənarlaşma diaqnostik dəyərə malikdir.

## Sərt data (hard data)

`_find_hard_data_matches()` hər hədəf hüceyrəni sərt data nöqtələri ilə
(tolerantlıq `hard_data_tolerance=1e-6`, defolt) MÜQAYİSƏ edir — üst-üstə
düşənlər HEÇ VAXT simulyasiya yoluna daxil edilmir, birbaşa müşahidə
dəyərini alır. Sınanıb: `test_hard_data_is_honored_exactly_regardless_of_seed`
(4 fərqli seed, hamısında sərt data dəyişməz).

## Fasiya nisbətləri

`FaciesProportions`: qlobal (məcburi) + region-əsaslı (istəyə görə) +
lay-əsaslı (istəyə görə), prioritet lay > region > qlobal. Hər üçü
AYRICA `validate_facies_proportions()` ilə yoxlanılır (cəm=1 ±1e-6,
mənfi/>[0,1] rədd edilir). **İSTİFADƏÇİNİN verdiyi nisbətlər HEÇ VAXT
avtomatik normallaşdırılmır** — cəmi 1 deyilsə `ValueError`, "özün
düzəlt" mesajı ilə. Müşahidə olunan (sərt datadan) və tələb olunan
nisbətlər arasında böyük fərq (>0.15, konfiqurasiya edilə bilər) YALNIZ
XƏBƏRDARLIQ doğurur, DƏYƏR DƏYİŞMİR.

## Anizotropluq

Ayrıca test edilib: (1) birbaşa kriging ehtimalı — güclü anizotropluqla
major-ox qonşusu minor-oxdan çox təsir edir
(`test_anisotropic_indicator_kriging_favours_major_axis_neighbour`); (2)
ANSAMBL statistikası — 60 realizasiya üzərində major-ox boyu hüceyrə
cütü minor-ox boyu cütdən DAHA TEZ-TEZ eyni fasiyada olur (198/200 vs
138/200 geniş sınaqda, bax test şərhi) —
`test_ensemble_spatial_continuity_is_stronger_along_major_axis`. Bu,
tapşırığın "TƏK hüceyrə deyil, ansambl statistikası" tələbinə düz
uyğundur (§11).

## Random seed / təkrarlana bilənlik

`np.random.default_rng(seed)` (kod bazasının hər yerdə işlətdiyi
konvensiya, bax `geology/cross_validation.k_fold`). EYNİ (data,
parametr, grid, seed) → EYNİ realizasiya (sınanıb,
`test_same_seed_produces_identical_realization`). FƏRQLİ seed → FƏRQLİ
realizasiya (`test_different_seeds_produce_different_realizations`).

## Çoxlu realizasiya memarlığı

`run_realizations(n, ..., seed=0)` → hər realizasiya
`seed + i*1000` (mövcud `application/scenarios.py` konvensiyası ilə
EYNİ). API `realization_id`/`seed`/`codes` daşıyan `FaciesRealization`
qaytarır — N=1-dən N=100+-a QƏDƏR memarlıq DƏYİŞMƏDƏN işləyir (sadə
Python siyahı-anlayışı, paralellik YOXDUR, bilərəkdən — bax aşağı).

`summarize_realized_proportions()` — tapşırıq §12: hər fasiya üçün
tələb olunan + reallaşan nisbətlərin mean/std/min/max-ı. **TƏK
realizasiyanın DƏQİQ tələb olunan nisbəti verməsi GÖZLƏNİLMİR** —
yalnız ORTALAMA yaxınlaşmalıdır (test tolerantlığı ±0.15, 30 realizasiya
üzərində).

## Mürəkkəblik (sənədləşdirilib, ƏVVƏLCƏDƏN optimallaşdırılmayıb)

n = simulyasiya ediləcək hüceyrə sayı, K = fasiya sayı, m =
`max_neighbors` (sabit, grid ölçüsündən ASILI DEYİL).

* Hər hüceyrə üçün: K dəfə `OrdinaryKriging.interpolate()` çağırılır.
* Hər çağırış: (a) məsafə matrisi — CƏMİ kondisioner nöqtə sayı ilə
  XƏTTİ, bu say hüceyrə-hüceyrə BÖYÜYÜR (0-dan n-ə qədər) → **O(n)**
  hər çağırışda; (b) ən yaxın `m` qonşunun seçilməsi (sıralama); (c)
  `(m+1)×(m+1)` xətti sistemin həlli → **O(m³)**, m SABİT olduğu üçün
  bu HƏMİŞƏ ucuzdur.
* **Ümumi**: `O(K · n · (n + m³)) ≈ O(K·n²)` böyük n üçün (məsafə
  axtarışı üstünlük təşkil edir, m³ sabit qalır).

**225 hüceyrəlik (15×15) bençmarkda** (`test_small_synthetic_sis_
benchmark_completes_and_is_self_consistent`) icra vaxtı sanıyələrlədir
(tam test faylı — 24 test, bir neçə yüz kriging çağırışı daxil olmaqla —
~3.5 saniyədə tamamlanır). **Bu, TƏSADÜFİ deyil** — `OrdinaryKriging`
HƏR çağırışda kondisioner nöqtələr üzərində YENİDƏN sıfırdan axtarış
aparır (keşləmə/inkremental yeniləmə YOXDUR, bax Phase 2-3 audit-i) —
minlərlə hüceyrəli TAM SAHƏ grid-ində (məs. 41×41×5 ≈ 8400 hüceyrə) bu
`O(K·n²)` termini nəzərəçarpacaq şəkildə yavaşlaya bilər. **Optimallaşdırma
BU FAZADA edilmədi** (tapşırıq "əvvəlcə düzgünlük" tələb etdi) — gələcək
optimallaşdırma yolu: k-d ağacı ilə İNKREMENTAL qonşu axtarışı (hər
addımda YENİ nöqtəni ağaca əlavə etmək, sıfırdan axtarmaq YOX), bu,
axtarış xərcini `O(n)`-dən `O(log n)`-ə endirər.

## Elmi çəkincə (tapşırıq §15)

Bu modul PROQRAM CƏHƏTDƏN düzgün SIS-i tətbiq edir — bu, SİMULYASİYA
NƏTİCƏSİNİN yerin gerçək fasiya paylanmasını əks etdirdiyi demək DEYİL.
Testlər yoxladığı: indikator çevirməsinin düzgünlüyü, ehtimal
normallaşdırmasının riyazi etibarlılığı, sərt-data hörməti, seed
təkrarlanabilənliyi, ANİZOTROP davamlılığın ANSAMBL SƏVİYYƏSİNDƏ
gözlənilən istiqamətdə olması. Testlərin YOXLAMADIĞI: real geoloji
"doğruluq" (bu, YALNIZ real quyu sıxlığı, ekspert geoloji şərh və
kənar doğrulama ilə qiymətləndirilə bilər — sintetik testlə YOXLANILA
BİLMƏZ).

## Sınaqlar

`tests/test_facies.py` — 24 test, tapşırıqda sadalanan 14 kateqoriyanın
HAMISI: indikator çevirməsi, nisbət yoxlanması (+ region/lay prioriteti),
indikator variogram (Phase 2-3 vasitəsilə), ehtimal normallaşdırması
(həm birbaşa, həm SIS daxilində), sərt data hörməti (4 seed), seed
təkrarlanabilənliyi, fərqli seed-lərin fərqli nəticə verməsi, anizotrop
davamlılıq (nöqtəvi + ansambl), çoxfasiyalı (3 fasiya), tək-fasiya kənar
halı, seyrək sərt data (ehtiyat evristika xəbərdarlığı ilə), çoxlu
realizasiya üzrə nisbət saxlanması (mean/std/min/max), etibarsız giriş
(3 ayrı hal), kiçik sintetik bençmark (15×15=225 hüceyrə).

## Qalan iş / NOT DONE

* **SIS-in `geology_service.py`/`GeologicalModel` iş axınına inteqrasiyası**
  — bax yuxarı, "Memarlıq" bölməsi. Hazırkı modul MÜSTƏQİL API-dir.
* Performans optimallaşdırması (k-d ağacı, inkremental axtarış) —
  bilərəkdən edilmədi (bax "Mürəkkəblik").
* SGS (Sequential Gaussian Simulation), fasiya-şərtli poroziya/keçiricilik
  paylanması, MPFA, tarixi uyğunlaşdırma, EOR, qeyri-müəyyənlik
  kəmiyyətləndirməsi — tapşırıqda AÇIQ İSTİSNA EDİLİB, bu fazada YOXDUR.
* Paralel realizasiya icrası YOXDUR (sadə ardıcıl Python siyahı-anlayışı)
  — tapşırıq bunu bu fazada tələb ETMƏDİ.
