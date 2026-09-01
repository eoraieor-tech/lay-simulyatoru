# SIS boru xətti inteqrasiyası + məkan axtarışı bünövrəsi (Phase 4.1)

## 1. Audit: FACIES sütunu əvvəllər HANSI yolla sistemə girirdi?

Tam yol (kod oxunaraq təsdiqlənib, TƏXMİN EDİLMƏYİB):

    CSV (`well_data_io.read_well_csv`)
        → `WellSample.values["FACIES"] = <float>`   (istənilən ədədi sütun kimi)
        → `WellDataset.property_names()` "FACIES"-i sadəcə BAŞQA sütun kimi qaytarır
        → `WellBasedGeologicalModelBuilder.build()`:
              rule = self.rules.get("FACIES", PropertyRule("FACIES"))   # DEFAULT_RULES-da YOX idi
              values = self._interpolate_volume(...)                    # KƏSİLMƏZ Kriging/IDW
              model.add_property(PropertyMap.from_array("FACIES", values, ...))

Yəni FACIES `interpolate_property()`-ə (Kriging/IDW) BİRBAŞA gedirdi —
0/1/2 kodları arasında "1.4" kimi mənasız ədəd ala bilərdi. Bu, HEÇ VAXT
düzəldilməmişdi (Phase 4 standalone SIS bunu HƏLL ETMİRDİ, sadəcə AYRI
bir imkan kimi mövcud idi). **Bu fazanın əsas işi məhz bu yolu
DƏYİŞDİRMƏKDİR** (bax aşağı, "Marşrutlaşdırma").

## 2. Fayllar

**Yeni:**
- `imex2d/geology/property_types.py` — CONTINUOUS/CATEGORICAL reyestri.
- `imex2d/geology/spatial_search.py` — `AnisotropicNeighborSearch` (sabit),
  `IncrementalAnisotropicSearch` (ARDICIL, SIS üçün).
- `imex2d/geology/hard_data.py` — hüceyrə xəritələnməsi + ziddiyyət aşkarlanması.
- `imex2d/domain/facies_field.py` — `FaciesField` (RegionSet-dən AYRI).
- `tests/test_property_types.py`, `test_spatial_search.py`, `test_hard_data.py`,
  `test_facies_integration.py`, `test_facies_performance.py`.

**Dəyişdirilmiş:**
- `imex2d/geology/facies.py` — `use_fast_search` (defolt `True`),
  `FaciesDiagnostics` (ayrıca sayğaclar).
- `imex2d/application/geology_service.py` — `FaciesBuildConfig`,
  kateqorik marşrutlaşdırma, `_simulate_categorical_field`.
- `imex2d/domain/geological_model.py` — `facies_fields` sahəsi.
- `imex2d/ui/panels.py` — `FaciesPanel` (minimal, bax §9).

## 3. Xassə növü arxitekturası

`geology/property_types.py`: `PropertyType.CONTINUOUS`/`CATEGORICAL`
reyestri, `classify_property(name, overrides=None)`. Naməlum ad →
`CONTINUOUS` (geriyə uyğun, mövcud davranış DƏYİŞMİR). Reyestr
GENİŞLƏNƏ BİLƏNDİR (`overrides` arqumenti) — sütun adı ilə səpələnmiş
`if name == "FACIES"` yoxlaması İSTİFADƏ OLUNMUR.

## 4. SIS inteqrasiya arxitekturası

`build()`-də hər sütun ÜÇÜN ƏVVƏLCƏ `classify_property()` çağırılır:

    CATEGORICAL → _simulate_categorical_field() → model.add_facies_field()
    CONTINUOUS  → _interpolate_volume() (DƏYİŞMƏYİB) → model.add_property()

Kateqorik sütun HEÇ VAXT `_interpolate_volume`-a ötürülmür (kod
səviyyəsində ayrılıb, "bəzən" yox). Bütün laylar BİR SIS çağırışında
(tam 3D X,Y,Z) işlənir — hər lay üçün ayrıca DEYİL, çünki şaquli
variogram (`range_v`) artıq laylar-arası davamlılığı modelləşdirir.

## 5. Sərt-data xəritələnməsi

`hard_data.map_samples_to_cells()` (`xy_to_ij`/`depth_to_k` — MÖVCUD,
TƏKRARLANMAYAN funksiyalar) hər nümunəni `(i,j,k)` hüceyrəsinə
xəritələyir. `detect_hard_data_conflicts()` eyni hüceyrəyə düşən
FƏRQLİ-kodlu nümunələri tapır. `resolve_hard_data(on_conflict=...)`:

    "raise" (defolt)  — HƏR ziddiyyətdə `HardDataConflictError`
    "majority"        — ən çox səs; BƏRABƏRLİKDƏ yenə atılır
    "keep_first/last" — sənədləşdirilmiş, deterministik seçim

**TAPILAN VƏ DÜZƏLDİLƏN HƏQİQİ SƏHV**: ilkin inteqrasiyada quyular öz
XAM koordinatları ilə kondisioner kimi verilirdi — SIS-in sərt-data
hörməti isə HƏDƏF massivi ilə DƏQİQ koordinat üst-üstə düşməsinə
əsaslanır (bax `facies._find_hard_data_matches`). Real quyu demək olar
HEÇ VAXT hüceyrə mərkəzində olmadığı üçün bu, sərt datanın SƏSSİZCƏ
HÖRMƏT EDİLMƏMƏSİNƏ səbəb olurdu (`test_hard_data_preserved_through_
full_pipeline` bunu TUTDU). Düzəliş: hər sərt nöqtə öz EV HÜCEYRƏSİNİN
mərkəzinə "sancılır" (snap) — indi hədəf massivi ilə DƏQİQ üst-üstə
düşür, sərt data HƏQİQƏTƏN qorunur.

## 6. Ziddiyyət idarəetməsi

Yuxarı bax — HEÇ VAXT sükutla "sonuncu qazanır". `resolve_hard_data`
strategiyası açıq seçilir, defolt HƏMİŞƏ sərt xəta (`raise`).

## 7. Məkan axtarışı bünövrəsi

`spatial_search.py`: `scipy.spatial.cKDTree` əsaslı. `AnisotropyParams.
transform()` İLƏ EYNİ fəzada indeksləyir (Phase 2-3-dən BİRBAŞA istifadə,
YENİDƏN YAZILMAYIB) — kriging-in özünün işlətdiyi anizotrop məsafə ilə
UYĞUN qonşuluq seçilir (tapşırıq §5 "IMPORTANT").

`IncrementalAnisotropicSearch` — ARDICIL (SIS-in bir-bir nöqtə əlavə
etdiyi) ssenari üçün: cKDTree DƏYİŞMƏZDİR, ona görə `rebuild_interval`
addımda bir YENİDƏN qurulur, aralıqda əlavə olunanlar KİÇİK bufer kimi
brute-force axtarılır və BİRLƏŞDİRİLİR — nəticə TAM brute-force ilə
EYNİDİR (approksimasiya DEYİL).

**Paritetin sübutu** (tapşırıq §6, "Only then replace"):
`tests/test_spatial_search.py` (10 test) hər ssenarini (izotrop/anizotrop,
radius/max_neighbors, ARDICIL əlavə) brute-force ilə MÜQAYİSƏ edir —
HAMISI dəqiq (element-be-element) uyğun gəlir. `tests/test_facies_
integration.py::test_fast_and_brute_force_search_agree_*` bunu TAM SIS
səviyyəsində TƏSDİQLƏYİR (eyni seed → eyni realizasiya, `use_fast_
search=True/False`).

## 8. Anizotropluq uyğunluğu

`simulate_sis` hər fasiya üçün `OrdinaryKriging`-in ÖZÜNÜN defolt
qaydası ilə `AnisotropyParams` qurur (`range_v`/`range_minor`/
`azimuth_deg` `None` olanda EYNİ defolt dəyərlər, bax `interpolation.py
._parameters`-in məntiqi) — axtarış VƏ kriging EYNİ transformasiya
fəzasında işləyir. Test: `test_fast_and_brute_force_agree_with_explicit_
anisotropy_and_radius`.

## 9. Ehtimal düzəliş statistikası

`FaciesDiagnostics` — AYRICA sayğaclar (heç biri digərinə yığılmır):

    negative_probability_events   — mənfi kriging ehtimalı olan hüceyrə sayı
    excess_probability_events     — kliplənmədən əvvəl cəm > 1 olan hüceyrə sayı
    nan_fallback_cells            — NaN səbəbindən qlobal nisbətə keçən hüceyrə
    zero_sum_fallback_cells       — sıfır-cəm səbəbindən qlobal nisbətə keçən hüceyrə

Hər biri üçün say + faiz (`rate()`), `correction_warn_threshold`-dan
(defolt 0.30) çox olanda GÜCLÜ XƏBƏRDARLIQ. Sınanıb:
`test_probability_correction_diagnostics_are_tracked_and_reported`
(dar axtarış radiusu ilə qəsdən yüksək NaN-geri-dönüş dərəcəsi yaradıb
xəbərdarlığın işə düşdüyünü təsdiqləyir).

## 10. Performans bənçmarkı (FAKTİKİ ölçmə, bu mühitdə — UYDURULMAMIŞ)

| Grid | Hüceyrə | sürətli (cKDTree) | brute-force | qeyd |
|---|---|---|---|---|
| 15×15 | 225 | 0.125 s | 0.091 s | brute-force hələ SÜRƏTLİ |
| 50×50 | 2500 | 1.33 s | 1.06 s | brute-force HƏLƏ DƏ rəqabətlidir |
| 100×100 | 10000 | 5.20 s | 7.41 s | sürətli yol İNDİ üstündür (~1.4×) |
| 20×20×5 (3D) | 2000 | 0.99 s | (ölçülmədi) | mötədil 3D hal, tamamlanır |

**DÜRÜST TAPINTI**: kiçik/orta şəbəkələrdə (bir neçə min hüceyrəyə qədər)
cKDTree-nin çağırış-başına sabit xərci (Python/C sərhədi + dövri
yenidən-qurma) numpy-nin vektorlaşdırılmış brute-force məsafə
hesablamasından BÖYÜK ola bilər — sürət ÜSTÜNLÜYÜ YALNIZ ~10000
hüceyrədən başlayaraq özünü göstərir (bu mühitdə ölçülüb). `rebuild_
interval` (64→1000) ilə tənzimləmə CÜZİ fərq yaratdı — sabit xərc əsasən
KDTree sorğusunun ÖZÜNDƏDİR, yenidənqurma tezliyində DEYİL. Bu, "əvvəlcə
düzgünlük, sonra performans" prinsipinə uyğun DÜRÜST hesabatdır —
"sürətli axtarış HƏR ZAMAN daha sürətlidir" İDDİASI EDİLMİR.

Mürəkkəblik: brute-force `O(K·n²)`; sürətli yol axtarış hissəsini
təxminən `O(K·n log n)`-ə endirir, AMMA sabit çağırış xərci `K·n`
sorğunun HƏR BİRİNƏ tətbiq olunur — n kiçik olanda bu sabit xərc üstünlük
təşkil edir. Böyük `n`-də (~10⁴+) asimptotik üstünlük qalib gəlir.

## 11. Qalan iş / NOT DONE

- `FaciesPanel` (UI) HƏLƏ `main_window.py`-nin saxla/yüklə (layihə
  serializasiyası) iş axınına BAĞLANMAYIB — dəyərləri oxumaq mümkündür,
  amma layihə faylında SAXLANMIR. **NOT DONE**.
- `WellPanel`/`GeologyPanel` (quyu CƏDVƏLLƏRİ) hələ fasiya sütunu
  seçimi TƏQDİM ETMİR — YALNIZ CSV idxalı ilə fasiya datası daxil edilə
  bilər. **NOT DONE**.
- Çoxlu realizasiyanın (`n_realizations`) `build()`-ə TAM inteqrasiyası
  yoxdur — `FaciesPanel.realization_count()` mövcuddur, amma `build()`-i
  N dəfə çağırıb N `GeologicalModel` yaratmaq HƏLƏ AVTOMATLAŞDIRILMAYIB
  (sadəcə N ayrı `build()` çağırışı, artıq mümkündür, sarğı YAZILMAYIB).
- Performans optimallaşdırması bu fazada YALNIZ SIS üçün edildi —
  `OrdinaryKriging`-in ÖZÜ (kontinual Kriging/IDW yolu) TOXUNULMADI.
- 100000+ hüceyrəli TAM SAHƏ grid performansı ÖLÇÜLMƏDİ (yalnız 10⁴-ə
  qədər) — **NOT VALIDATED YET**.
