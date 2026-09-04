# LAY-MƏLUMATLI (LAYER-AWARE) MODELLƏŞDİRMƏ

**Məqsəd:** quyu məlumatının hansı grid layında HƏQİQƏTƏN mövcud olduğunu
müəyyən etmək, YALNIZ məlumatla dəstəklənən laylarda şərtli interpolyasiya
aparmaq və məlumatsız layları SƏSSİZCƏ "interpolyasiya olunmuş" kimi
göstərməmək.

---

## 1. KÖK SƏBƏB (nə üçün bu iş lazım idi)

İki MÜSTƏQİL nöqsan var idi:

**(A) UI cədvəlində lay anlayışı ümumiyyətlə yox idi.**
`GeologicalWell`-də bir quyu üçün CƏMİ BİR `porosity`, BİR `permeability`
var idi. `geology_adapter.wells_to_dataset()` `WellSample.layer`-i HEÇ VAXT
doldurmurdu, ona görə `WellDataset.is_layered()` HƏMİŞƏ `False` idi və
`geology_service._interpolate_volume()` hər K üçün `layer=None` sorğusu
göndərirdi — TƏK dəyər BÜTÜN K-lara sinxron yayılırdı.

**(B) "Məlumat mövcudluğu" anlayışı SİSTEMDƏ YOX İDİ.**
Kod yalnız "bu layda nümunə varmı" sualını `_interpolate_volume` daxilində
verirdi və cavab "yox" olanda ya xəta atırdı, ya da (`allow_cross_layer_
fallback=True` ilə) bütün layların nöqtələrini hovuzlayırdı. Xassə-üzrə
("PORO L4-də var, PERMX yoxdur") ayrım MÜMKÜN DEYİLDİ.

**Nəticə:** `top=2000, bottom=2210` faktiki olaraq "L1–L5 üçün data var"
kimi işləyirdi, halbuki bu, sadəcə HƏNDƏSƏDİR.

---

## 2. ALTI ANLAYIŞ — QARIŞDIRILMASI QADAĞANDIR

| # | Anlayış | Harada yaşayır |
|---|---------|----------------|
| 1 | Quyunun FİZİKİ intervalı (`top`/`bottom`, m) | `domain/geology.GeologicalWell` |
| 2 | Grid təbəqə həndəsəsi (K sərhədləri, DZ) | `domain/geometry.CellGeometry`, `layer_edges()`, `interval_layers()` |
| 3 | MƏLUMAT MÖVCUDLUĞU (hansı layda ölçmə var) | `domain/data_availability.PropertyAvailability` |
| 4 | İNTERPOLYASİYA HƏDƏFİ (istifadəçi seçimi) | `geology_service.LayerInterpolationConfig.target_layers` |
| 5 | TAMAMLAMA (completion) üsulu | `geology_service.CompletionSpec` |
| 6 | QEYRİ-MÜƏYYƏNLİK / ETİBARLILIQ | `domain/properties.PropertyProvenance.confidence` |

**(1) ≠ (3).** Quyu beş layı kəsə bilər, amma ölçmə üç layda ola bilər.
**(3) ≠ (4).** Məlumat üç layda ola bilər, istifadəçi isə yalnız ikisini
interpolyasiya etmək istəyə bilər.

---

## 3. YENİ MƏLUMAT AXINI

```
    QUYU CƏDVƏLİ (top/bottom + "Data layları" sütunu)
        ↓  geology_adapter.wells_to_dataset(geometry, policy)
    WellDataset  — HƏR LAY ÜÇÜN AYRICA WellSample (layer=k, depth=hüceyrə mərkəzi)
        ↓  geology/layer_availability.compute_availability()
    ModelDataAvailability  — XASSƏ-SPESİFİK: PORO L1-L5, PERMX L1-L3 …
        ↓  LayerInterpolationConfig.targets_for()
    İNTERPOLYASİYA MASKASI (yalnız seçilmiş VƏ məlumatlı laylar)
        ↓  mövcud Kriging/SGS mühərriki (DƏYİŞMƏYİB)
    interpolated_field  → status INTERPOLATED / sərt-data hüceyrəsi MEASURED
        ↓  CompletionSpec (AÇIQ seçim)
    estimated_field     → ESTIMATED / EXTRAPOLATED / SIMULATED / PRESERVED
        ↓
    PropertyProvenance (original | interpolated | estimated | final +
                        status + method + confidence)
        ↓  GeologicalModel.completeness_issues()
    VALİDASİYA QAPISI  — MISSING qalıbsa ReservoirModelBuilder RƏDD EDİR
        ↓
    3D GEOLOJİ MODEL → REZERVUAR MODELİ → MPFA / TPFA (TOXUNULMAYIB)
```

---

## 4. STATUSLAR

`domain/data_availability.DataStatus`:

| Status | Mənası |
|--------|--------|
| `MEASURED` | quyuda ölçülüb (sərt-data hüceyrəsi) |
| `INTERPOLATED` | həmin layın ÖZ sərt datası ilə hesablanıb |
| `ESTIMATED` | məlumatsız lay, AÇIQ completion üsulu ilə |
| `EXTRAPOLATED` | məlumat zərfindən KƏNARDA qiymətləndirilib |
| `SIMULATED` | stoxastik realizasiya (SGS/SIS) |
| `PRESERVED` | orijinal sahədən toxunulmadan gətirilib |
| `MISSING` | məlumat yoxdur və tamamlanmayıb |

**QƏTİ QAYDA:** `ESTIMATED`/`SIMULATED`/`EXTRAPOLATED` nəticə HEÇ VAXT
`MEASURED` kimi qeyd edilmir.

---

## 5. MƏLUMAT MÖVCUDLUĞU SİYASƏTİ (`LayerDataPolicy`)

Lay etiketi OLMAYAN nümunə (`layer is None`, `depth is None`) nə deməkdir:

- **`BROADCAST`** (defolt, geriyə-uyğun) — bütün laylara aiddir. Köhnə
  davranış; `layer_config` verilmədikdə YEGANƏ rejimdir.
- **`STRICT`** (lay-məlumatlı rejimin defoltu) — heç bir laya aid deyil.
  Mövcudluq YALNIZ açıq `layer`/`depth`/"Data layları" bəyanından gəlir.
- **`INTERVAL`** — quyunun `top`/`bottom` intervalının kəsdiyi laylara
  aiddir. **AÇIQ istifadəçi seçimidir** və hesabatda "FƏRZİYYƏDİR"
  xəbərdarlığı ilə görünür.

İstisna: `WellSample.areal=True` (struktur səthlər — `TOP`/`BOTTOM`).
Bunlar TƏBİƏTƏN laydan asılı deyil, ona görə STRICT rejimdə də hər layda
keçərlidir. Bu, xassə ADINA görə deyil, ÖLÇMƏNİN TƏBİƏTİNƏ görə ayrımdır.

---

## 6. "DATA LAYLARI" SÜTUNU (UI)

Geologiya cədvəlində iki YENİ sütun var:

| Sütun | Növ | Mənası |
|-------|-----|--------|
| **Data layları** | redaktə olunan | hansı layda ÖLÇMƏ var |
| **Kəsdiyi laylar** | yalnız oxunan | `top`/`bottom` + grid həndəsəsinin NƏTİCƏSİ |

Format:

```
1-3                      → bütün xassələr üçün L1–L3
1,2,5                    → L1, L2, L5
PORO:1-5; PERMX:1-3      → XASSƏ-SPESİFİK (§4/§15)
*                        → bütün laylar
(boş)                    → bəyan edilməyib
```

Səhv mətn (`3-1`, `1-9` NZ=5 olanda, hərf) **AÇIQ xəta** verir — səssiz
düzəliş yoxdur.

---

## 7. TAMAMLAMA (COMPLETION) STRATEGİYALARI

`geology_service.CompletionMethod`:

| Üsul | Nəticə statusu | Nə edir |
|------|----------------|---------|
| `NONE` (defolt) | `MISSING` / `PRESERVED` | heç nə — lay toxunulmur |
| `PRESERVE_ORIGINAL` | `PRESERVED` | mövcud geoloji prior saxlanılır |
| `VERTICAL_TREND` | `ESTIMATED` / `EXTRAPOLATED` | məlumatlı layların lay-ortalarından dərinliyə görə xətti trend; laya YALNIZ ORTA verilir |
| `GEOSTATISTICAL_3D` | `ESTIMATED` / `EXTRAPOLATED` | bütün layların sərt datası ilə HƏQİQİ 3D (X,Y,Z) Kriging |
| `SGS` | `SIMULATED` | mövcud SGS mühərriki, YALNIZ məlumatsız layların hüceyrələri hədəflənir |
| `CONSTANT` | `ESTIMATED` | istifadəçinin AÇIQ verdiyi lay dəyəri |

**`VERTICAL_TREND` niyə lay daxilində SABİTDİR?** Şaquli trend YALNIZ
şaquli məlumat verir. Qonşu layın LATERAL xəritəsini kopyalayıb trendlə
miqyaslamaq "L4 = L3-ün surəti" ilə eyni elmi problemi yaradardı — ona görə
lateral struktur UYDURULMUR və bu, hesabatda açıq deyilir.

---

## 8. ETİBARLILIQ (CONFIDENCE)

`PropertyProvenance.confidence` — `[0,1]`, `confidence_kind =
"ordinal_support_score"`. **KALİBRLƏNMİŞ EHTİMAL DEYİL.**

| Hal | Mənbə |
|-----|-------|
| `MEASURED` hüceyrə | `1.0` |
| `INTERPOLATED` | mövcud `PropertyEstimate.confidence` (qonşu sayı + məsafə + nisbi kriginq variansı) → `high 0.90 / medium 0.60 / low 0.30 / extrapolated 0.10` |
| `VERTICAL_TREND` | `exp(-Δz / range_v)`, tavan `0.50` (zərf içi) / `0.30` (kənar) |
| `GEOSTATISTICAL_3D` | kriginq `Confidence` + şaquli məsafə cəzası |
| `SGS` | **`NaN`** — tək realizasiyadan etibarlılıq çıxarmaq olmaz |
| `CONSTANT` / `PRESERVED` | `NaN` (istifadəçi açıq verməyibsə) |
| `MISSING` | `NaN` |

**`range_v` (şaquli korrelyasiya radiusu) verilməyibsə etibarlılıq
HESABLANMIR (`NaN`)** — saxta rəqəm yaradılmır (tapşırıq §18).

---

## 9. VALİDASİYA QAPISI (simulyator inteqrasiyası)

`GeologicalModel.completeness_issues()` MISSING qalan hər məcburi xassə
üçün lay-adlı xəta mətni qaytarır. Bu, `GeologicalModel.validate()`-ə
daxildir, deməli:

- `ReservoirModelBuilder.build()` belə modeli **RƏDD EDİR** →
  simulyasiya işə düşmür.
- `build(layer_config=...)` isə modeli QAYTARIR (3D-də MISSING laylar
  GÖRÜNSÜN deyə), amma `report.blocking` doludur və UI xəbərdarlıq göstərir.

MISSING hüceyrə `NaN` daşıyır. Bu, "problemi NaN-ın arxasında gizlətmək"
DEYİL: status AÇIQ `MISSING`-dir, validasiya bloklayır, hesabat sayır.
`GeologicalModel._defined_values()` fiziki diapazon yoxlamasını yalnız
TƏYİN OLUNMUŞ hüceyrələrdə aparır — beləliklə MISSING-dən BAŞQA səbəbdən
yaranan NaN (əsl korlanma) HƏLƏ DƏ tutulur.

---

## 10. TƏSİR (IMPACT) ANALİZİ

```python
from imex2d.application.geology_service import compute_property_impact

impact = compute_property_impact(original_model, hypothetical_model, "PORO")
impact.delta            # hypothetical − original
impact.relative         # delta / |original|
impact.layer_mean_delta()
```

`ImpactResult` massivləri **KOPYALANIR** — nəticə üzərində aparılan heç bir
əməliyyat mənbə modelinə qayıtmır (tapşırıq §12).

---

## 11. 3D GÖRÜNTÜ

`rendering/renderers.py` DÖRD yeni açar prefiksi tanıyır:

```
"PROVENANCE:PORO"    → status kodu (kateqorik rəng şkalası)
"CONFIDENCE:PERMX"   → etibarlılıq balı [0,1]
"ORIGINAL:PORO"      → ORİJİNAL sahə (varsa)
"IMPACT:PORO"        → final − orijinal (simmetrik ± şkala)
```

Xassə adı AÇARIN İÇİNDƏDİR — rendering qatında heç bir xassə adı sabit
kodlanmayıb. `main_window._refresh_provenance_choices()` siyahını
`model.provenance`-dan qurur (ORIGINAL/IMPACT yalnız orijinal sahə
mövcud olanda təklif olunur). Rendering qatı HEÇ BİR geologiya hesablaması
aparmır, yalnız hazır massivləri göstərir; `IMPACT` fərqi HƏR ÇAĞIRIŞDA
yenidən hesablanır və `final` sahəni DƏYİŞMİR (§12).

**Status filtri.** 3D nəzarət sətrində iki combo var: hansı xassənin
mənşəyinə görə süzüləcəyi və hansı statusun göstəriləcəyi (`Hamısı`,
`ÖLÇÜLÜB`, `İNTERPOLYASİYA`, `QİYMƏTLƏNDİRİLİB`, `EKSTRAPOLYASİYA`,
`SİMULYASİYA`, `ORİJİNAL`, `MƏLUMAT YOXDUR`). Maska
`PropertyProvenance.mask()`-dan gəlir və `VolumeFilter.cell_mask` /
`VtkViewSettings.cell_mask` ilə HƏR İKİ motora (matplotlib və VTK) eyni
şəkildə tətbiq olunur. Mövcud `K` diapazon sürgüləri layları ayrıca
açıb-bağlamağa imkan verir (əvvəldən var idi, dəyişməyib).

---

## 12. NƏYƏ TOXUNULMAYIB

- **MPFA-O / TPFA riyazi nüvəsi** — bir sətir də dəyişməyib.
- **Kriging / variogram / anizotropluq / sərt-data honoring** — mövcud
  `interpolate_property_field()` DƏYİŞMƏDƏN çağırılır; yalnız HANSI layın
  HANSI nöqtələrlə hesablanacağı maskası əlavə olunub.
- **SGS mühərriki** (`geology/sgs.py`) — dəyişməyib; `_simulate_continuous_
  sgs_field()`-ə yalnız `layers=` (hədəf məhdudlaşdırması) parametri
  əlavə olunub, defolt `None` ilə davranış eynidir.
- **Kateqorik/fasiya yolu** (SIS / indikator kriginq) — dəyişməyib.
- **Köhnə iş axını** — `layer_config` verilmədikdə `build()` TAM ƏVVƏLKİ
  kimi işləyir (`provenance` boş, `availability is None`).

---

## 13. İSTİFADƏ (proqram səviyyəsində)

```python
from imex2d.application.geology_adapter import wells_to_dataset
from imex2d.application.geology_service import (
    CompletionMethod, CompletionSpec, GeologicalGridSpec,
    LayerInterpolationConfig, WellBasedGeologicalModelBuilder)
from imex2d.geology.layer_availability import LayerDataPolicy

dataset, skipped = wells_to_dataset(wells, "Kriging (adi)", geometry,
                                    LayerDataPolicy.STRICT)

config = LayerInterpolationConfig(
    policy=LayerDataPolicy.STRICT,
    target_layers=[0, 1, 2],                       # 0-ƏSASLI (UI 1-əsaslıdır)
    property_completion={
        "PORO": CompletionSpec(method=CompletionMethod.VERTICAL_TREND),
        "PERMX": CompletionSpec(method=CompletionMethod.NONE)},   # MISSING qalsın
    original_fields={"PORO": prior_poro})

model, report = builder.build(dataset, spec, layer_config=config)
print(report.as_text())            # lay-üzrə vəziyyət + xəbərdarlıq + bloklar
model.availability["PORO"].as_text()
model.provenance["PORO"].status_counts()
```

---

## 14. TESTLƏR

`tests/test_layer_data_availability.py` — 64 test:

- TEST A–O (tapşırıq §22) — hər biri ayrıca adlandırılıb
- Kənar hallar (§23): boş seçim, diapazondan kənar indeks, `top >= bottom`,
  grid-dən kənar interval, tək lay, lay sərhədi, `NZ = 1/2/5/8`, tək layda
  data, bütün laylarda data, data yoxdur, dublikat nümunə, NaN dəyər
- Performans (§24): `monkeypatch` ilə Kriging/SGS çağırış SAYI yoxlanılır
- Geriyə-uyğunluq (§25): `layer_config` verilmədikdə davranışın EYNİ qalması
- Serializasiya: `data_layers_text`-in saxlanılıb-yüklənməsi
- Simulyator inteqrasiyası (§20/§21): tamamlanmış modelin TPFA VƏ MPFA-O
  diskretizasiyalarına çatması, ucdan-uca simulyasiya, natamam modelin
  axın həlledicisinə HEÇ VAXT ÇATMAMASI
- UI (§6/§14): 1-əsaslı → 0-əsaslı lay çevirməsi, diapazon xətası,
  "Data layları" / "Kəsdiyi laylar" sütunlarının ayrılığı
