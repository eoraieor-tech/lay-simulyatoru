# Vahid çevirmə + fiziki yoxlama qatı (Phase 1)

## Audit nəticəsi: mühərrik hansı vahid sistemində işləyir?

Kod bazası **SI DEYİL** — "neft-sənayesi metrik" (Eclipse-in `METRIC`
vahid dəstinə uyğun): təzyiq **bar**, uzunluq **m**, keçiricilik **mD**,
lözlük **cP**, debit **m³/gün**, sıxlıq **kg/m³**, sıxılma **1/bar**.
Bu, `domain/units.py`-dəki `darcy_constant=0.008527`-in dəqiq törəməsidir
(audit ilə əl hesabı üst-üstə düşdü: `q[m³/gün]=0.008527·k[mD]·A[m²]·
Δp[bar]/(μ[cP]·L[m])`). `FIELD` vahid dəsti `units.py`-da TƏRİF EDİLİB,
amma heç yerdə İSTEHLAK EDİLMİR — istifadə edilsəydi (geometriya/PVT/
quyular metrik qalarkən) transmissibillik SƏSSİZCƏ səhv olardı. Bu, bu
fazanın həll etdiyi əsas risk idi.

## Qərar: kanonik daxili sistem SI-ya KEÇİRİLMƏDİ

Tapşırıq "üstünlük SI-dır" deyir, amma "672 mövcud testi pozma" və
"mövcud ədədi davranışı SƏSSİZCƏ dəyişmə" da tələb edir. Bu ikisi
ziddiyyətlidir: Nyuton qalığını (`residual.py`), Jakobianı (`jacobian.py`),
quyu modelini (`standard_well.py`) və TPFA-nı (`discretization.py`) Pa/m³
saniyə/Pa·s-ə köçürmək bütün validasiya edilmiş ədədi nüvəni yenidən
yazmaq deməkdir — bu, "sərt xəta" riskini "üslub təmizliyi" naminə qəbul
etmək olardı. Ona görə **BİLƏRƏKDƏN** qərar verildi:

* Mühərrikin daxili işlək vahidi **DƏYİŞMİR** (bar/m/mD/cP/m³/gün/kg/m³) —
  672 test bunu qoruyur.
* Yeni `domain/unit_conversions.py` HƏQİQİ SI baza vahidlərinə (Pa, m,
  m², Pa·s, m³, m³/san, kg/m³, K) PİVOT edən, elmi cəhətdən sərt,
  mühərrikdən müstəqil çevirmə təbəqəsidir — istənilən dəstəklənən vahid
  cütü arasında (o cümlədən mühərrikin öz vahidlərinə/vahidlərindən)
  çevirmə üçün.

Bu, tapşırığın "İDEAL" tələbindən (SI nüvə) BİLƏRƏKDƏN kənara çıxmadır —
səbəb yuxarıda izah edilib, "backward compatibility" tələbi "preferably
SI" tələbindən üstün tutulub.

## Çevirmə arxitekturası

```
GİRİŞ (istənilən vahid) → to_engine_units() → MÜHƏRRİK (bar/m/mD/cP/m³/gün)
MÜHƏRRİK → from_engine_units() → ÇIXIŞ (istənilən vahid)
```

Dəstəklənir: təzyiq (Pa/kPa/MPa/bar/psi), uzunluq (m/ft), sahə
(m²/ft²/acre), keçiricilik (m²/D/mD), lözlük (Pa·s/cP), həcm
(m³/ft³/bbl/stb/rb), debit (m³/san/m³/gün/bbl/gün/stb/gün), sıxlıq
(kg/m³/lb/ft³), temperatur (K/°C/°F, ofsetli — ayrıca funksiya), sıxılma
(TƏRS miqyaslı — ayrıca funksiya, aşağı bax).

Bütün sabit faktorlar BİR yerdədir (`PRESSURE_TO_PA`, `LENGTH_TO_M`,
və s.) — koda səpələnmiş magic number YOXDUR.

## Sıxılmanın TƏRS miqyaslanması — tapılan (və rədd edilən) "səhv"

Audit `simulation/pvt/correlations.py:76`-nı ("psi⁻¹ sıxılmanı
`* BAR_TO_PSI` ilə bar⁻¹-ə çevirir") şübhəli saydı ("vurma yox, bölmə
olmalıdır" ehtimalı ilə). **Əl hesabı ilə YOXLANILDI VƏ SƏHV
DEYİL**: `c = -(1/V)(dV/dP)` tərifindən, təzyiq vahidi dəyişəndə
sıxılmanın ədədi qiyməti TƏRS mütənasib dəyişir —

    c[1/bar] = c[1/psi] · (Pa-per-bar / Pa-per-psi) = c[1/psi] · 14.5037744

yəni **VURMA DOĞRUDUR** (bar daha "böyük" vahiddir, ona görə eyni fiziki
sıxılma bar⁻¹-də daha böyük ədədlə ifadə olunur — 4.5e-5 1/bar =
3.1e-6 1/psi, tipik süxur sıxılması). Bu, `convert_compressibility()`-nin
öz ilk versiyasında (TƏRS istiqamətdə yazılmışdı) `test_unit_conversions.py`-
dəki çarpaz-yoxlama testi ilə TUTULDU və düzəldildi — məhz bu səbəbdən
tapşırıqdakı "hər çevirmə üçün test yaz" tələbi vacibdir.
`correlations.py`-nin özü TOXUNULMADI (aşağı, "qalan risklər"ə bax).

## Yoxlama qatı (`domain/validation.py`)

Hər funksiya `(errors, warnings)` qaytarır — SƏRT xəta (fiziki qeyri-
mümkün) və XƏBƏRDARLIQ (qeyri-adi, amma mümkün) AYRILIR:

| Funksiya | Sərt xəta | Xəbərdarlıq (rədd edilmir) |
|---|---|---|
| `validate_porosity` | phi<0 və ya >=1 | phi>0.40 və ya <0.02 |
| `validate_saturation` | S∉[0,1] | — |
| `validate_permeability` | k<=0 | k>20000 mD və ya <0.01 mD |
| `validate_viscosity` | mu<=0 | mu>50000 cP və ya <0.05 cP |
| `validate_density` | rho<=0 | [200,1500] kg/m³ kənarında |
| `validate_compressibility` | c<=0 | [1e-6,1e-2] 1/bar kənarında |
| `validate_pressure` | p<=0 | çatlama qradiyentini aşır (`depth_m` verilsə) |
| `validate_thickness` | dz<=0 | — |
| `validate_grid_dimensions` | NX/NY/NZ<1, DX/DY<=0 | — |
| `validate_cell_volumes` | həcm<=0 (dejenerativ) | çox kiçik/böyük (vahid qarışıqlığı əlaməti) |
| `validate_well_rate` | debit<0 | debit=0 və ya >100000 m³/gün |

NaN/sonsuz HƏR yerdə ayrıca, aydın mesajla tutulur (əvvəllər bəzi
yerlərdə `np.diff(...) <= 0` kimi müqayisələr NaN üçün HƏMİŞƏ `False`
qaytardığı üçün SƏSSİZCƏ keçirdi — bax aşağı, "tapılan gizli boşluqlar").

## Harada bağlandı (mövcud `validate()` metodları DƏYİŞDİRİLMƏDİ, YALNIZ genişləndirildi)

* `CellGeometry.validate()` — **YENİ metod** (əvvəllər YOX idi — audit
  təsdiqlədi). `GeologicalModel.validate()`-ə bağlandı.
* `RockProperties.validate()` — `phi>=1` sərt xətası əlavə edildi
  (əvvəlki `phi<=0` yoxlaması qorunub). `validate_warnings()` YENİ.
* `FluidProperties` — **əvvəllər HEÇ BİR `validate()` yox idi**. İndi var
  (lözlük/sıxlıq/sıxılma/FVF sərt xətaları + `validate_warnings()`).
* `PVTTable.validate()` — NaN/sonsuz aydın mesajla əlavə edildi (əvvəllər
  gizli boşluq idi). `validate_warnings()`/`check_query_range()` YENİ.
* `SaturationTable.validate()` — eyni NaN/sonsuz boşluğu bağlandı.
  `check_query_range()` YENİ (ekstrapolyasiya xəbərdarlığı).
* `WellControl.validate()`/`validate_warnings()` — **YENİ** (əvvəllər
  quyu hədəfi heç vaxt yoxlanmırdı).
* `GeologicalModel.validate()` — geometriya + PORO/PERM sərt xətaları
  bağlandı (geriyə uyğunluq: 672+ testin heç biri pozulmadı, çünki
  bu yoxlamalar `geology_service.py`-nin artıq tətbiq etdiyi klip
  (`PORO∈[0.01,0.45]`) daxilindəki dəyərlərdə HEÇ VAXT işə düşmür).
* `GrdeclImporter` — deck `FIELD`/`LAB` bəyan edibsə (RUNSPEC-də) açıq
  xəbərdarlıq (əvvəllər SƏSSİZ idi — oxucu HƏMİŞƏ metrik kimi oxuyur).

## Mərkəzləşdirilən (dublikat) sabitlər

`equilibrium.py` və `residual.py`-da EYNİ `GRAVITY=9.80665`/
`PA_TO_BAR=1.0e-5` iki yerdə təkrarlanırdı. İndi hər ikisi
`unit_conversions.STANDARD_GRAVITY_M_S2`/`PRESSURE_TO_PA["bar"]`-dan
gəlir. **Bit-bə-bit eynidir** (yoxlanılıb: `1.0/1e5 == 1e-5` Python-da
`True`) — Nyuton/rezidual/inisiallaşdırma testlərinin hamısı YENİDƏN
İŞLƏDİLDİ, sıfır fərq.

`correlations.py`-dəki `BAR_TO_PSI=14.5037744` İSƏ **TOXUNULMADI** —
tam-dəqiq sabitlə (14.5037737730...) əvəz etmək 7-ci mənalı rəqəmdə
kiçik fərq yaradardı; bu fayl Pb-ətrafı Nyuton ossilyasiyası üçün
diqqətlə tənzimlənib (bax modulun öz şərhi), ona görə "səssizcə ədədi
davranışı dəyişmə" qaydasına görə TOXUNULMADI.

## Sınaqlar

* `tests/test_unit_conversions.py` — 47 test: konkret dəyərlər (psi↔Pa,
  bar↔Pa, ft↔m, mD↔m², cP↔Pa·s, bbl↔m³ və s.), dəyirmi-səyahət hər
  kəmiyyət üçün, sıxılmanın TƏRS istiqaməti `correlations.py` ilə
  çarpaz yoxlanıb, Darsi axını METRIC/FIELD girişindən EYNİ SI nəticəni
  verdiyi göstərilib.
* `tests/test_validation.py` — 18 test: hər funksiya üçün sərt-xəta vs.
  xəbərdarlıq ayrımı.
* `tests/test_domain_validation_wiring.py` — 22 test: yuxarıdakı
  bağlamaların doğru işlədiyi (CellGeometry/RockProperties/
  FluidProperties/PVTTable/SaturationTable/WellControl/GeologicalModel).
* `tests/test_unit_invariance_integration.py` — 2 test: **DÜRÜST ƏHATƏ**
  aşağı bax.
* `tests/test_eclipse_io.py` — 2 yeni test (`FIELD` bəyanının
  xəbərdarlıq yaratması).

### İnteqrasiya testinin dürüst əhatəsi

Tapşırıq "Case A psi/ft/mD/stb-day, Case B SI — nəticə fiziki eyni
olmalıdır" istəyir. **Mühərrikin giriş boru xətti (`ReservoirModelBuilder`,
`one_dimensional_model`) HƏLƏ vahid-etiketli giriş qəbul ETMİR** — yalnız
bar/m/mD/m³/gün ədədi gözləyir (bu, "qalan iş"dir, aşağı bax). Ona görə
inteqrasiya testi EYNİ ssenarini FIELD kimi ifadə edir, YENİ çevirmə
qatı ilə mühərrik vahidinə çevirir, və nəticənin birbaşa metrik girişlə
EYNİ olduğunu göstərir (OOIP/RF/addım sayı DƏQİQ, zaman seriyası
1e-4 nisbi tolerantlıqla — adaptiv addım seçimi son-bit fərqlərini
yığır, bu, çevirmənin YANLIŞLIĞI deyil, adaptiv Nyutonun həssaslığıdır).

## Qalan risklər (Phase 1-in ilk hissəsi, bax Phase 1b aşağı üçün YENİLƏNMƏ)

Aşağıdakı bəndlər BU bölmənin YAZILDIĞI andakı vəziyyəti əks etdirir —
Phase 1b (aşağı) bunların ƏKSƏRİYYƏTİNİ bağladı. Tarixi qeyd kimi
saxlanılıb, "əvvəl/sonra" fərqi görünsün deyə.

* ~~Giriş boru xətti hələ vahid-etiketli DEYİL~~ → **Phase 1b-də
  bağlandı**: CSV (`NAME[vahid]`), GRDECL (`FIELD`/`LAB` uzunluq
  çevirməsi), UI (seçilmiş sahələr), PVT (`PVTTable.from_values`).
* ~~`PropertyMap.unit` heç yerdə yoxlanmır~~ → **Phase 1b-də bağlandı**:
  `PROPERTY_QUANTITY` reyestri (keçiricilik/təzyiq) ilə `__post_init__`
  yoxlaması.
* ~~GRDECL FIELD avtomatik çevrilmir~~ → **Phase 1b-də QİSMƏN bağlandı**:
  uzunluq (DX/DY/DZ/TOPS/COORD/ZCORN) ÇEVRİLİR; PVT/təzyiq açar sözləri
  bu oxucuda hələ ÜMUMİYYƏTLƏ yoxdur (aşağı, "Phase 1b qalan risklər").
* Temperatur — DƏYİŞMƏDİ (bilərəkdən, bax aşağı).
* ~~`check_query_range` heç bir axına bağlanmayıb~~ → **Phase 1b-də
  bağlandı**: `ReservoirModel.diagnose()` (bax aşağı, "Ekstrapolyasiya
  darvazası").

---

# Phase 1b — Giriş boru xətti inteqrasiyası

Phase 1-in əsas çevirmə/yoxlama KİTABXANASI hazır olandan sonra bu
mərhələ onu FAKTİKİ giriş nöqtələrinə (CSV, GRDECL, UI, PVT/SCAL,
`ReservoirModelBuilder`) bağladı — məqsəd kitabxananın "istifadə
olunmayan utility" olaraq qalmamasıdır.

## Giriş API-si: `value + unit + quantity`

`unit_conversions.Quantity(value, unit, quantity).to_engine()` —
tapşırıqda tələb olunan üçlüyün rahatlıq örtüyü (məcburi deyil, bütün
daxili kod birbaşa `to_engine_units(value, unit, quantity)` çağırır,
ikisi eyni nəticəni verir).

## CSV (`geology/well_data_io.py`)

Sütun adının sonuna `[vahid]` əlavə edilə bilər: `PERMX[D]`, `x[ft]`.
Baza ad (mötərizədən əvvəlki hissə) DƏYİŞMİR — yalnız dəyər mühərrik
vahidinə çevrilir. **Vahid göstərilməyəndə (köhnə format) dəyər HEÇ
DƏYİŞMİR** — YALNIZ keçiricilik (PERMX/PERMY/PERMZ) üçün bu, açıq
xəbərdarlıqla (`WellDataset.warnings`, YENİ sahə) bildirilir, çünki
mD/Darsi qarışıqlığı real bir səhv mənbəyidir; X/Y/dərinlik üçün
xəbərdarlıq YOXDUR ("metr" bu kod bazasında hər yerdə sənədləşdirilmiş
mövcud konvensiyadır, xatırlatmaq səs-küy olardı). Tanınmayan xassə
üçün (`PORO[%]` kimi) vahid ETİBARSIZ sayılmır, sadəcə "bu xassə üçün
çevirmə dəstəklənmir" xəbərdarlığı ilə DƏYİŞMƏDƏN saxlanılır.

## GRDECL (`io/grdecl_import.py`)

RUNSPEC-də `FIELD`/`LAB` bəyan edilibsə, indi:

* **Uzunluq** (DX/DY/DZ/TOPS, həmçinin corner-point COORD/ZCORN) FIELD→ft,
  LAB→cm-dən mühərrik "m"-inə FAKTİKİ ÇEVRİLİR (əvvəllər sükutla METRIC
  sayılırdı — bu, audit-in tapdığı ƏN TƏHLÜKƏLİ boşluq idi).
* **Keçiricilik** (PERMX/PERMY/PERMZ) ÇEVRİLMİR — bu, çevirmənin
  YAZILMAMASI deyil, Eclipse-in ÜÇ vahid sistemində də (METRIC/FIELD/LAB)
  keçiriciliyin HƏMİŞƏ mD olması faktına görə EHTİYAC YOXLUĞUDUR (əl
  ilə yoxlanılıb, standart Eclipse konvensiyası).
* **Təzyiq/PVT açar sözləri bu oxucuda HƏLƏ ÜMUMİYYƏTLƏ YOXDUR** (nə
  METRIC, nə FIELD üçün) — gələcəkdə əlavə ediləndə FIELD/LAB üçün
  AYRICA çevirmə YAZILMALIDIR, bu hazırkı iş bunu ETMİR. **NOT VALIDATED
  YET** hər hansı gələcək təzyiq/PVT GRDECL importu üçün.

Xəbərdarlıq mətni indi NƏYİN çevrildiyini, NƏYİN YOX olduğunu AÇIQ
bildirir — "hər şey METRIC sayılır" ifadəsi TAMAMILƏ silindi.

## `ReservoirModelBuilder`

`rock_compressibility_unit: str = "bar"` YENİ parametr — defolt
mühərrik vahididir, ona görə DƏYİŞDİRİLMƏSƏ davranış ƏVVƏLKİ kimi qalır.
Digər skalyar girişlər (quyu hədəfləri, ilkin şərtlər) artıq özləri
`WellControl`/`InitialConditions` səviyyəsində qurulur (bax UI bölməsi)
— builder-in özü onları YENİDƏN çevirmir (təkrar çevirmənin qarşısı).

## UI (`ui/panels.py`)

Vahid seçicisi YALNIZ 3 sahəyə əlavə edildi (tapşırığın "hər sahəni
doldurma" xəbərdarlığına uyğun): `NumericalPanel.initial_pressure`
(bar/psi/MPa), `RockFluidPanel.permx` (mD/D/m²), `RockFluidPanel.mu_w`/
`mu_o` (cP/Pa·s, paylaşılan seçici). Defolt seçim HƏMİŞƏ mühərrik
vahididir — toxunulmayan panel ƏVVƏLKİ ƏDƏDİ dəyəri verir (funksional
Qt testi ilə yoxlanılıb, bax aşağı).

**Vahid dəyişəndə** spin-box-un DİAPAZONU və GÖSTƏRİLƏN ƏDƏD YENİDƏN
HESABLANIR (`_bind_unit_aware_spins`) — bunsuz diapazon ilk vahidin
miqyasında qalıb yeni vahiddə səssizcə kəsərdi (bu, YAZARKƏN öz
testimlə TUTULAN reqressiya idi — bax aşağı, "tapılan boşluqlar").

**Bağlanmadı (bilərəkdən)**: `WellPanel` (quyu BHP/RATE cədvəli) və
`GeologyPanel` (quyu-əsaslı keçiricilik cədvəli) — hər ikisi
`QTableWidget`-dir, sxemi dəyişmək serializasiya formatına toxunardı,
üç ayrı test faylına (`test_ui_wiring.py`, `test_ui_static.py`,
`test_serialization.py`) risk yaradardı. **NOT DONE** — bu sahələr
mühərrik vahidini (bar/m³-gün) qəbul etməyə davam edir.

## PropertyMap vahid reyestri

`PROPERTY_QUANTITY = {"PERMX": "permeability", ..., "PRESSURE": "pressure"}`
— YALNIZ real qarışıqlıq riski olan adlar. `unit` bu adlardan biri VƏ
qeyri-boşdursa, `known_units(quantity)`-ə uyğun olmalıdır, əks halda
`ValueError` (məs. `PERMX` + `"psi"`). Boş `unit` (defolt) HEÇ VAXT rədd
edilmir. Reyestrdə OLMAYAN ad (PORO, NTG, SW, REGION_ID) üçün `unit`
sərbəst mətn olaraq qalır — DƏYİŞMƏDİ.

**Phase 2 (tam tenzor permeabilite) əlavəsi**: `KXX/KYY/KZZ/KXY/KXZ/KYZ`
də reyestrə əlavə olundu — diaqonal `PERMX/PERMY/PERMZ` İLƏ EYNİ vahid-
etibarlılığı (mD/D/m²) qoruyur. **Məhdudiyyət**: bu, YALNIZ vahid-ETİKETİNİN
etibarlılığını yoxlayır (`PropertyMap.__post_init__`) — heç bir xassə üçün
(diaqonal DA daxil) konstruksiya zamanı AVTOMATİK ƏDƏDİ ÇEVİRMƏ YOXDUR,
tenzor da istisna deyil. `PermeabilityTensor.convert_units(from, to)`
BÜTÜN 6 komponentə EYNİ amili tətbiq edən hazır aləti təmin edir, AMMA
CSV/GRDECL idxal sərhədində bunu AVTOMATİK çağıran boru xətti HƏLƏ YOXDUR
— çünki tenzor K-nı fayldan oxuyan mexanizmin ÖZÜ hələ yoxdur (bax
`GEOSTATISTICS.md`/`ARCHITECTURE.md` §5.13, audit §13: "tensor components
currently need to be supplied explicitly").

## PVT (`domain/pvt.py`)

`PVTTable.from_values(..., pressure_unit="bar", viscosity_unit="cP",
solution_gor_unit="sm3/sm3", rock_compressibility_unit="bar")` — idxal
sərhədi. **Bo/Bw ÇEVRİLMİR** — bunlar ölçüsüz nisbətdir (rezervuar
həcmi/səth həcmi, hər ikisi eyni fiziki vahiddə ölçülürsə nisbət vahiddən
asılı deyil). Rs (`solution_gor`) `scf/stb`↔`sm3/sm3` arasında çevrilir
(`correlations.py`-dəki `SM3M3_TO_SCFSTB=5.61458` ilə EYNİ faktor, YENİ
`RS_TO_SM3_SM3` cədvəlində mərkəzləşdirilib — `correlations.py`-nin
özünə TOXUNULMADI). Defolt vahidlər mühərrik vahidləridir, ona görə
`from_values(...)` arqumentsiz-vahidlə birbaşa `PVTTable(...)` ilə
ƏDƏD-ƏDƏD eynidir.

`build_pvt_table(..., temperature_unit="C")` — YENİ, defolt "C"
(dəyişməzlik). Korrelyasiyaların ÖZÜ (`correlations.py`) TOXUNULMADI —
yalnız girişdəki temperatur `C`-yə çevrilir, sonra əvvəlki kimi işlədilir.

## SCAL

`SaturationTable.check_query_range()` artıq `ReservoirModel.diagnose()`-ə
bağlıdır (aşağı). Pc mühərrikin təzyiq vahidindədir (bar) — Sw/kr kimi
[0,1] ölçüsüz nisbətlər ÇEVRİLMİR.

## Ekstrapolyasiya darvazası — model qurularkən BİR DƏFƏ, Nyuton həlqəsində DEYİL

`ReservoirModel.diagnose()`-ə YENİ `_check_pvt_scal_ranges()` addımı
əlavə edildi. Yoxlanılan sorğu nöqtələri: ilkin təzyiq + bütün BHP quyu
hədəfləri (PVT cədvəli üçün), ilkin Sw (hər SCAL region cədvəli üçün).
YENİ `validate_query_range()` (bax `domain/validation.py`) DƏRƏCƏLİ
qərar verir:

    diapazon daxilində              -> problem yoxdur
    yüngül kənar (< 0.5×diapazon)    -> XƏBƏRDARLIQ (sərhədə kəsilir, davam edir)
    HƏDDİNDƏN ARTIQ kənar (>= 0.5×)  -> SƏRT XƏTA (vahid qarışıqlığı əlaməti sayılır)

Bu, `model.validate()`-i (yalnız SƏRT xətaları çıxaran mövcud metod)
YALNIZ HƏDDİNDƏN ARTIQ kənarlaşmada bloklayır — mövcud 763 test heç
birində PVT/SCAL diapazonu bu qədər kənar deyil (tam suitin YENİDƏN
işlədilməsi ilə TƏSDİQLƏNİB), ona görə geriyə uyğunluq qorunur.

## Tapılan boşluqlar (özüm-öz testimlə tutulan reqressiyalar)

* İlk `convert_compressibility()` versiyası TƏRS istiqamətdə yazılmışdı
  (bax əvvəlki Phase 1 bölməsi, indi ORADA YOX, BURADA təkrarlanmasın
  deyə silindi — əslində elə əvvəlki bölmədə izah olunub). Bu fazada
  TƏKRAR toxunulmadı, sadəcə `from_values()`-də istifadə edilərkən
  YENİDƏN çarpaz-yoxlanıldı (test #8).
* UI-da vahid seçici spin-box-un DİAPAZONUNU yeniləmirdisə (ilk versiya),
  istifadəçi "psi" seçib 3000 yazanda dəyər SƏSSİZCƏ 1200-ə (~83 bar-a)
  kəsilirdi — məhz YAZDIĞIM `test_ui_units.py` bunu TUTDU, `_bind_unit_
  aware_spins()` ilə düzəldildi.
* `convert()`-in ilk versiyası eyni-vahid (`from==to`) halında belə
  vurma/bölmə edirdi — `1000 mD -> mD` üçün `1000.0000000000001` verirdi.
  Tapşırığın öz test #3-ü ("1000 mD → mühərrik mD") bunu TUTDU; indi
  `from_unit == to_unit` DƏQİQ no-op-dur (bütün `convert*` funksiyalarında).

## Sınaqlar (Phase 1b, YENİ fayllar)

* `tests/test_csv_units.py` — 7 test (bracket sintaksisi, geriyə uyğunluq).
* `tests/test_ui_units.py` — 5 test, **real (offscreen) Qt widget-ləri
  ilə** (digər UI testlərindən fərqli olaraq AST deyil, faktiki widget
  qurulur — bax `conftest.py`).
* `tests/test_pvt_scal_units.py` — 10 test (FIELD/METRIC PVT ekvivalentliyi,
  temperatur vahidi, ekstrapolyasiya darvazası).
* `tests/test_phase1_pipeline_integration.py` — 14 test, tapşırıqda
  sadalanan 14 bəndin HƏR BİRİNİN birbaşa təmsili.
* `tests/test_eclipse_io.py` — 2 YENİ test (FIELD uzunluq FAKTİKİ
  çevrilməsi, keçiriciliyin çevrilməməsi).
* `tests/test_domain_validation_wiring.py` — 5 YENİ test (PropertyMap
  reyestri).
* `tests/test_unit_conversions.py` — 2 YENİ test (`ReservoirModelBuilder`
  sıxılma vahidi, `Quantity` örtüyü).

## Phase 1b qalan risklər (AÇIQ DEYİLİR)

* `WellPanel`/`GeologyPanel` (quyu cədvəlləri) vahid seçicisi YOXDUR —
  yuxarı bax, bilərəkdən bu mərhələdə edilmədi (serializasiya riski).
  **NOT DONE**.
* GRDECL-də təzyiq/PVT açar sözləri (SOLUTION bölməsi və s.) bu oxucuda
  ümumiyyətlə yoxdur — FIELD/LAB üçün bunların çevrilməsi **NOT VALIDATED
  YET** (kodun özü mövcud deyil).
* `PROPERTY_QUANTITY` reyestri yalnız keçiricilik+təzyiqi əhatə edir —
  həcm/sıxlıq kimi digər PropertyMap adları üçün genişləndirilə bilər,
  hazırda YOXDUR.
* Temperatur hələ də simulyasiya STATE-i deyil (bilərəkdən, tapşırıq
  bunu tələb etmirdi) — yalnız PVT cədvəl generatoruna giriş.
* Ekstrapolyasiya darvazasının `severe_factor=0.5` həddi EVRİSTİKADIR
  (diapazonun 50%-i) — real yataq məlumatı ilə tənzimlənməyib, YALNIZ
  daxili testlərlə (sintetik ssenari) yoxlanılıb.
