# İcra planı — sahibkarın hədəf proqramına çatmaq

**Tarix:** 10 sentyabr 2026
**Əsas:** Seans 3 təhlili ([ISH_HESABATI.md](ISH_HESABATI.md)) və
[QARARLAR.md](QARARLAR.md) → Q-08.

Bu sənəd "nə çatmır"dan "necə edirik"ə keçidi göstərir. Hər blok
müstəqil commit-lənə bilən, testlə bitən iş vahididir.

---

## 0. Hədəf tərifi — nə vaxt "bitdi" deyəcəyik

Sahibkarın orijinal tələbindən çıxarılan **qəbul meyarları**:

| # | Meyar | İndi |
|---|---|---|
| M1 | Kəşfiyyat quyusu faylı → 3D heterogen lay modeli, ekranda | ✅ |
| M2 | **MPFA-O** anizotrop tenzorla, istifadəçi tərəfindən seçilə bilən | ✅ **B1-də bağlandı** |
| M3 | 5-spot (1 vurma + 4 hasilat), günbəgün, **3 fazalı** | ✅ **B2-də bağlandı** |
| M4 | P < Psat → qazın ayrılması, GOR artımı | 🟡 **qoşuldu**, amma ilkin Rs sahəsi yoxdur + B3 bloklayır |
| M5 | **THP və BHP** hər quyu üçün, qrafikdə | 🟡 BHP ✅, THP ❌ |
| M6 | RF (%), Water Cut, GOR, orta təzyiq — günbəgün | ✅ |
| M7 | 3D-də cəbhənin hərəkəti + interaktiv kəsik | ❌ |
| M8 | Nəticələrin fayla ixracı | 🟡 PDF ✅, CSV ❌ |

---

## 1. Blok sırası və səbəbi

```
B1 (MPFA seçimi)  ─┐
                   ├─→ B3 (Nyuton möhkəmliyi) ─→ B7 (yekun doğrulama)
B2 (A7 qaytarılması)┘                              ↑
B4 (THP/VFP) ──────────────────────────────────────┤
B5 (Pcog + CSV) ───────────────────────────────────┤
B6 (3D animasiya + slice) ─────────────────────────┘
```

**Niyə bu sıra:**

- **B1 əvvəl** — ən kiçik iş (~1 seans), amma layihənin ƏSAS texniki
  fərqini (MPFA-O) istifadəçiyə açır. Kod hazırdır, yalnız düymə lazımdır.
- **B2 sonra** — M3/M4 bundan asılıdır. Ölçüldü: gözləniləndən **xeyli
  kiçikdir** (aşağı bax).
- **B3 B2-dən sonra** — qaz mühərriki qoşulmayınca rəqsi 41×41-də
  təkrar edib düzəltmək mümkün deyil.
- **B4, B5, B6 paraleldir** — bir-birindən asılı deyil, istənilən
  ardıcıllıqla görülə bilər.
- **B7 axırda** — hamısını bir 5-spot ssenarisində birləşdirən yekun
  doğrulama.

---

## B1 — MPFA-O-nun istifadəçiyə açılması ✅ BİTDİ

**Həcm:** kiçik (~1 seans) · **Risk:** aşağı · **Bitdi:** 10 sentyabr 2026 · 16 yeni test

> **İcra qeydi:** planlaşdırılan 5 addımın hamısı edildi. Planda
> nəzərdə tutulmayan, icra zamanı üzə çıxan bir məsələ:
> `MPFAODiscretization()`-un DEFOLT sərhəd bağlanışı
> (`DIRICHLET`) qalıq qatı tərəfindən hələ dəstəklənmir, ona
> görə servis onu `NEUMANN_ZERO` ilə qurur — bu, simulyatorun
> onsuz da tətbiq etdiyi axınsız sərhəddir, fizika DƏYİŞMİR.
> Təfərrüat: `ISH_HESABATI.md` → Seans 4.

MPFA-O artıq `FullyImplicitEngine`-ə qoşulub (Phase 5B-2, 145 test
keçir), amma `ModelAwareSimulationService.create_engine()`
`flux_discretization` ötürmür — yəni UI-dən seçmək mümkün deyil.

**Addımlar:**

1. `application/config.py` → `SimulationConfig`-ə `flux_scheme: str = "TPFA"`
   sahəsi (`"TPFA"` / `"MPFA-O"`), `validate()`-də dəyər yoxlaması.
2. `application/simulation_service.py` → `create_engine()`-də sxemə görə
   `TwoPointFluxDiscretization()` və ya `MPFAODiscretization(...)` qurulur
   və mühərrikə ötürülür.
3. **IMPES + MPFA-O qadağası** istifadəçi dilinə çevrilir: hazırda
   `_reject_multipoint_impes()` `NotImplementedError` atır — servis
   səviyyəsində tutulub aydın mesaja çevrilməlidir
   ("MPFA-O yalnız tam implicit mühərriklə işləyir").
4. `ui/panels.py` → `NumericalPanel`-ə "Axın diskretizasiyası" combobox-u
   + qısa izah etiketi (anizotropiya / qeyri-ortoqonal grid üçün MPFA-O).
5. `application/serialization.py` → `.imx` faylına açar; **köhnə fayllar
   açılanda defolt `"TPFA"`** (geriyə uyğunluq).

**Bitmə şərti (testlər):**
- `create_engine()` sxemə görə DÜZGÜN sinif qurur (2 test)
- IMPES + MPFA-O → aydın, dildə mesaj (1 test)
- `.imx` gedər-gələr dövrü; açarı olmayan köhnə fayl `"TPFA"` verir (2 test)
- UI: combobox dəyəri konfiqurasiyaya düşür (`test_ui_wiring.py` üslubunda)

---

## B2 — A7 qaz fazasının servis və UI-yə qaytarılması ✅ BİTDİ

**Həcm:** orta · **Risk:** orta · **Bitdi:** 10 sentyabr 2026

> **İcra qeydi — planın bir fərziyyəsi SƏHV idi.** Plan 4b-ni "ən çətin,
> ona görə birinci" sayırdı. Ölçüldü: `standard_well.py` və
> `coupled_newton.py` YALNIZ testlərdən çağırılır — heç bir mühərrik
> onları idxal etmir. 4b heç nəyi bloklamırdı; əsl iş addım 2 və 1 idi
> (~120 semantik sətir). Sıra dəyişdirildi: 2 → 1 → 3 → 4a/4b → 001cc12.
>
> Bərpa zamanı ÜÇ səhv tapılıb düzəldildi (faza sayından asılı kütlə
> balansı, ACTNUM-un dəyişən sayı, `self.reservoir` atributu) — heç biri
> statik köçürmə ilə görünməzdi.
>
> ⚠️ **Qalan iki məhdudiyyət:** (1) GOC verilməyəndə neft "ölü" başlayır
> (Rs=0) — domain-də ilkin Rs sahəsi YOXDUR, ayrıca iş lazımdır;
> (2) Pb > BHP rejimində üç fazalı mühərrik də yığılmır — **B3**.
>
> Təfərrüat: `ISH_HESABATI.md` → Seans 6.

### Ölçülmüş vəziyyət — gözləniləndən KİÇİKDİR

v69 silmə commit-lərinin diff-i nəhəng görünür, çünki patch skriptləri
faylların **sətir sonluğunu (CRLF↔LF) dəyişib**. Sətir sonluğu nəzərə
alınmayanda real ölçü:

| v69 addımı | Görünən | **Real semantik** | Sonrakı toqquşma |
|---|---|---|---|
| 1 — UI (`panels`, `main_window`) | 3 306 | **97 sətir** | ⚠️ `panels.py` sonra 8 dəfə dəyişib |
| 2 — application | 785 | **20 sətir** (+133 test) | ⚠️ 2 commit |
| 3 — rendering (`draw_gas`) | 547 | **20 sətir** | ✅ yoxdur |
| 4a — `well_state.py` | 556 | **16 sətir** | ✅ **fayl toxunulmayıb** |
| 4b — `standard_well` + `coupled_newton` | 443 | **443 sətir** | ⚠️ `coupled_newton` 2 commit |
| `001cc12` — versiya / `.imx` | 34 | **34 sətir** | seçici davranmaq lazımdır |

⚠️ **Düz `git revert` işləməyəcək** — yoxlanıldı, altı commit-in hamısı
toqquşur. Metod: hər addımı
`git diff --ignore-cr-at-eol <sha>^ <sha>` ilə oxuyub **əl ilə tətbiq
etmək**, sonra testlə təsdiqləmək.

### Addımlar (v69-un TƏRS sırası ilə)

1. **4b** — `standard_well.py` + `coupled_newton.py` üç fazalı formaya.
   `standard_well.py` v69-dan sonra **heç dəyişməyib** → `66deca8`-dakı
   versiyası birbaşa götürülə bilər. `coupled_newton.py` üçün əl ilə
   birləşdirmə (2 sonrakı commit-in dəyişikliyi saxlanılmalıdır).
   *Ən çətin addım budur — əvvəl edilir ki, qalanı üstünə otursun.*
2. **4a** — `well_state.py` üç fazalı vəziyyət (16 sətir, toqquşmasız).
3. **3** — `rendering/renderers.py::draw_gas` (20 sətir).
4. **2** — `simulation_service.create_engine()`-də **qaz budağı**:
   PVT-də qaz xassələri varsa `ThreePhaseSimulationEngine` seçilir;
   `model_builder`-də `gas_scal`; `serialization`-da `.imx` açarı;
   `tests/test_gas_ui_wiring.py` `berpa/`-dan `tests/`-ə qaytarılır.

   ⚠️ **Diqqət (A7_PLAN-da qeyd olunmuş, təkrarlanmamalı iki səhv):**
   - üç fazalı mühərrik `ScipyCgIluSolver`-i İŞLƏTMƏMƏLİDİR (2×2 CPR
     üçün tənzimlənib, 3×3 strukturu çökdürür) — öz `NewtonLinearSolver`
     nüsxəsini saxlayır;
   - PVT-də qaz aktiv, SCAL-da qaz əyriləri seçilməyibsə → defolt
     `GasCoreyParameters()` işlədilir (yoxsa Stone relperm gözləyən
     mühərrik yalnız Corey alır və çökür).
5. **1** — UI: PVT tabında "Qaz fazasını aktivləşdir", Ədədi
   parametrlərdə "Qaz papağı (GOC)" + dərinlik girişi, SCAL tabında
   qaz-neft əyrilərinin "(önizləmə)" statusunun götürülməsi (97 sətir).
6. **`001cc12` — SEÇİCİ**: versiya artımı və doyma təzyiqi xəbərdarlığı
   qaytarılır; `.imx` **geriyə uyğunluğu SAXLANILIR** (onu geri qaytarmaq
   köhnə layihə fayllarını sındırardı).

**Bitmə şərti:**
- `tests/test_gas_ui_wiring.py` daxil olmaqla mövcud 143 qaz testi keçir
- UI-dən qaz aktivləşdirilib **10×10 gridd**ə 5-spot işə düşür,
  `converged` və GOR əyrisi görünür
- Qaz söndürülmüş halda 2 fazalı nəticələr **bitə-bit dəyişmir**
  (`test_regression.py` 5-spot etalonu) — ən vacib qoruma budur

---

## B3 — Qaz mühərrikinin Nyuton möhkəmliyi

**Həcm:** orta-böyük (~2-3 seans) · **Risk:** yüksək (yeganə açıq
riyazi problem)

> ## ❗ YENİLƏNDİ (Seans 5) — B3 ARTIQ B2-ni GÖZLƏMİR
>
> Real istifadədə tapıldı ki, tam implicit mühərrik **2 fazalı** halda
> da, qaz fazası olmadan da yığılmır — kifayətdir ki, doyma təzyiqi
> istismarçının BHP-sindən yuxarı olsun (panelin DEFOLT dəyərləri:
> Pb=240 bar, ilkin 250 bar, BHP 150 bar).
>
> **Təkrarlanma 11×11 gridd 0.4 saniyə çəkir** — əvvəl bunun üçün A7
> qaz mühərriki (B2) lazım idi.
>
> | Pb | FIM | IMPES |
> |---|---|---|
> | 240 bar (defolt) | ❌ yığılmır | ✅ işləyir |
> | 140 / 100 / 50 bar | ✅ işləyir | ✅ |
>
> **Nəticə:** B3 indi B1-dən dərhal sonra görülə bilər və hər cəhd
> saniyələrlə ölçülür. Üstəlik bu, istifadəçinin FAKTİKİ qarşılaşdığı
> səhvdir — "PVT modelini işlət" defolt ayarlarla işləmir.
>
> Təfərrüat və ölçmələr: `ISH_HESABATI.md` → Seans 5.

**Məlum diaqnoz** (A7_PLAN, dəqiq ölçülüb): quyu öz BHP həddinə
yaxınlaşanda (P≈8 bar), istismarçının BİLAVASİTƏ QONŞUSU olan hüceyrədə
neft qalığı DÖVR-2 rəqs edir (8·10⁻⁴ ↔ 4.5·10⁻³). 41×41 defolt gridd
problem t≈0 ətrafında üzə çıxır.

**Artıq sınanmış və işə yaramayanlar** (təkrarlanmır): `mb_gas` bölmə
düzəlişi, faza-miqyaslı line search, upstream dondurma.

**Sınanacaq — ədəbiyyatda həlli olan sıra ilə:**

1. **Appleyard chopping** — Nyuton addımında doyma dəyişimini
   məhdudlaşdırmaq (|ΔS| ≤ 0.2), təzyiq dəyişimini də.
   *Ən sadə, ən çox ehtimalla kifayət edən.*
2. **Per-cell qəbul meyarı** — line search bütöv normla deyil, **ən pis
   hüceyrə** ilə qərar verir (CNV üslubu).
3. **Trust-region (Wang & Tchelepi)** — addımı axın funksiyasının
   (fractional flow) əyilmə nöqtəsindən keçməyə qoymamaq.
4. **Quyu rejim keçidində histerezis** — BHP↔RATE arasında "çırpınma"nı
   dayandırmaq.

**Bitmə şərti:**
- **2 fazalı hal (Seans 5):** 11×11, PVT korrelyasiyası, Pb=240 bar,
  FIM → `converged=True`. Reqressiya testi kimi yazılır.
- 5×5 gridd t≈6.7 gündəki qaz rəqs halı **reqressiya testi** kimi
  yazılır və keçir (B2-dən sonra)
- 41×41 defolt gridd 5-spot **`converged=True`** ilə sona çatır
- 2 fazalı nəticələr dəyişmir (yenə etalon testi)

---

## B4 — THP / VFP modulu (SIFIRDAN)

**Həcm:** orta (~2 seans) · **Risk:** aşağı · **Asılılıq: yoxdur**

Kodda THP/VFP **ümumiyyətlə yoxdur** — planın ən dəyərli bəndi (3.8).

**Arxitektura (heksaqonal quruluşa uyğun):**

| Qat | Fayl | Məzmun |
|---|---|---|
| domain | `domain/vfp.py` | `VFPTable` — Eclipse `VFPPROD` analoqu: oxlar (FLO, THP, WFR, GFR, ALQ), çoxxətti interpolyasiya → BHP |
| interfaces | `interfaces/providers.py` | `IWellHydraulicsProvider` müqaviləsi |
| simulation | `simulation/vfp_model.py` | **Analitik model** (cədvəl yoxdursa): hidrostatik sütun + sürtünmə itkisi |
| domain | `domain/wells.py` | `ControlMode.THP` + boru həndəsəsi (diametr, kələ-kötürlük, dərinlik) |
| io | `io/vfp_io.py` | `VFPPROD` cədvəlinin oxunması |

**Analitik model — I versiya (açıq sadələşdirmə ilə):**

```
BHP = THP + ρ_qarışıq·g·h·1e-5  +  ΔP_sürtünmə          [bar]

ρ_qarışıq   — no-slip qarışıq sıxlığı (slip HESABA ALINMIR — sənədləşir)
ΔP_sürtünmə — Darcy-Weisbach; sürtünmə əmsalı: Chen (1979) açıq düsturu
```

Slip və axın rejimi xəritəsi (Hagedorn-Brown, Beggs-Brill) **II
versiyaya** qalır — sənəddə ⏳ ilə açıq göstərilir, uydurulmur.

**Quyu ilə birləşmə:** `ControlMode.THP` seçiləndə addımın əvvəlində
debitlərə görə BHP hesablanır və Nyutona **sabit BHP kimi** verilir
(açıq / explicit birləşmə). Tam implicit THP birləşməsi ⏳ II versiya.

**Bitmə şərti:**
- Sürtünmə sıfır olanda nəticə **analitik hidrostatik sütunla** üst-üstə
  düşür (1 test)
- BHP→THP→BHP gedər-gələr dövrü (1 test)
- Monotonluq: debit artdıqca sürtünmə itkisi artır (1 test)
- Vahid invariantlığı (`test_unit_invariance_integration.py` üslubunda)
- `VFPPROD` cədvəl oxucusu + interpolyasiya (2 test)
- `SimulationResult`-da `well_thp` zaman sırası görünür

---

## B4b — İlkin Rs sahəsi (B2-də aşkarlandı, YENİ)

**Həcm:** kiçik-orta · **Risk:** aşağı · **Meyar:** M4

Domain modelində "ilkin həll olmuş qaz (Rs)" sahəsi YOXDUR. Ona görə
qaz papağı (GOC) verilmədikdə üç fazalı mühərrik nefti "ölü" kimi
başladır (Rs = 0) və **P < Pb olsa belə qaz ayrıla bilmir** — ayrılacaq
həll olmuş qaz yoxdur. OGIP = 0 çıxır.

Bu, v69-un sildiyi bir şey DEYİL — heç vaxt olmayıb (A7-nin açıq
sənədləşdirilmiş mühafizəkar defoltu).

**Addımlar:**
1. `domain/initial.py` → `InitialConditions.solution_gor` (və ya
   "doyma təzyiqindən hesabla" bayrağı).
2. `three_phase_engine._initial_state()` → GOC yoxdursa Rs-i
   `Rs(min(P, Pb))` ilə doldurmaq.
3. UI (PVT və ya Ədədi parametrlər tabı) + `.imx` açarı.
4. Testlər: OGIP > 0; P < Pb-də Sg artır; GOR yüksəlir.

**Bitmə şərti:** GOC OLMADAN da, sadəcə təzyiq Pb-dən aşağı düşəndə
qaz ayrılır və GOR əyrisi qalxır (M4-ün əsl tələbi).

---

## B5 — Kiçik boşluqlar: Pcog və CSV/JSON ixracı

**Həcm:** kiçik (~1 seans) · **Risk:** aşağı

1. **Pcog** — hazırda yalnız `BrooksCoreyCapillaryProvider.pcow` var.
   Qaz-neft kapilyar təzyiqi və `∂Pcog/∂Sg` (Jakobian üçün) əlavə olunur.
   *B2-dən sonra edilməlidir* — üç fazalı qalıq onu istifadə edəcək.
2. **CSV / JSON ixracı** — `TimeSeries`-in bütün sütunları (RF, WCT, GOR,
   orta təzyiq, kumulyativlər) + quyu üzrə debitlər, BHP, THP.
   `reporting/results_export.py`; UI-də "Nəticələri ixrac et" düyməsi.

**Bitmə şərti:** ixrac faylı yenidən oxunanda eyni ədədləri verir;
Pcog üçün monotonluq və son nöqtə testləri.

---

## B6 — 3D canlı vizualizasiya

**Həcm:** orta (~2 seans) · **Risk:** aşağı · **Yeni asılılıq: YOXDUR**

Q-08 təklifinə uyğun: **PyVista-ya keçmirik**, mövcud
`rendering/vtk_volume.py` (895 sətir, 49 test) genişləndirilir.

1. **İnteraktiv kəsik (slice plane)** — `vtkPlaneWidget` + `vtkCutter`,
   `VtkReservoirScene`-ə metod kimi; UI-də "Kəsik" açarı və ox seçimi.
2. **Cəbhə animasiyası** — `SimulationResult.snapshots` üzərində
   `QTimer` ilə oynatma; sürət tənzimi və "geri / irəli" idarəsi.
3. **GIF / PNG ixracı** — kadrlar `vtkWindowToImageFilter` ilə alınır,
   GIF **Pillow ilə** yazılır (Pillow 12.3.0 artıq quraşdırılıb —
   `imageio` / `ffmpeg` LAZIM DEYİL).
4. **Dashboard-a THP əlavəsi** — B4 bitəndən sonra
   `ProductionCurveRenderer`-ə THP/BHP oxu.

**Bitmə şərti:** offscreen rejimdə (test mühiti) kəsik və animasiya
kadrlarının yaradılması test edilir; GIF faylı yaranır və açılır.

---

## B7 — Yekun doğrulama və sənədləşdirmə

**Həcm:** kiçik-orta (~1-2 seans)

1. **Uc-uca ssenari:** 3D corner-point + anizotrop tenzor + **MPFA-O** +
   **üç fazalı** 5-spot + THP nəzarəti → RF / GOR / WCT / THP əyriləri +
   3D animasiya. Reqressiya testi kimi qeyd olunur (`tests/golden/`).
2. **SPE1 benchmark** (backlog-dan qaldırılır) — üç fazalı fizikanın
   sənaye etalonu ilə tutuşdurulması. *"Peşəkar" iddiasının ən güclü
   sübutu budur.*
3. Sənədlər: `ARCHITECTURE.md` (VFP qatı), `ROADMAP.md` statusları,
   `ISH_HESABATI.md`, `QARARLAR.md`.
4. **Test dəstinin sürətləndirilməsi** (audit §9): bütöv dəst 15
   dəqiqədir. `pytest-xdist` ilə paralelləşdirmə gündəlik işi
   rahatlaşdırır.

---

## 2. Nə ETMİRİK (Q-08 təklifi)

Bunlar plandan **qəsdən çıxarılıb**, səbəbi [QARARLAR.md](QARARLAR.md)
→ Q-08-də yazılıb:

- ❌ Python 3.11-ə enmək
- ❌ RBFInterpolator və PyKrige (öz kriginqimiz güclüdür)
- ❌ `spsolve` / SuperLU / AMG (CG+ILU + CPR var)
- ❌ PyVista-ya köçmək (eyni VTK-nın örtüyüdür, 895 sətir yenidən yazılardı)
- ❌ Plotly (masaüstü PyQt5 tətbiqi üçün uyğunsuzdur)
- ❌ Plandakı düz `src/` strukturuna keçmək (heksaqonal arxitektura geriyə getməz)

---

## 3. Ümumi mənzərə

| Blok | Həcm | Risk | Hansı meyarı bağlayır |
|---|---|---|---|
| ~~B1 MPFA seçimi~~ ✅ **BİTDİ** | kiçik | aşağı | **M2 bağlandı** |
| ~~B2 A7 qaytarılması~~ ✅ **BİTDİ** | orta | orta | **M3 ✅ · M4 qismən** |
| B3 Nyuton möhkəmliyi | orta-böyük | **yüksək** | M3 + **PVT səhvi (Seans 5)** |
| B4 THP/VFP | orta | aşağı | **M5** |
| B5 Pcog + CSV | kiçik | aşağı | **M8** |
| B6 3D animasiya + slice | orta | aşağı | **M7** |
| B7 yekun doğrulama | kiçik-orta | aşağı | hamısı |

**Yeganə yüksək riskli iş B3-dür** — riyazi problemdir, cəhd tələb edir.
Qalan hər şey mühəndislik işidir: yolu məlum, ölçüsü bilinir.

**Hər blokda dəyişməyən qayda:** 2 fazalı 5-spot etalon nəticəsi
(`test_regression.py`) pozulmamalıdır. Pozulursa — dəyişiklik səhvdir,
etalon deyil.

⏳ **Sahibkarın təsdiqi gözlənilir:** blok sırası və Q-08 texnologiya
qərarları.
