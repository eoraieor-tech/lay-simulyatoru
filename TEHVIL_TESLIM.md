# Təhvil-təslim — işi başqa kompüterdə davam etdirmək üçün

**Hazırlanıb:** 15 sentyabr 2026 (Seans 31) ·
**Yenilənib:** 16 sentyabr 2026 (Seans 37) · **Son kod commit-i:** G3b (bax `git log`)
**Növbəti iş:** **SPE1CASE2 modeli** (bax §5.4)

> Seans 31-dən bəri BEŞ mərhələ bitib (G6, G4, G1, G2, G3). §5-dəki təsvirlər
> tarixi kontekst kimi saxlanılıb, hər birinin üstündə cari vəziyyət yazılıb.

Bu sənəd işi davam etdirəcək şəxs (insan və ya AI köməkçisi) üçündür.
Tarixçənin tam təfərrüatı `ISH_HESABATI.md`, qərarların səbəbləri `QARARLAR.md`,
SPE1-in mənbəsi və boşluqları isə **`SPE1.md`** faylındadır — burada yalnız
davam etmək üçün LAZIM olanlar var.

Bu fayl 13–15 sentyabrdakı əvvəlki təhvil sənədinin yerinə yazılıb; köhnə
məzmun git tarixçəsindədir (`git show dd425be:TEHVIL_TESLIM.md`).

---

## 1 · Bir baxışda vəziyyət

| | |
|---|---|
| Budaq | `main` = `origin/main` (bu sənədin commit-i) |
| Test dəsti | **2646 keçdi, 1 buraxıldı, 1 xfailed** (~6–9 dəqiqə) |
| Son hesabat bölməsi | **Seans 37** → növbəti yazılacaq: **Seans 38** |
| Son qərar | **Q-30** → növbəti: **Q-31** |
| Aktiv blok | **B7 / SPE1** — boşluqlar bağlandı (G7/G8 istisna), növbəti SPE1CASE2 modeli |

Bitmiş son işlər (yenidən başlamağa ehtiyac yoxdur):

| Seans | İş | Commit |
|---|---|---|
| 27 | **Səssiz səhv:** RATE hədəfi hər perforasiyaya TAM yazılırdı (3 perf → 3× debit) — WI·λ payı ilə bölünür (Q-19) | `1ca084c` |
| 27 | **Səssiz səhv:** Eclipse ixracında RATE debiti GRAT sütununa düşürdü; qaz vurucusu WATER kimi yazılırdı (Q-20) | `fe09c78` |
| 27 | **B7 addım 2:** RATE quyusunda BHP limiti, rejim keçidi, histerezis (Q-21) | `a76c8b2` |
| 28 | **Səssiz səhv:** qaz vuran modeldə sahə GOR = 0 idi; vurulan qaz qrafikə/ixraca əlavə olundu (Q-22) | `769d245` |
| 29 | SPE1 parametrləri OPM deck-indən yoxlandı, hədəf SPE1CASE2, texniki borc siyahısı (Q-23) | `19301bb` |
| 30 | **B7 addım 3 / G5:** SƏTH debiti hədəfi (`ORAT`, qaz `RATE`) (Q-24) | `b9ed645` |
| 32 | **G6:** süxur sıxılmasının istinad təzyiqi ayrıca verilir (`None` → köhnə davranış) (Q-25) | `ef23c77` |
| 33 | **G4:** SGOF cədvəli ilə qaz relperm + **deck oxuyucusu** (Q-26) | `a4a9cf3` |
| 34 | **G1:** deck PVT oxuyucuları PVTW/PVDG/PVTO — itkisiz + açıq yaxınlaşdırma (Q-27) | `3638014` |
| 35 | **G2a:** doymamış neft özlülüyü μo(p, Rs) provider səviyyəsində (Q-28) | `32df7b4` |
| 36 | **G2b:** özlülük mühərriyə (qalıq + Jakobian) qoşuldu; **doyma nöqtəsi lövbəri düzəldildi** (Q-29) | `a07a627` |
| 37 | **G3:** c_o və n hər PVTO qolundan (`oil_branches`); sınaq deck-indəki səhv sətir düzəldildi (Q-30) | `git log` |

**Seans 36-nın ən vacib tapıntısı:** `Bo_sat(Pb)` və `μo_sat(Pb)` lövbərləri
cədvəldən adi interpolyasiya ilə götürülürdü, halbuki Pb adətən düyün deyil və
əyrinin orada sınığı var — μ lövbərində **2.27 % meyl** yaranırdı (B3-B-dən
qalan, indiyə qədər görünməyən qüsur). Düzəliş üç cəhddən sonra tapıldı;
təfərrüat və rəqəmlər `ISH_HESABATI.md` → Seans 36, qərar `QARARLAR.md` → Q-29.

---

## 2 · Başqa kompüterdə işə başlamaq

```bash
git fetch origin
git status                 # yerli dəyişiklik varsa ƏVVƏL onunla məşğul olun
git pull --ff-only origin main
git log --oneline -3       # ən üstdə təhvil commit-i, altında b9ed645
```

Virtual mühit: `run.bat` ardıcıllıqla `.venv`, `..\venv`, `venv` qovluqlarını
axtarır (Q-10, Q-13). Bu maşında **`.venv`, Python 3.12.10** işlədilib (əvvəlki
maşında `venv`, Python 3.14.3 idi — hər ikisində dəst keçib).

```bash
.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider -o addopts=""   # baza: 2626
.venv\Scripts\python.exe app.py                                           # proqram
```

**Başlamazdan əvvəl baza nəticəsini ÖZ maşınınızda alın.** 2626-dan fərqli
çıxarsa, işə başlamadan səbəbini tapın. (Seans 31-də baza 2549 idi; aradakı
fərq Seans 32–36-da yazılan testlərdir.)

`graphify` (CLAUDE.md-də qeyd olunub) bu maşında qurulmayıb, `graphify-out/`
yoxdur — sizdə varsa kod dəyişəndən sonra `graphify update .` işlədin.

---

## 3 · Layihə qaydaları (MÜTLƏQ — `CLAUDE.md`-dən)

1. **Hər iş seansı sənədlənir:** `ISH_HESABATI.md`-ə YENİ bölmə (köhnələr heç
   vaxt dəyişdirilmir), texniki seçim varsa `QARARLAR.md`-ə Q-N, statuslar
   `ROADMAP.md`, `ICRA_PLANI.md`, **`SPE1.md`** (§4 boşluq cədvəli, §5 plan),
   struktur dəyişibsə `ARCHITECTURE.md`.
2. **Dil:** sənədlər Azərbaycan dilində, kod və identifikatorlar ingiliscə.
3. **Uydurma məzmun yazılmır.** Bilinməyən yer ⏳. Rəqəm yazılırsa ÖLÇÜLMÜŞ olmalıdır.
4. **Sənədsiz push edilmir.** Hər mənalı addımdan sonra commit + push.
5. **Nömrə toqquşması:** iki maşın paralel işləyə bilər. Yeni Seans/Q nömrəsi
   yazmazdan əvvəl `git fetch` edib `origin/main`-dəki son nömrəni yoxlayın.
6. **İş axını:** qısa budaq → testlər → sənədlər → `git merge --ff-only` →
   budağı sil → push. Böyük dəyişikliyi bir commit-ə yığmayın.
7. **Qərarı sahibkara verin** nəticəni dəyişən hər seçimdə (Seans 27-də
   perforasiya səhvinin düzəldilməsi belə soruşulub).

---

## 4 · Mühərrikin əsas fəlsəfəsi və quyu qaydalarının yeri

* **Qalıq və Jakobiana mümkün qədər toxunmamaq.** Quyu qaydaları bağlantı
  obyektinin sahələrini (`mode`, `target`, `well_index`, `rate_share`) **Nyuton
  həlləri arasında** dəyişir. Hamısı `simulation/well_constraints.py`-dadır
  (THP isə `simulation/wellbore/thp_control.py`).
* **Yeni `ControlMode` əlavə etməyin** — `mode is ControlMode.BHP` yoxlamaları
  onu səssizcə RATE kimi işlədər (Q-15). Yeni mənanı ayrıca sahə ilə verin
  (`bhp_limit`, `rate_basis` belə edilib).
* **Addım dövrünün sırası** (`implicit/engine.py` və `three_phase_engine.py`
  `run()` — ikisi eyni saxlanılmalıdır):

  ```
  _update_rate_shares()          RATE payı (Q-19)
  surface_rate.predict(state)    səth hədəfi → lay həcmi (Q-24)
  time_stepper.advance(...)
  _thp_outer_loop                THP (Q-15, Q-17)
  _surface_rate_loop             səth debiti düzəlişi (Q-24)
  _bhp_limit_loop                BHP limiti (Q-21)
  record: thp_control / bhp_limit
  ```

  Təkrar həll həmişə `AdaptiveTimeStepper.resolve_step` ilədir (tarixçə şişmir).
* **Hər iddia ölçmə ilə.** Sahibkarın diaqnozu iki dəfə ölçmə ilə dəqiqləşdi
  (Seans 25: THP; Seans 28: "GOR PVT-dən sıfırdır" — əslində hesabat səhvi idi).
* **Səssiz yanlış davranış = ən pis səhv.** Dəstəklənməyən kombinasiya açıq xəta
  (IMPES + THP/BHP limiti/SƏTH bazası; ixracda THP/qaz vurucusu) və ya
  diaqnostika xəbərdarlığı verməlidir.
* **Qalıq/Jakobian dəyişibsə — sonlu fərq testi MÜTLƏQDİR.** Nümunələr:
  `tests/test_implicit_jacobian.py::_max_relative_error`,
  `tests/test_gas_injection.py::_jacobian_error`.

---

## 5 · NÖVBƏTİ İŞLƏR — SPE1CASE2 (sıra `SPE1.md` §5, Q-23 Qərar 4)

Mənbə deck-ləri repoda DEYİL (lisenziya yoxlanılmayıb). Endirmək:

```bash
curl -sSfLO https://raw.githubusercontent.com/OPM/opm-tests/master/spe1/SPE1CASE2.DATA
curl -sSfLO https://raw.githubusercontent.com/OPM/opm-tests/master/spe1/opm-simulation-reference/flow/SPE1CASE2.SMSPEC
curl -sSfLO https://raw.githubusercontent.com/OPM/opm-tests/master/spe1/opm-simulation-reference/flow/SPE1CASE2.UNSMRY
sha256sum SPE1CASE2.DATA   # f3de3d06ab5705381e14e6902a21c1af7bb6b37dbc9d6e4055fe146f50a4c249
python tools/eclipse_summary.py SPE1CASE2 "TIME,FOPR,FGOR,WBHP:PROD,WBHP:INJ,BPR:1,BPR:300" 12
```

Checksum fərqlidirsə OPM faylı dəyişib — `SPE1.md` §2/§3-dəki rəqəmləri yenidən
yoxlayın.

### 5.1 · G6 — süxur sıxılmasının istinad təzyiqi ✅ BİTDİ (Seans 32, Q-25, `ef23c77`)

> Domendə `compressibility_reference_pressure: Optional[float]` — `None` olanda
> datum işlədilir, yəni mövcud modellər dəyişməyib. IMPES-in öz sıxılma hesabı
> yoxlanıldı: orada istinad OXUNMUR ⏳ (backlog-a düşüb).
> Aşağıdakı təsvir tarixi kontekstdir.

* **Problem:** SPE1 `ROCK 14.7 3E-6` — istinad 14.7 psia. Bizdə istinad
  `model.initial_conditions.datum_pressure`-dır: `simulation/implicit/residual.py:71`,
  `simulation/implicit/three_phase_residual.py:182`; Jakobian törəmələri
  `derivatives.py` (`reference_pressure`, `jacobian.py:90`-da ötürülür).
  4800 psia-da məsamə həcmi fərqi `1 + 3e-6·(4800−14.7)` = ~1.4 %.
* **Tövsiyə:** domendə `Optional[float]` istinad sahəsi, `None` → datum (mövcud
  modellər bit-bit eyni). IMPES-in öz sıxılma hesabını da yoxlayın ⏳ (baxılmayıb).
* ⏳ Diqqət: `PVTTable.rock_compressibility` və `model.rock.compressibility`
  AYRI sahələrdir — hansının harada oxunduğunu əvvəl yoxlayın.
* **Qəbul:** `None` ilə tam dəst dəyişməz; istinad verilən modeldə PV testlə;
  akkumulyasiya Jakobianı sonlu fərqlə.

### 5.2 · G4 — SGOF cədvəli ilə qaz relperm ✅ BİTDİ (Seans 33, Q-26, `a4a9cf3`)

> `GasSaturationTable`/`GasSaturationTableSet` (`domain/scal_tables.py`) +
> deck oxuyucusu `io/scal_io.py::read_sgof`. Törəmələr parçalı xəttidir, yəni
> Jakobian qalıqla DƏQİQ uyğundur. Corey yolu bit-bit eyni qalıb.
> Aşağıdakı təsvir tarixi kontekstdir.

* **Problem:** qaz əyrisi yalnız `GasCoreyParameters` (`domain/scal.py`).
  `StoneRelativePermeabilityProvider` (`simulation/stone_relperm.py`) ondan
  `krg(sg, swc)`, `krog(sg, swc, kro_end)` və törəmələri çağırır; Jakobian
  `relperm.gas.krg_derivative(...)` işlədir (`three_phase_residual.py`),
  qaz vurulmasının son nöqtəsi `relperm.gas.krg_end`-dir (`three_phase_newton.py`).
* **Tövsiyə:** eyni interfeysi verən cədvəl sinfi (su-neft üçün
  `domain/scal_tables.py` + `simulation/scal_tables_provider.py` nümunədir —
  parçalı xətti interpolyasiya və interval meylləri).
* SGOF-un son sətri (0.88) OPM-in əlavəsidir ki, `Swc + Sg_maks = 1` olsun (deck şərhi).
* **Qəbul:** törəmələr sonlu fərqlə; Corey yolu bit-bit eyni.

### 5.3 · G1 ✅ + G2 ✅ + G3 ✅ — PVTO/PVDG/PVTW və doymamış Bo/özlülük

> **G1** ✅ Seans 34 (Q-27, `3638014`) — `io/pvt_io.py`: itkisiz `DeckPvt` +
> açıq `to_pvt_table` yaxınlaşdırması (tək təzyiq şəbəkəsinə köçürmə
> GİZLƏDİLMİR). Yeni vahidlər: `Mscf/stb`, `rb/Mscf`.
>
> **G2** ✅ Seans 35 (provider, Q-28) + Seans 36 (mühərrik + lövbər, Q-29) —
> `μo(p,Rs) = μo_sat(Pb(Rs))·(p/Pb)ⁿ`, `n` cədvəlin ÖZ doymamış sətirlərindən
> fit olunur (korrelyasiyadakı 0.278 SPE1 üçün yanlışdır).
>
> **G3** ✅ Seans 37 (Q-30) — `BlackOilPVTProvider(deck.to_pvt_table(),
> oil_branches=deck.oil.branches)`: c_o və n hər qoldan, Rs üzrə
> interpolyasiya. ⚠️ Aşağıdakı "0.3 %" ölçməsi YANLIŞ idi (sınaq
> deck-ində səhv sətir) — həqiqi fərq ~11 %, n 0.460 ↔ 0.580.
> ⏳ Qollar hələ `ReservoirModel`-ə/servisə ötürülmür — SPE1CASE2 modelində.
>
> *(tarixi mətn)* hazırda tək `c_o` var. Seans 34-də ÖLÇÜLÜB: SPE1
> deck-inin Rs = 1.27 və 1.618 qollarının c_o-su bir-birindən **0.3 %**
> fərqlənir. Yəni tək qiymət pis yaxınlaşma deyil, lakin açıq yazılmalıdır.
> Qərar sahibkara verilməlidir: hər Rs üçün öz c_o-su, yoxsa tək qiymət +
> sənəddə ölçülmüş fərq.

* **G1:** `PVTTable` (`domain/pvt.py`) BÜTÜN sütunlar üçün tək təzyiq şəbəkəsidir.
  PVDG 14.7…9014.7, PVTO-nun doymuş qolu isə 14.7…5014.7 psia-dır — ortaq
  şəbəkəyə köçürücü lazımdır. PVTW: Bw(p) sıxılma ilə, μw = 0.318 sabit.
* **G2 (fizika):** üç fazalı mühərrikdə `mu_o = pvt.oil_viscosity(pressure)` —
  yalnız doymuş əyri (`three_phase_newton.py:145`). SPE1 PVTO-da doymamış neftin
  özlülüyü 0.51 → 0.74 cP (Rs = 1.27, 4014.7 → 9014.7 psia). μo(p, Rs) qalığa
  (axın + quyu) və Jakobiana (∂μo/∂Rs doymamış hüceyrədə 3-cü dəyişəndir) girir
  — **iki ayrı commit tövsiyə olunur** (Seans 26 prinsipi), hər biri sonlu fərqlə.
* **G3:** doymamış Bo tək `c_o` ilə (`black_oil.py::_build_undersaturated_branch`,
  cədvəlin `bubble_point`-dən yuxarı hissəsindən fit). SPE1-də Rs = 1.27 və
  1.618 sətirlərinin sıxılması fərqlidir — əvvəl ÖLÇÜN, sonra qərar verin.
* **G7 ⏳:** etalonda vurucu hüceyrəsinin təzyiqi 7534 psia-ya qalxır, doymuş qol
  isə 5014.7-də bitir. `np.interp` sərhəddə saxlayır (Rs_sat = 1.618 plato).
  OPM-in bu haldakı qaydası mənbədən yoxlanılmayıb — uydurmayın.

### 5.4 · SPE1CASE2 modeli və etalonla müqayisə ⏳ ƏSAS QALAN İŞ

* **Vahidlər:** mühərrik METRIC-dir; `domain/unit_conversions.py`
  (`psi_to_bar`, `ft_to_m`, `stb_per_day_to_m3_per_day`, `convert(...)`).
  Qaz (Mscf), Rs (Mscf/STB), Bg (rb/Mscf), sıxlıq (lb/ft³) üçün hazır funksiya
  olub-olmadığını yoxlayın ⏳ — çevirmə əmsallarını əldən yazmazdan əvvəl
  modulla tutuşdurun.
* **Model:** grid, keçiricilik, quyular — `SPE1.md` §2. Perforasiyalar 0-dan
  indeksli: istismarçı `(9, 9, 2)`, vurucu `(0, 0, 0)`; rw = 0.25 ft.
  İstismarçı: `RATE`, `rate_basis=SURFACE`, 20 000 STB/gün, `bhp_limit` 1000 psia.
  Vurucu: `RATE`, qaz, `rate_basis=SURFACE`, 100 000 Mscf/gün, `bhp_limit` 9014 psia.
  İlkin: `use_equilibration`, 4800 psia @ 8400 ft, `solution_gor` = 1.27 Mscf/STB,
  Sw = 0.12, WOC 8450 / GOC 8300 (lay xaricində).
* **PVT (Seans 37):** cədvəli `DeckPvt.to_pvt_table()` ilə `reference_rs`
  OLMADAN qurun və provider-ə `oil_branches=deck.oil.branches` verin.
  `reference_rs = 1.27` Rs-i kəsir (Bo −9 %, μo +17 %), qolsuz defolt isə
  c_o-nu ehtiyat qiymətə salır (Bo −80 %) — ölçülüb. ⏳ Servis provider-i
  yalnız `model.pvt_table`-dan qurur: qolları `ReservoirModel`-ə (və layihə
  faylına) ötürmək bu mərhələnin ilk işidir.
* **Etalon nöqtələri** (`SPE1.md` §3): FOPR 1550-ci gündə 20 000-dən aşağı;
  3650-ci gündə FOPR 5732.65 STB/gün, FGOR 22.1403 Mscf/STB, WBHP PROD 1000 psia.
* **Müqayisə testi:** `tools/eclipse_summary.py`-ı `imex2d/io/`-ya köçürüb testlə
  örtün; reqressiya `tests/golden/` altında (bax `tools/golden.py`).
* ⏳ **Əvvəlcədən bilinən maneə:** `domain/validation.py::validate_well_rate`
  100 000 m³/gündən böyük debitə "qeyri-adi yüksək" xəbərdarlığı verir —
  SPE1 qaz vurması ≈ 2.83·10⁶ sm³/gündür. Yanıldıcı xəbərdarlığı modeldən ƏVVƏL
  həll edin (məs. səth qaz debiti üçün ayrı hədd).

### 5.5 · Sonra (müqayisə nəticəsinə görə)

* G7 (doymuş qolun ekstrapolyasiyası), CASE1 (`DRSDT 0` — Rs artımını qadağan
  edən qayda, qalığa toxunur).

---

## 6 · Texniki borc və açıq qalan backlog

**Texniki borc** (`ROADMAP.md` → «Texniki borc»):

* **TB-1:** üç fazalı RATE istismarçısının quyu Jakobianı — sonlu fərqə qarşı
  xəta 0.4893 (tək perforasiya). Sahibkarın qərarı: indi düzəldilmir.
* **TB-2:** üç fazalı AXIN Jakobianının TƏZYİQ sütunu güclü qarışıq vəziyyətdə
  — xəta 0.1936 (G2b-dən əvvəl 0.5672 idi, yəni YAXŞILAŞIB). Ən pis element
  quyusuz hüceyrənin qaz tənliyindədir, yəni mənbə quyu həddi deyil (Seans 36).
* **TB-3:** `black_oil.py::_saturation_pressure_slope` ən üst Rs düyünündə
  analitik 3.619, sonlu fərq 1.810 — tam 2 dəfə (interpolyasiya düyündən yuxarı
  sabit qalır). SPE1-in G7 boşluğu ilə eyni kökdəndir (Seans 36).

**Açıq qalan ⏳ (təcili deyil):**

* `stone_relperm.py` modul sənədində "Stone II … Eclipse-in defoltu" yazılıb —
  yanlışdır (Eclipse-də `STONE1`/`STONE2` açar sözü ilə seçilir), düzəldilməlidir.
* Başlanğıcdakı debit titrəməsi ("mişar dişi") — sahibkarın modelində görünüb,
  test modelində təkrarlanmayıb. Yoxlamaq üçün `layihe.imx` və ya o qaçışın CSV-si lazımdır.
* `well_control_mode` CSV/JSON ixracına və dashboard-a çıxarılmayıb.
* THP quyusu ilə BHP limitli quyu eyni modeldə: limit dövrəsinin təkrar həllindən
  sonra THP-nin BHP-si yenilənmir (ölçülməyib).
* Vurucu BHP limiti və iki fazalı SƏTH bazalı su vurucusu yalnız vahid səviyyədə sınanıb.
* BHP limitinin tərs düsturu BHP rejimindəki `min(q, 0)` kəsməsini nəzərə almır.
* `ui/main_window.py::export_results` (köhnə sadə CSV) qaz sütunlarını yazmır.
* Köhnə backlog: ilkin tarazlıqda Pcog yoxdur; Eclipse ixracı iki fazalıdır
  (SGOF/PVDG/PVTO yazılmır); cədvəldən Pcog; Vogel IPR; Beggs-Brill, VFPPROD,
  RATE quyusunda THP; üç fazalı mühərrikdə `soft_failure_*` tolerantlıqları.

---

## 7 · Bilinən tələlər (vaxt itirməmək üçün)

* **Windows konsolu cp1254-dür.** Azərbaycan hərfləri çap edən skript
  `UnicodeEncodeError` verir (Seans 27-də sənəd skripti bir `print`-də düşdü).
  Skriptin əvvəlində `sys.stdout = io.TextIOWrapper(sys.stdout.buffer,
  encoding="utf-8", errors="replace")` və ya çap etməyin.
* **Sənədlər CRLF sətir sonludur.** Proqramla dəyişəndə qoruyun (`file X.md`
  "with CRLF line terminators" göstərməlidir, "CRLF, LF" yox). Bu seansdakı sənəd
  skriptləri mətni `\r\n`-ə çevirib bayt kimi yazırdı.
* **Arxa planda gedən tam dəst** test modullarını başlanğıcda import edir: sonradan
  edilən kod dəyişikliyi həmin qaçışa təsir etmir, sonradan yaradılan test faylı
  isə ümumiyyətlə toplanmır. Nəticəni hansı kod vəziyyətinə aid olduğunu bilərək oxuyun.
* **Uzun Python heredoc-u bash-da qırıla bilər** — ayrıca `.py` faylına yazın.
* **PVT cədvəlində doyma təzyiqi DÜYÜN DEYİL.** `np.interp` ilə götürülən
  `Bo_sat(Pb)`/`μo_sat(Pb)` sınığın iki tərəfini qarışdırır. Lövbər həmişə
  doymamış qolun ÖZ fit-indən alınmalıdır (Q-29).
* **Testin keçməsi düzgünlüyün sübutu deyil.** Lövbər testi `atol=1e-5` ilə
  KEÇİRDİ, çünki lövbər və müqayisə sütunu EYNİ interpolyasiya səhvini
  bölüşürdü. İnvariant testi yazanda "hər iki tərəf eyni mənbədəndirmi?"
  sualını verin.
* **Jakobian xətasını sütun NÖVÜ üzrə ölçün** (p / Sw / 3-cü dəyişən). Ümumi
  maksimum təzyiq sütunundakı mövcud borcla (TB-2) üstələnir və yeni həddin
  təsirini gizlədir — nümunə:
  `tests/test_undersaturated_viscosity_engine.py::_column_kind_errors`.
* **`*.imx` artıq `.gitignore`-dadır** (Seans 36) — sahibkarın iş faylı
  təsadüfən repoya düşməsin.
* **`MainWindow()` testdə yaradılsa pytest çökür.** UI-ni mənbə yoxlaması ilə sınayın
  (`inspect.getsource(panels.WellPanel)`).
* **Proqramı arxa plan terminal əmri ilə açsanız seans bağlananda bağlanır** —
  PowerShell: `Start-Process -FilePath ".venv\Scripts\python.exe" -ArgumentList "app.py"`.
* **`layihe.imx` git-də DEYİL** (sahibkarın iş faylı). Seans 25-in 41×41×3 modeli
  odur; saxlanmış halda qaz fazası söndürülüb.
* **Ölçmə skriptləri repoda yoxdur** (müvəqqəti qovluqda idi); rəqəmlər
  `ISH_HESABATI.md`-də ölçmə şərtləri ilə yazılıb — təkrarlamaq üçün oradakı
  model təsvirindən istifadə edin. İstisna: `tools/eclipse_summary.py`.
* **Yalnız RATE quyulu model** servis doğrulamasından keçmir ("ən azı bir quyu
  BHP, THP və ya BHP limiti ilə idarə olunmalıdır") — testdə ya BHP limiti verin,
  ya mühərriki birbaşa qurun.

---

## 8 · Toxunmamalı olan

* **`b2a795e` commit-i** yalnız digər maşındadır və onu **sahibkar özü xilas
  edəcək**. İstənilən birləşmədən əvvəl `git log --all --oneline | grep b2a795e`
  ilə yoxlayın və sahibkarla razılaşdırın — özbaşına birləşdirməyin, silməyin.

---

*Bu fayl təhvil anının şəklidir. İş davam etdikcə cari vəziyyət `ISH_HESABATI.md`,
`ROADMAP.md` və `SPE1.md`-dədir; bu fayl köhnələrsə, həmin sənədlər əsasdır.*
