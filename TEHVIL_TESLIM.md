# Təhvil-təslim — işi başqa kompüterdə davam etdirmək üçün

**Hazırlanıb:** 16 sentyabr 2026 (Seans 38) · **Yenilənib:** 20 sentyabr 2026 (Seans 44)
**Son kod commit-i:** `647ac78` (`imex2d/`; Seans 44 yalnız `tools/` və sənəd əlavə etdi) · son sənəd commit-i bu faylın öz commit-idir
(`git log --oneline -3` ilə yoxlayın)
**Növbəti iş:** SPE1CASE2-nin QALAN fərqi üçün yeni namizəd tapmaq (§5) —
əsas səbəb (G7) tapılıb və bağlanıb, üç namizəd isə ölçülərək istisna olunub.

Bu sənəd işi öz kompüterində davam etdirəcək şəxs (insan və ya AI köməkçisi)
üçündür. Burada yalnız davam etmək üçün LAZIM olanlar var:

| Harada | Nə |
|---|---|
| `ISH_HESABATI.md` | hər seansın tam təfərrüatı və ölçmələri (son: **Seans 44**) |
| `docs/teqdimat/` | təqdimat materialları və onların linkləri (Seans 44) |
| `QARARLAR.md` | texniki qərarların səbəbləri (son: **Q-34**) |
| `SPE1.md` | SPE1-in mənbəsi, deck parametrləri, boşluq cədvəli, **§7 səbəb**, **§8 cari müqayisə** |
| `ROADMAP.md` | mərhələ statusu, texniki borc (TB-1…TB-3) |
| `CLAUDE.md` | layihə qaydaları (AI köməkçisi üçün) |

Əvvəlki təhvil sənədi (Seans 31–38 arası yenilənmiş variant):
`git show 94d04f4:TEHVIL_TESLIM.md`. Ondan əvvəlki: `git show dd425be:TEHVIL_TESLIM.md`.

---

## 1 · Bir baxışda vəziyyət

| | |
|---|---|
| Budaq | `main` = `origin/main`, açıq budaq yoxdur |
| Test dəsti | **2715 keçdi, 1 buraxıldı, 1 xfailed** (~6 dəqiqə boş maşında) |
| Son hesabat bölməsi | **Seans 44** → növbəti yazılacaq: **Seans 45** |
| Son qərar | **Q-34** → növbəti: **Q-35** |
| Aktiv blok | **B7 / SPE1** — əsas fərq bağlandı, qalan ~145 günlük fərqin səbəbi ⏳ |

### Seans 44 (19–20 sentyabr 2026) — təqdimat materialları

Kod dəyişmədi. Proqramın iş prinsipini təqdim etmək üçün slayd dəsti, interaktiv
xəritə, zaman addımı animasiyası, demo planı, icmal və lüğət hazırlandı —
bax `docs/teqdimat/README.md`. Yeni alətlər: `tools/module_graph.py`,
`tools/presentation_demo.py`. Növbəti iş dəyişmədi (SPE1 qalan fərqi).

### Seans 40–43-də bitənlər (17 sentyabr 2026)

| Seans | İş | Commit |
|---|---|---|
| 40 | **Fərqin SƏBƏBİ tapıldı** — etalonun 9 `BGSAT` bloku ilə qaz cəbhəsi ölçüldü; PVT/quyu indeksi/PERMZ deck ilə 0.00 % təsdiqləndi | `e0aadd1` |
| 41 | **G7** ✅ — doymuş Rs qolu OPM kimi xətti uzadılır (Q-32). WBHP PROD 1034-cü gündə −60.7 % → **−1.2 %** | `2073136` |
| 42 | **Vurucu mobilliyi** ✅ — hər iki mühərrikdə hüceyrənin tam mobilliyi (Q-33). WBHP INJ 1-ci gündə −34.8 % → **+1.1 %** | `3f9fd26` |
| 43 | **Doymamış qol XƏTTİ** (Q-34) — μo sapması 1.97 % → 0.000 %; zaman addımı və kro modeli ölçülərək istisna olundu | `647ac78` |

**Bu dörd seansın əsas dərsi:** hər qayda OPM-in MƏNBƏ KODUNDAN oxundu
(`LiveOilPvt.hpp`, `Tabulated1DFunction.hpp`, `StandardWell_impl.hpp`,
`EclDefaultMaterial.hpp`) — fərziyyə ilə kod yazılmadı.

### Daha əvvəl bitənlər (Seans 37–38)

| Seans | İş | Commit |
|---|---|---|
| 37 | Sınaq deck-ində səhv köçürülmüş PVTO sətri düzəldildi — "qollar 0.3 % fərqlənir" ölçməsi yanlış idi, həqiqi fərq ~11 % | `f56d95d` |
| 37 | **G3:** doymamış Bo/μo hər PVTO qolunun öz `c_o` və `n`-i ilə (`oil_branches`), Rs üzrə interpolyasiya (Q-30) | `1139e44`, `0b4e914` |
| 38 | **G9a:** cazibə həddində sıxlığın törəmələri üç fazalı Jakobiana | `67cb596` |
| 38 | **G9b:** canlı neftin sıxlığı `(ρo + Rs·ρg)/Bo` — SPE1-də +27 % (Q-31) | `a23b3b8` |
| 38 | **G9c:** ilkin tarazlıq eyni sıxlıq və Bo(p, Rs) ilə | `98f4400` |
| 38 | **Səssiz səhv:** layihə faylı PVT-nin qaz sütunlarını atırdı | `caa6e79` |
| 38 | PVT qolları modelə/servisə/layihə faylına (`OilBranch` → `domain/pvt.py`) | `344e5d7` |
| 38 | Qazın SƏTH debiti üçün ayrıca xəbərdarlıq həddi, `Mscf/day` vahidi | `fbc2134` |
| 38 | Eclipse summary oxuyucusu `imex2d/io/eclipse_summary.py` + testlər | `767756e` |
| 38 | **SPE1CASE2 qurucusu** `imex2d/benchmarks/spe1.py` + `tools/spe1_compare.py` | `097015b` |

Daha əvvəlki işlər (G1, G2, G4, G5, G6, BHP limiti, səth debiti və s.) —
`ROADMAP.md` və `ISH_HESABATI.md` Seans 27–36.

### SPE1 boşluqları — cari vəziyyət

| # | Nə | Status |
|---|---|---|
| G1–G6, G9 | PVT oxuyucuları, doymamış μo/Bo, SGOF, səth debiti, süxur istinadı, canlı neft sıxlığı | ✅ |
| G7 | Doymuş qol 5014.7 psia-dan yuxarı (Rs_sat plato) | ✅ Seans 41 (Q-32) |
| G8 | `DRSDT 0` (yalnız SPE1CASE1) | ⏳ hədəf CASE2 olduğu üçün təxirə salınıb |

---

## 2 · Yeni kompüterdə işə başlamaq

### 2.1 · Kod

```bash
git fetch origin
git status                    # yerli dəyişiklik varsa ƏVVƏL onunla məşğul olun
git pull --ff-only origin main
git log --oneline -3          # ən üstdə bu təhvil sənədinin commit-i, altında 94d04f4
```

### 2.2 · Python mühiti

`run.bat` ardıcıllıqla `.venv`, `..\venv`, `venv` qovluqlarını axtarır (Q-10, Q-13).
İşlədilmiş mühitlər — hər ikisində tam dəst keçib:

| Maşın | Mühit | Python |
|---|---|---|
| Seans 37–38 | `..\venv` (repo qovluğunun YANINDA) | 3.14.7 |
| Seans 32–36 | `.venv` | 3.12.10 |

Aşağıdakı əmrlərdə `PY` öz mühitinizin `python.exe` yoludur:

```bash
PY -m pytest -q -p no:cacheprovider -o addopts=""    # baza: 2715 keçdi, 1 buraxıldı, 1 xfailed
PY app.py                                            # proqram
```

**İşə başlamazdan ƏVVƏL bazanı öz maşınınızda alın.** 2715-dən fərqli çıxarsa,
səbəbini tapmadan dəyişiklik etməyin.

### 2.3 · SPE1 faylları (repoda DEYİL — lisenziya yoxlanılmayıb)

Repo qovluğundan KƏNARDA bir yerə endirin (təsadüfən commit olunmasın):

```bash
curl -sSfLO https://raw.githubusercontent.com/OPM/opm-tests/master/spe1/SPE1CASE2.DATA
curl -sSfLO https://raw.githubusercontent.com/OPM/opm-tests/master/spe1/opm-simulation-reference/flow/SPE1CASE2.SMSPEC
curl -sSfLO https://raw.githubusercontent.com/OPM/opm-tests/master/spe1/opm-simulation-reference/flow/SPE1CASE2.UNSMRY
sha256sum SPE1CASE2.DATA   # f3de3d06ab5705381e14e6902a21c1af7bb6b37dbc9d6e4055fe146f50a4c249
```

Checksum fərqlidirsə OPM faylı dəyişib — `SPE1.md` §2/§3-dəki rəqəmləri yenidən yoxlayın.

Yoxlama (etalonun oxunması və modelin tam qaçışı):

```bash
PY tools/eclipse_summary.py <qovluq>/SPE1CASE2 "TIME,FOPR,FGOR,WBHP:PROD,WBHP:INJ,BPR:1,BPR:300" 20
PY tools/spe1_compare.py <qovluq>/SPE1CASE2.DATA <qovluq>/SPE1CASE2    # ~80 san
```

İkinci əmr §4-dəki cədvəlin **«İNDİ» sütununu eyni rəqəmlərlə** verməlidir
(375 addım; Seans 43-də bir neçə dəfə işlədilib, nəticə təkrarlanıb). Fərqli
çıxarsa — işə başlamadan səbəbini tapın.

Alət cədvəldən əlavə **qaz cəbhəsi cədvəlini** (etalonun 9 `BGSAT` bloku)
və FOPR < 19 900 anını çap edir. **FGOR > 2 anını ÇAP ETMİR** — Seans 43-də
o, `imex2d.benchmarks.spe1.first_time_above(zaman, FGOR, 2.0)` ilə ayrıca
hesablandı. Zaman addımı həssaslığı: `--max-dt 10`.

Seans 40–43-ün maşınında etalon faylları `C:\Users\Acer\spe1-ref`
qovluğunda idi (repodan KƏNARDA, checksum yuxarıdakı ilə eyni çıxdı).

#### OPM mənbə faylları — qaydalar buradan oxundu

Seans 40–43-də hər qayda OPM-in mənbəyindən oxundu. Backlog-dakı N-1…N-4
(`SPE1.md` §8.2) üçün də eyni yol lazımdır:

| Fayl (`https://raw.githubusercontent.com/OPM/...`) | Nə üçün |
|---|---|
| `opm-common/master/opm/material/fluidsystems/blackoilpvt/LiveOilPvt.hpp` | G7 — `extrapolate=true` (sətir 512) |
| `opm-common/master/opm/material/common/Tabulated1DFunction.hpp` | xətti interpolyasiya və ekstrapolyasiya (266–283) |
| `opm-simulators/master/opm/simulators/wells/StandardWell_impl.hpp` | vurucu: `total_mob` (264–315) |
| `opm-simulators/master/opm/simulators/wells/WellInterface_impl.hpp` | `getMobility` (2261–2278) |
| `opm-common/master/opm/material/fluidmatrixinteractions/EclDefaultMaterial.hpp` | defolt kro (390–421, 447–458) |
| `opm-common/master/opm/material/fluidmatrixinteractions/EclMaterialLawManager.cpp` | kro modelinin seçimi (540–545) |

⚠️ **Tələlər:** `opm-simulators/opm/material/...` yolu **404** verir —
cari fayllar `opm-common`-dadır. `opm-material` deposu hələ də açılır, lakin
**köhnə versiyanı** verir (`LiveOilPvt.hpp` orada 706 sətir, `opm-common`-da
418) — onu İŞLƏTMƏYİN. Sətir nömrələri 17 sentyabr 2026-dakı `master`-ə
aiddir və dəyişə bilər: funksiya adı ilə axtarın.

### 2.4 · graphify

`CLAUDE.md` kod sualları üçün `graphify query` tələb edir. Seans 38-in
maşınında qurulu idi (`graphify-out/` `.gitignore`-dadır, repoya düşmür).
Sizdə yoxdursa bu addımı buraxa bilərsiniz; varsa kod dəyişəndən sonra
`graphify update .` işlədin.

### 2.5 · Proqramı işə salmaq və gözlənilməz bağlanmanı tutmaq

Arxa plan terminal əmri ilə açılan proqram seans bağlananda bağlanır (§8).
Seans 40–43-də proqram PowerShell ilə **stderr tutulması** ilə açıldı ki,
özbaşına bağlanarsa səbəbi itməsin:

```powershell
$env:PYTHONFAULTHANDLER = "1"; $env:PYTHONUNBUFFERED = "1"
Start-Process -FilePath "<PY>" -ArgumentList "app.py" -WorkingDirectory "<repo>" -RedirectStandardError "<qovluq>\app_stderr.log" -RedirectStandardOutput "<qovluq>\app_stdout.log"
```

`PYTHONFAULTHANDLER=1` C səviyyəli çökmədə (Qt, numpy) də yığın dökümünü
stderr-ə yazır. Proqram gözlənilmədən bağlananda bu ardıcıllıqla yoxlayın:

1. `app_stderr.log` — traceback və ya `Fatal Python error`;
2. `logs/imex2d.log` — son yazılar (proqram özü də jurnal aparır);
3. Windows hadisə jurnalı — `Get-WinEvent -FilterHashtable
   @{LogName='Application'; StartTime=(Get-Date).AddMinutes(-30)}`,
   `Application Error` / `Application Hang` qeydləri.

**Seans 40–43-də ölçülüb:** proqram iki dəfə "yox oldu", lakin nə jurnalda,
nə stderr-də, nə də Windows hadisə jurnalında çökmə izi tapılmadı — yəni
normal bağlanma idi. İkinci halda ölüm anı 3 saat 45 dəqiqəlik pəncərədə
idi, ona görə "qaçış bitəndən sonra çökür" fərziyyəsi TƏKZİB olundu.
Çıxışında `MainWindowTitle` boş görünən ikinci `python.exe` prosesi normaldır
— `main()` yalnız uşaq prosesdə işləyir (jurnalda bir "başladıldı" sətri).

---

## 3 · Layihə qaydaları (MÜTLƏQ — `CLAUDE.md`-dən)

1. **Hər iş seansı sənədlənir:** `ISH_HESABATI.md`-ə YENİ bölmə (köhnə bölmələr
   heç vaxt dəyişdirilmir), texniki seçim varsa `QARARLAR.md`-ə Q-N, statuslar
   `ROADMAP.md`, `ICRA_PLANI.md`, `SPE1.md`; struktur dəyişibsə `ARCHITECTURE.md`.
2. **Dil:** sənədlər Azərbaycan dilində, kod və identifikatorlar ingiliscə.
3. **Uydurma məzmun yazılmır.** Bilinməyən yer ⏳. Yazılan rəqəm ÖLÇÜLMÜŞ olmalıdır.
4. **Sənədsiz push edilmir.** Hər mənalı addımdan sonra commit + push.
5. **Nömrə toqquşması:** iki maşın paralel işləyə bilər. Yeni Seans/Q nömrəsi
   yazmazdan əvvəl `git fetch` edib `origin/main`-dəki son nömrəni yoxlayın.
6. **İş axını:** qısa budaq → testlər → sənədlər → `git merge --ff-only` →
   budağı sil → push. Böyük dəyişikliyi bir commit-ə yığmayın.
7. **Nəticəni dəyişən hər seçimi sahibkara verin** (G3 və G9 belə soruşulub).

---

## 4 · SPE1CASE2 — hazırkı nəticə

Model: [imex2d/benchmarks/spe1.py](imex2d/benchmarks/spe1.py). Cədvəllər deck-dən
oxunur, qalan parametrlər `SPE1.md` §2-də yoxlanılmış FIELD sabitləridir.
Qaçış: **375 addım**, orta Δt 9.7 gün, 146 təkrar (~2 dəqiqə).

| t, gün | Kəmiyyət | Seans 38-də | **İNDİ** | OPM Flow |
|---|---|---|---|---|
| 1 | WBHP INJ, psia | 5271 (−34.8 %) | **8173 (+1.1 %)** | 8082 |
| 304 | FOPR, STB/gün | 20 000 (0.0 %) | 20 000 (0.0 %) | 20 000 |
| 1034 | FGOR, Mscf/STB | 6.255 (+389 %) | **1.776 (+38.7 %)** | 1.280 |
| 1034 | WBHP PROD, psia | 1576 (−60.7 %) | **3965 (−1.2 %)** | 4013 |
| 1399 | FOPR, STB/gün | 14 720 (−26.4 %) | **19 450 (−2.8 %)** | 20 000 |
| 3650 | FOPR, STB/gün | 4980 (−13.1 %) | **5176 (−9.7 %)** | 5733 |
| 3650 | BPR (1,1,1) / (10,10,3) | −4.2 / −3.7 % | **−2.7 / −1.9 %** | 4101 / 3278 |

| Hadisə | Seans 38-də | **İNDİ** | OPM Flow |
|---|---|---|---|
| FGOR > 2 (qazın çatması) | 900 gün | **1132 gün** | 1276 gün |
| FOPR < 19 900 (limitə keçid) | 1120 gün | **1382 gün** | 1550 gün |
| Qaz cəbhəsi, blok 300 | — | **1161 gün** | 1307 gün |

Tam cədvəl və qaz cəbhəsi: `SPE1.md` §8;  alət:
`PY tools/spe1_compare.py <deck> <etalon>` (cəbhə cədvəlini də çap edir).

**Oxunuşu:** 1034-cü günə qədər FOPR dəqiq, təzyiqlər ≤ 5.7 % daxilindədir.
Qalan fərq 1399–1580-ci günlərdə toplanıb: qaz hələ ~145 gün tez gəlir.

**Etalon (golden) fayl YAZILMAYIB** — fərqlər izah olunmayana qədər yazılmamalıdır (Q-31 Qərar 7).

### 4.1 · Seans 40–43-ün MÖVCUD modellərə təsiri — köhnə nəticələrlə müqayisədə nəzərə alın

| Dəyişiklik | Hansı modellərə təsir edir | Ölçülmüş təsir |
|---|---|---|
| G7 (Q-32) — doymuş qolun xətti uzadılması | YALNIZ deck qolları (`pvt_oil_branches`) olan modellər | SPE1: WBHP PROD, 1034-cü gün −60.7 % → −1.2 % |
| **Q-33 — vurucuda tam mobillik** | **vurucusu olan BÜTÜN modellər, iki fazalı su vurulması da** | qazsız sınaq modelində RF 62.86 → 62.83 %, addım 31 → 29; vurmanın başlanğıcında BHP daha YÜKSƏK |
| Q-34 — doymamış qol xətti | YALNIZ deck qolları olan modellər | SPE1: FGOR > 2 1142 → 1132 gün |

**Sahibkarın iş faylı `layihe.imx` (Seans 43-də faylın özündə yoxlanıldı):**
PVT mənbəyi `correlation(API=32, γg=0.75, T=70°C)`, doyma təzyiqi 240 bar,
`pvt_oil_branches` açarı YOXDUR (fayl 14 sentyabrda, qollar layihə faylına
əlavə olunmazdan əvvəl saxlanılıb). Quyular: 1 su vurucusu (INJ) + 2
istismarçı. Deməli:

* G7 və Q-34 bu modelə **təsir ETMİR**;
* Q-33 **təsir EDİR** — köhnə qaçışlarla müqayisədə vurucunun BHP-si
  başlanğıcda daha yüksək, RF isə cüzi fərqli çıxacaq. Bu, səhv deyil,
  sahibkarın təsdiqlədiyi dəyişiklikdir.

`.imx` faylı **gzip ilə sıxılmış JSON-dur** — `grep` ilə oxunmur, Python
`gzip.open(..., "rt")` + `json.load` işlədin.

---

## 5 · NÖVBƏTİ İŞ — fərqin səbəblərini bir-bir ölçmək

> **DİQQƏT — bu bölmə köhnəlib.** §5.1 (vurucu mobilliyi) və §5.2 (G7)
> ARTIQ BİTİB (Seans 41–42). Aşağıdakı təsvirlər tarixi kontekstdir.
>
> **CARİ VƏZİYYƏT:** əsas səbəb tapılıb və bağlanıb. Qalan fərq ~145 gündür
> və ÜÇ namizəd ölçülərək İSTİSNA olunub (`SPE1.md` §8.1):
>
> | Namizəd | Ölçmə | Nəticə |
> |---|---|---|
> | Zaman addımı | `--max-dt` 31 → 10: cəbhə 13–28 gün DAHA TEZ | səbəb deyil |
> | Üç fazalı kro modeli | OPM defoltu ilə orta 0.174 %, istiqamət TƏRS | səbəb deyil |
> | Doymamış qolun forması | Q-34 ilə düzəldildi, təsir ±10 gün | səbəb deyil |
>
> **Ölçülərək TƏSDİQLƏNƏNLƏR** (bunlara qayıtmayın): Bg, μg, Bo, μo,
> Rs_sat, c_o deck ilə **0.00 %**; quyu indeksi analitik Peaceman ilə eyni;
> `PERMZ = PERMX`, `PORO 0.3` deck ilə eyni.
>
> **Növbəti namizədlər ⏳** (heç biri ölçülməyib): qaz-neft sistemində
> upstream çəkiləməsi; cazibə/şaquli axının cəbhəyə təsiri (G9-dan sonra
> yenidən); istismarçının qaz hasilatının bölünməsi; `DRSDT`-siz Rs
> artımının sürəti. Bir müşahidə: 1034-cü gündə quyudibi təzyiqimiz
> 3965 psia, OPM-də 4013 — hər ikisi doyma nöqtəsinin (4014.7) düz
> ətrafındadır və sistem orada həddindən həssasdır (48 psi → FGOR-da 38 %).

Aşağıdakı bölmələr TARİXİ kontekstdir (Seans 39-da yazılıb).

### 5.1 · Vurucu bağlantısının mobilliyi — tövsiyə olunan başlanğıc

* **Nə görünür:** 1-ci gündə vurucunun BHP-si bizdə 5271, OPM-də 8082 psia.
  Yəni bizim injektivlik xeyli yüksəkdir. Bu, ən birbaşa və ən ucuz yoxlanılan
  əlamətdir.
* **Bizdə:** vurulan fazanın **son nöqtə** mobilliyi (`krg_end/μg`, SGOF-da
  krg_end = 0.984) — `simulation/implicit/three_phase_residual.py::ThreePhaseWellModel`
  (sinfin sənədi: "Vurulan fazanın mobilliyi SON NÖQTƏ mobilliyidir").
* **Yoxlanılmalı:** OPM Flow vurucu bağlantısında hansı mobilliyi işlədir
  (`opm-simulators` mənbəsi, `StandardWell`/`WellInterface`). Qaydanı
  mənbədən oxumadan yazmayın.
* **Ölçmə təklifi:** qaydanı müvəqqəti dəyişib `tools/spe1_compare.py`-ı
  işlədin və 1-ci günün WBHP INJ-i ilə qazın çatma vaxtının necə dəyişdiyinə baxın.
* ⚠️ Quyu Jakobianı dəyişərsə sonlu fərq testi MÜTLƏQDİR (TB-1 bu sinifdədir).

### 5.2 · G7 — doymuş qolun 5014.7 psia-dan yuxarı davranışı

* Etalonda vurucu hüceyrəsinin təzyiqi 7534 psia-ya qalxır, PVTO-nun doymuş
  qolu isə 5014.7-də bitir. Bizdə `np.interp` sərhəddə saxlayır (Rs_sat = 1.618 plato).
* OPM-in bu haldakı qaydası (ekstrapolyasiya, yoxsa plato) mənbədən yoxlanılmayıb.
* **TB-3 eyni kökdəndir:** `black_oil.py::_saturation_pressure_slope` ən üst
  Rs düyünündə analitik 3.619, sonlu fərq 1.810.

### 5.3 · Üç fazalı nisbi keçiricilik

* Bizdə Stone II (`simulation/stone_relperm.py`); deck-də `STONE1`/`STONE2`
  açar sözü yoxdur, yəni Eclipse/OPM öz defolt modelini işlədir.
* Su hərəkətsiz qalanda (Sw ≈ Swc) hər ikisi kro = krog verir, lakin qaz
  zonasında fərq ola bilər — ölçülməyib.
* `stone_relperm.py`-nin modul sənədində "Stone II … Eclipse-in defoltu"
  yazılıb — bu, yanlışdır və düzəldilməlidir.

### 5.4 · Zaman addımı

* 197 təkrar həll; bir dəfə "Səth debiti: düzəldilmiş hədəflə addım yenidən
  yığılmadı" xəbərdarlığı. `--max-dt` dəyişərək (`tools/spe1_compare.py
  --max-dt 10`) nəticənin addımdan asılılığını ölçün.

### 5.5 · Fərqlər izah olunandan sonra

* `tests/golden/`-da SPE1CASE2 reqressiyası (bax `tools/golden.py`).
* CASE1 (`DRSDT 0`) — Rs artımını qadağan edən qayda, qalığa toxunur.

---

## 6 · Mühərrikin əsas qaydaları

* **Qalıq və Jakobiana mümkün qədər toxunmamaq.** Quyu qaydaları bağlantı
  obyektinin sahələrini (`mode`, `target`, `well_index`, `rate_share`) Nyuton
  həlləri arasında dəyişir — `simulation/well_constraints.py` (THP:
  `simulation/wellbore/thp_control.py`).
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

* **Qalıq/Jakobian dəyişibsə — sonlu fərq testi MÜTLƏQDİR**, xəta sütun NÖVÜ
  üzrə ölçülür (p / Sw / 3-cü dəyişən):
  `tests/test_undersaturated_viscosity_engine.py::_column_kind_errors`.
  Cazibəni yoxlamaq üçün ÇOX LAYLI model lazımdır —
  `tests/test_live_oil_density.py::_engine(nz=3)`.
* **Səssiz yanlış davranış = ən pis səhv.** Dəstəklənməyən kombinasiya açıq
  xəta və ya diaqnostika xəbərdarlığı verməlidir.
* **Hər iddia ölçmə ilə.** Seans 36-da "RF fizikaya görə düşüb" izahı, Seans
  38-də "TB-2 cazibədən gəlir" fərziyyəsi ölçmə ilə təkzib olundu.

---

## 7 · Texniki borc və açıq qalan backlog

**Texniki borc** (`ROADMAP.md` → «Texniki borc»):

* **TB-1:** üç fazalı RATE istismarçısının quyu Jakobianı — sonlu fərqə qarşı
  xəta 0.4893 (tək perforasiya). Sahibkarın qərarı: indi düzəldilmir.
* **TB-2:** üç fazalı axın Jakobianının təzyiq sütunu qarışıq vəziyyətdə —
  0.03–0.19. Seans 38-də ölçüldü: **cazibədən DEYİL** (cazibəsiz modeldə də
  var, G9a dəyişmədi). Mənbə ⏳.
* **TB-3:** `_saturation_pressure_slope` ən üst Rs düyünündə 2 dəfə fərq —
  **deck yolunda Seans 41-də (G7/Q-32) həll olundu**, korrelyasiya yolunda qalır.
* **Vurucu qaydası birləşmiş həllediciyə (`implicit/standard_well.py`)
  TƏTBİQ OLUNMAYIB** (Q-33 Qərar 4). O yol heç bir mühərrik tərəfindən
  işlədilmir (yatmış kod), lakin fayla toxunanda uyğunlaşdırılmalıdır.
* **OPM-in defolt üç fazalı kro modeli əlavə edilməyib** (Seans 43-də
  ölçüldü: SPE1-də fərq cəmi 0.174 %). Su HƏRƏKƏT EDƏN modellərdə fərq
  50 %-ə çatır, ona görə gələcəkdə əlavə edilməlidir. Düstur və mənbə
  sətirləri `SPE1.md` §8.1-də yazılıb.
* **Q-33-ün yığılmaya yan təsiri** (sahibkar bilərək qəbul etdi, Seans 43):
  vurma debiti doyuma bağlandığı üçün quyu ətrafındakı sərt ssenaridə CNV
  minimumu **0.001004 → 0.001322** (`test_implicit_newton.py`-nin
  oscillasiya ssenarisi; hər iki halda həmin addım yığılmır, mühərrik Δt-ni
  kəsir). Geri-izləmə hələ də işləyir — qalığın sıçrayışı 1.4, geri-izləməsiz
  9.6. Test mütləq həddən (0.0012) oscillasiyanı ölçən şərtə keçirilib.

**Açıq qalan ⏳:**

* **SPE1CASE2-nin qalan ~145 günlük fərqi** — dörd namizəd və hər birinin
  ÖLÇMƏ ÜSULU `SPE1.md` §8.2-də yazılıb (şaquli axın/cazibə, upstream
  çəkiləməsi, istismarçıda qazın bölünməsi, Rs artımının sürəti). Orada
  həm də "qovmayın" siyahısı var: lay təzyiqinin aşağı olması müstəqil
  səbəb deyil, FGOR-un yüksəkliyinin nəticəsidir.
* **Etalon (golden) fayl YAZILMAYIB** və fərqlər izah olunana qədər
  yazılmamalıdır (Q-31 Qərar 7, sahibkarın təsdiqi ilə Seans 43-də də
  saxlanıldı).

* THP hidravlikası (`simulation/wellbore/`) hələ ölü neft sıxlığı ilə hesablayır.
* İki fazalı Jakobianda cazibə sıxlığının törəməsi atılır — ölçülməyib.
* SWOF/SGOF-un Pc sütunu vahid çevirməsiz oxunur (SPE1-də 0).
* Quyu lüləsində hidrostatik hədd yoxdur (BHP bütün perforasiyalara eyni).
  SPE1-də problem deyil: hər quyuda bir perforasiya var.
* IMPES süxur sıxılmasının istinad təzyiqini (G6) oxumur və buna görə xəbərdarlıq vermir.
* Başlanğıcdakı debit titrəməsi ("mişar dişi") — sahibkarın modelində görünüb,
  təkrarlanmayıb; `layihe.imx` və ya o qaçışın CSV-si lazımdır.
* `well_control_mode` CSV/JSON ixracına və dashboard-a çıxarılmayıb;
  `ui/main_window.py::export_results` qaz sütunlarını yazmır.
* THP quyusu ilə BHP limitli quyu eyni modeldə — limit dövrəsindən sonra
  THP-nin BHP-si yenilənmir (ölçülməyib).
* Köhnə backlog: ilkin tarazlıqda Pcog yoxdur; Eclipse ixracı iki fazalıdır;
  cədvəldən Pcog; Vogel IPR; Beggs-Brill, VFPPROD; üç fazalı mühərrikdə
  `soft_failure_*` tolerantlıqları.

---

## 8 · Bilinən tələlər (vaxt itirməmək üçün)

* **`SPE1.md` CRLF-dir** (bu sənədin aşağısındakı sətir sonları cədvəli onu
  səhvən LF kimi göstərirdi — Seans 43-də `file` ilə yoxlanıldı). Hər fayl
  dəyişməzdən ƏVVƏL `file <ad>` işlədin, siyahıya güvənməyin.
* **Bash-da ikiqat dırnaq içində backtick ` KOMANDA KİMİ işləyir** —
  Azərbaycanca şərh yazanda mətn səssizcə itir (Seans 43-də bir şərh
  belə zədələndi). Heredoc (`<<'PYEOF'`) işlədin.
* **Ölçməni tətbiqdən ƏVVƏL iki dəfə yoxlayın.** Seans 43-də kro
  müqayisəsində `krow`-u səhv arqumentdə (sw, halbuki OPM sg+sw işlədir)
  hesabladım və sahibkara YANLIŞ rəqəm verdim (1.78 % ↔ həqiqi 0.174 %).
  Düsturu mənbədən oxuyanda HƏR arqumentin nə olduğunu yoxlayın.
* **Windows konsolu cp1254-dür.** Azərbaycan hərfləri çap edən skriptin
  əvvəlinə: `sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")`.
* **Sətir sonları qarışıqdır.** Bəzi fayllar CRLF-dir (`ISH_HESABATI.md`,
  `QARARLAR.md`, `ROADMAP.md`, `ARCHITECTURE.md`, `black_oil.py`,
  `serialization.py`…), bəziləri LF (`SPE1.md`, bu fayl, `three_phase_residual.py`,
  `pvt_io.py`, `benchmarks/spe1.py`…). Dəyişdikdən sonra `file <ad>`
  ilə yoxlayın. **Git Bash-dakı `sed -i` CRLF-i LF-ə çevirir** (Seans 37-də
  bütün fayl dəyişmiş göründü) — Python ilə bayt səviyyəsində yazın.
* **Arxa planda gedən tam dəst** modulları başlanğıcda import edir: sonrakı kod
  dəyişikliyi o qaçışa təsir etmir, sonradan yaradılan test faylı toplanmır.
  Paralel başqa iş gedəndə dəst 10–11 dəqiqəyə uzanır.
* **Test modellərinin çoxu tək laylıdır (`nz=1`) — orada cazibə həddi İŞLƏMİR.**
  Cazibə ilə bağlı hər ölçmə çox laylı modeldə aparılmalıdır.
* **`engine.initialization` mühərrik qurulanda seçilir.** Sonradan
  `initial_conditions.use_equilibration` dəyişsəniz tarazlıq işə düşmür — ya
  mühərriki yenidən qurun, ya provider-i əl ilə verin
  (`tests/test_live_oil_density.py::_equilibrated_initial_state`).
* **`DeckPvt.to_pvt_table()` qollar verilsə belə** "yalnız 1 düyün … c_o"
  xəbərdarlığı yazır — `oil_branches` provider-ə verilibsə gözləniləndir.
  `reference_rs` verməyin: doymuş qol kəsilir (ölçülüb: Bo −9 %, μo +17 %).
* **SPE1-də FOPR ilk günlərdə hədəfdən ~0.1 % aşağı olur** (19 981) — səth
  debiti tənzimləyicisinin toleransıdır, rejim keçidi deyil. Keçidi 0.5 %
  toleransla axtarın.
* **PVT cədvəlində doyma təzyiqi DÜYÜN DEYİL** — lövbər fit-dən alınır (Q-29).
* **Testin keçməsi düzgünlüyün sübutu deyil.** Seans 34-ün səhv rəqəmi
  testdəki səhv köçürülmüş deck sətrindən gəlirdi; Seans 36-da lövbər testi
  eyni səhvi bölüşən iki tərəfi müqayisə edirdi. Rəqəmi həmişə orijinal
  mənbə ilə tutuşdurun.
* **Yalnız RATE quyulu model** servis doğrulamasından keçmir — ən azı bir
  quyuda BHP limiti olmalıdır (SPE1-də hər ikisində var).
* **`MainWindow()` testdə yaradılsa pytest çökür** — UI-ni mənbə yoxlaması ilə sınayın.
* **Proqramı arxa plan terminal əmri ilə açsanız seans bağlananda bağlanır** —
  PowerShell: `Start-Process -FilePath "<PY>" -ArgumentList "app.py"`.
* **`*.imx` `.gitignore`-dadır**; `layihe.imx` (sahibkarın iş faylı) git-də deyil.
* **Ölçmə skriptləri repoda yoxdur** — rəqəmlər `ISH_HESABATI.md`-də ölçmə
  şərtləri ilə yazılıb. Repoda olan alətlər: `tools/eclipse_summary.py`,
  `tools/spe1_compare.py`, `tools/golden.py`.

---

## 9 · Toxunmamalı olan

* **`b2a795e` commit-i** yalnız digər maşındadır və onu **sahibkar özü xilas
  edəcək**. Birləşdirmədən əvvəl `git log --all --oneline | grep b2a795e` ilə
  yoxlayın və sahibkarla razılaşdırın — özbaşına birləşdirməyin, silməyin.

---

*Bu fayl təhvil anının şəklidir. İş davam etdikcə cari vəziyyət `ISH_HESABATI.md`,
`ROADMAP.md` və `SPE1.md`-dədir; bu fayl köhnələrsə, həmin sənədlər əsasdır.*
