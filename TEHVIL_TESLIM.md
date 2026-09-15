# Təhvil-təslim — işi başqa kompüterdə davam etdirmək üçün

**Hazırlanıb:** 15 sentyabr 2026 · **Son commit:** `6b40a88` (main, GitHub ilə sinxron)
**Növbəti iş:** B7 addım 2 — RATE rejimində BHP limiti (bax §5)

Bu sənəd işi yarımçıq qaldığı yerdən götürən şəxs (insan və ya AI köməkçisi)
üçündür. Tarixçənin tam təfərrüatı `ISH_HESABATI.md`, qərarların səbəbləri isə
`QARARLAR.md` faylındadır — burada yalnız davam etmək üçün LAZIM olanlar var.

---

## 1 · Bir baxışda vəziyyət

| | |
|---|---|
| Budaq | `main` = `origin/main` = `6b40a88` |
| Test dəsti | **2490 keçdi, 1 buraxıldı, 1 xfailed** (~5.5 dəqiqə) |
| Son hesabat bölməsi | **Seans 26** → növbəti yazılacaq: **Seans 27** |
| Son qərar | **Q-18** → növbəti: **Q-19** |
| Aktiv blok | **B7** (yekun doğrulama + SPE1 etalonu) |

Bitmiş son işlər (yenidən başlamağa ehtiyac yoxdur):

| Seans | İş | Commit |
|---|---|---|
| 23 | B4-B — quyunun THP ilə idarəsi | `3c67724` |
| 24 | B5-b — üç fazalı kapilyar təzyiq (Pcow qoşuldu, Pcog əlavə olundu) | `60a1393` |
| 25 | THP sabitliyi — nodal analiz (IPR ∩ VLP), yarı-implicit təkrar, avtomatik bağlanma | `372d4f4` |
| 26 | B7 addım 1 — qaz vuran quyu | `6b40a88` |

---

## 2 · Başqa kompüterdə işə başlamaq

```bash
git fetch origin
git status                 # yerli dəyişiklik varsa ƏVVƏL onunla məşğul olun
git pull --ff-only origin main
git log --oneline -1       # 6b40a88 və ya daha yeni olmalıdır
```

Virtual mühit: `run.bat` ardıcıllıqla `.venv`, `..\venv`, `venv` qovluqlarını
axtarır (bax `QARARLAR.md` → Q-10, Q-13). Hazırlandığı maşında `venv` işlədilib,
Python 3.14.3.

```bash
venv\Scripts\python.exe -m pytest -q -p no:cacheprovider   # baza: 2490 keçməlidir
venv\Scripts\python.exe app.py                              # proqramı açmaq
```

**Başlamazdan əvvəl baza test nəticəsini ÖZ maşınınızda alın.** 2490-dan fərqli
çıxarsa, işə başlamadan səbəbini tapın — mühit fərqi ola bilər.

---

## 3 · Layihə qaydaları (MÜTLƏQ — `CLAUDE.md`-dən)

1. **Hər iş seansı sənədlənir:** `ISH_HESABATI.md`-ə YENİ bölmə (köhnə bölmələr
   heç vaxt dəyişdirilmir, yalnız əlavə olunur), texniki seçim varsa
   `QARARLAR.md`-ə Q-N, statuslar `ROADMAP.md` və `ICRA_PLANI.md`-də,
   struktur dəyişibsə `ARCHITECTURE.md`.
2. **Dil:** sənədlər Azərbaycan dilində, kod və identifikatorlar ingiliscə.
3. **Uydurma məzmun yazılmır.** Bilinməyən yer ⏳ ilə açıq göstərilir.
   Rəqəm yazılırsa, ÖLÇÜLMÜŞ olmalıdır.
4. **Sənədsiz push edilmir.** Hər mənalı addımdan sonra commit + push.
5. **Nömrə toqquşması:** iki maşın paralel işləyib. Yeni Seans/Q nömrəsi
   yazmazdan əvvəl `git fetch` edin və `origin/main`-dəki son nömrəni yoxlayın.
   Toqquşma olarsa, ƏVVƏL dərc olunan nömrəsini saxlayır.
6. **İş axını:** qısa budaq → testlər → sənədlər → `git merge --ff-only` → budağı
   sil → push.

---

## 4 · Mühərrikin əsas fəlsəfəsi (son seansların dərsi)

Son dörd seansda eyni prinsip dəfələrlə təsdiqləndi — yeni işdə də gözlənilir:

* **Qalıq və Jakobiana mümkün qədər toxunmamaq.** Quyu idarəsi yeni rejimləri
  (THP, avtomatik bağlanma) bağlantı səviyyəsində, addımlar arasında həll olunur:
  `ThpController` bağlantının `target`/`well_index`-ini yerində dəyişir, qalıq
  isə sadəcə yeni ədədi görür. Bax Q-15, Q-17.
* **Hər iddia ölçmə ilə.** Sahibkarın ilkin diaqnozu ("açıq THP çökür") ölçmə
  ilə qismən təkzib olundu və əsl səbəb başqa çıxdı (Seans 25). Düzəlişdən
  əvvəl səbəbi variantları ayıraraq ölçün.
* **Səssiz yanlış davranış = ən pis səhv.** Bu seanslarda tapılan boşluqların
  hamısı səssiz idi: Pcow atılırdı, `injected_phase` oxunmurdu, THP ədədi debit
  kimi ixrac olunardı. Dəstəklənməyən kombinasiya **açıq xəta** verməlidir.
* **Jakobian dəyişibsə — sonlu fərq testi MÜTLƏQDİR.** Nümunələr:
  `tests/test_three_phase_capillary.py`, `tests/test_gas_injection.py`.

---

## 5 · NÖVBƏTİ İŞ — B7 addım 2: RATE rejimində BHP limiti

### 5.1 Niyə lazımdır

SPE1 etalonunda istismarçı **debitlə** idarə olunur, lakin quyu dibi təzyiqi
**minimal həddən aşağı düşə bilməz**; lay tükəndikcə quyu debiti saxlaya bilmir və
**BHP idarəsinə keçir**. Bizdə RATE rejimində heç bir limit və rejim keçidi yoxdur.

### 5.2 Mövcud kod — RATE budaqları harada

`connection.mode is ControlMode.BHP` yoxlamasının `else` budağı RATE-dir:

| Fayl | Yer |
|---|---|
| `simulation/implicit/three_phase_residual.py` | `ThreePhaseWellModel.well_rates` — istismarçı RATE budağı: `total = -abs(connection.target)`, su/neft mobillik payına görə bölünür |
| `simulation/implicit/three_phase_residual.py` | `ThreePhaseWellJacobian.blocks` — eyni budağın törəmələri |
| `simulation/implicit/residual.py` | `well_rates` (iki fazalı) |
| `simulation/implicit/jacobian.py` | `_wells` (iki fazalı) |
| `simulation/impes_engine.py` | bir neçə `c.mode is ControlMode.BHP` yoxlaması |
| `simulation/implicit/standard_well.py`, `well_state.py` | birləşmiş (coupled) Nyuton yolu |

⚠️ **TƏHLÜKƏ (B4-B-də tapılıb):** yeni `ControlMode` dəyəri əlavə etsəniz, yuxarıdakı
BÜTÜN `mode is BHP` yoxlamaları onu **səssizcə RATE kimi** işlədəcək. B4-B-də THP
bu səbəbdən bağlantıya BHP kimi ötürüldü (bax Q-15).

### 5.3 Tövsiyə olunan dizayn (qəbul etmək məcburi deyil, amma əsaslandırılıb)

THP-də işləyən yanaşmanın təkrarı — **Jakobiana toxunmadan**:

1. **Domen:** `WellControl`-a `bhp_limit: Optional[float] = None` (istismarçıda
   minimal, vurucuda maksimal). `None` → köhnə davranış, nəticə bit-bit eyni.
   Layihə faylı (`application/serialization.py`), UI sütunu, diaqnostika.
2. **Nəzarətçi** (məs. `simulation/wellbore/` altında və ya ayrıca modul):
   bağlantının `mode`/`target`-ini **addımlar arasında** dəyişir:
   * RATE quyusu: addımdan sonra debitə uyğun BHP qiymətləndirilir
     (Peaceman: `q = Σ WI·λ_t·(p_hüceyrə − BHP)` → BHP üçün həll). BHP limitdən
     aşağıdırsa → `mode = BHP`, `target = bhp_limit`.
   * Limitdə olan quyu: limit BHP-də debit tutumu hədəfdən böyükdürsə → RATE-ə qayıt.
   * **Histerezis mütləqdir** (THP-də `REOPEN_MARGIN_BAR` kimi) — yoxsa hər
     addımda rejim atılır və Nyuton Seans 25-dəki kimi rəqs edər.
3. **Gecikmə:** mühərriklərdə artıq `_thp_outer_loop` və
   `AdaptiveTimeStepper.resolve_step` var — addımı yeni rejimlə təkrar həll etmək
   üçün eyni mexanizmi işlədin (tarixçə şişmir).
4. Rejim yalnız **Nyuton həlləri ARASINDA** dəyişməlidir — iterasiya daxilində
   dəyişmək qalıqda sıçrayış yaradır.

### 5.4 Qəbul meyarları

* `bhp_limit = None` olan bütün mövcud modellər **bit-bit eyni** nəticə verir
  (tam dəst 2490 keçməlidir).
* Limitə çatan quyunun rejim keçidi testlə göstərilir; keçiddən sonra BHP limitdən
  aşağı düşmür.
* Histerezis testi: sərhəddə aç-qapa rəqsi yoxdur.
* Real qaçışda Δt davranışı ölçülür (Seans 25-dəki cədvəl formatında).
* Eclipse ixracı: `WCONPROD` BHP limit sahəsi (`io/eclipse_export.py` — hazırda
  `'LRAT' 2* ...` yazır) yoxlanılır.

---

## 6 · B7 addım 3 — səth debiti və SPE1 modeli

### 6.1 Səth debiti hədəfi

Bizdə RATE hədəfi **maye (su+neft), LAY HƏCMİ**-dir (`well_rates` RATE budağı
debiti `B`-yə bölməzdən əvvəl mobillik payına görə paylayır). SPE1 isə **neftin
SƏTH debitini** (STB/gün) təyin edir. Mövcud RATE modellərinin mənasını səssizcə
dəyişməmək üçün bunu AYRI rejim/parametr kimi etmək tövsiyə olunur — §5.2-dəki
TƏHLÜKƏNİ nəzərə alaraq.

### 6.2 SPE1 parametrləri — ⏳ İSTİFADƏDƏN ƏVVƏL MƏNBƏDƏN YOXLANILMALIDIR

Aşağıdakılar ədəbiyyatda geniş yayılmış dəyərlərdir, lakin bu repoda HƏLƏ
yoxlanılmayıb. Model qurmazdan əvvəl orijinal mənbələrlə tutuşdurun:

* Odeh, A.S. (1981), *Comparison of Solutions to a Three-Dimensional Black-Oil
  Reservoir Simulation Problem*, JPT;
* OPM layihəsinin `SPE1CASE1.DATA` deck faylı (açıq mənbə).

| Parametr | Yayılmış dəyər |
|---|---|
| Grid | 10 × 10 × 3, DX = DY = 1000 ft |
| Təbəqə qalınlıqları | 20 / 30 / 50 ft |
| Keçiricilik (təbəqələr) | 500 / 50 / 200 mD |
| Məsaməlilik | 0.3 |
| Vurucu | (1,1,1) hüceyrəsi, **qaz**, 100 MMscf/gün |
| İstismarçı | (10,10,3) hüceyrəsi, 20 000 STB/gün neft, minimal BHP 1000 psia |
| İlkin təzyiq | 4800 psia (datum 8400 ft) |
| İlkin su doymuşluğu | 0.12 |
| Doyma təzyiqi | ~4014.7 psia |
| Müddət | 10 il |

Mühərrik metrik vahidlərdədir (bar, m³/gün) — çevirmə üçün
`domain/unit_conversions.py` mövcuddur.

### 6.3 B7-nin qalan bəndləri (`ICRA_PLANI.md` → B7)

* uc-uca reqressiya ssenarisi `tests/golden/` altında;
* test dəstinin sürətləndirilməsi (`pytest-xdist`).

---

## 7 · Açıq qalan backlog (təcili deyil)

* İlkin tarazlıqda Pcog yoxdur — qaz-neft sərhədi kəskin qalır (Seans 24).
* Eclipse ixracında SGOF / Pcog sütunu yoxdur.
* Cədvəldən Pcog (`SaturationTableSet`) — yalnız analitik Brooks-Corey var.
* IPR düz xətdir; Vogel tipli əyri ⏳ (Seans 25).
* Beggs-Brill sürüşməsi, VFPPROD idxalı, RATE quyusunda THP.
* Üç fazalı mühərrikdə `soft_failure_*` tolerantlıqları verilmir
  (`three_phase_engine.py::_time_config` sənədində qeyd var).

**Mövcud Jakobian qeyri-dəqiqlikləri** (yeni işin səbəbi DEYİL — müqayisə üçün baza):

* qarışıq doyma vəziyyətində sonlu fərq xətası 2.7×10⁻⁵ (Seans 24);
* su vurucusunda 1.343141869382806×10⁻⁵ (Seans 26, köhnə kodda da eyni ölçülüb);
* köhnə qeyd: istismarçı hüceyrəsində qaz tənliyi ↔ Sw elementində 87 % xəta
  (`ISH_HESABATI.md`-də "87 %" axtarın).

---

## 8 · Bilinən tələlər (vaxt itirməmək üçün)

* **Windows konsolu cp1254-dür.** Azərbaycan hərfləri çap edən ad-hoc skriptlər
  `UnicodeEncodeError` verir. Skriptin əvvəlində:
  `sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")`.
  Bu sarğını iki dəfə (məs. iki modul import edəndə) qurmayın — birinci obyekt
  silinəndə axın bağlanır.
* **Sənədlər CRLF sətir sonludur.** Faylı proqramla dəyişəndə sətir sonlarını
  qoruyun, yoxsa `git diff` bütün faylı dəyişmiş göstərir.
* **Uzun Python heredoc-u bash-da qırıla bilər.** Yamağı ayrıca `.py` faylına
  yazıb işlədin.
* **`MainWindow()` testdə yaradılsa pytest çökür.** UI məntiqini ayrıca funksiyada
  və ya mənbə yoxlaması ilə sınayın.
* **Proqramı arxa plan terminal əmri ilə açsanız, seans bağlananda proqram da
  bağlanır.** Ayrıca proses kimi açın, məs. PowerShell:
  `Start-Process -FilePath "venv\Scripts\python.exe" -ArgumentList "app.py"`.
* **`layihe.imx` git-də DEYİL** (sahibkarın iş faylıdır, izlənmir). Seans 25-in
  41×41×3 modeli odur və saxlanmış halda **qaz fazası söndürülüb**. Həmin
  qaçışları təkrarlamaq lazımdırsa, faylı ayrıca köçürün və PVT-də qazı açın.
* **Ölçmə skriptləri repoda yoxdur** (müvəqqəti qovluqda idi). Seans 25-in
  nəticələrini təkrarlamaq üçün: `layihe.imx` yüklənir, PVT
  `build_pvt_table(..., include_gas=True)` ilə əvəz edilir, `engine.time_stepper.history`
  üzərindən Δt / kəsilmə sayılır.

---

## 9 · Toxunmamalı olan

* **`b2a795e` commit-i** yalnız digər maşındadır və onu **sahibkar özü xilas
  edəcək**. İstənilən birləşmədən əvvəl `git log --all --oneline | grep b2a795e`
  ilə yoxlayın və sahibkarla razılaşdırın — özbaşına birləşdirməyin, silməyin.

---

*Bu fayl təhvil anının şəklidir. İş davam etdikcə cari vəziyyət `ISH_HESABATI.md`
və `ROADMAP.md`-dədir; bu fayl köhnələrsə, həmin sənədlər əsasdır.*
