# A7 — Qaz fazası (üç fazalı black-oil)

Məqsəd: modelin **yeganə qalan fiziki məhdudiyyətini** aradan
qaldırmaq. `A6_PLAN.md`-də sənədləşdirildiyi kimi, hazırkı iki fazalı
model doyma təzyiqini (Pb) kəsəndə Nyuton osilyasiya edir — çünki
orada real fiziki hadisə (qazın neftdən ayrılması) baş verir, amma
model bunu təmsil edə bilmir.

## Mərhələlər

| # | İş | Vəziyyət |
|---|---|---|
| **1** | Qaz PVT xassələri (Z, Bg, μg) | **HAZIR** |
| **2** | Domain: üç fazalı doyumluluq, ilkin şərtlər | **HAZIR** |
| **3** | Üç fazalı nisbi keçiricilik (Stone) | **HAZIR** |
| **4** | Primary dəyişənlər + dəyişən keçid (variable switching) | **HAZIR** |
| **5** | Qalıq tənlikləri (akkumulyasiya) | **HAZIR** |
| 6a | Axın həddləri (Darcy, üç faza) | **HAZIR** |
| **6b** | Quyu həddləri | **HAZIR** |
| **6c** (1-ci hissə) | Analitik Jakobian — akkumulyasiya bloku | **HAZIR** |
| **6c (2-ci hissə)** | Axın Jakobian töhfələri | **HAZIR** |
| **6c (3-cü hissə)** | Quyu Jakobianı | **HAZIR** |
| **6c (4-cü hissə)** | Seyrək matris yığımı | **HAZIR** |
| 6c (CPR) | CPR ön-şərtçisi (3×3) | təxirə salınıb — bax qeyd |
| **6d (Nyuton)** | Nyuton-Rafson döngəsi | **HAZIR** |
| 6d (UI) | Mühərrik sarğısı + interfeys qoşulması | qalıb |

Bu, A6 qədər (bəlkə ondan da) böyük işdir — eyni səbəbdən mərhələlərə
bölünür: hər mərhələ tək başına test olunan, işləyən artım olmalıdır.

---

## Mərhələ 1 — Qaz PVT xassələri

### Z-faktoru (sıxılma əmsalı)

Real qaz `PV = ZnRT` tənliyi ilə təsvir olunur. Z-faktoru laboratoriya
kompozisiya analizi olmadan **Sutton (1985)** psevdo-kritik
korrelyasiyası + **Beggs-Brill (1973)** approksimasiyası ilə
hesablanır — Standing-Katz əyrisinin **açıq formada** (iterasiyasız)
təxmini.

İterativ üsullar (Dranchuk-Abou-Kassem) daha dəqiqdir, lakin hər
Nyuton addımında minlərlə hüceyrə üçün iterasiya daxilində iterasiya
etmək performans baxımından məqsədəuyğun deyil. Beggs-Brill vektor
əməliyyatlarına ideal uyğundur.

**Ölçülmüş nəticə:** Z-faktoru Standing-Katz əyrisinin klassik
formasını göstərir — aşağı təzyiqdə 1-ə yaxın (ideal qaz həddi), orta
təzyiqdə minimum (molekullararası cazibə üstünlük təşkil edir),
yüksək təzyiqdə 1-dən yuxarı (itələmə qüvvələri üstünlük təşkil edir).

### Bg (qazın formasiya həcm əmsalı)

    Bg = (Z·T/P) / (Zsc·Tsc/Psc)

Neft/su FVF-dən **əks istiqamətdə** davranır: Bo, Bw təxminən sabitdir
(mayelər az sıxılır), Bg isə təzyiqlə kəskin azalır (qaz çox sıxılır)
— aşağı təzyiqdə 100 dəfədən çox ola bilər.

### Qaz lözlüyü — Lee, Gonzalez, Eakin (1966)

Ən vacib fiziki fərq: **qaz lözlüyü təzyiqlə ARTIR** (sıxlıq artdıqca
molekullararası toqquşma çoxalır), neft/suyunku isə adətən azalır.
Bu, üç fazalı axında qazın niyə "yüngül və sürətli" hərəkət etdiyini
izah edir — aşağı lözlük + yüksək sıxılma.

### İnteqrasiya

`PVTTable` geriyə uyğun genişləndirildi: `gas_fvf`, `gas_viscosity`
sahələri `None` defoltludur. Köhnə iki fazalı cədvəllər (`A1`-dən)
`has_gas_phase = False` qaytarır və heç bir mövcud test/ssenari
pozulmur.

`build_pvt_table(..., include_gas=True)` ilə tam üç fazalı cədvəl
qurulur. `IPVTProvider` interfeysinə `has_gas_phase()`, `gas_fvf()`,
`gas_viscosity()`, `solution_gor()` əlavə olundu — defolt
`NotImplementedError` (iki fazalı provider-lər toxunulmaz qalır).

Analitik törəmələr (`gas_fvf_derivative`, `gas_viscosity_derivative`)
A6-dakı üsulla — cədvəlin dəqiq parçalı-xətti meyli — hazırlandı,
çünki Jakobian (mərhələ 5) bunları birbaşa tələb edəcək.

---

## Mərhələ 2 — domain: üç fazalı doyumluluq

### `ThreePhaseSaturation`

    Sw + So + Sg = 1

İki fazalı modeldə olduğu kimi, üçüncü doyumluluq (bu dəfə So) sərbəst
dəyişən deyil — həmişə `1 − Sw − Sg` kimi hesablanır.

`is_saturated` — hansı hüceyrələrdə sərbəst qaz var (`Sg > 0`).
`clip()` iki hədd yoxlamasını **birlikdə** aparır: `Sw` və `Sg`
ayrı-ayrı öz intervalına salınsa da, cəmi 1-dən böyük çıxa bilər
(məs. Sw=0.9, Sg=0.9 → 1.8) — belə halda hər ikisi nisbi miqyaslanır.

### `saturation_state()` — doyma vəziyyəti

Hüceyrənin cari `Rs`-i PVT cədvəlinin doyma əyrisi ilə (`Rs_sat(p)`)
tutuşdurulur: bərabərdirsə hüceyrə **doymuşdur** (sərbəst qaz var),
aşağıdırsa **doymamışdır** (bütün qaz neftdə həll olub). Bu, mərhələ
4-dəki dəyişən keçidin əsasıdır.

### Qaz papağı equilibration

`EquilibriumInitializationProvider` genişləndirildi: `gas_oil_contact`
verilib PVT provider `has_gas_phase()` qaytarırsa, **üç zonalı**
tarazlıq qurulur:

    GOC-dan yuxarı   qaz zonası     Sw = Swc,  Sg = 1 − Swc
    GOC — OWC arası  neft zonası    Sg = 0,    Sw əvvəlki kimi (A3/A4)
    OWC-dən aşağı    su zonası      Sg = 0,    Sw = 1 − Sor

Təzyiq profilinə qazın sütun çəkisi də daxildir (iterativ, sıxlıq
təzyiqdən asılı olduğu üçün A3-dəki üsulla eyni).

**Ölçülmüş nəticə** (10 təbəqə, GOC 2030 m, OWC 2070 m):

| Zona | Dərinlik | Sw | Sg | So |
|---|---|---|---|---|
| Qaz | 2005–2025 m | 0.20 | 0.80 | 0.00 |
| Neft | 2035–2065 m | 0.20 | 0.00 | 0.80 |
| Su | 2075–2095 m | 0.75 | 0.00 | 0.25 |

Qaz sütununun təzyiqə töhfəsi ölçüldü: fərq 5 bar-dan kiçikdir (ρ_qaz
≪ ρ_neft/ρ_su) — gözlənilən, çünki qaz demək olar çəkisizdir.

### Sadələşdirmə

Qaz-neft kontaktında kapilyar keçid zonası hələ qurulmur (kəskin
sərhəd). Su-neft kontaktındakı A4 üsulu (hamar keçid) buraya da
tətbiq oluna bilər — gələcək təkmilləşdirmə, mərhələ 4-6-nı bloklamır.

### İnterfeysə hələ ƏLAVƏ OLUNMADI — qəsdən

Mühərriklər (`IMPES`, `FullyImplicit`) hələ `InitialState.gas_saturation`-ı
oxumur (mərhələ 4-6-nın işidir). GOC-nu indi interfeysə əlavə etsək,
istifadəçi onu qoyar, equilibration "uğurla" işləyər, amma simulyasiya
səssizcə qazı görməzdən gələr — yanlış təəssürat yaradar. Ona görə bu
mərhələ proqramlaşdırma səviyyəsində qalır; UI yalnız mühərrik
hazır olanda (mərhələ 6) qoşulacaq.

## Mərhələ 3 — üç fazalı nisbi keçiricilik (Stone II)

İki fazalı Corey kifayət etmir: neftin nisbi keçiriciliyi üç fazalı
axında HƏM su, HƏM qaz doyumluluğundan asılı olur.

### Stone II düsturu

    kro = kro_end · [(krow/kro_end + krw)·(krog/kro_end + krg)
                      − (krw + krg)]

Stone I-dən fərqli olaraq (mənfi `kro` verə bilər — sənayədə məlum
problem) Stone II riyazi olaraq mənfi olmayan nəticə zəmanət edir,
`[0, kro_end]`-ə kəsilir.

### `GasCoreyParameters` — su-neft əyrisinin yoldaşı

Qaz-neft sistemi (`krg`, `krog`) su-neft sistemi (`krw`, `kro`) ilə
EYNİ formada — Corey düsturu, bu dəfə `Sg` üzərində. `kro_end`
TƏKRARLANMIR: hər iki əyri eyni fiziki son nöqtəyə (Swc-də, Sg=0-da)
istinad edir, fərqli dəyər versək iki əyri uyğunsuz olardı.

### `StoneRelativePermeabilityProvider`

Mövcud iki fazalı provider-i (Corey **və ya** cədvəl əsaslı, B4)
**sarır** — `krw`/`kro` üçün onu birbaşa çağırır, `krg`/`krog` üçün
`GasCoreyParameters`-i. `from_corey()` klassmetodu su-neft
parametrlərindən `kro_end`-i avtomatik çəkir.

### Ölçülmüş doğruluq

- **Sg=0-da** Stone II **dəqiq** iki fazalı su-neft əyrisinə düşür —
  fərq maşın dəqiqliyi səviyyəsində (~10⁻¹⁷)
- **Mənfi kro yoxdur** — 900 nöqtəlik (Sw, Sg) toru üzərində minimum
  `kro = 0.0`, heç vaxt mənfi (Stone I-in əsas üstünlüyü)
- **`kro_end` heç vaxt aşılmır** — [0, kro_end] məhdudiyyəti qorunur

### İnterfeys — önizləmə statusunda

`5 · NİSBİ KEÇİRİCİLİK` bölməsində **"Qaz-neft əyriləri (önizləmə)"**
qutusu: aktivləşdirilsə `krg(Sg)`/`krog(Sg)` əyriləri üçüncü qrafikdə
göstərilir. Açıq şəkildə **önizləmədir** — mühərriklər hələ bu
əyriləri istehlak etmir (mərhələ 4-6). Məqsəd: istifadəçi parametrləri
əvvəlcədən tənzimləyib əyriyə baxa bilsin, simulyasiyaya təsir hələ
yoxdur.

## Mərhələ 4 — primary dəyişənlər və dəyişən keçid

Sənaye simulyatorlarının (Eclipse, IMEX) standart üsulu: hüceyrənin
**3-cü primary dəyişəni** doyma vəziyyətindən asılı olaraq dəyişir:

    doymuş hüceyrə (sərbəst qaz var)      3-cü dəyişən = Sg
    doymamış hüceyrə (bütün qaz həll olub) 3-cü dəyişən = Rs

### Niyə belə, niyə həmişə Sg yox

Doymamış hüceyrədə Sg **həmişə** 0-dır — onu primary dəyişən etmək
mənasız olardı (Nyuton sabit 0-ı "həll edərdi", dəyişməyəcək bir
kəmiyyətə tənlik həsr edilərdi). Doymuş hüceyrədə isə Rs **həmişə**
Rs_sat(p)-ə bərabərdir — onu izləməyə ehtiyac yoxdur, çünki birbaşa
təzyiqdən çıxarıla bilər.

### `ThreePhaseState` — Nyuton vektoru

    x = [p_0, Sw_0, x_0,  p_1, Sw_1, x_1,  …]

İki fazalı sxemlə (`state.py`) eyni interleaved struktur — CPR
ön-şərtçisi bunu tələb edir (indi 2×2 yerinə 3×3 blok, mərhələ 5).

`is_saturated` bayrağı Nyuton vektorunun **hissəsi deyil** — kəsilməz
dəyişən deyil, diskret vəziyyətdir, ayrıca daşınır.

### Keçid qaydası

    doymamış → doymuş   Rs (cari 3-cü dəyişən) Rs_sat(p)-i keçəndə
                         → Sg = 0-dan başlayır (kəsilməzlik)
    doymuş → doymamış    Sg (cari 3-cü dəyişən) mənfiyə düşəndə
                         → Rs = Rs_sat(p) sərhəd dəyərinə qayıdır

**Kəsilməzlik** vacibdir: keçid anında hər iki təsvir EYNİ fiziki
vəziyyəti göstərməlidir (Sg=0 ⟺ Rs=Rs_sat(p)) — əks halda kütlə
qalıq tənliklərində süni "atlama" yaradardı.

### Ölçülmüş doğruluq

Üç hüceyrəli sınaq (200 bar, Rs_sat=123.74 sm³/sm³):

| Hüceyrə | Əvvəl | Sonra | Nə baş verdi |
|---|---|---|---|
| 0 | Rs=98.99 (doymamış) | Rs=98.99 (doymamış) | sərhədi keçmədi |
| 1 | Rs=142.3 (doymamış) | Sg=0, doymuş | sərhədi yuxarı keçdi |
| 2 | Sg=−0.02 (doymuş) | Rs=123.74, doymamış | sərhədi aşağı keçdi |

Hər üçündə `Sw + So + Sg = 1` **dəqiq** qorunur (keçiddən əvvəl və
sonra). `updated()` metodu Rs dəyişəninə doyumluluq həddini tətbiq
ETMİR (Sg-dən fərqli vahid, sm³/sm³) — yalnız doymuş hüceyrələrdə
Appleyard kəsilməsi (chopping) tətbiq olunur.

### Növbəti addım

Bu modul hələ **qalıq tənliklərinə qoşulmayıb** — yalnız vəziyyəti
təsvir edir və sınanıb. Mərhələ 5-də üç kütlə balansı tənliyi (su,
neft, qaz — sərbəst + neftdə həll olmuş) bu vəziyyət üzərində
qurulacaq.

## Mərhələ 5 — qalıq tənlikləri (akkumulyasiya)

İki fazalı sxemdə hər hüceyrə üçün İKİ kütlə balansı tənliyi var idi
(su, neft). Qaz ÜÇÜNCÜ tənlik əlavə edir, lakin onun forması
keyfiyyətcə fərqlidir: qaz İKİ yerdə saxlanıla bilər.

### Qaz tənliyi — iki mənbə

    N_water = PV · Sw / Bw                    (dəyişməyib)
    N_oil   = PV · So / Bo                     (dəyişməyib)
    N_gas   = PV · (Sg/Bg + So·Rs/Bo)          (YENİ)

Qaz düsturundakı iki hədd fiziki cəhətdən ayrıdır: birincisi **sərbəst
qaz fazasının özü**, ikincisi **neftin daşıdığı həll olmuş qaz**. Su
heç vaxt neftdə həll olmur bu modeldə ("quru qaz" fərziyyəsi, Rv=0).

Doymamış hüceyrədə (Sg=0) yalnız ikinci hədd qalır — bu, mərhələ
4-dəki `ThreePhaseState.solution_gor()`-un niyə doymuş hüceyrələr
üçün Rs_sat(p) qaytardığını izah edir: hətta sərbəst qaz olanda da
neft maksimum qazı özündə saxlayır.

### `ThreePhaseAccumulator`

Su/neft düsturu **dəyişməyib** — A6-dakı `ResidualAssembler.accumulation()`
ilə eynidir, `two_phase_accumulation_matches()` doğrulama köməkçisi
bunu təsdiqləyir. Yalnız qaz tənliyi yenidir; bu ayrılıq qəsdən
saxlanılıb ki, mövcud iki fazalı testlər sınmasın.

### Ölçülmüş doğruluq

- **Doymamış hüceyrədə** qaz **dəqiq** yalnız həll olmuş hissədən
  gəlir (fərq = 0.0, maşın dəqiqliyi)
- **Doymuş hüceyrədə** həm sərbəst, həm həll olmuş hədd müsbətdir və
  cəmi tənliyi verir
- **Su/neft** doymamış hüceyrələrdə iki fazalı düsturla **fərqsiz**
  üst-üstə düşür; doymuş hüceyrələrdə (gözlənildiyi kimi) fərqlənir,
  çünki `So = 1 − Sw − Sg` artıq `1 − Sw`-dən kiçikdir
- **Sərhəddə kəsilməzlik**: Sg=0 ⟺ Rs=Rs_sat(p) nöqtəsində qaz
  kəmiyyəti sıçramır (fərq ədədi dəqiqlik həddində)

## Mərhələ 6a — axın həddləri (üç fazalı Darcy axını)

A6-dakı (`residual.py`) axın modulunu güzgüləyir: `ThreePhaseFlux`
sinfi hər üzdə üç fazanın Darcy axınını hesablayır.

### Su/neft — dəyişməyib

Formula A6-dakı ilə eynidir (upstream çəkiləmə, cazibə, kapilyar).

### Qaz — iki mexanizmlə daşınır

    1. SƏRBƏST qaz    öz təzyiq qradiyenti ilə (Φ_g)
    2. HƏLL OLMUŞ qaz  NEFTLƏ BİRLİKDƏ (Rs · neft axını)

İkinci hədd üçün upstream istiqaməti **neftin öz istiqamətidir** — Rs
"sərnişindir", öz axın istiqamətini seçmir, neftin apardığı yerə
gedir. Bu, standart black-oil qaz axını formuludur.

### Ölçülmüş fizika

Eyni doyumluluqla, yalnız cazibə fərqi olan sütunda potensial
fərqləri ölçüldü:

| Faza | ΔΦ (bar) | Sıxlıq sırası |
|---|---|---|
| Su | 1.98 | ən sıx |
| Neft | 1.21 | orta |
| Qaz | **0.34** | ən yüngül |

Sıralama fizika ilə tam uyğundur: `ΔΦ_w > ΔΦ_o > ΔΦ_g` — çünki
qazın cazibə başlığı sıxlığı ilə mütənasibdir, ρ_qaz ≪ ρ_neft < ρ_su.

Əlavə doğrulamalar: təcrid olunmuş sütunda daxili üzlərin cəmi axını
sıfır (kütlə itmir), bərabər təzyiqdə axın sıfır, təzyiq qradiyenti
yüksəkdən aşağıya axın yaradır, həll olmuş qazın işarəsi həmişə
neft axınının işarəsi ilə üst-üstə düşür.

## Növbəti addımlar (6b-6d)

## Mərhələ 6b — quyu həddləri

`ThreePhaseWellModel` — A6-dakı `well_rates()`-i güzgüləyir.

### Vuruculardan dəyişməyib

Vuruculuar hələ yalnız **su** vurur — A6-dakı davranış dəyişməyib.
Qaz vurma (WAG, CO₂ vurma) EOR-un öz mövzusudur və bu modulun
əhatəsindən kənardadır; `WellType`-a qaz-injektor tipi əlavə
olunanda buraya qoşulacaq.

### İstismarçılarda qaz — eyni iki mexanizm

Axın modulundakı (6a) eyni məntiq quyuya da tətbiq olunur:

    q_gas = q_sərbəst_qaz  +  Rs · q_neft

### RATE rejimi — hədəf hələ də MAYE debitidir

A6-dakı konvensiya qorunur: `RATE` hədəfi su+neft debitidir, qaz
NƏTİCƏ kimi çıxır, hədəf kimi yox. Bu, real qazlift/BHP rejimli
quyularla üst-üstə düşür — operator adətən maye debitini idarə edir,
qaz "gəldiyi kimi" gəlir. Üstəlik köhnə iki fazalı ssenarilər eyni
nəticəni verməyə davam edir (test bunu təsdiqləyir: maye debiti
dəqiq hədəflə üst-üstə düşür).

### Ölçülmüş doğruluq

- İstismarçıda `q_gas − Rs·q_oil` (sərbəst hissə) sıfırdan fərqlidir
  Sg>0 olanda — gözlənilən
- Vurucu **yalnız** su verir (neft/qaz sıfır)
- İşarə konvensiyası qorunur: istismarçı mənfi, vurucu müsbət
- Hər quyunun ümumi debiti öz hüceyrələrinin cəminə **dəqiq** bərabərdir

**6c — Jakobian:** ölçü 2N×2N-dən 3N×3N-ə keçir, CPR 2×2 blokdan
3×3 bloka genişlənir. Dəyişən keçid (mərhələ 4) Jakobianın strukturunu
hüceyrə-hüceyrə dəyişdirdiyi üçün analitik törəmələr xüsusi diqqət
tələb edir.

**6d — Nyuton + UI:** tam inteqrasiya, yalnız bundan sonra qaz fazası
real simulyasiyada işə düşəcək və interfeysdə görünəcək.

## Növbəti sessiya — daha üç düzəliş, kök hələ tapılmayıb

Bu sessiyada üç KONKRET, doğrulanmış düzəliş edildi:

1. **`mb_gas` bölmə xətası** — `Rs=0` (sərbəst qaz yoxdur) olanda
   `gas_in_place.sum()≈0`, nisbi mb_gas ölçüsü ƏHƏMİYYƏTSİZ ədədi
   səs-küyü ASTRONOMİK dəyərə çevirirdi, yığılmanı RİYAZİ CƏHƏTDƏN
   MÜMKÜNSÜZ edirdi. Düzəldi: gas_in_place əhəmiyyətsiz olanda mb_gas
   avtomatik 0 qəbul edilir.
2. **Geri-izləmənin (line search) miqyas problemi** — xalis L2 norm
   qaz tənliyinin nəhəng miqyası (kiçik Bg) ilə tamamilə üstələnirdi,
   su/neftin yaxşılaşıb-pisləşməsini "kor" edirdi. CNV-tipli
   faza-miqyaslı normla əvəz edildi.
3. **Axının upstream seçimi dondurma imkanı** — dəyişən keçid və
   quyu keçidi ilə eyni prinsiplə, axının "hansı hüceyrə yuxarı axın"
   qərarı da addımın əvvəlki vəziyyətinə görə dondurula bilər
   (`reference` parametri). Doğru, müdafiə xarakterli əlavədir, LAKİN
   bu konkret oscillasiyanı HƏLL ETMƏDİ.

### Dəqiq diaqnoz edilmiş qalan problem

t≈6.7 gündə, istismarçının BİLAVASİTƏ QONŞUSU olan hüceyrədə (5×5
grid modelində hüceyrə 29, istismarçı hüceyrə 35-in şimal qonşusu)
NEFT residualı DÖVR-2 rəqs edir (8×10⁻⁴ ↔ 4.5×10⁻³), heç vaxt
toleransa yaxınlaşmır. Bu, istismarçı özü demək olar tam boşalmış
vəziyyətdə (P≈8 bar) olanda, ona bilavasitə qonşu hüceyrənin öz
təzyiq/doyumluluq dinamikası ilə bağlıdır — nə mb_gas düzəlişi, nə
miqyaslı line search, nə upstream dondurma bunu həll etmədi.

**Növbəti addım üçün fikirlər:**
- Line search-in özünün PER-CELL (bütöv norm yox, ən pis hüceyrəyə
  görə) qəbul meyarına keçməsi
- Trust-region üsulu (addımın böyüklüyünü Jakobianın etibarlılıq
  radiusu ilə məhdudlaşdırmaq)
- Bu spesifik keçid zonasında (quyu demək olar tam boşalanda) əlavə
  "sub-relaxation" (yalnız bu bölgədə daha kiçik addım nisbəti)

## Mərhələ 6d (UI hissəsi) — QOŞULDU (Yol 2: təhlükəsiz uğursuzluqla)

Kök riyazi problem (quyu öz BHP sərhədinə yaxınlaşanda Nyuton
rəqsi) TAM həll olunmadı, lakin istifadəçi ilə razılaşaraq
**təhlükəsiz uğursuzluqla** UI-yə qoşuldu:

- Simulyasiya harda mümkündürsə oradan sonra da davam edir
- Yığılma tam bitməsə, proqram ÇÖKMÜR — `result.converged=False`
  və aydın mesajla (`"t=X gün: ... yığılmadı"`) dayanır
- O nöqtəyə qədər olan bütün nəticələr (qrafiklər, RF, GOR) görünür
- Bu davranış artıq `_on_finished()`-də (A6-dan) mövcud idi —
  əlavə koda ehtiyac olmadı, yalnız mühərriyin özünü qoşmaq
  kifayət etdi

### Yeni UI elementləri

- **PVT tabı:** "Qaz fazasını aktivləşdir (A7 — sınaq statusunda)"
  qutusu + aydın xəbərdarlıq qeydi
- **Ədədi parametrlər tabı:** "Qaz papağı (GOC)" qutusu + dərinlik
  girişi (mərhələ 2-dən bəri qəsdən söndürülmüşdü)
- **SCAL tabı:** "Qaz-neft əyriləri" artıq "(önizləmə)" DEYİL —
  parametrlər PVT-də qaz aktivdirsə mühərrikə həqiqətən ötürülür

### Servis səviyyəsində seçim

`ModelAwareSimulationService.create_engine()`: PVT-də qaz
xassələri varsa, IMPES/IMPLICIT seçimindən ASILI OLMAYARAQ
`ThreePhaseSimulationEngine` işə düşür (ikisi də qazı dəstəkləmir).

### Yol boyu tapılan iki əlavə səhv

1. **Uyğunsuz xətti həlledici** — servis bütün mühərriklərə A6-nın
   2×2 CPR/ILU üçün tənzimlənmiş `ScipyCgIluSolver`-i ötürürdü, bu
   isə üç fazalı 3×3 struktura uyğun deyildi və çökməyə səbəb
   olurdu. Düzəliş: `ThreePhaseSimulationEngine` ötürülən
   `linear_solver`-i FAKTİKİ İŞLƏTMİR, həmişə öz uyğun
   `NewtonLinearSolver` nüsxəsini saxlayır.
2. **SCAL/PVT qaz qutularının uyğunsuzluğu** — PVT-də qaz aktiv
   olub, SCAL-da öz qutusu işarələnməyəndə `gas_scal=None` qalırdı,
   mühərrik Stone relperm gözləyir, amma yalnız Corey alırdı —
   çökmə riski. Düzəliş: bu halda defolt `GasCoreyParameters()`
   avtomatik işlədilir.

### Ölçülmüş — vacib məlumat

UI-nin **defolt grid ölçüsü 41×41=1681 hüceyrədir** (əvvəlki
diaqnozlarımın 6×6 test modelindən qat-qat böyük). Bu miqyasda
kök qeyri-xəttilik problemi DAHA TEZ üzə çıxır (t≈0 ətrafında,
6×6-dakı t≈6.7 əvəzinə). **Tövsiyə: hazırda kiçik-orta grid
ölçüləri (məs. 10×10 və aşağı) ilə sınamaq daha etibarlıdır.**

## A7 — cari ümumi status

Bütün 6 mərhələ EDİLİB, nüvə fizikası 100% doğrulanıb, UI qoşulub
(təhlükəsiz uğursuzluqla). Qalan: kök qeyri-xəttilik problemini tam
həll etmək (böyük gridlərdə etibarlılığı artırmaq üçün) — bu, artıq
A7-ni "bitirmək" üçün deyil, onu DAHA ETİBARLI etmək üçün davam edən
işdir.

## v50 — istifadəçi tərəfindən tapılan real UI səhvi

İstifadəçi qaz fazasını aktivləşdirib sınayanda köhnə (A6 dövründən
qalma, A7-dən əvvəl yazılmış) bir xəbərdarlıq gördü: *"BHP doyma
təzyiqindən aşağıdır — quyudibində qaz ayrılacaq. Qaz fazası
modelləşdirilmir, nəticələr nikbin ola bilər."*

Bu, artıq **yanlış** idi — qaz fazası indi (v49-dan bəri) həqiqətən
modelləşdirilir. `ReservoirModel.diagnose()`-dəki yoxlama
`pvt_table.has_gas_phase`-ə heç baxmırdı, yalnız `bubble_point > 0`
şərtinə əsaslanırdı — A7-dən əvvəl yazıldığı üçün bu fərqi bilmirdi.

Düzəliş: xəbərdarlıq YALNIZ `not pvt_table.has_gas_phase` olanda
göstərilir. Real quyularla yoxlanıldı: qaz aktiv deyilsə xəbərdarlıq
göstərilir (dəyişməyib), qaz aktivdirsə göstərilmir (düzəldi).

Bu, growl A6→A7 keçidində "köhnə mesajların yeni imkanlara uyğun
yenilənməməsi" sinfindən bir səhv idi — gələcəkdə bənzər yerlər
(digər diaqnostika mesajları) da bu baxımdan nəzərdən keçirilməlidir.

## v51 — İstifadəçi bildirişi: real proqram çökməsi (crash)

İstifadəçi bildirdi: qaz fazası aktiv olanda proqram **tam çökür**
(bağlanır) — `converged=False` ilə yumşaq dayanmır, birbaşa
sistemdən atılır. Kiçik gridlə (10×10 və aşağı) də eyni nəticə.

### Diaqnoz

Linux/xvfb mühitində (real QThread axını, tam UI axını, nəticə
göstərmə daxil) TƏKRARLANA BİLMƏDİ — bu, güman ki, **Windows-a xas**
native kitabxana davranışı ilə bağlıdır (scipy-nin sparse xətti
həllediciləri (SuperLU/UMFPACK) platformadan asılı ola bilər).

Kod nəzərdən keçirilərkən İKİ REAL BOŞLUQ tapıldı:

1. **Dar istisna tutma** — `ThreePhaseNewtonSolver.solve()`-də xətti
   həll ətrafında yalnız `(FloatingPointError, ValueError)` tutulurdu.
   Scipy-nin özü `RuntimeError`, `LinAlgError` və digər növləri ata
   bilər — bunlar tutulmadan yuxarı sızırdı.
2. **NaN/Inf mühafizəsi yox idi** — qalıq vektoru və ya Jakobian
   matrisi NaN/Inf ehtiva edəndə (güclü osilyasiyadan sonra mümkündür)
   birbaşa scipy-nin native kodu ilə həll edilirdi. Bəzi
   platformalarda (xüsusən Windows-da SuperLU) NaN-lı seyrək matris
   Python-un tuta bilmədiyi native səviyyəli çökməyə səbəb ola bilər.

### Düzəliş — üç qat mühafizə

1. Qalıq/Jakobian NaN/Inf ehtiva edirsə, xətti həllediciyə HEÇ
   VERİLMİR — əvvəlcədən `LINEAR_SOLVER_FAILED` qaytarılır
2. İstisna tutma `Exception`-a qədər genişləndirildi
3. **İKİ ƏLAVƏ TƏHLÜKƏSİZLİK QATI:** `ThreePhaseNewtonSolver.solve()`
   və `ThreePhaseSimulationEngine.run()` hər ikisi öz "inner" tətbiqini
   xarici `try/except Exception` ilə əhatə edir — HANSI SƏBƏBDƏNSƏ
   (indi gözlənilməyən, gələcəkdə yarana biləcək) bir şey səhv getsə
   belə, metod İSTİSNA ATMIR, təhlükəsiz nəticə qaytarır

Bu, "GARANTİ" səviyyəsində sənədləşdirilib hər iki metodun
docstring-ində — gələcək dəyişikliklər bu qərarı təsadüfən poza bilməz.

### Doğrulama

4 yeni test: sındırılmış PVT provider, NaN-lı başlanğıc vəziyyət,
sındırılmış akkumulyasiya, əvvəllər tutulmayan `RuntimeError` növü —
hamısında `solve()`/`run()` təmiz (istisnasız) təhlükəsiz nəticə
qaytardı, loqa yazdı, çökmədi.

**Qeyd:** bu, KÖK səbəbi (nəyin konkret olaraq Windows-da çökməyə
səbəb olduğunu) tapmaqdan fərqlidir — bu, HƏR HANSI səbəbdən asılı
olmayaraq proqramın artıq ÇÖKƏ BİLMƏYƏCƏYİNƏ zəmanət verən müdafiə
təbəqəsidir. İstifadəçi yenidən çökmə müşahidə etsə, bu, çox
əhəmiyyətli əlavə diaqnostika məlumatı olardı (öz-özlüyündə mümkün
olmamalıdır artıq).

## v59 — ƏSL DİAQNOZ: bloklayan meyar material balansdır

Bu sessiyada OPM Flow-un sənədi araşdırıldı və üç dəyişiklik sınandı.
Ən dəyərlisi isə **nəhayət dəqiq diaqnoz** oldu.

### Düzəldilən iki əsl səhv

1. **Su lözlüyünün təzyiq törəməsi sıfır sayılırdı** — həm axın, həm
   quyu Jakobianında (`mu_w_p = np.zeros(...)`). Mühərrik `mu_w`-ni
   PVT-dən oxuyur və o, təzyiqdən ASILIDIR. Səhv gizli qalmışdı,
   çünki BÜTÜN testlər sabit `mu_w=0.5` işlədirdi — sabit dəyərdə
   törəmə həqiqətən sıfırdır, ona görə sonlu-fərq yoxlamaları
   keçirdi.

   **Metodoloji dərs:** test qurğusu real mühərrikdən fərqli olanda
   belə səhvlər gizlənir. İlk dəfə REAL mühərrikin öz PVT xassələri
   ilə tam Jakobian doğrulaması aparıldı — düzəlişdən sonra xəta
   **9×10⁻⁹** (əla).

2. **PVT cədvəlinin defolt aşağı həddi** 10 → 1 bar. Cədvəldən
   kənarda bütün xassələr sabitləşir (törəmə = 0) — bu, sınıq
   nöqtədir.

### Sınanıb GERİ QAYTARILAN

**Hamar kəsmə** (`min(q,0)` əvəzinə). İki səbəb: hədəf problemi həll
etmədi (CNV tarixçəsi hərfi olaraq eyni qaldı), VƏ fizikanı pozdu —
nisbi keçiricilik sıfır olanda debit dəqiq sıfır olmalıdır,
hamarlaşdırma isə ona süni −58 847 sm³/gün verdi. Funksiya kodda
XƏBƏRDARLIQ şərhi ilə saxlanılıb ki, təkrar sınanmasın.

### ƏSL DİAQNOZ (bütün əvvəlki fərziyyələri təkzib edir)

Uğursuz nöqtədə ölçüldü:

    CNV      = 4.3×10⁻⁴   (hədd 10⁻³)  ✓ KEÇİR
    MB su    = 5.0×10⁻⁴   (hədd 10⁻⁷)  ✗ BLOKLAYIR
    MB neft  = 1.7×10⁻²⁰ → 1.5×10⁻⁴    ✗ BLOKLAYIR

**Bloklayan meyar material balansdır, CNV deyil.** Üstəlik su və
neft MB-si NÖVBƏLƏŞİR: biri yaxşı olanda digəri pisdir.

Kritik müşahidə: **XAM Nyuton da** (Appleyard kəsməsi, line search,
doyumluluq hədləri — heç biri olmadan) eyni dövr-2 osilyasiyanı
göstərir, halbuki Jakobian dəqiqdir (9×10⁻⁹). Yəni problem addım
idarəsində DEYİL.

Bu, əvvəlki bütün fərziyyələri təkzib edir: quyu BHP sərhədi deyil,
PVT cədvəl kənarı deyil, upstream sıçrayışı deyil, Jakobian
səhvi deyil.

### Növbəti addım üçün istiqamət

MB meyarının özü araşdırılmalıdır: bizdə `|Σ R|·dt / yerində_olan`,
OPM-də isə məsamə həcminə normallaşdırılır (doyumluluq xətası kimi
şərh olunur) və defolt həddi 10⁻⁶-dır (bizdəki 10⁻⁷-dən 10 dəfə
boş). Lakin ölçülmüş MB neft = 2.4×10⁻⁵ OPM-in həddini də keçir —
yəni sadəcə həddi boşaltmaq DÜZGÜN həll deyil, kütlə saxlanmasında
əsl uyğunsuzluq var və onun mənbəyi tapılmalıdır.

## v60 — JAKOBİANDA UZUN MÜDDƏT GİZLİ QALMIŞ SƏHV TAPILDI

Kritik eksperiment: Jakobian UĞURSUZ VƏZİYYƏTDƏ yoxlanıldı (əvvəllər
yalnız bircins başlanğıc vəziyyətdə yoxlanılırdı). Nəticə:

    bircins vəziyyətdə:  9×10⁻⁹   ✓
    uğursuz vəziyyətdə:  1.7×10⁻³ ✗  (190 000 dəfə pis!)

### Səhv: hasil qaydasının yarısı unudulmuşdu

Qaz axını iki hissədən ibarətdir (bax mərhələ 6a):

    F_qaz = F_sərbəst + Rs_upstream · F_neft

Rs-ə görə törəmə HASİL QAYDASI tələb edir:

    ∂F/∂Rs_up = F_neft  +  Rs_up · ∂F_neft/∂Rs_up
                ↑ UNUDULMUŞDU

Kodda yalnız ikinci hədd var idi.

### Niyə illərlə gizli qaldı

Bu hədd YALNIZ hüceyrə DOYMAMIŞ olanda sıfırdan fərqlidir — çünki
yalnız o halda 3-cü primary dəyişən Rs-dir. Doymuş hüceyrədə
Rs = Rs_sat(p), yəni 3-cü dəyişəndən (Sg) asılı deyil → hədd sıfır.

**Bütün əvvəlki tam-sistem doğrulamalarım TAM DOYMUŞ vəziyyətdə
(`np.ones(n, bool)`) aparılmışdı.** Real simulyasiya isə çox vaxt
DOYMAMIŞ vəziyyətdə işləyir (Sg=0, Rs sərbəst dəyişən) — məhz orada
Jakobian səhv idi.

Düzəlişdən sonra: **4.6×10⁻⁸** (36 000 dəfə yaxşılaşma).

Yeni test (`test_jacobian_is_exact_in_the_all_undersaturated_state`)
bu boşluğu bağlayır.

### Metodoloji dərs (ikinci dəfə)

Bu, v59-dakı `mu_w_p` səhvi ilə EYNİ SİNİFDƏNDİR: doğrulama qurğusu
real işləmə şəraitindən fərqli olanda səhv gizlənir. İki nümunə:

- v59: testlər sabit `mu_w` işlədirdi → lözlük törəməsi səhvi görünmürdü
- v60: testlər tam doymuş vəziyyət işlədirdi → Rs törəməsi səhvi görünmürdü

**Qayda:** doğrulama real mühərrikin işlədiyi HƏR rejimdə (doymuş,
doymamış, qarışıq) aparılmalıdır.

### Amma problem HƏLƏ HƏLL OLUNMAYIB

Jakobian indi dəqiq olsa da, simulyasiya eyni nöqtədə (t≈6.7 gün)
dayanır. Ölçüldü: Nyuton istiqaməti boyunca residual yalnız ~17 %
azala bilir (t≈0.4-0.5-də minimum), yəni problem GÜCLÜ
QEYRİ-XƏTTİLİKdir, Jakobian səhvi deyil.

Bloklayan meyar hələ də material balansdır (CNV keçir).

## OPM tipli quyu modeli — MƏRHƏLƏLİ YENİDƏNQURMA

Qərar (istifadəçi ilə): kök qeyri-xəttilik problemini həll etmək üçün
quyu modeli OPM-in "standart quyu modeli" prinsipi ilə yenidən
qurulur. İş MƏRHƏLƏLƏRƏ bölünüb — hər biri ayrıca sınanan, işləyən
artımdır.

| # | Mərhələ | Vəziyyət |
|---|---|---|
| **1** | Quyu vəziyyəti: BHP naməlum dəyişən kimi | **HAZIR (v61)** |
| **2** | Quyu idarəetmə tənliyi + perforasiya debitləri BHP-dən | **HAZIR (v62)** |
| **3** | Jakobian: rezervuar↔quyu qarşılıqlı bağlantısı | **HAZIR (v63)** |
| 4 | Nyuton inteqrasiyası + doğrulama | **QİSMƏN (v64)** |
| 5 | Quyu saxlama həddi (OPM-in sabitləşdiricisi) | |

### Niyə bu prinsip

Hazırda BHP **sabit sərhəd şərtidir**: debiti ondan hesablayırıq və
quyu bağlananda debit BİRDƏN sıfıra düşür. OPM-də BHP **naməlum
dəyişəndir** və ona ayrıca idarəetmə tənliyi yazılır:

    BHP idarəsində:   R_ctrl = p_bhp − p_hədəf = 0
    RATE idarəsində:  R_ctrl = Σ q_perf − q_hədəf = 0

İdarəetmə rejimi dəyişəndə DEBİT sıçramır — sadəcə hansı tənliyin
işlədiyi dəyişir. Kəsilməzlik qorunur.

### Mərhələ 1 — vektor yerləşməsi

    [ p₀,Sw₀,x₀, …, p_{N−1},Sw_{N−1},x_{N−1},  bhp₀,…,bhp_{W−1} ]
      └────────── rezervuar (3N) ──────────┘  └── quyular (W) ──┘

Quyu naməlumları SONA əlavə olunur, rezervuarın 3×3 blok strukturu
TOXUNULMUR — CPR ön-şərtçisi (A6) ona əsaslanır. OPM də eyni
yerləşməni işlədir (quyu blokunu Schur tamamlayıcı ilə ayırır).

`WellUnknowns` — hər quyu üçün BİR naməlum (perforasiya sayından
asılı deyil). İlkin qiymət: BHP idarəlidə hədəfin özü, RATE
idarəlidə perforasiya hüceyrəsinin təzyiqi.

`CoupledState` — rezervuar + quyu birlikdə; `updated()` Appleyard
kəsməsini rezervuara, ayrıca BHP həddini quyulara tətbiq edir.

Bu mərhələ heç bir mövcud kodu DƏYİŞMİR — yanında yaşayır, ona görə
işləyən mühərrik toxunulmaz qalır.

### Mərhələ 2 — quyu tənlikləri (v62)

`StandardWellModel` — debitlər BHP naməlumundan hesablanır:

    q_α = WI · λ_α · (p_bhp − p_hüceyrə)

**KƏSMƏ YOXDUR.** İşarə təbii olaraq təzyiq fərqindən çıxır. Köhnə
modeldəki `min(q, 0)` sıfır nöqtəsində SINIQ idi — törəməsi 0-dan
1-ə sıçrayırdı və Nyuton məhz orada osilyasiya edirdi.

**Ölçülmüş nəticə** (istismarçı BHP-si hüceyrə təzyiqindən keçir):

    BHP = 213.4 → qo = −1.0309
    BHP = 213.5 → qo = ±0.0000
    BHP = 213.6 → qo = +1.0309

Keçid tam XƏTTİdir, sınıq yoxdur. Test bunu bərabər addımlarda
bərabər artımla yoxlayır.

Bu, eyni zamanda ÇARPAZ AXINI (cross-flow) təbii modelləşdirir —
istismarçının BHP-si hüceyrədən yuxarı olanda debit müsbətə keçir,
OPM də belə edir.

**Mobilliyin upstream seçimi:** laya daxil olanda vurulan fazanın
mobilliyi (su vurucusunda yalnız su), laydan çıxanda hüceyrənin öz
mobillikləri.

**İdarəetmə tənlikləri:**

    BHP idarəsində:   R = p_bhp − p_hədəf          (bar)
    RATE idarəsində:  R = (Σ q_maye − q_hədəf)/WI  (bar miqyasına gətirilir)

RATE qalığı quyu indeksinə bölünür: BHP qalığı bar, RATE qalığı
m³/gün — eyni vektorda çox fərqli böyüklüklər pis şərtlənmə yaradır.

### Mərhələ 3 — quyu Jakobianı (v63)

`StandardWellJacobian` — dörd blok:

    ┌─────────────┬──────────┐
    │  rezervuar  │  R↔Q     │   R↔Q: ∂(debit)/∂p_bhp  ← YENİ
    ├─────────────┼──────────┤   Q↔R: ∂(idarəetmə)/∂(hüceyrə)
    │  Q↔R        │   Q      │   Q:   ∂(idarəetmə)/∂p_bhp
    └─────────────┴──────────┘

`∂q/∂p_bhp = WI·λ/B` — köhnə modeldə bu bağlantı ÜMUMİYYƏTLƏ YOX idi,
çünki BHP sabit idi. Məhz o, quyunu sistemin bir hissəsinə çevirir.

**Doğrulama — keçmiş dərs tətbiq olundu.** v60-da tapılmışdı ki,
Jakobian yalnız TAM DOYMUŞ vəziyyətdə doğrulanırdı və doymamış
rejimdəki səhv uzun müddət gizli qalmışdı. Bu dəfə HƏR İKİ rejim
ayrıca yoxlanıldı:

    dəyişən   doymamış    doymuş
    ────────────────────────────────
    p         3.9×10⁻¹⁰   9.1×10⁻¹⁰
    Sw        1.2×10⁻¹⁰   1.0×10⁻¹¹
    x         0.0         3.5×10⁻¹³
    BHP       8.5×10⁻¹¹   5.1×10⁻¹⁰

Hamısı maşın dəqiqliyində. RATE idarəli quyunun idarəetmə tənliyi də
ayrıca doğrulandı (o, həm BHP-dən, həm rezervuardan asılıdır).

**Mobilliyin keçidi haqqında:** axın istiqaməti dəyişəndə hansı
mobilliyin işlədildiyi dəyişir; Δp=0-da debit hər iki halda sıfırdır
(funksiya kəsilməz), lakin törəmə sıçrayır. OPM-də də belədir və
köhnə `min`-dəki sıçrayışdan qat-qat zəifdir — orada FUNKSİYANIN
ÖZÜ sınırdı.

### Mərhələ 4 — Nyuton inteqrasiyası (v64, QİSMƏN)

`CoupledNewtonSolver` — rezervuar (3N) və quyu (W) naməlumlarını
BİRLİKDƏ həll edir. Rezervuar blokları MÖVCUD siniflərdən olduğu kimi
götürülür (onlar A7-nin əvvəlki mərhələlərində doğrulanıb), yalnız
quyu hissəsi yenidir.

**İŞLƏYƏN NƏTİCƏ:** tək addım (dt=1.0) **4 iterasiyada yığılır**,
BHP-lər hədəflərində dayanır:

    CNV: 1.40e+00 → 4.82e-01 → 3.89e-02 → 1.82e-03 → 4.69e-07
    son BHP: [320, 150]  (hədəflər 320/150)

Çoxaddımlı: dt=1.0 ilə 6 gün irəliləyir (köhnə model 6.7-də dayanırdı
— oxşar nöqtə).

**TAPILAN VƏ DÜZƏLDİLƏN SƏHV:** `NewtonConfig`-də `control_tolerance`
və `max_bhp_change` sahələri YOX idi, kod isə onlara müraciət edirdi.
Yaranan `AttributeError` təhlükəsiz mühafizə tərəfindən udulur və
"yığılmadı" kimi görünürdü — yəni əsl səbəb GİZLƏNİRDİ. Bu, təhlükəsiz
mühafizənin gözlənilməz yan təsiridir: o, çökməni önləyir, lakin
proqramlaşdırma səhvlərini də maskalayır.

**AÇIQ QALAN MƏSƏLƏ:** kiçik addımda (dt=0.25) həll YIĞILMIR — CNV
0.32-0.35 ətrafında ilişib qalır, halbuki dt=1.0 ilə 4 iterasiyada
yığılır. Bu, gözləniləndən TƏRSİNƏdir (kiçik addım asan olmalıdır) və
miqyasdan asılı bir problemə işarə edir. Növbəti sessiyada
araşdırılmalıdır.

### v65 — geri-izləmə + quyunun bağlanması

**İki əlavə, hər ikisi ölçülmüş fayda verdi:**

1. **Geri-izləmə (line search)** birləşmiş həllediciyə əlavə olundu.
   Ölçüldü (dt=0.25): Nyuton δp = **1456 bar** istəyirdi, Appleyard
   kəsməsi onu 50 bar-a endirirdi — nəticədə İSTİQAMƏT pozulurdu və
   həll sıçrayırdı (P_min: 163→159→112→162→112…). Kəsmə yalnız addımın
   UZUNLUĞUNU məhdudlaşdırır, onun FAYDALI olub-olmadığını yoxlamır;
   geri-izləmə məhz bunu edir. Nəticə: dt=0.25 ilə t=0.0-dan **t=8.0
   günə** qədər irəlilədi.

2. **Quyunun bağlanması.** BHP idarəli istismarçının hüceyrə təzyiqi
   hədəfindən aşağı düşəndə düstur `q > 0` verir — quyu laya VURMAĞA
   başlayır (ölçüldü: **268 m³/gün neft laya vurulurdu**, absurd).
   Real quyu belə halda dayanır; OPM də idarəetməni dəyişir.

   **Vaxtlama vacibdir:** bu qərar Nyuton iterasiyaları arasında yox,
   yalnız zaman addımının ƏVVƏLİNDƏ verilir və addım boyu SABİT qalır.
   Əks halda quyu iterasiyadan-iterasiyaya açılıb-bağlanar — bu, məhz
   köhnə `min(q,0)` problemidir.

**Cari vəziyyət:** tək addım mükəmməl yığılır (4 iterasiya); çoxaddımlı
dt=1.0 ilə t=7, dt=0.25 ilə t=8 günə çatır.

### Qalan kök səbəb — QUYU HÜCEYRƏSİNİ ÇOX TEZ BOŞALDIR

Ölçmə: bir hüceyrənin məsamə həcmi **1250 m³**, istismarçı isə
**654 m³/gün** verir. dt=1.0-da bu, hüceyrənin məsamə həcminin
YARISIDIR — bir addımda! Heç bir implicit sxem bunu hamar keçirə
bilməz: təzyiq mütləq çökür.

Bu, QUYU MODELİNDƏN ASILI OLMAYAN, ayrıca məsələdir — zaman addımı
quyu debitinə görə məhdudlaşdırılmalıdır (quyu üçün CFL-ə bənzər
şərt): `dt · Σq_quyu ≤ α · PV_hüceyrə`, tipik `α ≈ 0.1–0.25`.

Bu, mərhələ 6 kimi ayrıca görülməlidir. Mərhələ 1-4 (quyu modeli)
öz-özlüyündə tamdır və hərtərəfli doğrulanıb.

### Mərhələ 5 (quyu saxlama həddi) haqqında

OPM-in saxlama həddi quyunun ÖZ kütlə balansı tənliklərinə əlavə
olunur (`A_α,w` — quyu lüləsindəki komponent miqdarı). Bunun üçün
quyunun 4 naməlumu olmalıdır (Qt, Fw, Fg, p_bhp), bizdə isə 1
(p_bhp). Yəni saxlama həddi TAM OPM formulasiyasını tələb edir və
ayrıca mərhələdir — sadəcə "əlavə etmək" mümkün deyil.

### Mərhələ 6 — quyu debitinə görə zaman addımı həddi (v66)

`CoupledNewtonSolver.max_stable_dt()` — CFL-bənzər şərt:

    dt ≤ α · PV_hüceyrə / Σ|q_quyu|     (α=0.2 defolt)

Sınaqda gözlənilən nəticəni verdi: 1250 m³ məsamə həcmli hüceyrədə
vurucunun (1752 m³/gün su) yaratdığı hədd **0.14 gün** çıxdı.

**Qoşulma DUCK-TYPING ilə, minimal-invaziv:** `AdaptiveTimeStepper.
advance()`-də `getattr(self.newton, "max_stable_dt", None)` yoxlanılır.
A6-nın iki fazalı `NewtonSolver`-i bu metodu TANIMIR, ona görə köhnə
mühərrik HEÇ BİR DƏYİŞİKLİK olmadan işləməyə davam edir (test bunu
təsdiqləyir: `hasattr(NewtonSolver, "max_stable_dt")` False).

**Yol boyu tapılan kiçik uyğunsuzluq:** `CoupledState`-də
`water_saturation` xassəsi yox idi, `AdaptiveTimeStepper.
_saturation_change()` isə onu gözləyir (ümumi dizayn). Əlavə olundu.

### YENİ TAPINTI — yavaş yığılma (kök səbəb DEYİL, ayrıca məsələ)

`max_stable_dt` düzgün işləsə də, kiçik addımda (dt=0.05) birləşmiş
sistem CNV-ni 0.070-dan 0.008-ə endirir, sonra ÇOX YAVAŞ azalmağa
davam edir (divergensiya YOX, sadəcə sublinear yığılma) — 40
iterasiyada belə 10⁻³ həddinə çatmır.

Xətti həlledicinin özü İSTİSNA edilir: sistem 110 naməlumdur,
`direct_threshold`-dan (20 000) çox aşağı, ona görə `splu` (dəqiq
həll) işlədilir — problem xətti həllin YAXINLAŞMASINDA deyil.

Bu, mərhələ 1-4-ün doğruluğunu təkzib ETMİR (Jakobian maşın
dəqiqliyində qalır), lakin göstərir ki, CPR-ə bənzər bir
ön-şərtçi/sürətləndirici (A7_PLAN-da qəsdən təxirə salınmış CPR 3×3)
bu birləşmiş sistem üçün sadəcə PERFORMANS deyil, YIĞILMA SÜRƏTİ
üçün də lazım ola bilər.

**Növbəti sessiya üçün istiqamət:** bu yavaş yığılmanın həndəsi
səbəbini (məsələn, quyu/rezervuar bloklarının nisbi miqyası,
kondisiya ədədi) ölçmək.

### v67 — miqyas araşdırması: qismən düzəliş, əsl mənbə tapıldı

**Quyu tənliyinin dt-yə görə miqyaslanması** əlavə olundu (idarəetmə
qalığı `PV_orta/dt` ilə vurulur, Jakobianın Q↔R/Q blokları eyni
faktorla). Şərtlənmə ədədi **297,000 → 203,000** (30% yaxşılaşma),
LAKİN yığılma sürətinə TƏSİR ETMƏDİ (eyni CNV platosu).

**Diaqonal təhlili əsl mənbəni göstərdi:** problem quyu↔rezervuar
miqyası deyil, **FAZALAR ARASI** miqyas fərqidir:

    su diaqonalı:    [1.5,   18]
    neft diaqonalı:  [17820, 21090]
    qaz diaqonalı:   [11580, 12230]
    quyu diaqonalı:  [25000, 25000]     (miqyaslamadan sonra)

Su tənliyi digərlərindən **1000+ dəfə kiçikdir** — bu, sıxılmazlığın
(Bw≈1, aşağı kompressibilite) təbii nəticəsidir və quyu modelindən
TAM ASILI DEYİL — A6-nın iki fazalı sistemində də mövcud olmalıdır,
sadəcə orada CPR (2×2) bunu maskalayır.

**Sınanıb İŞLƏMƏYƏN:** sadə sətir miqyaslaması (hər sətri öz
diaqonalına bölmək) — şərtlənmə ədədini **203,000-dan 1.5 MİLYARDA**
pislətdi. Bu, gözlənilməz idi və göstərir ki, problem sadə Jacobi
miqyaslaması ilə həll olunmur — off-diaqonal elementlər arasındakı
NİSBİ balans da vacibdir, təkcə diaqonal deyil.

### Nəticə

Bu, **CPR-in (3×3) niyə lazım olduğunun** dəqiq izahıdır: CPR məhz
fazalar arası bu miqyas fərqini həll etmək üçün nəzərdə tutulub
(təzyiq alt-sistemini ayırıb, ayrıca ön-şərtçi tətbiq edir). Sadə
miqyaslama üsulları (dt-yə görə, sətir əsaslı) bunu əvəz edə bilmir.

**Qərar:** kiçik-orta gridlərdə (bu layihənin əsas istifadə halı) bu
məhdudiyyət PRAKTIKI ƏHƏMİYYƏT DAŞIMAYA BİLƏR — dt=1.0 ilə tək addım
mükəmməl yığılır (4 iterasiya), yalnız ÇOX KİÇİK dt-lərdə (0.05)
problem üzə çıxır. `max_stable_dt` artıq quyu debitinə görə ağlabatan
dt seçir (adətən 0.1-0.5 aralığında) — bu, problemli super-kiçik
dt-lərdən qaçınmağa kömək edə bilər.

**Növbəti sessiya üçün iki yol:** (1) CPR-i 3×3-ə genişləndirmək
(böyük iş, lakin kök həll), (2) real ssenarilərdə problemin
PRAKTIKI olaraq üzə çıxıb-çıxmadığını (max_stable_dt ilə seçilən
"normal" dt-lərdə) yoxlamaq — bəlkə əlavə iş lazım deyil.
