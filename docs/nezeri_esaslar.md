# Nəzəri əsaslar — hansı model hansı nəzəriyyə ilə işləyir

**Yazılıb:** 25 sentyabr 2026 (Seans 47) · **Mənbə:** kodun özü (aşağıdakı
fayl istinadları) və `docs/teqdimat/slaydlar/slides/models.html` cədvəli.

Bu sənəd təqdimat slaydındakı «Hesablamanı hansı modellər aparır»
cədvəlinin GENİŞLƏNMİŞ formasıdır: hər sətir üçün (1) nəzəriyyə və ilkin
mənbə, (2) kodda FAKTİKİ yazılmış düstur, (3) fərziyyələr və hüdudlar,
(4) fayl. Heç bir düstur yaddaşdan yazılmayıb — hamısı göstərilən
fayldan oxunub.

---

## 0 · Bütün hesablamanın altındakı üç nəzəriyyə

Simulyatorun qalan hər şeyi bu üçünün üstündə qurulur:

| # | Nəzəriyyə | Nə deyir |
|---|---|---|
| 1 | **Darcy qanunu (1856)** | Məsaməli mühitdə faza axını potensial qradiyentinə düz, lözlüyə tərs mütənasibdir: `u_p = −(k·kr_p/μ_p)·∇Φ_p` |
| 2 | **Kütlənin saxlanması** | Hüceyrəyə girən − çıxan = yığılan. Hər faza üçün ayrı tənlik. |
| 3 | **Black-oil faza davranışı** | Yataqda üç komponent (su, neft, qaz) və üç faza var; qaz neftdə həll ola bilir (`Rs`), su isə heç nə həll etmir. |

Bu üçü birləşəndə hər hüceyrə üçün yazılan tənlik (kodda faktiki forma,
`simulation/implicit/residual.py`):

```
R_p,c = (PV_c / Δt) · [ (S_p/B_p)ⁿ⁺¹ − (S_p/B_p)ⁿ ]      ← akkumulyasiya
        − Σ_üzlər T · (λ_p/B_p)_upstream · ΔΦ_p            ← qonşu axını
        − q_p,c                                            ← quyu

Φ_o = p − ρ_o·g·D                 (neft potensialı)
Φ_w = p − Pc(Sw) − ρ_w·g·D        (su potensialı)
```

`R = 0` → tarazlıq. Tənliklər SƏTH həcmi vahidindədir (`B_p`-yə bölmə
kütləni məhz orada saxlayır).

**Qaz tənliyi fərqlidir** (`implicit/three_phase_residual.py`): qaz iki
yerdə saxlanılır — sərbəst faza kimi və neftdə həll olmuş kimi:

```
N_qaz = PV · ( Sg/Bg  +  So·Rs/Bo )
```

Bu, standart black-oil qaz balansıdır. Modeldə **Rv = 0** («quru qaz»
fərziyyəsi) — neft qaz fazasında buxarlanmır.

---

## 1 · Geostatistika — quyudan hüceyrəyə keçid

Nəzəri sual: 5 quyuda ölçülmüş keçiricilikdən 41×41×3 = 5043 hüceyrənin
dəyərini necə çıxarmaq olar?

### 1.1 Variogram nəzəriyyəsi (Matheron, «regionlaşmış dəyişənlər»)

Əsas ideya: iki nöqtə bir-birinə nə qədər yaxındırsa, dəyərləri bir o
qədər oxşardır. Bu oxşarlıq **yarım-dəyişkənlik** funksiyası ilə ölçülür:

```
γ(h) = ½ · E[ (Z(x) − Z(x+h))² ]
```

Kodda üç model qurulur (`geology/variogram.py`):

```
γ(h) = nugget + sill · g(h/a)        g(0) = 0,  g(∞) = 1

sferik:        g = 1.5t − 0.5t³  (t ≤ 1), sonra 1
eksponensial:  g = 1 − exp(−3t)
qauss:         g = 1 − exp(−3t²)
```

- `nugget` — sıfır məsafədə qalan sıçrayış (ölçmə xətası + çox kiçik
  miqyaslı dəyişkənlik).
- `sill` — QURULUŞLU hissə (nugget DAXİL DEYİL — layihənin konvensiyası).
- `range_` — **praktiki radius**: hər üç modeldə `γ(range_) ≈ 0.95·(nugget+sill)`.
  Ona görə modellər arasında birbaşa müqayisə oluna bilər.

Model quyu datasına **çəkili ən kiçik kvadratlarla** fit edilir
(`scipy.optimize.least_squares`), sonra doğrulanır — yararsız model
Kriging-ə ÇATMIR.

### 1.2 Anizotropluq (geometrik anizotropluq nəzəriyyəsi)

Çöküntü laylarında korrelyasiya istiqamətə görə fərqlidir (uzununa
100 m, eninə 30 m). Klassik həll: koordinatı çevirmək, sonra ADİ
Evklid məsafəsi işlətmək (`geology/anisotropy.py`):

```
u  = R₀ x              azimut dönməsi
u' = R_dip u           dip dönməsi
x''= S u'              S = diag(1, a_maj/a_min, a_maj/a_v)

d_ani(p,q) = ‖x''(p) − x''(q)‖₂
```

Layihədə bu, **yeganə** anizotrop-məsafə mənbəyidir — Kriging, qonşuluq
axtarışı, variogram, SGS və SIS hamısı eyni funksiyadan keçir.

### 1.3 Adi Kriging (BLUE — ən yaxşı xətti yansız qiymət)

Nəzəriyyə: qiymət qonşuların ÇƏKİLİ ortasıdır, çəkilər isə iki şərtlə
tapılır — (a) **yansızlıq**: `Σwᵢ = 1`, (b) **minimal qiymət variansı**.
Lagrange vuruğu ilə sistem:

```
[ Γ   1 ] [ w ]   [ γ₀ ]
[ 1ᵀ  0 ] [ μ ] = [ 1  ]

Γᵢⱼ = γ(d_ani(xᵢ, xⱼ)),   γ₀ᵢ = γ(d_ani(xᵢ, x₀))
```

Kriging-in İDW-dən üstünlüyü: variogramdan qonşuların bir-birini
«örtməsini» (clustering) nəzərə alır və hər hüceyrədə **kriging
variansı** verir — yəni qiymətin nə qədər etibarlı olduğunu.

Kod (`geology/interpolation.py`) hər hüceyrədə YERLİ (kiçik) sistem
qurur, ədədi problemi gizlətmir — həllin hansı yolla alındığı açıq
qeyd olunur: `direct` / `jitter` (diaqonal requlyarlaşdırma) / `lstsq`
(psevdo-tərs) / `renormalized` / `idw_fallback` / `exact_hard_data`.

### 1.4 IDW və ən yaxın qonşu

```
İDW:  Z(x₀) = Σ wᵢZᵢ / Σ wᵢ,   wᵢ = 1/dᵢᵖ
```

Nəzəri çatışmazlığı sənəddə açıq yazılıb: quyular arasında həmişə orta
dəyərə meyl edir («öküz gözü» effekti), məkan strukturunu tanımır.

### 1.5 Log-çevirmə

Keçiricilik **log-normal** paylanır (sənayedə qəbul olunmuş empirik
fakt), ona görə interpolyasiya `ln(k)` fəzasında aparılır
(`log_transform=True`), sonra `exp()` ilə geri qaytarılır. Kodda bunun
əsaslılığı ayrıca yoxlanılır (`geology/distribution_analysis.py`).

### 1.6 SGS — Sequential Gaussian Simulation (Deutsch & Journel, GSLIB)

Kriging **hamar** sahə verir — dəyişkənliyi süni azaldır. Şərti
simulyasiya bunu düzəldir:

1. **Normal-score çevirməsi** (`geology/gaussian_transform.py`) —
   sıra-əsaslı: Hazen mövqeyi `(rank − 0.5)/n` → `norm.ppf`. Nəticə
   cədvəldir, çevirmə onun xətti interpolyasiyasıdır, tam determinik.
2. Təsadüfi yol ilə hər hüceyrəyə gedilir, yerli kriging `mean` VƏ
   `variance` verir.
3. **Kriging qiyməti YAZILMIR** — `N(mean, variance)`-dən NÜMUNƏ çəkilir.
   (Bu, SGS-in əsas qaydasıdır; qiyməti yazmaq determinik hamar sahə
   verərdi.)
4. Nümunə kondisioner çoxluğa əlavə olunur — növbəti hüceyrə onu da görür.
5. Tərs normal-score → (varsa) `exp()` → fiziki hədd.

### 1.7 SIS — Sequential Indicator Simulation (fasiyalar üçün)

Kateqorik dəyəri (fasiya kodu) adi rəqəm kimi interpolyasiya etmək
mənasız «1.4» verir. Nəzəri həll — **indikator kriging**:

```
I_k(x) = 1  əgər fasiya(x) = k,  əks halda 0
```

Hər fasiya üçün ayrı indikator sahəsi Kriging edilir → nəticə həmin
fasiyanın EHTİMALIDIR → ehtimallar normallaşdırılır → təsadüfi nümunə
götürülür. `geology/facies.py`.

### 1.8 Elmi çəkincə (kodda yazılıb, burada təkrarlanır)

SGS və SIS **proqram cəhətdən** düzgündür (çevirmə doğru, varians
düzgün, sərt data hörmət olunur, təkrarlana bilir). Bu, nəticənin
YERİN gerçək paylanması olduğu demək DEYİL — o, yalnız kifayət qədər
təmsiledici quyu sıxlığı və DOĞRU variogram modeli ilə mümkündür.

---

## 2 · PVT — flüid xassələri

Bütün korrelyasiyalar `simulation/pvt/correlations.py`-dədir. Mühərrik
korrelyasiyaları GÖRMÜR: onlar yalnız `PVTTable` istehsal edir, mühərrik
isə həmişə cədvəllə (xətti interpolyasiya, `np.interp`) işləyir.

| Kəmiyyət | Korrelyasiya | Kodda faktiki forma |
|---|---|---|
| `Rs(p)`, `Pb` | **Standing (1947)** | `Rs = γg·[(p/18.2 + 1.4)·10^x]^1.2048`, `x = 0.0125·API − 0.00091·T` |
| `Bo(Rs)` | **Vazquez-Beggs (1980)** | API ≤ 30 və > 30 üçün AYRI əmsal dəsti |
| `co` (doymamış sıxılma) | **Vazquez-Beggs (1980)** | `co = (−1433 + 5Rs + 17.2T − 1180γg + 12.61API)/(10⁵·p)` |
| `μod` (ölü neft) | **Beggs-Robinson (1975)** | `μod = 10^x − 1`, `x = 10^(3.0324 − 0.02023·API)·T^−1.163` |
| `μo(Rs)` (doymuş) | **Beggs-Robinson (1975)** | `μo = A·μod^B`, `A = 10.715(Rs+100)^−0.515`, `B = 5.44(Rs+150)^−0.338` |
| `Bw`, `μw` | **McCain tipli / Meehan** | temperatur və duzluluqdan |
| `Tpc`, `Ppc` | **Sutton (1985)** | kompozisiya analizi olmadan, yalnız `γg`-dən psevdo-kritik nöqtə |
| `Z` | **Beggs-Brill (1973)** | Standing-Katz əyrisinin AÇIQ (iterasiyasız) approksimasiyası |
| `Bg` | **həqiqi qaz qanunu** | `Bg = (Z·T/p) / (Zsc·Tsc/psc)` |
| `μg` | **Lee-Gonzalez-Eakin (1966)** | `μg = K·10⁻⁴·exp(X·ρ^Y)` — qaz sıxlığından |

### 2.1 Niyə Beggs-Brill, Dranchuk-Abou-Kassem yox

DAK tənliyi hər çağırışda **iterasiya** tələb edir. Z hər Nyuton
iterasiyasında min hüceyrə üçün hesablanır — qapalı formalı Beggs-Brill
vektor əməliyyatına oturur. Bu, dəqiqlik–sürət mübadiləsinin AÇIQ
seçimidir.

### 2.2 Doyma nöqtəsində şaxələnmə (ən vacib PVT nəzəriyyəsi)

```
p < Pb  (DOYMUŞ — qaz ayrılır):
    Rs = Standing(p)                 p azaldıqca AZALIR
    Bo = Vazquez-Beggs(Rs(p))        p azaldıqca AZALIR
    μo = Beggs-Robinson(μod, Rs(p))  p azaldıqca ARTIR

p ≥ Pb  (DOYMAMIŞ — tək fazalı maye):
    Rs = Rsb = sabit
    Bo = Bob·exp(co·(Pb − p))        sadəcə mayenin sıxılması
    μo = μob·(p/Pb)^0.278
```

Hər iki qol EYNİ anchor-dan (Pb-də `Rs = Rsb`) qurulur, ona görə Pb-də
kəsilməzlik maşın dəqiqliyi ilə təminatlıdır. Kodda bu bölmənin
regressiya tarixçəsi də saxlanılıb: bir dəfə Bo/μo bütün diapazonda
doymamış düsturla hesablanmışdı (Nyuton ossilyasiyası üçün) — bu,
fizikanı dağıtdığı üçün geri qaytarıldı (μo(p→0) → 0 qeyri-fiziki idi),
ossilyasiya isə xətti axtarışla (line search) həll olundu.

### 2.3 Deck yolu (SPE1 kimi real modellərdə)

Laboratoriya `PVTO` cədvəli veriləndə korrelyasiya İŞLƏMİR: doymamış
sıxılma və özlülük üstəli hər qolun ÖZ sətirlərindən alınır və Rs üzrə
**xətti** interpolyasiya olunur (Q-34), doymuş qol deck düyünlərindən
yuxarı XƏTTİ uzadılır (Q-32) — OPM Flow ilə eyni qayda.

---

## 3 · Nisbi keçiricilik və kapilyar təzyiq (SCAL)

### 3.1 Corey modeli (Brooks-Corey ailəsi)

```
Sw_norm = (Sw − Swc) / (1 − Swc − Sor)

krw = krw_end · Sw_norm^nw
kro = kro_end · (1 − Sw_norm)^no
```

Nəzəri əsas: məsamə ölçüsü paylanmasının üstəl qanunu ilə təsviri
(Brooks & Corey, 1964). `nw`, `no` — əyrilik göstəriciləri; `Swc`, `Sor`
— hərəkətsiz doyumluluqlar.

### 3.2 Cədvəl yolu (SWOF / SGOF)

Real kern datası Corey-yə oturmur (asimmetrik, son nöqtələr ayrıca
ölçülür). Ona görə Eclipse `SWOF`/`SGOF` formatı birbaşa oxunur
(`domain/scal_tables.py`), `SATNUM` ilə region-region tətbiq olunur.
**Monotonluq sərt tələbdir** (`krw` artan, `kro` azalan) — pozulsa
diskretizasiya qeyri-stabil olur.

### 3.3 Stone II (1973) — üç fazada neftin keçiriciliyi

İki fazada `kro` tək dəyişəndən asılıdır. Üç fazada neft həm suyun, həm
qazın sıxışdırmasına məruz qalır. Stone-un ideyası: iki İKİFAZALI
ölçülmüş əyrini birləşdirmək:

```
kro = kro_end · [ (krow/kro_end + krw)·(krog/kro_end + krg) − (krw + krg) ]
```

Stone I (1970) mənfi `kro` verə bilirdi; Stone II bunu riyazi olaraq
aradan qaldırır və sənaye standartıdır (Eclipse defoltu).

**Ən güclü düzgünlük sınağı** — iki fazalı hala DƏQİQ reduksiya:
`Sg = 0` → `kro = krow(Sw)`; `Sw = Swc` → `kro = krog(Sg)`. Testlə
yoxlanılıb.

### 3.4 Brooks-Corey kapilyar təzyiqi

```
Pc(Sw) = Pe · Sw_norm^(−1/λ)
```

`Pe` — giriş (entry) təzyiqi, `λ` — məsamə ölçüsü paylanması indeksi.
`Sw → Swc` olduqda `Pc → ∞`, ona görə yuxarı kəsmə həddi qoyulur.
Törəmə analitikdir (Jakobian üçün). İki provider var, eyni ailə:
`Pcow = Po − Pw` (Sw-dən) və `Pcog = Pg − Po` (Sg-dən).

---

## 4 · Diskretizasiya — tənliyi hüceyrələrə necə yazırıq

### 4.1 Sonlu həcm metodu (finite volume)

Nəzəriyyə: tənlik hər hüceyrənin həcmi üzrə inteqrallanır, divergensiya
həddi Qauss teoremi ilə ÜZ axınlarının cəminə çevrilir. Nəticə: kütlə
**lokal olaraq** saxlanılır — bu, sonlu fərq və sonlu elementdən üstün
cəhətdir və layihədə açıq müqavilə kimi yazılıb:

```
Σ(üz axınları) + quyu/mənbə = akkumulyasiya
```

Şərt: `compute_flux()` üzün hər iki tərəfinə EYNİ ədədi qaytarmalıdır.

### 4.2 TPFA (two-point flux approximation) — defolt

Üz axını yalnız İKİ qonşu hüceyrənin təzyiqindən asılıdır:

```
T = C_darcy · A / (d_a/k_a + d_b/k_b)        ← harmonik orta
q_üz = T · (λ/B)_upstream · ΔΦ
```

Harmonik orta ardıcıl müqavimətlərin toplanmasından gəlir (az keçirici
hissə axını idarə edir). **Nəzəri hüdudu:** TPFA yalnız K-ortoqonal
şəbəkədə ardıcıldır (consistent) — tam tenzor `K` və ya əyri şəbəkədə
sistematik xəta verir. Kod bunu gizlətmir: tam tenzor aşkarlananda
xəbərdarlıq yazılır.

### 4.3 Upstream çəkilənmə

Mobillik `λ = kr/μ` **axının gəldiyi** hüceyrədən götürülür. Nəzəri
səbəb: doyumluluq tənliyi hiperbolikdir, mərkəzi çəkilənmə qeyri-fiziki
rəqs verir. Bu, birinci tərtib «upwind» sxemidir — dayanıqlıdır, amma
ədədi dispersiya (cəbhənin süni yayılması) gətirir.

### 4.4 MPFA-O — çoxnöqtəli axın approksimasiyası

Tam keçiricilik tenzoru və qeyri-ortoqonal şəbəkə üçün. Nəzəriyyə:
hüceyrə küncləri ətrafında **qarşılıqlı təsir zonası** (interaction
region) qurulur, orada təzyiqin kəsilməzliyi və axının kəsilməzliyi
eyni vaxtda tələb olunur, nəticədə üz axını N hüceyrənin təzyiqindən
asılı olur:

```
q_F = Σ_c T_cell[F,c] · p_c  +  Σ_β T_bnd[F,β] · π_β
```

`T_cell` yalnız HƏNDƏSƏ, TOPOLOGİYA və `K`-dan asılıdır — təzyiq,
doyumluluq, mobillik və PVT-dən ASILI DEYİL, ona görə bir dəfə qurulur
və Jakobianda törəməsi analitik alınır. Sənədlər:
`docs/mpfa_o_phase5a.md`, `…_validation.md`, `…_phase5b1.md`.

**Məhdudiyyət:** IMPES MPFA-O ilə İŞLƏMİR (açıq `NotImplementedError`)
— IMPES tək-üz skalyar transmissivlik tələb edir, MPFA-O-da belə
kəmiyyət yoxdur.

---

## 5 · Hesablama mühərrikləri

### 5.1 IMPES (Implicit Pressure, Explicit Saturation)

Nəzəriyyə: tənliklər elə birləşdirilir ki, doyumluluq həddi yox olsun və
yalnız təzyiq üçün xətti sistem qalsın; doyumluluq sonra açıq (explicit)
yenilənir.

- Üstünlüyü: hər addım ucuzdur (bir xətti sistem, Nyuton yoxdur).
- Nəzəri qiyməti: **şərti dayanıqlıdır** — CFL şərti:

```
Δt ≤ min_c  PV_c / (throughput_c · max|dfw/dSw|)
```

`fw = λw/(λw+λo)` — fraksional axın. PVT halında lözlüklər təzyiqdən
asılı olduğu üçün `max|dfw/dSw|` bütün təzyiq aralığında taranır.

- IMPES AÇIQ şəkildə imtina edir: MPFA-O, THP idarəsi, RATE quyusunun
  BHP limiti, səth bazalı debit. Səbəb hər yerdə eynidir — səssizcə
  yanlış fizika verməkdənsə, imtina etmək.

### 5.2 Tam implicit (fully implicit) + Nyuton-Rafson

Zamanda **geriyə Eyler** (backward Euler): bütün hədlər yeni qatda
qiymətləndirilir → şərtsiz dayanıqlı, amma qeyri-xətti sistem yaranır.
Həlli Nyuton-Rafson ilə:

```
J(x^k) · δ = −R(x^k)
x^{k+1} = x^k + chop(δ)
```

`chop` — dəyişmə həddi (saturation chopping): bir iterasiyada ΔSw-nin
böyük sıçrayışı kəsilir, əks halda Nyuton fiziki olmayan bölgəyə düşür.
Çətin addımlarda geri-izləmə (line search) işə düşür.

**İki konvergensiya meyarı (Eclipse/CMG standartı, hər ikisi ölçüsüz):**

```
CNV_p = max_c |R_p,c|·Δt / PV_c                    ← lokal: ən pis hüceyrə
MB_p  = |Σ_c R_p,c|·Δt / Σ_c (PV_c·S_p/B_p)        ← qlobal kütlə balansı
```

MB CNV-dən qat-qat sıx olmalıdır: lokal səhvlər bir-birini yeyə bilər,
ümumi kütlə itkisi isə yolverilməzdir. Ölçüsüz olduqları üçün eyni
tolerans bütün grid ölçülərində işləyir.

### 5.3 Üç fazalı black-oil və DƏYİŞƏN KEÇİD (variable switching)

Üçüncü primary dəyişənin MƏNASI hüceyrənin vəziyyətindən asılıdır —
Eclipse/IMEX-in standart üsulu:

```
doymuş hüceyrə    (Sg > 0)  →  3-cü dəyişən = Sg
doymamış hüceyrə  (Sg = 0)  →  3-cü dəyişən = Rs
```

Səbəb: doymamış hüceyrədə `Sg ≡ 0` (Nyuton sabiti «həll edərdi»),
doymuş hüceyrədə isə `Rs ≡ Rs_sat(p)` (təzyiqdən çıxarıla bilir).
Keçid şərtləri:

```
doymamış → doymuş   Rs > Rs_sat(p)   (neft doydu, qaz ayrılır)
doymuş → doymamış   Sg < 0           (bütün sərbəst qaz həll oldu)
```

Sərhəddə kəsilməzlik tələb olunur: `Sg = 0 ⟺ Rs = Rs_sat(p)`.

### 5.4 Adaptiv zaman addımı

```
Nyuton yığılmadı  → Δt kəs, addımı TƏKRARLA   (×0.5)
az iterasiya      → Δt böyüt                   (×1.5)
çox iterasiya     → Δt kiçilt (təkrar yox)
```

Əlavə məhdudiyyət: bir addımda doyumluluq dəyişməsi (ΔSw) həddi —
cəbhənin çox irəliləməsi həm Nyutonu çətinləşdirir, həm dəqiqliyi
azaldır. Ölçülüb: yığılma `Δt·q/PV` nisbətindən asılıdır, hədd ≈ 10–15.

### 5.5 Xətti həll: CPR ön-şərtçisi

Jakobian iki FƏRQLİ xarakterli tənliyi birləşdirir:

- **təzyiq** — elliptik, uzaqmənzilli, qlobal həll tələb edir;
- **doyumluluq** — hiperbolik, lokal, upstream istiqamətdə yayılır.

Tək ILU hər ikisini eyni cür emal edir və elliptik hissədə zəif qalır.
**CPR (Constrained Pressure Residual)** iki mərhələlidir: təzyiq
alt-sistemi ayrılıb dəqiq həll olunur, qalan qalıq tam sistemə ILU ilə
hamarlanır.

Dekuplinq (quasi-IMPES): hər hüceyrənin 2×2 blokundan çəkilər elə
seçilir ki, doyumluluq törəməsi yox olsun:

```
w = [ ∂R_o/∂Sw ,  −∂R_w/∂Sw ]
M⁻¹r = P·A_p⁻¹·(W·r) + M_ILU⁻¹·( r − A·P·A_p⁻¹·(W·r) )
```

Sabit `B` halında bu, `w = [Bw, Bo]` — yəni klassik həcm balansına
çevrilir.

Xətti həlledici: SciPy CG / BiCGStab + ILU ön-şərtçi, warm start,
uğursuzluqda birbaşa həllə qayıtma (`simulation/linear_solver.py`).

---

## 6 · Quyu modelləri

### 6.1 Peaceman quyu indeksi (1978, 1983)

Nəzəri problem: hüceyrə 20×20 m-dir, quyu radiusu 0.1 m. Hüceyrənin
ORTA təzyiqi quyu dibi təzyiqi DEYİL. Peaceman göstərdi ki, radial axın
həlli ilə uzlaşma üçün **ekvivalent radius** var:

```
WI = 2π·C·√(K1·K2)·h / ( ln(re/rw) + S )

re = 0.28·√( √(K2/K1)·d1² + √(K1/K2)·d2² ) / ( (K2/K1)^0.25 + (K1/K2)^0.25 )
```

Bu, anizotrop blokun (Peaceman 1983) formasıdır; izotrop kvadrat blokda
klassik `re ≈ 0.2·dx`-ə çevrilir. `S` — skin (lay zədəsi/stimullaşdırma).
Həndəsə həqiqi wellblock həndəsəsindən gəlir, ox-boyu sərhəd qutusundan
YOX (`domain/geometry.py::WellblockGeometry`).

Debit:
```
q_p = WI · (λ_p/B_p) · (p_hüceyrə − BHP)
```

### 6.2 Quyu qaydaları (Nyuton həlləri ARASINDA tətbiq olunur)

Bunlar qalıq/Jakobianı DƏYİŞMİR — addımlar arasında hədəfi dəyişirlər:

| Qayda | Nəzəri məzmun | Qərar |
|---|---|---|
| `assign_rate_shares` | RATE hədəfi perforasiyalara **mobilliyə görə** paylanır | Q-19 |
| `BhpLimitController` | RATE quyusu BHP həddini keçəndə rejim RATE → BHP olur | Q-21 |
| `SurfaceRateController` | SƏTH debiti `(1−f)/Bo`, `1/B` əmsalları ilə lay həcminə çevrilir | Q-24 |
| vurucu bağlantısı | vurucuda hüceyrənin **TAM mobilliyi** işlədilir (OPM qaydası) | Q-33 |

### 6.3 Quyu lüləsi hidravlikası — Darcy-Weisbach + Chen (1979)

Lay tənlikləri ilə lülənin içi QƏSDƏN ayrıdır; mühərrik yalnız BHP
tanıyır. Çoxseqmentli təzyiq traversi:

```
dp/dz = [ ρ_qrav·g + f·ρ_ns·v_m²/(2d) ] / 1e5      [bar/m]
THP = BHP − Σ( Δp_qravitasiya + Δp_sürtünmə )
```

Sürtünmə əmsalı — **Chen (1979)**, Colebrook-White-ın AÇIQ (iterasiyasız)
yaxınlaşması, turbulent aralıqda fərq < 0.5 %. Laminar aralıqda analitik
`f = 64/Re`. (Tam istinad koddadır: *Ind. Eng. Chem. Fundam.* 18(3),
296-297.)

**Nəzəri sadələşdirmələr (açıq yazılıb, uydurulmayıb):**

- **Sürətlənmə həddi** (`ρ·v·dv/dz`) DAXİL DEYİL — ümumi itkinin < 1 %-i, ⏳ V2.
- **Sürüşmə (slip) YOXDUR**: `H_L = λ_L` (no-slip holdup). Nəticə:
  qazlı quyuda hidrostatik sütun OLDUĞUNDAN YÜNGÜL çıxır, hesablanan
  THP həqiqi dəyərdən YÜKSƏK olur; fərq GOR artdıqca böyüyür.
  Hagedorn-Brown `CNL` və `ψ` QRAFİKLƏRİNİN rəqəmsallaşdırılmasını
  tələb edir — o cədvəl datası bizdə yoxdur və uydurulmayacaq.
  Beggs-Brill holdup-u qapalı düsturludur və sonradan eyni interfeysə
  oturur.
- Çox seqment lazımdır, çünki qaz yuxarı qalxdıqca `Bg` böyüyür, `ρ_m`
  azalır, `v_m` artır — tək seqment yalnız qazsız quyuda doğrudur.

### 6.4 Nodal analiz — IPR ∩ VLP

Klassik istehsal mühəndisliyi nəzəriyyəsi: quyu iki əyrinin kəsişməsində
işləyir.

```
IPR (lay):     q(BHP) = J · (p_lay − BHP)      düz xətt, J son addımdan
VLP (lülə):    THP(BHP, q(BHP))                çoxseqmentli traverse
```

`ThpController` iş nöqtəsini ikiqat bölmə (bisection) ilə tapır, addım
daxilində yarı-implicit təkrar aparır və axa bilməyən quyunu bağlayır
(`well_index = 0`).

**Birləşmə AÇIQDIR (explicit):** BHP bir addım əvvəlki debitlərə görə
hesablanır. Bilinən riski — gecikmə/rəqs; qarşısı relaksasiya (ω = 0.5)
və addım başına maksimal dəyişmə (25 bar) ilə alınır. Tam implicit THP
birləşməsi ⏳ sonraya.

---

## 7 · İlkin şərtlər — hidrostatik tarazlıq

```
p(z) = p_datum + ρ·g·(z − z_datum) / 1e5        [bar]
```

Sıxlıqlar LAY şəraitinə çevrilir: `ρ_lay = ρ_səth / B`.

**Kapilyar-cazibə tarazlığı** (keçid zonası nəzəriyyəsi): OWC-dən
yuxarıda kapilyar təzyiq sütunun çəkisi ilə tarazlaşır —

```
Pc(D) = (ρw − ρo)·g·(D_owc − D) / 1e5,     Sw = Pc⁻¹( Pc(D) )
```

Kapilyar provider verilmədikdə kontakt kəskin qalır (köhnə davranış).

**Qaz papağı:** GOC verilibsə üç zona qurulur (qaz / neft / su). Qaz-neft
kontaktında kapilyar keçid zonası ⏳ hələ qurulmur (kəskin sərhəd).

**Canlı neft (G9c):** PVT-də qaz varsa neftin lay sıxlığına həll olmuş
qazın kütləsi də daxildir:

```
ρo = (ρo_səth + Rs·ρg_səth) / Bo(p, Rs)
```

Bu vacibdir: əks halda ilkin təzyiq profili mühərrikin öz cazibə həddi
ilə TARAZ olmaz və ilk addımda süni axın yaranardı.

---

## 8 · Doğrulama nəzəriyyəsi

### 8.1 Buckley-Leverett (1942) — analitik etalon

Frontal advance nəzəriyyəsi: sıxılmayan, kapilyarsız, cazibəsiz 1D
sıxışdırmada doyumluluq **xarakteristikalar** boyunca yayılır:

```
x(Sw, t) = (q·t)/(φ·A) · dfw/dSw
```

Bu, çoxqiymətli profil verir — fiziki olmayan hissə **Welge (1952)
toxunan qaydası** ilə sıçrayışla (shock) əvəz olunur: `Swi`-dən
fraksional axın əyrisinə çəkilən toxunanın toxunma nöqtəsi `Sw_shock`-dur
(kodda: `chord` maksimumu).

Profil dörd hissədən qurulur: süpürülmüş zona → rarefaksiya yelpiyi →
şaquli sıçrayış → toxunulmamış zona.

**Doğruluq meyarı** — kütlə eyniliyi DƏQİQ ödənməlidir:

```
∫₀^∞ (Sw(x) − Swi) dx = q·t/(φ·A)
```

Ölçülmüş nəticə (Seans 44, 1D 120 hüceyrə, 250 gün, IMPES): RMS xəta
0.0231 Sw, cəbhə 190.9 m (analitik 182.6 m → 4.53 %), həcm balansı 2.71 %.
Fərqin nəzəri səbəbi bilinir: upstream sxemin **ədədi dispersiyası**.

### 8.2 SPE1 — sənaye etalonu

Odeh (1981) müqayisə məsələsi: qaz vurma ilə üç fazalı black-oil.
Müqayisə OPM Flow ilə aparılır (`SPE1.md`). Bu, korrelyasiya yolunu YOX,
**deck yolunu** (PVTO/PVDG/SWOF/SGOF cədvəlləri) yoxlayır.

---

## 9 · Uyğunlaşdırma və həssaslıq

### 9.1 Uyğunsuzluq funksiyası (ən kiçik kvadratlar nəzəriyyəsi)

```
SSE   = Σ w·(hesablanmış − ölçülmüş)²
RMSE  = √(SSE/n)
NRMSE = RMSE / (müşahidənin diapazonu)      ← ÖLÇÜSÜZ
```

Yekun — NRMSE-lərin çəkili ortası. Ölçüsüzlük vacibdir: təzyiq barla,
debit m³/günlə ölçülür; xam SSE-ləri toplamaq böyük ədədli kəmiyyəti
süni üstün edərdi. Model nəticəsi müşahidə vaxtlarına interpolyasiya
olunur — əksinə YOX, çünki müşahidə həqiqətdir.

### 9.2 Optimallaşdırma üsulları (hamısı **törəməsiz**)

| Üsul | Nəzəri xarakter |
|---|---|
| **Nelder-Mead** | simpleks; lokal, sürətli, başlanğıcdan asılı |
| **Powell** | ardıcıl istiqamətlər üzrə birölçülü axtarış; lokal |
| **Differential Evolution** | populyasiya əsaslı qlobal axtarış; bahalı, lokal minimumdan çıxa bilir |

Törəməsiz seçim təsadüfi deyil: hər qiymətləndirmə bir SİMULYASİYADIR,
analitik qradiyent mövcud deyil (adjoint metodu ⏳ yoxdur).

Axtarış `[0,1]` vahid fəzasında aparılır ki, müxtəlif vahidli
parametrlər bir-birinə mane olmasın. Yığılmayan model İSTİSNA ATMIR —
böyük, lakin SONLU cərimə (`1e6`) qaytarır; sonsuzluq simpleksi
sındırardı və DE seçim apara bilməzdi.

### 9.3 Həssaslıq — iki fərqli sual

```
Tornado (diapazon):  hər parametr öz TAM hədləri arasında (one-at-a-time)
Elastiklik (lokal):  (ΔÇıxış/Çıxış)/(ΔParametr/Parametr), mərkəzi fərq
```

Geniş hədli, lakin lokal olaraq zəif parametr Tornado-da yuxarıda,
elastiklikdə aşağıda görünə bilər — ikisi eyni sual DEYİL.

---

## 10 · Modelin nəzəri hüdudları (nə DAXİL DEYİL)

Bunlar qüsur deyil, **modelin sərhədidir** — təqdimatda soruşula bilər:

| Hüdud | İzah |
|---|---|
| **İzotermik** | temperatur girişdə sabitdir; istilik tənliyi yoxdur (buxar/SAGD modelləşdirilmir) |
| **Rv = 0** | quru qaz: neft qaz fazasında buxarlanmır (yaş qaz/kondensat deyil) |
| **Kompozisiya yoxdur** | black-oil üç komponentlə kifayətlənir; EOS (Peng-Robinson və s.) yoxdur |
| **Diffuziya/dispersiya yoxdur** | yalnız konveksiya (Darcy axını) |
| **Tək məsaməlilik** | dual-porosity / dual-permeability (təbii çatlı lay) yoxdur |
| **Süxur sadə sıxılır** | `c_r` sabit; geomexanika, çökmə, gərginlikdən asılı keçiricilik yoxdur |
| **Kimyəvi/EOR yoxdur** | polimer, səthi-aktiv maddə, qarışan (miscible) vurma modelləşdirilmir |
| **Sürüşmə yoxdur (lülə)** | no-slip holdup — bax §6.3 |
| **TPFA ardıcıllığı** | tam tenzor `K` və ya əyri şəbəkədə TPFA sistematik xəta verir → MPFA-O bunun üçündür |
| **Ədədi dispersiya** | birinci tərtib upstream sxem cəbhəni süni yayır (BL müqayisəsində ölçülüb) |

---

## 11 · Cədvəl: slayd sətri → nəzəriyyə → fayl

| Slayd sətri | Nəzəriyyə | Fayl |
|---|---|---|
| Geostatistika | Matheron variogram nəzəriyyəsi, BLUE-Kriging, Deutsch & Journel SGS/SIS | `geology/variogram.py`, `interpolation.py`, `sgs.py`, `facies.py`, `anisotropy.py`, `gaussian_transform.py` |
| PVT | Standing, Vazquez-Beggs, Beggs-Robinson, Sutton, Beggs-Brill, Lee-Gonzalez-Eakin | `simulation/pvt/correlations.py`, `black_oil.py` |
| Nisbi keçiricilik | Brooks-Corey üstəl qanunu, Stone II (1973), laboratoriya cədvəli | `domain/scal.py`, `domain/scal_tables.py`, `simulation/stone_relperm.py` |
| Kapilyar təzyiq | Brooks-Corey (1964) | `simulation/capillary.py` |
| Quyu | Peaceman (1978/1983), Darcy-Weisbach + Chen (1979), nodal analiz | `simulation/well_model.py`, `simulation/wellbore/*.py` |
| Diskretizasiya | Sonlu həcm, TPFA (harmonik orta), MPFA-O | `simulation/discretization.py`, `discretization/mpfa_o*.py` |
| Mühərrik | IMPES + CFL, backward Euler + Nyuton-Rafson, CNV/MB, CPR, dəyişən keçid | `simulation/impes_engine.py`, `simulation/implicit/*.py` |
| Uyğunlaşdırma | Ən kiçik kvadratlar (NRMSE), törəməsiz optimallaşdırma | `history/mismatch.py`, `optimizer.py`, `sensitivity.py` |
| Doğrulama | Buckley-Leverett (1942) + Welge (1952), SPE1 (Odeh 1981) | `simulation/analytical.py`, `benchmarks/spe1.py` |

---

## 12 · Açıq qalanlar (⏳)

- ⏳ Orijinal məqalələrin **səhifə/tənlik nömrələri** bu sənədə
  yazılmayıb — kodda əksər hallarda yalnız müəllif və il göstərilib
  (istisna: Chen 1979, tam istinadı `wellbore/friction.py`-dədir).
  Diplom/məqalə üçün biblioqrafiya ayrıca hazırlanmalıdır.
- ⏳ Üç fazalı Jakobianın analitik törəmələrinin tam çıxarılışı (dəyişən
  keçid halında) bu sənəddə yoxdur — `implicit/derivatives.py` və
  `three_phase_newton.py`-dədir; `docs/README.md`-də gözlənilən
  `riyaziyyat.md` bunun üçündür.
- ⏳ MPFA-O-nun lokal sisteminin (interaction region) tam riyazi
  çıxarılışı burada təkrarlanmır — `docs/mpfa_o_phase5a.md`-dədir.
