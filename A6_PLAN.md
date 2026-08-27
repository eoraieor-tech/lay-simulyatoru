# A6 — Fully implicit sxem

Məqsəd: CFL məhdudiyyətini aradan qaldırmaq. `PERFORMANCE.md`-də
ölçüldü ki, IMPES-in orta zaman addımı **0.34 gündə** ilişib qalır və
addım sayı müddətlə xətti artır. Fully implicit sxemdə Δt 20–30 gün
ola bilər.

## Mərhələlər

| # | İş | Vəziyyət |
|---|---|---|
| **1** | Qalıq (residual) vektoru + kütlə balansı yoxlaması | **HAZIR** |
| **2** | Analitik Jakobian | **HAZIR** |
| **3** | Nyuton döngəsi + konvergensiya meyarları | **HAZIR** |
| **4** | Adaptiv zaman addımı | **HAZIR** |
| **5** | CPR ön-şərtçisi | **HAZIR** |

---

## Mərhələ 1 — qalıq vektoru

### Primary dəyişənlər

    x = [p₀, Sw₀, p₁, Sw₁, …, p_{N−1}, Sw_{N−1}]

Hüceyrə üzrə növbələşdirmə (interleaved) qəsdən seçilib: Jakobian 2×2
blok-diaqonal struktura düşür və mərhələ 5-dəki CPR ön-şərtçisi məhz
bunu tələb edir.

### Tənliklər

Hər hüceyrə üçün iki kütlə balansı (səth həcmi vahidlərində):

    R_p,c = (PV(p)/Δt)·[(S_p/B_p)ⁿ⁺¹ − (S_p/B_p)ⁿ]
            − Σ_üzlər T·(λ_p/B_p)_upstream·ΔΦ_p
            − q_p,c

    Φ_o = p          − ρ_o·g·D
    Φ_w = p − Pc(Sw) − ρ_w·g·D

Bütün axın və quyu formulları IMPES mühərriki ilə **eynidir** — fərq
yalnız ondadır ki, burada onlar yeni zaman qatında qiymətləndirilir.

### Yoxlanılan xassələr

| Xassə | Nəticə |
|---|---|
| Quyusuz, bərabər vəziyyətdə R = 0 | maşın dəqiqliyi (< 1e-12) |
| Σ R_p = −Σ q_p (daxili axınlar yeyilir) | nisbi xəta < 1e-6 |
| İki hüceyrə arasında axın antisimmetrikdir | Σ = 0 |
| Akkumulyasiya = PV·S/B | dəqiq |
| IMPES addımının qalığı Δt ilə azalır | 1.16e-1 → 9.15e-4 |

Sonuncu ən vacibidir: iki sxemin **ardıcıllığını (consistency)**
təsdiqləyir.

### Konvergensiya ölçüsü

Xam qalıq m³/gün vahidindədir və hüceyrə həcmindən asılıdır. Ona görə
məsamə həcminə normallaşdırılır (CNV tipli ölçü):

    ‖R‖ = max_c |R_p,c| / (PV_c/Δt)

Test ilə yoxlanılıb: eyni fiziki vəziyyətdə 7×7 və 15×15 gridlərdə
norma eynidir.

---

## Mərhələ 2 — Jakobian

    J[i, j] = ∂R_i / ∂x_j,   ölçü 2N × 2N
    sətir  = 2c + faza       (0 = su, 1 = neft)
    sütun  = 2c + dəyişən    (0 = təzyiq, 1 = Sw)

Üç töhfə: akkumulyasiya (diaqonal), üzlər üzrə axın (diaqonal +
qonşu), quyular (diaqonal). Seyrəklik strukturu bir dəfə qurulur.

### Törəmələrin mənbəyi

| Törəmə | Üsul |
|---|---|
| `dkrw/dSw`, `dkro/dSw` | **analitik** — Corey düsturundan |
| `dPc/dSw` | **analitik** — A4-də yazılıb |
| `dB/dp`, `dμ/dp` (PVT) | **analitik** — parçalı xətti cədvəlin interval meyli |
| `dB/dp` (PVT-siz) | **analitik** — sıxılma düsturundan |
| Quyu debitləri | analitik, Peaceman düsturundan |

`DerivativeProvider` provider-də analitik metod varsa onu işlədir,
yoxsa sonlu fərqə keçir — interfeyslər dəyişmir.

### Yoxlama — sonlu fərqlə element-element müqayisə

| Konfiqurasiya | Maks. nisbi fərq |
|---|---|
| Baza (BHP quyular) | ~1e-10 |
| + kapilyar təzyiq | ~1e-10 |
| + PVT | ~1e-10 |
| + PVT + kapilyar | ~1e-10 |
| RATE idarəetmə | ~1e-10 |
| Cazibə olan modellər | ~1e-5 |

### Qəbul edilən sadələşdirmələr

**1. Cazibə üzvündə `∂ρ/∂p` buraxılıb.** `ρ = ρ_səth/B(p)` olduğuna
görə sıxlıq təzyiqdən asılıdır. Xəta yalnız dərinlik fərqi olan
modellərdə görünür və 1e-3-dən kiçikdir. Nyuton buna baxmayaraq
sürətlə yığılır, çünki konvergensiya **qalığa** görə yoxlanılır.

**2. Upstream seçiminin özü diferensiallaşdırılmır.** ΔΦ işarəsi
dəyişəndə funksiya kəsilməzdir, törəməsi isə sıçrayır. Sənaye
simulyatorlarında da belədir.

---

## Mərhələ 3 — Nyuton döngəsi

    x⁰ = xⁿ
    təkrarla:
        J(x^k) · δ = −R(x^k)
        x^{k+1} = x^k + chop(δ)
        konvergensiya yoxlanılır

### Xətti həlledici

Jakobian **simmetrik deyil** (upstream çəkilənmə və quyu üzvləri
simmetriyanı pozur), ona görə IMPES-in KQ həlledicisi yaramır.
`NewtonLinearSolver`: kiçik sistemlərdə birbaşa LU, böyüklərdə
BiCGStab + ILU. Ön-şərtçi hər Nyuton iterasiyasında yenilənir —
Jakobian bir addım ərzində kəskin dəyişir.

### İki konvergensiya meyarı

    CNV = max_c |R_p,c| · Δt / PV_c              < 1e-3
    MB  = |Σ_c R_p,c| · Δt / Σ_c (PV_c·S_p/B_p)  < 1e-7

Hər ikisi ölçüsüzdür — eyni tolerans bütün grid ölçülərində işləyir.
MB toleransı CNV-dən 10⁴ dəfə sıxdır: lokal səhvlər bir-birini yeyə
bilər, qlobal kütlə itkisi isə yolverilməzdir.

### Tapılan və düzəldilən üç problem

**1. Sıxılma buraxılmışdı.** İlk müqayisədə IMPES 14.22 % RF verirdi,
Nyuton isə 12.76 % — və fərq Δt kiçildikdə itmirdi. Səbəb: IMPES
sıxılmanı `ct` ilə modelləşdirir, implicit sxem isə PVT olmayanda
`B` sabit saxlayırdı. Hasilat yalnız vurulan həcmlə məhdudlaşırdı.

Həll: PVT olmasa da sıxılma modeli tətbiq olunur —
`B(p) = B_ref/(1 + c·Δp)` və `PV(p) = PV_ref·(1 + c_r·Δp)`.
Nəticə: **IMPES 14.22 %, Nyuton 14.10 %** — fərq 0.8 %, bu iki
diskretizasiyanın təbii fərqidir.

**2. Doyumluluq hədləri Nyutonu dayandırırdı.** Δt = 1 gündə CNV
1e-4-də ilişirdi. Ölçmə: qalığın **100 %-i** `Sw = Swc` həddində
ilişən 36 hüceyrədən gəlirdi, sərbəst hüceyrələrdə dəqiq sıfır idi.
Klassik bound-constrained məsələdir.

Həll: hədlərə ±0.01 ədədi zolaq. Fiziki nəticəyə təsir etmir (kr
onsuz da Corey düsturunda kəsilir); nəticə fiziki hədlərə qaytarılır.

**3. PVT törəməsi.** İlk yanaşmada `np.gradient` ilə hamar törəmə
işlədildi. Ölçmə: bu, Nyutonun iterasiya sayını dəyişmir (4/5/7 hər
iki halda), lakin Jakobianı sonlu fərqdən 10⁶ dəfə uzaqlaşdırır.
Ona görə parçalı xətti cədvəlin **dəqiq** interval meyli işlədilir.

### Nyutonun yığılma həddi

Yığılma `Δt·q/PV` nisbətindən asılıdır — bir addımda quyu
hüceyrəsindən neçə məsamə həcmi keçir:

| q/PV (1/gün) | Δt = 10 | 30 | 90 | 365 |
|---|---|---|---|---|
| 0.02 | 4 iter | 5 | 6 | 11 |
| 0.05 | 4 | 6 | 9 | yığılmır |
| 0.20 | 7 | 11 | yığılmır | yığılmır |

Praktik hədd: **Δt·q/PV ≲ 10–15**. IMPES-in CFL şərti eyni nisbəti
**0.45**-lə məhdudlaşdırır — implicit sxem təxminən **25–30 dəfə**
böyük addım ata bilir.

Quyusuz məsələ (yalnız daxili axın) **istənilən Δt-də bir
iterasiyada** yığılır.

### Ölçülmüş sürət (21×21, 900 gün, sabit Δt)

| Sxem | Addım | Vaxt | RF |
|---|---|---|---|
| IMPES | 942 | 0.71 san | 14.22 % |
| Nyuton, Δt = 10 gün | 90 | 0.40 san | 14.10 % |
| Nyuton, Δt = 30 gün | 30 | 0.29 san | 14.10 % |
| Nyuton, Δt = 90 gün | 10 | 0.18 san | 14.10 % |

Nyutonun nəticəsi Δt-dən **asılı deyil** — 10, 30, 90 gündə eyni RF.
Bu, ədədi diffuziyanın olmadığını göstərir.

### Məlum məhdudiyyət — doyma nöqtəsi

Model `Pb`-ni kəsəndə Nyuton osilyasiya edir: `∂Bo/∂p` orada işarə
dəyişir. Real black-oil simulyatorda həmin nöqtədə qaz fazası ayrılır
və primary dəyişən dəyişdirilir (variable switching) — bu, **A7**-nin
işidir. Hazırkı iki fazalı model bu rejimdə işləməməlidir;
`ReservoirModel.diagnose()` (V1) bu barədə xəbərdarlıq verir.

Məhdudiyyət test ilə sənədləşdirilib: davranış dəyişsə, test xəbər
verəcək.

---

## Mərhələ 4 — adaptiv zaman addımı

### Strategiya

    Nyuton yığılmadı   -> Δt-ni kəs, addımı TƏKRARLA
    az iterasiya       -> Δt-ni böyüt
    çox iterasiya      -> Δt-ni kiçilt (növbəti addım üçün)

Əlavə məhdudiyyət: bir addımda maksimal ΔSw. Bu, yığılma meyarı
deyil — **dəqiqlik** meyarıdır. Nyuton böyük ΔSw ilə də yığıla bilər,
lakin cəbhə bir addımda çox irəliləyəndə həll kobudlaşır.

### FullyImplicitEngine

`ISimulationEngine` implementasiyası: eyni `ReservoirModel`, eyni
`SimulationResult`, eyni provider-lər. IMPES mühərriki **silinmir** —
hər ikisi saxlanılır və interfeysdə "Hesablama sxemi" siyahısından
seçilir.

`SimulationService.with_engine()` metodu eyni provider dəsti ilə
başqa mühərrikli servis qaytarır; orijinal servis dəyişmir.

### Ölçülmüş nəticə

| Grid | Müddət | IMPES | Implicit | Sürət |
|---|---|---|---|---|
| 21×21 | 900 gün | 2 978 addım / 4.55 san | 38 / 0.38 san | **12x** |
| 31×31 | 900 gün | 2 734 / 7.30 san | 38 / 0.73 san | **10x** |
| 41×41 | 3 600 gün | 3 720 / 21.1 san | 127 / 2.85 san | **7.4x** |

Orta Δt: IMPES 0.30 gün, implicit **24 gün**.

Adaptiv nəzarətin işlədiyinin sübutu: `max_dt` 10 və 120 gün
verildikdə RF fərqi 0.3 %-dən azdır — istifadəçi Δt seçməli deyil.

Sərt hallarda (q/PV = 0.2) addım avtomatik kəsilir və təkrarlanır:
900 gün 22 addımda, cəmi 1 təkrarla.

## Mərhələ 5 — CPR ön-şərtçisi

CPR = Constrained Pressure Residual. İki mərhələli ön-şərtçi:

    M⁻¹r = P·A_p⁻¹·(W·r)  +  S⁻¹·(r − A·P·A_p⁻¹·(W·r))

Səbəb: Jakobian iki fərqli xarakterli tənliyi birləşdirir —

| Tənlik | Xarakter | Nə tələb edir |
|---|---|---|
| Təzyiq | elliptik | uzaqmənzilli, qlobal həll |
| Doyumluluq | hiperbolik | lokal, upstream istiqamətdə |

Tək ILU hər ikisini eyni cür emal edir və elliptik hissədə zəif qalır.

### Dekuplinq (quasi-IMPES)

Hər hüceyrənin 2×2 diaqonal blokundan çəkilər hesablanır ki, çəkili
cəmdə doyumluluq törəməsi **yox olsun**:

    w = [ ∂R_o/∂Sw,  −∂R_w/∂Sw ]

Sabit B halında bu, `w = [Bw, Bo]` — klassik həcm balansı.

Test ilə yoxlanılıb: qalıq doyumluluq törəməsi < 1e-10, təzyiq
matrisinin asimmetriyası < 10 % (yəni o, həqiqətən elliptik operatordur).

### İkinci mərhələ — hamarlayıcı

İlk versiyada yalnız ILU idi. Ölçmə göstərdi ki, **10 000 hüceyrədən
sonra qlobal ILU ümumiyyətlə qurula bilmir** və hamarlayıcısız CPR
yığılmır (500 iterasiya, 82 saniyə).

Həll: blok-Jakobi ehtiyatı — hər hüceyrənin 2×2 diaqonal bloku ayrıca
tərslənir. Həmişə qurulur, yaddaş tələb etmir.

### Ölçülmüş nəticə — dürüst mənzərə

**Sürətdə CPR üstün DEYİL.** Orta ölçülü sistemlərdə güclü ILU
1 iterasiyada həll edir:

| Hüceyrə | ILU | CPR (ILU hamarlayıcı) |
|---|---|---|
| 1 600 | 0.016 san, 1 iter | 0.044 san, 1 iter |
| 4 900 | 0.081 san, 1 iter | 0.217 san, 1 iter |
| 16 900 | 0.293 san, 1 iter | 0.637 san, 1 iter |

**Üstünlük yaddaşdadır.** Blok-Jakobi hamarlayıcı ilə:

| Hüceyrə | ILU nnz | CPR nnz | nisbət |
|---|---|---|---|
| 225 | 9 575 | 6 052 | 0.63 |
| 625 | 37 626 | 23 069 | 0.61 |
| 1 600 | 123 846 | 60 549 | 0.49 |
| 3 025 | 266 530 | 116 016 | **0.44** |

Nisbət ölçü ilə **yaxşılaşır**, çünki ILU-nun doldurulması superxətti
artır, CPR isə yarıya bölünmüş sistemdə faktorlaşdırma aparır.

Ekstrapolyasiya: 1 milyon hüceyrəli modeldə ILU ≈ 1.2 GB, CPR ≈ 480 MB.
Məhz orada ILU yaddaşa sığmır və CPR yeganə variant olur.

### Seçim qaydası

`NewtonLinearSolverConfig.preconditioner`:

| Dəyər | Davranış |
|---|---|
| `"ilu"` | həmişə ILU |
| `"cpr"` | həmişə CPR |
| `"auto"` (defolt) | sistem `cpr_threshold` (40 000) -dan böyükdürsə CPR |

Ön-şərtçi seçimi **həlli dəyişmir** — test ilə yoxlanılıb (nisbi fərq
< 1e-5, RF fərqi < 0.05 %). O, yalnız həllə gedən yola təsir edir.

### pyamg

`CprConfig.use_amg` mövcuddur: `pyamg` quraşdırılıbsa, təzyiq
alt-sistemi üçün AMG işlədilir. Bu mühitdə paket yoxdur, ona görə
təzyiq sistemi GMRES + ILU ilə həll olunur. AMG ilə iterasiya sayının
daha da azalması gözlənilir.


