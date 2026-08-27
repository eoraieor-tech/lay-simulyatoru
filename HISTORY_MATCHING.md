# Tarixçə uyğunluğu (C5)

Simulyator proqnoz verir, lakin proqnozun dəyəri keçmişi nə qədər
yaxşı təkrarladığından asılıdır. Bu modul həmin müqayisəni ölçür və
uyğunlaşdırma parametrlərini təyin edir.

| Mərhələ | Vəziyyət |
|---|---|
| 1 · Müşahidə məlumatı, uyğunsuzluq ölçüsü, qrafiklər | **HAZIR** |
| 2 · Parametr modifikatorları | **HAZIR** |
| 3 · Avtomatik optimallaşdırma | **HAZIR** |

---

## Mərhələ 1 — müşahidə və uyğunsuzluq

### Müşahidə faylı

Uzun (long) cədvəl formatı:

    time,well,quantity,value
    30,PROD-1,OIL_RATE,142.5
    60,PROD-1,OIL_RATE,138.0
    30,,AVERAGE_PRESSURE,247.8

Boş `well` sahəsi **yataq səviyyəsi** deməkdir.

| Kəmiyyət | Tanınan adlar |
|---|---|
| Neft debiti | `OIL_RATE`, `WOPR`, `FOPR`, `QO`, `NEFT` |
| Su debiti | `WATER_RATE`, `WWPR`, `FWPR`, `QW`, `SU` |
| Su vurma | `WATER_INJECTION`, `WWIR`, `FWIR`, `VURMA` |
| Sulaşma | `WATER_CUT`, `WWCT`, `FWCT`, `SULASMA` |
| Quyudibi təzyiq | `BHP`, `WBHP` |
| Orta lay təzyiqi | `AVERAGE_PRESSURE`, `FPR` |
| Kumulyativ neft | `CUM_OIL`, `WOPT`, `FOPT` |

Eclipse SUMMARY adları qəsdən dəstəklənir — real məlumat adətən
həmin formatda ixrac olunur.

### Uyğunsuzluq ölçüsü

    SSE   = Σ (hesablanmış − ölçülmüş)²
    RMSE  = √(SSE / n)                      ölçü vahidində
    NRMSE = RMSE / (müşahidənin diapazonu)  ÖLÇÜSÜZ

Yekun uyğunsuzluq — NRMSE-lərin çəkili ortası.

**Niyə ölçüsüz.** Təzyiq barla (200–250), debit m³/günlə (0–150)
ölçülür. Xam SSE-ləri toplasaydıq, böyük ədədli kəmiyyət yekun ölçüyə
hakim olardı və optimallaşdırma yalnız onu uyğunlaşdırardı.

Test bunu qoruyur: dəyərlər 1000 dəfə böyüdüləndə RMSE 1000 dəfə
artır, NRMSE isə **dəyişmir**.

### Çəkilər

| Kəmiyyət | Çəki |
|---|---|
| Kumulyativ neft | 1.5 |
| Neft/su debiti, sulaşma | 1.0 |
| Təzyiq | 0.7 |
| Su vurma | 0.5 |

Kumulyativ neft daha ağırdır: gündəlik debitdə səs-küy çox olur,
toplam isə ehtiyatın nə qədər çıxarıldığını birbaşa göstərir.

### Əlavə göstəricilər

| Göstərici | Nə deyir |
|---|---|
| **Meyl** (bias) | Sistematik səhv: müsbət = model çox verir |
| **Korrelyasiya** | Formanın uyğunluğu (miqyasdan asılı deyil) |

RMSE işarəni gizlədir — model həmişə 5 vahid çox verirsə və həmişə
5 vahid az verirsə, RMSE eynidir. Meyl bu ikisini ayırır.

### Zaman uyğunlaşdırması

Model nəticəsi müşahidə vaxtlarına **interpolyasiya olunur** — əksinə
yox, çünki müşahidə həqiqətdir və dəyişdirilməməlidir.

Simulyasiya dövrünə düşməyən müşahidələr **atlanır və hesabatda
göstərilir**. Susmaq təhlükəlidir: istifadəçi məlumatın nəzərə
alındığını düşünərdi.

### Qrafiklər

**Zaman qrafikləri** — ən pis dörd sıra. Müşahidə nöqtə ilə, model
xətt ilə (sənaye konvensiyası).

**Çarpaz qrafik** — ölçülmüş/hesablanmış. İdeal uyğunluqda bütün
nöqtələr 45° xəttindədir. Səpələnmə təsadüfi xətanı, xəttdən
sistematik meyl isə model səhvini göstərir; zaman qrafikində bu
ikisini ayırmaq çətindir.

### Ölçülmüş həssaslıq

| Model | Yekun uyğunsuzluq |
|---|---|
| Doğru | 0.068 |
| Keçiricilik 0.3× | 0.922 |

13 dəfə fərq — ölçü səhv modeli aydın seçir.

---

## Mərhələ 2 — uyğunlaşdırma parametrləri

| Parametr | Tip | Hədlər | Nə idarə edir |
|---|---|---|---|
| `PERM_MULT` | çarpan, **log** | 0.1 – 10 | Keçiriciliyin qlobal səviyyəsi |
| `PORO_MULT` | çarpan | 0.7 – 1.4 | Ehtiyat |
| `KV_KH` | mütləq, **log** | 0.001 – 1 | Şaquli axın, cazibə seqreqasiyası |
| `SOR` | mütləq | 0.05 – 0.45 | Qalıq neft, RF-in yuxarı həddi |
| `SWC` | mütləq | 0.05 – 0.40 | Bağlı su |
| `KRW_END` | mütləq | 0.05 – 1.0 | Su cəbhəsinin sürəti |
| `COREY_NW/NO` | mütləq | 1 – 6 | Cəbhənin forması |
| `MU_OIL` | mütləq, **log** | 0.3 – 50 cP | Lözlük (PVT olmadıqda) |
| `OWC` | mütləq | model dərinliyi | Su-neft kontaktı |

Dəst modeldən asılı olaraq qurulur: `PERMZ` yoxdursa `KV_KH`
əlavə olunmur, PVT cədvəli varsa `MU_OIL` çıxarılır.

### Çarpan və mütləq

Keçiricilik və məsaməlilik üçün **çarpan** işlədilir: geoloji model
sahənin FORMASINI verir (quyulardan interpolyasiya ilə), mütləq
səviyyəsi isə qeyri-müəyyəndir. Çarpan heterogenliyi qoruyur.

SCAL parametrləri skalyardır və birbaşa təyin olunur.

### Loqarifmik miqyas

Keçiricilik çarpanı üçün vacibdir: 0.5 və 2.0 fiziki cəhətdən
simmetrik dəyişikliklərdir, xətti miqyasda isə 1.0-dan məsafələri
fərqlidir (0.5 və 1.0) və optimallaşdırıcı azaltmağa üstünlük verərdi.

Bütün parametrlər optimallaşdırıcıya **[0, 1] fəzasında** verilir ki,
müxtəlif vahidlər və diapazonlar bir-birinə mane olmasın.

### Dəyişməzlik

`ModelModifier.apply()` baza modelini **heç vaxt dəyişdirmir**, kopiya
qaytarır. Əks halda təkrar tətbiqdə çarpanlar üst-üstə yığılardı:
1.4, sonra 1.96, sonra 2.74. Optimallaşdırma yüzlərlə dəfə tətbiq edir.

---

## Yol boyu tapılan iki qüsur

### 1. Provider modeldən oxunmurdu

`SimulationService` nisbi keçiricilik adapterini konstruktorda alıb
saxlayırdı. Optimallaşdırıcı `Sor`-u dəyişəndə model yenilənirdi,
lakin adapter köhnə `Sor` ilə qalırdı və **nəticə heç dəyişmirdi**.

Ölçmə: `Sor` 0.25-dən 0.40-a dəyişikliyi RF-ə **sıfır** təsir edirdi.
Düzəlişdən sonra RF 53.2 %-dən 41.3 %-ə düşür, su gəlişi 800-dən
582 günə çəkilir.

`ModelAwareSimulationService` artıq `application` qatındadır və
defolt davranışdır. Əvvəl yalnız `app.py`-də gizli alt-sinif idi,
ona görə testlər və gələcək modullar köhnə davranışı işlədirdi.

### 2. `Swc` dəyişəndə model yararsız olurdu

Bağlı su doyumluluğu artanda ilkin `Sw` hərəkətli intervaldan kənarda
qalır və model yoxlamadan keçmir. Optimallaşdırmada bu ölümcüldür,
axtarış ilk addımda dayanardı.

İndi ilkin doyumluluq avtomatik yeni intervala salınır. Bu, həm də
fiziki cəhətdən doğrudur: bağlı su doyumluluğu dəyişəndə lay şəraiti
də dəyişir.

---

## Parametrlərin ölçülmüş təsiri

1 500 günlük ssenari, su gəlişi 800-cü gündə:

| Dəyişiklik | RF | Su gəlişi |
|---|---|---|
| Baza | 53.2 % | 800 gün |
| `SOR` 0.25 → 0.40 | 41.3 % | 582 gün |
| `KRW_END` 0.35 → 0.90 | 46.6 % | 626 gün |
| `COREY_NW` 2.5 → 5.0 | 57.6 % | 986 gün |

**Vacib qeyd:** SCAL parametrləri yalnız **su gəlişindən sonra**
özünü göstərir. Qısa ssenaridə (400 gün, su gəlişi yox) onlar nəticəyə
təsir etmir. Bu, fiziki reallıqdır, qüsur deyil: uyğunlaşdırma dövrü
su gəlişini əhatə etməlidir, əks halda SCAL parametrləri təyin oluna
bilməz.

---

## Mərhələ 3 — avtomatik optimallaşdırma

    unit vektor -> model -> simulyasiya -> uyğunsuzluq -> skalyar

Axtarış [0, 1] fəzasında aparılır. Üç üsul: **Nelder-Mead** (defolt),
**Powell**, **Differential Evolution**.

### Twin experiment — əsas doğrulama

Modelin parametrləri məlum dəyərlərə qoyulur, nəticə "müşahidə" kimi
işlədilir, sonra optimallaşdırıcının həmin dəyərləri bərpa edib-etmədiyi
yoxlanılır.

Gizlədilmiş: `PERM_MULT = 2.5`, `SOR = 0.32`

| Üsul | Bərpa olunan | Uyğunsuzluq |
|---|---|---|
| Nelder-Mead | 2.501 / 0.320 | 0.830 → 0.0001 |
| Powell | 2.498 / 0.323 | 0.830 → 0.0019 |

### Uğursuz qiymətləndirmə

Optimallaşdırıcı hədlərin kənarını da sınayır və oradakı model
yığılmaya bilər. Belə hal **istisna atmamalıdır** — axtarış dayanardı.
Əvəzinə böyük, lakin **sonlu** cərimə qaytarılır.

Sonsuzluq yaramaz: Nelder-Mead simpleksi qura bilmir, differential
evolution isə seçim apara bilmir.

### Keşləmə

Nelder-Mead eyni nöqtəni təkrar sınaya bilir. Hər qiymətləndirmə bir
simulyasiyadır, ona görə nəticələr keşlənir.

### Başlanğıc simpleks

`scipy`-ın defolt simpleksi **nisbi** addım işlədir və [0, 1] fəzasında
sıfıra yaxın parametrlər üçün praktiki olaraq hərəkətsiz qalır. Ona
görə simpleks əl ilə qurulur: hər təpə bir parametr üzrə 0.15
sürüşdürülür.

---

## Təyin olunma qabiliyyəti (identifiability)

Ən vacib praktik nəticə: **müşahidə dəsti parametri təyin edə
bilməyə bilər.**

Səs-küylü məlumatla (5–7 %) eyni məsələ iki fərqli cavab verdi:

| Üsul | Tapılan `PERM_MULT` | Uyğunsuzluq |
|---|---|---|
| Nelder-Mead | 1.29 | 0.031 |
| Differential Evolution | 5.23 | 0.037 |

Dörd dəfə fərqli keçiricilik, demək olar eyni uyğunsuzluq. Bu, alqoritm
səhvi deyil — məsələnin özü belədir.

### Səbəb

`RATE` ilə idarə olunan vurucuda neft debiti vurulan həcmlə
müəyyənləşir və keçiricilikdən demək olar asılı deyil. Yalnız **təzyiq**
keçiriciliyə həssasdır.

Xəta funksiyasının kəskinliyi (`SOR` sabit, `PERM_MULT` süpürülür):

| `PERM_MULT` | Təzyiqsiz | Təzyiqli |
|---|---|---|
| 1.0 | 0.0047 | **0.6201** |
| 2.5 (həqiqi) | 0.0000 | 0.0000 |
| 6.0 | 0.0026 | 0.1968 |

Təzyiqsiz halda funksiya **yastıdır** və 3 %-lik səs-küy siqnalı
basdırır. Təzyiq əlavə olunanda 130 dəfə kəskinləşir.

Diqqət: siqnal **aşağı keçiricilik tərəfində** güclüdür. Yuxarıda
təzyiq düşməsi onsuz da kiçikdir, ona görə fərq azalır.

### Praktik nəticə

Uyğunlaşdırmaya başlamazdan əvvəl soruşulmalıdır: **bu müşahidə dəsti
seçdiyim parametrləri təyin edə bilirmi?** Cavab yoxdursa,
optimallaşdırıcı "yaxşı" uyğunluq tapacaq, lakin parametr dəyərləri
mənasız olacaq.

Eyni səbəbdən SCAL parametrləri yalnız su gəlişindən sonrakı dövrü
əhatə edən məlumatla təyin oluna bilər.

---

## İnterfeys

**"Uyğunlaşdırma"** tabı: üsul seçimi, büdcə, başlat/dayandır,
parametr cədvəli, konvergensiya əyrisi və parametr trayektoriyası.

Axtarış fon axınında gedir — interfeys bloklanmır. **"Nəticəni modelə
tətbiq et"** düyməsi tapılan dəyərləri panellərə yazır.

Konvergensiya qrafiki loqarifmik miqyasdadır; uğursuz qiymətləndirmələr
qırmızı × ilə ayrıca göstərilir.
