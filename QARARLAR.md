# Qərarlar jurnalı — Lay Simulyatoru

Burada layihədə verilən **texniki qərarlar** və onların **səbəbləri**
saxlanılır. Məqsəd: 6 ay sonra "biz niyə belə etmişdik?" sualına cavab olsun.

Hər qərarın formatı:

```
## Q-NN — Başlıq
**Tarix:** ...
**Vəziyyət:** Qəbul olunub / Əvəzlənib / Ləğv olunub
**Kontekst:** hansı problem vardı
**Qərar:** nə seçdik
**Səbəb:** niyə məhz bu
**Nəticələr:** bu seçimin qiyməti nədir (mənfi tərəflər də daxil)
```

---

## Q-01 — Sənədləşdirmə dili Azərbaycan dili olsun

**Tarix:** 10 sentyabr 2026
**Vəziyyət:** Qəbul olunub

**Kontekst:** Layihə sənədləri hansı dildə yazılsın?

**Qərar:** Bütün sənədlər (`.md` faylları) Azərbaycan dilində.
Kod, dəyişən/funksiya adları, commit mesajları texniki hissədə ingiliscə.

**Səbəb:** Komandanın iş dili Azərbaycan dilidir; sənəd oxunmursa,
yazılmasının mənası yoxdur.

**Nəticələr:** Beynəlxalq əməkdaşlıq lazım olsa, tərcümə yükü yaranacaq.
Bu qəbul olunan qiymətdir.

---

## Q-02 — Sənəd şablonları boş yaradıldı

**Tarix:** 10 sentyabr 2026
**Vəziyyət:** Qəbul olunub

**Kontekst:** Layihənin texniki tələbləri hələ verilməyib, amma struktur
lazımdır.

**Qərar:** Şablonlar başlıqlarla yaradıldı, məzmun ⏳ ilə boş saxlanıldı.

**Səbəb:** Uydurma məzmun (məsələn "Python istifadə olunur" — halbuki
seçilməyib) sonradan səhv qərarların əsasına çevrilir. Boş yer dürüstdür.

**Nəticələr:** Repo ilk baxışda "yarımçıq" görünür — amma bu, real
vəziyyəti düzgün əks etdirir.

---

## Q-03 — Python 3.12 seçildi (3.11 deyil, 3.14 deyil)

**Tarix:** 10 sentyabr 2026
**Vəziyyət:** Qəbul olunub

**Kontekst:** Sahibkarın planı "Python 3.11+" tələb edirdi. Maşında
3.11 **yoxdur** — yalnız 3.14.7 (aktiv) və 3.12 var.

**Qərar:** Yeni iş üçün Python 3.12.

**Səbəb:** 3.12 "3.11+" şərtini ödəyir və bütün elmi paketlərin hazır
wheel dəstəyi var. 3.14 çox yenidir — ölçüldü: **PyKrige 3.14 üçün
qurula bilmir** (Cython ilə kompilyasiya olunur, wheel yoxdur),
3.12-də isə 1.7.3 problemsiz quruldu.

**Nəticələr:** Planda yazılan "3.11" rəqəmi ilə fərq var; sənədlərdə
3.12 kimi qeyd olunur.

---

## Q-04 — Mövcud kod bazası atılmır, üzərində davam edilir

**Tarix:** 10 sentyabr 2026
**Vəziyyət:** Qəbul olunub

**Kontekst:** Sahibkar əvvəllər yazdığı proqramı "çox qarışıq və
səhvləri var" deyə təsvir etdi və sıfırdan başlamağa hazır idi.
Audit bunu yoxladı.

**Qərar:** Sıfırdan yazılmır. `lay-simulyatoru` kod bazası əsas
götürülür.

**Səbəb:** Ölçülmüş göstəricilər əksini sübut etdi:

- 65 175 sətir kod, bunun **28 509-u test**
- **1 993 test keçir** (real icra, 7 dəq 48 san)
- 0 `TODO`, 0 `FIXME`, 0 `HACK`, 0 çılpaq `except:`
- Təmiz hexagonal arxitektura, sənədləşdirilmiş refaktorinq
- Buckley-Leverett analitik doğrulaması + şəbəkə yaxınsaması sübutu
- **MPFA-O artıq yazılıb** (1 699 sətir) — planın baş tələbi
- Corner-point, Kriging, SGS, Peaceman, IMPES, FIM, CPR — hamısı var

Sıfırdan yazmaq bu doğrulanmış fizikanı itirmək demək olardı.

**Nəticələr:** Layihənin "yeni başlanğıc" xarakteri itir; əvəzində
mövcud kodun konvensiyalarına uyğunlaşmaq lazım gələcək.

---

## Q-05 — Bu layihə üçün yeni venv yaradılmır

**Tarix:** 10 sentyabr 2026
**Vəziyyət:** Qəbul olunub

**Kontekst:** Windows 11-də Smart App Control işlək vəziyyətdədir
(`VerifiedAndReputablePolicyState = 1`) və imzasız DLL-ləri bloklayır.

**Ölçülmüş fərq:**

| Mühit | SciPy | VTK |
|---|---|---|
| Yeni venv (3.12) | `sparse.linalg`, `interpolate`, `spatial` **bloklanır** | bloklanır |
| Sistem Python 3.14 | **hamısı işləyir** | bloklanır |

**Qərar:** Mövcud kod bazası ilə iş sistem Python 3.14-də aparılır.

**Səbəb:** Yeni yazılan imzasız fayllar SAC-da reputasiya qazanmayıb,
köhnələr qazanıb. Layihə əlavə paket tələb etmir — Kriging **öz içində**
yazılıb (`variogram.py` 940 sətir + `interpolation.py` 1 147 sətir),
PyKrige-a ehtiyac yoxdur.

**Nəticələr:** Asılılıqların izolyasiyası itir; paket versiyaları
qlobal mühitə bağlı qalır. Yeni paket əlavə etmək lazım gələrsə,
SAC yenidən problem çıxara bilər.

---

## Q-06 — Smart App Control söndürülmür

**Tarix:** 10 sentyabr 2026
**Vəziyyət:** Qəbul olunub

**Kontekst:** SAC-ı söndürmək VTK blokunu dərhal həll edərdi.

**Qərar:** Söndürülmür.

**Səbəb:** SAC-ın söndürülməsi **geri qaytarılmır** — yenidən
aktivləşdirmək üçün Windows-u tam yenidən qurmaq lazımdır. Bu maşın
ümumi iş kompüteridir, ayrılmış hesablama serveri deyil. Aylar sonra
lazım olacaq bir imkan üçün sistem müdafiəsini daimi zəiflətmək
düzgün deyil.

**Alternativlər:** vizualizasiya üçün matplotlib/Plotly (bloklanmır),
və ya WSL2 (SAC tətbiq olunmur, **geri qaytarıla bilir**).

**Nəticələr:** `rendering/vtk_volume.py` (895 sətir) bu maşında
işləmir. Vizualizasiya qatı üçün ayrıca qərar lazımdır.

⏳ *Yekun vizualizasiya qərarı hələ verilməyib.*

---

## Q-07 — VTK bloku öz-özünə həll oldu; Q-06 doğrulandı

**Tarix:** 10 sentyabr 2026
**Vəziyyət:** Qəbul olunub

**Kontekst:** Q-06-da Smart App Control-un söndürülməməsi qərara
alınmışdı. Həmin gün, bir neçə saat sonra VTK bloku **öz-özünə aradan
qalxdı**.

**Ölçülmüş:**

| Mühit | Əvvəl | Sonra |
|---|---|---|
| Sistem Python 3.14 — SciPy | ✅ | ✅ |
| Sistem Python 3.14 — **VTK** | ❌ bloklanırdı | ✅ **işləyir** |
| venv 3.12 — SciPy, PyKrige | ❌ | ✅ |
| venv 3.12 — VTK, PyVista | ❌ | ❌ |
| SAC siyasəti | `= 1` | `= 1` (**dəyişməyib**) |

Test dəsti: **2 042 keçdi** (əvvəl 1 993) — fərq məhz 49 VTK testidir.

**Qərar:** Vizualizasiya qatı olduğu kimi qalır — VTK + PyQt5.
Plotly-yə keçid, WSL2 və SAC-ın söndürülməsi variantları **ləğv edilir**.

**Səbəb:** SAC imzasız faylların reputasiyasını buludda yoxlayır; sorğu
tamamlanana qədər bloklayır. Problem daimi deyilmiş.

**Nəticələr:** Q-06 (SAC söndürülməsin) doğru çıxdı — geri qaytarılmayan
qərar verilsəydi, heç bir faydası olmayacaqdı. Yeni paket quraşdıranda
müvəqqəti blok yenidən görünə bilər; bu, gözlənilən davranışdır.

**Dərs:** mühit problemi diaqnoz edəndə "daimi maneə" nəticəsinə tələsik
gəlmək olmaz — ölçmə təkrarlanmalıdır.


---

## Q-08 — Sahibkarın plan stəkinin qiymətləndirilməsi

**Tarix:** 10 sentyabr 2026
**Vəziyyət:** ⏳ **Təklif — sahibkarın təsdiqi gözlənilir**

**Kontekst:** Sahibkarın verdiyi 3D 3-fazalı plan konkret texnologiyalar
adlandırır. Plan sətir-sətir kod bazası ilə tutuşduruldu (Seans 3).

**Təklif olunan qərarlar:**

| Planda | Təklif | Səbəb (ölçülmüş) |
|---|---|---|
| Python **3.11** | ❌ **qəbul edilməsin** | 3.14.7 işləyir, hər şey (VTK daxil) keçir; enmək geriyə addımdır |
| **RBFInterpolator** | ❌ **qəbul edilməsin** | öz kriginqimiz var (`geology/`, 8 253 sətir) və qiymətləndirmə **variance**-ı da verir — RBF vermir |
| **PyKrige** | ❌ **qəbul edilməsin** | kriginq öz içimizdə yazılıb; artıq xarici asılılıq SAC riski deməkdir |
| **`spsolve` / SuperLU / AMG** | ❌ **qəbul edilməsin** | bizdə CG+ILU və **CPR ön-şərtçisi** var; birbaşa həll böyük modeldə dayanır |
| **PyVista** | ⚠️ **saxlanılsın: mövcud VTK** | PyVista elə eyni VTK-nın örtüyüdür; `vtk_volume.py` (895 sətir, 49 test) işləyir. Yeganə qazanc — hazır slice/volume API; onu VTK ilə də yazmaq olar |
| **Plotly** | ❌ **qəbul edilməsin** | Plotly veb üçündür, bizim proqram masaüstü PyQt5 tətbiqidir (bax Q-07) |
| **Plandakı `src/` fayl strukturu** | ❌ **qəbul edilməsin** | düz modul yığınıdır; bizdə heksaqonal (onion) arxitektura var — modul ADLARI onsuz da uyğun gəlir |
| **THP / VFP modulu** | ✅ **qəbul edilsin** | kodda ümumiyyətlə yoxdur; planın ən dəyərli bəndi |
| Slice plane · animasiya · GIF · Pcog · CSV ixracı | ✅ **qəbul edilsin** | real boşluqlardır |

**Nəticələr:** Plan **yol xəritəsi kimi deyil, boşluq siyahısı kimi**
işlədilir. `ROADMAP.md`-dəki mərhələ cədvəlləri buna görə real fayl
yolları və ölçülmüş statuslarla dolduruldu.

⏳ *Sahibkar bu təklifi təsdiqləyənə qədər qərar qüvvəyə minmir.*

---

## Q-09 — Üç fazalı Bo doymamış qoldan hesablanır; ölü-neft düzəlişi TƏTBİQ EDİLMİR

**Tarix:** 11 sentyabr 2026 · **Kontekst:** B3-B (bax
[ISH_HESABATI.md](ISH_HESABATI.md) → Seans 10)

### Sual

Üç fazalı mühərrik Pb ≥ 240 bar-da yığılmırdı. İki yol vardı:

1. B3-A-nın **ölü-neft düzəlişini** qaz sütunlu cədvələ də tətbiq etmək;
2. **doymamış Bo qolunu** düzgün hesablamaq.

### Qərar — 2-ci yol

`Bo(p, Rs) = Bo_sat(Pb(Rs)) · exp(c_o · (Pb(Rs) − p))`

### Niyə 1-ci yol RƏDD EDİLDİ

Ölçüldü: məcburi ölü-neft düzəlişi də yığılmanı bərpa edir (Pb=240 →
23 addım). **Lakin nəticə fiziki olaraq YANLIŞ olardı** — ölü-neft Bo-su
"neftdən qaz ayrılmır, tərkib sabitdir" deməkdir. İki fazalı modeldə bu
qəbuledilən yamaqdır, çünki model ayrılan qazı onsuz da saymır. Üç
fazalı modeldə isə qaz tənliyi məhz o qazı sayır — Bo-nu ölü-neftə
çevirmək neftin şişməsini silib **qaz kütlə balansını pozardı**.

### Niyə 2-ci yol DÜZGÜNDÜR

Problem heç vaxt "Bo-nun düzəldilməli olması" deyildi — Bo **YANLIŞ
QOLDAN oxunurdu**. Doymamış hüceyrədə neftin tərkibi sabitdir, ona görə
Bo həmin Rs-in sıxılma qoluna aiddir. Düzgün qol:

* `dBo/dp < 0` hər yerdə → Jakobian kilidlənmir;
* `∂Bo/∂Rs ≠ 0` → neft tənliyi 3-cü dəyişənə bağlanır (kodda yazılmış
  "qaz tənliyi kompensasiya edir" fərziyyəsi məhz bu əlaqənin
  olmamasına görə SƏHV idi);
* `Rs` toxunulmur → **qaz kütlə balansı pozulmur**;
* doymuş hüceyrədə `Pb(Rs) = p` → keçid kəsilməz;
* Pb-dən yuxarı cədvəlin öz düsturudur → köhnə nəticələr dəyişmir.

Bu, Eclipse `PVTO`-nun və sənaye black-oil simulyatorlarının standart
davranışıdır — yeni bir şey icad edilmir.

### Nəticə

İki fazalı yol (`dead_oil_below_bubble_point`) **olduğu kimi qalır** —
o, ayrı problemin ayrı həllidir. Bax Q-08 və B3-A.

---

## Q-10 — `run.bat` əvvəlcə `.venv`-ə baxır

**Tarix:** 11 sentyabr 2026
**Vəziyyət:** Qəbul olunub

**Kontekst:** `run.bat` virtual mühiti iki yerdə axtarırdı: `..\venv`
və `venv`. Layihənin faktiki venv-i isə repo kökündə **`.venv`** adı
ilə yerləşir — yəni skript onu heç vaxt tapmırdı və hər dəfə
"XETA: venv tapilmadi" verirdi. Skript işlək deyildi.

**Qərar:** `.venv\Scripts\activate.bat` axtarış siyahısına **birinci**
sırada əlavə edildi. Köhnə iki yol olduğu kimi saxlanıldı.

**Səbəb:** `.venv` — Python alətlərinin (uv, VS Code, PyCharm) standart
adıdır və `.gitignore`-da artıq var. Birinci sırada olmalıdır ki, kökdə
həm `.venv`, həm də köhnə `venv` varsa, aktual olan seçilsin. Köhnə
yolların silinməməsi isə başqa maşınlarda qurulmuş mühitləri
sındırmamaq üçündür.

**Nəticələr:** `run.bat` indi standart quraşdırmada işləyir. Xəta
mesajındakı gözlənilən yerlər siyahısına da `.venv` əlavə olundu.

**Q-05 ilə münasibət:** Q-05 "bu layihə üçün yeni venv yaradılmır"
demişdi — səbəb SAC blokları idi. Q-07-də blokun aradan qalxdığı
qeydə alınıb və `.venv` (3.12) altında proqram problemsiz işləyir
(ölçüldü, 11 sentyabr). Q-05 **ləğv edilmir** — bu dəyişiklik yeni
venv yaratmır, yalnız mövcud olanı tanıyır.

**Qeyd:** skriptin özü `python app.py` çağırır — konsol pəncərəsi açıq
qalır və `errorlevel` yoxlanılır. Bu qəsdəndir: xəta olanda `pause`
mesajı göstərsin. Konsolsuz işə salma üçün `pythonw.exe` istifadə
olunur (bax `ISH_HESABATI.md` → Seans 11, Tapıntı 2).

---

## Q-11 — THP traversi çoxseqmentlidir; Hagedorn-Brown İŞLƏDİLMİR

**Tarix:** 11 sentyabr 2026 · **Kontekst:** B4, A variantı
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 12)

### Qərar 1 — tək seqment yox, N seqment (defolt 20)

`ICRA_PLANI.md`-dəki ilkin B4 eskizi təzyiq düşgüsünü **tək seqmentdə**
(orta təzyiqdə) hesablamağı nəzərdə tuturdu.

**Rədd edildi.** Səbəb ölçülüb: B4b-dən sonra modeldə real sərbəst qaz
var (maks S_g 0.11), qaz isə yuxarı qalxdıqca genişlənir — `Bg` böyüyür,
ρ azalır, v artır. Qazlı quyuda (GOR 150, 2000 m):

| Seqment | THP, bar | 200-dən fərq |
|---|---|---|
| 1 | 115.58 | **2.77 bar** |
| 5 | 118.20 | 0.15 bar |
| 20 | 118.34 | 0.006 bar |
| 200 | 118.35 | — |

Tək seqment yalnız qazsız quyuda düzgündür. Defolt **20** seçildi.

### Qərar 2 — Hagedorn-Brown DEYİL

Sahibkar HB-ni adı ilə soruşdu. Rədd səbəbi **texnikidir, fiziki deyil**:
HB korrelyasiyası `CNL` və `ψ` **qrafiklərinin** rəqəmsallaşdırılmasını
tələb edir. O cədvəl datası bizdə yoxdur; uyğunlaşdırılmış əyriləri
uydursaydım, nəticə "Hagedorn-Brown" adını daşıyıb HB OLMAZDI.

Sahibkar razılaşdı: v1-də **sürüşmə yoxdur** (`H_L = λ_L`).
Sürüşmə lazım olanda **Beggs-Brill** tövsiyə olunur — qapalı
düsturlarla yazılır, testlənə bilir və `IHoldupCorrelation`
interfeysinə çağıranı dəyişmədən oturur.

⚠️ **Sənədləşən nəticə:** no-slip qazlı quyuda hidrostatik sütunu
OLDUĞUNDAN YÜNGÜL göstərir, yəni hesablanan THP həqiqi dəyərdən
YÜKSƏKDİR. Fərq GOR artdıqca böyüyür.

### Qərar 3 — axa bilməyən quyuda `nan`, sıfır YOX

Sütunun çəkisi BHP-ni üstələyəndə traverse `nan` qaytarır. Sıfıra
"qısaldılmış" dəyər qrafikdə REAL ÖLÇMƏ kimi görünərdi — bu, səssiz
yanlış məlumat olardı. `nan` qrafikdə boşluq buraxır və mesajda səbəb
yazılır ("süni qaldırma lazımdır").

### Qərar 4 — THP post-prosesdir, mühərrik toxunulmur

A variantında quyu dibi təzyiqi onsuz da sabitdir (istifadəçinin
hədəfi), debit sıraları isə `SimulationResult`-da var. Ona görə THP
simulyasiyadan SONRA hesablanır — Nyutonun yığılmasına sıfır risk.
Test bunu kilidləyir: lüləli və lüləsiz qaçışın RF-i bitə-bit eynidir.

`ControlMode.THP` (B variantı) gələndə əlaqə addımın əvvəlinə köçəcək,
lakin `simulation/wellbore/` paketi YENİDƏN İŞLƏDİLƏCƏK — təkrar
yazılmayacaq.

---

## Q-12 — İşləməyən panel sahələri bozarılır, keş isə silinmir

**Tarix:** 11 sentyabr 2026 · **Kontekst:** "parametrlər təsir etmir"
bildirişi (bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 14)

### Sual

İstifadəçi paneldə dəyər dəyişirdi, nəticə isə dəyişmirdi. Üç səbəb
tapıldı və hər üçü **dizayn üzrə** idi:

* PVT modeli işlədiləndə μo `PVT cədvəlindən` gəlir;
* interpolyasiya olunmuş geologiya keşlənir;
* GRDECL idxalında φ/K faylın xəritələrindən gəlir.

Yəni fizika səhv deyildi — **davranış görünmürdü**.

### Qərar 1 — hesablama DEYİL, görünürlük düzəldilir

Mühərriyə toxunulmadı. Bunun əvəzinə işləməyən sahələr söndürülür
(`setEnabled(False)`) və səbəb yazılır.

**Alternativ rədd edildi:** "panel dəyəri PVT cədvəlini üstələsin".
Bu, black-oil modelini pozardı — cədvəldəki μo(p) təzyiqdən asılıdır,
paneldəki isə sabitdir; ikisini qarışdırmaq termodinamik olaraq
uyğunsuz nəticə verərdi.

### Qərar 2 — keş SİLİNMİR, yalnız köhnəlmiş işarələnir

`_geology_model_from_wells` keşinin öz məqsədi var: böyük gridi hər
klikdə yenidən interpolyasiya etməmək (funksiyanın sənədi bunu açıq
yazır). Avtomatik silmək həmin niyyəti pozardı.

Ona görə φ/K dəyişəndə yalnız **mövcud banner** ("Nəticə köhnəlib —
'İnterpolyasiya et' basın") işə düşür. Qərar istifadəçinindir.

### Qərar 3 — ayrıca `geology_changed` siqnalı

`RockFluidPanel.changed` bütün sahələr üçün atəş açır. Lözlük
dəyişəndə geologiyanı köhnəlmiş saymaq YANLIŞ siqnal olardı, ona görə
yalnız φ/K/heterogenlik sahələri ayrıca siqnal verir.

### Qərar 4 — bozarma dəyəri SİLMİR

Söndürülmüş widget öz dəyərini saxlayır və `fluids()` onu oxumağa
davam edir. Beləliklə PVT söndürüləndə istifadəçinin əvvəlki lözlük
dəyəri itmir. Ayrıca testlə kilidləndi.

---

## Q-13 — Əsas maşında `.venv` `--system-site-packages` ilə qurulur

**Tarix:** 12 sentyabr 2026 · **Kontekst:** sahibkar `run.bat`-ın bu
maşında işləməsi üçün `.venv` qurulmasını seçdi
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 15)

### Sual

`C:\Dev\LSM`-də `.venv` yox idi, ona görə `run.bat` (Q-10) venv
tapmayıb dayanırdı. Adi yol — `python -m venv .venv` + `pip install
-r requirements-dev.txt` — sınandı və **uğursuz oldu**.

### Tapıntı — Windows Application Control təzə DLL-ləri bloklayır

`.venv`-ə pip ilə endirilən təkərlərin imzasız yerli kitabxanaları
işə salınanda:

```
ImportError: DLL load failed while importing vtkCommonCore:
An Application Control policy has blocked this file.
```

Eyni paketin **sistem Python-undakı nüsxəsi işləyir** (2 322 test
keçir, VTK-nın 49 testi daxil). Yəni maneə paketin özündə deyil,
**yeni endirilmiş fayllarda**dır — bu, ROADMAP 0.7-də qeyd olunan
Smart App Control maneəsinin eyni ailəsidir.

Müşahidə: `scipy` bir neçə dəqiqə sonra keçdi, `vtk` isə bloklanmış
qaldı. Yəni siyasət vaxta görə dəyişir və **etibarlı deyil**.

### Qərar

`.venv` **öz paket kopyaları olmadan**, sistem paketlərinə baxan
nazik mühit kimi qurulur:

```
python -m venv --system-site-packages .venv
```

Beləliklə:

* `run.bat` venv tapır → Q-10 dəyişmədən işləyir;
* idxallar artıq **etibarlı sayılan** sistem DLL-lərinə düşür;
* heç bir versiya fərqi yaranmır.

Bu təhlükəsizdir, çünki sistem Python-u tələb olunan **18 paketin
hamısını dəqiq pinlənmiş versiyalarda** daşıyır — `requirements.txt`
və `requirements-dev.txt` ilə bir-bir tutuşduruldu, uyğunsuzluq 0.

### Rədd edilən alternativlər

| Alternativ | Niyə yox |
|---|---|
| Adi (izolyasiya olunmuş) `.venv` | Application Control bloklayır — proqram ümumiyyətlə açılmır |
| `Unblock-File` / siyasəti söndürmək | Maşının təhlükəsizlik siyasətinə müdaxilədir; sahibkardan belə icazə istənilməyib |
| `run.bat`-a "venv yoxdursa sistem Python-u" qolu | Sahibkar `.venv` variantını seçdi; həm də bu, Q-10-un məqsədini (mühitin müəyyənliyi) zəiflədərdi |

### Nəticə

`.venv` bu maşında `--system-site-packages` ilə qurulub və
`run.bat` ilə proqram AÇILIR (Seans 15-də gözlə təsdiqləndi).
Yeni asılılıq əlavə olunanda o, **sistem Python-una** qurulmalıdır
(`pip install --user` və ya sistem interpretatoru ilə) — yoxsa eyni
blok təkrarlanar.

**Diqqət:** bu qərar YALNIZ bu maşına aiddir. Application Control
maneəsi olmayan maşında adi izolyasiya olunmuş venv üstündür.

---

## Q-14 — İxrac formatının qaydaları: BOM, kvadrat mötərizə, boş xana

**Tarix:** 12 sentyabr 2026 · **Kontekst:** B5-a
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 16)

### Qərar 1 — vahid KVADRAT MÖTƏRİZƏDƏ, vergüllə deyil

İlk versiya `"t, gün"` yazırdı. Bu, `UNITS.md` prinsipinə ("vahid
başlıqda olsun") uyğun idi, LAKİN vergül CSV-nin öz ayırıcısıdır —
nəticədə hər başlıq dırnağa düşürdü və `awk`/`cut` kimi sadə alətlər
faylı sındırırdı.

İndi: `t [gün]`, `RF [%]`, `GOR [sm³/sm³]`.

⚠️ **Bu qüsuru test TAPMADI** — `csv` modulu dırnaqlanmış sahəni
düzgün oxuyur, ona görə gedər-gələr testi keçirdi. Qüsur **faylı açıb
gözlə baxanda** göründü. Dərs: format qərarlarında test kifayət deyil,
çıxışa baxmaq lazımdır.

### Qərar 2 — CSV `utf-8-sig` (BOM ilə)

Windows Excel BOM-suz UTF-8-i tanımır: `ə`, `ş`, `ğ` korlanır. BOM
əlavə olunur; `csv` və `pandas` onu özləri atır, yəni gedər-gələr
pozulmur.

**Alternativ rədd edildi:** "ASCII başlıqlar işlət" — sənədlərin
hamısı Azərbaycan dilindədir, ixracı ingiliscəyə çevirmək uyğunsuzluq
yaradardı.

### Qərar 3 — `nan` BOŞ xana, JSON-da `null`

CSV-də boş xana, JSON-da `null`. **`0` QƏTİYYƏN YOX.**

Səbəb B4-A-dan gəlir: `well_thp`-də quyunun səthə axa bilmədiyi
addımlar `nan`-dır və bu, qəsdən belədir (bax
[QARARLAR.md](QARARLAR.md) → Q-11). `0 bar` yazsaydıq, o, real ölçmə
kimi oxunardı — eyni səhv, sadəcə qrafik əvəzinə faylda.

`json.dump` defolt `NaN` yazır, lakin bu, JSON spesifikasiyasına
ziddir — `allow_nan=False` ilə qadağan olundu.

### Qərar 4 — metadata yalnız JSON-da

CSV TƏMİZ cədvəldir: başlıq + sətirlər, şərh sətri yoxdur. `#` ilə
başlayan sətirlər ciddi CSV oxuyucularını sındırır və gedər-gələr
müqaviləsini mənasızlaşdırır. Model adı, OOIP/OGIP, yığılma vəziyyəti
JSON-dadır.

### Qərar 5 — rəqəmlər tam dəqiqliklə

`repr(float)` işlədilir (`8.952477630417802e-05`), yuvarlaqlaşdırma
YOXDUR. İxrac faylı hesabat deyil, **məlumatdır**; yuvarlaqlaşdırma
"geri oxunanda eyni ədədlər" müqaviləsini pozardı. İnsan üçün
oxunaqlı təqdimat PDF hesabatın işidir.
