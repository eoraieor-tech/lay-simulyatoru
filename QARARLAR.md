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


## Q-15 — THP idarəsi açıq birləşmədir; THP quyusu mühərrikdə BHP bağlantısıdır

**Tarix:** 13 sentyabr 2026 · **Kontekst:** B4-B
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 23)

### Qərar 1 — `ControlMode.THP` bağlantıya ÖTÜRÜLMÜR

Mühərrikdə quyu rejimi `connection.mode is ControlMode.BHP` ilə
yoxlanılır — `else` budağı RATE-dir. Bu yoxlama altı faylda var. Yeni
rejim bağlantıya ötürülsəydi, hamısı onu **səssizcə RATE** kimi işlədərdi.

İndi THP quyusu `mode = BHP`, `thp_target = THP` ilə qurulur;
`ThpController` bağlantının `target`-ini addım-addım yeniləyir.

**Alternativ rədd edildi:** hər `is BHP` yoxlamasını `in (BHP, THP)` etmək —
altı faylda, qalıq və Jakobian daxil, dəyişiklik; biri unudulsa səhv
yenə səssiz olardı.

### Qərar 2 — açıq (explicit) birləşmə, relaksasiya + addım həddi ilə

BHP əvvəlki addımın debitlərindən tərs traverse ilə tapılır. Rəqsə qarşı
ω = 0.5 relaksasiya və addım başına 25 bar hədd. Ölçüldü (Seans 23): üç
qaçışda rəqs yoxdur, THP hədəfə 0.2 bar dəqiqliklə oturur, maksimal sapma
yalnız ilk addımlardadır.

**Alternativ təxirə salındı ⏳:** tam implicit birləşmə (THP tənliyi
Nyuton sisteminə) — qalıq/Jakobianı dəyişir, B3-B-dən sonra qazanılmış
yığılma sabitliyini riskə atardı.

### Qərar 3 — dəstəklənməyən yerlərdə AÇIQ imtina

IMPES və Eclipse ixracı THP-ni tanımır. Səssiz yaxınlaşma (IMPES-də
BHP = THP, ixracda LRAT) əvəzinə aydın xəta verilir.


## Q-16 — Qaz-neft kapilyar təzyiqi Brooks-Corey formasındadır; Pcow düzəlişi ayrı addımdır

**Tarix:** 14 sentyabr 2026 · **Kontekst:** B5-b
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 24)

### Qərar 1 — Pcog su-neft modelinin GÜZGÜSÜDÜR

    Pcog(Sg) = Pe · S_L,n^(−1/λ),   S_L,n = (1 − Sg − Swc − Sorg) / (1 − Swc − Sorg)

Səbəb: `BrooksCoreyCapillaryProvider` (A4) eyni düsturu su-neft üçün işlədir.
Qaz-neft üçün başqa korrelyasiya (məs. Leverett J) seçsəydik, iki əyri fərqli
ailədən olardı və istifadəçi eyni üç parametrlə (Pe, λ, yuxarı hədd) hər ikisini
idarə edə bilməzdi.

Sg = 0-da Pcog = Pe (sabit) — qazsız hüceyrələr arasında süni axın yaranmır,
çünki potensiala hər iki tərəfdə eyni dəyər əlavə olunur.

**Alternativ rədd edildi:** "Sg = 0-da Pcog = 0 olsun" (Pe çıxılmaqla). Bu,
giriş təzyiqi anlayışını pozardı və Brooks-Corey-nin öz formasından kənara
çıxardı.

### Qərar 2 — Pcow-un qoşulması AYRI addım kimi aparıldı

B5-b-yə başlayanda məlum oldu ki, üç fazalı yolda Pcow da işləmir
(`pc = None` sabit). Ölçüldü: Pc = 1 bar ilə və onsuz qaçış bit-bit eyni idi.

Əvvəl Pcow qoşuldu və sonlu fərqlə yoxlandı, sonra Pcog əlavə edildi. Səbəb:
ikisi bir yerdə ediləydi, Jakobian xətası çıxanda hansı həddin səhv olduğu
bilinməzdi.

### Qərar 3 — Pcog modeldən qurulur, mühərrik imzasına əlavə edilmir

`SimulationService.create_engine()` bütün mühərrikləri EYNİ açar sözlərlə
qurur. Yalnız üç fazalı mühərrikin oxuduğu bir provider üçün ortaq imzanı
genişlətmək iki fazalı mühərriyə mənasız parametr əlavə etmək olardı. Ona görə
Pcog provider-i `model.gas_capillary_parameters`-dan mühərrikin öz içində
qurulur (B4-B-də lülə həndəsəsi ilə eyni yanaşma).

### Qərar 4 — qaz fazası olmayan modeldə Pcog XƏBƏRDARLIQ verir

Səssizcə atmaq məhz bu seansda düzəldilən səhvin təkrarı olardı: istifadəçi
parametr verir, nəticədə heç nə dəyişmir və səbəbi heç yerdə yazılmır.


## Q-17 — THP iş nöqtəsi NODAL ANALİZLƏ tapılır (keçən addımın debiti ilə yox)

**Tarix:** 15 sentyabr 2026 · **Kontekst:** B4-B sabitliyi
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 25)

### Qərar 1 — iş nöqtəsi IPR ∩ VLP kəsişməsidir

B4-B-də BHP "keçən addımın debiti ilə tərs traverse" kimi hesablanırdı. Bu,
PRİNSİPCƏ səhvdir: tələb olunan BHP debitdən güclü asılıdır (ölçüldü — 5 m³/gün-də
71 bar, 1000-də 245 bar), debit isə BHP-dən. Nəzarətçi VLP əyrisi boyunca sıçrayır
və bistabil dövrəyə düşürdü.

İndi hər addımda `THP(BHP, q(BHP)) = THP_hədəf` tənliyi BHP üzrə həll olunur;
`q(BHP) = J·(p_lay − BHP)`, J son iş nöqtəsindən yenilənir. Funksiya monotondur
(BHP ↑ → debit ↓ → sürtünmə ↓ → THP ↑), ona görə kəsişmə yeganədir.

**Alternativ rədd edildi:** tam implicit THP birləşməsi (THP tənliyi Nyuton
sisteminə). Ölçmə göstərdi ki, problem gecikmə DEYİL — sabit BHP ilə qaçış da
qeyri-stabil idi, nodal düzəlişdən sonra isə açıq birləşmə kifayət etdi:
Δt medianı 38 dəfə yaxşılaşdı. Qalıq/Jakobianı riskə atmağa əsas qalmadı.

### Qərar 2 — quyu YALNIZ həqiqi axan tərkibə görə bağlanır

Bağlanma qərarı durğun (axınsız) sütuna əsaslana bilməz: o, ən ağır sütundur və
işləyən quyunu da "axa bilmir" kimi göstərir (ölçüldü — PROD-1 t = 0-da səhvən
bağlandı). Başlanğıcda BHP lay təzyiqinin `STARTUP_DRAWDOWN_BAR` = 5 bar altına
qoyulur ki, quyu işə düşsün və nəzarətçi həqiqi tərkibi öyrənsin.

### Qərar 3 — bağlanma = quyu indeksinin sıfırlanması

Alternativ (BHP-ni lay təzyiqində saxlamaq) sınandı və ÖLÇÜLDÜ: sıfıra yaxın
drawdown ədədi cəhətdən ən pis haldır, Δt 0.003 günə düşdü. Quyu indeksini
sıfırlamaq isə quyunu tənlikdən tamamilə çıxarır — qalıq və Jakobian kodu
dəyişmir, onlar sadəcə WI = 0 görür. Yenidən açılma 2 barlıq ehtiyatla
(histerezis) olur.


## Q-18 — Qaz vurulması: faza bağlantıya çıxarılır, qazsız modeldə isə XƏTADIR

**Tarix:** 15 sentyabr 2026 · **Kontekst:** B7 addım 1
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 26)

### Qərar 1 — vurulan faza `WellConnection`-a çıxarılır

Qalıq və Jakobian quyunun özünü deyil, yalnız bağlantını görür (B4-B-də THP
üçün seçilmiş eyni yanaşma). Ona görə `injected_phase` bağlantıya köçürüldü —
`well_rates` və `ThreePhaseWellJacobian` modelə müraciət etmir.

### Qərar 2 — vurulan fazanın mobilliyi SON NÖQTƏ mobilliyidir

`krg_end / μ_g` (qaz) və `krw_end / μ_w` (su). Səbəb: vurulan faza quyu dibini
öz doyma həddində doldurur, ona görə qarışığın nisbi keçiriciliyi deyil, həmin
fazanın son nöqtəsi işlədilir. Su vurulmasında bu, A6-dan bəri belədir —
qaz üçün eyni konvensiya saxlanıldı ki, iki yol uyğunsuz olmasın.

### Qərar 3 — qazsız modeldə qaz vurulması XƏTADIR, xəbərdarlıq deyil

Qaz fazası söndürüləndə iki fazalı mühərrik seçilir və o, vurulan fazanı
ümumiyyətlə oxumur — quyu səssizcə SU vurardı. Nəticə "parametr işləmədi"
deyil, "tamamilə başqa flüid vuruldu" olardı, ona görə model BLOKLANIR.


## Q-19 — RATE hədəfi perforasiyalara WI·λ nisbətində bölünür; λ əvvəlki addımdandır

**Tarix:** 15 sentyabr 2026 · **Kontekst:** B7 addım 2-yə hazırlıq
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 27)

### Qərar 1 — pay Peaceman nisbətidir

    pay_c = WI_c · λ_c / Σ_k WI_k · λ_k

Səbəb: BHP rejimində debit perforasiyalar arasında məhz belə paylanır (hamısı
eyni BHP-ni görür). Beləliklə RATE ↔ BHP keçidində (B7 addım 2) təbəqələrin
payı sıçramır.

**Alternativlər rədd edildi:** bərabər pay və yalnız WI ilə pay — sulanmış
(neft mobilliyi aşağı) təbəqəyə yenə eyni debit yazardı.

### Qərar 2 — λ addımın ƏVVƏLİNDƏKİ yığılmış vəziyyətdəndir

Cari iterasiyanın λ-sı ilə bir perforasiyanın debiti quyunun bütün digər
hüceyrələrinin doymuşluğundan asılı olardı — Jakobianda hüceyrələr arası yeni
törəmələr və yeni seyrəklik strukturu yaranardı. Pay Nyuton daxilində sabit
olanda RATE törəmələri sadəcə paya vurulur, struktur dəyişmir və Jakobian
dəqiq qalır (iki fazalıda ölçüldü: < 10⁻⁸).

Qiyməti: təbəqələr arası paylanma bir addım gecikir. Quyunun CƏMİ debiti isə
dəqiqdir.

**Alternativ təxirə salındı ⏳:** quyu BHP-sini Nyuton naməlumu etmək
(`StandardWellModel` / `CoupledNewtonSolver` kodu var, mühərrikə qoşulmayıb) —
paylanmanı təbii verir, lakin qalıq/Jakobian strukturunu dəyişir.

### Qərar 3 — tək perforasiyada pay dəqiq 1.0, yeniləmə yalnız lazım olanda

`needs_rate_allocation` çox perforasiyalı RATE quyusu yoxdursa mühərrikdə
yeniləməni tamamilə söndürür — mövcud modellər bit-bit eyni və əlavə xərcsizdir.

### Qərar 4 — çox perforasiyalı RATE modellərinin nəticəsi dəyişir

Sahibkarın seçimi (Seans 27): köhnə davranışı saxlamaq və ya bloklamaq yerinə
düzəltmək.


## Q-20 — Eclipse ixracında RATE quyusu RESV rejimi ilə yazılır; qaz vurucusu imtina edilir

**Tarix:** 15 sentyabr 2026 · **Kontekst:** B7 addım 2 (ixrac yoxlaması)
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 27 (davamı))

### Qərar 1 — `LRAT`/`RATE` əvəzinə `RESV`

Mühərrikdə RATE hədəfi lay həcmində maye (su+neft) debitidir. Eclipse-də bu
mənanı daşıyan rejim `RESV`-dir; iki fazalı deck-də RESV məhz su+neftin lay
həcmidir. `LRAT` səth debitidir — onu yazmaq deck-i modeldən `Bo` qədər fərqli
edərdi.

**Alternativ rədd edildi:** `LRAT`-ı düzgün (7-ci) sütunda saxlamaq — sütun
səhvi düzələrdi, lakin mənaca uyğunsuzluq qalardı.

### Qərar 2 — qaz vurucusu olan modeldə ixrac XƏTA verir

Deck iki fazalıdır. Qaz vurucusunu `WATER` kimi yazmaq Q-18-də bağlanan
səhvin ixracdakı təkrarı olardı; `GAS` yazmaq isə qaz fazası olmayan deck-də
etibarsızdır. THP-də (Q-15, Qərar 3) olduğu kimi açıq imtina seçildi.


## Q-21 — BHP limiti: rejim keçidi addımlar arasındadır, meyar Peaceman-ın tərsidir

**Tarix:** 15 sentyabr 2026 · **Kontekst:** B7 addım 2
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 27 (davamı): BHP limiti)

### Qərar 1 — nəzarətçi bağlantının `mode`/`target`-ini dəyişir, qalıq/Jakobian toxunulmur

`BhpLimitController` yığılmış addımdan sonra işləyir; limit pozulubsa bağlantı
`mode = BHP`, `target = limit` olur və addım `resolve_step` ilə yenidən həll
olunur. Qalıq və Jakobian yalnız mövcud BHP/RATE budaqlarını görür.

**Alternativlər rədd edildi:**
* yeni `ControlMode` dəyəri — Q-15-də sənədlənmiş TƏHLÜKƏ: bütün `mode is BHP`
  yoxlamaları onu səssizcə RATE kimi işlədərdi;
* Nyuton daxilində keçid (quyu tənliyi, `StandardWellModel`) — iterasiya
  daxilində rejim dəyişməsi qalıqda sıçrayış yaradır, Seans 25-dəki rəqs riski.

### Qərar 2 — keçid meyarı: hədəf debiti verən BHP, qalığın ÖZ mobilliyi ilə

`BHP = (Σ WI·λ·p ∓ q) / Σ WI·λ`, λ — `well_rates`-in BHP budağındakı eyni
ifadə (`connection_mobilities`). Nəticədə keçid anında cəmi debit sıçramır
(ölçüldü: < 10⁻⁹). Mobillik sıfırdırsa (quyu axa bilmir) quyu da limitə keçir —
RATE budağı belə halda hərəkətsiz hüceyrədən debit "yaradardı".

### Qərar 3 — histerezis 2 bar + addımda ən çox bir keçid

RATE-ə qayıdış yalnız hədəf limitdən `RATE_RESTORE_MARGIN_BAR` = 2 bar əlverişli
BHP ilə əldə olunanda (THP-dəki `REOPEN_MARGIN_BAR` ilə eyni qiymət). Bir
addımın təkrar həlləri arasında quyu ikinci dəfə keçə bilməz — dövr ən çox
`quyu sayı + 1` dəfə fırlanır və sonuncu qiymətləndirmə qəbul olunan həll
üzərindədir.

### Qərar 4 — təkrar həll yığılmasa əvvəlki həll saxlanılır

Rejim isə növbəti addım üçün qüvvədə qalır (limit növbəti addımda artıq
pozulmur). Rejimi geri qaytarmaq eyni addımın hər dəfə təkrar uğursuz olmasına
aparardı.

### Qərar 5 — limitli RATE quyusu təzyiq idarəsi sayılır; IMPES-də xəta, BHP/THP rejimində xəbərdarlıq

Bax ISH_HESABATI §5.


## Q-22 — Sahə qaz seriyaları quyular üzrə yığılır; qaz debiti ikinci oxda çəkilir

**Tarix:** 15 sentyabr 2026 · **Kontekst:** qaz vurulmasının hesabatı
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 28)

### Qərar 1 — hasilat və vurma quyu tipinə görə ayrılır, hüceyrə işarəsinə görə yox

Hüceyrə massivinin cəmi (`rates.gas.sum()`) vurucu və istismarçını qarışdırır.
Quyu lüğətləri (`per_well_gas`) quyu tipini açıq bilir, ona görə sahə qaz
debiti = istismarçıların cəmi, vurulan qaz = vurucuların cəmi.

**Alternativ rədd edildi:** hüceyrə işarəsinə görə maska (`rates.gas < 0`) —
su üçün belə edilir, lakin eyni hüceyrədə həm vurucu, həm istismarçı
perforasiyası olsa yenə qarışardı.

### Qərar 2 — qaz debitləri «Debitlər» panelində ikinci oxdadır

Qaz sm³/gün maye m³/gündən 2–3 tərtib böyükdür. Ayrıca panel 2×2 düzümü
pozardı (mövcud çağırışlar və testlər 2×2 gözləyir), eyni ox isə maye
xətlərini oxunmaz edərdi.


## Q-23 — SPE1 hədəfi CASE2-dir; etalon OPM Flow nəticəsidir, öz oxuyucumuzla oxunur

**Tarix:** 15 sentyabr 2026 · **Kontekst:** B7 addım 3
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 29, [SPE1.md](SPE1.md))

### Qərar 1 — əvvəl SPE1CASE2 (dəyişən doyma təzyiqi)

İki OPM deck-i yalnız `DRSDT 0` ilə fərqlənir. Mühərrik vurulan qazı həmişə
doymamış neftdə həll edir, yəni CASE2-nin fizikasını daşıyır. CASE1 üçün yeni
fizika (Rs-in artımını qadağan edən qayda — qalıq və Jakobiana toxunur) lazım
olardı.

**Alternativ təxirə salındı ⏳:** CASE1 — CASE2 müqayisəsi bitəndən sonra.

### Qərar 2 — etalon: OPM Flow-un `SPE1CASE2` summary faylı

Odeh (1981) məqaləsinin cədvəl rəqəmləri əlimizdə deyil; OPM Flow-un nəticəsi
açıq, rəqəmsaldır və məhz eyni deck-dən alınıb.

### Qərar 3 — `resdata` asılılığı əlavə edilmir

Summary formatı sadədir (big-endian Fortran qeydləri); oxumaq üçün ~90 sətir
kifayət etdi. Ağır native kitabxana (bax `OPM_IMPORT.md` — natamam faylda
SIGABRT) yalnız iki massiv üçün əsaslandırılmır.

### Qərar 4 — iş sırası: səth debiti → süxur istinadı → SGOF → PVT → model

Ən kiçik və ən az riskli boşluqlar əvvəl; fizikaya (qalıq/Jakobian) toxunan
G2 sona yaxın, hər biri ayrıca sonlu fərq testi ilə — Seans 26-dakı "iki yeni
xüsusiyyəti eyni anda qalığa salma" prinsipi.


## Q-24 — Səth debiti hədəfi ayrıca parametrdir və addımlar arasında lay həcminə çevrilir

**Tarix:** 15 sentyabr 2026 · **Kontekst:** B7 addım 3, SPE1 boşluğu G5
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 30)

### Qərar 1 — `ControlMode` deyil, `WellControl.rate_basis`

Yeni `ControlMode` dəyəri Q-15-də sənədlənmiş təhlükəni yaradardı (bütün
`mode is BHP` yoxlamaları onu səssizcə RATE kimi işlədərdi). Baza RATE-in
mənasını dəqiqləşdirir, rejimin özünü dəyişmir; defolt `RESERVOIR` — mövcud
modellər bit-bit eynidir.

### Qərar 2 — proqnoz + düzəliş, qalıq/Jakobian toxunulmur

Addımdan əvvəl yığılmış vəziyyətin `(1−f)/Bo`, `1/B` əmsalları ilə səth hədəfi
lay həcminə çevrilir; addımdan sonra əldə olunan səth debitinə görə miqyaslanıb
addım təkrarlanır. Qiymət ölçüldü (Seans 30): iki fazalıda 6, üç fazalı qaz
vurmada 21 əlavə həll (300 gün); dəqiqlik 10⁻³ daxilində.

**Alternativ rədd edildi:** neftin səth debitini birbaşa qalığa yazmaq
(`q_neft = −hədəf`, su `q_neft·(Bo/Bw)·λw/λo`). Dəqiqdir, lakin iki mühərrikin
RATE Jakobianını dəyişir, `λo → 0`-da təkil olur və üç fazalı RATE Jakobianı
onsuz da təqribidir (TB-1).

### Qərar 3 — BHP limiti cari lay hədəfini callback ilə alır

Səth bazalı quyuda lay hədəfi addım-addım dəyişir; limitdən RATE-ə qayıdan quyu
qurulma anındakı (səth) rəqəmi bərpa etsəydi hədəf ~Bo qədər səhv olardı.
`BhpLimitController(..., rate_target=SurfaceRateController.rate_target)`.

### Qərar 4 — IMPES-də xəta, RATE olmayan rejimdə xəbərdarlıq

IMPES-də səth rəqəmi səssizcə lay həcmi kimi işlənərdi. BHP/THP rejimində
baza sadəcə işləmir — Q-21, Qərar 5 ilə eyni prinsip.


## Q-25 — Süxur sıxılmasının istinad təzyiqi ayrıca sahədir; flüidinki toxunulmur

**Tarix:** 16 sentyabr 2026 · **Kontekst:** SPE1 boşluğu G6
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 32, [SPE1.md](SPE1.md) → §4)

### Qərar 1 — istəyə bağlı sahə, `None` → datum

`RockProperties.compressibility_reference_pressure` əlavə olundu. `None`
(defolt) olduqda istinad əvvəlki kimi datum təzyiqidir, yəni mövcud bütün
modellərdə nəticə **bit-bit eynidir** (testlə kilidlənib).

**Alternativ rədd edildi:** istinadı həmişə `ROCK`-dakı kimi ayrıca tələb etmək —
bu, mövcud modellərin nəticəsini səssizcə dəyişərdi.

### Qərar 2 — flüidin statik sıxılma istinadı DƏYİŞMİR

Kodda iki ayrı istinad var: məsamə həcmi üçün (G6 bunu ayırır) və PVT provider
olmayanda flüidin `B(p)` modeli üçün (`ResidualAssembler.reference_pressure`,
oradan `DerivativeProvider`-ə). İkincisi datum olaraq qaldı: süxur üçün verilən
rəqəmin flüidin sıxılmasını da dəyişməsi istifadəçinin gözlədiyi davranış deyil.

### Qərar 3 — Jakobian toxunulmur (riyazi əsas)

`d(PV)/dp = PV_ref·c_r` istinaddan asılı deyil. Kod dəyişmədi, iddia isə
istinad datumdan fərqli olan modeldə sonlu fərq testi ilə kilidləndi.

### Qərar 4 — IMPES genişləndirilmir

IMPES məsamə həcmini təzyiqlə miqyaslamır (ümumi sıxılma `ct` yanaşması) —
bu, G6-dan əvvəlki sadələşdirmədir. Yeni sahəni orada "bir az" tətbiq etmək
iki fərqli modeli qarışdırardı. Sahənin sənədində açıq yazılıb; istifadəçi
xəbərdarlığı ⏳ ayrıca iş kimi qeyd olundu.


## Q-26 — Qaz-neft cədvəli AYRI sinifdir; üç fazalı yol artıq SWOF cədvəlini atmır

**Tarix:** 16 sentyabr 2026 · **Kontekst:** SPE1 boşluğu G4
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 33, [SPE1.md](SPE1.md) → §4)

### Qərar 1 — `GasSaturationTable` ayrı sinifdir

Su-neft cədvəlinə sütun əlavə etmək əvəzinə ayrı sinif yazıldı: arqument `Sg`-dir,
`krg` artır, `krog` azalır — yoxlama qaydaları da tərsdir. Eyni sinifdə saxlamaq
sütun adlarını və monotonluq yoxlamasını yanıldıcı edərdi.

### Qərar 2 — modeldə AYRI sahə (`gas_scal_tables`)

Cədvəl `gas_scal_parameters` sahəsinə "sıxışdırılmadı" (o, Corey parametrləri
üçündür). Su-neft tərəfdə artıq `scal_parameters` + `scal_tables` cütü var —
qaz tərəf eyni konvensiyanı təkrarlayır, yəni oxuyan adam üçün sürpriz yoxdur.

### Qərar 3 — Stone ilə əlaqə duck-typing-lə qalır

`StoneRelativePermeabilityProvider` dəyişdirilmədi: o, qaz obyektindən yalnız bir
neçə metod istəyir və cədvəl həmin müqaviləni ödəyir. Beləliklə üç fazalı
qalıq/Jakobian G4-dən XƏBƏRSİZDİR.

### Qərar 4 — üç fazalı yol su-neft CƏDVƏLİNİ işlədir (səssiz səhvin düzəlişi)

`_create_three_phase_engine` su-neft provider-ini həmişə Corey-dən qururdu və
modeldəki SWOF cədvəli səssizcə atılırdı. İndi hər iki mühərrik eyni seçim
qaydasından (`_relative_permeability`) keçir. SPE1 həm SWOF, həm SGOF tələb edir,
ona görə bu, seçim deyil, zərurət idi.

### Qərar 5 — regionlar və layihə faylı

Qaz cədvəlində hazırda yalnız defolt region işlədilir (Stone müqaviləsində region
arqumenti yoxdur) ⏳. SCAL cədvəlləri (su-neft də daxil) layihə faylında
saxlanmır — bu, G4-dən əvvəlki vəziyyətdir və qaz üçün fərqli davranış icad
edilmədi; backlog-a yazıldı.


## Q-27 — Deck PVT-si İTKİSİZ oxunur, birləşdirmə isə AÇIQ təqribdir

**Tarix:** 16 sentyabr 2026 · **Kontekst:** SPE1 boşluğu G1
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 34, [SPE1.md](SPE1.md) → §4)

### Qərar 1 — oxuma `io` qatında, `scal_io` ilə eyni üslubda

Deck mətninin ayrışdırılması GİRİŞ sərhədidir: `io/pvt_io.py` açar söz axtarışı,
`--` şərhləri və `/` ayırıcısı üçün `scal_io.read_swof`-un EYNİ qaydalarını
işlədir. Domain qatı deck formatından xəbərsiz qalır.

### Qərar 2 — `DeckPvt` itkisizdir, `to_pvt_table` isə təqribdir

`PVTTable` bir doymamış qol daşıya bilir; deck-də bir neçə Rs qolu olur.
Ona görə oxuma bütün qolları saxlayır, birləşdirmə isə yalnız SEÇİLƏN qolu
cədvələ köçürür və neçə qolun köçürülmədiyini jurnala yazır. Beləliklə itki
SƏSSİZ deyil, ölçülə bilən və sənədlənmişdir.

**Alternativ rədd edildi:** `PVTTable`-ı çox qollu etmək — bu, domain
strukturunu və bütün mühərrik yollarını dəyişərdi; həmin iş G2/G3-ün öz
mövzusudur və ayrıca aparılmalıdır.

### Qərar 3 — Bg üçün YENİ vahid növü

Mühərrikdə Bg ölçüsüzdür (lay m³ / səth m³), deck-də isə `rb/Mscf`. Bu, gözdən
qaçan uyğunsuzluqdur, ona görə `gas_fvf` növü açıq şəkildə əlavə olundu.
Əmsallar mövcud həcm sabitlərindən TÖRƏDİLİR (`rb`, `ft3`) — əl ilə rəqəm
yazılmır ki, iki yerdə fərqli dəyər qalmasın.

### Qərar 4 — səssiz ehtiyat qiymət xəbərdarlığa çevrildi

Bir doymamış sətir olanda provider `c_o`-nu fit edə bilmir və ehtiyat qiymətə
düşür (ölçüldü: 15 dəfə böyük). Birləşdirmə indi bundan açıq xəbərdarlıq verir.
Düzgün həll — `c_o`-nu deck qolundan hesablamaq — G3-ə aiddir.


## Q-28 — Doymamış özlülük üstəli CƏDVƏLDƏN fit olunur; G2 iki commit-ə bölünür

**Tarix:** 16 sentyabr 2026 · **Kontekst:** SPE1 boşluğu G2
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 35)

### Qərar 1 — üstəl korrelyasiyadan GÖTÜRÜLMÜR

`μo = μo_b·(p/Pb)^n`-də `n` korrelyasiyada 0.278-dir, SPE1 deck-ində isə
0.4602 və 0.5085 (ölçüldü, Seans 34). Ona görə `n` cədvəlin öz doymamış
sətirlərindən fit olunur. Sabit qiymət işlətmək real deck-də özlülük artımını
təxminən iki dəfə az göstərərdi.

### Qərar 2 — fit mümkün olmayanda XƏBƏRDARLIQ

Pb-dən yuxarı iki sətir olmayanda korrelyasiya qiymətinə düşülür, lakin bu,
jurnalda açıq yazılır (Q-27-dəki eyni qayda — səssiz yanlış dəyər ən pis haldır).

### Qərar 3 — G2 İKİ commit-dir

(a) provider + testlər (mühərrik toxunulmur, tam dəst dəyişməz qalmalıdır);
(b) mühərriyə qoşulma + sonlu fərq testləri.

Səbəb Seans 26-dakı ilə eynidir: iki yeni şey eyni anda qalığa/Jakobiana
girsəydi, sonlu fərq xətası çıxanda mənbəyi ayırmaq çətinləşərdi.

### Qərar 4 — cədvəl diapazonundan kənarda törəmə SIFIRDIR

`dPb/dRs` cədvəlin Rs diapazonundan kənarda sıfırdır (Bo qolunda da belədir),
ona görə ∂μo/∂Rs də sıfır olur. Bu, EKSTRAPOLYASİYA ETMƏMƏK qərarıdır və sonlu
fərq də eyni nəticəni verir, yəni Jakobian cədvəlin öz davranışı ilə uyğundur.


## Q-29 — Özlülüyün Rs-asılılığı Jakobiana daxil edilir; doyma nöqtəsi lövbəri fit-dən alınır

**Tarix:** 16 sentyabr 2026 · **Kontekst:** SPE1 boşluğu G2b
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 36)

### Qərar 1 — μo flüid vəziyyəti vasitəsilə daşınır

`mu_o_p`/`mu_o_rs` sahələri əlavə olundu — Bo-dakı (B3-B) eyni naxış. Qalıq və
Jakobian provider-i tanımır, yalnız massivləri görür.

### Qərar 2 — mobilliyin Rs törəməsi İKİ hədddir

`∂mob_o/∂Rs = −mob_o·(Bo'_Rs/Bo + μo'_Rs/μo)`. Ölçüldü: ikinci hədd olmadan
3-cü sütun 3.6×10⁻¹¹-dən 4.2×10⁻¹-ə düşür.

### Qərar 3 — RATE quyusunda SU tənliyinə yeni element

`f` özlülükdən asılı olduğu üçün doymamış halda `blocks[c, 0, 2] ≠ 0`.

### Qərar 4 — LÖVBƏR cədvəl interpolyasiyasından DEYİL, fit-dən alınır

Pb adətən düyün deyil və əyrinin orada sınığı var; interpolyasiya μ lövbərini
2.27 %, Bo lövbərini 0.37 % şişirdirdi. Düzgün lövbər doymamış qolun öz
fit-inin Pb-dəki qiymətidir (fit onsuz da aparılır).

**Pb-dən yuxarı şəbəkə DOYMUŞ meylin davamıdır.** İki variant ölçülüb rədd
edildi: `np.where` keçidi (lövbərdə sıçrayış — sonlu fərq ∂/∂Rs 27-yə qalxırdı)
və doymamış qolun qiymətləri (meyl `−c_o` ilə kompensasiya olunub analitik
törəməni sıfıra endirirdi). Meyar: **törəmə uyğunluğu qiymət öz-özünə
uyğunluğundan vacibdir** — Nyuton törəmə ilə işləyir.

### Qərar 5 — invariant testlər yenilənir (sahibkarın qərarı ilə)

Köhnə `atol=1e-5` şərti yalnız ona görə keçirdi ki, lövbər və cədvəl sütunu
EYNİ interpolyasiya səhvini bölüşürdü. Yeni şərt: doymuş zonanın düyünlərində
tam dəqiqlik, qalan yerlərdə ölçülmüş hədd (0.0195 → 0.025).

### Qərar 6 — etalon hasilat əmsalı 62.86-da QALIR

Müvəqqəti 62.73 dəyişikliyi geri alındı: ölçmə göstərdi ki, fərq fizikadan yox,
lövbər qüsurundan gəlirdi. Lövbər düzəldiləndən sonra iki və üç fazalı
mühərriklərin qazsız rejimdəki uyğunluğu da bərpa olundu (62.8612).


## Q-30 — Doymamış neftin c_o və n-i HƏR PVTO qolundan; Rs üzrə interpolyasiya

**Tarix:** 16 sentyabr 2026 · **Kontekst:** SPE1 boşluğu G3
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 37)

### Qərar 1 — tək c_o/n YOX (sahibkarın qərarı)

Ölçüldü: SPE1CASE2-nin iki qolunda c_o 2.056e-4 ↔ 1.832e-4 1/bar (~11 %),
n 0.460 ↔ 0.580. Əvvəlki "0.3 %" rəqəmi sınaq deck-indəki səhv köçürülmüş
sətirdən gəlirdi. Üç variant sahibkara verildi; **çox qollu provider** seçildi.

### Qərar 2 — qollar provider-ə AYRICA verilir, `PVTTable` dəyişmir

`BlackOilPVTProvider(table, oil_branches=...)`. Cədvəl doymuş əyrini (bütün
qolların başları) daşıyır, qollar isə doymamış parametrləri. `PVTTable`-ın
quruluşuna (serializasiya, UI, IMPES) toxunulmadı. `simulation` qatı `io`-nu
import etmir — qollar duck-typing ilə qəbul olunur.

### Qərar 3 — qollar arasında Rs üzrə parçalı xətti interpolyasiya, kənarda sabit

Qolların özündə deck DƏQİQ təkrarlanır. Qollar arasında törəmə interval
meylidir (Jakobian qalıqla dəqiq uyğun), kənarda sıfırdır — ekstrapolyasiya
yoxdur (Q-28 Qərar 4). ⏳ OPM-in öz qaydası mənbədən yoxlanılmayıb; bu,
bizim fərziyyəmizdir və açıq yazılıb.

### Qərar 4 — lövbər deck-in doymuş başlarından

Deck-də Pb düyündür, ona görə lövbər fit-dən yox, birbaşa başlardan götürülür
(Q-29-un qolsuz yolu dəyişmir).

### Qərar 5 — uyğunsuzluq AÇIQ xətadır

Qollar cədvəlin doymuş əyrisinə uyğun gəlmirsə, doymamış sətir yoxdursa və ya
Bo təzyiqlə artırsa — `ValueError`. Cədvəl ən böyük qoldan əvvəl bitirsə
(Rs kəsiləcək) — xəbərdarlıq.


## Q-31 — Canlı neftin sıxlığı həll olmuş qazı daxil edir; SPE1 qurucusu ayrıca paketdədir

**Tarix:** 16 sentyabr 2026 · **Kontekst:** SPE1CASE2 hazırlığı
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 38)

### Qərar 1 — `ρo = (ρo_səth + Rs·ρg_səth)/Bo` (sahibkarın qərarı: modeldən ƏVVƏL)

Ölçüldü: SPE1-də +27 %. Üç fazalı axının cazibə həddi və ilkin tarazlıq
dəyişdi; iki fazalı yol və IMPES-də Rs = 0 olduğu üçün nəticə eynidir.

### Qərar 2 — cazibə həddində sıxlıq törəmələri Jakobiana daxildir (üç fazalı)

A6-dan qalan "sıxlığın təzyiqdən asılılığı nəzərə alınmır" sadələşdirməsi üç
fazalı yolda aradan qaldırıldı — Rs sıxlığa girəndən sonra 3-cü sütun onsuz
7.7×10⁻³ səhv olurdu. Upstream seçimi hələ də diferensiallaşdırılmır.

### Qərar 3 — ilkin tarazlıq mühərriklə EYNİ Rs və Bo qaydasını işlədir

Əks halda ilkin vəziyyət mühərrikin cazibə həddi ilə taraz olmur (ölçüldü:
5.7 m³/gün süni şaquli axın).

### Qərar 4 — `OilBranch` domain-dədir

Model və layihə faylı qolları daşıyır; `application` qatı `io`-nu import
etmir. `io/pvt_io.py` sinfi oradan import edir (köhnə importlar işləyir).

### Qərar 5 — etalon modellər `imex2d/benchmarks/`-dadır

Qurucu həm deck oxuyucularını (`io`), həm model qurucusunu (`application`)
işlədir — `app.py` kimi kompozisiya qatıdır; `domain`/`simulation` onu tanımır.

### Qərar 6 — qazın səth debiti üçün xəbərdarlıq həddi 1e7 sm³/gün

Evristikadır, fiziki sərhəd deyil; maye həddi (1e5 m³/gün) dəyişmədi.

### Qərar 7 — SPE1 nəticəsi etalon kimi YAZILMIR

İlk müqayisədə izah olunmamış fərqlər var (qaz 375 gün tez çatır). Golden
fayl yalnız fərqlərin səbəbi ölçüləndən sonra yazılacaq.


## Q-32 — Doymuş PVTO qolu cədvəldən yuxarı XƏTTİ uzadılır (yalnız deck yolunda)

**Tarix:** 17 sentyabr 2026 · **Kontekst:** SPE1 boşluğu G7
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 40 və 41)

### Problem

Deck-in PVTO cədvəli 5014.7 psia-da bitir, SPE1CASE2-də isə vurucunun
ətrafında təzyiq 7500 psia-dan yuxarı qalxır. Biz `np.interp` ilə Rs_sat-ı
platoda saxlayırdıq. Ölçüldü: eyni təzyiqdə bizdə sərbəst qaz OPM-dən çox
olur, cəbhə tez gedir və bütün sonrakı fərqlər bundan doğur.

### OPM-in qaydası — MƏNBƏDƏN oxundu

* `opm-common/opm/material/fluidsystems/blackoilpvt/LiveOilPvt.hpp:512` —
  `saturatedGasDissolutionFactor` cədvəli **`extrapolate=true`** ilə
  çağırılır;
* `opm-common/opm/material/common/Tabulated1DFunction.hpp:266-283` — bu
  bayraq qoyulanda son seqmentin xətti düsturu cədvəldən KƏNARDA da tətbiq
  olunur ("extended beyond its range by straight lines").

### Qərar 1 — uzantı YALNIZ deck yolunda

Korrelyasiya ilə qurulan cədvəldə Pb-dən yuxarı Rs platosu HƏQİQİ
fizikadır (neftin tərkibi sabitdir) — orada uzatmaq YANLIŞ olardı.
Deck-də isə plato sadəcə məlumatın bitməsidir. Ona görə uzantı
`oil_branches` verildikdə qurulur; korrelyasiya modelləri bit-bit eyni
qalır (testlə kilidlənib).

### Qərar 2 — meyl SON İKİ doymuş düyündən

OPM eyni qaydadır. Bo və μo da eyni meyllə uzanır, çünki onlar Rs ilə
uzlaşan doymuş vəziyyəti təsvir edir.

### Qərar 3 — törəmələr də uzanır

Uzantı qoyulub törəmə sıfır qalsaydı, Nyuton yanlış Jakobianla işləyərdi —
səssiz səhvin ən pis sinfi. `_slope` və `_saturation_pressure_slope`
uzantının öz meylini qaytarır; sonlu fərqlə yoxlanılıb.

### Qərar 4 — müsbətlik qoruyucusu

Xətti uzantı kifayət qədər yüksək təzyiqdə μo_sat-ı sıfırdan keçirər.
Hədd + xəbərdarlıq qoyuldu (SPE1-də işə düşmür). Səbəb: səssiz mənfi
özlülük fəlakət olardı.

### Ölçülmüş nəticə

İstismarçının BHP-si 1034-cü gündə −60.7 % → **−1.2 %**; FGOR +389 % →
**+38.7 %**; BHP limitinə keçid 1120 → **1381 gün** (OPM 1550). Addım
sayı 514 → 375.

### Açıq qalan

Vurucu bağlantısının mobilliyi (1-ci gündə −34.8 %) bu qərara DAXİL
DEYİL — ayrı məsələdir və OPM qaydası hələ mənbədən oxunmayıb ⏳.


## Q-33 — Vurucu bağlantısı hüceyrənin TAM mobilliyini işlədir (hər iki mühərrik)

**Tarix:** 17 sentyabr 2026 · **Kontekst:** SPE1CASE2-də vurucunun 1-ci
gündəki BHP-si −34.8 % (bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 42)

### Qayda MƏNBƏDƏN

`StandardWell_impl.hpp:264-315` — `total_mob = Σ mob[faza]`,
`cqt_i = − Tw · total_mob · Δp`; `WellInterface_impl.hpp:2261-2278` —
mobilliklər birbaşa hüceyrədən götürülür, vurucu üçün xüsusi hal yoxdur.

### Qərar 1 — üç fazalı mühərrikdə tam mobillik

`λw + λo + λg`. Köhnə qayda (vurulan fazanın son nöqtə mobilliyi) bloku
əvvəlcədən vurulan faza ilə dolmuş sayırdı. ÖLÇÜLDÜ: səhv məhz vurmanın
BAŞLANĞICINDA böyükdür (1-ci gün −34.8 %, 304-cü gün −1.2 %) — blok
həqiqətən dolandan sonra iki qayda üst-üstə düşür.

### Qərar 2 — Jakobian genişlənir

Mobillik doyumlardan asılı olduğu üçün vurucunun sətrində Sw (və üç
fazalıda 3-cü dəyişən) sütunları YARANIR. Sonlu fərqlə yoxlanılıb
(1.4×10⁻¹⁰).

### Qərar 3 — İKİ FAZALI mühərrikdə də eyni qayda (sahibkarın qərarı)

Yalnız üç fazalı yol dəyişəndə iki mühərrikin uyğunluq testi düşdü
(62.86 ↔ 62.83). Seçim sahibkara verildi; qərar: qayda hər iki yolda eyni
olsun. ÖLÇÜLMÜŞ TƏSİR: mövcud iki fazalı modellərdə RF 62.86 → 62.83 %
(0.03 pp), addım sayı 31 → 29. Bu, su vuran bütün iki fazalı modellərə
aiddir — vurmanın başlanğıcında BHP daha yüksək çıxır.

### Qərar 4 — birləşmiş həlledici (`standard_well.py`) TOXUNULMADI

Həmin yol heç bir mühərrik tərəfindən işlədilmir (yatmış kod) və Jakobian
quruluşu fərqlidir. Uyğunsuzluq AÇIQ yazılıb ki, səssiz qalmasın; həmin
fayla növbəti dəfə toxunanda uyğunlaşdırılmalıdır ⏳.

### Ölçülmüş nəticə

WBHP INJ, 1-ci gün: 5271 psia (−34.8 %) → **8173 (+1.1 %)**, etalon 8082.
SPE1-in axın nəticələri dəyişməyib (vurucu RATE rejimindədir).

### Ölçülmüş YAN TƏSİR (gizlədilmir)

Vurma debiti artıq doyumdan asılı olduğu üçün quyu ətrafındakı sərt
ssenari sərtləşdi: `test_implicit_newton.py`-nin oscillasiya ssenarisində
CNV minimumu 0.001004 → 0.001322 (hər iki halda addım onsuz da
yığılmırdı). Geri-izləmə isə hələ də işləyir — qalığın sıçrayışı 1.4,
geri-izləməsiz 9.6. Sahibkar bunu bilərək qərar verə bilsin deyə burada
yazılıb; qayda geri qaytarılarsa, yığılma bərpa olunar.


## Q-34 — Doymamış qol deck düyünləri arasında XƏTTİ interpolyasiya olunur

**Tarix:** 17 sentyabr 2026 · **Kontekst:** Seans 40-da ölçülmüş 1.97 %
μo fərqi (bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 43)

### Problem

OPM cədvəl düyünləri arasında XƏTTİdir (`Tabulated1DFunction.hpp:282`),
bizdə isə doymamış qol üstəl qanunla gedirdi (`Bo = Bo_b·exp(c_o·ΔP)`,
`μo = μ_b·(p/Pb)^n` — Q-28/Q-30). Düyünlərdə fərq sıfır, aralıqda isə
ölçüldü: **μo 1.97 %**, Bo 0.06 %.

### Qərar 1 — deck yolunda xətti, korrelyasiya yolunda üstəl

Deck qolları (`oil_branches`) verildikdə qolun forması xəttidir; qollar
yoxdursa (korrelyasiya cədvəli) üstəl qanun QALIR — orada cədvəl onsuz da
üstəl düsturla qurulur, dəyişiklik mövcud modelləri pozardı.

### Qərar 2 — LÖVBƏR dəyişmir

Qol Pb-dəki doymuş qiymətdən başlayır (Seans 36-nın lövbəri). Alternativ
— OPM kimi tam 2 ölçülü cədvəl — Pb-də kiçik kəsilmə yaradardı və
Nyutonun doymuş/doymamış keçidini pisləşdirərdi. İki düyünlü qolda (SPE1)
iki yanaşma EYNİ nəticə verir.

### Qərar 3 — fiziki olmayan qol açıq XƏTA verir

Doymamış Bo təzyiqlə artırsa və ya μo azalırsa, `ValueError` atılır
(səssiz qəbul yox). Testlə kilidlənib.

### Ölçülmüş nəticə

Qiymətlər deck-in xətti interpolyasiyası ilə **0.000 %** üst-üstə düşür
(əvvəl μo-da 1.97 %). ∂/∂p dəqiqdir; ∂/∂Rs qollar arasında sonlu fərqlə
0.00 %, qolun düyünündə isə birtərəflidir (parçalı xətti modelin təbii
xassəsi). SPE1CASE2-yə təsiri cüzidir (FGOR > 2: 1142 → 1132 gün).


## Q-35 — Günlük göstəricilər mühərrikə toxunmadan, addımlardan qurulur

**Tarix:** 21 sentyabr 2026 · **Kontekst:** sahibkar hər günün
göstəricilərini — yataq və hər quyu üzrə, simulyasiyanın bütün günləri
üçün — görmək istədi (bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 45).

### Problem

Mühərrik adaptiv addımla gedir: nümunə qaçışında 1500 gün 88 addımda
(orta Δt 17 gün), standart açılış modelində isə 4314 addımda keçilir.
Addımlar tam günlərə düşmür, yəni «N-ci gün» sətri nəticədə yoxdur.
Vurucu quyular üzrə isə ümumiyyətlə heç bir sıra saxlanılmırdı.

### Qərar 1 — mühərrik DƏYİŞMİR, günlük cədvəl nəticədən qurulur

Alternativ — mühərriki hər gün sonunda addımı kəsməyə məcbur etmək
(hesabat anları) — nəticəni dəyişərdi: addım sayı, deməli rəqəmlər də
dəyişərdi. Gündəlik dəqiq təzyiq lazım olsa, bu imkan onsuz da var:
«Maks. Δt» = 1 gün.

### Qərar 2 — günün debiti = günün HƏCMİ / günün uzunluğu

Implicit Eyler addım boyu debiti sabit götürür. Ona görə:

- debitlər: gün bir addımın içindədirsə, o addımın debiti; addım sərhədi
  günün içinə düşürsə, zamanla çəkili orta. Günlük həcmlərin cəmi
  mühərrikin kumulyativinə bərabərdir (ölçüldü: nisbi fərq < 1e-12);
- kumulyativ və RF: addım boyu xətti, interpolyasiya dəqiqdir;
- su kəsri, GOR: günün həcmlərindən;
- orta təzyiq: yalnız addım sonlarında məlumdur, xətti interpolyasiya —
  **təxminidir**, sütun adında `(interp.)` yazılır;
- BHP/THP: günün sonunu örtən addımın dəyəri; THP-nin `nan`-ı saxlanılır.

Rədd edilən variant: «günün sonunu örtən addımın debiti». Sadədir, lakin
addım sərhədi günün içinə düşəndə günlük həcmlərin cəmi kumulyativdən
fərqlənərdi.

### Qərar 3 — vurucular üzrə sıra hər üç mühərrikdə yazılır

`well_water_injection_rate` və `well_gas_injection_rate`. Vahid sahə
sırası (`series.water_injection_rate`) ilə eynidir, vurucuların cəmi
ona bərabərdir (testlə kilidlənib). Bu, yalnız YAZMADIR — heç bir tənliyə
girmir, nəticələr bit-bit eyni qalır.
## Q-36 — Layihə faylı hesablamanı TAM təkrarlamalıdır

**Tarix:** 21 sentyabr 2026 · **Kontekst:** sahibkar sabah eyni
simulyasiyaya əl ilə doldurmadan qayıtmaq istədi
(bax [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 46).

### Problem

Model hər qaçışda panellərdən qurulur, layihə faylı isə panellərin
yalnız bir hissəsini bərpa edirdi. Fayl özü də qaz SCAL-ı, SWOF/SGOF-u və
nəticənin qaz/BHP/THP sıralarını itirirdi (qazlı modeldə RF 64.75 → 64.92 %).

### Qərar 1 — panellərin vəziyyəti faylda, sahə siyahısı ƏL İLƏ YOX

`ui/panel_state.py` panelin ictimai sadə sahələrini (spin, seçici, qutu,
mətn) atribut adı ilə oxuyur. Panelə yeni sahə əlavə olunanda o da
avtomatik saxlanılır. Bərpa sırası: seçicilər → qutular → ədədlər (vahid
çevrilməsi üçün), iki keçid (asılı diapazonlar üçün). Cədvəllər (quyular,
geologiya, faultlar, SWOF) modeldən və layihədən bərpa olunur.

Rədd edilən variant: hər panelə əl ilə `get_state/set_state` yazmaq —
`_load_model_into_panels`-in taleyi göstərdi ki, belə siyahılar köhnəlir.

### Qərar 2 — modelin sadə dataclass-ları BÜTÜN sahələri ilə yazılır

`_all_fields` + `_known_fields`: yeni sahə avtomatik yazılır, köhnə
faylda olmayan sahə defolt alır, naməlum açar atılır.

### Qərar 3 — bərpa faylı AYRICA, istifadəçinin faylına avtomatik yazılmır

Avtomatik saxlama istifadəçinin qəsdən saxladığı variantı səssizcə
silə bilərdi. Bərpa faylı yalnız qəfil bağlanmaya qarşıdır və düzgün
bağlanmada silinir.

### Qərar 4 — son layihə açılışda avtomatik açılmır

Yeni model qurmaq istəyəndə mane olur, böyük fayl açılışı yavaşladır.
Əvəzində «Son layihələr» menyusu (5 fayl).
