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
