# İş hesabatı — Lay Simulyatoru

Bu fayl işin **gedişini** qeyd edir: nə edildi, nə vaxt, niyə, hansı
çətinliklə qarşılaşdıq və nə buraxıldı.

**Qayda:** hər iş seansından sonra yeni bölmə əlavə olunur — köhnə
qeydlər silinmir, dəyişdirilmir. Bu, layihənin yaddaşıdır.

---

## 10 sentyabr 2026 — Seans 1: Layihənin qurulması

### Görülən iş

1. **GitHub repo-su bağlandı**
   `https://github.com/eoraieor-tech/LAY-SIMULYATIR-MODELI-` (private, boş idi)
   Lokal qovluq: `C:\Users\LOQ\LaySimulyator`
   Branch: `main` · Remote: `origin`

2. **Sənədləşdirmə strukturu quruldu**
   - `README.md` — giriş nöqtəsi və naviqasiya
   - `ARCHITECTURE.md` — struktur (boş şablon, doldurulacaq)
   - `ROADMAP.md` — mərhələlər
   - `QARARLAR.md` — qərar jurnalı
   - `ISH_HESABATI.md` — bu fayl
   - `CLAUDE.md` — iş qaydaları
   - `docs/` — dərin mövzular üçün
   - `.gitignore`

### Qərarlar

- **Sənədlər Azərbaycan dilində yazılır.** Səbəb: komandanın iş dili budur.
  Kod və dəyişən adları ingiliscə qalacaq (standart praktika).
- **Şablonlar boş yaradıldı, uydurulmadı.** Layihənin texniki detalları
  hələ məlum deyil — uydurma məzmun sonradan səhv istiqamət verə bilər.
  Doldurulmamış yerlər ⏳ işarəsi ilə açıq göstərilib.

### Bilinməyənlər (növbəti seansda aydınlaşmalıdır)

| # | Sual |
|---|---|
| 1 | Simulyatorun məqsədi nədir — hansı prosesləri modelləşdirir? |
| 2 | Hansı proqramlaşdırma dili / texnologiya? |
| 3 | İnterfeys olacaq, yoxsa yalnız hesablama nüvəsi? |
| 4 | 2D, yoxsa 3D? Neçə faza (su/neft/qaz)? |
| 5 | Hazır layihədən (IMEX2D) nəsə götürülür, yoxsa tam təmiz başlanğıc? |

### Buraxılan iş

Heç nə — bu seansda yalnız struktur quruldu, kod yazılmadı.
Bu qəsdəndir: tələblər məlum olmadan kod yazmaq boş işdir.

---

## 10 sentyabr 2026 — Seans 2: Tələblər, mühit və mövcud kodun auditi

### Görülən iş

#### 1. Tələblər alındı

Sahibkar tam texniki plan verdi: **3D 3-fazalı (neft-su-qaz) hidrodinamik
lay simulyatoru + geoloji interpolyasiya platforması**, Python 3.11+ /
VS Code, **MPFA-O** diskretizasiyası (TPFA əvəzinə), 5-spot quyu şəbəkəsi,
RF/THP/BHP hesablanması, PyVista ilə 3D vizualizasiya. Altı mərhələlik
xəritə ilə birlikdə.

Bu, Faza 0-ın 1-4-cü açıq suallarını bağladı.

#### 2. Mühit yoxlandı — üç problem tapıldı

| Tapıntı | Nəticə |
|---|---|
| Python 3.11 **yoxdur** (3.14.7 və 3.12 var) | 3.12 seçildi (Q-03) |
| C++ alət zənciri yoxdur (MSVC, CMake, MinGW, LLVM) | C++ variantı praktiki deyil |
| **Smart App Control işlək** (`VerifiedAndReputablePolicyState = 1`) | imzasız DLL-lər bloklanır |

Smart App Control ölçmələri:

- Yeni yaradılan venv-lərdə (həm OneDrive-da, həm `C:\Dev`-də)
  `scipy.sparse.linalg` (`_superlu`), `scipy.interpolate` (`_specfun`),
  `scipy.spatial` (`_special_ufuncs`), PyKrige və VTK **bloklanır**
- **Sistem Python 3.14-də** isə SciPy-ın **bütün** modulları işləyir;
  yalnız VTK bloklanır
- Səbəb: köhnə fayllar reputasiya qazanıb, yeni yazılanlar yox

Yol boyu **iki səhv fərziyyəm oldu və düzəldildi**:

1. Əvvəlcə blokun səbəbini OneDrive qovluğunda gördüm — səhv idi.
   Yalnız `scipy.linalg` sınanmışdı; `scipy.sparse.linalg` sınananda
   məlum oldu ki, yeni yerdə də bloklanır.
2. Sonra "yalnız vizualizasiya bloklanır" dedim — bu da səhv idi;
   həqiqətdə xətti həlledici də bloklanırdı (yeni venv-də).

#### 3. Mövcud layihə tapıldı və auditdən keçirildi

Sahibkar bildirdi ki, proqramın bir hissəsi artıq yazılıb. `Downloads`
qovluğunda **14 zip versiyası** tapıldı (31 avqust → 10 sentyabr).
Ən yenisi (`(13).zip`) auditə verildi.

Bu, Faza 0-ın 5-ci sualını bağladı: **"IMEX2D" məhz bu layihədir.**

**Ölçülənlər:** 65 175 sətir Python (28 509-u test, 96 fayl), 27 sənəd,
təmiz hexagonal arxitektura, 0 TODO / 0 FIXME / 0 HACK / 0 çılpaq
`except:`.

**Testlər real icra edildi** (fərziyyə ilə deyil): **1 993 keçdi,
1 uğursuz, 1 ötürüldü, 1 xfail** — 7 dəq 48 san.

**Tapılan səhvlər:**

- **§8.3** — `pytest.importorskip("vtk")` pytest 8.2+ ilə işləmir
  (yalnız `ModuleNotFoundError`-da ötürür, burada `ImportError` gəlir).
  Nəticə: toplama dayanır, **1 995 testin heç biri işləmir**.
- **§8.4** — uğursuz test diaqnoz edildi: səhv **işlək kodda deyil,
  testin özündədir** (`_cell_centres`-ə `geometry` yerinə `spec`).
  Düzəliş tətbiq olunub yoxlanıldı → **23/23 keçdi**.
- **§8.1** — **MPFA-O yazılıb (1 699 sətir), test edilib, amma
  mühərrikə qoşulmayıb.** Yalnız testlər idxal edir; `impes_engine.py`
  MPFA-O ötürülsə qəsdən xəta atır; `jacobian.py` hər üz üçün dəqiq
  2 hüceyrə fərz edir.
- **§8.5** — THP / VFP modulu yoxdur (planda tələb olunur).

#### 4. Qaz fazası bərpa edildi

`A7_PLAN.md` göstərdi ki, üç fazalı black-oil işi demək olar
bitirilmişdi, sonra **v69-da tamamilə silinib** — çünki quyu öz BHP
həddinə yaxınlaşanda Nyuton yığılmırdı və kök tapılmamışdı. Bunun
üzərinə strateji dönüş edilib: fizika OPM Flow-a həvalə olunsun.

Kod 14 zip-in **heç birində yox idi**. Lakin GitHub reposu
(`eoraieor-tech/lay-simulyatoru`) public çıxdı və reponun **ilk
commit-i** məhz bu məqsədlə saxlanılıb:

```
66deca8   v67 - son 3 fazali versiya (qaz fazasi + matplotlib 3D)
```

Oradan **12 fayl / 4 118 sətir** çıxarıldı (2 165 mənbə + 1 953 test)
və `berpa/A7_qaz_fazasi/` qovluğunda saxlanıldı — Stone II, üç fazalı
qalıq tənlikləri, analitik Jakobian, Nyuton döngəsi, dəyişən keçidi.

### Qərarlar

- **Q-03** — Python 3.12 (3.11 yoxdur, 3.14 çox yenidir)
- **Q-04** — mövcud kod bazası atılmır, üzərində davam edilir
- **Q-05** — bu layihə üçün yeni venv yaradılmır; sistem Python 3.14
  işlədilir
- **Q-06** — Smart App Control söndürülmür

### Buraxılan iş

- Yeni repoda kod yazılmadı — audit nəticəsi strateji qərarı
  dəyişdirdiyi üçün qəsdən dayandırıldı
- `src/` altındakı boş modul skeletləri qaldı; hansı repoda davam
  ediləcəyi həll olunandan sonra taleyi müəyyənləşəcək
- VTK bloku həll olunmadı (qərar gözlənilir)
- MPFA-O qoşulmadı, qaz fazası əsas koda qaytarılmadı

### Açıq suallar (növbəti seansda həll olunmalıdır)

| # | Sual |
|---|---|
| 1 | **Öz fizikamız, yoxsa OPM Flow?** v69 qərarı yeni planla ziddiyyətdədir |
| 2 | **Hansı repoda davam edirik?** `lay-simulyatoru` (65k sətir) yoxsa `LAY-SIMULYATIR-MODELI-` (boş) |
| 3 | Qovluq təkrarı: `OneDrive\Desktop\L.S.M\...` və `C:\Dev\LSM` — hansı qalır |
| 4 | VTK bloku: Plotly/matplotlib, WSL2, yoxsa SAC söndürülsün |

Tam audit: [AUDIT_2026-09-10.md](AUDIT_2026-09-10.md)
