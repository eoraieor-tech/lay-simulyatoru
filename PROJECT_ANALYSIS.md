# IMEX-2D — Layihənin tam texniki analizi

**Tarix:** 26 sentyabr 2026 · **Commit:** `eacb029` (main) · **Analizi aparan:** Claude (AI köməkçi)
**Qısa xülasə:** [PROJECT_ANALYSIS_SUMMARY.md](PROJECT_ANALYSIS_SUMMARY.md) · **Diaqramlar:** [docs/diagrams/](docs/diagrams/README.md)

> **Oxuma qaydası.** Hər iddia konkret fayl və sətirlə əsaslandırılıb.
> Mətndə üç işarə var:
> **[FAKT]** — kodda və ya əmr nəticəsində birbaşa görünür;
> **[RƏY]** — analitikin qiymətləndirməsi, mübahisə oluna bilər;
> **[⏳ ƏMİN DEYİL]** — yoxlanmayıb və ya yoxlana bilmədi.
> Heç bir kod dəyişdirilməyib, heç bir fayl silinməyib.

---

## Mündəricat

0. [Analizin metodu: oxunan fayllar və işlədilən əmrlər](#0-analizin-metodu)
1. [Executive Summary](#1-executive-summary)
2. [Layihənin məqsədi](#2-layihənin-məqsədi)
3. [Qovluq və fayl strukturu](#3-qovluq-və-fayl-strukturu)
4. [Arxitektura](#4-arxitektura)
5. [Məlumat axını](#5-məlumat-axını)
6. [Əsas istifadə ssenariləri (addım-addım)](#6-əsas-istifadə-ssenariləri)
7. [Texnologiyalar və kitabxanalar](#7-texnologiyalar-və-kitabxanalar)
8. [İstifadə olunan modellər](#8-istifadə-olunan-modellər)
9. [Nəzəri əsaslar](#9-nəzəri-əsaslar)
10. [Əsas fayllar və funksiyalar](#10-əsas-fayllar-və-funksiyalar)
11. [Konfiqurasiya və idarəetmə](#11-konfiqurasiya-və-idarəetmə)
12. [Təhlükəsizlik analizi](#12-təhlükəsizlik-analizi)
13. [Test analizi](#13-test-analizi)
14. [Problemlər və risklər](#14-problemlər-və-risklər)
15. [Təkmilləşdirmə təklifləri](#15-təkmilləşdirmə-təklifləri)
16. [Layihəni öyrənmə planı](#16-layihəni-öyrənmə-planı)
17. [Terminlər lüğəti](#17-terminlər-lüğəti)
18. [Açıq qalan suallar](#18-açıq-qalan-suallar)
19. [Yekun nəticə](#19-yekun-nəticə)

---

## 0. Analizin metodu

### 0.1 İşlədilən əmrlər və nəticələri

| # | Əmr (qısaldılmış) | Nə üçün | Nəticə |
|---|---|---|---|
| 1 | `git ls-files` | Versiya nəzarətindəki bütün fayllar | 389 fayl |
| 2 | `ls -la` (kök) | İzlənməyən faylları görmək | `gunluk.csv` (izlənmir), `layihe.YOXLANIS 1.imx` (7.6 MB, `.gitignore`-da), `logs/`, `graphify-out/` |
| 3 | `wc -l` bütün `.py` faylları | Ölçü | `imex2d/` 46 240 sətir, `tests/` 38 863, `tools/` 1 199 |
| 4 | `graphify query "main data flow from UI to simulation engine"` | Layihənin bilik qrafiki ilə istiqamətlənmə (CLAUDE.md tələbi) | 994 node; `MainWindow → ModelAwareSimulationService → FullyImplicitEngine/ImpesEngine` zənciri |
| 5 | `head graphify-out/GRAPH_REPORT.md` | God node-lar | `default_scal()`, `OrdinaryKriging`, `BlackOilPVTProvider`, `SimulationConfig`, `ReservoirModel` … |
| 6 | `grep -nE "^(class \|def )"` bütün modullarda | Sinif/funksiya xəritəsi | 857 sətirlik siyahı (skretçpeddə) |
| 7 | Python AST skripti — paketlər arası `import` sayı | Arxitektura qaydasının yoxlanması | §4.3-də cədvəl; 2 qayda pozuntusu |
| 8 | `grep` — `eval`, `exec`, `pickle`, `subprocess`, `os.system`, `shell=True`, `yaml.load`, `urllib`, `requests`, `socket` | Təhlükəsizlik | `imex2d/`-də **heç biri yoxdur**; yalnız `run_tests.py`-də `importlib` |
| 9 | `grep` — `os.environ`, `QSettings`, `LOCALAPPDATA` | Konfiqurasiya mənbələri | `IMEX2D_DATA_DIR`, `LOCALAPPDATA`, `IMEX_SKIP_SLOW`, `QT_QPA_PLATFORM` |
| 10 | `ls pyproject.toml Dockerfile .env.example .github …` | Build/deploy faylları | **Heç biri yoxdur** |
| 11 | `../venv/Scripts/python.exe -m pip list` | Faktiki quraşdırılmış versiyalar | Python 3.14.7; `requirements.txt`-dəki versiyalarla üst-üstə düşür; `resdata`, `coverage`, `pip-audit` **yoxdur** |
| 12 | `diff -q berpa/A7_qaz_fazasi/imex2d/... imex2d/...` | Bərpa qovluğunun rolu | 3 fayl eyni, 3 fayl canlı kodda xeyli böyüyüb — `berpa/` arxivdir |
| 13 | `grep -c "ERROR" logs/imex2d.log` | Real istifadədə xətalar | 1 xəta: «VTK çəkilişi uğursuz — matplotlib-ə qayıdılır» |
| 14 | `python -m pytest -q -p no:cacheprovider --tb=line -o addopts=""` | Testlərin faktiki vəziyyəti | **2761 keçdi, 1 atlandı, 1 xfailed, 0 uğursuz** (8 dəq 47 san) — §13.4 |
| 15 | `python -m pytest --co -q` | Toplanan testlərin sayı | 2762 test |
| 16 | `python -m pytest -q -rsx tests/test_opm_import.py tests/test_implicit_newton.py` | Atlanan/xfail səbəbi | «resdata quraşdırılmayıb»; bilinən xfail (§13.4) |
| 17 | `grep` — `import csv`, `PIL`, `pandas` | Data kitabxanaları | `csv` 6 faylda, `PIL` `animation_export.py`-də, `pandas` yoxdur |

### 0.2 Oxunan fayllar

**Tam oxunub:** `app.py`, `requirements.txt`, `requirements-dev.txt`, `pytest.ini`,
`run.bat`, `test.bat`, `run_tests.py`, `.gitignore`, `AGENTS.md`, `CLAUDE.md`,
`README.md`, `ARCHITECTURE.md`, `ROADMAP.md`, `docs/is_axini_ardicilligi.md`,
`.claude/settings.json`, `.codex/hooks.json`,
`imex2d/application/{simulation_service,config,session,project,model_builder}.py`,
`imex2d/simulation/implicit/{engine,newton,time_stepping,linear}.py`,
`imex2d/simulation/linear_solver.py`, `imex2d/logging_setup.py`, `imex2d/version.py`,
`imex2d/ui/worker.py`, `imex2d/simulation/results.py`.

**Qismən oxunub (əsas hissələr):** `ISH_HESABATI.md` (Seans 46–47),
`docs/nezeri_esaslar.md` (§0–§2), `tests/README.md`,
`imex2d/application/{serialization,geology_service}.py`,
`imex2d/ui/main_window.py` (run/save/open/worker bağlantıları),
`imex2d/simulation/impes_engine.py`, `imex2d/simulation/implicit/three_phase_engine.py`,
`imex2d/io/{grdecl,grdecl_import}.py`, `imex2d/history/optimizer.py`,
`imex2d/domain/{wells,reservoir_model}.py`, `imex2d/simulation/pvt/correlations.py`,
`imex2d/simulation/scal_adapter.py`, `imex2d/simulation/well_model.py`,
`tools/*.py` (başlıqlar), `berpa/A7_qaz_fazasi/MENBE.md`, nümunə CSV-lər.

**Yalnız strukturu (sinif/funksiya siyahısı) oxunub:** qalan bütün `imex2d/**/*.py` modulları.

> **[⏳ ƏMİN DEYİL]** 86 000+ sətirlik kod bazasının hər sətri oxunmayıb.
> Oxunmamış modullar haqqında deyilənlər onların docstring-lərinə, sinif/funksiya
> adlarına və layihənin öz sənədlərinə (`ARCHITECTURE.md`, `docs/nezeri_esaslar.md`)
> əsaslanır. Sənəddəki düsturlar üç yerdə kodla nümunəvi tutuşduruldu
> (Standing — `correlations.py:36-43`, Corey törəməsi — `scal_adapter.py`,
> Peaceman — `well_model.py:184-186`) və üçü də uyğun gəldi.

---

## 1. Executive Summary

**IMEX-2D** — neft-qaz yatağının (rezervuarın) **masaüstü simulyatorudur**.
Python + PyQt5 ilə yazılıb, Windows-da işləyir, internetə, verilənlər bazasına
və ya süni intellekt modelinə qoşulmur. Kəşfiyyat quyularının ölçmələrindən
3D geoloji model qurur (Kriging, SGS, SIS), onu hidrodinamik modelə çevirir
və neft/su/qaz axınını zamana görə hesablayır (black-oil, IMPES və tam implicit
Nyuton), nəticəni qrafik, 3D görüntü, günlük cədvəl, CSV/JSON/PDF kimi verir.

**Güclü tərəflər [RƏY, faktlarla]:**

- Aydın qatlı arxitektura (`domain` heç nədən asılı deyil — §4.3-də AST ilə
  təsdiqləndi), dependency injection (`app.py` composition root).
- Fizika ciddi: analitik Jakobian, CNV + kütlə balansı meyarları, adaptiv Δt,
  CPR ön-şərtçisi, MPFA-O (tam tenzor), corner-point həndəsə, SPE1 etalonu ilə
  müqayisə aləti.
- Çox geniş test dəsti — 137 test faylı, ~2 500 test funksiyası (parametrləşdirmə
  ilə 2 700+ test), analitik həllərlə (Buckley-Leverett) və qızıl etalonlarla
  (`tests/golden/*.json`).
- Hər qərarın səbəbi yazılıb (`QARARLAR.md`, `ISH_HESABATI.md` — 5 400+ sətir).

**Əsas risklər [RƏY]:**

1. **UI monoliti** — `ui/main_window.py` 3 155 sətir, 121 metod; pytest-də qurula
   bilmir (VTK səhnəsi qaçırıcını çökdürür, `application/session.py:3-5`), yalnız
   AST ilə statik yoxlanılır.
2. **Layihə faylı atomik yazılmır** — `serialization.py:375` birbaşa hədəf fayla
   yazır; yazma zamanı çökmə istifadəçinin yeganə `.imx` faylını korlaya bilər.
3. **Üç fazalı Jakobianda bilinən qeyri-dəqiqliklər** (`ROADMAP.md` TB-1: nisbi
   xəta 0.49; TB-2: 0.19) və SPE1CASE2-də OPM Flow-dan ~145 günlük qaz cəbhəsi fərqi.
4. **Konfiqurasiyanın bir hissəsi səssizcə təsirsizdir** — `LinearSolverConfig`
   tam implicit mühərriklərdə işlədilmir (`implicit/engine.py:60,84-86`).
5. **Build/CI yoxdur** — `pyproject.toml`, CI, coverage ölçməsi, asılılıq
   zəiflik yoxlaması yoxdur; tam test dəsti 9–13 dəqiqə çəkir (bu analizdə: 2761 keçdi, 0 uğursuz).

**Kritik (dərhal təhlükəli) problem tapılmadı.** Ən ciddi 4 problem «High» kimi
qiymətləndirildi (§14).

---

## 2. Layihənin məqsədi

| Sual | Cavab | Əsas |
|---|---|---|
| **Əsas məqsəd nədir?** | Quyu məlumatından 3D heterogen lay modeli qurub orada neft/su/qazın hərəkətini günbəgün simulyasiya etmək; RF, debit, BHP/THP, lay təzyiqini hesablamaq və göstərmək. | `README.md` §1; `version.py:15-36` |
| **Hansı problemi həll edir?** | Yataq mühəndisinin sualı: «bu quyu sxemi ilə nə qədər neft çıxara bilərik və nə vaxt su gələcək?» Kommersiya simulyatorları (Eclipse, CMG IMEX) bahalı və qapalıdır; bu layihə şəffaf, öyrənilə bilən və sahibkarın öz nəzarətində olan alternativdir. **[RƏY]** — ikinci hissə sənədlərdən çıxarılıb, açıq yazılmayıb. | `README.md` §1; `QARARLAR.md` |
| **Əsas funksiyalar** | (1) quyu CSV-dən və ya sintetik geologiya; (2) GRDECL/Eclipse deck idxalı; (3) Kriging/IDW/SGS/SIS; (4) PVT korrelyasiyaları və ya deck cədvəlləri; (5) Corey/SWOF/SGOF, Stone II; (6) IMPES, tam implicit 2 və 3 fazalı mühərriklər; (7) TPFA və MPFA-O; (8) BHP/RATE/THP quyu idarəsi, BHP limiti, səth debiti; (9) 2D/3D görüntü, animasiya, GIF; (10) günlük göstəricilər; (11) CSV/JSON/PDF/Eclipse ixracı; (12) history matching və həssaslıq; (13) layihə faylı və qəza bərpası. | `version.py:15-36`, `EXPECTED_TABS` (`version.py:38-43`) |
| **Kim üçün?** | Rezervuar mühəndisi, geoloq, tələbə/tədqiqatçı (sənədlərdə «diplom» qeyd olunur — `ISH_HESABATI.md` Seans 47 §4). | `ISH_HESABATI.md` |
| **Hansı ssenarilərdə?** | Su vurma ilə hasilat proqnozu (5-spot), qaz papağı/həll olmuş qazla üç fazalı proqnoz, quyu sxemlərinin müqayisəsi, tarixi hasilata uyğunlaşdırma, parametr həssaslığı, tədris/təqdimat. | `application/scenarios.py:34-70`; `history/` |
| **Giriş məlumatları** | Quyu CSV (`well,x,y,k,PORO,PERMX,NTG`), GRDECL, Eclipse PVT/SCAL deck, SCAL CSV, fay CSV, müşahidə CSV (`time,well,quantity,value`), UI panellərindəki parametrlər, `.imx` layihə faylı. | `numune_*.csv`; `io/`; `geology/well_data_io.py:90` |
| **Çıxış məlumatları** | `SimulationResult` (zaman sıraları, quyu sıraları, 3D anlar), qrafiklər, 3D səhnə, günlük cədvəl, CSV/JSON, PDF hesabat, Eclipse `.DATA`, PNG/GIF, `.imx`, `logs/imex2d.log`. | `simulation/results.py:43`; `reporting/`; `io/eclipse_export.py` |
| **Əsas iş prinsipi** | Hər hüceyrə üçün kütlə balansı tənliyi yazılır (girən − çıxan = yığılan), axın Darcy qanunu ilə hesablanır, qeyri-xətti sistem hər zaman addımında Nyuton üsulu ilə həll olunur. | `simulation/implicit/residual.py`, `newton.py:1-23` |

---

## 3. Qovluq və fayl strukturu

```
IMEX2D/
├── app.py                  Giriş nöqtəsi — composition root (61 sətir)
├── run.bat / test.bat      Windows başladıcıları (venv aktivləşdirir)
├── run_tests.py            pytest olmayanda ehtiyat test qaçırıcısı
├── requirements*.txt       Asılılıqlar (dəqiq versiyalarla)
├── pytest.ini              testpaths=tests, marker: performance
├── numune_*.csv/.GRDECL    Nümunə giriş faylları (5 ədəd)
├── *.md (29 ədəd)          Layihə sənədləri — AZ dilində
├── imex2d/                 ƏSAS PAKET — 46 240 sətir, 145 modul
│   ├── domain/        27 fayl  7 059 s.  Model obyektləri (grid, həndəsə, xassələr, quyu, PVT, SCAL)
│   ├── interfaces/     6       401     ABC müqavilələri (IPVTProvider, ILinearSolver …)
│   ├── application/   10     4 398     Servislər, konfiqurasiya, layihə, serializasiya, geologiya boru xətti
│   ├── geology/       19     8 253     Geostatistika: variogram, kriging, SGS, SIS, QC, CV
│   ├── simulation/    42    11 046     Mühərriklər, PVT, SCAL, quyu, lülə hidravlikası, ilkin şərt
│   ├── discretization/ 5     1 699     MPFA-O (çoxnöqtəli axın)
│   ├── io/             9     2 293     GRDECL, PVT/SCAL deck, fay, Eclipse ixrac/summary, OPM
│   ├── history/        6     1 283     Uyğunsuzluq, optimallaşdırma, həssaslıq
│   ├── rendering/      6     2 610     matplotlib renderer-lər, VTK 3D, animasiya
│   ├── reporting/      4       842     PDF, CSV/JSON, günlük cədvəl
│   ├── ui/             9     5 913     PyQt5 pəncərə, panellər, fon işçiləri
│   ├── benchmarks/     2       288     SPE1CASE2 qurucusu
│   ├── logging_setup.py, version.py
├── tests/                  137 test faylı, 38 863 sətir; golden/*.json etalonlar
├── tools/                  8 skript: benchmark, SPE1 müqayisə, modul qrafiki, təqdimat demo
├── docs/                   Dərin mövzular (MPFA-O, nəzəriyyə, iş axını), teqdimat/ (HTML slaydlar, PDF)
├── berpa/A7_qaz_fazasi/    Git tarixçəsindən bərpa olunmuş köhnə qaz kodu — ARXİV, işlədilmir
├── graphify-out/           Bilik qrafiki (generasiya olunur, .gitignore-da)
├── logs/                   imex2d.log (.gitignore-da)
├── .claude/ .codex/        AI köməkçi alətlərinin hook-ları (graphify)
```

**[FAKT]** `berpa/` `pytest.ini`-nin `testpaths`-ına daxil deyil və heç bir
canlı modul onu idxal etmir; `three_phase_state.py`, `stone_relperm.py`,
`domain/three_phase.py` canlı kodla eynidir, qalan 3 fayl canlı kodda xeyli
böyüyüb (məs. `three_phase_residual.py` 1 066 → 1 401 sətir).

---

## 4. Arxitektura

### 4.1 Arxitektura modeli

| Model | Uyğundurmu? | Əsas |
|---|---|---|
| **Monolitik** | **Bəli** — tək proses, tək paket, tək masaüstü tətbiq | `app.py` |
| **Layered (qatlı)** | **Bəli, əsas model** — `domain → interfaces → simulation/geology → application → rendering/ui` | `ARCHITECTURE.md` §2; §4.3 |
| **Modular** | **Bəli** — hər mövzu ayrı paketdə, müqavilələr `interfaces/`-də | `imex2d/interfaces/providers.py` |
| **Ports & Adapters / DI** | **Bəli** — provider-lər konstruktorla ötürülür, konkret siniflər yalnız `app.py`-də və `ModelAwareSimulationService`-də seçilir | `simulation_service.py:54-73`, `app.py:47-54` |
| **Pipeline** | **Bəli, məntiq səviyyəsində** — geologiya boru xətti (QC → çevirmə → variogram → kriging → CV) və «geologiya → model → mühərrik → nəticə» | `geology_service.py:441`; `docs/is_axini_ardicilligi.md` |
| **Event-driven** | **Qismən** — yalnız UI-də Qt siqnal/slot (`progress`, `finished_ok`, `failed`) | `ui/worker.py:103-105` |
| **MVC** | **Qismən / qeyri-rəsmi** — Model = `domain/`, View = `rendering/` + Qt vidjetlər, Controller = `MainWindow`. Lakin `MainWindow` həm View, həm Controller rolunu daşıyır. **[RƏY]** | `ui/main_window.py` |
| Microservice, Agent-based | **Xeyr** | Şəbəkə kodu yoxdur (§0.1 əmr 8) |

### 4.2 Modullar və məsuliyyətlər

| Paket | Məsuliyyət | Rol |
|---|---|---|
| `domain/` | Saf məlumat obyektləri: `CartesianGrid`, `CellGeometry`/`CornerPointGeometry`, `PropertyMap`, `PermeabilityTensor`, `ReservoirModel`, `GeologicalModel`, `Well`, `PVTTable`, `CoreyParameters`, vahid çevirmələri, validasiya | «Verilənlər modeli» (DB yoxdur — hər şey yaddaşda) |
| `interfaces/` | ABC müqavilələri: `IPVTProvider`, `IRelativePermeabilityProvider`, `ICapillaryPressureProvider`, `IInitializationProvider`, `IWellHydraulicsProvider`, `ILinearSolver`, `IProgressReporter`, `ISimulationEngine`, `IFluxDiscretization`, `IPropertyInterpolator` | Portlar |
| `simulation/` | Hesablama nüvəsi: 3 mühərrik, PVT provider və korrelyasiyalar, SCAL provider-ləri, kapilyar, Peaceman, quyu qaydaları, lülə hidravlikası, ilkin tarazlıq, analitik B-L | **«Backend» (hesablama)** |
| `discretization/` | MPFA-O lokal və qlobal operatorları | Hesablama |
| `geology/` | Geostatistika | Hesablama |
| `application/` | İş axını: `SimulationService`, `WellBasedGeologicalModelBuilder`, `ReservoirModelBuilder`, `Project`, `ProjectSerializer`, `SimulationConfig`, seans | **Servis qatı** |
| `history/` | Parametr dəyişdirici, uyğunsuzluq, optimallaşdırıcı, həssaslıq | Servis |
| `io/` | Fayl formatları (oxu/yaz) | Adapterlər |
| `rendering/`, `reporting/` | Qt-siz çəkmə və ixrac | Çıxış adapterləri |
| `ui/` | PyQt5 pəncərə, 10 panel, 3 fon işçisi | **«Frontend»** |
| `app.py` | Asılılıqları bağlayır, loglamanı qurur, pəncərəni açır | Composition root |

**API, verilənlər bazası, xarici servis — yoxdur [FAKT]** (§0.1 əmr 8, 10).
Yeganə «API» daxili Python API-dir — `ARCHITECTURE.md` §7 skript rejimini göstərir.

### 4.3 Faktiki asılılıqlar (AST ilə sayılıb)

Diaqram: [docs/diagrams/01_qat_arxitekturasi.md](docs/diagrams/01_qat_arxitekturasi.md).

| Haradan → hara | import sayı | Qaydaya uyğun? |
|---|---|---|
| `application → domain / geology / simulation / interfaces` | 52 / 17 / 13 / 3 | ✅ |
| `simulation → domain / interfaces` | 36 / 19 | ✅ |
| `ui → application / rendering / domain / io / history / reporting` | 14 / 9 / 22 / 6 / 5 / 4 | ✅ |
| `domain → (heç nə)` | 0 | ✅ tam qorunur |
| **`simulation → application`** | **4** (`impes_engine.py:17`, `linear_solver.py:14`, `implicit/engine.py:18`, üç fazalı mühərrik) | ❌ `config.py` idxalı — `application ↔ simulation` dövri asılılığı |
| **`ui → simulation`** | **6** (`main_window.py:78-80` və s.) | ⚠️ mühərrik sinfi UI-də seçilir (`main_window.py:2120-2123`, `1412-1414`, `1491-1494` — eyni kod 3 dəfə) |
| `rendering → history`, `reporting → history` | 1 / 3 | ✅ (history servis qatıdır) |

### 4.4 Modullar necə əlaqə qurur

1. **Konstruktor inyeksiyası** — `app.py:47-54` servisə Corey adapterini və
   CG+ILU həlledicini, pəncərəyə servisi, geologiya və model qurucusunu verir.
2. **Model-yönlü provider seçimi** — `ModelAwareSimulationService.create_engine`
   (`simulation_service.py:262-271`) hər qaçışda provider-ləri modelin özündən qurur:
   SCAL cədvəli varsa `TableRelativePermeabilityProvider`, yoxsa Corey
   (`:401-408`); PVT cədvəli varsa `BlackOilPVTProvider` (`:344-367`); PVT-də qaz
   varsa **avtomatik** üç fazalı mühərrik + Stone II (`:269-341`).
3. **Fon axını** — `SimulationWorker(QThread)` servisin `run()`-unu çağırır,
   proqres `QtProgressReporter` vasitəsilə Qt siqnalına çevrilir (`worker.py:14-21`).
4. **Qt-siz çıxış** — renderer-lər hazır `matplotlib.Axes` alır (`ARCHITECTURE.md` §5.10),
   ona görə eyni kod PDF-də və testdə işləyir.

---

## 5. Məlumat axını

Diaqramlar: [02_melumat_axini.md](docs/diagrams/02_melumat_axini.md),
[03_simulyasiya_ardicilligi.md](docs/diagrams/03_simulyasiya_ardicilligi.md),
[04_zaman_addimi_nyuton.md](docs/diagrams/04_zaman_addimi_nyuton.md).

```
GİRİŞ (CSV · GRDECL · deck · UI panelləri · .imx)
  ↓ bir dəfə
GEOLOGİYA   QC → çevirmə → variogram → anizotropluq → qonşuluq → kriging/SGS/SIS → CV
  ↓
MODEL       GeologicalModel → ReservoirModelBuilder.build() → ReservoirModel.validate()
  ↓ bir dəfə (mühərrik qurulanda)
ÖN-HESAB    T_ij (TPFA/MPFA-O) · PV · WI (Peaceman) · ilkin p, Sw, Sg
  ↓
ZAMAN       hər addım: RATE payları → Nyuton (PVT→kr→Pc→Φ→upstream→axın→R→J→δ)
            → THP → səth debiti → BHP limiti → qəbul → sıralar → adaptiv Δt
  ↓
NƏTİCƏ      THP post-prosesi → SimulationResult → qrafik / 3D / günlük / CSV / PDF / .imx
```

**İstifadəçi sorğusu hansı mərhələlərdən keçir** («İşə sal» düyməsi, `main_window.py:2115-2160`):

1. `rebuild_model()` — panellərdən `ReservoirModel` yığılır.
2. `numerical_panel.simulation_config()` — `SimulationConfig`.
3. Mühərrik sinfi seçilir və `service.with_engine(factory)` servisin **surətini** yaradır.
4. `reservoir_model.diagnose()` — xəta varsa dayanır, xəbərdarlıq varsa soruşur.
5. `service.create_engine()` — **ön yoxlama** (model + konfiqurasiya validasiyası, uyğunsuz kombinasiyalar: IMPES+MPFA-O, IMPES+THP, IMPES+BHP limiti, IMPES+səth debiti, 3 faza+MPFA-O).
6. `project.new_run()` → `RUN-00N`.
7. `SimulationWorker.start()` → fon axınında `service.run()` → mühərrik **yenidən** qurulur və işə düşür.
8. Mühərrik hər `progress_every_n_steps` addımda proqres bildirir; istifadəçi dayandıra bilər.
9. Nəticə `finished_ok` siqnalı ilə UI-yə qayıdır → bütün görüntülər yenilənir → bərpa faylı yazılır.

---

## 6. Əsas istifadə ssenariləri

### 6.1 Ssenari A — Simulyasiyanı işə salmaq (əsas ssenari)

| Addım | Nə baş verir | Fayl / funksiya |
|---|---|---|
| 1. Giriş | İstifadəçi panellərdə grid, süxur, PVT, SCAL, quyu, ədədi parametrləri doldurur | `ui/panels.py` — `GridGeometryPanel:115`, `RockFluidPanel:1313`, `ScalPanel:1463`, `PvtPanel:1574`, `WellPanel:1694`, `NumericalPanel:2069` |
| 2. Ötürmə | «İşə sal» → `MainWindow.run_simulation()` | `main_window.py:2115` |
| 3. Preprocessing | Panellərdən domain obyektləri (`ScalPanel.values() → CoreyParameters` və s.); vahidlər mühərrik vahidinə (bar, m, mD, cP, gün) çevrilir | `main_window.py:1704` `rebuild_model`; `domain/unit_conversions.py:238` `to_engine_units` |
| 4. Validasiya | (a) `ReservoirModel.diagnose()` — xəta/xəbərdarlıq; (b) `ReservoirModel.validate() + SimulationConfig.validate()`; (c) uyğunsuz seçimlərin istifadəçi dilində rəddi | `reservoir_model.py:179,183`; `config.py:76-89`; `simulation_service.py:75-167` |
| 5. Provider seçimi | Model məlumatına görə Corey/cədvəl, PVT, kapilyar, ilkin şərt; qaz varsa 3 fazalı mühərrik | `simulation_service.py:262-423` |
| 6. Mühərrik qurulması | Diskretizasiya (`TwoPointFluxDiscretization.build` və ya `MPFAODiscretization`), Peaceman bağlantıları, qalıq/Jakobian, Nyuton, adaptiv addım, quyu nəzarətçiləri, ilkin vəziyyət | `implicit/engine.py:56-122` |
| 7. Emal | Zaman dövrəsi (§5; diaqram 04) | `implicit/engine.py:169-290` |
| 8. Nəticə | `SimulationResult` — sıralar, quyu sıraları, anlar, mesaj | `simulation/results.py:43` |
| 9. Post-proses | Lülə həndəsəsi olan quyular üçün THP traversi | `simulation_service.py:211-227`; `wellbore/hydraulics.py` |
| 10. Göstərmə/saxlama | Qrafiklər, xəritə, 3D, müqayisə; bərpa faylı `%LOCALAPPDATA%\IMEX2D\berpa.imx` | `main_window.py:2195-2218`, `:2766` |
| **Xəta davranışı** | Validasiya xətası → `QMessageBox` və dayanma. Nyuton minimal Δt-də də yığılmasa → `result.converged=False`, **istisna atılmır**, o ana qədərki nəticə saxlanılır (`engine.py:199-204`). Gözlənilməz istisna → `worker.failed` → tam traceback jurnala, istifadəçiyə yalnız son sətir (`main_window.py:2186-2193`). |

### 6.2 Ssenari B — Quyu məlumatından geoloji model

| Addım | Nə baş verir | Fayl / funksiya |
|---|---|---|
| 1 | CSV oxunur; sütun adları sinonimlərlə tanınır, vahid şəkilçiləri ayrılır (`PERMX[mD]`) | `geology/well_data_io.py:90` `read_well_csv`, `:60` `_split_unit_suffix` |
| 2 | `dataset.validate()` — yararsız data variogramı korlamasın | `geology_service.py` `build()` başı |
| 3 | QC: dublikat koordinat siyasəti, kənar dəyər (robust miqyas) | `geology/data_quality.py:390` `run_quality_control` |
| 4 | Struktur səthlər (TOP/BOTTOM) **əvvəl** — qalan xassələr hüceyrə mərkəzinin Z-ini işlədir | `geology_service.py`, `domain/structural_grid.py` |
| 5 | Kateqorik sahələr (fasiya) SIS ilə | `geology/facies.py:237` `simulate_sis` |
| 6 | Kəsilməz xassələr: strategiya (`PERMX → log`, `PORO → logit`) → variogram fit → anizotropluq → qonşuluq → kriging (və ya SGS) → geri çevirmə | `geology/property_config.py:345`, `variogram.py:702`, `anisotropy.py`, `spatial_search.py:429`, `interpolation.py:332`, `sgs.py:158` |
| 7 | Çarpaz-doğrulama (LOO / k-fold) → RMSE, bias, etibar kateqoriyası | `geology/cross_validation.py:378` |
| 8 | Nəticə `GeologicalModel` + `InterpolationReport` | `geology_service.py:360` |
| **Xəta** | Yararsız data → `ValueError`; məlumatsız lay (lay-məlumatlı rejim) → `MISSING` + `report.blocking`, simulyasiyaya **buraxılmır** (`model_builder.py:61-63`) |

### 6.3 Ssenari C — Layihəni saxlamaq və açmaq

| Addım | Fayl / funksiya |
|---|---|
| Panellərin vəziyyəti toplanır (`panel_state.capture`) | `main_window.py:2539` `_collect_project_state`, `ui/panel_state.py:47` |
| `ProjectSerializer.save` — gzip + JSON, `FORMAT_VERSION = 3` | `serialization.py:367-377` |
| Açma: gzip → uğursuzsa adi JSON → versiya yoxlaması (gələcək versiya rədd edilir) | `serialization.py:379-393` |
| Panellər bərpa olunur, son nəticə göstərilir | `main_window.py:2576-2640` |
| «Son layihələr» (5 ədəd) `QSettings("IMEX2D","IMEX-2D")` — Windows reyestri `HKCU\Software\IMEX2D` | `main_window.py:174`, `application/session.py:31` |
| **Xəta** | `ProjectFileError` → `QMessageBox «Fayl açılmadı»`; qalan istisnalar «Gözlənilməz xəta» (`main_window.py:2583-2590`) |

### 6.4 Ssenari D — History matching (tarixi uyğunlaşdırma)

1. Müşahidə CSV oxunur — `history/observation_io.py:81`.
2. `standard_parameters(model)` — k çarpanı, Kv/Kh, SCAL, μo, kontakt, PV (`history/parameters.py:158`).
3. `HistoryMatchingService.run(method)` — Nelder-Mead / Powell / Differential Evolution (`scipy.optimize`), axtarış `[0,1]ⁿ` fəzasında (`optimizer.py:116-298`).
4. Hər qiymətləndirmə = **tam simulyasiya**; nəticə müşahidə vaxtlarına interpolyasiya → NRMSE (`history/mismatch.py:130`).
5. Yığılmayan model → `FAILURE_PENALTY = 1e6` (istisna atılmır), nəticələr keşlənir (`optimizer.py:40-45`).
6. Fon axını: `MatchingWorker` (`ui/worker.py:24`).

### 6.5 Ssenari E — GRDECL idxalı

`io/grdecl.py:132` `read_grdecl` — tokenləşdirmə, `n*value` genişlənməsi, `INCLUDE`
dəstəklənmir (xəbərdarlıq) → `io/grdecl_import.py:37` `GrdeclImporter` — COORD/ZCORN
varsa **həqiqi** corner-point həndəsə (`CornerPointGeometry.from_grdecl`, `grdecl_import.py:124-152`).

---

## 7. Texnologiyalar və kitabxanalar

Versiyalar `requirements.txt`-dən və venv-dəki `pip list`-dən (üst-üstə düşür).

### 7.1 Proqramlaşdırma dilləri

| Texnologiya | Versiya | Məqsəd | Fayllar | Əhəmiyyəti |
|---|---|---|---|---|
| Python | 3.14.7 (venv) | Bütün kod | hamısı | Əsas |
| Windows Batch | — | Başladıcılar | `run.bat`, `test.bat` | Köməkçi |
| HTML/CSS/JS | — | Təqdimat materialları (proqramın hissəsi deyil) | `docs/teqdimat/` | Sənəd |

### 7.2 Frontend (masaüstü UI)

| Kitabxana | Versiya | Məqsəd | Fayllar | Əhəmiyyəti |
|---|---|---|---|---|
| PyQt5 (+PyQt5-Qt5 5.15.2, PyQt5_sip 12.19.0) | 5.15.11 | Pəncərə, panellər, QThread, QSettings | `ui/*`, `app.py` | Kritik |
| matplotlib (Qt5Agg backend) | 3.11.1 | 2D qrafiklər, xəritələr, psevdo-3D, PDF | `rendering/renderers.py`, `rendering/volume.py`, `reporting/report.py` | Kritik |
| VTK | 9.7.0 | ResInsight tipli 3D səhnə, kəsik, PNG kadr | `rendering/vtk_volume.py`, `rendering/animation_export.py` | Yüksək (yoxdursa matplotlib-ə qayıdır — `logs/imex2d.log`) |
| Pillow | 12.3.0 | Kadrları `PIL.Image`-ə çevirmək, PNG/GIF yazmaq (funksiya daxilində idxal) | `rendering/animation_export.py:36,61` | Orta |

### 7.3 Backend / hesablama

| Kitabxana | Versiya | Məqsəd | Fayllar | Əhəmiyyəti |
|---|---|---|---|---|
| NumPy | 2.5.2 | Bütün massiv hesablamaları | hər yerdə | Kritik |
| SciPy | 1.18.1 | `scipy.sparse` (Jakobian, MPFA), `sparse.linalg` (CG, BiCGStab, spilu, spsolve), `optimize` (least_squares, minimize, differential_evolution), `spatial.cKDTree`, `stats.norm` | `simulation/implicit/linear.py`, `cpr.py`, `linear_solver.py`, `geology/variogram.py`, `spatial_search.py`, `gaussian_transform.py`, `history/optimizer.py`, `discretization/` | Kritik |
| Backend framework (Flask, FastAPI …) | — | **Yoxdur** | — | — |

### 7.4 Database, API, Auth, ML/AI

| Kateqoriya | Vəziyyət |
|---|---|
| Database | **Yoxdur.** Davamlı saxlama — `.imx` (gzip+JSON) faylı və `QSettings` |
| API (REST/gRPC) | **Yoxdur** |
| Authentication/authorization | **Yoxdur** (tək istifadəçili masaüstü tətbiq) |
| Machine learning / AI / LLM | **Yoxdur.** «Model» sözü layihədə fiziki və statistik modelləri bildirir (§8) |

### 7.5 Data processing

NumPy/SciPy (yuxarıda). Pandas **işlədilmir** [FAKT] — CSV-lər standart `csv` modulu ilə
oxunur/yazılır: `geology/well_data_io.py:34`, `io/scal_io.py:41`, `io/fault_io.py:31`,
`history/observation_io.py:19`, `reporting/results_export.py:21`, `ui/main_window.py:9`.

### 7.6 Logging və monitoring

| Texnologiya | Məqsəd | Fayl |
|---|---|---|
| `logging` (stdlib) + `RotatingFileHandler` (2 MB × 3) | Konsol + `logs/imex2d.log` | `logging_setup.py:39-75` |
| `QtLogHandler` | Eyni axın «Jurnal» tabında | `ui/main_window.py:90` |
| Monitoring (Prometheus, Sentry …) | **Yoxdur** | — |

### 7.7 Testing

| Kitabxana | Versiya | Məqsəd |
|---|---|---|
| pytest | 9.1.1 | Əsas test çərçivəsi (`pytest.ini`) |
| pypdf | 6.16.2 | PDF hesabatın məzmununu yoxlamaq (`requirements-dev.txt`) |
| `run_tests.py` | — | pytest olmadan ehtiyat qaçırıcı |
| vulture, pyflakes | 2.16 / 3.4.0 | venv-də quraşdırılıb, `requirements-dev.txt`-də **yoxdur** — ölü kod/lint üçün əl ilə işlədildiyi güman edilir **[⏳ ƏMİN DEYİL]** |
| coverage / pytest-cov | — | **Quraşdırılmayıb** |

### 7.8 Build, deployment, konfiqurasiya, asılılıq idarəsi

| Mövzu | Vəziyyət |
|---|---|
| `requirements.txt` | Var — 16 paket, **dəqiq versiyalarla** (`==`) |
| `requirements-dev.txt` | Var — `-r requirements.txt` + pytest + pypdf |
| `pyproject.toml`, `setup.py`, `Pipfile`, `poetry.lock` | **Yoxdur** — paket quraşdırılmır, `app.py` kök qovluqdan işlədilir |
| `Dockerfile`, `docker-compose.yml` | **Yoxdur** (Qt masaüstü tətbiqi üçün adi haldır) |
| `.env`, `.env.example` | **Yoxdur**; yalnız 4 mühit dəyişəni (§11.3) |
| CI (GitHub Actions və s.) | **Yoxdur** (`.github/` yoxdur) |
| Paketləmə (PyInstaller, MSI) | **Yoxdur** — istifadəçi venv yaradıb `run.bat` işlədir |

### 7.9 Xarici servislər və inteqrasiyalar

| İnteqrasiya | Necə | Fayl |
|---|---|---|
| OPM Flow nəticələri (`.EGRID`, `.UNRST`) | İstəyə bağlı `resdata` kitabxanası ilə — **venv-də yoxdur, `requirements.txt`-də yoxdur** | `io/opm_import.py:104` `_require_resdata` |
| Eclipse summary (`.SMSPEC`/`.UNSMRY`) | Öz binar oxuyucusu, `resdata`-sız | `io/eclipse_summary.py` |
| Eclipse deck | Öz oxuyucu/yazıcı | `io/grdecl.py`, `io/pvt_io.py`, `io/scal_io.py`, `io/eclipse_export.py` |
| graphify | Yalnız inkişaf aləti (AI köməkçinin hook-u), tətbiqin hissəsi deyil | `.claude/settings.json`, `.codex/hooks.json` |

---

## 8. İstifadə olunan modellər

### 8.1 Əvvəlcə vacib aydınlıq

**Layihədə süni intellekt, LLM, embedding, vektor bazası, RAG, neyron şəbəkə,
classification, clustering, computer vision, NLP və ya recommendation sistemi
YOXDUR [FAKT]** (asılılıqlarda heç bir ML kitabxanası yoxdur — §7; şəbəkə kodu yoxdur — §0.1).

Layihədəki «modellər» iki növdür:
- **Fiziki/empirik modellər** — neft-qaz mühəndisliyinin düsturları (Darcy, PVT korrelyasiyaları, Corey, Peaceman …).
- **Statistik modellər** — geostatistika (variogram, Kriging, SGS, SIS) və optimallaşdırma (Nelder-Mead, DE). Bunlar klassik statistikadır, «machine learning» deyil; lakin **regressiya** (variogram fit — `least_squares`) və **optimallaşdırma** anlayışları burada həqiqətən işlədilir.

Hər modelin nəzəri izahı və düsturu: `docs/nezeri_esaslar.md`. Aşağıda tələb olunan formatda xülasə.

### 8.2 Fiziki modellər

| Model | Növ | Rol | Giriş → Çıxış | Harada çağırılır | Parametrlər | Niyə seçilib / alternativ | Məhdudiyyət, səhv ehtimalı | Performans | Fayl |
|---|---|---|---|---|---|---|---|---|---|
| **Black-oil (3 komponent)** | PDE sistemi | Bütün simulyasiyanın əsası | p, Sw, Sg/Rs → qalıq R | Hər Nyuton iterasiyası | — | Sənaye standartı; alt. — kompozisiya (EOS) modeli (yoxdur) | Rv = 0, izotermik, tək məsaməlilik (`nezeri_esaslar.md` §10) | Əsas xərc | `implicit/residual.py`, `three_phase_residual.py` |
| **Darcy qanunu + TPFA** | Diskretizasiya | Hüceyrələr arası axın | K, həndəsə → T_ij | Mühərrik qurulanda | harmonik orta | Sadə, sürətli; alt. MPFA-O | K-ortoqonal olmayan gridlərdə səhv verir | O(N) | `simulation/discretization.py:85` |
| **MPFA-O** | Çoxnöqtəli diskretizasiya | Anizotrop/qeyri-ortoqonal gridlərdə dəqiq axın | tam K tenzoru + təpələr → T_cell | `flux_scheme="MPFA-O"` | `NEUMANN_ZERO` sərhəd | Layihənin əsas fərqi (`README.md`); alt. MPFA-L | Yalnız 2 fazalı FIM; fay/ACTNUM ilə rədd edilir | Stensil ≤18, lokal ≤12×12 | `discretization/*.py` |
| **Standing (1947)** | Empirik korrelyasiya | Rs(p), Pb | p, API, γg, T → Rs | PVT cədvəli qurulanda | — | Klassik; alt. Glaso, Vasquez-Beggs Rs | ±10-15 % tipik xəta **[RƏY — ədəbiyyat, kodda yoxdur]** | Bir dəfə | `pvt/correlations.py:36-53` |
| **Vazquez-Beggs (1980)** | Empirik | Bo, co | Rs, API, γg, T → Bo | eyni | API ≤ 30 / > 30 əmsal dəstləri | — | — | Bir dəfə | `correlations.py:56-78` |
| **Beggs-Robinson (1975)** | Empirik | μod, μo | API, T, Rs → μ | eyni | — | — | — | Bir dəfə | `correlations.py:81-96` |
| **Sutton + Beggs-Brill Z + Lee-Gonzalez-Eakin** | Empirik | Qaz: Tpc/Ppc, Z, Bg, μg | γg, p, T | eyni | — | Beggs-Brill iterasiyasız (DAK-dan sürətli) — `nezeri_esaslar.md` §2.1 | Z approksimasiyası | Bir dəfə | `correlations.py:178-256` |
| **Deck PVT (PVTO/PVDG/PVTW)** | Cədvəl | Laboratoriya məlumatı | cədvəl → np.interp | Mühərrikdə | Q-32, Q-34 qaydaları | OPM Flow ilə eyni | TB-3 (korrelyasiya yolunda) | `np.interp` | `io/pvt_io.py`, `pvt/black_oil.py` |
| **Corey** | Analitik relperm | krw, kro, törəmələr | Sw → kr | Hər iterasiya | nw, no, krw_end, kro_end, Swc, Sor | Sadə; alt. SWOF cədvəli | Laboratoriya əyrisini təxmini verir | Vektor | `scal_adapter.py:20` |
| **SWOF/SGOF cədvəlləri** | Cədvəl relperm (+Pc) | Region əsaslı SCAL | Sw/Sg → kr, Pc | Hər iterasiya | SATNUM | — | — | Vektor | `scal_tables_provider.py`, `domain/scal_tables.py` |
| **Stone II (1973)** | 3 fazalı relperm | kro(Sw, Sg) | krow, krog → kro | 3 fazalı | — | Sənayedə geniş; alt. OPM defolt modeli | Su hərəkətli modellərdə OPM-dən 50 %-ə qədər fərq (`ROADMAP.md` backlog) | Vektor | `simulation/stone_relperm.py:37` |
| **Brooks-Corey Pc** | Kapilyar | Pcow, Pcog | S → Pc | Hər iterasiya | entry_pressure, λ | — | — | Vektor | `simulation/capillary.py` |
| **Peaceman (1978/1983)** | Quyu indeksi | Quyu–hüceyrə bağlantısı | K, həndəsə, rw, skin → WI | Bir dəfə | `re = 0.28·√…` | Sənaye standartı | Tək hüceyrə ətrafı; `ln(re/rw)` 1.01-də kəsilir (`well_model.py:186`) | Bir dəfə | `simulation/well_model.py:53` |
| **Darcy-Weisbach + Chen (1979) + no-slip holdup** | Lülə hidravlikası | BHP ↔ THP | debit, həndəsə → Δp | Post-proses və THP idarəsi | seqment sayı | Sadə; alt. Hagedorn-Brown, Beggs-Brill çoxfazalı | No-slip — sürüşmə yoxdur; sürətlənmə həddi yoxdur | Seqment başına | `simulation/wellbore/` |
| **Hidrostatik tarazlıq** | İlkin şərt | p(z), Sw(z) | datum, OWC/GOC, Pc | Bir dəfə | — | — | — | Bir dəfə | `initialization/equilibrium.py:61` |
| **Buckley-Leverett + Welge** | Analitik həll | Doğrulama etalonu | SCAL, μ → cəbhə profili | Testlər, «Validasiya» tabı | — | — | 1D, cazibəsiz | — | `simulation/analytical.py:36` |

### 8.3 Statistik modellər

| Model | Növ | Rol | Giriş → Çıxış | Parametrlər | Məhdudiyyət | Fayl |
|---|---|---|---|---|---|---|
| **Variogram (sferik/eksponensial/qauss)** | Regressiya (çəkili ən kiçik kvadratlar) | Məkan korrelyasiyası | nöqtələr, dəyərlər → nugget, sill, range | model="auto" → ən yaxşısı | Az quyuda qeyri-sabit fit | `geology/variogram.py:702, 741` |
| **Adi Kriging (BLUE)** | Məkan interpolyasiyası | Hüceyrə dəyəri + variansı | quyular → Ẑ, σ² | qonşu sayı, radius, anizotropluq | Hamar sahə verir; həll yolu açıq qeyd olunur (`direct/jitter/lstsq/idw_fallback`) | `geology/interpolation.py:332` |
| **IDW, ən yaxın qonşu** | Deterministik | Sadə alternativ | — | güc p | «Öküz gözü» effekti | `interpolation.py:119,131` |
| **SGS** | Stoxastik simulyasiya | Realistik dəyişkənlik, çoxlu realizasiya | quyular + seed → sahə | seed, realizasiya sayı | Nəticə variogramın doğruluğundan asılıdır | `geology/sgs.py:158`, `sgs_ensemble.py` |
| **SIS (indikator kriging)** | Kateqorik stoxastik | Fasiya xəritəsi | fasiya kodları → sahə | nisbətlər, variogram | — | `geology/facies.py:237` |
| **Normal-score / log / logit çevirmələri** | Verilən çevirməsi | Paylanmanı normallaşdırmaq | — | offset, hədlər | — | `geology/transforms.py`, `gaussian_transform.py` |
| **Çarpaz-doğrulama (LOO, k-fold, məkan blokları)** | Model qiymətləndirmə | RMSE, bias, model seçimi | — | fold sayı | — | `geology/cross_validation.py` |
| **Nelder-Mead / Powell / Differential Evolution** | Törəməsiz optimallaşdırma | History matching | parametr vektoru → NRMSE | büdcə | Hər qiymətləndirmə tam simulyasiya — çox bahalı | `history/optimizer.py` |
| **Tornado (OAT) və yerli elastiklik** | Həssaslıq | Parametr təsiri | — | addım payı | Qarşılıqlı təsirləri görmür | `history/sensitivity.py` |

**Qiymət təsiri:** heç bir model pullu API çağırmır — xərc yalnız CPU vaxtıdır.

---

## 9. Nəzəri əsaslar

Yalnız layihədə **faktiki işlədilən** nəzəriyyələr. Hər birinin sadə izahı,
sonra harada işlədiyi. Riyazi detallar: `docs/nezeri_esaslar.md`.

### 9.1 Fizika və rəqəmsal üsullar

| Anlayış | Sadə dillə | Layihədə harada |
|---|---|---|
| **Darcy qanunu** | Su borudan necə axırsa, məsaməli daşda da maye təzyiq fərqinə görə axır; daş nə qədər «keçirici», maye nə qədər «duru»dursa, axın o qədər güclüdür. | `residual.py` — `compute_flux(ΔΦ)·λ/B` |
| **Kütlənin saxlanması** | Hər hüceyrə bir «qab»dır: girən − çıxan = içində artan. Bunun pozulma ölçüsü «qalıq» (residual) adlanır. | `residual.py`; `newton.py:150-165` (MB meyarı) |
| **Sonlu həcm metodu** | Yatağı minlərlə kiçik kubiklərə bölüb hər biri üçün yuxarıdakı balansı yazmaq. | `domain/grid.py`, `simulation/discretization.py` |
| **Upstream çəkilənmə** | Maye hansı tərəfdən gəlirsə, onun xassəsini götür — əks halda hesab «dalğalanır». | `residual.py`, `impes_engine.py` |
| **IMPES** | Əvvəl təzyiqi bir dəfəyə hesabla, sonra doyumluluğu köhnə dəyərlərlə «irəli at». Sürətli, amma addım kiçik olmalıdır (CFL). | `simulation/impes_engine.py:94` |
| **Tam implicit + Nyuton-Rafson** | Bütün naməlumları eyni anda tap; tənlik qeyri-xəttidir, ona görə «təxmin et → səhvi ölç → düzəlt» dövrəsi (Nyuton) aparılır. Addım böyük ola bilər. | `implicit/engine.py`, `newton.py` |
| **Jakobian** | «Hər naməlumu bir az dəyişsəm, hər tənlik nə qədər dəyişər?» cədvəli. Nyuton bu cədvəllə düzəlişi hesablayır. Burada analitik (düsturla) qurulur. | `implicit/jacobian.py`, `three_phase_residual.py:607-1401` |
| **Seyrək xətti cəbr, ILU, CPR, BiCGStab** | Milyonlarla sıfırlı nəhəng tənlik sistemi; yalnız sıfır olmayanlar saxlanılır. ILU/CPR — həlli sürətləndirən «ipucu» (ön-şərtçi). | `implicit/linear.py`, `cpr.py`, `linear_solver.py` |
| **Adaptiv zaman addımı** | Hesab asan gedirsə addımı böyüt, çətin gedirsə kiçilt, alınmırsa yarıya böl və təkrarla. | `implicit/time_stepping.py` |
| **Dəyişən keçid (variable switching)** | Neftdə qaz ya tam həll olub (Rs naməlumdur), ya da ayrılıb (Sg naməlumdur) — hüceyrənin vəziyyətinə görə naməlum dəyişdirilir. | `three_phase_state.py`, `three_phase_newton.py` |

### 9.2 Geostatistika

| Anlayış | Sadə dillə | Harada |
|---|---|---|
| **Variogram** | «Bir-birindən h metr aralı iki nöqtə orta hesabla nə qədər fərqlənir?» qrafiki. | `geology/variogram.py` |
| **Kriging** | Qonşu quyuların çəkili ortası, amma çəkilər variogramdan və quyuların bir-birinə yaxınlığından «ağıllı» seçilir; üstəlik «nə qədər əminik» (varians) verir. | `geology/interpolation.py` |
| **SGS / SIS** | Kriging-in hamar nəticəsi əvəzinə təsadüfi, amma quyulara hörmət edən «mümkün reallıqlar» yaratmaq. | `geology/sgs.py`, `facies.py` |
| **Anizotropluq** | Laylarda korrelyasiya uzununa və eninə fərqlidir; koordinatları dartıb adi məsafə işlətmək. | `geology/anisotropy.py` |
| **Çarpaz-doğrulama** | Bir quyunu gizlət, qalanlarından onu təxmin et, səhvi ölç — metodun etibarlılığı. | `geology/cross_validation.py` |

### 9.3 Optimallaşdırma və regressiya

- **Ən kiçik kvadratlar (regressiya):** variogram modelini ölçülmüş nöqtələrə
  «ən yaxşı» uyğunlaşdırmaq — `scipy.optimize.least_squares` (`variogram.py:593`).
- **Törəməsiz optimallaşdırma:** history matching-də simulyatorun törəməsi yoxdur,
  ona görə Nelder-Mead (simpleks), Powell və DE (təkamül) işlədilir (`optimizer.py:129-160`).

### 9.4 Proqram mühəndisliyi prinsipləri

| Prinsip | Sadə dillə | Layihədə nümunə |
|---|---|---|
| **OOP** | Məlumat və onunla işləyən funksiyalar bir «sinif»də. | Bütün `domain/`, mühərriklər |
| **Layered architecture** | Hər qat yalnız aşağıdakını tanıyır. | §4.3 |
| **Dependency Inversion (SOLID-in D-si)** | Mühərrik «hansısa relperm verən»dən asılıdır, konkret Corey-dən yox. | `interfaces/providers.py`; `ARCHITECTURE.md` §5.5 |
| **Open/Closed (SOLID-in O-su)** | Yeni PVT əlavə etmək üçün mühərriki açmaq lazım deyil. | `ARCHITECTURE.md` §6 |
| **Single Responsibility (S)** | Hər sinfin bir işi. **Qismən pozulur** — `MainWindow`. | §14 P-01 |
| **Liskov (L)** | `CornerPointGeometry` `CellGeometry`-nin yerinə keçir, çağıranlar dəyişmir. | `domain/corner_point_geometry.py:354` |
| **Interface Segregation (I)** | Kiçik, məqsədli ABC-lər (`ILinearSolver` 1 metod). | `interfaces/services.py` |
| **Composition root** | Bütün «new» bir yerdə. | `app.py` |
| **Strategy pattern** | Dəyişdirilə bilən alqoritm obyekti: diskretizasiya, interpolyator, holdup korrelyasiyası. | `IFluxDiscretization`, `IPropertyInterpolator`, `IHoldupCorrelation` |
| **Adapter pattern** | Mövcud düsturu interfeysə bağlayan nazik təbəqə. | `CoreyRelativePermeabilityAdapter` |
| **Factory** | `engine_factory` parametri — hansı mühərrikin qurulacağını çağıran seçir. | `simulation_service.py:60` |
| **Observer (siqnal/slot)** | Fon axını hadisə göndərir, UI dinləyir. | `ui/worker.py` |
| **Dataclass-lar (dəyər obyektləri)** | Konfiqurasiya və nəticə — sadə məlumat qabları. | `config.py`, `results.py` |
| **Funksional proqramlaşdırma** | Saf funksiyalar (yan təsirsiz): korrelyasiyalar, `daily.py`, `session.py`. Paradiqma kimi deyil, üslub kimi. | `pvt/correlations.py`, `reporting/daily.py` |

### 9.5 Concurrency, xəta idarəsi, serializasiya

- **Concurrency:** yalnız `QThread` — hər uzun iş (simulyasiya, uyğunlaşdırma, həssaslıq) bir fon axınında; hesablamanın özü tək axınlıdır. `async`/multiprocessing **yoxdur**.
- **Xəta idarəsi:** domain yoxlamaları istisna atmır, **siyahı qaytarır** (`validate() -> list`); servis onları `ModelValidationError`-a çevirir; ədədi uğursuzluq istisna yox, `converged=False` ilə bildirilir («təhlükəsiz uğursuzluq»).
- **Serializasiya:** gzip + JSON, versiyalı format (`FORMAT_VERSION = 3`), köhnə fayllar oxunur, gələcək versiya rədd edilir.

### 9.6 Layihədə İŞLƏDİLMƏYƏN anlayışlar

Tələb siyahısındakı bu mövzular layihəyə **aid deyil** və fakt kimi təqdim edilmir:
deep learning, neural networks, transformers, attention, embeddings, vector similarity,
RAG, prompt engineering, fine-tuning, classification, clustering, NLP, REST API,
authentication/authorization, database normalization, event-driven arxitektura
(UI siqnallarından başqa), asinxron proqramlaşdırma, caching (history matching-in
yaddaş keşindən — `optimizer.py:21` — başqa).

---

## 10. Əsas fayllar və funksiyalar

### 10.1 `app.py` — composition root

- **Məqsəd:** loglamanı qurmaq, `QApplication` yaratmaq, asılılıqları bağlamaq, pəncərəni açmaq.
- **Funksiya:** `main()` (`app.py:36-57`).
- **Çağırır:** `logging_setup.configure`, `ModelAwareSimulationService`, `CoreyRelativePermeabilityAdapter`, `ScipyCgIluSolver`, `MainWindow`, `SyntheticGeologicalModelBuilder`, `ReservoirModelBuilder`.
- **Kim çağırır:** `run.bat`, `python app.py`.
- **Dəyişməsinin təsiri:** bütün tətbiqin başlanğıc konfiqurasiyası.

### 10.2 `imex2d/application/simulation_service.py`

| Sinif/funksiya | Giriş | Çıxış | Alqoritm | Xəta halları | Yan təsir |
|---|---|---|---|---|---|
| `SimulationService.create_engine(model, config)` `:75` | `ReservoirModel`, `SimulationConfig` | mühərrik obyekti | validate → uyğunsuzluq yoxlamaları → `engine_factory(...)` | `ModelValidationError` (validasiya, MPFA/THP/limit uyğunsuzluğu, `NotImplementedError`-ın çevrilməsi) | həlledicini `reset()` edir |
| `SimulationService.run(...)` `:205` | model, config, reporter | `SimulationResult` | `create_engine().run()` + THP post-prosesi | yuxarıdakılar | — |
| `ModelAwareSimulationService.create_engine` `:262` | eyni | eyni | provider-ləri modeldən qurur; qaz varsa 3 fazalı mühərrik | eyni + «MPFA-O 3 fazada yoxdur» | **`self.*_provider` sahələrini dəyişir** (servis vəziyyətlidir) |

- **Kim çağırır:** `ui/main_window.py` (run, sensitivity, matching), `history/optimizer.py`, `history/sensitivity.py`, `benchmarks/spe1.py`, `tools/`, testlər.
- **Dəyişməsinin təsiri:** hər simulyasiya yolu — çox yüksək risk.

### 10.3 `imex2d/simulation/implicit/engine.py` — `FullyImplicitEngine`

- **`__init__` (`:56-122`)** — diskretizasiya, Peaceman bağlantıları, `ResidualAssembler`, `JacobianAssembler`, `NewtonSolver`, `AdaptiveTimeStepper`, 3 quyu nəzarətçisi, ilkin vəziyyət. **Diqqət:** `linear_solver` parametri qəbul edilir, amma **işlədilmir** — Nyuton həmişə `NewtonLinearSolver()` defoltunu alır (`:84-86`); `growth_factor` konfiqurasiyadakından **+0.35** böyük götürülür (`:145`).
- **`run(reporter)` (`:169-290`)** — zaman dövrəsi. Giriş: reporter; çıxış: `SimulationResult`. Xəta: istisna atmır, `converged=False` + mesaj (`:199-204`); istifadəçi dayandırsa mesaj yazılır (`:269-271`). Performans: hər addımda 1+ Nyuton həlli, THP/səth/BHP dövrələri əlavə həllər edə bilər.
- **Kim çağırır:** `SimulationService.create_engine` (engine_factory), `MainWindow` (sinif seçimi).

### 10.4 `imex2d/simulation/implicit/newton.py` — `NewtonSolver.solve`

- **Giriş:** əvvəlki vəziyyət, `dt`, başlanğıc təxmin. **Çıxış:** `NewtonResult` (status, vəziyyət, iterasiya, CNV tarixçəsi, MB, debitlər).
- **Alqoritm (`:181-248`):** qalıq → CNV/MB → yığılma yoxlaması → divergensiya yoxlaması (CNV > 10×ilk) → Jakobian → ACTNUM daraltması → xətti həll → Appleyard kəsməsi → xətti axtarış (`:260-294`, maks. 10 dəfə yarıya bölmə).
- **Parametrlər (`NewtonConfig` `:50-100`):** `max_iterations=20`, `cnv_tolerance=1e-3`, `material_balance_tolerance=1e-7`, `max_pressure_change=50 bar`, `max_saturation_change=0.2`.
- **Xəta halları:** `LINEAR_SOLVER_FAILED` (FloatingPointError/ValueError/NaN), `DIVERGED`, `MAX_ITERATIONS`.
- **Performans:** hər iterasiyada 1 Jakobian + 1 xətti həll + 1-11 qalıq qiymətləndirməsi (xətti axtarış).

### 10.5 `imex2d/simulation/implicit/linear.py` — `NewtonLinearSolver.solve`

- `n ≤ 20 000` → birbaşa `spsolve`; böyükdürsə BiCGStab + ön-şərtçi (`n ≥ 40 000` → CPR, yoxsa ILU); BiCGStab yığılmasa birbaşa həllə qayıdır; o da alınmasa `FloatingPointError` (`:59-90`).

### 10.6 `imex2d/application/geology_service.py` — `WellBasedGeologicalModelBuilder.build`

- **Giriş:** `WellDataset`, `GeologicalGridSpec`, anizotropluq nisbətləri, fasiya/SGS/lay konfiqurasiyaları. **Çıxış:** `GeologicalModel` + `InterpolationReport` (`:441`).
- 2 369 sətirlik faylın mərkəzi; geologiya boru xəttinin bütün ardıcıllığı (§6.2) buradadır.
- **Dəyişməsinin təsiri:** hər quyu-əsaslı model, çarpaz-doğrulama, keyfiyyət hesabatı.

### 10.7 `imex2d/application/serialization.py` — `ProjectSerializer`

- `save(project, path, include_snapshots)` (`:367`) — gzip+JSON, **birbaşa hədəf fayla**.
- `load(path)` (`:379`) — gzip və ya adi JSON; versiya yoxlaması.
- 855 sətir — hər domain obyekti üçün `*_to_dict` / `*_from_dict`. Yeni sahə əlavə edəndə həm burada, həm testdə dəyişiklik lazımdır.

### 10.8 `imex2d/ui/main_window.py` — `MainWindow`

- 3 155 sətir, 121 metod, 13 tab (`version.py:38-43`).
- Əsas metodlar: `rebuild_model` `:1704`, `run_simulation` `:2115`, `_on_finished` `:2195`, `save_project` `:2525`, `open_project` `:2576`, `_write_recovery` `:2766`, `closeEvent` `:2830`.
- **Test edilməsi çətindir** — pytest-də qurulmur (§13).

### 10.9 Digər vacib fayllar (qısa)

| Fayl | Nə edir |
|---|---|
| `domain/reservoir_model.py` | «Yeganə həqiqət mənbəyi» — simulyasiyaya hazır model, `validate()`, `diagnose()` |
| `domain/grid.py` | `CartesianGrid`, `Connections`, `ActiveMap` (ACTNUM konvensiyası) |
| `domain/geometry.py`, `corner_point_geometry.py` | Hüceyrə həndəsəsi (Kartezian və corner-point) |
| `simulation/implicit/residual.py` | 2 fazalı qalıq — hər Nyuton iterasiyasının əsası |
| `simulation/implicit/three_phase_residual.py` | 3 fazalı qalıq və Jakobian (1 401 sətir, ən böyük fizika faylı) |
| `simulation/pvt/black_oil.py` | PVT provider — cədvəl interpolyasiyası və törəmələr |
| `simulation/well_constraints.py` | RATE payı, BHP limiti, səth debiti nəzarətçiləri |
| `simulation/wellbore/thp_control.py` | THP → BHP, nodal analiz |
| `discretization/mpfa_o.py` | MPFA-O əmsalları və `MPFAODiscretization` |
| `geology/interpolation.py` | `OrdinaryKriging` (god node, 203 əlaqə) |

---

## 11. Konfiqurasiya və idarəetmə

### 11.1 Proqram necə başladılır

```bat
:: Windows — run.bat venv-i bu üç yerdə axtarır: .venv\, ..\venv\, venv\
run.bat
```
```bash
# və ya əl ilə
..\venv\Scripts\activate
python app.py
```
İlk quraşdırma: `python -m venv ..\venv` → `pip install -r requirements-dev.txt`
(`README.md` §3 bu hissəni hələ ⏳ kimi saxlayır; `ARCHITECTURE.md` §7-də qısa təlimat var).

### 11.2 Əmrlər

| Əmr | Nə edir |
|---|---|
| `python app.py` | Tətbiqi açır |
| `pytest` | Bütün testlər (`pytest.ini`: `-v --tb=short`) |
| `pytest -m "not performance"` | Uzun benchmark testləri olmadan |
| `IMEX_SKIP_SLOW=1 pytest` | Yavaş testləri keç (`tests/helpers.py:26`) |
| `python run_tests.py [-q]` | pytest olmadan |
| `test.bat` | `..\venv` aktivləşdirir, `pytest` |
| `python tools/benchmark.py [--profile 41]` | Performans ölçməsi |
| `python tools/spe1_compare.py <DATA> <etalon>` | SPE1CASE2 — OPM Flow müqayisəsi |
| `python tools/presentation_demo.py --out <qovluq>` | Təqdimat rəqəmlərini yenidən yaratmaq |
| `python tools/module_graph.py` | Modul asılılıq qrafiki |
| `graphify update .` | Bilik qrafikini yeniləmək (inkişaf aləti) |

### 11.3 Mühit dəyişənləri

| Dəyişən | Təsiri | Fayl |
|---|---|---|
| `IMEX2D_DATA_DIR` | Bərpa faylının qovluğunu dəyişir (test/xüsusi quraşdırma) | `application/session.py:23,50` |
| `LOCALAPPDATA` | Defolt məlumat qovluğu (`%LOCALAPPDATA%\IMEX2D`) | `session.py:53` |
| `IMEX_SKIP_SLOW` | Yavaş testləri keçir | `tests/helpers.py:26` |
| `QT_QPA_PLATFORM` | Testlərdə `offscreen` (başsız Qt) | `tests/conftest.py:13` |

API key, parol, token **yoxdur** — idarə olunmalı secret yoxdur.

### 11.4 Konfiqurasiya harada saxlanılır

| Nə | Harada |
|---|---|
| Ədədi sabitlər (tolerans, Δt, CFL, snapshot sayı, axın sxemi) | `application/config.py` dataclass-ları — UI-dən `NumericalPanel` doldurur |
| Nyuton/xətti həll sabitləri | `NewtonConfig` (`newton.py:50`), `NewtonLinearSolverConfig` (`linear.py:27`), `AdaptiveTimeStepConfig` (`time_stepping.py:31`) — **UI-dən dəyişdirilmir** |
| Layihə parametrləri | `.imx` faylı (`Project.ui_state` daxil) |
| Son layihələr, bərpa mənbəyi | `QSettings` → `HKCU\Software\IMEX2D` |
| Log | `logs/imex2d.log` (tətbiq qovluğunda), 2 MB × 3 rotasiya |

### 11.5 Mühitlər (development / testing / production)

**Ayrıca mühit anlayışı yoxdur [FAKT].** Fərqlər yalnız bunlardır: testlərdə
`QT_QPA_PLATFORM=offscreen` və `IMEX_SKIP_SLOW`; real istifadədə `run.bat`.
«Production» — istifadəçinin öz kompüterindəki venv.

### 11.6 Log-lara baxmaq və xəta araşdırmaq

1. Proqramda **«Jurnal»** tabı — eyni log axını (`QtLogHandler`).
2. Fayl: `logs/imex2d.log` (UTF-8). Hər qaçışın başında `Mühərrik qurulur: … | TPFA/MPFA-O` sətri var (`simulation_service.py:84-86`) — hansı sxemin işlədiyini göstərir.
3. Hesablama xətası → tam traceback `LOG.error` ilə yazılır (`main_window.py:2190`).
4. Nyuton problemləri `DEBUG` səviyyəsindədir (`newton.py:210`, `linear.py:84`) — `app.py:39`-da səviyyə `INFO`-dur, ona görə defoltda görünmür. Ətraflı diaqnoz üçün səviyyəni `logging.DEBUG` etmək lazımdır (kod dəyişikliyi).
5. Yumşaq qəbul edilmiş addımlar `WARNING` kimi yazılır və nəticə mesajında sayılır (`engine.py:277-281`).

### 11.7 Yeni funksiya əlavə etmək — hansı fayllar

| Nə əlavə edilir | Dəyişən fayllar |
|---|---|
| Yeni PVT/SCAL/kapilyar modeli | yeni provider sinfi (`simulation/`) → `interfaces/providers.py` müqaviləsinə uyğun → `ModelAwareSimulationService`-də seçim qaydası → (lazımdırsa) `domain/` parametr sinfi + `serialization.py` + UI paneli |
| Yeni quyu idarə rejimi | `domain/wells.py` (`ControlMode`), `simulation/well_model.py`, qalıq/Jakobian və ya `well_constraints.py`, hər iki FIM mühərriki (kod təkrarı — P-08), IMPES-də rədd, `serialization.py`, `ui/panels.py` `WellPanel` |
| Yeni nəticə sırası | `simulation/results.py` → 3 mühərrikin hamısı → `serialization.py` → `results_export.py` → `renderers.py` → `reporting/daily.py` |
| Yeni fayl formatı | `io/` yeni modul → `main_window.py` menyu əmri |
| Yeni tab | `main_window.py` + `version.py` `EXPECTED_TABS` (açılışda yoxlanılır) |

### 11.8 Deployment

**Formal deployment prosesi yoxdur [FAKT].** Repo klonlanır, venv yaradılır,
`pip install -r requirements.txt`, `run.bat`. Docker, installer, avtomatik
yeniləmə, CI/CD **yoxdur**. `version.py` versiyanı pəncərə başlığında göstərir.

---

## 12. Təhlükəsizlik analizi

**Kontekst:** oflayn, tək istifadəçili masaüstü tətbiq. Şəbəkə, verilənlər
bazası, autentifikasiya, veb interfeys yoxdur. Buna görə veb tipli risklərin
çoxu (SQL injection, XSS, CSRF, rate limiting) **tətbiq olunmur**.

| Risk | Harada | Ciddilik | İstismar ssenarisi | Həll |
|---|---|---|---|---|
| API key / secret sızması | — | **Yoxdur** | Secret yoxdur | — |
| Authentication / authorization | — | **Tətbiq olunmur** | Tək istifadəçi; OS istifadəçi hüquqları | — |
| SQL injection | — | **Tətbiq olunmur** | DB yoxdur | — |
| Command injection | — | **Yoxdur** | `subprocess`/`os.system`/`eval`/`exec` yoxdur (§0.1 əmr 8) | — |
| XSS / CSRF | — | **Tətbiq olunmur** | Veb yoxdur. `docs/teqdimat/*.html` statik sənəddir | — |
| **Etibarsız faylın deserializasiyası** | `serialization.py:379-393` | **Aşağı** | `.imx` JSON-dur, `pickle` deyil — kod icrası mümkün deyil ✅. Lakin gzip «bombası» (kiçik fayl → nəhəng JSON) yaddaşı doldura bilər | Oxunan ölçüyə hədd (`gzip` axınından məhdud oxu) |
| **GRDECL `n*value` genişlənməsi** | `io/grdecl.py:110-117` | **Aşağı** | Səhv və ya zərərli faylda `999999999*0.1` → milyardlarla element → yaddaş tükənir, proqram donur | `count`-u `nx·ny·nz` (və ya ağlabatan hədd) ilə məhdudlaşdırmaq |
| **Input validation** | `domain/validation.py`, `config.py:76`, `reservoir_model.py:179` | **Yaxşı** ✅ | Porozite, doyumluluq, keçiricilik, təzyiq, debit, grid ölçüsü, tenzorun müsbət-müəyyənliyi yoxlanılır | — |
| **Atomik olmayan yazma (data bütövlüyü)** | `serialization.py:375` | **Orta** | Yazma zamanı elektrik kəsilməsi/disk dolması → istifadəçinin `.imx` faylı yarımçıq qalır, köhnə nüsxə yoxdur | Müvəqqəti fayla yaz → `os.replace` |
| **Bərpa faylında məlumat** | `%LOCALAPPDATA%\IMEX2D\berpa.imx` (`session.py:58`) | **Aşağı** | Layihənin tam surəti (yataq məlumatı — kommersiya baxımından həssas ola bilər) şifrəsiz saxlanılır; düzgün bağlanmada silinir | Sahibkarın qərarı — adətən qəbul edilə bilən |
| **Log-da həssas məlumat** | `logs/imex2d.log` | **Aşağı** | Fayl yolları (istifadəçi adı daxil) və model xülasəsi yazılır; parol/şəxsi məlumat yoxdur | — |
| **Fayl upload** | — | **Tətbiq olunmur** | Yalnız istifadəçinin öz seçdiyi lokal fayllar | — |
| **Prompt injection** | — | **Tətbiq olunmur** | LLM yoxdur | — |
| **Data leakage** | — | **Yoxdur** | Şəbəkə kodu yoxdur | — |
| **Asılılıq zəiflikləri** | `requirements.txt` | **⏳ Yoxlanmayıb** | `pip-audit`/`safety` quraşdırılmayıb; versiyalar sabitlənib (yaxşı), amma avtomatik zəiflik yoxlaması yoxdur | `pip-audit -r requirements.txt` dövri işlətmək |
| **İcazələr** | `logs/` tətbiq qovluğuna yazılır (`app.py:37-38`) | **Aşağı** | Proqram `Program Files` kimi yazılmaz qovluğa quraşdırılsa, log yazılmaz → `configure()` `os.makedirs`-də istisna ata bilər **[⏳ yoxlanmayıb]** | Logu `data_dir()`-ə köçürmək |
| **İnkişaf alətinin hook-ları** | `.claude/settings.json`, `.codex/hooks.json` | **Məlumat** | Hər AI alət çağırışında `C:/Users/LOQ/.local/bin/graphify.EXE` işlədilir — mütləq yol, başqa kompüterdə işləməz; tətbiqə təsiri yoxdur | — |

**Yekun [RƏY]:** tətbiqin hücum səthi çox kiçikdir. Əsas praktiki risk
təhlükəsizlikdən çox **məlumat bütövlüyüdür** (atomik olmayan saxlama).

---

## 13. Test analizi

### 13.1 Mövcud testlər

- **137 test faylı**, 38 863 sətir, **~2 489 test funksiyası** (parametrləşdirmə ilə daha çox).
- `tests/golden/{bl_1d,five_spot,five_spot_small}.json` — qızıl etalonlar.
- `tests/helpers.py` — ortaq qurucular (`default_scal()` — qrafikdə 380 əlaqə ilə ən mərkəzi node).
- `tests/conftest.py` — Qt-ni başsız rejimə keçirir.
- `pytest.ini` — `performance` markeri (7 yerdə).

### 13.2 Test növləri

| Növ | Varmı? | Nümunə |
|---|---|---|
| **Unit** | ✅ çox geniş | `test_pvt.py`, `test_variogram.py`, `test_domain.py`, `test_unit_conversions.py` |
| **Riyazi/ədədi doğrulama** | ✅ | Analitik Jakobian vs sonlu fərq (`test_implicit_jacobian.py`), MPFA-O manufactured həllər (`test_mpfa_o.py`), Buckley-Leverett (`test_analytical_bl.py`) |
| **Integration** | ✅ | `test_phase1_pipeline_integration.py`, `test_phase_b_production_integration.py`, `test_phase_d_mpfa_integration.py`, `test_facies_integration.py`, `test_sgs_integration.py` |
| **Reqressiya** | ✅ | `test_regression.py` (RF 16.840 %, 4314 addım), `test_numerical_benchmarks.py` |
| **Round-trip (saxla/aç)** | ✅ | `test_serialization.py`, `test_project_restore.py` |
| **UI** | ⚠️ qismən | `test_ui_static.py` (AST), `test_ui_wiring.py`, `test_panel_state.py`, `test_playback_controls.py` — `MainWindow`-un özü pytest-də **qurulmur** |
| **End-to-end (real pəncərə)** | ❌ avtomatik yox | Seans 46-da əl ilə «iki ayrı proses» sınağı aparılıb (`ISH_HESABATI.md`) |
| **Performans** | ✅ | `@pytest.mark.performance`, `test_sgs_performance.py`, `test_facies_performance.py` |

### 13.3 Coverage

**Ölçülməyib [FAKT]** — `coverage`/`pytest-cov` quraşdırılmayıb. Evristik
yoxlama (hər modulun ən azı bir test faylında birbaşa və ya paket vasitəsilə
idxal olunması) göstərdi ki, heç bir test faylının idxal **etmədiyi** modullar:
`ui/worker.py`, `ui/geology_map.py`, `ui/style.py`, `rendering/theme.py`,
`interfaces/geometry.py`, `interfaces/interpolation.py`. (`simulation/wellbore/*`
əvvəlcə siyahıya düşdü, lakin `test_wellbore_thp.py` və `test_thp_control.py`
onları paket vasitəsilə test edir.) **Sətir səviyyəsində coverage faizi
məlum deyil [⏳].**

### 13.4 Testlərin işə salınması — FAKTİKİ nəticə

```
Əmr:   ../venv/Scripts/python.exe -m pytest -q -p no:cacheprovider --tb=line -o addopts=""
Tarix: 26 sentyabr 2026, commit eacb029, Python 3.14.7, Windows 11
```

**Nəticə [FAKT]:**

```
2761 passed, 1 skipped, 1 xfailed in 527.80s (0:08:47)
EXIT=0
```

| Nə | Səbəb (əlavə əmrlə: `pytest -rsx tests/test_opm_import.py tests/test_implicit_newton.py`) |
|---|---|
| **1 skipped** | `tests/test_opm_import.py:6` — «resdata quraşdırılmayıb». Modul səviyyəsində atlanır, ona görə içindəki testlər toplanmır da (`pytest --co` → 2762 test). |
| **1 xfailed** | `tests/test_implicit_newton.py::test_crossing_the_bubble_point_now_converges` — gözlənilən uğursuzluq: «İki fazalı (neft-su) mühərrikdə sərbəst qaz fazası yoxdur, ona görə doymuş qolda mənfi görünən neft sıxılması Nyutonu dayandırır». Bu, bilinən və sənədləşdirilmiş məhdudiyyətdir, səhv deyil. |
| **0 failed** | Heç bir test uğursuz olmadı. |

**Müqayisə:** Seans 46 — «2773 keçdi, 1 xfailed» (`ISH_HESABATI.md`).
Fərq (12 test) tam olaraq `test_opm_import.py`-dəki testlərin sayıdır
(`tests/README.md`: 12) — yəni reqressiya yoxdur, sadəcə bu venv-də `resdata`
yoxdur. Müddət 13 dəq 18 san-dan 8 dəq 47 san-a düşüb (`-v` olmadan və fərqli
maşın yükü ilə — **[RƏY]** müqayisə ciddi performans ölçməsi deyil).

### 13.5 Test olunmayan və ya zəif test olunan vacib hissələr

1. **`MainWindow` (3 155 sətir)** — orkestrasiyanın özü: `run_simulation`, `_on_finished`, `open_project`, `closeEvent`, bərpa axını. Yalnız statik AST testi və əl ilə sınaq.
2. **`ui/worker.py`** — fon axınları, dayandırma (`request_stop`) davranışı.
3. **Atomik olmayan saxlama / yarımçıq fayl** — korlanmış `.imx`-in açılması halında davranış **[⏳ yoxlanmayıb]**.
4. **`resdata` asılı OPM idxalı** — venv-də `resdata` yoxdur; `test_opm_import.py` bu mühitdə **modul səviyyəsində atlanır** (12 test işləmir) — OPM idxalı bu venv-də yoxlanmır.
5. **Bilinən Jakobian qeyri-dəqiqlikləri** (TB-1, TB-2) — ölçülüb, amma düzəlişi yoxdur; test onları «qəbul edilmiş» kimi qıfıllamır **[⏳ ƏMİN DEYİL]**.

### 13.6 Mövcud testlərdə potensial problemlər

- **Uzunluq:** tam dəst 9–13 dəqiqə — hər dəyişiklikdən sonra işlədilməsi çətindir, CI yoxdur.
- **`tests/README.md` köhnəlib:** «Cəmi 1841 test» yazır, faktiki 2 700+.
- **Mühitdən asılılıq:** VTK və Qt offscreen rejimi tələb olunur; `resdata` olmadan bəzi testlər işləmir.
- **`berpa/A7_qaz_fazasi/tests/`** — `testpaths`-dan kənardadır (düzgündür), amma graphify qrafikində canlı testlərlə qarışır.

### 13.7 Əlavə edilməli testlər [RƏY]

1. `ProjectSerializer.save` üçün «yazma yarıda kəsildi» testi (atomik yazma əlavə ediləndən sonra).
2. `read_grdecl` üçün nəhəng `n*` sayı ilə qoruma testi.
3. `FullyImplicitEngine`-in `LinearSolverConfig`-i nəzərə aldığını yoxlayan test (düzəlişdən sonra).
4. `SimulationWorker` üçün `QThread` + `request_stop` inteqrasiya testi (offscreen).
5. `MainWindow.run_simulation`-u VTK olmadan qurmağa imkan verən «başsız» test fikstürü.

---

## 14. Problemlər və risklər

**Critical: tapılmadı.** Hər problem koda istinadla.

### High

| ID | Problem | Yer | Təsir | Həll | Prioritet | Mürəkkəblik |
|---|---|---|---|---|---|---|
| **P-01** | **UI monoliti və test boşluğu.** `MainWindow` 3 155 sətir / 121 metod; `panels.py` 2 216; `geology_service.py` 2 369. `MainWindow` pytest-də qurulmur (`session.py:3-5`). | `ui/main_window.py` | Hər UI dəyişikliyi reqressiya riski daşıyır; Seans 46-da `triggered(bool)` səhvi məhz bu boşluqdan uzun müddət gizli qalıb (`ISH_HESABATI.md`) | Tab başına ayrı controller siniflərinə bölmək; VTK səhnəsini inyeksiya edilə bilən etmək ki, `MainWindow` başsız qurulsun | Yüksək | Yüksək |
| **P-02** | **Layihə faylı atomik yazılmır.** | `application/serialization.py:375` | Yazma zamanı çökmə/disk dolması → yeganə `.imx` faylı korlanır | `path + ".tmp"`-yə yaz, `flush`+`fsync`, sonra `os.replace` | Yüksək | Aşağı |
| **P-03** | **Üç fazalı Jakobianın bilinən xətaları** — TB-1 (RATE istismarçısı, nisbi xəta 0.4893), TB-2 (axın təzyiq sütunu, 0.1936, mənbə ⏳). | `implicit/three_phase_residual.py`; `ROADMAP.md` «Texniki borc» | Nyuton daha çox iterasiya/kəsmə edir, çətin modellərdə yığılmama və «yumşaq qəbul» riski. Yığılmış həllin özü qalığa görə düzgündür **[RƏY]** | Sonlu fərq ilə müqayisəni hüceyrə/tənlik səviyyəsində lokallaşdırmaq (TB-2 üçün `ROADMAP`-da artıq plan var) | Yüksək | Yüksək |
| **P-04** | **SPE1CASE2 etalonundan ~145 gün qaz cəbhəsi fərqi.** | `benchmarks/spe1.py`, `SPE1.md` §8 | Üç fazalı proqnozların etibarlılığı sənaye etalonu ilə tam təsdiqlənməyib | `SPE1.md` §8.2-dəki 4 namizədin ölçülməsi | Yüksək | Orta-yüksək |

### Medium

| ID | Problem | Yer | Təsir | Həll | Prioritet | Mürəkkəblik |
|---|---|---|---|---|---|---|
| **P-05** | **Xətti həlledici konfiqurasiyası tam implicit mühərriklərdə işləmir.** `linear_solver` parametri qəbul edilir, saxlanılmır; `NewtonLinearSolver()` defoltla yaradılır. `config.linear_solver` (tolerans, ILU parametrləri) FIM-ə çatmır. | `implicit/engine.py:60, 84-86`; `three_phase_engine.py:55-62, 146-148` | İstifadəçi/skript `LinearSolverConfig`-i dəyişir, nəticə dəyişmir — səssiz təsirsizlik | `NewtonLinearSolverConfig`-i `SimulationConfig`-dən qurmaq və ya sahəni «yalnız IMPES» kimi açıq adlandırmaq | Orta | Aşağı |
| **P-06** | **Mühərrik hər qaçışda iki dəfə qurulur** — UI ön yoxlama üçün `create_engine`, sonra işçi `service.run` yenidən. | `main_window.py:2139` + `simulation_service.py:207` | Qurulma xərci (MPFA-O əmsalları, ilkin tarazlıq, Peaceman) ikiqat; böyük modellərdə UI donur, çünki birinci qurulma **əsas Qt axınında** gedir | Yalnız validasiya edən `service.validate(model, config)` metodu və ya qurulmuş mühərriki işçiyə ötürmək | Orta | Aşağı |
| **P-07** | **Qat qaydası iki yerdə pozulur**: `simulation → application.config` (4 idxal, dövri asılılıq), `ui → simulation` mühərrik sinifləri (eyni seçim kodu 3 dəfə). | §4.3; `main_window.py:2120, 1412, 1491` | Arxitektura sənədi ilə kod ayrılır; mühərrik seçimi composition root-dan kənarda | Konfiqurasiya dataclass-larını `domain/` və ya `interfaces/`-ə köçürmək; mühərrik seçimini servisə (`engine="IMPLICIT"` sətri ilə) vermək | Orta | Orta |
| **P-08** | **Kod təkrarı: iki FIM mühərriki.** `engine.py` və `three_phase_engine.py` eyni metod dəstini (`_update_rate_shares`, `_connection_mobilities`, `_surface_factors`, `_surface_rate_loop`, `_bhp_limit_loop`, `_thp_outer_loop`, `_record_snapshot`, run dövrəsi) ayrıca saxlayır. | `implicit/engine.py:292-396`; `three_phase_engine.py:440-533` | Quyu qaydasına düzəliş iki yerdə edilməlidir; biri unudulsa mühərriklər fərqli davranır (tarixçədə belə hallar olub — Q-33: «birləşmiş həlledici hələ köhnə qaydadadır») | Ortaq baza sinif və ya «quyu nəzarət dövrəsi» komponenti | Orta | Orta |
| **P-09** | **Gizli sabit**: FIM-də `growth_factor = config + 0.35` (1.15 → 1.5). | `implicit/engine.py:145` | `TimeSteppingConfig.growth_factor`-un mənası mühərrikə görə fərqlidir; sənədsiz | Ayrıca `implicit_growth_factor` sahəsi və ya şərh + sənəd | Orta | Aşağı |
| **P-10** | **Build/CI/coverage yoxdur**; `resdata` isteğe bağlı asılılıq heç yerdə elan olunmayıb; tam dəst 9–13 dəq. | kök qovluq | Reqressiyalar yalnız əl ilə işlədiləndə tutulur | `pyproject.toml` (extras: `opm = ["resdata"]`), GitHub Actions (Windows runner, `-m "not performance"`), `pytest-cov` | Orta | Aşağı-orta |
| **P-11** | **Vurma sırasının vahidi mühərrikə görə fərqlidir**: FIM-də səth, IMPES-də lay həcmi. | `simulation/results.py` (`well_water_injection_rate` şərhi) | «Müqayisə» tabında IMPES və FIM qaçışlarını yan-yana qoyanda yanlış nəticə **[⏳ UI-nin bunu göstərib-göstərmədiyi yoxlanmayıb]** | Nəticəyə vahid sahəsi əlavə etmək və ya IMPES-də səthə çevirmək | Orta | Aşağı |

### Low

| ID | Problem | Yer | Həll |
|---|---|---|---|
| **P-12** | **Köhnəlmiş sənədlər**: `README.md` — «Status: Başlanğıc mərhələsi», «Quraşdırma ⏳»; `io/grdecl.py:16-26` — corner-point «bərabər bloka approksimasiya olunur» (əslində `grdecl_import.py:124` həqiqi həndəsə qurur); `version.py` — `VERSION="69"`, `RELEASE_DATE="2026-08-28"`, halbuki Seans 45 (21 sentyabr) funksiyaları siyahıdadır; `tests/README.md` — «1841 test». | Sənədləri yeniləmək |
| **P-13** | **Yetim docstring**: `saturation_bound_relaxation`-un izahı `max_bhp_change`-dən sonra ayrıca sətir kimi qalıb. | `implicit/newton.py:90-100` | Sahənin yanına köçürmək |
| **P-14** | **GRDECL `n*` genişlənməsi hədsizdir** (§12). | `io/grdecl.py:110-117` | Hədd əlavə etmək |
| **P-15** | **Log-da tarix yoxdur** — `DATE_FORMAT = "%H:%M:%S"`; rotasiya olunan faylda günlər qarışır. | `logging_setup.py:25` | `"%Y-%m-%d %H:%M:%S"` |
| **P-16** | **Geniş `except Exception`** — `main_window.py`-də 17 dəfə; istifadəçiyə yalnız traceback-in son sətri göstərilir. `except Exception: pass` — 7 yerdə (`cross_validation.py:428`, `interpolation.py:109`, `opm_import.py:162,221`, `renderers.py:180,537`, `main_window.py:104`). | Səssiz udulan istisnaları ən azı `LOG.debug` ilə yazmaq |
| **P-17** | **`berpa/` arxiv kodu** graphify qrafikinə və axtarışa düşür (qrafik hesabatı canlı testi `berpa/`-dakı sinfə bağlayır). | `.graphifyignore` və ya qovluğu `docs/`-a/zip-ə köçürmək (sahibkarın qərarı — silinmir) |
| **P-18** | **`ProjectSerializer.load`** — `payload["project"]` yoxdursa `KeyError`, `ProjectFileError`-a çevrilmir (UI ümumi `except`-lə tutur). | `serialization.py:393` | `ProjectFileError`-a bükmək |
| **P-19** | **Vəziyyətli servis** — `ModelAwareSimulationService.create_engine` öz sahələrini dəyişir; paralel istifadədə təhlükəli olardı. İndi hər işçi `with_engine` ilə surət alır, ona görə real risk yoxdur. | `simulation_service.py:262-267` | Provider-ləri lokal dəyişənlərdə saxlamaq |
| **P-20** | **Log faylı tətbiq qovluğunda** — yazılmaz qovluğa quraşdırılanda problem. | `app.py:37-38` | `session.data_dir()` altına |

### Kateqoriyalar üzrə qısa qiymət [RƏY]

| Mövzu | Qiymət | Əsas |
|---|---|---|
| Kod keyfiyyəti | Yaxşı (nüvə) / zəif (UI) | Sənədləşmə sıx, adlar aydın; UI faylları çox böyük |
| Arxitektura | Yaxşı, 2 pozuntu | §4.3 |
| Təkrarlanan kod | Orta | P-08; mühərrik seçimi 3 dəfə |
| Performans | Ölçülür, sənədləşir | `PERFORMANCE.md`, `tools/benchmark*.py`; P-06 |
| Təhlükəsizlik | Yaxşı (kiçik səth) | §12 |
| Error handling | Düşünülmüş | Validasiya siyahıları, təhlükəsiz uğursuzluq; P-16 |
| Scalability | Tək nüvə, tək proses | `README.md` scope: «paralel/GPU yoxdur»; CPR böyük modellər üçün |
| Maintainability | Orta | Böyük fayllar, təkrar, 17 000+ sətir sənəd |
| Test coverage | Geniş, ölçülməyib | §13 |
| Documentation | Çox geniş, qismən köhnəlib | P-12 |
| Dependency | Sabitlənib, elan edilməmiş isteğe bağlı asılılıq | P-10 |
| İstifadə rahatlığı | Yaxşı | Azərbaycanca mesajlar, «nə etməli» göstərişi (`simulation_service.py:124-167`), bərpa faylı |

---

## 15. Təkmilləşdirmə təklifləri

Prioritet sırası ilə [RƏY]:

1. **Atomik saxlama** (P-02) — 10-15 sətirlik dəyişiklik, məlumat itkisi riskini aradan qaldırır.
2. **`LinearSolverConfig`-i FIM-ə ötürmək və ya adını dəqiqləşdirmək** (P-05), **gizli `+0.35`-i açıq etmək** (P-09).
3. **İkiqat mühərrik qurulmasını aradan qaldırmaq** (P-06) — həm sürət, həm UI donmasının qarşısı.
4. **`pyproject.toml` + CI** (P-10): Windows runner, `pytest -m "not performance"`, `IMEX_SKIP_SLOW=1`; `pytest-cov` ilə coverage hesabatı; `pip-audit`.
5. **FIM mühərriklərinin ortaq hissəsini çıxarmaq** (P-08).
6. **`MainWindow`-u bölmək** (P-01): hər tab üçün controller sinfi; VTK səhnəsini fabrik vasitəsilə inyeksiya edib başsız test.
7. **TB-1/TB-2 Jakobian xətaları** (P-03) və **SPE1 fərqi** (P-04) — `ROADMAP.md`-dakı planla.
8. **Sənədlərin yenilənməsi** (P-12): `README.md` §3–4, `io/grdecl.py` docstring, `version.py`, `tests/README.md`.
9. **Kiçik düzəlişlər**: log tarixi (P-15), GRDECL həddi (P-14), `KeyError` (P-18), səssiz `except` (P-16).

---

## 16. Layihəni öyrənmə planı

### 16.1 Əvvəlcə hansı faylları oxumalı (1-ci həftə)

1. `README.md` → `ARCHITECTURE.md` §2–4 → bu hesabatın §4–5 və diaqramları.
2. `app.py` (61 sətir) — hər şeyin necə bağlandığı.
3. `imex2d/application/config.py` (89 sətir) — bütün ədədi parametrlər.
4. `imex2d/domain/grid.py`, `domain/wells.py`, `domain/reservoir_model.py` — model nədir.
5. `imex2d/application/scenarios.py` (131 sətir) — sintetik 5-spot modeli necə qurulur.
6. `ARCHITECTURE.md` §7-dəki skript nümunəsi — UI-siz simulyasiya.

### 16.2 Modulların öyrənmə ardıcıllığı

```
domain  →  interfaces  →  application/model_builder + config
        →  simulation/discretization.py (TPFA) + well_model.py (Peaceman)
        →  simulation/impes_engine.py        (ən sadə mühərrik)
        →  simulation/implicit/{state, residual, jacobian, newton, time_stepping, engine}.py
        →  simulation/pvt/, scal_adapter.py, capillary.py
        →  application/simulation_service.py (hamısını birləşdirən)
        →  geology/{variogram, interpolation, sgs}.py + application/geology_service.py
        →  simulation/implicit/three_phase_*.py (ən mürəkkəb)
        →  discretization/ (MPFA-O)  →  ui/  →  history/
```

### 16.3 Əsas funksiyalar (bunları dərindən başa düş)

| Funksiya | Niyə |
|---|---|
| `ReservoirModelBuilder.build` (`model_builder.py:30`) | Geologiya → simulyasiya keçidi |
| `ModelAwareSimulationService.create_engine` (`simulation_service.py:262`) | Hansı fizika və mühərrikin seçildiyi |
| `TwoPointFluxDiscretization.build` (`discretization.py:85`) | Transmissivlik |
| `ResidualAssembler.residual` (`implicit/residual.py:76`) | Fizikanın ürəyi |
| `NewtonSolver.solve` (`newton.py:181`) | Qeyri-xətti həll |
| `AdaptiveTimeStepper.advance` (`time_stepping.py:118`) | Zaman addımı məntiqi |
| `FullyImplicitEngine.run` (`engine.py:169`) | Bütün dövrə |
| `OrdinaryKriging` (`geology/interpolation.py:332`) | Geostatistikanın ürəyi |

### 16.4 Ayrıca öyrənilməli texnologiyalar

1. **Python** — dataclass, ABC, type hints, `from __future__ import annotations`.
2. **NumPy** — vektorlaşdırma, fancy indexing, `np.add.at`.
3. **SciPy** — `scipy.sparse` (CSR/CSC), `sparse.linalg` (spsolve, spilu, cg, bicgstab), `optimize`.
4. **PyQt5** — QMainWindow, layout-lar, siqnal/slot, QThread.
5. **matplotlib** — Figure/Axes API (pyplot yox).
6. **pytest** — fikstür, `parametrize`, marker.
7. **Rezervuar mühəndisliyi əsasları** — Darcy, porozite, keçiricilik, doyumluluq, nisbi keçiricilik, PVT (Bo, Rs, Pb), black-oil. Tövsiyə: Ahmed «Reservoir Engineering Handbook»; Aziz & Settari «Petroleum Reservoir Simulation» **[RƏY — layihədə istinad yoxdur]**.
8. **Ədədi üsullar** — sonlu həcm, Nyuton-Rafson, seyrək xətti sistemlər.
9. **Geostatistika** — variogram, kriging (Deutsch & Journel «GSLIB» — kodda istinad var).

### 16.5 Layihəni özün dəyişmək üçün lazım olan biliklər

Python + NumPy (mütləq), fizikanın əsasları (qalıq tənliyini oxuya bilmək),
pytest (hər dəyişiklikdən sonra testlər), git, və layihənin öz qaydaları
(`CLAUDE.md`: hər iş `ISH_HESABATI.md`-a, hər texniki qərar `QARARLAR.md`-a yazılır).

### 16.6 Kiçik praktiki tapşırıqlar

| # | Tapşırıq | Nəyi öyrədir |
|---|---|---|
| 1 | `ARCHITECTURE.md` §7-dəki skripti işlət, `result.final_recovery_factor`-u çap et | Bütün boru xəttini UI-siz |
| 2 | Eyni skriptdə `SimulationConfig(end_time=1500, flux_scheme="MPFA-O")` + `FullyImplicitEngine` — RF-i TPFA ilə müqayisə et | Diskretizasiya seçimi |
| 3 | `CoreyParameters(no=3.0)` ilə RF necə dəyişir? | SCAL-ın təsiri |
| 4 | `NewtonConfig(cnv_tolerance=1e-4)` ötür, addım/iterasiya sayını müqayisə et | Konvergensiya meyarları |
| 5 | `tests/test_analytical_bl.py`-ni oxu və işlət; sonra bir parametri dəyiş | Doğrulama necə aparılır |
| 6 | `numune_quyular.csv`-yə bir quyu əlavə et, UI-də Kriging modeli qur, 3D-də bax | Geologiya boru xətti |
| 7 | `logging_setup.DATE_FORMAT`-a tarix əlavə et (P-15), test yaz, `ISH_HESABATI.md`-a qeyd yaz | Layihənin iş qaydası (kiçik, təhlükəsiz dəyişiklik) |
| 8 | `serialization.save`-i atomik et (P-02) + test | Real faydalı ilk töhfə |

---

## 17. Terminlər lüğəti

| Termin | İzah |
|---|---|
| **Rezervuar (lay)** | Yeraltı məsaməli süxur — neft/qaz/su onun boşluqlarındadır |
| **Porozite (φ, PORO)** | Süxurun boşluq payı (0–1) |
| **Keçiricilik (k, PERMX/Y/Z)** | Süxurun maye buraxma qabiliyyəti, mD (millidarsi) |
| **NTG** | Net-to-gross — layın məhsuldar hissəsinin payı |
| **Doyumluluq (Sw, So, Sg)** | Boşluğun neçə faizini su/neft/qaz tutur; cəmi 1 |
| **Nisbi keçiricilik (kr)** | Bir fazanın, digərləri olanda, nə qədər asan axdığı (0–1) |
| **SCAL** | Special Core Analysis — laboratoriya kr və Pc ölçmələri |
| **Kapilyar təzyiq (Pc)** | İki faza arasındakı təzyiq fərqi |
| **PVT** | Pressure-Volume-Temperature — flüid xassələri təzyiqə görə |
| **Bo, Bg, Bw** | Həcm əmsalı: lay şəraitindəki həcm / səth həcmi |
| **Rs** | Həll olmuş qaz-neft nisbəti (sm³/sm³) |
| **Pb** | Doyma (qabarcıq) təzyiqi — bundan aşağı qaz ayrılır |
| **Black-oil** | Neft, su, qazın sadələşdirilmiş 3 komponentli modeli |
| **OOIP / STOIIP, OGIP** | Yataqdakı ilkin neft / qaz həcmi |
| **RF** | Recovery Factor — çıxarılmış neftin OOIP-ə nisbəti, % |
| **WCT (su kəsri)** | Hasilatda suyun payı, % |
| **GOR** | Qaz-neft nisbəti |
| **BHP / THP** | Quyu dibi / quyu başı təzyiqi |
| **Peaceman WI** | Quyu indeksi — quyu ilə hüceyrə arasındakı «keçiricilik» |
| **Skin** | Quyu ətrafındakı zədələnmə/stimulyasiya əmsalı |
| **Nodal analiz, IPR, VLP** | Laydan axın (IPR) ilə lülə axınının (VLP) kəsişməsi — iş nöqtəsi |
| **5-spot** | 4 küncdə vurucu, mərkəzdə istismarçı (və ya tərsinə) quyu sxemi |
| **Grid, hüceyrə** | Yatağın kiçik bloklara bölünməsi |
| **Corner-point grid** | Hər hüceyrənin 8 küncü ayrıca verilən (əyri, maili) grid; COORD/ZCORN |
| **ACTNUM** | Hansı hüceyrənin aktiv olduğunu bildirən massiv |
| **Transmissivlik (T)** | İki hüceyrə arasındakı axın əmsalı |
| **TPFA / MPFA-O** | İki-nöqtəli / çoxnöqtəli axın approksimasiyası |
| **Tenzor keçiricilik** | İstiqamətə görə fərqli keçiricilik (3×3 matris) |
| **IMPES** | Implicit Pressure, Explicit Saturation |
| **FIM** | Fully Implicit Method — tam implicit |
| **Qalıq (residual)** | Tənliyin nə qədər pozulduğu; 0 = həll |
| **Jakobian** | Qalığın naməlumlara görə törəmələri matrisi |
| **CNV, MB** | Lokal və qlobal (kütlə balansı) yığılma meyarları |
| **CFL** | IMPES-də stabil addımın hədd şərti |
| **ILU, CPR, BiCGStab, CG** | Seyrək xətti sistemlərin həlli üçün üsullar/ön-şərtçilər |
| **Appleyard kəsməsi** | Nyuton addımının ölçüsünü məhdudlaşdırmaq |
| **Variogram** | Məsafəyə görə dəyişkənlik funksiyası (nugget, sill, range) |
| **Kriging** | Variogram əsaslı optimal çəkili interpolyasiya |
| **SGS / SIS** | Sequential Gaussian / Indicator Simulation |
| **Fasiya** | Süxur növü (qum, gil və s.) — kateqorik xassə |
| **History matching** | Model parametrlərini ölçülmüş hasilata uyğunlaşdırmaq |
| **NRMSE** | Normallaşdırılmış orta kvadratik xəta |
| **GRDECL, .DATA, SWOF, SGOF, PVTO, PVDG, PVTW** | Eclipse simulyatorunun fayl formatları və açar sözləri |
| **SPE1** | Society of Petroleum Engineers-in 1-ci etalon məsələsi |
| **OPM Flow** | Açıq mənbəli sənaye simulyatoru — etalon kimi işlədilir |
| **Composition root** | Bütün asılılıqların bağlandığı yeganə yer (`app.py`) |
| **Dependency injection** | Obyektə lazım olanı xaricdən vermək |
| **Provider** | Müəyyən xassəni (PVT, kr, Pc) hesablayan dəyişdirilə bilən obyekt |
| **`.imx`** | IMEX-2D layihə faylı (gzip + JSON) |

---

## 18. Açıq qalan suallar

1. ⏳ **Coverage faizi** — `pytest-cov` quraşdırılmayıb, ölçülməyib.
2. ⏳ **Asılılıq zəiflikləri** — `pip-audit` işlədilməyib.
3. ✅ ~~Pillow və `csv` idxal yerləri~~ — analiz zamanı yoxlanıldı (§7.2, §7.5).
4. ⏳ **P-11** — «Müqayisə» tabının IMPES və FIM vurma sıralarını eyni qrafikdə göstərib-göstərmədiyi.
5. ⏳ **Yazılmaz qovluqda** (`Program Files`) başladılanda `logging_setup.configure`-un davranışı.
6. ⏳ **TB-2-nin mənbəyi** — layihənin öz açıq sualıdır (`ROADMAP.md`).
7. ⏳ **SPE1CASE2 qalıq fərqi** — layihənin öz açıq sualıdır (`SPE1.md` §8.2).
8. ⏳ **`coupled_newton.py` / `standard_well.py`** — «defolt axında işləmir» (`docs/is_axini_ardicilligi.md` §10); onların gələcəyi (silinmə və ya qoşulma) sahibkarın qərarıdır.
9. ⏳ **`berpa/` qovluğunun saxlanma forması** — arxiv kimi qalsın, yoxsa zip/ayrı branch.
10. ⏳ **Layihənin məqsəd auditoriyası** — kommersiya məhsulu, yoxsa diplom/tədqiqat? Bu, P-01 və P-10-un prioritetini dəyişir.

---

## 19. Yekun nəticə

**[RƏY, faktlara əsaslanan]** IMEX-2D ciddi mühəndislik işidir: fizika nüvəsi
(black-oil, analitik Jakobian, MPFA-O, corner-point, adaptiv Δt, CPR) sənaye
simulyatorlarının istifadə etdiyi üsullarla qurulub, hər mühüm addım test
və sənədlə qorunur, qərarların səbəbləri yazılıb. Qatlı arxitektura əsasən
qorunur — `domain` tam müstəqildir, provider-lər inyeksiya olunur.

Zəif tərəflər nüvədə deyil, **ətrafdadır**: çox böyük və başsız test oluna bilməyən
UI sinfi, atomik olmayan fayl yazma, iki mühərrik arasında kod təkrarı, səssizcə
təsirsiz qalan konfiqurasiya sahəsi, build/CI infrastrukturunun olmaması.
Fizikada isə ən vacib açıq iş üç fazalı Jakobian xətaları və SPE1 etalonundan
qalan fərqdir — bunlar layihənin özü tərəfindən artıq ölçülüb və qeydə alınıb.

İlk addım kimi ən az xərclə ən çox fayda verən düzəlişlər: **P-02** (atomik
saxlama), **P-05/P-09** (konfiqurasiyanın dürüstlüyü), **P-06** (ikiqat qurulma)
və **P-10** (CI + coverage).
</content>
</invoke>
