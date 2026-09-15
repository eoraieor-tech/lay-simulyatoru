# Yol xəritəsi — Lay Simulyatoru

Layihənin mərhələləri. Hər mərhələ bitəndə burada işarələnir və
[ISH_HESABATI.md](ISH_HESABATI.md)-na təfərrüatlı yazılır.

**İşarələr:** ✅ bitib · 🟡 gedir · ⬜ gözləyir · ❌ ləğv olunub

**Son yenilənmə:** 10 sentyabr 2026

---

> ## ℹ️ Bu xəritə 10 sentyabr 2026-da kodla tutuşduruldu
>
> Sahibkarın 3D 3-fazalı planı **sətir-sətir kod bazası ilə yoxlanıldı**
> (Seans 3). Nəticə: planın **~80 %-i artıq mövcuddur** və testlərlə
> doğrulanıb. Aşağıdakı cədvəllərdəki statuslar həmin ölçmənin nəticəsidir
> — `Modul` sütunu artıq plandakı fərzi adları yox, **REAL fayl yollarını**
> göstərir.
>
> Tam təhlil: [ISH_HESABATI.md](ISH_HESABATI.md) → Seans 3 ·
> [AUDIT_2026-09-10.md](AUDIT_2026-09-10.md)
>
> **Həqiqətən qalan iş:** THP/VFP · ~~qaz mühərrikinin servisə
> qaytarılması~~ ✅ B2 · ~~qaz Nyuton rəqsi~~ ✅ B3-A + B3-B ·
> ~~MPFA-nın UI-dən seçilməsi~~ ✅ B1 · animasiya/slice · Pcog ·
> ~~CSV ixracı~~ ✅ B5-a.
>
> **11 sentyabr 2026:** B3-B bağlandı — üç fazalı mühərrik artıq
> Pb = 300 bar-a qədər yığılır (doymamış Bo qolu düzəldildi, bax
> [QARARLAR.md](QARARLAR.md) → Q-09).
>
> **11 sentyabr 2026 (davamı):** **B4-A** bitdi — quyu başı təzyiqi
> (THP) çoxseqmentli şaquli traversdən hesablanır və BHP ilə yanaşı
> dashboard-a düşür (Q-11). Yol boyu iki əsassız status düzəldildi:
> `well_bhp` sahəsi ümumiyyətlə yox idi, GOR isə heç bir qrafikdə
> çəkilmirdi — hər ikisi indi var.
>
> Qalan bloklar: **B4-B** (`ControlMode.THP`), **B5-b** (Pcog),
> **B6** (3D animasiya), **B7** (yekun doğrulama).
>
> **13 sentyabr 2026:** **B6** (kəsik, oynatma, GIF) və **B4-B** bitdi —
> quyu artıq birbaşa **THP ilə idarə olunur** (açıq birləşmə, Q-15).
>
> **14 sentyabr 2026:** **B5-b** bitdi — üç fazalı mühərrikdə kapilyar
> təzyiq: **Pcow qoşuldu** (əvvəl səssizcə atılırdı!) və **Pcog** əlavə
> olundu (Q-16).
>
> **15 sentyabr 2026:** THP idarəsi **nodal analizə** keçirildi (Q-17) —
> real modeldə Δt medianı 38 dəfə yaxşılaşdı, quyular isə süni bağlanmadan
> hasilat verir.
>
> **B7 başladı:** SPE1 üçün üç boşluq tapıldı; **qaz vurulması** ✅ bitdi
> (Q-18 — əvvəl bütün vurucular su vururdu). Qalan: RATE + BHP limiti,
> səth debiti hədəfi, SPE1 modelinin qurulması.
>
> **Seans 27:** BHP limitinə hazırlıqda **səssiz səhv** tapılıb düzəldildi —
> çox perforasiyalı RATE quyusu hədəfin N qatını hasil edirdi (Q-19).
>
> **B7 addım 2 ✅ bitdi:** RATE quyusunda **BHP limiti** — quyu limitə çatanda
> BHP idarəsinə keçir, histerezislə RATE-ə qayıdır (Q-21). Yol boyu Eclipse
> ixracında iki səssiz səhv düzəldildi (Q-20). Qalan: səth debiti hədəfi, SPE1.
>
> **Seans 28:** qaz vuran qaçışda sahə **GOR = 0** səhvi (vurulan qaz hasilatla
> cəmlənib sıfırlanırdı) düzəldildi; **vurulan qaz** qrafikə və ixraca əlavə
> olundu (Q-22).
>
> 📋 **Bu işlərin addım-addım icra planı:** [ICRA_PLANI.md](ICRA_PLANI.md)
> (B1–B7 blokları, qəbul meyarları M1–M8).

---

## Faza 0 — Hazırlıq

| # | Tapşırıq | Status |
|---|---|---|
| 0.1 | Repo qurulması və GitHub-a bağlanması | ✅ |
| 0.2 | Sənədləşdirmə strukturunun yaradılması | ✅ |
| 0.3 | Layihənin məqsədi və əhatə dairəsinin dəqiqləşdirilməsi | ✅ |
| 0.4 | Texnologiya yığınının seçilməsi | ✅ |
| 0.5 | Qovluq strukturunun qurulması | ✅ |
| 0.6 | Python 3.12 venv + asılılıqların quraşdırılması | ✅ |
| 0.7 | Mühit maneələrinin aşkarlanması (Smart App Control) | ✅ |
| 0.8 | **Mövcud kod bazasının auditi** | ✅ |
| 0.9 | **A7 qaz fazasının git tarixçəsindən bərpası** | ✅ |
| 0.10 | Strateji qərar: öz fizikamız, yoxsa OPM Flow? | ✅ öz fizikamız |
| 0.11 | Strateji qərar: hansı repoda davam edilir? | ✅ birləşdirildi (`C:\Dev\LSM`) |
| 0.12 | **Sahibkarın planının kod bazası ilə tutuşdurulması** | ✅ |

---

## Mövcud kodda NƏ VAR (audit nəticəsi)

Aşağıdakı mərhələ cədvəllərini oxumazdan əvvəl bunu nəzərə alın:

| Komponent | Vəziyyət |
|---|---|
| 3D corner-point grid (COORD/ZCORN) | ✅ hazır |
| **MPFA-O** (1 699 sətir) | ✅ **tam işlək** — FIM-ə qoşulub (Phase 5B-2) və UI-dən seçilir (B1) |
| Kriging + variogram + SGS + fasiya | ✅ hazır (RBF-dən güclü) |
| Peaceman quyu indeksi | ✅ hazır |
| IMPES + FIM (Nyuton, analitik Jakobian, CPR) | ✅ hazır |
| CFL ilə adaptiv `dt` | ✅ hazır |
| Neft-su PVT və nisbi keçiricilik | ✅ hazır |
| **Qaz fazası + Stone II** (4 323 sətir) | ⚠️ **kod bazasına qaytarılıb** (143 test keçir) — servis/UI qoşulması qalır |
| GRDECL / Eclipse `.DATA` giriş-çıxışı | ✅ hazır |
| History matching + həssaslıq analizi | ✅ hazır |
| 5-spot reqressiya testi | ✅ hazır |
| **THP / VFP** | ❌ yoxdur |
| PyVista / Plotly | ❌ yoxdur — birbaşa **VTK 9.7.0 işləyir** (49 test keçir) |

---

## Mərhələ 1 — Pre-processing və 3D geologiya

**Məqsəd:** kəşfiyyat quyularının səpələnmiş ölçmələrindən tam 3D heterogen
lay modeli qurmaq və onu ekranda görmək.

| # | Tapşırıq | Modul | Status |
|---|---|---|---|
| 1.1 | Giriş məlumat sxemi (`.json` / `.csv`) və oxuyucu | `geology/well_data_io.py`, `io/` | ✅ |
| 1.2 | 3D structured grid generatoru (nx, ny, nz, DX/DY/DZ, TOPS) | `domain/grid.py` | ✅ |
| 1.3 | Hüceyrə həndəsəsi: mərkəz, həcm, üz sahəsi və normalı | `domain/general_grid_geometry.py` | ✅ |
| 1.4 | Corner-point həndəsə (COORD / ZCORN) və ACTNUM | `domain/corner_point_geometry.py` | ✅ |
| 1.5 | Phi və K tenzorunun 3D yayılması | `geology/` | ✅ **kriginq/SGS ilə** (RBF yox — bax Q-08) |
| 1.6 | Kriging (Ordinary / Universal) | `geology/interpolation.py` | ✅ **öz kodumuz** (PyKrige lazım deyil) |
| 1.7 | Keçiricilik tenzorunun müsbət-müəyyənliyinin yoxlanması | `domain/properties.py` | ✅ |
| 1.8 | 3D lay mesh-i + quyuların renderi | `rendering/vtk_volume.py` | ✅ **VTK ilə** (PyVista yox) |
| 1.9 | Testlər: grid həcm balansı, interpolyasiya nöqtə dəqiqliyi | `tests/` | ✅ |

**Bitmə şərti:** kəşfiyyat quyusu faylı verilir → ekranda rəngli 3D heterogen
lay modeli görünür, quyu trayektoriyaları üstündə.

---

## Mərhələ 2 — Fiziki və PVT modulları

**Məqsəd:** flüidin təzyiqə reaksiyasını və faza axıcılığını modelləşdirmək.

| # | Tapşırıq | Modul | Status |
|---|---|---|---|
| 2.1 | PVT cədvəl strukturu | `domain/pvt.py` | ✅ |
| 2.2 | Təzyiqə görə Bo, Bg, Bw, Rs, µ interpolyasiyası | `simulation/pvt/` | ✅ |
| 2.3 | Doyma təzyiqi (Psat) və saturated/undersaturated ayrımı | `simulation/pvt/` | ✅ |
| 2.4 | Sıxılma əmsalları və törəmələr (∂b/∂p) | `domain/pvt.py` | ✅ |
| 2.5 | 2-fazalı nisbi keçiricilik əyriləri | `domain/scal*.py` | ✅ |
| 2.6 | Stone II modeli ilə 3-fazalı Kro | `simulation/stone_relperm.py` | ✅ analitik törəmələri ilə |
| 2.7 | Kapilyar təzyiq | `simulation/capillary.py` | 🟡 **Pcow ✅ · Pcog ❌** |
| 2.8 | Testlər: monotonluq, son nöqtə, fiziki sərhədlər | `tests/` | ✅ |

**Bitmə şərti:** istənilən (P, Sw, Sg) üçün bütün flüid və axın xassələri
fiziki cəhətdən düzgün qaytarılır, testlərlə təsdiqlənir.

---

## Mərhələ 3 — MPFA-O mühərriki və Peaceman quyu modeli

**Məqsəd:** layihənin əsas texniki fərqi — dəqiq axın diskretizasiyası.

| # | Tapşırıq | Modul | Status |
|---|---|---|---|
| 3.1 | Interaction volume-ların qurulması | `discretization/mpfa_o_interaction.py` | ✅ |
| 3.2 | Yerli lokal sistem | `discretization/mpfa_o_local_system.py` | ✅ |
| 3.3 | MPFA-O transmissibilite əmsalları | `discretization/mpfa_o.py` | ✅ |
| 3.4 | Qlobal `scipy.sparse` matrisinin yığılması | `discretization/mpfa_global.py` | ✅ |
| 3.4b | **MPFA-O-nun konfiqurasiya/UI-dən seçilməsi (B1)** | `application/`, `ui/panels.py` | ✅ |
| 3.5 | Doğrulama: TPFA ilə müqayisə | `tests/` | ✅ |
| 3.6 | Doğrulama: anizotrop tenzorda | `tests/` | ✅ |
| 3.7 | Peaceman 3D quyu indeksi (WI), anizotrop r_e | `simulation/well_model.py` | ✅ |
| **3.8** | **VFP: hidrostatik sütun + sürtünmə itkisi → BHP ↔ THP** | `simulation/wellbore/` | 🟡 **BHP → THP ✅** (B4-A: Chen sürtünməsi + çoxseqmentli traverse) · **tərs istiqamət (`ControlMode.THP`) ✅** (B4-B, açıq birləşmə) · VFPPROD ⏳ |
| 3.9 | Quyu rejimləri: sabit BHP, sabit debet | `domain/wells.py` | ✅ |
| 3.10 | Testlər: axın balansı, WI analitik yoxlama | `tests/` | ✅ |

**Bitmə şərti:** MPFA-O matrisi qurulur, axın balansı maşın dəqiqliyində
saxlanılır, anizotrop testdə TPFA-dan üstünlüyü ədədlə göstərilir.

---

## Mərhələ 4 — Simulyasiya mühərriki (IMPES → FIM)

**Məqsəd:** zaman üzrə həqiqi simulyasiya.

| # | Tapşırıq | Modul | Status |
|---|---|---|---|
| 4.1 | İlkin şərtlər: hidrostatik tarazlıq | `simulation/initialization/` | ✅ |
| 4.2 | Təzyiq tənliyinin qurulması və həlli | `simulation/impes_engine.py` | ✅ |
| 4.3 | Xətti həlledici | `simulation/linear_solver.py`, `implicit/cpr.py` | ✅ **CG+ILU və CPR** (SuperLU/AMG yox — bax Q-08) |
| 4.4 | Upwind sxemi ilə faza axınları | `implicit/residual.py` | ✅ |
| 4.5 | Doymaların yenilənməsi (Sw, So, Sg) | `implicit/three_phase_state.py` | ✅ |
| 4.6 | Qazın ayrılması, faza keçidi (variable switching) | `implicit/three_phase_*.py` | ✅ kod var, **servisə qoşulmayıb** |
| 4.7 | CFL şərtinə görə dinamik `dt` | `implicit/time_stepping.py` | ✅ |
| 4.8 | Quyuların həlledici ilə birləşdirilməsi | `implicit/coupled_newton.py` | ✅ |
| 4.9 | Kütlə balansı yoxlaması hər addımda | `implicit/newton.py` | ✅ |
| 4.10 | FIM — Newton-Raphson + analitik Jakobian | `implicit/` | ✅ |
| 4.11 | 5-spot ssenarisi | `application/scenarios.py` | ✅ reqressiya testi ilə |

**Bitmə şərti:** 5-spot modeli günbəgün işləyir, kütlə balansı qorunur,
su cəbhəsi fiziki cəhətdən düzgün irəliləyir.

---

## Mərhələ 5 — Analitika və Recovery Factor

| # | Tapşırıq | Modul | Status |
|---|---|---|---|
| 5.1 | STOIIP (bizdə adı `ooip`) | `simulation/results.py` | ✅ |
| 5.2 | Kumulyativ Neft / Su / Qaz hasilatı | `simulation/results.py` | ✅ |
| 5.3 | Water Cut (%) və GOR | `simulation/results.py` | ✅ |
| 5.4 | Layın orta təzyiqi | `simulation/results.py` | ✅ |
| 5.5 | Günbəgün RF (%) və hasilat profili | `simulation/results.py` | ✅ |
| 5.6 | Nəticələrin fayla yazılması | `reporting/report.py`, `reporting/results_export.py` | ✅ PDF · **CSV/JSON** (B5-a: BHP/THP daxil) |

**Bitmə şərti:** simulyasiyadan sonra tam hasilat hesabatı avtomatik çıxır.

---

## Mərhələ 6 — Post-processing və 3D canlı vizualizasiya

| # | Tapşırıq | Modul | Status |
|---|---|---|---|
| 6.1 | 3D render: doyma və təzyiq sahələri | `rendering/vtk_volume.py` | ✅ **VTK ilə** |
| 6.2 | Su/qaz cəbhəsinin irəliləmə animasiyası | `rendering/` | ❌ |
| 6.3 | İnteraktiv kəsik (slice plane) aləti | `rendering/` | ❌ (statik kəsik ✅) |
| 6.4 | Volumetric rendering | `rendering/volume.py` | 🟡 psevdo-3D var, əsl volume yox |
| 6.5 | Dashboard: debet, BHP, RF | `rendering/renderers.py` | ✅ matplotlib · THP/BHP paneli (B4-A, 3×2 düzüm) |
| 6.6 | Nəticə animasiyasının video/GIF ixracı | `rendering/` | ❌ |

**Bitmə şərti:** istifadəçi simulyasiyanı işə salır, 3D-də cəbhənin
hərəkətini izləyir və dashboard-da bütün göstəriciləri görür.

---

## Texniki borc

Bilinən, lakin hazırda proqramı çökdürməyən qeyri-dəqiqliklər — sahibkarın
qərarı ilə indi düzəldilmir, unudulmasın deyə burada saxlanılır.

| # | Nə | Ölçülüb | Mənbə |
|---|---|---|---|
| TB-1 | Üç fazalı RATE istismarçısının quyu Jakobianı: sərbəst qazın təzyiq/Sg törəmələri sıfır qoyulub | sonlu fərqə qarşı nisbi xəta **0.4893** (tək perforasiya), 0.3239 (iki) | Seans 27, `three_phase_residual.py` RATE budağı |

---

## Gələcək ideyalar (backlog)

Hələ mərhələyə salınmayan, amma unudulmaması lazım olan fikirlər:

- Qırılma (fault) modeli və NNC (non-neighbor connection) dəstəyi
- ECLIPSE `.DATA` deck oxuyucusu — sənaye modelləri ilə uyğunluq
- SPE benchmark testləri (SPE1, SPE9) ilə doğrulama
- Paralel hesablama (çox nüvəli / GPU)
- History matching (tarixi uyğunlaşdırma) modulu
- Qeyri-struktur (PEBI) grid dəstəyi
