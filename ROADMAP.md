# Yol xəritəsi — Lay Simulyatoru

Layihənin mərhələləri. Hər mərhələ bitəndə burada işarələnir və
[ISH_HESABATI.md](ISH_HESABATI.md)-na təfərrüatlı yazılır.

**İşarələr:** ✅ bitib · 🟡 gedir · ⬜ gözləyir · ❌ ləğv olunub

**Son yenilənmə:** 10 sentyabr 2026

---

> ## ⚠️ DİQQƏT — bu xəritə yenidən baxılmalıdır
>
> 10 sentyabr 2026 auditi mənzərəni dəyişdi: sahibkarın əvvəllər
> yazdığı **65 175 sətirlik işlək kod bazası** tapıldı
> (`lay-simulyatoru`, 1 993 keçən test). Aşağıdakı mərhələlərin
> **böyük hissəsi artıq mövcuddur** — sıfırdan yazılmayacaq.
>
> Tam mənzərə: [AUDIT_2026-09-10.md](AUDIT_2026-09-10.md)
>
> ⏳ **Xəritə iki strateji qərardan sonra yenidən yazılacaq:**
> (1) öz fizikamız, yoxsa OPM Flow? (2) hansı repoda davam edirik?

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
| 0.10 | Strateji qərar: öz fizikamız, yoxsa OPM Flow? | ⬜ |
| 0.11 | Strateji qərar: hansı repoda davam edilir? | ⬜ |

---

## Mövcud kodda NƏ VAR (audit nəticəsi)

Aşağıdakı mərhələ cədvəllərini oxumazdan əvvəl bunu nəzərə alın:

| Komponent | Vəziyyət |
|---|---|
| 3D corner-point grid (COORD/ZCORN) | ✅ hazır |
| **MPFA-O** (1 699 sətir) | ⚠️ yazılıb, **mühərriyə qoşulmayıb** |
| Kriging + variogram + SGS + fasiya | ✅ hazır (RBF-dən güclü) |
| Peaceman quyu indeksi | ✅ hazır |
| IMPES + FIM (Nyuton, analitik Jakobian, CPR) | ✅ hazır |
| CFL ilə adaptiv `dt` | ✅ hazır |
| Neft-su PVT və nisbi keçiricilik | ✅ hazır |
| **Qaz fazası + Stone II** (4 118 sətir) | ⚠️ **bərpa edilib**, qaytarılmayıb |
| GRDECL / Eclipse `.DATA` giriş-çıxışı | ✅ hazır |
| History matching + həssaslıq analizi | ✅ hazır |
| 5-spot reqressiya testi | ✅ hazır |
| **THP / VFP** | ❌ yoxdur |
| PyVista / Plotly | ❌ yoxdur (birbaşa VTK — bloklanır) |

---

## Mərhələ 1 — Pre-processing və 3D geologiya

**Məqsəd:** kəşfiyyat quyularının səpələnmiş ölçmələrindən tam 3D heterogen
lay modeli qurmaq və onu ekranda görmək.

| # | Tapşırıq | Modul | Status |
|---|---|---|---|
| 1.1 | Giriş məlumat sxemi (`.json` / `.csv`) və oxuyucu | `data/`, `geostats/` | ⬜ |
| 1.2 | 3D structured grid generatoru (nx, ny, nz, DX/DY/DZ, TOPS) | `grid/` | ⬜ |
| 1.3 | Hüceyrə həndəsəsi: mərkəz, həcm, üz sahəsi və normalı | `grid/` | ⬜ |
| 1.4 | Corner-point həndəsə (COORD / ZCORN) və ACTNUM | `grid/` | ⬜ |
| 1.5 | RBF interpolyasiya ilə Phi və K tenzorunun 3D yayılması | `geostats/` | ⬜ |
| 1.6 | Kriging (Ordinary / Universal) alternativi — PyKrige | `geostats/` | ⬜ |
| 1.7 | Keçiricilik tenzorunun müsbət-müəyyənliyinin yoxlanması | `geostats/` | ⬜ |
| 1.8 | PyVista ilə 3D lay mesh-i + kəşfiyyat quyularının renderi | `viz/` | ⬜ |
| 1.9 | Testlər: grid həcm balansı, interpolyasiya nöqtə dəqiqliyi | `tests/` | ⬜ |

**Bitmə şərti:** kəşfiyyat quyusu faylı verilir → ekranda rəngli 3D heterogen
lay modeli görünür, quyu trayektoriyaları üstündə.

---

## Mərhələ 2 — Fiziki və PVT modulları

**Məqsəd:** flüidin təzyiqə reaksiyasını və faza axıcılığını modelləşdirmək.

| # | Tapşırıq | Modul | Status |
|---|---|---|---|
| 2.1 | PVT cədvəl strukturu (PVTO, PVTG, PVTW analoqu) | `pvt/` | ⬜ |
| 2.2 | Təzyiqə görə Bo, Bg, Bw, Rs, µ interpolyasiyası | `pvt/` | ⬜ |
| 2.3 | Doyma təzyiqi (Psat) və saturated/undersaturated ayrımı | `pvt/` | ⬜ |
| 2.4 | Sıxılma əmsalları və törəmələr (∂b/∂p — FIM üçün lazım) | `pvt/` | ⬜ |
| 2.5 | 2-fazalı nisbi keçiricilik əyriləri (Krw-Kro, Krg-Krog) | `relperm/` | ⬜ |
| 2.6 | Stone II modeli ilə 3-fazalı Kro | `relperm/` | ⬜ |
| 2.7 | Kapilyar təzyiq (Pcow, Pcog) | `relperm/` | ⬜ |
| 2.8 | Testlər: monotonluq, son nöqtə, fiziki sərhədlər | `tests/` | ⬜ |

**Bitmə şərti:** istənilən (P, Sw, Sg) üçün bütün flüid və axın xassələri
fiziki cəhətdən düzgün qaytarılır, testlərlə təsdiqlənir.

---

## Mərhələ 3 — MPFA-O mühərriki və Peaceman quyu modeli

**Məqsəd:** layihənin əsas texniki fərqi — dəqiq axın diskretizasiyası.

| # | Tapşırıq | Modul | Status |
|---|---|---|---|
| 3.1 | Interaction volume-ların (təpə nöqtəsi ətrafı) qurulması | `numerical/` | ⬜ |
| 3.2 | Yerli lokal sistem: təzyiq davamlılığı + axın davamlılığı | `numerical/` | ⬜ |
| 3.3 | MPFA-O transmissibilite əmsallarının çıxarılması | `numerical/` | ⬜ |
| 3.4 | Qlobal `scipy.sparse` matrisinin yığılması | `numerical/` | ⬜ |
| 3.5 | Doğrulama: TPFA ilə müqayisə (ortoqonal gridd nəticə üst-üstə düşməli) | `tests/` | ⬜ |
| 3.6 | Doğrulama: anizotrop tenzorda analitik həll ilə müqayisə | `tests/` | ⬜ |
| 3.7 | Peaceman 3D quyu indeksi (WI), anizotrop r_e düsturu | `wells/` | ⬜ |
| 3.8 | VFP: hidrostatik sütun + sürtünmə itkisi → BHP ↔ THP | `wells/` | ⬜ |
| 3.9 | Quyu rejimləri: sabit BHP, sabit debet, rejim keçidi | `wells/` | ⬜ |
| 3.10 | Testlər: axın balansı (flux conservation), WI analitik yoxlama | `tests/` | ⬜ |

**Bitmə şərti:** MPFA-O matrisi qurulur, axın balansı maşın dəqiqliyində
saxlanılır, anizotrop testdə TPFA-dan üstünlüyü ədədlə göstərilir.

---

## Mərhələ 4 — Simulyasiya mühərriki (IMPES → FIM)

**Məqsəd:** zaman üzrə həqiqi simulyasiya.

| # | Tapşırıq | Modul | Status |
|---|---|---|---|
| 4.1 | İlkin şərtlər: hidrostatik tarazlıq (equilibration) | `engine/` | ⬜ |
| 4.2 | Təzyiq tənliyinin qurulması və həlli (A·P = B) | `engine/` | ⬜ |
| 4.3 | Xətti həlledici: SuperLU → BiCGStab/AMG (böyük modellər üçün) | `engine/` | ⬜ |
| 4.4 | Upwind sxemi ilə faza axınları | `engine/` | ⬜ |
| 4.5 | Doymaların yenilənməsi (Sw, So, Sg) və normalizasiya | `engine/` | ⬜ |
| 4.6 | Qazın ayrılması: P < Psat olduqda Rs nəzarəti və faza keçidi | `engine/` | ⬜ |
| 4.7 | CFL şərtinə görə dinamik zaman addımı (dt) | `engine/` | ⬜ |
| 4.8 | Quyuların həlledici ilə birləşdirilməsi (well coupling) | `engine/` | ⬜ |
| 4.9 | Kütlə balansı yoxlaması hər addımda | `engine/` | ⬜ |
| 4.10 | FIM (Fully Implicit) sxemi — Newton-Raphson + Jacobian | `engine/` | ⬜ |
| 4.11 | 5-spot ssenarisinin işə salınması (1 vurma + 4 hasilat) | `main.py` | ⬜ |

**Bitmə şərti:** 5-spot modeli günbəgün işləyir, kütlə balansı qorunur,
su cəbhəsi fiziki cəhətdən düzgün irəliləyir.

---

## Mərhələ 5 — Analitika və Recovery Factor

| # | Tapşırıq | Modul | Status |
|---|---|---|---|
| 5.1 | STOIIP hesablanması (həcm üsulu ilə) | `analytics/` | ⬜ |
| 5.2 | Kumulyativ Neft / Su / Qaz hasilatı | `analytics/` | ⬜ |
| 5.3 | Water Cut (%) və GOR | `analytics/` | ⬜ |
| 5.4 | Layın orta təzyiqi (həcm-çəkili) | `analytics/` | ⬜ |
| 5.5 | Günbəgün RF (%) və hasilat profili | `analytics/` | ⬜ |
| 5.6 | Nəticələrin fayla yazılması (CSV / JSON) | `analytics/` | ⬜ |

**Bitmə şərti:** simulyasiyadan sonra tam hasilat hesabatı avtomatik çıxır.

---

## Mərhələ 6 — Post-processing və 3D canlı vizualizasiya

| # | Tapşırıq | Modul | Status |
|---|---|---|---|
| 6.1 | PyVista 3D render: doyma və təzyiq sahələri | `viz/` | ⬜ |
| 6.2 | Su/qaz cəbhəsinin irəliləmə animasiyası | `viz/` | ⬜ |
| 6.3 | İnteraktiv kəsik (slice plane) aləti | `viz/` | ⬜ |
| 6.4 | Volumetric rendering | `viz/` | ⬜ |
| 6.5 | Matplotlib / Plotly dashboard: debet, THP/BHP, RF | `viz/` | ⬜ |
| 6.6 | Nəticə animasiyasının video/GIF kimi ixracı | `viz/` | ⬜ |

**Bitmə şərti:** istifadəçi simulyasiyanı işə salır, 3D-də cəbhənin
hərəkətini izləyir və dashboard-da bütün göstəriciləri görür.

---

## Gələcək ideyalar (backlog)

Hələ mərhələyə salınmayan, amma unudulmaması lazım olan fikirlər:

- Qırılma (fault) modeli və NNC (non-neighbor connection) dəstəyi
- ECLIPSE `.DATA` deck oxuyucusu — sənaye modelləri ilə uyğunluq
- SPE benchmark testləri (SPE1, SPE9) ilə doğrulama
- Paralel hesablama (çox nüvəli / GPU)
- History matching (tarixi uyğunlaşdırma) modulu
- Qeyri-struktur (PEBI) grid dəstəyi
