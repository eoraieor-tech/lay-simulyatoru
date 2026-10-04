# Lay Simulyatoru (IMEX-2D)

> Neft-qaz layının (rezervuarın) masaüstü simulyatoru.
> Bu repo layihənin **bütün detallarını** — kodu, qərarları, mərhələləri və iş
> jurnalını — bir yerdə saxlayır.

**Status:** 🟡 İcra planının B1–B6 blokları bitib; **B7 (yekun doğrulama, SPE1)** gedir.
Cari vəziyyət: [ROADMAP.md](ROADMAP.md) · qalan işlərin siyahısı:
[ISH_HESABATI.md](ISH_HESABATI.md) → Seans 50 · açıq xətalar (X-2, X-3;
X-1 çökməsi Seans 53-də düzəldildi): [ROADMAP.md](ROADMAP.md) → «Açıq xətalar».
**Bu repoda sənədləşdirmə:** 10 sentyabr 2026-dan (kod bazası ondan əvvəl mövcud idi —
bax [AUDIT_2026-09-10.md](AUDIT_2026-09-10.md)).
**Son yenilənmə:** 4 oktyabr 2026 (Seans 53)

---

## 1. Layihə nədir

Lay Simulyatoru — neft-qaz yataqlarında lay (rezervuar) proseslərinin
riyazi modelləşdirilməsi və simulyasiyası üçün proqram təminatıdır
(Python + PyQt5, Windows).

**Məqsəd:** kəşfiyyat quyularının 3D koordinat və geoloji məlumatlarından
(məsaməlik, keçiricilik tenzoru) başlayaraq 3D heterogen lay modelini quran,
5-spot quyu şəbəkəsində su ilə sıxışdırma və qazın neftdən ayrılmasını (Rs)
günbəgün simulyasiya edən, RF / THP / BHP / lay təzyiqini hesablayan və
nəticələri interaktiv 3D-də göstərən platforma.

**Əsas texniki fərq:** klassik TPFA əvəzinə **MPFA-O** (Multi-Point Flux
Approximation, O-sxemi) — anizotrop keçiricilik tenzorlarında və
qeyri-ortoqonal / corner-point gridlərdə axın dəqiqliyini qorumaq üçün.

**Əhatə dairəsi (scope):**

| Var | Yoxdur |
|---|---|
| 3D corner-point grid (COORD/ZCORN), ACTNUM | Kompozisiya (EOS) modeli |
| 3 faza: neft, su, qaz (black-oil, Rv = 0) | Termal (buxar) proseslər |
| TPFA və MPFA-O diskretizasiyası (MPFA-O yalnız iki fazalı tam implicit mühərrikdə) | Qeyri-struktur (PEBI) grid |
| Geostatistik interpolyasiya (Kriging, SGS, SIS) | Paralel / GPU hesablama |
| IMPES + tam implicit (FIM) həlledici, CPR | Eclipse `VFPPROD` cədvəl idxalı |
| Peaceman quyuları; BHP / RATE / THP idarəsi, BHP limiti, səth debiti | Lülədə sürüşmə (slip) modeli |
| 3D vizualizasiya (VTK), dashboard, günlük göstəricilər | |
| History matching və həssaslıq analizi | |
| GRDECL / Eclipse deck oxuma, `.DATA` ixracı, CSV / JSON / PDF | |

Layihə sıfırdan yazılmayıb: 10 sentyabr 2026 auditi göstərdi ki, işin böyük
hissəsi sahibkarın əvvəlki `lay-simulyatoru` kod bazasında artıq mövcud idi və
onun üzərində davam edildi (bax [QARARLAR.md](QARARLAR.md) → Q-04).

---

## 2. Sənədlərin xəritəsi

Layihənin hər aspekti ayrıca sənəddə qeyd olunur. Bir şey axtarırsansa,
başlanğıc nöqtəsi budur:

**Vəziyyət və tarixçə**

| Sənəd | Nə var içində |
|---|---|
| [ROADMAP.md](ROADMAP.md) | Mərhələlər, nə bitib / nə qalıb, texniki borc |
| [ISH_HESABATI.md](ISH_HESABATI.md) | İş jurnalı — hər seans, hər ölçmə, tarixlə (yalnız əlavə olunur) |
| [QARARLAR.md](QARARLAR.md) | Texniki qərarlar (Q-01 … Q-36) və **niyə** məhz belə seçildi |
| [ICRA_PLANI.md](ICRA_PLANI.md) | B1–B7 blokları və M1–M8 qəbul meyarları |
| [TEHVIL_TESLIM.md](TEHVIL_TESLIM.md) | İşi başqa kompüterdə davam etdirmək üçün: mühit, tələlər, növbəti iş |
| [SPE1.md](SPE1.md) | SPE1CASE2 etalonu: parametrlər, OPM Flow ilə müqayisə, qalan fərq |

**Quruluş və nəzəriyyə**

| Sənəd | Nə var içində |
|---|---|
| [ARCHITECTURE.md](ARCHITECTURE.md) | Qatlar, modullar, dizayn qərarları, skript rejimi |
| [PROJECT_ANALYSIS.md](PROJECT_ANALYSIS.md) · [xülasə](PROJECT_ANALYSIS_SUMMARY.md) | 26 sentyabr 2026 texniki analizi, problemlər P-01…P-20 |
| [docs/](docs/README.md) | Nəzəri əsaslar, icra ardıcıllığı, MPFA-O spesifikasiyası, diaqramlar, təqdimat |
| [UNITS.md](UNITS.md) | Vahid sistemi və çevirmə qatı |

**Modul sənədləri**

| Mövzu | Sənəd |
|---|---|
| Geostatistika | [INTERPOLATION_CORE.md](INTERPOLATION_CORE.md) · [PROPERTY_ENGINE.md](PROPERTY_ENGINE.md) · [GEOSTATISTICS.md](GEOSTATISTICS.md) · [SGS.md](SGS.md) · [FACIES.md](FACIES.md) · [FACIES_INTEGRATION.md](FACIES_INTEGRATION.md) · [LAYER_AWARE_MODELING.md](LAYER_AWARE_MODELING.md) |
| Fizika və mühərrik | [A6_PLAN.md](A6_PLAN.md) (tam implicit) · [A7_PLAN.md](A7_PLAN.md) (qaz fazası, tarixi) · [SCAL.md](SCAL.md) · [FAULTS.md](FAULTS.md) · [PERFORMANCE.md](PERFORMANCE.md) |
| Giriş / çıxış | [ECLIPSE_IO.md](ECLIPSE_IO.md) · [OPM_IMPORT.md](OPM_IMPORT.md) · [REPORTING.md](REPORTING.md) · [VISUALIZATION.md](VISUALIZATION.md) |
| Uyğunlaşdırma | [HISTORY_MATCHING.md](HISTORY_MATCHING.md) · [SENSITIVITY.md](SENSITIVITY.md) |
| Testlər | [tests/README.md](tests/README.md) |
| Tarixi | [AUDIT_2026-09-10.md](AUDIT_2026-09-10.md) · [berpa/](berpa/A7_qaz_fazasi/MENBE.md) (git tarixçəsindən bərpa olunmuş köhnə qaz kodu — arxiv) |

Modul sənədləri yazıldıqları anın vəziyyətini təsvir edir; ziddiyyət olanda
`ROADMAP.md` və `ISH_HESABATI.md`-nin son bölmələri əsasdır.

---

## 3. Quraşdırma

Tələb: Windows, Python 3.12 və ya 3.14 (hər ikisində tam test dəsti keçib).
Asılılıqlar dəqiq versiyalarla sabitlənib ([requirements.txt](requirements.txt)):
NumPy 2.5.2, SciPy 1.18.1, PyQt5 5.15.11, matplotlib 3.11.1, VTK 9.7.0, Pillow 12.3.0.

```bat
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements-dev.txt
```

`run.bat` virtual mühiti bu sıra ilə axtarır: `.venv\`, `..\venv\`, `venv\`.

⚠️ **Smart App Control işlək olan Windows-da** təzə endirilmiş imzasız DLL-lər
(xüsusən VTK) bloklana bilər. Həll yolu və səbəbi: [QARARLAR.md](QARARLAR.md) →
Q-07, Q-13.

`pyproject.toml`, quraşdırıcı (installer) və CI yoxdur — proqram repo
qovluğundan işlədilir.

---

## 4. İstifadə

```bat
run.bat                 :: və ya: python app.py
python -m pytest -q     :: testlər (2762 test; tam dəst ~9–13 dəqiqə)
```

Proqramda iş axını:

1. Sol paneldə bölmələr doldurulur: grid, geologiya (quyu cədvəli və ya
   sintetik), süxur/flüid, faultlar, nisbi keçiricilik, PVT, quyular, ədədi
   parametrlər.
2. **«MODELİ İŞƏ SAL»** — hesablama fon axınında gedir.
3. Nəticələr 13 tabda: Layihə, Model, Nəticələr, Günlük göstəricilər, Nisbi
   keçiricilik, 3D görüntü, PVT, Validasiya (B-L), Müqayisə, Tarixçə,
   Uyğunlaşdırma, Həssaslıq, Jurnal.
4. «Layihə» menyusu: layihəni saxla/aç (`.imx`), son layihələr, GRDECL oxu,
   Eclipse deck yaz, CSV/JSON ixracı, PDF hesabat.

Nümunə giriş faylları kök qovluqdadır: `numune_quyular.csv`,
`numune_grid.GRDECL`, `numune_scal.csv`, `numune_faultlar.csv`,
`numune_musahide.csv`.

İnterfeyssiz (skript) rejim: [ARCHITECTURE.md](ARCHITECTURE.md) §7.
Köməkçi alətlər (`tools/`): `benchmark.py`, `golden.py`, `spe1_compare.py`,
`eclipse_summary.py`, `module_graph.py`, `presentation_demo.py`.

Jurnal faylı: `logs/imex2d.log`.

---

## 5. Bu repo necə işləyir

Qayda sadədir: **layihədə baş verən hər şey burada qeyd olunur.**

- Yeni tələb gəldi → [ROADMAP.md](ROADMAP.md)-ə mərhələ kimi yazılır
- Texniki seçim edildi → [QARARLAR.md](QARARLAR.md)-a səbəbi ilə yazılır
- İş görüldü → [ISH_HESABATI.md](ISH_HESABATI.md)-na tarixlə yazılır
- Struktur dəyişdi → [ARCHITECTURE.md](ARCHITECTURE.md) yenilənir
- Hər mərhələdən sonra → `git commit` + `git push`

Detallı qaydalar: [CLAUDE.md](CLAUDE.md)

---

## 6. Əlaqə

**Müəllif:** eoraieor-tech
**Repo:** https://github.com/eoraieor-tech/lay-simulyatoru
