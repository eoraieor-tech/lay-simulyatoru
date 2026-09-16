# SPE1 etalonu (B7 addım 3) — mənbə, doğrulanmış parametrlər, boşluqlar

**Tarix:** 15 sentyabr 2026 · **Status:** parametrlər mənbədən yoxlanıldı,
model HƏLƏ qurulmayıb · bax `ISH_HESABATI.md` → Seans 29, `QARARLAR.md` → Q-23

---

## 1 · Mənbə

| | |
|---|---|
| Orijinal məqalə | Odeh, A.S. (1981), *Comparison of Solutions to a Three-Dimensional Black-Oil Reservoir Simulation Problem*, JPT, yanvar 1981 |
| Deck faylları | OPM layihəsi, `opm-tests` repozitoriyası, `spe1/` qovluğu |
| `SPE1CASE1.DATA` | https://raw.githubusercontent.com/OPM/opm-tests/master/spe1/SPE1CASE1.DATA — 436 sətir, sha256 `7e30d000d61aa6c0d9cc9bc5d9841b5f59a17c0dee628f69d8d4af8d7fcccb1e` |
| `SPE1CASE2.DATA` | https://raw.githubusercontent.com/OPM/opm-tests/master/spe1/SPE1CASE2.DATA — sha256 `f3de3d06ab5705381e14e6902a21c1af7bb6b37dbc9d6e4055fe146f50a4c249` |
| Etalon nəticə | `spe1/opm-simulation-reference/flow/SPE1CASE2.SMSPEC` + `.UNSMRY` (OPM Flow) |

Fayllar 15 sentyabr 2026-da endirilib. Üçüncü tərəf faylları olduğu üçün (lisenziya
yoxlanılmayıb ⏳) repoya **əlavə edilməyib** — yuxarıdakı URL və checksum ilə
bərpa olunur.

**CASE1 ↔ CASE2 fərqi** (şərhlər çıxarılmaqla `diff` ilə yoxlanıldı): yalnız başlıq
və `DRSDT 0` sətri. CASE1-də sərbəst qaz doymamış neftdə həll ola bilmir
(sabit doyma təzyiqi), CASE2-də həll olur (dəyişən doyma təzyiqi). Bütün
cədvəl və quyu parametrləri eynidir.

## 2 · Doğrulanmış parametrlər (deck-dən, FIELD vahidləri)

| Parametr | Deck | Təhvil sənədindəki "yayılmış" dəyər |
|---|---|---|
| Grid | `DIMENS 10 10 3`, `DX = DY = 300*1000` ft | ✅ eyni |
| Qalınlıq | `DZ 100*20 100*30 100*50` ft | ✅ eyni |
| Üst dərinlik | `TOPS 100*8325` ft | — (sənəddə yox idi) |
| Keçiricilik | `PERMX = PERMY = PERMZ = 500 / 50 / 200` mD | ✅ (PERMZ də eynidir) |
| Məsaməlilik | `PORO 300*0.3` | ✅ |
| Fazalar | `OIL GAS WATER DISGAS` (VAPOIL yox) | — |
| Süxur | `ROCK 14.7 3E-6` — istinad 14.7 psia, 3·10⁻⁶ 1/psi | — |
| Su | `PVTW 4017.55 1.038 3.22E-6 0.318 0.0` | — |
| Sıxlıq (səth) | `DENSITY 53.66 64.49 0.0533` lb/ft³ (neft, su, qaz) | — |
| İlkin tarazlıq | `EQUIL 8400 4800 8450 0 8300 0` — 4800 psia @ 8400 ft, WOC 8450 ft, GOC 8300 ft | ✅ 4800 psia @ 8400 ft |
| İlkin Rs | `RSVD` 8300–8450 ft: 1.270 Mscf/STB (sabit) | — |
| Doyma təzyiqi | PVTO-da Rs = 1.270 → 4014.7 psia | ✅ ~4014.7 |
| İlkin Sw | SWOF-un ilk sətri 0.12 | ✅ |
| İstismarçı | `(10,10,3)`, `WCONPROD 'ORAT' 20000 4* 1000` — neft 20 000 STB/gün, min BHP 1000 psia; BHP istinadı 8400 ft | ✅ |
| Vurucu | `(1,1,1)`, `WCONINJE 'GAS' 'RATE' 100000 1* 9014` — qaz 100 000 Mscf/gün, maks BHP 9014 psia; istinad 8335 ft | ✅ (100 MMscf/gün) |
| Quyu radiusu | `COMPDAT` diametr 0.5 ft → rw = 0.25 ft | — |
| Müddət | `TSTEP` 10 il, aylıq, cəmi 3650 gün | ✅ |

Qeyd: WOC (8450) layın altında (8425), GOC (8300) üstündədir — lay tamamilə
neftlidir, ilkin sərbəst qaz yoxdur. BHP istinad dərinlikləri perforasiya
edilmiş hüceyrələrin mərkəzi ilə dəqiq üst-üstə düşür (8325+20+30+25 = 8400,
8325+10 = 8335) — dərinlik düzəlişi lazım deyil.

Cədvəllər (SWOF 15 sətir, SGOF 15 sətir, PVDG 10 sətir, PVTO 7 doymuş + 2
doymamış sətir) deck-də olduğu kimidir; model qurularkən birbaşa oradan
köçürüləcək.

## 3 · Etalon nəticə — OPM Flow, SPE1CASE2 (öz oxuyucumuzla oxundu)

| t, gün | FOPR, STB/gün | FGOR, Mscf/STB | WBHP PROD, psia | WBHP INJ, psia | BPR (1,1,1), psia |
|---|---|---|---|---|---|
| 1 | 20000 | 1.27 | 2905.09 | 8081.86 | 5063.13 |
| 304 | 20000 | 1.28182 | 2653.21 | 6429.19 | 6146.66 |
| 1034 | 20000 | 1.27988 | 4013.49 | 7458.84 | 7207.69 |
| 1399 | 20000 | 5.9606 | 2459.59 | 7632.99 | 7388.97 |
| 1580 | 18865.9 | 8.41719 | 1000 | 6982.48 | — |
| 2129 | 11624.9 | 11.4945 | 1000 | 5542.88 | 5318.34 |
| 3650 | 5732.65 | 22.1403 | 1000 | 4333.43 | 4101.32 |

* FOPR ilk dəfə 1550-ci gündə 20 000-dən aşağı düşür — istismarçı BHP limitinə keçir.
* Diapazonlar: BPR(1,1,1) 4101.32 – 7534.16 psia; BPR(10,10,3) 3278.11 – 6233.51 psia;
  WBHP INJ 4333.43 – 8081.86 psia. Hesabatda 123 addım.

## 4 · Mühərrik ilə fərqlər (gap cədvəli)

| # | SPE1 tələbi | Mühərrikdə indi | Qiymət |
|---|---|---|---|
| G1 | PVT cədvəlləri PVTO/PVDG/PVTW (FIELD) | ✅ `io/pvt_io.py` — itkisiz oxuma + açıq birləşdirmə (Seans 34, Q-27) | bağlandı; çox qollu Bo/μo — G2/G3-də bağlandı |
| G2 | Doymamış neft özlülüyü μo(p, Rs) — PVTO-da 0.51 → 0.74 cP | ✅ provider (Seans 35) + mühərrik (Seans 36, Q-29) | bağlandı; 3-cü sütun 3.6×10⁻¹¹, lövbər qüsuru da düzəldildi |
| G3 | Doymamış Bo hər Rs üçün öz sıxılması ilə | ✅ `BlackOilPVTProvider(..., oil_branches=...)` — c_o(Rs), n(Rs) qollardan (Seans 37, Q-30) | bağlandı; qollar arasında 11 % fərq ölçüldü; modelə qoşulma ⏳ SPE1CASE2 ilə |
| G4 | SGOF cədvəli (krg, krog) | ✅ `GasSaturationTable` + `read_sgof` — Seans 33, Q-26 | bağlandı (SWOF-un səssiz atılması da düzəldildi) |
| G5 | Neftin SƏTH debiti (`ORAT`), qazın SƏTH vurma debiti (`RATE`) | ✅ `RateBasis.SURFACE` — Seans 30, Q-24 | bağlandı |
| G6 | Süxur sıxılmasının istinad təzyiqi 14.7 psia | ✅ `rock.compressibility_reference_pressure` — Seans 32, Q-25 | bağlandı (ölçüldü: 4800 psia-da 1.44 % fərq) |
| G7 | Doymuş qol 5014.7 psia-dan yuxarı (etalonda hüceyrə təzyiqi 7534 psia-ya qalxır) | `np.interp` sərhəddə saxlayır (Rs_sat = 1.618 plato) | ✅ **Seans 41 (Q-32)** — OPM-in qaydası (xətti uzantı) tətbiq olundu; ölçülmüş təsir §8 |
| G9 | Canlı neftin sıxlığı `(ρo + Rs·ρg)/Bo` (cazibə, ilkin tarazlıq) | ✅ Seans 38, Q-31 — əvvəl `ρo/Bo` idi | bağlandı; SPE1-də +27 %; THP hidravlikası ⏳ |
| G8 | `DRSDT 0` (yalnız CASE1) | dəstəklənmir — qaz həmişə yenidən həll olur | CASE2 hədəf seçildi (Q-23) |
| — | Üç fazalı kro | Stone II. SPE1-də su hərəkətsizdir (Sw = Swc = 0.12); orada Stone II və Eclipse defolt modeli eyni `kro = krog` verir | fərq gözlənilmir |
| — | Etalonun oxunması | `resdata` qurulmayıb | öz oxuyucumuz işləyir (Seans 29) |

## 5 · İş planı (hər biri ayrıca commit, sonlu fərq testi ilə)

1. **G5 — səth debiti hədəfi** (istismarçıda neft, vurucuda qaz). ✅ Seans 30 (Q-24).
2. **G6 — süxur sıxılmasının istinad təzyiqi** ✅ Seans 32 (Q-25) — mövcud
   modellər bit-bit eyni qaldı; flüidin öz istinadı toxunulmadı.
3. **G4 — SGOF cədvəli** ilə qaz relperm provider-i. ✅ Seans 33 (Q-26).
4. **G1** ✅ Seans 34 (Q-27) — PVTO/PVDG/PVTW oxuyucuları.
   **G2** ✅ Seans 35 (provider) + Seans 36 (mühərrik + lövbər), Q-28/Q-29.
   **G3** ✅ Seans 37 (Q-30) — c_o və n hər qoldan, Rs üzrə interpolyasiya.
   Ölçüldü (real deck): c_o 2.056e-4 ↔ 1.832e-4 1/bar (~11 %), n 0.460 ↔
   0.580. Əvvəl yazılmış "0.3 %, 0.46/0.51" sınaq deck-indəki səhv
   sətirdən gəlirdi. Korrelyasiyadakı 0.278 SPE1 üçün YANLIŞDIR.
5. **SPE1CASE2 modeli** ✅ qurulub (Seans 38, `imex2d/benchmarks/spe1.py`,
   `tools/spe1_compare.py`). İlk müqayisə etalondan fərqlənir — bax §6.
   `tests/golden/` reqressiyası fərqlər izah olunandan SONRA (Q-31 Qərar 7).
6. ⏳ G7 və CASE1 (`DRSDT 0`) — müqayisə nəticəsinə görə.


## 6 · İlk müqayisə — Seans 38 (G1–G6, G9 bağlı; G7 açıq)

`python tools/spe1_compare.py SPE1CASE2.DATA SPE1CASE2` — 514 addım, ~80 san.

| t, gün | kəmiyyət | bizdə | OPM Flow | fərq |
|---|---|---|---|---|
| 1 | WBHP INJ, psia | 5271 | 8082 | −34.8 % |
| 304 | FGOR, Mscf/STB | 1.359 | 1.282 | +6.0 % |
| 1034 | FGOR | 6.255 | 1.280 | +389 % |
| 1034 | WBHP PROD | 1576 | 4013 | −60.7 % |
| 1034 | BPR (10,10,3) | 4821 | 5806 | −17.0 % |
| 1399 | FOPR, STB/gün | 14 720 | 20 000 | −26.4 % |
| 1399 | BPR (1,1,1) | 5842 | 7389 | −20.9 % |
| 2129 | FOPR | 9735 | 11 625 | −16.3 % |
| 3650 | FOPR | 4980 | 5733 | −13.1 % |
| 3650 | FGOR | 24.55 | 22.14 | +10.9 % |
| 3650 | WBHP INJ | 4142 | 4333 | −4.4 % |
| 3650 | BPR (1,1,1) / (10,10,3) | 3931 / 3158 | 4101 / 3278 | −4.2 / −3.7 % |

| Hadisə (0.5 % toleransla) | bizdə | OPM Flow |
|---|---|---|
| FOPR 19 900-dən aşağı (istismarçı BHP limitinə keçir) | 1120 gün | 1550 gün |
| FGOR > 2 Mscf/STB (qazın çatması) | 900 gün | 1276 gün |

Fərqin səbəbi **Seans 40-da ölçülərək tapıldı — §7**. Qısaca: G7
(Rs_sat platosu). Vurucunun 1-ci gündəki BHP fərqi (5271 ↔ 8082) AYRI
məsələdir və cəbhənin sürətini izah etmir — §7.2.


## 7 · Fərqin səbəbi — Seans 40 (ölçülmüş)

### 7.1 · Əsas səbəb: G7 (doymuş Rs qolunun platosu)

OPM mənbəsindən oxundu: `LiveOilPvt.hpp:512` cədvəli `extrapolate=true`
ilə çağırır, `Tabulated1DFunction.hpp:266–283` isə bu halda son seqmentin
xəttini cədvəldən kənarda davam etdirir. Yəni OPM-in nefti 5014.7 psia-dan
yuxarı daha çox qaz həll edir:

| p, psia | Rs_sat OPM | bizdə | fərq |
|---|---|---|---|
| 6 150 | 2.013 | 1.618 | +24 % |
| 7 459 | 2.469 | 1.618 | +53 % |

Ölçülmüş izi: eyni təzyiqdə bizdə sərbəst qaz çoxdur (304-cü gün, blok
201: **0.195 ↔ 0.169**), ona görə cəbhə hər istiqamətdə tez gedir
(blok 10: **944 ↔ 1611 gün**) və istismarçı sütununa ~365 gün tez çatır.
Qalan bütün fərqlər (kro-nun çökməsi, BHP limitinə tez keçid, təzyiqin
ayrılması, FGOR) bunun nəticəsidir.

### 7.2 · İkinci, AYRI məsələ: vurucu bağlantısının mobilliyi

Vaxt imzası göstərir ki, bu, cəbhə ilə bağlı deyil:

| gün | WBHP INJ fərqi |
|---|---|
| 1 | **−34.8 %** (ΔP bizdə 449, OPM 3019 psi) |
| 304 | −1.2 % |

Bizim qayda vurulan fazanın SON NÖQTƏ mobilliyini işlədir, yəni bloku
əvvəlcədən qazla dolmuş sayır. Blok həqiqətən qazla dolanda (≈300 gün)
fərq itir. Debit sabit olduğu üçün bu, cəbhənin sürətinə təsir etmir.
OPM-in vurucu bağlantısındakı qaydası hələ mənbədən oxunmayıb ⏳.

### 7.3 · Ölçülərək TƏKZİB olunanlar

* Bg, μg, Bo, μo, Rs_sat, c_o — deck ilə **0.00 %** fərq.
* Quyu indeksi — analitik Peaceman ilə eyni (10.610 ↔ 10.608).
* `PERMZ`, `PORO` — deck ilə eyni.
* İstismarçının blokundakı erkən qaz ayrılması — **etalonda da var**
  (`BGSAT:300`: 0.0131 @ 31 g, 0.0222 @ 59 g, 0 @ 212 g).

### 7.4 · Kiçik, ayrıca fərq

Doymamış qolda OPM düyünlər arasında xəttidir, biz üstəl qanunla gedirik:
Bo-da 0.06 %, **μo-da 1.97 %** (5500–6500 psia).


## 8 · G7-dən sonrakı müqayisə — Seans 41

Eyni qaçış, eyni etalon; yeganə dəyişiklik doymuş qolun xətti uzadılmasıdır.

| Kəmiyyət | G7-dən əvvəl | **sonra** | OPM Flow |
|---|---|---|---|
| WBHP PROD, 1034 | 1576 (−60.7 %) | **3965 (−1.2 %)** | 4013 |
| FGOR, 1034 | 6.255 (+389 %) | **1.776 (+38.7 %)** | 1.280 |
| FOPR, 1399 | 14 720 (−26.4 %) | **19 450 (−2.8 %)** | 20 000 |
| FOPR, 3650 | 4980 (−13.1 %) | **5176 (−9.7 %)** | 5733 |
| FGOR, 3650 | 24.55 (+10.9 %) | **24.29 (+9.7 %)** | 22.14 |
| BHP limitinə keçid | 1120 gün | **1381 gün** | 1550 gün |

Qaz cəbhəsi (Sg > 0.05):

| blok | əvvəl | **sonra** | OPM |
|---|---|---|---|
| 10 (10,1,1) | 944 | **1533** | 1611 |
| 200 (10,10,2) | — | **1132** | 1246 |
| 300 (10,10,3) | — | **1157** | 1307 |

**Qalan fərqlər:** vurucunun 1-ci gündəki BHP-si (−34.8 %, §7.2 — ayrı
məsələ), cəbhənin 110–150 gün tez gəlməsi, FGOR-un 1034-cü gündə +38.7 %
olması. Etalon (golden) fayl hələ YAZILMAYIB.
