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
| G1 | PVT cədvəlləri PVTO/PVDG/PVTW (FIELD) | `PVTTable` tək təzyiq şəbəkəsidir; PVDG (9014.7-yə qədər) və PVTO doymuş qolu (5014.7-yə qədər) fərqli təzyiqlərdədir | köçürücü lazımdır |
| G2 | Doymamış neft özlülüyü μo(p, Rs) — PVTO-da 0.51 → 0.74 cP | `mu_o = pvt.oil_viscosity(p)` — yalnız doymuş əyri (`three_phase_newton.py:145`) | **fizika boşluğu** (qalıq + Jakobian) |
| G3 | Doymamış Bo hər Rs üçün öz sıxılması ilə | tək `c_o`, cədvəlin Pb-dən yuxarı hissəsindən (`black_oil.py::_build_undersaturated_branch`) | təqribi — ölçülməlidir |
| G4 | SGOF cədvəli (krg, krog) | qaz əyrisi yalnız Corey (`GasCoreyParameters`) | **provider boşluğu** |
| G5 | Neftin SƏTH debiti (`ORAT`), qazın SƏTH vurma debiti (`RATE`) | ✅ `RateBasis.SURFACE` — Seans 30, Q-24 | bağlandı |
| G6 | Süxur sıxılmasının istinad təzyiqi 14.7 psia | ✅ `rock.compressibility_reference_pressure` — Seans 32, Q-25 | bağlandı (ölçüldü: 4800 psia-da 1.44 % fərq) |
| G7 | Doymuş qol 5014.7 psia-dan yuxarı (etalonda hüceyrə təzyiqi 7534 psia-ya qalxır) | `np.interp` sərhəddə saxlayır (Rs_sat = 1.618 plato) | OPM-in ekstrapolyasiya qaydası mənbədən yoxlanılmalıdır ⏳ |
| G8 | `DRSDT 0` (yalnız CASE1) | dəstəklənmir — qaz həmişə yenidən həll olur | CASE2 hədəf seçildi (Q-23) |
| — | Üç fazalı kro | Stone II. SPE1-də su hərəkətsizdir (Sw = Swc = 0.12); orada Stone II və Eclipse defolt modeli eyni `kro = krog` verir | fərq gözlənilmir |
| — | Etalonun oxunması | `resdata` qurulmayıb | öz oxuyucumuz işləyir (Seans 29) |

## 5 · İş planı (hər biri ayrıca commit, sonlu fərq testi ilə)

1. **G5 — səth debiti hədəfi** (istismarçıda neft, vurucuda qaz). ✅ Seans 30 (Q-24).
2. **G6 — süxur sıxılmasının istinad təzyiqi** ✅ Seans 32 (Q-25) — mövcud
   modellər bit-bit eyni qaldı; flüidin öz istinadı toxunulmadı.
3. **G4 — SGOF cədvəli** ilə qaz relperm provider-i.
4. **G1 + G2 (+G3 ölçmə)** — PVTO/PVDG köçürücüsü və doymamış özlülük.
5. **SPE1CASE2 modeli** + etalonla müqayisə (FOPR, FGOR, WBHP, BPR), `tests/golden/`-da reqressiya.
6. ⏳ G7 və CASE1 (`DRSDT 0`) — müqayisə nəticəsinə görə.
