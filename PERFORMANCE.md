# Performans ölçməsi (C3)

Ölçmə aləti: `python tools/benchmark.py`

```bash
python tools/benchmark.py --sizes 21,31,41 --days 300
python tools/benchmark.py --profile 41            # cProfile hesabatı
python tools/benchmark.py --sizes 21,31 --nz 4    # 3D
```

## 1. Profil — vaxt harada gedir

41×41, 300 gün, optimallaşdırmadan ƏVVƏL:

| Funksiya | Pay |
|---|---|
| `SuperLU.solve` (ILU tətbiqi, KQ daxilində) | 22 % |
| `gstrf` (ILU faktorlaşdırması) | 13 % |
| `coo_tocsr` + `csr_sort_indices` + `sum_duplicates` | **12 %** |
| `_update_saturation` | 7 % |
| `_solve_pressure` (matris yığımı) | 7 % |
| KQ döngəsi və qalanı | 39 % |

Nəticə: vaxtın **35 %-i xətti həlledicidə**, **12 %-i isə hər addımda
matrisin sıfırdan qurulmasında** gedirdi.

## 2. Tətbiq edilən optimallaşdırma

### 2.1 CSR strukturunun keşlənməsi

Matrisin seyrəklik strukturu simulyasiya boyunca **dəyişmir** — yalnız
dəyərlər dəyişir. Əvvəl hər addımda COO matris yaradılıb CSR-ə
çevrilirdi (sıralama + dublikatların toplanması daxil).

İndi struktur bir dəfə qurulur, COO girişlərinin CSR `data`
massivindəki mövqeyi (`_data_index`) yadda saxlanılır, hər addımda
yalnız `np.bincount` ilə toplama qalır.

| Ölçü | Əvvəl | Sonra | Qazanc |
|---|---|---|---|
| 21×21 | 2.79 ms/addım | 1.82 ms/addım | **35 %** |
| 41×41 | 5.29 ms/addım | 3.78 ms/addım | **29 %** |
| Test dəsti (173 test) | 146 san | 66 san | **55 %** |

### 2.2 Ön-şərtçinin yenilənmə tezliyi

`preconditioner_refresh_steps`: 25 → **50**.

| Tezlik | ms/addım (41×41) |
|---|---|
| 25 | 3.88 |
| **50** | **3.52** |
| 100 | 3.53 |
| 400 | 4.66 |
| heç vaxt | 7.71 |

### 2.3 Yoxlanılıb və RƏDD edilib: zəif ILU

`ilu_drop_tolerance` böyüdüləndə faktorlaşdırma ucuzlaşır, amma KQ
iterasiyaları partlayır:

| drop_tol | fill | ms/addım |
|---|---|---|
| 1e-4 | 10 | 4.00 |
| 1e-3 | 5 | 120.2 |
| 1e-2 | 3 | 109.4 |

Güclü ILU bu məsələ üçün məcburidir — parametr azaldılmamalıdır.

## 3. Miqyaslanma

### 3.1 Grid ölçüsü üzrə (200 gün)

| Grid | Hüceyrə | Addım | ms/addım | µs/hüceyrə·addım |
|---|---|---|---|---|
| 21×21 | 441 | 641 | 1.79 | 4.06 |
| 31×31 | 961 | 582 | 2.46 | 2.56 |
| 41×41 | 1 681 | 546 | 3.29 | 1.96 |
| 51×51 | 2 601 | 521 | 6.06 | 2.33 |

Hüceyrə sayı **5.9 dəfə** artanda ümumi vaxt **2.7 dəfə** artır —
həlledici yaxşı miqyaslanır. Hüceyrəyə düşən xərc azalır, çünki
addım başına sabit Python yükü amortizasiya olunur.

### 3.2 Müddət üzrə (41×41) — ƏSAS MƏHDUDİYYƏT

| Müddət | Addım | Addım/gün |
|---|---|---|
| 300 gün | 826 | 2.75 |
| 1 000 gün | 2 840 | 2.84 |
| 3 000 gün | 8 830 | **2.94** |

Addım sayı müddətlə **xətti** artır və orta zaman addımı ≈ **0.34 gün**
səviyyəsində qalır. Bu, IMPES-in CFL şərtidir və optimallaşdırma ilə
aradan qalxmır.

## 4. A6 (fully implicit) üçün əsaslandırma

Real işlək model — 100×100×20 = **200 000 hüceyrə**, 30 illik proqnoz
(10 950 gün) — üçün ekstrapolyasiya:

| | IMPES (hazırkı) | Fully implicit (gözlənilən) |
|---|---|---|
| Orta Δt | 0.34 gün | 20–30 gün |
| Addım sayı | ~32 000 | ~400 |
| Addım xərci | ~0.46 san | ~2.5 san (Nyuton iterasiyaları) |
| **Ümumi** | **~4 saat** | **~17 dəqiqə** |

Yəni A6-nın gözlənilən qazancı təxminən **15 dəfədir** və o, hüceyrə
sayından deyil, **zaman addımının uzunluğundan** gəlir.

Vacib qeyd: kiçik modellərdə (< 5 000 hüceyrə, < 1 000 gün) IMPES
daha sürətli qalacaq, çünki hər addımı ucuzdur və Nyuton iterasiyası
yoxdur. Ona görə hər iki mühərrik `ISimulationEngine` interfeysi
altında saxlanılmalıdır.

## 5. A6 nəticəsi — ölçülmüş qazanc

Fully implicit mühərrik (A6, mərhələ 1-4) tamamlandıqdan sonra
`tools/benchmark.py --engine both` ilə ölçülən:

| Grid | Müddət | IMPES addım/vaxt | Implicit addım/vaxt | Sürət |
|---|---|---|---|---|
| 21×21 | 900 gün | 2 978 / 4.55 san | 38 / 0.38 san | **12x** |
| 31×31 | 900 gün | 2 734 / 7.30 san | 38 / 0.73 san | **10x** |
| 41×41 | 3 600 gün | 3 720 / 21.1 san | 127 / 2.85 san | **7.4x** |

Diqqətəlayiq: implicit addımın özü **6-13 dəfə bahalıdır**
(10-19 ms vs 1.5-2.7 ms), çünki hər addımda 2-4 Nyuton iterasiyası
və hər iterasiyada Jakobian + xətti həll var. Qazanc tamamilə
**addım sayından** gəlir: 2 978 → 38.

Orta zaman addımı: IMPES 0.30 gün, implicit **24 gün** — 80 dəfə fərq.

Nəticə fərqi 0.8 % səviyyəsindədir (RF 14.22 % vs 14.10 %) və bu,
iki diskretizasiyanın təbii fərqidir, xəta deyil.

### Hansı mühərriki seçmək

| Hal | Tövsiyə |
|---|---|
| Kiçik model, qısa proqnoz, sürətli sınaq | IMPES |
| Uzun proqnoz (illər) | Fully implicit |
| Böyük grid | Fully implicit |
| Cəbhənin dəqiq yeri vacibdir | IMPES (kiçik Δt ilə) |

Hər ikisi `ISimulationEngine` interfeysi altındadır və interfeysdə
"Hesablama sxemi" siyahısından seçilir.

### CPR ön-şərtçisi (A6/5)

Xətti həlledicidə ikinci seçim. Dürüst nəticə: **sürətdə ILU-dan
üstün deyil**, çünki orta ölçülü sistemlərdə güclü ILU onsuz da
1 iterasiyada həll edir. Üstünlük **yaddaşdadır** — 3 025 hüceyrədə
CPR ILU-nun 44 %-i qədər yer tutur və nisbət ölçü ilə yaxşılaşır.

Praktik məna: 100 000+ hüceyrəli modellərdə ILU yaddaşa sığmır,
CPR isə işləyir. Ətraflı: `A6_PLAN.md`, mərhələ 5.

## 6. Növbəti optimallaşdırma namizədləri

| Namizəd | Gözlənilən qazanc | Qeyd |
|---|---|---|
| CPR ön-şərtçisi | böyük (implicit üçün) | Təzyiq və doyumluluq bloklarını ayırır |
| `_update_saturation`-da Numba | 5–7 % | `np.add.at` çağırışları |
| Nisbi keçiricilik keşi | 2–3 % | Corey düsturları hər addımda hesablanır |

`np.add.at` və Corey hesablamaları ümumi vaxtın yalnız ~10 %-ni
tutur — Numba-ya keçid A6-dan ƏVVƏL sərfəli deyil.
