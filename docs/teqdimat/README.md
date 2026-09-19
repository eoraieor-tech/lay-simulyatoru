# Təqdimat materialları — proqramın iş prinsipi

**Hazırlanıb:** 19–20 sentyabr 2026 (Seans 44) · **Kod commit-i:** `dcaa74d`
(materiallar hazırlananda `imex2d/` bu vəziyyətdə idi)

Sahibkarın tələbi: proqramın necə işlədiyini təqdim etmək üçün bir neçə izah
yolu. Yeddi variant təklif olundu, hamısı hazırlandı.

## 1 · Materiallar

| # | Material | Link (claude.ai artifakt) | Mənbə bu qovluqda |
|---|---|---|---|
| — | İnteraktiv xəritə: iş axını + modul qrafiki | https://claude.ai/artifact/8DmUT7Qwhu1Zfi3dk6AifE | `xerite.html` |
| 1, 3, 4 | Slayd dəsti, 18 slayd (nümunə hekayəsi, doğrulama, müqayisə) | https://claude.ai/artifact/8GC5g2DRyDasc4eWUdXR8T | `slaydlar/` |
| 2 | «Zaman addımının içi» — Nyuton döngüsünün animasiyası | https://claude.ai/artifact/EPBLhf7YB6Tx1Dkvr97oFL | `zaman-addimi.html` |
| 5 | Demo planı: 7 dəqiqəlik canlı ssenari, ehtiyat planı, suallar | https://claude.ai/artifact/KEaSKHUKS2D21dRqV64Rd5 | `demo-plani.html` |
| 6, 7 | İki səhifəlik icmal + terminlər lüğəti | https://claude.ai/artifact/3xGmWTWXQ3ttYsUvAJXdme | `icmal.html`, `IMEX-2D_icmal.pdf` |

⚠️ Artifaktlar **şəxsidir** — başqaları linki yalnız səhifənin «Share» menyusundan
paylaşıldıqdan sonra aça bilər.

`slaydlar/slides/*.html` fayllarındakı şəkillər (`/_blob/…`) artifaktın öz
fayl anbarındadır — bu HTML-lər artifaktdan kənarda şəkilsiz görünür. Şəkilləri
yenidən yaratmaq üçün bax §3.

## 2 · Ölçülmüş rəqəmlər (slaydlarda və səhifələrdə işlədilənlər)

Heç bir rəqəm əl ilə yazılmayıb — hamısı aşağıdakı qaçışlardan və layihə
sənədlərindən götürülüb.

**Nümunə qaçışı** — `numune_quyular.csv` (5 quyu × 3 qat) → adi Kriging →
41 × 41 × 3 grid (hüceyrə 20 × 20 × 5 m) → five-spot (INJ-1 BHP 320 bar,
PROD-1 BHP 150 bar) → `FullyImplicitEngine`, TPFA, 1500 gün:

| Kəmiyyət | Dəyər |
|---|---|
| Keçiricilik (Kriging nəticəsi) | 178–1119 mD, geologiya 0.18 s |
| İlkin neft ehtiyatı (OOIP) | 896 963 m³ |
| Neftvermə əmsalı, 1500-cü gün | **57.80 %** |
| Suyun istismarçıya çatması | 557-ci gün |
| Sulaşma, 1500-cü gün | 93.6 % |
| Zaman addımı / Nyuton iterasiyası | 88 addım, 268 iterasiya, orta Δt 17.0 gün, 2 təkrar |
| Nyuton cəhdləri (rədd olunanlar daxil) | 90 — ilk iki cəhd yığılıb, amma ΔSw > 0.2 olduğu üçün rədd edilib |
| Hesablama vaxtı | 39 s (birinci qaçış), 27 s (təkrar qaçış) — maşının yükündən asılıdır |

**Buckley-Leverett** — 1D, 120 hüceyrə, 250 gün, IMPES; metrikalar
`MainWindow.run_validation` ilə eynidir:

| Metrika | Dəyər |
|---|---|
| RMS xəta | 0.0231 Sw |
| Cəbhə (orta nöqtə) | 190.9 m, analitik 182.6 m → **4.53 %** |
| Həcm balansı | 2.71 % |

**SPE1CASE2** — yenidən işlədilmədi; rəqəmlər `SPE1.md` §8.3-dən (Seans 43).

**Kod bazası** — `tools/module_graph.py`: 130 modul, 491 idxal əlaqəsi,
45 098 sətir.

## 3 · Yenidən qurmaq

```bash
PY tools/presentation_demo.py --out <qovluq>   # demo.json, bl_metrics.json, figs/*.png (~1 dəq)
PY tools/module_graph.py --out <qovluq>/graph.json
```

`presentation_demo.py` Seans 44-də iki dəfə işlədildi, vaxtdan başqa bütün
rəqəmlər eyni çıxdı.

## 4 · Açıq qalanlar

- ⏳ Demo planının 4-cü addımı: standart açılış modelində (41 × 41 × 1)
  «MODELİ İŞƏ SAL»-ın vaxtı ölçülməyib — demo planında boş sahə var,
  sahibkar öz maşınında doldurur.
- ⏳ Standart açılışda geologiya cədvəlində quyu xassələri yoxdur (yalnız
  mövqelər). Kriging-i canlı göstərmək üçün xassələri əvvəlcədən doldurub
  layihə faylı kimi saxlamaq lazımdır.
- Slayd 15-dəki kommersiya sütunu (Eclipse, CMG, tNavigator) ümumi imkanları
  göstərir; konkret versiya və lisenziya ilə yoxlanılmayıb.
