# A7 — qaz fazası kodunun bərpası

**Bərpa tarixi:** 10 sentyabr 2026

## Haradan gəldi

| | |
|---|---|
| Repo | `github.com/eoraieor-tech/lay-simulyatoru` (public) |
| Commit | `66deca8295617244761092945f42fc896c2aa2aa` |
| Commit mesajı | `v67 - son 3 fazali versiya (qaz fazasi + matplotlib 3D)` |
| Tarix | 27 avqust 2026 |

Bu, reponun **ilk commit-idir** — qaz fazası silinməzdən əvvəl qəsdən
saxlanılmış snapshot.

## Niyə burada saxlanılır

Kod `v69`-da əsas kod bazasından tamamilə çıxarılıb. Yüklənmiş 14 zip
arxivinin **heç birində yoxdur** (silinmə mərhələli olub, bütün
arxivlər ondan sonrakıdır). Yeganə mənbə git tarixçəsidir.

Repo private edilsə və ya silinsə, bu kod itərdi. Ona görə burada
nüsxə saxlanılır.

## Nə var

**Mənbə kodu — 2 165 sətir**

| Fayl | Sətir | Nə edir |
|---|---|---|
| `three_phase_residual.py` | 1 066 | Qalıq tənlikləri + analitik Jakobian |
| `three_phase_newton.py` | 385 | Nyuton-Rafson döngəsi (3×3 blok) |
| `three_phase_engine.py` | 261 | Mühərrik sarğısı |
| `stone_relperm.py` | 172 | **Stone II** 3-fazalı nisbi keçiricilik |
| `three_phase_state.py` | 170 | Primary dəyişənlər + dəyişən keçidi |
| `three_phase.py` | 111 | Üç fazalı doyumluluq (domain) |

**Testlər — 2 158 sətir**

`test_three_phase_residual.py` (964), `test_three_phase_newton.py` (257),
`test_three_phase.py` (215), `test_variable_switching.py` (205),
`test_gas_pvt.py` (192), `test_stone_relperm.py` (192),
`test_gas_ui_wiring.py` (133).

(İlk bərpada `test_variable_switching.py` buraxılmışdı — sonradan əlavə
olundu; cəmi 13 fayl, 4 323 sətir. Bax `AUDIT_2026-09-10.md` §11.)

## Vəziyyət

`A7_PLAN.md`-ə görə bütün mərhələlər (qaz PVT → Stone II → dəyişən
keçidi → qalıq → Jakobian → Nyuton → UI) **HAZIR** idi.

Həll olunmamış **yeganə** problem: quyu öz BHP həddinə yaxınlaşanda
Nyutonun yığılmaması (t≈6.7 gündə istismarçının qonşu hüceyrəsində
neft residualı DÖVR-2 rəqs edir).

Bu, tanınmış problemdir; həll yolları: Appleyard chopping, doyma
dəyişiminin məhdudlaşdırılması, trust-region Nyuton (Wang & Tchelepi),
per-cell qəbul meyarı.

**Yenilənmə (1 oktyabr 2026):** yuxarıdakı vəziyyət bərpa GÜNÜNƏ aiddir.
Sahibkar «öz fizikamız» seçdi və kod həmin gün əsas kod bazasına (`imex2d/`,
`tests/`) qaytarıldı (B2 — `ISH_HESABATI.md` → Seans 6). Yığılmama problemi
B3-A və B3-B-də həll olundu (Seans 7 və 10, `QARARLAR.md` → Q-09).

Bu qovluq indi yalnız **ARXİVDİR**: canlı kod onu idxal etmir, `pytest.ini`-nin
`testpaths`-ına daxil deyil. Canlı fayllar o vaxtdan xeyli dəyişib (məs.
`three_phase_residual.py` 1 066 → 1 401 sətir) — buradakı nüsxələr v67
vəziyyətidir və işlək kod kimi götürülməməlidir.

Ətraflı: [`../../AUDIT_2026-09-10.md`](../../AUDIT_2026-09-10.md) §8.2
