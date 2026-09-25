# IMEX-2D — Analizin qısa xülasəsi

**Tarix:** 26 sentyabr 2026 · **Commit:** `eacb029` · **Tam hesabat:** [PROJECT_ANALYSIS.md](PROJECT_ANALYSIS.md) · **Diaqramlar:** [docs/diagrams/](docs/diagrams/README.md)

## Nədir

Neft-qaz yatağının **masaüstü simulyatoru** (Python 3.14 + PyQt5, Windows).
Quyu ölçmələrindən 3D geoloji model qurur (Kriging / SGS / SIS), neft-su-qaz
axınını zamana görə hesablayır (black-oil; IMPES və tam implicit Nyuton;
TPFA və MPFA-O), nəticəni qrafik, 3D, günlük cədvəl, CSV/JSON/PDF kimi verir.
History matching və həssaslıq analizi var.

**Yoxdur:** süni intellekt/LLM/ML, verilənlər bazası, şəbəkə/API, autentifikasiya,
Docker, CI. «Model» sözü burada fiziki (Darcy, PVT korrelyasiyaları, Corey,
Stone II, Peaceman) və statistik (variogram, kriging, optimallaşdırma) modelləri bildirir.

## Rəqəmlər

| | |
|---|---|
| Əsas kod (`imex2d/`) | 145 modul, 46 240 sətir, 12 paket |
| Testlər | 137 fayl, 2 762 test |
| **Test nəticəsi (bu analizdə)** | **2761 keçdi · 1 atlandı (`resdata` yoxdur) · 1 xfail (bilinən) · 0 uğursuz** — 8 dəq 47 san |
| Sənədlər | 39 `.md` faylı, 17 000+ sətir |
| Asılılıqlar | numpy 2.5.2, scipy 1.18.1, PyQt5 5.15.11, matplotlib 3.11.1, vtk 9.7.0 (dəqiq sabitlənib) |

## Arxitektura bir cümlə ilə

Qatlı monolit: `domain → interfaces → simulation/geology/discretization →
application → rendering/reporting → ui`, asılılıqlar `app.py`-də inyeksiya olunur.
`domain` tam müstəqildir; iki pozuntu var: `simulation → application.config`
və `ui → simulation` (mühərrik sinfi UI-də seçilir).

## Ən vacib problemlər

| ID | Səviyyə | Problem | Yer |
|---|---|---|---|
| P-01 | High | `MainWindow` 3 155 sətir, pytest-də qurulmur — UI orkestrasiyası avtomatik test olunmur | `ui/main_window.py` |
| P-02 | High | Layihə faylı atomik yazılmır — çökmə `.imx`-i korlaya bilər | `application/serialization.py:375` |
| P-03 | High | Üç fazalı Jakobianda bilinən xətalar (TB-1: 0.49, TB-2: 0.19) | `implicit/three_phase_residual.py` |
| P-04 | High | SPE1CASE2-də OPM Flow-dan ~145 günlük qaz cəbhəsi fərqi | `SPE1.md` §8 |
| P-05 | Medium | `LinearSolverConfig` tam implicit mühərriklərdə səssizcə işlədilmir | `implicit/engine.py:60,84-86` |
| P-06 | Medium | Mühərrik hər qaçışda iki dəfə qurulur (biri UI axınında) | `main_window.py:2139` |
| P-08 | Medium | İki FIM mühərriki arasında quyu-dövrəsi kodunun təkrarı | `engine.py` / `three_phase_engine.py` |
| P-10 | Medium | `pyproject.toml`, CI, coverage, zəiflik yoxlaması yoxdur | kök |

Critical problem tapılmadı. Təhlükəsizlik səthi kiçikdir (oflayn, `pickle`/`eval`/`subprocess` yoxdur);
əsas risk məlumat bütövlüyüdür (P-02).

## İlk 5 öyrənmə addımı

1. `README.md` + `ARCHITECTURE.md` §2–4 + [diaqram 1–2](docs/diagrams/README.md).
2. `app.py` → `application/config.py` → `domain/reservoir_model.py`.
3. `ARCHITECTURE.md` §7-dəki skripti işlət (UI-siz simulyasiya), RF-i çap et.
4. `simulation/impes_engine.py`, sonra `implicit/{residual,newton,time_stepping,engine}.py` ([diaqram 4](docs/diagrams/04_zaman_addimi_nyuton.md)).
5. `application/simulation_service.py` — provider və mühərrik seçimi; sonra `tests/test_analytical_bl.py` ilə doğrulamanı gör.
</content>
</invoke>
