# Diaqram 1 — Qatlı arxitektura və FAKTİKİ asılılıqlar

**Mənbə:** `ARCHITECTURE.md` §2 (nəzərdə tutulan qaydalar) və 2026-09-26
tarixində `imex2d/` paketinin bütün `import` sətirlərinin AST ilə sayılması
(bax `PROJECT_ANALYSIS.md` §4.3). Oxların üstündəki rəqəm — neçə import
sətri olduğu.

Qırmızı oxlar sənədləşdirilmiş istiqaməti («həmişə içəriyə doğru») POZUR.

```mermaid
flowchart TB
    app["app.py<br/>(composition root)"]

    subgraph UI["ui/ — PyQt5"]
        mw["main_window.py<br/>MainWindow (3155 sətir)"]
        panels["panels.py<br/>10 panel"]
        worker["worker.py<br/>QThread işçiləri"]
    end

    subgraph OUT["Çıxış qatları (Qt-siz)"]
        rendering["rendering/<br/>matplotlib + VTK"]
        reporting["reporting/<br/>PDF, CSV/JSON, günlük"]
    end

    subgraph APP["application/"]
        svc["simulation_service.py"]
        geosvc["geology_service.py"]
        ser["serialization.py (.imx)"]
        cfg["config.py"]
    end

    history["history/<br/>uyğunlaşdırma, həssaslıq"]
    geology["geology/<br/>variogram, kriging, SGS, SIS"]
    io["io/<br/>GRDECL, PVT/SCAL deck, Eclipse, OPM"]
    benchmarks["benchmarks/<br/>SPE1"]

    subgraph SIM["simulation/"]
        engines["impes_engine.py<br/>implicit/engine.py<br/>implicit/three_phase_engine.py"]
        physics["pvt/, scal, capillary,<br/>well_model, wellbore/"]
    end

    disc["discretization/<br/>MPFA-O"]
    interfaces["interfaces/<br/>ABC müqavilələri"]
    domain["domain/<br/>model obyektləri (heç kimdən asılı deyil)"]

    app --> mw
    app --> svc
    mw --> panels
    mw --> worker
    UI -->|9| rendering
    UI -->|4| reporting
    UI -->|14| APP
    UI -->|5| history
    UI -->|6| io
    UI -.->|"6 — mühərrik sinfini UI seçir"| SIM
    rendering -->|1| history
    reporting -->|3| history
    history -->|4| APP
    benchmarks --> APP
    benchmarks --> io
    APP -->|17| geology
    APP -->|13| SIM
    APP -->|1| disc
    SIM -.->|"4 — config.py idxalı"| cfg
    SIM -->|19| interfaces
    disc --> interfaces
    geology --> interfaces
    SIM -->|36| domain
    APP -->|52| domain
    geology --> domain
    io -->|20| domain
    disc --> domain
    interfaces --> domain

    linkStyle 9 stroke:#c0392b,stroke-width:2px
    linkStyle 18 stroke:#c0392b,stroke-width:2px
```

## Qeydlər

- **`simulation → application` (4 sətir):** `impes_engine.py:17`,
  `linear_solver.py:14`, `implicit/engine.py:18` və üç fazalı mühərrik
  `application/config.py`-dən `SimulationConfig`/`LinearSolverConfig`
  idxal edir. `application` isə `simulation`-dan idxal edir — iki paket
  arasında dövri asılılıq yaranır. Python bunu işlədir (modul səviyyəsində
  dövr yoxdur), lakin sənəddəki qayda formal pozulur.
- **`ui → simulation` (6 sətir):** `main_window.py` `FullyImplicitEngine` və
  `ImpesEngine` siniflərini birbaşa idxal edib `service.with_engine(...)`-ə
  ötürür (`main_window.py:2120-2123`). Mühərrik seçimi composition root-da
  deyil, UI-dədir.
- `domain/` heç bir yuxarı qatı idxal etmir — sənədləşdirilmiş qayda burada
  tam qorunur.
</content>
</invoke>
