# Diaqram 3 — «Simulyasiyanı işə sal» düyməsindən nəticəyə (ardıcıllıq)

**Mənbə:** `imex2d/ui/main_window.py:2115-2218`,
`imex2d/ui/worker.py:100-124`,
`imex2d/application/simulation_service.py:75-110, 205-227, 262-341`.

```mermaid
sequenceDiagram
    actor U as İstifadəçi
    participant MW as MainWindow (Qt əsas axın)
    participant SVC as ModelAwareSimulationService
    participant W as SimulationWorker (QThread)
    participant E as Mühərrik
    participant H as WellboreHydraulics

    U->>MW: «İşə sal»
    MW->>MW: rebuild_model() — panellərdən ReservoirModel
    MW->>MW: engine_choice() → FullyImplicitEngine | ImpesEngine
    MW->>SVC: with_engine(factory) — servisin SURƏTİ
    MW->>MW: reservoir_model.diagnose()
    alt xəta var
        MW-->>U: QMessageBox «Model işə salına bilməz»
    else yalnız xəbərdarlıq
        MW-->>U: «Davam edilsin?» (Yes/Cancel)
    end
    MW->>SVC: create_engine(model, config) — ön yoxlama
    Note over MW,SVC: Mühərrik burada BİR DƏFƏ tam qurulur<br/>və atılır (bax PROJECT_ANALYSIS.md, P-06)
    alt ModelValidationError / NotImplementedError
        SVC-->>MW: xəta
        MW-->>U: QMessageBox
    end
    MW->>MW: project.new_run() → RUN-00N, status RUNNING
    MW->>W: start()
    W->>SVC: run(model, config, QtProgressReporter)
    SVC->>SVC: create_engine() — İKİNCİ dəfə
    SVC->>SVC: provider-lər modeldən<br/>qaz PVT varsa → ThreePhaseSimulationEngine
    SVC->>E: engine.run(reporter)
    loop hər addım
        E->>E: Nyuton / IMPES addımı
        E-->>W: reporter.report(%, mesaj)
        W-->>MW: progress siqnalı
        MW-->>U: proqres zolağı, status sətri
        opt İstifadəçi «Dayandır»
            U->>MW: stop
            MW->>W: request_stop()
            W-->>E: report() False qaytarır → dövrə dayanır
        end
    end
    E-->>SVC: SimulationResult
    SVC->>H: annotate() — yalnız lülə həndəsəsi olan quyular
    SVC-->>W: result
    alt istisna
        W-->>MW: failed(traceback)
        MW-->>U: Jurnal tabı + QMessageBox (traceback-in son sətri)
    else uğur
        W-->>MW: finished_ok(result)
        MW->>MW: qrafiklər, xəritə, 3D, müqayisə yenilənir
        MW->>MW: _write_recovery() → %LOCALAPPDATA%\IMEX2D\berpa.imx
        MW-->>U: «Nəticələr» tabı
    end
```
</content>
</invoke>
