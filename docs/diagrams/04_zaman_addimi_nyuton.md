# Diaqram 4 — Tam implicit mühərrikdə bir zaman addımı və Nyuton dövrəsi

**Mənbə:** `imex2d/simulation/implicit/engine.py:192-271`,
`imex2d/simulation/implicit/time_stepping.py:118-201`,
`imex2d/simulation/implicit/newton.py:181-248`,
`imex2d/simulation/implicit/linear.py:59-90`.
Ətraflı mətn izahı: `docs/is_axini_ardicilligi.md` §4.

```mermaid
flowchart TD
    start(["Addımın başı: state, t"]) --> shares["RATE payları<br/>assign_rate_shares"]
    shares --> surf["Səth debiti proqnozu<br/>SurfaceRateController.predict"]
    surf --> adv["AdaptiveTimeStepper.advance<br/>Δt = min(dt, max_dt, qalan)"]

    adv --> newton

    subgraph newton["NewtonSolver.solve — maks. 20 iterasiya"]
        r["Qalıq R(x)<br/>PVT → kr → Pc → Φ → upstream → axın<br/>+ akkumulyasiya + quyu"]
        conv{"CNV < 1e-3<br/>və MB < 1e-7?"}
        div{"CNV > 10 × ilk CNV?"}
        j["Analitik Jakobian J"]
        red["ACTNUM daraltması"]
        lin["J·δ = −R<br/>n ≤ 20 000: splu<br/>böyük: BiCGStab + ILU/CPR"]
        chop["Appleyard kəsməsi<br/>|Δp| ≤ 50 bar, |ΔSw| ≤ 0.2"]
        ls["Xətti axtarış<br/>qalıq azalana qədər ½"]
        r --> conv
        conv -- yox --> div
        div -- yox --> j --> red --> lin --> chop --> ls --> r
    end

    conv -- bəli --> ok{"ΔSw ≤ 0.2?"}
    div -- bəli --> cut
    lin -- "xəta / NaN" --> cut
    ok -- yox --> cut
    cut{"Δt ≤ min_dt?"} -- yox --> half["Δt ← Δt × 0.5<br/>addımı təkrarla"] --> newton
    cut -- bəli --> soft{"Yumşaq qəbul?<br/>CNV < 1e-2, MB < 1e-4,<br/>ardıcıl < 20"}
    soft -- bəli --> accept
    soft -- yox --> fail(["dt = 0 → result.converged = False<br/>simulyasiya dayanır"])

    ok -- bəli --> accept["Qəbul"]
    accept --> thp["THP dövrəsi<br/>(BHP dəyişdi > tolerans → resolve_step)"]
    thp --> srl["Səth debiti dövrəsi"]
    srl --> bhp["BHP limiti dövrəsi"]
    bhp --> rec["state ← new_state, t ← t + Δt<br/>sıralar: q_o, q_w, WCT, RF, p̄"]
    rec --> snap{"Snapshot vaxtı?"}
    snap --> nextdt["Növbəti Δt:<br/>≤ 6 iter → × 1.5<br/>≤ 9 iter → × 1.0<br/>> 9 iter → × 0.5"]
    nextdt --> start
```
</content>
</invoke>
