# Diaqram 2 — Ümumi məlumat axını: girişdən nəticəyə

**Mənbə:** `imex2d/application/geology_service.py:441`,
`imex2d/application/model_builder.py:30`,
`imex2d/application/simulation_service.py:262-271`,
`imex2d/application/serialization.py:367-393`, `docs/is_axini_ardicilligi.md` §0.

```mermaid
flowchart LR
    subgraph IN["GİRİŞ"]
        csv["Quyu CSV<br/>well,x,y,k,PORO,PERMX,NTG<br/>(geology/well_data_io.py)"]
        table["UI quyu cədvəli<br/>(GeologyPanel)"]
        synth["Sintetik geologiya<br/>(application/scenarios.py)"]
        grdecl["GRDECL<br/>SPECGRID, COORD/ZCORN, PORO…<br/>(io/grdecl.py)"]
        deck["Eclipse deck<br/>PVTO/PVDG/PVTW, SWOF/SGOF<br/>(io/pvt_io.py, io/scal_io.py)"]
        faults["Fay CSV / FAULTS<br/>(io/fault_io.py)"]
        panels["UI panelləri<br/>PVT, SCAL, quyular, ədədi"]
        obs["Müşahidə CSV<br/>(history/observation_io.py)"]
        imx[".imx layihə faylı<br/>gzip + JSON"]
    end

    subgraph GEO["GEOLOGİYA (bir dəfə)"]
        qc["QC: dublikat, kənar dəyər"]
        tr["Çevirmə: ln, logit, normal-score"]
        vg["Variogram fit"]
        kr["Kriging / IDW / SGS / SIS"]
        cv["Çarpaz-doğrulama"]
        qc --> tr --> vg --> kr --> cv
    end

    geomodel[("GeologicalModel")]
    resmodel[("ReservoirModel<br/>validate()")]

    subgraph ENG["MÜHƏRRİK"]
        pre["Bir dəfə: T (TPFA/MPFA-O),<br/>PV, WI (Peaceman), ilkin p/Sw/Sg"]
        loop["Zaman dövrəsi:<br/>IMPES və ya Nyuton (2/3 faza)"]
        pre --> loop
    end

    result[("SimulationResult<br/>sıralar + 3D anlar")]

    subgraph OUT["ÇIXIŞ"]
        charts["Qrafiklər (matplotlib)"]
        vtk["3D (VTK / psevdo-3D)"]
        daily["Günlük cədvəl + CSV"]
        exp["CSV / JSON ixracı"]
        pdf["PDF hesabat"]
        ecl["Eclipse .DATA ixracı"]
        gif["PNG / GIF animasiya"]
        save[".imx saxlama + bərpa faylı"]
    end

    hm["History matching / həssaslıq<br/>(bütün dövrəni N dəfə təkrarlayır)"]

    csv --> GEO
    table --> GEO
    synth --> geomodel
    grdecl --> geomodel
    GEO --> geomodel
    geomodel -->|ReservoirModelBuilder| resmodel
    deck --> resmodel
    faults --> resmodel
    panels --> resmodel
    imx --> resmodel
    resmodel -->|SimulationService| ENG
    ENG --> result
    result --> charts & vtk & daily & exp & pdf & gif & save
    resmodel --> ecl
    obs --> hm
    result --> hm
    hm -->|parametr dəyişir| resmodel
```
</content>
</invoke>
