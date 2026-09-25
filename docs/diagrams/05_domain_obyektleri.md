# Diaqram 5 — Əsas domain obyektləri (sinif diaqramı, sadələşdirilmiş)

**Mənbə:** `imex2d/application/project.py`, `imex2d/domain/reservoir_model.py:32`,
`imex2d/domain/geological_model.py:63`, `imex2d/domain/wells.py`,
`imex2d/simulation/results.py`. Yalnız əsas sahələr göstərilib.

```mermaid
classDiagram
    class Project {
        name
        geological_models: dict
        reservoir_models: dict
        runs: dict
        geology_wells: list
        ui_state: dict
        new_run(model_name, config)
    }
    class GeologicalModel {
        grid: CartesianGrid
        geometry: CellGeometry
        property_maps: dict
        facies_fields
        provenance
        validate()
    }
    class ReservoirModel {
        grid, geometry
        rock: RockProperties
        fluids: FluidProperties
        wells: list~Well~
        initial_conditions
        scal_parameters / scal_tables
        gas_scal_parameters / gas_scal_tables
        pvt_table / pvt_oil_branches
        fault_references
        validate()
        diagnose()
    }
    class SimulationRun {
        run_id
        config: SimulationConfig
        result: SimulationResult
        status
    }
    class SimulationConfig {
        end_time
        time_stepping
        linear_solver
        output
        flux_scheme: TPFA | MPFA-O
    }
    class Well {
        name
        well_type: PROD | INJ
        control: WellControl
        perforations
        tubing
    }
    class WellControl {
        mode: BHP | RATE | THP
        target
        bhp_limit
        rate_basis
        injected_phase
    }
    class SimulationResult {
        series: TimeSeries
        snapshots: list~Snapshot~
        well_oil_rate, well_bhp, well_thp ...
        ooip, ogip
        converged, message
    }
    class CartesianGrid
    class CellGeometry
    class CornerPointGeometry

    Project "1" o-- "*" GeologicalModel
    Project "1" o-- "*" ReservoirModel
    Project "1" o-- "*" SimulationRun
    GeologicalModel ..> ReservoirModel : ReservoirModelBuilder.build()
    ReservoirModel "1" o-- "*" Well
    Well "1" *-- "1" WellControl
    SimulationRun --> SimulationConfig
    SimulationRun --> SimulationResult
    ReservoirModel --> CartesianGrid
    ReservoirModel --> CellGeometry
    CellGeometry <|-- CornerPointGeometry
```
</content>
</invoke>
