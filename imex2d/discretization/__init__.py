"""MPFA-O diskretizasiya nüvəsi (Phase 5A).

Bax `docs/mpfa_o_phase5a.md` — TAM riyazi spesifikasiya.

Qat diaqramı (bax `ARCHITECTURE.md` §5.16):

    Geometry (`domain/general_grid_geometry.py::GeneralGridGeometry`)
        ↓
    Interaction Region (`mpfa_o_interaction.py::MPFAOInteractionRegion`)
        ↓
    Lokal MPFA-O riyaziyyatı (`mpfa_o_local_system.py::MPFAOLocalSystem`)
        ↓
    Axın əmsalları (`mpfa_o.py::MPFAOCoefficients` / `MPFAODiscretization`)

Bu paket NYUTON/PVT/QUYU/DOYMA/ZAMAN ADDIMI anlayışlarını TANIMIR
(tapşırıq §21) — girişi YALNIZ `{həndəsə, topologiya, K, Γ}`-dır.
TPFA (`imex2d.simulation.discretization`) TOXUNULMAZDIR.
"""

from .mpfa_o import (MPFAOCoefficients, MPFAODiscretization,
                     MPFADiscretizedGrid, build_mpfa_o_coefficients)
from .mpfa_o_interaction import (MPFAOInteractionRegion, MPFAOSubFace,
                                 build_interaction_regions)
from .mpfa_o_local_system import (MPFAOBoundaryClosure, MPFAOLocalSystem,
                                  MPFAORegionDiagnostics, MPFAOSingularSystemError,
                                  MPFAOTensorError, validate_permeability_matrices)

__all__ = [
    "MPFAOCoefficients", "MPFAODiscretization", "MPFADiscretizedGrid",
    "build_mpfa_o_coefficients",
    "MPFAOInteractionRegion", "MPFAOSubFace", "build_interaction_regions",
    "MPFAOBoundaryClosure", "MPFAOLocalSystem", "MPFAORegionDiagnostics",
    "MPFAOSingularSystemError", "MPFAOTensorError",
    "validate_permeability_matrices",
]
