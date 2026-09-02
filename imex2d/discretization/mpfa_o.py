"""MPFA-O diskretizasiyası — qlobal əmsal yığımı (Phase 5A).

Bax `docs/mpfa_o_phase5a.md` §8/§14/§15.

    Geometry  →  Interaction Region  →  Lokal MPFA-O  →  ƏMSALLAR (BU FAYL)

Bu fayl lokal `T_cell`/`T_bnd` matrislərini QLOBAL `(nface × ncell)` və
`(nface × n_bnd_dof)` əmsal matrislərinə yığır:

    q_F = Σ_{v ∈ künclər(F)} q_(F,v)
        = Σ_c T[F,c] p_c  +  Σ_β T_bnd[F,β] π_β

**TPFA TOXUNULMAZDIR** (tapşırıq §22): bu modul
`imex2d.simulation.discretization`-dan HEÇ NƏ İDXAL ETMİR və onun
davranışına təsir edə BİLMİR.

**Phase 5A HÜDUDU** (tapşırıq §36): `compute_flux(d_phi)` müqaviləsi
(bax `interfaces/discretization.py`) BURADA TƏTBİQ EDİLMİR — MPFA axını
`ΔΦ`-nin funksiyası DEYİL, TAM təzyiq vektorunun funksiyasıdır. Saxta
uyğunluq yaratmaq əvəzinə metod AÇIQ `NotImplementedError` verir və
Phase 5B-yə istinad edir.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy import sparse

from ..domain.general_grid_geometry import (GeneralGridGeometry,
                                            hexahedral_vertices_from_cartesian)
from ..domain.grid import CartesianGrid, Connections
from ..domain.properties import PermeabilityTensor
from ..interfaces.discretization import IFluxDiscretization
from .mpfa_global import MPFAGlobalOperator
from .mpfa_o_interaction import (MPFAOInteractionRegion, build_interaction_regions,
                                 validate_interaction_regions)
from .mpfa_o_local_system import (MPFAOBoundaryClosure, MPFAOLocalSystem,
                                  MPFAORegionDiagnostics, MPFAOSingularSystemError,
                                  validate_permeability_matrices)


@dataclass(frozen=True)
class MPFAOCoefficients:
    """Qlobal MPFA-O axın əmsalları — HƏNDƏSƏ/TOPOLOGİYA/K-dan ASILI,
    təzyiq/doyma/mobilite/PVT-dən ASILI DEYİL (tapşırıq §24).

    `T_cell` `(nface, ncell)` — `q_F = Σ_c T_cell[F,c] p_c + ...`
    `T_bnd`  `(nface, n_bnd_dof)` — sərhəd kəsilməzlik təzyiqlərinin payı.
    `boundary_dofs` — `(üz, node)` cütlərinin sırası (sərhəd DOF-ları).

    Hər ikisi `scipy.sparse.csr_matrix`-dir: sətir başına sıfırdan
    fərqli elementlərin sayı SABİT-məhduddur (§14: ≤18), ona görə
    yaddaş `O(N)`-dir. SIX matris `O(nface·ncell) = O(N²)` olardı və
    tapşırıq §31-i pozardı.

    İŞARƏ: `q_F > 0` ⟺ axın `face_owner[F]` → `face_neighbor[F]`
    (bax `docs/mpfa_o_phase5a.md` §6).
    """
    T_cell: sparse.csr_matrix
    T_bnd: sparse.csr_matrix
    boundary_dofs: List[Tuple[int, int]]
    boundary_points: np.ndarray                 #: (n_bnd_dof, 3) kəsilməzlik nöqtələri
    regions: List[MPFAOInteractionRegion]
    local_systems: List[MPFAOLocalSystem]
    diagnostics: List[MPFAORegionDiagnostics]
    eta: float
    closure: MPFAOBoundaryClosure
    darcy_constant: float
    build_seconds: float = 0.0
    #: `(üz, node) -> sərhəd DOF indeksi` (daxili axtarış cədvəli).
    _boundary_dof: Dict[Tuple[int, int], int] = field(default_factory=dict)

    @property
    def n_face(self) -> int:
        return self.T_cell.shape[0]

    @property
    def n_cell(self) -> int:
        return self.T_cell.shape[1]

    # ── axın hesablanması ────────────────────────────────────────────────
    def face_fluxes(self, pressures: np.ndarray,
                    boundary_pressures: Optional[np.ndarray] = None) -> np.ndarray:
        """`(nface,)` — hər üz üçün owner→neighbor axını.

        `boundary_pressures` `(n_bnd_dof,)` — `DIRICHLET` bağlanışında
        TƏLƏB OLUNUR (sərhəd DOF varsa). Bu, "icad edilmiş sərhəd şərti"
        DEYİL — çağıran qat verir (bax `docs/mpfa_o_phase5a.md` §10).
        """
        p = np.asarray(pressures, float)
        if p.shape != (self.n_cell,):
            raise ValueError(f"pressures ({self.n_cell},) olmalıdır, alındı {p.shape}")
        flux = np.asarray(self.T_cell @ p).ravel()
        if self.T_bnd.shape[1]:
            if boundary_pressures is None:
                raise ValueError(
                    f"{self.T_bnd.shape[1]} sərhəd DOF var və bağlanış "
                    f"{self.closure.value}-dır — `boundary_pressures` TƏLƏB OLUNUR. "
                    "`boundary_pressures_from(fn)` köməkçisinə bax.")
            pi_b = np.asarray(boundary_pressures, float)
            if pi_b.shape != (self.T_bnd.shape[1],):
                raise ValueError(f"boundary_pressures ({self.T_bnd.shape[1]},) "
                                 f"olmalıdır, alındı {pi_b.shape}")
            flux = flux + np.asarray(self.T_bnd @ pi_b).ravel()
        return flux

    def boundary_pressures_from(self, pressure_function: Callable[[np.ndarray], float]
                                ) -> np.ndarray:
        """Analitik təzyiq funksiyasından sərhəd `π` vektoru — manufactured
        solution testləri üçün (tapşırıq §14). Fiziki BC İCAD ETMİR:
        çağıran ARTIQ bildiyi həlli ötürür."""
        if not self.boundary_dofs:
            return np.zeros(0)
        return np.array([float(pressure_function(x)) for x in self.boundary_points])

    # ── stensil müayinəsi (tapşırıq §10/§29/§F) ──────────────────────────
    def face_stencil(self, face: int, tolerance: float = 0.0) -> Dict[int, float]:
        """`{hüceyrə: əmsal}` — `face` üzündəki axının HANSI hüceyrə
        təzyiqlərindən asılı olduğu. `len(...) > 2` ⟺ HƏQİQİ çoxnöqtəli
        bağlantı."""
        row = self.T_cell.getrow(face).toarray().ravel()
        return {int(c): float(row[c])
                for c in np.flatnonzero(np.abs(row) > tolerance)}

    def stencil_sizes(self, tolerance: float = 1e-12) -> np.ndarray:
        """`(nface,)` — hər üzün stensil ölçüsü (NİSBİ tolerantlıqla —
        sətrin ən böyük əmsalına görə)."""
        matrix = self.T_cell.tocsr()
        sizes = np.zeros(matrix.shape[0], dtype=int)
        for face in range(matrix.shape[0]):
            values = np.abs(matrix.data[matrix.indptr[face]:matrix.indptr[face + 1]])
            if values.size:
                sizes[face] = int(np.sum(values > tolerance * max(values.max(), 1e-300)))
        return sizes

    # ── diaqnostika (tapşırıq §19) ───────────────────────────────────────
    def condition_numbers(self) -> np.ndarray:
        return np.array([d.condition_number for d in self.diagnostics])

    def ill_conditioned_regions(self) -> List[MPFAORegionDiagnostics]:
        return [d for d in self.diagnostics if d.ill_conditioned]

    def region_by_node(self, node_id: int) -> MPFAOLocalSystem:
        for system in self.local_systems:
            if system.region.node_id == node_id:
                return system
        raise KeyError(f"node_id={node_id} üçün bölgə yoxdur.")

    def conservation_residual(self, pressures: np.ndarray,
                              boundary_pressures: Optional[np.ndarray] = None) -> float:
        """BÜTÜN bölgələr üzrə `max |q_{o,σ} + q_{n,σ}|` (tapşırıq §15).

        Hər bölgədə hər tərəf MÜSTƏQİL qradiyentdən hesablanır — bax
        `MPFAOLocalSystem.conservation_residual`.
        """
        p = np.asarray(pressures, float)
        pi_b = None if boundary_pressures is None else np.asarray(boundary_pressures, float)
        worst = 0.0
        for system in self.local_systems:
            local_p = p[system.region.cells]
            local_b = None
            if system.known_boundary_sub_faces:
                if pi_b is None:
                    raise ValueError("Dirichlet bağlanışı üçün boundary_pressures lazımdır.")
                local_b = pi_b[[self._boundary_dof[(system.region.sub_faces[s].face_index,
                                                    system.region.node_id)]
                                for s in system.known_boundary_sub_faces]]
            worst = max(worst, system.conservation_residual(local_p, local_b))
        return worst


def build_mpfa_o_coefficients(grid: CartesianGrid, geometry: GeneralGridGeometry,
                              k_matrices: np.ndarray, darcy_constant: float,
                              eta: float = 1.0,
                              closure: MPFAOBoundaryClosure = MPFAOBoundaryClosure.DIRICHLET,
                              condition_warning_threshold: float = 1e12,
                              ) -> MPFAOCoefficients:
    """MPFA-O nüvəsinin SAF girişi — `ReservoirModel`-siz (tapşırıq §21).

    Giriş YALNIZ: topologiya (`grid`), həndəsə (`geometry`), tam tenzor
    permeabilite (`k_matrices`, `(ncell,3,3)`) və vahid sabiti. Nyuton/
    PVT/quyu/doyma YOXDUR.

    Mürəkkəblik `O(N)`: bölgə sayı `(nx+1)(ny+1)(nz+1)`, hər bölgədə
    SABİT ölçülü (≤12×12) lokal həll (tapşırıq §31).
    """
    started = time.perf_counter()
    k_matrices = validate_permeability_matrices(k_matrices)
    if k_matrices.shape[0] != geometry.ncell:
        raise ValueError(f"K massivi {k_matrices.shape[0]} hüceyrə üçündür, "
                         f"həndəsədə {geometry.ncell} hüceyrə var.")

    regions = build_interaction_regions(grid, geometry, eta=eta)
    issues = validate_interaction_regions(regions)
    if issues:
        raise ValueError("MPFA-O bölgələri etibarsızdır:\n  - " + "\n  - ".join(issues[:10]))

    centroids = geometry.cell_centroids
    n_face, n_cell = len(geometry.faces), geometry.ncell

    # ── sərhəd DOF-larının qlobal nömrələnməsi (§12) ─────────────────────
    boundary_dofs: List[Tuple[int, int]] = []
    boundary_points: List[np.ndarray] = []
    boundary_dof: Dict[Tuple[int, int], int] = {}
    dirichlet = closure is MPFAOBoundaryClosure.DIRICHLET
    if dirichlet:
        for region in regions:
            for s in region.boundary_sub_faces:
                sub = region.sub_faces[s]
                key = (sub.face_index, region.node_id)
                boundary_dof[key] = len(boundary_dofs)
                boundary_dofs.append(key)
                boundary_points.append(sub.continuity_point)

    # COO üçlükləri — SIX matris `O(nface·ncell)` yaddaş tutardı (§31).
    cell_rows: List[np.ndarray] = []
    cell_cols: List[np.ndarray] = []
    cell_values: List[np.ndarray] = []
    bnd_rows: List[np.ndarray] = []
    bnd_cols: List[np.ndarray] = []
    bnd_values: List[np.ndarray] = []
    systems: List[MPFAOLocalSystem] = []
    diagnostics: List[MPFAORegionDiagnostics] = []

    for region in regions:
        system = MPFAOLocalSystem(
            region=region, cell_centroids=centroids, k_matrices=k_matrices,
            darcy_constant=darcy_constant, closure=closure,
            condition_warning_threshold=condition_warning_threshold)
        systems.append(system)
        diagnostics.append(system.diagnostics())

        faces = np.array([sub.face_index for sub in region.sub_faces], dtype=int)
        cells = np.array(region.cells, dtype=int)
        cell_rows.append(np.repeat(faces, cells.size))
        cell_cols.append(np.tile(cells, faces.size))
        cell_values.append(system.T_cell.ravel())

        if system.T_bnd.size:
            dofs = np.array([boundary_dof[(region.sub_faces[s].face_index,
                                           region.node_id)]
                             for s in system.known_boundary_sub_faces], dtype=int)
            bnd_rows.append(np.repeat(faces, dofs.size))
            bnd_cols.append(np.tile(dofs, faces.size))
            bnd_values.append(system.T_bnd.ravel())

    T_cell = sparse.coo_matrix(
        (np.concatenate(cell_values),
         (np.concatenate(cell_rows), np.concatenate(cell_cols))),
        shape=(n_face, n_cell)).tocsr()
    T_bnd = sparse.coo_matrix(
        (np.concatenate(bnd_values) if bnd_values else np.zeros(0),
         (np.concatenate(bnd_rows) if bnd_rows else np.zeros(0, int),
          np.concatenate(bnd_cols) if bnd_cols else np.zeros(0, int))),
        shape=(n_face, len(boundary_dofs))).tocsr()
    # `tocsr()` təkrarlanan (sətir, sütun) girişlərini AVTOMATİK CƏMLƏYİR
    # — bir üzün 4 küncündən gələn paylar məhz belə toplanır (§8).
    # Lokal `T_cell` blokları SIX olduğu üçün struktur SIFIRLAR da gəlir;
    # onlar atılır (ƏDƏDİ nəticə DƏYİŞMİR, yalnız yaddaş azalır — bax
    # `docs/mpfa_o_phase5b1.md` §9).
    T_cell.eliminate_zeros()
    T_bnd.eliminate_zeros()

    coefficients = MPFAOCoefficients(
        T_cell=T_cell, T_bnd=T_bnd, boundary_dofs=boundary_dofs,
        boundary_points=(np.array(boundary_points) if boundary_points
                         else np.zeros((0, 3))),
        regions=regions, local_systems=systems, diagnostics=diagnostics,
        eta=eta, closure=closure, darcy_constant=float(darcy_constant),
        build_seconds=time.perf_counter() - started,
        _boundary_dof=boundary_dof)
    return coefficients


@dataclass(frozen=True)
class MPFADiscretizedGrid:
    """`MPFAODiscretization.build()`-in qaytardığı obyekt.

    `DiscretizedGrid` (TPFA) ilə EYNİ ORTAQ sahələri daşıyır
    (`connections`, `pore_volume`, `cell_volume`, `warnings`) — ona görə
    həndəsə/həcm istifadə edən mövcud kod işləyir — AMMA `compute_flux
    (d_phi)` müqaviləsini TƏTBİQ ETMİR (bax modul docstring-i, Phase 5A
    hüdudu).
    """
    connections: Connections
    pore_volume: np.ndarray
    cell_volume: np.ndarray
    geometry: GeneralGridGeometry
    coefficients: MPFAOCoefficients
    #: `Connections` sırası ilə uyğun QLOBAL üz indeksləri — TPFA ilə
    #: MÜQAYİSƏ üçün YEGANƏ düzgün körpü (bax `GeneralGridGeometry.
    #: connection_faces`).
    connection_faces: np.ndarray = field(default_factory=lambda: np.zeros(0, int))
    warnings: List[str] = field(default_factory=list)
    #: MPFA-nın HƏLƏ DƏSTƏKLƏMƏDİYİ, amma modeldə MÖVCUD olan
    #: xüsusiyyətlər (məs. fay). Qalıq qatı bu siyahı BOŞ DEYİLSƏ
    #: hesablamanı RƏDD EDİR — səssiz, yanıldıcı nəticə QADAĞANDIR
    #: (bax `docs/mpfa_o_phase5b1.md` §8, tapşırıq §26).
    unsupported_features: List[str] = field(default_factory=list)
    #: Phase 5B-1 QLOBAL operatoru (`MPFAGlobalOperator`) — `None` yalnız
    #: `Connections` verilmədikdə.
    global_operator: Optional[MPFAGlobalOperator] = None

    def supports_multipoint_stencil(self) -> bool:
        """`ResidualAssembler` bu metodla MPFA yolunu seçir (duck-typed;
        `DiscretizedGrid`-də bu metod YOXDUR → TPFA yolu, bax
        `docs/mpfa_o_phase5b1.md` §5)."""
        return True

    # ── Phase 5B-1: potensial → axın (qalıq qatının işlətdiyi giriş) ────
    def connection_fluxes_from_potential(self, potential: np.ndarray,
                                         boundary_potential=None) -> np.ndarray:
        """`(nconn,)` — HÜCEYRƏ üzrə skalyar potensialdan (`Φ_α`) baza
        Darsi axını, `cell_a → cell_b`, mobilitəsiz.

        Bu, `IFluxDiscretization` müqaviləsinin ÇOXNÖQTƏLİ qarşılığıdır
        (tapşırıq §4): giriş TAM vektordur, tək `Δp` deyil."""
        return self._require_operator().connection_fluxes(potential,
                                                          boundary_potential)

    def upstream_cells(self, flux: np.ndarray) -> np.ndarray:
        """HƏQİQİ çoxnöqtəli axının işarəsindən upstream hüceyrəsi (§10)."""
        return self._require_operator().upstream_cells(flux)

    def global_stencil_pattern(self):
        """HƏQİQİ MPFA bağlantı naxışı (§22)."""
        return self._require_operator().global_stencil_pattern()

    def _require_operator(self) -> MPFAGlobalOperator:
        if self.global_operator is None:
            raise RuntimeError(
                "Qlobal MPFA operatoru qurulmayıb — `MPFAODiscretization.build()` "
                "`Connections` olan model tələb edir.")
        return self.global_operator

    def compute_flux(self, d_phi: np.ndarray) -> np.ndarray:
        """**Phase 5A-da QƏSDƏN TƏTBİQ EDİLMƏYİB.**

        MPFA-O axını `q_F = Σ_c T[F,c] p_c` — TAM təzyiq vektorunun
        funksiyasıdır, üz-başına `ΔΦ`-nin DEYİL. `d_phi`-dən "təxmini"
        bir cavab uydurmaq metodun ÇOXNÖQTƏLİ mahiyyətini gizlədərdi
        (yəni TPFA-ya çevirərdi) — bu, tapşırıq §1/§22 ilə QADAĞANDIR.

        Residual/Jacobian inteqrasiyası Phase 5B-nin işidir (bax
        `docs/mpfa_o_phase5a.md` §17).
        """
        raise NotImplementedError(
            "MPFA-O `compute_flux(d_phi)` müqaviləsini TƏTBİQ ETMİR: çoxnöqtəli "
            "axın tək üzün ΔΦ-sindən çıxarıla BİLMƏZ. `compute_flux_from_pressure"
            "(p, boundary_pressures)` işlədin. Residual/Nyuton inteqrasiyası — "
            "Phase 5B (bax docs/mpfa_o_phase5a.md §17).")

    def compute_flux_from_pressure(self, pressures: np.ndarray,
                                   boundary_pressures: Optional[np.ndarray] = None
                                   ) -> np.ndarray:
        """`(nface,)` — MPFA-O baza Darsi axını (mobilite/upstream OLMADAN)."""
        return self.coefficients.face_fluxes(pressures, boundary_pressures)

    def connection_fluxes(self, pressures: np.ndarray,
                          boundary_pressures: Optional[np.ndarray] = None) -> np.ndarray:
        """`(nconn,)` — `Connections` sırasında, `cell_a`→`cell_b` axını.

        `Connections`-da `cell_a` HƏMİŞƏ aşağı indeksdir və
        `GeneralGridGeometry` həmin üzü məhz `cell_a` owner-i ilə qurur,
        ona görə işarə çevirməsi LAZIM DEYİL."""
        flux = self.compute_flux_from_pressure(pressures, boundary_pressures)
        return flux[self.connection_faces]


class MPFAODiscretization(IFluxDiscretization):
    """TAM tenzorlu, çoxnöqtəli MPFA-O diskretizasiyası (Phase 5A).

    `TwoPointFluxDiscretization`-ın YANINDA, ONU ƏVƏZ ETMƏDƏN yaşayır —
    `default_flux_discretization()` HƏLƏ DƏ TPFA qaytarır (tapşırıq §22).
    Bu sinifdən istifadə AÇIQ seçim tələb edir.

    Parametrlər
    -----------
    `eta` — kəsilməzlik nöqtəsi parametri η ∈ (0,1], defolt `1.0`
        (bax `docs/mpfa_o_phase5a.md` §4/§16).
    `closure` — sərhəd bağlanışı, defolt `DIRICHLET` (§10).
    `condition_warning_threshold` — bundan böyük şərt ədədi
        diaqnostikada `ill_conditioned` kimi işarələnir (§11);
        sinqulyarlıq HƏR HALDA AÇIQ xəta verir.
    """

    def __init__(self, eta: float = 1.0,
                 closure: MPFAOBoundaryClosure = MPFAOBoundaryClosure.DIRICHLET,
                 condition_warning_threshold: float = 1e12):
        self.eta = float(eta)
        self.closure = closure
        self.condition_warning_threshold = float(condition_warning_threshold)

    # ── SAF nüvə girişi (ReservoirModel-siz) ─────────────────────────────
    def build_coefficients(self, grid: CartesianGrid, geometry: GeneralGridGeometry,
                           k_matrices: np.ndarray,
                           darcy_constant: float) -> MPFAOCoefficients:
        return build_mpfa_o_coefficients(
            grid, geometry, k_matrices, darcy_constant, eta=self.eta,
            closure=self.closure,
            condition_warning_threshold=self.condition_warning_threshold)

    # ── IFluxDiscretization ──────────────────────────────────────────────
    def supports_multipoint_stencil(self) -> bool:
        return True

    def build(self, model) -> MPFADiscretizedGrid:
        """`ReservoirModel` → MPFA-O diskretizasiyası.

        Həndəsə: `CellGeometry`-dən hekzahedral təpələr qurulur
        (`hexahedral_vertices_from_cartesian`) — Kartezian modellər üçün
        DƏQİQ. Qeyri-ortoqonal/corner-point təpə mənbəyi HƏLƏ YOXDUR
        (Phase 5D); belə həndəsə ilə işləmək üçün `build_coefficients`-i
        BİRBAŞA öz `GeneralGridGeometry`-nizlə çağırın.
        """
        conn = model.connections()
        vertices = hexahedral_vertices_from_cartesian(model.grid, model.geometry)
        geometry = GeneralGridGeometry(vertices, conn)
        k_matrices = permeability_matrices(model)

        warnings: List[str] = []
        unsupported: List[str] = []
        if any(f.has_geometry for f in model.fault_references):
            unsupported.append(
                "fay (fault) transmissivlik çarpanları — MPFA-da HƏLƏ "
                "implement edilməyib (Phase 5D)")
            warnings.append(
                "Modeldə fay (fault) var, AMMA MPFA-O (Phase 5A) fay "
                "transmissivlik çarpanlarını TƏTBİQ ETMİR — bu, TPFA-da "
                "işləyir (`TwoPointFluxDiscretization._apply_fault_multipliers`), "
                "MPFA-da isə Phase 5D işidir. Alınan MPFA axınları fayları "
                "NƏZƏRƏ ALMIR.")
        if self.closure is not MPFAOBoundaryClosure.NEUMANN_ZERO:
            unsupported.append(
                f"'{self.closure.value}' sərhəd bağlanışı — qalıq (residual) qatı "
                "hələ sərhəd π dəyərlərini ötürmür (Phase 5B-2). Qalıq yolu üçün "
                "`MPFAOBoundaryClosure.NEUMANN_ZERO` (axınsız sərhəd) seçin — "
                "bu, mövcud simulyatorun ARTIQ tətbiq etdiyi şərtdir")
        if model.rock.permeability_tensor is None:
            warnings.append(
                "`rock.permeability_tensor` verilməyib — MPFA-O DİAQONAL "
                "K = diag(PERMX, PERMY, PERMZ) işlədir. Bu, riyazi cəhətdən "
                "etibarlıdır, amma off-diaqonal anizotropluq MODELDƏ YOXDUR.")

        coefficients = self.build_coefficients(
            model.grid, geometry, k_matrices, model.units.darcy_constant)

        connection_faces = geometry.connection_faces()
        operator = MPFAGlobalOperator(
            coefficients=coefficients, connections=conn,
            connection_faces=connection_faces, face_owner=geometry.face_owner,
            face_neighbor=geometry.face_neighbor)

        return MPFADiscretizedGrid(
            connections=conn, pore_volume=model.pore_volume(),
            cell_volume=model.geometry.volumes(), geometry=geometry,
            coefficients=coefficients, connection_faces=connection_faces,
            warnings=warnings, unsupported_features=unsupported,
            global_operator=operator)


def permeability_matrices(model) -> np.ndarray:
    """`(ncell,3,3)` tam tenzor massivi — `ReservoirModel`-dən.

    `rock.permeability_tensor` VARSA o işlədilir (`as_matrices()` — tam
    6 komponent). YOXDURSA diaqonal `diag(PERMX, PERMY, PERMZ)` qurulur
    (PERMZ yoxdursa PERMX). Heç bir dəyər "təmir edilmir"/klipənmir —
    etibarsız tenzor `validate_permeability_matrices` ilə RƏDD EDİLİR.
    """
    tensor: Optional[PermeabilityTensor] = model.rock.permeability_tensor
    if tensor is not None:
        return tensor.as_matrices()
    kx = model.rock.permx.values
    ky = model.rock.permy.values
    kz = model.rock.permz.values if model.rock.permz is not None else kx
    matrices = np.zeros((kx.size, 3, 3))
    matrices[:, 0, 0] = kx
    matrices[:, 1, 1] = ky
    matrices[:, 2, 2] = kz
    return matrices
