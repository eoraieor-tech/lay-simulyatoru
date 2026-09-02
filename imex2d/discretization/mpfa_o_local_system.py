"""MPFA-O LOKAL riyazi sistemi — Phase 5A nüvəsi.

Bax `docs/mpfa_o_phase5a.md` §4–§11 — bu modul HƏMİN düsturların
implementasiyasıdır, sətir-sətir.

Riyaziyyat (təkrar, koda yaxın formada)
---------------------------------------
Sub-cell (c,v) daxilində təzyiq XƏTTİDİR:

    p(x) = p_c + g_{c,v} · (x − x_c)

3 kəsilməzlik nöqtəsi (sub-cell-in DƏQİQ 3 sub-üzü var):

    D_{c,v} g_{c,v} = π_S − p_c 1,     D sətirləri (x_σ − x_c)ᵀ
    ⇒ g_{c,v} = D⁻¹ (π_S − p_c 1)

`c`-dən `σ` vasitəsilə ÇIXAN axın (TAM tenzor K, heç bir nᵀKn YOXDUR):

    q_{c,σ} = −Γ a_σ^(c)ᵀ K_c g_{c,v}
            = −Σ_u W[t,u] π_u + (Σ_u W[t,u]) p_c
    W = Γ · A_c · K_c · D_{c,v}⁻¹        (A_c sətirləri a_σ^(c)ᵀ)

Lokal sistem (§7) və axın bərpası (§8):

    C π_unk = D p + E π_bnd
    q_sub   = (F C⁻¹ D + G) p + (F C⁻¹ E + H) π_bnd
            =      T_cell   p +      T_bnd    π_bnd

BU MODULDA YOXDUR (tapşırıq §21/§25): mobilite, doyma, PVT, upstream,
Nyuton, Jacobian, sonlu-fərq pertürbasiyası.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Tuple

import numpy as np

from .mpfa_o_interaction import MPFAOInteractionRegion


class MPFAOBoundaryClosure(Enum):
    """Sərhəd sub-üzlərinin bağlanışı — bax `docs/mpfa_o_phase5a.md` §10.

    `DIRICHLET` (DEFOLT) — sərhəd `π_σ` NAMƏLUM DEYİL, XARİCDƏN verilən
        giriş kimi qəbul edilir. Nüvə `T_bnd` əmsallarını qaytarır ki,
        gələcək BC qatı dəyəri ötürsün. Bu, FİZİKİ sərhəd şərti İCAD
        ETMƏK DEYİL (tapşırıq §20) — sadəcə riyazi strukturdur.
    `NEUMANN_ZERO` — `q_{o,σ} = 0` (axınsız sərhəd). Bu, HƏQİQİ fiziki
        şərtdir, ona görə AÇIQ seçim tələb edir, defolt DEYİL.
    """
    DIRICHLET = "dirichlet"
    NEUMANN_ZERO = "neumann_zero"


class MPFAOTensorError(ValueError):
    """Etibarsız permeabilite tenzoru — NaN/Inf, qeyri-simmetrik, və ya
    qeyri-müsbət-müəyyən. Bax tapşırıq §18: SƏSSİZ təmir QADAĞANDIR."""


class MPFAOSingularSystemError(RuntimeError):
    """Lokal MPFA sistemi (və ya sub-cell `D` matrisi) sinqulyardır.

    Bax tapşırıq §19: səssiz TPFA-ya keçid / səssiz requlyarizasiya
    QADAĞANDIR — AÇIQ xəta verilir, diaqnostika `diagnostics` sahəsindədir.
    """

    def __init__(self, message: str, diagnostics: Optional[dict] = None):
        super().__init__(message)
        self.diagnostics = diagnostics or {}


@dataclass(frozen=True)
class MPFAORegionDiagnostics:
    """Bir lokal sistemin şərtlənmə diaqnostikası — tapşırıq §19."""
    region_id: int
    node_ijk: Tuple[int, int, int]
    n_cells: int
    n_sub_faces: int
    n_unknowns: int
    is_boundary_region: bool
    condition_number: float
    rank: int
    determinant: float
    singular: bool
    ill_conditioned: bool
    #: Sub-cell `D_{c,v}` matrislərinin ƏN PİS şərt ədədi.
    max_sub_cell_condition: float

    def as_dict(self) -> dict:
        return {
            "region_id": self.region_id, "node_ijk": self.node_ijk,
            "n_cells": self.n_cells, "n_sub_faces": self.n_sub_faces,
            "n_unknowns": self.n_unknowns,
            "is_boundary_region": self.is_boundary_region,
            "condition_number": self.condition_number, "rank": self.rank,
            "determinant": self.determinant, "singular": self.singular,
            "ill_conditioned": self.ill_conditioned,
            "max_sub_cell_condition": self.max_sub_cell_condition,
        }


def validate_permeability_matrices(k_matrices: np.ndarray,
                                   symmetry_tolerance: float = 1e-10) -> np.ndarray:
    """`(ncell,3,3)` tenzor massivini FİZİKİ etibarlılıq üçün yoxlayır.

    Bax tapşırıq §18 — `PermeabilityTensor.validate()` AVTORİTETDİR, bu
    funksiya EYNİ meyarları XAM matris massivinə tətbiq edir (MPFA nüvəsi
    `PermeabilityTensor`-suz, birbaşa matrislərlə də çağırıla bilər):

      1. NaN/Inf → rədd
      2. qeyri-simmetrik → rədd (K = Kᵀ RİYAZİ tələbdir)
      3. λ_min ≤ 0 (qeyri-SPD, sinqulyar daxil) → rədd

    HEÇ NƏ DÜZƏLDİLMİR: eigenvalue klipləmə YOX, diaqonala ε əlavəsi
    YOX, simmetrikləşdirmə YOX.
    """
    k_matrices = np.asarray(k_matrices, float)
    if k_matrices.ndim != 3 or k_matrices.shape[1:] != (3, 3):
        raise MPFAOTensorError(
            f"K massivi (ncell,3,3) olmalıdır, alındı {k_matrices.shape}")

    finite = np.isfinite(k_matrices).all(axis=(1, 2))
    if not np.all(finite):
        bad = np.flatnonzero(~finite)
        raise MPFAOTensorError(
            f"{bad.size} hüceyrədə NaN/Inf permeabilite komponenti var "
            f"(ilk hüceyrələr: {bad[:5].tolist()}) — MPFA-O bunu QƏBUL ETMİR.")

    asymmetry = np.max(np.abs(k_matrices - np.transpose(k_matrices, (0, 2, 1))),
                       axis=(1, 2))
    scale = np.maximum(np.max(np.abs(k_matrices), axis=(1, 2)), 1e-300)
    bad_sym = np.flatnonzero(asymmetry > symmetry_tolerance * scale)
    if bad_sym.size:
        raise MPFAOTensorError(
            f"{bad_sym.size} hüceyrədə K simmetrik DEYİL (ilk: {bad_sym[:5].tolist()}, "
            f"maks. asimmetriya {asymmetry.max():.3g}) — MPFA-O K=Kᵀ tələb edir "
            "və tenzoru SƏSSİZCƏ simmetrikləşdirMİR.")

    eig = np.linalg.eigvalsh(k_matrices)
    bad_spd = np.flatnonzero(eig[:, 0] <= 0.0)
    if bad_spd.size:
        raise MPFAOTensorError(
            f"{bad_spd.size} hüceyrədə K müsbət-müəyyən DEYİL (λ_min = "
            f"{eig[:, 0].min():.6g} ≤ 0, ilk: {bad_spd[:5].tolist()}) — fiziki "
            "cəhətdən etibarsız permeabilite tenzoru; MPFA-O onu TƏMİR ETMİR.")
    return k_matrices


class MPFAOLocalSystem:
    """TƏK bir qarşılıqlı təsir bölgəsinin lokal MPFA-O sistemi.

    Bütün matrislər AÇIQ atributdur (tapşırıq §8: "Do not hide the
    mathematical system behind an opaque helper"):

        C  (n_u × n_u)   sətir = naməlum π-nin bağlanış tənliyi
                         sütun = naməlum π
        D  (n_u × m)     sütun = bölgədəki hüceyrə təzyiqi p_c
        E  (n_u × n_b)   sütun = XARİCDƏN verilən sərhəd π
        F  (n_σ × n_u)   sub-üz axınının naməlum π-lərdən asılılığı
        G  (n_σ × m)     ... hüceyrə təzyiqlərindən
        H  (n_σ × n_b)   ... sərhəd π-lərindən
        T_cell (n_σ × m) = F C⁻¹ D + G      ← ÇOXNÖQTƏLİ transmissivlik
        T_bnd  (n_σ × n_b) = F C⁻¹ E + H

    `sub_cell_weights[a]` — `W = Γ A_a K_a D_a⁻¹` (3×3), `sub_cell_D[a]`
    — `D_{c,v}` (3×3). Hər ikisi audit/əl ilə yoxlama üçün saxlanılır
    (tapşırıq §30).
    """
    def __init__(self, region: MPFAOInteractionRegion, cell_centroids: np.ndarray,
                 k_matrices: np.ndarray, darcy_constant: float,
                 closure: MPFAOBoundaryClosure = MPFAOBoundaryClosure.DIRICHLET,
                 condition_warning_threshold: float = 1e12,
                 singular_tolerance: float = 1e-12):
        """`cell_centroids` — QLOBAL `(ncell,3)`; `k_matrices` — QLOBAL
        `(ncell,3,3)`. Bölgəyə aid sətirlər burada seçilir.

        Konstruktor DƏRHAL `assemble()` çağırır — obyekt ya TAM qurulur,
        ya da AÇIQ xəta ilə (`MPFAOTensorError`/`MPFAOSingularSystemError`)
        çökür; "yarım qurulmuş" vəziyyət YOXDUR.
        """
        self.region = region
        self._centroids = np.asarray(cell_centroids, float)
        k_all = np.asarray(k_matrices, float)
        if k_all.ndim != 3 or k_all.shape[1:] != (3, 3):
            raise MPFAOTensorError(f"K massivi (ncell,3,3) olmalıdır, alındı {k_all.shape}")
        self.k_matrices: np.ndarray = k_all[region.cells]
        self.darcy_constant = float(darcy_constant)
        self.closure = closure
        self.condition_warning_threshold = float(condition_warning_threshold)
        self.singular_tolerance = float(singular_tolerance)

        self.unknown_sub_faces: List[int] = []
        self.known_boundary_sub_faces: List[int] = []
        self.sub_cell_D: List[np.ndarray] = []
        self.sub_cell_D_inv: List[np.ndarray] = []
        self.sub_cell_weights: List[np.ndarray] = []
        self.C = self.D = self.E = self.F = self.G = self.H = None
        self.T_cell = self.T_bnd = None
        self._diagnostics: Optional[MPFAORegionDiagnostics] = None

        validate_permeability_matrices(self.k_matrices)
        self.assemble()

    # ── §6: sub-cell çəkiləri ────────────────────────────────────────────
    def _build_sub_cell_weights(self) -> None:
        region = self.region
        for a, cell in enumerate(region.cells):
            sub_indices = region.cell_sub_faces[a]
            centroid = self._cell_centroid(cell)
            d_matrix = np.array([region.sub_faces[t].continuity_point - centroid
                                 for t in sub_indices])
            area_rows = np.array([region.sub_faces[t].outward_area_vector(cell)
                                  for t in sub_indices])

            condition = _safe_condition(d_matrix)
            if not np.isfinite(condition) or condition > 1.0 / self.singular_tolerance:
                raise MPFAOSingularSystemError(
                    f"Bölgə node={region.node_id}, hüceyrə {cell}: sub-cell D matrisi "
                    f"sinqulyardır (şərt ədədi {condition:.3g}, det "
                    f"{np.linalg.det(d_matrix):.3g}). Kəsilməzlik nöqtələri hüceyrə "
                    "mərkəzi ilə eyni müstəvidədir — degenerativ/yastılanmış hüceyrə. "
                    "MPFA-O bunu REQULYARİZASİYA ETMİR (bax docs/mpfa_o_phase5a.md §11).",
                    diagnostics={"region_id": region.node_id, "cell": cell,
                                 "sub_cell_condition": condition,
                                 "D": d_matrix.tolist()})

            d_inv = np.linalg.inv(d_matrix)
            self.sub_cell_D.append(d_matrix)
            self.sub_cell_D_inv.append(d_inv)
            # W = Γ · A · K · D⁻¹ ; W[t,u] = π_u-nun q_{c,σ_t}-dəki çəkisi
            self.sub_cell_weights.append(
                self.darcy_constant * area_rows @ self.k_matrices[a] @ d_inv)

    def _cell_centroid(self, cell: int) -> np.ndarray:
        return self._centroids[cell]

    # ── §7/§8: lokal sistemin yığılması və həlli ────────────────────────
    def assemble(self) -> None:
        region = self.region
        self._build_sub_cell_weights()

        n_sub = region.n_sub_faces
        m = region.n_cells
        dirichlet = self.closure is MPFAOBoundaryClosure.DIRICHLET

        unknown = [s.local_index for s in region.sub_faces
                   if (not s.is_boundary) or (not dirichlet)]
        known = [s.local_index for s in region.sub_faces
                 if s.is_boundary and dirichlet]
        self.unknown_sub_faces, self.known_boundary_sub_faces = unknown, known
        col_unk = {s: t for t, s in enumerate(unknown)}
        col_bnd = {s: t for t, s in enumerate(known)}
        n_u, n_b = len(unknown), len(known)

        self.C = np.zeros((n_u, n_u))
        self.D = np.zeros((n_u, m))
        self.E = np.zeros((n_u, n_b))
        self.F = np.zeros((n_sub, n_u))
        self.G = np.zeros((n_sub, m))
        self.H = np.zeros((n_sub, n_b))

        # ── kəsilməzlik tənlikləri: hər NAMƏLUM π üçün BİR sətir ───────
        for row, s in enumerate(unknown):
            sub = region.sub_faces[s]
            for cell in sub.cells():
                a = region.cell_local[cell]
                t = region.cell_sub_faces[a].index(s)
                weights = self.sub_cell_weights[a][t]           # (3,)
                #   Σ_c Σ_u W π_u = Σ_c (Σ_u W) p_c
                for u, other in enumerate(region.cell_sub_faces[a]):
                    if other in col_unk:
                        self.C[row, col_unk[other]] += weights[u]
                    else:
                        self.E[row, col_bnd[other]] -= weights[u]
                self.D[row, a] += float(weights.sum())

        # ── axın bərpası: hər sub-üz üçün OWNER tərəfindən ─────────────
        for s, sub in enumerate(region.sub_faces):
            a = region.cell_local[sub.owner]
            t = region.cell_sub_faces[a].index(s)
            weights = self.sub_cell_weights[a][t]
            for u, other in enumerate(region.cell_sub_faces[a]):
                if other in col_unk:
                    self.F[s, col_unk[other]] -= weights[u]
                else:
                    self.H[s, col_bnd[other]] -= weights[u]
            self.G[s, a] += float(weights.sum())

        self._solve()

    def _solve(self) -> None:
        region = self.region
        n_u = self.C.shape[0]
        max_sub_cond = max((_safe_condition(d) for d in self.sub_cell_D), default=0.0)

        if n_u == 0:
            # Bütün sub-üzlər sərhəddir və Dirichlet bağlanışı seçilib —
            # naməlum yoxdur, axın BİRBAŞA giriş dəyərlərindən çıxır.
            self.T_cell, self.T_bnd = self.G.copy(), self.H.copy()
            self._diagnostics = MPFAORegionDiagnostics(
                region_id=region.node_id, node_ijk=region.node_ijk,
                n_cells=region.n_cells, n_sub_faces=region.n_sub_faces,
                n_unknowns=0, is_boundary_region=region.is_boundary_region,
                condition_number=1.0, rank=0, determinant=1.0, singular=False,
                ill_conditioned=False, max_sub_cell_condition=max_sub_cond)
            return

        condition = _safe_condition(self.C)
        rank = int(np.linalg.matrix_rank(self.C))
        determinant = float(np.linalg.det(self.C))
        singular = (rank < n_u) or (not np.isfinite(condition)) \
            or condition > 1.0 / self.singular_tolerance
        ill = (not singular) and condition > self.condition_warning_threshold

        self._diagnostics = MPFAORegionDiagnostics(
            region_id=region.node_id, node_ijk=region.node_ijk,
            n_cells=region.n_cells, n_sub_faces=region.n_sub_faces,
            n_unknowns=n_u, is_boundary_region=region.is_boundary_region,
            condition_number=float(condition), rank=rank, determinant=determinant,
            singular=singular, ill_conditioned=ill,
            max_sub_cell_condition=max_sub_cond)

        if singular:
            raise MPFAOSingularSystemError(
                f"Bölgə node={region.node_id} {region.node_ijk}: lokal MPFA-O "
                f"sistemi sinqulyardır (rank {rank}/{n_u}, şərt ədədi "
                f"{condition:.3g}, det {determinant:.3g}). MPFA-O bunu SƏSSİZCƏ "
                "requlyarizasiya ETMİR və TPFA-ya KEÇMİR "
                "(bax docs/mpfa_o_phase5a.md §11).",
                diagnostics=self._diagnostics.as_dict())

        rhs = np.hstack([self.D, self.E]) if self.E.size else self.D
        solved = np.linalg.solve(self.C, rhs)
        c_inv_d = solved[:, :self.D.shape[1]]
        c_inv_e = solved[:, self.D.shape[1]:]

        self.T_cell = self.F @ c_inv_d + self.G
        self.T_bnd = (self.F @ c_inv_e + self.H) if self.H.size else self.H.copy()

    # ── nəticələr ────────────────────────────────────────────────────────
    def diagnostics(self) -> MPFAORegionDiagnostics:
        return self._diagnostics

    def sub_face_fluxes(self, cell_pressures: np.ndarray,
                        boundary_pressures: Optional[np.ndarray] = None) -> np.ndarray:
        """`q_σ` — HƏR sub-üz üçün OWNER-dən ÇIXAN axın (§6 işarəsi).

        `cell_pressures` — bölgə sırası ilə `(m,)`; `boundary_pressures`
        — `known_boundary_sub_faces` sırası ilə `(n_b,)`.
        """
        p = np.asarray(cell_pressures, float)
        q = self.T_cell @ p
        if self.T_bnd.size:
            if boundary_pressures is None:
                raise ValueError(
                    f"Bölgə node={self.region.node_id}: Dirichlet bağlanışında "
                    f"{self.T_bnd.shape[1]} sərhəd π dəyəri TƏLƏB OLUNUR "
                    "(boundary_pressures=None verildi).")
            q = q + self.T_bnd @ np.asarray(boundary_pressures, float)
        return q

    def continuity_pressures(self, cell_pressures: np.ndarray,
                             boundary_pressures: Optional[np.ndarray] = None
                             ) -> np.ndarray:
        """Naməlum `π` vektorunun həlli — `unknown_sub_faces` sırası ilə."""
        p = np.asarray(cell_pressures, float)
        if self.C.shape[0] == 0:
            return np.zeros(0)
        rhs = self.D @ p
        if self.E.size:
            if boundary_pressures is None:
                raise ValueError("Dirichlet bağlanışı üçün boundary_pressures lazımdır.")
            rhs = rhs + self.E @ np.asarray(boundary_pressures, float)
        return np.linalg.solve(self.C, rhs)

    def half_fluxes(self, cell_pressures: np.ndarray,
                    boundary_pressures: Optional[np.ndarray] = None
                    ) -> Dict[Tuple[int, int], float]:
        """`q_{c,σ}` — HƏR (hüceyrə, sub-üz) cütü üçün AYRI-AYRILIQDA,
        hər tərəfin ÖZ `g_{c,v}` qradiyentindən hesablanmış yarım-axın.

        Bu, `sub_face_fluxes`-dən MÜSTƏQİL yoldur (o, yalnız owner
        tərəfini işlədir) — məhz buna görə `conservation_residual` lokal
        konservasiyanı HƏQİQƏTƏN yoxlaya bilir (tapşırıq §15).
        Açar: `(qlobal hüceyrə, sub-üzün bölgə-yerli indeksi)`.
        """
        region = self.region
        p = np.asarray(cell_pressures, float)
        pi_unknown = self.continuity_pressures(p, boundary_pressures)
        pi = np.zeros(region.n_sub_faces)
        for t, s in enumerate(self.unknown_sub_faces):
            pi[s] = pi_unknown[t]
        if self.known_boundary_sub_faces:
            pi_b = np.asarray(boundary_pressures, float)
            for t, s in enumerate(self.known_boundary_sub_faces):
                pi[s] = pi_b[t]

        result: Dict[Tuple[int, int], float] = {}
        for a, cell in enumerate(region.cells):
            sub_indices = region.cell_sub_faces[a]
            local_pi = np.array([pi[s] for s in sub_indices])
            for t, s in enumerate(sub_indices):
                weights = self.sub_cell_weights[a][t]
                result[(cell, s)] = float(-weights @ local_pi + weights.sum() * p[a])
        return result

    def sub_cell_gradients(self, cell_pressures: np.ndarray,
                           boundary_pressures: Optional[np.ndarray] = None
                           ) -> List[np.ndarray]:
        """`g_{c,v}` — bölgə sırası ilə hər sub-cell-in qradiyenti (§4)."""
        region = self.region
        p = np.asarray(cell_pressures, float)
        pi_unknown = self.continuity_pressures(p, boundary_pressures)
        pi = np.zeros(region.n_sub_faces)
        for t, s in enumerate(self.unknown_sub_faces):
            pi[s] = pi_unknown[t]
        if self.known_boundary_sub_faces:
            pi_b = np.asarray(boundary_pressures, float)
            for t, s in enumerate(self.known_boundary_sub_faces):
                pi[s] = pi_b[t]
        gradients = []
        for a in range(region.n_cells):
            local_pi = np.array([pi[s] for s in region.cell_sub_faces[a]])
            gradients.append(self.sub_cell_D_inv[a] @ (local_pi - p[a]))
        return gradients

    def conservation_residual(self, cell_pressures: np.ndarray,
                              boundary_pressures: Optional[np.ndarray] = None
                              ) -> float:
        """`max_σ |q_{o,σ} + q_{n,σ}|` DAXİLİ sub-üzlər üzrə (§9/§15).

        Hər tərəf MÜSTƏQİL hesablanır (`half_fluxes`), ona görə bu,
        lokal sistemin HƏQİQƏTƏN həll olunduğunun sübutudur, tavtologiya
        deyil.
        """
        halves = self.half_fluxes(cell_pressures, boundary_pressures)
        worst = 0.0
        for sub in self.region.sub_faces:
            if sub.is_boundary:
                continue
            residual = abs(halves[(sub.owner, sub.local_index)]
                           + halves[(sub.neighbor, sub.local_index)])
            worst = max(worst, residual)
        return worst

    def describe(self) -> str:
        """Əl ilə yoxlana bilən TAM dump — tapşırıq §30."""
        lines = [self.region.describe(),
                 f"  bağlanış = {self.closure.value}, Γ = {self.darcy_constant}",
                 f"  naməlum π ({len(self.unknown_sub_faces)}): "
                 f"{[self.region.sub_faces[s].face_index for s in self.unknown_sub_faces]}",
                 f"  verilən sərhəd π ({len(self.known_boundary_sub_faces)}): "
                 f"{[self.region.sub_faces[s].face_index for s in self.known_boundary_sub_faces]}"]
        for a, cell in enumerate(self.region.cells):
            lines.append(f"  hüceyrə {cell}: K =\n"
                         f"{np.array2string(self.k_matrices[a], precision=6, prefix='    ')}")
            lines.append(f"    D_(c,v) =\n"
                         f"{np.array2string(self.sub_cell_D[a], precision=6, prefix='    ')}")
            lines.append(f"    W = Γ·A·K·D⁻¹ =\n"
                         f"{np.array2string(self.sub_cell_weights[a], precision=6, prefix='    ')}")
        for name in ("C", "D", "E", "F", "G", "H", "T_cell", "T_bnd"):
            matrix = getattr(self, name)
            lines.append(f"  {name} {matrix.shape} =\n"
                         f"{np.array2string(matrix, precision=6, prefix='    ')}")
        lines.append(f"  diaqnostika: {self._diagnostics.as_dict()}")
        return "\n".join(lines)


def _safe_condition(matrix: np.ndarray) -> float:
    """`np.linalg.cond`, NaN/Inf girişində ÇÖKMƏDƏN `inf` qaytarır."""
    if not np.all(np.isfinite(matrix)):
        return float("inf")
    if matrix.size == 0:
        return 1.0
    try:
        singular_values = np.linalg.svd(matrix, compute_uv=False)
    except np.linalg.LinAlgError:
        return float("inf")
    if singular_values[-1] <= 0.0:
        return float("inf")
    return float(singular_values[0] / singular_values[-1])
