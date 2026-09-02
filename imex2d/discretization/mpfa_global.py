"""MPFA-O QLOBAL operatoru — Phase 5B-1.

Bax `docs/mpfa_o_phase5b1.md` §3/§9 — bu modul HƏMİN spesifikasiyanın
implementasiyasıdır.

Nə edir
-------
Phase 5A `MPFAOCoefficients.T_cell` `(nface × ncell)` seyrək matrisini
QALIQ qatının işlədə biləcəyi formaya gətirir:

    T_face  (nface × ncell)   ← Phase 5A, OLDUĞU KİMİ saxlanılır
        ↓  sətir seçimi (connection_faces — DETERMİNİSTİK biyeksiya)
    T_conn  (nconn × ncell)   ← `Connections` sırasında, cell_a → cell_b
        ↓  divergensiya D (nconn → ncell scatter)
    A = D · T_conn (ncell × ncell)   ← SABİT mobilitəli XƏTTİ operator
                                        + HƏQİQİ bağlantı naxışı

Nə ETMİR
--------
* Phase 5A əmsallarını YENİDƏN HESABLAMIR/yığmır (§3 qadağası).
* Mobilitə/doyma/PVT/Nyuton/quyu TANIMIR (§11 ayrılığı) — girişi
  YALNIZ skalyar potensial vektorudur.
* `A`-nı qeyri-xətti Jacobian kimi TƏQDİM ETMİR (§23):
  `xətti MPFA operatoru ≠ tam rezervuar Jacobian-ı`.
* HEÇ BİR yerdə `np.zeros((ncell, ncell))` (§5).

İşarə konvensiyası (mövcud simulyatorla EYNİ)
---------------------------------------------
`connection_fluxes(...)[k] > 0`  ⟺  axın `cell_a[k] → cell_b[k]`.
`Connections.cell_a` HƏMİŞƏ aşağı indeksdir və `GeneralGridGeometry`
həmin üzü məhz `cell_a` owner-i ilə qurur, ona görə üz `owner→neighbor`
konvensiyası ilə `cell_a→cell_b` konvensiyası ÜST-ÜSTƏ DÜŞÜR —
işarə çevirməsi LAZIM DEYİL (`_validate_orientation` bunu yoxlayır).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, List, Optional, Tuple

import numpy as np
from scipy import sparse

from ..domain.grid import Connections

if TYPE_CHECKING:                      # dairəvi idxalın qarşısını alır —
    from .mpfa_o import MPFAOCoefficients   # `mpfa_o` BU modulu idxal edir


class MPFAStateError(ValueError):
    """Qalıq qatından gələn dövlət vektoru etibarsızdır (NaN/Inf, yanlış
    ölçü) — bax tapşırıq §27: SƏSSİZ təmir QADAĞANDIR."""


class MPFAGlobalOperator:
    """Phase 5A əmsallarından qurulmuş QLOBAL seyrək MPFA operatoru.

    Bu obyekt DÖVLƏTDƏN (təzyiq/doyma/mobilite) TAM ASILI DEYİL — bir
    dəfə qurulur və bütün Nyuton/qalıq qiymətləndirmələrində TƏKRAR
    İSTİFADƏ olunur (tapşırıq §11/§28).
    """

    def __init__(self, coefficients: "MPFAOCoefficients", connections: Connections,
                 connection_faces: np.ndarray, face_owner: np.ndarray,
                 face_neighbor: np.ndarray):
        self.coefficients = coefficients
        self.connections = connections
        self.connection_faces = np.asarray(connection_faces, dtype=int)
        self.ncell = int(coefficients.n_cell)
        self.nconn = int(connections.count)

        if self.connection_faces.size != self.nconn:
            raise MPFAStateError(
                f"connection_faces ölçüsü ({self.connection_faces.size}) "
                f"əlaqə sayı ({self.nconn}) ilə uyğun gəlmir.")
        self._validate_orientation(face_owner, face_neighbor)

        #: (nconn × ncell) — Phase 5A sətirlərinin SEÇİMİ, YENİDƏN
        #: HESABLAMA YOX (§3).
        self.T_conn: sparse.csr_matrix = coefficients.T_cell[self.connection_faces]
        self._divergence = self._build_divergence()
        self._cell_operator: Optional[sparse.csr_matrix] = None
        self._boundary_index = coefficients._boundary_dof

    # ────────────────────────────────────────────────── qurulma/yoxlama
    def _validate_orientation(self, face_owner: np.ndarray,
                              face_neighbor: np.ndarray) -> None:
        """`Connections` ↔ fiziki üz biyeksiyasının İSTİQAMƏTİ.

        Bu, sükutla yanlış işarə verə biləcək YEGANƏ xəritələmə
        səhvidir, ona görə AÇIQ yoxlanılır (tapşırıq §17/§27)."""
        if self.nconn == 0:
            return
        owner = np.asarray(face_owner, int)[self.connection_faces]
        neighbor = np.asarray(face_neighbor, int)[self.connection_faces]
        bad = np.flatnonzero((owner != self.connections.cell_a)
                             | (neighbor != self.connections.cell_b))
        if bad.size:
            k = int(bad[0])
            raise MPFAStateError(
                f"Üz↔əlaqə xəritələməsi pozulub: əlaqə {k} "
                f"(cell_a={self.connections.cell_a[k]}, "
                f"cell_b={self.connections.cell_b[k]}) üzü {self.connection_faces[k]}-dir, "
                f"amma həmin üzün owner/neighbor-u {owner[k]}/{neighbor[k]}-dir. "
                "MPFA axınlarının işarəsi etibarsız olardı.")

    def _build_divergence(self) -> sparse.csr_matrix:
        """`D` (ncell × nconn): `cell_a`-ya `−1`, `cell_b`-yə `+1`.

        Bu, `ResidualAssembler.net_influx`-un `np.add.at` məntiqinin
        MATRİS formasıdır — EYNİ işarə, EYNİ sahiblik (§8)."""
        conn = self.connections
        rows = np.concatenate([conn.cell_a, conn.cell_b])
        cols = np.concatenate([np.arange(self.nconn), np.arange(self.nconn)])
        values = np.concatenate([-np.ones(self.nconn), np.ones(self.nconn)])
        return sparse.coo_matrix((values, (rows, cols)),
                                 shape=(self.ncell, self.nconn)).tocsr()

    # ──────────────────────────────────────────────── axın hesablanması
    @property
    def n_boundary_dof(self) -> int:
        return int(self.coefficients.T_bnd.shape[1])

    def _boundary_contribution(self, matrix, boundary_potential) -> np.ndarray:
        """Sərhəd `π` payı. Sərhəd DOF varsa VƏ dəyər verilməyibsə
        AÇIQ xəta — SƏSSİZCƏ SIFIR SAYMAQ QADAĞANDIR (tapşırıq §13/§27):
        `T_bnd` payını buraxmaq axını fiziki cəhətdən YANLIŞ edərdi."""
        if not self.n_boundary_dof:
            return 0.0
        if boundary_potential is None:
            raise MPFAStateError(
                f"Bu operator {self.n_boundary_dof} sərhəd DOF-lu "
                f"`{self.coefficients.closure.value}` bağlanışı ilə qurulub — "
                "`boundary_potential` TƏLƏB OLUNUR. Onu buraxmaq `T_bnd` payını "
                "SƏSSİZCƏ atmaq və yanlış axın vermək olardı. Qalıq qatı üçün "
                "`MPFAOBoundaryClosure.NEUMANN_ZERO` işlədin "
                "(bax docs/mpfa_o_phase5b1.md §7).")
        values = np.asarray(boundary_potential, float)
        if values.shape != (self.n_boundary_dof,):
            raise MPFAStateError(f"boundary_potential ({self.n_boundary_dof},) "
                                 f"olmalıdır, alındı {values.shape}")
        return np.asarray(matrix @ values).ravel()

    def connection_fluxes(self, potential: np.ndarray,
                          boundary_potential: Optional[np.ndarray] = None
                          ) -> np.ndarray:
        """`(nconn,)` — baza (mobilitəsiz) Darsi axını, `cell_a → cell_b`.

        `potential` — HÜCEYRƏ üzrə skalyar sahə: təzyiq VƏ YA faza
        potensialı `Φ_α` (bax `docs/mpfa_o_phase5b1.md` §4). MPFA
        operatoru xətti olduğu üçün ikisi də eyni şəkildə tətbiq olunur.
        """
        phi = self._validate_field(potential, "potensial")
        flux = np.asarray(self.T_conn @ phi).ravel()
        contribution = self._boundary_contribution(
            self.coefficients.T_bnd[self.connection_faces], boundary_potential)
        return flux + contribution

    def face_fluxes(self, potential: np.ndarray,
                    boundary_potential: Optional[np.ndarray] = None) -> np.ndarray:
        """`(nface,)` — BÜTÜN fiziki üzlər (sərhəd daxil). Qalıq qatı
        bunu işlətmir (sərhəd üzlərinin axını `NEUMANN_ZERO`-da ≡ 0),
        amma konservasiya/diaqnostika testləri üçün lazımdır."""
        phi = self._validate_field(potential, "potensial")
        flux = np.asarray(self.coefficients.T_cell @ phi).ravel()
        return flux + self._boundary_contribution(self.coefficients.T_bnd,
                                                  boundary_potential)

    def net_influx(self, potential: np.ndarray,
                   mobility: Optional[np.ndarray] = None,
                   boundary_potential: Optional[np.ndarray] = None) -> np.ndarray:
        """`(ncell,)` — hüceyrəyə DAXİL OLAN xalis axın.

        `mobility` verilibsə (əlaqə üzrə, artıq upstream seçilmiş),
        axın ona vurulur. `D @ (M ⊙ q)` — `ResidualAssembler.net_influx`
        ilə RİYAZİ olaraq EYNİ (§8 sahiblik konvensiyası)."""
        flux = self.connection_fluxes(potential, boundary_potential)
        if mobility is not None:
            mobility = np.asarray(mobility, float)
            if mobility.shape != (self.nconn,):
                raise MPFAStateError(f"mobility ({self.nconn},) olmalıdır, "
                                     f"alındı {mobility.shape}")
            flux = flux * mobility
        return np.asarray(self._divergence @ flux).ravel()

    def upstream_cells(self, flux: np.ndarray) -> np.ndarray:
        """`(nconn,)` — HƏQİQİ çoxnöqtəli axının işarəsindən upstream
        hüceyrəsi (tapşırıq §10).

        `q ≥ 0` → `cell_a` (owner), `q < 0` → `cell_b` (neighbor).
        İki-nöqtəli limitdə `q = T·ΔΦ`, `T > 0` olduğu üçün bu, TPFA-nın
        `np.where(ΔΦ >= 0, cell_a, cell_b)` qaydası ilə BİRƏBİR eynidir.
        """
        flux = np.asarray(flux, float)
        return np.where(flux >= 0.0, self.connections.cell_a, self.connections.cell_b)

    # ──────────────────────────────── qlobal operator / naxış (§16/§22)
    def cell_operator(self) -> sparse.csr_matrix:
        """`A = D · T_conn` `(ncell × ncell)` — SABİT vahid mobilitəli
        halda `net_influx = A · Φ`.

        **BU, QEYRİ-XƏTTİ REZERVUAR JACOBIAN-I DEYİL** (tapşırıq §23):
        mobilitə, upstream seçimi, akkumulyasiya, quyu törəmələri
        BURADA YOXDUR. Yalnız MPFA-nın XƏTTİ məkan operatorudur.
        Keşlənir (dövlətdən asılı deyil)."""
        if self._cell_operator is None:
            self._cell_operator = (self._divergence @ self.T_conn).tocsr()
            self._cell_operator.sum_duplicates()
        return self._cell_operator

    def global_stencil_pattern(self) -> sparse.csr_matrix:
        """HƏQİQİ MPFA bağlantı naxışı `(ncell × ncell)`, `bool`.

        `Connections` qonşuluğundan GENİŞDİR — MPFA bölgəsi diaqonal
        (üz paylaşmayan) hüceyrələri də bağlayır (tapşırıq §22).
        `JacobianAssembler` Phase 5B-2-də məhz BU naxışa keçməlidir."""
        matrix = self.cell_operator().copy()
        matrix.data = np.abs(matrix.data) > 0.0
        matrix.eliminate_zeros()
        return matrix.astype(bool)

    def face_stencil(self, connection: int, tolerance: float = 0.0) -> dict:
        """`{hüceyrə: əmsal}` — `connection` üzündəki axının HANSI hüceyrə
        potensiallarından asılı olduğu (§16 auditi üçün)."""
        row = self.T_conn.getrow(connection).toarray().ravel()
        return {int(c): float(row[c])
                for c in np.flatnonzero(np.abs(row) > tolerance)}

    def cell_stencil(self, cell: int, tolerance: float = 0.0) -> dict:
        """`{hüceyrə: əmsal}` — `A` matrisinin `cell` sətri (§30 auditi)."""
        row = self.cell_operator().getrow(cell).toarray().ravel()
        return {int(c): float(row[c])
                for c in np.flatnonzero(np.abs(row) > tolerance)}

    def stencil_sizes(self, tolerance: float = 1e-12) -> np.ndarray:
        """`(nconn,)` — hər daxili üzün qlobal stensil ölçüsü."""
        matrix = self.T_conn.tocsr()
        sizes = np.zeros(matrix.shape[0], dtype=int)
        for row in range(matrix.shape[0]):
            values = np.abs(matrix.data[matrix.indptr[row]:matrix.indptr[row + 1]])
            if values.size:
                sizes[row] = int(np.sum(values > tolerance * max(values.max(), 1e-300)))
        return sizes

    # ───────────────────────────────────────────── konservasiya (§7/§8)
    def conservation_report(self, potential: np.ndarray,
                            boundary_potential: Optional[np.ndarray] = None) -> dict:
        """Daxili üz və qlobal balans qalıqları (tapşırıq §7/§D).

        `max_internal_face_error` — MÜSTƏQİL yolla hesablanır: hər
        qarşılıqlı təsir bölgəsində sub-üzün HƏR İKİ tərəfi ÖZ
        `g_(c,v)` qradiyentindən ayrıca qurulur (Phase 5A
        `MPFAOLocalSystem.conservation_residual`) — yəni bu, "axını bir
        dəfə hesabladıq, deməli qorunur" tavtologiyası DEYİL.
        `max_boundary_flux` — `NEUMANN_ZERO`-da sıfır OLMALIDIR; sıfırdan
        fərqlidirsə sərhəd payı qlobal balansdan İTİRİLİR (§8).
        `global_imbalance` — `Σ_i net_influx_i`, qapalı sistemdə ≡ 0.
        """
        phi = self._validate_field(potential, "potensial")
        face_flux = self.face_fluxes(phi, boundary_potential)

        internal = np.zeros(face_flux.size, dtype=bool)
        internal[self.connection_faces] = True
        return {
            "max_internal_face_error": _region_conservation(
                self.coefficients, phi, boundary_potential, self._boundary_index),
            "max_boundary_flux": (float(np.abs(face_flux[~internal]).max())
                                  if np.any(~internal) else 0.0),
            "global_imbalance": float(np.abs(np.sum(
                self.net_influx(phi, boundary_potential=boundary_potential)))),
            "max_face_flux": float(np.abs(face_flux).max()) if face_flux.size else 0.0,
        }

    # ───────────────────────────────────────────────────── validasiya
    def _validate_field(self, field: np.ndarray, label: str) -> np.ndarray:
        """Tapşırıq §27 — NaN/Inf və ölçü uyğunsuzluğu AÇIQ rədd edilir,
        problemli HÜCEYRƏLƏR göstərilir. SƏSSİZ təmir YOXDUR."""
        array = np.asarray(field, float)
        if array.shape != (self.ncell,):
            raise MPFAStateError(
                f"{label} vektoru ({self.ncell},) olmalıdır, alındı {array.shape} "
                "— MPFA stensil ölçüsü ilə uyğunsuzluq.")
        bad = np.flatnonzero(~np.isfinite(array))
        if bad.size:
            raise MPFAStateError(
                f"{label} vektorunda {bad.size} hüceyrədə NaN/Inf var "
                f"(ilk hüceyrələr: {bad[:5].tolist()}) — MPFA bunu TƏMİR ETMİR.")
        return array


def _region_conservation(coefficients: "MPFAOCoefficients", potential: np.ndarray,
                         boundary_potential: Optional[np.ndarray],
                         boundary_index) -> float:
    """Bölgə-daxili konservasiya qalığı — hər tərəf ÖZ `g_(c,v)`
    qradiyentindən MÜSTƏQİL hesablanır (Phase 5A
    `MPFAOLocalSystem.conservation_residual`)."""
    worst = 0.0
    for system in coefficients.local_systems:
        local_boundary = None
        if system.known_boundary_sub_faces:
            if boundary_potential is None:
                raise MPFAStateError(
                    "Dirichlet bağlanışında konservasiya hesabatı üçün "
                    "`boundary_potential` TƏLƏB OLUNUR.")
            local_boundary = np.asarray(boundary_potential, float)[
                [boundary_index[(system.region.sub_faces[s].face_index,
                                 system.region.node_id)]
                 for s in system.known_boundary_sub_faces]]
        worst = max(worst, system.conservation_residual(
            potential[system.region.cells], local_boundary))
    return float(worst)
