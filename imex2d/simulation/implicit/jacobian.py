"""Analitik Jakobian — A6, mərhələ 2.

J[i, j] = ∂R_i / ∂x_j,  ölçü 2N × 2N.

Sətir indeksi:  2c + faza      (faza: 0 = su, 1 = neft)
Sütun indeksi:  2c + dəyişən   (dəyişən: 0 = təzyiq, 1 = Sw)

Üç töhfə mənbəyi var:

1. AKKUMULYASİYA — yalnız diaqonal bloka
       ∂/∂p  [PV·Sw/Bw]   = −PV·Sw·B'w / B²w
       ∂/∂Sw [PV·Sw/Bw]   =  PV / Bw
       ∂/∂p  [PV·So/Bo]   = −PV·So·B'o / B²o
       ∂/∂Sw [PV·So/Bo]   = −PV / Bo

2. ÜZLƏR ÜZRƏ AXIN — həm diaqonal, həm qonşu blokuna
       F_p = T · M_p,upstream · ΔΦ_p,      M_p = kr_p / (μ_p·B_p)
       ∂F/∂x = T · [M · ∂ΔΦ/∂x  +  ΔΦ · ∂M/∂x]
   İkinci üzv YALNIZ upstream hüceyrəsinə görə sıfırdan fərqlidir.

3. QUYULAR — yalnız diaqonal bloka

STANDART SADƏLƏŞDİRMƏLƏR (sənaye simulyatorlarında da belədir):
  · Upstream seçiminin özü diferensiallaşdırılmır (ΔΦ işarəsi
    dəyişəndə funksiya kəsilməzdir, törəməsi isə sıçrayır).
  · Cazibə üzvündə sıxlığın təzyiqdən asılılığı nəzərə alınmır.
Hər ikisi Nyutonun yığılma sürətinə cüzi təsir edir, doğruluğuna yox —
çünki konvergensiya qalığa görə yoxlanılır, Jakobiana görə yox.

PHASE 5B-2 — MPFA-O ANALİTİK JACOBIAN (bu fayl artıq DƏYİŞDİRİLİB):

`_flux()`/`_build_pattern()` (TPFA, YUXARIDA, BİR SƏTİR BELƏ DƏYİŞMƏDİ)
`self.R.transmissibility`-ni BİRBAŞA oxuyur, ∂ΔΦ/∂p_a=+1, ∂ΔΦ/∂p_b=−1
tək-cüt fərziyyəsi ilə. `ResidualAssembler._multipoint=True` (MPFA-O)
olanda ƏVƏZİNƏ `_flux_multipoint()`/`_build_pattern_diag_only()` işə
düşür (bax aşağıda, "PHASE 5B-2" bölməsi) — bu, `docs/mpfa_o_phase5b1.md`
§11-in "qalan" siyahısını tamamlayır: `∂M_α/∂p`, `∂M_α/∂S_w`, `∂Φ/∂S_w`
(kapilyar) İNDİ çoxnöqtəli stensillə hesablanır, `T_conn`-un (Phase 5B-1,
`MPFAGlobalOperator`) HƏNDƏSƏ/K-dan asılı, dövlətdən ASILI OLMAYAN seyrək
matrisindən İSTİFADƏ edərək — heç bir lokal MPFA sistemi TƏKRAR
QURULMUR, heç bir sonlu-fərq QISAYOLU yoxdur (tapşırıq §23/§26 — analitik
düstur aşağıda İKİ-NÖQTƏLİ HƏDDƏ endirilib və TPFA-nın öz `_flux()`-u ilə
ƏL İLƏ yoxlanıb, bax `tests/test_mpfa_jacobian.py`).

Riyazi düstur (su fazası üçün, neft analojidir, `Pc` üzvü YOXDUR):

    Φ_w,j = p_j − Pc_j(Sw_j) − ρ_w,j(p_j)·g·D_j·PA_TO_BAR    (HÜCEYRƏ, `cell_potentials`)
    q_pot,f = Σ_j T_conn[f,j]·Φ_w,j                            (T_conn — Phase 5B-1, DÖVLƏTSİZ)
    u(f) = upstream(q_pot,f)                                   (TƏK hüceyrə, DİFERENSİALLAŞMIR — TPFA-dakı EYNİ sadələşdirmə)
    q_w,f = q_pot,f · M_w[u(f)],    M_w = λ_w/B_w

    ∂q_w,f/∂p_j  = M_w[u(f)]·T_conn[f,j]  +  q_pot,f·∂M_w/∂p[u(f)]·[j=u(f)]
    ∂q_w,f/∂Sw_j = M_w[u(f)]·T_conn[f,j]·(−∂Pc_j/∂Sw_j)  +  q_pot,f·∂M_w/∂Sw[u(f)]·[j=u(f)]

Birinci həd (BAZA, `T·M`) — bütün stensil sütunlarına; ikinci həd
(UPSTREAM mobilitə törəməsi, TPFA-dakı `ΔΦ·∂M/∂p[a]` üzvünün analoqu) —
YALNIZ `u(f)` sütununa. İkiqat oxşarlıq (`D_R[a,f]=+1,D_R[b,f]=−1`
işarəsi ilə) İKİ-NÖQTƏLİ HƏDDƏ (`T_conn[f,a]=T,T_conn[f,b]=−T`) qoyulanda
TPFA-nın mövcud `_flux()` düsturunu BİRƏBİR (4 element: aa/ab/ba/bb)
bərpa edir — `tests/test_mpfa_jacobian.py::test_two_point_limit_matches_tpfa_jacobian_exactly`.
"""

from __future__ import annotations
from typing import Optional

import numpy as np
import scipy.sparse as sp

from ...domain.wells import ControlMode
from .derivatives import DerivativeProvider
from .residual import OIL, WATER, FluidState, ResidualAssembler
from .state import (PRESSURE, VARIABLES_PER_CELL, WATER_SATURATION,
                    ReservoirState)


class JacobianAssembler:
    """Qalıq vektoru ilə eyni diskretizasiyadan Jakobian qurur."""

    def __init__(self, residual_assembler: ResidualAssembler,
                 derivatives: Optional[DerivativeProvider] = None):
        self.R = residual_assembler
        self.model = residual_assembler.model
        self.ncell = residual_assembler.ncell
        self.size = self.ncell * VARIABLES_PER_CELL
        self.derivatives = derivatives or DerivativeProvider(
            residual_assembler.relperm,
            pvt=residual_assembler.pvt,
            capillary=residual_assembler.capillary,
            fluids=self.model.fluids,
            reference_pressure=residual_assembler.reference_pressure)
        #: ÇOXNÖQTƏLİ (MPFA-O) yolu — `ResidualAssembler`-in ÖZ (duck-typed)
        #: bayrağı ilə EYNİ (bax `residual.py`). Phase 5B-2: bu ARTIQ AÇIQ
        #: rədd EDİLMİR — `_build_pattern_diag_only()` + `_flux_multipoint()`
        #: işə düşür (aşağıda). TPFA yolu (`_multipoint=False`) BİR SƏTİR
        #: BELƏ dəyişmədi.
        self._multipoint = bool(getattr(residual_assembler, "_multipoint", False))
        if self._multipoint:
            self._build_pattern_diag_only()
        else:
            self._build_pattern()

    # ═════════════════════════════════════════════ struktur (bir dəfə)
    def _build_pattern(self):
        """Seyrəklik strukturu simulyasiya boyu sabitdir — bir dəfə qurulur."""
        conn = self.R.connections
        cells = np.arange(self.ncell)
        a, b = conn.cell_a, conn.cell_b

        rows, cols = [], []

        def block(row_cells, col_cells):
            """Hər (hüceyrə cütü) üçün 2×2 blokun dörd mövqeyi."""
            for phase in (WATER, OIL):
                for variable in (PRESSURE, WATER_SATURATION):
                    rows.append(row_cells * VARIABLES_PER_CELL + phase)
                    cols.append(col_cells * VARIABLES_PER_CELL + variable)

        block(cells, cells)      # diaqonal: akkumulyasiya + quyular + axın
        block(a, a)              # üz: R_a ↔ x_a
        block(a, b)              # üz: R_a ↔ x_b
        block(b, a)              # üz: R_b ↔ x_a
        block(b, b)              # üz: R_b ↔ x_b

        self._rows = np.concatenate(rows)
        self._cols = np.concatenate(cols)
        self._values = np.zeros(self._rows.size)

        entries = np.arange(self._rows.size, dtype=np.int64)
        pattern = sp.coo_matrix((entries + 1.0, (self._rows, self._cols)),
                                shape=(self.size, self.size)).tocsr()
        pattern.sum_duplicates()
        self._matrix = sp.csr_matrix(
            (np.zeros(pattern.nnz), pattern.indices.copy(),
             pattern.indptr.copy()), shape=(self.size, self.size))

        locator = sp.csr_matrix(
            (np.arange(pattern.nnz, dtype=np.int64), pattern.indices.copy(),
             pattern.indptr.copy()), shape=(self.size, self.size))
        self._data_index = np.asarray(
            locator[self._rows, self._cols]).ravel().astype(np.int64)
        self._nnz = pattern.nnz

        # blokların `_values` massivindəki başlanğıc mövqeləri
        n_cell_entries = self.ncell
        n_face_entries = conn.count
        self._slices = {}
        offset = 0
        for name, count in (("diag", n_cell_entries), ("aa", n_face_entries),
                            ("ab", n_face_entries), ("ba", n_face_entries),
                            ("bb", n_face_entries)):
            for phase in (WATER, OIL):
                for variable in (PRESSURE, WATER_SATURATION):
                    self._slices[(name, phase, variable)] = slice(offset,
                                                                  offset + count)
                    offset += count

    # ═════════════════════════════ struktur (bir dəfə) — PHASE 5B-2: MPFA-O
    def _build_pattern_diag_only(self):
        """`_build_pattern()`-in DİAQONAL-YALNIZ variantı — akkumulyasiya
        VƏ quyular (`_accumulation`/`_wells`, BİR SƏTİR BELƏ DƏYİŞMƏDƏN,
        MPFA-da da hüceyrə-lokaldır) üçün EYNİ `_set`/`_add` mexanizmini
        işlədir. Axın (`_flux_multipoint`) BU naxışdan KƏNARDA, AYRI seyrək
        matris kimi qurulur (çoxnöqtəli stensil — sabit 2-hüceyrəli blok
        fərziyyəsi YOXDUR) və `assemble()`-də bu matrisə ƏLAVƏ OLUNUR."""
        cells = np.arange(self.ncell)
        rows, cols = [], []
        for phase in (WATER, OIL):
            for variable in (PRESSURE, WATER_SATURATION):
                rows.append(cells * VARIABLES_PER_CELL + phase)
                cols.append(cells * VARIABLES_PER_CELL + variable)
        self._rows = np.concatenate(rows)
        self._cols = np.concatenate(cols)
        self._values = np.zeros(self._rows.size)

        entries = np.arange(self._rows.size, dtype=np.int64)
        pattern = sp.coo_matrix((entries + 1.0, (self._rows, self._cols)),
                                shape=(self.size, self.size)).tocsr()
        pattern.sum_duplicates()
        self._matrix = sp.csr_matrix(
            (np.zeros(pattern.nnz), pattern.indices.copy(),
             pattern.indptr.copy()), shape=(self.size, self.size))
        locator = sp.csr_matrix(
            (np.arange(pattern.nnz, dtype=np.int64), pattern.indices.copy(),
             pattern.indptr.copy()), shape=(self.size, self.size))
        self._data_index = np.asarray(
            locator[self._rows, self._cols]).ravel().astype(np.int64)
        self._nnz = pattern.nnz

        self._slices = {}
        offset = 0
        for phase in (WATER, OIL):
            for variable in (PRESSURE, WATER_SATURATION):
                self._slices[("diag", phase, variable)] = slice(offset, offset + self.ncell)
                offset += self.ncell

    def _set(self, name, phase, variable, values):
        self._values[self._slices[(name, phase, variable)]] = values

    def _add(self, name, phase, variable, values):
        self._values[self._slices[(name, phase, variable)]] += values

    # ═══════════════════════════════════════════════════════ yığım
    def assemble(self, state: ReservoirState, fluid: FluidState,
                 dt: float) -> sp.csr_matrix:
        if self._multipoint:
            return self._assemble_multipoint(state, fluid, dt)

        self._values[:] = 0.0
        self._accumulation(state, fluid, dt)
        self._flux(state, fluid)
        self._wells(state, fluid)

        self._matrix.data[:] = np.bincount(self._data_index,
                                           weights=self._values,
                                           minlength=self._nnz)
        return self._matrix

    # ═══════════════════════════════════ yığım (PHASE 5B-2: MPFA-O) ─────
    def _assemble_multipoint(self, state: ReservoirState, fluid: FluidState,
                             dt: float) -> sp.csr_matrix:
        """Akkumulyasiya/quyular (diaqonal, `_accumulation`/`_wells`
        DƏYİŞMƏDƏN) + çoxnöqtəli axın Jacobian-ı (`_flux_multipoint`,
        AYRI seyrək matris) — ƏLAVƏ olaraq CSR toplanır."""
        self._values[:] = 0.0
        self._accumulation(state, fluid, dt)
        self._wells(state, fluid)
        self._matrix.data[:] = np.bincount(self._data_index,
                                           weights=self._values,
                                           minlength=self._nnz)
        flux_matrix = self._flux_multipoint(state, fluid)
        return (self._matrix + flux_matrix).tocsr()

    # ─────────────────────────────────────────────── akkumulyasiya
    def _accumulation(self, state: ReservoirState, fluid: FluidState,
                      dt: float) -> None:
        """∂/∂x [PV(p)·S_p/B_p(p)] / Δt.

        Həm məsamə həcmi, həm formasiya həcm əmsalı təzyiqdən asılıdır:

            ∂/∂p [PV·S/B] = S · [PV'·B − PV·B'] / B²
            ∂/∂Sw[PV·Sw/Bw] = +PV/Bw,   ∂/∂Sw[PV·So/Bo] = −PV/Bo
        """
        pore_volume = self.R.pore_volume_at(state.pressure) / dt
        d_pore_volume = (self.R.pore_volume
                         * self.model.rock.compressibility / dt)
        sw = state.water_saturation
        so = 1.0 - sw
        dbw = self.derivatives.dbw_dp(state.pressure)
        dbo = self.derivatives.dbo_dp(state.pressure)

        self._add("diag", WATER, PRESSURE,
                  sw * (d_pore_volume * fluid.bw - pore_volume * dbw)
                  / fluid.bw ** 2)
        self._add("diag", WATER, WATER_SATURATION, pore_volume / fluid.bw)
        self._add("diag", OIL, PRESSURE,
                  so * (d_pore_volume * fluid.bo - pore_volume * dbo)
                  / fluid.bo ** 2)
        self._add("diag", OIL, WATER_SATURATION, -pore_volume / fluid.bo)

    # ────────────────────────────────────────────────────── axın
    def _flux(self, state: ReservoirState, fluid: FluidState) -> None:
        conn = self.R.connections
        a, b = conn.cell_a, conn.cell_b
        transmissibility = self.R.transmissibility

        d_phi_w, d_phi_o = self.R.potentials(state, fluid)
        upstream_w = np.where(d_phi_w >= 0, a, b)
        upstream_o = np.where(d_phi_o >= 0, a, b)

        mobility_w = fluid.lam_w / fluid.bw
        mobility_o = fluid.lam_o / fluid.bo
        dmw_dp, dmw_dsw = self.derivatives.water_transport_derivatives(
            state.pressure, state.water_saturation, fluid)
        dmo_dp, dmo_dsw = self.derivatives.oil_transport_derivatives(
            state.pressure, state.water_saturation, fluid)
        dpc = self.derivatives.dpc_dsw(state.water_saturation)

        upstream_is_a_w = upstream_w == a
        upstream_is_a_o = upstream_o == a

        # ---- su fazası -------------------------------------------------
        m_w = mobility_w[upstream_w]
        # ∂ΔΦw/∂p_a = +1, ∂ΔΦw/∂p_b = −1
        dfw_dpa = transmissibility * (m_w + d_phi_w * np.where(
            upstream_is_a_w, dmw_dp[a], 0.0))
        dfw_dpb = transmissibility * (-m_w + d_phi_w * np.where(
            upstream_is_a_w, 0.0, dmw_dp[b]))
        # ∂ΔΦw/∂Sw_a = −dPc_a/dSw,  ∂ΔΦw/∂Sw_b = +dPc_b/dSw
        dfw_dsa = transmissibility * (-m_w * dpc[a] + d_phi_w * np.where(
            upstream_is_a_w, dmw_dsw[a], 0.0))
        dfw_dsb = transmissibility * (m_w * dpc[b] + d_phi_w * np.where(
            upstream_is_a_w, 0.0, dmw_dsw[b]))

        # ---- neft fazası -----------------------------------------------
        m_o = mobility_o[upstream_o]
        dfo_dpa = transmissibility * (m_o + d_phi_o * np.where(
            upstream_is_a_o, dmo_dp[a], 0.0))
        dfo_dpb = transmissibility * (-m_o + d_phi_o * np.where(
            upstream_is_a_o, 0.0, dmo_dp[b]))
        dfo_dsa = transmissibility * (d_phi_o * np.where(
            upstream_is_a_o, dmo_dsw[a], 0.0))
        dfo_dsb = transmissibility * (d_phi_o * np.where(
            upstream_is_a_o, 0.0, dmo_dsw[b]))

        # R_a = … + F,  R_b = … − F
        for phase, dpa, dpb, dsa, dsb in (
                (WATER, dfw_dpa, dfw_dpb, dfw_dsa, dfw_dsb),
                (OIL, dfo_dpa, dfo_dpb, dfo_dsa, dfo_dsb)):
            self._set("aa", phase, PRESSURE, dpa)
            self._set("aa", phase, WATER_SATURATION, dsa)
            self._set("ab", phase, PRESSURE, dpb)
            self._set("ab", phase, WATER_SATURATION, dsb)
            self._set("ba", phase, PRESSURE, -dpa)
            self._set("ba", phase, WATER_SATURATION, -dsa)
            self._set("bb", phase, PRESSURE, -dpb)
            self._set("bb", phase, WATER_SATURATION, -dsb)

    # ─────────────────────────────── axın (PHASE 5B-2: MPFA-O, ÇOXNÖQTƏLİ)
    def _flux_multipoint(self, state: ReservoirState, fluid: FluidState) -> sp.csr_matrix:
        """MPFA-O axınının `(2N × 2N)` Jacobian töhfəsi — bax modul
        docstring-i (düstur + iki-nöqtəli həddə TPFA-ya bərabərliyin
        izahı). `T_conn` (Phase 5B-1, `MPFAGlobalOperator`) HƏNDƏSƏ/K-dan
        asılı, DÖVLƏTDƏN asılı DEYİL — burada YALNIZ oxunur, YENİDƏN
        QURULMUR."""
        grid = self.R.grid
        operator = grid.global_operator
        conn = self.R.connections
        a, b = conn.cell_a, conn.cell_b
        nconn = conn.count
        T_conn = operator.T_conn

        phi_w, phi_o = self.R.cell_potentials(state, fluid)
        q_pot_w = np.asarray(T_conn @ phi_w).ravel()
        q_pot_o = np.asarray(T_conn @ phi_o).ravel()
        u_w = grid.upstream_cells(q_pot_w)
        u_o = grid.upstream_cells(q_pot_o)

        mobility_w = fluid.lam_w / fluid.bw
        mobility_o = fluid.lam_o / fluid.bo
        dmw_dp, dmw_dsw = self.derivatives.water_transport_derivatives(
            state.pressure, state.water_saturation, fluid)
        dmo_dp, dmo_dsw = self.derivatives.oil_transport_derivatives(
            state.pressure, state.water_saturation, fluid)
        dpc = self.derivatives.dpc_dsw(state.water_saturation)

        # D_R: (ncell × nconn) — `cell_a`: +1, `cell_b`: −1. Bu, `_flux()`-un
        # ("aa"=+dpa, "ba"=−dpa) mövcud işarə konvensiyası ilə EYNİDİR
        # (bax modul docstring-i, iki-nöqtəli həddə yoxlanılıb) — `Residual
        # Assembler.net_influx`-un `np.add.at` işarəsi (`cell_a`: −, `cell_b`:
        # +) ilə QARIŞDIRILMASIN: `R = acc/dt − influx − q`, ona görə
        # `∂R/∂flux` = `−∂influx/∂flux`, işarə burada ARTIQ HESABA ALINIB.
        idx = np.arange(nconn)
        d_rows = np.concatenate([a, b])
        d_cols = np.concatenate([idx, idx])
        d_vals = np.concatenate([np.ones(nconn), -np.ones(nconn)])
        D_R = sp.coo_matrix((d_vals, (d_rows, d_cols)),
                            shape=(self.ncell, nconn)).tocsr()

        def base(mobility_at_upstream: np.ndarray,
                column_scale: Optional[np.ndarray] = None) -> sp.csr_matrix:
            """`D_R · diag(M[u]) · T_conn[· diag(column_scale)]` — BAZA
            (`T·M`) həddi, bütün stensil sütunlarına."""
            matrix = T_conn if column_scale is None else T_conn @ sp.diags(column_scale)
            scaled = sp.diags(mobility_at_upstream) @ matrix
            return (D_R @ scaled).tocsr()

        def extra(q_pot: np.ndarray, upstream: np.ndarray,
                 d_mobility_at_upstream: np.ndarray) -> sp.csr_matrix:
            """UPSTREAM mobilitə-törəməsi həddi — YALNIZ `u(f)` sütununa
            (TPFA-dakı `ΔΦ·∂M/∂x[upstream]` üzvünün çoxnöqtəli analoqu)."""
            weighted = q_pot * d_mobility_at_upstream
            rows = np.concatenate([a, b])
            cols = np.concatenate([upstream, upstream])
            values = np.concatenate([weighted, -weighted])
            return sp.coo_matrix((values, (rows, cols)),
                                 shape=(self.ncell, self.ncell)).tocsr()

        j_water_p = base(mobility_w[u_w]) + extra(q_pot_w, u_w, dmw_dp[u_w])
        j_water_s = (base(mobility_w[u_w], column_scale=-dpc)
                    + extra(q_pot_w, u_w, dmw_dsw[u_w]))
        j_oil_p = base(mobility_o[u_o]) + extra(q_pot_o, u_o, dmo_dp[u_o])
        j_oil_s = extra(q_pot_o, u_o, dmo_dsw[u_o])   # Φ_o Sw-dan asılı deyil → baza yoxdur

        def place(block: sp.csr_matrix, phase: int, variable: int) -> sp.coo_matrix:
            coo = block.tocoo()
            rows = coo.row * VARIABLES_PER_CELL + phase
            cols = coo.col * VARIABLES_PER_CELL + variable
            return sp.coo_matrix((coo.data, (rows, cols)), shape=(self.size, self.size))

        total = (place(j_water_p, WATER, PRESSURE) + place(j_water_s, WATER, WATER_SATURATION)
                + place(j_oil_p, OIL, PRESSURE) + place(j_oil_s, OIL, WATER_SATURATION))
        return total.tocsr()

    # ───────────────────────────────────────────────────── quyular
    def _wells(self, state: ReservoirState, fluid: FluidState) -> None:
        """R = … − q,  ona görə ∂R/∂x = −∂q/∂x."""
        pressure = state.pressure
        sw = state.water_saturation
        endpoint = self.R.relperm.endpoint_water_mobility(1.0)

        dmuw = self.derivatives.dmuw_dp(pressure)
        dmuo = self.derivatives.dmuo_dp(pressure)
        dbw = self.derivatives.dbw_dp(pressure)
        dbo = self.derivatives.dbo_dp(pressure)
        dmw_dp, dmw_dsw = self.derivatives.water_transport_derivatives(
            pressure, sw, fluid)
        dmo_dp, dmo_dsw = self.derivatives.oil_transport_derivatives(
            pressure, sw, fluid)

        water_p = np.zeros(self.ncell)
        water_s = np.zeros(self.ncell)
        oil_p = np.zeros(self.ncell)
        oil_s = np.zeros(self.ncell)

        for connection in self.R.wells:
            c = connection.cell
            wi = connection.well_index

            if connection.is_injector:
                # q = WI · (krw_end/μw) · (target − p) / Bw,  q = max(q, 0)
                transport = endpoint / (fluid.mu_w[c] * fluid.bw[c])
                drawdown = connection.target - pressure[c]
                if connection.mode is ControlMode.BHP:
                    if wi * transport * drawdown <= 0.0:
                        continue                      # kəsilmiş (q = 0)
                    d_transport = -endpoint * (
                        dmuw[c] * fluid.bw[c] + fluid.mu_w[c] * dbw[c]) / (
                        (fluid.mu_w[c] * fluid.bw[c]) ** 2)
                    water_p[c] += wi * (-transport + drawdown * d_transport)
                else:
                    rate = abs(connection.target)
                    water_p[c] += -rate * dbw[c] / fluid.bw[c] ** 2
                continue

            # hasilat
            if connection.mode is ControlMode.BHP:
                drawdown = connection.target - pressure[c]
                mobility_w = fluid.lam_w[c] / fluid.bw[c]
                mobility_o = fluid.lam_o[c] / fluid.bo[c]
                if wi * mobility_w * drawdown < 0.0:
                    water_p[c] += wi * (-mobility_w + drawdown * dmw_dp[c])
                    water_s[c] += wi * drawdown * dmw_dsw[c]
                if wi * mobility_o * drawdown < 0.0:
                    oil_p[c] += wi * (-mobility_o + drawdown * dmo_dp[c])
                    oil_s[c] += wi * drawdown * dmo_dsw[c]
            else:
                # q_total sabit, fazalara fraksional axınla bölünür:
                #     f = λw / (λw + λo),   λp = kr_p / μ_p(p)
                # f həm Sw-dan, həm də TƏZYİQDƏN asılıdır, çünki lözlüklər
                # PVT ilə təzyiqə bağlıdır.
                total = -abs(connection.target)
                lam_w, lam_o = fluid.lam_w[c], fluid.lam_o[c]
                lam_t = max(lam_w + lam_o, 1e-30)
                fraction = lam_w / lam_t

                krw = float(self.R.relperm.krw(np.array([sw[c]]))[0])
                kro = float(self.R.relperm.kro(np.array([sw[c]]))[0])

                dlw_dsw = self.derivatives.dkrw_dsw(np.array([sw[c]]))[0] / fluid.mu_w[c]
                dlo_dsw = self.derivatives.dkro_dsw(np.array([sw[c]]))[0] / fluid.mu_o[c]
                dfraction_dsw = (dlw_dsw * lam_o - lam_w * dlo_dsw) / lam_t ** 2

                dlw_dp = -krw * dmuw[c] / fluid.mu_w[c] ** 2
                dlo_dp = -kro * dmuo[c] / fluid.mu_o[c] ** 2
                dfraction_dp = (dlw_dp * lam_o - lam_w * dlo_dp) / lam_t ** 2

                water_s[c] += total * dfraction_dsw / fluid.bw[c]
                oil_s[c] += -total * dfraction_dsw / fluid.bo[c]
                water_p[c] += total * (dfraction_dp / fluid.bw[c]
                                       - fraction * dbw[c] / fluid.bw[c] ** 2)
                oil_p[c] += total * (-dfraction_dp / fluid.bo[c]
                                     - (1.0 - fraction) * dbo[c] / fluid.bo[c] ** 2)

        self._add("diag", WATER, PRESSURE, -water_p)
        self._add("diag", WATER, WATER_SATURATION, -water_s)
        self._add("diag", OIL, PRESSURE, -oil_p)
        self._add("diag", OIL, WATER_SATURATION, -oil_s)

    # ═════════════════════════════════════════ yoxlama (sonlu fərq)
    def numerical(self, state: ReservoirState, previous: ReservoirState,
                  dt: float, step: float = 1e-6) -> np.ndarray:
        """Sonlu fərqlə Jakobian — YALNIZ testlər üçün, O(N²) bahalıdır."""
        base = state.to_vector()
        result = np.zeros((self.size, self.size))
        for column in range(self.size):
            scale = step * max(abs(base[column]), 1.0)
            forward, backward = base.copy(), base.copy()
            forward[column] += scale
            backward[column] -= scale
            r_forward, _, _ = self.R.residual(
                ReservoirState.from_vector(forward), previous, dt)
            r_backward, _, _ = self.R.residual(
                ReservoirState.from_vector(backward), previous, dt)
            result[:, column] = (r_forward - r_backward) / (2.0 * scale)
        return result
