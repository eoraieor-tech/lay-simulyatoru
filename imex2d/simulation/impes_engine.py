"""IMPES hesablama mühərriki.

DƏYİŞƏN: mühərrik artıq öz modelini QURMUR. Konstruktora hazır
ReservoirModel, konfiqurasiya və provider-lər verilir.

DƏYİŞMƏYƏN: riyaziyyat. Təzyiq implicit, doyumluluq explicit, upstream
çəkilənmə, CFL nəzarəti — hamısı əvvəlki nüvədən olduğu kimi köçürülüb.
Bu, refaktorinqin nəticəni dəyişmədiyini yoxlamağa imkan verir.
"""

from __future__ import annotations
from typing import List, Optional

import numpy as np
import scipy.sparse as sp

from ..application.config import SimulationConfig
from ..domain.reservoir_model import ReservoirModel
from ..domain.wells import ControlMode
from ..interfaces.providers import (ICapillaryPressureProvider,
                                    IInitializationProvider, IPVTProvider,
                                    IRelativePermeabilityProvider)
from ..interfaces.discretization import IFluxDiscretization
from ..logging_setup import get_logger
from ..interfaces.services import (ILinearSolver, IProgressReporter,
                                   ISimulationEngine, NullProgressReporter)
from .discretization import default_flux_discretization

LOG = get_logger(__name__)

GRAVITY = 9.80665
PA_TO_BAR = 1.0e-5
from .results import SimulationResult, Snapshot
from .well_model import PeacemanWellModel, WellConnection


def _reject_multipoint_impes(discretization: IFluxDiscretization) -> None:
    """ÇOXNÖQTƏLİ (MPFA-O) diskretizasiya ilə AÇIQ imtina — Phase 5B-2.

    IMPES `_discretization.transmissibility`-ni BİRBAŞA oxuyur (tək-üz
    skalyar). MPFA-O-da belə kəmiyyət YOXDUR; saxta bir dəyər uydurmaq
    metodu TPFA-ya çevirmək olardı (tapşırıq §24/§32).
    """
    if getattr(discretization, "supports_multipoint_stencil", lambda: False)():
        raise NotImplementedError(
            "ImpesEngine çoxnöqtəli diskretizasiya (MPFA-O) ilə HƏLƏ İŞLƏMİR "
            "(Phase 5B-2): IMPES tək-üz transmissivliyi tələb edir, MPFA-O-da "
            "isə belə kəmiyyət yoxdur. Bax docs/mpfa_o_phase5b1.md §11.")


class ImpesEngine(ISimulationEngine):

    def __init__(self,
                 model: ReservoirModel,
                 config: SimulationConfig,
                 relperm: IRelativePermeabilityProvider,
                 linear_solver: ILinearSolver,
                 pvt: Optional[IPVTProvider] = None,
                 capillary: Optional[ICapillaryPressureProvider] = None,
                 initialization: Optional[IInitializationProvider] = None,
                 flux_discretization: Optional[IFluxDiscretization] = None):
        self.model = model
        self.config = config
        self.relperm = relperm
        self.linear_solver = linear_solver
        self.pvt = pvt
        self.capillary = capillary
        self.initialization = initialization

        #: DEFOLT = TPFA (bax audit tapşırığı §4, `FullyImplicitEngine`-də
        #: EYNİ pattern). IMPES-in öz daxili axın hesablaması (aşağıda,
        #: `self._trans`-a görə) BU FAZADA DƏYİŞMİR — yalnız GRID-in necə
        #: QURULDUĞU pluggable edilib.
        self.flux_discretization = flux_discretization or default_flux_discretization()
        _reject_multipoint_impes(self.flux_discretization)
        self._discretization = self.flux_discretization.build(model)
        self._connections = self._discretization.connections
        self._trans = self._discretization.transmissibility
        self._pv = self._discretization.pore_volume
        self._well_conn: List[WellConnection] = PeacemanWellModel().build_connections(model)

        self._setup_fluid_model()
        self._setup_gravity()
        self._sw_min, self._sw_max = relperm.saturation_limits()

        self._producer_names = sorted({c.well_name for c in self._well_conn
                                       if not c.is_injector})
        self._prealloc()
        self._apply_initial_conditions()

    # ------------------------------------------------------------ cazibə
    def _setup_gravity(self):
        """Üzlər üzrə dərinlik fərqi. Düz layda sıfırdır — cazibə üzvü itir
        və hesablama əvvəlki ilə eyni qalır (reqressiya testi bunu qoruyur)."""
        conn = self._connections
        depths = self.model.geometry.cell_depths()
        self._depth_difference = depths[conn.cell_a] - depths[conn.cell_b]
        self._has_gravity = bool(np.any(np.abs(self._depth_difference) > 1e-12))

    def _phase_densities(self, bw, bo):
        """Lay şəraitində sıxlıq, kg/m3."""
        fluids = self.model.fluids
        return (fluids.water_density / np.maximum(bw, 1e-9),
                fluids.oil_density / np.maximum(bo, 1e-9))

    def _face_potential_terms(self, pressure, sw, lam_w, lam_o, bw, bo):
        """(dPhi_o, dPhi_w) — üzlər üzrə faza potensialı fərqləri, bar.

            Φ_o = p            − ρ_o·g·D
            Φ_w = p − Pc(Sw)   − ρ_w·g·D
        """
        conn = self._connections
        dp = pressure[conn.cell_a] - pressure[conn.cell_b]
        d_phi_o = dp
        d_phi_w = dp

        if self._has_gravity:
            rho_w, rho_o = self._phase_densities(bw, bo)
            face_rho_w = 0.5 * (rho_w[conn.cell_a] + rho_w[conn.cell_b])
            face_rho_o = 0.5 * (rho_o[conn.cell_a] + rho_o[conn.cell_b])
            head = GRAVITY * self._depth_difference * PA_TO_BAR
            d_phi_o = d_phi_o - face_rho_o * head
            d_phi_w = d_phi_w - face_rho_w * head

        if self.capillary is not None:
            pc = self.capillary.pcow(sw)
            d_phi_w = d_phi_w - (pc[conn.cell_a] - pc[conn.cell_b])

        return d_phi_o, d_phi_w

    # ------------------------------------------------------- flüid modeli
    def _setup_fluid_model(self):
        """Statik və PVT yollarını bir interfeys altında birləşdirir.

        PVT provider verilmədikdə modelin sabit dəyərləri massivə çevrilir —
        hesablama tam eyni qalır (reqressiya testi bunu qoruyur).
        """
        model = self.model
        n = model.ncell
        fl = model.fluids
        ic = model.initial_conditions

        if self.pvt is None:
            ct = max(fl.water_compressibility * ic.water_saturation
                     + fl.oil_compressibility * (1.0 - ic.water_saturation)
                     + model.rock.compressibility, 1e-6)
            self._static = dict(
                mu_w=np.full(n, fl.water_viscosity),
                mu_o=np.full(n, fl.oil_viscosity),
                bw=np.full(n, fl.water_fvf),
                bo=np.full(n, fl.oil_fvf),
                ct=np.full(n, ct),
            )
            self._dfw_max = max(self.relperm.max_fractional_flow_derivative(
                fl.water_viscosity, fl.oil_viscosity), 1e-6)
        else:
            self._static = None
            self._dfw_max = self._worst_case_fractional_flow_derivative()

    def _worst_case_fractional_flow_derivative(self) -> float:
        """CFL limiti üçün ən pis (ən böyük) dfw/dSw.

        PVT halında lözlüklər təzyiqdən asılıdır, ona görə cədvəlin bütün
        təzyiq aralığı taranır və maksimum götürülür. Bu, yeni fizika deyil —
        mövcud CFL şərtinin dəyişkən lözlüyə uyğunlaşdırılmasıdır.
        """
        table = self.pvt.table if hasattr(self.pvt, "table") else None
        if table is None:
            ic = self.model.initial_conditions
            return max(self.relperm.max_fractional_flow_derivative(
                float(self.pvt.water_viscosity(ic.datum_pressure)),
                float(self.pvt.oil_viscosity(ic.datum_pressure))), 1e-6)
        worst = 0.0
        for pressure in table.pressure:
            worst = max(worst, self.relperm.max_fractional_flow_derivative(
                float(self.pvt.water_viscosity(pressure)),
                float(self.pvt.oil_viscosity(pressure))))
        return max(worst, 1e-6)

    def _fluid_state(self, pressure, sw):
        """(mu_w, mu_o, bw, bo, ct) — hüceyrə üzrə massivlər."""
        if self._static is not None:
            s = self._static
            return s["mu_w"], s["mu_o"], s["bw"], s["bo"], s["ct"]
        return (self.pvt.water_viscosity(pressure),
                self.pvt.oil_viscosity(pressure),
                self.pvt.water_fvf(pressure),
                self.pvt.oil_fvf(pressure),
                self.pvt.total_compressibility(pressure, sw))

    def _injection_mobility(self, mu_w):
        krw_end = self.relperm.endpoint_water_mobility(1.0)
        return krw_end / mu_w

    # ------------------------------------------------------------ setup
    def _apply_initial_conditions(self):
        n = self.model.ncell
        ic = self.model.initial_conditions
        if self.initialization is not None:
            state = self.initialization.initialize(self.model)
            self.pressure = np.asarray(state.pressure, float).copy()
            self.sw = np.asarray(state.water_saturation, float).copy()
        else:
            self.pressure = np.full(n, ic.datum_pressure, float)
            self.sw = np.full(n, ic.water_saturation, float)
        self.sw = np.clip(self.sw, self._sw_min, self._sw_max)

    def _prealloc(self):
        """Matris strukturu bir dəfə qurulur, hər addımda yalnız dəyərlər dəyişir.

        Əvvəl hər zaman addımında COO matris yaradılıb CSR-ə çevrilirdi:
        `coo_tocsr` + `csr_sort_indices` + `csr_sum_duplicates` ümumi vaxtın
        təxminən 12 %-ni yeyirdi. Struktur sabit olduğuna görə bir dəfə
        hesablanır və COO girişlərinin CSR `data` massivindəki mövqeyi
        (`_data_index`) yadda saxlanılır. Sonrakı addımlarda yalnız
        `bincount` ilə toplama qalır.
        """
        n = self.model.ncell
        conn = self._connections
        diag = np.arange(n)
        bhp_cells = np.array([c.cell for c in self._well_conn
                              if c.mode is ControlMode.BHP], dtype=int)
        self._nface = conn.count
        self._rows = np.concatenate([conn.cell_a, conn.cell_b, conn.cell_a,
                                     conn.cell_b, diag, bhp_cells])
        self._cols = np.concatenate([conn.cell_a, conn.cell_b, conn.cell_b,
                                     conn.cell_a, diag, bhp_cells])
        self._vals = np.zeros(self._rows.size)

        # COO -> CSR uyğunluğu: hər COO girişi hansı CSR mövqeyinə düşür
        entries = np.arange(self._rows.size, dtype=np.int64)
        pattern = sp.coo_matrix((entries + 1.0, (self._rows, self._cols)),
                                shape=(n, n)).tocsr()
        pattern.sum_duplicates()
        self._matrix = sp.csr_matrix(
            (np.zeros(pattern.nnz), pattern.indices.copy(),
             pattern.indptr.copy()), shape=(n, n))

        locator = sp.csr_matrix(
            (np.arange(pattern.nnz, dtype=np.int64),
             pattern.indices.copy(), pattern.indptr.copy()), shape=(n, n))
        self._data_index = np.asarray(
            locator[self._rows, self._cols]).ravel().astype(np.int64)
        self._nnz = pattern.nnz

    # -------------------------------------------------------- mobility
    def _mobilities(self, sw, mu_w, mu_o):
        lam_w = self.relperm.krw(sw) / mu_w
        lam_o = self.relperm.kro(sw) / mu_o
        return lam_w, lam_o

    # -------------------------------------------------------- pressure
    def _solve_pressure(self, dt: float):
        model = self.model
        conn = self._connections
        n = model.ncell
        mu_w, mu_o, bw, bo, ct = self._fluid_state(self.pressure, self.sw)
        lam_w, lam_o = self._mobilities(self.sw, mu_w, mu_o)
        lam_t = lam_w + lam_o
        inj_mobility = self._injection_mobility(mu_w)

        d_phi_o, d_phi_w = self._face_potential_terms(
            self.pressure, self.sw, lam_w, lam_o, bw, bo)
        up_o = np.where(d_phi_o >= 0, conn.cell_a, conn.cell_b)
        up_w = np.where(d_phi_w >= 0, conn.cell_a, conn.cell_b)
        t_o = self._trans * lam_o[up_o]
        t_w = self._trans * lam_w[up_w]
        tt = t_o + t_w

        # cazibə və kapilyar üzvləri implicit deyil — sağ tərəfə keçir
        dp_face = self.pressure[conn.cell_a] - self.pressure[conn.cell_b]
        gravity_capillary_flux = (t_o * (d_phi_o - dp_face)
                                  + t_w * (d_phi_w - dp_face))

        acc = self._pv * ct / dt
        nf = self._nface
        v = self._vals
        v[0:nf] = tt
        v[nf:2 * nf] = tt
        v[2 * nf:3 * nf] = -tt
        v[3 * nf:4 * nf] = -tt
        v[4 * nf:4 * nf + n] = acc
        rhs = acc * self.pressure
        if self._has_gravity or self.capillary is not None:
            np.add.at(rhs, conn.cell_a, -gravity_capillary_flux)
            np.add.at(rhs, conn.cell_b, +gravity_capillary_flux)

        k = 4 * nf + n
        for c in self._well_conn:
            lam = inj_mobility[c.cell] if c.is_injector else lam_t[c.cell]
            if c.mode is ControlMode.BHP:
                a = c.well_index * lam
                v[k] = a
                k += 1
                rhs[c.cell] += a * c.target
            else:
                rhs[c.cell] += abs(c.target) if c.is_injector else -abs(c.target)

        self._matrix.data[:] = np.bincount(self._data_index, weights=v,
                                           minlength=self._nnz)
        pressure = self.linear_solver.solve(self._matrix, rhs, x0=self.pressure)
        return pressure, lam_w, lam_o, lam_t, bw, bo, inj_mobility

    # ------------------------------------------------------ saturation
    def _update_saturation(self, pressure, lam_w, lam_o, lam_t, bw, bo,
                           inj_mobility, dt):
        conn = self._connections
        n = self.model.ncell
        d_phi_o, d_phi_w = self._face_potential_terms(
            pressure, self.sw, lam_w, lam_o, bw, bo)
        up_o = np.where(d_phi_o >= 0, conn.cell_a, conn.cell_b)
        up_w = np.where(d_phi_w >= 0, conn.cell_a, conn.cell_b)

        water_flux = self._trans * lam_w[up_w] * d_phi_w
        oil_flux = self._trans * lam_o[up_o] * d_phi_o
        total_flux = water_flux + oil_flux

        net_water = np.zeros(n)
        np.add.at(net_water, conn.cell_a, -water_flux)
        np.add.at(net_water, conn.cell_b, +water_flux)

        throughput = np.zeros(n)
        np.add.at(throughput, conn.cell_a, np.abs(total_flux))
        np.add.at(throughput, conn.cell_b, np.abs(total_flux))

        qo_total = qw_total = qwi_total = 0.0
        well_oil = {name: 0.0 for name in self._producer_names}
        well_water = {name: 0.0 for name in self._producer_names}

        for c in self._well_conn:
            cell = c.cell
            if c.is_injector:
                if c.mode is ControlMode.BHP:
                    q = c.well_index * inj_mobility[cell] * (c.target - pressure[cell])
                else:
                    q = abs(c.target)
                q = max(q, 0.0)
                net_water[cell] += q
                throughput[cell] += q
                qwi_total += q
            else:
                if c.mode is ControlMode.BHP:
                    dpw = c.target - pressure[cell]
                    qw = c.well_index * lam_w[cell] * dpw
                    qo = c.well_index * (lam_t[cell] - lam_w[cell]) * dpw
                else:
                    qt = -abs(c.target)
                    frac = lam_w[cell] / max(lam_t[cell], 1e-30)
                    qw, qo = qt * frac, qt * (1.0 - frac)
                qw, qo = min(qw, 0.0), min(qo, 0.0)
                net_water[cell] += qw
                throughput[cell] += abs(qw) + abs(qo)
                qw_total += -qw / bw[cell]
                qo_total += -qo / bo[cell]
                well_oil[c.well_name] += -qo / bo[cell]
                well_water[c.well_name] += -qw / bw[cell]

        sw_new = self.sw + dt * net_water / np.maximum(self._pv, 1e-12)
        sw_new = np.clip(sw_new, self._sw_min, self._sw_max)

        with np.errstate(divide="ignore", invalid="ignore"):
            dt_cfl = np.nanmin(self._pv / np.maximum(throughput * self._dfw_max, 1e-12))

        return (sw_new, dt_cfl, qo_total, qw_total,
                qwi_total, well_oil, well_water)

    # ------------------------------------------------------------ ooip
    def original_oil_in_place(self) -> float:
        _, _, _, bo, _ = self._fluid_state(self.pressure, self.sw)
        return float(np.sum(self._pv * (1.0 - self.sw) / bo))

    # ------------------------------------------------------------- run
    def run(self, reporter: Optional[IProgressReporter] = None) -> SimulationResult:
        reporter = reporter or NullProgressReporter()
        cfg = self.config
        ts, out = cfg.time_stepping, cfg.output

        result = SimulationResult(model_name=self.model.name,
                                  grid_shape=self.model.grid.shape)
        result.ooip = self.original_oil_in_place()
        result.well_oil_rate = {n: [] for n in self._producer_names}
        result.well_water_rate = {n: [] for n in self._producer_names}
        s = result.series

        t = 0.0
        dt = ts.initial_dt
        cum_o = cum_w = 0.0
        snap_every = max(cfg.end_time / max(out.snapshot_count, 1), 1e-9)
        next_snap = 0.0
        step = 0

        self._record_snapshot(result, t)
        next_snap += snap_every

        while t < cfg.end_time - 1e-9 and step < ts.max_steps:
            dt = min(dt, ts.max_dt, cfg.end_time - t)
            try:
                (pressure, lam_w, lam_o, lam_t, bw, bo,
                 inj_mob) = self._solve_pressure(dt)
                (sw_new, dt_cfl, qo, qw, qwi,
                 well_oil, well_water) = self._update_saturation(
                    pressure, lam_w, lam_o, lam_t, bw, bo, inj_mob, dt)
            except (FloatingPointError, RuntimeError) as exc:
                result.converged = False
                result.message = f"Addım {step}, t={t:.1f} gün: {exc}"
                LOG.error("Divergensiya: %s", result.message)
                break

            if dt > ts.cfl_factor * dt_cfl and dt > ts.min_dt:
                dt = max(ts.cfl_factor * dt_cfl, ts.min_dt)
                continue

            self.pressure, self.sw = pressure, sw_new
            t += dt
            step += 1
            cum_o += qo * dt
            cum_w += qw * dt

            s.time.append(t)
            s.oil_rate.append(qo)
            s.water_rate.append(qw)
            s.water_injection_rate.append(qwi)
            s.cumulative_oil.append(cum_o)
            s.cumulative_water.append(cum_w)
            s.water_cut.append(qw / max(qo + qw, 1e-12) * 100.0)
            s.average_pressure.append(float(np.mean(self.pressure)))
            s.recovery_factor.append(cum_o / max(result.ooip, 1e-12) * 100.0)
            if out.record_well_rates:
                for name in self._producer_names:
                    result.well_oil_rate[name].append(well_oil[name])
                    result.well_water_rate[name].append(well_water[name])

            if t >= next_snap - 1e-9:
                self._record_snapshot(result, t)
                next_snap += snap_every

            dt = min(dt * ts.growth_factor, ts.cfl_factor * dt_cfl, ts.max_dt)

            if step % out.progress_every_n_steps == 0:
                msg = (f"t = {t:8.1f} gün | RF = {s.recovery_factor[-1]:5.2f} % "
                       f"| dt = {dt:.3f}")
                if not reporter.report(t / cfg.end_time * 100.0, msg):
                    result.message = "İstifadəçi tərəfindən dayandırıldı."
                    break

        if result.snapshots and result.snapshots[-1].time < t - 1e-9:
            self._record_snapshot(result, t)
        result.steps = step
        if not result.message:
            result.message = f"Tamamlandı: {step} addım, t = {t:.1f} gün."
        LOG.info("%s  RF = %.2f %%  OOIP = %.0f m3", result.message,
                 result.final_recovery_factor, result.ooip)
        return result

    def _record_snapshot(self, result: SimulationResult, t: float):
        shape = self.model.grid.shape
        result.snapshots.append(Snapshot(
            time=t,
            pressure=self.pressure.reshape(shape).copy(),
            water_saturation=self.sw.reshape(shape).copy(),
        ))
