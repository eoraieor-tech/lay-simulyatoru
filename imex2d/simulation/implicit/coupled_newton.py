"""Birləşmiş Nyuton həlledicisi — OPM tipli quyu modeli, MƏRHƏLƏ 4.

v69 addım 4b: modul İKİ FAZALI (neft-su) oldu. Rezervuar hissəsi artıq
əsas mühərrikin `ResidualAssembler` / `JacobianAssembler` sinifləri ilə
qurulur (quyusuz, `wells=[]` ilə), quyu töhfələri isə
`StandardWellModel` / `StandardWellJacobian`-dan gəlir. Beləliklə
rezervuar fizikası TƏK yerdə qalır — kod təkrarlanmır və A6 ilə
avtomatik uyğun olur.

Mərhələ 1-3-də hazırlanan hissələr burada BİRLƏŞİR:

    mərhələ 1  →  `CoupledState` (rezervuar + quyu naməlumları)
    mərhələ 2  →  `StandardWellModel` (debitlər, idarəetmə tənlikləri)
    mərhələ 3  →  `StandardWellJacobian` (dörd blok)

Nəticədə tam sistem (2N + W) ölçüsündə həll olunur:

    ┌─────────────┬──────────┐  ┌────┐   ┌────┐
    │  rezervuar  │  R↔Q     │  │ δx │   │ R  │
    ├─────────────┼──────────┤  ├────┤ = ├────┤
    │  Q↔R        │   Q      │  │δbhp│   │ Rc │
    └─────────────┴──────────┘  └────┘   └────┘

MÖVCUD MÜHƏRRİKƏ TƏSİR ETMİR

Bu sinif `NewtonSolver`-in YANINDA yaşayır. Mühərrik seçimi
mərhələ 5-də (və ya doğrulama uğurlu olandan sonra) ediləcək;
ondan əvvəl işləyən kod tam toxunulmaz qalır.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from ...logging_setup import get_logger
from .jacobian import JacobianAssembler
from .newton import NewtonConfig, NewtonStatus
from .residual import FluidState, OIL, ResidualAssembler, WATER
from .standard_well import StandardWellJacobian, StandardWellModel
from .state import VARIABLES_PER_CELL, ReservoirState
from .well_state import CoupledState, WellUnknowns

LOG = get_logger(__name__)


@dataclass
class CoupledNewtonResult:
    status: NewtonStatus
    state: CoupledState
    iterations: int
    history: List[float] = field(default_factory=list)
    material_balance: tuple = (0.0, 0.0)
    fluid: Optional[FluidState] = None
    rates: Optional[object] = None

    @property
    def converged(self) -> bool:
        return self.status is NewtonStatus.CONVERGED


class CoupledNewtonSolver:
    """Rezervuar + quyu naməlumlarını BİRLİKDƏ həll edir."""

    def __init__(self, model, relperm, pvt, linear_solver, grid, connections,
                 config: Optional[NewtonConfig] = None,
                 endpoint_water_mobility: float = 0.35,
                 capillary=None):
        self.model = model
        self.relperm = relperm
        self.pvt = pvt
        self.linear_solver = linear_solver
        self.grid = grid

        # Rezervuar hissəsi QUYUSUZ qurulur — quyu töhfələri
        # `StandardWell*`-dan gəlir (kod təkrarlanmır).
        self.reservoir = ResidualAssembler(model, grid, [], relperm, pvt,
                                           capillary)
        self.reservoir_jacobian = JacobianAssembler(self.reservoir)

        self.well_model = StandardWellModel(connections, model.ncell,
                                            endpoint_water_mobility)
        self.well_jacobian = StandardWellJacobian(
            self.well_model, pvt, relperm,
            derivatives=self.reservoir_jacobian.derivatives)

        # ── quyu tənliklərinin miqyaslanması ────────────────────────
        # Ölçülüb (A7_PLAN.md): rezervuarın akkumulyasiya diaqonalı
        # (PV/dt) kiçik dt-də (məs. 0.05) 25 000-ə çatır, quyu
        # tənliyinin diaqonalı isə HƏMİŞƏ 1 qalır (`R=p_bhp−hədəf`,
        # bar). Sabit çarpanla düzəlmir — dt-dən ASILI miqyaslama
        # lazımdır.
        self._reference_pore_volume = float(
            np.median(grid.pore_volume[grid.pore_volume > 0]))

        self.config = config or NewtonConfig(max_iterations=40)
        # Birləşmiş sistem (rezervuar+quyu) A6-nın xalis rezervuar
        # sisteminə görə bir qədər YAVAŞ yığılır (əlavə W tənlik,
        # əlavə qoşulma) — 25 bəzən kifayət etmir (CNV yavaş-yavaş
        # azalır, divergensiya YOX). 40-a qaldırıldı.
        sw_min, sw_max = relperm.saturation_limits()
        relaxation = self.config.saturation_bound_relaxation
        self.sw_min = sw_min - relaxation
        self.sw_max = sw_max + relaxation
        self.physical_sw_min = sw_min
        self.physical_sw_max = sw_max

    # ── flüid xassələri ────────────────────────────────────────────
    def build_fluid(self, reservoir: ReservoirState) -> FluidState:
        return self.reservoir.fluid_state(reservoir)

    # ── qalıq ──────────────────────────────────────────────────────
    def compute_residual(self, state: CoupledState, previous: CoupledState,
                         previous_fluid: FluidState, dt: float):
        """Birləşmiş qalıq: rezervuar (2N) + idarəetmə (W)."""
        reservoir_part, fluid, _ = self.reservoir.residual(
            state.reservoir, previous.reservoir, dt, previous_fluid)

        rates = self.well_model.perforation_rates(state.reservoir, fluid,
                                                  state.wells)
        # `R = … − q` (quyu töhfəsi mənfi işarə ilə daxil olur)
        reservoir_part[WATER::VARIABLES_PER_CELL] -= rates.water
        reservoir_part[OIL::VARIABLES_PER_CELL] -= rates.oil

        residual = np.empty(state.size)
        residual[:state.well_offset] = reservoir_part
        # bax `__init__`-in şərhi: dt-yə görə miqyaslanır ki,
        # Jakobianla (bax `assemble_jacobian`) uyğun qalsın
        control_scale = self._reference_pore_volume / dt
        residual[state.well_offset:] = (
            self.well_model.control_residuals(rates, state.wells)
            * control_scale)
        return residual, fluid, rates

    # ── Jakobian ───────────────────────────────────────────────────
    def assemble_jacobian(self, state: CoupledState, fluid: FluidState,
                          dt: float):
        import scipy.sparse as sp

        reservoir = state.reservoir
        size = state.size
        offset = state.well_offset

        # rezervuar bloku (quyusuz) — mövcud yığıcı
        reservoir_matrix = self.reservoir_jacobian.assemble(reservoir, fluid,
                                                            dt)

        blocks = self.well_jacobian.blocks(reservoir, fluid, state.wells)

        rows: List[np.ndarray] = []
        cols: List[np.ndarray] = []
        values: List[np.ndarray] = []

        # quyu debitlərinin rezervuar dəyişənlərinə görə törəməsi.
        # `R = … − q`, ona görə İŞARƏ MƏNFİdir.
        cell_index = np.arange(reservoir.ncell)
        well_diagonal = -blocks.rate_wrt_reservoir
        for r in range(VARIABLES_PER_CELL):
            for c in range(VARIABLES_PER_CELL):
                entries = well_diagonal[:, r, c]
                nonzero = entries != 0.0
                if not np.any(nonzero):
                    continue
                rows.append(cell_index[nonzero] * VARIABLES_PER_CELL + r)
                cols.append(cell_index[nonzero] * VARIABLES_PER_CELL + c)
                values.append(entries[nonzero])

        # R↔Q: rezervuar tənliklərinin BHP-yə görə törəməsi
        for cell, column in blocks.rate_wrt_bhp.items():
            well_position = blocks.rate_bhp_owner[cell]
            for phase in range(VARIABLES_PER_CELL):
                if column[phase] == 0.0:
                    continue
                rows.append(np.array([cell * VARIABLES_PER_CELL + phase]))
                cols.append(np.array([offset + well_position]))
                values.append(np.array([-column[phase]]))

        # Q↔R VƏ Q: idarəetmə tənliyinin BÜTÖV sətri eyni faktorla
        # miqyaslanır ki, rezervuarla UYĞUN qalsın (bax `__init__`-in
        # şərhi) — bütöv sətir eyni tənliyin hissəsidir, ona görə
        # HƏR İKİ blok eyni `control_scale` ilə vurulur.
        control_scale = self._reference_pore_volume / dt
        for well_position, coupling in blocks.control_wrt_reservoir.items():
            for cell, vector in coupling.items():
                for variable in range(VARIABLES_PER_CELL):
                    if vector[variable] == 0.0:
                        continue
                    rows.append(np.array([offset + well_position]))
                    cols.append(np.array([cell * VARIABLES_PER_CELL + variable]))
                    values.append(np.array([vector[variable] * control_scale]))

        # Q: idarəetmə tənliklərinin BHP-yə görə törəməsi (miqyaslı)
        for well_position, value in enumerate(blocks.control_wrt_bhp):
            if value == 0.0:
                continue
            rows.append(np.array([offset + well_position]))
            cols.append(np.array([offset + well_position]))
            values.append(np.array([value * control_scale]))

        extra = sp.coo_matrix(
            (np.concatenate(values) if values else np.array([]),
             (np.concatenate(rows) if rows else np.array([], dtype=int),
              np.concatenate(cols) if cols else np.array([], dtype=int))),
            shape=(size, size)).tocsr()

        reservoir_coo = reservoir_matrix.tocoo()
        expanded = sp.coo_matrix(
            (reservoir_coo.data, (reservoir_coo.row, reservoir_coo.col)),
            shape=(size, size)).tocsr()
        return (expanded + extra).tocsr()

    # ── yığılma meyarları ──────────────────────────────────────────
    def convergence_measures(self, residual: np.ndarray, state: CoupledState,
                             fluid: FluidState, dt: float):
        reservoir = state.reservoir
        pore_volume = self.reservoir.pore_volume_at(reservoir.pressure)
        scale = np.maximum(pore_volume / dt, 1e-30)

        reservoir_part = residual[:state.well_offset]
        water = reservoir_part[WATER::VARIABLES_PER_CELL]
        oil = reservoir_part[OIL::VARIABLES_PER_CELL]
        cnv = float(max((np.abs(water) / scale).max(),
                        (np.abs(oil) / scale).max()))

        # idarəetmə qalığı da yığılma meyarına daxildir — quyu tənliyi
        # ödənilməsə həll natamamdır. `residual`-in quyu hissəsi artıq
        # `control_scale`-lə vurulub (bax `compute_residual`) —
        # yığılma yoxlaması üçün ORİJİNAL (bar) miqyasına GERİ
        # qaytarırıq ki, `control_tolerance` meyarı mənalı qalsın.
        control_scale = self._reference_pore_volume / dt
        control = (float(np.abs(residual[state.well_offset:]).max())
                   / control_scale) if state.wells.count else 0.0

        in_place = self.reservoir.accumulation(reservoir, fluid)
        total_pore_volume = float(pore_volume.sum())
        balances = []
        for phase_residual, phase_in_place in zip((water, oil), in_place):
            reference = phase_in_place.sum()
            if reference < 1e-6 * max(total_pore_volume, 1e-30):
                balances.append(0.0)
            else:
                balances.append(float(abs(phase_residual.sum()) * dt
                                      / reference))
        return cnv, control, tuple(balances)

    def _is_converged(self, cnv, control, balances) -> bool:
        tolerance = self.config.material_balance_tolerance
        return (cnv < self.config.cnv_tolerance
                and control < self.config.control_tolerance      # bar
                and all(value < tolerance for value in balances))

    # ── döngə ──────────────────────────────────────────────────────
    def solve(self, previous: CoupledState, dt: float) -> CoupledNewtonResult:
        """Bu metod HEÇ VAXT istisna atmır (mövcud həlledici ilə eyni
        zəmanət — bax `NewtonSolver.solve`)."""
        try:
            return self._solve_inner(previous, dt)
        except Exception as error:
            LOG.exception("Birləşmiş Nyuton gözlənilməz istisna (%s)",
                          type(error).__name__)
            return CoupledNewtonResult(NewtonStatus.LINEAR_SOLVER_FAILED,
                                       previous.copy(), 0, [], (0.0, 0.0))

    def _solve_inner(self, previous: CoupledState, dt: float
                     ) -> CoupledNewtonResult:
        config = self.config
        state = previous.copy()
        previous_fluid = self.build_fluid(previous.reservoir)
        # Quyunun açıq/bağlı vəziyyəti ADDIMIN ƏVVƏLİNDƏ bir dəfə
        # qərara alınır və bütün iterasiyalar boyu SABİT qalır —
        # bax `StandardWellModel.update_shut_wells()`.
        self.well_model.update_shut_wells(previous.reservoir, previous_fluid,
                                          previous.wells)
        history: List[float] = []
        balances = (0.0, 0.0)
        fluid = rates = None

        for iteration in range(config.max_iterations + 1):
            residual, fluid, rates = self.compute_residual(
                state, previous, previous_fluid, dt)
            cnv, control, balances = self.convergence_measures(
                residual, state, fluid, dt)
            history.append(cnv)

            if self._is_converged(cnv, control, balances):
                final = CoupledState(self._clamp(state.reservoir),
                                     state.wells.copy())
                return CoupledNewtonResult(NewtonStatus.CONVERGED, final,
                                           iteration, history, balances,
                                           fluid, rates)

            if iteration == config.max_iterations:
                break
            if not np.all(np.isfinite(residual)):
                return CoupledNewtonResult(NewtonStatus.LINEAR_SOLVER_FAILED,
                                           state, iteration, history,
                                           balances, fluid, rates)

            try:
                matrix = self.assemble_jacobian(state, fluid, dt)
                delta = self.linear_solver.solve(matrix, -residual)
            except Exception as error:
                LOG.debug("Xətti həll uğursuz (%s)", type(error).__name__)
                return CoupledNewtonResult(NewtonStatus.LINEAR_SOLVER_FAILED,
                                           state, iteration, history,
                                           balances, fluid, rates)
            if not np.all(np.isfinite(delta)):
                return CoupledNewtonResult(NewtonStatus.LINEAR_SOLVER_FAILED,
                                           state, iteration, history,
                                           balances, fluid, rates)

            pressure_limit, saturation_limit = config.limits_at(iteration)
            trial = state.updated(
                delta, self.sw_min, self.sw_max,
                max_pressure_change=pressure_limit,
                max_saturation_change=saturation_limit,
                max_bhp_change=config.max_bhp_change)

            # ── GERİ-İZLƏMƏ (line search) ───────────────────────────
            #
            # Kəsmə (Appleyard) TƏK BAŞINA kifayət etmir. Ölçüldü
            # (dt=0.25): Nyuton δp = 1456 bar istəyir, kəsmə onu 50
            # bar-a endirir — nəticədə İSTİQAMƏT pozulur və həll
            # sıçrayır, CNV isə 0.32-0.35-də ilişib qalır.
            #
            # Kəsmə yalnız addımın UZUNLUĞUNU məhdudlaşdırır, lakin
            # onun FAYDALI olub-olmadığını yoxlamır. Geri-izləmə məhz
            # bunu edir: qalığı azaltmayan addım qəbul edilmir, addım
            # yarıya bölünüb yenidən sınanılır.
            trial_residual, _, _ = self.compute_residual(
                trial, previous, previous_fluid, dt)
            current_norm = self._scaled_norm(residual, state, dt)
            trial_norm = self._scaled_norm(trial_residual, trial, dt)

            scale = 1.0
            for _ in range(10):
                if trial_norm <= current_norm * 0.999 or scale < 1.0 / 512.0:
                    break
                scale *= 0.5
                trial = state.updated(
                    delta * scale, self.sw_min, self.sw_max,
                    max_pressure_change=pressure_limit,
                    max_saturation_change=saturation_limit,
                    max_bhp_change=config.max_bhp_change)
                trial_residual, _, _ = self.compute_residual(
                    trial, previous, previous_fluid, dt)
                trial_norm = self._scaled_norm(trial_residual, trial, dt)
            state = trial

        return CoupledNewtonResult(NewtonStatus.MAX_ITERATIONS, state,
                                   config.max_iterations, history, balances,
                                   fluid, rates)

    def _scaled_norm(self, residual: np.ndarray, state: CoupledState,
                     dt: float) -> float:
        """Faza-miqyaslı qalıq normu — geri-izləmənin qəbul meyarı.

        Rezervuar hissəsi məsamə həcminə görə normallaşdırılır (CNV
        ilə eyni prinsip), quyu hissəsi isə onsuz da bar
        vahidindədir və olduğu kimi götürülür.
        """
        offset = state.well_offset
        pore_volume = self.reservoir.pore_volume_at(state.reservoir.pressure)
        scale = np.maximum(pore_volume / dt, 1e-30)
        water = residual[WATER:offset:VARIABLES_PER_CELL] / scale
        oil = residual[OIL:offset:VARIABLES_PER_CELL] / scale
        control = residual[offset:]
        return float(np.sqrt(np.mean(water ** 2 + oil ** 2)
                             + np.sum(control ** 2)))

    def _clamp(self, reservoir: ReservoirState) -> ReservoirState:
        water = np.clip(reservoir.water_saturation, self.physical_sw_min,
                        self.physical_sw_max)
        return ReservoirState(reservoir.pressure, water)

    # ── MƏRHƏLƏ 6: zaman addımının quyu debitinə görə hədd ─────────
    def max_stable_dt(self, state: CoupledState, alpha: float = 0.2) -> float:
        """Quyu debitinə görə "CFL-bənzər" zaman addımı həddi, gün.

        KÖK SƏBƏB (bax `A7_PLAN.md`): bir hüceyrənin məsamə həcmi
        (PV) 1250 m³, istismarçı isə 654 m³/gün verə bilir — dt=1.0-da
        bu, hüceyrənin YARISIdır BİR ADDIMDA. Heç bir implicit sxem
        bunu hamar keçirə bilməz, quyu modelindən ASILI OLMAYARAQ.

        Şərt (klassik CFL prinsipi, quyu axınına tətbiq olunub):

            dt · Σ|q_quyu| ≤ α · PV_hüceyrə     hər perforasiya üçün
            dt ≤ α · PV_hüceyrə / Σ|q_quyu|

        `alpha` (defolt 0.2): bir addımda hüceyrənin ən çox 20 %-i
        çıxarılsın/vurulsun. OPM-də bənzər fikir "quyu CFL ədədi"
        adı ilə tanınır.

        Debit ƏVVƏLKİ (converged) vəziyyətdən hesablanır — bu, dt-ni
        seçərkən hələ YENİ vəziyyət məlum olmadığı üçün TƏBİİDİR (A6
        özü də eyni prinsiplə addım BÖYÜMƏSİNİ keçmiş nəticəyə görə
        seçir).
        """
        fluid = self.build_fluid(state.reservoir)
        rates = self.well_model.perforation_rates(state.reservoir, fluid,
                                                  state.wells)
        pore_volume = self.reservoir.pore_volume_at(state.reservoir.pressure)

        # səth həcmi debitlərini TƏXMİNİ rezervuar həcminə gətiririk
        # (FVF ilə) ki, müqayisə mənalı olsun
        withdrawal = (np.abs(rates.water) * fluid.bw
                      + np.abs(rates.oil) * fluid.bo)
        active = withdrawal > 1e-9
        if not np.any(active):
            return float("inf")           # heç bir quyu axını yoxdur

        limits = alpha * pore_volume[active] / withdrawal[active]
        return float(limits.min())
