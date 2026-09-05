"""Nyuton-Rafson döngəsi — A6, mərhələ 3.

    x⁰ = xⁿ
    təkrarla:
        J(x^k) · δ = −R(x^k)
        x^{k+1} = x^k + chop(δ)
        konvergensiya yoxlanılır

İKİ KONVERGENSİYA MEYARI (sənaye standartı, Eclipse/CMG kimi):

  CNV — lokal qalıq. Ən pis hüceyrədə tənlik nə qədər pozulub:
        CNV_p = max_c |R_p,c| · Δt / PV_c
        Bu, hər hüceyrənin ayrıca həll olunduğunu təmin edir.

  MB  — qlobal kütlə balansı. Bütün model üzrə cəmdə nə qədər
        flüid itir və ya yaranır:
        MB_p = |Σ_c R_p,c| · Δt / Σ_c (PV_c·S_p/B_p)
        Bu, CNV-dən qat-qat sıx olmalıdır — lokal səhvlər bir-birini
        yeyə bilər, amma ümumi kütlə itkisi yolverilməzdir.

Hər iki meyar ÖLÇÜSÜZDÜR, ona görə eyni tolerans bütün grid
ölçülərində və zaman addımlarında işləyir.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional

import numpy as np

from ...logging_setup import get_logger
from .active_reduction import ActiveDofReduction
from .jacobian import JacobianAssembler
from .linear import NewtonLinearSolver
from .residual import OIL, WATER, ResidualAssembler
from .state import ReservoirState

LOG = get_logger(__name__)


class NewtonStatus(Enum):
    CONVERGED = "yığıldı"
    MAX_ITERATIONS = "iterasiya limiti"
    DIVERGED = "divergensiya"
    LINEAR_SOLVER_FAILED = "xətti həlledici uğursuz"


@dataclass
class NewtonConfig:
    max_iterations: int = 20
    """Geri-izləmə əlavə olunandan sonra ölçüldü: çətin addımlar CNV-ni
    MONOTON (dövr etmədən) azaldır, lakin tolerantlığa çatmaq üçün
    12-dən çox — tipik 14-18 — iterasiya tələb edə bilir. 20 bu
    marjını verir, xərci isə əhəmiyyətsizdir (asan addımlar 1-4
    iterasiyada, bu hədəd ora heç toxunmur)."""
    cnv_tolerance: float = 1e-3
    material_balance_tolerance: float = 1e-7
    max_pressure_change: float = 50.0        # bar, bir iterasiyada
    max_saturation_change: float = 0.2       # bir iterasiyada
    relaxation_start: int = 0
    """Chopping-in boşaldılmağa başladığı iterasiya (0 = heç vaxt).

    Nəzəri olaraq chopping-i boşaltmaq Nyutonun kvadratik yığılmasına
    imkan verməlidir. Ölçmə isə əksini göstərdi: sərt hədlər sərt
    məsələlərdə daha yaxşı işləyir, çünki addımın kəsilməsi həllin
    fiziki intervalda qalmasını təmin edir. Ona görə defolt = söndürülü;
    parametr xüsusi hallar üçün saxlanılıb.
    """
    relaxation_factor: float = 2.0
    divergence_factor: float = 10.0
    """Qalıq başlanğıcdan bu qədər dəfə böyüsə, divergensiya sayılır."""
    saturation_bound_relaxation: float = 0.01

    # ── OPM tipli quyu modeli üçün (bax `coupled_newton.py`) ─────────
    control_tolerance: float = 1.0e-4
    """Quyu idarəetmə qalığının həddi, BAR.

    BHP idarəsində qalıq `p_bhp − p_hədəf`-dir, yəni birbaşa bar
    vahidindədir — 10⁻⁴ bar praktiki olaraq dəqiq deməkdir.
    """

    max_bhp_change: float = 50.0
    """Quyu BHP-sinin bir iterasiyada maksimum dəyişməsi, bar.

    Rezervuar təzyiqinin Appleyard kəsməsi ilə eyni məqsəd: quyu
    təzyiqi bir iterasiyada həddindən çox sıçrasa, perforasiya
    debitləri qeyri-real dəyərlərə gedər.
    """
    """Doyumluluq hədlərinin ədədi zolağı.

    Sərt hədd (Sw ≥ Swc) məsələni bound-constrained edir: hədə ilişən
    hüceyrələrdə qalıq sıfır ola bilmir və Nyoton dayanır. Ölçmə
    göstərdi ki, belə hallarda qalığın 100 %-i məhz ilişən hüceyrələrdən
    gəlir, sərbəstlərdə isə dəqiq sıfırdır.

    Kiçik zolaq bu problemi aradan qaldırır və fiziki nəticəyə təsir
    etmir: nisbi keçiricilik onsuz da Corey düsturunda kəsilir, yəni
    Sw = Swc − 0.001 halında krw hələ də sıfırdır.
    """

    def limits_at(self, iteration: int) -> tuple:
        """(Δp_maks, ΔSw_maks) — iterasiya nömrəsinə görə."""
        if self.relaxation_start <= 0 or iteration < self.relaxation_start:
            return self.max_pressure_change, self.max_saturation_change
        factor = self.relaxation_factor ** (iteration - self.relaxation_start + 1)
        return (self.max_pressure_change * factor,
                min(self.max_saturation_change * factor, 1.0))


@dataclass
class NewtonResult:
    status: NewtonStatus
    state: ReservoirState
    iterations: int = 0
    cnv_history: List[float] = field(default_factory=list)
    material_balance: tuple = (0.0, 0.0)
    fluid: object = None
    rates: object = None
    linear_iterations: int = 0

    @property
    def converged(self) -> bool:
        return self.status is NewtonStatus.CONVERGED


class NewtonSolver:
    """Bir zaman addımını implicit həll edir."""

    def __init__(self, residual: ResidualAssembler,
                 jacobian: Optional[JacobianAssembler] = None,
                 linear_solver: Optional[NewtonLinearSolver] = None,
                 config: Optional[NewtonConfig] = None):
        self.R = residual
        self.J = jacobian or JacobianAssembler(residual)
        self.linear_solver = linear_solver or NewtonLinearSolver()
        self.config = config or NewtonConfig()
        #: ACTNUM (tapşırıq §1): xətti sistem `2·n_active` naməlumla həll
        #: olunur — bax `active_reduction.py`. Bütün hüceyrələr aktivdirsə
        #: bu obyekt matrisə TOXUNMUR.
        self.reduction = ActiveDofReduction(residual.model.grid.active)
        sw_min, sw_max = residual.relperm.saturation_limits()
        relaxation = self.config.saturation_bound_relaxation
        self.sw_min = sw_min - relaxation
        self.sw_max = sw_max + relaxation
        self.physical_sw_min = sw_min
        self.physical_sw_max = sw_max

    # ═══════════════════════════════════════════════ konvergensiya
    def convergence_measures(self, residual_vector: np.ndarray,
                             state: ReservoirState, fluid, dt: float):
        """(CNV, MB_su, MB_neft) — üçü də ölçüsüzdür."""
        pore_volume = self.R.pore_volume
        scale = np.maximum(pore_volume / dt, 1e-30)

        water = residual_vector[WATER::2]
        oil = residual_vector[OIL::2]
        cnv = float(max((np.abs(water) / scale).max(),
                        (np.abs(oil) / scale).max()))

        water_in_place, oil_in_place = self.R.accumulation(state, fluid)
        mb_water = float(abs(water.sum()) * dt
                         / max(water_in_place.sum(), 1e-30))
        mb_oil = float(abs(oil.sum()) * dt / max(oil_in_place.sum(), 1e-30))
        return cnv, mb_water, mb_oil

    def _is_converged(self, cnv: float, mb_water: float, mb_oil: float) -> bool:
        config = self.config
        return (cnv < config.cnv_tolerance
                and mb_water < config.material_balance_tolerance
                and mb_oil < config.material_balance_tolerance)

    def _clamp_to_physical(self, state: ReservoirState) -> ReservoirState:
        """Ədədi zolaq yalnız həll müddətindədir — nəticə fiziki hədlərdədir."""
        return ReservoirState(
            state.pressure,
            np.clip(state.water_saturation,
                    self.physical_sw_min, self.physical_sw_max))

    # ═════════════════════════════════════════════════════ döngə
    def solve(self, previous: ReservoirState, dt: float,
              initial_guess: Optional[ReservoirState] = None) -> NewtonResult:
        config = self.config
        state = (initial_guess or previous).copy()
        previous_fluid = self.R.fluid_state(previous)

        history: List[float] = []
        first_cnv = None
        linear_iterations = 0
        fluid = rates = None
        material_balance = (0.0, 0.0)

        for iteration in range(config.max_iterations + 1):
            residual_vector, fluid, rates = self.R.residual(
                state, previous, dt, previous_fluid)
            cnv, mb_water, mb_oil = self.convergence_measures(
                residual_vector, state, fluid, dt)
            history.append(cnv)
            material_balance = (mb_water, mb_oil)

            if self._is_converged(cnv, mb_water, mb_oil):
                return NewtonResult(NewtonStatus.CONVERGED,
                                    self._clamp_to_physical(state), iteration,
                                    history, material_balance, fluid, rates,
                                    linear_iterations)

            if first_cnv is None:
                first_cnv = max(cnv, 1e-30)
            elif cnv > first_cnv * config.divergence_factor:
                LOG.debug("Nyuton divergensiya etdi: CNV %.3e -> %.3e",
                          first_cnv, cnv)
                return NewtonResult(NewtonStatus.DIVERGED, state, iteration,
                                    history, material_balance, fluid, rates,
                                    linear_iterations)

            if iteration == config.max_iterations:
                break

            try:
                matrix = self.J.assemble(state, fluid, dt)
                reduced_matrix, reduced_rhs = self.reduction.restrict(
                    matrix, -residual_vector)
                delta = self.reduction.expand(
                    self.linear_solver.solve(reduced_matrix, reduced_rhs))
                linear_iterations += self.linear_solver.last_iterations
            except (FloatingPointError, ValueError) as error:
                LOG.debug("Xətti həlledici uğursuz: %s", error)
                return NewtonResult(NewtonStatus.LINEAR_SOLVER_FAILED, state,
                                    iteration, history, material_balance,
                                    fluid, rates, linear_iterations)

            if not np.all(np.isfinite(delta)):
                return NewtonResult(NewtonStatus.LINEAR_SOLVER_FAILED, state,
                                    iteration, history, material_balance,
                                    fluid, rates, linear_iterations)

            pressure_limit, saturation_limit = config.limits_at(iteration)
            trial = state.updated(
                delta, self.sw_min, self.sw_max,
                max_pressure_change=pressure_limit,
                max_saturation_change=saturation_limit)
            state = self._line_search(state, trial, delta, previous, dt,
                                      previous_fluid, residual_vector,
                                      pressure_limit, saturation_limit)

        return NewtonResult(NewtonStatus.MAX_ITERATIONS, state,
                            config.max_iterations, history, material_balance,
                            fluid, rates, linear_iterations)

    # ═══════════════════════════════════════════════════ geri-izləmə
    def _scaled_norm(self, residual_vector: np.ndarray,
                     state: ReservoirState, dt: float) -> float:
        """Geri-izləmənin qəbul meyarı — CNV ilə eyni miqyaslı RMS qalıq."""
        pore_volume = self.R.pore_volume_at(state.pressure)
        scale = np.maximum(pore_volume / dt, 1e-30)
        water = residual_vector[WATER::2] / scale
        oil = residual_vector[OIL::2] / scale
        return float(np.sqrt(np.mean(water ** 2 + oil ** 2)))

    def _line_search(self, state: ReservoirState, trial: ReservoirState,
                     delta: np.ndarray, previous: ReservoirState, dt: float,
                     previous_fluid, residual_vector: np.ndarray,
                     pressure_limit: float, saturation_limit: float
                     ) -> ReservoirState:
        """Kəsmə (Appleyard) tək başına kifayət etmir.

        O, addımın YALNIZ UZUNLUĞUNU məhdudlaşdırır — istiqamətin
        FAYDALI olub-olmadığını (qalığı azaldıb-azaltmadığını) yoxlamır.
        Quyu hüceyrəsi ətrafında axının yuxarı axın (upstream) istiqaməti
        iterasiyadan-iterasiyaya dəyişəndə bu, Nyutonu bir neçə vəziyyət
        arasında SONSUZ DÖVRƏYƏ sala bilir (ölçüldü: CNV heç vaxt
        yığılmadan eyni 3 qiymət arasında rəqs edir).

        `CoupledNewtonSolver`-də tapılan və doğrulanmış həllin eynisi:
        addımı qalığı AZALDANA qədər yarıya böl (bax onun sənədləşməsi).
        """
        trial_residual, _, _ = self.R.residual(trial, previous, dt,
                                               previous_fluid)
        current_norm = self._scaled_norm(residual_vector, state, dt)
        trial_norm = self._scaled_norm(trial_residual, trial, dt)

        scale = 1.0
        for _ in range(10):
            if trial_norm <= current_norm * 0.999 or scale < 1.0 / 512.0:
                break
            scale *= 0.5
            trial = state.updated(
                delta * scale, self.sw_min, self.sw_max,
                max_pressure_change=pressure_limit,
                max_saturation_change=saturation_limit)
            trial_residual, _, _ = self.R.residual(trial, previous, dt,
                                                   previous_fluid)
            trial_norm = self._scaled_norm(trial_residual, trial, dt)
        return trial
