"""Avtomatik uyğunlaşdırma — C5, mərhələ 3.

Axtarış [0, 1] fəzasında aparılır (bax `parameters.py`): müxtəlif
vahidli və diapazonlu parametrlər bir-birinə mane olmasın.

    unit vektor -> model -> simulyasiya -> uyğunsuzluq -> skalyar

UĞURSUZ QİYMƏTLƏNDİRMƏ
Optimallaşdırıcı hədlərin kənarını da sınayır və oradakı model
yığılmaya bilər. Belə hal İSTİSNA ATMAMALIDIR — axtarış dayanardı.
Əvəzinə böyük cərimə qaytarılır: nöqtə "pis" sayılır və axtarış davam
edir. Uğursuzluqların sayı hesabatda göstərilir.

KEŞLƏMƏ
Nelder-Mead eyni nöqtəni təkrar sınaya bilir. Hər qiymətləndirmə bir
simulyasiyadır, ona görə nəticələr keşlənir.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional

import numpy as np
from scipy import optimize

from ..application.config import SimulationConfig
from ..application.simulation_service import (ModelValidationError,
                                              SimulationService)
from ..domain.observations import ObservationSet
from ..domain.reservoir_model import ReservoirModel
from ..logging_setup import get_logger
from ..simulation.results import SimulationResult
from .mismatch import MismatchCalculator, MismatchReport
from .parameters import ModelModifier, ParameterSet

LOG = get_logger(__name__)

FAILURE_PENALTY = 1.0e6
"""Yığılmayan və ya yararsız model üçün cərimə.

Sonsuzluq YARAMAZ: Nelder-Mead simpleksi qura bilmir və
differential evolution seçim apara bilmir. Böyük, lakin sonlu
ədəd axtarışı həmin bölgədən uzaqlaşdırır.
"""


@dataclass
class Evaluation:
    """Bir qiymətləndirmənin qeydi."""
    iteration: int
    unit_values: np.ndarray
    values: np.ndarray
    mismatch: float
    succeeded: bool
    seconds: float
    message: str = ""


@dataclass
class MatchResult:
    parameters: ParameterSet
    best_values: np.ndarray
    best_mismatch: float
    best_report: Optional[MismatchReport] = None
    best_result: Optional[SimulationResult] = None
    history: List[Evaluation] = field(default_factory=list)
    initial_mismatch: float = float("nan")
    method: str = ""
    stopped_early: bool = False

    @property
    def evaluations(self) -> int:
        return len(self.history)

    @property
    def failures(self) -> int:
        return sum(1 for item in self.history if not item.succeeded)

    @property
    def improvement(self) -> float:
        """Uyğunsuzluğun neçə faiz yaxşılaşdığı."""
        if not np.isfinite(self.initial_mismatch) or self.initial_mismatch <= 0:
            return 0.0
        return (1.0 - self.best_mismatch / self.initial_mismatch) * 100.0

    @property
    def convergence_curve(self) -> np.ndarray:
        """Hər addımda o ana qədərki ən yaxşı nəticə."""
        best = float("inf")
        curve = []
        for item in self.history:
            best = min(best, item.mismatch)
            curve.append(best)
        return np.asarray(curve, float)

    def as_dict(self) -> Dict[str, float]:
        return self.parameters.as_dict(self.best_values)

    def summary(self) -> str:
        lines = [
            f"Üsul: {self.method}",
            f"Qiymətləndirmə: {self.evaluations}  "
            f"(uğursuz: {self.failures})",
            f"Uyğunsuzluq: {self.initial_mismatch:.5f} -> "
            f"{self.best_mismatch:.5f}   "
            f"({self.improvement:+.1f} % yaxşılaşma)",
            "Ən yaxşı parametrlər:",
        ]
        for name, value in self.as_dict().items():
            lines.append(f"  {name:<11} {value:>12.5f}")
        return "\n".join(lines)


class HistoryMatchingService:
    """Uyğunsuzluğu minimuma endirən axtarış."""

    METHODS = ("Nelder-Mead", "Differential Evolution", "Powell")

    def __init__(self,
                 base_model: ReservoirModel,
                 parameters: ParameterSet,
                 observations: ObservationSet,
                 simulation_service: SimulationService,
                 config: SimulationConfig,
                 calculator: Optional[MismatchCalculator] = None):
        self.modifier = ModelModifier(base_model, parameters)
        self.parameters = parameters
        self.observations = observations
        self.service = simulation_service
        self.config = config
        self.calculator = calculator or MismatchCalculator()

        self.history: List[Evaluation] = []
        self._cache: Dict[tuple, float] = {}
        self._best_mismatch = float("inf")
        self._best_unit: Optional[np.ndarray] = None
        self._cancelled = False
        self._progress: Optional[Callable[[Evaluation], bool]] = None

    # ═══════════════════════════════════════════ qiymətləndirmə
    def evaluate(self, unit_values: np.ndarray) -> float:
        """[0,1] vektor -> uyğunsuzluq. Optimallaşdırıcı bunu çağırır."""
        if self._cancelled:
            return FAILURE_PENALTY

        unit_values = np.clip(np.asarray(unit_values, float), 0.0, 1.0)
        key = tuple(np.round(unit_values, 8))
        if key in self._cache:
            return self._cache[key]

        values = self.parameters.from_unit(unit_values)
        started = time.perf_counter()
        message = ""
        succeeded = True

        try:
            model = self.modifier.apply(values)
            result = self.service.run(model, self.config)
            if not result.converged:
                succeeded = False
                message = result.message
                mismatch = FAILURE_PENALTY
            else:
                report = self.calculator.evaluate(result, self.observations)
                mismatch = report.total
                if not np.isfinite(mismatch):
                    succeeded = False
                    message = "Müşahidə ilə örtüşmə yoxdur."
                    mismatch = FAILURE_PENALTY
                elif mismatch < self._best_mismatch:
                    self._best_mismatch = mismatch
                    self._best_unit = unit_values.copy()
        except (ModelValidationError, ValueError, FloatingPointError,
                RuntimeError) as error:
            # Optimallaşdırıcı hədlərin kənarını sınayır — orada model
            # yararsız ola bilər. İstisna axtarışı dayandırmamalıdır.
            succeeded = False
            message = str(error)
            mismatch = FAILURE_PENALTY

        evaluation = Evaluation(
            iteration=len(self.history) + 1,
            unit_values=unit_values.copy(),
            values=values.copy(),
            mismatch=float(mismatch),
            succeeded=succeeded,
            seconds=time.perf_counter() - started,
            message=message)
        self.history.append(evaluation)
        self._cache[key] = float(mismatch)

        if self._progress is not None and not self._progress(evaluation):
            self._cancelled = True
        return float(mismatch)

    def cancel(self) -> None:
        self._cancelled = True

    # ═══════════════════════════════════════════════════ axtarış
    def run(self, method: str = "Nelder-Mead",
            max_evaluations: int = 60,
            progress: Optional[Callable[[Evaluation], bool]] = None,
            seed: int = 0) -> MatchResult:
        self.history.clear()
        self._cache.clear()
        self._best_mismatch = float("inf")
        self._best_unit = None
        self._cancelled = False
        self._progress = progress

        start_unit = self.parameters.to_unit(self.parameters.initial_values)
        initial = self.evaluate(start_unit)
        LOG.info("Uyğunlaşdırma başladı: %s, başlanğıc uyğunsuzluq %.5f",
                 method, initial)

        try:
            self._search(method, start_unit, max_evaluations, seed)
        except Exception as error:                      # optimallaşdırıcının özü
            LOG.warning("Optimallaşdırıcı dayandı: %s", error)

        best_unit = (self._best_unit if self._best_unit is not None
                     else start_unit)
        best_values = self.parameters.from_unit(best_unit)
        report, result = self._final_evaluation(best_values)

        match = MatchResult(
            parameters=self.parameters,
            best_values=best_values,
            best_mismatch=self._best_mismatch,
            best_report=report,
            best_result=result,
            history=list(self.history),
            initial_mismatch=initial,
            method=method,
            stopped_early=self._cancelled)
        LOG.info("Uyğunlaşdırma bitdi: %.5f -> %.5f (%d qiymətləndirmə, "
                 "%d uğursuz)", initial, match.best_mismatch,
                 match.evaluations, match.failures)
        self._progress = None
        return match

    def _search(self, method: str, start_unit: np.ndarray,
                max_evaluations: int, seed: int) -> None:
        dimensions = len(self.parameters)
        bounds = [(0.0, 1.0)] * dimensions

        if method == "Differential Evolution":
            # populyasiya × nəsil ≈ büdcə
            population = max(4, min(15, max_evaluations // (dimensions * 3)))
            generations = max(1, max_evaluations
                              // max(population * dimensions, 1))
            optimize.differential_evolution(
                self.evaluate, bounds, maxiter=generations,
                popsize=population, seed=seed, polish=False,
                init="sobol" if dimensions <= 21 else "latinhypercube",
                tol=1e-6)
            return

        if method == "Powell":
            optimize.minimize(
                self.evaluate, start_unit, method="Powell", bounds=bounds,
                options={"maxfev": max_evaluations, "xtol": 1e-3,
                         "ftol": 1e-4})
            return

        # Nelder-Mead: başlanğıc simpleksi hədlərin daxilində qururuq
        simplex = self._initial_simplex(start_unit)
        optimize.minimize(
            self.evaluate, start_unit, method="Nelder-Mead", bounds=bounds,
            options={"maxfev": max_evaluations, "initial_simplex": simplex,
                     "xatol": 1e-3, "fatol": 1e-4, "adaptive": True})

    @staticmethod
    def _initial_simplex(start_unit: np.ndarray,
                         step: float = 0.15) -> np.ndarray:
        """n+1 təpəli simpleks — hər təpə bir parametr üzrə sürüşdürülüb.

        Defolt simpleks nisbi addım işlədir və [0,1] fəzasında sıfıra
        yaxın parametrlər üçün praktiki olaraq hərəkətsiz qalır.
        """
        dimensions = start_unit.size
        simplex = np.repeat(start_unit[None, :], dimensions + 1, axis=0)
        for index in range(dimensions):
            offset = step if start_unit[index] <= 0.5 else -step
            simplex[index + 1, index] = np.clip(
                start_unit[index] + offset, 0.0, 1.0)
        return simplex

    def _final_evaluation(self, values: np.ndarray):
        """Ən yaxşı variantın tam nəticəsi — qrafiklər üçün."""
        try:
            result = self.service.run(self.modifier.apply(values), self.config)
            return self.calculator.evaluate(result, self.observations), result
        except Exception as error:
            LOG.warning("Ən yaxşı variant təkrar hesablana bilmədi: %s", error)
            return None, None
