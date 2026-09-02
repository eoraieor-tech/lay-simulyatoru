"""Həssaslıq analizi — C6.

Uyğunlaşdırmadan (C5) əvvəl soruşulmalı sual: hansı parametr nəticəyə
nə qədər təsir edir? Cavab olmadan optimallaşdırıcı bütün parametrləri
eyni ciddiyyətlə axtarır, halbuki bəziləri praktik olaraq əhəmiyyətsizdir
— bu, axtarış fəzasını lüzumsuz genişləndirir və yığılmanı çətinləşdirir.

İKİ ÜSUL

    Tornado (diapazon əsaslı)
        Hər parametr öz TAM hədləri arasında dəyişdirilir, digərləri
        baza dəyərində saxlanılır (one-at-a-time). Nəticədəki
        "yayılma" parametrlərin nisbi əhəmiyyətini göstərir.
        Sual: "hansı parametr öz mümkün diapazonunda ən çox təsir edir?"

    Yerli elastiklik (nöqtə əsaslı)
        Baza nöqtəsi ətrafında kiçik ±addım, mərkəzi fərqlə törəmə.
        Elastiklik ölçüsüzdür: (ΔÇıxış/Çıxış) / (ΔParametr/Parametr).
        Sual: "hazırkı modeldə kiçik dəyişiklik nəyə təsir edir?"

Fərq vacibdir: geniş hədli, lakin lokal olaraq az təsirli parametr
Tornado-da yuxarıda, elastiklikdə aşağıda görünə bilər (və əksinə).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional

import numpy as np

from ..application.config import SimulationConfig
from ..application.simulation_service import SimulationService
from ..domain.reservoir_model import ReservoirModel
from ..logging_setup import get_logger
from ..simulation.results import SimulationResult
from .parameters import ModelModifier, ParameterSet

LOG = get_logger(__name__)


# ══════════════════════════════════════════════════════ çıxış ölçüləri

def _final_recovery_factor(result: SimulationResult) -> float:
    return result.final_recovery_factor


def _cumulative_oil(result: SimulationResult) -> float:
    return result.series.cumulative_oil[-1] if result.series.cumulative_oil else 0.0


def _final_water_cut(result: SimulationResult) -> float:
    return result.series.water_cut[-1] if result.series.water_cut else 0.0


def _breakthrough_time(result: SimulationResult) -> float:
    """Su gəlişi baş verməyibsə simulyasiya müddəti qaytarılır.

    `None` qaytarmaq həssaslıq hesablamasını sındırardı — "hələ gəlməyib"
    "sonsuz gecikmə" mənasında ən böyük dəyər kimi işlədilir.
    """
    value = result.breakthrough_time
    if value is not None:
        return float(value)
    return float(result.series.time[-1]) if result.series.time else 0.0


def _average_pressure(result: SimulationResult) -> float:
    values = result.series.average_pressure
    return values[-1] if values else 0.0


OUTPUT_METRICS: Dict[str, Callable[[SimulationResult], float]] = {
    "RF (%)": _final_recovery_factor,
    "Kumulyativ neft (m³)": _cumulative_oil,
    "Son sulaşma (%)": _final_water_cut,
    "Su gəlişi vaxtı (gün)": _breakthrough_time,
    "Son orta təzyiq (bar)": _average_pressure,
}


# ══════════════════════════════════════════════════════════ nəticələr

@dataclass
class ParameterSensitivity:
    """Bir parametrin bir çıxış ölçüsünə təsiri."""
    name: str
    description: str
    baseline_value: float
    baseline_output: float
    low_value: float
    low_output: float
    high_value: float
    high_output: float
    failed_low: bool = False
    failed_high: bool = False

    @property
    def swing(self) -> float:
        """Tornado eninin xam qiyməti — çıxışın diapazonu."""
        return abs(self.high_output - self.low_output)

    @property
    def direction_reversed(self) -> bool:
        """Parametr artanda çıxış AZALIR."""
        return self.high_output < self.low_output

    def elasticity_at(self, low_output: float, high_output: float,
                      step_fraction: float) -> float:
        centre = self.baseline_output
        if abs(centre) < 1e-12 or step_fraction <= 0:
            return 0.0
        relative_output = (high_output - low_output) / abs(centre)
        return relative_output / (2.0 * step_fraction)


@dataclass
class SensitivityReport:
    metric_name: str
    baseline_output: float
    items: List[ParameterSensitivity] = field(default_factory=list)
    failures: int = 0

    def sorted_by_swing(self) -> List[ParameterSensitivity]:
        return sorted(self.items, key=lambda item: -item.swing)

    def as_text(self) -> str:
        lines = [f"Çıxış: {self.metric_name}   (baza = {self.baseline_output:.4g})"]
        for item in self.sorted_by_swing():
            lines.append(
                f"  {item.name:<11} yayılma {item.swing:10.4g}   "
                f"[{item.low_output:.4g} .. {item.high_output:.4g}]"
                f"{'  (tərs)' if item.direction_reversed else ''}"
                f"{'  UĞURSUZ' if item.failed_low or item.failed_high else ''}")
        if self.failures:
            lines.append(f"  ({self.failures} qiymətləndirmə yığılmadı)")
        return "\n".join(lines)


# ══════════════════════════════════════════════════════════ analizator

class SensitivityAnalyzer:
    """Parametrləri bir-bir dəyişdirib çıxışa təsiri ölçür."""

    def __init__(self, base_model: ReservoirModel, parameters: ParameterSet,
                 simulation_service: SimulationService,
                 config: SimulationConfig):
        self.modifier = ModelModifier(base_model, parameters)
        self.parameters = parameters
        self.service = simulation_service
        self.config = config

    # ─────────────────────────────────────────────────── işə salma
    def _run(self, values: np.ndarray) -> Optional[SimulationResult]:
        try:
            model = self.modifier.apply(values)
            result = self.service.run(model, self.config)
            return result if result.converged else None
        except Exception as error:               # noqa: BLE001
            LOG.debug("Həssaslıq qiymətləndirməsi uğursuz: %s", error)
            return None

    # ═══════════════════════════════════════════════════════ tornado
    def run_tornado(self, metric: str = "RF (%)",
                    progress: Optional[Callable[[int, int], bool]] = None
                    ) -> SensitivityReport:
        """Hər parametr öz TAM hədləri arasında, digərləri baza dəyərində."""
        extractor = OUTPUT_METRICS[metric]
        baseline_values = self.parameters.initial_values
        baseline_result = self._run(baseline_values)
        if baseline_result is None:
            raise RuntimeError("Baza model yığılmadı — həssaslıq təhlili "
                              "üçün etibarlı başlanğıc nöqtəsi lazımdır.")
        baseline_output = extractor(baseline_result)

        report = SensitivityReport(metric_name=metric,
                                   baseline_output=baseline_output)
        total = len(self.parameters) * 2

        for index, definition in enumerate(self.parameters.definitions):
            low_values = baseline_values.copy()
            high_values = baseline_values.copy()
            low_values[index] = definition.minimum
            high_values[index] = definition.maximum

            low_result = self._run(low_values)
            if progress and not progress(index * 2 + 1, total):
                break
            high_result = self._run(high_values)
            if progress and not progress(index * 2 + 2, total):
                break

            failed_low = low_result is None
            failed_high = high_result is None
            if failed_low:
                report.failures += 1
            if failed_high:
                report.failures += 1

            low_output = (extractor(low_result) if not failed_low
                         else baseline_output)
            high_output = (extractor(high_result) if not failed_high
                          else baseline_output)

            report.items.append(ParameterSensitivity(
                name=definition.name, description=definition.description,
                baseline_value=baseline_values[index],
                baseline_output=baseline_output,
                low_value=definition.minimum, low_output=low_output,
                high_value=definition.maximum, high_output=high_output,
                failed_low=failed_low, failed_high=failed_high))

        return report

    # ═══════════════════════════════════════════════ yerli elastiklik
    def run_local(self, metric: str = "RF (%)", step_fraction: float = 0.1,
                  progress: Optional[Callable[[int, int], bool]] = None
                  ) -> SensitivityReport:
        """Baza nöqtəsi ətrafında kiçik ±addım — elastiklik.

        `step_fraction` [0,1] fəzasında normallaşdırılmış addımdır
        (bax `ParameterSet.to_unit`), fiziki vahiddə deyil — beləliklə
        müxtəlif vahidli parametrlər müqayisə oluna bilir.
        """
        extractor = OUTPUT_METRICS[metric]
        baseline_unit = self.parameters.to_unit(self.parameters.initial_values)
        baseline_result = self._run(self.parameters.initial_values)
        if baseline_result is None:
            raise RuntimeError("Baza model yığılmadı — həssaslıq təhlili "
                              "üçün etibarlı başlanğıc nöqtəsi lazımdır.")
        baseline_output = extractor(baseline_result)

        report = SensitivityReport(metric_name=metric,
                                   baseline_output=baseline_output)
        total = len(self.parameters) * 2

        for index, definition in enumerate(self.parameters.definitions):
            low_unit = baseline_unit.copy()
            high_unit = baseline_unit.copy()
            low_unit[index] = np.clip(baseline_unit[index] - step_fraction,
                                      0.0, 1.0)
            high_unit[index] = np.clip(baseline_unit[index] + step_fraction,
                                       0.0, 1.0)
            low_values = self.parameters.from_unit(low_unit)
            high_values = self.parameters.from_unit(high_unit)

            low_result = self._run(low_values)
            if progress and not progress(index * 2 + 1, total):
                break
            high_result = self._run(high_values)
            if progress and not progress(index * 2 + 2, total):
                break

            failed_low = low_result is None
            failed_high = high_result is None
            if failed_low:
                report.failures += 1
            if failed_high:
                report.failures += 1

            low_output = (extractor(low_result) if not failed_low
                         else baseline_output)
            high_output = (extractor(high_result) if not failed_high
                          else baseline_output)

            report.items.append(ParameterSensitivity(
                name=definition.name, description=definition.description,
                baseline_value=self.parameters.initial_values[index],
                baseline_output=baseline_output,
                low_value=low_values[index], low_output=low_output,
                high_value=high_values[index], high_output=high_output,
                failed_low=failed_low, failed_high=failed_high))

        return report
