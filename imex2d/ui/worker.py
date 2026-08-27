"""Qt fon axını və proqres adapteri.

QtProgressReporter — IProgressReporter interfeysinin Qt implementasiyası.
Mühərrik Qt-ni tanımır; yalnız interfeysi tanıyır.
"""

from __future__ import annotations

from PyQt5.QtCore import QThread, pyqtSignal

from ..interfaces.services import IProgressReporter


class QtProgressReporter(IProgressReporter):

    def __init__(self, worker: "SimulationWorker"):
        self._worker = worker

    def report(self, fraction: float, message: str) -> bool:
        self._worker.progress.emit(float(fraction), message)
        return not self._worker.stop_requested


class MatchingWorker(QThread):
    """Uyğunlaşdırma axtarışı fon axınında.

    Hər qiymətləndirmə bir simulyasiyadır — onlarla dəqiqə çəkə bilər,
    ona görə interfeys bloklanmamalıdır.
    """

    progress = pyqtSignal(object)
    finished_ok = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, service, method: str, max_evaluations: int):
        super().__init__()
        self.service = service
        self.method = method
        self.max_evaluations = max_evaluations
        self.stop_requested = False

    def request_stop(self):
        self.stop_requested = True
        self.service.cancel()

    def run(self):
        import traceback
        try:
            result = self.service.run(
                method=self.method,
                max_evaluations=self.max_evaluations,
                progress=self._report)
            self.finished_ok.emit(result)
        except Exception:
            self.failed.emit(traceback.format_exc())

    def _report(self, evaluation) -> bool:
        self.progress.emit(evaluation)
        return not self.stop_requested


class SensitivityWorker(QThread):
    """Həssaslıq skanı fon axınında — hər addım bir simulyasiyadır."""

    progress = pyqtSignal(int, int)
    finished_ok = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, analyzer, method: str, metric: str,
                 step_fraction: float = 0.1):
        super().__init__()
        self.analyzer = analyzer
        self.method = method
        self.metric = metric
        self.step_fraction = step_fraction
        self.stop_requested = False

    def request_stop(self):
        self.stop_requested = True

    def run(self):
        import traceback
        try:
            if self.method == "Tornado":
                result = self.analyzer.run_tornado(metric=self.metric,
                                                   progress=self._report)
            else:
                result = self.analyzer.run_local(
                    metric=self.metric, step_fraction=self.step_fraction,
                    progress=self._report)
            self.finished_ok.emit(result)
        except Exception:
            self.failed.emit(traceback.format_exc())

    def _report(self, done: int, total: int) -> bool:
        self.progress.emit(done, total)
        return not self.stop_requested


class SimulationWorker(QThread):
    """Servis + model + konfiqurasiya qəbul edir, hesablamanı fon axınında aparır."""

    progress = pyqtSignal(float, str)
    finished_ok = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, service, model, config):
        super().__init__()
        self.service = service
        self.model = model
        self.config = config
        self.stop_requested = False

    def request_stop(self):
        self.stop_requested = True

    def run(self):
        import traceback
        try:
            result = self.service.run(self.model, self.config,
                                      QtProgressReporter(self))
            self.finished_ok.emit(result)
        except Exception:
            self.failed.emit(traceback.format_exc())
