"""PROJECT — iş axınının kök obyekti.

Project → GeologicalModel(-lər) → ReservoirModel(-lər) → SimulationRun(-lar)

Bir geoloji modeldən bir neçə rezervuar modeli (fərqli quyu sxemi,
fərqli SCAL) yaradıla bilər; bir rezervuar modelindən bir neçə
simulyasiya (fərqli müddət, fərqli idarəetmə) işə salına bilər.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional

from ..domain.geological_model import GeologicalModel
from ..domain.geology import GeologicalWell
from ..domain.reservoir_model import ReservoirModel
from .config import SimulationConfig


@dataclass
class SimulationRun:
    """Bir işə salınmanın qeydi: hansı model, hansı konfiqurasiya, nə nəticə."""
    run_id: str
    reservoir_model_name: str
    config: SimulationConfig
    result: object = None
    created_at: str = field(default_factory=lambda: datetime.now().isoformat(timespec="seconds"))
    status: str = "CREATED"


@dataclass
class Project:
    name: str = "Adsız layihə"
    geological_models: Dict[str, GeologicalModel] = field(default_factory=dict)
    reservoir_models: Dict[str, ReservoirModel] = field(default_factory=dict)
    runs: Dict[str, SimulationRun] = field(default_factory=dict)
    # Geologiya bölməsinin (2 ·) redaktə olunan quyu cədvəli — CSV-nin
    # əvəzidir. Cədvəl boşdursa mənbə avtomatik "synthetic" sayılır (bax
    # `geology_source`), ona görə ayrıca saxlanmır.
    geology_wells: List[GeologicalWell] = field(default_factory=list)
    geology_method: str = "Kriging (adi)"
    geology_params: Dict[str, float] = field(default_factory=dict)
    geology_defaults: Dict[str, float] = field(default_factory=lambda: {"sw": 0.30})
    _counter: int = 0

    @property
    def geology_source(self) -> str:
        return "wells" if self.geology_wells else "synthetic"

    def add_geological_model(self, model: GeologicalModel) -> GeologicalModel:
        self.geological_models[model.name] = model
        return model

    def add_reservoir_model(self, model: ReservoirModel) -> ReservoirModel:
        self.reservoir_models[model.name] = model
        return model

    def new_run(self, model_name: str, config: SimulationConfig) -> SimulationRun:
        if model_name not in self.reservoir_models:
            raise KeyError(f"Layihədə '{model_name}' rezervuar modeli yoxdur")
        self._counter += 1
        run = SimulationRun(f"RUN-{self._counter:03d}", model_name, config)
        self.runs[run.run_id] = run
        return run

    def latest_run(self) -> Optional[SimulationRun]:
        if not self.runs:
            return None
        return self.runs[sorted(self.runs)[-1]]

    def tree(self) -> List[str]:
        """UI-də layihə ağacını göstərmək üçün sadə mətn təsviri."""
        lines = [f"{self.name}"]
        for g in self.geological_models.values():
            lines.append(f"  Geoloji model: {g.name}")
        for r in self.reservoir_models.values():
            lines.append(f"  Rezervuar modeli: {r.name}  ({r.ncell} hüceyrə)")
        for run in self.runs.values():
            lines.append(f"    {run.run_id}: {run.reservoir_model_name} [{run.status}]")
        return lines
