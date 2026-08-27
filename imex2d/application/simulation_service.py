"""SimulationService — iş axınının "Simulyasiya" addımı.

Məsuliyyəti: modeli yoxlamaq, provider-ləri toplamaq, mühərriki qurmaq,
işə salmaq və nəticəni layihəyə yazmaq. Bu qat olmasaydı, həmin
məntiq yenidən UI-yə düşərdi.

DEPENDENCY INJECTION: provider-lər konstruktorda verilir. Servis heç bir
konkret sinfi özü yaratmır — hansı SCAL və ya PVT modulunun işlədiyini
tətbiqin giriş nöqtəsi (composition root) həll edir.
"""

from __future__ import annotations

import copy
from typing import Optional

from ..domain.reservoir_model import ReservoirModel
from ..logging_setup import get_logger
from ..interfaces.providers import (ICapillaryPressureProvider,
                                    IInitializationProvider, IPVTProvider,
                                    IRelativePermeabilityProvider)
from ..interfaces.services import ILinearSolver, IProgressReporter
from ..simulation.impes_engine import ImpesEngine
from ..simulation.implicit.engine import FullyImplicitEngine
from ..simulation.capillary import BrooksCoreyCapillaryProvider
from ..simulation.initialization.equilibrium import (
    EquilibriumInitializationProvider)
from ..simulation.linear_solver import ScipyCgIluSolver
from ..simulation.pvt.black_oil import BlackOilPVTProvider
from ..simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from ..simulation.scal_tables_provider import (
    TableCapillaryPressureProvider, TableRelativePermeabilityProvider)
from ..simulation.results import SimulationResult
from .config import SimulationConfig
from .project import Project, SimulationRun


LOG = get_logger(__name__)


class ModelValidationError(ValueError):
    def __init__(self, issues):
        super().__init__("; ".join(issues))
        self.issues = list(issues)


class SimulationService:

    def __init__(self,
                 relperm_provider: IRelativePermeabilityProvider,
                 linear_solver: Optional[ILinearSolver] = None,
                 pvt_provider: Optional[IPVTProvider] = None,
                 capillary_provider: Optional[ICapillaryPressureProvider] = None,
                 initialization_provider: Optional[IInitializationProvider] = None,
                 engine_factory=ImpesEngine):
        """`engine_factory` mühərrik seçimidir.

        İki mühərrik eyni interfeysi paylaşır və hər ikisi saxlanılır:
        kiçik modellərdə IMPES daha sürətlidir (hər addımı ucuzdur),
        böyüklərdə və uzun proqnozlarda isə fully implicit qazanır,
        çünki zaman addımı CFL ilə məhdud deyil.
        """
        self.relperm_provider = relperm_provider
        self.linear_solver = linear_solver
        self.pvt_provider = pvt_provider
        self.capillary_provider = capillary_provider
        self.initialization_provider = initialization_provider
        self.engine_factory = engine_factory

    def create_engine(self, model: ReservoirModel, config: SimulationConfig):
        issues = model.validate() + config.validate()
        if issues:
            LOG.error("Model yoxlamadan keçmədi: %s", "; ".join(issues))
            raise ModelValidationError(issues)
        LOG.info("Mühərrik qurulur: %s | %d hüceyrə | %d quyu | %.0f gün",
                 model.name, model.ncell, len(model.active_wells()),
                 config.end_time)
        solver = self.linear_solver or ScipyCgIluSolver(config.linear_solver)
        solver.reset()
        return self.engine_factory(
            model=model,
            config=config,
            relperm=self.relperm_provider,
            linear_solver=solver,
            pvt=self.pvt_provider,
            capillary=self.capillary_provider,
            initialization=self.initialization_provider,
        )

    def with_engine(self, engine_factory) -> "SimulationService":
        """Eyni provider-lərlə, başqa mühərriklə yeni servis."""
        service = copy.copy(self)
        service.engine_factory = engine_factory
        return service

    def run(self, model: ReservoirModel, config: SimulationConfig,
            reporter: Optional[IProgressReporter] = None) -> SimulationResult:
        return self.create_engine(model, config).run(reporter)

    def run_in_project(self, project: Project, model_name: str,
                       config: SimulationConfig,
                       reporter: Optional[IProgressReporter] = None) -> SimulationRun:
        run = project.new_run(model_name, config)
        run.status = "RUNNING"
        try:
            run.result = self.run(project.reservoir_models[model_name],
                                  config, reporter)
            run.status = "FINISHED" if run.result.converged else "FAILED"
        except Exception:
            run.status = "FAILED"
            run.result = None
            raise
        return run



class ModelAwareSimulationService(SimulationService):
    """Provider-ləri HƏR DƏFƏ modeldən qurur.

    Baza sinif provider-ləri konstruktorda alır və saxlayır. Bu, model
    dəyişmədikdə düzgündür, lakin modelin SCAL və PVT parametrləri
    dəyişəndə köhnə provider işlədilir və nəticə səssizcə yanlış olur.

    Problem history matching-də üzə çıxdı: optimallaşdırıcı `Sor`-u
    dəyişir, model yenilənir, amma nisbi keçiricilik adapteri köhnə
    `Sor` ilə qalır — nəticə heç dəyişmir və axtarış mənasız olur.

    Ona görə bu davranış artıq DEFOLTDUR (`app.py`-dəki alt-sinifdən
    buraya köçürüldü). Real region əsaslı SCAL modulu yazılanda
    (B4) provider modeli özü oxuyacaq və bu sinif silinəcək.
    """

    def create_engine(self, model, config):
        self.relperm_provider = self._relative_permeability(model)
        self.pvt_provider = (BlackOilPVTProvider(model.pvt_table)
                             if model.pvt_table is not None else None)
        self.capillary_provider = self._capillary(model)
        self.initialization_provider = (
            EquilibriumInitializationProvider(self.pvt_provider,
                                              self.capillary_provider)
            if model.initial_conditions.use_equilibration else None)

        return super().create_engine(model, config)

    @staticmethod
    def _relative_permeability(model):
        """Modeldə SCAL cədvəli varsa onu, yoxsa Corey düsturunu işlədir."""
        tables = getattr(model, "scal_tables", None)
        if tables is not None and len(tables):
            regions = (model.regions.region_id.values
                       if model.regions is not None else None)
            return TableRelativePermeabilityProvider(tables, regions)
        return CoreyRelativePermeabilityAdapter(model.scal_parameters)

    @staticmethod
    def _capillary(model):
        tables = getattr(model, "scal_tables", None)
        if tables is not None and len(tables):
            regions = (model.regions.region_id.values
                       if model.regions is not None else None)
            provider = TableCapillaryPressureProvider(tables, regions)
            if provider.has_capillary_pressure():
                return provider
            return None
        if model.capillary_parameters.enabled:
            return BrooksCoreyCapillaryProvider(model.capillary_parameters,
                                                model.scal_parameters)
        return None
