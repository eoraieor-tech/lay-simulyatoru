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
from ..simulation.capillary import BrooksCoreyCapillaryProvider
from ..simulation.initialization.equilibrium import (
    EquilibriumInitializationProvider)
from ..simulation.initialization.saturation_map import (
    SaturationMapInitializationProvider, SaturationMapOverride,
    water_saturation_from_map)
from ..simulation.linear_solver import ScipyCgIluSolver
from ..simulation.wellbore.hydraulics import WellboreHydraulics
from ..simulation.pvt.black_oil import BlackOilPVTProvider
from ..simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from ..simulation.scal_tables_provider import (
    TableCapillaryPressureProvider, TableRelativePermeabilityProvider)
from ..simulation.results import SimulationResult
from .config import MPFA_O, SimulationConfig
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
        # Axın sxemi də jurnala yazılır: eyni model iki sxemlə fərqli
        # nəticə verir, ona görə jurnal onsuz OXUNMAZ olur (B1-dən sonra
        # aşkarlandı — RUN-002 diaqnozunda hansı sxemin işlədildiyi
        # bilinmirdi).
        LOG.info("Mühərrik qurulur: %s | %d hüceyrə | %d quyu | %.0f gün | %s",
                 model.name, model.ncell, len(model.active_wells()),
                 config.end_time, config.flux_scheme)
        solver = self.linear_solver or ScipyCgIluSolver(config.linear_solver)
        solver.reset()
        self._reject_incompatible_engine(config)
        self._reject_thp_with_impes(model)
        try:
            return self.engine_factory(
                model=model,
                config=config,
                relperm=self.relperm_provider,
                linear_solver=solver,
                pvt=self.pvt_provider,
                capillary=self.capillary_provider,
                initialization=self.initialization_provider,
                flux_discretization=self._flux_discretization(config),
            )
        except NotImplementedError as exc:
            # MPFA-O-nun HƏLƏ dəstəkləmədiyi model xüsusiyyəti (fay, NNC,
            # qeyri-aktiv hüceyrə — bax `residual.py::_reject_unsupported`).
            # Mesaj ORİJİNAL SAXLANILIR: orada məhz hansı xüsusiyyətin
            # maneə olduğu yazılıb, onu öz sözümüzlə əvəz etsək
            # istifadəçi səbəbi İTİRƏRDİ. Yalnız çevrilir ki, proqram
            # çökməsin (UI `ModelValidationError` gözləyir).
            LOG.error("Diskretizasiya modeli dəstəkləmir: %s", exc)
            raise ModelValidationError([str(exc)]) from exc

    def _reject_incompatible_engine(self, config: SimulationConfig) -> None:
        """IMPES + MPFA-O — İSTİFADƏÇİ DİLİNDƏ imtina.

        `impes_engine.py::_reject_multipoint_impes()` bunu onsuz da
        rədd edir, amma `NotImplementedError` texniki mesajdır. Burada
        seçim mühərrik qurulmazdan ƏVVƏL yoxlanılır ki, istifadəçi nə
        etməli olduğunu göstərən mesaj alsın.
        """
        if not config.uses_multipoint_flux:
            return
        if self.engine_factory is not ImpesEngine:
            return
        raise ModelValidationError([
            "MPFA-O yalnız tam implicit (Nyuton) mühərriklə işləyir. "
            "Ədədi parametrlər tabında ya hesablama sxemini "
            "«Fully implicit» edin, ya da diskretizasiyanı TPFA seçin."])

    def _reject_thp_with_impes(self, model) -> None:
        """IMPES + THP — İSTİFADƏÇİ DİLİNDƏ imtina (B4-B).

        `impes_engine._reject_thp_impes()` onsuz da rədd edir, amma
        texniki dildə. Burada quyuların adı ilə nə etməli olduğu deyilir.
        """
        if self.engine_factory is not ImpesEngine:
            return
        from ..domain.wells import ControlMode
        names = [well.name for well in model.active_wells()
                 if well.control.mode is ControlMode.THP]
        if names:
            raise ModelValidationError([
                "THP idarəsi yalnız tam implicit (Nyuton) mühərriklə işləyir "
                f"({', '.join(names)}). Ədədi parametrlər tabında hesablama "
                "sxemini «Fully implicit» edin və ya quyunu BHP rejiminə "
                "qaytarın."])
        # B7 addım 2: RATE quyusunun BHP limiti də yalnız tam implicit
        # mühərrikdədir (bax `impes_engine._reject_bhp_limit_impes`)
        limited = [well.name for well in model.active_wells()
                   if well.control.mode is ControlMode.RATE
                   and well.control.bhp_limit is not None]
        if limited:
            raise ModelValidationError([
                "RATE quyusunun BHP limiti yalnız tam implicit (Nyuton) "
                f"mühərriklə işləyir ({', '.join(limited)}). Ədədi parametrlər "
                "tabında hesablama sxemini «Fully implicit» edin və ya limiti "
                "silin."])

    @staticmethod
    def _flux_discretization(config: SimulationConfig):
        """Konfiqurasiyadan axın diskretizasiyası.

        MPFA-O nüvəsi (`imex2d/discretization/`) və onun tam implicit
        mühərrikə qoşulması (Phase 5B-2) ARTIQ mövcud idi — çatışmayan
        yeganə halqa məhz bu seçim idi: servis `flux_discretization`
        ötürmədiyi üçün istifadəçi MPFA-O-nu heç cür işə sala bilmirdi
        (bax `ISH_HESABATI.md` → Seans 3, `ICRA_PLANI.md` → B1).

        SƏRHƏD BAĞLANIŞI — `NEUMANN_ZERO`, defolt `DIRICHLET` DEYİL.
        Səbəb ikiqatdır:
          1. Simulyatorun özü onsuz da AXINSIZ (no-flow) xarici sərhəd
             tətbiq edir — TPFA yolunda da belədir, yəni bu seçim
             fizikanı DƏYİŞMİR, mövcud şərti təkrarlayır.
          2. Qalıq qatı Dirichlet üçün sərhəd π dəyərlərini hələ
             ötürmür (Phase 5B-2 məhdudiyyəti, `residual.py`) — defolt
             ilə qurulsaydı, hər MPFA-O seçimi dərhal xəta verərdi.
        Eyni seçim doğrulama testlərində də işlədilir
        (`tests/test_phase_d_mpfa_integration.py::NEUMANN`).

        İdxal QƏSDƏN funksiyanın içindədir: TPFA yolu MPFA-O modulunu
        ümumiyyətlə yükləmir, yəni köhnə davranışın başlanğıc xərci
        dəyişmir.
        """
        if config.flux_scheme != MPFA_O:
            return None            # → mühərrik `default_flux_discretization()` işlədir
        from ..discretization import MPFAOBoundaryClosure, MPFAODiscretization
        return MPFAODiscretization(closure=MPFAOBoundaryClosure.NEUMANN_ZERO)

    def with_engine(self, engine_factory) -> "SimulationService":
        """Eyni provider-lərlə, başqa mühərriklə yeni servis."""
        service = copy.copy(self)
        service.engine_factory = engine_factory
        return service

    def run(self, model: ReservoirModel, config: SimulationConfig,
            reporter: Optional[IProgressReporter] = None) -> SimulationResult:
        result = self.create_engine(model, config).run(reporter)
        self._annotate_wellbore(model, result)
        return result

    def _annotate_wellbore(self, model: ReservoirModel,
                           result: SimulationResult) -> None:
        """B4 — quyu başı təzyiqi (THP) POST-PROSES kimi hesablanır.

        Mühərrik toxunulmur: BHP rejimli quyuda quyu dibi təzyiqi
        onsuz da sabitdir, debit sıraları isə nəticədə artıq var.
        Bax `simulation/wellbore/hydraulics.py`.

        HEÇ BİR quyuda lülə həndəsəsi yoxdursa dərhal qayıdır — yəni
        mövcud modellərin davranışı və sürəti DƏYİŞMİR.
        """
        if not any(getattr(well, "tubing", None) is not None
                   for well in model.wells):
            return
        hydraulics = WellboreHydraulics(
            pvt=getattr(self, "pvt_provider", None), fluids=model.fluids)
        hydraulics.annotate(model, result)

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
        self.pvt_provider = self._build_pvt_provider(model)
        self.capillary_provider = self._capillary(model)
        self.initialization_provider = self._initialization(
            model, self.pvt_provider, self.capillary_provider)

        if self._gas_phase_requested():
            return self._create_three_phase_engine(model, config)
        return super().create_engine(model, config)

    def _gas_phase_requested(self) -> bool:
        return (self.pvt_provider is not None
                and self.pvt_provider.has_gas_phase())

    def _create_three_phase_engine(self, model, config):
        """A7 — PVT-də qaz xassələri varsa ÜÇ FAZALI mühərrik.

        İstifadəçi "Qaz fazasını aktivləşdir" qutusunu işarələyibsə,
        IMPES/IMPLICIT seçimindən ASILI OLMAYARAQ üç fazalı mühərrik
        işlədilir — hər ikisi qazı dəstəkləmir.

        Bax `A7_PLAN.md`: nüvə doğrulanıb, lakin quyu öz BHP hədəfinə
        çox yaxınlaşan hallarda simulyasiya vaxtından əvvəl (yığılmadan)
        dayana bilər. Bu, xəta kimi ATILMIR —
        `SimulationResult.converged=False` və aydın mesajla qaytarılır,
        son yığılmış nöqtəyə qədərki nəticələr saxlanılır.

        İKİ FƏRQ (v69-dan əvvəlki orijinal koda nisbətən):

        1. `engine_factory` DAİMİ dəyişdirilmir. Orijinal kod
           `self.engine_factory = ThreePhaseSimulationEngine` yazırdı —
           bu, servis TƏKRAR işlədiləndə (history matching, həssaslıq)
           qaz söndürülsə belə üç fazalı mühərriki seçili saxlayardı və
           `ValueError` verərdi. Burada seçim yalnız BU çağırışa aiddir.
        2. `flux_discretization` ötürülmür — üç fazalı mühərrik B1-dən
           sonra əlavə olunan bu açar sözü qəbul edir, lakin yalnız
           TPFA-nı dəstəkləyir (bax aşağıdakı yoxlama).
        """
        if config.uses_multipoint_flux:
            raise ModelValidationError([
                "MPFA-O hələ üç fazalı (qazlı) mühərriklə işləmir. "
                "Ya PVT-də qaz fazasını söndürün, ya da Ədədi "
                "parametrlər tabında diskretizasiyanı TPFA seçin."])

        from ..simulation.implicit.three_phase_engine import (
            ThreePhaseSimulationEngine)
        from ..simulation.stone_relperm import StoneRelativePermeabilityProvider

        gas_scal = getattr(model, "gas_scal_parameters", None)
        if gas_scal is None:
            # PVT-də qaz aktiv, SCAL-da qaz əyriləri seçilməyib —
            # mühərrik Stone relperm gözləyir, yalnız Corey alsaydı
            # ÇÖKƏRDİ (A7_PLAN-da qeyd olunmuş səhv). Defolt işlədilir.
            from ..domain.scal import GasCoreyParameters
            gas_scal = GasCoreyParameters()
            LOG.info("Qaz-neft SCAL verilməyib — defolt "
                     "GasCoreyParameters() işlədilir.")
        self.relperm_provider = StoneRelativePermeabilityProvider.from_corey(
            model.scal_parameters, gas_scal)

        previous_factory = self.engine_factory
        self.engine_factory = ThreePhaseSimulationEngine
        try:
            return super().create_engine(model, config)
        finally:
            self.engine_factory = previous_factory

    @staticmethod
    def _build_pvt_provider(model):
        """PVT provider — qaz modelləşdirilmirsə ölü-neft düzəlişi ilə.

        `dead_oil_below_bubble_point=True` HƏMİŞƏ verilir, çünki
        provider onu YALNIZ qaz sütunları OLMAYAN cədvələ tətbiq edir
        (bax `BlackOilPVTProvider._dead_oil_fvf`). Yəni:

          * qaz sütunlu cədvəl → üç fazalı mühərrik → cədvəl
            TOXUNULMUR (doymuş budaq orada DÜZGÜNDÜR);
          * qaz sütunsuz cədvəl → iki fazalı mühərrik → doyma
            təzyiqindən aşağı Bo ölü-neft budağı ilə əvəz olunur,
            əks halda Jakobian təkləşir və Nyuton donur (B3-A).

        Bu seçim QƏSDƏN burada, `application` qatındadır: "qaz
        modelləşdirilirmi" sualı iş axınının qərarıdır, cədvəl
        interpolyasiyasının deyil.
        """
        if model.pvt_table is None:
            return None
        return BlackOilPVTProvider(model.pvt_table,
                                   dead_oil_below_bubble_point=True)

    @staticmethod
    def _initialization(model, pvt_provider=None, capillary_provider=None):
        """İki bayraqdan ilkin şərt provider-i (bax `domain/initial.py`).

            use_equilibration  use_saturation_map  Provider
            ─────────────────  ──────────────────  ──────────────────────────────
            False              False               None (mühərrik skalyar işlədir)
            False              True                SaturationMapInitializationProvider
            True               False               EquilibriumInitializationProvider
            True               True                SaturationMapOverride(Equilibrium)

        Xəritə TƏLƏB OLUNUB, amma oxunmursa — SƏSSİZCƏ skalyara qayıtmırıq:
        istifadəçi bayraq qaldırıbsa, o, xəritə gözləyir; skalyar nəticə
        onun bilmədiyi başqa bir hesabdır.
        """
        ic = model.initial_conditions
        inner = (EquilibriumInitializationProvider(pvt_provider, capillary_provider)
                 if ic.use_equilibration else None)
        if not ic.use_saturation_map:
            return inner

        try:
            water_saturation_from_map(model)
        except ValueError as exc:
            LOG.error("İlkin Sw xəritəsi oxunmadı: %s", exc)
            raise ModelValidationError([str(exc)]) from exc

        if inner is None:
            return SaturationMapInitializationProvider()
        return SaturationMapOverride(inner)

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
