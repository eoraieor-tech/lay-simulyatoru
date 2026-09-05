"""Layihə faylı (.imx) — saxlama və oxuma.

Format: gzip ilə sıxılmış JSON. Mətn formatı seçilməsinin səbəbi:
fayl açılıb baxıla bilir, versiyalar arasında fərq izlənə bilir və
xarici alətlə emal oluna bilir. Sıxılma anlıq şəkillərin (snapshot)
həcmini idarə edir.

Bu modul YALNIZ domain obyektlərini tanıyır — Qt, matplotlib və
hesablama mühərriki ilə əlaqəsi yoxdur.
"""

from __future__ import annotations

import gzip
import json
from datetime import datetime
from typing import Optional

import numpy as np

from ..domain.corner_point_geometry import CornerPointGeometry
from ..domain.facies_field import FaciesField
from ..domain.geological_model import Fault, GeologicalModel, Horizon
from ..domain.geology import GeologicalWell
from ..domain.geometry import CellGeometry
from ..domain.grid import CartesianGrid
from ..domain.initial import InitialConditions
from ..domain.properties import (FluidProperties, PermeabilityTensor, PropertyMap,
                                 RockProperties)
from ..domain.pvt import PVTTable
from ..domain.reservoir_model import ReservoirModel
from ..domain.scal import CapillaryParameters, CoreyParameters
from ..domain.structure import FaultReference, HorizonReference, RegionSet
from ..domain.units import FIELD, METRIC
from ..domain.wells import (ControlMode, Perforation, Phase, Well, WellControl,
                            WellType)
from ..geology.facies import FaciesVariogramParams
from ..geology.sgs import (DEFAULT_MIN_HARD_DATA_FOR_OWN_MODEL, FaciesPropertyConfig,
                           PropertyVariogramParams)
from ..simulation.results import SimulationResult, Snapshot, TimeSeries
from .config import (LinearSolverConfig, OutputConfig, SimulationConfig,
                     TimeSteppingConfig)
from .geology_service import ContinuousSGSConfig, FaciesBuildConfig
from .project import Project, SimulationRun

FORMAT_VERSION = 3
FILE_EXTENSION = ".imx"

# v1 → v2: geologiya quyu cədvəli (`Project.geology_wells` və
# `Well.perf_top`/`perf_bottom`) əlavə olundu. Hər ikisi köhnə fayllarda
# yoxdursa boş/None ilə doldurulur (aşağıda `.get()`), ona görə v1 faylı
# dəyişikliksiz açılır — versiya yalnız GƏLƏCƏK (bu proqramın anlamadığı)
# faylları rədd etmək üçün yoxlanılır, dəqiq bərabərliklə yox.
#
# v2 → v3: `GeologicalWell.data_layers_text` ("Data layları" sütunu, lay-
# məlumatlı rejim) əlavə olundu. v1/v2 faylında bu açar YOXDUR və boş
# mətnlə doldurulur (`GeologicalWell.from_dict`) — yəni köhnə layihə
# EYNİ (lay-məlumatsız) davranışla açılır. Versiya ona görə qaldırıldı ki,
# lay bəyanı OLAN faylı bu sahədən xəbərsiz KÖHNƏ build səssizcə açıb
# saxlayanda bəyanları İTİRMƏSİN — köhnə build v3 faylını AÇIQ rədd edir.

_UNIT_SYSTEMS = {"METRIC": METRIC, "FIELD": FIELD}


class ProjectFileError(Exception):
    """Fayl oxuna bilmədikdə və ya format uyğun gəlmədikdə."""


# ══════════════════════════════════════════════════════════ köməkçilər

def _array(values) -> list:
    return np.asarray(values, dtype=float).ravel().tolist()


def _property_map(prop: Optional[PropertyMap]) -> Optional[dict]:
    if prop is None:
        return None
    return {"name": prop.name, "unit": prop.unit, "values": _array(prop.values)}


def _property_map_from(data: Optional[dict]) -> Optional[PropertyMap]:
    if data is None:
        return None
    return PropertyMap(data["name"], np.asarray(data["values"], dtype=float),
                       data.get("unit", ""))


def _permeability_tensor_to_dict(tensor: Optional[PermeabilityTensor]) -> Optional[dict]:
    """Phase 2 — bütün 6 komponent (3 diaqonal + 3 off-diaqonal) ayrıca
    saxlanılır, heç biri itirilmir/təxmin edilmir."""
    if tensor is None:
        return None
    return {
        "kxx": _property_map(tensor.kxx),
        "kyy": _property_map(tensor.kyy),
        "kzz": _property_map(tensor.kzz),
        "kxy": _property_map(tensor.kxy),
        "kxz": _property_map(tensor.kxz),
        "kyz": _property_map(tensor.kyz),
    }


def _permeability_tensor_from_dict(data: Optional[dict]) -> Optional[PermeabilityTensor]:
    if data is None:
        return None
    return PermeabilityTensor(
        kxx=_property_map_from(data["kxx"]),
        kyy=_property_map_from(data["kyy"]),
        kzz=_property_map_from(data["kzz"]),
        kxy=_property_map_from(data.get("kxy")),
        kxz=_property_map_from(data.get("kxz")),
        kyz=_property_map_from(data.get("kyz")))


def _dataclass_to_dict(obj, fields) -> dict:
    return {field: getattr(obj, field) for field in fields}


def _facies_field_to_dict(facies: FaciesField) -> dict:
    """`FaciesField` (SIS realizasiyası) tam saxlanılır — kod, kateqoriya
    adları, tələb/reallaşan nisbətlər, variogram/anizotropluq metadatası,
    kondisioner statistikası, xəbərdarlıqlar. RNG-in ÖZÜ (transient) yox,
    yalnız NƏTİCƏ (`.codes`) və onu YARADAN parametrlər (`.seed` daxil,
    təkrarlana bilmə üçün) saxlanılır — bax modul §4 tələbi."""
    return {
        "name": facies.name,
        "codes": [int(c) for c in facies.codes.tolist()],
        "category_names": {str(k): v for k, v in facies.category_names.items()},
        "realization_id": facies.realization_id,
        "seed": facies.seed,
        "requested_proportions": {str(k): v for k, v in facies.requested_proportions.items()},
        "realized_proportions": {str(k): v for k, v in facies.realized_proportions.items()},
        "variogram_metadata": {str(k): v for k, v in facies.variogram_metadata.items()},
        "anisotropy_metadata": {str(k): v for k, v in facies.anisotropy_metadata.items()},
        "conditioning_data_stats": dict(facies.conditioning_data_stats),
        "warnings": list(facies.warnings),
    }


def _facies_field_from_dict(data: dict) -> FaciesField:
    return FaciesField(
        name=data["name"],
        codes=np.asarray(data["codes"], dtype=int),
        category_names={int(k): v for k, v in data.get("category_names", {}).items()},
        realization_id=data.get("realization_id", 0),
        seed=data.get("seed", 0),
        requested_proportions={int(k): v
                               for k, v in data.get("requested_proportions", {}).items()},
        realized_proportions={int(k): v
                              for k, v in data.get("realized_proportions", {}).items()},
        variogram_metadata={int(k): v for k, v in data.get("variogram_metadata", {}).items()},
        anisotropy_metadata={int(k): v for k, v in data.get("anisotropy_metadata", {}).items()},
        conditioning_data_stats=dict(data.get("conditioning_data_stats", {})),
        warnings=list(data.get("warnings", [])))


def _variogram_params_to_dict(vp) -> Optional[dict]:
    """`FaciesVariogramParams` VƏ `PropertyVariogramParams`-ın İKİSİ üçün
    də işləyir — sahə çoxluğu eynidir (bax `sgs.py`/`facies.py`)."""
    if vp is None:
        return None
    return {"model": vp.model, "nugget": vp.nugget, "sill": vp.sill, "range_": vp.range_,
            "range_v": vp.range_v, "azimuth_deg": vp.azimuth_deg,
            "range_minor": vp.range_minor}


def _facies_build_config_to_dict(config: FaciesBuildConfig) -> dict:
    return {
        "proportions": ({str(k): v for k, v in config.proportions.items()}
                       if config.proportions is not None else None),
        "category_names": ({str(k): v for k, v in config.category_names.items()}
                           if config.category_names is not None else None),
        "variograms": ({str(k): _variogram_params_to_dict(v)
                        for k, v in config.variograms.items()}
                       if config.variograms is not None else None),
        "seed": config.seed,
        "realization_id": config.realization_id,
        "search_radius": config.search_radius,
        "max_neighbors": config.max_neighbors,
        "min_neighbors": config.min_neighbors,
        "on_conflict": config.on_conflict,
    }


def _facies_build_config_from_dict(data: dict) -> FaciesBuildConfig:
    proportions = data.get("proportions")
    category_names = data.get("category_names")
    variograms = data.get("variograms")
    return FaciesBuildConfig(
        proportions=({int(k): v for k, v in proportions.items()}
                    if proportions is not None else None),
        category_names=({int(k): v for k, v in category_names.items()}
                        if category_names is not None else None),
        variograms=({int(k): FaciesVariogramParams(**v) for k, v in variograms.items()}
                   if variograms is not None else None),
        seed=data.get("seed", 0),
        realization_id=data.get("realization_id", 0),
        search_radius=data.get("search_radius"),
        max_neighbors=data.get("max_neighbors", 24),
        min_neighbors=data.get("min_neighbors", 1),
        on_conflict=data.get("on_conflict", "raise"))


def _facies_property_config_to_dict(fpc: FaciesPropertyConfig) -> dict:
    return {
        "variogram": _variogram_params_to_dict(fpc.variogram),
        "log_space": fpc.log_space,
        "bounds": (list(fpc.bounds) if fpc.bounds is not None else None),
    }


def _facies_property_config_from_dict(data: dict) -> FaciesPropertyConfig:
    variogram = data.get("variogram")
    bounds = data.get("bounds")
    return FaciesPropertyConfig(
        variogram=(PropertyVariogramParams(**variogram) if variogram is not None else None),
        log_space=data.get("log_space"),
        bounds=(tuple(bounds) if bounds is not None else None))


def _continuous_sgs_config_to_dict(config: ContinuousSGSConfig) -> dict:
    return {
        "variogram": _variogram_params_to_dict(config.variogram),
        "log_space": config.log_space,
        "bounds": (list(config.bounds) if config.bounds is not None else None),
        "seed": config.seed,
        "realization_id": config.realization_id,
        "search_radius": config.search_radius,
        "max_neighbors": config.max_neighbors,
        "min_neighbors": config.min_neighbors,
        "on_conflict": config.on_conflict,
        "conflict_tolerance": config.conflict_tolerance,
        "facies_field_name": config.facies_field_name,
        "facies_configs": ({str(k): _facies_property_config_to_dict(v)
                            for k, v in config.facies_configs.items()}
                           if config.facies_configs is not None else None),
        "min_hard_data_for_own_model": config.min_hard_data_for_own_model,
    }


def _continuous_sgs_config_from_dict(data: dict) -> ContinuousSGSConfig:
    variogram = data.get("variogram")
    bounds = data.get("bounds")
    facies_configs = data.get("facies_configs")
    return ContinuousSGSConfig(
        variogram=(PropertyVariogramParams(**variogram) if variogram is not None else None),
        log_space=data.get("log_space"),
        bounds=(tuple(bounds) if bounds is not None else None),
        seed=data.get("seed", 0),
        realization_id=data.get("realization_id", 0),
        search_radius=data.get("search_radius"),
        max_neighbors=data.get("max_neighbors", 24),
        min_neighbors=data.get("min_neighbors", 1),
        on_conflict=data.get("on_conflict", "raise"),
        conflict_tolerance=data.get("conflict_tolerance", 0.0),
        facies_field_name=data.get("facies_field_name"),
        facies_configs=({int(k): _facies_property_config_from_dict(v)
                        for k, v in facies_configs.items()}
                       if facies_configs is not None else None),
        min_hard_data_for_own_model=data.get("min_hard_data_for_own_model",
                                             DEFAULT_MIN_HARD_DATA_FOR_OWN_MODEL))


def _phase_or_water(value) -> Phase:
    """Naməlum faza adını (məs. v67-nin "GAS"-ı) suya endirir.

    Qaz fazası v69-da silindi, ona görə `Phase("GAS")` artıq
    `ValueError` atardı və köhnə layihə faylı açılmazdı.
    """
    try:
        return Phase(value) if value else Phase.WATER
    except ValueError:
        return Phase.WATER


def _known_fields(cls, data: dict) -> dict:
    """Yalnız dataclass-ın TANIDIĞI açarları buraxır.

    NİYƏ: köhnə layihə fayllarında artıq mövcud olmayan sahələr ola
    bilər (məs. v67-nin `gas_oil_contact` açarı). Süzgəc olmasa
    `Cls(**data)` `TypeError: unexpected keyword argument` atır və
    istifadəçinin köhnə `.imx` faylı ÜMUMİYYƏTLƏ açılmır. Naməlum
    açarı sükutla atmaq geriyə uyğunluğu saxlayır.
    """
    import dataclasses
    names = {f.name for f in dataclasses.fields(cls)}
    return {key: value for key, value in data.items() if key in names}


# ══════════════════════════════════════════════════════════ serializer

class ProjectSerializer:

    # ---------------------------------------------------------- public
    def save(self, project: Project, path: str,
             include_snapshots: bool = True) -> str:
        payload = {
            "format": "IMEX-2D project",
            "version": FORMAT_VERSION,
            "saved_at": datetime.now().isoformat(timespec="seconds"),
            "project": self.project_to_dict(project, include_snapshots),
        }
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False)
        return path

    def load(self, path: str) -> Project:
        try:
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                payload = json.load(handle)
        except OSError:
            # sıxılmamış fayl da qəbul edilir
            try:
                with open(path, "r", encoding="utf-8") as handle:
                    payload = json.load(handle)
            except Exception as exc:
                raise ProjectFileError(f"Fayl oxuna bilmədi: {exc}") from exc
        except json.JSONDecodeError as exc:
            raise ProjectFileError(f"Fayl formatı pozulub: {exc}") from exc

        version = payload.get("version")
        if not isinstance(version, int) or version > FORMAT_VERSION:
            raise ProjectFileError(
                f"Fayl versiyası {version}, bu proqram {FORMAT_VERSION} versiyasını oxuyur.")
        return self.project_from_dict(payload["project"])

    # -------------------------------------------------------- project
    def project_to_dict(self, project: Project, include_snapshots=True) -> dict:
        return {
            "name": project.name,
            "counter": project._counter,
            "geological_models": [self.geological_model_to_dict(m)
                                  for m in project.geological_models.values()],
            "reservoir_models": [self.reservoir_model_to_dict(m)
                                 for m in project.reservoir_models.values()],
            "runs": [self.run_to_dict(r, include_snapshots)
                     for r in project.runs.values()],
            "geology_wells": [w.to_dict() for w in project.geology_wells],
            "geology_method": project.geology_method,
            "geology_params": dict(project.geology_params),
            "geology_defaults": dict(project.geology_defaults),
            "geology_facies_configs": {name: _facies_build_config_to_dict(cfg)
                                       for name, cfg in project.geology_facies_configs.items()},
            "geology_sgs_configs": {name: _continuous_sgs_config_to_dict(cfg)
                                    for name, cfg in project.geology_sgs_configs.items()},
        }

    def project_from_dict(self, data: dict) -> Project:
        project = Project(name=data.get("name", "Layihə"))
        project._counter = int(data.get("counter", 0))
        for item in data.get("geological_models", []):
            project.add_geological_model(self.geological_model_from_dict(item))
        for item in data.get("reservoir_models", []):
            project.add_reservoir_model(self.reservoir_model_from_dict(item))
        for item in data.get("runs", []):
            run = self.run_from_dict(item)
            project.runs[run.run_id] = run

        if "geology_wells" in data:
            project.geology_wells = [GeologicalWell.from_dict(item)
                                     for item in data["geology_wells"]]
        else:
            project.geology_wells = self._migrate_geology_wells(project)
        project.geology_method = data.get("geology_method", project.geology_method)
        project.geology_params = data.get("geology_params", project.geology_params)
        project.geology_defaults = data.get("geology_defaults", project.geology_defaults)
        project.geology_facies_configs = {
            name: _facies_build_config_from_dict(item)
            for name, item in data.get("geology_facies_configs", {}).items()}
        project.geology_sgs_configs = {
            name: _continuous_sgs_config_from_dict(item)
            for name, item in data.get("geology_sgs_configs", {}).items()}
        return project

    @staticmethod
    def _migrate_geology_wells(project: Project) -> list:
        """v1 (`geology_wells` yoxdur) → v2 miqrasiyası.

        Köhnə layihələrdə heç vaxt geologiya quyu cədvəli olmayıb (yalnız
        CSV və ya sintetik generator). Son rezervuar modelinin quyuları
        varsa, onların i/j indeksindən TƏXMİNİ X/Y hesablanıb geologiya
        cədvəlinə bir dəfə köçürülür ki, istifadəçi boş cədvəllə
        qarşılaşmasın. Petrofiziki sahələr `None` qalır (bunlar əvvəllər
        heç vaxt quyu-səviyyəsində saxlanılmayıb, yalnız grid xəritəsi kimi).
        """
        if not project.reservoir_models:
            return []
        model = list(project.reservoir_models.values())[-1]
        # Perforasiya (i, j) → X/Y: hüceyrənin HƏQİQİ mərkəzindən
        # (Phase 5E). Nominal `(i+0.5)*dx` corner-point modeldə quyunu
        # yanlış koordinata miqrasiya edərdi.
        centroids = model.geometry.cell_centroid()
        wells = []
        for well in model.wells:
            if not well.perforations:
                continue
            p = well.perforations[0]
            wells.append(GeologicalWell(
                name=well.name, in_model=True,
                x=float(centroids[model.grid.index(p.i, p.j, 0)][0]),
                y=float(centroids[model.grid.index(p.i, p.j, 0)][1]),
                note="köhnə layihədən miqrasiya edilib"))
        return wells

    # ---------------------------------------------------- grid/geometry
    def _geometry_to_dict(self, geometry: CellGeometry) -> dict:
        """Kartezian sahələr + (corner-point olduqda) TƏPƏLƏR.

        `nodes` OLMADAN saxlanan corner-point modeli geri oxunanda
        SÜKUTLA bərabər bloka çevrilərdi — yəni layihəni saxlayıb açmaq
        həqiqi həndəsəni İTİRƏRDİ. Ona görə təpələr də yazılır; köhnə
        (nodes-suz) fayllar əvvəlki kimi `CellGeometry` kimi oxunur.
        """
        data = {
            "dx": geometry.dx, "dy": geometry.dy, "dz": geometry.dz.tolist(),
            "top_depth": geometry.top_depth,
            "top_depth_map": (None if geometry.top_depth_map is None
                              else _array(geometry.top_depth_map)),
        }
        nodes = getattr(geometry, "nodes", None)
        if nodes is not None:
            data["nodes"] = _array(nodes)          # (ncell·8·3,) düz siyahı
        return data

    def _geometry_from_dict(self, grid: CartesianGrid, data: dict) -> CellGeometry:
        surface = data.get("top_depth_map")
        surface = None if surface is None else np.asarray(surface, float)
        nodes = data.get("nodes")
        if nodes is not None:
            return CornerPointGeometry(
                grid=grid, dx=data["dx"], dy=data["dy"], dz=data["dz"],
                top_depth=data.get("top_depth", 0.0), top_depth_map=surface,
                nodes=np.asarray(nodes, float).reshape(grid.ncell, 8, 3))
        return CellGeometry(
            grid=grid, dx=data["dx"], dy=data["dy"], dz=data["dz"],
            top_depth=data.get("top_depth", 0.0), top_depth_map=surface)

    # -------------------------------------------------- geological model
    def geological_model_to_dict(self, model: GeologicalModel) -> dict:
        return {
            "name": model.name,
            "grid": {"nx": model.grid.nx, "ny": model.grid.ny, "nz": model.grid.nz},
            "geometry": self._geometry_to_dict(model.geometry),
            "property_maps": [_property_map(p) for p in model.property_maps.values()],
            "regions": {"region_id": _property_map(model.regions.region_id),
                        "names": {str(k): v for k, v in model.regions.names.items()}},
            "horizons": [{"name": h.name} for h in model.horizons],
            "faults": [{"name": f.name, "throw": f.throw, "dip": f.dip}
                       for f in model.faults],
            "coordinate_system": model.coordinate_system,
            "facies_fields": [_facies_field_to_dict(f) for f in model.facies_fields.values()],
        }

    def geological_model_from_dict(self, data: dict) -> GeologicalModel:
        grid = CartesianGrid(**data["grid"])
        geometry = self._geometry_from_dict(grid, data["geometry"])
        regions = RegionSet(
            _property_map_from(data["regions"]["region_id"]),
            {int(k): v for k, v in data["regions"].get("names", {}).items()})
        model = GeologicalModel(
            name=data["name"], grid=grid, geometry=geometry, regions=regions,
            horizons=[Horizon(h["name"]) for h in data.get("horizons", [])],
            faults=[Fault(f["name"], None, f.get("throw", 0.0), f.get("dip", 90.0))
                    for f in data.get("faults", [])],
            coordinate_system=data.get("coordinate_system", "LOCAL"))
        for item in data.get("property_maps", []):
            model.add_property(_property_map_from(item))
        for item in data.get("facies_fields", []):
            model.add_facies_field(_facies_field_from_dict(item))
        return model

    # --------------------------------------------------- reservoir model
    def reservoir_model_to_dict(self, model: ReservoirModel) -> dict:
        rock = model.rock
        return {
            "name": model.name,
            "grid": {"nx": model.grid.nx, "ny": model.grid.ny, "nz": model.grid.nz},
            "geometry": self._geometry_to_dict(model.geometry),
            "rock": {
                "porosity": _property_map(rock.porosity),
                "permx": _property_map(rock.permx),
                "permy": _property_map(rock.permy),
                "permz": _property_map(rock.permz),
                "net_to_gross": _property_map(rock.net_to_gross),
                "compressibility": rock.compressibility,
                "permeability_tensor": _permeability_tensor_to_dict(rock.permeability_tensor),
            },
            "fluids": _dataclass_to_dict(model.fluids, [
                "water_viscosity", "oil_viscosity", "water_fvf", "oil_fvf",
                "water_compressibility", "oil_compressibility",
                "water_density", "oil_density"]),
            "property_maps": [_property_map(p) for p in model.property_maps.values()],
            "regions": {"region_id": _property_map(model.regions.region_id),
                        "names": {str(k): v for k, v in model.regions.names.items()}},
            "faults": [_dataclass_to_dict(f, ["name", "source_id",
                                              "transmissibility_multiplier", "sealing",
                                              "axis", "plane_index", "range_a", "range_b"])
                       for f in model.fault_references],
            "horizons": [_dataclass_to_dict(h, ["name", "source_id", "role"])
                         for h in model.horizon_references],
            "wells": [self._well_to_dict(w) for w in model.wells],
            "initial_conditions": _dataclass_to_dict(model.initial_conditions, [
                "datum_depth", "datum_pressure", "water_saturation",
                "oil_water_contact", "equilibration_region",
                "use_equilibration"]),
            "scal": _dataclass_to_dict(model.scal_parameters, [
                "swc", "sor", "krw_end", "kro_end", "nw", "no"]),
            "capillary": _dataclass_to_dict(model.capillary_parameters, [
                "entry_pressure", "lambda_exponent", "max_pressure"]),
            "pvt": self._pvt_to_dict(model.pvt_table),
            "units": model.units.name,
            "source_geological_model": model.source_geological_model,
        }

    def reservoir_model_from_dict(self, data: dict) -> ReservoirModel:
        grid = CartesianGrid(**data["grid"])
        geometry = self._geometry_from_dict(grid, data["geometry"])
        rock_data = data["rock"]
        rock = RockProperties(
            porosity=_property_map_from(rock_data["porosity"]),
            permx=_property_map_from(rock_data["permx"]),
            permy=_property_map_from(rock_data["permy"]),
            permz=_property_map_from(rock_data.get("permz")),
            net_to_gross=_property_map_from(rock_data.get("net_to_gross")),
            compressibility=rock_data.get("compressibility", 4.5e-5),
            permeability_tensor=_permeability_tensor_from_dict(
                rock_data.get("permeability_tensor")))
        maps = {}
        for item in data.get("property_maps", []):
            prop = _property_map_from(item)
            maps[prop.name] = prop
        regions = RegionSet(
            _property_map_from(data["regions"]["region_id"]),
            {int(k): v for k, v in data["regions"].get("names", {}).items()})
        return ReservoirModel(
            name=data["name"], grid=grid, geometry=geometry, rock=rock,
            fluids=FluidProperties(**data["fluids"]),
            property_maps=maps, regions=regions,
            fault_references=[FaultReference(**f) for f in data.get("faults", [])],
            horizon_references=[HorizonReference(**h) for h in data.get("horizons", [])],
            wells=[self._well_from_dict(w) for w in data.get("wells", [])],
            initial_conditions=InitialConditions(
                **_known_fields(InitialConditions,
                                data["initial_conditions"])),
            scal_parameters=CoreyParameters(**data["scal"]),
            capillary_parameters=CapillaryParameters(**data.get("capillary", {})),
            pvt_table=self._pvt_from_dict(data.get("pvt")),
            units=_UNIT_SYSTEMS.get(data.get("units", "METRIC"), METRIC),
            source_geological_model=data.get("source_geological_model", ""))

    # ------------------------------------------------------------ wells
    @staticmethod
    def _well_to_dict(well: Well) -> dict:
        return {
            "name": well.name,
            "well_type": well.well_type.value,
            "control": {"mode": well.control.mode.value,
                        "target": well.control.target,
                        "injected_phase": well.control.injected_phase.value},
            "perforations": [{"i": p.i, "j": p.j, "k": p.k,
                              "open": p.open, "skin": p.skin,
                              "direction": p.direction}
                             for p in well.perforations],
            "radius": well.radius,
            "active": well.active,
            "perf_top": well.perf_top,
            "perf_bottom": well.perf_bottom,
        }

    @staticmethod
    def _well_from_dict(data: dict) -> Well:
        control = data["control"]
        return Well(
            name=data["name"],
            well_type=WellType(data["well_type"]),
            control=WellControl(ControlMode(control["mode"]), control["target"],
                                _phase_or_water(control.get("injected_phase"))),
            perforations=[Perforation(**p) for p in data.get("perforations", [])],
            radius=data.get("radius", 0.1),
            active=data.get("active", True),
            perf_top=data.get("perf_top"),
            perf_bottom=data.get("perf_bottom"))

    # -------------------------------------------------------------- pvt
    @staticmethod
    def _pvt_to_dict(table: Optional[PVTTable]) -> Optional[dict]:
        if table is None:
            return None
        return {
            "pressure": _array(table.pressure),
            "oil_fvf": _array(table.oil_fvf),
            "oil_viscosity": _array(table.oil_viscosity),
            "solution_gor": _array(table.solution_gor),
            "water_fvf": _array(table.water_fvf),
            "water_viscosity": _array(table.water_viscosity),
            "bubble_point": table.bubble_point,
            "rock_compressibility": table.rock_compressibility,
            "source": table.source,
        }

    @staticmethod
    def _pvt_from_dict(data: Optional[dict]) -> Optional[PVTTable]:
        if data is None:
            return None
        return PVTTable(**data)

    # -------------------------------------------------------------- run
    def run_to_dict(self, run: SimulationRun, include_snapshots=True) -> dict:
        return {
            "run_id": run.run_id,
            "reservoir_model_name": run.reservoir_model_name,
            "created_at": run.created_at,
            "status": run.status,
            "config": self._config_to_dict(run.config),
            "result": (None if run.result is None
                       else self._result_to_dict(run.result, include_snapshots)),
        }

    def run_from_dict(self, data: dict) -> SimulationRun:
        run = SimulationRun(
            run_id=data["run_id"],
            reservoir_model_name=data["reservoir_model_name"],
            config=self._config_from_dict(data["config"]),
            created_at=data.get("created_at", ""),
            status=data.get("status", "FINISHED"))
        if data.get("result") is not None:
            run.result = self._result_from_dict(data["result"])
        return run

    @staticmethod
    def _config_to_dict(config: SimulationConfig) -> dict:
        return {
            "end_time": config.end_time,
            "time_stepping": _dataclass_to_dict(config.time_stepping, [
                "initial_dt", "max_dt", "min_dt", "cfl_factor",
                "growth_factor", "max_steps"]),
            "linear_solver": _dataclass_to_dict(config.linear_solver, [
                "tolerance", "max_iterations", "preconditioner_refresh_steps",
                "ilu_drop_tolerance", "ilu_fill_factor", "fallback_to_direct"]),
            "output": _dataclass_to_dict(config.output, [
                "snapshot_count", "record_well_rates", "progress_every_n_steps"]),
        }

    @staticmethod
    def _config_from_dict(data: dict) -> SimulationConfig:
        return SimulationConfig(
            end_time=data["end_time"],
            time_stepping=TimeSteppingConfig(**data["time_stepping"]),
            linear_solver=LinearSolverConfig(**data["linear_solver"]),
            output=OutputConfig(**data["output"]))

    @staticmethod
    def _result_to_dict(result: SimulationResult, include_snapshots=True) -> dict:
        series = result.series
        return {
            "model_name": result.model_name,
            "grid_shape": list(result.grid_shape),
            "ooip": result.ooip,
            "steps": result.steps,
            "converged": result.converged,
            "message": result.message,
            "series": {name: list(getattr(series, name)) for name in (
                "time", "oil_rate", "water_rate", "water_injection_rate",
                "cumulative_oil", "cumulative_water", "water_cut",
                "average_pressure", "recovery_factor")},
            "well_oil_rate": {k: list(v) for k, v in result.well_oil_rate.items()},
            "well_water_rate": {k: list(v) for k, v in result.well_water_rate.items()},
            "snapshots": ([{"time": s.time,
                            "pressure": _array(s.pressure),
                            "water_saturation": _array(s.water_saturation)}
                           for s in result.snapshots] if include_snapshots else []),
        }

    @staticmethod
    def _result_from_dict(data: dict) -> SimulationResult:
        series = TimeSeries()
        for name, values in data["series"].items():
            setattr(series, name, list(values))
        shape = tuple(data["grid_shape"])
        result = SimulationResult(
            model_name=data.get("model_name", ""),
            grid_shape=shape,
            series=series,
            ooip=data.get("ooip", 0.0),
            steps=data.get("steps", 0),
            converged=data.get("converged", True),
            message=data.get("message", ""))
        result.well_oil_rate = {k: list(v) for k, v in data.get("well_oil_rate", {}).items()}
        result.well_water_rate = {k: list(v) for k, v in data.get("well_water_rate", {}).items()}
        result.snapshots = [
            Snapshot(time=s["time"],
                     pressure=np.asarray(s["pressure"], float).reshape(shape),
                     water_saturation=np.asarray(s["water_saturation"], float).reshape(shape))
            for s in data.get("snapshots", [])]
        return result
