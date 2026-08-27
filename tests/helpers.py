"""Testlər üçün ümumi model qurucuları.

Qəsdən pytest fixture-ları yox, adi funksiyalardır — bu sayədə testlər
həm pytest ilə, həm də run_tests.py ehtiyat runner-i ilə işləyir.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from imex2d.application.config import (OutputConfig, SimulationConfig,
                                       TimeSteppingConfig)
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder, five_spot
from imex2d.application.simulation_service import (
    ModelAwareSimulationService, SimulationService)
from imex2d.domain.initial import InitialConditions
from imex2d.domain.scal import CoreyParameters
from imex2d.domain.wells import (ControlMode, Perforation, Well, WellControl,
                                 WellType)
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter

SKIP_SLOW = bool(os.environ.get("IMEX_SKIP_SLOW"))

# ── Etalon dəyərlər: refaktorinqdən ƏVVƏLKİ core.py-nin nəticələri ──────
REFERENCE_FIVE_SPOT = {
    "recovery_factor": 16.840,   # %
    "steps": 4314,
    "ooip": 1029064.3,           # m3
}


def default_scal() -> CoreyParameters:
    return CoreyParameters()


def make_service(scal: CoreyParameters = None) -> SimulationService:
    """Provider-ləri modeldən quran servis — proqramın əsl davranışı.

    Sabit provider işlətsək, model SCAL/PVT parametrləri dəyişəndə
    nəticə səssizcə köhnə dəyərlərlə hesablanardı.
    """
    scal = scal or default_scal()
    return ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(scal))


def five_spot_model(nx=41, ny=41, dx=20.0, dy=20.0, dz=10.0,
                    porosity=0.22, permeability=150.0, scal=None):
    """Etalon five-spot modeli — reqressiya testinin bazası."""
    scal = scal or default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=dx, dy=dy, dz=dz,
        porosity=porosity, permx_base=permeability)
    return ReservoirModelBuilder().build(
        geological_model=geology,
        wells=five_spot(geology.grid),
        scal=scal,
        name="Test five-spot")


def one_dimensional_model(nx=120, dx=8.0, dy=100.0, dz=10.0,
                          porosity=0.20, permeability=200.0,
                          injection_rate=60.0, scal=None):
    """1D xətti sıxışdırma — Bukley-Leverett müqayisəsi üçün."""
    scal = scal or default_scal()
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=1, dx=dx, dy=dy, dz=dz,
        porosity=porosity, permx_base=permeability)
    wells = [
        Well("INJ", WellType.INJECTOR,
             WellControl(ControlMode.RATE, injection_rate), [Perforation(0, 0, 0)]),
        Well("PROD", WellType.PRODUCER,
             WellControl(ControlMode.BHP, 200.0), [Perforation(nx - 1, 0, 0)]),
    ]
    return ReservoirModelBuilder().build(
        geological_model=geology, wells=wells, scal=scal,
        initial=InitialConditions(datum_pressure=200.0,
                                  water_saturation=scal.swc),
        name="1D BL model")


def short_config(end_time=300.0, snapshots=10) -> SimulationConfig:
    return SimulationConfig(end_time=end_time,
                            output=OutputConfig(snapshot_count=snapshots))


def bl_config(end_time=250.0) -> SimulationConfig:
    return SimulationConfig(
        end_time=end_time,
        time_stepping=TimeSteppingConfig(max_dt=2.0, cfl_factor=0.4),
        output=OutputConfig(snapshot_count=2))
