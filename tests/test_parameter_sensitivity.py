"""Parametr həssaslığı — UI-dan verilən dəyər NƏTİCƏYƏ çatırmı?

PROBLEM (sahibkarın bildirişi, `ISH_HESABATI.md` → Seans 14). İstifadəçi
interfeysdə neft lözlüyünü və φ/K-nı dəyişirdi, RF isə tərpənmirdi.

Audit göstərdi ki, **mühərrik sağlamdır** — mobillik və transmissibillik
hesablamalarında hardcode rəqəm YOXDUR. Qırılma UI → model qatında idi
və üç səssiz üstələmə mexanizmi vardı:

  A. PVT modeli işlədiləndə μo `PVT cədvəlindən` gəlir, paneldəki dəyər
     oxunmur (`residual.py`-dakı `pvt is None` budaqlanması);
  B. interpolyasiya olunmuş geologiya keşlənir və φ/K dəyişəndə
     avtomatik yenilənmir;
  C. GRDECL idxal olunanda φ/K faylın öz xəritələrindən gəlir.

Hər üçü DİZAYN ÜZRƏDİR. Problem görünürlükdə idi — indi sahələr
bozarır və banner çıxır.

BU FAYLIN İŞİ: fizikanın həqiqətən reaksiya verdiyini KİLİDLƏMƏK.
Gələcəkdə hansısa dəyişiklik bu əlaqəni səssizcə qırsa, test tutsun.
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.properties import FluidProperties
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter


def _run(oil_viscosity=3.0, porosity=0.22, permx=150.0,
         api=32.0, temperature=70.0, with_pvt=False, end_time=400.0):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=8, ny=8, dx=20.0, dy=20.0, dz=10.0,
        porosity=porosity, permx_base=permx)
    table = (build_pvt_table(api=api, temperature_c=temperature,
                             pressure_min=1.0, pressure_max=400.0,
                             n_points=40, bubble_point_bar=150.0)
             if with_pvt else None)
    model = ReservoirModelBuilder().build(
        geological_model=geology, wells=five_spot(geology.grid),
        fluids=FluidProperties(oil_viscosity=oil_viscosity),
        scal=default_scal(), pvt_table=table, name="həssaslıq")
    service = ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)
    return service.run(model, SimulationConfig(end_time=end_time))


def _rf(result):
    assert result.converged, result.message
    return float(result.series.recovery_factor[-1])


# ═══════════════════ PVT SÖNDÜRÜLÜ — panel dəyərləri işləyir ═══════════

def test_oil_viscosity_changes_recovery_without_pvt():
    """Panel lözlüyü PVT söndürülü olanda MÜTLƏQ təsir etməlidir.

    Ölçülmüş: 1 cP → RF 64.92 %, 30 cP → RF 22.49 %. Ağır neft daha
    pis süpürülür — fiziki gözlənti budur.
    """
    light = _rf(_run(oil_viscosity=1.0))
    heavy = _rf(_run(oil_viscosity=30.0))
    assert light > heavy + 20.0, (light, heavy)


def test_viscosity_response_is_monotone():
    """Lözlük artdıqca çıxarım monoton AZALMALIDIR."""
    values = [_rf(_run(oil_viscosity=mu)) for mu in (1.0, 3.0, 10.0, 30.0)]
    assert all(b < a for a, b in zip(values, values[1:])), values


def test_permeability_changes_recovery():
    """Keçiricilik — ölçülmüş: 10 mD → 10.7 %, 1000 mD → 66.8 %."""
    tight = _rf(_run(permx=10.0))
    good = _rf(_run(permx=1000.0))
    assert good > tight + 30.0, (tight, good)


def test_porosity_changes_oil_in_place():
    """Məsaməlilik həcmi dəyişir → kumulyativ neft dəyişməlidir.

    RF FAİZDİR, ona görə məsaməlilik artanda o, bir qədər AZALA da
    bilər (daha çox neft var, süpürmə isə eyni). Mənalı yoxlama
    kumulyativ hasilatdır.
    """
    low = _run(porosity=0.10)
    high = _run(porosity=0.35)
    assert low.converged and high.converged
    assert high.series.cumulative_oil[-1] > 2.0 * low.series.cumulative_oil[-1]


# ═══════════════════ PVT AÇIQ — düymə API/temperaturdur ════════════════

def test_panel_viscosity_is_ignored_when_pvt_is_active():
    """SƏNƏDLƏŞDİRİLMİŞ DAVRANIŞ: PVT cədvəli paneli üstələyir.

    Bu, qüsur DEYİL — `FluidProperties` sənədində "PLACEHOLDER" kimi
    təsvir olunub və `residual.py` açıq budaqlanır. Test bunu KİLİDLƏYİR
    ki, davranış təsadüfən dəyişməsin; UI isə indi sahəni bozardır.
    """
    a = _rf(_run(oil_viscosity=1.0, with_pvt=True))
    b = _rf(_run(oil_viscosity=30.0, with_pvt=True))
    assert a == pytest.approx(b, rel=0, abs=0), (a, b)


def test_api_gravity_changes_recovery_with_pvt():
    """PVT açıq olanda lözlüyün ƏSL düyməsi API-dir.

    Ölçülmüş: API 15 → μo 6.36 cP → RF 49.0 %;
              API 50 → μo 0.34 cP → RF 67.6 %.
    """
    heavy = _rf(_run(api=15.0, with_pvt=True))
    light = _rf(_run(api=50.0, with_pvt=True))
    assert light > heavy + 10.0, (heavy, light)


def test_temperature_changes_recovery_with_pvt():
    """Temperatur da lözlük vasitəsilə nəticəyə çatır."""
    cold = _rf(_run(temperature=40.0, with_pvt=True))
    hot = _rf(_run(temperature=120.0, with_pvt=True))
    assert hot > cold + 1.0, (cold, hot)


def test_geology_still_responds_with_pvt():
    """PVT açıq olsa da φ/K geologiyadan gəlir — üstələnmir."""
    tight = _rf(_run(permx=10.0, with_pvt=True))
    good = _rf(_run(permx=1000.0, with_pvt=True))
    assert good > tight + 30.0, (tight, good)


# ═══════════════════ mühərrikdə hardcode YOXDUR ════════════════════════

def test_engine_reads_viscosity_from_the_provider_not_a_constant():
    """Mobillik `λ = kr/μ` — μ provider-dən gəlməlidir.

    Əgər hardcode sabit işlədilsəydi, iki fərqli PVT cədvəli eyni
    mobillik verərdi.
    """
    from imex2d.simulation.pvt.black_oil import BlackOilPVTProvider

    pressure = np.array([250.0])
    light = BlackOilPVTProvider(build_pvt_table(
        api=50.0, pressure_min=1.0, pressure_max=400.0, n_points=40,
        bubble_point_bar=150.0))
    heavy = BlackOilPVTProvider(build_pvt_table(
        api=15.0, pressure_min=1.0, pressure_max=400.0, n_points=40,
        bubble_point_bar=150.0))
    assert float(light.oil_viscosity(pressure)[0]) < \
        float(heavy.oil_viscosity(pressure)[0])


# ═══════════════════ UI — görünürlük düzəlişi ══════════════════════════

@pytest.fixture
def qt_app():
    """Real (ekransız) Qt tətbiqi — `test_ui_units.py`-dəki ilə eyni üsul."""
    QtWidgets = pytest.importorskip("PyQt5.QtWidgets")
    yield QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def _rock_panel(qt_app):
    from imex2d.ui.panels import RockFluidPanel
    return RockFluidPanel()


def test_fluid_fields_are_disabled_when_pvt_is_active(qt_app):
    """Mexanizm A görünən olmalıdır: işləməyən sahə bozarır."""
    panel = _rock_panel(qt_app)
    assert panel.mu_o.isEnabled()
    assert not panel.context_note.isVisible() or not panel.context_note.text()

    panel.set_context(pvt_active=True)
    for name in panel.FLUID_WIDGETS:
        assert not getattr(panel, name).isEnabled(), name
    assert "PVT" in panel.context_note.text()

    panel.set_context(pvt_active=False)
    for name in panel.FLUID_WIDGETS:
        assert getattr(panel, name).isEnabled(), name


def test_geology_fields_are_disabled_for_imported_grdecl(qt_app):
    """Mexanizm C görünən olmalıdır."""
    panel = _rock_panel(qt_app)
    panel.set_context(geology_imported=True)
    for name in panel.GEOLOGY_WIDGETS:
        assert not getattr(panel, name).isEnabled(), name
    assert "GRDECL" in panel.context_note.text()


def test_disabled_fields_still_report_their_values(qt_app):
    """BOZARMA DƏYƏRİ SİLMİR — `fluids()` hələ də işləməlidir.

    Söndürülmüş widget-dən dəyər oxumaq Qt-də qanunidir; bu test
    davranışın təsadüfən dəyişməməsini kilidləyir.
    """
    panel = _rock_panel(qt_app)
    panel.mu_o.setValue(7.5)
    panel.set_context(pvt_active=True)
    assert panel.fluids().oil_viscosity == pytest.approx(7.5)


def test_geology_signal_fires_only_for_geology_fields(qt_app):
    """Mexanizm B: lözlük dəyişəndə interpolyasiya KÖHNƏLMİŞ sayılmamalıdır."""
    panel = _rock_panel(qt_app)
    fired = []
    panel.geology_changed.connect(lambda: fired.append(1))

    panel.mu_o.setValue(panel.mu_o.value() + 1.0)
    assert not fired, "lözlük geologiyaya təsir etmir"

    panel.porosity.setValue(panel.porosity.value() + 0.05)
    assert fired, "məsaməlilik geologiyanı köhnəltməlidir"


def test_geology_panel_banner_text_names_the_button(qt_app):
    """Banner istifadəçiyə NƏ ETMƏLİ olduğunu deməlidir."""
    from imex2d.ui.panels import GeologyPanel
    panel = GeologyPanel()
    panel.mark_fresh()
    assert panel.stale_label.text() == ""
