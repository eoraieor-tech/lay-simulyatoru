"""Qaz doyumluluğu (Sg) görüntüdə — Model və 3D tabının "Xassə" siyahısı.

BİLDİRİŞ (sahibkar, `ISH_HESABATI.md` → Seans 20). Üç fazalı mühərrik
işləyir, GOR hesablanır, amma "Xassə" siyahısında Sg YOX idi — qazın
cəbhəsini ekranda izləmək mümkün deyildi. Snapshot-lar Sg-ni onsuz da
saxlayırdı (`three_phase_engine._record_snapshot`); çatışmayan yalnız
görüntü qatı idi.

GİZLİ TƏLƏ. Hər iki seçicinin (`_select_volume`, `_select`) sonunda
tanınmayan açar üçün SƏSSİZ defolt var — məsaməlilik qaytarılır. Yəni
siyahıya sadəcə "Sg" əlavə etmək ekranda Sg adı altında φ göstərərdi.
`test_selector_does_not_fall_through_to_porosity` bunu kilidləyir.
"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from helpers import default_scal
from imex2d.application.config import SimulationConfig
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.scenarios import (SyntheticGeologicalModelBuilder,
                                          five_spot)
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.scal import GasCoreyParameters
from imex2d.rendering import renderers as R
from imex2d.simulation.implicit.engine import FullyImplicitEngine
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.pvt.correlations import build_pvt_table
from imex2d.simulation.results import Snapshot
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter


def _model(with_gas=True, nx=8, ny=8, bubble_point=240.0):
    geology = SyntheticGeologicalModelBuilder().build(
        nx=nx, ny=ny, dx=20.0, dy=20.0, dz=10.0,
        porosity=0.22, permx_base=150.0)
    return ReservoirModelBuilder().build(
        geological_model=geology, wells=five_spot(geology.grid),
        scal=default_scal(),
        gas_scal=GasCoreyParameters() if with_gas else None,
        pvt_table=build_pvt_table(pressure_min=1.0, pressure_max=400.0,
                                  n_points=40, bubble_point_bar=bubble_point,
                                  include_gas=with_gas),
        name="Sg görüntü sınağı")


def _snapshot(model, gas=None):
    n = model.grid.ncell
    return Snapshot(time=10.0, pressure=np.full(n, 200.0),
                    water_saturation=np.full(n, 0.3), gas_saturation=gas)


# ═══════════════════════ etiket və açar ══════════════════════════════

def test_gas_saturation_has_a_label():
    assert R.PROPERTY_LABELS[R.GAS_SATURATION] == "Qaz doyumluluğu (Sg)"
    assert R.property_label(R.GAS_SATURATION) == "Qaz doyumluluğu (Sg)"


def test_key_does_not_collide_with_existing_keys():
    others = {R.SATURATION, R.PRESSURE, R.PERMEABILITY, R.POROSITY, R.DEPTH}
    assert R.GAS_SATURATION not in others


# ═══════════════════════ seçicilər — dəyərlər ════════════════════════

def test_three_phase_values_are_passed_through_unchanged():
    model = _model()
    gas = np.linspace(0.0, 0.2, model.grid.ncell)
    values, cmap, low, high = R.MapRenderer()._select_volume(
        model, R.GAS_SATURATION, _snapshot(model, gas))
    assert values.shape == model.grid.shape
    assert np.allclose(values.ravel(), gas)
    assert cmap is R.GAS_SATURATION_CMAP


def test_two_phase_snapshot_shows_zero_gas():
    """İki fazalı modeldə Sg ≡ 0 — fiziki olaraq dəqiq, uydurma deyil."""
    model = _model(with_gas=False)
    values, *_ = R.MapRenderer()._select_volume(
        model, R.GAS_SATURATION, _snapshot(model, gas=None))
    assert np.all(values == 0.0)


def test_no_snapshot_gives_nan_not_a_made_up_zero():
    """Simulyasiyadan əvvəl ilkin qaz papağı ola bilər — 0 UYDURULMUR."""
    model = _model()
    values, *_ = R.MapRenderer()._select_volume(model, R.GAS_SATURATION, None)
    assert np.all(np.isnan(values))


def test_selector_does_not_fall_through_to_porosity():
    """GİZLİ TƏLƏ: son `return` tanınmayan açarı φ kimi göstərir."""
    model = _model()
    gas = np.full(model.grid.ncell, 0.123)
    values, cmap, *_ = R.MapRenderer()._select_volume(
        model, R.GAS_SATURATION, _snapshot(model, gas))
    porosity = model.rock.porosity.values
    assert not np.allclose(values.ravel(), porosity)
    assert cmap is not R.POROSITY_CMAP


def test_layer_selector_matches_volume_selector():
    """Model tabının təbəqə kəsiyi 3D massivin eyni təbəqəsi olmalıdır."""
    model = _model()
    gas = np.random.default_rng(3).random(model.grid.ncell) * 0.2
    snap = _snapshot(model, gas)
    volume, *_ = R.MapRenderer()._select_volume(model, R.GAS_SATURATION, snap)
    layer, cmap, low, high = R.MapRenderer()._select(
        model, R.GAS_SATURATION, snap, layer=0)
    assert np.allclose(layer, volume[0])
    assert cmap is R.GAS_SATURATION_CMAP


# ═══════════════════════ rəng şkalası ════════════════════════════════

def test_colour_scale_is_fixed_across_frames():
    """Animasiya boyu şkala DƏYİŞMƏMƏLİDİR — yoxsa cəbhə gözə çarpmaz."""
    model = _model()
    early = _snapshot(model, np.full(model.grid.ncell, 0.001))
    late = _snapshot(model, np.full(model.grid.ncell, 0.15))
    _, _, low1, high1 = R.MapRenderer()._select_volume(model, R.GAS_SATURATION, early)
    _, _, low2, high2 = R.MapRenderer()._select_volume(model, R.GAS_SATURATION, late)
    assert (low1, high1) == (low2, high2)
    assert low1 == 0.0
    assert high1 == pytest.approx(1.0 - model.scal_parameters.swc)


# ═══════════════════════ UI siyahıları ═══════════════════════════════

def test_both_property_dropdowns_list_gas_saturation():
    """Model və 3D tabının siyahısı — pəncərə QURULMADAN yoxlanılır.

    `MainWindow()` qurmaq testdə `pytest`-i çökdürür (bax
    `test_playback_controls.py`), ona görə mənbə mətni yoxlanılır.
    """
    pytest.importorskip("PyQt5.QtWidgets")
    from imex2d.ui import main_window as MW
    source = inspect.getsource(MW.MainWindow)
    listing = "for key in (R.SATURATION, R.GAS_SATURATION, R.PRESSURE"
    assert source.count(listing) == 2, "Model və 3D siyahısının İKİSİ də"


# ═══════════════════════ real üç fazalı qaçış ════════════════════════

def test_real_three_phase_run_shows_a_gas_front():
    """Uc-uca: Pb = 240 bar-da qaz ayrılır → görüntüdə Sg > 0 hüceyrələr var."""
    model = _model(bubble_point=240.0)
    service = ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(default_scal()),
        linear_solver=ScipyCgIluSolver(), engine_factory=FullyImplicitEngine)
    result = service.run(model, SimulationConfig(end_time=400.0))
    assert result.converged, result.message

    first, *_ = R.MapRenderer()._select_volume(
        model, R.GAS_SATURATION, result.snapshots[0])
    last, *_ = R.MapRenderer()._select_volume(
        model, R.GAS_SATURATION, result.snapshots[-1])
    assert float(np.nanmax(last)) > 1e-3, "son kadrda sərbəst qaz görünməlidir"
    assert np.count_nonzero(last > 1e-3) >= np.count_nonzero(first > 1e-3), \
        "qaz zonası zamanla kiçilməməlidir"
