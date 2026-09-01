"""Phase 1 (giriş boru xətti) — panel vahid seçicilərinin FUNKSIONAL testi.

Digər UI testlərindən (`test_ui_static.py`/`test_ui_wiring.py`) fərqli
olaraq bu, real (ekransız `offscreen` platforma ilə) Qt widget-ləri
qurur — çünki burada YALNIZ struktur deyil, faktiki ƏDƏDİ çevirmə
yoxlanılır (bax `conftest.py`, `QT_QPA_PLATFORM=offscreen`).
"""

from __future__ import annotations

import pytest

QtWidgets = pytest.importorskip("PyQt5.QtWidgets")


@pytest.fixture(scope="module")
def qapp():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield app


def test_rock_fluid_panel_default_units_reproduce_previous_numeric_behaviour(qapp):
    """Toxunulmayan panel (defolt vahidlər) ƏVVƏLKİ kimi ədədi olaraq
    dəyişməz — yəni `geology_values()['permx_base']` == vidjetin özü."""
    from imex2d.ui.panels import RockFluidPanel
    panel = RockFluidPanel()
    assert panel.permx_unit.currentText() == "mD"
    assert panel.viscosity_unit.currentText() == "cP"
    panel.permx.setValue(150.0)
    panel.mu_w.setValue(0.5)
    panel.mu_o.setValue(3.0)

    values = panel.geology_values()
    fluids = panel.fluids()
    assert values["permx_base"] == pytest.approx(150.0)
    assert fluids.water_viscosity == pytest.approx(0.5)
    assert fluids.oil_viscosity == pytest.approx(3.0)


def test_rock_fluid_panel_converts_darcy_permeability_to_engine_md(qapp):
    """Realist UI axını: ƏVVƏLCƏ vahid seçilir, SONRA dəyər yazılır —
    bax `_bind_unit_aware_spins`-in niyə vahid dəyişəndə mövcud dəyəri
    YENİDƏN HESABLADIĞI (əks halda diapazon köhnə vahiddə qalıb səssizcə
    kəsərdi, bax modulun öz şərhi)."""
    from imex2d.domain import unit_conversions as uc
    from imex2d.ui.panels import RockFluidPanel
    panel = RockFluidPanel()
    panel.permx_unit.setCurrentText("D")
    panel.permx.setValue(1.5)          # 1.5 D

    values = panel.geology_values()
    assert values["permx_base"] == pytest.approx(uc.convert(1.5, "D", "mD", "permeability"))
    assert values["permx_base"] == pytest.approx(1500.0)   # 1 D = 1000 mD


def test_rock_fluid_panel_converts_pas_viscosity_to_engine_cp(qapp):
    from imex2d.domain import unit_conversions as uc
    from imex2d.ui.panels import RockFluidPanel
    panel = RockFluidPanel()
    panel.viscosity_unit.setCurrentText("Pa.s")
    panel.mu_o.setValue(0.003)           # 0.003 Pa.s

    fluids = panel.fluids()
    assert fluids.oil_viscosity == pytest.approx(uc.pas_to_cp(0.003))
    assert fluids.oil_viscosity == pytest.approx(3.0)


def test_numerical_panel_default_unit_reproduces_previous_numeric_behaviour(qapp):
    from imex2d.ui.panels import NumericalPanel
    panel = NumericalPanel()
    assert panel.initial_pressure_unit.currentText() == "bar"
    panel.initial_pressure.setValue(250.0)
    ic = panel.initial_conditions()
    assert ic.datum_pressure == pytest.approx(250.0)


def test_facies_panel_parses_proportions_and_builds_config(qapp):
    from imex2d.ui.panels import FaciesPanel
    panel = FaciesPanel()
    assert panel.column_name_value() == "FACIES"
    assert panel.parse_proportions() is None   # boş -> müşahidə olunana keçid

    panel.proportions_text.setText("0:0.6, 1:0.4")
    assert panel.parse_proportions() == {0: 0.6, 1: 0.4}

    panel.seed.setValue(7)
    config = panel.build_config(realization_id=2, seed_offset=1000)
    assert config.proportions == {0: 0.6, 1: 0.4}
    assert config.seed == 1007
    assert config.realization_id == 2


def test_facies_panel_rejects_malformed_proportions(qapp):
    from imex2d.ui.panels import FaciesPanel
    panel = FaciesPanel()
    panel.proportions_text.setText("not-a-valid-format")
    with pytest.raises(ValueError):
        panel.parse_proportions()


def test_numerical_panel_converts_psi_pressure_to_engine_bar(qapp):
    from imex2d.domain import unit_conversions as uc
    from imex2d.ui.panels import NumericalPanel
    panel = NumericalPanel()
    panel.initial_pressure_unit.setCurrentText("psi")
    panel.initial_pressure.setValue(3000.0)
    ic = panel.initial_conditions()
    assert ic.datum_pressure == pytest.approx(uc.psi_to_bar(3000.0), rel=1e-9)
    assert ic.datum_pressure == pytest.approx(206.8427, rel=1e-4)
