"""Seans 46 — panellərin vəziyyəti saxlanılıb bərpa olunur (`ui/panel_state.py`).

Hər paneldə BÜTÜN sadə sahələr defoltdan fərqli dəyərə çəkilir, vəziyyət
götürülür və TƏZƏ panelə bərpa olunur. İki şey yoxlanılır:

  1. bərpadan sonra hər sahə eynidir (heç biri defolta qayıtmır);
  2. panelin qurduğu domain obyektləri (konfiqurasiya, SCAL, flüidlər,
     ilkin şərtlər) eynidir — model məhz bunlardan qurulur.
"""

from __future__ import annotations

import dataclasses
import json

import numpy as np
import pytest
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QLineEdit,
                             QSpinBox)

from imex2d.ui import panels as P
from imex2d.ui.panel_state import _widgets, capture, restore

PANELS = (P.GridGeometryPanel, P.GeologyPanel, P.RockFluidPanel, P.ScalPanel,
          P.ScalSourcePanel, P.PvtPanel, P.WellPanel, P.NumericalPanel)


@pytest.fixture(scope="module")
def qapp():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield app


def _scramble(panel) -> None:
    """Hər sahəni defoltdan fərqli, lakin icazəli dəyərə çəkir."""
    widgets = list(_widgets(panel))
    for _, widget in widgets:
        if isinstance(widget, QComboBox) and widget.count() > 1:
            widget.setCurrentIndex(widget.count() - 1)
    for _, widget in widgets:
        if isinstance(widget, QCheckBox):
            widget.setChecked(not widget.isChecked())
    for _, widget in widgets:
        if isinstance(widget, QDoubleSpinBox):
            low = widget.minimum()
            high = min(widget.maximum(), low + 1000.0)
            widget.setValue(round(low + 0.37 * (high - low), widget.decimals()))
        elif isinstance(widget, QSpinBox):
            low = widget.minimum()
            high = min(widget.maximum(), low + 50)
            widget.setValue(low + int(0.37 * (high - low)) + 1)
        elif isinstance(widget, QLineEdit):
            widget.setText(widget.text() + "7")


def _restored(cls):
    original = cls()
    _scramble(original)
    copy = cls()
    restore(copy, json.loads(json.dumps(capture(original))))   # fayldan keçmiş kimi
    return original, copy


@pytest.mark.parametrize("cls", PANELS, ids=lambda cls: cls.__name__)
def test_every_field_survives_capture_and_restore(qapp, cls):
    original, copy = _restored(cls)
    state = capture(original)
    assert state, "panelin heç bir sahəsi tapılmadı"
    changed = [name for name, value in state.items()
               if value != capture(cls()).get(name)]
    assert len(changed) >= len(state) - 1, "sınaq sahələri dəyişdirə bilmədi"
    assert capture(copy) == state


def test_restored_panels_build_the_same_domain_objects(qapp):
    rock, rock_copy = _restored(P.RockFluidPanel)
    assert rock_copy.fluids() == rock.fluids()
    assert rock_copy.geology_values() == rock.geology_values()

    scal, scal_copy = _restored(P.ScalPanel)
    assert scal_copy.values() == scal.values()
    assert scal_copy.gas_values() == scal.gas_values()
    assert scal_copy.capillary_values() == scal.capillary_values()

    numerical, numerical_copy = _restored(P.NumericalPanel)
    assert numerical_copy.simulation_config() == numerical.simulation_config()
    assert numerical_copy.initial_conditions() == numerical.initial_conditions()
    assert numerical_copy.engine_choice() == numerical.engine_choice()

    grid, grid_copy = _restored(P.GridGeometryPanel)
    assert grid_copy.values() == grid.values()
    assert grid_copy.structure_values() == grid.structure_values()


def test_restored_pvt_panel_builds_the_same_table(qapp):
    pvt, pvt_copy = _restored(P.PvtPanel)
    first, second = pvt.values(), pvt_copy.values()
    assert (first is None) == (second is None)
    assert first is not None, "sınaq PVT-ni aktivləşdirməli idi"
    for field in dataclasses.fields(first):
        a, b = getattr(first, field.name), getattr(second, field.name)
        if isinstance(a, np.ndarray):
            np.testing.assert_array_equal(a, b, err_msg=field.name)
        else:
            assert a == b, field.name


def test_unit_selector_restored_before_value(qapp):
    """psi seçilib 3000 yazılıb: bərpa bar-da 3000 YOX, psi-də 3000 verməlidir."""
    panel = P.NumericalPanel()
    panel.initial_pressure_unit.setCurrentText("psi")
    panel.initial_pressure.setValue(3000.0)
    copy = P.NumericalPanel()
    restore(copy, capture(panel))
    assert copy.initial_pressure_unit.currentText() == "psi"
    assert copy.initial_pressure.value() == pytest.approx(3000.0)
    assert (copy.initial_conditions().datum_pressure
            == pytest.approx(panel.initial_conditions().datum_pressure))


def test_unknown_and_missing_keys_are_ignored(qapp):
    panel = P.ScalPanel()
    before = capture(panel)
    restore(panel, {"yoxdur": 5, "swc": "səhv tip", "nw": True})
    assert capture(panel) == before
    restore(panel, {})
    assert capture(panel) == before
