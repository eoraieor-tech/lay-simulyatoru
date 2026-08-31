"""Giriş panelləri.

Hər panelin yeganə işi: istifadəçi girişini DOMAIN obyektinə çevirmək.
Panellər hesablama aparmır, yoxlama etmir və simulyatordan xəbərsizdir.
Yoxlama ReservoirModel.validate() və SimulationConfig.validate()
metodlarındadır.
"""

from __future__ import annotations
from typing import List, Optional

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QBrush, QColor
from PyQt5.QtWidgets import (QCheckBox, QComboBox, QDialog, QDoubleSpinBox,
                             QFileDialog, QFormLayout, QHBoxLayout,
                             QHeaderView, QLabel, QLineEdit, QMessageBox,
                             QPushButton, QSpinBox, QTableWidget,
                             QTableWidgetItem, QTextEdit, QVBoxLayout, QWidget)

from ..application.config import (LinearSolverConfig, OutputConfig,
                                  SimulationConfig, TimeSteppingConfig)
import os

from ..application.scenarios import WELL_PATTERNS
from ..geology.interpolation import (INTERPOLATORS, InverseDistance,
                                     NearestNeighbour, OrdinaryKriging)
from ..domain.geology import GeologicalWell, validate_wells
from ..domain.geometry import depth_to_k, xy_to_ij
from ..domain.structure import FaultReference
from ..io.fault_io import (FaultFormatError, read_eclipse_faults,
                          read_faults_csv)
from ..io.scal_io import ScalFormatError, read_scal_csv, read_swof
from ..domain.grid import CartesianGrid
from ..domain.initial import InitialConditions
from ..domain.properties import FluidProperties
from ..domain.pvt import PVTTable
from ..domain.scal import (CapillaryParameters, CoreyParameters)
from ..domain.wells import ControlMode, Well, WellControl, WellType, Perforation
from ..rendering.theme import PALETTE
from .geology_map import GeologyMapWidget


def _spin(value, lo, hi, decimals=2, step=1.0, suffix=""):
    box = QDoubleSpinBox()
    box.setRange(lo, hi)
    box.setDecimals(decimals)
    box.setSingleStep(step)
    box.setValue(value)
    if suffix:
        box.setSuffix(f"  {suffix}")
    return box


def _ispin(value, lo, hi):
    box = QSpinBox()
    box.setRange(lo, hi)
    box.setValue(value)
    return box


class GridGeometryPanel(QWidget):
    changed = pyqtSignal()

    def __init__(self):
        super().__init__()
        form = QFormLayout(self)
        self.nx = _ispin(41, 3, 300)
        self.ny = _ispin(41, 1, 300)
        self.dx = _spin(20, 0.5, 1000, 1, 5, "m")
        self.dy = _spin(20, 0.5, 1000, 1, 5, "m")
        self.dz = _spin(10, 0.1, 500, 1, 1, "m")
        self.nz = _ispin(1, 1, 60)
        self.top_depth = _spin(2000.0, 0.0, 8000.0, 1, 50.0, "m")
        self.base_depth = _spin(2010.0, 0.1, 9000.0, 1, 50.0, "m")
        self.thickness_mode = QComboBox()
        self.thickness_mode.addItem("Təbəqə qalınlığı (DZ) ilə", "DZ")
        self.thickness_mode.addItem("Baza dərinliyi ilə", "BASE")
        self.thickness_mode.currentIndexChanged.connect(self._on_mode_changed)
        self.thickness_mode.currentIndexChanged.connect(self.changed)
        self.per_layer = QCheckBox("Hər təbəqə üçün ayrı qalınlıq")
        self.per_layer.toggled.connect(self._on_per_layer_toggled)
        self.per_layer.toggled.connect(self.changed)
        self.dz_table = QTableWidget(0, 1)
        self.dz_table.setHorizontalHeaderLabels(["Qalınlıq (m)"])
        self.dz_table.horizontalHeader().setStretchLastSection(True)
        self.dz_table.verticalHeader().setDefaultSectionSize(24)
        self.dz_table.setMaximumHeight(160)
        self.dz_table.setVisible(False)
        self.dip_x = _spin(0.0, -20.0, 20.0, 2, 0.5, "m/hüc")
        self.dip_y = _spin(0.0, -20.0, 20.0, 2, 0.5, "m/hüc")

        top_rows = [("NX", self.nx), ("NY", self.ny), ("DX", self.dx),
                    ("DY", self.dy), ("NZ (təbəqə sayı)", self.nz),
                    ("Qalınlıq necə verilir", self.thickness_mode),
                    ("Təbəqə qalınlığı DZ", self.dz)]
        bottom_rows = [("Tavan dərinliyi", self.top_depth),
                       ("Baza dərinliyi", self.base_depth),
                       ("Maillik X üzrə", self.dip_x),
                       ("Maillik Y üzrə", self.dip_y)]
        for label, widget in top_rows:
            form.addRow(label, widget)
            signal = getattr(widget, "valueChanged", None)
            if signal is not None:
                signal.connect(self.changed)
        form.addRow(self.per_layer)
        form.addRow("Təbəqə qalınlıqları", self.dz_table)
        for label, widget in bottom_rows:
            form.addRow(label, widget)
            signal = getattr(widget, "valueChanged", None)
            if signal is not None:
                signal.connect(self.changed)

        self._sync_table_rows(self.nz.value())

        self.info = QLabel()
        self.info.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        form.addRow(self.info)
        self.changed.connect(self._refresh_info)
        # `_on_mode_changed` məlumat sətrini yeniləyir, ona görə YALNIZ
        # `self.info` yaradıldıqdan sonra çağırıla bilər.
        self._on_mode_changed()

    def _on_mode_changed(self):
        """Baza dərinliyi ilə DZ eyni kəmiyyəti təyin edir.

        Üç kəmiyyətdən (tavan, qalınlıq, baza) yalnız ikisi müstəqildir.
        Hər üçü sərbəst olsaydı, ziddiyyət yaranardı. Ona görə istifadəçi
        seçir: qalınlığı birbaşa verir, yoxsa baza dərinliyindən
        hesablatdırır — geoloji işdə karotajdan məhz tavan və daban
        oxunur, qalınlıq isə onlardan çıxır. Baza dərinliyi rejimi tək
        orta qalınlıq təyin etdiyi üçün hər-təbəqə cədvəli ilə birgə
        işlədilmir.
        """
        by_base = self.thickness_mode.currentData() == "BASE"
        self.per_layer.setEnabled(not by_base)
        if by_base and self.per_layer.isChecked():
            self.per_layer.setChecked(False)
        self.dz.setEnabled(not by_base and not self.per_layer.isChecked())
        self.base_depth.setEnabled(by_base)
        self._refresh_info()

    def _on_per_layer_toggled(self, checked: bool) -> None:
        if checked:
            self._sync_table_rows(self.nz.value())
        by_base = self.thickness_mode.currentData() == "BASE"
        self.dz.setEnabled(not checked and not by_base)
        self.dz_table.setVisible(checked)

    def _sync_table_rows(self, n: int) -> None:
        """Cədvəlin sətir sayını NZ-yə uyğunlaşdırır, mövcud dəyərləri saxlayır."""
        current = self.dz_table.rowCount()
        if n == current:
            return
        if n > current:
            default = (self.dz_table.cellWidget(current - 1, 0).value()
                      if current > 0 else self.dz.value())
            self.dz_table.setRowCount(n)
            for row in range(current, n):
                spin = _spin(default, 0.1, 500, 2, 1, "m")
                spin.valueChanged.connect(self.changed)
                self.dz_table.setCellWidget(row, 0, spin)
                self.dz_table.setVerticalHeaderItem(
                    row, QTableWidgetItem(f"Təbəqə {row + 1}"))
        else:
            self.dz_table.setRowCount(n)

    def layer_thickness(self) -> float:
        """Bir təbəqənin qalınlığı — seçilmiş rejimdən asılı olaraq (uniform)."""
        if self.thickness_mode.currentData() != "BASE":
            return self.dz.value()
        span = self.base_depth.value() - self.top_depth.value()
        if span <= 0.0:
            return self.dz.minimum()
        return max(span / max(self.nz.value(), 1), self.dz.minimum())

    def _per_layer_active(self) -> bool:
        return (self.per_layer.isChecked()
                and self.thickness_mode.currentData() != "BASE")

    def layer_thicknesses(self) -> list:
        """Hər təbəqənin qalınlığı, uzunluq NZ — uniform yoxsa cədvəldən."""
        n = self.nz.value()
        if self._per_layer_active():
            self._sync_table_rows(n)
            return [self.dz_table.cellWidget(row, 0).value() for row in range(n)]
        return [self.layer_thickness()] * n

    def set_layer_thicknesses(self, dz) -> None:
        """Fayldan/modeldən gələn qalınlığı panelə yazır (skalyar və ya massiv)."""
        try:
            values = [float(v) for v in dz]
        except TypeError:
            values = [float(dz)]
        if not values:
            return
        self.nz.setValue(len(values))
        uniform = all(abs(v - values[0]) < 1e-9 for v in values)
        if uniform:
            self.dz.setValue(values[0])
            if self.per_layer.isChecked():
                self.per_layer.setChecked(False)
        else:
            self._sync_table_rows(len(values))
            for row, value in enumerate(values):
                self.dz_table.cellWidget(row, 0).setValue(value)
            if not self.per_layer.isChecked():
                self.per_layer.setChecked(True)
        self._refresh_info()

    def _refresh_info(self):
        # Qoruyucu: siqnal panel tam qurulmamış da gələ bilər.
        if not hasattr(self, "info"):
            return
        lx = self.nx.value() * self.dx.value()
        ly = self.ny.value() * self.dy.value()
        n = self.nx.value() * self.ny.value() * self.nz.value()
        per_layer_active = self._per_layer_active()

        if per_layer_active:
            thicknesses = self.layer_thicknesses()
            thickness = sum(thicknesses)
            if len(set(thicknesses)) <= 1:
                dz_text = f"{thicknesses[0]:.2f}"
            else:
                dz_text = f"{min(thicknesses):.2f}–{max(thicknesses):.2f}"
        else:
            dz = self.layer_thickness()
            thickness = self.nz.value() * dz
            dz_text = f"{dz:.2f}"

        top = self.top_depth.value()
        dip = ((self.nx.value() - 1) * self.dip_x.value()
               + (self.ny.value() - 1) * self.dip_y.value())
        base = top + max(dip, 0.0) + thickness

        if not per_layer_active:
            if self.thickness_mode.currentData() == "BASE":
                self.dz.blockSignals(True)
                self.dz.setValue(self.layer_thickness())
                self.dz.blockSignals(False)
            else:
                self.base_depth.blockSignals(True)
                self.base_depth.setValue(top + thickness)
                self.base_depth.blockSignals(False)

        self.info.setText(
            f"Hüceyrə sayı: {n}     Sahə: {lx:.0f} × {ly:.0f} m     "
            f"DZ: {dz_text} m     Ümumi qalınlıq: {thickness:.1f} m\n"
            f"Dərinlik: {top:.0f} – {base:.0f} m     "
            f"Həcm: {lx * ly * thickness / 1e6:.2f} mln m³")

    def values(self) -> dict:
        thicknesses = self.layer_thicknesses()
        dz_value = thicknesses if self._per_layer_active() else thicknesses[0]
        return dict(nx=self.nx.value(), ny=self.ny.value(), dx=self.dx.value(),
                    dy=self.dy.value(), dz=dz_value,
                    nz=self.nz.value(), top_depth=self.top_depth.value(),
                    dip_x=self.dip_x.value(), dip_y=self.dip_y.value())

    def depth_range(self) -> tuple:
        """Layın dərinlik intervalı — OWC seçimində istifadəçiyə göstərilir."""
        top = self.top_depth.value()
        dipped = top + ((self.nx.value() - 1) * self.dip_x.value()
                        + (self.ny.value() - 1) * self.dip_y.value())
        thickness = sum(self.layer_thicknesses())
        return min(top, dipped), max(top, dipped) + thickness

    def grid(self) -> CartesianGrid:
        return CartesianGrid(self.nx.value(), self.ny.value(), self.nz.value())


class GeologyPanel(QWidget):
    """Quyu cədvəli (2 ·) → interpolyasiya parametrləri.

    CSV yükləməsinin əvəzidir: istifadəçi quyuları birbaşa cədvəldə
    redaktə edir. Panel heç bir hesablama aparmır — yalnız
    `list[GeologicalWell]` istehsal/qəbul edir; interpolyasiyanı
    `İnterpolyasiya et` düyməsi ilə application qatı işə salır
    (`MainWindow._interpolate_geology`).

    Cədvəl dəyişəndə interpolyasiya AVTOMATİK işə düşmür (böyük gridə
    yavaşdır) — yalnız `changed` siqnalı ilə "nəticə köhnəlib" bildirilir.
    """

    changed = pyqtSignal()
    interpolate_requested = pyqtSignal()

    COLUMNS = ["Ad", "Modeldə", "X, m", "Y, m", "(i, j)", "Lay üstü, m",
              "Lay altı, m", "φ", "k, mD", "Sw", "Qeyd"]
    COL_NAME = 0
    COL_IN_MODEL = 1
    COL_X = 2
    COL_Y = 3
    COL_IJ = 4
    COL_TOP = 5
    COL_BOTTOM = 6
    COL_PORO = 7
    COL_PERM = 8
    COL_SW = 9
    COL_NOTE = 10
    _NUMERIC_COLUMNS = {COL_X: "x", COL_Y: "y", COL_TOP: "top",
                        COL_BOTTOM: "bottom", COL_PORO: "porosity",
                        COL_PERM: "permeability", COL_SW: "water_saturation"}

    def __init__(self):
        super().__init__()
        self._geometry = None          # CellGeometry, grid qurulanda gəlir
        self._stale = False
        self._well_counter = 0
        layout = QVBoxLayout(self)

        toolbar = QHBoxLayout()
        self.add_button = QPushButton("Quyu əlavə et")
        self.duplicate_button = QPushButton("Dublikat")
        self.delete_button = QPushButton("Sil")
        self.centre_button = QPushButton("Grid mərkəzinə at")
        for button in (self.add_button, self.duplicate_button,
                      self.delete_button, self.centre_button):
            toolbar.addWidget(button)
        self.add_button.clicked.connect(self.add_row)
        self.duplicate_button.clicked.connect(self._duplicate_selected)
        self.delete_button.clicked.connect(self._delete_selected)
        self.centre_button.clicked.connect(self._centre_selected)
        layout.addLayout(toolbar)

        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.verticalHeader().setVisible(False)
        self.table.setMinimumHeight(160)
        self.table.itemChanged.connect(self._on_item_changed)
        layout.addWidget(self.table)

        self.map_widget = GeologyMapWidget()
        layout.addWidget(self.map_widget)

        self.validation_view = QTextEdit()
        self.validation_view.setReadOnly(True)
        self.validation_view.setMaximumHeight(90)
        self.validation_view.setStyleSheet("font-family:monospace;font-size:11px")
        layout.addWidget(self.validation_view)

        form = QFormLayout()
        self.method = QComboBox()
        self.method.addItems(list(INTERPOLATORS.keys()))
        self.method.setCurrentText("Kriging (adi)")
        self.method.currentIndexChanged.connect(self._on_method_changed)
        self.method.currentIndexChanged.connect(self._on_table_edited)
        form.addRow("Üsul", self.method)

        self.power = _spin(2.0, 0.5, 8.0, 2, 0.5)
        self.search_radius = _spin(0.0, 0.0, 100000.0, 0, 50.0, "m")
        self.range_ = _spin(0.0, 0.0, 100000.0, 0, 50.0, "m")
        self.sill = _spin(0.0, 0.0, 1e6, 4, 0.01)
        self.nugget = _spin(0.0, 0.0, 1e6, 4, 0.01)
        rows = [("IDW dərəcəsi p", self.power),
                ("Axtarış radiusu (0 = limitsiz)", self.search_radius),
                ("Kriging radiusu a (0 = avto)", self.range_),
                ("Sill c (0 = avto)", self.sill),
                ("Nugget c₀", self.nugget)]
        for label, widget in rows:
            form.addRow(label, widget)
            widget.valueChanged.connect(self._on_table_edited)
        layout.addLayout(form)

        action_row = QHBoxLayout()
        self.interpolate_button = QPushButton("İnterpolyasiya et")
        self.interpolate_button.clicked.connect(self.interpolate_requested)
        action_row.addWidget(self.interpolate_button)
        self.stale_label = QLabel("")
        self.stale_label.setStyleSheet("color:#e0a020;font-size:11px")
        action_row.addWidget(self.stale_label, 1)
        layout.addLayout(action_row)

        self.report = QTextEdit()
        self.report.setReadOnly(True)
        self.report.setMinimumHeight(110)
        self.report.setStyleSheet("font-family:monospace;font-size:11px")
        layout.addWidget(self.report, 1)

        note = QLabel("Boş xana = məlumat yoxdur (sıfır DEYİL). Cədvəl boşdursa "
                      "sintetik model işlədilir.")
        note.setWordWrap(True)
        note.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        layout.addWidget(note)
        self._on_method_changed()
        self._refresh_map_and_validation()

    # ------------------------------------------------------------ slots
    def _on_method_changed(self):
        method = self.method.currentText()
        self.power.setEnabled("IDW" in method)
        self.search_radius.setEnabled("IDW" in method)
        for widget in (self.range_, self.sill, self.nugget):
            widget.setEnabled("Kriging" in method)

    def _on_item_changed(self, item: QTableWidgetItem):
        if item.column() == self.COL_IN_MODEL:
            pass   # checkbox dəyişikliyi də buradan gəlir, əlavə iş lazım deyil
        self._on_table_edited()

    def _on_table_edited(self):
        self._recompute_indices()
        self._refresh_map_and_validation()
        self.mark_stale()
        self.changed.emit()

    def _selected_row(self) -> Optional[int]:
        row = self.table.currentRow()
        return row if row >= 0 else None

    def _duplicate_selected(self):
        row = self._selected_row()
        if row is None:
            return
        well = self._well_from_row(row)
        if well is None:
            return
        well.name = self._unique_name(well.name + "-kopya")
        self.add_row(well)

    def _delete_selected(self):
        row = self._selected_row()
        if row is not None:
            self.table.removeRow(row)
            self._on_table_edited()

    def _centre_selected(self):
        row = self._selected_row()
        if row is None or self._geometry is None:
            return
        x_max, y_max = self._geometry.areal_extent()
        self.table.blockSignals(True)
        self.table.setItem(row, self.COL_X, QTableWidgetItem(f"{x_max / 2.0:g}"))
        self.table.setItem(row, self.COL_Y, QTableWidgetItem(f"{y_max / 2.0:g}"))
        self.table.blockSignals(False)
        self._on_table_edited()

    # ----------------------------------------------------------- public
    def set_geometry(self, geometry) -> None:
        """Grid həndəsəsi (varsa) — (i, j) sütunu və sərhəd yoxlaması üçün."""
        self._geometry = geometry
        self._recompute_indices()
        self._refresh_map_and_validation()

    def method_text(self) -> str:
        return self.method.currentText()

    def interpolator(self):
        method = self.method.currentText()
        if "IDW" in method:
            radius = self.search_radius.value()
            return InverseDistance(power=self.power.value(),
                                   search_radius=radius if radius > 0 else None)
        if "Kriging" in method:
            return OrdinaryKriging(
                range_=self.range_.value() or None,
                sill=self.sill.value() or None,
                nugget=self.nugget.value())
        return NearestNeighbour()

    def set_report(self, text: str):
        self.report.setPlainText(text)

    def mark_stale(self):
        self._stale = True
        self.stale_label.setText("Nəticə köhnəlib — 'İnterpolyasiya et' basın."
                                 if self.wells() else "")

    def mark_fresh(self):
        self._stale = False
        self.stale_label.setText("")

    @property
    def is_stale(self) -> bool:
        return self._stale

    def set_validation(self, issues) -> None:
        if not issues:
            self.validation_view.setPlainText("Xəta/xəbərdarlıq yoxdur.")
            return
        prefixes = {"error": "[XƏTA] ", "warning": "[XƏBƏRDARLIQ] ", "info": "[MƏLUMAT] "}
        lines = [prefixes.get(issue.level, "") + issue.message for issue in issues]
        self.validation_view.setPlainText("\n".join(lines))

    def has_blocking_errors(self) -> bool:
        issues = validate_wells(self.wells(), self._geometry, self.method_text())
        return any(issue.level == "error" for issue in issues)

    def add_row(self, well: Optional[GeologicalWell] = None):
        if well is None:
            self._well_counter += 1
            well = GeologicalWell(name=self._unique_name(f"W-{self._well_counter}"),
                                  in_model=True, x=0.0, y=0.0)
        r = self.table.rowCount()
        self.table.blockSignals(True)
        self.table.insertRow(r)
        self.table.setItem(r, self.COL_NAME, QTableWidgetItem(well.name))

        check_item = QTableWidgetItem()
        check_item.setFlags(check_item.flags() | Qt.ItemIsUserCheckable)
        check_item.setCheckState(Qt.Checked if well.in_model else Qt.Unchecked)
        self.table.setItem(r, self.COL_IN_MODEL, check_item)

        self.table.setItem(r, self.COL_X, QTableWidgetItem(f"{well.x:g}"))
        self.table.setItem(r, self.COL_Y, QTableWidgetItem(f"{well.y:g}"))

        ij_item = QTableWidgetItem("—")
        ij_item.setFlags(ij_item.flags() & ~Qt.ItemIsEditable)
        ij_item.setForeground(QBrush(QColor(Qt.gray)))
        self.table.setItem(r, self.COL_IJ, ij_item)

        optional_columns = [(self.COL_TOP, "top"), (self.COL_BOTTOM, "bottom"),
                           (self.COL_PORO, "porosity"), (self.COL_PERM, "permeability"),
                           (self.COL_SW, "water_saturation")]
        for column, attr in optional_columns:
            value = getattr(well, attr)
            self.table.setItem(r, column, QTableWidgetItem(
                "" if value is None else f"{value:g}"))
        self.table.setItem(r, self.COL_NOTE, QTableWidgetItem(well.note))
        self.table.blockSignals(False)
        self._on_table_edited()

    def load(self, wells: List[GeologicalWell]):
        self.table.blockSignals(True)
        self.table.setRowCount(0)
        self.table.blockSignals(False)
        for well in wells:
            self.add_row(well)
        if not wells:
            self._on_table_edited()

    def wells(self) -> List[GeologicalWell]:
        result = []
        for row in range(self.table.rowCount()):
            well = self._well_from_row(row)
            if well is not None:
                result.append(well)
        return result

    # -------------------------------------------------------- internal
    def _unique_name(self, base: str) -> str:
        existing = set()
        for row in range(self.table.rowCount()):
            item = self.table.item(row, self.COL_NAME)
            if item is not None:
                existing.add(item.text().strip())
        name, suffix = base, 1
        while name in existing:
            suffix += 1
            name = f"{base}-{suffix}"
        return name

    def _well_from_row(self, row: int) -> Optional[GeologicalWell]:
        name_item = self.table.item(row, self.COL_NAME)
        if name_item is None:
            return None
        check_item = self.table.item(row, self.COL_IN_MODEL)
        in_model = check_item is not None and check_item.checkState() == Qt.Checked
        values = {}
        for column, attr in self._NUMERIC_COLUMNS.items():
            item = self.table.item(row, column)
            text = (item.text().strip() if item is not None else "")
            values[attr] = self._to_float(text)
        note_item = self.table.item(row, self.COL_NOTE)
        return GeologicalWell(
            name=name_item.text().strip(),
            in_model=in_model,
            x=values["x"] or 0.0, y=values["y"] or 0.0,
            top=values["top"], bottom=values["bottom"],
            porosity=values["porosity"], permeability=values["permeability"],
            water_saturation=values["water_saturation"],
            note=note_item.text() if note_item is not None else "")

    @staticmethod
    def _to_float(text: str) -> Optional[float]:
        text = text.strip().replace(",", ".")
        if not text:
            return None
        try:
            return float(text)
        except ValueError:
            return None

    def _recompute_indices(self):
        self.table.blockSignals(True)
        try:
            for row in range(self.table.rowCount()):
                item = self.table.item(row, self.COL_IJ)
                if item is None:
                    continue
                if self._geometry is None:
                    item.setText("grid qurulduqdan sonra")
                    continue
                well = self._well_from_row(row)
                if well is None:
                    continue
                x_max, y_max = self._geometry.areal_extent()
                if not (0.0 <= well.x <= x_max and 0.0 <= well.y <= y_max):
                    item.setText("kənar")
                    continue
                i, j = xy_to_ij(well.x, well.y, self._geometry)
                item.setText(f"({i}, {j})")
        finally:
            self.table.blockSignals(False)

    def _refresh_map_and_validation(self):
        wells = self.wells()
        if self._geometry is not None:
            x_max, y_max = self._geometry.areal_extent()
        else:
            x_max = max((w.x for w in wells), default=1.0) or 1.0
            y_max = max((w.y for w in wells), default=1.0) or 1.0
        selected = None
        row = self._selected_row()
        if row is not None:
            item = self.table.item(row, self.COL_NAME)
            selected = item.text().strip() if item is not None else None
        self.map_widget.set_data(wells, x_max, y_max, selected)
        issues = validate_wells(wells, self._geometry, self.method_text())
        self.set_validation(issues)


class FaultPanel(QWidget):
    """Fault siyahısı: CSV, Eclipse FAULTS/MULTFLT, ya da əl ilə.

    Fault tam əl ilə həndəsə düzəltmək (I/J/K müstəvi, diapazon)
    sürgü ilə edilə bilməz — buna görə bu panel siyahını sadəcə
    GÖSTƏRİR; fault yaratmaq faylla (CSV/Eclipse) və ya "Fault əlavə
    et" pəncərəsi ilədir.
    """

    changed = pyqtSignal()

    def __init__(self):
        super().__init__()
        self.faults: List = []
        layout = QVBoxLayout(self)

        buttons = QHBoxLayout()
        self.load_csv = QPushButton("CSV yüklə…")
        self.load_eclipse = QPushButton("Eclipse FAULTS yüklə…")
        self.add_manual = QPushButton("Fault əlavə et…")
        self.clear_button = QPushButton("Təmizlə")
        self.load_csv.clicked.connect(self._load_csv)
        self.load_eclipse.clicked.connect(self._load_eclipse)
        self.add_manual.clicked.connect(self._add_manual)
        self.clear_button.clicked.connect(self._clear)
        for button in (self.load_csv, self.load_eclipse, self.add_manual,
                      self.clear_button):
            buttons.addWidget(button)
        layout.addLayout(buttons)

        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(
            ["Ad", "Ox", "Müstəvi", "A diapazonu", "B diapazonu", "Çarpan"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.verticalHeader().setVisible(False)
        self.table.setMaximumHeight(160)
        layout.addWidget(self.table)

        hint = QLabel("Müstəvi 0-based-dir: 'plane_index=10' i=10 ilə i=11 "
                      "arasındakı sərhəddir. Boş diapazon bütün grid-i "
                      "əhatə edir. Çarpan 0 = tam bağlı fay.")
        hint.setWordWrap(True)
        hint.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        layout.addWidget(hint)

    def _refresh_table(self):
        self.table.setRowCount(0)
        for fault in self.faults:
            row = self.table.rowCount()
            self.table.insertRow(row)
            span = lambda b: "hamısı" if b is None else f"{b[0]}-{b[1]}"
            cells = [fault.name, fault.axis or "—",
                    str(fault.plane_index) if fault.plane_index is not None else "—",
                    span(fault.range_a), span(fault.range_b),
                    "BAĞLI" if fault.sealing else f"{fault.transmissibility_multiplier:g}"]
            for column, text in enumerate(cells):
                self.table.setItem(row, column, QTableWidgetItem(text))

    def _load(self, reader, title, filter_text):
        path, _ = QFileDialog.getOpenFileName(self, title, "", filter_text)
        if not path:
            return
        try:
            loaded = reader(path)
        except FaultFormatError as error:
            QMessageBox.warning(self, "Fault faylı oxunmadı", str(error))
            return
        except Exception as error:
            QMessageBox.critical(self, "Fault faylı oxunmadı",
                                 f"Gözlənilməz xəta: {error}")
            return
        self.faults = loaded
        self._refresh_table()
        self.changed.emit()

    def _load_csv(self):
        self._load(read_faults_csv, "Fault CSV",
                  "CSV (*.csv *.txt);;Bütün fayllar (*)")

    def _load_eclipse(self):
        self._load(read_eclipse_faults, "Eclipse deck (FAULTS)",
                  "Eclipse (*.DATA *.data *.inc);;Bütün fayllar (*)")

    def _add_manual(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Fault əlavə et")
        form = QFormLayout(dialog)

        name = QLineEdit(f"F{len(self.faults) + 1}")
        axis = QComboBox(); axis.addItems(["I", "J", "K"])
        plane = QSpinBox(); plane.setRange(0, 9999)
        multiplier = QDoubleSpinBox(); multiplier.setRange(0.0, 1.0)
        multiplier.setDecimals(3); multiplier.setValue(0.1)
        sealing = QCheckBox("Tam bağlı (çarpanı görməzdən gəlir)")
        for label, widget in (("Ad", name), ("Ox", axis),
                              ("Müstəvi indeksi", plane),
                              ("Çarpan", multiplier), ("", sealing)):
            form.addRow(label, widget)

        buttons = QHBoxLayout()
        ok = QPushButton("Əlavə et"); cancel = QPushButton("Ləğv et")
        ok.clicked.connect(dialog.accept)
        cancel.clicked.connect(dialog.reject)
        buttons.addWidget(ok); buttons.addWidget(cancel)
        form.addRow(buttons)

        if dialog.exec_() != QDialog.Accepted:
            return
        try:
            fault = FaultReference(
                name=name.text().strip() or f"F{len(self.faults) + 1}",
                source_id=name.text().strip(), axis=axis.currentText(),
                plane_index=plane.value(),
                transmissibility_multiplier=multiplier.value(),
                sealing=sealing.isChecked())
        except ValueError as error:
            QMessageBox.warning(self, "Yanlış dəyər", str(error))
            return
        self.faults.append(fault)
        self._refresh_table()
        self.changed.emit()

    def _clear(self):
        self.faults = []
        self._refresh_table()
        self.changed.emit()

    def values(self) -> Optional[List]:
        return list(self.faults) if self.faults else None


class ScalSourcePanel(QWidget):
    """SCAL mənbəyi: Corey düsturu, yoxsa laboratoriya cədvəli.

    Corey sadə modellər üçün kifayətdir. Real kern məlumatı olanda
    cədvəl işlədilməlidir — əyrilər asimmetrik olur və düsturla
    ifadə edilmir.
    """

    changed = pyqtSignal()

    def __init__(self):
        super().__init__()
        self.tables = None
        layout = QVBoxLayout(self)

        self.mode = QComboBox()
        self.mode.addItem("Corey düsturu (aşağıdakı parametrlər)", "COREY")
        self.mode.addItem("Laboratoriya cədvəli", "TABLE")
        self.mode.currentIndexChanged.connect(self._on_mode_changed)
        self.mode.currentIndexChanged.connect(self.changed)
        layout.addWidget(self.mode)

        buttons = QHBoxLayout()
        self.load_csv = QPushButton("CSV yüklə…")
        self.load_swof = QPushButton("Eclipse SWOF yüklə…")
        self.load_csv.clicked.connect(self._load_csv)
        self.load_swof.clicked.connect(self._load_swof)
        buttons.addWidget(self.load_csv)
        buttons.addWidget(self.load_swof)
        layout.addLayout(buttons)

        self.info = QLabel("Cədvəl yüklənməyib.")
        self.info.setWordWrap(True)
        self.info.setStyleSheet(
            f"background:{PALETTE.panel_alt};border:1px solid {PALETTE.line};"
            f"border-radius:3px;padding:6px;font-family:monospace;"
            f"font-size:11px;color:{PALETTE.text}")
        layout.addWidget(self.info)

        hint = QLabel("CSV sütunları: region, sw, krw, kro [, pc]. "
                      "'region' olmasa hamısı bir zonaya düşür. "
                      "Region nömrələri GRDECL-dəki SATNUM ilə uyğun gəlməlidir.")
        hint.setWordWrap(True)
        hint.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        layout.addWidget(hint)
        self._on_mode_changed()

    def _on_mode_changed(self):
        by_table = self.mode.currentData() == "TABLE"
        self.load_csv.setEnabled(by_table)
        self.load_swof.setEnabled(by_table)

    def _load(self, reader, title, filter_text):
        path, _ = QFileDialog.getOpenFileName(self, title, "", filter_text)
        if not path:
            return
        try:
            self.tables = reader(path)
        except ScalFormatError as error:
            QMessageBox.warning(self, "Cədvəl oxunmadı", str(error))
            return
        except Exception as error:
            QMessageBox.critical(self, "Cədvəl oxunmadı",
                                 f"Gözlənilməz xəta: {error}")
            return
        self.info.setText(f"{os.path.basename(path)}\n"
                          f"{len(self.tables)} region\n"
                          + self.tables.summary())
        self.mode.setCurrentIndex(1)
        self.changed.emit()

    def _load_csv(self):
        self._load(read_scal_csv, "SCAL cədvəli",
                   "CSV (*.csv *.txt);;Bütün fayllar (*)")

    def _load_swof(self):
        self._load(read_swof, "Eclipse deck (SWOF)",
                   "Eclipse (*.DATA *.data *.inc);;Bütün fayllar (*)")

    def is_enabled(self) -> bool:
        return self.mode.currentData() == "TABLE" and self.tables is not None


class RockFluidPanel(QWidget):
    changed = pyqtSignal()

    def __init__(self):
        super().__init__()
        form = QFormLayout(self)
        self.porosity = _spin(0.22, 0.01, 0.45, 3, 0.01)
        self.permx = _spin(150, 0.01, 20000, 1, 10, "mD")
        self.ky_over_kx = _spin(1.0, 0.01, 10, 2, 0.1)
        self.kv_over_kh = _spin(0.10, 0.001, 1.0, 3, 0.01)
        self.heterogeneity = QComboBox()
        self.heterogeneity.addItems(["Homogen", "Təsadüfi (log-normal)"])
        self.sigma = _spin(0.5, 0.05, 2.0, 2, 0.05)
        self.seed = _ispin(7, 0, 9999)
        self.mu_w = _spin(0.5, 0.05, 50, 2, 0.05, "cP")
        self.mu_o = _spin(3.0, 0.05, 5000, 2, 0.5, "cP")
        self.bo = _spin(1.15, 1.0, 3.0, 3, 0.01)
        self.rock_compressibility = _spin(4.5e-5, 1e-6, 1e-3, 7, 1e-5, "1/bar")
        rows = [("Məsaməlilik φ", self.porosity), ("Keçiricilik Kx", self.permx),
                ("Ky/Kx", self.ky_over_kx), ("Kv/Kh (şaquli)", self.kv_over_kh),
                ("Heterogenlik", self.heterogeneity),
                ("σ (log-normal)", self.sigma), ("Seed", self.seed),
                ("Su lözlüyü μw", self.mu_w), ("Neft lözlüyü μo", self.mu_o),
                ("Bo", self.bo), ("Süxur sıxılması", self.rock_compressibility)]
        for label, widget in rows:
            form.addRow(label, widget)
            sig = getattr(widget, "valueChanged", None) or widget.currentIndexChanged
            sig.connect(self.changed)

    def geology_values(self) -> dict:
        return dict(porosity=self.porosity.value(), permx_base=self.permx.value(),
                    ky_over_kx=self.ky_over_kx.value(),
                    kv_over_kh=self.kv_over_kh.value(),
                    heterogeneous=self.heterogeneity.currentIndex() == 1,
                    sigma=self.sigma.value(), seed=self.seed.value())

    def fluids(self) -> FluidProperties:
        return FluidProperties(water_viscosity=self.mu_w.value(),
                               oil_viscosity=self.mu_o.value(),
                               oil_fvf=self.bo.value())

    def rock_compressibility_value(self) -> float:
        return self.rock_compressibility.value()


class ScalPanel(QWidget):
    changed = pyqtSignal()

    def __init__(self):
        super().__init__()
        form = QFormLayout(self)
        self.swc = _spin(0.20, 0.0, 0.6, 3, 0.01)
        self.sor = _spin(0.25, 0.0, 0.6, 3, 0.01)
        self.krw_end = _spin(0.35, 0.01, 1.0, 3, 0.05)
        self.kro_end = _spin(0.90, 0.01, 1.0, 3, 0.05)
        self.nw = _spin(2.5, 1.0, 6.0, 2, 0.1)
        self.no = _spin(2.0, 1.0, 6.0, 2, 0.1)
        self.pc_entry = _spin(0.0, 0.0, 20.0, 3, 0.05, "bar")
        self.pc_lambda = _spin(2.0, 0.2, 10.0, 2, 0.1)
        self.pc_max = _spin(5.0, 0.1, 100.0, 2, 1.0, "bar")
        for label, widget in [("Swc (bağlı su)", self.swc), ("Sor (qalıq neft)", self.sor),
                              ("krw @ 1-Sor", self.krw_end), ("kro @ Swc", self.kro_end),
                              ("Corey nw", self.nw), ("Corey no", self.no),
                              ("Pc giriş təzyiqi Pe", self.pc_entry),
                              ("Brooks-Corey λ", self.pc_lambda),
                              ("Pc yuxarı həddi", self.pc_max)]:
            form.addRow(label, widget)
            widget.valueChanged.connect(self.changed)
        note = QLabel("Pe = 0 → kapilyar təzyiq söndürülür (köhnə davranış).")
        note.setWordWrap(True)
        note.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        form.addRow(note)

    def values(self) -> CoreyParameters:
        return CoreyParameters(self.swc.value(), self.sor.value(),
                               self.krw_end.value(), self.kro_end.value(),
                               self.nw.value(), self.no.value())

    def capillary_values(self) -> CapillaryParameters:
        return CapillaryParameters(entry_pressure=self.pc_entry.value(),
                                   lambda_exponent=self.pc_lambda.value(),
                                   max_pressure=self.pc_max.value())


class PvtPanel(QWidget):
    """PVT cədvəli — korrelyasiyalardan qurulur.

    Panel yalnız PVTTable istehsal edir; provider-i application qatı yaradır.
    Söndürüləndə None qaytarır və mühərrik statik dəyərlərlə işləyir.
    """

    changed = pyqtSignal()

    def __init__(self):
        super().__init__()
        form = QFormLayout(self)
        self.enabled = QCheckBox("PVT modelini işlət (təzyiqdən asılı flüid)")
        self.enabled.setChecked(False)
        form.addRow(self.enabled)

        self.api = _spin(32.0, 5.0, 60.0, 1, 1.0, "°API")
        self.gas_gravity = _spin(0.75, 0.55, 1.5, 3, 0.01)
        self.temperature = _spin(70.0, 5.0, 250.0, 1, 5.0, "°C")
        self.salinity = _spin(30000.0, 0.0, 300000.0, 0, 5000.0, "ppm")
        self.bubble_point = _spin(240.0, 5.0, 1000.0, 1, 10.0, "bar")
        # Defolt 1 bar (əvvəl 10 idi) — CƏDVƏLDƏN KƏNAR
        # EKSTRAPOLYASİYA ciddi problem yaradır: cədvəldən aşağıda
        # bütün xassələr SABİTLƏŞİR (törəmə = 0), yuxarıda isə yox —
        # bu, SINIQ nöqtədir və Nyuton onun ətrafında osilyasiya edir.
        # Ölçülüb: qaz fazasının yığılmama probleminin kökü məhz bu
        # idi — problemli hüceyrənin təzyiqi 10.01 bar, yəni köhnə
        # sərhədin dəqiq üstündə qalmışdı.
        self.pressure_min = _spin(1.0, 0.5, 500.0, 1, 10.0, "bar")
        self.pressure_max = _spin(400.0, 20.0, 1200.0, 1, 20.0, "bar")
        self.points = _ispin(40, 5, 200)

        rows = [("Neftin sıxlığı (API)", self.api),
                ("Qaz sıxlığı γg", self.gas_gravity),
                ("Lay temperaturu", self.temperature),
                ("Su duzluluğu", self.salinity),
                ("Doyma təzyiqi Pb", self.bubble_point),
                ("Cədvəl: min təzyiq", self.pressure_min),
                ("Cədvəl: maks təzyiq", self.pressure_max),
                ("Nöqtə sayı", self.points)]
        for label, widget in rows:
            form.addRow(label, widget)
            widget.valueChanged.connect(self.changed)
        self.enabled.stateChanged.connect(self.changed)

        note = QLabel("Söndürülübsə, sabit Bo və μ dəyərləri işlədilir "
                      "(2-ci paneldəki qiymətlər).")
        note.setWordWrap(True)
        note.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        form.addRow(note)

    def is_enabled(self) -> bool:
        return self.enabled.isChecked()

    def values(self):
        """PVTTable və ya None."""
        if not self.enabled.isChecked():
            return None
        from ..simulation.pvt.correlations import build_pvt_table
        return build_pvt_table(
            api=self.api.value(),
            gas_gravity=self.gas_gravity.value(),
            temperature_c=self.temperature.value(),
            salinity_ppm=self.salinity.value(),
            pressure_min=self.pressure_min.value(),
            pressure_max=max(self.pressure_max.value(),
                             self.pressure_min.value() + 10.0),
            n_points=self.points.value(),
            bubble_point_bar=self.bubble_point.value())


class WellPanel(QWidget):
    """Quyu rejimi (7 ·) — geologiya cədvəlinə (2 ·) bağlıdır.

    `in_model = True` olan hər geologiya quyusu buraya avtomatik sətir kimi
    düşür (`set_geology_context`). İstifadəçi burada YALNIZ rejimi (Tip,
    İdarə, Qiymət, rw) və perforasiya intervalını (METRLƏ) təyin edir —
    `Ad`/`i`/`j`/`k` geologiya cədvəlindən və qrid həndəsəsindən HESABLANIR,
    redaktə olunmur.

    Bayraq (`in_model`) söndürüləndə sətir cədvəldən yoxa çıxır, AMMA rejim
    məlumatı `_retained`-də saxlanılır — istifadəçi yenidən işarələsə,
    əvvəlki BHP/rate geri qayıdır (bax `_sync_rows`).
    """

    changed = pyqtSignal()

    COLUMNS = ["Ad", "i", "j", "Perf üst, m", "Perf alt, m", "k",
              "Tip", "İdarə", "Qiymət", "rw"]
    COL_NAME = 0
    COL_I = 1
    COL_J = 2
    COL_PERF_TOP = 3
    COL_PERF_BOTTOM = 4
    COL_K = 5
    COL_TYPE = 6
    COL_MODE = 7
    COL_TARGET = 8
    COL_RW = 9

    def __init__(self):
        super().__init__()
        self._geology_by_name: dict = {}
        self._geometry = None
        self._retained: dict = {}
        self._last_in_model_names: set = set()
        layout = QVBoxLayout(self)
        row = QHBoxLayout()
        self.pattern = QComboBox()
        self.pattern.addItems(list(WELL_PATTERNS.keys()))
        self.apply_button = QPushButton("Tətbiq et")
        row.addWidget(self.pattern, 1)
        row.addWidget(self.apply_button)
        layout.addLayout(row)

        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.verticalHeader().setVisible(False)
        self.table.setMinimumHeight(210)
        self.table.itemChanged.connect(self._on_item_changed)
        layout.addWidget(self.table)

        self.warning_label = QLabel("")
        self.warning_label.setWordWrap(True)
        self.warning_label.setStyleSheet("color:#e0a020;font-size:11px")
        layout.addWidget(self.warning_label)

        hint = QLabel("BHP → bar,   RATE → m³/gün (rezervuar həcmi)\n"
                      "Ad/i/j/k geologiya cədvəlindən avtomatik gəlir. "
                      "Perf üst/alt boşdursa bütün lay perforasiya olunur.")
        hint.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        layout.addWidget(hint)

    # ------------------------------------------------------------ slots
    def _on_item_changed(self, item: QTableWidgetItem):
        if item.column() in (self.COL_PERF_TOP, self.COL_PERF_BOTTOM):
            self._recompute_ij_k()
        self.changed.emit()

    # ----------------------------------------------------------- public
    def set_geology_context(self, wells: List[GeologicalWell], geometry) -> None:
        """Geologiya cədvəli və ya grid dəyişəndə çağırılır."""
        self._geology_by_name = {w.name: w for w in wells}
        self._geometry = geometry
        self._sync_rows(wells)
        self._recompute_ij_k()
        self.changed.emit()

    def load(self, wells: List[Well]):
        """Fayldan bərpa: rejimi (və perf metrlərini) birbaşa yazır."""
        self.table.blockSignals(True)
        self.table.setRowCount(0)
        try:
            for well in wells:
                self._add_row(well.name, {
                    "kind": well.well_type.value, "mode": well.control.mode.value,
                    "target": well.control.target, "rw": well.radius,
                    "perf_top": well.perf_top, "perf_bottom": well.perf_bottom})
        finally:
            self.table.blockSignals(False)
        self._recompute_ij_k()

    def values(self) -> List[Well]:
        wells: List[Well] = []
        for row in range(self.table.rowCount()):
            name_item = self.table.item(row, self.COL_NAME)
            if name_item is None:
                continue
            name = name_item.text().strip()
            perf_top = self._to_float(self.table.item(row, self.COL_PERF_TOP).text())
            perf_bottom = self._to_float(self.table.item(row, self.COL_PERF_BOTTOM).text())
            try:
                kind = self.table.cellWidget(row, self.COL_TYPE).currentText()
                mode = self.table.cellWidget(row, self.COL_MODE).currentText()
            except AttributeError:
                continue
            target = self._to_float(self.table.item(row, self.COL_TARGET).text()) or 0.0
            rw = self._to_float(self.table.item(row, self.COL_RW).text()) or 0.1

            i, j, first, last = self._resolve_ijk(name, perf_top, perf_bottom)
            wells.append(Well(
                name=name, well_type=WellType(kind),
                control=WellControl(ControlMode(mode), target),
                perforations=[Perforation(i, j, k) for k in range(first, last + 1)],
                radius=rw, perf_top=perf_top, perf_bottom=perf_bottom))
        return wells

    # -------------------------------------------------------- internal
    @staticmethod
    def _to_float(text: str) -> Optional[float]:
        text = (text or "").strip().replace(",", ".")
        if not text:
            return None
        try:
            return float(text)
        except ValueError:
            return None

    @staticmethod
    def _default_state() -> dict:
        return {"kind": "PROD", "mode": "BHP", "target": 150.0, "rw": 0.1,
                "perf_top": None, "perf_bottom": None}

    def _find_row(self, name: str) -> Optional[int]:
        for row in range(self.table.rowCount()):
            item = self.table.item(row, self.COL_NAME)
            if item is not None and item.text().strip() == name:
                return row
        return None

    def _row_state(self, row: int) -> dict:
        return {
            "kind": self.table.cellWidget(row, self.COL_TYPE).currentText(),
            "mode": self.table.cellWidget(row, self.COL_MODE).currentText(),
            "target": self._to_float(self.table.item(row, self.COL_TARGET).text()),
            "rw": self._to_float(self.table.item(row, self.COL_RW).text()),
            "perf_top": self._to_float(self.table.item(row, self.COL_PERF_TOP).text()),
            "perf_bottom": self._to_float(self.table.item(row, self.COL_PERF_BOTTOM).text()),
        }

    def _add_row(self, name: str, state: dict):
        r = self.table.rowCount()
        self.table.insertRow(r)
        name_item = QTableWidgetItem(name)
        name_item.setFlags(name_item.flags() & ~Qt.ItemIsEditable)
        self.table.setItem(r, self.COL_NAME, name_item)

        for column in (self.COL_I, self.COL_J, self.COL_K):
            item = QTableWidgetItem("—")
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            item.setForeground(QBrush(QColor(Qt.gray)))
            self.table.setItem(r, column, item)

        for column, key in ((self.COL_PERF_TOP, "perf_top"),
                            (self.COL_PERF_BOTTOM, "perf_bottom")):
            value = state.get(key)
            self.table.setItem(r, column, QTableWidgetItem(
                "" if value is None else f"{value:g}"))

        type_box = QComboBox()
        type_box.addItems(["PROD", "INJ"])
        type_box.setCurrentText(state.get("kind", "PROD"))
        type_box.currentIndexChanged.connect(lambda *_: self.changed.emit())
        self.table.setCellWidget(r, self.COL_TYPE, type_box)

        mode_box = QComboBox()
        mode_box.addItems(["BHP", "RATE"])
        mode_box.setCurrentText(state.get("mode", "BHP"))
        mode_box.currentIndexChanged.connect(lambda *_: self.changed.emit())
        self.table.setCellWidget(r, self.COL_MODE, mode_box)

        self.table.setItem(r, self.COL_TARGET,
                           QTableWidgetItem(f"{state.get('target') or 150.0:g}"))
        self.table.setItem(r, self.COL_RW,
                           QTableWidgetItem(f"{state.get('rw') or 0.1:g}"))

    def _sync_rows(self, wells: List[GeologicalWell]):
        in_model_names = {w.name for w in wells if w.in_model}
        newly_off = self._last_in_model_names - in_model_names
        newly_on = in_model_names - self._last_in_model_names
        self.table.blockSignals(True)
        try:
            for name in newly_off:
                row = self._find_row(name)
                if row is not None:
                    self._retained[name] = self._row_state(row)
                    self.table.removeRow(row)
            for name in newly_on:
                if self._find_row(name) is None:
                    state = self._retained.pop(name, None) or self._default_state()
                    self._add_row(name, state)
        finally:
            self.table.blockSignals(False)
        self._last_in_model_names = in_model_names

    def _resolve_ijk(self, name: str, perf_top: Optional[float],
                     perf_bottom: Optional[float]):
        geo = self._geology_by_name.get(name)
        if geo is None or self._geometry is None:
            return 0, 0, 0, 0
        i, j = xy_to_ij(geo.x, geo.y, self._geometry)
        nz = self._geometry.grid.nz
        k_top = (depth_to_k(geo.x, geo.y, perf_top, self._geometry)
                if perf_top is not None else 0)
        k_bottom = (depth_to_k(geo.x, geo.y, perf_bottom, self._geometry)
                   if perf_bottom is not None else nz - 1)
        k_top = 0 if k_top is None else k_top
        k_bottom = nz - 1 if k_bottom is None else k_bottom
        first, last = (k_top, k_bottom) if k_top <= k_bottom else (k_bottom, k_top)
        return i, j, first, last

    def _recompute_ij_k(self):
        self.table.blockSignals(True)
        warnings = []
        try:
            for row in range(self.table.rowCount()):
                name_item = self.table.item(row, self.COL_NAME)
                if name_item is None:
                    continue
                name = name_item.text().strip()
                i_item = self.table.item(row, self.COL_I)
                j_item = self.table.item(row, self.COL_J)
                k_item = self.table.item(row, self.COL_K)
                if self._geometry is None:
                    i_item.setText("—")
                    j_item.setText("—")
                    k_item.setText("grid qurulduqdan sonra")
                    continue
                geo = self._geology_by_name.get(name)
                if geo is None:
                    i_item.setText("?")
                    j_item.setText("?")
                    k_item.setText("geologiyada yoxdur")
                    continue
                perf_top = self._to_float(self.table.item(row, self.COL_PERF_TOP).text())
                perf_bottom = self._to_float(self.table.item(row, self.COL_PERF_BOTTOM).text())
                i, j, first, last = self._resolve_ijk(name, perf_top, perf_bottom)
                i_item.setText(str(i))
                j_item.setText(str(j))
                k_item.setText(f"{first + 1}–{last + 1}")
                if (perf_top is not None
                        and depth_to_k(geo.x, geo.y, perf_top, self._geometry) is None):
                    warnings.append(f"'{name}': perforasiya üstü lay qalınlığından kənardadır.")
                if (perf_bottom is not None
                        and depth_to_k(geo.x, geo.y, perf_bottom, self._geometry) is None):
                    warnings.append(f"'{name}': perforasiya altı lay qalınlığından kənardadır.")
        finally:
            self.table.blockSignals(False)
        self.warning_label.setText("\n".join(warnings))


class NumericalPanel(QWidget):
    changed = pyqtSignal()

    def __init__(self):
        super().__init__()
        form = QFormLayout(self)
        self.initial_pressure = _spin(250, 1, 1200, 1, 10, "bar")
        self.initial_sw = _spin(0.20, 0.0, 1.0, 3, 0.01)
        self.end_time = _spin(1500, 1, 40000, 0, 100, "gün")
        self.max_dt = _spin(20, 0.01, 365, 2, 1, "gün")
        self.cfl = _spin(0.45, 0.05, 0.95, 2, 0.05)
        self.snapshots = _ispin(60, 5, 400)
        self.engine = QComboBox()
        self.engine.addItem("IMPES (explicit doyumluluq)", "IMPES")
        self.engine.addItem("Fully implicit (Nyuton)", "IMPLICIT")
        self.engine.currentIndexChanged.connect(self.changed)
        form.addRow("Hesablama sxemi", self.engine)
        self.use_equilibration = QCheckBox("Equilibration (dərinlikdən asılı ilkin şərtlər)")
        self.datum_depth = _spin(2000.0, 0.0, 8000.0, 1, 50.0, "m")
        self.owc = _spin(2050.0, 0.0, 8000.0, 1, 10.0, "m")
        form.addRow(self.use_equilibration)
        for label, widget in [("Başlanğıc təzyiq", self.initial_pressure),
                              ("Başlanğıc Sw", self.initial_sw),
                              ("Datum dərinliyi", self.datum_depth),
                              ("Su-neft kontaktı (OWC)", self.owc),
                              ("Simulyasiya müddəti", self.end_time),
                              ("Maks. Δt", self.max_dt),
                              ("CFL əmsalı", self.cfl),
                              ("Yaddaş anlarının sayı", self.snapshots)]:
            form.addRow(label, widget)
            widget.valueChanged.connect(self.changed)
        self.use_equilibration.stateChanged.connect(self.changed)
        note = QLabel("Söndürülübsə, bütün hüceyrələrdə eyni təzyiq və Sw "
                      "işlədilir (köhnə davranış).")
        note.setWordWrap(True)
        note.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        form.addRow(note)

    def initial_conditions(self) -> InitialConditions:
        equilibrate = self.use_equilibration.isChecked()
        return InitialConditions(
            datum_depth=self.datum_depth.value(),
            datum_pressure=self.initial_pressure.value(),
            water_saturation=self.initial_sw.value(),
            oil_water_contact=self.owc.value() if equilibrate else None,
            use_equilibration=equilibrate)

    def engine_choice(self) -> str:
        return self.engine.currentData()

    def simulation_config(self) -> SimulationConfig:
        return SimulationConfig(
            end_time=self.end_time.value(),
            time_stepping=TimeSteppingConfig(max_dt=self.max_dt.value(),
                                             cfl_factor=self.cfl.value()),
            linear_solver=LinearSolverConfig(),
            output=OutputConfig(snapshot_count=self.snapshots.value()),
        )
