"""Giriş panelləri.

Hər panelin yeganə işi: istifadəçi girişini DOMAIN obyektinə çevirmək.
Panellər hesablama aparmır, yoxlama etmir və simulyatordan xəbərsizdir.
Yoxlama ReservoirModel.validate() və SimulationConfig.validate()
metodlarındadır.
"""

from __future__ import annotations
from typing import List, Optional

from PyQt5.QtCore import pyqtSignal
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
from ..geology.well_data_io import (WellDataFormatError, read_well_csv,
                                    write_example_csv)
from ..domain.structure import FaultReference
from ..io.fault_io import (FaultFormatError, read_eclipse_faults,
                          read_faults_csv)
from ..io.scal_io import ScalFormatError, read_scal_csv, read_swof
from ..domain.grid import CartesianGrid
from ..domain.initial import InitialConditions
from ..domain.properties import FluidProperties
from ..domain.pvt import PVTTable
from ..domain.scal import (CapillaryParameters, CoreyParameters,
                           GasCoreyParameters)
from ..domain.wells import ControlMode, Well, WellControl, WellType, Perforation
from ..rendering.theme import PALETTE


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
        self.dip_x = _spin(0.0, -20.0, 20.0, 2, 0.5, "m/hüc")
        self.dip_y = _spin(0.0, -20.0, 20.0, 2, 0.5, "m/hüc")
        rows = [("NX", self.nx), ("NY", self.ny), ("DX", self.dx),
                ("DY", self.dy), ("NZ (təbəqə sayı)", self.nz),
                ("Qalınlıq necə verilir", self.thickness_mode),
                ("Təbəqə qalınlığı DZ", self.dz),
                ("Tavan dərinliyi", self.top_depth),
                ("Baza dərinliyi", self.base_depth),
                ("Maillik X üzrə", self.dip_x),
                ("Maillik Y üzrə", self.dip_y)]
        for label, widget in rows:
            form.addRow(label, widget)
            signal = getattr(widget, "valueChanged", None)
            if signal is not None:
                signal.connect(self.changed)
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
        oxunur, qalınlıq isə onlardan çıxır.
        """
        by_base = self.thickness_mode.currentData() == "BASE"
        self.dz.setEnabled(not by_base)
        self.base_depth.setEnabled(by_base)
        self._refresh_info()

    def layer_thickness(self) -> float:
        """Bir təbəqənin qalınlığı — seçilmiş rejimdən asılı olaraq."""
        if self.thickness_mode.currentData() != "BASE":
            return self.dz.value()
        span = self.base_depth.value() - self.top_depth.value()
        if span <= 0.0:
            return self.dz.minimum()
        return max(span / max(self.nz.value(), 1), self.dz.minimum())

    def _refresh_info(self):
        # Qoruyucu: siqnal panel tam qurulmamış da gələ bilər.
        if not hasattr(self, "info"):
            return
        lx = self.nx.value() * self.dx.value()
        ly = self.ny.value() * self.dy.value()
        n = self.nx.value() * self.ny.value() * self.nz.value()
        dz = self.layer_thickness()
        thickness = self.nz.value() * dz

        top = self.top_depth.value()
        dip = ((self.nx.value() - 1) * self.dip_x.value()
               + (self.ny.value() - 1) * self.dip_y.value())
        base = top + max(dip, 0.0) + thickness

        if self.thickness_mode.currentData() == "BASE":
            self.dz.blockSignals(True)
            self.dz.setValue(dz)
            self.dz.blockSignals(False)
        else:
            self.base_depth.blockSignals(True)
            self.base_depth.setValue(top + thickness)
            self.base_depth.blockSignals(False)

        self.info.setText(
            f"Hüceyrə sayı: {n}     Sahə: {lx:.0f} × {ly:.0f} m     "
            f"DZ: {dz:.2f} m     Ümumi qalınlıq: {thickness:.1f} m\n"
            f"Dərinlik: {top:.0f} – {base:.0f} m     "
            f"Həcm: {lx * ly * thickness / 1e6:.2f} mln m³")

    def values(self) -> dict:
        return dict(nx=self.nx.value(), ny=self.ny.value(), dx=self.dx.value(),
                    dy=self.dy.value(), dz=self.layer_thickness(),
                    nz=self.nz.value(), top_depth=self.top_depth.value(),
                    dip_x=self.dip_x.value(), dip_y=self.dip_y.value())

    def depth_range(self) -> tuple:
        """Layın dərinlik intervalı — OWC seçimində istifadəçiyə göstərilir."""
        top = self.top_depth.value()
        dipped = top + ((self.nx.value() - 1) * self.dip_x.value()
                        + (self.ny.value() - 1) * self.dip_y.value())
        thickness = self.nz.value() * self.layer_thickness()
        return min(top, dipped), max(top, dipped) + thickness

    def grid(self) -> CartesianGrid:
        return CartesianGrid(self.nx.value(), self.ny.value(), self.nz.value())


class WellDataPanel(QWidget):
    """Quyu məlumatı (CSV) → interpolyasiya parametrləri.

    Panel yalnız məlumatı yükləyir və interpolyator qurur; grid xassələrini
    hesablamaq application qatının işidir (WellBasedGeologicalModelBuilder).
    """

    changed = pyqtSignal()

    def __init__(self):
        super().__init__()
        self.dataset = None
        layout = QVBoxLayout(self)

        self.enabled = QCheckBox("Quyu məlumatından geoloji model qur")
        self.enabled.setEnabled(False)
        self.enabled.stateChanged.connect(self.changed)
        layout.addWidget(self.enabled)

        buttons = QHBoxLayout()
        self.load_button = QPushButton("CSV yüklə…")
        self.example_button = QPushButton("Nümunə fayl yarat…")
        self.load_button.clicked.connect(self.load_csv)
        self.example_button.clicked.connect(self.create_example)
        buttons.addWidget(self.load_button)
        buttons.addWidget(self.example_button)
        layout.addLayout(buttons)

        self.summary = QLabel("Məlumat yüklənməyib.")
        self.summary.setWordWrap(True)
        self.summary.setStyleSheet(
            f"background:{PALETTE.panel_alt};border:1px solid {PALETTE.line};"
            f"border-radius:3px;padding:6px;font-size:11px;color:{PALETTE.text}")
        layout.addWidget(self.summary)

        form = QFormLayout()
        self.method = QComboBox()
        self.method.addItems(list(INTERPOLATORS.keys()))
        self.method.setCurrentText("Kriging (adi)")
        self.method.currentIndexChanged.connect(self._on_method_changed)
        self.method.currentIndexChanged.connect(self.changed)
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
            widget.valueChanged.connect(self.changed)
        layout.addLayout(form)

        self.report = QTextEdit()
        self.report.setReadOnly(True)
        self.report.setMinimumHeight(150)
        self.report.setStyleSheet("font-family:monospace;font-size:11px")
        layout.addWidget(self.report, 1)

        note = QLabel("Sütunlar: well, x, y [, k, depth], PORO, PERMX [, NTG …]. "
                      "Söndürülübsə sintetik model işlədilir.")
        note.setWordWrap(True)
        note.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        layout.addWidget(note)
        self._on_method_changed()

    # ------------------------------------------------------------ slots
    def _on_method_changed(self):
        method = self.method.currentText()
        self.power.setEnabled("IDW" in method)
        self.search_radius.setEnabled("IDW" in method)
        for widget in (self.range_, self.sill, self.nugget):
            widget.setEnabled("Kriging" in method)

    def load_csv(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Quyu məlumatı", "", "CSV (*.csv *.txt);;Bütün fayllar (*)")
        if not path:
            return
        try:
            self.dataset = read_well_csv(path)
        except WellDataFormatError as error:
            QMessageBox.warning(self, "Fayl oxunmadı", str(error))
            return
        except Exception as error:
            QMessageBox.critical(self, "Fayl oxunmadı", f"Gözlənilməz xəta: {error}")
            return
        info = self.dataset.summary()
        self.summary.setText(
            f"{os.path.basename(path)}\n"
            f"{info['quyu']} quyu · {info['nöqtə']} nöqtə · "
            f"{'təbəqəli' if info['təbəqəli'] else 'təbəqəsiz'}\n"
            f"Xassələr: {info['xassə']}\n{info['sahə']}")
        self.enabled.setEnabled(True)
        self.enabled.setChecked(True)
        self.changed.emit()

    def create_example(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Nümunə fayl", "quyular.csv", "CSV (*.csv)")
        if not path:
            return
        try:
            write_example_csv(path, nz=3)
        except Exception as error:
            QMessageBox.critical(self, "Yazılmadı", str(error))
            return
        QMessageBox.information(
            self, "Nümunə hazırdır",
            f"{os.path.basename(path)} yaradıldı.\n\n"
            "Onu 'CSV yüklə…' ilə açıb formatı görə bilərsən.")

    # ----------------------------------------------------------- public
    def is_enabled(self) -> bool:
        return self.enabled.isChecked() and self.dataset is not None

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

        # ── Qaz-neft SCAL (A7) ──────────────────────────────────────
        self.gas_enabled = QCheckBox("Qaz-neft əyriləri")
        self.gas_enabled.setChecked(False)
        self.gas_enabled.stateChanged.connect(self._on_gas_toggled)
        self.gas_enabled.stateChanged.connect(self.changed)
        form.addRow(self.gas_enabled)

        self.sgc = _spin(0.05, 0.0, 0.4, 3, 0.01)
        self.sorg = _spin(0.10, 0.0, 0.4, 3, 0.01)
        self.krg_end = _spin(0.80, 0.01, 1.0, 3, 0.05)
        self.ng = _spin(2.0, 1.0, 6.0, 2, 0.1)
        self.nog = _spin(2.0, 1.0, 6.0, 2, 0.1)
        self._gas_rows = [("Sgc (bağlı qaz)", self.sgc),
                          ("Sorg (qaza qarşı qalıq neft)", self.sorg),
                          ("krg @ 1-Swc-Sorg", self.krg_end),
                          ("Corey ng", self.ng), ("Corey nog", self.nog)]
        for label, widget in self._gas_rows:
            form.addRow(label, widget)
            widget.valueChanged.connect(self.changed)

        gas_note = QLabel(
            "Bu parametrlər simulyasiyaya TƏTBİQ OLUNUR — PVT tabında "
            "\"Qaz fazasını aktivləşdir\" işarələnibsə. SINAQ STATUSU: "
            "quyu öz BHP hədəfinə çox yaxınlaşan hallarda simulyasiya "
            "vaxtından əvvəl dayana bilər (bax PVT tabındakı qeyd) — "
            "proqram çökmür, son yığılmış nöqtəyə qədər nəticələr qalır.")
        gas_note.setWordWrap(True)
        gas_note.setStyleSheet(f"color:{PALETTE.oil};font-size:11px")
        form.addRow(gas_note)
        self._on_gas_toggled()

    def _on_gas_toggled(self):
        active = self.gas_enabled.isChecked()
        for _, widget in self._gas_rows:
            widget.setEnabled(active)

    def values(self) -> CoreyParameters:
        return CoreyParameters(self.swc.value(), self.sor.value(),
                               self.krw_end.value(), self.kro_end.value(),
                               self.nw.value(), self.no.value())

    def capillary_values(self) -> CapillaryParameters:
        return CapillaryParameters(entry_pressure=self.pc_entry.value(),
                                   lambda_exponent=self.pc_lambda.value(),
                                   max_pressure=self.pc_max.value())

    def gas_values(self) -> Optional[GasCoreyParameters]:
        """`None` — söndürülübsə (defolt).

        PVT tabında qaz aktivdirsə, bu dəyərlər mühərrikə ötürülür.
        Söndürülübsə, mühərrik (qaz aktiv olduğu halda) defolt Corey
        parametrlərinə keçir — bax `MainWindow.rebuild_model()`.
        """
        if not self.gas_enabled.isChecked():
            return None
        return GasCoreyParameters(self.sgc.value(), self.sorg.value(),
                                  self.krg_end.value(), self.ng.value(),
                                  self.nog.value())


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

        self.gas_phase_enabled = QCheckBox(
            "Qaz fazasını aktivləşdir (A7 — sınaq statusunda)")
        self.gas_phase_enabled.setChecked(False)
        self.gas_phase_enabled.stateChanged.connect(self.changed)

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

        form.addRow(self.gas_phase_enabled)
        gas_note = QLabel(
            "Üç fazalı mühərrik istifadə olunur. SINAQ STATUSU: quyu öz "
            "BHP hədəfinə çox yaxınlaşan hallarda simulyasiya vaxtından "
            "əvvəl (yığılmadan) dayana bilər — bu halda son yığılmış "
            "nöqtəyə qədər olan nəticələr göstərilir, proqram çökmür.")
        gas_note.setWordWrap(True)
        gas_note.setStyleSheet(f"color:{PALETTE.oil};font-size:11px")
        form.addRow(gas_note)

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
            bubble_point_bar=self.bubble_point.value(),
            include_gas=self.gas_phase_enabled.isChecked())

    def gas_phase_active(self) -> bool:
        return self.enabled.isChecked() and self.gas_phase_enabled.isChecked()


class WellPanel(QWidget):
    changed = pyqtSignal()

    COLUMNS = ["Ad", "i", "j", "K üst", "K alt", "Tip", "İdarə", "Qiymət", "rw"]

    def __init__(self):
        super().__init__()
        self.layer_count = 1
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
        self.table.itemChanged.connect(lambda *_: self.changed.emit())
        layout.addWidget(self.table)

        row = QHBoxLayout()
        self.add_button = QPushButton("+ Quyu")
        self.remove_button = QPushButton("− Sil")
        self.add_button.clicked.connect(lambda: self.add_row())
        self.remove_button.clicked.connect(self.remove_selected)
        row.addWidget(self.add_button)
        row.addWidget(self.remove_button)
        row.addStretch()
        layout.addLayout(row)

        hint = QLabel("BHP → bar,   RATE → m³/gün (rezervuar həcmi)\n"
                      "K üst / K alt → perforasiya intervalı (1-dən başlayır). "
                      "Su zonasındakı təbəqələri bağlamaq üçün K alt-ı azalt.")
        hint.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        layout.addWidget(hint)

    def set_layer_count(self, layer_count: int):
        """Grid dəyişəndə perforasiya intervalının yuxarı həddi yenilənir."""
        self.layer_count = max(int(layer_count), 1)

    def clamp_to_grid(self, nx: int, ny: int, nz: int) -> int:
        """Grid kiçiləndə quyu indekslərini avtomatik hüdud daxilinə salır.

        Səbəb: NX 41-dən 21-ə düşəndə cədvəldəki i = 40 mövcud olmayan
        hüceyrəyə işarə edir və model bloklanır. İstifadəçinin hər dəfə
        quyu sxemini yenidən tətbiq etməsi lazım gəlirdi.

        Qaytarır: düzəldilmiş xana sayı.
        """
        self.set_layer_count(nz)
        limits = {1: nx - 1, 2: ny - 1, 3: nz, 4: nz}
        minimums = {1: 0, 2: 0, 3: 1, 4: 1}
        corrected = 0

        self.table.blockSignals(True)
        try:
            for row in range(self.table.rowCount()):
                for column, maximum in limits.items():
                    item = self.table.item(row, column)
                    if item is None:
                        continue
                    try:
                        value = int(float(item.text()))
                    except ValueError:
                        continue
                    clamped = min(max(value, minimums[column]), maximum)
                    if clamped != value:
                        item.setText(str(clamped))
                        corrected += 1
                # K üst > K alt olarsa yerlərini dəyiş
                top, bottom = self.table.item(row, 3), self.table.item(row, 4)
                if top is not None and bottom is not None:
                    try:
                        if int(float(top.text())) > int(float(bottom.text())):
                            top.setText(bottom.text())
                            corrected += 1
                    except ValueError:
                        pass
        finally:
            self.table.blockSignals(False)
        return corrected

    def add_row(self, name="WELL", i=0, j=0, kind="PROD", mode="BHP",
                target=150.0, radius=0.1, k_top=1, k_bottom=None):
        k_bottom = self.layer_count if k_bottom is None else k_bottom
        r = self.table.rowCount()
        self.table.blockSignals(True)
        self.table.insertRow(r)
        for c, text in enumerate([name, str(i), str(j),
                                  str(int(k_top)), str(int(k_bottom))]):
            self.table.setItem(r, c, QTableWidgetItem(text))
        type_box = QComboBox()
        type_box.addItems(["PROD", "INJ"])
        type_box.setCurrentText(kind)
        type_box.currentIndexChanged.connect(lambda *_: self.changed.emit())
        self.table.setCellWidget(r, 5, type_box)
        mode_box = QComboBox()
        mode_box.addItems(["BHP", "RATE"])
        mode_box.setCurrentText(mode)
        self.table.setCellWidget(r, 6, mode_box)
        self.table.setItem(r, 7, QTableWidgetItem(f"{target:g}"))
        self.table.setItem(r, 8, QTableWidgetItem(f"{radius:g}"))
        self.table.blockSignals(False)
        self.changed.emit()

    def remove_selected(self):
        r = self.table.currentRow()
        if r >= 0:
            self.table.removeRow(r)
            self.changed.emit()

    def load(self, wells: List[Well]):
        self.table.blockSignals(True)
        self.table.setRowCount(0)
        self.table.blockSignals(False)
        for well in wells:
            perforations = well.perforations or [Perforation(0, 0, 0)]
            layers = [p.k for p in perforations]
            self.add_row(well.name, perforations[0].i, perforations[0].j,
                         well.well_type.value, well.control.mode.value,
                         well.control.target, well.radius,
                         k_top=min(layers) + 1, k_bottom=max(layers) + 1)

    def values(self) -> List[Well]:
        wells: List[Well] = []
        for r in range(self.table.rowCount()):
            try:
                name = self.table.item(r, 0).text().strip() or f"W{r + 1}"
                i = int(float(self.table.item(r, 1).text()))
                j = int(float(self.table.item(r, 2).text()))
                k_top = int(float(self.table.item(r, 3).text()))
                k_bottom = int(float(self.table.item(r, 4).text()))
                kind = self.table.cellWidget(r, 5).currentText()
                mode = self.table.cellWidget(r, 6).currentText()
                target = float(self.table.item(r, 7).text())
                radius = float(self.table.item(r, 8).text())
            except (AttributeError, ValueError):
                continue

            first = max(min(k_top, k_bottom), 1)
            last = min(max(k_top, k_bottom), self.layer_count)
            if last < first:
                last = first
            wells.append(Well(
                name=name,
                well_type=WellType(kind),
                control=WellControl(ControlMode(mode), target),
                perforations=[Perforation(i, j, k - 1)
                              for k in range(first, last + 1)],
                radius=radius,
            ))
        return wells


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
        self.use_goc = QCheckBox(
            "Qaz papağı (GOC) — yalnız qaz fazası aktivdirsə təsir edir")
        self.goc = _spin(2010.0, 0.0, 8000.0, 1, 10.0, "m")
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
        form.addRow(self.use_goc)
        form.addRow("Qaz-neft kontaktı (GOC)", self.goc)
        self.use_goc.stateChanged.connect(self.changed)
        self.goc.valueChanged.connect(self.changed)
        note = QLabel("Söndürülübsə, bütün hüceyrələrdə eyni təzyiq və Sw "
                      "işlədilir (köhnə davranış).")
        note.setWordWrap(True)
        note.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        form.addRow(note)

    def initial_conditions(self) -> InitialConditions:
        equilibrate = self.use_equilibration.isChecked()
        use_goc = equilibrate and self.use_goc.isChecked()
        return InitialConditions(
            datum_depth=self.datum_depth.value(),
            datum_pressure=self.initial_pressure.value(),
            water_saturation=self.initial_sw.value(),
            oil_water_contact=self.owc.value() if equilibrate else None,
            gas_oil_contact=self.goc.value() if use_goc else None,
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
