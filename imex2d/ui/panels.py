"""Giriş panelləri.

Hər panelin yeganə işi: istifadəçi girişini DOMAIN obyektinə çevirmək.
Panellər hesablama aparmır, yoxlama etmir və simulyatordan xəbərsizdir.
Yoxlama ReservoirModel.validate() və SimulationConfig.validate()
metodlarındadır.
"""

from __future__ import annotations
from typing import Dict, List, Optional

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QBrush, QColor
from PyQt5.QtWidgets import (QCheckBox, QComboBox, QDialog, QDoubleSpinBox,
                             QGroupBox,
                             QFileDialog, QFormLayout, QHBoxLayout,
                             QHeaderView, QLabel, QLineEdit, QMessageBox,
                             QPushButton, QSpinBox, QTableWidget,
                             QTableWidgetItem, QTextEdit, QVBoxLayout, QWidget)

from ..application.config import (MPFA_O, TPFA, LinearSolverConfig,
                                  OutputConfig, SimulationConfig,
                                  TimeSteppingConfig)
import os

from ..application.scenarios import WELL_PATTERNS
from ..geology.interpolation import (INTERPOLATORS, InverseDistance,
                                     NearestNeighbour, OrdinaryKriging)
from ..application.geology_adapter import well_layer_summary
from ..application.geology_service import (CompletionMethod, CompletionSpec,
                                           LayerInterpolationConfig)
from ..domain.data_availability import format_layers, parse_layers
from ..domain.geology import GeologicalWell, validate_wells, well_effective_layers
from ..domain.geometry import depth_to_k, xy_to_ij
from ..geology.layer_availability import LayerDataPolicy
from ..domain.structure import FaultReference
from ..io.fault_io import (FaultFormatError, read_eclipse_faults,
                          read_faults_csv)
from ..io.scal_io import (ScalFormatError, read_scal_csv, read_sgof,
                          read_swof)
from ..domain.grid import CartesianGrid
from ..domain.initial import InitialConditions
from ..domain.properties import FluidProperties
from ..domain.scal import (CapillaryParameters, CoreyParameters,
                           GasCoreyParameters)
from ..domain.unit_conversions import convert, to_engine_units
from ..domain.tubing import TubingGeometry
from ..domain.wells import (ControlMode, Perforation, Phase, RateBasis, Well,
                            WellControl, WellType)
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


def _unit_combo(options, default: str) -> QComboBox:
    """Vahid seçici — DEFOLT HƏMİŞƏ mühərrik vahididir (`options[0]`),
    ona görə toxunulmayan panel əvvəlki kimi (mühərrik vahidində)
    davranır — 'Mövcud istifadəçilər əvvəlki ədədi defoltları görməlidir'
    tələbi budur."""
    box = QComboBox()
    box.addItems(options)
    box.setCurrentText(default)
    return box


def _bind_unit_aware_spins(spins, combo: QComboBox, quantity: str) -> None:
    """Vahid seçici dəyişəndə spin-box-ların DİAPAZONUNU və GÖSTƏRİLƏN
    ƏDƏDİ yeni vahidə YENİDƏN HESABLAYIR (fiziki kəmiyyəti SAXLAYIR).

    Bunsuz: diapazon (`setRange`) HƏMİŞƏ ilk (defolt) vahidin miqyasında
    qalırdı — məs. təzyiq sahəsi bar üçün 1–1200 diapazonlu idi; istifadəçi
    "psi" seçib 3000 yazanda dəyər SƏSSİZCƏ 1200-ə (~83 bar-a) kəsilirdi
    (bax `test_ui_units.py`-də tutulan reqressiya). Bu funksiya hər spin-un
    ilkin (`__init__`-də verilmiş) min/maks-ını REFERANS kimi saxlayıb hər
    vahid dəyişimində yeni vahidə çevirir.
    """
    spins = list(spins)
    base_unit = combo.currentText()
    base_ranges = [(spin.minimum(), spin.maximum()) for spin in spins]

    def on_unit_changed():
        new_unit = combo.currentText()
        for spin, (base_lo, base_hi) in zip(spins, base_ranges):
            old_value = spin.value()
            new_lo = convert(base_lo, base_unit, new_unit, quantity)
            new_hi = convert(base_hi, base_unit, new_unit, quantity)
            new_value = convert(old_value, on_unit_changed.previous_unit, new_unit, quantity)
            spin.blockSignals(True)
            spin.setRange(min(new_lo, new_hi), max(new_lo, new_hi))
            spin.setValue(new_value)
            spin.blockSignals(False)
        on_unit_changed.previous_unit = new_unit

    on_unit_changed.previous_unit = base_unit
    combo.currentIndexChanged.connect(on_unit_changed)


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

        # ── A2: struktur QUYULARDAN ─────────────────────────────────
        # İşarələnəndə grid həndəsəsi geologiya cədvəlindəki «lay üstü»/
        # «lay altı» dərinliklərinin interpolyasiyasından qurulur; tavan
        # dərinliyi, maillik və (mənbə "quyulardan" olanda) DZ ARTIQ
        # İŞTİRAK ETMİR, ona görə həmin sahələr SÖNDÜRÜLÜR — istifadəçi
        # işləməyən dəyər doldurmasın (bax `_on_structure_toggled`).
        self.structure_from_wells = QCheckBox(
            "Struktur quyulardan (lay üstü/altı)")
        self.structure_from_wells.setToolTip(
            "Grid həndəsəsi geologiya cədvəlindəki lay üstü/altı dərinliklərindən "
            "qurulur (əyri səth, sütundan sütuna dəyişən qalınlıq).\n"
            "Yalnız 'İnterpolyasiya et' düyməsi ilə qurulan modelə aiddir.")
        self.thickness_source = QComboBox()
        self.thickness_source.addItem("Quyulardan (üst + alt)", "wells")
        self.thickness_source.addItem("Sabit (aşağıdakı DZ)", "constant")
        self.thickness_source.setToolTip(
            "Quyulardan: qalınlıq = lay altı − lay üstü.\n"
            "Sabit: tavan quyulardan (əyri səth), qalınlıq isə verilmiş DZ cəmi — "
            "praktikada çox vaxt yalnız lay tavanı məlum olur.")
        self.structure_from_wells.toggled.connect(self._on_structure_toggled)
        self.structure_from_wells.toggled.connect(self.changed)
        self.thickness_source.currentIndexChanged.connect(self._on_structure_toggled)
        self.thickness_source.currentIndexChanged.connect(self.changed)

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

        form.addRow(self.structure_from_wells)
        form.addRow("Qalınlıq mənbəyi (struktur)", self.thickness_source)

        self._sync_table_rows(self.nz.value())

        self.info = QLabel()
        self.info.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        form.addRow(self.info)
        self.changed.connect(self._refresh_info)
        # `_on_mode_changed` məlumat sətrini yeniləyir, ona görə YALNIZ
        # `self.info` yaradıldıqdan sonra çağırıla bilər.
        self._on_mode_changed()
        self._on_structure_toggled()

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

    #: Struktur rejimində İŞLƏMƏYƏN sahələrin tooltip mətni — səssiz
    #: "parametr təsir etmir" halının qarşısını alır (bax A2, §3.5).
    _UNUSED_TOOLTIP = ("Bu dəyər struktur rejimində istifadə olunmur — "
                       "həndəsə quyuların lay üstü/altı dərinliklərindən qurulur.")

    def _on_structure_toggled(self, *_args) -> None:
        """Struktur rejimində iştirak etməyən sahələri SÖNDÜRÜR.

        `top_depth`/`dip_x`/`dip_y` HƏR İKİ qalınlıq mənbəyində
        iştirak etmir; `dz` (və onunla bağlı `base_depth`/`per_layer`/
        cədvəl) yalnız mənbə "quyulardan" olanda söndürülür — "sabit"
        rejimində qalınlıq məhz DZ-dən gəlir.
        """
        structural = self.structure_from_wells.isChecked()
        self.thickness_source.setEnabled(structural)
        thickness_from_wells = (structural
                                and self.thickness_source.currentData() == "wells")

        for widget in (self.top_depth, self.dip_x, self.dip_y):
            widget.setEnabled(not structural)
            widget.setToolTip(self._UNUSED_TOOLTIP if structural else "")
        for widget in (self.thickness_mode, self.base_depth, self.per_layer,
                       self.dz_table, self.dz):
            widget.setToolTip(self._UNUSED_TOOLTIP if thickness_from_wells else "")

        if thickness_from_wells:
            for widget in (self.thickness_mode, self.base_depth, self.per_layer,
                           self.dz_table, self.dz):
                widget.setEnabled(False)
        else:
            # Söndürülməmiş vəziyyəti öz normal məntiqi bərpa etsin
            # (rejim/`per_layer` asılılıqları burada TƏKRARLANMIR).
            self.thickness_mode.setEnabled(True)
            self.dz_table.setEnabled(True)
            self._on_mode_changed()
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

    def structure_values(self) -> dict:
        """A2 struktur seçimləri — `values()`-dən AYRICA saxlanılır,
        çünki `values()` sintetik model qurucusuna (`geology_builder.
        build(**grid_values)`) BİRBAŞA açılır və ona bu açarlar
        naməlumdur."""
        return dict(structure_from_wells=self.structure_from_wells.isChecked(),
                    thickness_source=self.thickness_source.currentData())

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
    cross_validate_requested = pyqtSignal()

    # "Lay üstü/altı" (FİZİKİ interval) və "Data layları" (MƏLUMAT
    # MÖVCUDLUĞU) AYRI sütunlardır və AYRI mənaları var (tapşırıq §1/§14).
    # "Kəsdiyi laylar" YALNIZ OXUNAN, hesablanan sütundur — istifadəçi
    # onu redaktə etmir, çünki o, top/bottom + grid həndəsəsinin
    # NƏTİCƏSİDİR, müstəqil giriş deyil.
    COLUMNS = ["Ad", "Modeldə", "X, m", "Y, m", "(i, j)", "Lay üstü, m",
              "Lay altı, m", "Data layları", "Kəsdiyi laylar",
              "φ", "k, mD", "Sw", "Qeyd"]
    COL_NAME = 0
    COL_IN_MODEL = 1
    COL_X = 2
    COL_Y = 3
    COL_IJ = 4
    COL_TOP = 5
    COL_BOTTOM = 6
    COL_DATA_LAYERS = 7
    COL_EFFECTIVE = 8
    COL_PORO = 9
    COL_PERM = 10
    COL_SW = 11
    COL_NOTE = 12
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
        self.add_button.clicked.connect(self._on_add_clicked)
        self.duplicate_button.clicked.connect(self._duplicate_selected)
        self.delete_button.clicked.connect(self._delete_selected)
        self.centre_button.clicked.connect(self._centre_selected)
        layout.addLayout(toolbar)

        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        # 13 sütun `Stretch` rejimində başlıqları KƏSİRDİ ("Data layları"
        # → "ata layl"), yəni yeni sütunların ADI oxunmurdu. Həll: `Stretch`
        # SAXLANILIR (geniş pəncərədə sütunlar boşluğu BƏRABƏR bölür, tək
        # bir sütun onu udmur), amma minimum sütun eni ƏN UZUN BAŞLIĞA görə
        # təyin olunur. Ölçü font metrikasından hesablanır — sabit piksel
        # DEYİL, ona görə fərqli DPI/şriftdə də, gələcəkdə yeni sütun
        # əlavə olunanda da başlıqlar tam görünür; pəncərə çox darsa
        # (məs. sol alət qutusu) üfüqi sürüşdürmə çubuğu çıxır.
        header = self.table.horizontalHeader()
        metrics = header.fontMetrics()
        widest = max(metrics.horizontalAdvance(text) for text in self.COLUMNS)
        header.setMinimumSectionSize(widest + 22)
        header.setSectionResizeMode(QHeaderView.Stretch)
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
        self.range_v = _spin(0.0, 0.0, 100000.0, 0, 5.0, "m")
        self.sill = _spin(0.0, 0.0, 1e6, 4, 0.01)
        self.nugget = _spin(0.0, 0.0, 1e6, 4, 0.01)
        self.min_neighbors = _ispin(1, 1, 999)
        self.max_neighbors = _ispin(0, 0, 999)
        rows = [("IDW dərəcəsi p", self.power),
                ("Axtarış radiusu (0 = limitsiz)", self.search_radius),
                ("Kriging radiusu a — üfüqi (0 = avto)", self.range_),
                ("Kriging radiusu — şaquli (0 = üfüqi ilə eyni)", self.range_v),
                ("Sill c (0 = avto)", self.sill),
                ("Nugget c₀", self.nugget),
                ("Min qonşu sayı (yerli axtarışda)", self.min_neighbors),
                ("Maks qonşu sayı (0 = limitsiz)", self.max_neighbors)]
        for label, widget in rows:
            form.addRow(label, widget)
            widget.valueChanged.connect(self._on_table_edited)

        self.allow_cross_layer_fallback = QCheckBox(
            "Boş laya digər laylardan borc verməyə icazə ver (ekstrapolyasiya)")
        self.allow_cross_layer_fallback.setToolTip(
            "Söndürülübsə (defolt): heç bir quyu nöqtəsi olmayan lay üçün "
            "açıq xəta verilir — başqa layların dəyəri sükutla işlədilmir.\n"
            "Yandırılsa: 3D/anizotrop Kriging seçilibsə yaxın laylar uzaq "
            "laylardan çox təsir edir (kor-koranə bərabər hovuzlama yox).")
        self.allow_cross_layer_fallback.toggled.connect(self._on_table_edited)
        form.addRow("", self.allow_cross_layer_fallback)

        # ── LAY-MƏLUMATLI (layer-aware) rejim ─────────────────────────
        # SÖNDÜRÜLÜ olanda proqram TAM ƏVVƏLKİ kimi işləyir (§25).
        self.layer_aware = QCheckBox(
            "Lay-məlumatlı rejim (hansı layda HƏQİQƏTƏN data var)")
        self.layer_aware.setToolTip(
            "Yandırılanda quyunun 'Data layları' sütunu oxunur və hər lay "
            "YALNIZ ÖZ məlumatı ilə interpolyasiya olunur.\n"
            "Məlumatı olmayan lay sükutla doldurulmur — MISSING qalır və "
            "simulyasiyadan əvvəl açıq xəta verir (tamamlama üsulu "
            "seçilməyibsə).\n"
            "Söndürülü (defolt): əvvəlki davranış — bir quyu dəyəri bütün "
            "K-lara yayılır.")
        self.layer_aware.toggled.connect(self._on_layer_mode_changed)
        self.layer_aware.toggled.connect(self._on_table_edited)
        form.addRow("", self.layer_aware)

        self.data_policy = QComboBox()
        self.data_policy.addItem("Ciddi — yalnız bəyan olunmuş laylar", "strict")
        self.data_policy.addItem("İnterval — top/bottom-un kəsdiyi laylar (FƏRZİYYƏ)",
                                 "interval")
        self.data_policy.setToolTip(
            "Ciddi (tövsiyə olunan): 'Data layları' sütunu boş olan quyu heç "
            "bir laya məlumat vermir.\n"
            "İnterval: bəyan olmayan quyunun fiziki intervalı məlumat sayılır — "
            "bu, AÇIQ FƏRZİYYƏDİR və hesabatda xəbərdarlıq kimi görünür.")
        self.data_policy.currentIndexChanged.connect(self._on_table_edited)
        form.addRow("Məlumat mövcudluğu siyasəti", self.data_policy)

        self.target_layers = QLineEdit("")
        self.target_layers.setPlaceholderText("boş = məlumatı olan bütün laylar (məs. 1-3)")
        self.target_layers.setToolTip(
            "İnterpolyasiya HƏDƏFİ — 1-əsaslı lay nömrələri ('1-3' və ya '1,2,5').\n"
            "Bu, məlumat mövcudluğundan AYRI anlayışdır: seçilmiş, lakin "
            "məlumatı olmayan lay interpolyasiya OLUNMUR.")
        self.target_layers.textChanged.connect(self._on_table_edited)
        form.addRow("İnterpolyasiya layları", self.target_layers)

        self.completion_method = QComboBox()
        for label, value in (
                ("Yoxdur — MISSING qalsın (defolt)", CompletionMethod.NONE.value),
                ("Orijinal sahəni saxla", CompletionMethod.PRESERVE_ORIGINAL.value),
                ("Şaquli trend (ESTIMATED)", CompletionMethod.VERTICAL_TREND.value),
                ("3D geostatistik qiymətləndirmə (ESTIMATED)",
                 CompletionMethod.GEOSTATISTICAL_3D.value),
                ("SGS realizasiyası (SIMULATED)", CompletionMethod.SGS.value),
                ("Sabit lay dəyəri (ESTIMATED)", CompletionMethod.CONSTANT.value)):
            self.completion_method.addItem(label, value)
        self.completion_method.setToolTip(
            "Məlumatı OLMAYAN laylar necə tamamlansın. Heç bir variant "
            "nəticəni 'ölçülmüş' kimi qeyd etmir — mənşə (provenance) "
            "modeldə saxlanılır.")
        self.completion_method.currentIndexChanged.connect(self._on_completion_changed)
        self.completion_method.currentIndexChanged.connect(self._on_table_edited)
        form.addRow("Məlumatsız layların tamamlanması", self.completion_method)

        self.completion_value = _spin(0.0, -1e9, 1e9, 5, 0.01)
        self.completion_value.setToolTip(
            "Yalnız 'Sabit lay dəyəri' üsulu üçün — həmin laylara veriləcək "
            "dəyər (xassənin fiziki vahidində).")
        self.completion_value.valueChanged.connect(self._on_table_edited)
        form.addRow("Sabit tamamlama dəyəri", self.completion_value)
        layout.addLayout(form)

        self.layer_summary = QTextEdit()
        self.layer_summary.setReadOnly(True)
        self.layer_summary.setMaximumHeight(90)
        self.layer_summary.setStyleSheet("font-family:monospace;font-size:11px")
        self.layer_summary.setToolTip(
            "Quyu üzrə: kəsdiyi laylar / məlumatı olan laylar — "
            "interpolyasiya və MISSING laylar isə 'İnterpolyasiya et'-dən "
            "sonra hesabatda görünür.")
        layout.addWidget(self.layer_summary)

        action_row = QHBoxLayout()
        self.interpolate_button = QPushButton("İnterpolyasiya et")
        self.interpolate_button.clicked.connect(self.interpolate_requested)
        action_row.addWidget(self.interpolate_button)
        self.cross_validate_button = QPushButton("Cross-validation et")
        self.cross_validate_button.setToolTip(
            "Quyu nöqtələrindən bəzilərini müvəqqəti gizlədib qalanları ilə "
            "proqnozlaşdırır, real RMSE/MAE/R²/MAPE göstərir. '100% dəqiq' "
            "vəd etmir — az nöqtə ilə R² mənfi ola bilər.")
        self.cross_validate_button.clicked.connect(self.cross_validate_requested)
        action_row.addWidget(self.cross_validate_button)
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
        self._on_layer_mode_changed()
        self._on_completion_changed()
        self._refresh_map_and_validation()

    # ------------------------------------------------------------ slots
    def _on_layer_mode_changed(self):
        enabled = self.layer_aware.isChecked()
        for widget in (self.data_policy, self.target_layers,
                       self.completion_method, self.completion_value):
            widget.setEnabled(enabled)
        if enabled:
            self._on_completion_changed()
        # Köhnə "boş laya borc ver" seçimi lay-məlumatlı rejimin
        # tamamlama mexanizmi ilə ƏVƏZLƏNİR — ikisi eyni anda mənasızdır.
        self.allow_cross_layer_fallback.setEnabled(
            not enabled and "Kriging" in self.method.currentText())

    def _on_completion_changed(self):
        self.completion_value.setEnabled(
            self.layer_aware.isChecked()
            and self.completion_method.currentData() == CompletionMethod.CONSTANT.value)

    def _on_method_changed(self):
        method = self.method.currentText()
        is_kriging = "Kriging" in method
        self.power.setEnabled("IDW" in method)
        self.search_radius.setEnabled("IDW" in method or is_kriging)
        for widget in (self.range_, self.range_v, self.sill, self.nugget,
                      self.min_neighbors, self.max_neighbors,
                      self.allow_cross_layer_fallback):
            widget.setEnabled(is_kriging)
        if self.layer_aware.isChecked():
            self.allow_cross_layer_fallback.setEnabled(False)

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

    def _on_add_clicked(self):
        """`add_button.clicked` `bool checked` göndərir — `add_row(well=None)`
        birbaşa qoşulsaydı, bu bool `well` parametrinə düşərdi."""
        self.add_row()

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
            radius = self.search_radius.value()
            max_n = self.max_neighbors.value()
            return OrdinaryKriging(
                range_=self.range_.value() or None,
                range_v=self.range_v.value() or None,
                sill=self.sill.value() or None,
                nugget=self.nugget.value(),
                search_radius=radius if radius > 0 else None,
                min_neighbors=self.min_neighbors.value(),
                max_neighbors=max_n if max_n > 0 else None)
        return NearestNeighbour()

    def cross_layer_fallback_allowed(self) -> bool:
        """Lay-məlumatlı rejim AÇIQ olanda HƏMİŞƏ `False` — köhnə "boş laya
        borc ver" yolu ilə yeni tamamlama mexanizmi eyni anda işləsəydi,
        məlumatsız lay iki müxtəlif məntiqlə doldurulardı."""
        if self.layer_aware.isChecked():
            return False
        return self.allow_cross_layer_fallback.isChecked()

    # ─────────────────────────────────────── lay-məlumatlı rejim (public)
    def layer_aware_enabled(self) -> bool:
        return self.layer_aware.isChecked()

    def layer_data_policy(self) -> LayerDataPolicy:
        if not self.layer_aware.isChecked():
            return LayerDataPolicy.BROADCAST
        return (LayerDataPolicy.INTERVAL
                if self.data_policy.currentData() == "interval"
                else LayerDataPolicy.STRICT)

    def layer_config(self, nz: int) -> Optional[LayerInterpolationConfig]:
        """Paneldəki seçimlərdən `LayerInterpolationConfig` qurur.

        Rejim söndürülübsə `None` — `build()` ƏVVƏLKİ yolu ilə gedir.
        Mətn səhvdirsə `ValueError` (çağıran istifadəçiyə göstərir) —
        SƏSSİZ DÜZƏLİŞ YOXDUR.
        """
        if not self.layer_aware.isChecked():
            return None
        targets = parse_layers(self.target_layers.text(), nz)
        method = CompletionMethod(self.completion_method.currentData())
        spec = CompletionSpec(
            method=method,
            value=(self.completion_value.value()
                   if method is CompletionMethod.CONSTANT else None))
        return LayerInterpolationConfig(
            policy=self.layer_data_policy(),
            target_layers=targets or None,
            completion=spec)

    def refresh_layer_summary(self) -> None:
        """Quyu-üzrə "kəsdiyi / məlumatı olan laylar" xülasəsi (§14).

        HESABLAMA APARMIR — yalnız bəyanları oxuyur (`well_layer_summary`),
        ona görə hər cədvəl redaktəsində çağırıla bilər.
        """
        if self._geometry is None:
            self.layer_summary.setPlainText("Grid qurulduqdan sonra göstərilir.")
            return
        wells = self.wells()
        if not wells:
            self.layer_summary.setPlainText("Quyu yoxdur.")
            return
        try:
            summary = well_layer_summary(wells, self._geometry,
                                         policy=self.layer_data_policy())
        except ValueError as error:
            self.layer_summary.setPlainText(f"Lay bəyanı oxunmadı: {error}")
            return
        lines = []
        for name, entry in summary.items():
            lines.append(
                f"{name}: kəsir {format_layers(entry.get('effective', []))} | "
                f"data {format_layers(entry.get('data', []))} | "
                + " ".join(f"{key} {format_layers(entry[key])}"
                           for key in ("PORO", "PERMX", "SW") if key in entry))
        self.layer_summary.setPlainText("\n".join(lines))

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
        try:
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

            self.table.setItem(r, self.COL_DATA_LAYERS,
                               QTableWidgetItem(well.data_layers_text))
            effective_item = QTableWidgetItem("—")
            effective_item.setFlags(effective_item.flags() & ~Qt.ItemIsEditable)
            effective_item.setForeground(QBrush(QColor(Qt.gray)))
            self.table.setItem(r, self.COL_EFFECTIVE, effective_item)

            optional_columns = [(self.COL_TOP, "top"), (self.COL_BOTTOM, "bottom"),
                               (self.COL_PORO, "porosity"), (self.COL_PERM, "permeability"),
                               (self.COL_SW, "water_saturation")]
            for column, attr in optional_columns:
                value = getattr(well, attr)
                self.table.setItem(r, column, QTableWidgetItem(
                    "" if value is None else f"{value:g}"))
            self.table.setItem(r, self.COL_NOTE, QTableWidgetItem(well.note))
        finally:
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
        layers_item = self.table.item(row, self.COL_DATA_LAYERS)
        return GeologicalWell(
            name=name_item.text().strip(),
            in_model=in_model,
            x=values["x"] or 0.0, y=values["y"] or 0.0,
            top=values["top"], bottom=values["bottom"],
            porosity=values["porosity"], permeability=values["permeability"],
            water_saturation=values["water_saturation"],
            note=note_item.text() if note_item is not None else "",
            data_layers_text=(layers_item.text().strip()
                              if layers_item is not None else ""))

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
        """(i, j) VƏ "Kəsdiyi laylar" sütunlarını yeniləyir.

        "Kəsdiyi laylar" YALNIZ `top`/`bottom` + grid həndəsəsindən
        hesablanır — "Data layları" sütununa TƏSİR ETMİR və ondan
        ASILI DEYİL (tapşırıq §1: interval ≠ məlumat)."""
        self.table.blockSignals(True)
        try:
            for row in range(self.table.rowCount()):
                item = self.table.item(row, self.COL_IJ)
                effective_item = self.table.item(row, self.COL_EFFECTIVE)
                if item is None:
                    continue
                if self._geometry is None:
                    item.setText("grid qurulduqdan sonra")
                    if effective_item is not None:
                        effective_item.setText("grid qurulduqdan sonra")
                    continue
                well = self._well_from_row(row)
                if well is None:
                    continue
                x_max, y_max = self._geometry.areal_extent()
                if not (0.0 <= well.x <= x_max and 0.0 <= well.y <= y_max):
                    item.setText("kənar")
                    if effective_item is not None:
                        effective_item.setText("kənar")
                    continue
                i, j = xy_to_ij(well.x, well.y, self._geometry)
                item.setText(f"({i}, {j})")
                if effective_item is not None:
                    try:
                        effective_item.setText(
                            format_layers(well_effective_layers(well, self._geometry)))
                    except ValueError:
                        effective_item.setText("interval səhvdir")
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
        self.refresh_layer_summary()


class FaciesPanel(QWidget):
    """Fasiya/SIS konfiqurasiyası (Phase 4.1 §9) — YALNIZ minimal
    parametrlər: sütun adı, seed, realizasiya sayı, nisbətlər. Qabaqcıl
    parametrlər (variogram, anizotropluq, qonşuluq axtarışı) BU PANELDƏ
    YOXDUR — proqram səviyyəsində (`geology_service.FaciesBuildConfig`)
    əlçatandır, "Do NOT expose every advanced parameter" qaydasına görə.

    QEYD (bilərəkdən məhdud əhatə): bu panel HƏLƏ `main_window.py`-nin
    əsas saxla/yüklə iş axınına BAĞLANMAYIB — bu fazanın fokusu backend
    inteqrasiyasıdır (bax `FACIES.md`/hesabat "Qalan iş"). `build_config()`
    dəyərləri oxumaq üçün hazırdır, layihə serializasiyası onları HƏLƏ
    saxlamır.
    """
    changed = pyqtSignal()

    def __init__(self):
        super().__init__()
        form = QFormLayout(self)
        self.column_name = QLineEdit("FACIES")
        self.seed = _ispin(0, 0, 999999)
        self.n_realizations = _ispin(1, 1, 500)
        self.proportions_text = QLineEdit("")
        hint = QLabel("Nisbətlər formatı: kod:nisbət, kod:nisbət, … (cəm 1.0 olmalıdır). "
                      "Boş buraxılsa quyu datasından müşahidə olunan nisbətlər işlədilir "
                      "(bax hesabatda xəbərdarlıq).")
        hint.setWordWrap(True)
        hint.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        for label, widget in [("Fasiya sütunu", self.column_name),
                              ("Seed", self.seed),
                              ("Realizasiya sayı", self.n_realizations),
                              ("Nisbətlər (kod:nisbət)", self.proportions_text)]:
            form.addRow(label, widget)
            sig = getattr(widget, "textChanged", None) or widget.valueChanged
            sig.connect(self.changed)
        form.addRow(hint)

    def parse_proportions(self) -> Optional[dict]:
        """Boş mətn -> `None` (müşahidə olunan nisbətlərə keçid, bax
        `geology_service._simulate_categorical_field`). Səhv formatda
        `ValueError` — SƏSSİZCƏ boş/yanlış nisbətə keçilmir."""
        text = self.proportions_text.text().strip()
        if not text:
            return None
        proportions: Dict[int, float] = {}
        for token in text.split(","):
            token = token.strip()
            if not token:
                continue
            code_text, sep, value_text = token.partition(":")
            if not sep:
                raise ValueError(f"Nisbət formatı yanlışdır: {token!r} (gözlənilən 'kod:nisbət').")
            proportions[int(code_text.strip())] = float(value_text.strip())
        return proportions

    def build_config(self, realization_id: int = 0, seed_offset: int = 0):
        from ..application.geology_service import FaciesBuildConfig
        return FaciesBuildConfig(proportions=self.parse_proportions(),
                                 seed=self.seed.value() + seed_offset,
                                 realization_id=realization_id)

    def column_name_value(self) -> str:
        return self.column_name.text().strip()

    def realization_count(self) -> int:
        return self.n_realizations.value()


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
        #: Qaz-neft cədvəlləri (SGOF) — G4; ayrıca saxlanılır, çünki
        #: su-neft cədvəli olmadan da yüklənə bilər.
        self.gas_tables = None
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
        self.load_sgof = QPushButton("Eclipse SGOF yüklə…")
        self.load_csv.clicked.connect(self._load_csv)
        self.load_swof.clicked.connect(self._load_swof)
        self.load_sgof.clicked.connect(self._load_sgof)
        buttons.addWidget(self.load_csv)
        buttons.addWidget(self.load_swof)
        buttons.addWidget(self.load_sgof)
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
        self.load_sgof.setEnabled(by_table)

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

    def _load_sgof(self):
        """Qaz-neft cədvəli (G4) — AYRI saxlanılır, `self.tables`-ı əvəz
        etmir: model hər ikisini eyni vaxtda işlədə bilər."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Eclipse deck (SGOF)", "",
            "Eclipse (*.DATA *.data *.inc);;Bütün fayllar (*)")
        if not path:
            return
        try:
            self.gas_tables = read_sgof(path)
        except ScalFormatError as error:
            QMessageBox.warning(self, "Cədvəl oxunmadı", str(error))
            return
        except Exception as error:
            QMessageBox.critical(self, "Cədvəl oxunmadı",
                                 f"Gözlənilməz xəta: {error}")
            return
        self.info.setText(f"{os.path.basename(path)} (SGOF)\n"
                          f"{len(self.gas_tables)} region\n"
                          + self.gas_tables.summary())
        self.mode.setCurrentIndex(1)
        self.changed.emit()

    def is_enabled(self) -> bool:
        return self.mode.currentData() == "TABLE" and self.tables is not None

    def gas_tables_enabled(self) -> bool:
        return (self.mode.currentData() == "TABLE"
                and self.gas_tables is not None)


class RockFluidPanel(QWidget):
    changed = pyqtSignal()

    #: YALNIZ geoloji modeli dəyişdirən sahələr (φ, K, heterogenlik).
    #: `changed`-dən ayrıdır, çünki lözlük dəyişəndə interpolyasiyanı
    #: köhnəlmiş saymaq YANLIŞ olardı — flüid geologiyaya təsir etmir.
    geology_changed = pyqtSignal()

    def __init__(self):
        super().__init__()
        form = QFormLayout(self)
        self.porosity = _spin(0.22, 0.01, 0.45, 3, 0.01)
        # decimals=4 (mD üçün lazım olandan çox) — vahid "D"-yə keçəndə
        # (1 D = 1000 mD) tipik dəyərlər (0.001-20 D) hələ də mənalı
        # göstərilsin deyə. "m2" tipik reservoir dəyərləri (~1e-16 m²)
        # ÜÇÜN İSƏ HEÇ BİR spin-box dəqiqliyi kifayət deyil — bax UNITS.md.
        self.permx = _spin(150, 0.01, 20000, 4, 10, "mD")
        self.permx_unit = _unit_combo(["mD", "D", "m2"], "mD")
        _bind_unit_aware_spins([self.permx], self.permx_unit, "permeability")
        self.ky_over_kx = _spin(1.0, 0.01, 10, 2, 0.1)
        self.kv_over_kh = _spin(0.10, 0.001, 1.0, 3, 0.01)
        self.heterogeneity = QComboBox()
        self.heterogeneity.addItems(["Homogen", "Təsadüfi (log-normal)"])
        self.sigma = _spin(0.5, 0.05, 2.0, 2, 0.05)
        self.seed = _ispin(7, 0, 9999)
        # decimals=6 — "Pa.s" seçiləndə (1 cP = 0.001 Pa·s) tipik neft
        # lözlüyü (0.5-5000 cP -> 0.0005-5 Pa·s) mənalı göstərilsin.
        self.mu_w = _spin(0.5, 0.05, 50, 6, 0.05, "cP")
        self.mu_o = _spin(3.0, 0.05, 5000, 6, 0.5, "cP")
        self.viscosity_unit = _unit_combo(["cP", "Pa.s"], "cP")
        _bind_unit_aware_spins([self.mu_w, self.mu_o], self.viscosity_unit, "viscosity")
        self.bo = _spin(1.15, 1.0, 3.0, 3, 0.01)
        self.rock_compressibility = _spin(4.5e-5, 1e-6, 1e-3, 7, 1e-5, "1/bar")
        rows = [("Məsaməlilik φ", self.porosity), ("Keçiricilik Kx", self.permx),
                ("Kx vahidi", self.permx_unit),
                ("Ky/Kx", self.ky_over_kx), ("Kv/Kh (şaquli)", self.kv_over_kh),
                ("Heterogenlik", self.heterogeneity),
                ("σ (log-normal)", self.sigma), ("Seed", self.seed),
                ("Su lözlüyü μw", self.mu_w), ("Neft lözlüyü μo", self.mu_o),
                ("Lözlük vahidi", self.viscosity_unit),
                ("Bo", self.bo), ("Süxur sıxılması", self.rock_compressibility)]
        for label, widget in rows:
            form.addRow(label, widget)
            sig = getattr(widget, "valueChanged", None) or widget.currentIndexChanged
            sig.connect(self.changed)

        # Geologiyaya təsir edən sahələr ayrıca siqnal verir — bax
        # `geology_changed`.
        for widget in (self.porosity, self.permx, self.permx_unit,
                       self.ky_over_kx, self.kv_over_kh,
                       self.heterogeneity, self.sigma, self.seed):
            sig = getattr(widget, "valueChanged", None) or widget.currentIndexChanged
            sig.connect(self.geology_changed)

        # G6 — süxur sıxılmasının istinad təzyiqi (Eclipse `ROCK`).
        # Söndürülü olanda istinad datum təzyiqidir (köhnə davranış).
        self.rock_reference_enabled = QCheckBox(
            "Süxur sıxılmasının istinad təzyiqini ayrıca ver")
        self.rock_reference_enabled.setToolTip(
            "Söndürülü: istinad = ilkin (datum) təzyiq — mövcud modellərdəki "
            "davranış. Açıq: Eclipse ROCK açar sözündəki kimi ayrıca təzyiq.")
        self.rock_reference_pressure = _spin(1.0, 0.01, 2000.0, 3, 1.0, "bar")
        self.rock_reference_enabled.stateChanged.connect(
            self._on_rock_reference_toggled)
        self.rock_reference_enabled.stateChanged.connect(self.changed)
        self.rock_reference_pressure.valueChanged.connect(self.changed)
        form.addRow(self.rock_reference_enabled)
        form.addRow("Süxur istinad təzyiqi", self.rock_reference_pressure)
        self._on_rock_reference_toggled()

        self.context_note = QLabel("")
        self.context_note.setWordWrap(True)
        self.context_note.setStyleSheet("color:#e0a020;font-size:11px")
        self.context_note.setVisible(False)
        form.addRow(self.context_note)

    def geology_values(self) -> dict:
        permx_engine = to_engine_units(self.permx.value(), self.permx_unit.currentText(),
                                       "permeability")
        return dict(porosity=self.porosity.value(), permx_base=permx_engine,
                    ky_over_kx=self.ky_over_kx.value(),
                    kv_over_kh=self.kv_over_kh.value(),
                    heterogeneous=self.heterogeneity.currentIndex() == 1,
                    sigma=self.sigma.value(), seed=self.seed.value())

    def fluids(self) -> FluidProperties:
        unit = self.viscosity_unit.currentText()
        return FluidProperties(
            water_viscosity=to_engine_units(self.mu_w.value(), unit, "viscosity"),
            oil_viscosity=to_engine_units(self.mu_o.value(), unit, "viscosity"),
            oil_fvf=self.bo.value())

    def rock_compressibility_value(self) -> float:
        return self.rock_compressibility.value()

    def _on_rock_reference_toggled(self):
        self.rock_reference_pressure.setEnabled(
            self.rock_reference_enabled.isChecked())

    def rock_compressibility_reference_value(self):
        """İstinad təzyiqi, bar — söndürülübsə `None` (G6)."""
        if not self.rock_reference_enabled.isChecked():
            return None
        return self.rock_reference_pressure.value()

    # ───────────────────────────────────── kontekst (görünürlük düzəlişi)
    #: Geologiya sahələri — GRDECL idxal olunanda bunlar İŞLƏMİR.
    GEOLOGY_WIDGETS = ("porosity", "permx", "permx_unit", "ky_over_kx",
                       "kv_over_kh", "heterogeneity", "sigma", "seed")
    #: Flüid sahələri — PVT modeli işlədiləndə bunlar İŞLƏMİR.
    FLUID_WIDGETS = ("mu_w", "mu_o", "viscosity_unit", "bo")

    def set_context(self, pvt_active: bool = False,
                    geology_imported: bool = False) -> None:
        """Hansı sahələrin FAKTİKİ təsiri olduğunu görünən edir.

        PROBLEM (ölçülüb, `ISH_HESABATI.md` → Seans 14). Bu paneldəki
        bəzi sahələr müəyyən şəraitdə mühərriyə ÜMUMİYYƏTLƏ çatmırdı,
        lakin redaktə edilə bilən qalırdı — istifadəçi dəyəri dəyişir,
        nəticə isə dəyişmir və simulyator sınmış kimi görünür.

          * PVT modeli işlədiləndə μw/μo/Bo `PVT cədvəlindən` gəlir
            (bax `simulation/implicit/residual.py` — `pvt is None`
            budaqlanması). Paneldəki dəyərlər oxunmur.
          * GRDECL modeli idxal olunanda φ/K faylın öz xəritələrindən
            gəlir (bax `main_window._build_geological_model`).

        Bu metod HEÇ BİR hesablamanı dəyişmir — yalnız işləməyən
        sahələri söndürür və səbəbini yazır.
        """
        notes = []
        for name in self.FLUID_WIDGETS:
            getattr(self, name).setEnabled(not pvt_active)
        if pvt_active:
            notes.append(
                "Lözlük və Bo PVT tabından gəlir — buradakı dəyərlər "
                "işlədilmir. Onları dəyişmək üçün PVT tabındakı API, "
                "temperatur və doyma təzyiqini redaktə edin.")

        for name in self.GEOLOGY_WIDGETS:
            getattr(self, name).setEnabled(not geology_imported)
        if geology_imported:
            notes.append(
                "Məsaməlilik və keçiricilik GRDECL faylından gəlir — "
                "buradakı dəyərlər işlədilmir.")

        self.context_note.setText("  ".join(notes))
        self.context_note.setVisible(bool(notes))


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
        # Qaz-neft kapilyar təzyiqi (B5-b) — su-neft ilə eyni üç parametr
        self.pcog_entry = _spin(0.0, 0.0, 20.0, 3, 0.05, "bar")
        self.pcog_lambda = _spin(2.0, 0.2, 10.0, 2, 0.1)
        self.pcog_max = _spin(5.0, 0.1, 100.0, 2, 1.0, "bar")
        self._gas_rows = [("Sgc (bağlı qaz)", self.sgc),
                          ("Sorg (qaza qarşı qalıq neft)", self.sorg),
                          ("krg @ 1-Swc-Sorg", self.krg_end),
                          ("Corey ng", self.ng), ("Corey nog", self.nog),
                          ("Pcog giriş təzyiqi Pe", self.pcog_entry),
                          ("Pcog Brooks-Corey λ", self.pcog_lambda),
                          ("Pcog yuxarı həddi", self.pcog_max)]
        for label, widget in self._gas_rows:
            form.addRow(label, widget)
            widget.valueChanged.connect(self.changed)

        pcog_note = QLabel(
            "Pcog Pe = 0 → qaz-neft kapilyar təzyiqi söndürülür. Pcog yalnız "
            "PVT-də qaz fazası aktiv olanda işlədilir.")
        pcog_note.setWordWrap(True)
        pcog_note.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        form.addRow(pcog_note)

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

    def gas_capillary_values(self) -> CapillaryParameters:
        """Qaz-neft Pc parametrləri (B5-b) — qaz əyriləri söndürülübsə
        söndürülmüş (Pe = 0) dəyər qaytarılır."""
        if not self.gas_enabled.isChecked():
            return CapillaryParameters()
        return CapillaryParameters(entry_pressure=self.pcog_entry.value(),
                                   lambda_exponent=self.pcog_lambda.value(),
                                   max_pressure=self.pcog_max.value())

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

        self.gas_phase_enabled = QCheckBox(
            "Qaz fazasını aktivləşdir (A7 — sınaq statusunda)")
        self.gas_phase_enabled.setChecked(False)
        self.gas_phase_enabled.stateChanged.connect(self.changed)

        # ── İlkin həll olmuş qaz (B4b) ──────────────────────────────
        self.manual_rs = QCheckBox("İlkin Rs-i əl ilə ver")
        self.manual_rs.setChecked(False)
        self.solution_gor = _spin(120.0, 0.0, 600.0, 1, 5.0, "sm³/sm³")
        self.manual_rs.stateChanged.connect(self._on_manual_rs_toggled)
        self.manual_rs.stateChanged.connect(self.changed)
        self.solution_gor.valueChanged.connect(self.changed)

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
        form.addRow(self.manual_rs)
        form.addRow("İlkin Rs", self.solution_gor)
        rs_note = QLabel(
            "Söndürülübsə (tövsiyə olunan), ilkin həll olmuş qaz PVT "
            "cədvəlindən çıxarılır: Rs = Rs_sat(min(P, Pb)) — yəni neft "
            "öz doyma təzyiqinə uyğun qədər qaz saxlayır. Bu dəyər "
            "olmadan neft «ölü» başlayır və təzyiq doyma təzyiqindən "
            "aşağı düşsə belə qaz ayrılmır.")
        rs_note.setWordWrap(True)
        rs_note.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        form.addRow(rs_note)
        self._on_manual_rs_toggled()

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

    def _on_manual_rs_toggled(self):
        self.solution_gor.setEnabled(self.manual_rs.isChecked())

    def initial_solution_gor(self):
        """`None` — PVT cədvəlindən çıxarılsın (defolt)."""
        if not self.manual_rs.isChecked():
            return None
        return self.solution_gor.value()


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
              "Tip", "İdarə", "Qiymət", "rw", "Vurulan faza", "BHP limiti",
              "Debit bazası"]
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
    COL_PHASE = 10
    #: RATE quyusunun BHP həddi, bar (B7 addım 2) — boş = hədd yoxdur
    COL_BHP_LIMIT = 11
    #: RATE hədəfinin həcm bazası (B7 addım 3): LAY həcmi / SƏTH debiti
    COL_RATE_BASIS = 12

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

        layout.addWidget(self._build_tubing_group())

        hint = QLabel("BHP → bar,   RATE → m³/gün (rezervuar həcmi),   "
                      "THP → bar (aşağıdakı lülə qrupu açıq olmalıdır)\n"
                      "BHP limiti → bar, yalnız RATE: istismarçıda minimal, "
                      "vurucuda maksimal BHP (boş = limit yoxdur).\n"
                      "Debit bazası (RATE): LAY = rezervuar həcmi, SƏTH = "
                      "istismarçıda neftin, vurucuda vurulan fazanın səth debiti.\n"
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
                    "phase": ("QAZ" if well.control.injected_phase is Phase.GAS
                              else "SU"),
                    "bhp_limit": well.control.bhp_limit,
                    "basis": ("SƏTH" if well.control.rate_basis is RateBasis.SURFACE
                              else "LAY"),
                    "perf_top": well.perf_top, "perf_bottom": well.perf_bottom})
        finally:
            self.table.blockSignals(False)
        self._recompute_ij_k()

    def _build_tubing_group(self) -> QGroupBox:
        """Lülə həndəsəsi — quyu başı təzyiqi (THP) hesabatı üçün (B4).

        V1-də parametrlər BÜTÜN istismarçılara eyni tətbiq olunur.
        Quyu cədvəli onsuz da 10 sütundur; hər quyuya ayrıca lülə
        vermək VFP cədvəlləri gələndə mənalı olacaq (⏳ B4 v2).
        """
        group = QGroupBox("Lülə — quyu başı təzyiqi (THP)")
        form = QFormLayout(group)

        self.thp_enabled = QCheckBox("THP hesabla")
        self.thp_enabled.setChecked(False)
        self.thp_enabled.setToolTip(
            "Quyu dibi təzyiqindən (BHP) lülə boyunca yuxarı təzyiq "
            "düşgüsü hesablanır və quyu başı təzyiqi tapılır. "
            "V1: yalnız BHP rejimli istismarçılar, tam şaquli lülə, "
            "sürüşmə (slip) nəzərə alınmır.")
        form.addRow("", self.thp_enabled)

        self.tubing_diameter = _spin(62.0, 10.0, 400.0, 1, 1.0, "mm")
        self.tubing_diameter.setToolTip("Borunun DAXİLİ diametri")
        form.addRow("Daxili diametr", self.tubing_diameter)

        self.tubing_roughness = _spin(0.06, 0.0, 5.0, 3, 0.01, "mm")
        self.tubing_roughness.setToolTip(
            "Mütləq kələ-kötürlük ε. Yeni polad boru üçün ~0.06 mm.")
        form.addRow("Kələ-kötürlük", self.tubing_roughness)

        self.wellhead_depth = _spin(0.0, 0.0, 3000.0, 1, 10.0, "m")
        self.wellhead_depth.setToolTip(
            "Quyu başının dərinliyi. Quruda 0, su altında müsbət.")
        form.addRow("Quyu başı dərinliyi", self.wellhead_depth)

        self.tubing_segments = _ispin(20, 1, 500)
        self.tubing_segments.setToolTip(
            "Traverse neçə hissəyə bölünsün. Qazlı quyuda az seqment "
            "kobud nəticə verir — 20 tövsiyə olunur.")
        form.addRow("Seqment sayı", self.tubing_segments)

        for widget in (self.thp_enabled, self.tubing_diameter,
                       self.tubing_roughness, self.wellhead_depth,
                       self.tubing_segments):
            signal = getattr(widget, "stateChanged", None) or widget.valueChanged
            signal.connect(lambda *_: self.changed.emit())
        return group

    def tubing_values(self) -> Optional[TubingGeometry]:
        """Paneldəki lülə həndəsəsi — THP söndürülübsə `None`."""
        if not self.thp_enabled.isChecked():
            return None
        return TubingGeometry(
            diameter=self.tubing_diameter.value() / 1000.0,
            roughness=self.tubing_roughness.value() / 1000.0,
            wellhead_depth=self.wellhead_depth.value(),
            segments=self.tubing_segments.value())

    def values(self) -> List[Well]:
        wells: List[Well] = []
        tubing = self.tubing_values()
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
            phase_widget = self.table.cellWidget(row, self.COL_PHASE)
            phase = (Phase.GAS
                     if phase_widget is not None
                     and phase_widget.currentText() == "QAZ" else Phase.WATER)
            limit_item = self.table.item(row, self.COL_BHP_LIMIT)
            bhp_limit = (self._to_float(limit_item.text())
                         if limit_item is not None else None)
            basis_widget = self.table.cellWidget(row, self.COL_RATE_BASIS)
            rate_basis = (RateBasis.SURFACE
                          if basis_widget is not None
                          and basis_widget.currentText() == "SƏTH"
                          else RateBasis.RESERVOIR)

            i, j, first, last = self._resolve_ijk(name, perf_top, perf_bottom)
            wells.append(Well(
                name=name, well_type=WellType(kind),
                control=WellControl(ControlMode(mode), target,
                                    injected_phase=phase, bhp_limit=bhp_limit,
                                    rate_basis=rate_basis),
                perforations=[Perforation(i, j, k) for k in range(first, last + 1)],
                radius=rw, perf_top=perf_top, perf_bottom=perf_bottom,
                tubing=(tubing if WellType(kind) is WellType.PRODUCER else None)))
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
                "phase": "SU", "bhp_limit": None, "basis": "LAY",
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
            "phase": (self.table.cellWidget(row, self.COL_PHASE).currentText()
                      if self.table.cellWidget(row, self.COL_PHASE) else "SU"),
            "bhp_limit": (self._to_float(self.table.item(row, self.COL_BHP_LIMIT).text())
                          if self.table.item(row, self.COL_BHP_LIMIT) else None),
            "basis": (self.table.cellWidget(row, self.COL_RATE_BASIS).currentText()
                      if self.table.cellWidget(row, self.COL_RATE_BASIS) else "LAY"),
            "perf_top":self._to_float(self.table.item(row, self.COL_PERF_TOP).text()),
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
        mode_box.addItems(["BHP", "RATE", "THP"])
        mode_box.setCurrentText(state.get("mode", "BHP"))
        mode_box.currentIndexChanged.connect(lambda *_: self.changed.emit())
        self.table.setCellWidget(r, self.COL_MODE, mode_box)

        # Vurulan faza — YALNIZ vurucu quyuda mənalıdır (B7)
        phase_box = QComboBox()
        phase_box.addItems(["SU", "QAZ"])
        phase_box.setCurrentText(state.get("phase", "SU"))
        phase_box.currentIndexChanged.connect(lambda *_: self.changed.emit())
        self.table.setCellWidget(r, self.COL_PHASE, phase_box)

        self.table.setItem(r, self.COL_TARGET,
                           QTableWidgetItem(f"{state.get('target') or 150.0:g}"))
        self.table.setItem(r, self.COL_RW,
                           QTableWidgetItem(f"{state.get('rw') or 0.1:g}"))
        limit = state.get("bhp_limit")
        self.table.setItem(r, self.COL_BHP_LIMIT, QTableWidgetItem(
            "" if limit is None else f"{limit:g}"))

        # Debit bazası — YALNIZ RATE rejimində mənalıdır (B7 addım 3)
        basis_box = QComboBox()
        basis_box.addItems(["LAY", "SƏTH"])
        basis_box.setCurrentText(state.get("basis", "LAY"))
        basis_box.currentIndexChanged.connect(lambda *_: self.changed.emit())
        self.table.setCellWidget(r, self.COL_RATE_BASIS, basis_box)

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
        self.initial_pressure = _spin(250, 1, 1200, 3, 10, "bar")
        self.initial_pressure_unit = _unit_combo(["bar", "psi", "MPa"], "bar")
        _bind_unit_aware_spins([self.initial_pressure], self.initial_pressure_unit, "pressure")
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
        self.flux_scheme = QComboBox()
        self.flux_scheme.addItem("TPFA (iki-nöqtəli)", TPFA)
        self.flux_scheme.addItem("MPFA-O (çoxnöqtəli, tam tenzor)", MPFA_O)
        self.flux_scheme.currentIndexChanged.connect(self.changed)
        form.addRow("Axın diskretizasiyası", self.flux_scheme)
        flux_note = QLabel(
            "TPFA — defolt, ortoqonal gridd və diaqonal keçiricilikdə "
            "kifayətdir. MPFA-O — anizotrop tenzor (Kxy≠0) və "
            "qeyri-ortoqonal corner-point gridd axını düzgün verir; "
            "YALNIZ «Fully implicit» mühərriklə işləyir.")
        flux_note.setWordWrap(True)
        flux_note.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        form.addRow(flux_note)
        self.use_equilibration = QCheckBox("Equilibration (dərinlikdən asılı ilkin şərtlər)")
        self.use_saturation_map = QCheckBox("İlkin Sw geologiya xəritəsindən (SW)")
        self.saturation_map_info = QLabel()
        self.saturation_map_info.setStyleSheet(
            f"color:{PALETTE.text_dim};font-size:11px")
        self.datum_depth = _spin(2000.0, 0.0, 8000.0, 1, 50.0, "m")
        self.owc = _spin(2050.0, 0.0, 8000.0, 1, 10.0, "m")
        self.use_goc = QCheckBox(
            "Qaz papağı (GOC) — yalnız qaz fazası aktivdirsə təsir edir")
        self.goc = _spin(2010.0, 0.0, 8000.0, 1, 10.0, "m")
        form.addRow(self.use_equilibration)
        form.addRow(self.use_saturation_map)
        form.addRow(self.saturation_map_info)
        for label, widget in [("Başlanğıc təzyiq", self.initial_pressure),
                              ("Başlanğıc təzyiq vahidi", self.initial_pressure_unit),
                              ("Başlanğıc Sw", self.initial_sw),
                              ("Datum dərinliyi", self.datum_depth),
                              ("Su-neft kontaktı (OWC)", self.owc),
                              ("Simulyasiya müddəti", self.end_time),
                              ("Maks. Δt", self.max_dt),
                              ("CFL əmsalı", self.cfl),
                              ("Yaddaş anlarının sayı", self.snapshots)]:
            form.addRow(label, widget)
            sig = getattr(widget, "valueChanged", None) or widget.currentIndexChanged
            sig.connect(self.changed)
        form.addRow(self.use_goc)
        form.addRow("Qaz-neft kontaktı (GOC)", self.goc)
        self.use_goc.stateChanged.connect(self.changed)
        self.goc.valueChanged.connect(self.changed)
        self.use_equilibration.stateChanged.connect(self.changed)
        self.use_saturation_map.stateChanged.connect(self._saturation_source_changed)
        self.set_saturation_map(None)
        note = QLabel("Söndürülübsə, bütün hüceyrələrdə eyni təzyiq və Sw "
                      "işlədilir (köhnə davranış).")
        note.setWordWrap(True)
        note.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        form.addRow(note)

    def _saturation_source_changed(self):
        """Xəritə rejimində skalyar sahə İŞLƏMİR — söndürülür ki,
        istifadəçi nəticəyə təsir etməyən qutunu doldurmasın."""
        self.initial_sw.setEnabled(not self.use_saturation_map.isChecked())
        self.changed.emit()

    def set_saturation_map(self, stats: Optional[dict]) -> None:
        """`SW` xəritəsinin mövcudluğunu bildirir (`PropertyMap.stats()`).

        Xəritə yoxdursa seçim SEÇİLƏ BİLMİR — amma ARTIQ seçilibsə
        (məs. `.imx` faylından belə açılıb) söndürülmür: əks halda
        istifadəçi qutunu geri qaldıra bilməzdi. Bu halda vəziyyət
        etiketdə AÇIQ yazılır və model diaqnostikası xəta verir.
        """
        available = stats is not None
        self.use_saturation_map.setEnabled(
            available or self.use_saturation_map.isChecked())
        if available:
            self.saturation_map_info.setText(
                f"SW: min {stats['min']:.3f} · orta {stats['mean']:.3f} · "
                f"maks {stats['max']:.3f}")
            self.use_saturation_map.setToolTip("")
        else:
            self.saturation_map_info.setText("SW xəritəsi yoxdur.")
            self.use_saturation_map.setToolTip(
                "Əvvəlcə geologiya cədvəlini interpolyasiya edin.")

    def initial_conditions(self) -> InitialConditions:
        equilibrate = self.use_equilibration.isChecked()
        pressure_bar = to_engine_units(self.initial_pressure.value(),
                                       self.initial_pressure_unit.currentText(), "pressure")
        return InitialConditions(
            datum_depth=self.datum_depth.value(),
            datum_pressure=pressure_bar,
            water_saturation=self.initial_sw.value(),
            oil_water_contact=self.owc.value() if equilibrate else None,
            # GOC yalnız equilibration açıq olanda mənalıdır — hidrostatik
            # tarazlıq qurulmadan kontakt dərinliyinin təsiri yoxdur.
            gas_oil_contact=(self.goc.value()
                             if equilibrate and self.use_goc.isChecked()
                             else None),
            solution_gor=getattr(self, "_solution_gor", None),
            use_equilibration=equilibrate,
            use_saturation_map=self.use_saturation_map.isChecked())

    def engine_choice(self) -> str:
        return self.engine.currentData()

    def flux_scheme_choice(self) -> str:
        return self.flux_scheme.currentData()

    def set_solution_gor(self, value) -> None:
        """PVT panelindən gələn ilkin Rs (B4b) — `None` = PVT-dən çıxar.

        Dəyər `InitialConditions`-a burada qoşulur, çünki `initial_conditions()`
        məhz bu paneldədir; mənbəyi isə PVT panelidir (`MainWindow` ötürür).
        """
        self._solution_gor = value

    def set_flux_scheme(self, scheme: str) -> None:
        """`.imx` faylından və ya proqramlı olaraq sxemi seçir.

        Naməlum dəyər gələndə TPFA-ya qayıdır — köhnə fayl açanda
        istifadəçi qarşısına xəta çıxmasın (dəyər onsuz da
        `SimulationConfig.validate()`-də yoxlanılır).
        """
        index = self.flux_scheme.findData(scheme)
        self.flux_scheme.setCurrentIndex(index if index >= 0 else 0)

    def simulation_config(self) -> SimulationConfig:
        return SimulationConfig(
            end_time=self.end_time.value(),
            time_stepping=TimeSteppingConfig(max_dt=self.max_dt.value(),
                                             cfl_factor=self.cfl.value()),
            linear_solver=LinearSolverConfig(),
            output=OutputConfig(snapshot_count=self.snapshots.value()),
            flux_scheme=self.flux_scheme_choice(),
        )
