"""Əsas pəncərə — YALNIZ orkestrasiya.

Bu fayl heç bir fiziki hesablama aparmır və heç bir massiv qurmur.
İşi: panellərdən domain obyektlərini toplamaq, application servislərini
çağırmaq, nəticəni rendering qatına ötürmək.
"""

from __future__ import annotations
import csv
import logging
import os
from typing import Optional

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavToolbar
from matplotlib.figure import Figure

from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import (QAction, QCheckBox, QComboBox, QDoubleSpinBox,
                             QFileDialog,
                             QHBoxLayout, QHeaderView, QLabel, QTableWidget,
                             QTableWidgetItem,
                             QMainWindow, QMessageBox, QProgressBar, QPushButton,
                             QSlider, QSpinBox, QSplitter, QStackedWidget,
                             QTabWidget, QTextEdit,
                             QToolBox, QTreeWidget, QTreeWidgetItem, QVBoxLayout,
                             QWidget)

from ..application.config import (OutputConfig, SimulationConfig,
                                  TimeSteppingConfig)
from ..application.model_builder import ReservoirModelBuilder
from ..application.project import Project
from ..application.serialization import (FILE_EXTENSION, ProjectFileError,
                                         ProjectSerializer)
from ..application.geology_adapter import wells_to_dataset
from ..application.geology_service import (GeologicalGridSpec,
                                          WellBasedGeologicalModelBuilder,
                                          format_cross_validation_report)
from ..application.scenarios import WELL_PATTERNS, SyntheticGeologicalModelBuilder
from ..application.simulation_service import (ModelValidationError,
                                              SimulationService)
from ..domain.data_availability import STATUS_LABEL, DataStatus
from ..domain.diagnostics import DiagnosticReport, Severity
from ..domain.geology import GeologicalWell, validate_wells
from ..domain.grid import CartesianGrid
from ..domain.geometry import CellGeometry
from ..domain.reservoir_model import ReservoirModel
from ..domain.wells import ControlMode, Perforation, Well, WellControl, WellType
from ..rendering import renderers as R
from ..rendering.volume import (VIEW_ANGLES, VolumeFilter, VolumeRenderer,
                                apply_zoom)

from ..rendering.theme import PALETTE
from ..history.mismatch import MismatchCalculator
from ..history.optimizer import HistoryMatchingService
from ..history.sensitivity import OUTPUT_METRICS, SensitivityAnalyzer
from ..history.parameters import ParameterSet, standard_parameters
from ..history.observation_io import (ObservationFormatError,
                                       read_observations_csv)
from ..io.eclipse_export import EclipseDeckWriter
from ..reporting.report import ReportContext, ReportGenerator
from ..io.grdecl import GrdeclError, read_grdecl
from ..io.grdecl_import import GrdeclImporter
from ..logging_setup import add_handler, get_logger
from ..version import EXPECTED_TABS, VERSION, summary, title
from ..simulation.analytical import buckley_leverett
from ..simulation.impes_engine import ImpesEngine
from ..simulation.implicit.engine import FullyImplicitEngine
from .panels import (FaultPanel, GeologyPanel, GridGeometryPanel,
                     NumericalPanel, PvtPanel, RockFluidPanel, ScalPanel,
                     ScalSourcePanel, WellPanel)
from .worker import MatchingWorker, SensitivityWorker, SimulationWorker


LOG = get_logger(__name__)


class QtLogHandler(logging.Handler):
    """Log mesajlarını jurnal tabına yönləndirir.

    Beləliklə istifadəçinin gördüyü mətn ilə fayldakı log eyni olur.
    """

    def __init__(self, sink):
        super().__init__()
        self._sink = sink

    def emit(self, record):
        try:
            self._sink(self.format(record))
        except Exception:
            pass


def _figure(nrows=1, ncols=1):
    fig = Figure(facecolor=PALETTE.background, tight_layout=True)
    canvas = FigureCanvas(fig)
    return fig, canvas, fig.subplots(nrows, ncols)


class MainWindow(QMainWindow):

    def __init__(self, project: Project, service: SimulationService,
                 geology_builder: SyntheticGeologicalModelBuilder,
                 model_builder: ReservoirModelBuilder):
        super().__init__()
        self.project = project
        self.service = service
        self.geology_builder = geology_builder
        self.model_builder = model_builder

        self.reservoir_model: Optional[ReservoirModel] = None
        self.imported_geology = None
        self._geology_model_from_wells = None
        self._dirty = False
        self.result = None
        self.worker: Optional[SimulationWorker] = None
        self.colorbar = None

        self.map_renderer = R.MapRenderer()
        self.section_renderer = R.CrossSectionRenderer()
        self.volume_renderer = VolumeRenderer()
        self.volume_colorbar = None
        self.curve_renderer = R.ProductionCurveRenderer()
        self.scal_renderer = R.ScalRenderer()
        self.pvt_renderer = R.PvtRenderer()
        self.comparison_renderer = R.RunComparisonRenderer()
        self.history_renderer = R.HistoryMatchRenderer()
        self.crossplot_renderer = R.CrossPlotRenderer()
        self.optimisation_renderer = R.OptimisationRenderer()
        self.tornado_renderer = R.TornadoRenderer()
        self.sensitivity_worker = None
        self.sensitivity_report = None
        self.match_worker = None
        self.match_result = None
        self.observations = None
        self.mismatch_report = None
        self.serializer = ProjectSerializer()
        self.project_path = None
        self.validation_renderer = R.ValidationRenderer()

        self._player = QTimer(self)
        self._player.timeout.connect(self._next_frame)

        self.resize(1560, 940)
        self._build_menu()
        self._build_ui()
        self._refresh_title()
        self._ready = True
        default_wells = WELL_PATTERNS["Five-spot (1/4)"](self.grid_panel.grid())
        self.well_panel.load(default_wells)
        self.geology_panel.load(
            self._wells_to_geology_rows(default_wells, self._current_geometry()))
        self._sync_geology_geometry()
        self.geology_panel.mark_fresh()
        self._mark_clean()
        self.rebuild_model()
        self._verify_build()

    # ---------------------------------------------------- pəncərə başlığı
    def _refresh_title(self):
        """`*` — geologiya cədvəlində yadda saxlanılmamış dəyişiklik var."""
        self.setWindowTitle(title() + (" *" if self._dirty else ""))

    def _mark_dirty(self):
        self._dirty = True
        self._refresh_title()

    def _mark_clean(self):
        self._dirty = False
        self._refresh_title()

    # ------------------------------------------------------------ menyu
    def _build_menu(self):
        file_menu = self.menuBar().addMenu("&Layihə")
        for text, slot in [("Layihəni aç…", self.open_project),
                           ("Layihəni yadda saxla…", self.save_project),
                           ("Layihəni saxla (nəticəsiz)…",
                            lambda: self.save_project(include_snapshots=False))]:
            action = QAction(text, self)
            action.triggered.connect(slot)
            file_menu.addAction(action)
        file_menu.addSeparator()
        for text, slot in [("GRDECL grid oxu…", self.import_grdecl),
                           ("Eclipse deck yaz (.DATA)…", self.export_eclipse),
                           ("PDF hesabat yaz…", self.export_pdf_report),
                           ("OPM Flow nəticəsini yüklə (.EGRID+.UNRST)…",
                            self.import_opm_case)]:
            action = QAction(text, self)
            action.triggered.connect(slot)
            file_menu.addAction(action)
        file_menu.addSeparator()
        for text, slot in [("Nəticələri CSV kimi yaz…", self.export_results),
                           ("Grid anını CSV kimi yaz…", self.export_snapshot)]:
            action = QAction(text, self)
            action.triggered.connect(slot)
            file_menu.addAction(action)
        file_menu.addSeparator()
        quit_action = QAction("Çıxış", self)
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)

        help_menu = self.menuBar().addMenu("&Kömək")
        version_action = QAction(f"Versiya {VERSION} …", self)
        version_action.triggered.connect(self.show_version)
        help_menu.addAction(version_action)
        about = QAction("Arxitektura haqqında", self)
        about.triggered.connect(self.show_about)
        help_menu.addAction(about)

    # --------------------------------------------------------------- ui
    def _build_ui(self):
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self._build_input_panel())
        splitter.addWidget(self._build_view_panel())
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([430, 1130])
        self.setCentralWidget(splitter)

        self.progress = QProgressBar()
        self.progress.setMaximumWidth(280)
        self.statusBar().addPermanentWidget(self.progress)
        self.statusBar().showMessage("Hazır.")

    def _build_input_panel(self):
        self.toolbox = QToolBox()
        self.toolbox.setMinimumWidth(400)

        self.grid_panel = GridGeometryPanel()
        self.geology_panel = GeologyPanel()
        self.rock_panel = RockFluidPanel()
        self.scal_panel = ScalPanel()
        self.scal_source_panel = ScalSourcePanel()
        self.fault_panel = FaultPanel()
        self.pvt_panel = PvtPanel()
        self.well_panel = WellPanel()
        self.numerical_panel = NumericalPanel()

        self.toolbox.addItem(self.grid_panel, "1 · GRID (geoloji model)")
        self.toolbox.addItem(self.geology_panel, "2 · GEOLOJİ MƏLUMAT (quyular)")
        self.toolbox.addItem(self.rock_panel, "3 · SÜXUR & FLÜİD")
        self.toolbox.addItem(self.fault_panel, "4 · FAULTS")
        scal_page = QWidget()
        scal_layout = QVBoxLayout(scal_page)
        scal_layout.setContentsMargins(0, 0, 0, 0)
        scal_layout.addWidget(self.scal_source_panel)
        scal_layout.addWidget(self.scal_panel, 1)
        self.toolbox.addItem(scal_page, "5 · NİSBİ KEÇİRİCİLİK")
        self.toolbox.addItem(self.pvt_panel, "6 · PVT (flüid modeli)")
        self.toolbox.addItem(self.well_panel, "7 · QUYULAR (rezervuar modeli)")
        self.toolbox.addItem(self.numerical_panel, "8 · ƏDƏDİ PARAMETRLƏR")

        # `grid_panel.changed`: `_sync_geology_geometry` BİRİNCİ qoşulur ki,
        # (i,j,k) təzələnsin, SONRA `rebuild_model` köhnəlməmiş indekslərlə
        # işləsin (bağlantı sırası = çağırılma sırası).
        self.grid_panel.changed.connect(self._sync_geology_geometry)
        # `geology_panel` bilərəkdən aşağıdakı dövrədə YOXDUR: cədvəl
        # dəyişəndə model avtomatik yenidən qurulmur (böyük gridə yavaşdır)
        # — yalnız "İnterpolyasiya et" düyməsi (`_interpolate_geology`) qurur.
        for panel in (self.grid_panel, self.rock_panel,
                      self.fault_panel, self.scal_panel, self.scal_source_panel,
                      self.pvt_panel, self.well_panel, self.numerical_panel):
            panel.changed.connect(self.rebuild_model)
        self.geology_panel.changed.connect(self._on_geology_table_changed)
        self.geology_panel.interpolate_requested.connect(self._interpolate_geology)
        self.geology_panel.cross_validate_requested.connect(self._cross_validate_geology)
        self.well_panel.apply_button.clicked.connect(self._apply_pattern)

        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(self.toolbox, 1)

        self.run_button = QPushButton("MODELİ İŞƏ SAL")
        self.run_button.setObjectName("run")
        self.run_button.clicked.connect(self.run_simulation)
        self.stop_button = QPushButton("Dayandır")
        self.stop_button.setObjectName("stop")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self.stop_simulation)
        layout.addWidget(self.run_button)
        layout.addWidget(self.stop_button)
        return container

    def _build_view_panel(self):
        self.tabs = QTabWidget()

        # Layihə ağacı
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Layihə strukturu"])
        self.tabs.addTab(self.tree, "Layihə")

        # Model xəritəsi
        page = QWidget(); layout = QVBoxLayout(page)
        bar = QHBoxLayout()
        self.property_box = QComboBox()
        for key in (R.SATURATION, R.PRESSURE, R.PERMEABILITY, R.POROSITY,
                    R.DEPTH):
            self.property_box.addItem(R.PROPERTY_LABELS[key], key)
        self.property_box.currentIndexChanged.connect(self.update_map)
        self.view_mode = QComboBox()
        self.view_mode.addItem("Areal xəritə (təbəqə)", "AREAL")
        self.view_mode.addItem("Şaquli kəsik — X üzrə", "SECTION_J")
        self.view_mode.addItem("Şaquli kəsik — Y üzrə", "SECTION_I")
        self.view_mode.currentIndexChanged.connect(self.update_map)
        self.layer_spin = QSpinBox()
        self.layer_spin.setRange(1, 1)
        self.layer_spin.setPrefix("K = ")
        self.layer_spin.valueChanged.connect(self.update_map)
        self.play_button = QPushButton("▶  Oynat")
        self.play_button.clicked.connect(self.toggle_play)
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setEnabled(False)
        self.slider.valueChanged.connect(self.update_map)
        self.time_label = QLabel("t = 0 gün")
        self.time_label.setMinimumWidth(150)
        self.time_label.setStyleSheet(f"color:{PALETTE.accent};font-family:monospace")
        bar.addWidget(QLabel("Xassə:")); bar.addWidget(self.property_box, 1)
        bar.addWidget(self.view_mode, 1); bar.addWidget(self.layer_spin)
        bar.addWidget(self.play_button); bar.addWidget(self.slider, 2)
        bar.addWidget(self.time_label)
        layout.addLayout(bar)
        self.map_fig = Figure(facecolor=PALETTE.background)
        self.map_canvas = FigureCanvas(self.map_fig)
        gridspec = self.map_fig.add_gridspec(1, 2, width_ratios=[40, 1],
                                             wspace=0.04, left=0.08, right=0.94,
                                             top=0.92, bottom=0.10)
        self.map_ax = self.map_fig.add_subplot(gridspec[0, 0])
        self.map_cax = self.map_fig.add_subplot(gridspec[0, 1])
        layout.addWidget(NavToolbar(self.map_canvas, self))
        layout.addWidget(self.map_canvas, 1)
        self.tabs.addTab(page, "Model")

        # Nəticələr
        page = QWidget(); layout = QVBoxLayout(page)
        self.kpi = QLabel("Model hələ işə salınmayıb.")
        self.kpi.setStyleSheet(
            f"background:{PALETTE.panel_alt};border:1px solid {PALETTE.line};"
            f"border-radius:4px;padding:9px;font-family:monospace;"
            f"font-size:12px;color:{PALETTE.text}")
        self.result_fig, self.result_canvas, self.result_axes = _figure(2, 2)
        layout.addWidget(self.kpi)
        layout.addWidget(self.result_canvas, 1)
        self.tabs.addTab(page, "Nəticələr")

        # SCAL
        page = QWidget(); layout = QVBoxLayout(page)
        self.scal_fig, self.scal_canvas, scal_axes = _figure(1, 2)
        self.scal_axes = scal_axes[:2]
        layout.addWidget(self.scal_canvas)
        self.tabs.addTab(page, "Nisbi keçiricilik")

        # 3D görüntü
        self.tabs.addTab(self._build_volume_tab(), "3D görüntü")

        # PVT
        page = QWidget(); layout = QVBoxLayout(page)
        self.pvt_info = QLabel("PVT söndürülüb — sabit flüid xassələri işlədilir.")
        self.pvt_info.setStyleSheet(f"color:{PALETTE.text_dim}")
        self.pvt_fig, self.pvt_canvas, self.pvt_axes = _figure(2, 2)
        layout.addWidget(self.pvt_info)
        layout.addWidget(self.pvt_canvas, 1)
        self.tabs.addTab(page, "PVT")

        # Validasiya
        page = QWidget(); layout = QVBoxLayout(page)
        top = QHBoxLayout()
        button = QPushButton("1D Bukley-Leverett testini işə sal")
        button.clicked.connect(self.run_validation)
        self.validation_label = QLabel("Ədədi həlli analitik həll ilə müqayisə edir.")
        self.validation_label.setStyleSheet(f"color:{PALETTE.text_dim}")
        top.addWidget(button); top.addWidget(self.validation_label, 1)
        layout.addLayout(top)
        self.validation_fig, self.validation_canvas, self.validation_ax = _figure()
        layout.addWidget(self.validation_canvas, 1)
        self.tabs.addTab(page, "Validasiya (B-L)")

        # --- Müqayisə
        page = QWidget(); layout = QVBoxLayout(page)
        top = QHBoxLayout()
        refresh = QPushButton("Yenilə")
        refresh.clicked.connect(self.update_comparison)
        clear = QPushButton("İşə salınmaları təmizlə")
        clear.clicked.connect(self.clear_runs)
        self.comparison_info = QLabel("Ən azı iki dəfə modeli işə salın.")
        self.comparison_info.setStyleSheet(f"color:{PALETTE.text_dim}")
        top.addWidget(refresh); top.addWidget(clear)
        top.addWidget(self.comparison_info, 1)
        layout.addLayout(top)

        self.comparison_table = QTableWidget(0, 7)
        self.comparison_table.setHorizontalHeaderLabels(
            ["Run", "Model", "Müddət, gün", "RF, %", "Su gəlişi, gün",
             "Kum. neft, min m³", "Son WC, %"])
        self.comparison_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch)
        self.comparison_table.verticalHeader().setVisible(False)
        self.comparison_table.setMaximumHeight(190)
        layout.addWidget(self.comparison_table)

        self.comparison_fig, self.comparison_canvas, self.comparison_axes = _figure(2, 2)
        layout.addWidget(self.comparison_canvas, 1)
        self.tabs.addTab(page, "Müqayisə")

        # --- Tarixçə uyğunluğu
        self.tabs.addTab(self._build_history_tab(), "Tarixçə")

        # --- Avtomatik uyğunlaşdırma
        self.tabs.addTab(self._build_matching_tab(), "Uyğunlaşdırma")

        # --- Həssaslıq analizi
        self.tabs.addTab(self._build_sensitivity_tab(), "Həssaslıq")

        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setFont(QFont("monospace", 9))
        self.tabs.addTab(self.log, "Jurnal")

        self.log_handler = QtLogHandler(self._append_log)
        add_handler(self.log_handler)
        return self.tabs

    def _build_volume_tab(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)

        top = QHBoxLayout()

        # ── motor seçimi (VTK / matplotlib) ─────────────────────────
        # VTK MƏCBURİ ASILILIQ DEYİL: quraşdırılmayıbsa siyahıda
        # görünmür və proqram tam olduğu kimi (matplotlib ilə) işləyir.
        from ..rendering import vtk_volume
        self._vtk_available = vtk_volume.available()
        self.volume_engine = QComboBox()
        if self._vtk_available:
            self.volume_engine.addItem("VTK (sürətli)", "vtk")
        self.volume_engine.addItem("matplotlib", "matplotlib")
        self.volume_engine.setToolTip(
            "VTK — ResInsight tipli sürətli 3D (rəvan fırlatma, real "
            "işıqlandırma). matplotlib — köhnə, hər yerdə işləyən motor."
            if self._vtk_available else
            "VTK quraşdırılmayıb (pip install vtk) — yalnız matplotlib "
            "mövcuddur.")
        self.volume_engine.currentIndexChanged.connect(self._on_engine_changed)

        self.volume_property = QComboBox()
        for key in (R.SATURATION, R.PRESSURE, R.PERMEABILITY, R.POROSITY,
                    R.DEPTH):
            self.volume_property.addItem(R.PROPERTY_LABELS[key], key)
        self.volume_property.currentIndexChanged.connect(self.update_volume)

        self.volume_time = QSlider(Qt.Horizontal)
        self.volume_time.setEnabled(False)
        self.volume_time.valueChanged.connect(self.update_volume)
        self.volume_time_label = QLabel("t = 0 gün")
        self.volume_time_label.setMinimumWidth(130)
        self.volume_time_label.setStyleSheet(
            f"color:{PALETTE.accent};font-family:monospace")

        top.addWidget(QLabel("Motor:"))
        top.addWidget(self.volume_engine)
        top.addWidget(QLabel("Xassə:"))
        top.addWidget(self.volume_property, 1)
        top.addWidget(self.volume_time, 2)
        top.addWidget(self.volume_time_label)
        layout.addLayout(top)

        controls = QHBoxLayout()
        self.volume_threshold = QSlider(Qt.Horizontal)
        self.volume_threshold.setRange(0, 100)
        self.volume_threshold.setValue(0)
        self.volume_threshold.valueChanged.connect(self.update_volume)
        self.volume_threshold_label = QLabel("hamısı")
        self.volume_threshold_label.setMinimumWidth(110)
        self.volume_threshold_label.setStyleSheet(
            f"color:{PALETTE.text_dim};font-size:11px")

        self.volume_k_from = QSpinBox(); self.volume_k_from.setPrefix("K ")
        self.volume_k_to = QSpinBox(); self.volume_k_to.setPrefix("… ")
        for widget in (self.volume_k_from, self.volume_k_to):
            widget.setRange(1, 1)
            widget.valueChanged.connect(self.update_volume)

        self.volume_zoom = QSpinBox()
        self.volume_zoom.setRange(25, 500)
        self.volume_zoom.setSingleStep(25)
        self.volume_zoom.setValue(100)
        self.volume_zoom.setSuffix(" %")
        self.volume_zoom.setToolTip(
            "Ekranı yaxınlaşdırma. 100 % = tam model çərçivəyə sığır; "
            "artırdıqca modelə yaxından baxılır.")
        self.volume_zoom.valueChanged.connect(self.update_volume)

        self.volume_wells = QCheckBox("Quyular")
        self.volume_wells.setChecked(True)
        self.volume_wells.stateChanged.connect(self.update_volume)

        self.volume_faults = QCheckBox("Faults")
        self.volume_faults.setChecked(True)
        self.volume_faults.stateChanged.connect(self.update_volume)

        self.volume_edges = QCheckBox("Kənarlar")
        self.volume_edges.setChecked(True)
        self.volume_edges.stateChanged.connect(self.update_volume)

        # ── status filtri (lay-məlumatlı rejim, §13) ─────────────────
        # Siyahı `DataStatus`-dan qurulur — burada heç bir status adı
        # SABİT KODLANMIR. Model provenance daşımırsa filtr sadəcə
        # təsirsizdir (heç nə gizlədilmir).
        self.volume_status_property = QComboBox()
        self.volume_status_property.addItem("— status filtri yoxdur —", "")
        self.volume_status_property.setToolTip(
            "Hansı xassənin MƏNŞƏ statusuna görə süzüləcəyi.")
        self.volume_status_property.currentIndexChanged.connect(self.update_volume)

        self.volume_status = QComboBox()
        self.volume_status.addItem("Hamısı", "")
        for status in DataStatus:
            self.volume_status.addItem(STATUS_LABEL[status.value], status.value)
        self.volume_status.setToolTip(
            "Yalnız seçilmiş mənşəli hüceyrələr göstərilir — ölçülmüş, "
            "interpolyasiya olunmuş, qiymətləndirilmiş, simulyasiya edilmiş "
            "və ya məlumatsız (MISSING) hüceyrələr ayrıca görünə bilər.")
        self.volume_status.currentIndexChanged.connect(self.update_volume)

        controls.addWidget(QLabel("Kəsim həddi:"))
        controls.addWidget(self.volume_threshold, 2)
        controls.addWidget(self.volume_threshold_label)
        controls.addWidget(self.volume_k_from)
        controls.addWidget(self.volume_k_to)
        controls.addWidget(QLabel("Yaxınlaşdırma:"))
        controls.addWidget(self.volume_zoom)
        controls.addWidget(self.volume_wells)
        controls.addWidget(self.volume_faults)
        controls.addWidget(self.volume_edges)
        controls.addWidget(QLabel("Status:"))
        controls.addWidget(self.volume_status_property)
        controls.addWidget(self.volume_status)
        layout.addLayout(controls)

        appearance = QHBoxLayout()
        self.volume_view = QComboBox()
        self.volume_view.addItem("Sərbəst", None)
        for name in VIEW_ANGLES:
            self.volume_view.addItem(name, name)
        self.volume_view.currentIndexChanged.connect(self.update_volume)

        self.volume_shading = QSlider(Qt.Horizontal)
        self.volume_shading.setRange(0, 100)
        self.volume_shading.setValue(45)
        self.volume_shading.setMaximumWidth(140)
        self.volume_shading.valueChanged.connect(self.update_volume)

        self.volume_opacity = QSlider(Qt.Horizontal)
        self.volume_opacity.setRange(20, 100)
        self.volume_opacity.setValue(100)
        self.volume_opacity.setMaximumWidth(140)
        self.volume_opacity.valueChanged.connect(self.update_volume)

        self.volume_reset = QPushButton("Görünüşü sıfırla")
        self.volume_reset.clicked.connect(self.reset_volume_view)

        appearance.addWidget(QLabel("Baxış:"))
        appearance.addWidget(self.volume_view)
        appearance.addWidget(QLabel("İşıq:"))
        appearance.addWidget(self.volume_shading)
        appearance.addWidget(QLabel("Şəffaflıq:"))
        appearance.addWidget(self.volume_opacity)
        appearance.addWidget(self.volume_reset)
        appearance.addStretch()
        layout.addLayout(appearance)

        self.volume_fig = Figure(facecolor=PALETTE.background)
        self.volume_canvas = FigureCanvas(self.volume_fig)
        spec = self.volume_fig.add_gridspec(1, 2, width_ratios=[40, 1],
                                            wspace=0.02, left=0.02, right=0.93,
                                            top=0.95, bottom=0.05)
        self.volume_ax = self.volume_fig.add_subplot(spec[0, 0],
                                                     projection="3d")
        self.volume_cax = self.volume_fig.add_subplot(spec[0, 1])
        # matplotlib alət paneli 3D-də ZƏRƏRLİDİR: "pan"/"zoom" rejimi
        # aktiv olanda `Axes3D._on_move` fırlatmanı tamamilə bloklayır
        # (`get_navigate_mode() is not None -> return`). Ona görə burada
        # yalnız şəkli saxlamaq düyməsi qalır.
        toolbar = QHBoxLayout()
        save_button = QPushButton("Şəkli saxla…")
        save_button.clicked.connect(self.save_volume_image)
        toolbar.addWidget(save_button)
        toolbar.addStretch()
        layout.addLayout(toolbar)

        # ── iki motor eyni yerdə, biri görünür ──────────────────────
        # `QStackedWidget` seçilib (hər dəfə widget yaratmaq/silmək
        # əvəzinə): VTK səhnəsinin qurulması bahalıdır, motorlar arası
        # keçid ANİ olmalıdır. Hər ikisi yaddaşda qalır.
        self.volume_stack = QStackedWidget()
        self.volume_stack.addWidget(self.volume_canvas)

        self.vtk_widget = None
        self.vtk_scene = None
        if self._vtk_available:
            try:
                from vtkmodules.qt.QVTKRenderWindowInteractor import (
                    QVTKRenderWindowInteractor)
                self.vtk_widget = QVTKRenderWindowInteractor(page)
                self.volume_stack.addWidget(self.vtk_widget)
            except Exception:
                # VTK modulu var, lakin Qt körpüsü qurula bilmədi
                # (OpenGL sürücüsü, uzaq masaüstü və s.) — səssizcə
                # matplotlib-ə qayıdırıq, proqram işləməyə davam edir.
                LOG.exception("VTK Qt widget qurulmadı — matplotlib işlədilir")
                self._vtk_available = False
                self.vtk_widget = None
                index = self.volume_engine.findData("vtk")
                if index >= 0:
                    self.volume_engine.removeItem(index)

        layout.addWidget(self.volume_stack, 1)
        self.volume_canvas.setFocusPolicy(Qt.WheelFocus)
        self.volume_canvas.mpl_connect("scroll_event", self._on_volume_scroll)

        hint = QLabel("Siçanla sürüşdürərək modeli fırlat. Yaxınlaşdırma "
                      "100 % = tam model çərçivədə. Kəsim həddi aşağı "
                      "dəyərləri gizlədir — daxili struktur görünür. "
                      "İşıq sürgüsü formanı, şəffaflıq isə daxili qatları "
                      "aydınlaşdırır.")
        hint.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        layout.addWidget(hint)

        # `addItem("VTK …")` yuxarıda combo-nu boşdan ilk elementə keçirir
        # və bunu edərkən `currentIndexChanged`-i DƏRHAL atəşləyir — bu isə
        # `connect()`-dən ƏVVƏL baş verir (VTK ilk element olduğu üçün).
        # Nəticədə ilk seçim heç vaxt `_on_engine_changed`-ə çatmır: yığın
        # matplotlib kətanında qalır və VTK interaktoru işə salınmır, halbuki
        # combo "VTK" göstərir. `Initialize()` isə yalnız pəncərə göründükdən
        # sonra təhlükəsizdir, ona görə sinxronizasiyanı bir dəfə, hadisə
        # dövrü pəncərəni göstərdikdən sonra həyata keçiririk.
        QTimer.singleShot(0, self._on_engine_changed)
        return page

    def _on_engine_changed(self):
        """Motor dəyişəndə uyğun widget-ə keç və yenidən çək."""
        use_vtk = (self.volume_engine.currentData() == "vtk"
                   and self.vtk_widget is not None)
        self.volume_stack.setCurrentWidget(
            self.vtk_widget if use_vtk else self.volume_canvas)
        if use_vtk and self.vtk_widget is not None:
            # interaktor yalnız GÖRÜNƏN olduqdan sonra işə salına bilər
            self.vtk_widget.Initialize()
            self.vtk_widget.Start()
        self.update_volume()

    def _update_volume_vtk(self, values, key, colormap, limits, volume_filter,
                           view):
        """VTK səhnəsini yenilə — həndəsə YALNIZ model dəyişəndə qurulur.

        Bu, matplotlib motorundan əsas fərqdir: orada hər yeniləmədə
        bütün üzlər yenidən hesablanırdı. Burada səhnə keşlənir, zaman
        slider-i sürüşdürəndə yalnız DƏYƏRLƏR yenilənir.
        """
        from ..rendering.vtk_volume import VtkReservoirScene, VtkViewSettings

        settings = VtkViewSettings(
            colormap=colormap or "viridis",
            value_limits=limits,
            k_range=volume_filter.k_range,
            value_min=volume_filter.value_min,
            # status filtri HƏR İKİ motorda EYNİ maskadan gəlir — VTK və
            # matplotlib görüntüsü bir-birindən FƏRQLƏNMƏMƏLİDİR.
            cell_mask=volume_filter.cell_mask,
            show_edges=self.volume_edges.isChecked(),
            opacity=self.volume_opacity.value() / 100.0,
            vertical_exaggeration=1.0,
            show_wells=self.volume_wells.isChecked(),
            show_faults=self.volume_faults.isChecked(),
            shading=self.volume_shading.value() / 100.0,
            zoom=self.volume_zoom.value() / 100.0)

        rebuild = (self.vtk_scene is None
                   or self.vtk_scene.model is not self.reservoir_model)
        if rebuild:
            self.vtk_scene = VtkReservoirScene(self.reservoir_model, settings)
            window = self.vtk_widget.GetRenderWindow()
            window.RemoveRenderer(window.GetRenderers().GetFirstRenderer()) \
                if window.GetRenderers().GetNumberOfItems() else None
            window.AddRenderer(self.vtk_scene.renderer)
        else:
            self.vtk_scene.settings = settings

        # istiqamət oxu (sağ aşağı künc) — interaktor tələb edir, ona
        # görə səhnə qurulanda yox, burada bağlanır
        self.vtk_scene.attach_orientation_marker(
            self.vtk_widget.GetRenderWindow().GetInteractor())

        self.vtk_scene.update_values(values, R.property_label(key))
        if rebuild or view is not None:
            # `reset_camera()` istifadəçinin fırlatdığı bucağı SIFIRLAYIR
            # — ona görə yalnız model dəyişəndə və ya hazır baxış
            # seçiləndə çağırılır (matplotlib motorunda eyni qayda).
            self.vtk_scene.reset_camera(view)
        self.vtk_widget.GetRenderWindow().Render()

    def update_volume(self):
        if not hasattr(self, "volume_ax") or self.reservoir_model is None:
            return
        model = self.reservoir_model
        grid = model.grid

        for widget in (self.volume_k_from, self.volume_k_to):
            widget.blockSignals(True)
            widget.setMaximum(grid.nz)
            widget.blockSignals(False)
        if self.volume_k_to.value() < self.volume_k_from.value():
            self.volume_k_to.setValue(grid.nz)

        snapshot = None
        if self.result and self.result.snapshots:
            snapshot = self.result.snapshots[
                min(self.volume_time.value(), len(self.result.snapshots) - 1)]

        key = self.volume_property.currentData()
        values, colormap, low, high = self.map_renderer._select_volume(
            model, key, snapshot)
        values = np.asarray(values, float).ravel()

        percentile = self.volume_threshold.value()
        if percentile > 0:
            cut = float(np.percentile(values, percentile))
            self.volume_threshold_label.setText(f"≥ {cut:.4g}")
        else:
            cut = None
            self.volume_threshold_label.setText("hamısı")

        volume_filter = VolumeFilter(
            value_min=cut,
            k_range=(self.volume_k_from.value() - 1,
                     self.volume_k_to.value() - 1),
            cell_mask=self._status_cell_mask())

        time = snapshot.time if snapshot else 0.0
        self.volume_time_label.setText(f"t = {time:8.0f} gün")

        # istifadəçinin fırlatdığı bucaq qorunur: yalnız hazır baxış
        # seçiləndə dəyişdirilir
        view = self.volume_view.currentData()

        if (self.volume_engine.currentData() == "vtk"
                and self.vtk_widget is not None):
            try:
                self._update_volume_vtk(
                    values, key, colormap,
                    (low, high) if low is not None and high is not None
                    else None,
                    volume_filter, view)
            except Exception:
                # VTK çəkilişi uğursuz oldusa (sürücü problemi və s.),
                # istifadəçini boş ekranla qoymuruq — matplotlib-ə
                # qayıdırıq və səbəbi loqa yazırıq.
                LOG.exception("VTK çəkilişi uğursuz — matplotlib-ə qayıdılır")
                index = self.volume_engine.findData("matplotlib")
                if index >= 0:
                    self.volume_engine.setCurrentIndex(index)
            return
        if view is None:
            elevation = self.volume_ax.elev
            azimuth = self.volume_ax.azim
        else:
            elevation = azimuth = None

        self.volume_colorbar = self.volume_renderer.draw(
            self.volume_ax, self.volume_fig, model, values,
            label=R.property_label(key), colormap=colormap,
            value_limits=(low, high) if low is not None and high is not None
            else None,
            volume_filter=volume_filter,
            show_wells=self.volume_wells.isChecked(),
            show_faults=self.volume_faults.isChecked(),
            edge_width=0.15 if self.volume_edges.isChecked() else 0.0,
            cax=self.volume_cax,
            shading=self.volume_shading.value() / 100.0,
            opacity=self.volume_opacity.value() / 100.0,
            zoom=self.volume_zoom.value() / 100.0,
            view=view)
        if view is None and elevation is not None:
            self.volume_ax.view_init(elev=elevation, azim=azimuth)
        self.volume_canvas.draw_idle()

    def _on_volume_scroll(self, event):
        """Siçan çarxı ilə yaxınlaşdırma.

        Səhnə yenidən çəkilmir — yalnız kamera nisbəti dəyişir, ona görə
        böyük modellərdə də hərəkət hamardır.
        """
        step = 1.15 if event.button == "up" else 1.0 / 1.15
        value = int(round(self.volume_zoom.value() * step))
        value = max(self.volume_zoom.minimum(),
                    min(self.volume_zoom.maximum(), value))
        if value == self.volume_zoom.value():
            return

        self.volume_zoom.blockSignals(True)
        self.volume_zoom.setValue(value)
        self.volume_zoom.blockSignals(False)

        if apply_zoom(self.volume_ax, value / 100.0):
            self.volume_canvas.draw_idle()
        else:
            self.update_volume()

    def save_volume_image(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "3D görüntünü saxla", "model_3d.png",
            "PNG (*.png);;PDF (*.pdf)")
        if not path:
            return
        try:
            self.volume_fig.savefig(path, dpi=200,
                                    facecolor=self.volume_fig.get_facecolor())
        except Exception as error:
            QMessageBox.critical(self, "Yazılmadı", str(error))
            return
        self.statusBar().showMessage(f"Saxlanıldı: {os.path.basename(path)}")

    def reset_volume_view(self):
        self.volume_view.setCurrentIndex(
            self.volume_view.findData("İzometrik"))
        self.volume_shading.setValue(45)
        self.volume_opacity.setValue(100)
        self.volume_zoom.setValue(100)
        self.update_volume()

    def _verify_build(self) -> None:
        """Faktiki tablar gözlənilənlə uyğun gəlirmi?

        Köhnə `main_window.py` başqa qovluqdan qarışanda proqram
        işləyir, lakin bəzi tablar yox olur və səbəbi görünmür.
        Bu yoxlama fərqi dərhal jurnala yazır.
        """
        actual = [self.tabs.tabText(index) for index in range(self.tabs.count())]
        missing = [name for name in EXPECTED_TABS if name not in actual]
        LOG.info("IMEX-2D v%s  ·  %d tab: %s", VERSION, len(actual),
                 ", ".join(actual))
        if missing:
            LOG.error("DİQQƏT: gözlənilən tablar yoxdur: %s. "
                      "Böyük ehtimalla köhnə fayl qarışıb — arxivi "
                      "təmiz qovluğa aç.", ", ".join(missing))
            self.statusBar().showMessage(
                f"Diqqət: {', '.join(missing)} tabı yoxdur — köhnə fayl "
                f"qarışmış ola bilər (bax: Jurnal).")

    def _build_history_tab(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)

        top = QHBoxLayout()
        load = QPushButton("Müşahidə CSV yüklə…")
        load.clicked.connect(self.load_observations)
        clear = QPushButton("Təmizlə")
        clear.clicked.connect(self.clear_observations)
        self.history_info = QLabel("Müşahidə məlumatı yüklənməyib.")
        self.history_info.setWordWrap(True)
        self.history_info.setStyleSheet(
            f"background:{PALETTE.panel_alt};border:1px solid {PALETTE.line};"
            f"border-radius:3px;padding:6px;font-family:monospace;"
            f"font-size:11px;color:{PALETTE.text}")
        top.addWidget(load)
        top.addWidget(clear)
        top.addWidget(self.history_info, 1)
        layout.addLayout(top)

        self.history_table = QTableWidget(0, 5)
        self.history_table.setHorizontalHeaderLabels(
            ["Sıra", "NRMSE", "RMSE", "Meyl", "Korrelyasiya"])
        self.history_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch)
        self.history_table.verticalHeader().setVisible(False)
        self.history_table.setMaximumHeight(170)
        layout.addWidget(self.history_table)

        self.history_fig, self.history_canvas, history_axes = _figure(2, 2)
        self.history_axes = list(np.ravel(history_axes))
        layout.addWidget(self.history_canvas, 2)

        self.crossplot_fig, self.crossplot_canvas, self.crossplot_ax = _figure()
        layout.addWidget(self.crossplot_canvas, 1)

        hint = QLabel("CSV sütunları: time, well, quantity, value. "
                      "Boş 'well' yataq səviyyəsi deməkdir. "
                      "Kəmiyyət adları: OIL_RATE / WOPR, WATER_CUT / WWCT, "
                      "FPR, CUM_OIL …")
        hint.setWordWrap(True)
        hint.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        layout.addWidget(hint)
        return page

    def _build_matching_tab(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)

        controls = QHBoxLayout()
        self.match_method = QComboBox()
        self.match_method.addItems(list(HistoryMatchingService.METHODS))
        self.match_budget = QSpinBox()
        self.match_budget.setRange(10, 2000)
        self.match_budget.setValue(60)
        self.match_budget.setPrefix("büdcə ")
        self.match_start = QPushButton("Uyğunlaşdırmanı başlat")
        self.match_start.setObjectName("run")
        self.match_start.clicked.connect(self.start_matching)
        self.match_stop = QPushButton("Dayandır")
        self.match_stop.setObjectName("stop")
        self.match_stop.setEnabled(False)
        self.match_stop.clicked.connect(self.stop_matching)
        self.match_apply = QPushButton("Nəticəni modelə tətbiq et")
        self.match_apply.setEnabled(False)
        self.match_apply.clicked.connect(self.apply_match_result)

        controls.addWidget(QLabel("Üsul:"))
        controls.addWidget(self.match_method)
        controls.addWidget(self.match_budget)
        controls.addWidget(self.match_start)
        controls.addWidget(self.match_stop)
        controls.addWidget(self.match_apply)
        controls.addStretch()
        layout.addLayout(controls)

        self.match_info = QLabel(
            "Müşahidə yükləyin (Tarixçə tabı) və modeli işə salın.")
        self.match_info.setWordWrap(True)
        self.match_info.setStyleSheet(
            f"background:{PALETTE.panel_alt};border:1px solid {PALETTE.line};"
            f"border-radius:3px;padding:6px;font-family:monospace;"
            f"font-size:11px;color:{PALETTE.text}")
        layout.addWidget(self.match_info)

        self.match_table = QTableWidget(0, 5)
        self.match_table.setHorizontalHeaderLabels(
            ["Parametr", "Başlanğıc", "Nəticə", "Hədlər", "İzah"])
        self.match_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch)
        self.match_table.verticalHeader().setVisible(False)
        self.match_table.setMaximumHeight(220)
        layout.addWidget(self.match_table)

        self.match_fig, self.match_canvas, match_axes = _figure(1, 2)
        self.match_axes = list(np.ravel(match_axes))
        layout.addWidget(self.match_canvas, 1)

        hint = QLabel("Hər qiymətləndirmə bir simulyasiyadır — büdcəni "
                      "model ölçüsünə görə seçin. Axtarış [0, 1] fəzasında "
                      "aparılır; uğursuz variantlar cərimə alır və axtarışı "
                      "dayandırmır.")
        hint.setWordWrap(True)
        hint.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        layout.addWidget(hint)
        return page

    def _build_sensitivity_tab(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)

        controls = QHBoxLayout()
        self.sensitivity_method = QComboBox()
        self.sensitivity_method.addItems(["Tornado", "Yerli elastiklik"])
        self.sensitivity_method.currentIndexChanged.connect(
            self._on_sensitivity_method_changed)

        self.sensitivity_metric = QComboBox()
        self.sensitivity_metric.addItems(list(OUTPUT_METRICS))

        self.sensitivity_step = QDoubleSpinBox()
        self.sensitivity_step.setRange(0.01, 0.5)
        self.sensitivity_step.setSingleStep(0.01)
        self.sensitivity_step.setValue(0.10)
        self.sensitivity_step.setPrefix("addım ")
        self.sensitivity_step.setToolTip(
            "Baza nöqtəsindən [0,1] fəzasında ± addım (yalnız yerli "
            "elastiklik üçün).")

        self.sensitivity_start = QPushButton("Təhlili başlat")
        self.sensitivity_start.setObjectName("run")
        self.sensitivity_start.clicked.connect(self.start_sensitivity)
        self.sensitivity_stop = QPushButton("Dayandır")
        self.sensitivity_stop.setObjectName("stop")
        self.sensitivity_stop.setEnabled(False)
        self.sensitivity_stop.clicked.connect(self.stop_sensitivity)

        controls.addWidget(QLabel("Üsul:"))
        controls.addWidget(self.sensitivity_method)
        controls.addWidget(QLabel("Çıxış:"))
        controls.addWidget(self.sensitivity_metric)
        controls.addWidget(self.sensitivity_step)
        controls.addWidget(self.sensitivity_start)
        controls.addWidget(self.sensitivity_stop)
        controls.addStretch()
        layout.addLayout(controls)

        self.sensitivity_info = QLabel(
            "Modeli işə salın, sonra təhlili başladın. Hər parametr üçün "
            "iki əlavə simulyasiya aparılır (aşağı və yuxarı hədd).")
        self.sensitivity_info.setWordWrap(True)
        self.sensitivity_info.setStyleSheet(
            f"background:{PALETTE.panel_alt};border:1px solid {PALETTE.line};"
            f"border-radius:3px;padding:6px;font-family:monospace;"
            f"font-size:11px;color:{PALETTE.text}")
        layout.addWidget(self.sensitivity_info)

        self.sensitivity_table = QTableWidget(0, 5)
        self.sensitivity_table.setHorizontalHeaderLabels(
            ["Parametr", "Yayılma", "Aşağı", "Yuxarı", "İstiqamət"])
        self.sensitivity_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch)
        self.sensitivity_table.verticalHeader().setVisible(False)
        self.sensitivity_table.setMaximumHeight(200)
        layout.addWidget(self.sensitivity_table)

        self.sensitivity_fig, self.sensitivity_canvas, self.sensitivity_ax = \
            _figure()
        layout.addWidget(self.sensitivity_canvas, 1)

        hint = QLabel("Tornado: hər parametr öz TAM hədləri arasında, "
                      "digərləri baza dəyərində — 'hansı parametr öz "
                      "diapazonunda ən çox təsir edir?'. Yerli elastiklik: "
                      "baza nöqtəsi ətrafında kiçik addım — 'hazırkı "
                      "modeldə kiçik dəyişiklik nəyə təsir edir?'. İkisi "
                      "fərqli sual cavablandırır və fərqli sıralama verə bilər.")
        hint.setWordWrap(True)
        hint.setStyleSheet(f"color:{PALETTE.text_dim};font-size:11px")
        layout.addWidget(hint)
        self._on_sensitivity_method_changed()
        return page

    # ═══════════════════════════════════════════════ həssaslıq analizi
    def _on_sensitivity_method_changed(self):
        self.sensitivity_step.setEnabled(
            self.sensitivity_method.currentText() == "Yerli elastiklik")

    def start_sensitivity(self):
        self.rebuild_model()
        if self.reservoir_model is None:
            return

        parameters = ParameterSet(standard_parameters(self.reservoir_model))
        config = self.numerical_panel.simulation_config()
        engine_factory = (FullyImplicitEngine
                          if self.numerical_panel.engine_choice() == "IMPLICIT"
                          else ImpesEngine)
        analyzer = SensitivityAnalyzer(
            self.reservoir_model, parameters,
            self.service.with_engine(engine_factory), config)

        self.sensitivity_start.setEnabled(False)
        self.sensitivity_stop.setEnabled(True)
        self.sensitivity_table.setRowCount(0)
        self.sensitivity_info.setText(
            f"{self.sensitivity_method.currentText()} · "
            f"{self.sensitivity_metric.currentText()} · "
            f"{len(parameters)} parametr\nTəhlil başladı…")

        self.sensitivity_worker = SensitivityWorker(
            analyzer, self.sensitivity_method.currentText(),
            self.sensitivity_metric.currentText(),
            self.sensitivity_step.value())
        self.sensitivity_worker.progress.connect(self._on_sensitivity_progress)
        self.sensitivity_worker.finished_ok.connect(self._on_sensitivity_finished)
        self.sensitivity_worker.failed.connect(self._on_sensitivity_failed)
        self.sensitivity_worker.start()

    def stop_sensitivity(self):
        if self.sensitivity_worker:
            self.sensitivity_worker.request_stop()
            self.statusBar().showMessage("Həssaslıq təhlili dayandırılır…")

    def _on_sensitivity_progress(self, done, total):
        self.statusBar().showMessage(
            f"Həssaslıq analizi: {done} / {total} qiymətləndirmə")

    def _on_sensitivity_failed(self, traceback_text):
        self.sensitivity_start.setEnabled(True)
        self.sensitivity_stop.setEnabled(False)
        LOG.error("Həssaslıq analizi xətası:\n%s", traceback_text)
        QMessageBox.critical(self, "Xəta", traceback_text.splitlines()[-1])

    def _on_sensitivity_finished(self, report):
        self.sensitivity_report = report
        self.sensitivity_start.setEnabled(True)
        self.sensitivity_stop.setEnabled(False)

        self.sensitivity_table.setRowCount(0)
        for item in report.sorted_by_swing():
            row = self.sensitivity_table.rowCount()
            self.sensitivity_table.insertRow(row)
            direction = "tərs" if item.direction_reversed else "düz"
            if item.failed_low or item.failed_high:
                direction += "  (uğursuz hədd)"
            for column, text in enumerate([
                    item.name, f"{item.swing:.4g}", f"{item.low_output:.4g}",
                    f"{item.high_output:.4g}", direction]):
                self.sensitivity_table.setItem(row, column,
                                               QTableWidgetItem(text))

        self.tornado_renderer.draw(self.sensitivity_ax, report)
        self.sensitivity_canvas.draw_idle()

        self.sensitivity_info.setText(
            f"{self.sensitivity_method.currentText()} · {report.metric_name}"
            f"   ·   baza = {report.baseline_output:.4g}"
            f"{'   ·   ' + str(report.failures) + ' uğursuz' if report.failures else ''}")
        LOG.info("%s", report.as_text())

    # ═══════════════════════════════════════ avtomatik uyğunlaşdırma
    def start_matching(self):
        if self.observations is None:
            QMessageBox.information(
                self, "Müşahidə yoxdur",
                "Əvvəlcə 'Tarixçə' tabında müşahidə faylı yükləyin.")
            return
        self.rebuild_model()
        if self.reservoir_model is None:
            return

        parameters = ParameterSet(standard_parameters(self.reservoir_model))
        config = self.numerical_panel.simulation_config()
        engine_factory = (FullyImplicitEngine
                          if self.numerical_panel.engine_choice() == "IMPLICIT"
                          else ImpesEngine)
        service = HistoryMatchingService(
            self.reservoir_model, parameters, self.observations,
            self.service.with_engine(engine_factory), config)

        self._fill_match_table(parameters)
        self.match_start.setEnabled(False)
        self.match_stop.setEnabled(True)
        self.match_apply.setEnabled(False)
        self.match_info.setText(
            f"{self.match_method.currentText()} · büdcə "
            f"{self.match_budget.value()} · {len(parameters)} parametr\n"
            f"Axtarış başladı…")

        self.match_worker = MatchingWorker(service,
                                           self.match_method.currentText(),
                                           self.match_budget.value())
        self.match_worker.progress.connect(self._on_match_progress)
        self.match_worker.finished_ok.connect(self._on_match_finished)
        self.match_worker.failed.connect(self._on_match_failed)
        self.match_worker.start()

    def stop_matching(self):
        if self.match_worker:
            self.match_worker.request_stop()
            self.statusBar().showMessage("Uyğunlaşdırma dayandırılır…")

    def _fill_match_table(self, parameters, best_values=None):
        self.match_table.setRowCount(0)
        for index, definition in enumerate(parameters.definitions):
            row = self.match_table.rowCount()
            self.match_table.insertRow(row)
            best = ("—" if best_values is None
                    else f"{best_values[index]:.4f}")
            cells = [definition.name, f"{definition.initial:.4f}", best,
                     f"{definition.minimum:g} … {definition.maximum:g}"
                     + (" (log)" if definition.log_scale else ""),
                     definition.description]
            for column, text in enumerate(cells):
                self.match_table.setItem(row, column, QTableWidgetItem(text))

    def _on_match_progress(self, evaluation):
        status = "uğursuz" if not evaluation.succeeded else \
            f"{evaluation.mismatch:.5f}"
        self.statusBar().showMessage(
            f"Uyğunlaşdırma: {evaluation.iteration} qiymətləndirmə  ·  "
            f"cari {status}  ·  {evaluation.seconds:.1f} san")

    def _on_match_failed(self, traceback_text):
        self.match_start.setEnabled(True)
        self.match_stop.setEnabled(False)
        LOG.error("Uyğunlaşdırma xətası:\n%s", traceback_text)
        QMessageBox.critical(self, "Uyğunlaşdırma xətası",
                             traceback_text.splitlines()[-1])

    def _on_match_finished(self, result):
        self.match_result = result
        self.match_start.setEnabled(True)
        self.match_stop.setEnabled(False)
        self.match_apply.setEnabled(True)

        self._fill_match_table(result.parameters, result.best_values)
        self.optimisation_renderer.draw(self.match_axes, result)
        self.match_canvas.draw_idle()

        text = (f"{result.method} · {result.evaluations} qiymətləndirmə "
                f"({result.failures} uğursuz)"
                f"{' · dayandırıldı' if result.stopped_early else ''}\n"
                f"Uyğunsuzluq: {result.initial_mismatch:.5f} → "
                f"{result.best_mismatch:.5f}   "
                f"({result.improvement:+.1f} % yaxşılaşma)")
        self.match_info.setText(text)
        LOG.info("%s", result.summary())

        if result.best_report is not None:
            self.mismatch_report = result.best_report
            self.history_renderer.draw(self.history_fig, self.history_axes,
                                       result.best_report)
            self.crossplot_renderer.draw(self.crossplot_ax,
                                         result.best_report)
            self.history_canvas.draw_idle()
            self.crossplot_canvas.draw_idle()

    def apply_match_result(self):
        """Tapılan parametrləri interfeys sahələrinə yazır."""
        if self.match_result is None:
            return
        values = self.match_result.as_dict()
        self._ready = False
        try:
            scal = self.scal_panel
            for name, widget in (("SOR", scal.sor), ("SWC", scal.swc),
                                 ("KRW_END", scal.krw_end),
                                 ("COREY_NW", scal.nw),
                                 ("COREY_NO", scal.no)):
                if name in values:
                    widget.setValue(values[name])
            if "MU_OIL" in values:
                self.rock_panel.mu_o.setValue(values["MU_OIL"])
            if "PERM_MULT" in values:
                self.rock_panel.permx.setValue(
                    self.rock_panel.permx.value() * values["PERM_MULT"])
            if "PORO_MULT" in values:
                self.rock_panel.porosity.setValue(
                    self.rock_panel.porosity.value() * values["PORO_MULT"])
            if "KV_KH" in values:
                self.rock_panel.kv_over_kh.setValue(values["KV_KH"])
            if "OWC" in values:
                self.numerical_panel.owc.setValue(values["OWC"])
        finally:
            self._ready = True
        self.rebuild_model()
        self.statusBar().showMessage(
            "Uyğunlaşdırma nəticəsi panellərə tətbiq edildi.")

    # ═══════════════════════════════════════════ tarixçə uyğunluğu
    def load_observations(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Müşahidə məlumatı", "", "CSV (*.csv *.txt);;Hamısı (*)")
        if not path:
            return
        try:
            self.observations = read_observations_csv(path)
        except ObservationFormatError as error:
            QMessageBox.warning(self, "Fayl oxunmadı", str(error))
            return
        except Exception as error:
            QMessageBox.critical(self, "Fayl oxunmadı",
                                 f"Gözlənilməz xəta: {error}")
            LOG.exception("Müşahidə faylı oxunmadı")
            return
        LOG.info("Müşahidə yükləndi: %s — %s", os.path.basename(path),
                 self.observations.summary())
        self.update_history_match()

    def clear_observations(self):
        self.observations = None
        self.mismatch_report = None
        self.update_history_match()

    def update_history_match(self):
        if not hasattr(self, "history_axes"):
            return
        if self.observations is None:
            self.history_info.setText("Müşahidə məlumatı yüklənməyib.")
            self.history_table.setRowCount(0)
            self.history_renderer.draw(self.history_fig, self.history_axes,
                                       R.MismatchReportPlaceholder())
            self.crossplot_renderer.draw(self.crossplot_ax,
                                         R.MismatchReportPlaceholder())
            self.history_canvas.draw_idle()
            self.crossplot_canvas.draw_idle()
            return

        summary = self.observations.summary()
        if self.result is None:
            self.history_info.setText(
                f"{summary['sıra']} sıra · {summary['nöqtə']} nöqtə · "
                f"{summary['müddət']}\nModeli işə salın — uyğunluq sonra "
                f"hesablanacaq.")
            return

        report = MismatchCalculator().evaluate(self.result, self.observations)
        self.mismatch_report = report

        self.history_table.setRowCount(0)
        for label, rmse, nrmse, bias, correlation, _ in report.as_rows():
            row = self.history_table.rowCount()
            self.history_table.insertRow(row)
            for column, text in enumerate([label, f"{nrmse:.4f}",
                                           f"{rmse:.3f}", f"{bias:+.3f}",
                                           f"{correlation:.3f}"]):
                self.history_table.setItem(row, column,
                                           QTableWidgetItem(text))

        worst = report.worst
        text = (f"Yekun uyğunsuzluq: {report.total:.4f}   "
                f"({len(report.series)} sıra · {summary['nöqtə']} nöqtə)")
        if worst is not None:
            text += f"\nƏn pis uyğunluq: {worst.label} (NRMSE {worst.nrmse:.3f})"
        if report.skipped:
            text += ("\nAtlanan: " + ", ".join(report.skipped)
                     + "  — modeldə qarşılığı yoxdur")
        self.history_info.setText(text)

        self.history_renderer.draw(self.history_fig, self.history_axes, report)
        self.crossplot_renderer.draw(self.crossplot_ax, report)
        self.history_canvas.draw_idle()
        self.crossplot_canvas.draw_idle()
        LOG.info("Tarixçə uyğunluğu: %.4f", report.total)

    def show_tab(self, title: str) -> None:
        """Tabı ADLA seçir.

        Sərt indekslər kövrəkdir: yeni tab əlavə olunanda hamısı sürüşür
        və proqram səhv səhifəyə keçir (3D tab əlavə olunanda məhz bu
        baş verdi). Ad dəyişməzdir.
        """
        for index in range(self.tabs.count()):
            if self.tabs.tabText(index) == title:
                self.tabs.setCurrentIndex(index)
                return
        LOG.debug("Tab tapılmadı: %s", title)

    def _append_log(self, text: str):
        if self.log.document().blockCount() > 5000:
            self.log.clear()
        self.log.append(text)

    # ══════════════════════════════════ iş axını: geoloji → rezervuar
    def rebuild_model(self):
        if not getattr(self, "_ready", False):
            return
        if self.imported_geology is not None:
            imported = self.imported_geology.grid
            if (self.grid_panel.nx.value(), self.grid_panel.ny.value(),
                    self.grid_panel.nz.value()) != (imported.nx, imported.ny,
                                                    imported.nz):
                LOG.info("Grid ölçüsü dəyişdirildi — GRDECL modeli ləğv olundu.")
                self.imported_geology = None
        try:
            geology = self._build_geological_model()
            self.project.add_geological_model(geology)
            model = self.model_builder.build(
                geological_model=geology,
                wells=self.well_panel.values(),
                fluids=self.rock_panel.fluids(),
                scal=self.scal_panel.values(),
                capillary=self.scal_panel.capillary_values(),
                initial=self.numerical_panel.initial_conditions(),
                pvt_table=self.pvt_panel.values(),
                scal_tables=self.scal_tables(),
                fault_references=self.fault_panel.values(),
                rock_compressibility=self.rock_panel.rock_compressibility_value(),
                name="Aktiv rezervuar modeli")
            self.project.add_reservoir_model(model)
            self.reservoir_model = model
        except Exception as exc:
            self.statusBar().showMessage(f"Model qurula bilmədi: {exc}")
            return
        self.refresh_tree()
        self._refresh_provenance_choices()
        self.update_scal_plot()
        self.update_pvt_plot()
        self.update_map()
        self.update_volume()

    def _status_cell_mask(self):
        """Status filtrinin bool maskası (`None` — filtr yoxdur).

        MASKA MODELİN `provenance`-INDAN OXUNUR — UI heç bir geologiya
        hesablaması aparmır (§13), yalnız hazır statusları müqayisə edir.
        """
        if self.reservoir_model is None:
            return None
        name = self.volume_status_property.currentData()
        wanted = self.volume_status.currentData()
        if not name or not wanted:
            return None
        entry = (getattr(self.reservoir_model, "provenance", {}) or {}).get(name)
        if entry is None:
            return None
        return entry.mask(wanted)

    def _refresh_provenance_choices(self):
        """Modeldə provenance varsa xəritə/3D combo-larına "Mənşə/status"
        və "Etibarlılıq balı" seçimlərini əlavə edir (§13).

        Siyahı MODELDƏN qurulur — heç bir xassə adı burada SABİT
        KODLANMIR. Provenance yoxdursa (köhnə rejim) heç nə əlavə
        olunmur, combo ƏVVƏLKİ kimi qalır.
        """
        provenance = getattr(self.reservoir_model, "provenance", {}) or {}
        wanted = []
        for name in sorted(provenance):
            keys_for_name = [R.provenance_key(name), R.confidence_key(name)]
            # ORİJİNAL və TƏSİR yalnız orijinal sahə MÖVCUD olanda mənalıdır
            # — əks halda boş (NaN) görüntü təklif etmirik.
            if provenance[name].original is not None:
                keys_for_name += [R.original_key(name), R.impact_key(name)]
            for key in keys_for_name:
                wanted.append((key, R.property_label(key)))
        keys = {key for key, _label in wanted}

        # status filtri: hansı xassənin mənşəyinə görə süzüləcəyi
        selected = self.volume_status_property.currentData()
        self.volume_status_property.blockSignals(True)
        self.volume_status_property.clear()
        self.volume_status_property.addItem("— status filtri yoxdur —", "")
        for name in sorted(provenance):
            self.volume_status_property.addItem(name, name)
        index = self.volume_status_property.findData(selected or "")
        self.volume_status_property.setCurrentIndex(max(index, 0))
        self.volume_status_property.blockSignals(False)
        for combo in (self.property_box, self.volume_property):
            for index in range(combo.count() - 1, -1, -1):
                data = str(combo.itemData(index) or "")
                if data.startswith(R.PROVENANCE_PREFIXES) and data not in keys:
                    combo.removeItem(index)
            existing = {combo.itemData(i) for i in range(combo.count())}
            for key, label in wanted:
                if key not in existing:
                    combo.addItem(label, key)

    def scal_tables(self):
        """Yüklənmiş SCAL cədvəlləri (yoxdursa None)."""
        return (self.scal_source_panel.tables
                if self.scal_source_panel.is_enabled() else None)

    def _build_geological_model(self):
        """Geoloji modelin mənbəyi: quyu cədvəli və ya sintetik generator.

        `İnterpolyasiya et` düyməsi ilə hesablanan model `_geology_model_from_wells`-
        də keşlənir — cədvəl dəyişəndə (`geology_panel.changed`) bura TOXUNULMUR,
        yalnız düymə basılanda (`_interpolate_geology`) yenilənir. Beləliklə
        cədvəldəki hər redaktə böyük gridi yenidən interpolyasiya etmir.
        """
        grid_values = self.grid_panel.values()
        rock_values = self.rock_panel.geology_values()

        if self.imported_geology is not None:
            self.geology_panel.set_report(
                f"GRDECL faylından oxunmuş model işlədilir: "
                f"{self.imported_geology.name}\n"
                f"Xassələr: {', '.join(sorted(self.imported_geology.property_maps))}\n"
                f"(Panel parametrləri bu modelə təsir etmir. Sintetik modelə "
                f"qayıtmaq üçün grid ölçüsünü dəyişin.)")
            return self.imported_geology

        if not self.geology_panel.wells():
            self.geology_panel.set_report(
                "Sintetik model işlədilir — geologiya cədvəli boşdur.")
            return self.geology_builder.build(**grid_values, **rock_values)

        if self._geology_model_from_wells is None:
            self.geology_panel.set_report(
                "Quyu cədvəli dolduruldu, amma hələ interpolyasiya edilməyib.\n"
                "'İnterpolyasiya et' düyməsini basın. Hazırda sintetik model işlədilir.")
            return self.geology_builder.build(**grid_values, **rock_values)

        return self._geology_model_from_wells

    def _wells_to_geology_rows(self, wells, geometry: CellGeometry) -> list:
        """İndeks-əsaslı quyuları (ssenari generatoru) metrə çevirib
        geologiya sətirlərinə çevirir — koordinat/rejim birlikdə.

        Yalnız MÖVQE köçürülür (petrofizika YOX) — ssenari generatoru
        heç vaxt φ/k/Sw dəyəri verməyib, bunu istifadəçi əl ilə doldurur.
        """
        rows = []
        for well in wells:
            if not well.perforations:
                continue
            p = well.perforations[0]
            rows.append(GeologicalWell(
                name=well.name, in_model=True,
                x=(p.i + 0.5) * geometry.dx, y=(p.j + 0.5) * geometry.dy))
        return rows

    def _current_geometry(self) -> CellGeometry:
        grid_values = self.grid_panel.values()
        grid = CartesianGrid(grid_values["nx"], grid_values["ny"], grid_values["nz"])
        return CellGeometry(grid, grid_values["dx"], grid_values["dy"],
                            grid_values["dz"], top_depth=grid_values["top_depth"])

    def _sync_geology_geometry(self):
        """Grid dəyişəndə (2·) və (7·) bölmələrinin (i,j,k) sütunları yenilənir."""
        geometry = self._current_geometry()
        self.geology_panel.set_geometry(geometry)
        self.well_panel.set_geology_context(self.geology_panel.wells(), geometry)

    def _on_geology_table_changed(self):
        if not getattr(self, "_ready", False):
            return
        self._mark_dirty()
        self.well_panel.set_geology_context(self.geology_panel.wells(), self._current_geometry())

    def _interpolate_geology(self):
        """'İnterpolyasiya et' düyməsi — yeganə yer ki, geologiya cədvəli
        grid xassələrinə çevrilir. `geology_service.py` TOXUNULMUR."""
        grid_values = self.grid_panel.values()
        rock_values = self.rock_panel.geology_values()
        geometry = self._current_geometry()
        wells = self.geology_panel.wells()
        method = self.geology_panel.method_text()
        issues = validate_wells(wells, geometry, method,
                                reservoir_well_names=[w.name for w in self.well_panel.values()])
        self.geology_panel.set_validation(issues)
        if any(issue.level == "error" for issue in issues):
            QMessageBox.warning(
                self, "İnterpolyasiya edilmədi",
                "Geologiya cədvəlində xəta var — aşağıdakı yoxlama panelinə baxın.")
            return

        spec = GeologicalGridSpec(
            nx=grid_values["nx"], ny=grid_values["ny"], nz=grid_values["nz"],
            dx=grid_values["dx"], dy=grid_values["dy"], dz=grid_values["dz"],
            top_depth=grid_values["top_depth"],
            dip_x=grid_values["dip_x"], dip_y=grid_values["dip_y"])
        # LAY-MƏLUMATLI rejim: `geometry`/`policy` verilir ki, quyunun
        # "Data layları" bəyanı HƏR LAY ÜÇÜN AYRICA nümunəyə çevrilsin.
        # Rejim söndürülübsə hər ikisi defolt qalır → ƏVVƏLKİ davranış.
        policy = self.geology_panel.layer_data_policy()
        try:
            layer_config = self.geology_panel.layer_config(grid_values["nz"])
            dataset, skipped = wells_to_dataset(
                wells, method,
                geometry=geometry if layer_config is not None else None,
                policy=policy)
        except ValueError as exc:
            QMessageBox.critical(self, "İnterpolyasiya edilmədi",
                                 f"Lay seçimi oxunmadı:\n{exc}")
            return
        builder = WellBasedGeologicalModelBuilder(self.geology_panel.interpolator())
        try:
            geology, report = builder.build(
                dataset, spec, ky_over_kx=rock_values["ky_over_kx"],
                kv_over_kh=rock_values["kv_over_kh"],
                allow_cross_layer_fallback=self.geology_panel.cross_layer_fallback_allowed(),
                name="Quyu cədvəlindən geoloji model",
                layer_config=layer_config)
        except ValueError as exc:
            QMessageBox.critical(self, "İnterpolyasiya edilmədi", str(exc))
            return

        text = report.as_text()
        if skipped:
            text += "\n\nBuraxılan xassələr:\n" + "\n".join(
                f"  {message}" for message in skipped.values())
        self.geology_panel.set_report(text)
        if report.has_blocking:
            # Model QAYTARILIR (3D-də MISSING laylar görünsün deyə), amma
            # rezervuar modeli qurulmur — istifadəçi bunu AÇIQ görməlidir.
            QMessageBox.warning(
                self, "Model natamamdır",
                "Aşağıdakı laylar üçün məlumat yoxdur və tamamlama üsulu "
                "seçilməyib — simulyasiya bu modellə işə düşməyəcək:\n\n"
                + "\n\n".join(report.blocking))
        self._geology_model_from_wells = geology
        self.geology_panel.mark_fresh()
        self.rebuild_model()

    def _cross_validate_geology(self):
        """'Cross-validation et' düyməsi — M4. Nəticəni QURULAN modelə
        DEYİL, birbaşa quyu cədvəlinə tətbiq edir (interpolyasiya
        edilməmiş olsa belə işləyir) — real dəqiqliyi göstərir, "100%
        dəqiq" vəd etmir."""
        wells = self.geology_panel.wells()
        method = self.geology_panel.method_text()
        geometry = self._current_geometry()
        layer_aware = self.geology_panel.layer_aware_enabled()
        try:
            dataset, _skipped = wells_to_dataset(
                wells, method, geometry=geometry if layer_aware else None,
                policy=self.geology_panel.layer_data_policy())
        except ValueError as exc:
            QMessageBox.critical(self, "Cross-validation", f"Lay bəyanı oxunmadı:\n{exc}")
            return
        if len(dataset) < 3:
            QMessageBox.information(
                self, "Cross-validation",
                "Cross-validation üçün ən azı 3 quyu nöqtəsi lazımdır.")
            return
        builder = WellBasedGeologicalModelBuilder(self.geology_panel.interpolator())
        # `nz` verilir ki, DOĞRULAMA MƏLUMATI OLMAYAN laylar da hesabatda
        # AÇIQ görünsün — "RMSE = 0" kimi saxta uğur yaranmasın (§19).
        all_results = builder.cross_validate_all(
            dataset, nz=geometry.grid.nz if layer_aware else None)
        text = format_cross_validation_report(all_results)
        QMessageBox.information(self, "Cross-validation nəticəsi", text)

    def _apply_pattern(self):
        """Ssenari daxildə İNDEKSLƏ işləyir (`five_spot` və s. dəyişmir) —
        yalnız tətbiq ediləndə metrə çevrilib geologiya cədvəlini doldurur."""
        existing = self.geology_panel.wells()
        if existing:
            reply = QMessageBox.question(
                self, "Ssenari tətbiq edilsin?",
                f"Mövcud {len(existing)} quyu əvəz olunacaq. Davam edilsin?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply != QMessageBox.Yes:
                return
        pattern = WELL_PATTERNS[self.well_panel.pattern.currentText()]
        wells = pattern(self.grid_panel.grid())
        geometry = self._current_geometry()
        self.well_panel.load(wells)
        self.geology_panel.load(self._wells_to_geology_rows(wells, geometry))
        self._sync_geology_geometry()

    def refresh_tree(self):
        self.tree.clear()
        root = QTreeWidgetItem([self.project.name])
        self.tree.addTopLevelItem(root)
        for geology in self.project.geological_models.values():
            node = QTreeWidgetItem([f"Geoloji model: {geology.name}"])
            node.addChild(QTreeWidgetItem(
                [f"Grid {geology.grid.nx}×{geology.grid.ny}×{geology.grid.nz}"
                 f" = {geology.grid.ncell} hüceyrə"]))
            node.addChild(QTreeWidgetItem(
                [f"Xassə xəritələri: {', '.join(geology.property_maps)}"]))
            node.addChild(QTreeWidgetItem([f"Horizontlar: {len(geology.horizons)}"]))
            node.addChild(QTreeWidgetItem([f"Faylar: {len(geology.faults)}"]))
            depths = geology.geometry.cell_depths()
            node.addChild(QTreeWidgetItem(
                [f"Dərinlik: {depths.min():.0f} – {depths.max():.0f} m"]))
            root.addChild(node)
        for model in self.project.reservoir_models.values():
            node = QTreeWidgetItem([f"Rezervuar modeli: {model.name}"])
            for key, value in model.summary().items():
                node.addChild(QTreeWidgetItem([f"{key}: {value}"]))
            for well in model.active_wells():
                layers = sorted({p.k + 1 for p in well.open_perforations()})
                node.addChild(QTreeWidgetItem(
                    [f"{well.name}: K {layers[0]}–{layers[-1]} "
                     f"({len(layers)} təbəqə)"]))
            ic = model.initial_conditions
            node.addChild(QTreeWidgetItem(
                [f"initialization: equilibration, OWC = {ic.oil_water_contact:.0f} m"
                 if ic.use_equilibration and ic.oil_water_contact is not None
                 else "initialization: bərabər paylanma"]))
            root.addChild(node)
        for run in self.project.runs.values():
            root.addChild(QTreeWidgetItem(
                [f"{run.run_id} — {run.status} ({run.config.end_time:.0f} gün)"]))
        self.tree.expandAll()

    # ══════════════════════════════════════════════════════ simulyasiya
    def run_simulation(self):
        self.rebuild_model()
        if self.reservoir_model is None:
            return
        config = self.numerical_panel.simulation_config()
        engine_factory = (FullyImplicitEngine
                          if self.numerical_panel.engine_choice() == "IMPLICIT"
                          else ImpesEngine)
        service = self.service.with_engine(engine_factory)

        report = self.reservoir_model.diagnose()
        for item in report.items:
            LOG.warning("%s", item) if item.severity is not Severity.INFO \
                else LOG.info("%s", item)
        if report.has_errors:
            QMessageBox.warning(
                self, "Model yoxlaması",
                "Model işə salına bilməz:\n\n"
                + "\n".join(f"• {d.message}" for d in report.errors))
            return
        if report.warnings and not self._confirm_warnings(report):
            return

        try:
            service.create_engine(self.reservoir_model, config)
        except ModelValidationError as exc:
            QMessageBox.warning(self, "Model yoxlaması", "\n".join(exc.issues))
            return
        except NotImplementedError as exc:
            QMessageBox.information(self, "Hələ hazır deyil", str(exc))
            return

        self.run = self.project.new_run(self.reservoir_model.name, config)
        self.run.status = "RUNNING"
        self.log.clear()
        self._log(f"{self.run.run_id}: {self.reservoir_model.summary()}")

        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.slider.setEnabled(False)
        self.worker = SimulationWorker(service, self.reservoir_model, config)
        self.worker.progress.connect(self._on_progress)
        self.worker.finished_ok.connect(self._on_finished)
        self.worker.failed.connect(self._on_failed)
        self.worker.start()

    def _confirm_warnings(self, report) -> bool:
        """Xəbərdarlıqlar bloklamır, amma istifadəçi onları görməlidir."""
        lines = []
        for diagnostic in report.warnings:
            lines.append(f"• {diagnostic.message}")
            if diagnostic.hint:
                lines.append(f"    {diagnostic.hint}")
        answer = QMessageBox.warning(
            self, "Diqqət",
            "Model işə salına bilər, amma aşağıdakılar çox güman səhvdir:\n\n"
            + "\n".join(lines) + "\n\nDavam edilsin?",
            QMessageBox.Yes | QMessageBox.Cancel, QMessageBox.Cancel)
        return answer == QMessageBox.Yes

    def stop_simulation(self):
        if self.worker:
            self.worker.request_stop()
            self.statusBar().showMessage("Dayandırılır…")

    def _on_progress(self, fraction, message):
        self.progress.setValue(int(fraction))
        self.statusBar().showMessage(message)
        LOG.debug("%s", message)

    def _on_failed(self, traceback_text):
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.run.status = "FAILED"
        LOG.error("Hesablama xətası:\n%s", traceback_text)
        self.show_tab("Jurnal")
        self.refresh_tree()
        QMessageBox.critical(self, "Hesablama xətası", traceback_text.splitlines()[-1])

    def _on_finished(self, result):
        self.result = result
        self.run.result = result
        self.run.status = "FINISHED" if result.converged else "FAILED"
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress.setValue(100)
        self._log(result.message)
        self.slider.setEnabled(True)
        self.slider.setRange(0, len(result.snapshots) - 1)
        self.slider.setValue(len(result.snapshots) - 1)
        self.statusBar().showMessage(result.message)
        self.volume_time.setEnabled(True)
        self.volume_time.setRange(0, len(result.snapshots) - 1)
        self.volume_time.setValue(len(result.snapshots) - 1)
        self.update_results()
        self.update_map()
        self.update_volume()
        self.refresh_tree()
        self.update_comparison()
        self.update_history_match()
        self.show_tab("Nəticələr")

    def _log(self, text):
        LOG.info("%s", text)

    # ══════════════════════════════════════════════════════════ görüntü
    def update_map(self):
        if not getattr(self, "_ready", False) or self.reservoir_model is None:
            return
        snapshot = None
        if self.result and self.result.snapshots:
            snapshot = self.result.snapshots[self.slider.value()]

        grid = self.reservoir_model.grid
        mode = self.view_mode.currentData()
        self.layer_spin.setEnabled(grid.nz > 1 or mode != "AREAL")
        if mode == "AREAL":
            self.layer_spin.setPrefix("K = ")
            self.layer_spin.setMaximum(max(grid.nz, 1))
            self.colorbar = self.map_renderer.draw(
                self.map_ax, self.map_fig, self.reservoir_model,
                self.property_box.currentData(), snapshot,
                layer=self.layer_spin.value() - 1, cax=self.map_cax)
        else:
            axis = "J" if mode == "SECTION_J" else "I"
            self.layer_spin.setPrefix(f"{axis} = ")
            self.layer_spin.setMaximum(grid.ny if axis == "J" else grid.nx)
            self.colorbar = self.section_renderer.draw(
                self.map_ax, self.map_fig, self.reservoir_model,
                self.property_box.currentData(), snapshot,
                axis=axis, index=self.layer_spin.value() - 1, cax=self.map_cax)
        self.time_label.setText(f"t = {(snapshot.time if snapshot else 0.0):8.0f} gün")
        self.map_canvas.draw_idle()

    def update_results(self):
        if not self.result:
            return
        self.curve_renderer.draw(self.result_axes, self.result)
        self.result_canvas.draw_idle()
        s = self.result.series
        bt = self.result.breakthrough_time
        self.kpi.setText(
            f"OOIP {self.result.ooip / 1e3:10.1f} min m³      "
            f"Kum. neft {s.cumulative_oil[-1] / 1e3:8.1f} min m³      "
            f"RF {s.recovery_factor[-1]:6.2f} %      "
            f"Water cut {s.water_cut[-1]:6.1f} %      "
            f"Orta P {s.average_pressure[-1]:6.1f} bar      "
            f"Su gəlişi {('%.0f gün' % bt) if bt else 'baş verməyib'}")

    def update_scal_plot(self):
        if self.reservoir_model is None:
            return
        fluids = self.reservoir_model.fluids
        capillary = None
        if self.reservoir_model.capillary_parameters.enabled:
            from ..simulation.capillary import BrooksCoreyCapillaryProvider
            capillary = BrooksCoreyCapillaryProvider(
                self.reservoir_model.capillary_parameters,
                self.reservoir_model.scal_parameters)
        self.scal_renderer.draw(self.scal_axes, self.reservoir_model.scal_parameters,
                                fluids.water_viscosity, fluids.oil_viscosity,
                                capillary)
        self.scal_canvas.draw_idle()

    def update_pvt_plot(self):
        if not hasattr(self, "pvt_axes"):
            return
        table = self.reservoir_model.pvt_table if self.reservoir_model else None
        if table is None:
            for row in self.pvt_axes:
                for ax in row:
                    ax.clear()
            self.pvt_info.setText("PVT söndürülüb — sabit flüid xassələri işlədilir "
                                  "(nəticələr köhnə versiya ilə eyni olacaq).")
        else:
            self.pvt_renderer.draw(self.pvt_axes, table)
            self.pvt_info.setText(
                f"Aktiv: {table.source}  ·  {table.size} nöqtə  ·  "
                f"Pb = {table.bubble_point:.0f} bar  ·  "
                f"Bo = {table.oil_fvf.min():.3f}–{table.oil_fvf.max():.3f}  ·  "
                f"μo = {table.oil_viscosity.min():.2f}–{table.oil_viscosity.max():.2f} cP")
        self.pvt_canvas.draw_idle()

    def toggle_play(self):
        if self._player.isActive():
            self._player.stop()
            self.play_button.setText("▶  Oynat")
        elif self.result:
            self._player.start(140)
            self.play_button.setText("❚❚  Dayan")

    def _next_frame(self):
        if not self.result:
            self._player.stop()
            return
        value = self.slider.value() + 1
        self.slider.setValue(0 if value > self.slider.maximum() else value)

    # ═════════════════════════════════════════════════ B-L validasiyası
    def run_validation(self):
        if self.reservoir_model is None:
            return
        base = self.reservoir_model
        nx, dx, dy, dz = 120, 8.0, 100.0, 10.0
        rate, end_time = 60.0, 250.0
        porosity = float(base.rock.porosity.values.mean())
        permeability = float(base.rock.permx.values.mean())

        geology = self.geology_builder.build(nx=nx, ny=1, dx=dx, dy=dy, dz=dz,
                                             porosity=porosity,
                                             permx_base=permeability,
                                             name="B-L validasiya modeli")
        wells = [
            Well("INJ", WellType.INJECTOR, WellControl(ControlMode.RATE, rate),
                 [Perforation(0, 0, 0)]),
            Well("PROD", WellType.PRODUCER, WellControl(ControlMode.BHP, 200.0),
                 [Perforation(nx - 1, 0, 0)]),
        ]
        scal = base.scal_parameters
        from ..domain.initial import InitialConditions
        model = self.model_builder.build(
            geology, wells, fluids=base.fluids, scal=scal,
            initial=InitialConditions(datum_pressure=200.0,
                                      water_saturation=scal.swc),
            name="B-L validasiya modeli")
        config = SimulationConfig(
            end_time=end_time,
            time_stepping=TimeSteppingConfig(max_dt=2.0, cfl_factor=0.4),
            output=OutputConfig(snapshot_count=2))
        result = self.service.run(model, config)

        analytical = buckley_leverett(scal, base.fluids.water_viscosity,
                                      base.fluids.oil_viscosity, porosity,
                                      rate, dy * dz, end_time)
        x_cells = (np.arange(nx) + 0.5) * dx
        sw_numeric = result.snapshots[-1].water_saturation.ravel()
        self.validation_renderer.draw(self.validation_ax, analytical, x_cells,
                                      sw_numeric, end_time, nx)
        self.validation_canvas.draw_idle()

        mask = sw_numeric > scal.swc + 0.01
        x_numeric = x_cells[mask][-1] if mask.any() else 0.0
        error = abs(x_numeric - analytical.front_position) / max(analytical.front_position, 1e-9) * 100
        self.validation_label.setText(
            f"Analitik front {analytical.front_position:.1f} m  ·  "
            f"ədədi front {x_numeric:.1f} m  ·  fərq {error:.1f} %  ·  "
            f"shock Sw = {analytical.shock_saturation:.3f}")
        self.show_tab("Validasiya (B-L)")

    # ═══════════════════════════════════════════════════════ müqayisə
    def finished_runs(self):
        return [run for run in self.project.runs.values()
                if run.result is not None and run.result.series.time]

    def update_comparison(self):
        if not hasattr(self, "comparison_axes"):
            return
        runs = self.finished_runs()
        labels = [(f"{run.run_id} · {run.config.end_time:.0f} gün", run.result)
                  for run in runs]
        self.comparison_renderer.draw(self.comparison_axes, labels)
        self.comparison_canvas.draw_idle()

        self.comparison_table.setRowCount(0)
        for run in runs:
            result = run.result
            breakthrough = result.breakthrough_time
            row = self.comparison_table.rowCount()
            self.comparison_table.insertRow(row)
            values = [run.run_id, run.reservoir_model_name,
                      f"{run.config.end_time:.0f}",
                      f"{result.final_recovery_factor:.2f}",
                      f"{breakthrough:.0f}" if breakthrough else "—",
                      f"{result.series.cumulative_oil[-1] / 1e3:.1f}",
                      f"{result.series.water_cut[-1]:.1f}"]
            for column, text in enumerate(values):
                self.comparison_table.setItem(row, column, QTableWidgetItem(text))

        if len(runs) < 2:
            self.comparison_info.setText(
                "Ən azı iki dəfə modeli işə salın — sonra fərq burada görünəcək.")
        else:
            best = max(runs, key=lambda r: r.result.final_recovery_factor)
            worst = min(runs, key=lambda r: r.result.final_recovery_factor)
            delta = (best.result.final_recovery_factor
                     - worst.result.final_recovery_factor)
            self.comparison_info.setText(
                f"{len(runs)} işə salınma  ·  ən yüksək RF: {best.run_id} "
                f"({best.result.final_recovery_factor:.2f} %)  ·  fərq {delta:.2f} %")

    def clear_runs(self):
        self.project.runs.clear()
        self.project._counter = 0
        self.update_comparison()
        self.refresh_tree()
        self.statusBar().showMessage("İşə salınmalar siyahısı təmizləndi.")

    # ═══════════════════════════════════════════════ layihə faylı (.imx)
    def save_project(self, include_snapshots: bool = True):
        suggested = self.project_path or f"layihe{FILE_EXTENSION}"
        path, _ = QFileDialog.getSaveFileName(
            self, "Layihəni yadda saxla", suggested,
            f"IMEX-2D layihəsi (*{FILE_EXTENSION})")
        if not path:
            return
        if not path.endswith(FILE_EXTENSION):
            path += FILE_EXTENSION
        self.project.geology_wells = self.geology_panel.wells()
        self.project.geology_method = self.geology_panel.method_text()
        self.project.geology_params = dict(
            power=self.geology_panel.power.value(),
            search_radius=self.geology_panel.search_radius.value(),
            range_a=self.geology_panel.range_.value(),
            sill_c=self.geology_panel.sill.value(),
            nugget_c0=self.geology_panel.nugget.value())
        try:
            self.serializer.save(self.project, path, include_snapshots)
        except Exception as exc:
            QMessageBox.critical(self, "Yadda saxlanmadı", str(exc))
            return
        self.project_path = path
        self._mark_clean()
        size = os.path.getsize(path) / 1024.0
        self.statusBar().showMessage(
            f"Yazıldı: {os.path.basename(path)}  ({size:.0f} KB, "
            f"{len(self.project.runs)} işə salınma)")

    def open_project(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Layihəni aç", "", f"IMEX-2D layihəsi (*{FILE_EXTENSION})")
        if not path:
            return
        try:
            project = self.serializer.load(path)
        except ProjectFileError as exc:
            QMessageBox.critical(self, "Fayl açılmadı", str(exc))
            return
        except Exception as exc:
            QMessageBox.critical(self, "Fayl açılmadı", f"Gözlənilməz xəta: {exc}")
            return

        self.project = project
        self.project_path = path
        self._ready = False
        try:
            self.geology_panel.load(project.geology_wells)
            if project.geology_method:
                self.geology_panel.method.setCurrentText(project.geology_method)
            params = project.geology_params
            if params:
                self.geology_panel.power.setValue(params.get("power", self.geology_panel.power.value()))
                self.geology_panel.search_radius.setValue(
                    params.get("search_radius", self.geology_panel.search_radius.value()))
                self.geology_panel.range_.setValue(params.get("range_a", self.geology_panel.range_.value()))
                self.geology_panel.sill.setValue(params.get("sill_c", self.geology_panel.sill.value()))
                self.geology_panel.nugget.setValue(params.get("nugget_c0", self.geology_panel.nugget.value()))
        finally:
            self._ready = True
        self._sync_geology_geometry()
        self._geology_model_from_wells = project.geological_models.get(
            "Quyu cədvəlindən geoloji model")
        self.geology_panel.mark_fresh()
        self._mark_clean()
        if project.reservoir_models:
            self.reservoir_model = list(project.reservoir_models.values())[-1]
            self._load_model_into_panels(self.reservoir_model)
            # `_load_model_into_panels` `well_panel.load(model.wells)` ilə
            # cədvəli TAM sıfırlayır (saxlanılmış rejim üçün) — bundan
            # sonra geologiya ilə yenidən uzlaşdırmaq lazımdır (in_model
            # olub faylda quyusu olmayanlar da sətir kimi görünsün).
            self._sync_geology_geometry()
        latest = project.latest_run()
        self.result = latest.result if latest else None
        if self.result and self.result.snapshots:
            self.slider.setEnabled(True)
            self.slider.setRange(0, len(self.result.snapshots) - 1)
            self.slider.setValue(len(self.result.snapshots) - 1)
            self.update_results()
        self.refresh_tree()
        self.update_scal_plot()
        self.update_pvt_plot()
        self.update_map()
        self.update_volume()
        self.update_comparison()
        self.statusBar().showMessage(
            f"Açıldı: {os.path.basename(path)}  ·  "
            f"{len(project.reservoir_models)} model, {len(project.runs)} işə salınma")

    def _load_model_into_panels(self, model):
        """Fayldan gələn modeli interfeys sahələrinə yazır.

        Yalnız panellərin dəstəklədiyi sadə parametrlər bərpa olunur;
        heterogen xassə xəritələri modeldə qalır, amma paneldəki
        "Homogen/Təsadüfi" seçimi ilə yenidən qurulmur.
        """
        self._ready = False
        try:
            grid, geometry = model.grid, model.geometry
            self.grid_panel.nx.setValue(grid.nx)
            self.grid_panel.ny.setValue(grid.ny)
            self.grid_panel.nz.setValue(grid.nz)
            self.grid_panel.dx.setValue(geometry.dx)
            self.grid_panel.dy.setValue(geometry.dy)
            self.grid_panel.thickness_mode.setCurrentIndex(0)   # DZ rejimi
            self.grid_panel.set_layer_thicknesses(geometry.dz)
            self.grid_panel.top_depth.setValue(geometry.top_depth)

            self.rock_panel.porosity.setValue(float(model.rock.porosity.values.mean()))
            # `model.rock.permx`/`model.fluids.*` HƏMİŞƏ mühərrik vahidindədir
            # (mD/cP) — vahid seçicini dəyər ilə BİRLİKDƏ defolta qaytarırıq
            # ki, panel köhnə (vahidsiz) davranışı ilə eyni ədədi göstərsin.
            self.rock_panel.permx_unit.setCurrentText("mD")
            self.rock_panel.permx.setValue(float(model.rock.permx.values.mean()))
            self.rock_panel.viscosity_unit.setCurrentText("cP")
            self.rock_panel.mu_w.setValue(model.fluids.water_viscosity)
            self.rock_panel.mu_o.setValue(model.fluids.oil_viscosity)
            self.rock_panel.bo.setValue(model.fluids.oil_fvf)

            scal = model.scal_parameters
            for widget, value in ((self.scal_panel.swc, scal.swc),
                                  (self.scal_panel.sor, scal.sor),
                                  (self.scal_panel.krw_end, scal.krw_end),
                                  (self.scal_panel.kro_end, scal.kro_end),
                                  (self.scal_panel.nw, scal.nw),
                                  (self.scal_panel.no, scal.no),
                                  (self.scal_panel.pc_entry,
                                   model.capillary_parameters.entry_pressure),
                                  (self.scal_panel.pc_lambda,
                                   model.capillary_parameters.lambda_exponent)):
                widget.setValue(value)

            ic = model.initial_conditions
            self.numerical_panel.initial_pressure_unit.setCurrentText("bar")
            self.numerical_panel.initial_pressure.setValue(ic.datum_pressure)
            self.numerical_panel.initial_sw.setValue(ic.water_saturation)
            self.numerical_panel.datum_depth.setValue(ic.datum_depth)
            self.numerical_panel.use_equilibration.setChecked(ic.use_equilibration)
            if ic.oil_water_contact is not None:
                self.numerical_panel.owc.setValue(ic.oil_water_contact)

            self.pvt_panel.enabled.setChecked(model.pvt_table is not None)
            self.well_panel.load(model.wells)
        finally:
            self._ready = True

    # ══════════════════════════════════════════════════════════ eksport
    def export_results(self):
        if not self.result or not self.result.series.time:
            QMessageBox.information(self, "Nəticə yoxdur", "Əvvəlcə modeli işə salın.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Nəticələri yaz", "results.csv",
                                              "CSV (*.csv)")
        if not path:
            return
        s = self.result.series
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["time_day", "qo_m3d", "qw_m3d", "qwinj_m3d",
                             "cum_oil_m3", "cum_wat_m3", "watercut_pct",
                             "avg_p_bar", "RF_pct"])
            for i in range(len(s.time)):
                writer.writerow([f"{s.time[i]:.4f}", f"{s.oil_rate[i]:.4f}",
                                 f"{s.water_rate[i]:.4f}",
                                 f"{s.water_injection_rate[i]:.4f}",
                                 f"{s.cumulative_oil[i]:.2f}",
                                 f"{s.cumulative_water[i]:.2f}",
                                 f"{s.water_cut[i]:.3f}",
                                 f"{s.average_pressure[i]:.3f}",
                                 f"{s.recovery_factor[i]:.4f}"])
        self.statusBar().showMessage(f"Yazıldı: {os.path.basename(path)}")

    def export_snapshot(self):
        if not self.result or self.reservoir_model is None:
            QMessageBox.information(self, "Nəticə yoxdur", "Əvvəlcə modeli işə salın.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Grid anını yaz", "grid.csv",
                                              "CSV (*.csv)")
        if not path:
            return
        model = self.reservoir_model
        snapshot = self.result.snapshots[self.slider.value()]
        sw = snapshot.water_saturation.ravel()
        pressure = snapshot.pressure.ravel()
        poro = model.rock.porosity.values
        permx = model.rock.permx.values
        region = model.regions.region_id.values
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["i", "j", "k", "x_m", "y_m", "poro", "permx_mD",
                             "region", "Sw", "P_bar"])
            for cell in range(model.ncell):
                i, j, k = model.grid.ijk(cell)
                writer.writerow([i, j, k, (i + .5) * model.geometry.dx,
                                 (j + .5) * model.geometry.dy,
                                 f"{poro[cell]:.4f}", f"{permx[cell]:.3f}",
                                 int(region[cell]), f"{sw[cell]:.5f}",
                                 f"{pressure[cell]:.3f}"])
        self.statusBar().showMessage(
            f"Yazıldı: {os.path.basename(path)}  (t = {snapshot.time:.0f} gün)")

    # ═══════════════════════════════════════════ Eclipse mübadiləsi
    def import_grdecl(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "GRDECL grid faylı", "",
            "Eclipse grid (*.grdecl *.GRDECL *.data *.DATA *.inc);;"
            "Bütün fayllar (*)")
        if not path:
            return

        report = DiagnosticReport()
        try:
            deck = read_grdecl(path, report)
            geology = GrdeclImporter().build(
                deck, report, name=os.path.basename(path))
        except GrdeclError as error:
            QMessageBox.critical(self, "Fayl oxunmadı", str(error))
            LOG.error("GRDECL oxunmadı: %s", error)
            return
        except Exception as error:
            QMessageBox.critical(self, "Fayl oxunmadı",
                                 f"Gözlənilməz xəta: {error}")
            LOG.exception("GRDECL oxunmadı")
            return

        for item in report.items:
            LOG.info("%s", item)

        self.imported_geology = geology
        self.project.add_geological_model(geology)
        self._apply_imported_geometry(geology)
        self.rebuild_model()

        summary = deck.summary()
        message = (f"{os.path.basename(path)}\n\n"
                   f"Grid: {summary['ölçü']}  ·  {summary['hüceyrə']} hüceyrə\n"
                   f"Massivlər: {summary['massivlər']}")
        if report.warnings:
            message += "\n\nDiqqət:\n" + "\n".join(
                f"• {item.message}" for item in report.warnings)
        QMessageBox.information(self, "GRDECL oxundu", message)
        self.statusBar().showMessage(
            f"Oxundu: {os.path.basename(path)} — {summary['hüceyrə']} hüceyrə")

    def _apply_imported_geometry(self, geology):
        """Fayldan gələn grid ölçülərini panelə yazır.

        Xassə xəritələri modeldə qalır; panel yalnız ölçüləri göstərir,
        çünki heterogen massivlər sürgülərlə ifadə oluna bilməz.
        """
        self._ready = False
        try:
            grid, geometry = geology.grid, geology.geometry
            self.grid_panel.nx.setValue(grid.nx)
            self.grid_panel.ny.setValue(grid.ny)
            self.grid_panel.nz.setValue(grid.nz)
            self.grid_panel.dx.setValue(geometry.dx)
            self.grid_panel.dy.setValue(geometry.dy)
            self.grid_panel.thickness_mode.setCurrentIndex(0)   # DZ rejimi
            self.grid_panel.set_layer_thicknesses(geometry.dz)
            self.grid_panel.top_depth.setValue(geometry.top_depth)
            default_wells = WELL_PATTERNS["Five-spot (1/4)"](grid)
            self.well_panel.load(default_wells)
            self.geology_panel.load(self._wells_to_geology_rows(default_wells, geometry))
        finally:
            self._ready = True
        self._sync_geology_geometry()

    def import_opm_case(self):
        """OPM Flow (Eclipse formatlı) halının nəticələrini idxal edir
        və öz 3D görüntümüzdə göstərir.

        Strateji qərar (bax A7_PLAN.md/jurnal): A7-nin öz üç fazalı
        Nyuton həlledicisi hələ açıq bir davamlılıq problemi daşıyır.
        Bunun əvəzinə: FİZİKANI OPM Flow-a həvalə edirik, öz güclü
        tərəfimizi (3D görüntü/analiz) onun NƏTİCƏLƏRİNİ göstərmək
        üçün işlədirik. Mövcud "3D görüntü" tabındakı bütün idarəetmə
        (vaxt slider-i, xassə seçimi) DƏYİŞMƏDƏN işləyir — yalnız
        `self.reservoir_model`/`self.result` doldurulur.
        """
        path, _ = QFileDialog.getOpenFileName(
            self, "OPM Flow halı (.EGRID seç)", "", "Eclipse grid (*.EGRID)")
        if not path:
            return
        case_root = path[:-6] if path.upper().endswith(".EGRID") else path

        try:
            from ..io.opm_import import build_display_model, load_opm_case
            case = load_opm_case(case_root)
        except Exception as error:
            QMessageBox.critical(self, "OPM halı yüklənmədi", str(error))
            LOG.exception("OPM Flow idxalı uğursuz oldu")
            return

        from ..simulation.results import SimulationResult, Snapshot
        self.reservoir_model = build_display_model(case)
        self.result = SimulationResult(
            model_name=self.reservoir_model.name,
            grid_shape=self.reservoir_model.grid.shape,
            snapshots=[Snapshot(time=s.time, pressure=s.pressure,
                               water_saturation=s.water_saturation)
                      for s in case.snapshots])
        self.volume_time.setMaximum(max(len(case.snapshots) - 1, 0))
        self.volume_time.setValue(len(case.snapshots) - 1)

        message = f"{len(case.snapshots)} addım idxal olundu."
        if case.warnings:
            message += "\n\nXəbərdarlıqlar:\n" + "\n".join(case.warnings)
        QMessageBox.information(self, "OPM halı idxal olundu", message)

        self.show_tab("3D görüntü")
        self.update_volume()

    def export_pdf_report(self):
        if self.reservoir_model is None:
            QMessageBox.information(self, "Model yoxdur",
                                    "Əvvəlcə model qurun.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "PDF hesabat yaz", "hesabat.pdf", "PDF (*.pdf)")
        if not path:
            return
        context = ReportContext(
            model=self.reservoir_model, result=self.result,
            mismatch=self.mismatch_report)
        try:
            ReportGenerator().write(context, path)
        except Exception as error:
            QMessageBox.critical(self, "Yazılmadı", str(error))
            LOG.exception("PDF hesabat yazılmadı")
            return
        size = os.path.getsize(path) / 1024.0
        QMessageBox.information(
            self, "Hesabat hazırdır",
            f"{os.path.basename(path)}  ({size:.0f} KB)")
        self.statusBar().showMessage(f"Yazıldı: {os.path.basename(path)}")

    def export_eclipse(self):
        if self.reservoir_model is None:
            QMessageBox.information(self, "Model yoxdur",
                                    "Əvvəlcə model qurun.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Eclipse deck yaz", "model.DATA",
            "Eclipse deck (*.DATA *.data)")
        if not path:
            return
        config = self.numerical_panel.simulation_config()
        try:
            EclipseDeckWriter(end_time=config.end_time,
                              report_steps=20).write(self.reservoir_model, path)
        except Exception as error:
            QMessageBox.critical(self, "Yazılmadı", str(error))
            LOG.exception("Eclipse deck yazılmadı")
            return
        size = os.path.getsize(path) / 1024.0
        QMessageBox.information(
            self, "Eclipse deck hazırdır",
            f"{os.path.basename(path)}  ({size:.0f} KB)\n\n"
            "Deck avtomatik yaradılıb. Hər simulyatorun öz "
            "xüsusiyyətləri var — işə salmazdan əvvəl yoxlayın.")
        self.statusBar().showMessage(f"Yazıldı: {os.path.basename(path)}")

    def show_version(self):
        QMessageBox.information(self, f"IMEX-2D v{VERSION}", summary())

    def show_about(self):
        QMessageBox.information(
            self, "Arxitektura",
            "IMEX-2D — qatlara ayrılmış rezervuar modelləşdirmə platforması.\n\n"
            "Qatlar:\n"
            "  domain       — model obyektləri (asılılıqsız)\n"
            "  interfaces   — provider müqavilələri\n"
            "  application  — layihə, iş axını, servislər\n"
            "  simulation   — hesablama mühərriki\n"
            "  rendering    — çəkmə (Qt-siz)\n"
            "  ui           — Qt interfeysi\n\n"
            "İş axını: Layihə → Geoloji model → Rezervuar modeli →\n"
            "Simulyasiya → Nəticələr\n\n"
            "PVT, kapilyar təzyiq və initialization üçün interfeyslər hazırdır,\n"
            "implementasiya növbəti mərhələdədir.")
