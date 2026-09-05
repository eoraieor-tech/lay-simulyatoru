"""Giriş nöqtəsi — COMPOSITION ROOT.

Bütün asılılıqlar YALNIZ burada bağlanır. Aşağıdakı qatların heç biri
konkret implementasiyanı özü seçmir. PVT modulu yazılanda dəyişəcək
yeganə fayl budur.

    pip install PyQt5 matplotlib numpy scipy
    python app.py
"""

from __future__ import annotations

import logging
import os
import sys

import matplotlib
matplotlib.use("Qt5Agg")

from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import QApplication

from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.project import Project
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder
from imex2d.application.simulation_service import ModelAwareSimulationService
from imex2d.domain.scal import CoreyParameters
from imex2d.logging_setup import configure as configure_logging
from imex2d.logging_setup import get_logger
from imex2d.simulation.linear_solver import ScipyCgIluSolver
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter
from imex2d.ui.main_window import MainWindow
from imex2d.ui.style import stylesheet


def main():
    log_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "logs", "imex2d.log")
    configure_logging(level=logging.INFO, log_file=log_file)
    logger = get_logger()
    from imex2d.version import VERSION
    logger.info("IMEX-2D v%s başladıldı  ·  log faylı: %s", VERSION, log_file)

    app = QApplication(sys.argv)
    app.setStyleSheet(stylesheet())
    app.setFont(QFont("Segoe UI" if sys.platform.startswith("win") else "Sans", 9))

    service = ModelAwareSimulationService(
        relperm_provider=CoreyRelativePermeabilityAdapter(CoreyParameters()),
        linear_solver=ScipyCgIluSolver(),
    )
    window = MainWindow(project=Project("IMEX-2D layihəsi"),
                        service=service,
                        geology_builder=SyntheticGeologicalModelBuilder(),
                        model_builder=ReservoirModelBuilder())
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
