"""pytest üçün yol qurulması — layihə kökü sys.path-a əlavə olunur."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# PyQt5 widget-lərini ekran/displeysiz mühitdə (CI, bu mühit) qurmağa
# imkan verir — YALNIZ artıq təyin edilməyibsə (istifadəçinin öz mühiti
# üstünlük təşkil edir). Yalnız `test_ui_units.py` real QApplication
# yaradır; qalan UI testləri (`test_ui_static.py`/`test_ui_wiring.py`)
# AST təhlili aparır və Qt-ni heç işə salmır — bu dəyişiklikdən təsirlənmir.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
