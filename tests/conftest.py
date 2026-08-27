"""pytest üçün yol qurulması — layihə kökü sys.path-a əlavə olunur."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
