"""Geologiya bölməsinin kiçik 2D xəritəsi.

Yeganə işi: grid düzbucaqlısı + quyu nöqtələrini çəkmək. Heç bir hesablama
aparmır — yalnız `imex2d.domain.geology.GeologicalWell` siyahısını göstərir.
Koordinat səhvini (məs. X və Y qarışdırılıb) rəqəmə baxaraq tutmaq çətindir;
nöqtə sərhəddən kənarda görünəndə dərhal aydın olur.
"""

from __future__ import annotations

from typing import List, Optional

from PyQt5.QtCore import QRectF, Qt
from PyQt5.QtGui import QBrush, QColor, QPainter, QPen
from PyQt5.QtWidgets import QWidget

from ..domain.geology import GeologicalWell
from ..rendering.theme import PALETTE

_MARGIN = 14
_RADIUS = 4.5

_COLOR_IN_MODEL = QColor("#3fa7ff")
_COLOR_DATA_ONLY = QColor("#9aa0a6")
_COLOR_OUT_OF_BOUNDS = QColor("#e05252")
_COLOR_SELECTED = QColor("#ffb020")


class GeologyMapWidget(QWidget):
    """Grid sərhədi + quyu nöqtələri. Statik göstərici, siçan girişi yoxdur."""

    def __init__(self):
        super().__init__()
        self.setMinimumHeight(170)
        self._wells: List[GeologicalWell] = []
        self._x_max = 0.0
        self._y_max = 0.0
        self._selected: Optional[str] = None

    def set_data(self, wells: List[GeologicalWell], x_max: float, y_max: float,
                selected: Optional[str] = None):
        self._wells = list(wells)
        self._x_max = max(x_max, 1e-9)
        self._y_max = max(y_max, 1e-9)
        self._selected = selected
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.fillRect(self.rect(), Qt.transparent)

        w, h = self.width() - 2 * _MARGIN, self.height() - 2 * _MARGIN
        if w <= 0 or h <= 0:
            return
        scale = min(w / self._x_max, h / self._y_max)
        draw_w, draw_h = self._x_max * scale, self._y_max * scale
        ox = _MARGIN + (w - draw_w) / 2.0
        oy = _MARGIN + (h - draw_h) / 2.0

        def to_screen(x, y):
            # Y ekranda aşağı-yuxarı əksdir (0 aşağıda olsun)
            return ox + x * scale, oy + draw_h - y * scale

        painter.setPen(QPen(Qt.gray, 1))
        painter.setBrush(QBrush(Qt.NoBrush))
        painter.drawRect(QRectF(ox, oy, draw_w, draw_h))

        for well in self._wells:
            in_bounds = 0.0 <= well.x <= self._x_max and 0.0 <= well.y <= self._y_max
            px = min(max(well.x, 0.0), self._x_max)
            py = min(max(well.y, 0.0), self._y_max)
            sx, sy = to_screen(px, py)

            if not in_bounds:
                color = _COLOR_OUT_OF_BOUNDS
            elif well.in_model:
                color = _COLOR_IN_MODEL
            else:
                color = _COLOR_DATA_ONLY

            pen_color = _COLOR_SELECTED if well.name == self._selected else color
            painter.setPen(QPen(pen_color, 1.5))
            painter.setBrush(QBrush(color if well.in_model or not in_bounds
                                    else Qt.transparent))
            painter.drawEllipse(sx - _RADIUS, sy - _RADIUS, 2 * _RADIUS, 2 * _RADIUS)
            painter.setPen(QPen(QColor(PALETTE.text), 1))
            painter.drawText(sx + _RADIUS + 2, sy + 3, well.name)
