"""Günlük göstəricilər cədvəlinin Qt modeli — Seans 45.

NİYƏ `QTableWidget` YOX. 1500 günlük qaçış yataq üçün ~14 sütunla
~21 000 xana deməkdir, SPE1 kimi 10 illik qaçış isə 3650 sətir. Hər
xana üçün `QTableWidgetItem` yaratmaq obyekti seçəndə hiss olunan
gecikmə verir; `QAbstractTableModel` isə yalnız ekranda GÖRÜNƏN xanaları
soruşur. Rəqəmlər `reporting/daily.py`-dən gəlir — burada hesablama yoxdur.
"""

from __future__ import annotations

import math
from typing import List, Optional

from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt

from ..reporting import daily
from ..reporting.daily import DailyTable

#: Boş dəyər (THP-nin `nan`-ı — quyu axmır) — sıfır YAZILMIR.
BLANK = "—"


def format_value(column: str, value: float) -> str:
    """Sütunun mənasına görə dəqiqlik: debit/təzyiq 2, həcm 1 rəqəm."""
    if value is None or math.isnan(value):
        return BLANK
    if column == daily.DAY:
        return f"{value:.0f}" if abs(value - round(value)) < 1e-9 else f"{value:.2f}"
    if column.startswith("kum_"):
        return f"{value:,.1f}".replace(",", " ")
    if column == daily.GOR:
        return f"{value:.1f}"
    return f"{value:.2f}"


def format_row(row: dict) -> str:
    """Seçilmiş günün qısa xülasəsi — cədvəlin üstündəki sətir üçün."""
    parts = [f"{name} = {format_value(name, value)}"
             for name, value in row.items() if name != daily.DAY]
    return "     ".join(parts)


class DailyTableModel(QAbstractTableModel):
    """Bir obyektin (yataq və ya quyu) günlük cədvəli."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._table: DailyTable = DailyTable()
        self._target: str = daily.FIELD
        self._headers: List[str] = []

    # ------------------------------------------------------------ məlumat
    def set_table(self, table: DailyTable, target: Optional[str] = None) -> None:
        self.beginResetModel()
        self._table = table
        targets = table.targets()
        self._target = target if target in targets else daily.FIELD
        self._headers = [daily.DAY] + list(table.columns(self._target))
        self.endResetModel()

    @property
    def target(self) -> str:
        return self._target

    @property
    def headers(self) -> List[str]:
        return list(self._headers)

    def _value(self, row: int, column: int) -> float:
        name = self._headers[column]
        if name == daily.DAY:
            return float(self._table.days[row])
        return float(self._table.columns(self._target)[name][row])

    # ----------------------------------------------------------- Qt API
    def rowCount(self, parent=QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._table)

    def columnCount(self, parent=QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._headers)

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        if role == Qt.DisplayRole:
            name = self._headers[index.column()]
            return format_value(name, self._value(index.row(), index.column()))
        if role == Qt.TextAlignmentRole:
            return int(Qt.AlignRight | Qt.AlignVCenter)
        return None

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return None
        if orientation == Qt.Horizontal and 0 <= section < len(self._headers):
            # "q_neft [m³/gün]" → iki sətir: ad və vahid (sütunlar dar qalsın)
            return self._headers[section].replace(" [", "\n[", 1)
        return None
