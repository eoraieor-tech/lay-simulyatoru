"""Qt stil vərəqi — rendering palitrasından qurulur (tək rəng mənbəyi)."""

from __future__ import annotations

from ..rendering.theme import PALETTE as P


def stylesheet() -> str:
    return f"""
QMainWindow, QWidget {{ background: {P.background}; color: {P.text}; }}
QLabel {{ color: {P.text}; }}
QGroupBox {{ border: 1px solid {P.line}; border-radius: 4px; margin-top: 14px;
    padding: 10px 8px 8px 8px; background: {P.panel}; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 5px;
    color: {P.text_dim}; font-size: 10px; letter-spacing: 1px; }}
QToolBox::tab {{ background: {P.panel_alt}; border: 1px solid {P.line};
    border-radius: 3px; color: {P.text}; padding: 7px; font-weight: 600; }}
QToolBox::tab:selected {{ background: {P.accent}; color: #0B1116; }}
QSpinBox, QDoubleSpinBox, QComboBox, QTextEdit, QTableWidget, QLineEdit {{
    background: {P.panel_alt}; border: 1px solid {P.line}; border-radius: 3px;
    padding: 3px 5px; color: {P.text}; selection-background-color: {P.accent}; }}
QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{ border: 1px solid {P.accent}; }}
QComboBox QAbstractItemView {{ background: {P.panel_alt}; color: {P.text};
    selection-background-color: {P.accent}; }}
QPushButton {{ background: {P.panel_alt}; border: 1px solid {P.line};
    border-radius: 3px; padding: 6px 12px; color: {P.text}; }}
QPushButton:hover {{ border: 1px solid {P.accent}; }}
QPushButton:disabled {{ color: {P.text_dim}; border-color: #2A343D; }}
QPushButton#run {{ background: {P.accent}; color: #0B1116; font-weight: 700;
    letter-spacing: 1px; border: none; padding: 10px; }}
QPushButton#run:hover {{ background: #63B2E3; }}
QPushButton#stop {{ background: {P.danger}; color: #FFF; border: none; font-weight: 700; }}
QTabWidget::pane {{ border: 1px solid {P.line}; background: {P.panel}; }}
QTabBar::tab {{ background: {P.background}; color: {P.text_dim}; padding: 8px 18px;
    border: 1px solid {P.line}; border-bottom: none; }}
QTabBar::tab:selected {{ background: {P.panel}; color: {P.text};
    border-top: 2px solid {P.accent}; }}
QHeaderView::section {{ background: {P.panel_alt}; color: {P.text_dim}; border: none;
    border-right: 1px solid {P.line}; border-bottom: 1px solid {P.line}; padding: 5px; }}
QTableWidget {{ gridline-color: {P.line}; }}
QTreeWidget {{ background: {P.panel}; border: 1px solid {P.line}; color: {P.text}; }}
QProgressBar {{ background: {P.panel_alt}; border: 1px solid {P.line};
    border-radius: 3px; text-align: center; color: {P.text}; height: 18px; }}
QProgressBar::chunk {{ background: {P.accent}; }}
QSlider::groove:horizontal {{ background: {P.panel_alt}; height: 5px; border-radius: 2px; }}
QSlider::handle:horizontal {{ background: {P.accent}; width: 13px; margin: -5px 0;
    border-radius: 6px; }}
QStatusBar {{ background: {P.panel}; color: {P.text_dim}; border-top: 1px solid {P.line}; }}
QMenuBar {{ background: {P.panel}; color: {P.text}; }}
QMenuBar::item:selected {{ background: {P.accent}; color: #0B1116; }}
QMenu {{ background: {P.panel_alt}; color: {P.text}; border: 1px solid {P.line}; }}
QMenu::item:selected {{ background: {P.accent}; color: #0B1116; }}
"""
