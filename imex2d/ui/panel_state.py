"""Panellərin vəziyyətinin saxlanması və bərpası — Seans 46.

NİYƏ LAZIMDIR. «MODELİ İŞƏ SAL» modeli hər dəfə PANELLƏRDƏN yenidən
qurur (`MainWindow.rebuild_model`). Layihə açılanda isə panellərin
yalnız bir hissəsi doldurulurdu (`_load_model_into_panels`): simulyasiya
müddəti, Maks. Δt, mühərrik, PVT parametrləri, qaz SCAL və s. defoltda
qalırdı — sabah açılan layihə dünənki hesablamanı TƏKRARLAMIRDI.

NECƏ. Hər panelin sadə sahələri (spin, seçici, qutu, mətn) atribut adı
ilə oxunur və lüğətə yazılır — sahə siyahısı ƏL İLƏ yazılmır, ona görə
panelə yeni sahə əlavə olunanda o da avtomatik saxlanılır. Cədvəllər
(quyular, geologiya, faultlar, SWOF) bura DAXİL DEYİL: onlar modeldə və
layihədə artıq saxlanılır və öz `load()` metodları ilə bərpa olunur.

BƏRPA SIRASI VACİBDİR: əvvəl seçicilər, sonra qutular, sonra ədədlər.
Vahid seçicisi dəyişəndə spin-in göstərdiyi ədəd yeni vahidə çevrilir
(`panels._bind_unit_aware_spins`) — ədəd seçicidən SONRA yazılmalıdır.
Bəzi sahələrin diapazonu başqa sahədən asılıdır, ona görə hamısı İKİ
dəfə keçilir.
"""

from __future__ import annotations

from typing import Any, Dict

from PyQt5.QtWidgets import (QAbstractSpinBox, QCheckBox, QComboBox,
                             QDoubleSpinBox, QLineEdit, QSpinBox)

#: Bərpa neçə dəfə keçilir — asılı diapazonlar üçün (bax modulun şərhi).
RESTORE_PASSES = 2


def _json_scalar(value) -> bool:
    return value is None or isinstance(value, (str, int, float, bool))


def _widgets(panel):
    """Panelin İCTİMAİ sadə sahələri — `(ad, widget)` cütləri."""
    for name, widget in vars(panel).items():
        if name.startswith("_"):
            continue
        if isinstance(widget, (QComboBox, QCheckBox, QSpinBox,
                               QDoubleSpinBox, QLineEdit)):
            yield name, widget


def capture(panel) -> Dict[str, Any]:
    """Panelin sadə sahələrinin dəyərləri — JSON-a yazıla bilən lüğət."""
    state: Dict[str, Any] = {}
    for name, widget in _widgets(panel):
        if isinstance(widget, QComboBox):
            data = widget.currentData()
            state[name] = {"text": widget.currentText(),
                           "data": data if _json_scalar(data) else None}
        elif isinstance(widget, QCheckBox):
            state[name] = widget.isChecked()
        elif isinstance(widget, QAbstractSpinBox):
            state[name] = widget.value()
        elif isinstance(widget, QLineEdit):
            state[name] = widget.text()
    return state


def _set_combo(widget: QComboBox, value) -> None:
    if not isinstance(value, dict):
        return
    index = -1
    if value.get("data") is not None:
        index = widget.findData(value["data"])
    if index < 0:
        index = widget.findText(value.get("text", ""))
    if index >= 0 and index != widget.currentIndex():
        widget.setCurrentIndex(index)


def restore(panel, state: Dict[str, Any]) -> None:
    """`capture`-in tərsi. Naməlum açar (köhnə/yeni versiya) atılır.

    Siqnallar BLOKLANMIR: panelin daxili vəziyyəti (sahələrin aktivliyi,
    qat cədvəlinin sətirləri, vahid çevirməsi) məhz onlarla yenilənir.
    """
    if not state:
        return
    widgets = dict(_widgets(panel))
    for _ in range(RESTORE_PASSES):
        for kinds in ((QComboBox,), (QCheckBox,),
                      (QSpinBox, QDoubleSpinBox, QLineEdit)):
            for name, value in state.items():
                widget = widgets.get(name)
                if widget is None or not isinstance(widget, kinds):
                    continue
                if isinstance(widget, QComboBox):
                    _set_combo(widget, value)
                elif isinstance(widget, QCheckBox):
                    if isinstance(value, bool) and widget.isChecked() != value:
                        widget.setChecked(value)
                elif isinstance(widget, QAbstractSpinBox):
                    if isinstance(value, (int, float)) and not isinstance(value, bool):
                        widget.setValue(value)
                elif isinstance(widget, QLineEdit) and isinstance(value, str):
                    widget.setText(value)
