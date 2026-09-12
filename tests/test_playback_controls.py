"""B6-b — oynatma idarəsi: sürət, kadr-kadr irəli/geri, dövrə.

ƏVVƏLKİ VƏZİYYƏT. Zaman slider-i və `▶ Oynat` düyməsi ARTIQ vardı
(`main_window`-da `QTimer` 140 ms ilə). Slider həm 2D xəritəni, həm də
3D VTK görüntüsünü yeniləyir — yəni cəbhənin hərəkəti onsuz da
izlənilirdi. Çatışmayan: **sürət tənzimi**, **kadr-kadr addımlama** və
sonda **dayanma seçimi**.

NİYƏ MƏNTİQ AYRICA MODULDADIR (`ui/playback.py`). İlk versiyada testlər
`MainWindow`-u qururdu və `pytest` FATAL XƏTA ilə çökürdü — pəncərənin
qurulması VTK səhnəsi və bir neçə matplotlib kanvası yaradır. Mövcud
testlərin heç biri təsadüfən `MainWindow` qurmur. Ona görə saf məntiq
Qt-dən ayrıldı; burada məhz o yoxlanılır.

GERİYƏ UYĞUNLUQ ŞƏRTİ: defolt (1× sürət, dövrə açıq) davranış
B6-b-dən ƏVVƏLKİ ilə eyni olmalıdır — 140 ms və sonda əvvələ qayıtma.
"""

from __future__ import annotations

import pytest

from imex2d.ui.playback import (MIN_INTERVAL_MS, PLAYBACK_BASE_MS,
                                PLAYBACK_DEFAULT_INDEX, PLAYBACK_SPEEDS,
                                advance, interval_ms, step_value)


# ═══════════════════════ sürət cədvəli ═══════════════════════════════

def test_default_speed_is_one_times():
    """Defolt seçim 1× olmalıdır — köhnə davranış dəyişmir."""
    label, factor = PLAYBACK_SPEEDS[PLAYBACK_DEFAULT_INDEX]
    assert label == "1×"
    assert factor == 1.0


def test_base_interval_is_unchanged():
    """140 ms — B6-b-dən əvvəlki `QTimer.start(140)` dəyəri."""
    assert PLAYBACK_BASE_MS == 140


def test_speed_table_is_ordered_and_positive():
    factors = [factor for _, factor in PLAYBACK_SPEEDS]
    assert all(f > 0 for f in factors)
    assert factors == sorted(factors), factors
    assert len(set(factors)) == len(factors), "təkrarlanan sürət"


# ═══════════════════════ interval hesabı ═════════════════════════════

def test_one_times_gives_exactly_the_old_interval():
    assert interval_ms(1.0) == 140


def test_faster_speed_means_shorter_interval():
    intervals = [interval_ms(factor) for _, factor in PLAYBACK_SPEEDS]
    assert intervals == sorted(intervals, reverse=True), intervals


def test_double_speed_halves_the_interval():
    assert interval_ms(2.0) == 70
    assert interval_ms(0.5) == 280


def test_interval_never_drops_below_the_floor():
    """Çox yüksək sürətdə Qt hadisə növbəsi boğulmamalıdır."""
    assert interval_ms(1000.0) == MIN_INTERVAL_MS
    assert interval_ms(1e9) >= MIN_INTERVAL_MS


@pytest.mark.parametrize("bad", [0.0, -2.0, None, "sürətli"])
def test_invalid_speed_falls_back_to_one_times(bad):
    """Yararsız dəyər taymeri partlatmamalıdır."""
    assert interval_ms(bad) == 140


# ═══════════════════════ kadr-kadr addımlama ═════════════════════════

def test_step_moves_one_frame_each_way():
    assert step_value(0, +1, 0, 4) == 1
    assert step_value(3, +1, 0, 4) == 4
    assert step_value(3, -1, 0, 4) == 2


def test_step_clamps_at_both_ends():
    """Addımlama DÖVRƏ ETMİR — sərhəddə dayanır.

    Dövrə yalnız avtomatik oynatmaya aiddir; əl ilə addımlayan
    istifadəçi son kadrdan birdən başlanğıca atılmağı gözləmir.
    """
    assert step_value(4, +1, 0, 4) == 4
    assert step_value(0, -1, 0, 4) == 0
    assert step_value(4, +10, 0, 4) == 4
    assert step_value(0, -10, 0, 4) == 0


def test_step_handles_an_empty_range():
    """Nəticə bir kadrlıdırsa (və ya boşdursa) çökməməlidir."""
    assert step_value(0, +1, 0, 0) == 0
    assert step_value(0, +1, 0, -1) == 0


# ═══════════════════════ dövrə / sonda dayanma ═══════════════════════

def test_advance_moves_forward_in_the_middle():
    assert advance(1, 4, loop=True) == (2, True)
    assert advance(1, 4, loop=False) == (2, True)


def test_loop_wraps_to_the_first_frame():
    """DEFOLT davranış — B6-b-dən əvvəlki ilə eyni."""
    assert advance(4, 4, loop=True) == (0, True)


def test_without_loop_it_stops_on_the_last_frame():
    value, keep_playing = advance(4, 4, loop=False)
    assert value == 4, "son kadrda qalmalıdır, əvvələ atılmamalıdır"
    assert keep_playing is False


def test_single_frame_result_does_not_spin():
    """Tək kadrlı nəticədə dövrəsiz rejim dərhal dayanmalıdır."""
    assert advance(0, 0, loop=False) == (0, False)
    assert advance(0, 0, loop=True) == (0, True)


# ═══════════════════════ UI bağlantısı (yüngül) ══════════════════════

def test_main_window_uses_the_shared_logic():
    """`main_window` öz nüsxəsini saxlamamalıdır — TƏKRAR KOD OLMASIN.

    Pəncərənin ÖZÜ qurulmur: `MainWindow()` qurmaq VTK səhnəsi və
    matplotlib kanvasları yaradır və `pytest`-i çökdürür (ölçüldü).
    Burada yalnız modulun düzgün mənbədən idxal etdiyi yoxlanılır.
    """
    pytest.importorskip("PyQt5.QtWidgets")
    from imex2d.ui import main_window as MW
    from imex2d.ui import playback

    assert MW.PLAYBACK_SPEEDS is playback.PLAYBACK_SPEEDS
    assert MW.interval_ms is playback.interval_ms
    assert MW.advance is playback.advance
    assert MW.step_value is playback.step_value


def test_main_window_exposes_the_new_controls():
    """Yeni widget-lər sinifdə QURULUR — mənbə mətnindən yoxlanılır."""
    pytest.importorskip("PyQt5.QtWidgets")
    import inspect

    from imex2d.ui import main_window as MW

    source = inspect.getsource(MW.MainWindow)
    for name in ("step_back_button", "step_forward_button",
                 "speed_box", "loop_box", "_set_playback_enabled",
                 "step_frame", "_set_playing"):
        assert name in source, name
