"""Oynatma məntiqi — B6-b. **Qt-dən ASILI DEYİL.**

NİYƏ AYRICA MODUL. Bu məntiq əvvəl birbaşa `main_window`-un içində idi
və orada test oluna bilmirdi: `MainWindow`-un qurulması VTK səhnəsi və
bir neçə matplotlib kanvası yaradır, testdə isə bu, qaçırıcını
ÇÖKDÜRÜR (ölçüldü — `pytest` fatal xəta ilə dayanır). Layihənin öz
qaydası da bunu tələb edir: "UI yalnız orkestrasiyadır" (bax
`ARCHITECTURE.md` → 1.3).

Burada YALNIZ saf funksiyalar var: kadr arifmetikası və interval
hesabı. `main_window` onları çağırır, taymeri isə özü idarə edir.
"""

from __future__ import annotations

from typing import Sequence, Tuple

#: Oynatma sürətləri — `(etiket, əmsal)`. Əmsal kadr intervalını BÖLÜR,
#: yəni 2× iki dəfə tez oynadır.
PLAYBACK_SPEEDS: Sequence[Tuple[str, float]] = (
    ("0.25×", 0.25), ("0.5×", 0.5), ("1×", 1.0), ("2×", 2.0), ("4×", 4.0))

#: Defolt seçim — 1×. B6-b-dən ƏVVƏLKİ davranışı olduğu kimi saxlayır.
PLAYBACK_DEFAULT_INDEX = 2

#: 1× sürətdə bir kadrın müddəti, ms. Bu dəyər DƏYİŞMƏDİ.
PLAYBACK_BASE_MS = 140

#: Taymerin aşağı həddi. Çox yüksək sürətdə interval sıfıra yaxınlaşsa,
#: Qt hadisə növbəsi boğulur və interfeys cavab verməz olur.
MIN_INTERVAL_MS = 10


def interval_ms(factor: float) -> int:
    """Sürət əmsalından taymer intervalı, ms.

    `factor = 1.0` → tam `PLAYBACK_BASE_MS` (140 ms), yəni köhnə
    davranış bitə-bit qorunur.
    """
    try:
        value = float(factor)
    except (TypeError, ValueError):
        value = 1.0
    if value <= 0.0:
        value = 1.0
    return max(int(round(PLAYBACK_BASE_MS / value)), MIN_INTERVAL_MS)


def step_value(current: int, delta: int, minimum: int, maximum: int) -> int:
    """Kadr-kadr addımlama — sərhəddə DAYANIR, dövrə etmir.

    Dövrə yalnız avtomatik oynatmaya aiddir (`advance`). Əl ilə
    addımlayan istifadəçi son kadrdan birdən başlanğıca atılmağı
    gözləmir.
    """
    if maximum < minimum:
        return minimum
    return min(max(int(current) + int(delta), int(minimum)), int(maximum))


def advance(current: int, maximum: int, loop: bool) -> Tuple[int, bool]:
    """Avtomatik oynatmanın növbəti kadrı.

    Qaytarır: `(yeni_kadr, oynatma_davam_etsin)`.

    Sonda:
      * `loop=True`  → əvvələ qayıdır, oynatma davam edir (DEFOLT,
        B6-b-dən əvvəlki davranış);
      * `loop=False` → son kadrda qalır və oynatma dayanır.
    """
    value = int(current) + 1
    if value > int(maximum):
        if loop:
            return 0, True
        return int(maximum), False
    return value, True
