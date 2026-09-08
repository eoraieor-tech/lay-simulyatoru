"""Analitik həllər — validasiya üçün. Mövcud koddan köçürülüb.

Simulyasiya mühərrikindən ayrı saxlanılır: bu, verifikasiya alətidir,
hesablama zəncirinin hissəsi deyil.

Bu modul ÖLÇÜ CİHAZIDIR — mühərriki onunla ölçürük. Ona görə burada
"təxminən düzgün" qəbul edilmir: sınıq etalon fizika səhvindən daha
təhlükəlidir, çünki gələcək səhvləri maskalayır. Doğruluğun meyarı
Bakli-Leverett həllində DƏQİQ ödənən kütlə eyniliyidir:

    ∫₀^∞ (Sw(x) − Swi) dx  =  q·t / (φ·A)  =  velocity · t

`tests/test_analytical_bl.py` məhz bunu yoxlayır.
"""

from __future__ import annotations
from dataclasses import dataclass

import numpy as np

from ..domain.scal import CoreyParameters


@dataclass
class BuckleyLeverettSolution:
    distance: np.ndarray
    water_saturation: np.ndarray
    shock_saturation: float
    front_position: float
    length: float = 0.0
    breakthrough: bool = False
    monotonicity_corrections: int = 0
    notes: tuple = ()


def buckley_leverett(scal: CoreyParameters, mu_w: float, mu_o: float,
                     porosity: float, total_rate: float, area: float,
                     time: float, sw_initial: float = None,
                     length: float = None) -> BuckleyLeverettSolution:
    """Bakli-Leverett sıxışdırma profili.

    `length` — modelin fiziki uzunluğu, m. Verilməzsə `2·x_front`
    işlədilir (YALNIZ göstərmə üçün; müqayisə həmişə `length` ilə
    aparılmalıdır).

    Qaytarılan profil dörd hissədən ibarətdir və ciddi artan `distance`
    ilə ARTIQ SIRALANMIŞ qurulur:

        1) x = 0 … x(1−Sor)   : Sw = 1 − Sor   (tam süpürülmüş zona)
        2) rarefaksiya yelpiyi: x(s) = velocity·(df/ds)(s)·t
        3) x_front            : ŞAQULİ sıçrayış sw_shock → swi
        4) x_front⁺ … length  : Sw = swi       (toxunulmamış zona)
    """
    if not time > 0:
        raise ValueError(f"Zaman müsbət olmalıdır (verildi: {time}).")
    if not total_rate > 0:
        raise ValueError(f"Vurulma debiti müsbət olmalıdır (verildi: {total_rate}).")
    if not area > 0:
        raise ValueError(f"En kəsik sahəsi müsbət olmalıdır (verildi: {area}).")
    if not porosity > 0:
        raise ValueError(f"Məsaməlik müsbət olmalıdır (verildi: {porosity}).")
    if not mu_w > 0 or not mu_o > 0:
        raise ValueError(
            f"Lözlüklər müsbət olmalıdır (su: {mu_w}, neft: {mu_o}).")
    if length is not None and not length > 0:
        raise ValueError(f"Model uzunluğu müsbət olmalıdır (verildi: {length}).")

    swi = scal.swc if sw_initial is None else float(sw_initial)
    if not (scal.swc <= swi < 1.0 - scal.sor):
        raise ValueError(
            f"İlkin su doyumluluğu {swi} intervaldan kənardadır: "
            f"Swc = {scal.swc} ≤ Swi < 1 − Sor = {1.0 - scal.sor}.")

    s = np.linspace(swi + 1e-6, 1.0 - scal.sor - 1e-6, 2000)

    lw = scal.krw(s) / mu_w
    lo = scal.kro(s) / mu_o
    f = lw / np.maximum(lw + lo, 1e-30)
    lw_i = scal.krw(np.array([swi]))[0] / mu_w
    lo_i = scal.kro(np.array([swi]))[0] / mu_o
    f_i = lw_i / max(lw_i + lo_i, 1e-30)

    chord = (f - f_i) / (s - swi)
    k = int(np.argmax(chord))
    sw_shock, slope_shock = s[k], chord[k]

    dfds = np.gradient(f, s)
    velocity = total_rate / (area * porosity)
    x_front = velocity * slope_shock * time

    fan = s >= sw_shock
    s_prof = s[fan]
    x_raw = velocity * dfds[fan] * time

    # Rarefaksiyada df/ds fiziki olaraq monotondur, yəni bu sətir HEÇ NƏ
    # etməməlidir. Səssiz düzəliş qadağandır: neçə nöqtənin zorla
    # dəyişdirildiyini sayırıq və nəticə ilə birlikdə qaytarırıq —
    # sıfırdan böyükdürsə, SCAL parametrləri qeyri-adidir.
    x_prof = np.maximum.accumulate(x_raw[::-1])[::-1]
    corrections = int(np.count_nonzero(x_prof != x_raw))
    notes = []
    if corrections:
        notes.append(
            f"Rarefaksiyada {corrections} nöqtə monotonluq üçün zorla "
            "düzəldildi — SCAL parametrləri qeyri-adidir (df/ds monoton deyil).")

    # ── (2) Yelpik: s-ə görə artan massiv x-ə görə AZALANDIR → çeviririk.
    x_fan, s_fan = x_prof[::-1], s_prof[::-1]

    # Şok nöqtəsi (3) dəqiq `x_front` ilə ayrıca qoyulur; yelpiyin ondan
    # kənara düşən son nöqtələri ədədi `np.gradient` ilə Welge chord-u
    # arasındakı fərqdən yaranır və atılır. `x = 0` nöqtəsi də (1)-də
    # ayrıca verilir.
    inside = (x_fan > 0.0) & (x_fan < x_front)
    x_fan, s_fan = x_fan[inside], s_fan[inside]

    # Ciddi artanlıq: yelpikdə eyni x-ə düşən qonşu nöqtələr yalnız float
    # yuvarlaqlaşmasından yaranır — birincisini saxlayırıq.
    if x_fan.size:
        keep = np.concatenate([[True], np.diff(x_fan) > 0.0])
        x_fan, s_fan = x_fan[keep], s_fan[keep]

    # ── (3) ŞOK: fiziki qalınlığı SIFIR olan şaquli sıçrayış.
    # `np.nextafter` növbəti təmsil oluna bilən float-u verir: massiv
    # ciddi artan qalır (np.interp və matplotlib bunu tələb edir), amma
    # `x_front · 1.001` kimi saxta keçid zonası YARANMIR.
    x_shock_plus = np.nextafter(x_front, np.inf)

    # ── (4) Toxunulmamış zona. Uzunluq verilməyibsə 2·x_front (göstərmə).
    model_length = 2.0 * x_front if length is None else float(length)
    x_tail = max(model_length, float(np.nextafter(x_shock_plus, np.inf)))

    x = np.concatenate([[0.0], x_fan, [x_front, x_shock_plus, x_tail]])
    sw = np.concatenate([[1.0 - scal.sor], s_fan, [sw_shock, swi, swi]])

    breakthrough = x_front >= model_length
    if breakthrough:
        notes.append(
            f"Cəbhə modeldən çıxıb: x_front = {x_front:.1f} m ≥ "
            f"uzunluq = {model_length:.1f} m — cəbhə mövqeyinə görə "
            "müqayisə mənasızdır.")

    return BuckleyLeverettSolution(
        distance=x, water_saturation=sw,
        shock_saturation=float(sw_shock), front_position=float(x_front),
        length=float(model_length), breakthrough=bool(breakthrough),
        monotonicity_corrections=corrections, notes=tuple(notes))
