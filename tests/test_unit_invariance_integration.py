"""Phase 1 — inteqrasiya: vahid çevirməsi simulyasiya nəticəsini dəyişmir.

Tapşırıq: "Case A: psi/ft/mD/stb-day, Case B: SI/metrik — nəticələr FİZİKİ
CƏHƏTDƏN EYNİ olmalıdır."

DÜRÜST ƏHATƏ: mühərrikin giriş boru xətti (`ReservoirModelBuilder`,
`one_dimensional_model`) HƏLƏ vahid-etiketli GİRİŞ qəbul etmir (bax
Phase 1 hesabatı, "qalan iş") — yalnız bar/m/mD/m³/gün ədədini gözləyir.
Ona görə bu test EYNİ fiziki ssenarini psi/ft/stb-day kimi İFADƏ EDİR,
YENİ `unit_conversions.py` ilə mühərrikin vahidlərinə ÇEVİRİR, və
nəticənin (A) birbaşa mühərrik-vahidli girişlə qurulan modeldən FƏRQLİ
OLMADIĞINI göstərir — yəni çevirmə qatı sədaqətli round-trip verir və
bunun simulyasiya nəticəsinə TƏSİRİ YOXDUR (üzən nöqtə itkisindən başqa).
"""

from __future__ import annotations

import numpy as np
from helpers import default_scal, make_service, one_dimensional_model, short_config

from imex2d.domain import unit_conversions as uc


def _run(dx, permeability, injection_rate):
    scal = default_scal()
    model = one_dimensional_model(nx=20, dx=dx, permeability=permeability,
                                  injection_rate=injection_rate, scal=scal)
    return make_service(scal).run(model, short_config(end_time=60.0, snapshots=4))


def test_field_unit_inputs_converted_to_engine_units_reproduce_metric_result():
    # Case A — mühərrikin öz vahidləri (m, mD, m3/gün) birbaşa
    dx_m, k_md, rate_m3day = 8.0, 200.0, 60.0
    result_a = _run(dx_m, k_md, rate_m3day)

    # Case B — EYNİ fiziki ssenari FIELD vahidlərində ifadə olunub
    # (ft, mD, stb/gün), sonra YENİ çevirmə qatı ilə mühərrik vahidinə
    # geri çevrilib (istifadəçi FIELD daxil etsəydi baş verəcək yol)
    dx_ft = uc.m_to_ft(dx_m)
    rate_stbday = uc.m3_per_day_to_stb_per_day(rate_m3day)

    dx_m_roundtrip = uc.to_engine_units(dx_ft, "ft", "length")
    rate_m3day_roundtrip = uc.to_engine_units(rate_stbday, "stb/day", "rate")

    result_b = _run(dx_m_roundtrip, k_md, rate_m3day_roundtrip)

    assert abs(result_a.ooip - result_b.ooip) / result_a.ooip < 1e-9
    assert abs(result_a.final_recovery_factor - result_b.final_recovery_factor) < 1e-6
    assert result_a.steps == result_b.steps
    # Qeyd: `dx`-in son-bit fərqi (dəyirmi-səyahətdən) adaptiv addım
    # ölçüsü seçimini CÜzi (~1e-12) dəyişdirir, bu da 43 addım boyunca
    # yığılaraq ~1e-6 nisbi fərqə çatır — ƏSAS fiziki nəticələr (OOIP,
    # RF, addım sayı, yuxarıda) İSƏ dəqiq eynidir. Ona görə seriya
    # müqayisəsi realist (sərt deyil, amma mənalı) tolerantlıqla aparılır.
    assert np.allclose(result_a.series.oil_rate, result_b.series.oil_rate, rtol=1e-4)
    assert np.allclose(result_a.series.average_pressure, result_b.series.average_pressure,
                       rtol=1e-4)


def test_pressure_unit_round_trip_preserves_darcy_constant_inputs():
    """Təzyiq fərqi (psi <-> bar) Darsi axınına bir dəyişən kimi girir —
    dəyirmi-səyahətdən sonra ədədi giriş DƏYİŞMƏMƏLİDİR."""
    dp_bar = 35.0
    dp_psi = uc.bar_to_psi(dp_bar)
    dp_bar_roundtrip = uc.psi_to_bar(dp_psi)
    assert dp_bar_roundtrip == dp_bar or abs(dp_bar_roundtrip - dp_bar) < 1e-9
