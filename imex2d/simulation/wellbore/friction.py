"""Darcy-Weisbach sürtünmə əmsalı — Chen (1979) açıq düsturu.

NİYƏ CHEN, COLEBROOK YOX. Colebrook-White tənliyi `f`-ə görə qapalı
formada həll olunmur — hər çağırışda iterasiya tələb edir. Traverse
hər quyu üçün hər zaman addımında 20 seqment hesablayır, ona görə
iterasiyasız yaxınlaşma seçildi. Chen (1979) Colebrook-dan turbulent
aralıqda **< 0.5 %** fərqlənir (orijinal məqalənin öz ölçməsi) —
bu, PVT və quyu modelinin öz qeyri-müəyyənliyindən xeyli kiçikdir.

Mənbə: Chen, N.H. (1979), "An Explicit Equation for Friction Factor in
Pipe", Ind. Eng. Chem. Fundam. 18(3), 296-297.
"""

from __future__ import annotations
import math

#: Laminar → turbulent keçid. Aşağısı analitik (`64/Re`), yuxarısı Chen.
LAMINAR_LIMIT = 2000.0
TURBULENT_LIMIT = 4000.0

#: Re sıfıra çox yaxın olanda (quyu bağlıdır, debit ≈ 0) `64/Re` partlayır.
#: Bu həddən aşağı axın sürtünməsiz sayılır — fiziki olaraq da doğrudur,
#: çünki dayanmış mayedə sürtünmə itkisi yoxdur.
MIN_REYNOLDS = 1e-9


def reynolds(density: float, velocity: float, diameter: float,
             viscosity_cp: float) -> float:
    """Reynolds ədədi — ÖLÇÜSÜZ.

    `viscosity_cp` sentipuazdadır (layihənin daxili vahidi), daxildə
    Pa·s-ə çevrilir: 1 cP = 1e-3 Pa·s.
    """
    mu = max(float(viscosity_cp), 1e-12) * 1e-3
    return abs(float(density) * float(velocity) * float(diameter)) / mu


def chen_friction_factor(reynolds_number: float,
                         relative_roughness: float) -> float:
    """Darcy sürtünmə əmsalı `f` (Fanning DEYİL — Fanning-in 4 misli).

    Darcy seçildi, çünki təzyiq düşgüsü düsturu `Δp = f·ρ·v²·L/(2d)`
    bu formada yazılır; Fanning işlətsək düsturda əlavə 4 əmsalı
    unudula bilərdi.

    Keçid aralığında (2000 < Re < 4000) laminar və turbulent qiymətlər
    arasında XƏTTİ keçid tətbiq olunur. Real axın orada qeyri-sabitdir və
    heç bir düstur onu dəqiq vermir; xətti keçid süni SIÇRAYIŞIN
    qarşısını alır — sıçrayış Nyutonun yığılmasını pozardı.
    """
    re = float(reynolds_number)
    if re <= MIN_REYNOLDS:
        return 0.0
    if re < LAMINAR_LIMIT:
        return 64.0 / re

    turbulent = _chen_turbulent(max(re, LAMINAR_LIMIT), relative_roughness)
    if re >= TURBULENT_LIMIT:
        return turbulent

    laminar = 64.0 / LAMINAR_LIMIT
    turbulent_at_limit = _chen_turbulent(TURBULENT_LIMIT, relative_roughness)
    weight = (re - LAMINAR_LIMIT) / (TURBULENT_LIMIT - LAMINAR_LIMIT)
    return laminar + weight * (turbulent_at_limit - laminar)


def _chen_turbulent(reynolds_number: float, relative_roughness: float) -> float:
    """Chen (1979) — yalnız turbulent budaq.

        1/√f = −2·log₁₀[ ε/(3.7065·d)
                         − (5.0452/Re)·log₁₀( (ε/d)^1.1098/2.8257
                                              + (7.149/Re)^0.8981 ) ]
    """
    eps = max(float(relative_roughness), 0.0)
    re = float(reynolds_number)

    inner = (eps ** 1.1098) / 2.8257 + (7.149 / re) ** 0.8981
    bracket = eps / 3.7065 - (5.0452 / re) * math.log10(inner)

    # Hamar boruda (ε=0) və çox böyük Re-də mötərizə sıfıra yaxınlaşa
    # bilər; log10(≤0) riyazi xətadır. Fiziki olaraq bu, sürtünməsiz
    # hədddir — `f`-i Blasius hamar-boru qiyməti ilə əvəz edirik.
    if bracket <= 0.0:
        return 0.3164 * re ** -0.25

    return (-2.0 * math.log10(bracket)) ** -2.0
