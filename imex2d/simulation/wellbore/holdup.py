"""Maye tutumu (liquid holdup) — sürüşmə (slip) modelinin nöqtəsi.

NƏ ÜÇÜNDÜR. İki fazalı şaquli axında qaz mayedən SÜRÜŞÜR — daha sürətli
qalxır. Ona görə borunun en kəsiyində mayenin tutduğu HƏQİQİ pay `H_L`
onun həcm payından `λ_L` BÖYÜKDÜR. Qravitasiya həddi məhz `H_L` ilə
çəkilir:

    ρ_slip = ρ_maye·H_L + ρ_qaz·(1 − H_L)

V1-də SÜRÜŞMƏ NƏZƏRƏ ALINMIR (`H_L = λ_L`) və bu, sənədləşir. Nəticə:
qazlı quyuda hidrostatik sütun OLDUĞUNDAN YÜNGÜL çıxır, yəni hesablanan
THP həqiqi dəyərdən YÜKSƏK olur. Fərq GOR artdıqca böyüyür.

NİYƏ BELƏ QƏRAR VERİLDİ. Hagedorn-Brown korrelyasiyası `CNL` və `ψ`
QRAFİKLƏRİNİN rəqəmsallaşdırılmasını tələb edir — o cədvəl datası bizdə
yoxdur və uydurulmayacaq. Beggs-Brill qapalı düsturlarla yazılır və
sonradan bu interfeysə TAM OTURUR: çağıran tərəf dəyişmir.
"""

from __future__ import annotations
from abc import ABC, abstractmethod


class IHoldupCorrelation(ABC):
    """Sürüşmə modelinin müqaviləsi.

    Bir metod, çünki traverse-in ehtiyacı budur: en kəsikdəki maye payı.
    Sürtünmə həddi `ρ_no-slip` ilə hesablanır və o, holdup-dan ASILI
    DEYİL — ona görə bu interfeysə girmir.

    Beggs-Brill əlavə olunanda `superficial_liquid_velocity`,
    `superficial_gas_velocity`, `diameter` və maili bucaq arqumentləri
    onun axın rejimi xəritəsi üçün kifayətdir — imza GENİŞLƏNMİR.
    """

    #: İnsan oxuya bilən ad — jurnalda və UI-də göstərilir.
    name: str = "holdup"

    @abstractmethod
    def liquid_holdup(self, no_slip_holdup: float,
                      superficial_liquid_velocity: float,
                      superficial_gas_velocity: float,
                      diameter: float,
                      inclination_rad: float = 0.0) -> float:
        """En kəsikdə mayenin HƏQİQİ payı `H_L` ∈ [0, 1].

        `no_slip_holdup` — λ_L, yəni sürüşməsiz həcm payı.
        `superficial_*_velocity` — faza öz-özünə bütün en kəsiyi
        tutsaydı alacağı sürət, m/s.
        `inclination_rad` — şaquliyə görə bucaq; v1-də həmişə 0
        (tam şaquli lülə), Beggs-Brill üçün lazım olacaq.
        """


class NoSlipHoldup(IHoldupCorrelation):
    """Sürüşməsiz: `H_L = λ_L`.

    Ən sadə və ən şəffaf model. Bir fazalı axında (λ_L = 0 və ya 1)
    DƏQİQDİR — sürüşmə yalnız iki faza birlikdə axanda mövcuddur.
    Ona görə qazsız neft-su quyusunda bu seçim heç bir xəta vermir.
    """

    name = "sürüşməsiz (no-slip)"

    def liquid_holdup(self, no_slip_holdup: float,
                      superficial_liquid_velocity: float,
                      superficial_gas_velocity: float,
                      diameter: float,
                      inclination_rad: float = 0.0) -> float:
        return min(max(float(no_slip_holdup), 0.0), 1.0)
