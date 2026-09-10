"""Üç fazalı nisbi keçiricilik — Stone II modeli (A7, mərhələ 3).

İki fazalı sistemdə neftin keçiriciliyi tək dəyişəndən (Sw) asılıdır.
Üç fazada neft həm suyun, həm qazın "sıxışdırmasına" məruz qalır və
tək əyri kifayət etmir — Stone (1973) İKİ İKİFAZALI əyrini (su-neft
və qaz-neft) BİRLƏŞDİRİR:

    kro = kro_end · [(krow/kro_end + krw) · (krog/kro_end + krg)
                     − (krw + krg)]

    krow(Sw) — su-neft sistemində neft (mövcud iki fazalı əyri)
    krog(Sg) — qaz-neft sistemində neft (GasCoreyParameters)
    krw(Sw)  — mövcud iki fazalı əyri
    krg(Sg)  — GasCoreyParameters

Stone I (1970) orijinal versiyada mənfi `kro` verə bilirdi; Stone II
(1973) bunu riyazi olaraq aradan qaldırır və sənayedə standart
seçimdir (Eclipse-in defoltu).

TUTARLILIQ YOXLAMASI (ən vacib fiziki xassə)
    Sg = 0-da:   kro(Sw, 0) = krow(Sw)     — iki fazalı su-neft
    Sw = Swc-də: kro(Swc, Sg) = krog(Sg)   — iki fazalı qaz-neft

Yəni Stone modeli iki fazalı hallara DƏQIQ REDUKSİYA olunur — bu,
düzgünlüyün ən güclü sınağıdır və testlə yoxlanılıb.
"""

from __future__ import annotations
from typing import Optional

import numpy as np

from ..domain.scal import CoreyParameters, GasCoreyParameters
from ..interfaces.providers import IRelativePermeabilityProvider


class StoneRelativePermeabilityProvider(IRelativePermeabilityProvider):
    """İki fazalı provider-i (Corey və ya cədvəl) qaz əyrisi ilə birləşdirir.

    Mövcud su-neft provider-ə TOXUNMUR — onu sarğı (wrap) edir. Bütün
    iki fazalı metodlar (`krw`, `kro`, `saturation_limits`, CFL) birbaşa
    ona ötürülür, ona görə mövcud IMPES/implicit mühərrik dəyişmədən
    işləməyə davam edir (qaz hələ residual/Jacobian-a bağlanmayıb).
    """

    def __init__(self, water_oil: IRelativePermeabilityProvider,
                 gas: GasCoreyParameters, swc: float, kro_end: float):
        self.water_oil = water_oil
        self.gas = gas
        self.swc = float(swc)
        self.kro_end = float(kro_end)
        issues = gas.validate(swc)
        if issues:
            raise ValueError("Qaz SCAL parametrləri yararsızdır: "
                             + "; ".join(issues))

    @classmethod
    def from_corey(cls, water_oil: CoreyParameters,
                   gas: GasCoreyParameters) -> "StoneRelativePermeabilityProvider":
        """Ən çox işlədilən yol: iki Corey dəstindən Stone provider qurur.

        `kro_end` avtomatik su-neft parametrlərindən götürülür — əl ilə
        təkrar yazılmır, iki əyrinin uyğunsuz olma riski aradan qalxır.
        """
        from ..simulation.scal_adapter import CoreyRelativePermeabilityAdapter
        return cls(CoreyRelativePermeabilityAdapter(water_oil), gas,
                  water_oil.swc, water_oil.kro_end)

    # ─────────────────────────────────────────── iki fazalı (ötürmə)
    def krw(self, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self.water_oil.krw(sw, region)

    def krw_derivative(self, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        """dkrw/dSw — sarılmış su-neft provider-ə ötürülür (dəyişməyib,
        Stone II yalnız kro-nu dəyişir, krw-yə toxunmur)."""
        analytic = getattr(self.water_oil, "krw_derivative", None)
        if analytic is not None:
            return np.asarray(analytic(sw, region), float)
        step = 1e-6
        sw = np.asarray(sw, float)
        return (np.asarray(self.water_oil.krw(sw + step, region), float)
               - np.asarray(self.water_oil.krw(sw - step, region), float)) / (2 * step)

    def kro(self, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        """İki fazalı `kro(Sw)` — Sg=0 fərziyyəsi ilə (geriyə uyğunluq)."""
        return self.water_oil.kro(sw, region)

    def saturation_limits(self, region: Optional[int] = None) -> tuple:
        return self.water_oil.saturation_limits(region)

    def endpoint_water_mobility(self, water_viscosity: float,
                                region: Optional[int] = None) -> float:
        return self.water_oil.endpoint_water_mobility(water_viscosity, region)

    def max_fractional_flow_derivative(self, water_viscosity: float,
                                       oil_viscosity: float,
                                       region: Optional[int] = None) -> float:
        return self.water_oil.max_fractional_flow_derivative(
            water_viscosity, oil_viscosity, region)

    # ─────────────────────────────────────────────────────── üç fazalı
    def has_gas_phase(self, region: Optional[int] = None) -> bool:
        return True

    def krg(self, sg, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self.gas.krg(sg, self.swc)

    def krog(self, sg, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self.gas.krog(sg, self.swc, self.kro_end)

    def kro_three_phase(self, sw, sg,
                        region: Optional[np.ndarray] = None) -> np.ndarray:
        """Stone II — iki iki-fazalı əyridən üç fazalı kro qurur."""
        sw = np.asarray(sw, float)
        sg = np.asarray(sg, float)
        krow = np.asarray(self.water_oil.kro(sw, region), float)
        krow_g = self.krog(sg, region)
        krw = np.asarray(self.water_oil.krw(sw, region), float)
        krg = self.krg(sg, region)

        end = max(self.kro_end, 1e-9)
        result = end * ((krow / end + krw) * (krow_g / end + krg)
                        - (krw + krg))
        return np.clip(result, 0.0, self.kro_end)

    def kro_three_phase_derivatives(self, sw, sg,
                                    region: Optional[np.ndarray] = None):
        """(∂kro/∂Sw, ∂kro/∂Sg) — Stone II analitik zəncirvari qaydası.

        Stone II sadə hasil+cəm formasında olduğu üçün (krow, krw yalnız
        Sw-dən; krog, krg yalnız Sg-dən asılıdır) tam analitik törəmə
        mümkündür — sonlu fərqə ehtiyac yoxdur:

            kro = end·[(krow/end+krw)(krog/end+krg) − (krw+krg)]

            ∂kro/∂Sw = end·[(dkrow/dSw/end + dkrw/dSw)·(krog/end+krg)
                            − dkrw/dSw]
            ∂kro/∂Sg = end·[(krow/end+krw)·(dkrog/dSg/end + dkrg/dSg)
                            − dkrg/dSg]

        Kəsilmə (clip) tətbiq olunan nöqtələrdə (nəticə [0,kro_end]
        xaricinə çıxan yerlərdə) törəmə sıfırdır — orada kro sabitdir.
        """
        sw = np.asarray(sw, float)
        sg = np.asarray(sg, float)
        end = max(self.kro_end, 1e-9)

        krow = np.asarray(self.water_oil.kro(sw, region), float)
        krw = np.asarray(self.water_oil.krw(sw, region), float)
        krog = self.krog(sg, region)
        krg = self.krg(sg, region)

        dkrow_dsw = getattr(self.water_oil, "kro_derivative",
                            lambda s, r=None: np.zeros_like(sw))(sw, region)
        dkrw_dsw = getattr(self.water_oil, "krw_derivative",
                           lambda s, r=None: np.zeros_like(sw))(sw, region)
        dkrog_dsg = self.gas.krog_derivative(sg, self.swc, self.kro_end)
        dkrg_dsg = self.gas.krg_derivative(sg, self.swc)

        d_dsw = end * ((dkrow_dsw / end + dkrw_dsw) * (krog / end + krg)
                       - dkrw_dsw)
        d_dsg = end * ((krow / end + krw) * (dkrog_dsg / end + dkrg_dsg)
                       - dkrg_dsg)

        unclipped = end * ((krow / end + krw) * (krog / end + krg)
                           - (krw + krg))
        active = (unclipped > 0.0) & (unclipped < self.kro_end)
        return np.where(active, d_dsw, 0.0), np.where(active, d_dsg, 0.0)

    def gas_saturation_limits(self, region: Optional[int] = None) -> tuple:
        low, high = self.gas.sgc, 1.0 - self.swc - self.gas.sorg
        return (low, max(high, low))
