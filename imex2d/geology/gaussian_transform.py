"""Normal-score (Gauss) çevirməsi — SGS-in tələb etdiyi ikitərəfli çevirmə.

SGS (Sequential Gaussian Simulation) STANDART normal N(0,1) fəzasında
kriging aparır (çünki şərti Gauss paylanması yalnız bu fəzada analitik
asandır: `N(mean_kriging, variance_kriging)`). Real geoloji xassələr
(PORO, PERMX...) isə demək olar HEÇ VAXT Gauss deyil (çarpıq, hədlərlə
məhdud). Bu modul ORİJİNAL paylanma ↔ Gauss fəzası arasında SIRA-ƏSASLI
(rank-based) çevirmə edir — Deutsch & Journel-in "normal score transform"-
unun standart metodu:

    1. Dəyərləri sırala, HƏR dəyərin empirik kvantilini tap (Hazen
       plotting position: `(rank - 0.5) / n`, TIES üçün ORTA RANK —
       `scipy.stats.rankdata(method="average")`).
    2. Bu kvantilə uyğun standart normal kvantili tap
       (`scipy.stats.norm.ppf`).
    3. Nəticə: (orijinal dəyər, Gauss dəyəri) CÜTLƏRİ CƏDVƏLİ — irəli/
       tərs çevirmə bu CƏDVƏLİN XƏTTİ İNTERPOLYASİYASIDIR.

Bu, TAMAMILƏ DETERMİNİSTİKDİR (heç bir təsadüfilik yoxdur) — eyni giriş
HƏMİŞƏ eyni cədvəli verir, ona görə TƏKRARLANA BİLƏNDİR.

**Hədlərdən kənar sorğu (ekstrapolyasiya)**: cədvəlin diapazonundan kənar
dəyər SƏRHƏD dəyərinə KƏSİLİR (clamp) — YENİ "quyruq" UYDURULMUR. Bu,
`OrdinaryKriging`/`BlackOilPVTProvider`-in də işlətdiyi eyni prinsipdir
(bax UNITS.md/GEOSTATISTICS.md) — sənədləşdirilmiş, AÇIQ seçimdir.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import norm, rankdata

#: Kvantil [0,1] sərhədində DƏQİQ 0/1 olanda `norm.ppf` ±sonsuzluq
#: qaytarır — bunun qarşısını almaq üçün kiçik kəsilmə.
_QUANTILE_EPS = 1e-6


@dataclass
class NormalScoreTransform:
    """Fit edilmiş (dəyər ↔ Gauss) axtarış cədvəli."""
    sorted_values: np.ndarray
    sorted_gaussian: np.ndarray
    is_constant: bool

    @classmethod
    def fit(cls, values) -> "NormalScoreTransform":
        values = np.asarray(values, float).ravel()
        if values.size == 0:
            raise ValueError("Normal-score çevirməsi boş massiv üçün qurula bilməz.")
        if np.any(~np.isfinite(values)):
            raise ValueError("Normal-score çevirməsi NaN/sonsuz dəyər qəbul etmir.")

        order = np.argsort(values, kind="stable")
        sorted_values = values[order]

        if np.ptp(sorted_values) < 1e-12:
            # SABİT XASSƏ: Gauss çevirməsi RİYAZİ CƏHƏTDƏN DEGENERATİVDİR
            # (bütün rütbələr bərabərdir) — irəli/tərs çevirmə sabit dəyəri
            # SAXLAYIR, UYDURULMUŞ dəyişkənlik ƏLAVƏ EDİLMİR.
            return cls(sorted_values=sorted_values, sorted_gaussian=np.zeros_like(sorted_values),
                      is_constant=True)

        n = values.size
        ranks = rankdata(values, method="average")
        quantile = np.clip((ranks - 0.5) / n, _QUANTILE_EPS, 1.0 - _QUANTILE_EPS)
        gaussian_all = norm.ppf(quantile)
        sorted_gaussian = gaussian_all[order]
        return cls(sorted_values=sorted_values, sorted_gaussian=sorted_gaussian, is_constant=False)

    def forward(self, values) -> np.ndarray:
        """Orijinal fəza → Gauss fəzası. Diapazondan kənar → sərhədə kəsilir."""
        values = np.atleast_1d(np.asarray(values, float))
        if self.is_constant:
            return np.zeros_like(values)
        return np.interp(values, self.sorted_values, self.sorted_gaussian,
                         left=self.sorted_gaussian[0], right=self.sorted_gaussian[-1])

    def inverse(self, gaussian_values) -> np.ndarray:
        """Gauss fəzası → orijinal fəza. Diapazondan kənar → sərhədə kəsilir
        (YENİ "quyruq" UYDURULMUR — bax modul docstring-i)."""
        gaussian_values = np.atleast_1d(np.asarray(gaussian_values, float))
        if self.is_constant:
            return np.full_like(gaussian_values, self.sorted_values[0])
        return np.interp(gaussian_values, self.sorted_gaussian, self.sorted_values,
                         left=self.sorted_values[0], right=self.sorted_values[-1])
