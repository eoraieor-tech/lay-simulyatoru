"""Cədvəl əsaslı SCAL provider — müvəqqəti Corey adapterinin əvəzi.

`CoreyRelativePermeabilityAdapter` refaktorinqdə körpü kimi yazılmışdı
(bax onun sənədləşməsi). Bu modul onun yerini tutur və iki şey əlavə edir:

    · laboratoriya cədvəlləri (Corey düsturu deyil)
    · REGION üzrə fərqli əyrilər (SATNUM)

Corey adapteri SİLİNMİR — sadə modellər üçün hələ də faydalıdır və
bütün mövcud testlər ondan asılıdır. İkisi eyni interfeysi paylaşır.

REGION İNDEKSLƏŞMƏSİ
Mühərrik `krw(sw)` çağıranda bütün hüceyrələr üçün bir massiv gözləyir.
Regionlar fərqli olduqda hər region üçün ayrıca interpolyasiya aparılır
və nəticə birləşdirilir. Tək region halında bu, sadə interpolyasiyaya
düşür — əlavə xərc yoxdur.
"""

from __future__ import annotations
from typing import Optional

import numpy as np

from ..domain.scal_tables import SaturationTable, SaturationTableSet
from ..interfaces.providers import (ICapillaryPressureProvider,
                                    IRelativePermeabilityProvider)


class TableRelativePermeabilityProvider(IRelativePermeabilityProvider):
    """Cədvəldən nisbi keçiricilik, region dəstəyi ilə."""

    def __init__(self, tables: SaturationTableSet,
                 region_ids: Optional[np.ndarray] = None):
        issues = tables.validate()
        if issues:
            raise ValueError("SCAL cədvəlləri yararsızdır: "
                             + "; ".join(issues))
        self.tables = tables
        self.region_ids = (None if region_ids is None
                           else np.asarray(region_ids, int).ravel())
        self._single = (self.region_ids is None
                        or np.unique(self.region_ids).size <= 1)

    # ─────────────────────────────────────────── region köməkçisi
    def _evaluate(self, sw, extractor, region: Optional[np.ndarray] = None):
        sw = np.asarray(sw, float)
        if self._single or region is not None:
            index = None
            if region is not None and np.ndim(region) == 0:
                index = int(region)
            elif self.region_ids is not None and self.region_ids.size:
                index = int(self.region_ids[0])
            return extractor(self.tables.get(index), sw)

        result = np.empty_like(sw)
        ids = self.region_ids
        if ids.size != sw.size:
            # ölçü uyğun gəlmirsə (məsələn tək nöqtə) defolt cədvəl
            return extractor(self.tables.get(), sw)
        for region_id in np.unique(ids):
            mask = ids == region_id
            result[mask] = extractor(self.tables.get(int(region_id)), sw[mask])
        return result

    # ─────────────────────────────────────────── interfeys
    def krw(self, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._evaluate(sw, lambda t, s: t.interpolate_krw(s), region)

    def kro(self, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._evaluate(sw, lambda t, s: t.interpolate_kro(s), region)

    def krw_derivative(self, sw, region: Optional[np.ndarray] = None):
        return self._evaluate(sw, lambda t, s: t.slope(t.krw, s), region)

    def kro_derivative(self, sw, region: Optional[np.ndarray] = None):
        return self._evaluate(sw, lambda t, s: t.slope(t.kro, s), region)

    def saturation_limits(self, region: Optional[int] = None) -> tuple:
        """Bütün regionların ƏN DAR ortaq intervalı.

        Fərqli regionlarda Swc/Sor fərqli olur. Mühərrik doyumluluğu
        bir hədd cütü ilə kəsir, ona görə ən məhdudlaşdırıcı interval
        götürülür — əks halda bəzi hüceyrələrdə kr cədvəldən kənara
        çıxardı.
        """
        if region is not None:
            table = self.tables.get(region)
            return (table.swc, 1.0 - table.sor)
        low = max(table.swc for table in self.tables.tables.values())
        high = min(1.0 - table.sor for table in self.tables.tables.values())
        return (low, high)

    def endpoint_water_mobility(self, sw=1.0,
                                region: Optional[int] = None) -> float:
        return float(self.tables.get(region).krw_end)

    def max_fractional_flow_derivative(self, water_viscosity: float,
                                       oil_viscosity: float,
                                       region: Optional[int] = None) -> float:
        """max |dfw/dSw| — IMPES-in CFL limiti bundan asılıdır.

        Bütün regionlar üzrə ƏN BÖYÜK dəyər götürülür: zaman addımı
        modelin ən sərt hüceyrəsi ilə məhdudlaşmalıdır, orta ilə yox.
        Əks halda gilli zonada həll qeyri-stabil olardı.
        """
        candidates = ([self.tables.get(region)] if region is not None
                      else list(self.tables.tables.values()))
        worst = 0.0
        for table in candidates:
            low, high = table.swc, 1.0 - table.sor
            if high <= low:
                continue
            sw = np.linspace(low, high, 400)
            mobility_w = table.interpolate_krw(sw) / water_viscosity
            mobility_o = table.interpolate_kro(sw) / oil_viscosity
            fractional = mobility_w / np.maximum(mobility_w + mobility_o, 1e-30)
            worst = max(worst,
                        float(np.nanmax(np.abs(np.gradient(fractional, sw)))))
        return worst


class TableCapillaryPressureProvider(ICapillaryPressureProvider):
    """SCAL cədvəlinin Pc sütunundan kapilyar təzyiq.

    Brooks-Corey analitik modelindən fərqli olaraq burada ölçülmüş
    əyri işlədilir — laboratoriya məlumatı olan hallarda daha doğrudur.
    """

    def __init__(self, tables: SaturationTableSet,
                 region_ids: Optional[np.ndarray] = None):
        self._provider = TableRelativePermeabilityProvider(tables, region_ids)
        self.tables = tables

    def pcow(self, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._provider._evaluate(
            sw, lambda t, s: t.interpolate_pc(s), region)

    def dpcow_dsw(self, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._provider._evaluate(
            sw, lambda t, s: (t.slope(t.pc, s) if t.pc is not None
                              else np.zeros_like(np.atleast_1d(s))), region)

    def has_capillary_pressure(self) -> bool:
        return any(table.has_capillary
                   for table in self.tables.tables.values())
