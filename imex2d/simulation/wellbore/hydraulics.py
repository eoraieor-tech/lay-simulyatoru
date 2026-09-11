"""`IWellHydraulicsProvider` implementasiyası + nəticəyə THP yazılması.

V1 (A variantı) QƏSDƏN POST-PROSESDİR. Mühərrik toxunulmur: BHP rejimli
quyuda quyu dibi təzyiqi onsuz da SABİTDİR (istifadəçinin verdiyi hədəf),
debit sıraları isə `SimulationResult`-da artıq var. Ona görə THP
simulyasiya bitəndən sonra hesablanır.

NƏ QAZANIRIQ: Nyutonun yığılmasına SIFIR risk. `ControlMode.THP`
(B variantı) gələndə əlaqə addımın əvvəlinə köçəcək — o zaman bu sinif
YENİDƏN İŞLƏDİLİR, təkrar yazılmır.

MƏHDUDİYYƏTLƏR (v1, sənədləşir — uydurulmur):
  * yalnız **BHP rejimli** quyular — RATE rejimində mühərrik BHP-ni
    ümumiyyətlə hesablamır (Peaceman-ın tərsi ⏳ sonraya);
  * yalnız **istismarçılar** — vurucuda axın aşağıdır, traversin işarə
    konvensiyası fərqlidir ⏳;
  * lülə **tam şaquli** (TVD = MD);
  * **sürüşmə yoxdur** (no-slip) — bax `holdup.py`.
"""

from __future__ import annotations
from typing import Dict, List, Optional

import numpy as np

from ...domain.tubing import TubingGeometry
from ...domain.wells import ControlMode, WellType
from ...interfaces.providers import IWellHydraulicsProvider
from ...logging_setup import get_logger
from .holdup import IHoldupCorrelation, NoSlipHoldup
from .traverse import WellStream, pressure_traverse

LOG = get_logger(__name__)


class WellboreHydraulics(IWellHydraulicsProvider):
    """Analitik traversə əsaslanan hidravlika provider-i."""

    def __init__(self, pvt=None, fluids=None,
                 holdup: Optional[IHoldupCorrelation] = None):
        self.pvt = pvt
        self.fluids = fluids
        self.holdup = holdup or NoSlipHoldup()

    # ─────────────────────────────────────────────── müqavilə metodu
    def tubing_head_pressure(self, bhp: float, perforation_depth: float,
                             stream: WellStream,
                             tubing: TubingGeometry) -> float:
        result = pressure_traverse(bhp, perforation_depth, tubing, stream,
                                   pvt=self.pvt, fluids=self.fluids,
                                   holdup=self.holdup)
        return result.thp

    # ────────────────────────────────────── nəticənin zənginləşdirilməsi
    def annotate(self, model, result) -> None:
        """`result`-a `well_bhp` və `well_thp` zaman sıralarını yazır.

        Şərtləri ödəməyən quyu SƏSSİZCƏ atılmır — jurnala səbəbi ilə
        yazılır, çünki istifadəçi UI-də boş qrafik görəndə niyəsini
        bilməlidir.
        """
        depths = self._perforation_depths(model)
        for well in model.wells:
            reason = self._skip_reason(well, depths)
            if reason is not None:
                LOG.info("THP hesablanmadı — quyu '%s': %s", well.name, reason)
                continue

            tubing = well.tubing
            depth = depths[well.name]
            bhp = float(well.control.target)

            oil = result.well_oil_rate.get(well.name, [])
            water = result.well_water_rate.get(well.name, [])
            gas = result.well_gas_rate.get(well.name, [])
            steps = len(oil)
            if steps == 0:
                LOG.info("THP hesablanmadı — quyu '%s': debit sırası boşdur "
                         "(`record_well_rates` söndürülüb?)", well.name)
                continue

            thp_series: List[float] = []
            for i in range(steps):
                stream = WellStream(
                    oil=float(oil[i]),
                    water=float(water[i]) if i < len(water) else 0.0,
                    gas=float(gas[i]) if i < len(gas) else 0.0,
                    oil_density=model.fluids.oil_density,
                    water_density=model.fluids.water_density,
                    gas_density=model.fluids.gas_density)
                thp_series.append(self.tubing_head_pressure(
                    bhp, depth, stream, tubing))

            result.well_thp[well.name] = thp_series
            result.well_bhp[well.name] = [bhp] * steps

        self._log_summary(result)

    # ─────────────────────────────────────────────────────── köməkçilər
    @staticmethod
    def _skip_reason(well, depths: Dict[str, float]) -> Optional[str]:
        if not well.active:
            return "quyu söndürülüb"
        if well.tubing is None:
            return "lülə həndəsəsi verilməyib"
        if well.well_type is not WellType.PRODUCER:
            return "v1-də yalnız istismarçılar dəstəklənir"
        if well.control.mode is not ControlMode.BHP:
            return ("v1-də yalnız BHP rejimi dəstəklənir — RATE rejimində "
                    "mühərrik quyu dibi təzyiqini hesablamır")
        if well.name not in depths:
            return "aktiv perforasiya yoxdur"
        return None

    @staticmethod
    def _perforation_depths(model) -> Dict[str, float]:
        """Hər quyu üçün BHP istinad dərinliyi, m.

        Çoxperforasiyalı quyuda WI ilə ÇƏKİLMİŞ orta dərinlik götürülür —
        quyu dibi təzyiqinin fiziki istinad nöqtəsi budur. Mühərrikin
        özü BHP-ni bütün perforasiyalara eyni tətbiq edir (dərinlik
        düzəlişi YOXDUR), ona görə bu, mövcud sadələşdirmə ilə uyğundur.
        """
        from ..well_model import PeacemanWellModel

        try:
            connections = PeacemanWellModel().build_connections(model)
        except Exception as error:            # model natamam ola bilər
            LOG.info("Perforasiya dərinlikləri alınmadı: %s", error)
            return {}

        cell_depths = model.geometry.cell_depths()
        weighted: Dict[str, float] = {}
        weights: Dict[str, float] = {}
        for connection in connections:
            w = max(float(connection.well_index), 0.0)
            depth = float(cell_depths[connection.cell])
            weighted[connection.well_name] = weighted.get(connection.well_name, 0.0) + w * depth
            weights[connection.well_name] = weights.get(connection.well_name, 0.0) + w

        out: Dict[str, float] = {}
        for name, total in weights.items():
            if total > 0.0:
                out[name] = weighted[name] / total
        return out

    @staticmethod
    def _log_summary(result) -> None:
        if not result.well_thp:
            return
        parts = []
        for name, series in result.well_thp.items():
            finite = [v for v in series if np.isfinite(v)]
            if finite:
                parts.append(f"{name}: {finite[-1]:.1f} bar")
        if parts:
            LOG.info("THP hesablandı (son addım) — %s", ", ".join(parts))
