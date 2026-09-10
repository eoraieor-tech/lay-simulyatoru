"""BlackOilPVTProvider — IPVTProvider-in ilk implementasiyası.

Yeganə hesablama üsulu: cədvəl üzrə xətti interpolyasiya (np.interp).
Cədvəldən kənarda sərhəd dəyəri saxlanılır (np.interp-in defolt davranışı) —
bu, ekstrapolyasiyadan daha təhlükəsizdir.

Sıxılma cədvəldən ədədi törəmə ilə alınır və bir dəfə əvvəlcədən
hesablanır ki, hər zaman addımında yenidən hesablanmasın.
"""

from __future__ import annotations
from typing import Optional

import numpy as np

from ...domain.pvt import PVTTable
from ...interfaces.providers import IPVTProvider
from ...logging_setup import get_logger

LOG = get_logger(__name__)


class BlackOilPVTProvider(IPVTProvider):

    def __init__(self, table: PVTTable,
                 dead_oil_below_bubble_point: bool = False):
        """`dead_oil_below_bubble_point` — bax `_dead_oil_fvf`.

        DEFOLT `False`: bu sinif cədvəlin TƏMİZ interpolyatorudur və
        düyünlərdə cədvəl dəyərlərini olduğu kimi qaytarır
        (`test_pvt.py::test_provider_reproduces_table_values_at_nodes`).

        `True` veriləndə (bunu `application` qatı seçir — bax
        `simulation_service._build_pvt_provider`) və cədvəldə QAZ
        SÜTUNLARI YOXDURSA, doyma təzyiqindən aşağı Bo ölü-neft budağı
        ilə əvəz olunur. Qaz sütunları varsa bayraq TƏSİRSİZDİR.
        """
        issues = table.validate()
        if issues:
            raise ValueError("PVT cədvəli yararsızdır: " + "; ".join(issues))
        self.table = table
        #: Doyma təzyiqindən aşağı Bo düzəlişi — bax `_dead_oil_fvf`.
        #: Qaz modelləşdirilməyəndə cədvəlin doymuş budağı Nyutonu
        #: pozur; bu massiv düzəldilmiş (və ya toxunulmamış) Bo-dur.
        self._oil_fvf, self.dead_oil_corrected = (
            self._dead_oil_fvf(table) if dead_oil_below_bubble_point
            else (table.oil_fvf, False))
        self._co = table.compressibility("oil_fvf")
        self._cw = table.compressibility("water_fvf")
        self._cr = table.rock_compressibility
        # Törəmələr bir dəfə hesablanır. Cədvəl PARÇALI XƏTTİ olduğuna
        # görə törəmə hər intervalda sabitdir — `np.gradient`-in verdiyi
        # hamar qiymət deyil, məhz interval meyli. Bu seçim ölçmə ilə
        # təsdiqləndi: hamar törəmə Nyutonun iterasiya sayını
        # dəyişmir, lakin Jakobianı sonlu fərqdən 10⁶ dəfə uzaqlaşdırır.
        pressure = table.pressure
        columns = [("oil_fvf", self._oil_fvf), ("water_fvf", table.water_fvf),
                  ("oil_viscosity", table.oil_viscosity),
                  ("water_viscosity", table.water_viscosity),
                  ("solution_gor", table.solution_gor)]
        if table.has_gas_phase:
            columns += [("gas_fvf", table.gas_fvf),
                       ("gas_viscosity", table.gas_viscosity)]
        self._slopes = {name: np.diff(values) / np.diff(pressure)
                        for name, values in columns}

    # ═══════════════════════ doyma təzyiqindən aşağı Bo düzəlişi (B3-A)
    @staticmethod
    def _dead_oil_fvf(table: PVTTable):
        """Qaz MODELLƏŞDİRİLMƏYƏNDƏ Bo-nun doymuş budağını əvəz edir.

        PROBLEM (ölçülüb, `ISH_HESABATI.md` → Seans 7). Korrelyasiya
        cədvəlində Bo doyma təzyiqinə qədər ARTIR (qaz həll olur, neft
        şişir), ondan sonra AZALIR (sıxılma). Yəni `dBo/dp` işarə
        dəyişir və yolda SIFIRDAN KEÇİR.

        Neft tənliyinin təzyiq üzrə diaqonal törəməsi
        `−So·B'o/Bo²`-yə mütənasibdir. İki fazalı (neft-su) modeldə bu
        o deməkdir ki, həmin keçid nöqtəsində diaqonal SIFIRA düşür və
        Jakobian təkləşir: kiçik qalıq nəhəng Nyuton addımına çevrilir
        (ölçüldü: 656 bar, halbuki bütün lay 176–299 bar aralığında
        idi), qoruyucular onu kəsir və iterasiya DONUR.

        NİYƏ BU, ƏDƏDİ QÜSUR DEYİL — MODEL NATAMAMLIĞIDIR. Doyma
        təzyiqindən aşağı neftdən qaz ayrılır; doymuş budaq məhz o
        ayrılmanı təsvir edir. İki fazalı model ayrılan qazı NƏZƏRƏ
        ALMIR, ona görə tənliklərdə "mənfi neft sıxılması" kimi görünür.
        Üç fazalı mühərrikdə isə qaz tənliyi bunu kompensasiya edir və
        doymuş budaq DÜZGÜNDÜR — ona görə düzəliş YALNIZ qaz sütunları
        olmayan cədvələ tətbiq olunur.

        DÜZƏLİŞ: doyma təzyiqindən aşağı undersaturated (sıxılma)
        meyli uzadılır — yəni "neftdən qaz ayrılmır, tərkibi sabit
        qalır" (ölü neft) fərziyyəsi. Belədə `dBo/dp < 0` HƏR YERDƏ
        qalır, sıxılma müsbətdir, Jakobian təkləşmir.

        Ölçülmüş nəticə (11×11, ilkin 250 bar, istismarçı BHP 150 bar):

            Pb     xam cədvəl        düzəlişlə
            150    ✅ 28 addım       ✅ 28 addım  (DƏYİŞMİR)
            200    ✅ 31 addım       ✅ 31 addım
            220    ❌ yığılmır       ✅ 32 addım, RF 62.66 %
            240    ❌ yığılmır       ✅ 34 addım, RF 63.21 %
            300    ❌ yığılmır       ✅ 36 addım, RF 64.36 %

        Heç bir hüceyrə Pb-dən aşağı düşmürsə nəticə BİTƏ-BİT eynidir
        (Pb=150 sətri) — köhnə modellər təsirlənmir.

        ⏳ AÇIQ SUAL: özlülük də doymuş budaqda qalır (Pb-dən aşağı
        ARTIR, sanki qaz ayrılıb). Onu da ölü-neft budağına keçirmək
        olar; ölçüldü — yığılmaya TƏSİRİ YOXDUR, yalnız RF dəyişir
        (Pb=300: 62.92 % → 64.36 %). Minimal müdaxilə prinsipi ilə
        TOXUNULMADI; sahibkar qərar verməlidir.

        Qaytarır: (bo_massivi, düzəliş_edildimi).
        """
        if table.has_gas_phase:
            return table.oil_fvf, False        # üç fazalı — cədvəl DÜZGÜNDÜR

        bo = np.asarray(table.oil_fvf, float)
        pressure = np.asarray(table.pressure, float)
        if bo.size < 3:
            return table.oil_fvf, False

        peak = int(np.argmax(bo))
        # Zirvə ya kənardadırsa, ya da artan budaq yoxdursa — düzəliş lazım deyil
        if peak == 0 or peak + 1 >= bo.size:
            return table.oil_fvf, False

        slope = (bo[peak + 1] - bo[peak]) / (pressure[peak + 1] - pressure[peak])
        if slope >= 0.0:
            return table.oil_fvf, False        # undersaturated budaq yoxdur

        corrected = bo.copy()
        below = np.arange(peak)
        corrected[below] = bo[peak] + slope * (pressure[below] - pressure[peak])
        LOG.info("PVT: qaz modelləşdirilmir — doyma təzyiqindən (%.0f bar) "
                 "aşağı Bo ölü-neft budağı ilə əvəz edildi (%d düyün). "
                 "Səbəb: doymuş budaq iki fazalı Jakobianı təkləşdirir.",
                 float(table.bubble_point), int(peak))
        return corrected, True

    # ------------------------------------------------------ interpolyasiya
    def _interp(self, values: np.ndarray, pressure) -> np.ndarray:
        return np.interp(np.asarray(pressure, float), self.table.pressure, values)

    def oil_fvf(self, pressure, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._interp(self._oil_fvf, pressure)

    def oil_viscosity(self, pressure, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._interp(self.table.oil_viscosity, pressure)

    def water_fvf(self, pressure, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._interp(self.table.water_fvf, pressure)

    def water_viscosity(self, pressure, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._interp(self.table.water_viscosity, pressure)

    def solution_gor(self, pressure, region: Optional[np.ndarray] = None) -> np.ndarray:
        return self._interp(self.table.solution_gor, pressure)

    def total_compressibility(self, pressure, sw, region: Optional[np.ndarray] = None) -> np.ndarray:
        sw = np.asarray(sw, float)
        co = self._interp(self._co, pressure)
        cw = self._interp(self._cw, pressure)
        return cw * sw + co * (1.0 - sw) + self._cr

    def bubble_point(self, region: Optional[np.ndarray] = None) -> float:
        return self.table.bubble_point

    def has_gas_phase(self, region: Optional[np.ndarray] = None) -> bool:
        return self.table.has_gas_phase

    def gas_fvf(self, pressure, region: Optional[np.ndarray] = None) -> np.ndarray:
        if not self.table.has_gas_phase:
            raise NotImplementedError(
                "Bu PVT cədvəlində qaz xassələri yoxdur "
                "(build_pvt_table(..., include_gas=True) işlədin).")
        return self._interp(self.table.gas_fvf, pressure)

    def gas_viscosity(self, pressure, region: Optional[np.ndarray] = None) -> np.ndarray:
        if not self.table.has_gas_phase:
            raise NotImplementedError(
                "Bu PVT cədvəlində qaz xassələri yoxdur "
                "(build_pvt_table(..., include_gas=True) işlədin).")
        return self._interp(self.table.gas_viscosity, pressure)

    # ─────────────────────────────────────────── analitik törəmələr
    def _slope(self, name: str, pressure) -> np.ndarray:
        """Parçalı xətti cədvəlin dəqiq törəməsi.

        Cədvəldən kənarda sıfırdır, çünki `np.interp` orada sərhəd
        dəyərini saxlayır (funksiya sabitdir).
        """
        pressure = np.atleast_1d(np.asarray(pressure, float))
        nodes = self.table.pressure
        index = np.clip(np.searchsorted(nodes, pressure, side="right") - 1,
                        0, nodes.size - 2)
        slopes = self._slopes[name][index]
        outside = (pressure < nodes[0]) | (pressure > nodes[-1])
        return np.where(outside, 0.0, slopes)

    def oil_fvf_derivative(self, pressure, region=None) -> np.ndarray:
        return self._slope("oil_fvf", pressure)

    def water_fvf_derivative(self, pressure, region=None) -> np.ndarray:
        return self._slope("water_fvf", pressure)

    def oil_viscosity_derivative(self, pressure, region=None) -> np.ndarray:
        return self._slope("oil_viscosity", pressure)

    def water_viscosity_derivative(self, pressure, region=None) -> np.ndarray:
        return self._slope("water_viscosity", pressure)

    def gas_fvf_derivative(self, pressure, region=None) -> np.ndarray:
        return self._slope("gas_fvf", pressure)

    def gas_viscosity_derivative(self, pressure, region=None) -> np.ndarray:
        return self._slope("gas_viscosity", pressure)

    def solution_gor_derivative(self, pressure, region=None) -> np.ndarray:
        """dRs_sat/dp — Jakobianda doymuş hüceyrələr üçün lazımdır (A7/6c).

        Doyma təzyiqindən yuxarıda Rs sabitdir, ona görə bu, sıfıra
        düşür — cədvəlin özündə bu sabitlik artıq mövcuddur, əlavə
        şərtə ehtiyac yoxdur.
        """
        return self._slope("solution_gor", pressure)
