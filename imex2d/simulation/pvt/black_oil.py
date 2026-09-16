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


#: Doymamış özlülük üstəlinin EHTİYAT qiyməti — `correlations.
#: saturated_undersaturated_oil_properties`-dəki `μ = μ_b·(p/Pb)^0.278`
#: ilə EYNİ ədəd. Cədvəldə doymamış sətir olmayanda işlədilir.
#:
#: ⚠️ ÖLÇÜLDÜ (Seans 34): SPE1 deck-inin qolları 0.4602 və 0.5085 verir,
#: yəni bu ehtiyat qiymət REAL deck üçün yanlışdır — ona görə işlədiləndə
#: xəbərdarlıq yazılır.
CORRELATION_VISCOSITY_EXPONENT = 0.278


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
        self._build_undersaturated_branch()

    # ═════════════════════════ doymamış qol — Bo(p, Rs) (B3-B)
    def _build_undersaturated_branch(self):
        """Rs → Pb tərsini və doymamış sıxılmanı bir dəfə hazırlayır.

        Bax `oil_fvf_undersaturated` — niyə lazım olduğu orada.
        """
        pressure = np.asarray(self.table.pressure, float)
        rs = np.asarray(self.table.solution_gor, float)
        pb = float(self.table.bubble_point)

        # ── Rs_sat-ın tərsi: yalnız Rs-in ARTAN hissəsi (Pb-yə qədər).
        # Pb-dən yuxarı Rs plato olduğu üçün oraya `np.interp` ilə
        # müraciət birmənalı olmazdı.
        increasing = np.concatenate(([True], np.diff(rs) > 0.0))
        self._rs_nodes = rs[increasing]
        self._pb_nodes = pressure[increasing]
        if self._rs_nodes.size >= 2:
            self._dpb_drs = (np.diff(self._pb_nodes) / np.diff(self._rs_nodes))
        else:                       # Rs sabitdir (ölü neft) — tərs yoxdur
            self._dpb_drs = np.zeros(0)

        # ── doymamış sıxılma c_o: cədvəlin ÖZ doymamış qolundan.
        # Cədvəl `Bo = Bo_b·exp(c_o·(Pb − p))` ilə qurulur (bax
        # `correlations.saturated_undersaturated_oil_properties`), yəni
        # ln(Bo) orada p-yə görə XƏTTİDİR və meyli tam `−c_o`-dur.
        above = pressure > pb
        bo = np.asarray(self._oil_fvf, float)
        if int(np.count_nonzero(above)) >= 2 and np.all(bo[above] > 0.0):
            slope = np.polyfit(pressure[above], np.log(bo[above]), 1)[0]
            self._undersaturated_co = float(max(-slope, 1e-6))
        else:
            # Cədvəldə doymamış qol yoxdur (Pb ≥ p_maks). Süxur-flüid
            # sıxılmasının neft hissəsi ən yaxın fiziki əvəzdir.
            fallback = float(np.max(self._co)) if np.size(self._co) else 0.0
            self._undersaturated_co = max(fallback, 1e-6)

        # ── doymamış özlülük üstəli n: μ = μ_b·(p/Pb)^n (G2)
        # ln μ ilə ln p arasındakı meyl məhz `n`-dir.
        mu = np.asarray(self.table.oil_viscosity, float)
        usable = above & (pressure > 0.0) & (mu > 0.0)
        if int(np.count_nonzero(usable)) >= 2:
            self._undersaturated_viscosity_exponent = float(np.polyfit(
                np.log(pressure[usable]), np.log(mu[usable]), 1)[0])
            self._viscosity_exponent_fitted = True
        else:
            self._undersaturated_viscosity_exponent = \
                CORRELATION_VISCOSITY_EXPONENT
            self._viscosity_exponent_fitted = False
            LOG.warning(
                "PVT: cədvəldə doyma təzyiqindən yuxarı iki sətir yoxdur — "
                "doymamış özlülük üstəli korrelyasiya qiymətinə (%.3f) "
                "düşdü. Ölçüldü: real deck-lərdə bu üstəl 0.46-0.51 ola "
                "bilir, yəni nəticə OLDUĞUNDAN ZƏİF özlülük artımı verir.",
                CORRELATION_VISCOSITY_EXPONENT)

    def saturation_pressure(self, rs) -> np.ndarray:
        """Pb(Rs) — verilmiş həll olmuş qazın doyma təzyiqi.

        `Rs_sat(p)` artan olduğu üçün tərsi parçalı xətti interpolyasiya
        ilə DƏQİQ alınır. Doymuş hüceyrədə `Rs = Rs_sat(p)` olduğundan
        bu funksiya elə `p`-nin özünü qaytarır — doymuş/doymamış keçid
        beləcə KƏSİLMƏZ olur.
        """
        rs = np.asarray(rs, float)
        if self._rs_nodes.size < 2:
            return np.full(rs.shape, float(self.table.bubble_point))
        return np.interp(rs, self._rs_nodes, self._pb_nodes)

    def _saturation_pressure_slope(self, rs) -> np.ndarray:
        """dPb/dRs — `saturation_pressure`-in parçalı meyli.

        Tərsin ÖZ meylindən çıxarılır (`1/Rs_sat'` kimi ayrıca
        hesablanmır) ki, sonlu fərqlə bitə-bit uyğun gəlsin.
        """
        rs = np.atleast_1d(np.asarray(rs, float))
        if self._dpb_drs.size == 0:
            return np.zeros_like(rs)
        index = np.clip(np.searchsorted(self._rs_nodes, rs, side="right") - 1,
                        0, self._dpb_drs.size - 1)
        slope = self._dpb_drs[index]
        outside = (rs < self._rs_nodes[0]) | (rs > self._rs_nodes[-1])
        return np.where(outside, 0.0, slope)

    def oil_fvf_undersaturated(self, pressure, rs) -> np.ndarray:
        """Bo(p, Rs) — doymamış qol, sənaye standartı (Eclipse `PVTO`).

        NİYƏ LAZIMDIR (ölçülüb, `ISH_HESABATI.md` → Seans 10). Üç fazalı
        mühərrik Bo-nu `oil_fvf(p)` ilə, yəni cədvəlin DOYMUŞ qolundan
        oxuyurdu — hüceyrənin doymamış olub-olmamasından ASILI OLMAYARAQ.
        Bu, termodinamik olaraq yanlışdır: doymamış hüceyrədə neftin
        tərkibi sabitdir (`Rs` sərbəst dəyişəndir və `Rs < Rs_sat(p)`),
        ona görə Bo həmin `Rs`-in doyma təzyiqindən başlayan SIXILMA
        qoluna aiddir, `Rs_sat(p)`-in doymuş qoluna yox.

        İki nəticəsi vardı:

          1. `dBo/dp` doymuş qolda Pb-də işarə dəyişir və SIFIRDAN keçir
             → neft tənliyinin təzyiq diaqonalı (`−So·B'o/Bo²`) itir,
             Jakobian kilidlənir, Nyuton donur (B3-A-nın üç fazalı
             analoqu — Pb ≥ 240 bar-da yığılma alınmırdı).
          2. Neft tənliyi 3-cü dəyişəndən (doymamışda `Rs`) HEÇ asılı
             olmurdu (`∂N_o/∂Rs = 0`) — ona görə qaz tənliyi 1-ci
             problemi kompensasiya EDƏ BİLMİRDİ.

        DÜZGÜN QOL:

            Bo(p, Rs) = Bo_sat(Pb(Rs)) · exp(c_o · (Pb(Rs) − p))

        Bu, cədvəlin öz doymamış qolu ilə EYNİ düsturdur — orada
        `Rs = Rs_sat(Pb_cədvəl)` xüsusi halıdır. Yəni Pb-dən yuxarıda
        nəticə DƏYİŞMİR; dəyişən yalnız Pb-dən aşağı DOYMAMIŞ
        hüceyrələrdir, harada ki əvvəl yanlış qol işlədilirdi.

        Fayda: `dBo/dp = −c_o·Bo < 0` HƏR YERDƏ (kilidlənmə yoxdur) və
        `dBo/dRs ≠ 0` (neft tənliyi 3-cü dəyişənə bağlanır). `Rs`-in
        ÖZÜ toxunulmur — qaz kütlə balansı POZULMUR.
        """
        pressure = np.asarray(pressure, float)
        pb = self.saturation_pressure(rs)
        bo_at_pb = self._interp(self._oil_fvf, pb)
        return bo_at_pb * np.exp(self._undersaturated_co * (pb - pressure))

    def oil_fvf_undersaturated_derivatives(self, pressure, rs):
        """`(∂Bo/∂p, ∂Bo/∂Rs)` — `oil_fvf_undersaturated`-in törəmələri.

        ∂Bo/∂p  = −c_o·Bo
        ∂Bo/∂Rs = (dPb/dRs)·Bo·(B'o_sat(Pb)/Bo_sat(Pb) + c_o)
        """
        pressure = np.asarray(pressure, float)
        pb = self.saturation_pressure(rs)
        bo_at_pb = self._interp(self._oil_fvf, pb)
        bo = bo_at_pb * np.exp(self._undersaturated_co * (pb - pressure))

        d_dp = -self._undersaturated_co * bo
        safe = np.where(bo_at_pb > 0.0, bo_at_pb, 1.0)
        d_drs = (self._saturation_pressure_slope(rs) * bo
                 * (self._slope("oil_fvf", pb) / safe
                    + self._undersaturated_co))
        return d_dp, d_drs

    # ═══════════════════════ doymamış qol — μo(p, Rs) (G2)
    def oil_viscosity_undersaturated(self, pressure, rs) -> np.ndarray:
        """μo(p, Rs) — doymamış qol (Eclipse `PVTO`-nun davam sətirləri).

        Bo qolunun (`oil_fvf_undersaturated`) GÜZGÜSÜDÜR:

            μo(p, Rs) = μo_sat(Pb(Rs)) · (p / Pb(Rs))^n

        NİYƏ LAZIMDIR: doymamış hüceyrədə neftin tərkibi sabitdir (Rs
        sərbəst dəyişəndir), ona görə özlülük həmin Rs-in doyma
        təzyiqindən başlayan SIXILMA qoluna aiddir — cədvəlin doymuş
        qolundakı `μo_sat(p)`-yə YOX. Doymuş qoldan oxumaq doymamış
        neftin təzyiqlə QATILAŞMASINI tamamilə itirir (SPE1-də 0.51 →
        0.74 cP, yəni 45 %).

        `Rs = Rs_sat(p)` (doymuş hüceyrə) olduqda `Pb(Rs) = p` və nəticə
        elə `μo_sat(p)`-nin özüdür — keçid KƏSİLMƏZDİR.
        """
        pressure = np.asarray(pressure, float)
        pb = self.saturation_pressure(rs)
        mu_at_pb = self._interp(self.table.oil_viscosity, pb)
        ratio = np.maximum(pressure / np.maximum(pb, 1e-12), 1e-12)
        return mu_at_pb * ratio ** self._undersaturated_viscosity_exponent

    def oil_viscosity_undersaturated_derivatives(self, pressure, rs):
        """`(∂μo/∂p, ∂μo/∂Rs)` — yuxarıdakı düsturun törəmələri.

            ∂μo/∂p  = n·μo / p
            ∂μo/∂Rs = (dPb/dRs) · [ μ'o_sat(Pb)·(p/Pb)^n − n·μo/Pb ]

        İkinci hədd Bo-dakı ilə eyni quruluşdadır: Pb dəyişəndə HƏM
        anchor (μo_sat(Pb)), HƏM də nisbət (p/Pb) dəyişir.
        """
        pressure = np.asarray(pressure, float)
        pb = self.saturation_pressure(rs)
        safe_pb = np.maximum(pb, 1e-12)
        exponent = self._undersaturated_viscosity_exponent
        mu_at_pb = self._interp(self.table.oil_viscosity, pb)
        ratio = np.maximum(pressure / safe_pb, 1e-12)
        mu = mu_at_pb * ratio ** exponent

        d_dp = exponent * mu / np.maximum(pressure, 1e-12)
        d_drs = self._saturation_pressure_slope(rs) * (
            self._slope("oil_viscosity", pb) * ratio ** exponent
            - exponent * mu / safe_pb)
        return d_dp, d_drs

    @property
    def undersaturated_viscosity_exponent(self) -> float:
        """Doymamış özlülük üstəli — cədvəldən fit olunub, yoxsa ehtiyat."""
        return self._undersaturated_viscosity_exponent

    @property
    def viscosity_exponent_fitted(self) -> bool:
        """`True` — üstəl cədvəlin ÖZ doymamış sətirlərindən alınıb."""
        return self._viscosity_exponent_fitted

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
