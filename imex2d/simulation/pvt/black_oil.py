"""BlackOilPVTProvider — IPVTProvider-in ilk implementasiyası.

Yeganə hesablama üsulu: cədvəl üzrə xətti interpolyasiya (np.interp).
Cədvəldən kənarda sərhəd dəyəri saxlanılır (np.interp-in defolt davranışı) —
bu, ekstrapolyasiyadan daha təhlükəsizdir.

Sıxılma cədvəldən ədədi törəmə ilə alınır və bir dəfə əvvəlcədən
hesablanır ki, hər zaman addımında yenidən hesablanmasın.
"""

from __future__ import annotations
from typing import Optional, Sequence

import numpy as np

from ...domain.pvt import PVTTable
from ...interfaces.providers import IPVTProvider
from ...logging_setup import get_logger

LOG = get_logger(__name__)


#: Doymamış özlülük üstəlinin EHTİYAT qiyməti — `correlations.
#: saturated_undersaturated_oil_properties`-dəki `μ = μ_b·(p/Pb)^0.278`
#: ilə EYNİ ədəd. Cədvəldə doymamış sətir olmayanda işlədilir.
#:
#: ⚠️ ÖLÇÜLDÜ (Seans 37): SPE1 deck-inin qolları 0.4602 və 0.5802 verir,
#: yəni bu ehtiyat qiymət REAL deck üçün yanlışdır — ona görə işlədiləndə
#: xəbərdarlıq yazılır.
CORRELATION_VISCOSITY_EXPONENT = 0.278


class BlackOilPVTProvider(IPVTProvider):

    def __init__(self, table: PVTTable,
                 dead_oil_below_bubble_point: bool = False,
                 oil_branches: Optional[Sequence] = None):
        """`dead_oil_below_bubble_point` — bax `_dead_oil_fvf`.

        `oil_branches` (G3) — deck-in `PVTO` qolları (məs.
        `io.pvt_io.OilBranch`; duck-typing: `solution_gor`, `pressure`,
        `formation_volume_factor`, `viscosity`). Verilibsə doymamış
        sıxılma `c_o` və özlülük üstəli `n` cədvəlin TƏK doymamış
        hissəsindən yox, HƏR QOLUN ÖZ sətirlərindən alınır və Rs üzrə
        interpolyasiya olunur — bax `_build_branch_parameters`.
        Cədvəl deck-in BÜTÜN doymuş qolunu daşımalıdır
        (`DeckPvt.to_pvt_table()` `reference_rs` olmadan).

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
        self._oil_branches = oil_branches
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

        self._branch_rs = None
        if self._oil_branches is not None:
            self._build_branch_parameters(pb)
            self._build_saturated_grid(pb)
            return

        # ── doymamış sıxılma c_o: cədvəlin ÖZ doymamış qolundan.
        # Cədvəl `Bo = Bo_b·exp(c_o·(Pb − p))` ilə qurulur (bax
        # `correlations.saturated_undersaturated_oil_properties`), yəni
        # ln(Bo) orada p-yə görə XƏTTİDİR və meyli tam `−c_o`-dur.
        above = pressure > pb
        bo = np.asarray(self._oil_fvf, float)
        self._bo_anchor = None
        self._mu_anchor = None
        if int(np.count_nonzero(above)) >= 2 and np.all(bo[above] > 0.0):
            slope, intercept = np.polyfit(pressure[above], np.log(bo[above]), 1)
            self._undersaturated_co = float(max(-slope, 1e-6))
            # LÖVBƏR (ölçülmüş düzəliş, Seans 36): `Bo_sat(Pb)`-ni cədvəldən
            # interpolyasiya ilə GÖTÜRMƏK OLMAZ — Pb adətən düyün deyil və
            # əyrinin orada SINIĞI var (Bo aşağıda artır, yuxarıda azalır).
            # Ölçüldü: iki tərəfli interpolyasiya Bo_b-ni 0.37 %, μ_b-ni
            # 2.27 % səhv verir. Doğru lövbər doymamış qolun ÖZ fit-inin
            # Pb-dəki qiymətidir — həmin fit onsuz da aparılır.
            self._bo_anchor = float(np.exp(intercept + slope * pb))
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
            slope_mu, intercept_mu = np.polyfit(
                np.log(pressure[usable]), np.log(mu[usable]), 1)
            self._undersaturated_viscosity_exponent = float(slope_mu)
            self._mu_anchor = float(np.exp(intercept_mu
                                           + slope_mu * np.log(max(pb, 1e-12))))
            self._viscosity_exponent_fitted = True
        else:
            self._undersaturated_viscosity_exponent = \
                CORRELATION_VISCOSITY_EXPONENT
            self._viscosity_exponent_fitted = False
            LOG.warning(
                "PVT: cədvəldə doyma təzyiqindən yuxarı iki sətir yoxdur — "
                "doymamış özlülük üstəli korrelyasiya qiymətinə (%.3f) "
                "düşdü. Ölçüldü: real deck-lərdə bu üstəl 0.46-0.58 ola "
                "bilir, yəni nəticə OLDUĞUNDAN ZƏİF özlülük artımı verir.",
                CORRELATION_VISCOSITY_EXPONENT)

        self._build_saturated_grid(pb)

    def _build_branch_parameters(self, pb: float):
        """G3 — `c_o(Rs)` və `n(Rs)` deck-in HƏR qolundan.

        PROBLEM (ölçüldü, Seans 37, SPE1CASE2.DATA): qolların parametrləri
        FƏRQLİDİR — c_o 2.056e-4 ↔ 1.832e-4 1/bar (~11 %), n 0.460 ↔ 0.580.
        Tək doymamış hissəli cədvəl isə ya Rs-i 1.27-də kəsirdi (qaz vurulan
        hüceyrədə Bo −9 %, μo +17 %), ya da c_o ehtiyat qiymətə düşürdü
        (Bo −80 %).

        HƏR QOL ÜÇÜN (doymuş başı `pressure[0]`-dan keçən ən kiçik kvadratlar;
        iki sətirli qolda deck sətrini DƏQİQ təkrarlayır):

            ln(Bo/Bo_b) = −c_o·(p − Pb)        ln(μ/μ_b) = n·ln(p/Pb)

        QOLLAR ARASINDA — Rs üzrə parçalı xətti interpolyasiya; ən kənar
        doymamış qollardan kənarda sabit (törəmə sıfır — ekstrapolyasiya
        YOX, Q-28 Qərar 4 ilə eyni qayda).

        ⏳ FƏRZİYYƏ: doymamış sətri olmayan qollar üçün (SPE1-də Rs < 1.27)
        ən yaxın doymamış qolun parametri götürülür. OPM/Eclipse-in bu
        haldakı qaydası mənbədən YOXLANILMAYIB.
        """
        rs_values, co_values, n_values = [], [], []
        heads = []
        for branch in self._oil_branches:
            p = np.asarray(branch.pressure, float).ravel()
            bo = np.asarray(branch.formation_volume_factor, float).ravel()
            mu = np.asarray(branch.viscosity, float).ravel()
            rs_branch = float(branch.solution_gor)
            heads.append((rs_branch, p[0], bo[0], mu[0]))
            if p.size < 2:
                continue
            dp = p[1:] - p[0]
            log_bo = np.log(bo[1:] / bo[0])
            co = -float(np.dot(dp, log_bo) / np.dot(dp, dp))
            if not co > 0.0:
                raise ValueError(
                    f"PVTO qolu (Rs = {rs_branch:.4g}): doymamış Bo təzyiqlə "
                    f"azalmır (c_o = {co:.3g}) — fiziki deyil.")
            log_p = np.log(p[1:] / p[0])
            n = float(np.dot(log_p, np.log(mu[1:] / mu[0]))
                      / np.dot(log_p, log_p))
            rs_values.append(rs_branch)
            co_values.append(co)
            n_values.append(n)

        if not rs_values:
            raise ValueError("PVTO: heç bir qolda doymamış sətir yoxdur — "
                             "c_o və n qollardan alına bilməz.")
        order = np.argsort(rs_values)
        self._branch_rs = np.asarray(rs_values, float)[order]
        if np.any(np.diff(self._branch_rs) <= 0.0):
            raise ValueError("PVTO: doymamış qolların Rs-i təkrarlanır.")
        self._branch_co = np.asarray(co_values, float)[order]
        self._branch_n = np.asarray(n_values, float)[order]
        self._branch_co_slope = (np.diff(self._branch_co)
                                 / np.diff(self._branch_rs))
        self._branch_n_slope = (np.diff(self._branch_n)
                                / np.diff(self._branch_rs))

        # Qolların doymuş başları cədvəlin doymuş əyrisində olmalıdır —
        # əks halda qollar BAŞQA cədvələ aiddir (səssiz qarışıqlıq olmasın).
        pressure = np.asarray(self.table.pressure, float)
        rs_table = np.asarray(self.table.solution_gor, float)
        head_rs, head_pb, head_bo, head_mu = (np.array(column) for column in
                                              zip(*sorted(heads)))
        inside = head_pb <= pb * (1.0 + 1e-12)
        mismatch = np.abs(np.interp(head_pb[inside], pressure, rs_table)
                          - head_rs[inside])
        if np.any(mismatch > 1e-6 * np.maximum(head_rs[inside], 1.0)):
            raise ValueError("PVTO qolları cədvəlin doymuş əyrisi ilə uyğun "
                             "gəlmir (Rs_sat(Pb) ≠ qolun Rs-i).")
        if np.any(~inside):
            LOG.warning(
                "PVT: cədvəlin doymuş qolu %.2f bar-da bitir, deck-də isə "
                "Rs = %.4g-ə qədər qol var — Rs həmin nöqtədə KƏSİLƏCƏK. "
                "Cədvəli `DeckPvt.to_pvt_table()` ilə `reference_rs` "
                "olmadan qurun.", pb, float(head_rs[-1]))

        # Lövbər — deck-in DƏQİQ doymuş başlarından (Pb burada düyündür,
        # interpolyasiya sınığı kəsmir; bax Seans 36).
        self._bo_anchor = float(np.interp(pb, head_pb, head_bo))
        self._mu_anchor = float(np.interp(pb, head_pb, head_mu))
        top_rs = float(self._rs_nodes[-1]) if self._rs_nodes.size else pb
        self._undersaturated_co = float(np.interp(top_rs, self._branch_rs,
                                                  self._branch_co))
        self._undersaturated_viscosity_exponent = float(
            np.interp(top_rs, self._branch_rs, self._branch_n))
        self._viscosity_exponent_fitted = True

    def _branch_parameter(self, values: np.ndarray, slopes: np.ndarray, rs):
        """`(qiymət, d/dRs)` — qollar üzrə parçalı xətti interpolyasiya."""
        rs = np.atleast_1d(np.asarray(rs, float))
        value = np.interp(rs, self._branch_rs, values)
        if slopes.size == 0:
            return value, np.zeros_like(rs)
        index = np.clip(np.searchsorted(self._branch_rs, rs, side="right") - 1,
                        0, slopes.size - 1)
        outside = (rs < self._branch_rs[0]) | (rs > self._branch_rs[-1])
        return value, np.where(outside, 0.0, slopes[index])

    def _compressibility(self, rs):
        """`(c_o, dc_o/dRs)` — qollar yoxdursa tək qiymət və `None`."""
        if self._branch_rs is None:
            return self._undersaturated_co, None
        return self._branch_parameter(self._branch_co, self._branch_co_slope, rs)

    def _viscosity_exponent(self, rs):
        """`(n, dn/dRs)` — qollar yoxdursa tək qiymət və `None`."""
        if self._branch_rs is None:
            return self._undersaturated_viscosity_exponent, None
        return self._branch_parameter(self._branch_n, self._branch_n_slope, rs)

    def _build_saturated_grid(self, pb: float):
        """DOYMUŞ qol üçün düzəldilmiş şəbəkə — lövbərin mənbəyi.

        PROBLEM (ölçüldü, Seans 36): `Bo_sat(Pb)`/`μo_sat(Pb)` cədvəldən
        adi interpolyasiya ilə götürüləndə Pb ADƏTƏN DÜYÜN OLMUR və
        interpolyasiya əyrinin SINIĞINI kəsir (μo Pb-dən aşağı azalır,
        yuxarı artır). Nəticə: μ lövbəri 2.27 %, Bo lövbəri 0.37 % səhv.

        İLK CƏHD SƏHV İDİ: lövbəri `np.where` ilə "Pb-də fit, aşağıda
        interpolyasiya" kimi seçmək SIÇRAYIŞ yaradırdı — sonlu fərq
        ∂/∂Rs üçün 27 kimi qiymətlər verirdi (kəsilməzlik Nyuton üçün
        2 %-lik meyldən DAHA VACİBDİR).

        DÜZGÜN HƏLL: şəbəkəni Pb-də KƏSMƏK və oraya fit-dən gələn dəqiq
        nöqtəni qoymaq. Onda lövbər həm sınığı kəsmir, həm də `pb` üzrə
        kəsilməz qalır; törəməsi isə elə bu əyrinin öz meylidir.
        """
        pressure = np.asarray(self.table.pressure, float)
        tolerance = 1e-9 * max(abs(pb), 1.0)

        # ŞƏBƏKƏ Pb-DƏ KƏSİLMİR (2-ci düzəliş, ölçülmüş səbəb): əvvəlki
        # variantda şəbəkə Pb-də bitirdi və `np.interp` ondan yuxarı sabit
        # qalırdı. Doymamış hüceyrələrin ƏKSƏRİYYƏTİ məhz Pb-də oturur, ona
        # görə sonlu fərq Rs-i artıranda lövbər "donurdu" və akkumulyasiya
        # Jakobianı 4.4 dəfə səhv çıxırdı (ölçüldü). İndi şəbəkə tam
        # diapazondadır: Pb-dən AŞAĞI cədvəlin öz düyünləri, Pb-dən YUXARI
        # isə fit əyrisinin öz qiymətləri — yəni hər iki tərəfdə hamardır.
        grid = np.unique(np.concatenate([pressure, [pb]]))
        below = grid < pb - tolerance
        above = grid > pb + tolerance

        def augmented(values, anchor):
            """DOYMUŞ əyri `pb`-nin funksiyası kimi.

            Pb-yə qədər — cədvəlin öz düyünləri (sınıq kəsilmir, çünki
            şəbəkə orada bitir). Pb-də — fit-dən gələn dəqiq lövbər.
            Pb-dən YUXARI — həmin DOYMUŞ meylin davamı.

            NİYƏ DOYMUŞ MEYL (2-ci cəhdin səhvi): oraya doymamış qolun
            qiymətlərini yazmaq lövbərin meylini `−c_o` edirdi və törəmə
            düsturundakı `+c_o` ilə tam kompensasiya olunurdu — analitik
            ∂/∂Rs sıfıra düşürdü (ölçüldü: SF 4.2e-3, analitik 7.7e-7).
            Lövbər DOYMUŞ qiymətdir, ona görə meyli də doymuş meyl olmalıdır.
            Yuxarı davam yalnız ona görə lazımdır ki, `pb` sərhəddə olanda
            mərkəzi sonlu fərq simmetrik qalsın (şəbəkə "donmasın").
            """
            values = np.asarray(values, float)
            merged = np.interp(grid, pressure, values)
            if anchor is None:          # fit yoxdur — köhnə davranış
                return merged
            merged[np.isclose(grid, pb, atol=tolerance)] = float(anchor)
            if np.any(below) and np.any(above):
                last = grid[below][-1]
                slope = (float(anchor) - merged[below][-1]) / max(pb - last, 1e-12)
                merged[above] = float(anchor) + slope * (grid[above] - pb)
            return merged

        self._sat_pressure = grid
        self._sat_oil_fvf = augmented(self._oil_fvf, self._bo_anchor)
        self._sat_oil_viscosity = augmented(self.table.oil_viscosity,
                                            self._mu_anchor)

    def _anchor_at(self, values: np.ndarray, pb) -> np.ndarray:
        """Düzəldilmiş doymuş şəbəkə üzrə lövbər (bax `_build_saturated_grid`)."""
        return np.interp(np.asarray(pb, float), self._sat_pressure, values)

    def _anchor_slope(self, values: np.ndarray, pb) -> np.ndarray:
        """`_anchor_at`-in DƏQİQ törəməsi — eyni şəbəkənin interval meyli."""
        pb = np.atleast_1d(np.asarray(pb, float))
        nodes = self._sat_pressure
        if nodes.size < 2:
            return np.zeros_like(pb)
        slopes = np.diff(values) / np.diff(nodes)
        index = np.clip(np.searchsorted(nodes, pb, side="right") - 1,
                        0, slopes.size - 1)
        result = slopes[index]
        outside = (pb < nodes[0]) | (pb > nodes[-1])
        return np.where(outside, 0.0, result)

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
        bo_at_pb = self._anchor_at(self._sat_oil_fvf, pb)
        co, _ = self._compressibility(rs)
        return bo_at_pb * np.exp(co * (pb - pressure))

    def oil_fvf_undersaturated_derivatives(self, pressure, rs):
        """`(∂Bo/∂p, ∂Bo/∂Rs)` — `oil_fvf_undersaturated`-in törəmələri.

        ∂Bo/∂p  = −c_o·Bo
        ∂Bo/∂Rs = (dPb/dRs)·Bo·(B'o_sat(Pb)/Bo_sat(Pb) + c_o)
                  + Bo·(dc_o/dRs)·(Pb − p)          ← yalnız deck qolları ilə (G3)
        """
        pressure = np.asarray(pressure, float)
        pb = self.saturation_pressure(rs)
        bo_at_pb = self._anchor_at(self._sat_oil_fvf, pb)
        co, co_rs = self._compressibility(rs)
        bo = bo_at_pb * np.exp(co * (pb - pressure))

        d_dp = -co * bo
        safe = np.where(bo_at_pb > 0.0, bo_at_pb, 1.0)
        d_drs = (self._saturation_pressure_slope(rs) * bo
                 * (self._anchor_slope(self._sat_oil_fvf, pb) / safe + co))
        if co_rs is not None:
            d_drs = d_drs + bo * co_rs * (pb - pressure)
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
        mu_at_pb = self._anchor_at(self._sat_oil_viscosity, pb)
        ratio = np.maximum(pressure / np.maximum(pb, 1e-12), 1e-12)
        exponent, _ = self._viscosity_exponent(rs)
        return mu_at_pb * ratio ** exponent

    def oil_viscosity_undersaturated_derivatives(self, pressure, rs):
        """`(∂μo/∂p, ∂μo/∂Rs)` — yuxarıdakı düsturun törəmələri.

            ∂μo/∂p  = n·μo / p
            ∂μo/∂Rs = (dPb/dRs) · [ μ'o_sat(Pb)·(p/Pb)^n − n·μo/Pb ]
                      + μo·(dn/dRs)·ln(p/Pb)        ← yalnız deck qolları ilə (G3)

        İkinci hədd Bo-dakı ilə eyni quruluşdadır: Pb dəyişəndə HƏM
        anchor (μo_sat(Pb)), HƏM də nisbət (p/Pb) dəyişir.
        """
        pressure = np.asarray(pressure, float)
        pb = self.saturation_pressure(rs)
        safe_pb = np.maximum(pb, 1e-12)
        exponent, exponent_rs = self._viscosity_exponent(rs)
        mu_at_pb = self._anchor_at(self._sat_oil_viscosity, pb)
        ratio = np.maximum(pressure / safe_pb, 1e-12)
        mu = mu_at_pb * ratio ** exponent

        d_dp = exponent * mu / np.maximum(pressure, 1e-12)
        d_drs = self._saturation_pressure_slope(rs) * (
            self._anchor_slope(self._sat_oil_viscosity, pb) * ratio ** exponent
            - exponent * mu / safe_pb)
        if exponent_rs is not None:
            d_drs = d_drs + mu * exponent_rs * np.log(ratio)
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
