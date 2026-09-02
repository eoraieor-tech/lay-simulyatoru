"""Axın diskretizasiyası interfeysi — TPFA və gələcək MPFA-O üçün ORTAQ müqavilə.

Niyə: Nyuton/qalıq (residual) qatı KONKRET diskretizasiya sxemindən (TPFA
və ya gələcək MPFA-O) ASILI OLMAMALIDIR — yalnız bu müqavilədən. Hazırda
YEGANƏ implementasiya `imex2d.simulation.discretization.
TwoPointFluxDiscretization`-dır; gələcək MPFA-O eyni bu ABC-ni tətbiq
edəcək. `ResidualAssembler` (bax `simulation/implicit/residual.py`) bu
obyektin hansı sinifdən gəldiyini heç vaxt bilmir — yalnız `build()`-in
qaytardığı obyektin `connections`/`pore_volume`/`cell_volume`/
`compute_flux()` müqaviləsinə güvənir (bax `simulation/discretization.py`
modul docstring-i, tam mərhələ diaqramı üçün).

Bu fayl YALNIZ müqavilə təyin edir — `imex2d/interfaces/providers.py`-in
eyni fəlsəfəsi (Dependency Inversion). MPFA-O-nun ÖZÜ bu fazada
İMPLEMENTASİYA EDİLMİR.
"""

from __future__ import annotations
from abc import ABC, abstractmethod


class IFluxDiscretization(ABC):
    """`ReservoirModel`-dən axın diskretizasiyası qurur.

    `build()`-in qaytardığı obyekt (bax `DiscretizedGrid`) ən azı bunları
    daşımalıdır:

        connections   — `domain.grid.Connections` (üzlərin hüceyrə cütləri)
        pore_volume   — (ncell,) hüceyrə üzrə məsamə həcmi
        cell_volume   — (ncell,) hüceyrə üzrə tam həcm
        compute_flux(d_phi) — potensial fərqindən (və ya gələcək MPFA
            üçün — çoxnöqtəli stensil ağırlıqlı cəmindən) BAZA Darcy
            axınını hesablayır, HƏLƏ mobilite/upstream çəkiləndirmə
            OLMADAN (bax `ResidualAssembler.face_fluxes`-in bunu necə
            işlətdiyi — mobilite YENƏ DƏ residual qatında qalır, çünki
            bu, diskretizasiya sxemindən ASILI OLMAYAN fizikadır).
    """

    @abstractmethod
    def build(self, model) -> "object":
        """`ReservoirModel` -> diskretizasiya olunmuş grid obyekti.

        Qaytarılan obyektin dəqiq sinfi diskretizasiya sxeminə görə
        fərqli ola bilər (TPFA üçün `DiscretizedGrid`), AMMA yuxarıdakı
        müqaviləyə (duck-typing) əməl etməlidir ki, `ResidualAssembler`
        dəyişməsin.
        """

    # ── Phase 5A ƏLAVƏSİ (GERİYƏ UYĞUN) ─────────────────────────────────
    def supports_multipoint_stencil(self) -> bool:
        """Bu sxem ÇOXNÖQTƏLİ stensil (bir üzdə 2-dən çox hüceyrə)
        qururmu?

        DEFOLT `False` — ona görə MÖVCUD implementasiyalar (`TwoPoint
        FluxDiscretization`) HEÇ DƏYİŞMƏDƏN doğru cavab verir; bu metod
        ABSTRAKT DEYİL, geriyə uyğunluq POZULMUR (tapşırıq §23).

        `True` qaytaran sxem (`MPFAODiscretization`, bax
        `imex2d/discretization/mpfa_o.py`) üçün:
          - `compute_flux(d_phi)` müqaviləsi TƏTBİQ EDİLMİR (çoxnöqtəli
            axın tək üzün ΔΦ-sindən çıxarıla bilməz);
          - `build()`-in qaytardığı obyekt bunun ƏVƏZİNƏ TAM təzyiq
            vektoru qəbul edən bir metod verir (`compute_flux_from_
            pressure`), və `JacobianAssembler._build_pattern`-in 2-
            hüceyrəli blok fərziyyəsi Phase 5B-də genişləndirilməlidir
            (bax `simulation/discretization.py` modul docstring-i).

        `ResidualAssembler`/`JacobianAssembler` BU FAZADA bu metodu
        ÇAĞIRMIR — o, Phase 5B inteqrasiyasının AÇIQ giriş nöqtəsidir.
        """
        return False
