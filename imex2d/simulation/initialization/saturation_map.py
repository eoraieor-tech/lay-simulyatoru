"""İlkin doyumluluq GEOLOGİYA XƏRİTƏSİNDƏN — `property_maps["SW"]`.

Geologiya cədvəlindəki `Sw` sütunu interpolyasiya olunub modelə
`SW` xəritəsi kimi düşür (`geology_adapter` → `geology_service` →
`model_builder`), lakin ƏVVƏLLƏR onu HEÇ BİR mühərrik oxumurdu:
ilkin doyumluluq həmişə `InitialConditions.water_saturation`
SKALYARINDAN qurulurdu. Nəticədə OOIP (`Σ pv·(1−Sw)/Bo`) və RF
istifadəçinin verdiyi Sw məlumatını əks etdirmirdi.

Bu modul həmin boşluğu MÜHƏRRİKLƏRƏ TOXUNMADAN doldurur: hər iki
mühərrik onsuz da `IInitializationProvider`-ə hörmət edir, ona görə
dəyişiklik yalnız provider qatındadır.

Bayraq: `InitialConditions.use_saturation_map` (defolt `False` —
köhnə davranış). Prioritet cədvəli üçün bax `domain/initial.py`.
"""

from __future__ import annotations

import numpy as np

from ...domain.reservoir_model import ReservoirModel
from ...interfaces.providers import IInitializationProvider, InitialState
from ...logging_setup import get_logger

LOG = get_logger(__name__)

#: Geologiya cədvəlindəki `Sw` sütununun xəritə adı
#: (`application/geology_adapter.py::PROPERTY_MAP`).
SATURATION_MAP_NAME = "SW"


def water_saturation_from_map(model: ReservoirModel) -> np.ndarray:
    """`ReservoirModel.property_maps["SW"]`-dən ilkin doyumluluq massivi.

    Xəta halları (hamısı AÇIQ `ValueError`, səssiz düzəliş YOXDUR):
      · "SW" xəritəsi yoxdur
      · ölçü != model.ncell
      · massivdə NaN/sonsuz dəyər var (AKTİV hüceyrələrdə)
      · dəyərlər [0, 1] intervalından kənardadır

    QEYRİ-AKTİV hüceyrələr istisnadır: orada NaN normaldır (bax
    `domain/grid.py`, "KONVENSİYA"), ona görə yoxlanmır və nəticədə
    `ic.water_saturation` skalyarı ilə doldurulur — massiv QLOBAL
    (`ncell`) formatda qalsın deyə.

    SCAL hədlərinə (`swc … 1−Sor`) KƏSMƏ BURADA EDİLMİR — mühərriklər
    onu onsuz da `relperm.saturation_limits()` ilə edir. Neçə hüceyrənin
    kəsiləcəyini saymaq üçün bax `count_outside_scal_limits`.
    """
    prop = (model.property_maps or {}).get(SATURATION_MAP_NAME)
    if prop is None:
        raise ValueError(
            f"İlkin doyumluluq xəritəsi ({SATURATION_MAP_NAME}) modeldə yoxdur — "
            "`use_saturation_map` açıqdır, amma geologiya cədvəlindəki Sw sütunu "
            "interpolyasiya edilməyib.")

    values = np.asarray(prop.values, dtype=float).ravel()
    if values.size != model.ncell:
        raise ValueError(
            f"{SATURATION_MAP_NAME} xəritəsinin ölçüsü grid ilə uyğun gəlmir: "
            f"{values.size} != {model.ncell}.")

    active = model.grid.active.mask
    finite = np.isfinite(values)
    bad = active & ~finite
    if bad.any():
        raise ValueError(
            f"{SATURATION_MAP_NAME} xəritəsində {int(bad.sum())} AKTİV hüceyrədə "
            "NaN/sonsuz dəyər var — interpolyasiya bu hüceyrələri örtməyib. "
            "Səssiz doldurma edilmir: ya geologiya məlumatını genişləndirin, "
            "ya da həmin hüceyrələri ACTNUM ilə söndürün.")

    out_of_range = active & ((values < 0.0) | (values > 1.0))
    if out_of_range.any():
        worst = values[out_of_range]
        raise ValueError(
            f"{SATURATION_MAP_NAME} xəritəsində {int(out_of_range.sum())} hüceyrə "
            f"[0, 1] intervalından kənardadır (min {worst.min():.3f}, "
            f"maks {worst.max():.3f}).")

    sw = values.copy()
    if not active.all():
        sw[~active] = float(model.initial_conditions.water_saturation)
    return sw


def count_outside_scal_limits(sw: np.ndarray, model: ReservoirModel) -> int:
    """Mühərrikin SCAL hədlərinə KƏSƏCƏYİ aktiv hüceyrələrin sayı.

    Kəsmənin özü mühərrikdədir (`relperm.saturation_limits()`); burada
    yalnız sayılır ki, istifadəçi xəritəsinin nə qədərinin dəyişəcəyini
    jurnalda görsün.
    """
    scal = model.scal_parameters
    low, high = scal.swc, 1.0 - scal.sor
    active = model.grid.active.mask
    return int((active & ((sw < low - 1e-12) | (sw > high + 1e-12))).sum())


def _log_map_statistics(model: ReservoirModel, sw: np.ndarray, source: str) -> None:
    active = model.grid.active.mask
    values = sw[active]
    clipped = count_outside_scal_limits(sw, model)
    LOG.info("İlkin Sw %s xəritəsindən (%s): min %.3f · orta %.3f · maks %.3f"
             " | SCAL hədlərindən kənarda %d hüceyrə (mühərrik kəsəcək)",
             SATURATION_MAP_NAME, source, float(values.min()),
             float(values.mean()), float(values.max()), clipped)


class SaturationMapInitializationProvider(IInitializationProvider):
    """Bərabər paylanmış təzyiq + `SW` xəritəsindən doyumluluq.

    `use_equilibration` SÖNDÜRÜLÜB, `use_saturation_map` AÇIQ olanda
    işlədilir: təzyiq əvvəlki kimi `datum_pressure` ilə bərabərdir,
    yalnız doyumluluq mənbəyi dəyişir.
    """

    def initialize(self, model: ReservoirModel) -> InitialState:
        sw = water_saturation_from_map(model)
        _log_map_statistics(model, sw, "bərabər təzyiq")
        pressure = np.full(model.ncell,
                           float(model.initial_conditions.datum_pressure))
        return InitialState(pressure=pressure, water_saturation=sw)


class SaturationMapOverride(IInitializationProvider):
    """Başqa provider-i bükür: TƏZYİQİ ondan alır, DOYUMLULUĞU xəritədən.

    Hər iki bayraq açıq olanda işlədilir. Təzyiq fiziki MODELDƏN
    (hidrostatika) gəlir, doyumluluq isə ÖLÇÜLMÜŞ məlumatdır — ölçmə
    modelin üstündədir, ona görə xəritə bükülmüş provider-in
    kontakt/keçid zonası hesabını ƏVƏZ EDİR. İstifadəçi keçid zonasının
    söndürüldüyünü bilməlidir, ona görə jurnala açıq INFO yazılır.
    """

    def __init__(self, inner: IInitializationProvider):
        self.inner = inner

    def initialize(self, model: ReservoirModel) -> InitialState:
        inner_state = self.inner.initialize(model)
        sw = water_saturation_from_map(model)
        _log_map_statistics(model, sw, "hidrostatik təzyiq")
        LOG.info("Kontakt/keçid zonası hesabı ƏVƏZ OLUNDU: təzyiq %s "
                 "provider-indən, doyumluluq %s xəritəsindən.",
                 type(self.inner).__name__, SATURATION_MAP_NAME)
        return InitialState(pressure=inner_state.pressure, water_saturation=sw)
