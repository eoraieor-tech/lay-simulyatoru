"""Uyğunlaşdırma parametrləri — C5, mərhələ 2.

History matching-in mahiyyəti: modelin naməlum parametrlərini elə
seçmək ki, hesablanmış tarixçə ölçülmüşə uyğun gəlsin. Bu modul
"nəyi dəyişdiririk" sualına cavab verir.

İKİ NÖV PARAMETR

    ÇARPAN (multiplier)  — mövcud sahəyə vurulur: PERMX × 1.4
                           Heterogenlik qorunur, yalnız səviyyə dəyişir.
    ƏVƏZLƏMƏ (absolute)  — dəyər birbaşa təyin olunur: Sor = 0.22

Keçiricilik və məsaməlilik üçün çarpan işlədilir, çünki geoloji model
sahənin FORMASINI verir (quyulardan interpolyasiya ilə), mütləq
səviyyəsi isə qeyri-müəyyəndir. SCAL parametrləri isə skalyardır.

ƏSAS PRİNSİP — dəyişməzlik

`apply()` baza modelini HEÇ VAXT dəyişdirmir, kopiya qaytarır. Əks
halda təkrar tətbiqdə çarpanlar üst-üstə yığılardı: 1.4 → 1.96 → 2.74.
Optimallaşdırma yüzlərlə dəfə tətbiq edir, ona görə bu, kritikdir.

LOQARİFMİK MİQYAS

Keçiricilik çarpanı üçün `log_scale=True` tövsiyə olunur: 0.5 və 2.0
fiziki cəhətdən simmetrik dəyişikliklərdir, xətti miqyasda isə
1.0-dan məsafələri fərqlidir (0.5 vs 1.0). Optimallaşdırıcı log
fəzasında daha balanslı hərəkət edir.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from enum import Enum
from typing import Callable, Dict, List, Optional

import numpy as np

from ..domain.reservoir_model import ReservoirModel


class ParameterKind(Enum):
    MULTIPLIER = "çarpan"
    ABSOLUTE = "mütləq"


@dataclass
class ParameterDefinition:
    """Bir uyğunlaşdırma parametri."""
    name: str
    apply_to: Callable[[ReservoirModel, float], None]
    minimum: float
    maximum: float
    initial: float = 1.0
    kind: ParameterKind = ParameterKind.MULTIPLIER
    log_scale: bool = False
    unit: str = ""
    description: str = ""

    def __post_init__(self):
        if self.minimum >= self.maximum:
            raise ValueError(f"{self.name}: minimum maksimumdan kiçik olmalıdır.")
        if self.log_scale and self.minimum <= 0.0:
            raise ValueError(f"{self.name}: log miqyasda minimum müsbət olmalıdır.")
        self.initial = float(np.clip(self.initial, self.minimum, self.maximum))

    # ─────────────────────────────── normallaşdırılmış [0, 1] fəza
    def to_unit(self, value: float) -> float:
        """Fiziki dəyər -> [0, 1]. Optimallaşdırıcı bu fəzada işləyir."""
        value = float(np.clip(value, self.minimum, self.maximum))
        if self.log_scale:
            return float((np.log(value) - np.log(self.minimum))
                         / (np.log(self.maximum) - np.log(self.minimum)))
        return float((value - self.minimum) / (self.maximum - self.minimum))

    def from_unit(self, unit_value: float) -> float:
        """[0, 1] -> fiziki dəyər."""
        unit_value = float(np.clip(unit_value, 0.0, 1.0))
        if self.log_scale:
            return float(np.exp(np.log(self.minimum) + unit_value
                                * (np.log(self.maximum) - np.log(self.minimum))))
        return float(self.minimum + unit_value * (self.maximum - self.minimum))

    def clip(self, value: float) -> float:
        return float(np.clip(value, self.minimum, self.maximum))


# ══════════════════════════════════════════════ hazır parametrlər

def _scale_property(attribute: str):
    def apply(model: ReservoirModel, value: float) -> None:
        prop = getattr(model.rock, attribute)
        if prop is not None:
            prop.values[:] = prop.values * value
    return apply


def _scale_horizontal_permeability(model: ReservoirModel, value: float) -> None:
    """PERMX və PERMY birlikdə — anizotropluq nisbəti qorunur."""
    model.rock.permx.values[:] = model.rock.permx.values * value
    model.rock.permy.values[:] = model.rock.permy.values * value


def _set_vertical_ratio(model: ReservoirModel, value: float) -> None:
    """Kv/Kh — PERMZ üfüqi keçiricilikdən yenidən qurulur."""
    if model.rock.permz is None:
        return
    horizontal = np.sqrt(model.rock.permx.values * model.rock.permy.values)
    model.rock.permz.values[:] = horizontal * value


def _set_scal(attribute: str):
    """SCAL parametrini təyin edir və ilkin doyumluluğu uyğunlaşdırır.

    Swc və ya Sor dəyişəndə hərəkətli doyumluluq intervalı sürüşür.
    İlkin Sw həmin intervaldan kənarda qalarsa, model YOXLAMADAN
    KEÇMİR və simulyasiya xəta atır.

    Optimallaşdırma zamanı bu ölümcüldür: optimallaşdırıcı Swc-ni
    sınayır, model rədd olunur, axtarış dayanır. Ona görə ilkin
    doyumluluq avtomatik olaraq yeni intervala salınır — bu, həm də
    fiziki cəhətdən doğrudur: bağlı su doyumluluğu dəyişəndə lay
    şəraiti də dəyişir.
    """
    def apply(model: ReservoirModel, value: float) -> None:
        setattr(model.scal_parameters, attribute, value)
        _reconcile_initial_saturation(model)
    return apply


def _reconcile_initial_saturation(model: ReservoirModel) -> None:
    scal = model.scal_parameters
    low, high = scal.swc, 1.0 - scal.sor
    if high <= low:                       # yararsız SCAL — yoxlama tutacaq
        return
    initial = model.initial_conditions
    initial.water_saturation = float(
        np.clip(initial.water_saturation, low, high))


def _set_capillary(attribute: str):
    def apply(model: ReservoirModel, value: float) -> None:
        setattr(model.capillary_parameters, attribute, value)
    return apply


def _set_oil_viscosity(model: ReservoirModel, value: float) -> None:
    model.fluids.oil_viscosity = value


def _set_contact(model: ReservoirModel, value: float) -> None:
    model.initial_conditions.oil_water_contact = value


def _scale_pore_volume(model: ReservoirModel, value: float) -> None:
    """Məsaməlilik çarpanı — ehtiyata birbaşa təsir edir."""
    model.rock.porosity.values[:] = np.clip(
        model.rock.porosity.values * value, 1e-4, 0.95)


def standard_parameters(model: ReservoirModel) -> List[ParameterDefinition]:
    """Tipik uyğunlaşdırma dəsti — modeldən asılı olaraq.

    Sıralama təsadüfi deyil: praktikada keçiricilik və məsamə həcmi
    ən böyük təsirə malikdir, SCAL isə cəbhənin formasını dəqiqləşdirir.
    """
    definitions = [
        ParameterDefinition(
            "PERM_MULT", _scale_horizontal_permeability,
            minimum=0.1, maximum=10.0, initial=1.0, log_scale=True,
            description="Üfüqi keçiriciliyin qlobal çarpanı"),
        ParameterDefinition(
            "PORO_MULT", _scale_pore_volume,
            minimum=0.7, maximum=1.4, initial=1.0,
            description="Məsaməliliyin çarpanı (ehtiyata təsir edir)"),
    ]

    if model.rock.permz is not None:
        definitions.append(ParameterDefinition(
            "KV_KH", _set_vertical_ratio,
            minimum=0.001, maximum=1.0, initial=0.1, log_scale=True,
            kind=ParameterKind.ABSOLUTE,
            description="Şaquli/üfüqi keçiricilik nisbəti"))

    scal = model.scal_parameters
    definitions += [
        ParameterDefinition(
            "SOR", _set_scal("sor"),
            minimum=0.05, maximum=0.45, initial=scal.sor,
            kind=ParameterKind.ABSOLUTE,
            description="Qalıq neft doyumluluğu"),
        ParameterDefinition(
            "SWC", _set_scal("swc"),
            minimum=0.05, maximum=0.40, initial=scal.swc,
            kind=ParameterKind.ABSOLUTE,
            description="Bağlı su doyumluluğu"),
        ParameterDefinition(
            "KRW_END", _set_scal("krw_end"),
            minimum=0.05, maximum=1.0, initial=scal.krw_end,
            kind=ParameterKind.ABSOLUTE,
            description="krw son nöqtəsi — su cəbhəsinin sürətini idarə edir"),
        ParameterDefinition(
            "COREY_NW", _set_scal("nw"),
            minimum=1.0, maximum=6.0, initial=scal.nw,
            kind=ParameterKind.ABSOLUTE,
            description="Corey su göstəricisi"),
        ParameterDefinition(
            "COREY_NO", _set_scal("no"),
            minimum=1.0, maximum=6.0, initial=scal.no,
            kind=ParameterKind.ABSOLUTE,
            description="Corey neft göstəricisi"),
    ]

    if model.pvt_table is None:
        definitions.append(ParameterDefinition(
            "MU_OIL", _set_oil_viscosity,
            minimum=0.3, maximum=50.0, initial=model.fluids.oil_viscosity,
            log_scale=True, kind=ParameterKind.ABSOLUTE, unit="cP",
            description="Neftin lözlüyü (PVT cədvəli olmadıqda)"))

    if model.initial_conditions.oil_water_contact is not None:
        depths = model.geometry.cell_depths()
        low, high = float(depths.min()), float(depths.max())
        # Düz layda bütün hüceyrələr eyni dərinlikdədir — kontaktı
        # dəyişdirmək mənasızdır (ya hamısı su, ya hamısı neft zonası).
        # Belə halda hədləri bir hüceyrə qalınlığı qədər genişləndiririk
        # ki, parametr heç olmasa "kontakt var / yox" seçimini versin.
        mean_dz = float(np.mean(model.geometry.dz))
        if high - low < mean_dz:
            margin = mean_dz
            low, high = low - margin, high + margin
        definitions.append(ParameterDefinition(
            "OWC", _set_contact, minimum=low, maximum=high,
            initial=float(np.clip(
                model.initial_conditions.oil_water_contact, low, high)),
            kind=ParameterKind.ABSOLUTE, unit="m",
            description="Su-neft kontaktının dərinliyi"))

    return definitions


# ══════════════════════════════════════════════════ tətbiq

@dataclass
class ParameterSet:
    """Parametrlərin dəsti + cari dəyərlər."""
    definitions: List[ParameterDefinition] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.definitions)

    @property
    def names(self) -> List[str]:
        return [item.name for item in self.definitions]

    @property
    def initial_values(self) -> np.ndarray:
        return np.array([item.initial for item in self.definitions], float)

    @property
    def bounds(self) -> List[tuple]:
        return [(item.minimum, item.maximum) for item in self.definitions]

    def to_unit(self, values: np.ndarray) -> np.ndarray:
        return np.array([definition.to_unit(value)
                         for definition, value in zip(self.definitions, values)])

    def from_unit(self, unit_values: np.ndarray) -> np.ndarray:
        return np.array([definition.from_unit(value)
                         for definition, value in zip(self.definitions,
                                                      unit_values)])

    def clip(self, values: np.ndarray) -> np.ndarray:
        return np.array([definition.clip(value)
                         for definition, value in zip(self.definitions, values)])

    def as_dict(self, values: np.ndarray) -> Dict[str, float]:
        return {definition.name: float(value)
                for definition, value in zip(self.definitions, values)}

    def describe(self, values: Optional[np.ndarray] = None) -> str:
        values = self.initial_values if values is None else values
        lines = []
        for definition, value in zip(self.definitions, values):
            lines.append(
                f"  {definition.name:<10} {value:>10.4f}  "
                f"[{definition.minimum:g} … {definition.maximum:g}]"
                f"{' log' if definition.log_scale else ''}"
                f"   {definition.description}")
        return "\n".join(lines)


class ModelModifier:
    """Parametr dəyərlərini modelə tətbiq edir — KOPİYA üzərində."""

    def __init__(self, base_model: ReservoirModel,
                 parameters: ParameterSet):
        self.base_model = base_model
        self.parameters = parameters

    def apply(self, values: np.ndarray) -> ReservoirModel:
        """Yeni model qaytarır; baza model toxunulmaz qalır."""
        if len(values) != len(self.parameters):
            raise ValueError(
                f"{len(values)} dəyər verildi, {len(self.parameters)} gözlənilir.")

        model = self._copy(self.base_model)
        for definition, value in zip(self.parameters.definitions, values):
            definition.apply_to(model, definition.clip(float(value)))
        return model

    def apply_unit(self, unit_values: np.ndarray) -> ReservoirModel:
        return self.apply(self.parameters.from_unit(unit_values))

    @staticmethod
    def _copy(model: ReservoirModel) -> ReservoirModel:
        """Dərin kopiya — massivlər paylaşılmamalıdır.

        Qonşuluq qrafı (`_connections`) yenidən qurulur, çünki o,
        keşdir və kopiyada köhnə obyektə istinad etməməlidir.
        """
        clone = copy.deepcopy(model)
        clone._connections = None
        return clone
