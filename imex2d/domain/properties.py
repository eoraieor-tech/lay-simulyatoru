"""Xassə xəritələri və süxur/flüid xassələri."""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from .unit_conversions import convert, known_units
from .validation import (ValidationResult, validate_compressibility, validate_density,
                         validate_permeability, validate_porosity, validate_viscosity)

#: Xassə adı -> gözlənilən kəmiyyət növü (bax `unit_conversions.py`).
#: YALNIZ real vahid-qarışıqlığı riski olan xassələr üçün (Phase 1
#: audit: keçiricilik mD/Darcy/m² arasında, təzyiq bar/psi arasında
#: real çaşqınlıq mənbəyidir). Reyestrdə OLMAYAN ad (PORO, NTG, SW,
#: REGION_ID...) üçün `unit` sərbəst mətn olaraq qalır — DƏYİŞMİR.
#: KXX/KYY/KZZ/KXY/KXZ/KYZ — Phase 2 (tam tenzor permeabilite): DİAQONAL
#: PERMX/PERMY/PERMZ İLƏ EYNİ vahid-etibarlılığı qoruyur (bax
#: `PermeabilityTensor`) — off-diaqonal komponentlər üçün AYRICA/zəif
#: yoxlama qalmasın deyə.
PROPERTY_QUANTITY = {
    "PERMX": "permeability", "PERMY": "permeability", "PERMZ": "permeability",
    "KXX": "permeability", "KYY": "permeability", "KZZ": "permeability",
    "KXY": "permeability", "KXZ": "permeability", "KYZ": "permeability",
    "PRESSURE": "pressure",
}


@dataclass
class PropertyMap:
    """Adlandırılmış hüceyrə massivi. Hər xassə öz adını və vahidini daşıyır.

    `unit` boş deyilsə (istifadəçi/idxal AÇIQ vahid göstərib) VƏ `name`
    `PROPERTY_QUANTITY`-də qeydiyyatdan keçibsə, vahid bu xassənin
    kəmiyyət növünə (məs. keçiricilik üçün mD/D/m²) uyğun OLMALIDIR —
    əks halda `ValueError` (məs. `PropertyMap("PERMX", ..., "psi")`).
    Boş `unit` (defolt) HEÇ VAXT rədd edilmir — bu, "vahid göstərilməyib"
    halıdır, mövcud mühərrik vahidi kimi qəbul edilir (bax `UNITS.md`).
    """
    name: str
    values: np.ndarray
    unit: str = ""

    def __post_init__(self):
        quantity = PROPERTY_QUANTITY.get(self.name)
        if quantity is not None and self.unit:
            allowed = known_units(quantity)
            if self.unit not in allowed:
                raise ValueError(
                    f"{self.name}: vahid {self.unit!r} bu xassə üçün etibarsızdır "
                    f"(gözlənilən kəmiyyət: {quantity}, dəstəklənən vahidlər: {allowed}).")

    @classmethod
    def uniform(cls, name: str, value: float, ncell: int, unit: str = "") -> "PropertyMap":
        return cls(name, np.full(ncell, float(value)), unit)

    @classmethod
    def from_array(cls, name: str, arr, ncell: int, unit: str = "") -> "PropertyMap":
        a = np.asarray(arr, dtype=float)
        if a.ndim == 0:
            return cls.uniform(name, float(a), ncell, unit)
        if a.size != ncell:
            raise ValueError(f"{name}: {a.size} dəyər, gözlənilən {ncell}")
        return cls(name, a.ravel().copy(), unit)

    def as_grid(self, shape) -> np.ndarray:
        return self.values.reshape(shape)

    def stats(self) -> dict:
        return {"min": float(self.values.min()), "max": float(self.values.max()),
                "mean": float(self.values.mean())}


@dataclass
class PropertyUncertainty:
    """Kəsilməz xassənin qeyri-müəyyənlik/keyfiyyət diaqnostikası (grid üzrə).

    `geology.property_interpolation.PropertyEstimate`-in DOMAIN-səviyyəli,
    yüngül əksidir — `domain` qatı `geology` alqoritm qatından ASILI OLA
    BİLMƏZ (qat sərhədi), ona görə zəngin obyekt deyil, sadə massivlər
    daşınır. Çevirən tərəf `application/geology_service.py`-dir.
    """
    name: str
    variance: np.ndarray
    std: np.ndarray
    confidence: np.ndarray          #: object massiv — "high"/"medium"/"low"/"extrapolated"
    support: np.ndarray             #: object massiv — Phase A `SUPPORT_*` həndəsi təsnifat
    neighbor_count: np.ndarray
    nearest_distance: np.ndarray
    data_density: np.ndarray
    extrapolated: np.ndarray
    variance_kind: str = ""
    warnings: List[str] = field(default_factory=list)

    @property
    def ncell(self) -> int:
        return int(self.variance.size)


@dataclass
class PropertyProvenance:
    """BİR xassənin hüceyrə-üzrə MƏNŞƏYİ — "bu ədəd haradan gəldi".

    `PropertyMap` (yalnız DƏYƏR) ilə YANAŞI saxlanılır, onu ƏVƏZ ETMİR:
    mövcud oxucular (`RockProperties`, TPFA/MPFA, Eclipse ixracı) heç nə
    bilmədən ƏVVƏLKİ kimi işləyir; mənşəyə ehtiyacı olan (validasiya,
    3D görüntü, hesabat) buradan oxuyur.

    SAHƏLƏR (tapşırıq §11 — ayrı-ayrı saxlanılır, birləşdirilmir):

        original      — dəyişdirilməmiş ilkin sahə (varsa)
        interpolated  — YALNIZ interpolyasiya olunmuş hüceyrələr, qalanı NaN
        estimated     — YALNIZ completion ilə doldurulanlar, qalanı NaN
        final         — modelə/simulyatora GEDƏN sahə
        status        — hər hüceyrə üçün `DataStatus` dəyəri (mətn)
        method        — hansı üsulla (məs. "kriging", "vertical_trend", "sgs")
        confidence    — `[0,1]` ORDİNAL dəstək balı, hesablanmayanda NaN

    `confidence_kind` AÇIQ şəkildə "ordinal_support_score"-dur:
    KALİBRLƏNMİŞ EHTİMAL DEYİL (bax `geology/property_interpolation.
    Confidence` docstring-i — eyni qayda).
    """

    name: str
    status: np.ndarray               #: object massiv — `DataStatus` dəyərləri (str)
    method: np.ndarray               #: object massiv — üsul adı (boş ola bilər)
    confidence: np.ndarray           #: float massiv — hesablanmayan hüceyrə NaN
    final: np.ndarray
    original: Optional[np.ndarray] = None
    interpolated: Optional[np.ndarray] = None
    estimated: Optional[np.ndarray] = None
    confidence_kind: str = "ordinal_support_score"
    layer_methods: Dict[int, str] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    def __post_init__(self):
        n = int(np.asarray(self.final).size)
        for label, array in (("status", self.status), ("method", self.method),
                             ("confidence", self.confidence)):
            if np.asarray(array).size != n:
                raise ValueError(
                    f"{self.name}: '{label}' ölçüsü final sahə ilə uyğun gəlmir "
                    f"({np.asarray(array).size} != {n})")

    @property
    def ncell(self) -> int:
        return int(np.asarray(self.final).size)

    def mask(self, *statuses: str) -> np.ndarray:
        """Verilmiş statuslu hüceyrələrin bool maskası."""
        wanted = {str(s) for s in statuses}
        return np.asarray([str(s) in wanted for s in self.status], dtype=bool)

    def status_counts(self) -> Dict[str, int]:
        names, counts = np.unique(np.asarray(self.status, dtype=object).astype(str),
                                  return_counts=True)
        return {str(name): int(count) for name, count in zip(names, counts)}


@dataclass
class CategoricalUncertainty:
    """Kateqorik xassənin qeyri-müəyyənlik diaqnostikası (grid üzrə).

    `geology.property_interpolation.CategoricalEstimate`-in DOMAIN-
    səviyyəli əksi — bax `PropertyUncertainty` docstring-i (eyni qat
    sərhədi qaydası)."""
    name: str
    categories: np.ndarray
    probabilities: np.ndarray       #: (ncell, k)
    entropy: np.ndarray
    normalized_entropy: np.ndarray
    max_probability: np.ndarray
    confidence: np.ndarray
    support: np.ndarray              #: object massiv — Phase A `SUPPORT_*` həndəsi təsnifat
    neighbor_count: np.ndarray
    nearest_distance: np.ndarray
    extrapolated: np.ndarray
    n_probability_corrections: int = 0
    warnings: List[str] = field(default_factory=list)

    @property
    def ncell(self) -> int:
        return int(self.max_probability.size)


@dataclass
class PermeabilityTensor:
    """Tam simmetrik permeabilite tenzoru — birinci-dərəcəli fiziki xassə
    (Phase 2: "Full Tensor Permeability Implementation"), gələcək MPFA-O
    üçün HAZIRLIQ (bax `imex2d/simulation/discretization.py` modul
    docstring-i).

        K = [[Kxx, Kxy, Kxz],
             [Kxy, Kyy, Kyz],
             [Kxz, Kyz, Kzz]]     (simmetrik fərz edilir — Kyx=Kxy, Kzx=Kxz, Kzy=Kyz)

    YALNIZ 6 MÜSTƏQİL komponent (Kxx,Kyy,Kzz,Kxy,Kxz,Kyz) SAXLANILIR —
    simmetrik cütlər (Kyx və s.) AYRICA sahə kimi TƏKRARLANMIR, `as_
    matrices()` onları hər çağırışda RİYAZİ olaraq bərpa edir. Hər
    komponent `PropertyMap`-dir, ona görə hüceyrə-hüceyrə DƏYİŞƏ bilər
    (heterojen rezervuar, bax audit §8).

    Öz-qiymətlər (`eigenvalues`), müsbət-müəyyənlik (`validate`),
    fırlanma (`rotate`) və vahid çevirməsi (`convert_units`) təmin
    edilir — hamısı VEKTORLAŞDIRILIB (`np.linalg.eigvalsh` bütün
    hüceyrələr üçün TƏK çağırış, Python dövrü YOXDUR).

    BU FAZADA DA HEÇ BİR HƏLLEDİCİ (TPFA) bunu İSTİFADƏ ETMİR —
    `RockProperties.permx/permy/permz` (diaqonal) YEGANƏ TPFA-nın
    oxuduğu mənbədir, DƏYİŞMİR. Bu sinif YALNIZ off-diaqonal
    (Kxy/Kxz/Kyz) anizotropluq məlumatını İTİRMƏDƏN daşımaq/doğrulamaq/
    saxlamaq üçündür ki, gələcək MPFA-O onu birbaşa istifadə edə bilsin
    — TPFA-nın bunu SƏSSİZCƏ diaqonala "yumşaltması" (scalarize)
    QADAĞANDIR (bax `has_off_diagonal`/`TwoPointFluxDiscretization.build`
    xəbərdarlığı). DOĞRU İFADƏ: "tam tenzor permeabilite TƏMSİL OLUNUR
    VƏ DOĞRULANIR, amma tam tenzor axın diskretizasiyası MPFA TƏLƏB
    EDİR" — "simulyator tərəfindən tam dəstəklənir" DEYİL.
    """
    kxx: PropertyMap
    kyy: PropertyMap
    kzz: PropertyMap
    kxy: Optional[PropertyMap] = None
    kxz: Optional[PropertyMap] = None
    kyz: Optional[PropertyMap] = None

    def __post_init__(self):
        n = self.kxx.values.size
        for label, component in (("Kyy", self.kyy), ("Kzz", self.kzz),
                                 ("Kxy", self.kxy), ("Kxz", self.kxz), ("Kyz", self.kyz)):
            if component is not None and component.values.size != n:
                raise ValueError(
                    f"{label}: hüceyrə sayı Kxx ilə uyğun gəlmir "
                    f"({component.values.size} != {n})")

    @property
    def ncell(self) -> int:
        return int(self.kxx.values.size)

    def has_off_diagonal(self, tol: float = 1e-12) -> bool:
        """TPFA-nın DÜZGÜN HƏLL EDƏ BİLMƏDİYİ off-diaqonal komponent varmı
        (yəni real, sıfırdan fərqli anizotropluq bucağı)."""
        for component in (self.kxy, self.kxz, self.kyz):
            if component is not None and np.any(np.abs(component.values) > tol):
                return True
        return False

    def _off_diagonal_values(self, component: Optional[PropertyMap]) -> np.ndarray:
        return component.values if component is not None else np.zeros(self.ncell)

    def as_matrices(self) -> np.ndarray:
        """(ncell, 3, 3) — hər hüceyrə üçün tam simmetrik matris.

        Yalnız `Kxx,Kyy,Kzz,Kxy,Kxz,Kyz` (6 komponent) SAXLANILIR (bax
        audit §3/§5: "avoid storing redundant symmetric components") —
        `Kyx=Kxy` və s. HƏR ÇAĞIRIŞDA bu metodla RİYAZİ olaraq bərpa
        olunur, ayrıca sahə kimi YOX.
        """
        n = self.ncell
        kxy = self._off_diagonal_values(self.kxy)
        kxz = self._off_diagonal_values(self.kxz)
        kyz = self._off_diagonal_values(self.kyz)
        matrix = np.zeros((n, 3, 3))
        matrix[:, 0, 0] = self.kxx.values
        matrix[:, 1, 1] = self.kyy.values
        matrix[:, 2, 2] = self.kzz.values
        matrix[:, 0, 1] = matrix[:, 1, 0] = kxy
        matrix[:, 0, 2] = matrix[:, 2, 0] = kxz
        matrix[:, 1, 2] = matrix[:, 2, 1] = kyz
        return matrix

    # ── invariantlar (bax audit §10) — HEÇ BİRİ tenzoru DƏYİŞMİR ─────────
    def eigenvalues(self) -> np.ndarray:
        """(ncell, 3) — artan sırada (`np.linalg.eigvalsh` konvensiyası).

        Vektorlaşdırılıb — `np.linalg.eigvalsh` bütün hüceyrələr üçün TƏK
        LAPACK çağırışında işləyir (bax audit §17: "validation should
        scale approximately O(N)... avoid Python loops")."""
        return np.linalg.eigvalsh(self.as_matrices())

    def min_eigenvalue(self) -> np.ndarray:
        return self.eigenvalues()[:, 0]

    def max_eigenvalue(self) -> np.ndarray:
        return self.eigenvalues()[:, -1]

    def determinant(self) -> np.ndarray:
        return np.linalg.det(self.as_matrices())

    def trace(self) -> np.ndarray:
        return self.kxx.values + self.kyy.values + self.kzz.values

    def anisotropy_ratio(self) -> np.ndarray:
        """λ_max / λ_min (hər hüceyrə üzrə) — Kmin<=0 olan hüceyrələr
        üçün `inf` (belə hüceyrə artıq `validate()`-də XƏTA kimi
        tutulur, bu, YALNIZ bu metodun sıfıra bölünmə ilə ÇÖKMƏMƏSİ
        üçündür)."""
        eig = self.eigenvalues()
        lo, hi = eig[:, 0], eig[:, -1]
        return np.divide(hi, lo, out=np.full_like(hi, np.inf), where=lo > 0.0)

    # ── validasiya (bax audit §4/§6) ──────────────────────────────────────
    def validate(self, label: str = "Permeabilite tenzoru") -> ValidationResult:
        """Simmetriya KONSTRUKSİYACA təmin olunur (bax `as_matrices` —
        `Kyx`/`Kzx`/`Kzy` AYRICA SAXLANILMIR, ona görə "simmetriya
        pozulması" strukturca MÜMKÜN DEYİL, yoxlanmasına EHTİYAC yoxdur).

        Bu metod YALNIZ FİZİKİ etibarlılığı yoxlayır:
          1. NaN/sonsuz komponent (hər hansı).
          2. Müsbət-müəyyənlik: λ_min(K) > 0 — TƏKCƏ diaqonalın müsbət
             olması KİFAYƏT DEYİL (bax audit §4: "a matrix can have
             positive diagonal entries while still being non-positive-
             definite") — ona görə HƏQİQİ məxsusi qiymət hesablanır,
             `Kxx>0 and Kyy>0 and Kzz>0` kimi SƏTHİ yoxlama İŞLƏDİLMİR.

        Heç bir dəyər DÜZƏLDİLMİR/"təmir edilmir" (bax audit §6: "do not
        repair the tensor silently") — yalnız AÇIQ xəta bildirilir.
        """
        result = ValidationResult()
        matrices = self.as_matrices()
        finite = np.isfinite(matrices).all(axis=(1, 2))
        if not np.all(finite):
            n_bad = int(np.sum(~finite))
            result.errors.append(
                f"{label}: {n_bad} hüceyrədə NaN/sonsuz komponent var — "
                "eigenvalue hesablamaq mümkün deyil.")
            return result   # NaN-lı matrisdə eigvalsh ƏLAVƏ NaN yayar, davam etmə

        eig = np.linalg.eigvalsh(matrices)
        min_eig = eig[:, 0]
        invalid = min_eig <= 0.0
        if np.any(invalid):
            n_bad = int(np.sum(invalid))
            worst = float(min_eig.min())
            result.errors.append(
                f"{label}: {n_bad} hüceyrə müsbət-müəyyən DEYİL "
                f"(minimum məxsusi qiymət {worst:.6g} <= 0) — fiziki cəhətdən "
                "etibarsız permeabilite tenzoru (λ_min(K) > 0 tələb olunur).")
        return result

    # ── fırlanma (bax audit §9) ────────────────────────────────────────────
    def rotate(self, rotation: np.ndarray) -> "PermeabilityTensor":
        """K_rotated = R · K · Rᵀ — bütün hüceyrələrə EYNİ (qlobal)
        fırlanma tətbiq edir.

        `rotation` (3,3) ORTOQONAL olmalıdır (R·Rᵀ=I) — yoxlanılır, əks
        halda `ValueError`: qeyri-ortoqonal "fırlanma" məxsusi
        qiymətləri DƏYİŞDİRƏR (bax audit: "must preserve positive
        definiteness under orthogonal rotation" — bu YALNIZ HƏQİQİ
        ortoqonal matrislər üçün riyazi olaraq DOĞRUDUR).
        """
        rotation = np.asarray(rotation, float)
        if rotation.shape != (3, 3):
            raise ValueError(f"rotation (3,3) olmalıdır, alındı {rotation.shape}")
        if not np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-8):
            raise ValueError(
                "rotation ORTOQONAL olmalıdır (R @ R.T == I) — əks halda fırlanma "
                "fiziki cəhətdən düzgün deyil və məxsusi qiymətləri (deməli "
                "müsbət-müəyyənliyi) poza bilər.")
        matrices = self.as_matrices()
        rotated = np.einsum("ij,njk,lk->nil", rotation, matrices, rotation)
        return PermeabilityTensor(
            kxx=PropertyMap("KXX", rotated[:, 0, 0].copy()),
            kyy=PropertyMap("KYY", rotated[:, 1, 1].copy()),
            kzz=PropertyMap("KZZ", rotated[:, 2, 2].copy()),
            kxy=PropertyMap("KXY", rotated[:, 0, 1].copy()),
            kxz=PropertyMap("KXZ", rotated[:, 0, 2].copy()),
            kyz=PropertyMap("KYZ", rotated[:, 1, 2].copy()))

    # ── vahid çevirməsi (bax audit §7) ─────────────────────────────────────
    def convert_units(self, from_unit: str, to_unit: str) -> "PermeabilityTensor":
        """BÜTÜN 6 komponentə EYNİ çevirməni tətbiq edir (bax audit §7:
        "must apply to every component consistently... Do not convert
        only the diagonal"). Mövcud `unit_conversions.convert`-i işlədir
        — YENİ çevirmə düsturu İCAD EDİLMİR."""
        def _convert(component: Optional[PropertyMap]) -> Optional[PropertyMap]:
            if component is None:
                return None
            return PropertyMap(component.name,
                               convert(component.values, from_unit, to_unit, "permeability"),
                               unit=to_unit)
        return PermeabilityTensor(
            kxx=_convert(self.kxx), kyy=_convert(self.kyy), kzz=_convert(self.kzz),
            kxy=_convert(self.kxy), kxz=_convert(self.kxz), kyz=_convert(self.kyz))


@dataclass
class RockProperties:
    """Statik süxur xassələri — geoloji modeldən gəlir."""
    porosity: PropertyMap
    permx: PropertyMap
    permy: PropertyMap
    permz: Optional[PropertyMap] = None
    net_to_gross: Optional[PropertyMap] = None
    compressibility: float = 4.5e-5
    #: HƏLƏ HEÇ BİR HƏLLEDİCİ TƏRƏFİNDƏN İSTİFADƏ OLUNMUR (bax
    #: `PermeabilityTensor` docstring-i) — yalnız gələcək MPFA-O üçün
    #: opt-in verilənlər daşıyıcısı. `None` (defolt) — mövcud bütün
    #: modellər ÜÇÜN DAVRANIŞ TAM EYNİDİR.
    permeability_tensor: Optional[PermeabilityTensor] = None

    def validate(self) -> list:
        """Sərt fiziki xətalar. `validate_warnings()` — qeyri-adi (amma
        mümkün) diapazon xəbərdarlıqları üçün, bax Phase 1 hesabatı.

        HƏQİQİ SƏHV (bax tapşırıq: NaN/inf idarəetməsi auditi): əvvəllər
        bu metod `self.porosity.values <= 0` kimi XAM müqayisələr
        işlədirdi — NaN dəyər üçün `NaN <= 0` HƏMİŞƏ `False` qaytarır,
        ona görə NaN/sonsuz PORO/PERMX/PERMY SƏSSİZCƏ bu yoxlamadan
        keçib `ReservoirModel.validate()` (bax `reservoir_model.py`)
        vasitəsilə simulyasiyaya buraxıla bilərdi. İndi eyni fayldan
        onsuz da idxal edilən (`validate_warnings()`-in artıq işlətdiyi)
        `validate_porosity`/`validate_permeability` istifadə olunur —
        bunlar NaN/sonsuzu AÇIQ xəta kimi tuturlar (bax `validation.py`
        `_finite_issue`)."""
        issues = []
        issues += validate_porosity(self.porosity.values, "PORO").errors
        issues += validate_permeability(self.permx.values, "PERMX").errors
        issues += validate_permeability(self.permy.values, "PERMY").errors
        if self.permeability_tensor is not None:
            issues += self.permeability_tensor.validate().errors
        return issues

    def validate_warnings(self) -> list:
        """Rədd edilməyən, amma qeyri-adi diapazon xəbərdarlıqları."""
        warnings = validate_porosity(self.porosity.values, "PORO").warnings
        warnings += validate_permeability(self.permx.values, "PERMX").warnings
        warnings += validate_permeability(self.permy.values, "PERMY").warnings
        if self.permz is not None:
            warnings += validate_permeability(self.permz.values, "PERMZ").warnings
        return warnings


@dataclass
class FluidProperties:
    """PLACEHOLDER — sabit flüid xassələri.

    Bu sinif müvəqqətidir. Real black-oil davranışı IPVTProvider vasitəsilə
    gələcək; PVT provider inject edilmədikdə mühərrik bu dəyərləri oxuyur.
    """
    water_viscosity: float = 0.5
    oil_viscosity: float = 3.0
    water_fvf: float = 1.0
    oil_fvf: float = 1.15
    water_compressibility: float = 4.4e-5
    oil_compressibility: float = 1.4e-4
    water_density: float = 1010.0   # səth sıxlığı, kg/m3
    oil_density: float = 850.0      # səth sıxlığı, kg/m3

    def validate(self) -> list:
        """Sərt fiziki xətalar (əvvəllər bu sinifdə HEÇ BİR yoxlama yox
        idi — bax Phase 1 audit)."""
        issues = []
        issues += validate_viscosity(self.water_viscosity, "su lözlüyü").errors
        issues += validate_viscosity(self.oil_viscosity, "neft lözlüyü").errors
        issues += validate_density(self.water_density, "su sıxlığı").errors
        issues += validate_density(self.oil_density, "neft sıxlığı").errors
        issues += validate_compressibility(self.water_compressibility, "su sıxılması").errors
        issues += validate_compressibility(self.oil_compressibility, "neft sıxılması").errors
        if self.water_fvf <= 0 or self.oil_fvf <= 0:
            issues.append("Formasiya həcm əmsalı (Bw/Bo) müsbət olmalıdır.")
        return issues

    def validate_warnings(self) -> list:
        warnings = []
        warnings += validate_viscosity(self.water_viscosity, "su lözlüyü").warnings
        warnings += validate_viscosity(self.oil_viscosity, "neft lözlüyü").warnings
        warnings += validate_density(self.water_density, "su sıxlığı").warnings
        warnings += validate_density(self.oil_density, "neft sıxlığı").warnings
        warnings += validate_compressibility(self.water_compressibility, "su sıxılması").warnings
        warnings += validate_compressibility(self.oil_compressibility, "neft sıxılması").warnings
        return warnings
