"""Həqiqi GEOMETRİK ANİZOTROPLUQ — VAHİD (tək) transformasiya yolu (A4).

Bu modul layihədəki YEGANƏ anizotrop-məsafə mənbəyidir. Kriging matrisi
(`geology/interpolation.py`), qonşuluq axtarışı (`geology/spatial_
search.py`), variogram qiymətləndirməsi (`geology/variogram.py`) və
SGS/SIS (`geology/sgs.py`, `geology/facies.py`) HAMISI eyni
`AnisotropyParams.transform()`-dan keçir — məhz buna görə "A2 bir
həndəsə, A1 başqa həndəsə işlədir" uyğunsuzluğu MÜMKÜN DEYİL.

Riyaziyyat
----------
Xam koordinat ``x = [X, Y, Z]``. Əvvəlcə DÖNMƏ (rotasiya), sonra
MİQYASLANMA::

    u  = R₀ x          (azimut: X,Y → major, minor)
    u' = R_dip u       (dip: major, vertical müstəvisində dönmə)
    x''= S u'          (S = diag(1, a_maj/a_min, a_maj/a_v))

Nəticədə ``x''`` fəzasında ADİ Evklid məsafəsi birbaşa `range_major`
ilə müqayisə edilə bilər::

    d_ani(p, q) = ‖x''(p) − x''(q)‖₂

və variogram `γ(d_ani)`-ni major radiusla qiymətləndirir. Yəni::

    d_ani² = (u'_maj)² + (u'_min·a_maj/a_min)² + (u'_v·a_maj/a_v)²

bu, tapşırıqdakı ``sqrt((x'/a_x)² + (y'/a_y)² + (z'/a_z)²)`` düsturunun
`a_maj` ilə vurulmuş (ölçüsü saxlanılmış) formasıdır — nisbətlər
eynidir, yalnız məsafə `range_major` vahidində DEYİL, uzunluq vahidində
qalır ki, `range_`-i olduğu kimi variograma ötürmək mümkün olsun.

Konvensiyalar
-------------
* `azimuth_deg` — major oxun istiqaməti, +Y (Şimal) oxundan SAAT ƏQRƏBİ
  istiqamətində (0°=Şimal, 90°=Şərq). Geoloji proqramlarda adət olunan.
* `dip_deg` — major oxun ÜFÜQİ müstəvidən çıxma bucağı; müsbət dip
  major oxu azimut istiqamətində Z-nin ARTAN tərəfinə (dərinliyə) əyir.
  `dip_deg == 0` (defolt) olanda `transform()` ƏDƏD-ƏDƏD əvvəlki (M2)
  davranışı təkrarlayır — qısa qapanma ilə dönmə heç tətbiq edilmir.
* Z oxu layihənin qalanında DƏRİNLİKDİR (aşağı müsbət, bax
  `domain/geometry.cell_depths`) — transform bunu bilmir, sadəcə üçüncü
  oxu miqyaslayır, ona görə hər iki konvensiya ilə düzgün işləyir.

Tenzor-uyğunluq (A4.5)
----------------------
`matrix()` və `metric_tensor()` eyni həndəsəni XƏTTİ CƏBR formasında
verir: ``d² = Δxᵀ G Δx``, ``G = MᵀM``. Gələcək tam keçiricilik-tenzoru
işi (Part D) üçün genişləndirmə nöqtəsi budur — burada NƏ tam tenzor
interpolyasiyası, NƏ də səhv bir tenzor modeli implementasiya OLUNMUR,
yalnız mövcud həndəsə tenzor dilində AÇIQ şəkildə ifadə edilir.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

#: Sıfıra bölünmənin qarşısını alan minimal radius. Radiuslar uzunluq
#: vahidindədir (m/ft), ona görə 1e-9 istənilən real şəbəkədən kiçikdir.
_MIN_RANGE = 1e-9


class AnisotropyError(ValueError):
    """Etibarsız anizotropluq parametrləri — SƏSSİZ düzəliş EDİLMİR."""


@dataclass
class AnisotropyParams:
    """Geometrik anizotropluq transformu (bax modul docstring-i).

    `range_major` ≥ `range_minor` şərtini bu sinif MƏCBUR ETMİR (kiçik
    quyu çoxluğunda fit tərsini verə bilər) — `validate()` bunu
    XƏBƏRDARLIQ kimi qaytarır, çünki riyazi cəhətdən nəticə hələ də
    etibarlıdır (sadəcə "major" adı yanıltıcı olur).
    """

    azimuth_deg: float = 0.0
    range_major: float = 1.0
    range_minor: float = 1.0
    range_vertical: float = 1.0
    dip_deg: float = 0.0

    # ── əsas yol: nöqtə transformu ─────────────────────────────────────
    def transform(self, points_xyz: np.ndarray) -> np.ndarray:
        """Nöqtələri anizotrop fəzadan izotrop (vahid `range_major`
        radiuslu) fəzaya çevirir.

        `range_minor == range_major`, `azimuth_deg == 0`, `dip_deg == 0`
        olanda transform (X,Y)-i sadəcə YERDƏYİŞDİRİR (major=Y, minor=X)
        — Evklid normuna təsir etmir, ona görə əvvəlki "yalnız Z
        miqyaslama" davranışı ilə ƏDƏD-ƏDƏD eynidir (bax
        `tests/test_kriging_3d_anisotropy.py::
        test_default_parameters_reproduce_pre_m2_2d_behaviour`).

        Hesablama qəsdən `matrix()` ilə matris hasili KİMİ APARILMIR:
        ardıcıl (dön → miqyasla) forma `dip_deg == 0` halında əvvəlki
        versiya ilə BİT-BİT eyni nəticə verir, matris hasili isə
        yuvarlaqlaşdırmanı dəyişərdi. `matrix()` eyni həndəsənin cəbri
        ifadəsidir (tenzor genişlənməsi üçün), ədədi ETALON isə BUDUR.
        """
        pts = np.asarray(points_xyz, float)
        if pts.ndim == 1:
            pts = pts.reshape(1, -1)
        if pts.ndim != 2 or pts.shape[1] < 2:
            raise AnisotropyError(
                f"Nöqtələr (n,2) və ya (n,3) olmalıdır, alındı: {pts.shape}")
        x, y = pts[:, 0], pts[:, 1]
        z = pts[:, 2] if pts.shape[1] > 2 else np.zeros(pts.shape[0])

        theta = np.radians(self.azimuth_deg)
        c, s = np.cos(theta), np.sin(theta)
        major = x * s + y * c
        minor = x * c - y * s
        vertical = z

        if self.dip_deg:
            # QISA QAPANMA: dip=0 (defolt) olanda heç bir əməliyyat
            # aparılmır → M2 ədədləri toxunulmaz qalır.
            phi = np.radians(self.dip_deg)
            cd, sd = np.cos(phi), np.sin(phi)
            major, vertical = major * cd + vertical * sd, -major * sd + vertical * cd

        ref = max(self.range_major, _MIN_RANGE)
        return np.column_stack([
            major,
            minor * (ref / max(self.range_minor, _MIN_RANGE)),
            vertical * (ref / max(self.range_vertical, _MIN_RANGE)),
        ])

    # ── törəmə həndəsə ─────────────────────────────────────────────────
    def rotation_matrix(self) -> np.ndarray:
        """(3,3) dönmə matrisi ``R = R_dip · R₀`` — miqyaslanma YOX."""
        theta = np.radians(self.azimuth_deg)
        c, s = np.cos(theta), np.sin(theta)
        r0 = np.array([[s, c, 0.0],
                       [c, -s, 0.0],
                       [0.0, 0.0, 1.0]])
        if not self.dip_deg:
            return r0
        phi = np.radians(self.dip_deg)
        cd, sd = np.cos(phi), np.sin(phi)
        r_dip = np.array([[cd, 0.0, sd],
                          [0.0, 1.0, 0.0],
                          [-sd, 0.0, cd]])
        return r_dip @ r0

    def scale_matrix(self) -> np.ndarray:
        """(3,3) diaqonal miqyaslanma ``S = diag(1, a_maj/a_min, a_maj/a_v)``."""
        ref = max(self.range_major, _MIN_RANGE)
        return np.diag([1.0,
                        ref / max(self.range_minor, _MIN_RANGE),
                        ref / max(self.range_vertical, _MIN_RANGE)])

    def matrix(self) -> np.ndarray:
        """Tam transform matrisi ``M = S · R`` (``x'' = M x``).

        `transform()` ilə RİYAZİ olaraq eynidir; ədədi etalon
        `transform()`-dur (bax onun docstring-i).
        """
        return self.scale_matrix() @ self.rotation_matrix()

    def metric_tensor(self) -> np.ndarray:
        """``G = Mᵀ M`` — ``d_ani² = Δxᵀ G Δx``.

        A4.5 genişlənmə nöqtəsi: diaqonal anizotropluq, dönmüş baş
        oxlar və (gələcəkdə) tam tenzor eyni `G` ilə ifadə olunur.
        """
        m = self.matrix()
        return m.T @ m

    def principal_axes(self) -> List[Tuple[str, float, np.ndarray]]:
        """`[(ad, radius, vahid vektor)]` — major/minor/vertical baş
        oxlar XAM (transformasiya edilməmiş) koordinatlarda."""
        rot = self.rotation_matrix()
        return [
            ("major", float(self.range_major), rot[0].copy()),
            ("minor", float(self.range_minor), rot[1].copy()),
            ("vertical", float(self.range_vertical), rot[2].copy()),
        ]

    # ── məsafə ─────────────────────────────────────────────────────────
    def distance(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """Anizotrop məsafə (uzunluq vahidində, `range_major` ilə
        müqayisə edilə bilən) — `a` (n,2|3), `b` (m,2|3) → (n,m)."""
        ta = self.transform(a)
        tb = self.transform(b)
        diff = ta[:, None, :] - tb[None, :, :]
        return np.sqrt(np.sum(diff * diff, axis=-1))

    def pairwise_distance(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """Sətir-sətir (element-üzrə) anizotrop məsafə — `a`, `b` eyni
        sətir sayında, nəticə (n,)."""
        ta = self.transform(a)
        tb = self.transform(b)
        return np.sqrt(np.sum((ta - tb) ** 2, axis=-1))

    # ── introspeksiya / doğrulama ──────────────────────────────────────
    @property
    def is_isotropic(self) -> bool:
        """Üç radius da (nisbi 1e-12 dəqiqliklə) bərabərdirmi."""
        ref = max(abs(self.range_major), _MIN_RANGE)
        return (abs(self.range_major - self.range_minor) <= 1e-12 * ref
                and abs(self.range_major - self.range_vertical) <= 1e-12 * ref)

    @property
    def horizontal_ratio(self) -> float:
        """`range_minor / range_major` — 1.0 = üfüqi izotrop."""
        return float(max(self.range_minor, _MIN_RANGE) / max(self.range_major, _MIN_RANGE))

    def validate(self) -> List[str]:
        """Etibarsız parametrdə `AnisotropyError` atır, ŞÜBHƏLİ (amma
        riyazi olaraq etibarlı) hal üçün xəbərdarlıq siyahısı qaytarır.

        Səssiz düzəliş EDİLMİR — mənfi/sıfır/qeyri-sonlu radius Kriging
        matrisini pozar, ona görə AÇIQ xəta atılır (A3.7/A4 qaydası).
        """
        for name in ("range_major", "range_minor", "range_vertical"):
            value = float(getattr(self, name))
            if not np.isfinite(value):
                raise AnisotropyError(f"{name} sonlu olmalıdır, alındı: {value!r}")
            if value <= 0.0:
                raise AnisotropyError(f"{name} müsbət olmalıdır, alındı: {value!r}")
        for name in ("azimuth_deg", "dip_deg"):
            value = float(getattr(self, name))
            if not np.isfinite(value):
                raise AnisotropyError(f"{name} sonlu olmalıdır, alındı: {value!r}")
        if abs(float(self.dip_deg)) > 90.0:
            raise AnisotropyError(
                f"dip_deg [-90, 90] aralığında olmalıdır, alındı: {self.dip_deg!r}")

        warnings: List[str] = []
        if self.range_minor > self.range_major * (1.0 + 1e-12):
            warnings.append(
                f"range_minor ({self.range_minor:.4g}) > range_major "
                f"({self.range_major:.4g}) — 'major' adı yanıltıcıdır; həndəsə "
                "riyazi olaraq etibarlıdır, amma azimut minor oxu göstərir.")
        ratio = self.horizontal_ratio
        if ratio < 1e-4 or ratio > 1e4:
            warnings.append(
                f"Üfüqi anizotropluq nisbəti həddindən artıq kəskindir "
                f"(minor/major = {ratio:.3g}) — Kriging matrisi pis şərtlənə bilər.")
        return warnings

    @classmethod
    def from_ranges(cls, range_major: float, range_minor: Optional[float] = None,
                    range_vertical: Optional[float] = None, azimuth_deg: float = 0.0,
                    dip_deg: float = 0.0) -> "AnisotropyParams":
        """Verilməyən radiusları `range_major`-a bərabər sayan konstruktor
        (A4.4: major/minor/vertical + azimut baş istiqamətləri)."""
        major = float(range_major)
        return cls(azimuth_deg=float(azimuth_deg),
                   range_major=major,
                   range_minor=major if range_minor is None else float(range_minor),
                   range_vertical=major if range_vertical is None else float(range_vertical),
                   dip_deg=float(dip_deg))


#: İzotrop (heç bir təsir göstərməyən) transform — `None` yoxlamasını
#: təkrarlamamaq üçün.
ISOTROPIC = AnisotropyParams()


def transform_points(points_xyz: np.ndarray,
                     anisotropy: Optional[AnisotropyParams]) -> np.ndarray:
    """`anisotropy is None` olanda nöqtələri (n,3)-ə PADDİNQ edib olduğu
    kimi qaytarır, əks halda `transform()` tətbiq edir.

    Bütün istifadəçi modullar bu funksiyanı çağırır ki, "None = xam
    koordinat" qaydası BİR yerdə yazılsın."""
    pts = np.asarray(points_xyz, float)
    if pts.ndim == 1:
        pts = pts.reshape(1, -1)
    if anisotropy is not None:
        return anisotropy.transform(pts)
    if pts.shape[1] == 2:
        return np.column_stack([pts, np.zeros(pts.shape[0])])
    return pts
