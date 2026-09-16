"""Eclipse deck-dən PVT cədvəllərinin oxunması — G1 (SPE1 üçün).

Üç açar söz oxunur (deck-in `PROPS` bölməsi):

    PVTW   — su: istinad təzyiqi, Bw, sıxılma, μw, özlülük əmsalı (bir sətir)
    PVDG   — quru qaz: p, Bg, μg
    PVTO   — neft: hər Rs üçün DOYMUŞ sətir + ondan yuxarı DOYMAMIŞ sətirlər

    PVTO
    --  Rs      Pb       Bo       μo
        0.165   400.0    1.012    1.17
                800.0    1.008    1.20      ← doymamış davam (eyni Rs)
    /
        1.270   4014.7   1.695    0.51
                9014.7   1.579    0.74      ← doymamış davam
    /
    /

NİYƏ AYRI MODUL: `PVTTable` (domain) TƏK təzyiq şəbəkəsidir — bütün
sütunlar eyni `pressure` massivinə bağlıdır. Deck isə üç fərqli şəbəkə
verir (PVDG 14.7…9014.7, PVTO doymuş qolu 14.7…5014.7 psia) və üstəlik
PVTO-da HƏR Rs üçün ayrıca doymamış qol var. Bu modul deck-i İTKİSİZ
oxuyur (`DeckPvt`), birləşdirmə isə AYRI addımdır (`DeckPvt.to_pvt_table`).

⏳ MƏHDUDİYYƏT (açıq yazılır, gizlədilmir): `PVTTable` bir doymamış qol
daşıya bilir. Deck-də bir neçə Rs qolu olanda `to_pvt_table` yalnız SEÇİLƏN
qolu cədvələ köçürür; qalanları `DeckPvt`-də QALIR və μo(p, Rs) / Bo(p, Rs)
işi (SPE1 boşluqları G2/G3) onları oradan oxuyacaq.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from ..domain.pvt import PVTTable
from ..domain.unit_conversions import convert, to_engine_units
from ..logging_setup import get_logger

LOG = get_logger(__name__)


class PvtFormatError(Exception):
    """Deck-in PVT bölməsi gözləniləndən fərqlidir."""


# ═══════════════════════════════ deck mətninin ayrılması ══════════════

def _keyword_body(text: str, keyword: str) -> str:
    """Açar sözdən sonrakı hissə — növbəti açar sözə qədər.

    `scal_io.read_swof` ilə EYNİ qayda: açar söz öz sətrində durur,
    növbəti BÖYÜK HƏRFLİ açar söz bölməni bitirir.
    """
    match = re.search(rf"^\s*{keyword}\s*$", text, re.MULTILINE | re.IGNORECASE)
    if match is None:
        raise PvtFormatError(f"Faylda {keyword} açar sözü tapılmadı.")
    body = text[match.end():]
    stop = re.search(r"^\s*[A-Z][A-Z0-9_]{2,}\s*$", body, re.MULTILINE)
    return body if stop is None else body[:stop.start()]


def _numbers(line: str) -> List[float]:
    return [float(token) for token in
            re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)]


def _clean(raw_line: str) -> str:
    return re.sub(r"--.*", "", raw_line).strip()


def _read_text(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        return handle.read()


# ═══════════════════════════════ PVTW ════════════════════════════════

@dataclass
class WaterPvt:
    """`PVTW` — suyun xassələri istinad təzyiqinə görə (mühərrik vahidləri).

    Eclipse düsturu (X = c_w·(p − p_ref)):

        Bw(p) = Bw_ref / (1 + X + X²/2)
        μw(p) = μw_ref / (1 + Y + Y²/2),   Y = −c_μ·(p − p_ref)

    SPE1-də `c_μ = 0`, yəni μw sabitdir.
    """

    reference_pressure: float          # bar
    formation_volume_factor: float     # ölçüsüz
    compressibility: float             # 1/bar
    viscosity: float                   # cP
    viscosibility: float = 0.0         # 1/bar

    def water_fvf(self, pressure) -> np.ndarray:
        x = self.compressibility * (np.asarray(pressure, float)
                                    - self.reference_pressure)
        return self.formation_volume_factor / (1.0 + x + 0.5 * x * x)

    def water_viscosity(self, pressure) -> np.ndarray:
        pressure = np.asarray(pressure, float)
        if self.viscosibility == 0.0:
            return np.full(pressure.shape, float(self.viscosity))
        y = -self.viscosibility * (pressure - self.reference_pressure)
        return self.viscosity / (1.0 + y + 0.5 * y * y)


def read_pvtw(path: str, pressure_unit: str = "psi",
              viscosity_unit: str = "cP") -> WaterPvt:
    """`PVTW` sətrini oxuyur (FIELD deck-i üçün defolt vahidlər)."""
    for raw_line in _keyword_body(_read_text(path), "PVTW").splitlines():
        line = _clean(raw_line)
        if not line or line.startswith("/"):
            continue
        values = _numbers(line.rstrip("/"))
        if len(values) < 4:
            continue
        reference = to_engine_units(values[0], pressure_unit, "pressure")
        return WaterPvt(
            reference_pressure=reference,
            formation_volume_factor=values[1],
            # sıxılma "1/[təzyiq]" formasındadır — bax `to_engine_units`
            compressibility=to_engine_units(values[2], pressure_unit,
                                            "compressibility"),
            viscosity=to_engine_units(values[3], viscosity_unit, "viscosity"),
            viscosibility=(to_engine_units(values[4], pressure_unit,
                                           "compressibility")
                           if len(values) > 4 else 0.0))
    raise PvtFormatError("PVTW bölməsində rəqəm tapılmadı.")


# ═══════════════════════════════ PVDG ════════════════════════════════

@dataclass
class GasPvt:
    """`PVDG` — quru qaz (mühərrik vahidləri: bar, m³/sm³, cP)."""

    pressure: np.ndarray
    formation_volume_factor: np.ndarray
    viscosity: np.ndarray

    def __post_init__(self):
        for name in ("pressure", "formation_volume_factor", "viscosity"):
            setattr(self, name, np.asarray(getattr(self, name), float).ravel())

    def validate(self) -> List[str]:
        issues = []
        if self.pressure.size < 2:
            issues.append("PVDG: ən azı iki sətir lazımdır.")
            return issues
        if not (self.pressure.size == self.formation_volume_factor.size
                == self.viscosity.size):
            issues.append("PVDG: sütun uzunluqları fərqlidir.")
            return issues
        if np.any(~np.isfinite(self.pressure)) or np.any(~np.isfinite(
                self.formation_volume_factor)) or np.any(~np.isfinite(self.viscosity)):
            issues.append("PVDG: NaN/sonsuz dəyər var.")
        if np.any(np.diff(self.pressure) <= 0):
            issues.append("PVDG: təzyiq artan sıralı olmalıdır.")
        if np.any(np.diff(self.formation_volume_factor) > 1e-12):
            issues.append("PVDG: Bg təzyiqlə artır — azalmalıdır.")
        return issues


def read_pvdg(path: str, pressure_unit: str = "psi",
              gas_fvf_unit: str = "rb/Mscf",
              viscosity_unit: str = "cP") -> GasPvt:
    """`PVDG` cədvəlini oxuyur."""
    pressure, bg, mu = [], [], []
    for raw_line in _keyword_body(_read_text(path), "PVDG").splitlines():
        line = _clean(raw_line)
        if not line or line.startswith("/"):
            continue
        values = _numbers(line.rstrip("/"))
        if len(values) >= 3:
            pressure.append(values[0])
            bg.append(values[1])
            mu.append(values[2])
    if not pressure:
        raise PvtFormatError("PVDG bölməsində rəqəm tapılmadı.")

    table = GasPvt(
        pressure=to_engine_units(np.asarray(pressure, float), pressure_unit,
                                 "pressure"),
        formation_volume_factor=convert(np.asarray(bg, float), gas_fvf_unit,
                                        "m3/sm3", "gas_fvf"),
        viscosity=to_engine_units(np.asarray(mu, float), viscosity_unit,
                                  "viscosity"))
    issues = table.validate()
    if issues:
        raise PvtFormatError("; ".join(issues))
    return table


# ═══════════════════════════════ PVTO ════════════════════════════════

@dataclass
class OilBranch:
    """Bir `Rs` üçün doymuş nöqtə + ondan yuxarı doymamış sətirlər.

    `pressure[0]` doyma təzyiqidir (Pb); qalan sətirlər DOYMAMIŞ qoldur
    (eyni Rs, artan təzyiq).
    """

    solution_gor: float                # sm³/sm³
    pressure: np.ndarray               # bar
    formation_volume_factor: np.ndarray
    viscosity: np.ndarray              # cP

    def __post_init__(self):
        for name in ("pressure", "formation_volume_factor", "viscosity"):
            setattr(self, name, np.asarray(getattr(self, name), float).ravel())

    @property
    def bubble_point(self) -> float:
        return float(self.pressure[0])

    @property
    def has_undersaturated(self) -> bool:
        return self.pressure.size > 1


@dataclass
class OilPvt:
    """`PVTO` — Rs qolları (artan Rs sırası ilə)."""

    branches: List[OilBranch] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.branches)

    @property
    def saturated(self):
        """(Rs, Pb, Bo, μo) — hər qolun DOYMUŞ nöqtəsi."""
        rs = np.array([branch.solution_gor for branch in self.branches])
        pressure = np.array([branch.bubble_point for branch in self.branches])
        bo = np.array([branch.formation_volume_factor[0] for branch in self.branches])
        mu = np.array([branch.viscosity[0] for branch in self.branches])
        return rs, pressure, bo, mu

    def branch_for(self, solution_gor: float) -> OilBranch:
        """Verilmiş Rs-ə ƏN YAXIN qol (dəqiq uyğunluq tələb olunmur)."""
        if not self.branches:
            raise PvtFormatError("PVTO-da heç bir Rs qolu yoxdur.")
        rs = np.array([branch.solution_gor for branch in self.branches])
        return self.branches[int(np.argmin(np.abs(rs - float(solution_gor))))]

    def validate(self) -> List[str]:
        issues = []
        if len(self.branches) < 2:
            issues.append("PVTO: ən azı iki Rs qolu lazımdır.")
            return issues
        rs, pressure, bo, _ = self.saturated
        if np.any(np.diff(rs) <= 0):
            issues.append("PVTO: Rs artan sıralı olmalıdır.")
        if np.any(np.diff(pressure) <= 0):
            issues.append("PVTO: doyma təzyiqi Rs ilə artmalıdır.")
        for branch in self.branches:
            if np.any(np.diff(branch.pressure) <= 0):
                issues.append(f"PVTO (Rs = {branch.solution_gor:.4g}): "
                              f"təzyiq artan sıralı olmalıdır.")
            if np.any(~np.isfinite(branch.formation_volume_factor)):
                issues.append(f"PVTO (Rs = {branch.solution_gor:.4g}): "
                              f"Bo sütununda NaN/sonsuz var.")
        return issues


def read_pvto(path: str, pressure_unit: str = "psi",
              rs_unit: str = "Mscf/stb",
              viscosity_unit: str = "cP") -> OilPvt:
    """`PVTO` cədvəlini oxuyur — hər `/` bir Rs qolunu bitirir."""
    branches: List[OilBranch] = []
    current_rs: Optional[float] = None
    pressure: List[float] = []
    bo: List[float] = []
    mu: List[float] = []

    def close():
        nonlocal current_rs, pressure, bo, mu
        if current_rs is not None and pressure:
            branches.append(OilBranch(
                solution_gor=convert(current_rs, rs_unit, "sm3/sm3",
                                     "solution_gor"),
                pressure=to_engine_units(np.asarray(pressure, float),
                                         pressure_unit, "pressure"),
                formation_volume_factor=np.asarray(bo, float),
                viscosity=to_engine_units(np.asarray(mu, float),
                                          viscosity_unit, "viscosity")))
        current_rs, pressure, bo, mu = None, [], [], []

    for raw_line in _keyword_body(_read_text(path), "PVTO").splitlines():
        line = _clean(raw_line)
        if not line:
            continue
        if line.startswith("/"):
            close()
            continue
        terminated = line.endswith("/")
        values = _numbers(line.rstrip("/"))
        if len(values) >= 4:
            # Rs ilə başlayan sətir = YENİ qolun doymuş nöqtəsi
            close()
            current_rs = values[0]
            pressure, bo, mu = [values[1]], [values[2]], [values[3]]
        elif len(values) >= 3 and current_rs is not None:
            # davam sətri = eyni Rs-in doymamış nöqtəsi
            pressure.append(values[0])
            bo.append(values[1])
            mu.append(values[2])
        if terminated:
            close()
    close()

    table = OilPvt(branches)
    issues = table.validate()
    if issues:
        raise PvtFormatError("; ".join(issues))
    return table


# ═══════════════════════════ hamısı birlikdə ═════════════════════════

@dataclass
class DeckPvt:
    """Deck-dən oxunan BÜTÜN PVT məlumatı (itkisiz, mühərrik vahidlərində)."""

    water: WaterPvt
    gas: GasPvt
    oil: OilPvt

    def to_pvt_table(self, reference_rs: Optional[float] = None,
                     source: str = "deck") -> PVTTable:
        """Deck-i mühərrikin TƏK ŞƏBƏKƏLİ `PVTTable`-ına köçürür.

        Şəbəkə = PVTO-nun doymuş təzyiqləri ∪ PVDG təzyiqləri.

        `reference_rs` — hansı doymamış qolun cədvələ köçürüləcəyi
        (verilməsə ƏN BÖYÜK Rs qolu). Doyma təzyiqindən YUXARI sətirlər
        həmin qoldan gəlir, çünki `PVTTable` yalnız BİR doymamış qol
        daşıya bilir.

        ⏳ Bu, TƏQRİBDİR: deck-də hər Rs-in öz doymamış qolu var və
        `DeckPvt`-də hamısı SAXLANILIR. Həqiqi `Bo(p, Rs)` / `μo(p, Rs)`
        SPE1 boşluqları G2/G3-ün işidir.
        """
        rs_nodes, pb_nodes, bo_nodes, mu_nodes = self.oil.saturated
        branch = (self.oil.branches[-1] if reference_rs is None
                  else self.oil.branch_for(reference_rs))

        grid = np.unique(np.concatenate([pb_nodes, self.gas.pressure,
                                         branch.pressure]))

        saturated_mask = grid <= branch.bubble_point + 1e-12
        oil_fvf = np.empty_like(grid)
        oil_viscosity = np.empty_like(grid)
        solution_gor = np.empty_like(grid)

        # doymuş hissə — PVTO-nun qol başlarından
        oil_fvf[saturated_mask] = np.interp(grid[saturated_mask], pb_nodes, bo_nodes)
        oil_viscosity[saturated_mask] = np.interp(grid[saturated_mask], pb_nodes,
                                                  mu_nodes)
        solution_gor[saturated_mask] = np.interp(grid[saturated_mask], pb_nodes,
                                                 rs_nodes)

        # doymamış hissə — SEÇİLMİŞ qoldan; Rs plato qalır
        above = ~saturated_mask
        if np.any(above):
            if branch.has_undersaturated:
                oil_fvf[above] = np.interp(grid[above], branch.pressure,
                                           branch.formation_volume_factor)
                oil_viscosity[above] = np.interp(grid[above], branch.pressure,
                                                 branch.viscosity)
            else:
                LOG.warning("PVTO: seçilmiş Rs = %.4g qolunda doymamış sətir "
                            "yoxdur — Bo/μo doyma nöqtəsində saxlanılır.",
                            branch.solution_gor)
                oil_fvf[above] = branch.formation_volume_factor[0]
                oil_viscosity[above] = branch.viscosity[0]
            solution_gor[above] = branch.solution_gor

        # ÖLÇÜLMÜŞ TƏHLÜKƏ (G1): `BlackOilPVTProvider` doymamış sıxılmanı
        # (c_o) cədvəlin Pb-dən YUXARI düyünlərindən fit edir və ƏN AZI
        # İKİ düyün tələb edir. Deck-də hər Rs qolunda çox vaxt CƏMİ BİR
        # doymamış sətir olur (SPE1 belədir) — o zaman provider səssizcə
        # ehtiyat qiymətə düşür. Ölçüldü: həqiqi c_o = 2.06e-4 1/bar,
        # ehtiyat qiymət 3.12e-3 1/bar, yəni 15 DƏFƏ böyük.
        if int(np.count_nonzero(above)) < 2:
            LOG.warning(
                "PVTO: seçilmiş Rs = %.4g qolu cədvələ doyma təzyiqindən "
                "yuxarı yalnız %d düyün verir — doymamış sıxılma (c_o) "
                "etibarlı fit oluna bilməz və provider ehtiyat qiymətə "
                "düşəcək. Başqa `reference_rs` seçin və ya deck-in öz "
                "qolundan hesablayın (SPE1 boşluqları G2/G3).",
                branch.solution_gor, int(np.count_nonzero(above)))

        dropped = sum(1 for other in self.oil.branches
                      if other is not branch and other.has_undersaturated)
        if dropped:
            LOG.info("PVTO: %d doymamış qol cədvələ köçürülmədi (cədvəl BİR qol "
                     "daşıyır) — məlumat `DeckPvt`-də saxlanılır, bax G2/G3.",
                     dropped)

        return PVTTable(
            pressure=grid,
            oil_fvf=oil_fvf,
            oil_viscosity=oil_viscosity,
            solution_gor=solution_gor,
            water_fvf=self.water.water_fvf(grid),
            water_viscosity=self.water.water_viscosity(grid),
            bubble_point=branch.bubble_point,
            gas_fvf=np.interp(grid, self.gas.pressure,
                              self.gas.formation_volume_factor),
            gas_viscosity=np.interp(grid, self.gas.pressure, self.gas.viscosity),
            source=source)


def read_deck_pvt(path: str, **kwargs) -> DeckPvt:
    """`PVTW` + `PVDG` + `PVTO` — üçünü birlikdə oxuyur."""
    return DeckPvt(water=read_pvtw(path, **{k: v for k, v in kwargs.items()
                                            if k in ("pressure_unit",
                                                     "viscosity_unit")}),
                   gas=read_pvdg(path, **{k: v for k, v in kwargs.items()
                                          if k in ("pressure_unit", "gas_fvf_unit",
                                                   "viscosity_unit")}),
                   oil=read_pvto(path, **{k: v for k, v in kwargs.items()
                                          if k in ("pressure_unit", "rs_unit",
                                                   "viscosity_unit")}))
