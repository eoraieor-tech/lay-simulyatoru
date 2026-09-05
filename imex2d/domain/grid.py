"""Grid topologiyası — YALNIZ hüceyrələrin bir-birinə necə bağlandığı.

Ölçülər (həcm, sahə, məsafə) burada deyil, CellGeometry-dədir.
Bu ayrılıq gələcəkdə corner-point grid-ə keçidi mümkün edir: topologiya
qalır, yalnız geometriya dəyişir.

AKTİV/QEYRİ-AKTİV HÜCEYRƏLƏR (ACTNUM) — bax `ActiveMap`
========================================================
Eclipse-in `ACTNUM` açar sözü hansı hüceyrənin simulyasiyada İŞTİRAK
ETDİYİNİ göstərir (0 = qeyri-aktiv). Əvvəllər bu massiv YALNIZ
XƏBƏRDARLIQ üçün oxunurdu — model bütün hüceyrələri aktiv sayırdı, bu isə
üç yerdə SƏSSİZ YANLIŞ nəticə verirdi:

  1. həcm/ehtiyat (PV, OOIP) qeyri-aktiv hüceyrələri də sayırdı;
  2. TPFA aktiv↔qeyri-aktiv üzlər qurur, axın "boşluğa" sızırdı;
  3. qeyri-aktiv hüceyrələr xətti sistemdə naməlum kimi qalırdı.

İndi `ACTNUM` grid topologiyasının BİR HİSSƏSİDİR: `build_connections()`
qeyri-aktiv hüceyrəyə toxunan hər üzü ATIR, `ActiveMap` isə
qlobal↔aktiv indeks çevirməsini (`global_to_active`/`active_to_global`)
təmin edir ki, xətti sistem YALNIZ `n_active` naməlumla qurulsun.

KONVENSİYA: massivlərin (təzyiq, doyumluluq, PORO/PERM, PropertyMap …)
SAXLAMA formatı QLOBAL qalır (`ncell` uzunluqda) — 3D görüntü, hesabat,
serializasiya və faylların hamısı qlobal indeksləmə üzərində qurulub.
REDUKSİYA yalnız XƏTTİ SİSTEM SƏRHƏDİNDƏ baş verir (bax
`impes_engine._prealloc` və `implicit/newton.py::_solve_active`).
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Union, Sequence

import numpy as np


@dataclass(frozen=True)
class ActiveMap:
    """Qlobal ↔ aktiv hüceyrə indeks çevirməsi (ACTNUM).

    `actnum`            — (ncell,) int8, 0/1; qlobal indekslə.
    `global_to_active`  — (ncell,) int64; qeyri-aktiv hüceyrə üçün **−1**.
    `active_to_global`  — (n_active,) int64; hər aktiv hüceyrənin qlobal indeksi.

    İkisi bir-birinin tərsidir:
        `global_to_active[active_to_global[a]] == a`  (hər a üçün)
        `active_to_global[global_to_active[g]] == g`  (aktiv g üçün)
    """
    actnum: np.ndarray
    global_to_active: np.ndarray
    active_to_global: np.ndarray

    @classmethod
    def all_active(cls, ncell: int) -> "ActiveMap":
        """ACTNUM verilməyəndə — HƏR hüceyrə aktivdir (köhnə davranış)."""
        ncell = int(ncell)
        identity = np.arange(ncell, dtype=np.int64)
        return cls(np.ones(ncell, dtype=np.int8), identity, identity.copy())

    @classmethod
    def from_actnum(cls, actnum, ncell: int) -> "ActiveMap":
        """`ACTNUM`-dan iki tərəfli xəritə qurur.

        `actnum` istənilən ədədi tipdə ola bilər (Eclipse onu float kimi
        də yaza bilir); `NaN` **aktiv** sayılır — deck-də `n*` defoltu ilə
        yazılmış, oxunmamış hüceyrəni səssizcə "yox etmək" təhlükəlidir.
        `> 0` aktivdir (Eclipse konvensiyası: 0 = qeyri-aktiv).
        """
        ncell = int(ncell)
        values = np.asarray(actnum, dtype=float).ravel()
        if values.size != ncell:
            raise ValueError(
                f"ACTNUM ölçüsü grid ilə uyğun gəlmir: {values.size} != {ncell}")
        mask = np.nan_to_num(values, nan=1.0) > 0.5
        global_to_active = np.full(ncell, -1, dtype=np.int64)
        active_to_global = np.flatnonzero(mask).astype(np.int64)
        global_to_active[active_to_global] = np.arange(
            active_to_global.size, dtype=np.int64)
        return cls(mask.astype(np.int8), global_to_active, active_to_global)

    # ────────────────────────────────────────────────────────── ölçülər
    @property
    def n_global(self) -> int:
        return int(self.actnum.size)

    @property
    def n_active(self) -> int:
        return int(self.active_to_global.size)

    @property
    def n_inactive(self) -> int:
        return self.n_global - self.n_active

    @property
    def mask(self) -> np.ndarray:
        """(ncell,) bool — aktivdirmi (qlobal indekslə)."""
        return self.actnum > 0

    @property
    def all_cells_active(self) -> bool:
        return self.n_active == self.n_global

    def is_active(self, cells) -> np.ndarray:
        """Verilmiş QLOBAL indeks(lər) aktivdirmi."""
        return self.actnum[np.asarray(cells, dtype=np.int64)] > 0

    # ──────────────────────────────────────────────────────── çevirmə
    def to_active(self, values: np.ndarray) -> np.ndarray:
        """Qlobal (ncell,…) massivdən aktiv (n_active,…) altmassiv."""
        return np.asarray(values)[self.active_to_global]

    def to_global(self, values: np.ndarray, fill: float = 0.0) -> np.ndarray:
        """Aktiv (n_active,) massivi qlobal (ncell,) massivə yayır;
        qeyri-aktiv hüceyrələr `fill` ilə doldurulur."""
        values = np.asarray(values)
        out = np.full(self.n_global, fill, dtype=values.dtype)
        out[self.active_to_global] = values
        return out

    def active_face_mask(self, cell_a: np.ndarray, cell_b: np.ndarray
                         ) -> np.ndarray:
        """Üz maskası: YALNIZ hər iki tərəfi aktiv olan üzlər üçün `True`
        (tapşırıq §3 — `ACTNUM[A] > 0 AND ACTNUM[B] > 0`)."""
        return self.is_active(cell_a) & self.is_active(cell_b)


@dataclass(frozen=True)
class Connections:
    """Hüceyrələr arası üzlərin siyahısı (qonşuluq qrafı)."""
    cell_a: np.ndarray
    cell_b: np.ndarray
    axis: np.ndarray

    @property
    def count(self) -> int:
        return int(self.cell_a.size)


@dataclass(frozen=True)
class CartesianGrid:
    """Struktur (kartezian) grid topologiyası.

    `actnum` — İSTƏYƏ BAĞLI (ncell,) ACTNUM massivi. Verilməyəndə (defolt
    `None`) BÜTÜN hüceyrələr aktivdir və bu sinif ƏVVƏLKİ KİMİ davranır —
    heç bir mövcud model/test dəyişmir. Verildikdə:

      * `build_connections()` qeyri-aktiv hüceyrəyə toxunan üzləri ATIR;
      * `active` (`ActiveMap`) qlobal↔aktiv indeks çevirməsini verir;
      * `ReservoirModel.pore_volume()` qeyri-aktiv hüceyrədə SIFIRDIR.

    Sahə `compare=False, repr=False` ilə elan olunub: `CartesianGrid`
    frozen dataclass-dır və `numpy` massivi ilə avtomatik `__eq__`/`__hash__`
    POZULARDI (massivin həqiqət dəyəri qeyri-müəyyəndir, həm də hash
    edilə bilmir). Bərabərlik/hash əvvəlki kimi YALNIZ (nx, ny, nz)
    üzərindədir.
    """
    nx: int
    ny: int
    nz: int = 1
    actnum: Optional[np.ndarray] = field(default=None, compare=False,
                                         repr=False)

    def __post_init__(self):
        if self.actnum is None:
            object.__setattr__(self, "_active", None)
            return
        active = ActiveMap.from_actnum(self.actnum, self.ncell)
        # normallaşdırılmış (int8, 0/1) forma saxlanılır — girişin tipi
        # (float, list, …) nə olursa olsun, sonrakı kod eynidir.
        object.__setattr__(self, "actnum", active.actnum)
        object.__setattr__(self, "_active", active)

    # ═════════════════════════════════════════════ aktiv hüceyrə xəritəsi
    @property
    def active(self) -> ActiveMap:
        """`ActiveMap` — ACTNUM verilməyibsə "hamısı aktiv" xəritəsi.

        Xəritə İLK sorğuda qurulur və saxlanılır (grid frozen-dir, ona
        görə `object.__setattr__`).
        """
        cached = getattr(self, "_active", None)
        if cached is None:
            cached = ActiveMap.all_active(self.ncell)
            object.__setattr__(self, "_active", cached)
        return cached

    @property
    def n_active(self) -> int:
        """Simulyasiyada İŞTİRAK EDƏN hüceyrə sayı (ACTNUM > 0)."""
        return self.active.n_active

    @property
    def has_inactive_cells(self) -> bool:
        return self.actnum is not None and not self.active.all_cells_active

    def with_actnum(self, actnum) -> "CartesianGrid":
        """Eyni topologiya, ACTNUM əlavə edilmiş YENİ grid (sinif frozen-dir).

        `actnum=None` verildikdə ACTNUM olmayan grid qaytarılır."""
        return CartesianGrid(self.nx, self.ny, self.nz, actnum)

    @property
    def ncell(self) -> int:
        return self.nx * self.ny * self.nz

    @property
    def shape(self) -> tuple:
        return (self.nz, self.ny, self.nx)

    def index(self, i: int, j: int, k: int = 0) -> int:
        return (k * self.ny + j) * self.nx + i

    def ijk(self, c: int) -> tuple:
        i = c % self.nx
        j = (c // self.nx) % self.ny
        k = c // (self.nx * self.ny)
        return i, j, k

    def ijk_array(self, cells: np.ndarray) -> tuple:
        """`ijk`-nin vektor forması — massiv üçün bir dəfə hesablanır."""
        cells = np.asarray(cells, dtype=np.int64)
        i = cells % self.nx
        j = (cells // self.nx) % self.ny
        k = cells // (self.nx * self.ny)
        return i, j, k

    def build_connections(self) -> Connections:
        idx = np.arange(self.ncell).reshape(self.shape)
        a_list, b_list, ax_list = [], [], []
        if self.nx > 1:
            a_list.append(idx[:, :, :-1].ravel())
            b_list.append(idx[:, :, 1:].ravel())
            ax_list.append(np.zeros(idx[:, :, :-1].size, dtype=np.int8))
        if self.ny > 1:
            a_list.append(idx[:, :-1, :].ravel())
            b_list.append(idx[:, 1:, :].ravel())
            ax_list.append(np.ones(idx[:, :-1, :].size, dtype=np.int8))
        if self.nz > 1:
            a_list.append(idx[:-1, :, :].ravel())
            b_list.append(idx[1:, :, :].ravel())
            ax_list.append(np.full(idx[:-1, :, :].size, 2, dtype=np.int8))
        if not a_list:
            e = np.zeros(0, dtype=int)
            return Connections(e, e, np.zeros(0, dtype=np.int8))
        cell_a = np.concatenate(a_list)
        cell_b = np.concatenate(b_list)
        axis = np.concatenate(ax_list)

        # ACTNUM (tapşırıq §3): üz YALNIZ hər İKİ tərəf aktiv olanda
        # yaradılır. Aktiv↔qeyri-aktiv üz HEÇ VAXT qurulmur — deməli
        # nə transmissivlik hesablanır, nə də seyrək matrisin
        # topologiyasına düşür (T_AB = 0 yazmaq kifayət deyildi: sıfır
        # element yenə də strukturda qalır və qeyri-aktiv hüceyrəni
        # sistemə naməlum kimi daxil edirdi).
        if self.has_inactive_cells:
            keep = self.active.active_face_mask(cell_a, cell_b)
            cell_a, cell_b, axis = cell_a[keep], cell_b[keep], axis[keep]

        return Connections(cell_a, cell_b, axis)
