# -*- coding: utf-8 -*-
"""Nyuton sistemini AKTİV hüceyrələrə (ACTNUM > 0) daraldan köməkçi.

NİYƏ LAZIMDIR
=============
`ACTNUM = 0` olan hüceyrə simulyasiyadan TAM çıxarılıb (bax
`domain/grid.py` modul sənədləşməsi):

    * məsamə həcmi SIFIRDIR      → akkumulyasiya həddi ≡ 0
    * qonşuluq bağlantısı YOXDUR → axın həddi ≡ 0
    * perforasiyası söndürülüb   → quyu həddi ≡ 0

Deməli həmin hüceyrənin Jakobiandakı İKİ sətri (təzyiq və doyumluluq)
TAMAMİLƏ SIFIRDIR — sistem SİNQULYARDIR və `spsolve`/BiCGStab NaN
qaytarır. Yəni "qeyri-aktiv hüceyrəni fizikadan çıxarmaq" işi YARIMÇIQ
qalsaydı (yalnız PV = 0 + bağlantı filtri), simulyasiya ÜMUMİYYƏTLƏ
işləməzdi.

NƏ EDİR
=======
Sistem AKTİV naməlumlara daraldılır — `A_aa · δ_a = −r_a`, sonra həll
qlobal vektora geri yayılır (qeyri-aktiv naməlumların dəyişməsi 0-dır).
Bu, sadəcə "sinqulyarlıqdan qaçmaq" deyil: xətti sistemin ÖLÇÜSÜ
`2·n_active` olur (tapşırıq §1), yəni qeyri-aktiv hüceyrələr nə yaddaş
tutur, nə də ön-şərtçi (ILU/CPR) vaxtı yeyir.

Daraltma DÜZGÜNDÜR (yaxınlaşma DEYİL): qeyri-aktiv sətir/sütunlar
tamamilə sıfır olduğuna görə `A_ai = A_ia = 0`, yəni tam sistem BLOK-
DİAQONALDIR və aktiv blok qalan hər şeydən ASILI DEYİL.

Bütün hüceyrələr aktiv olanda (`ActiveMap.all_active` — ACTNUM verilməmiş
HƏR mövcud model) `is_identity` doğrudur və bu sinif matrisə HEÇ
TOXUNMUR: nə kəsmə, nə kopyalama olur — davranış və performans
əvvəlkinin BİRƏBİR eynisidir.
"""

from __future__ import annotations
from typing import Optional, Tuple

import numpy as np
import scipy.sparse as sp

from ...domain.grid import ActiveMap
from .state import VARIABLES_PER_CELL


class ActiveDofReduction:
    """Qlobal ↔ aktiv SƏRBƏSTLİK DƏRƏCƏSİ (DOF) çevirməsi.

    `extra_dofs` — rezervuar naməlumlarından SONRA gələn, hüceyrə ilə
    bağlı OLMAYAN naməlumlar (`CoupledNewtonSolver`-də quyu BHP-ləri).
    Onlar HƏMİŞƏ saxlanılır: quyu hüceyrəyə aid deyil, aktiv hüceyrələrlə
    əlaqəlidir (qeyri-aktiv hüceyrədəki perforasiya artıq
    `PeacemanWellModel`-də atılıb).
    """

    def __init__(self, active: ActiveMap, extra_dofs: int = 0,
                 variables_per_cell: int = VARIABLES_PER_CELL):
        self.active = active
        self.variables_per_cell = int(variables_per_cell)
        self.extra_dofs = int(extra_dofs)
        self.reservoir_dofs = active.n_global * self.variables_per_cell
        self.size = self.reservoir_dofs + self.extra_dofs

        if active.all_cells_active:
            self._kept: Optional[np.ndarray] = None
            self.reduced_size = self.size
            return

        cell_dofs = (active.active_to_global[:, None] * self.variables_per_cell
                     + np.arange(self.variables_per_cell)[None, :]).ravel()
        if self.extra_dofs:
            cell_dofs = np.concatenate([
                cell_dofs,
                np.arange(self.reservoir_dofs, self.size, dtype=cell_dofs.dtype)])
        self._kept = cell_dofs.astype(np.int64)
        self.reduced_size = int(self._kept.size)

    # ─────────────────────────────────────────────────────────────────
    @property
    def is_identity(self) -> bool:
        """Daraltmağa EHTİYAC yoxdur (bütün hüceyrələr aktivdir)."""
        return self._kept is None

    @property
    def kept_dofs(self) -> np.ndarray:
        """Saxlanılan qlobal DOF indeksləri (identity halda `arange`)."""
        if self._kept is None:
            return np.arange(self.size, dtype=np.int64)
        return self._kept

    def restrict(self, matrix: sp.spmatrix, rhs: np.ndarray
                 ) -> Tuple[sp.spmatrix, np.ndarray]:
        """`(A[aktiv, aktiv], b[aktiv])` — daraldılmış sistem."""
        if self._kept is None:
            return matrix, rhs
        kept = self._kept
        reduced = matrix.tocsr()[kept][:, kept]
        return reduced, np.asarray(rhs)[kept]

    def expand(self, solution: np.ndarray) -> np.ndarray:
        """Daraldılmış həlli qlobal DOF vektoruna yayır (qeyri-aktiv
        naməlumların dəyişməsi SIFIRDIR — onlar həll olunmur)."""
        if self._kept is None:
            return solution
        full = np.zeros(self.size, dtype=float)
        full[self._kept] = solution
        return full
