"""Qeyri-simmetrik xətti həlledici — Nyuton sistemi üçün.

IMPES-in təzyiq matrisi simmetrik idi, ona görə KQ (conjugate gradient)
işləyirdi. Jakobian isə simmetrik DEYİL — upstream çəkilənmə və quyu
üzvləri simmetriyanı pozur. Ona görə BiCGStab (və ya GMRES) lazımdır.

Ön-şərtçi olaraq ILU işlədilir. Jakobian bir Nyuton addımı ərzində
kəskin dəyişdiyi üçün ön-şərtçi HƏR NYUTON İTERASİYASINDA yenilənir —
IMPES-dəki 50 addımlıq gecikdirmə burada yaramır.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ...logging_setup import get_logger
from .cpr import CprConfig, CprPreconditioner

LOG = get_logger(__name__)


@dataclass
class NewtonLinearSolverConfig:
    tolerance: float = 1e-8
    max_iterations: int = 300
    ilu_drop_tolerance: float = 1e-4
    ilu_fill_factor: float = 10.0
    fallback_to_direct: bool = True
    direct_threshold: int = 20_000
    """Bu ölçüdən kiçik sistemlərdə birbaşa həll (`splu`) daha sürətlidir."""

    preconditioner: str = "auto"
    """"ilu", "cpr" və ya "auto".

    Ölçmə (bax `A6_PLAN.md`, mərhələ 5): orta ölçülü sistemlərdə güclü
    ILU 1 iterasiyada həll edir və CPR-dən sürətlidir. CPR-in üstünlüyü
    YADDAŞDADIR — blok-Jakobi hamarlayıcı ilə o, 2× fill tutur, ILU isə
    5×. Böyük modellərdə ILU yaddaşa sığmır.

    "auto" — sistem `cpr_threshold`-dan böyükdürsə CPR, yoxsa ILU."""

    cpr_threshold: int = 40_000
    cpr_config: Optional[CprConfig] = None


class NewtonLinearSolver:
    """BiCGStab + ILU, kiçik sistemlərdə birbaşa LU."""

    def __init__(self, config: Optional[NewtonLinearSolverConfig] = None):
        self.config = config or NewtonLinearSolverConfig()
        self.last_iterations = 0
        self.last_method = ""
        self.last_preconditioner = ""

    def solve(self, matrix: sp.csr_matrix, rhs: np.ndarray,
              x0: Optional[np.ndarray] = None) -> np.ndarray:
        config = self.config
        size = rhs.size

        if size <= config.direct_threshold:
            solution = self._direct(matrix, rhs)
            if solution is not None:
                self.last_method = "LU"
                self.last_iterations = 0
                return solution

        preconditioner = self._preconditioner(matrix)
        counter = {"count": 0}

        def callback(_):
            counter["count"] += 1

        solution, info = spla.bicgstab(
            matrix, rhs, x0=x0, rtol=config.tolerance, atol=0.0,
            maxiter=config.max_iterations, M=preconditioner, callback=callback)
        self.last_iterations = counter["count"]
        self.last_method = "BiCGStab"

        if info != 0 or not np.all(np.isfinite(solution)):
            LOG.debug("BiCGStab yığılmadı (info=%s), birbaşa həllə keçilir", info)
            direct = self._direct(matrix, rhs)
            if direct is None:
                raise FloatingPointError("Nyuton sistemi həll oluna bilmədi.")
            self.last_method = "LU (ehtiyat)"
            return direct
        return solution

    # ------------------------------------------------------------ daxili
    @staticmethod
    def _direct(matrix: sp.csr_matrix, rhs: np.ndarray) -> Optional[np.ndarray]:
        try:
            solution = spla.spsolve(matrix.tocsc(), rhs)
        except (RuntimeError, ValueError, spla.MatrixRankWarning):
            return None
        return solution if np.all(np.isfinite(solution)) else None

    def _preconditioner(self, matrix: sp.csr_matrix):
        config = self.config
        use_cpr = (config.preconditioner == "cpr"
                   or (config.preconditioner == "auto"
                       and matrix.shape[0] >= config.cpr_threshold))

        if use_cpr:
            try:
                preconditioner = CprPreconditioner(matrix, config.cpr_config)
                self.last_preconditioner = f"CPR ({preconditioner.smoother_kind})"
                return preconditioner
            except (ValueError, RuntimeError) as error:
                LOG.debug("CPR qurula bilmədi (%s), ILU-ya keçilir", error)

        try:
            ilu = spla.spilu(matrix.tocsc(),
                             drop_tol=config.ilu_drop_tolerance,
                             fill_factor=config.ilu_fill_factor)
            self.last_preconditioner = "ILU"
            return spla.LinearOperator(matrix.shape, matvec=ilu.solve)
        except RuntimeError:
            LOG.debug("ILU qurula bilmədi — ön-şərtçisiz davam edilir")
            self.last_preconditioner = "yoxdur"
            return None
