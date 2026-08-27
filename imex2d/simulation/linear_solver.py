"""SciPy əsaslı xətti həlledici — ILinearSolver implementasiyası.

Mövcud koddan köçürülüb (KQ + ILU ön-şərtçi, warm start, birbaşa həllə
qayıtma). Ayrıca sinfə çıxarılmasının səbəbi: həlledicini dəyişmək üçün
mühərrikə toxunmaq lazım gəlməsin.
"""

from __future__ import annotations
from typing import Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ..application.config import LinearSolverConfig
from ..interfaces.services import ILinearSolver


class ScipyCgIluSolver(ILinearSolver):

    def __init__(self, config: Optional[LinearSolverConfig] = None):
        self.config = config or LinearSolverConfig()
        self._preconditioner = None
        self._age = 0

    def reset(self) -> None:
        self._preconditioner = None
        self._age = 0

    def solve(self, matrix, rhs: np.ndarray, x0: Optional[np.ndarray] = None) -> np.ndarray:
        cfg = self.config
        n = rhs.size
        if (self._preconditioner is None
                or self._age >= cfg.preconditioner_refresh_steps):
            try:
                ilu = spla.spilu(matrix.tocsc(),
                                 drop_tol=cfg.ilu_drop_tolerance,
                                 fill_factor=cfg.ilu_fill_factor)
                self._preconditioner = spla.LinearOperator((n, n), matvec=ilu.solve)
                self._age = 0
            except RuntimeError:
                self._preconditioner = None
        self._age += 1

        x, info = spla.cg(matrix, rhs, x0=x0, rtol=cfg.tolerance, atol=0.0,
                          maxiter=cfg.max_iterations, M=self._preconditioner)
        if (info != 0 or not np.all(np.isfinite(x))) and cfg.fallback_to_direct:
            self._preconditioner = None
            x = spla.spsolve(matrix.tocsr(), rhs)
        if not np.all(np.isfinite(x)):
            raise FloatingPointError("Xətti sistem həll oluna bilmədi.")
        return x
