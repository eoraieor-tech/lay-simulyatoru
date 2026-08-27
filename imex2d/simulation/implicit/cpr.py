"""CPR ön-şərtçisi — A6, mərhələ 5.

CPR = Constrained Pressure Residual. İki mərhələli ön-şərtçi:

    1) Təzyiq alt-sistemi ayrılır və dəqiq həll olunur
    2) Qalan qalıq tam sistemə ILU ilə hamarlanır

Səbəb: Jakobian iki fərqli xarakterli tənliyi birləşdirir —

    təzyiq        elliptikdir, uzaqmənzilli, qlobal həll tələb edir
    doyumluluq    hiperbolikdir, lokal, upstream istiqamətdə yayılır

Tək ILU hər ikisini eyni cür emal edir və elliptik hissədə zəif qalır.
CPR təzyiqi ayırıb ona uyğun həlledici verir.

DEKUPLİNQ (quasi-IMPES)
Hər hüceyrənin 2×2 diaqonal bloku:

    D_c = [[∂R_w/∂p, ∂R_w/∂Sw],
           [∂R_o/∂p, ∂R_o/∂Sw]]

Təzyiq tənliyi iki tənliyin çəkili cəmidir; çəkilər elə seçilir ki,
doyumluluq törəməsi YOX OLSUN:

    w = [ ∂R_o/∂Sw,  −∂R_w/∂Sw ]
    w · [∂R_w/∂Sw, ∂R_o/∂Sw]ᵀ = 0

Sabit B halında bu, `w = [Bw, Bo]` — yəni klassik həcm balansıdır.

TƏTBİQ

    M⁻¹r = P·A_p⁻¹·(W·r)  +  M_ILU⁻¹·(r − A·P·A_p⁻¹·(W·r))

burada W — çəkili restriction (N×2N), P — prolongation (2N×N).
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ...logging_setup import get_logger
from .state import VARIABLES_PER_CELL

LOG = get_logger(__name__)


@dataclass
class CprConfig:
    pressure_tolerance: float = 1e-3
    """Təzyiq alt-sistemi DƏQİQ həll olunmamalıdır — o, yalnız
    ön-şərtçinin bir hissəsidir. Sıx tolerans vaxt itkisidir."""

    pressure_max_iterations: int = 30
    pressure_drop_tolerance: float = 1e-4
    pressure_fill_factor: float = 8.0
    global_drop_tolerance: float = 1e-4
    global_fill_factor: float = 5.0
    smoother: str = "auto"
    """İkinci mərhələnin hamarlayıcısı: "ilu", "block_jacobi" və ya "auto".

    "auto" — əvvəlcə ILU sınanır, uğursuz olarsa blok-Jakobiyə keçilir.
    Ölçmə göstərdi ki, böyük sistemlərdə (10 000+ hüceyrə) ILU
    ümumiyyətlə qurula bilmir; hamarlayıcısız CPR isə yığılmır
    (500 iterasiya). Blok-Jakobi HƏMİŞƏ qurulur və ucuzdur —
    hər hüceyrənin 2×2 diaqonal bloku ayrıca tərslənir."""

    use_amg: bool = True
    """pyamg mövcuddursa təzyiq sistemi üçün AMG işlədilsin."""


class CprPreconditioner(spla.LinearOperator):
    """Jakobian üçün iki mərhələli ön-şərtçi."""

    def __init__(self, matrix: sp.csr_matrix, config: Optional[CprConfig] = None):
        self.config = config or CprConfig()
        self.matrix = matrix.tocsr()
        size = matrix.shape[0]
        if size % VARIABLES_PER_CELL:
            raise ValueError("Sistem ölçüsü hüceyrə başına dəyişən sayına "
                             "bölünmür.")
        self.ncell = size // VARIABLES_PER_CELL
        super().__init__(dtype=np.float64, shape=(size, size))

        self._build_transfer_operators()
        self._build_pressure_solver()
        self._build_smoother()
        self.pressure_iterations = 0

    # ═══════════════════════════════════════════ dekuplinq və ötürmə
    def _build_transfer_operators(self) -> None:
        """Çəkili restriction W (N×2N) və prolongation P (2N×N)."""
        cells = np.arange(self.ncell)
        water_rows = cells * VARIABLES_PER_CELL
        oil_rows = water_rows + 1
        saturation_columns = water_rows + 1

        # diaqonal blokun doyumluluq sütunu
        d_water = np.asarray(
            self.matrix[water_rows, saturation_columns]).ravel()
        d_oil = np.asarray(
            self.matrix[oil_rows, saturation_columns]).ravel()

        weight_water = d_oil
        weight_oil = -d_water
        scale = np.maximum(np.abs(weight_water) + np.abs(weight_oil), 1e-30)
        weight_water = weight_water / scale
        weight_oil = weight_oil / scale

        # dekuplinq mümkün olmayan hüceyrələrdə sadə cəmə qayıdırıq
        degenerate = np.abs(d_water) + np.abs(d_oil) < 1e-30
        weight_water[degenerate] = 1.0
        weight_oil[degenerate] = 1.0

        self.restriction = sp.csr_matrix(
            (np.concatenate([weight_water, weight_oil]),
             (np.concatenate([cells, cells]),
              np.concatenate([water_rows, oil_rows]))),
            shape=(self.ncell, self.shape[0]))

        self.prolongation = sp.csr_matrix(
            (np.ones(self.ncell), (water_rows, cells)),
            shape=(self.shape[0], self.ncell))

    # ═══════════════════════════════════════════ təzyiq alt-sistemi
    def _build_pressure_solver(self) -> None:
        config = self.config
        self.pressure_matrix = (self.restriction @ self.matrix
                                @ self.prolongation).tocsr()

        self._amg = None
        if config.use_amg:
            try:
                import pyamg
                self._amg = pyamg.smoothed_aggregation_solver(
                    self.pressure_matrix.tocsr())
                LOG.debug("CPR: təzyiq sistemi üçün AMG quruldu")
                return
            except ImportError:
                pass
            except Exception as error:
                LOG.debug("CPR: AMG qurula bilmədi (%s)", error)

        try:
            ilu = spla.spilu(self.pressure_matrix.tocsc(),
                             drop_tol=config.pressure_drop_tolerance,
                             fill_factor=config.pressure_fill_factor)
            self._pressure_preconditioner = spla.LinearOperator(
                self.pressure_matrix.shape, matvec=ilu.solve)
        except RuntimeError:
            self._pressure_preconditioner = None

    def _solve_pressure(self, residual: np.ndarray) -> np.ndarray:
        if self._amg is not None:
            return self._amg.solve(residual,
                                   tol=self.config.pressure_tolerance,
                                   maxiter=self.config.pressure_max_iterations)
        counter = {"n": 0}

        def callback(_):
            counter["n"] += 1

        solution, info = spla.gmres(
            self.pressure_matrix, residual,
            rtol=self.config.pressure_tolerance, atol=0.0,
            maxiter=self.config.pressure_max_iterations,
            M=self._pressure_preconditioner, callback=callback,
            callback_type="pr_norm")
        self.pressure_iterations += counter["n"]
        if not np.all(np.isfinite(solution)):
            return np.zeros_like(residual)
        return solution

    # ═══════════════════════════════════════════════ hamarlayıcı
    def _build_smoother(self) -> None:
        config = self.config
        self.smoother_kind = "yoxdur"

        if config.smoother in ("ilu", "auto"):
            try:
                ilu = spla.spilu(self.matrix.tocsc(),
                                 drop_tol=config.global_drop_tolerance,
                                 fill_factor=config.global_fill_factor)
                self._smoother = ilu.solve
                self.smoother_kind = "ILU"
                return
            except RuntimeError:
                if config.smoother == "ilu":
                    LOG.debug("CPR: qlobal ILU qurula bilmədi")
                    self._smoother = None
                    return
                LOG.debug("CPR: ILU alınmadı, blok-Jakobiyə keçilir")

        self._build_block_jacobi()

    def _build_block_jacobi(self) -> None:
        """Hər hüceyrənin 2×2 diaqonal bloku ayrıca tərslənir.

        Bu, həmişə qurulur (yalnız 4 ədəd oxunuş və 2×2 tərsləmə) və
        yaddaş tələb etmir. Cəbri baxımdan zəifdir, lakin CPR-in birinci
        mərhələsi elliptik hissəni onsuz da həll edib — ikinci mərhələyə
        yalnız lokal (hiperbolik) qalığı təmizləmək qalır.
        """
        cells = np.arange(self.ncell)
        rows = cells * VARIABLES_PER_CELL
        a = np.asarray(self.matrix[rows, rows]).ravel()
        b = np.asarray(self.matrix[rows, rows + 1]).ravel()
        c = np.asarray(self.matrix[rows + 1, rows]).ravel()
        d = np.asarray(self.matrix[rows + 1, rows + 1]).ravel()

        determinant = a * d - b * c
        singular = np.abs(determinant) < 1e-30
        determinant = np.where(singular, 1.0, determinant)

        self._inverse = (np.where(singular, 1.0, d / determinant),
                         np.where(singular, 0.0, -b / determinant),
                         np.where(singular, 0.0, -c / determinant),
                         np.where(singular, 1.0, a / determinant))
        self._smoother = self._apply_block_jacobi
        self.smoother_kind = "blok-Jakobi"

    def _apply_block_jacobi(self, residual: np.ndarray) -> np.ndarray:
        inv_a, inv_b, inv_c, inv_d = self._inverse
        water = residual[0::VARIABLES_PER_CELL]
        oil = residual[1::VARIABLES_PER_CELL]
        result = np.empty_like(residual)
        result[0::VARIABLES_PER_CELL] = inv_a * water + inv_b * oil
        result[1::VARIABLES_PER_CELL] = inv_c * water + inv_d * oil
        return result

    # ═══════════════════════════════════════════════════ tətbiq
    def _matvec(self, residual: np.ndarray) -> np.ndarray:
        residual = np.asarray(residual, float).ravel()

        # 1-ci mərhələ: təzyiq
        pressure_correction = self.prolongation @ self._solve_pressure(
            self.restriction @ residual)

        # 2-ci mərhələ: qalan qalığın hamarlanması
        if self._smoother is None:
            return pressure_correction
        remainder = residual - self.matrix @ pressure_correction
        return pressure_correction + self._smoother(remainder)
