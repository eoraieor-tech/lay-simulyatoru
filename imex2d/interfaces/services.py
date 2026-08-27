"""Texniki xidmət interfeysləri — xətti həlledici, proqres, mühərrik."""

from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Optional

import numpy as np


class ILinearSolver(ABC):
    """Xətti sistem həlledicisi.

    Ayrılma səbəbi: gələcəkdə AMG, PETSc və ya GPU həlledicisi əlavə
    ediləndə mühərrikin kodu dəyişməməlidir.
    """

    @abstractmethod
    def solve(self, matrix, rhs: np.ndarray, x0: Optional[np.ndarray] = None) -> np.ndarray: ...

    def reset(self) -> None:
        """Ön-şərtçi keşini təmizləyir."""


class IProgressReporter(ABC):
    """Gedişat bildirişi. Mühərrik UI-dən deyil, bu interfeysdən asılıdır."""

    @abstractmethod
    def report(self, fraction: float, message: str) -> bool:
        """False qaytarsa, mühərrik dayanır."""


class NullProgressReporter(IProgressReporter):
    """Skript rejimi üçün — heç nə etmir."""

    def report(self, fraction: float, message: str) -> bool:
        return True


class ISimulationEngine(ABC):
    """Hesablama mühərriki. ReservoirModel QƏBUL EDİR, yaratmır."""

    @abstractmethod
    def run(self, reporter: Optional[IProgressReporter] = None): ...
