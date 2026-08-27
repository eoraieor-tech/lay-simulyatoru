"""Diaqnostika — səviyyəli yoxlama nəticələri.

Əvvəllər `validate()` sadəcə mətn siyahısı qaytarırdı və hər şey
"xəta" sayılırdı. Praktikada iki fərqli hal var:

    ERROR    — model işə salına bilməz (grid-dən kənar perforasiya)
    WARNING  — işə salına bilər, amma çox güman səhvdir
               (vurucu quyunun BHP-si lay təzyiqindən aşağıdır —
               quyu heç nə vurmayacaq, amma proqram susacaq)

Xəbərdarlıqlar bloklamır, çünki bəzən qəsdən belə edilir. Lakin
istifadəçi onları İŞƏ SALMAZDAN ƏVVƏL görməlidir.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional


class Severity(Enum):
    ERROR = "XƏTA"
    WARNING = "XƏBƏRDARLIQ"
    INFO = "MƏLUMAT"


@dataclass(frozen=True)
class Diagnostic:
    severity: Severity
    message: str
    source: str = ""          # hansı obyektə aiddir (quyu adı, panel...)
    hint: Optional[str] = None

    def __str__(self) -> str:
        prefix = f"[{self.severity.value}]"
        location = f" {self.source}:" if self.source else ""
        suffix = f"  → {self.hint}" if self.hint else ""
        return f"{prefix}{location} {self.message}{suffix}"


@dataclass
class DiagnosticReport:
    items: List[Diagnostic] = field(default_factory=list)

    def add(self, severity: Severity, message: str, source: str = "",
            hint: Optional[str] = None) -> None:
        self.items.append(Diagnostic(severity, message, source, hint))

    def error(self, message: str, source: str = "", hint: str = None) -> None:
        self.add(Severity.ERROR, message, source, hint)

    def warning(self, message: str, source: str = "", hint: str = None) -> None:
        self.add(Severity.WARNING, message, source, hint)

    def info(self, message: str, source: str = "", hint: str = None) -> None:
        self.add(Severity.INFO, message, source, hint)

    def extend(self, other: "DiagnosticReport") -> None:
        self.items.extend(other.items)

    def of(self, severity: Severity) -> List[Diagnostic]:
        return [item for item in self.items if item.severity is severity]

    @property
    def errors(self) -> List[Diagnostic]:
        return self.of(Severity.ERROR)

    @property
    def warnings(self) -> List[Diagnostic]:
        return self.of(Severity.WARNING)

    @property
    def has_errors(self) -> bool:
        return bool(self.errors)

    def messages(self, severity: Optional[Severity] = None) -> List[str]:
        items = self.items if severity is None else self.of(severity)
        return [item.message for item in items]

    def as_text(self, severity: Optional[Severity] = None) -> str:
        items = self.items if severity is None else self.of(severity)
        return "\n".join(str(item) for item in items)

    def __len__(self) -> int:
        return len(self.items)

    def __bool__(self) -> bool:
        return bool(self.items)
