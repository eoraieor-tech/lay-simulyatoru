"""FaciesField — kateqorik fasiya sahəsi (Phase 4.1 §8).

BU, `structure.RegionSet` (SATNUM/PVTNUM/SCAL region) İLƏ
QARIŞDIRILMAMALIDIR:

    RegionSet   — SCAL/PVT cədvəlinin SEÇİLMƏSİ üçün sabit, adətən
                  Eclipse-dən (SATNUM) gələn tam ədəd kod. Stoxastik
                  DEYİL, "seed"/"realizasiya" anlayışı YOXDUR.
    FaciesField — geoloji HETEROGENLİYİ təmsil edən, SIS ilə stoxastik
                  yaradılan kateqorik kod. Öz seed-i, variogram/
                  anizotropluq metadatası, realizasiya nömrəsi var.

`PropertyMap`-dan da FƏRQLİDİR — sadə ədədi massiv deyil, kod→ad lüğəti
VƏ realizasiyanın necə yaradıldığının tam qeydini (§9 "audit/report"
tələbinin fasiya analoqu) daşıyır.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np


@dataclass
class FaciesField:
    """Bir SIS realizasiyasının GeologicalModel-ə bağlana bilən forması."""
    name: str
    codes: np.ndarray                                    # (ncell,) tam ədəd fasiya kodları
    category_names: Dict[int, str] = field(default_factory=dict)
    realization_id: int = 0
    seed: int = 0
    requested_proportions: Dict[int, float] = field(default_factory=dict)
    realized_proportions: Dict[int, float] = field(default_factory=dict)
    variogram_metadata: Dict[int, dict] = field(default_factory=dict)
    anisotropy_metadata: Dict[int, dict] = field(default_factory=dict)
    conditioning_data_stats: Dict[str, object] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    def __post_init__(self):
        self.codes = np.asarray(self.codes, dtype=int).ravel()

    @property
    def ncell(self) -> int:
        return int(self.codes.size)

    @property
    def categories(self) -> List[int]:
        return sorted(self.requested_proportions) or sorted(set(np.unique(self.codes).tolist()))

    def label(self, code: int) -> str:
        return self.category_names.get(code, str(code))

    def summary(self) -> str:
        lines = [f"{self.name}: realizasiya #{self.realization_id} (seed={self.seed}), "
                f"{self.ncell} hüceyrə"]
        for code in self.categories:
            requested = self.requested_proportions.get(code)
            realized = self.realized_proportions.get(code)
            req_text = f"{requested:.3f}" if requested is not None else "—"
            real_text = f"{realized:.3f}" if realized is not None else "—"
            lines.append(f"  {self.label(code)} (kod {code}): tələb {req_text}, "
                        f"reallaşan {real_text}")
        if self.warnings:
            lines.append(f"  {len(self.warnings)} xəbərdarlıq (bax `.warnings`)")
        return "\n".join(lines)
