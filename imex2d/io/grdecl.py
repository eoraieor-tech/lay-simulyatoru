"""GRDECL oxuyucusu — Eclipse/Petrel grid formatı.

Format sadədir: açar söz, sonra dəyərlər, sonra `/`.

    SPECGRID
      41 41 5 1 F /
    PORO
      8405*0.22 /
    PERMX
      100 120 3*150 200 /

TƏKRAR SİNTAKSİSİ: `n*value` — dəyər n dəfə təkrarlanır. Bu, real
fayllarda mütləqdir: 200 000 hüceyrəli modeldə hər dəyəri ayrıca
yazmaq faylı yüzlərlə meqabayt edərdi.

DƏSTƏKLƏNƏN HƏNDƏSƏ

    Block-centered (DX/DY/DZ/TOPS)   tam dəstək
    Corner-point (COORD/ZCORN)       oxunur, bərabər bloka
                                     APPROKSİMASİYA olunur

İkinci hal dürüst adlandırılmalıdır: `CellGeometry` hazırda yalnız
bərabər ölçülü bloklar saxlayır (bax `ARCHITECTURE.md`, 5.1). Corner-point
faylı oxunanda orta hüceyrə ölçüsü hesablanır və istifadəçi
xəbərdarlıq alır — nəticələr orijinal həndəsə ilə tam üst-üstə
düşməyəcək.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from ..domain.diagnostics import DiagnosticReport
from ..logging_setup import get_logger

LOG = get_logger(__name__)

# şərh: `--` sətrin sonuna qədər
_COMMENT = re.compile(r"--.*")
_REPEAT = re.compile(r"^(\d+)\*(.*)$")

# hüceyrə üzrə massiv gözlənilən açar sözlər
CELL_ARRAYS = {
    "PORO", "PERMX", "PERMY", "PERMZ", "NTG", "ACTNUM", "SATNUM",
    "PVTNUM", "EQLNUM", "FIPNUM", "DX", "DY", "DZ", "TOPS", "MULTX",
    "MULTY", "MULTZ", "SWAT", "PRESSURE",
}

# yalnız bir sətir ölçü məlumatı
SCALAR_KEYWORDS = {"SPECGRID", "DIMENS"}

# Bölmə başlıqları və bayraq açar sözləri `/` TƏLƏB ETMİR.
# Bu, oxumada incə, lakin ağır səhv mənbəyidir: `GRID` açar sözünü
# adi massiv kimi qəbul etsək, oxucu növbəti `/` işarəsinə qədər hər
# şeyi udur və ilk massiv (adətən DX) itir. Nəticədə model səssizcə
# defolt ölçülərlə qurulur.
STANDALONE_KEYWORDS = {
    "RUNSPEC", "GRID", "EDIT", "PROPS", "REGIONS", "SOLUTION", "SUMMARY",
    "SCHEDULE", "END", "ENDINC", "INIT", "NOECHO", "ECHO",
    "OIL", "WATER", "GAS", "DISGAS", "VAPOIL",
    "METRIC", "FIELD", "LAB", "NOSIM", "UNIFOUT", "UNIFIN",
}


class GrdeclError(Exception):
    """Fayl oxuna bilmədikdə və ya struktur uyğun gəlmədikdə."""


@dataclass
class GrdeclDeck:
    """Oxunmuş fayl — hələ domain obyekti deyil."""
    dimensions: Optional[tuple] = None          # (nx, ny, nz)
    arrays: Dict[str, np.ndarray] = field(default_factory=dict)
    keywords_seen: List[str] = field(default_factory=list)
    source: str = ""

    @property
    def ncell(self) -> int:
        if self.dimensions is None:
            return 0
        nx, ny, nz = self.dimensions
        return nx * ny * nz

    def has(self, keyword: str) -> bool:
        return keyword in self.arrays

    def get(self, keyword: str) -> Optional[np.ndarray]:
        return self.arrays.get(keyword)

    def summary(self) -> dict:
        return {
            "ölçü": ("×".join(str(n) for n in self.dimensions)
                     if self.dimensions else "—"),
            "hüceyrə": self.ncell,
            "massivlər": ", ".join(sorted(self.arrays)) or "—",
            "açar sözlər": len(self.keywords_seen),
        }


# ══════════════════════════════════════════════════════════ oxuma

def _expand(tokens: List[str]) -> List[float]:
    """`3*0.25` -> [0.25, 0.25, 0.25]."""
    values: List[float] = []
    for token in tokens:
        match = _REPEAT.match(token)
        if match:
            count = int(match.group(1))
            body = match.group(2)
            if body == "":                 # `3*` — Eclipse-də "defolt"
                values.extend([np.nan] * count)
                continue
            values.extend([float(body)] * count)
        else:
            values.append(float(token))
    return values


def _tokenise(text: str) -> List[str]:
    lines = []
    for line in text.splitlines():
        line = _COMMENT.sub("", line)
        if line.strip():
            lines.append(line)
    return " ".join(lines).replace("/", " / ").split()


def read_grdecl(path: str, report: Optional[DiagnosticReport] = None
                ) -> GrdeclDeck:
    """GRDECL faylını oxuyur. `INCLUDE` dəstəklənmir (xəbərdarlıq verilir)."""
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            text = handle.read()
    except OSError as error:
        raise GrdeclError(f"Fayl açıla bilmədi: {error}") from error

    deck = GrdeclDeck(source=path)
    tokens = _tokenise(text)
    index = 0

    while index < len(tokens):
        token = tokens[index].upper()
        index += 1

        if token == "/":
            continue
        if token in STANDALONE_KEYWORDS:
            deck.keywords_seen.append(token)
            continue
        if not token.isalpha() and not token.replace("_", "").isalnum():
            continue
        if token[0].isdigit():
            continue

        if token == "INCLUDE":
            if report is not None:
                report.warning(
                    "INCLUDE açar sözü dəstəklənmir — həmin fayl "
                    "oxunmayacaq.", "GRDECL",
                    "Faylları bir sənəddə birləşdir")
            while index < len(tokens) and tokens[index] != "/":
                index += 1
            continue

        deck.keywords_seen.append(token)

        # açar sözdən `/` işarəsinə qədər olan hissə
        body: List[str] = []
        while index < len(tokens) and tokens[index] != "/":
            body.append(tokens[index])
            index += 1
        index += 1                          # `/` işarəsini keç

        if token in SCALAR_KEYWORDS:
            numbers = [t for t in body if t.replace(".", "").isdigit()]
            if len(numbers) >= 3:
                deck.dimensions = tuple(int(float(n)) for n in numbers[:3])
            continue

        if token in CELL_ARRAYS or token in ("COORD", "ZCORN"):
            try:
                deck.arrays[token] = np.asarray(_expand(body), dtype=float)
            except ValueError as error:
                raise GrdeclError(
                    f"'{token}' massivi oxuna bilmədi: {error}") from error

    if deck.dimensions is None:
        raise GrdeclError("Faylda SPECGRID/DIMENS açar sözü tapılmadı.")
    _validate_sizes(deck, report)
    return deck


def _validate_sizes(deck: GrdeclDeck, report: Optional[DiagnosticReport]
                    ) -> None:
    """Massiv uzunluqlarını grid ölçüsü ilə tutuşdurur."""
    nx, ny, nz = deck.dimensions
    expected = nx * ny * nz

    for keyword, values in list(deck.arrays.items()):
        if keyword in ("COORD", "ZCORN"):
            continue
        if keyword == "TOPS" and values.size == nx * ny:
            continue                        # yalnız üst təbəqə — qanunidir
        if values.size == expected:
            continue
        message = (f"'{keyword}' massivində {values.size} dəyər var, "
                   f"gözlənilən {expected}.")
        if report is not None:
            report.error(message, "GRDECL")
        raise GrdeclError(message)
