"""Quyu lüləsi (tubing) həndəsəsi — YALNIZ məlumat, hesablama yoxdur.

`domain/wells.py`-dakı `Well` ilə eyni prinsip: bu modul quyu lüləsinin
NƏ olduğunu təsvir edir, təzyiq düşgüsünü NECƏ hesablamağı yox. Traverse
hesabatı `simulation/wellbore/` altındadır.

V1 MƏHDUDİYYƏTİ — TAM ŞAQULİ LÜLƏ. Ölçülmüş dərinlik (MD) ilə şaquli
dərinlik (TVD) EYNİ sayılır. Layihədəki bütün quyular `Well.vertical()`
ilə qurulur, ona görə bu, hazırkı modellər üçün dəqiqdir. Əyri (deviated)
quyu üçün `TubingGeometry`-yə MD→TVD profili əlavə olunmalıdır —
⏳ hələ YOXDUR, uydurulmur.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List

#: Kommersiya boru kələ-kötürlüyü (mütləq), m. Yeni karbon-polad boru üçün
#: sənaye standartı ~0.0018 hüceyrə (0.06 mm) — Beggs & Brill, "Two-Phase
#: Flow in Pipes" cədvəlindəki `commercial steel` dəyəri.
DEFAULT_ROUGHNESS_M = 6.0e-5

#: 2⅞ düym (tubing) daxili diametri, m — sahədə ən yayılmış ölçü.
DEFAULT_DIAMETER_M = 0.062


@dataclass
class TubingGeometry:
    """Bir quyunun lülə borusu.

    `diameter` — borunun DAXİLİ diametri (m). Sürtünmə və sürət bundan
    çıxır, ona görə xarici diametr YARAMAZ.

    `roughness` — mütləq kələ-kötürlük ε (m), nisbi deyil. Chen (1979)
    düsturu ε/d nisbətini özü qurur.

    `wellhead_depth` — quyu başının dərinliyi (m). Defolt 0.0, yəni
    səthdə. Dəniz platforması / su altı quyu başı üçün müsbət dəyər
    verilir; traverse yalnız `wellhead_depth`-dən perforasiyaya qədər
    olan intervalı hesablayır.

    `segments` — traverse neçə hissəyə bölünsün. Qaz yuxarı qalxdıqca
    genişlənir, ona görə ρ və v lülə boyunca kəskin dəyişir; tək
    seqmentli hesabat yalnız qazsız quyuda düzgündür. Bax
    `simulation/wellbore/traverse.py`.
    """

    diameter: float = DEFAULT_DIAMETER_M
    roughness: float = DEFAULT_ROUGHNESS_M
    wellhead_depth: float = 0.0
    segments: int = 20

    # ───────────────────────────────────────────────── törəmə kəmiyyətlər
    @property
    def area(self) -> float:
        """Axın en kəsiyi, m²."""
        import math
        return math.pi * self.diameter ** 2 / 4.0

    @property
    def relative_roughness(self) -> float:
        """ε/d — Chen düsturunun arqumenti."""
        return self.roughness / self.diameter

    def length_to(self, perforation_depth: float) -> float:
        """Quyu başından perforasiyaya qədər lülə uzunluğu, m.

        V1-də şaquli olduğu üçün bu, həm MD, həm TVD-dir.
        """
        return float(perforation_depth) - float(self.wellhead_depth)

    # ─────────────────────────────────────────────────────── yoxlamalar
    def validate(self) -> List[str]:
        """Sərt fiziki xətalar — modelin qurulmasını dayandırır."""
        issues: List[str] = []
        if not (self.diameter > 0.0):
            issues.append("lülə diametri müsbət olmalıdır")
        elif self.diameter > 1.0:
            issues.append(f"lülə diametri qeyri-real böyükdür "
                          f"({self.diameter:g} m) — dəyər metrlədirmi?")
        if self.roughness < 0.0:
            issues.append("kələ-kötürlük mənfi ola bilməz")
        elif self.diameter > 0.0 and self.roughness >= self.diameter:
            issues.append("kələ-kötürlük diametrdən kiçik olmalıdır")
        if self.wellhead_depth < 0.0:
            issues.append("quyu başının dərinliyi mənfi ola bilməz")
        if self.segments < 1:
            issues.append("seqment sayı ən azı 1 olmalıdır")
        return issues

    def validate_warnings(self) -> List[str]:
        """Şübhəli, lakin qanuni dəyərlər."""
        warnings: List[str] = []
        if self.diameter > 0.0 and self.relative_roughness > 0.05:
            warnings.append(
                f"nisbi kələ-kötürlük çox böyükdür "
                f"(ε/d = {self.relative_roughness:.3f}) — Chen düsturu "
                f"ε/d ≤ 0.05 aralığında doğrulanıb")
        if 1 <= self.segments < 5:
            warnings.append(
                f"seqment sayı azdır ({self.segments}) — sərbəst qazı olan "
                f"quyuda sıxlıq lülə boyunca kəskin dəyişir, nəticə kobud olar")
        return warnings
