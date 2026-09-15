"""Quyu məhdudiyyətləri — Nyuton həlləri ARASINDA tətbiq olunan qaydalar (B7).

Bu moduldakı hər şey `ThpController` (Q-15, Q-17) ilə eyni fəlsəfədədir:
QALIQ və JAKOBİAN yeni rejim tanımır. Qayda bağlantı obyektinin sahəsini
addımlar arasında yerində dəyişir, qalıq isə sadəcə yeni ədədi görür.

RATE HƏDƏFİNİN PERFORASİYALARA BÖLÜNMƏSİ (Seans 27)

ƏVVƏLKİ SƏHV (ölçüldü): qalıq, Jakobian və IMPES RATE hədəfini HƏR
perforasiyaya TAM yazırdı — `total = -abs(connection.target)` hər bağlantı
üçün ayrıca. Nəticədə 3 perforasiyalı quyu hədəfin 3 QATINI hasil edirdi
(50 m³/gün hədəf → 150 m³/gün). Heç bir xəbərdarlıq yox idi.

İNDİ hədəf quyunun bağlantıları arasında Peaceman nisbətində bölünür:

    pay_c = WI_c · λ_c / Σ_k WI_k · λ_k,        Σ pay = 1

Bu, BHP rejimində debitin təbii paylanmasıdır (`q_c = WI_c·λ_c·Δp`,
bütün perforasiyalarda eyni BHP).

NİYƏ λ ƏVVƏLKİ ADDIMDANDIR. Pay cari iterasiyanın λ-sı ilə hesablansaydı,
bir perforasiyanın debiti quyunun BÜTÜN digər hüceyrələrinin doymuşluğundan
asılı olardı — Jakobianda hüceyrələr arası yeni törəmələr (yeni seyrəklik
strukturu) yaranardı. Pay addımın ƏVVƏLİNDƏ (yığılmış vəziyyətdən) verilir və
Nyuton daxilində SABİTDİR, ona görə Jakobianda RATE törəmələri sadəcə paya
vurulur — struktur dəyişmir. Quyunun CƏMİ debiti dəqiqdir; yalnız təbəqələr
arası paylanma bir addım gecikir.

GERİYƏ UYĞUNLUQ: tək perforasiyalı quyuda pay DƏQİQ 1.0-dır (`x/x`), yəni
`target * 1.0` — nəticə bit-bit eyni.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence

from ..domain.wells import ControlMode


def _group_by_well(connections: Sequence) -> Dict[str, List[int]]:
    groups: Dict[str, List[int]] = {}
    for position, connection in enumerate(connections):
        groups.setdefault(connection.well_name, []).append(position)
    return groups


def needs_rate_allocation(connections: Sequence) -> bool:
    """Bir neçə perforasiyalı RATE quyusu varmı — yoxdursa mühərrik payı
    heç vaxt yeniləmir (əlavə hesab xərci yoxdur)."""
    connections = list(connections)
    return any(len(positions) > 1
               and connections[positions[0]].mode is ControlMode.RATE
               for positions in _group_by_well(connections).values())


def assign_rate_shares(connections: Sequence,
                       mobility: Optional[Sequence[float]] = None) -> None:
    """Hər bağlantının `rate_share`-ini yerində yazır.

    `mobility` — bağlantılarla EYNİ sırada LAY HƏCMİ mobilliyi (istismarçı
    `λw + λo`, vurucu vurulan fazanın son nöqtə mobilliyi). `None` olanda
    pay yalnız WI ilə verilir (mühərrik ilk addımda λ ilə yeniləyənə qədər).

    Bütün çəkilər sıfırdırsa (məs. hər perforasiyada λ = 0), WI nisbətinə,
    o da sıfırdırsa bərabər paya qayıdılır — cəm HƏMİŞƏ 1 qalır.
    """
    connections = list(connections)
    for positions in _group_by_well(connections).values():
        if len(positions) == 1:
            connections[positions[0]].rate_share = 1.0
            continue
        weights = [_weight(connections[p],
                           None if mobility is None else mobility[p])
                   for p in positions]
        total = sum(weights)
        if not total > 0.0:
            weights = [_weight(connections[p], None) for p in positions]
            total = sum(weights)
        for position, weight in zip(positions, weights):
            connections[position].rate_share = (
                weight / total if total > 0.0 else 1.0 / len(positions))


def _weight(connection, mobility: Optional[float]) -> float:
    value = float(connection.well_index)
    if mobility is not None:
        value *= float(mobility)
    return value if math.isfinite(value) and value > 0.0 else 0.0
