"""QUYULARDAN GƏLƏN STRUKTUR — interpolyasiya edilmiş lay üstü/altı
səthlərindən corner-point təpələri.

PROBLEM (bax `ISH_HESABATI.md`, A2). Geologiya cədvəlində hər quyu üçün
«lay üstü» (TOP) və «lay altı» (BOTTOM) dərinlikləri var; onlar
`geology_adapter.AREAL_TARGETS` ilə areal nümunə kimi toplanır və
`geology_service` tərəfindən interpolyasiya olunur, YƏNİ
`model.property_maps["TOP"]/["BOTTOM"]` HAZIRDIR. Amma grid həndəsəsi
onlardan QURULMURDU: `CellGeometry` `top_depth + i·dip_x + j·dip_y` düz
maili müstəvidir və `dz` LAY ÜZRƏDİR (`(nz,)`), yəni sütundan sütuna
dəyişən qalınlığı ifadə EDƏ BİLMİR. Nəticədə istifadəçinin doldurduğu
struktur dərinlikləri məsamə həcminə, cazibəyə, equilibration-a və
perforasiya→K uyğunlaşdırmasına HEÇ CÜR TƏSİR ETMİRDİ — SƏSSİZCƏ.

BU MODULUN İŞİ. Saf həndəsə çevirməsidir: `(nx·ny,)` TOP və BOTTOM
səthlərindən `(ncell, 8, 3)` təpə massivi qurur. Həndəsə MOTORU burada
YAZILMIR — `CornerPointGeometry` (bax `corner_point_geometry.py`) artıq
dəqiq çoxüzlü həcm, sütun-üzrə lay sərhədləri, əyri üz sahəsi və
Peaceman effektiv həndəsəsi verir; bu modul yalnız onun
`from_nodes()` konstruktoruna GİRİŞ hazırlayır.

Modul QƏSDƏN UI-dən və application qatından ASILI DEYİL (yalnız numpy +
`CartesianGrid`) — testi də, istifadəsi də tək başına mümkündür.

MƏHDUDİYYƏT (gizlədilmir): PİLLARLAR ŞAQULİDİR
==============================================
X/Y hələ də `dx`/`dy` düzbucaqlı şəbəkəsidir — yalnız Z dəyişir. Yəni
bu modul MAİLİ PİLLAR, FAY ATIMI (sıçrayış) və PINCH-OUT (sıfır
qalınlığa yığılan lay) QURMUR. Belə həndəsə lazımdırsa mənbə `.GRDECL`
COORD/ZCORN-dur (`CornerPointGeometry.from_grdecl`), quyu cədvəli deyil:
dörd quyudan interpolyasiya edilmiş iki səth fay atımını təyin edə
bilməz.

LAYLARA BÖLGÜ: `proportional`
=============================
Hər `(i, j)` sütunu üçün `H = BOTTOM − TOP`, sonra `nz` BƏRABƏR hissəyə:

    z_k = TOP + H · k / nz,      k = 0 … nz

Yəni layların NİSBİ qalınlığı bütün sahədə sabit, MÜTLƏQ qalınlığı isə
sütundan sütuna dəyişir. Bu, Eclipse-in "proportional" bölgüsüdür və
`spec.dz`-i ƏVƏZ EDİR (bu rejimdə `dz` ARTIQ İŞTİRAK ETMİR — çağıran
bunu istifadəçiyə AÇIQ bildirməlidir, bax `geology_service.
_build_structural_geometry`).
"""

from __future__ import annotations

from typing import List

import numpy as np

from .grid import CartesianGrid

#: Dəstəklənən lay bölgüsü qaydaları. Hazırda yalnız biri var; siyahı
#: AÇIQ saxlanılır ki, yeni qayda (məs. `equal` — sabit mütləq qalınlıq)
#: əlavə olunanda çağıranların yoxlaması bir yerdə qalsın.
LAYERING_RULES = ("proportional",)

#: Xəta mesajlarında sadalanan sütunların ƏN ÇOX sayı — problem 8000
#: sütunda olanda mesajı oxunmaz etməmək üçün. Ümumi say HƏMİŞƏ
#: göstərilir, yəni heç nə gizlədilmir, yalnız SİYAHI kəsilir.
_MAX_LISTED = 10


def _as_column_array(values, grid: CartesianGrid, label: str) -> np.ndarray:
    """`(nx·ny,)` areal massiv — `(ncell,)` verilsə İLK lay götürülür.

    `geology_service` xassələri həmişə `(ncell,)` kimi saxlayır (TOP/
    BOTTOM hər K-da eyni areal dəyəri daşıyır, çünki onlar areal
    nümunələrdən gəlir), ona görə hər iki forma qəbul edilir. Başqa
    ölçü SƏSSİZCƏ kəsilmir — açıq xəta atılır.
    """
    array = np.asarray(values, float).ravel()
    areal = grid.nx * grid.ny
    if array.size == areal:
        return array.copy()
    if array.size == grid.ncell:
        return array[:areal].copy()
    raise ValueError(
        f"{label}: {array.size} dəyər — gözlənilən {areal} (NX·NY) "
        f"və ya {grid.ncell} (NX·NY·NZ).")


def validate_structure(top_column, bottom_column, min_thickness: float,
                       grid: CartesianGrid = None) -> List[str]:
    """Struktur səthlərinin fiziki cəhətdən mümkünlüyünü yoxlayır.

    Aşkar edilənlər:
      · TOP və ya BOTTOM-da NaN/sonsuz (interpolyasiya alınmayıb);
      · `BOTTOM <= TOP` — səthlər KƏSİŞİR (dabanın tavandan yuxarı
        olması geoloji cəhətdən mümkünsüzdür);
      · `BOTTOM − TOP < min_thickness` — həddindən nazik lay (sıfıra
        yaxın həcm ədədi olaraq həll edilə bilməz).

    Hər problem üçün NEÇƏ sütunun təsirləndiyi və İLK {_MAX_LISTED}
    sütunun `(i, j)` indeksi göstərilir.

    HEÇ NƏ DÜZƏLDİLMİR — bu funksiya yalnız MƏSƏLƏNİ TAPIR. Nə ediləcəyi
    (xəta / clamp / deaktivasiya) çağıranın AÇIQ siyasətidir, bax
    `geology_service.GeologicalGridSpec.on_zero_thickness`.

    `grid` verilməyəndə `(i, j)` əvəzinə düz sütun indeksi göstərilir.
    """
    top = np.asarray(top_column, float).ravel()
    bottom = np.asarray(bottom_column, float).ravel()
    if top.size != bottom.size:
        return [f"TOP ({top.size}) və BOTTOM ({bottom.size}) səthlərinin "
                "ölçüsü uyğun gəlmir."]

    issues: List[str] = []
    bad_top = ~np.isfinite(top)
    bad_bottom = ~np.isfinite(bottom)
    for mask, name in ((bad_top, "TOP"), (bad_bottom, "BOTTOM")):
        if np.any(mask):
            issues.append(
                f"'{name}' səthində {int(mask.sum())} sütunda NaN/sonsuz var "
                f"— interpolyasiya bu sütunlara dəyər verməyib "
                f"({_describe(mask, grid)}).")

    finite = ~(bad_top | bad_bottom)
    if not np.any(finite):
        return issues

    thickness = np.where(finite, bottom - top, np.nan)
    crossing = finite & (thickness <= 0.0)
    if np.any(crossing):
        worst = float(np.nanmin(thickness[crossing]))
        issues.append(
            f"Səthlər KƏSİŞİR: {int(crossing.sum())} sütunda lay altı ≤ lay üstü "
            f"(ən pis fərq {worst:.3f} m) — {_describe(crossing, grid)}.")

    thin = finite & (thickness > 0.0) & (thickness < float(min_thickness))
    if np.any(thin):
        issues.append(
            f"{int(thin.sum())} sütunda qalınlıq minimumdan "
            f"({float(min_thickness):g} m) azdır (ən nazik "
            f"{float(np.nanmin(thickness[thin])):.3f} m) — {_describe(thin, grid)}.")
    return issues


def _describe(mask: np.ndarray, grid: CartesianGrid = None) -> str:
    """Problemli sütunların İLK bir neçəsini oxunaqlı sadalayır."""
    flat = np.flatnonzero(mask)
    listed = flat[:_MAX_LISTED]
    if grid is None:
        text = ", ".join(str(int(index)) for index in listed)
    else:
        text = ", ".join(f"({int(index) % grid.nx}, {int(index) // grid.nx})"
                         for index in listed)
    more = "" if flat.size <= _MAX_LISTED else f", … (+{flat.size - _MAX_LISTED})"
    return f"sütunlar: {text}{more}"


def structural_nodes(grid: CartesianGrid, top_column, bottom_column,
                     dx: float, dy: float,
                     layering: str = "proportional") -> np.ndarray:
    """`(nx·ny,)` TOP/BOTTOM səthlərindən `(ncell, 8, 3)` təpələr.

    Təpə sırası `corner_point_geometry.HEX_FACE_VERTEX_INDICES`
    konvensiyası ilə EYNİDİR: `v0..v3` = TAVAN (kiçik z), `v4..v7` =
    DABAN, hər ikisi `(i−,j−) → (i+,j−) → (i+,j+) → (i−,j+)` sırası ilə.
    Hüceyrə indeksi `(k·ny + j)·nx + i` — `CartesianGrid.index` ilə eyni.

    Sabit TOP və sabit BOTTOM verilsə nəticə `cartesian_nodes()`-un
    qaytardığı massivlə BİRƏBİR eynidir (bax
    `tests/test_structural_grid.py`) — yəni köhnə Kartezian həndəsə bu
    qurmanın XÜSUSİ HALIDIR.

    `top_column`/`bottom_column` sütun-üzrə (hüceyrə mərkəzinin XY-i
    üçün) dərinliklərdir; hüceyrənin DÖRD künc pillar-ı da EYNİ dərinliyi
    alır (bax modul docstring-i — pillarlar şaquli, künc dərinlikləri
    sütun daxilində sabit). Bu, `column_layer_edges()`-in qaytardığı
    (üzün MƏRKƏZ dərinliyi) qiymətlərin məhz interpolyasiya edilmiş
    səth qiymətləri olmasını təmin edir.
    """
    if layering not in LAYERING_RULES:
        raise ValueError(
            f"layering '{layering}' tanınmır — dəstəklənən: "
            f"{', '.join(LAYERING_RULES)}.")

    top = _as_column_array(top_column, grid, "TOP")
    bottom = _as_column_array(bottom_column, grid, "BOTTOM")

    nx, ny, nz = grid.nx, grid.ny, grid.nz
    i, j, k = grid.ijk_array(np.arange(grid.ncell))
    column = j * nx + i                                  # (ncell,) areal indeks

    # `proportional`: sütunun öz qalınlığı nz bərabər hissəyə bölünür.
    # `H` sütundan sütuna dəyişir, ona görə MÜTLƏQ lay qalınlığı da
    # dəyişir — məhz `CellGeometry.dz`-in (nz,) formasının ifadə EDƏ
    # BİLMƏDİYİ şey budur.
    thickness = bottom[column] - top[column]
    z0 = top[column] + thickness * (k / nz)
    z1 = top[column] + thickness * ((k + 1) / nz)

    x0, x1 = i * float(dx), (i + 1) * float(dx)
    y0, y1 = j * float(dy), (j + 1) * float(dy)
    corners = (
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    )
    return np.stack([np.stack(corner, axis=-1) for corner in corners], axis=1)


def structure_statistics(top_column, bottom_column,
                         grid: CartesianGrid = None) -> dict:
    """Struktur haqqında hesabata yazılan rəqəmlər (min/orta/maks
    qalınlıq, TOP dərinliyi diapazonu). NaN-a DÖZÜMLÜDÜR — natamam
    səthdə də statistika verilir, `nan_columns` isə neçə sütunun
    hesabdan kənarda qaldığını AÇIQ göstərir."""
    top = np.asarray(top_column, float).ravel()
    bottom = np.asarray(bottom_column, float).ravel()
    if grid is not None:
        top = _as_column_array(top, grid, "TOP")
        bottom = _as_column_array(bottom, grid, "BOTTOM")
    finite = np.isfinite(top) & np.isfinite(bottom)
    thickness = bottom - top
    if not np.any(finite):
        nan = float("nan")
        return {"columns": int(top.size), "nan_columns": int(top.size),
                "min_thickness": nan, "mean_thickness": nan, "max_thickness": nan,
                "min_top": nan, "max_top": nan, "relief": nan}
    return {
        "columns": int(top.size),
        "nan_columns": int((~finite).sum()),
        "min_thickness": float(thickness[finite].min()),
        "mean_thickness": float(thickness[finite].mean()),
        "max_thickness": float(thickness[finite].max()),
        "min_top": float(top[finite].min()),
        "max_top": float(top[finite].max()),
        "relief": float(top[finite].max() - top[finite].min()),
    }
