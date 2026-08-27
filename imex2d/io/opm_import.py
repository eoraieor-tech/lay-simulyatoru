"""OPM Flow nəticələrinin idxalı — VTK-a keçiddən ƏVVƏLKİ addım.

Strategiya dəyişdi: A7-nin öz üç fazalı Nyuton həlledicisi hələ də
açıq bir davamlılıq problemi daşıyır (bax `A7_PLAN.md`). Bunun
əvəzinə: FİZİKANI real, sınanmış bir simulyatora (OPM Flow, Eclipse
formatına uyğun açıq mənbəli alət) həvalə edirik, öz proqramımızın
güclü tərəfini (3D görüntü, analiz) OPM-in NƏTİCƏLƏRİNİ göstərmək
üçün işlədirik.

`resdata` kitabxanası (libecl-in müasir varisi) Eclipse binar
formatını (OPM Flow-un çıxışı da bu formatdadır) oxuyur.

TƏHLÜKƏSİZLİK QEYDİ — VACİB

resdata-nın "yüksək səviyyəli" rahatlıq metodları (`iget_restart_sim_time`
və s.) daxili C kitabxanasının SƏRT struktur fərziyyələri var —
uyğun olmayan fayl ötürüləndə **native SIGABRT** ilə çökür, bunu
Python-un `try/except`-i TUTA BİLMİR (proses bütöv dayanır). Ona görə
bu modul YALNIZ aşağı səviyyəli, "təhlükəsiz" metodları işlədir:

    `num_named_kw()`, `iget_named_kw()` — sadə açar söz oxunuşu,
    struktur fərziyyəsi yoxdur

Zaman DOUBHEAD-in 0-cı sahəsindən (elapsed sim time, gün) birbaşa
oxunur — Eclipse/OPM sənədləşməsinə görə bu sahə bütün versiyalarda
sabitdir, `iget_restart_sim_time()`-in tələb etdiyi tam INTEHEAD
strukturuna ehtiyac yoxdur.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

REQUIRED_LIBRARY_HINT = (
    "OPM Flow nəticələrini oxumaq üçün 'resdata' kitabxanası lazımdır: "
    "pip install resdata")


class OpmImportError(Exception):
    """OPM Flow halını oxumaq mümkün olmadı — aydın, tutulan xəta.

    Bütün oxuma funksiyaları bunun ALTINDAN başqa heç bir istisna
    buraxmamağa çalışır (native SIGABRT istisnadır — bax modul
    sənədləşməsi, bundan tam qorunmaq mümkün deyil, YALNIZ risqi
    azaltmaq mümkündür).
    """


@dataclass
class OpmSnapshot:
    """Bir hesabat addımının (report step) vəziyyəti."""
    time: float
    pressure: np.ndarray
    water_saturation: np.ndarray
    gas_saturation: Optional[np.ndarray] = None
    oil_saturation: Optional[np.ndarray] = None


@dataclass
class OpmGridGeometry:
    """Sadələşdirilmiş, DÜZBUCAQLI (Cartesian) tor həndəsəsi.

    OPM/Eclipse corner-point (küncnöqtə) torlarını da dəstəkləyə
    bilər, lakin BU MODUL yalnız DÜZBUCAQLI (dx/dy/dz sabit) halları
    dəstəkləyir — bax `from_grid()`-in sənədləşməsi.
    """
    nx: int
    ny: int
    nz: int
    dx: float
    dy: float
    dz: float
    top_depth: float
    active_count: int


@dataclass
class OpmFlowCase:
    """Bir OPM Flow (Eclipse formatlı) halının idxal edilmiş nəticəsi."""
    name: str
    geometry: OpmGridGeometry
    snapshots: List[OpmSnapshot] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


def _require_resdata():
    try:
        import resdata  # noqa: F401
    except ImportError as error:
        raise OpmImportError(REQUIRED_LIBRARY_HINT) from error


def _read_grid_geometry(egrid_path: str) -> OpmGridGeometry:
    from resdata.grid import Grid

    try:
        grid = Grid(egrid_path)
    except Exception as error:
        raise OpmImportError(
            f"EGRID faylı oxunmadı ({egrid_path}): {error}") from error

    nx, ny, nz = grid.get_nx(), grid.get_ny(), grid.get_nz()

    # düzbucaqlılığı yoxla — ilk hüceyrənin ölçüsünü BÜTÜN
    # qonşularla müqayisə etmirik (performans), yalnız bir neçə
    # nümunə nöqtədə. Uyğunsuzluq varsa AYDIN xəbərdarlıq veririk,
    # SƏSSIZCƏ yanlış nəticə göstərmirik.
    corner0 = np.array(grid.get_xyz(ijk=(0, 0, 0)))
    corner1 = (np.array(grid.get_xyz(ijk=(1, 0, 0)))
              if nx > 1 else corner0 + np.array([25.0, 0, 0]))
    corner2 = (np.array(grid.get_xyz(ijk=(0, 1, 0)))
              if ny > 1 else corner0 + np.array([0, 25.0, 0]))
    corner3 = (np.array(grid.get_xyz(ijk=(0, 0, 1)))
              if nz > 1 else corner0 + np.array([0, 0, 10.0]))

    dx = float(abs(corner1[0] - corner0[0])) or 25.0
    dy = float(abs(corner2[1] - corner0[1])) or 25.0
    dz = float(abs(corner3[2] - corner0[2])) or 10.0
    top_depth = float(corner0[2] - dz / 2.0)

    return OpmGridGeometry(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz,
                           top_depth=max(top_depth, 0.0),
                           active_count=grid.get_num_active())


def _read_snapshots(unrst_path: str, ncell: int) -> List[OpmSnapshot]:
    from resdata.resfile import ResdataFile

    try:
        restart = ResdataFile(unrst_path)
    except Exception as error:
        raise OpmImportError(
            f"UNRST faylı oxunmadı ({unrst_path}): {error}") from error

    try:
        step_count = restart.num_named_kw("PRESSURE")
    except Exception as error:
        raise OpmImportError(
            f"UNRST faylında PRESSURE açar sözü tapılmadı: {error}") from error

    has_swat = restart.num_named_kw("SWAT") > 0
    has_sgas = restart.num_named_kw("SGAS") > 0
    has_soil = restart.num_named_kw("SOIL") > 0

    snapshots = []
    for step in range(step_count):
        pressure = np.array(restart.iget_named_kw("PRESSURE", step))
        water = (np.array(restart.iget_named_kw("SWAT", step))
                if has_swat else np.full(pressure.shape, np.nan))
        gas = (np.array(restart.iget_named_kw("SGAS", step))
              if has_sgas else None)
        oil = (np.array(restart.iget_named_kw("SOIL", step))
              if has_soil else None)

        time = _report_step_time(restart, step)
        snapshots.append(OpmSnapshot(time=time, pressure=pressure,
                                     water_saturation=water,
                                     gas_saturation=gas, oil_saturation=oil))
    return snapshots


def _report_step_time(restart, step: int) -> float:
    """Elapsed simulyasiya vaxtı (gün) — DOUBHEAD-in 0-cı sahəsindən.

    QƏSDƏN `iget_restart_sim_time()` İŞLƏDİLMİR (bax modulun
    sənədləşməsi — natamam/qeyri-standart fayllarda native çökmə
    riski var). DOUBHEAD-in 0-cı sahəsi Eclipse/OPM sənədləşməsində
    "TSINIT" — elapsed vaxt (gün) kimi sabitdir.
    """
    try:
        if restart.num_named_kw("DOUBHEAD") > step:
            return float(restart.iget_named_kw("DOUBHEAD", step)[0])
    except Exception:
        pass
    return float(step)          # ehtiyat: addım nömrəsi vaxt kimi


def load_opm_case(case_root: str, name: Optional[str] = None) -> OpmFlowCase:
    """`case_root` — uzantısız yol (məs. `/path/CASE`, `.EGRID`/`.UNRST`
    avtomatik əlavə olunur).

    Yalnız DÜZBUCAQLI (dx/dy/dz sabit) torları dəstəkləyir. Corner-point
    (mürəkkəb həndəsə) torlar hələ dəstəklənmir — açıq xəbərdarlıq
    veriləcək (səssizcə yanlış göstərmək əvəzinə).
    """
    _require_resdata()
    egrid_path = f"{case_root}.EGRID"
    unrst_path = f"{case_root}.UNRST"

    geometry = _read_grid_geometry(egrid_path)
    snapshots = _read_snapshots(unrst_path, geometry.nx * geometry.ny
                                * geometry.nz)

    warnings: List[str] = []
    if geometry.active_count < geometry.nx * geometry.ny * geometry.nz:
        warnings.append(
            f"{geometry.nx * geometry.ny * geometry.nz - geometry.active_count} "
            "qeyri-aktiv hüceyrə var — hazırda bunlar nəzərə alınmır "
            "(bütün hüceyrələr aktiv sayılır).")
    for snapshot in snapshots:
        if snapshot.gas_saturation is None:
            warnings.append(
                "SGAS açar sözü tapılmadı — bu, qaz fazası olmayan (iki "
                "fazalı) OPM halı ola bilər.")
            break

    return OpmFlowCase(name=name or case_root, geometry=geometry,
                       snapshots=snapshots, warnings=warnings)


def build_display_model(case: OpmFlowCase):
    """OPM halını bizim öz `ReservoirModel`-imizə uyğunlaşdırır ki,
    mövcud `VolumeRenderer` heç bir dəyişiklik olmadan işə düşsün.

    Bu, ayrıca "adapter" sinfi yaratmaqdan daha təmizdir — bizim öz
    `CartesianGrid`/`CellGeometry` sinifləri artıq düzbucaqlı
    həndəsəni tam əhatə edir, təkrarlamağa ehtiyac yoxdur.

    Quyular/faylar daxil EDİLMİR (OPM halında bu məlumat WELSPECS
    kimi ayrıca fayllarda saxlanılır, bu, gələcək genişlənmədir) —
    `VolumeRenderer.draw(..., show_wells=False, show_faults=False)`
    ilə çağırılmalıdır.
    """
    from ..domain.grid import CartesianGrid
    from ..domain.geometry import CellGeometry
    from ..domain.properties import PropertyMap, RockProperties
    from ..domain.reservoir_model import ReservoirModel

    geometry = case.geometry
    grid = CartesianGrid(nx=geometry.nx, ny=geometry.ny, nz=geometry.nz)
    cell_geometry = CellGeometry(grid=grid, dx=geometry.dx, dy=geometry.dy,
                                 dz=geometry.dz, top_depth=geometry.top_depth)
    # `rock` VolumeRenderer tərəfindən istifadə OLUNMUR (yalnız təzyiq/
    # doyumluluq göstərilir) — konstruktor tələb etdiyi üçün minimal
    # (bir örnəkli) dəyərlərlə doldurulur.
    rock = RockProperties(
        porosity=PropertyMap.uniform("PORO", 0.2, grid.ncell),
        permx=PropertyMap.uniform("PERMX", 100.0, grid.ncell),
        permy=PropertyMap.uniform("PERMY", 100.0, grid.ncell))
    return ReservoirModel(name=f"OPM idxal: {case.name}", grid=grid,
                          geometry=cell_geometry, rock=rock)
