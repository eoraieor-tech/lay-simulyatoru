"""Axın diskretizasiyası — TPFA (indiki) + gələcək MPFA-O üçün ARXİTEKTURA.

Qat diaqramı (bax `ARCHITECTURE.md` — "Flow Solver" bölməsi):

    Flow Solver (FullyImplicitEngine / ImpesEngine)
        ↓
    Residual Assembly (`implicit/residual.py::ResidualAssembler`)
        ↓  yalnız `grid.connections` / `grid.pore_volume` / `grid.compute_flux()`-a güvənir
    Flux Discretization Interface (`interfaces/discretization.py::IFluxDiscretization`)
        ↓                                              ↓
    TPFA (`TwoPointFluxDiscretization`, İNDİKİ)   Future MPFA-O (HƏLƏ YOXDUR)

`ResidualAssembler` VƏ `JacobianAssembler` HANSI SİNİFDƏN `DiscretizedGrid`
gəldiyini BİLMİR — yalnız bu faylın müqaviləsinə (`connections`,
`pore_volume`, `cell_volume`, `compute_flux()`) güvənirlər. Bu, TPFA-nın
riyaziyyatını DƏYİŞMİR — mövcud nüvədən köçürülüb (harmonik orta), YALNIZ
əlavə bir METOD (`compute_flux`) vasitəsilə çağırılır ki, residual qatı
`transmissibility`-ni BİRBAŞA "T · ΔΦ" kimi işlətməsin (bax
`ResidualAssembler.face_fluxes`).

**Jacobian inteqrasiya nöqtəsi (MPFA-O üçün, HƏLƏ DƏYİŞDİRİLMƏYİB)** — bax
`implicit/jacobian.py::JacobianAssembler`:
  - `_build_pattern()` (sətir ~62-114) HƏR üz üçün DƏQİQ 2 hüceyrəli
    blok (`cell_a`↔`cell_a`, `cell_a`↔`cell_b`, və s.) fərz edir — MPFA-O
    interaction-region stensili BİR ÜZDƏ 2-dən ÇOX hüceyrəni bağlayacaq,
    ona görə bu metod GƏLƏCƏKDƏ ümumiləşdirilməli olacaq (2-hüceyrəli
    "block()" çağırışları N-hüceyrəli stensil siyahısına çevrilməli).
  - `_flux()` (sətir ~163-218) `self.R.transmissibility`-ni BİRBAŞA oxuyur
    və ∂ΔΦ/∂p_a=+1, ∂ΔΦ/∂p_b=−1 TƏK-CÜT fərziyyəsi ilə törəmə qurur —
    MPFA-O gələndə bu, `grid.compute_flux()`-un öz JACOBIAN-ını (çoxnöqtəli
    stensilin hər üzvünə görə qismən törəmə) qaytaran YENİ bir metodla
    (məs. `compute_flux_jacobian(d_phi)`) əvəz olunmalıdır.
  - Bu fayl bu iki nöqtəni SADƏCƏ SƏNƏDLƏŞDİRİR — `jacobian.py`-in özü
    BU FAZADA DƏYİŞMİR (tapşırıq: "Do not rewrite the entire Jacobian").

Konservasiya müqaviləsi (bax audit tapşırığı §8): istənilən
`IFluxDiscretization` implementasiyası üçün, hər hüceyrədə
`Σ(üz axınları) + quyu/mənbə = akkumulyasiya qalığı` LOKAL SAXLANMASI
POZULMAMALIDIR — `ResidualAssembler.net_influx`-un `np.add.at` ilə hər üz
axınını EYNİ ölçüdə (+) `cell_b`-yə, (−) `cell_a`-ya yazması (bax
`residual.py`) bunu AVTOMATİK təmin edir, bir şərtlə: `compute_flux()`
İKİ TƏRƏFDƏN EYNİ ƏDƏDİ qaytarsın (TPFA-da qaçılmazdır, MPFA-O-da da
qorunmalıdır — bax `test_flux_discretization_conserves_locally` testi).
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List

import numpy as np

from ..domain.grid import Connections
from ..domain.reservoir_model import ReservoirModel
from ..interfaces.discretization import IFluxDiscretization


@dataclass(frozen=True)
class DiscretizedGrid:
    """TPFA-nın (və gələcək MPFA-O-nun) `IFluxDiscretization.build()`-dan
    qaytardığı obyekt — bax modul docstring-i, `ResidualAssembler`-in
    güvəndiyi TAM müqavilə budur."""
    connections: Connections
    transmissibility: np.ndarray
    pore_volume: np.ndarray
    cell_volume: np.ndarray
    #: Diskretizasiya qurulan zaman yaranan, istifadəçiyə göstərilə bilən
    #: xəbərdarlıqlar (məs. tam tenzor K aşkarlanıb, TPFA onu düzgün həll
    #: EDƏ BİLMİR) — defolt boş siyahı, mövcud DAVRANIŞ DƏYİŞMİR.
    warnings: List[str] = field(default_factory=list)

    def compute_flux(self, d_phi: np.ndarray) -> np.ndarray:
        """Potensial fərqindən BAZA Darcy axını (mobilite/upstream
        çəkiləndirmə OLMADAN) — bax modul docstring-i.

        TPFA üçün bu, sadəcə `T · ΔΦ`-dir (xətti, tək-nöqtəli). Gələcək
        MPFA-O bu metodu ÖZ çoxnöqtəli stensil cəmi ilə əvəz edəcək —
        `ResidualAssembler.face_fluxes`-in özü DƏYİŞMƏYƏCƏK, çünki o
        yalnız bu metodu çağırır, `transmissibility`-ni birbaşa YOX.
        """
        return self.transmissibility * d_phi


class TwoPointFluxDiscretization(IFluxDiscretization):
    """İki nöqtəli axın approksimasiyası (TPFA) — DEFOLT diskretizasiya.

    Bax `interfaces/discretization.py::IFluxDiscretization` — gələcək
    `MPFAODiscretization` bu bazanı EYNİ ŞƏKİLDƏ tətbiq edəcək.
    """

    def build(self, model: ReservoirModel) -> DiscretizedGrid:
        conn = model.connections()
        self._assert_active_connections(model, conn)
        geom = model.geometry
        area = geom.face_areas(conn)
        half_a, half_b = geom.face_half_distances(conn)
        perm = self._directional_permeability(model, conn)

        k_a = np.maximum(perm[0], 1e-9)
        k_b = np.maximum(perm[1], 1e-9)
        trans = model.units.darcy_constant * area / (half_a / k_a + half_b / k_b)
        trans = self._apply_fault_multipliers(model, conn, trans)

        warnings = self._tensor_warnings(model)
        warnings += self._inactive_cell_warnings(model)

        return DiscretizedGrid(
            connections=conn,
            transmissibility=trans,
            pore_volume=model.pore_volume(),
            cell_volume=model.bulk_volume(),
            warnings=warnings,
        )

    @staticmethod
    def _assert_active_connections(model: ReservoirModel,
                                   conn: Connections) -> None:
        """MÜQAVİLƏ: bura çatan hər üzün HƏR İKİ tərəfi aktivdir
        (tapşırıq §3).

        Filtr `CartesianGrid.build_connections()`-dədir — üz orada
        ÜMUMİYYƏTLƏ yaradılmır, ona görə burada `T = 0` yazmağa ehtiyaq
        qalmır (sıfır transmissivlik seyrək matrisin TOPOLOGİYASINDA
        yenə qalardı və qeyri-aktiv hüceyrəni sistemə naməlum kimi
        gətirərdi — məhz bundan qaçılır).

        Bu isə HƏMİN müqavilənin POZULMADIĞINI yoxlayan ucuz (O(nconn),
        yalnız qeyri-aktiv hüceyrə olanda işləyən) qoruyucudur: kimsə
        `Connections`-ı əl ilə qurub ötürsə, səssiz yanlış nəticə əvəzinə
        AÇIQ xəta alsın.
        """
        if not model.grid.has_inactive_cells or conn.count == 0:
            return
        active = model.grid.active
        bad = ~active.active_face_mask(conn.cell_a, conn.cell_b)
        if np.any(bad):
            count = int(np.count_nonzero(bad))
            first = int(np.flatnonzero(bad)[0])
            raise ValueError(
                f"Diskretizasiyaya qeyri-aktiv hüceyrəyə toxunan {count} üz "
                f"gəldi (məs. hüceyrə {int(conn.cell_a[first])}↔"
                f"{int(conn.cell_b[first])}). Aktiv↔qeyri-aktiv bağlantı "
                f"QURULA BİLMƏZ — bax `CartesianGrid.build_connections()`.")

    @staticmethod
    def _inactive_cell_warnings(model: ReservoirModel) -> List[str]:
        """Qeyri-aktiv hüceyrələr NƏZƏRƏ ALINDIĞINI açıq bildirir.

        Əvvəllər buradakı xəbərdarlıq "ACTNUM oxundu, AMMA nəzərə
        alınmadı" idi (bax `grdecl_import._check_inactive_cells`-in köhnə
        mətni) — indi əksi doğrudur, istifadəçi bunu görməlidir."""
        grid = model.grid
        if not grid.has_inactive_cells:
            return []
        active = grid.active
        return [
            f"ACTNUM: {active.n_inactive} hüceyrə qeyri-aktivdir "
            f"({active.n_active}/{active.n_global} aktiv). Bu hüceyrələr "
            f"simulyasiyadan TAM ÇIXARILIB: PV = 0, bağlantı qurulmur, "
            f"xətti sistemdə naməlumları yoxdur."]

    @staticmethod
    def _tensor_warnings(model: ReservoirModel) -> List[str]:
        """Tam tenzor K aşkarlansa AÇIQ, DETERMİNİSTİK xəbərdarlıq verir
        — TPFA onu SƏSSİZCƏ diaqonala "yumşaltmır" (scalarize etmir),
        sadəcə İSTİFADƏ ETMİR və bunu, HANSI komponentlərin (Kxy/Kxz/Kyz)
        və neçə hüceyrənin təsirləndiyini AÇIQ göstərərək bildirir (bax
        audit Phase 2 §11: "the existing warning mechanism can be
        improved, make it clear and deterministic").
        """
        tensor = model.rock.permeability_tensor
        if tensor is None or not tensor.has_off_diagonal():
            return []
        tol = 1e-12
        affected: List[str] = []
        for label, component in (("Kxy", tensor.kxy), ("Kxz", tensor.kxz),
                                 ("Kyz", tensor.kyz)):
            if component is not None and np.any(np.abs(component.values) > tol):
                n_nonzero = int(np.sum(np.abs(component.values) > tol))
                affected.append(f"{label} ({n_nonzero} hüceyrə)")
        return [
            "Tam permeabilite tenzoru aşkarlanıb, sıfırdan fərqli off-diaqonal "
            f"komponent(lər): {', '.join(affected)}. TPFA (TwoPointFluxDiscretization) "
            "YALNIZ diaqonal Kx/Ky/Kz-i işlədir — off-diaqonal anizotropluq bu "
            "simulyasiyada NƏZƏRƏ ALINMIR (səssizcə diaqonala yumşaldılmır, "
            "sadəcə istifadə edilmir). Düzgün nəticə üçün MPFA-O lazımdır "
            "(hələ implement edilməyib)."]

    @staticmethod
    def _directional_permeability(model: ReservoirModel, conn: Connections):
        kx = model.rock.permx.values
        ky = model.rock.permy.values
        kz = model.rock.permz.values if model.rock.permz is not None else kx
        per_axis = (kx, ky, kz)
        k_a = np.empty(conn.count)
        k_b = np.empty(conn.count)
        for axis in (0, 1, 2):
            m = conn.axis == axis
            if np.any(m):
                k_a[m] = per_axis[axis][conn.cell_a[m]]
                k_b[m] = per_axis[axis][conn.cell_b[m]]
        return k_a, k_b

    @staticmethod
    def _apply_fault_multipliers(model: ReservoirModel, conn: Connections,
                                 trans: np.ndarray) -> np.ndarray:
        """Fay transmissivlik çarpanları (B3).

        Hər fay bir müstəvidir: yalnız həmin müstəvidəki üzlərin
        transmissivliyi çarpanla dəyişir, qalan bütün üzlər toxunulmaz
        qalır. Bir üz birdən çox faya düşərsə, çarpanlar VURULUR —
        iki qismən keçirici fay eyni yerdə üst-üstə düşəndə axının
        daha da azalması fiziki cəhətdən doğrudur.
        """
        faults = [f for f in model.fault_references if f.has_geometry]
        if not faults:
            return trans

        grid = model.grid
        i_a, j_a, k_a = grid.ijk_array(conn.cell_a)
        coordinates = {0: (j_a, k_a), 1: (i_a, k_a), 2: (i_a, j_a)}
        axis_code = {"I": 0, "J": 1, "K": 2}
        boundary = {0: i_a, 1: j_a, 2: k_a}

        trans = trans.copy()
        for fault in faults:
            code = axis_code[fault.axis.upper()]
            coordinate_a, coordinate_b = coordinates[code]
            mask = ((conn.axis == code)
                    & (boundary[code] == fault.plane_index)
                    & fault.matches(code, coordinate_a, coordinate_b))
            if np.any(mask):
                trans[mask] *= fault.effective_multiplier
        return trans


def default_flux_discretization() -> IFluxDiscretization:
    """DEFOLT axın diskretizasiyası — bax audit tapşırığı §4: "default
    discretization = TPFA". `FullyImplicitEngine`/`ImpesEngine`
    konstruktoruna `flux_discretization=None` veriləndə BUNDAN istifadə
    edilir. Gələcəkdə `MPFAODiscretization` əlavə olunanda, defolt
    DƏYİŞMİR — istifadəçi AÇIQ seçim etməli olacaq."""
    return TwoPointFluxDiscretization()
