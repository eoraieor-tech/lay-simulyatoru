"""VTK əsaslı 3D rezervuar görüntüsü — ResInsight tipli.

NİYƏ VTK

Mövcud `volume.py` matplotlib-in `Poly3DCollection`-ına əsaslanır.
Matplotlib 3D üçün nəzərdə tutulmayıb — o, 2D çəkici olub 3D-ni
"təqlid" edir:

    · dərinlik sıralaması dəqiq deyil (z-buffer yoxdur, mərkəz
      nöqtəsinə görə təxmini sıralama) — bax `volume.py`-dəki fault
      görüntüsündə buna görə dolğu əvəzinə kontur işlədilib
    · hər fırlatmada BÜTÜN üzlər yenidən sıralanır və çəkilir —
      böyük gridlərdə (100k+ hüceyrə) dözülməz ləng
    · işıqlandırma primitivdir (əl ilə kölgə hesablanır)

VTK isə əsl 3D qrafika kitabxanasıdır (OpenGL üzərində) — ResInsight,
ParaView, Petrel kimi peşəkar alətlər onu işlədir. Rezervuar
görüntüsündə istifadəçinin gözlədiyi hər şey (rəvan fırlatma, real
işıqlandırma, kəsim müstəviləri, böyük gridlərdə sürət) buradan gəlir.

GERİYƏ UYĞUNLUQ

Bu modul mövcud `VolumeRenderer`-i ƏVƏZ ETMİR — onun yanında yaşayır.
İstifadəçi interfeysdə hansı motoru işlədəcəyini seçir; VTK
quraşdırılmayıbsa, `is_available()` False qaytarır və proqram
matplotlib motoruna qayıdır (bax `available()` funksiyası).
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


def available() -> bool:
    """VTK quraşdırılıbmı — proqramın qalan hissəsi buna görə qərar verir.

    VTK məcburi asılılıq DEYİL: onsuz proqram tam işləyir (matplotlib
    motoru ilə). Bu ayrılıq sayəsində `pip install vtk` etməyən
    istifadəçi heç nə itirmir.
    """
    try:
        import vtk  # noqa: F401
    except ImportError:
        return False
    return True


@dataclass
class VtkViewSettings:
    """Görüntü parametrləri — matplotlib motoru ilə EYNİ məna daşıyır.

    Eyni adlar qəsdən saxlanılıb ki, interfeys eyni idarəetmə
    elementlərini hər iki motora ötürə bilsin.
    """
    colormap: object = "viridis"
    """Ad (str) və ya matplotlib `Colormap` obyekti — bax `_build_lookup_table`."""
    value_limits: Optional[Tuple[float, float]] = None
    k_range: Optional[Tuple[int, int]] = None
    value_min: Optional[float] = None
    """Kəsim həddi — bu dəyərdən aşağı hüceyrələr gizlədilir."""
    cell_mask: Optional[object] = None
    """Hazır (ncell,) bool maska — `volume.VolumeFilter.cell_mask` ilə EYNİ
    məna (məs. status filtri). Bu qat onu HESABLAMIR, yalnız tətbiq edir."""
    show_edges: bool = True
    opacity: float = 1.0
    vertical_exaggeration: float = 1.0
    show_wells: bool = True
    show_faults: bool = True
    shading: float = 0.45
    """İşıq gücü (0..1) — matplotlib motorundakı "İşıq" sürgüsü ilə eyni."""
    zoom: float = 1.0
    """Yaxınlaşdırma (1.0 = tam model çərçivədə)."""


# ── rəng xəritələri ──────────────────────────────────────────────────
# ResInsight/Petrel-də standart olan rəng ardıcıllıqları. matplotlib
# adları ilə eyni açarlar işlədilir ki, interfeys dəyişməsin.
_COLORMAPS = {
    "viridis": [(0.267, 0.005, 0.329), (0.128, 0.567, 0.551),
                (0.993, 0.906, 0.144)],
    "plasma": [(0.050, 0.030, 0.528), (0.798, 0.280, 0.470),
               (0.940, 0.975, 0.131)],
    "coolwarm": [(0.230, 0.299, 0.754), (0.865, 0.865, 0.865),
                 (0.706, 0.016, 0.150)],
    "jet": [(0.0, 0.0, 0.5), (0.0, 1.0, 1.0), (1.0, 1.0, 0.0),
            (0.5, 0.0, 0.0)],
    "turbo": [(0.190, 0.072, 0.232), (0.129, 0.756, 0.769),
              (0.994, 0.907, 0.144), (0.479, 0.010, 0.013)],
}


def _build_lookup_table(colormap, low: float, high: float):
    """Rəng xəritəsi — VTK-nın öz `vtkLookupTable`-i.

    `colormap` HƏM ad (str), HƏM matplotlib `Colormap` obyekti ola
    bilər: interfeys `MapRenderer._select_volume()`-dən OBYEKT alır,
    lakin bu modul öz-özlüyündə matplotlib-dən asılı olmamalıdır.
    Obyekt gələndə matplotlib-in ÖZ rənglərini birbaşa nümunələyirik
    — bu, iki motor arasında rəng fərqini tam aradan qaldırır
    (`_COLORMAPS` yalnız ad gələndə, ehtiyat kimi işlədilir).
    """
    import vtk

    sampler = getattr(colormap, "__call__", None)
    if sampler is not None and not isinstance(colormap, str):
        table = vtk.vtkLookupTable()
        table.SetNumberOfTableValues(256)
        table.SetRange(low, high)
        for index in range(256):
            red, green, blue, _ = colormap(index / 255.0)
            table.SetTableValue(index, red, green, blue, 1.0)
        table.Build()
        return table

    name = colormap if isinstance(colormap, str) else "viridis"
    anchors = _COLORMAPS.get(name, _COLORMAPS["viridis"])
    table = vtk.vtkLookupTable()
    table.SetNumberOfTableValues(256)
    table.SetRange(low, high)

    segments = len(anchors) - 1
    for index in range(256):
        position = index / 255.0 * segments
        segment = min(int(position), segments - 1)
        weight = position - segment
        start = anchors[segment]
        end = anchors[segment + 1]
        colour = [start[c] + (end[c] - start[c]) * weight for c in range(3)]
        table.SetTableValue(index, colour[0], colour[1], colour[2], 1.0)
    table.Build()
    return table


class VtkReservoirScene:
    """Rezervuar modelinin VTK səhnəsi.

    Səhnə BİR DƏFƏ qurulur, sonra yalnız DƏYƏRLƏR yenilənir
    (`update_values()`). Bu, matplotlib motorundan əsas fərqdir —
    orada hər yeniləmədə bütün həndəsə yenidən qurulurdu. Zaman
    slider-ini sürüşdürəndə fərq dərhal hiss olunur.
    """

    def __init__(self, model, settings: Optional[VtkViewSettings] = None):
        import vtk

        self.model = model
        self.settings = settings or VtkViewSettings()
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(0.078, 0.102, 0.122)   # PALETTE.background
        self.renderer.SetBackground2(0.110, 0.145, 0.173)
        self.renderer.GradientBackgroundOn()

        self._grid = None
        self._mapper = None
        self._actor = None
        self._scalar_bar = None
        self._well_actors = []
        self._fault_actors = []
        self._axes_actor = None
        self._orientation_widget = None
        self._build_grid()

    # ── həndəsə ────────────────────────────────────────────────────
    def _cell_corner_points(self):
        """KARTEZİAN grid düyün koordinatları — (nx+1)·(ny+1)·(nz+1).

        YALNIZ Kartezian yol. Model corner-point həndəsəsi daşıyırsa bu
        metod ÇAĞIRILMIR — orada `x = i·dx`, `y = j·dy` uydurma olardı;
        `_corner_point_cell_nodes()` deck-in ÖZ təpələrini işlədir
        (bax `_build_grid`).

        `vtkStructuredGrid` DÜYÜNLƏRLƏ işləyir, hüceyrə mərkəzləri ilə
        yox. Hüceyrə mərkəz dərinliklərindən (`cell_depths()`) düyün
        dərinliklərini çıxarırıq: hər düyün ona toxunan hüceyrələrin
        ORTALAMASI kimi hesablanır — bu, maili/qırışıqlı laylarda
        (`top_depth_map`) hamar səth verir.

        Vektorlaşdırılıb: əvvəlki (üç iç-içə döngəli) versiya böyük
        gridlərdə çox ləng idi VƏ kənar düyünləri səhv yerləşdirirdi
        (`min(i, nx-1)` sıxılması hüceyrə eninə bərabər çıxıntı
        yaradırdı — ilk render sınağında görüldü).
        """
        grid = self.model.grid
        geometry = self.model.geometry
        nx, ny, nz = grid.nx, grid.ny, grid.nz
        exaggeration = max(self.settings.vertical_exaggeration, 1e-6)

        centres = geometry.cell_depths().reshape((nz, ny, nx))

        # hüceyrə mərkəzlərini areal olaraq düyün şəbəkəsinə genişləndir:
        # kənarda ekstrapolyasiya YOX, sadəcə ən yaxın hüceyrənin dəyəri
        # (rezervuar kənarı düz kəsilir — fiziki cəhətdən doğru)
        areal = np.empty((nz, ny + 1, nx + 1))
        padded = np.pad(centres, ((0, 0), (0, 1), (0, 1)), mode="edge")
        areal[:, :ny + 1, :nx + 1] = padded

        # düyün səviyyələri: hüceyrə mərkəzindən ±dz/2 (hər təbəqənin
        # öz qalınlığı ilə — dz artıq (nz,) massivdir)
        half_dz = geometry.dz / 2.0
        levels = np.empty((nz + 1, ny + 1, nx + 1))
        levels[:nz] = areal - half_dz[:, None, None]
        levels[nz] = areal[nz - 1] + half_dz[-1]

        x = np.arange(nx + 1) * geometry.dx
        y = np.arange(ny + 1) * geometry.dy
        xx, yy = np.meshgrid(x, y)                    # (ny+1, nx+1)

        points = np.empty(((nz + 1) * (ny + 1) * (nx + 1), 3))
        points[:, 0] = np.tile(xx.ravel(), nz + 1)
        points[:, 1] = np.tile(yy.ravel(), nz + 1)
        points[:, 2] = -levels.ravel() * exaggeration
        return points

    # ── HƏQİQİ corner-point həndəsəsi ──────────────────────────────
    def uses_corner_point_geometry(self) -> bool:
        """Model HƏQİQİ corner-point həndəsəsi daşıyırmı.

        Yoxlama TİPƏ görə deyil, MÜQAVİLƏYƏ görədir (`nodes` (ncell,8,3)):
        bu qat `domain`-dən yalnız həndəsə oxuyur, ona görə hansı sinifin
        həmin təpələri verdiyi onun işi deyil.
        """
        nodes = getattr(self.model.geometry, "nodes", None)
        if nodes is None:
            return False
        nodes = np.asarray(nodes)
        return nodes.ndim == 3 and nodes.shape == (self.model.grid.ncell, 8, 3)

    def _corner_point_cell_nodes(self) -> np.ndarray:
        """`(ncell, 8, 3)` — hər hüceyrənin HƏQİQİ 8 təpəsi, VTK
        koordinat sistemində.

            COORD + ZCORN  →  8 təpə / hüceyrə  →  VTK həndəsəsi

        Domendə Z DƏRİNLİKDİR (aşağı müsbət), VTK-da isə Z YUXARI
        müsbətdir — ona görə işarə çevrilir və şaquli şişirtmə həmin
        addımda tətbiq olunur (X/Y TOXUNULMUR, bax
        `test_vertical_exaggeration_scales_depth_only`).

        Təpə SIRASI da çevrilir: domen konvensiyası (bax
        `polyhedral_geometry.HEX_FACE_VERTEX_INDICES`) `0..3`-ü DAHA AZ
        DƏRİNLİKDƏ (yəni VTK-da DAHA YUXARIDA) saxlayır, `VTK_HEXAHEDRON`
        isə ƏVVƏLCƏ aşağı üzü gözləyir. Sıra tərs verilsəydi hüceyrənin
        yakobianı mənfi çıxar, üz normalları içəri baxar və işıqlandırma
        tərsinə düşərdi.
        """
        nodes = np.asarray(self.model.geometry.nodes, float)
        exaggeration = max(self.settings.vertical_exaggeration, 1e-6)
        vtk_nodes = np.empty_like(nodes)
        vtk_nodes[:, 0:4] = nodes[:, 4:8]        # daban (dərin) → VTK-da aşağı
        vtk_nodes[:, 4:8] = nodes[:, 0:4]        # tavan (dayaz) → VTK-da yuxarı
        vtk_nodes[..., 2] = -vtk_nodes[..., 2] * exaggeration
        return vtk_nodes

    def _build_grid(self):
        if self.uses_corner_point_geometry():
            self._grid = self._build_corner_point_grid()
        else:
            self._grid = self._build_structured_grid()

    def _build_corner_point_grid(self):
        """`vtkUnstructuredGrid` — hər hüceyrə ÖZ 8 təpəsi ilə.

        NİYƏ STRUKTUR ŞƏBƏKƏ DEYİL: `vtkStructuredGrid` düyünləri
        QONŞULAR ARASINDA PAYLAŞIR, yəni bir düyünə yalnız BİR koordinat
        düşür. Corner-point grid isə məhz bunu POZUR — fay atımında,
        pinch-out-da və qeyri-konformal mesh-də eyni pillar düyünü qonşu
        hüceyrələrdə FƏRQLİ dərinlikdədir. Paylaşılan düyünlə belə
        həndəsəni göstərmək üçün onu ortalamaq lazım gəlir və məhz bu
        ortalama əvvəlki versiyada həqiqi həndəsəni itirirdi.

        Təpələr QƏSDƏN paylaşılmır (hüceyrə başına 8 nöqtə): yaddaş
        artımı kiçikdir, əvəzində hər hüceyrə deck-dəki KOORDİNATI
        HƏRFİ OLARAQ saxlayır.

        Hüceyrə SIRASI dəyişmir (`i` ən sürətli, sonra `j`, sonra `k`) —
        `update_values()` hüceyrə skalyarlarını elə həmin sıra ilə
        bağlayır.
        """
        import vtk
        from vtkmodules.util import numpy_support

        ncell = self.model.grid.ncell
        coordinates = self._corner_point_cell_nodes().reshape(ncell * 8, 3)

        points = vtk.vtkPoints()
        points.SetData(numpy_support.numpy_to_vtk(
            np.ascontiguousarray(coordinates), deep=True))

        # VTK hüceyrə massivi — `offsets` + `connectivity` cütü (VTK 9-un
        # öz daxili təmsili; köhnə `SetCells()` 9.6-dan etibarən
        # köhnəlib). Vektorlaşdırılıb: 100k+ hüceyrədə Python dövrü
        # nəzərəçarpan gecikmə verirdi.
        offsets = np.arange(ncell + 1, dtype=np.int64) * 8
        connectivity = np.arange(ncell * 8, dtype=np.int64)
        cells = vtk.vtkCellArray()
        cells.SetData(numpy_support.numpy_to_vtkIdTypeArray(offsets, deep=True),
                      numpy_support.numpy_to_vtkIdTypeArray(connectivity, deep=True))

        unstructured = vtk.vtkUnstructuredGrid()
        unstructured.SetPoints(points)
        unstructured.SetCells(vtk.VTK_HEXAHEDRON, cells)
        return unstructured

    def _build_structured_grid(self):
        """Kartezian model — DƏYİŞMƏYİB (paylaşılan düyünlü struktur
        şəbəkə: bərabər bloklarda o, həm daha yığcamdır, həm də həqiqi
        həndəsəni tam təsvir edir)."""
        import vtk
        from vtkmodules.util import numpy_support

        grid = self.model.grid
        nx, ny, nz = grid.nx, grid.ny, grid.nz

        structured = vtk.vtkStructuredGrid()
        structured.SetDimensions(nx + 1, ny + 1, nz + 1)

        points = vtk.vtkPoints()
        coordinates = self._cell_corner_points()
        points.SetData(numpy_support.numpy_to_vtk(coordinates, deep=True))
        structured.SetPoints(points)
        return structured

    # ── dəyərlər ───────────────────────────────────────────────────
    def update_values(self, values: np.ndarray, label: str = ""):
        """Hüceyrə dəyərlərini yeniləyir — həndəsə TOXUNULMUR."""
        import vtk
        from vtkmodules.util import numpy_support

        values = np.asarray(values, float).ravel()
        settings = self.settings

        array = numpy_support.numpy_to_vtk(values, deep=True)
        array.SetName(label or "dəyər")
        self._grid.GetCellData().SetScalars(array)

        visible = self._visibility_mask(values)
        self._apply_visibility(visible)

        low, high = settings.value_limits or self._auto_limits(values, visible)
        table = _build_lookup_table(settings.colormap, low, high)

        if self._mapper is None:
            self._mapper = vtk.vtkDataSetMapper()
            self._actor = vtk.vtkActor()
            self._actor.SetMapper(self._mapper)
            self.renderer.AddActor(self._actor)

        self._mapper.SetInputData(self._grid)
        self._mapper.SetLookupTable(table)
        self._mapper.SetScalarRange(low, high)
        self._mapper.SetScalarModeToUseCellData()
        self._mapper.ScalarVisibilityOn()

        properties = self._actor.GetProperty()
        properties.SetOpacity(settings.opacity)
        # İşıq sürgüsü: matplotlib-də əl ilə kölgə hesablanırdı, VTK-da
        # bu, materialın əsl işıq xassələridir (diffuse/ambient nisbəti).
        shading = min(max(settings.shading, 0.0), 1.0)
        properties.SetAmbient(1.0 - shading * 0.75)
        properties.SetDiffuse(0.25 + shading * 0.75)
        properties.SetSpecular(shading * 0.12)
        if settings.show_edges:
            properties.EdgeVisibilityOn()
            properties.SetEdgeColor(0.15, 0.18, 0.21)
            properties.SetLineWidth(0.5)
        else:
            properties.EdgeVisibilityOff()

        self._update_scalar_bar(table, label)
        self.update_wells()
        self.update_faults()
        self.update_axes()

    def _visibility_mask(self, values: np.ndarray) -> np.ndarray:
        """Kəsim həddi + K-təbəqə filtri — matplotlib motoru ilə eyni məntiq."""
        grid = self.model.grid
        visible = np.ones(values.size, dtype=bool)

        if self.settings.value_min is not None:
            visible &= values >= self.settings.value_min

        if self.settings.k_range is not None:
            k_from, k_to = self.settings.k_range
            layer = np.repeat(np.arange(grid.nz), grid.nx * grid.ny)
            visible &= (layer >= k_from) & (layer <= k_to)

        if self.settings.cell_mask is not None:
            visible &= np.asarray(self.settings.cell_mask, bool).ravel()
        return visible

    def _apply_visibility(self, visible: np.ndarray):
        """Gizli hüceyrələri VTK-nın `BlankCell` mexanizmi ilə söndürür.

        Hüceyrələri SİLMƏK əvəzinə "blank" edirik — həndəsə toxunulmaz
        qalır, ona görə kəsim həddini sürüşdürəndə yenidən qurma
        aparılmır (matplotlib motorunda hər dəfə yenidən qurulurdu).
        """
        import vtk

        self._grid.AllocateCellGhostArray()
        ghosts = self._grid.GetCellGhostArray()
        hidden = vtk.vtkDataSetAttributes.HIDDENCELL
        for index, is_visible in enumerate(visible):
            ghosts.SetValue(index, 0 if is_visible else hidden)
        ghosts.Modified()

    @staticmethod
    def _auto_limits(values: np.ndarray, visible: np.ndarray):
        shown = values[visible]
        if shown.size == 0:
            return 0.0, 1.0
        low, high = float(shown.min()), float(shown.max())
        return (low, high) if high > low else (low, low + 1.0)

    # ═══════════════════════════════════════════ koordinat şəbəkəsi
    def update_axes(self):
        """Modelin ətrafında ölçü şəbəkəsi — ResInsight-dakı kimi.

        `vtkCubeAxesActor` X/Y/Z oxlarını rəqəmlərlə çəkir. Bu, modelə
        MİQYAS HİSSİ verir: onsuz model "havada asılı" görünür və
        istifadəçi ölçüləri qiymətləndirə bilmir.

        Dərinlik oxu (Z) MƏNFİ işarə ilə saxlanılır (bizdə dərinlik
        aşağı müsbətdir, VTK-da isə Z yuxarı müsbətdir) — etiketlərdə
        istifadəçi ƏSL dərinliyi (müsbət, metr) görməlidir, ona görə
        Z aralığı işarə dəyişdirilərək verilir.
        """
        import vtk

        if self._axes_actor is not None:
            self.renderer.RemoveViewProp(self._axes_actor)

        bounds = self._grid.GetBounds()
        axes = vtk.vtkCubeAxesActor()
        axes.SetBounds(bounds)
        axes.SetCamera(self.renderer.GetActiveCamera())

        exaggeration = max(self.settings.vertical_exaggeration, 1e-6)
        top_depth = -bounds[5] / exaggeration
        base_depth = -bounds[4] / exaggeration

        # DƏRİNLİK ETİKETLƏRİ — MÜTLƏQ yox, NİSBİ göstərilir.
        #
        # Mütləq dərinlik (2000, 2004, 2008 … m) tipik rezervuarda
        # 4 rəqəmlidir, aralıq isə dardır (onlarla metr). Nəticədə
        # etiketlər ÜST-ÜSTƏ DÜŞÜR — iki fərqli üsul sınandı və hər
        # ikisi kifayət etmədi: kiçik font (etiketlər hələ sıx) və
        # `SetScreenSize` (VTK etiket SAYINA təsir etmir).
        #
        # Nisbi dərinlik (0, 10, 20 … m — tavandan aşağı) rəqəmləri
        # 1-2 rəqəmə endirir və həmişə oxunur. Mütləq dərinlik oxun
        # BAŞLIĞINDA göstərilir ki, məlumat itməsin.
        axes.SetZAxisRange(0.0, base_depth - top_depth)

        # BAŞLIQLAR QƏSDƏN BOŞDUR — yalnız rəqəmlər göstərilir.
        #
        # İki səbəb:
        #  1. Sağ aşağıdakı istiqamət oxu (X/Y/Z) onsuz da hansı oxun
        #     hansı olduğunu göstərir — başlıq təkrardır
        #  2. VTK-nın defolt fontu Azərbaycan hərflərini (ə, ı, ü)
        #     DƏSTƏKLƏMİR: "Dərinlik" ekranda "Dinlik" kimi görünürdü
        #     (istifadəçi şəkildə göstərdi). Başlıqsız bu problem
        #     tamamilə aradan qalxır.
        # BOŞ SƏTİR YOX, BOŞLUQ SİMVOLU (" ") verilir.
        #
        # `SetXTitle("")` — VTK-nın daxili `vtkVectorText` filtri boş
        # mətni QƏBUL ETMİR: hər kadrda "Text is not set!" xətası atır
        # və Windows-da fasiləsiz `vtkOutputWindow` pəncərəsi açılır
        # (istifadəçi şəkildə göstərdi). Boşluq simvolu keçərli mətndir
        # — xəta yoxdur, ekranda isə heç nə görünmür.
        axes.SetXTitle(" ")
        axes.SetYTitle(" ")
        axes.SetZTitle(" ")

        # Etiket formatı — VTK-nın defoltu ("%-#6.3g") dərinlik üçün
        # YARARSIZDIR: 2000–2032 m aralığı "2.0 2.0 2.0" kimi üst-üstə
        # düşən etiketlərə çevrilirdi (ilk render sınağında göründü),
        # çünki 3 əhəmiyyətli rəqəm bu aralığı ayırd edə bilmir.
        # Tam ədəd formatı dərinlik üçün həmişə oxunaqlıdır.
        axes.SetZLabelFormat("%.0f")
        axes.SetXLabelFormat("%.0f")
        axes.SetYLabelFormat("%.0f")
        # `SetLabelScaling(False, ...)` — VTK-nın AVTOMATİK üstlü
        # miqyaslamasını (etiketlərin "×10³" ilə sıxılması) söndürür.
        # O olmadan format ayarı təsirsizdir: dərinlik "2 2 2 2" kimi
        # görünürdü, çünki 2000–2032 aralığı ×10³-ə bölünəndə fərqlər
        # itirdi.
        axes.SetLabelScaling(False, 0, 0, 0)

        grid_colour = (0.35, 0.42, 0.48)
        text_colour = (0.72, 0.78, 0.82)
        for index in range(3):
            title = axes.GetTitleTextProperty(index)
            label = axes.GetLabelTextProperty(index)
            title.SetColor(*text_colour)
            label.SetColor(*text_colour)
            title.ShadowOff()
            label.ShadowOff()
            title.ItalicOff()
            label.ItalicOff()
            # Kiçik font — dərinlik oxu qısa olduğu üçün etiketlər
            # üst-üstə düşürdü (ilk sınaqda göründü). VTK etiket
            # SAYINI birbaşa idarə etmir, ona görə ölçü ilə həll olunur.
            # ETİKET ÖLÇÜSÜ SABİT EKRAN PİKSELİNDƏ saxlanılır.
            # VTK-nın defoltu mətni MODEL ölçüsünə görə miqyaslayır —
            # nəticədə rəqəmlər nəhəng olub bütün ekranı tuturdu
            # (istifadəçi şəkildə göstərdi). `SetFontSize` yalnız
            # `SetScreenSize`-la birlikdə işləyir.
            label.SetFontSize(10)
            title.SetFontSize(10)

        # kiçik bölgülər sıxlığı artırır — söndürülür
        axes.XAxisMinorTickVisibilityOff()
        axes.YAxisMinorTickVisibilityOff()
        axes.ZAxisMinorTickVisibilityOff()


        # ETİKET SIXLIĞI — VTK etiket SAYINI birbaşa idarə etmir
        # (`vtkCubeAxesActor`-da belə metod yoxdur, yoxlanılıb). Onun
        # əvəzinə `SetScreenSize` işlədilir: bu, etiketlərin EKRANDA
        # tutduğu ölçünü təyin edir və VTK ona uyğun sıxlığı özü
        # seçir. Defolt 10.0 dərinlik oxu üçün çox sıxdır — ölçülüb:
        # 2000–2032 m aralığında etiketlər tam üst-üstə düşürdü.
        #
        # Dərinlik aralığı areal ölçüdən nə qədər dardırsa, etiketlər
        # bir o qədər seyrək olmalıdır.
        areal_span = max(bounds[1] - bounds[0], bounds[3] - bounds[2], 1.0)
        depth_span = max(bounds[5] - bounds[4], 1e-6)
        ratio = depth_span / areal_span
        # `SetScreenSize` — etiketin EKRANDA tutduğu ölçü (piksel).
        # Kiçik dəyər həm mətn ölçüsünü, həm sıxlığı azaldır.
        # Ölçülüb: defolt 10.0 belə çox böyükdür, çünki VTK bunu
        # model ölçüsü ilə birlikdə miqyaslayır.
        axes.SetScreenSize(9.0 if ratio >= 0.15 else 11.0)

        # NAZİK modeldə (qalınlıq areal ölçünün kiçik hissəsi) dərinlik
        # oxu ekranda o qədər qısa olur ki, etiketlər bir-birinin
        # üstünə yığılır və oxunmur (istifadəçi şəkildə göstərdi).
        # Belə halda oxun ETİKETLƏRİ gizlədilir — areal oxlar (X/Y)
        # onsuz da miqyası verir, dərinlik isə rəng legendindən və
        # modelin öz formasından anlaşılır.
        #
        # Hədd ölçülüb: 41×41×5 model (820 m areal, 50 m qalınlıq,
        # nisbət 0.06) hələ sıx görünürdü; 0.12 real modellərdə
        # etibarlı işləyir.
        if ratio < 0.12:
            axes.ZAxisLabelVisibilityOff()
            axes.ZAxisTickVisibilityOff()
        for prop in (axes.GetXAxesLinesProperty(), axes.GetYAxesLinesProperty(),
                     axes.GetZAxesLinesProperty()):
            prop.SetColor(*grid_colour)
        for prop in (axes.GetXAxesGridlinesProperty(),
                     axes.GetYAxesGridlinesProperty(),
                     axes.GetZAxesGridlinesProperty()):
            prop.SetColor(*grid_colour)
            prop.SetOpacity(0.35)

        # ResInsight-dakı kimi: şəbəkə xətləri XARİCİ divarlarda
        axes.DrawXGridlinesOn()
        axes.DrawYGridlinesOn()
        axes.DrawZGridlinesOn()
        axes.SetGridLineLocation(axes.VTK_GRID_LINES_FURTHEST)
        axes.SetFlyModeToOuterEdges()

        self._axes_actor = axes
        self.renderer.AddViewProp(axes)

    # ══════════════════════════════════════════════ istiqamət oxu
    def attach_orientation_marker(self, interactor):
        """Sağ aşağı küncdə X/Y/Z istiqamət oxu — model ilə fırlanır.

        `vtkOrientationMarkerWidget` İNTERAKTOR tələb edir (adi aktyor
        deyil, öz kiçik görüntü sahəsi var) — ona görə səhnə qurulanda
        yox, interfeys interaktoru hazırlayandan SONRA çağırılır.
        Offscreen render zamanı (test/şəkil saxlama) interaktor
        olmadığı üçün bu, sadəcə atlanılır.
        """
        import vtk

        if interactor is None or self._orientation_widget is not None:
            return

        axes = vtk.vtkAxesActor()
        axes.SetXAxisLabelText("X")
        axes.SetYAxisLabelText("Y")
        axes.SetZAxisLabelText("Z")
        axes.SetTotalLength(1.0, 1.0, 1.0)
        for caption in (axes.GetXAxisCaptionActor2D(),
                       axes.GetYAxisCaptionActor2D(),
                       axes.GetZAxisCaptionActor2D()):
            caption.GetCaptionTextProperty().SetColor(0.86, 0.89, 0.92)
            caption.GetCaptionTextProperty().ShadowOff()
            caption.GetCaptionTextProperty().ItalicOff()

        widget = vtk.vtkOrientationMarkerWidget()
        widget.SetOrientationMarker(axes)
        widget.SetInteractor(interactor)
        widget.SetViewport(0.82, 0.02, 0.99, 0.22)   # sağ aşağı künc
        widget.SetEnabled(1)
        widget.InteractiveOff()          # istifadəçi sürüşdürə bilməsin
        self._orientation_widget = widget

    # ═══════════════════════════════════════════════════════ quyular
    def update_wells(self):
        """Quyu lülələri + perforasiya nöqtələri + adlar.

        matplotlib motorundakı `_draw_wells()` ilə EYNİ məntiq: lülə
        səthdən ən dərin perforasiyaya qədər, INJ mavi / PROD narıncı,
        ad lülənin başında.

        VTK-da lülə əsl SİLİNDR-dir (matplotlib-də sadəcə qalın xətt)
        — dərinlik hissi düzgün verilir, model fırlananda lülə də
        həqiqi 3D obyekt kimi davranır.
        """
        from ..domain.wells import WellType

        for actor in self._well_actors:
            self.renderer.RemoveViewProp(actor)
        self._well_actors = []

        if not self.settings.show_wells:
            return

        model = self.model
        geometry = model.geometry
        exaggeration = max(self.settings.vertical_exaggeration, 1e-6)
        depths = geometry.cell_depths().reshape(model.grid.shape)
        # Quyu lüləsi hüceyrənin HƏQİQİ mərkəzindən keçir. Əvvəl X/Y
        # `(i+0.5)·dx` kimi hesablanırdı — corner-point modeldə bu,
        # nominal (ortalama) ölçüdür, ona görə maili/əyri sütunda lülə
        # perforasiya etdiyi hüceyrələrin YANINDA qalırdı. Kartezian
        # modeldə `cell_centroid()` məhz `(i+0.5)·dx` verir, yəni
        # mövcud görüntü DƏYİŞMİR.
        centroids = geometry.cell_centroid().reshape(model.grid.shape + (3,))
        # Lülə modelin SƏTHİNDƏN YUXARIDA başlayır — əks halda tamamilə
        # hüceyrələrin içində gizlənərdi (ilk render sınağında məhz bu
        # baş verdi: yalnız adlar görünürdü). Yuxarı uzantı modelin
        # ümumi qalınlığının 25 %-i qədərdir ki, istənilən qalınlıqda
        # nisbətli görünsün.
        model_thickness = max(float(depths.max() - depths.min()),
                              float(geometry.dz.mean()))
        surface = float(depths.min()) - float(geometry.dz[0]) * 0.5
        top = surface - model_thickness * 0.25

        grid = model.grid
        for well in model.active_wells():
            perforations = [p for p in well.open_perforations()
                            if 0 <= p.i < grid.nx and 0 <= p.j < grid.ny
                            and 0 <= p.k < grid.nz]
            if not perforations:
                continue
            head = perforations[0]
            x = float(centroids[head.k, head.j, head.i, 0])
            y = float(centroids[head.k, head.j, head.i, 1])
            z_values = [depths[p.k, p.j, p.i] for p in perforations]
            colour = ((0.165, 0.655, 0.627) if well.well_type is WellType.INJECTOR
                      else (0.851, 0.557, 0.169))     # PALETTE.water / .oil

            deepest = max(z_values)
            self._well_actors.append(self._well_bore(
                x, y, -top * exaggeration, -deepest * exaggeration, colour,
                geometry))
            for depth in z_values:
                self._well_actors.append(self._perforation_marker(
                    x, y, -depth * exaggeration, colour, geometry))
            self._well_actors.append(self._well_label(
                well.name, x, y, -top * exaggeration, geometry))

        for actor in self._well_actors:
            self.renderer.AddViewProp(actor)

    def _well_bore(self, x, y, z_top, z_bottom, colour, geometry):
        import vtk

        source = vtk.vtkCylinderSource()
        source.SetRadius(min(geometry.dx, geometry.dy) * 0.09)
        source.SetHeight(abs(z_top - z_bottom))
        source.SetResolution(16)
        source.CappingOn()

        # vtkCylinderSource Y oxu boyunca qurur — Z-ə çevirmək lazımdır
        transform = vtk.vtkTransform()
        transform.Translate(x, y, (z_top + z_bottom) / 2.0)
        transform.RotateX(90)
        # `vtkTransformFilter` — `vtkTransformPolyDataFilter` VTK 9.7-də
        # köhnəlib (deprecated)
        filter_ = vtk.vtkTransformFilter()
        filter_.SetTransform(transform)
        filter_.SetInputConnection(source.GetOutputPort())

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(filter_.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(*colour)
        return actor

    def _perforation_marker(self, x, y, z, colour, geometry):
        import vtk

        source = vtk.vtkSphereSource()
        source.SetCenter(x, y, z)
        source.SetRadius(min(geometry.dx, geometry.dy) * 0.16)
        source.SetThetaResolution(14)
        source.SetPhiResolution(14)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(*colour)
        return actor

    def _well_label(self, name, x, y, z, geometry):
        """Ad — `vtkBillboardTextActor3D`: model fırlananda mətn HƏMİŞƏ
        kameraya baxır (adi 3D mətn fırlanıb oxunmaz olardı)."""
        import vtk

        actor = vtk.vtkBillboardTextActor3D()
        actor.SetInput(f" {name}")
        actor.SetPosition(x, y, z + float(geometry.dz.mean()) * 0.6)
        properties = actor.GetTextProperty()
        properties.SetColor(1.0, 1.0, 1.0)
        properties.SetFontSize(15)
        properties.BoldOn()
        properties.ShadowOn()
        return actor

    # ═══════════════════════════════════════════════════════ faultlar
    def update_faults(self):
        """Fault müstəviləri — VTK-da ƏSL DOLĞU (matplotlib-dən fərqli).

        matplotlib motorunda dolğu işlətmək mümkün olmamışdı: onun 3D
        dərinlik sıralaması dəqiq olmadığı üçün yarı-şəffaf müstəvi
        hüceyrə üzləri ilə kəsişəndə "cırıq" görünürdü, ona görə orada
        yalnız KONTUR çəkilir (bax `FAULTS.md`). VTK-da əsl z-buffer
        var — dolğu düzgün görünür.

        Şəffaflıq çarpandan asılıdır: tam bağlı (sealing) fault daha
        qeyri-şəffaf, şəffaf fault daha açıq.
        """
        import vtk

        for actor in self._fault_actors:
            self.renderer.RemoveViewProp(actor)
        self._fault_actors = []

        if not self.settings.show_faults:
            return
        if not getattr(self.model, "fault_references", None):
            return

        polygons, multipliers = self._fault_planes()
        for polygon, multiplier in zip(polygons, multipliers):
            points = vtk.vtkPoints()
            for corner in polygon:
                points.InsertNextPoint(*corner)
            quad = vtk.vtkQuad()
            for index in range(4):
                quad.GetPointIds().SetId(index, index)
            cells = vtk.vtkCellArray()
            cells.InsertNextCell(quad)
            data = vtk.vtkPolyData()
            data.SetPoints(points)
            data.SetPolys(cells)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(data)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            properties = actor.GetProperty()
            properties.SetColor(0.753, 0.341, 0.294)      # PALETTE.danger
            properties.SetOpacity(0.25 + 0.55 * (1.0 - min(max(multiplier, 0.0), 1.0)))
            properties.EdgeVisibilityOn()
            properties.SetEdgeColor(0.851, 0.400, 0.345)
            properties.SetLineWidth(2.0)
            self._fault_actors.append(actor)
            self.renderer.AddViewProp(actor)

    def _fault_planes(self):
        """Fault müstəvilərinin künc nöqtələri — `volume.py`-dəki
        `_fault_polygons()` ilə EYNİ həndəsə (hədd qutusu yanaşması)."""
        model = self.model
        grid, geometry = model.grid, model.geometry
        nz, ny, nx = grid.shape
        exaggeration = max(self.settings.vertical_exaggeration, 1e-6)
        depths = geometry.cell_depths().reshape(grid.shape)
        half = float(geometry.dz.mean()) * 0.5

        polygons, multipliers = [], []
        for fault in model.fault_references:
            if not fault.has_geometry:
                continue
            axis = fault.axis.upper()
            plane = fault.plane_index

            if axis == "I" and 0 <= plane < nx - 1:
                j_low, j_high = fault.range_a or (0, ny - 1)
                k_low, k_high = fault.range_b or (0, nz - 1)
                j_low, j_high = max(j_low, 0), min(j_high, ny - 1)
                k_low, k_high = max(k_low, 0), min(k_high, nz - 1)
                x = (plane + 1) * geometry.dx
                column = depths[k_low:k_high + 1, j_low:j_high + 1,
                                plane:plane + 2]
                z_top = -(float(column.min()) - half) * exaggeration
                z_bottom = -(float(column.max()) + half) * exaggeration
                y0, y1 = j_low * geometry.dy, (j_high + 1) * geometry.dy
                polygons.append([(x, y0, z_top), (x, y1, z_top),
                                (x, y1, z_bottom), (x, y0, z_bottom)])
            elif axis == "J" and 0 <= plane < ny - 1:
                i_low, i_high = fault.range_a or (0, nx - 1)
                k_low, k_high = fault.range_b or (0, nz - 1)
                i_low, i_high = max(i_low, 0), min(i_high, nx - 1)
                k_low, k_high = max(k_low, 0), min(k_high, nz - 1)
                y = (plane + 1) * geometry.dy
                column = depths[k_low:k_high + 1, plane:plane + 2,
                                i_low:i_high + 1]
                z_top = -(float(column.min()) - half) * exaggeration
                z_bottom = -(float(column.max()) + half) * exaggeration
                x0, x1 = i_low * geometry.dx, (i_high + 1) * geometry.dx
                polygons.append([(x0, y, z_top), (x1, y, z_top),
                                (x1, y, z_bottom), (x0, y, z_bottom)])
            elif axis == "K" and 0 <= plane < nz - 1:
                i_low, i_high = fault.range_a or (0, nx - 1)
                j_low, j_high = fault.range_b or (0, ny - 1)
                i_low, i_high = max(i_low, 0), min(i_high, nx - 1)
                j_low, j_high = max(j_low, 0), min(j_high, ny - 1)
                column = depths[plane:plane + 2, j_low:j_high + 1,
                                i_low:i_high + 1]
                z = -float(column.mean()) * exaggeration
                x0, x1 = i_low * geometry.dx, (i_high + 1) * geometry.dx
                y0, y1 = j_low * geometry.dy, (j_high + 1) * geometry.dy
                polygons.append([(x0, y0, z), (x1, y0, z),
                                (x1, y1, z), (x0, y1, z)])
            else:
                continue
            multipliers.append(fault.effective_multiplier)
        return polygons, multipliers

    def _update_scalar_bar(self, table, label: str):
        """Rəng legendi — ResInsight-dakı kimi sağ tərəfdə, şaquli."""
        import vtk

        if self._scalar_bar is None:
            self._scalar_bar = vtk.vtkScalarBarActor()
            self._scalar_bar.SetOrientationToVertical()
            self._scalar_bar.SetPosition(0.90, 0.15)
            self._scalar_bar.SetWidth(0.07)
            self._scalar_bar.SetHeight(0.70)
            self._scalar_bar.SetNumberOfLabels(6)
            for text in (self._scalar_bar.GetLabelTextProperty(),
                        self._scalar_bar.GetTitleTextProperty()):
                text.SetColor(0.863, 0.894, 0.918)      # PALETTE.text
                text.SetFontSize(12)
                text.ItalicOff()
                text.BoldOff()
                text.ShadowOff()
            # `AddViewProp` — VTK 9.7-də `AddActor2D` artıq yoxdur
            self.renderer.AddViewProp(self._scalar_bar)
        self._scalar_bar.SetLookupTable(table)
        self._scalar_bar.SetTitle(label or "")

    # ── kamera ─────────────────────────────────────────────────────
    def reset_camera(self, view: Optional[str] = None):
        """Baxış bucağı — matplotlib motorundakı adlarla eyni."""
        camera = self.renderer.GetActiveCamera()
        bounds = self._grid.GetBounds()
        centre = ((bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2,
                  (bounds[4] + bounds[5]) / 2)
        span = max(bounds[1] - bounds[0], bounds[3] - bounds[2],
                   bounds[5] - bounds[4], 1.0)

        presets = {
            "Yuxarıdan": (centre[0], centre[1], centre[2] + span * 2.0),
            "Yandan (X)": (centre[0] + span * 2.0, centre[1], centre[2]),
            "Yandan (Y)": (centre[0], centre[1] + span * 2.0, centre[2]),
            "İzometrik": (centre[0] + span, centre[1] - span,
                          centre[2] + span * 0.8),
        }
        position = presets.get(view, presets["İzometrik"])
        camera.SetFocalPoint(*centre)
        camera.SetPosition(*position)
        camera.SetViewUp(0, 0, 1)
        self.renderer.ResetCamera()
        self.apply_zoom()
        self.renderer.ResetCameraClippingRange()

    def apply_zoom(self):
        """Yaxınlaşdırma — `ResetCamera()`-dan SONRA tətbiq olunmalıdır
        (o, kameranı tam modelə sığdırır, sonra biz miqyaslayırıq)."""
        zoom = self.settings.zoom
        if zoom and abs(zoom - 1.0) > 1e-6:
            self.renderer.GetActiveCamera().Zoom(zoom)
