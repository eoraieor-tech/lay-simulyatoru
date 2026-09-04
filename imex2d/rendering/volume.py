"""3D həcm görüntüsü.

Əsas texniki məsələ — performans. Naiv yanaşmada hər hüceyrə üçün
6 üz çəkilir: 10 000 hüceyrə → 60 000 çoxbucaqlı, matplotlib bunu
dəqiqələrlə çəkir.

Həll: yalnız GÖRÜNƏN üzlər çəkilir. Üz görünəndir əgər:
  · modelin sərhədindədirsə, YAXUD
  · qonşusu filtr ilə gizlədilibsə

10 000 hüceyrəli modeldə bu, 60 000 üzü ~2 400-ə endirir — 25 dəfə az.

Filtr və kəsim də məhz bu mexanizmlə işləyir: hüceyrə gizlədiləndə
onun qonşularının üzləri "açılır" və daxili struktur görünür.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from matplotlib.colors import Normalize
from matplotlib.ticker import FixedLocator, MaxNLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from ..domain.reservoir_model import ReservoirModel
from ..domain.wells import WellType
from .theme import PALETTE


def _hex_to_rgb(colour: str):
    colour = colour.lstrip("#")
    return tuple(int(colour[i:i + 2], 16) / 255.0 for i in (0, 2, 4))

# hüceyrənin 6 üzü: (ox, istiqamət) və künclərin lokal indeksləri
# küncler: (i,j,k) ofsetləri ilə 0..7
_CORNERS = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                     [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)

_FACES = {
    (0, -1): (0, 3, 7, 4),   # i−
    (0, +1): (1, 2, 6, 5),   # i+
    (1, -1): (0, 1, 5, 4),   # j−
    (1, +1): (3, 2, 6, 7),   # j+
    (2, -1): (0, 1, 2, 3),   # k− (yuxarı, dərinlik azalır)
    (2, +1): (4, 5, 6, 7),   # k+
}

# üzlərin xarici normalları — işıqlandırma üçün
_NORMALS = {
    (0, -1): np.array([-1.0, 0.0, 0.0]),
    (0, +1): np.array([+1.0, 0.0, 0.0]),
    (1, -1): np.array([0.0, -1.0, 0.0]),
    (1, +1): np.array([0.0, +1.0, 0.0]),
    (2, -1): np.array([0.0, 0.0, -1.0]),
    (2, +1): np.array([0.0, 0.0, +1.0]),
}

# 100 % yaxınlaşdırmada modelin çərçivəni doldurması üçün baza əmsalı.
# matplotlib-in 3D oxu defolt olaraq geniş boş kənar buraxır.
BASE_FIT = 1.35

# işıq mənbəyinin istiqaməti (yuxarı-sol-öndən)
_LIGHT = np.array([-0.4, -0.6, -0.7])
_LIGHT = _LIGHT / np.linalg.norm(_LIGHT)

# hazır baxış bucaqları: (elevation, azimuth)
VIEW_ANGLES = {
    "İzometrik": (22.0, -60.0),
    "Üstdən": (89.0, -90.0),
    "Yandan (X)": (4.0, -90.0),
    "Yandan (Y)": (4.0, 0.0),
    "Künc": (35.0, -135.0),
}


@dataclass
class VolumeFilter:
    """Hansı hüceyrələrin göstəriləcəyi."""
    value_min: Optional[float] = None
    value_max: Optional[float] = None
    i_range: Optional[Tuple[int, int]] = None
    j_range: Optional[Tuple[int, int]] = None
    k_range: Optional[Tuple[int, int]] = None

    def mask(self, values: np.ndarray, shape: tuple) -> np.ndarray:
        """(nz, ny, nx) formasında bool maska."""
        visible = np.ones(shape, dtype=bool)
        grid_values = values.reshape(shape)

        if self.value_min is not None:
            visible &= grid_values >= self.value_min
        if self.value_max is not None:
            visible &= grid_values <= self.value_max

        nz, ny, nx = shape
        for axis, limits, size in ((2, self.i_range, nx),
                                   (1, self.j_range, ny),
                                   (0, self.k_range, nz)):
            if limits is None:
                continue
            low, high = max(limits[0], 0), min(limits[1], size - 1)
            index = np.zeros(size, dtype=bool)
            index[low:high + 1] = True
            visible &= np.expand_dims(
                index, axis=tuple(a for a in range(3) if a != axis))
        return visible


def apply_zoom(ax, zoom: float) -> bool:
    """Yaxınlaşmanı SƏHNƏNİ YENİDƏN ÇƏKMƏDƏN tətbiq edir.

    Siçan çarxı ilə yaxınlaşdırma hər addımda bütün üzləri yenidən
    hesablasaydı, 50 000 hüceyrəli modeldə hərəkət kəsikli olardı.
    Burada yalnız kamera nisbəti dəyişir.
    """
    aspect = getattr(ax, "_imex_aspect", None)
    if aspect is None:
        return False
    compensation = getattr(ax, "_imex_compensation", 1.0)
    try:
        ax.set_box_aspect(aspect,
                          zoom=compensation * max(float(zoom), 0.05))
    except (TypeError, AttributeError):
        return False
    return True


class VolumeRenderer:
    """3D hüceyrə görüntüsü — yalnız görünən üzlər."""

    def draw(self, ax, figure, model: ReservoirModel, values: np.ndarray,
             label: str = "", colormap="viridis",
             value_limits: Optional[Tuple[float, float]] = None,
             volume_filter: Optional[VolumeFilter] = None,
             show_wells: bool = True, show_faults: bool = True,
             vertical_exaggeration: float = 1.0,
             edge_width: float = 0.15, cax=None,
             shading: float = 0.45, view: Optional[str] = None,
             opacity: float = 1.0, zoom: float = 1.0):
        ax.clear()
        if cax is not None:
            cax.clear()

        grid = model.grid
        shape = grid.shape
        values = np.asarray(values, float).ravel()
        visible = (volume_filter or VolumeFilter()).mask(values, shape)

        polygons, face_values, shades = self._visible_faces(
            model, values, visible, vertical_exaggeration)

        if value_limits is None:
            shown = values.reshape(shape)[visible]
            value_limits = ((float(shown.min()), float(shown.max()))
                            if shown.size else (0.0, 1.0))
        norm = Normalize(*value_limits)

        collection = Poly3DCollection(
            polygons, linewidths=edge_width,
            edgecolors=(0, 0, 0, 0.25) if edge_width > 0 else "none")
        collection.set_cmap(colormap)
        collection.set_norm(norm)
        collection.set_array(np.asarray(face_values))

        if shading > 0.0 and len(polygons):
            colours = self._shaded_colours(collection, face_values, norm,
                                           shades, shading, opacity)
            collection.set_facecolors(colours)
        elif opacity < 1.0:
            collection.set_alpha(opacity)
        ax.add_collection3d(collection)

        if show_wells:
            self._draw_wells(ax, model, vertical_exaggeration)

        if show_faults:
            self._draw_faults(ax, model, vertical_exaggeration)

        self._style(ax, model, vertical_exaggeration, label, zoom)
        if view in VIEW_ANGLES:
            elevation, azimuth = VIEW_ANGLES[view]
            ax.view_init(elev=elevation, azim=azimuth)

        bar = None
        if cax is not None:
            bar = figure.colorbar(collection, cax=cax)
            bar.set_label(label, color=PALETTE.text_dim, fontsize=9)
            bar.ax.tick_params(colors=PALETTE.text_dim, labelsize=8)
            bar.outline.set_edgecolor(PALETTE.line)
        return bar

    @staticmethod
    def _shaded_colours(collection, face_values, norm, shades, strength,
                        opacity):
        """Rəng xəritəsinin üzərinə Lambert parlaqlığı.

        `strength` = 0 -> düz rəng (rəqəm oxumaq üçün ən dəqiq)
        `strength` = 1 -> güclü kölgə (forma ən yaxşı görünür)
        """
        colours = collection.get_cmap()(norm(np.asarray(face_values)))
        factor = 1.0 - strength + strength * np.clip(shades, 0.0, 1.0)
        colours[:, :3] *= factor[:, None]
        colours[:, 3] = opacity
        return np.clip(colours, 0.0, 1.0)

    # ═══════════════════════════════════════════ görünən üz çıxarışı
    @staticmethod
    def _visible_faces(model: ReservoirModel, values: np.ndarray,
                       visible: np.ndarray, exaggeration: float):
        """Yalnız sərhəd və gizlədilmiş qonşuya baxan üzləri qaytarır."""
        grid = model.grid
        geometry = model.geometry
        nz, ny, nx = grid.shape
        depths = geometry.cell_depths().reshape(grid.shape)
        grid_values = values.reshape(grid.shape)

        # Koordinatlar HƏQİQİ dərinlikdədir. Şaquli mübaliğə yalnız
        # görüntü nisbətinə tətbiq olunur (bax `_style`) — əks halda
        # hüceyrələr uzadılmış, ox etiketləri isə həqiqi qalır və ikisi
        # bir-birinə zidd olur.
        dz_grid = geometry.dz_per_cell().reshape(grid.shape)
        polygons, face_values, shades = [], [], []

        for (axis, direction), corner_indices in _FACES.items():
            neighbour = np.roll(visible, -direction, axis=axis)
            # sərhəd təbəqəsində qonşu yoxdur -> üz həmişə görünür
            edge = [slice(None)] * 3
            edge[axis] = -1 if direction > 0 else 0
            neighbour = neighbour.copy()
            neighbour[tuple(edge)] = False

            exposed = visible & ~neighbour
            if not exposed.any():
                continue

            k, j, i = np.nonzero(exposed)
            dz_k = dz_grid[k, j, i]
            origin = np.column_stack([i * geometry.dx, j * geometry.dy,
                                      depths[k, j, i] - dz_k * 0.5])
            # təbəqə qalınlığı dəyişkən ola bilər -> hər üz öz dz-i ilə
            # miqyaslanır (x/y ölçüsü sabit, z ölçüsü hüceyrəyə görə)
            cell_size = np.column_stack([
                np.full_like(dz_k, geometry.dx), np.full_like(dz_k, geometry.dy), dz_k])
            corners = _CORNERS[list(corner_indices)]
            quad = origin[:, None, :] + corners[None, :, :] * cell_size[:, None, :]
            polygons.extend(quad)
            face_values.extend(grid_values[k, j, i])

            # Lambert işıqlandırması: parlaqlıq üzün işığa baxma bucağından.
            # Bütün üzlər eyni rəngdə olanda model yastı görünür və
            # dərinlik hiss olunmur — bu, 3D-də ən çox itən məlumatdır.
            brightness = float(np.dot(_NORMALS[(axis, direction)], -_LIGHT))
            shades.extend([brightness] * len(quad))

        return polygons, face_values, np.asarray(shades)

    # ═══════════════════════════════════════════════════════ faultlar
    @staticmethod
    def _fault_polygons(model: ReservoirModel, exaggeration: float):
        """Fault müstəviləri — hər fault üçün BİR həmvar (flat) quad.

        Hər hüceyrə cütünə ayrıca quad (dip-i dəqiq izləyən) yaradıla
        bilərdi, lakin çoxlu nazik, demək olar üst-üstə düşən üzlər
        matplotlib-in 3D dərinlik sıralamasında "cırıq" görünür (məlum
        məhdudiyyət). Ona görə fault öz tam diapazonunun HƏDD
        qutusunu (bounding box) əhatə edən tək düz səth kimi çəkilir —
        Petrel/CMG kimi vasitələr də sadələşdirilmiş fault
        görüntüsündə eyni yanaşmanı işlədir.
        """
        grid, geometry = model.grid, model.geometry
        depths = geometry.cell_depths().reshape(grid.shape)
        nz, ny, nx = grid.shape
        half_thickness = float(geometry.dz.mean()) * 0.5

        polygons, multipliers, labels = [], [], []
        for fault in model.fault_references:
            if not fault.has_geometry:
                continue
            axis = fault.axis.upper()
            plane = fault.plane_index

            if axis == "I":
                if not (0 <= plane < nx - 1):
                    continue
                j_low, j_high = fault.range_a or (0, ny - 1)
                k_low, k_high = fault.range_b or (0, nz - 1)
                j_low, j_high = max(j_low, 0), min(j_high, ny - 1)
                k_low, k_high = max(k_low, 0), min(k_high, nz - 1)
                x = (plane + 1) * geometry.dx
                column = depths[k_low:k_high + 1, j_low:j_high + 1,
                                plane:plane + 2]
                z_top = float(column.min()) - half_thickness
                z_bottom = float(column.max()) + half_thickness
                y0, y1 = j_low * geometry.dy, (j_high + 1) * geometry.dy
                polygons.append(np.array([
                    [x, y0, z_top], [x, y1, z_top],
                    [x, y1, z_bottom], [x, y0, z_bottom]]))
                labels.append((fault.name, x, 0.5 * (y0 + y1),
                               0.5 * (z_top + z_bottom)))

            elif axis == "J":
                if not (0 <= plane < ny - 1):
                    continue
                i_low, i_high = fault.range_a or (0, nx - 1)
                k_low, k_high = fault.range_b or (0, nz - 1)
                i_low, i_high = max(i_low, 0), min(i_high, nx - 1)
                k_low, k_high = max(k_low, 0), min(k_high, nz - 1)
                y = (plane + 1) * geometry.dy
                column = depths[k_low:k_high + 1, plane:plane + 2,
                                i_low:i_high + 1]
                z_top = float(column.min()) - half_thickness
                z_bottom = float(column.max()) + half_thickness
                x0, x1 = i_low * geometry.dx, (i_high + 1) * geometry.dx
                polygons.append(np.array([
                    [x0, y, z_top], [x1, y, z_top],
                    [x1, y, z_bottom], [x0, y, z_bottom]]))
                labels.append((fault.name, 0.5 * (x0 + x1), y,
                               0.5 * (z_top + z_bottom)))

            else:  # axis == "K" — üfüqi (təbəqələr arası) fault
                if not (0 <= plane < nz - 1):
                    continue
                i_low, i_high = fault.range_a or (0, nx - 1)
                j_low, j_high = fault.range_b or (0, ny - 1)
                i_low, i_high = max(i_low, 0), min(i_high, nx - 1)
                j_low, j_high = max(j_low, 0), min(j_high, ny - 1)
                column = depths[plane:plane + 2, j_low:j_high + 1,
                                i_low:i_high + 1]
                z = float(column.mean())
                x0, x1 = i_low * geometry.dx, (i_high + 1) * geometry.dx
                y0, y1 = j_low * geometry.dy, (j_high + 1) * geometry.dy
                polygons.append(np.array([
                    [x0, y0, z], [x1, y0, z], [x1, y1, z], [x0, y1, z]]))
                labels.append((fault.name, 0.5 * (x0 + x1), 0.5 * (y0 + y1), z))

            multipliers.append(fault.effective_multiplier)

        return polygons, multipliers, labels

    @classmethod
    def _draw_faults(cls, ax, model: ReservoirModel, exaggeration: float) -> None:
        """Fault müstəvilərini QALIN KONTUR kimi çəkir.

        Dolğu (Poly3DCollection) yerinə kontur seçilib, çünki
        matplotlib-in 3D dərinlik-sıralaması iri, yarı-şəffaf poliqon
        hüceyrə üzləri ilə kəsişəndə "cırıq"/qara zolaq artefaktı
        yaradır (mplot3d-nin məlum məhdudiyyəti — z-buffer deyil,
        mərkəz nöqtəsinə görə təxmini sıralama işlədir). Xətlər bu
        problemə məruz qalmır və nəticə daha aydın oxunur.

        Xəttin qalınlığı çarpandan asılıdır: tam bağlı (sealing,
        çarpan 0) fault qalın, şəffaf fault (çarpan 1) nazik xətlə
        çəkilir — istifadəçi cədvələ baxmadan hansı faultun axını nə
        qədər məhdudlaşdırdığını görür.
        """
        polygons, multipliers, labels = cls._fault_polygons(model, exaggeration)
        if not polygons:
            return

        for polygon, multiplier in zip(polygons, multipliers):
            # Dolğu qəsdən YOXDUR: hətta çox aşağı alfa ilə belə,
            # matplotlib-in dərinlik-sıralama artefaktı (yuxarıdakı
            # sənədləşmə) yüngül zolaqlar kimi qalır. Xalis kontur
            # bu problemsiz və daha oxunaqlıdır.
            width = 1.2 + 2.3 * (1.0 - np.clip(multiplier, 0.0, 1.0))
            loop = np.vstack([polygon, polygon[:1]])
            ax.plot(loop[:, 0], loop[:, 1], loop[:, 2],
                   color=PALETTE.danger, lw=width, zorder=12)

        seen = set()
        for name, x, y, z in labels:
            if name in seen:
                continue
            seen.add(name)
            ax.text(x, y, z, f" {name}", color=PALETTE.danger, fontsize=8,
                   weight="bold", zorder=13)

    # ═══════════════════════════════════════════════════════ quyular
    @staticmethod
    def _draw_wells(ax, model: ReservoirModel, exaggeration: float) -> None:
        geometry = model.geometry
        depths = geometry.cell_depths().reshape(model.grid.shape)
        top = float(depths.min()) - float(geometry.dz[0]) * 0.5

        grid = model.grid
        for well in model.active_wells():
            perforations = [p for p in well.open_perforations()
                            if 0 <= p.i < grid.nx and 0 <= p.j < grid.ny
                            and 0 <= p.k < grid.nz]
            if not perforations:
                continue
            x = (perforations[0].i + 0.5) * geometry.dx
            y = (perforations[0].j + 0.5) * geometry.dy
            z_values = [depths[p.k, p.j, p.i] for p in perforations]
            colour = (PALETTE.water if well.well_type is WellType.INJECTOR
                      else PALETTE.oil)

            # lülə: səthdən perforasiyanın altına qədər
            ax.plot([x, x], [y, y], [top, max(z_values)],
                    color=colour, lw=2.5, zorder=10)
            ax.scatter([x] * len(z_values), [y] * len(z_values), z_values,
                       s=18, c=colour, edgecolors="white", linewidths=0.6,
                       depthshade=False, zorder=11)
            ax.text(x, y, top, f" {well.name}", color="white", fontsize=8,
                    weight="bold", zorder=12)

    # ═══════════════════════════════════════════════════ görünüş
    @staticmethod
    def _style(ax, model: ReservoirModel, exaggeration: float,
               label: str, zoom: float = 1.0) -> None:
        grid, geometry = model.grid, model.geometry
        depths = geometry.cell_depths()
        length_x = grid.nx * geometry.dx
        length_y = grid.ny * geometry.dy

        # modelin HƏQİQİ dərinlik intervalı: hüceyrə mərkəzləri ± yarım DZ
        top = float(depths.min()) - float(geometry.dz[0]) * 0.5
        base = float(depths.max()) + float(geometry.dz[-1]) * 0.5
        thickness = base - top

        ax.set_xlim(0, length_x)
        ax.set_ylim(0, length_y)
        margin = thickness * 0.05
        ax.set_zlim(base + margin, top - margin)

        # ─── dərinlik oxu TƏBƏQƏ SƏRHƏDLƏRİNDƏ işarələnir ───────────
        # İxtiyari addımlar (2050, 2150 …) təbəqələrlə üst-üstə düşmür
        # və hansı hüceyrənin harada bitdiyini görmək olmur.
        boundaries = top + np.concatenate(([0.0], np.cumsum(geometry.dz)))
        if boundaries.size <= 12:
            ax.zaxis.set_major_locator(FixedLocator(boundaries))
        else:
            step = int(np.ceil(boundaries.size / 10))
            ax.zaxis.set_major_locator(FixedLocator(boundaries[::step]))

        # Şaquli miqyas qutunu uzadanda plan sıxılır və X/Y etiketləri
        # bir-birinə qarışır — ona görə onların sayı məhdudlaşdırılır.
        tick_count = 6 if exaggeration <= 3.0 else 4
        ax.xaxis.set_major_locator(MaxNLocator(tick_count, integer=False))
        ax.yaxis.set_major_locator(MaxNLocator(tick_count, integer=False))

        # ─── şaquli miqyas ─────────────────────────────────────────
        # Mübaliğə YALNIZ qutunun nisbətinə tətbiq olunur, koordinatlara
        # yox — ona görə ox etiketləri həqiqi dərinliyi göstərir.
        #
        #     ekran nisbəti = (qalınlıq / plan) × mübaliğə
        #
        # `zoom` olmadan matplotlib qutunu ən böyük ölçüyə görə
        # normallaşdırır: hündürlük artanda plan kiçilir və model nazik
        # qülləyə çevrilir. `zoom` bunu kompensasiya edir.
        try:
            plan = max(length_x, length_y)
            height = max(thickness * exaggeration / plan, 0.02)
            aspect = (length_x / plan, length_y / plan, height)
            ax._imex_aspect = aspect
            ax._imex_compensation = 1.0 + 0.20 * min(max(height - 1.0, 0.0), 3.0)
            # Üç töhfə:
            #   BASE_FIT — matplotlib 3D defolt olaraq çox boş kənar
            #              buraxır; 100 % model çərçivəni doldurmalıdır
            #   compensation — hündür qutu normallaşdırmada kiçilir
            #   zoom — istifadəçinin seçimi
            compensation = 1.0 + 0.20 * min(max(height - 1.0, 0.0), 3.0)
            try:
                ax.set_box_aspect(
                    aspect, zoom=BASE_FIT * compensation * max(zoom, 0.05))
            except TypeError:              # matplotlib < 3.6
                ax.set_box_aspect(aspect)
        except AttributeError:
            pass

        ax.set_facecolor(PALETTE.background)
        ax.set_xlabel("X, m", color=PALETTE.text_dim, fontsize=9)
        ax.set_ylabel("Y, m", color=PALETTE.text_dim, fontsize=9)
        ax.set_zlabel("Dərinlik, m", color=PALETTE.text_dim, fontsize=9)
        ax.tick_params(colors=PALETTE.text_dim, labelsize=8)
        for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
            pane.pane.set_facecolor(PALETTE.panel)
            pane.pane.set_edgecolor(PALETTE.line)
            pane.pane.set_alpha(0.6)
        if label:
            ax.set_title(label, color=PALETTE.text, fontsize=10, pad=6)
