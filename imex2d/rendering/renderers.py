"""Çəkmə (rendering) qatı — matplotlib, Qt YOXDUR.

Hər renderer hazır ox (Axes) qəbul edir və üzərinə çəkir. Bu sayədə
eyni kod həm interfeysdə, həm PDF hesabatda, həm də testdə işləyir.
Rendering heç vaxt hesablama aparmır — yalnız hazır nəticəni göstərir.
"""

from __future__ import annotations
from typing import Optional

import numpy as np

from ..domain.reservoir_model import ReservoirModel
from ..domain.scal import CoreyParameters
from ..domain.wells import WellType
from ..simulation.results import SimulationResult, Snapshot
from .theme import (PALETTE, PERMEABILITY_CMAP, POROSITY_CMAP, PRESSURE_CMAP,
                    SATURATION_CMAP, legend, style_axes)

SATURATION = "SW"
PRESSURE = "PRESSURE"
PERMEABILITY = "PERMX"
POROSITY = "PORO"
DEPTH = "DEPTH"

PROPERTY_LABELS = {
    SATURATION: "Su doyumluluğu (Sw)",
    PRESSURE: "Təzyiq (bar)",
    PERMEABILITY: "Keçiricilik Kx (mD)",
    POROSITY: "Məsaməlilik φ",
    DEPTH: "Dərinlik (m)",
}

# ── LAY-MƏLUMATLI görüntü açarları (tapşırıq §13) ─────────────────────
# Açar formatı `"PROVENANCE:PORO"` / `"CONFIDENCE:PERMX"` — XASSƏ ADI
# AÇARIN İÇİNDƏDİR, ona görə burada HEÇ BİR xassə adı SABİT KODLANMIR
# (§15). UI mövcud `model.provenance` açarlarından siyahını ÖZÜ qurur.
# Bu qat HESABLAMA APARMIR — yalnız artıq hazır `PropertyProvenance`
# massivlərini göstərir (§13: "Visualization core geology hesablamalarını
# özündə daşımamalıdır").
PROVENANCE_PREFIX = "PROVENANCE:"
CONFIDENCE_PREFIX = "CONFIDENCE:"
ORIGINAL_PREFIX = "ORIGINAL:"
IMPACT_PREFIX = "IMPACT:"

#: Bu qatın tanıdığı BÜTÜN lay-məlumatlı prefikslər.
PROVENANCE_PREFIXES = (PROVENANCE_PREFIX, CONFIDENCE_PREFIX,
                       ORIGINAL_PREFIX, IMPACT_PREFIX)


def provenance_key(name: str) -> str:
    return f"{PROVENANCE_PREFIX}{name}"


def confidence_key(name: str) -> str:
    return f"{CONFIDENCE_PREFIX}{name}"


def original_key(name: str) -> str:
    return f"{ORIGINAL_PREFIX}{name}"


def impact_key(name: str) -> str:
    return f"{IMPACT_PREFIX}{name}"


def property_label(key: str) -> str:
    """Combo/başlıq üçün etiket — naməlum açar öz adını qaytarır."""
    if key.startswith(PROVENANCE_PREFIX):
        return f"Mənşə/status — {key[len(PROVENANCE_PREFIX):]}"
    if key.startswith(CONFIDENCE_PREFIX):
        return f"Etibarlılıq balı — {key[len(CONFIDENCE_PREFIX):]}"
    if key.startswith(ORIGINAL_PREFIX):
        return f"Orijinal sahə — {key[len(ORIGINAL_PREFIX):]}"
    if key.startswith(IMPACT_PREFIX):
        return f"Təsir (final − orijinal) — {key[len(IMPACT_PREFIX):]}"
    return PROPERTY_LABELS.get(key, key)


def _provenance_arrays(model, key: str):
    """`(dəyərlər, colormap, vmin, vmax)` — status/etibar/orijinal/təsir.

    Provenance (və ya orijinal sahə) yoxdursa HAMISI `NaN` qaytarılır —
    uydurma "hər şey ölçülüb"/"təsir sıfırdır" mənzərəsi YARADILMIR.

    QEYD (§12): TƏSİR burada YALNIZ GÖSTƏRİLİR — `final` sahəsi
    DƏYİŞDİRİLMİR, fərq hər çağırışda yenidən hesablanır.
    """
    from ..domain.data_availability import STATUS_CODE, STATUS_ORDER
    name = key.split(":", 1)[1]
    entry = getattr(model, "provenance", {}).get(name)
    ncell = model.grid.ncell
    if key.startswith(PROVENANCE_PREFIX):
        if entry is None:
            return np.full(ncell, np.nan), "tab10", 0.0, float(len(STATUS_ORDER) - 1)
        codes = np.asarray([STATUS_CODE.get(str(s), np.nan) for s in entry.status], float)
        return codes, "tab10", 0.0, float(len(STATUS_ORDER) - 1)
    if key.startswith(ORIGINAL_PREFIX):
        if entry is None or entry.original is None:
            return np.full(ncell, np.nan), POROSITY_CMAP, None, None
        return np.asarray(entry.original, float), POROSITY_CMAP, None, None
    if key.startswith(IMPACT_PREFIX):
        if entry is None or entry.original is None:
            return np.full(ncell, np.nan), "coolwarm", None, None
        delta = np.asarray(entry.final, float) - np.asarray(entry.original, float)
        limit = float(np.nanmax(np.abs(delta))) if np.isfinite(delta).any() else 0.0
        limit = limit or 1.0                      # sıfır təsir → dejenerativ şkala olmasın
        return delta, "coolwarm", -limit, limit
    if entry is None:
        return np.full(ncell, np.nan), "viridis", 0.0, 1.0
    return np.asarray(entry.confidence, float), "viridis", 0.0, 1.0


def _colorbar(figure, image, ax, cax, label):
    """Rəng şkalası.

    `cax` verildikdə şkala ÖNCƏDƏN ayrılmış oxda çəkilir. Bu vacibdir:
    `figure.colorbar(..., ax=ax)` hər çağırışda əsas oxdan yer alır və
    şkala silinəndə həmin yeri geri qaytarmır — nəticədə təkrar
    çəkilişlərdə qrafik tədricən daralır.
    """
    if cax is not None:
        bar = figure.colorbar(image, cax=cax)
    else:
        bar = figure.colorbar(image, ax=ax, fraction=0.04, pad=0.02)
    bar.set_label(label, color=PALETTE.text_dim, fontsize=9)
    bar.ax.tick_params(colors=PALETTE.text_dim, labelsize=8)
    bar.outline.set_edgecolor(PALETTE.line)
    return bar


class MapRenderer:
    """Areal xəritə + quyu simvolları."""

    def draw(self, ax, figure, model: ReservoirModel, property_key: str,
             snapshot: Optional[Snapshot] = None, colorbar=None, layer: int = 0,
             cax=None):
        ax.clear()
        if cax is not None:
            cax.clear()
        elif colorbar is not None:
            try:
                colorbar.remove()
            except Exception:
                pass

        data, cmap, vmin, vmax = self._select(model, property_key, snapshot, layer)
        lx, ly = model.geometry.areal_extent()
        image = ax.imshow(data, origin="lower", extent=[0, lx, 0, ly],
                          aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                          interpolation="nearest")

        bar = _colorbar(figure, image, ax, cax,
                        property_label(property_key))

        self._draw_wells(ax, model, layer)
        t = snapshot.time if snapshot else 0.0
        suffix = f"    |    təbəqə K = {layer + 1}" if model.grid.nz > 1 else ""
        style_axes(ax, f"{property_label(property_key)}"
                       f"    |    t = {t:.0f} gün{suffix}", "X, m", "Y, m")
        return bar

    def _select_volume(self, model, key, snapshot):
        """(nz, ny, nx) massivi — kəsik renderer-i üçün."""
        shape3d = model.grid.shape
        scal = model.scal_parameters
        if key.startswith(PROVENANCE_PREFIXES):
            values, cmap, low, high = _provenance_arrays(model, key)
            return values.reshape(shape3d), cmap, low, high
        if key == SATURATION:
            data = (np.asarray(snapshot.water_saturation).reshape(shape3d) if snapshot
                    else np.full(shape3d, model.initial_conditions.water_saturation))
            return data, SATURATION_CMAP, scal.swc, 1.0 - scal.sor
        if key == PRESSURE:
            data = (np.asarray(snapshot.pressure).reshape(shape3d) if snapshot
                    else np.full(shape3d, model.initial_conditions.datum_pressure))
            return data, PRESSURE_CMAP, None, None
        if key == PERMEABILITY:
            return model.rock.permx.values.reshape(shape3d), PERMEABILITY_CMAP, None, None
        if key == DEPTH:
            return model.geometry.cell_depths().reshape(shape3d), "magma_r", None, None
        return model.rock.porosity.values.reshape(shape3d), POROSITY_CMAP, None, None

    def _select(self, model, key, snapshot, layer=0):
        shape3d = model.grid.shape
        shape2d = (model.grid.ny, model.grid.nx)
        k = int(np.clip(layer, 0, model.grid.nz - 1))
        scal = model.scal_parameters

        def slice_of(flat_or_grid, default=None):
            if flat_or_grid is None:
                return np.full(shape2d, default)
            return np.asarray(flat_or_grid).reshape(shape3d)[k]

        if key.startswith(PROVENANCE_PREFIXES):
            values, cmap, low, high = _provenance_arrays(model, key)
            return slice_of(values), cmap, low, high
        if key == SATURATION:
            data = (slice_of(snapshot.water_saturation) if snapshot
                    else np.full(shape2d, model.initial_conditions.water_saturation))
            return data, SATURATION_CMAP, scal.swc, 1.0 - scal.sor
        if key == PRESSURE:
            data = (slice_of(snapshot.pressure) if snapshot
                    else np.full(shape2d, model.initial_conditions.datum_pressure))
            return data, PRESSURE_CMAP, None, None
        if key == PERMEABILITY:
            return slice_of(model.rock.permx.values), PERMEABILITY_CMAP, None, None
        if key == DEPTH:
            return slice_of(model.geometry.cell_depths()), "magma_r", None, None
        return slice_of(model.rock.porosity.values), POROSITY_CMAP, None, None

    @staticmethod
    def _draw_wells(ax, model: ReservoirModel, layer: int = 0):
        dx, dy = model.geometry.dx, model.geometry.dy
        for well in model.active_wells():
            for perf in well.open_perforations():
                if perf.k != layer:
                    continue
                x, y = (perf.i + 0.5) * dx, (perf.j + 0.5) * dy
                injector = well.well_type is WellType.INJECTOR
                ax.plot(x, y, "v" if injector else "^", ms=12,
                        mfc=PALETTE.water if injector else PALETTE.oil,
                        mec="white", mew=1.4, zorder=5)
                ax.annotate(well.name, (x, y), textcoords="offset points",
                            xytext=(9, 7), color="white", fontsize=8,
                            weight="bold", zorder=6)


class ProductionCurveRenderer:
    """Debit, sulaşma, kumulyativ hasilat və RF qrafikləri."""

    def draw(self, axes, result: SimulationResult):
        (a1, a2), (a3, a4) = axes
        for ax in (a1, a2, a3, a4):
            ax.clear()
        s = result.series
        if not s.time:
            return
        t = np.array(s.time)

        a1.plot(t, s.oil_rate, color=PALETTE.oil, lw=1.8, label="Neft")
        a1.plot(t, s.water_rate, color=PALETTE.water, lw=1.8, label="Su")
        a1.plot(t, s.water_injection_rate, color=PALETTE.accent, lw=1.2,
                ls="--", label="Vurulan su")
        style_axes(a1, "Debitlər", "Zaman, gün", "q, m³/gün")
        legend(a1)

        a2.plot(t, s.water_cut, color=PALETTE.water, lw=2)
        a2.set_ylim(0, 100)
        style_axes(a2, "Sulaşma (water cut)", "Zaman, gün", "WC, %")

        a3.plot(t, np.array(s.cumulative_oil) / 1e3, color=PALETTE.oil, lw=2,
                label="Kum. neft")
        a3.plot(t, np.array(s.cumulative_water) / 1e3, color=PALETTE.water,
                lw=1.5, label="Kum. su")
        style_axes(a3, "Kumulyativ hasilat", "Zaman, gün", "min m³")
        legend(a3)

        a4.plot(t, s.recovery_factor, color=PALETTE.oil, lw=2.2)
        twin = a4.twinx()
        twin.plot(t, s.average_pressure, color=PALETTE.text_dim, lw=1.2, ls=":")
        twin.set_ylabel("Orta P, bar", color=PALETTE.text_dim, fontsize=9)
        twin.tick_params(colors=PALETTE.text_dim, labelsize=8)
        for spine in twin.spines.values():
            spine.set_color(PALETTE.line)
        style_axes(a4, "Recovery Factor", "Zaman, gün", "RF, %")


class ScalRenderer:
    """Nisbi keçiricilik və fraksional axın əyriləri."""

    def draw(self, axes, scal: CoreyParameters, mu_w: float, mu_o: float,
             capillary=None):
        a1, a2 = axes
        a1.clear(); a2.clear()
        s = np.linspace(scal.swc, 1.0 - scal.sor, 250)

        a1.plot(s, scal.krw(s), color=PALETTE.water, lw=2, label="krw")
        a1.plot(s, scal.kro(s), color=PALETTE.oil, lw=2, label="kro")
        style_axes(a1, "Corey nisbi keçiriciliyi", "Sw", "kr")
        legend(a1)

        if capillary is not None:
            twin = a1.twinx()
            twin.plot(s, capillary.pcow(s), color=PALETTE.accent, lw=1.4, ls="-.")
            twin.set_ylabel("Pc, bar", color=PALETTE.accent, fontsize=9)
            twin.tick_params(colors=PALETTE.accent, labelsize=8)
            for spine in twin.spines.values():
                spine.set_color(PALETTE.line)

        lw = scal.krw(s) / mu_w
        lo = scal.kro(s) / mu_o
        fw = lw / np.maximum(lw + lo, 1e-30)
        a2.plot(s, fw, color=PALETTE.accent, lw=2)
        chord = (fw - fw[0]) / np.maximum(s - s[0], 1e-9)
        k = int(np.argmax(chord[1:])) + 1
        a2.plot([s[0], s[k]], [fw[0], fw[k]], color=PALETTE.oil, lw=1.2, ls="--")
        a2.plot(s[k], fw[k], "o", color=PALETTE.oil, ms=7)
        a2.annotate(f"Swf = {s[k]:.3f}", (s[k], fw[k]), xytext=(8, -14),
                    textcoords="offset points", color=PALETTE.oil, fontsize=9)
        style_axes(a2, "Fraksional axın + Welge toxunanı", "Sw", "fw")


class CrossSectionRenderer:
    """Şaquli kəsik — X (və ya Y) oxu boyunca dərinlik üzrə görüntü.

    3D modeldə cazibə seqreqasiyası və şaquli süpürmə YALNIZ burada
    görünür: areal xəritə hər təbəqəni ayrı-ayrı göstərir, kəsik isə
    təbəqələr arasındakı fərqi bir şəkildə birləşdirir.
    """

    def draw(self, ax, figure, model: ReservoirModel, property_key: str,
             snapshot: Optional[Snapshot] = None, colorbar=None,
             axis: str = "J", index: int = 0, cax=None):
        ax.clear()
        if cax is not None:
            cax.clear()
        elif colorbar is not None:
            try:
                colorbar.remove()
            except Exception:
                pass

        grid = model.grid
        data, cmap, vmin, vmax = MapRenderer()._select_volume(
            model, property_key, snapshot)

        if axis.upper() == "J":
            j = int(np.clip(index, 0, grid.ny - 1))
            section = data[:, j, :]
            horizontal = grid.nx * model.geometry.dx
            xlabel = f"X, m    (J = {j + 1})"
        else:
            i = int(np.clip(index, 0, grid.nx - 1))
            section = data[:, :, i]
            horizontal = grid.ny * model.geometry.dy
            xlabel = f"Y, m    (I = {i + 1})"

        depths = model.geometry.cell_depths().reshape(grid.shape)
        top = float(depths.min() - model.geometry.dz[0] * 0.5)
        base = float(depths.max() + model.geometry.dz[-1] * 0.5)

        image = ax.imshow(section, origin="upper", aspect="auto", cmap=cmap,
                          vmin=vmin, vmax=vmax, interpolation="nearest",
                          extent=[0, horizontal, base, top])
        bar = _colorbar(figure, image, ax, cax,
                        property_label(property_key))

        self._draw_well_tracks(ax, model, axis, index)
        t = snapshot.time if snapshot else 0.0
        style_axes(ax, f"Şaquli kəsik — {property_label(property_key)}"
                       f"    |    t = {t:.0f} gün", xlabel, "Dərinlik, m")
        return bar

    @staticmethod
    def _draw_well_tracks(ax, model: ReservoirModel, axis: str, index: int):
        grid = model.grid
        depths = model.geometry.cell_depths().reshape(grid.shape)
        for well in model.active_wells():
            perfs = well.open_perforations()
            if not perfs:
                continue
            if axis.upper() == "J":
                if perfs[0].j != int(np.clip(index, 0, grid.ny - 1)):
                    continue
                x = (perfs[0].i + 0.5) * model.geometry.dx
            else:
                if perfs[0].i != int(np.clip(index, 0, grid.nx - 1)):
                    continue
                x = (perfs[0].j + 0.5) * model.geometry.dy
            zs = [depths[p.k, p.j, p.i] for p in perfs]
            colour = (PALETTE.water if well.well_type is WellType.INJECTOR
                      else PALETTE.oil)
            ax.plot([x, x], [min(zs), max(zs)], color=colour, lw=3,
                    solid_capstyle="butt", zorder=5)
            ax.plot([x] * len(zs), zs, "o", ms=4, mfc="white", mec=colour, zorder=6)
            ax.annotate(well.name, (x, min(zs)), textcoords="offset points",
                        xytext=(6, -12), color="white", fontsize=8,
                        weight="bold", zorder=7)


class RunComparisonRenderer:
    """Bir neçə simulyasiya nəticəsinin üst-üstə müqayisəsi.

    B1-in əsas məqsədi budur: "bu dəyişiklik neft hasilatına necə təsir
    etdi" sualına cavab vermək üçün iki işə salınmanı yan-yana qoymaq.
    """

    PALETTE_CYCLE = (PALETTE.oil, PALETTE.water, PALETTE.accent,
                     "#B06FD0", "#8FBF3F", PALETTE.danger)

    def draw(self, axes, runs):
        """runs — (etiket, SimulationResult) cütlərinin siyahısı."""
        (a1, a2), (a3, a4) = axes
        for ax in (a1, a2, a3, a4):
            ax.clear()

        for index, (label, result) in enumerate(runs):
            series = result.series
            if not series.time:
                continue
            colour = self.PALETTE_CYCLE[index % len(self.PALETTE_CYCLE)]
            t = np.array(series.time)
            a1.plot(t, series.recovery_factor, color=colour, lw=2, label=label)
            a2.plot(t, series.water_cut, color=colour, lw=1.8, label=label)
            a3.plot(t, np.array(series.cumulative_oil) / 1e3, color=colour,
                    lw=1.8, label=label)
            a4.plot(t, series.oil_rate, color=colour, lw=1.5, label=label)

        style_axes(a1, "Recovery Factor", "Zaman, gün", "RF, %")
        style_axes(a2, "Sulaşma", "Zaman, gün", "WC, %")
        a2.set_ylim(0, 100)
        style_axes(a3, "Kumulyativ neft", "Zaman, gün", "min m³")
        style_axes(a4, "Neft debiti", "Zaman, gün", "qo, m³/gün")
        if runs:
            legend(a1)


class TornadoRenderer:
    """Tornado diaqramı — parametrlərin nəticəyə təsiri, azalan sırayla.

    Hər parametr üçün üfüqi zolaq: sol ucu aşağı hədddəki, sağ ucu
    yuxarı hədddəki çıxış dəyəri. Zolağın eni parametrin əhəmiyyətini
    göstərir — ən əhəmiyyətlisi yuxarıda (klassik tornado tərtibatı).
    """

    def draw(self, ax, report):
        ax.clear()
        items = report.sorted_by_swing() if report.items else []
        if not items:
            style_axes(ax, "Həssaslıq təhlili aparılmayıb", "", "")
            return

        names = [item.name for item in items]
        positions = np.arange(len(items))[::-1]
        baseline = report.baseline_output

        for position, item in zip(positions, items):
            left = min(item.low_output, item.high_output)
            width = item.swing
            colour = PALETTE.water if not item.direction_reversed else PALETTE.oil
            ax.barh(position, width, left=left, height=0.6,
                   color=colour, alpha=0.85,
                   edgecolor=(0, 0, 0, 0.3) if item.failed_low or item.failed_high
                   else "none",
                   hatch="//" if (item.failed_low or item.failed_high) else None)
            low_label = "✗" if item.failed_low else f"{item.low_output:.3g}"
            high_label = "✗" if item.failed_high else f"{item.high_output:.3g}"
            ax.text(left - width * 0.02, position, low_label, ha="right",
                   va="center", fontsize=8, color=PALETTE.text_dim)
            ax.text(left + width * 1.02, position, high_label, ha="left",
                   va="center", fontsize=8, color=PALETTE.text_dim)

        ax.axvline(baseline, color=PALETTE.text_dim, lw=1.2, ls="--",
                  alpha=0.7)
        ax.set_yticks(positions)
        ax.set_yticklabels(names, fontsize=9)
        style_axes(ax, f"Baza = {baseline:.4g}", report.metric_name, "")
        ax.set_ylim(-0.7, len(items) - 0.3)


class OptimisationRenderer:
    """Axtarışın gedişatı: konvergensiya əyrisi və parametr trayektoriyası."""

    def draw(self, axes_pair, match_result):
        left, right = axes_pair
        left.clear(); right.clear()
        history = getattr(match_result, "history", None)
        if not history:
            style_axes(left, "Axtarış aparılmayıb", "Qiymətləndirmə", "")
            style_axes(right, "", "", "")
            return

        iterations = np.arange(1, len(history) + 1)
        scores = np.array([item.mismatch for item in history])
        best = match_result.convergence_curve

        # cərimə dəyərləri qrafiki sıxır — onları ayrıca göstəririk
        failed = np.array([not item.succeeded for item in history])
        finite = scores[~failed]
        ceiling = finite.max() * 1.5 if finite.size else 1.0

        left.plot(iterations[~failed], np.minimum(scores[~failed], ceiling),
                  "o", ms=3.5, color=PALETTE.text_dim, alpha=0.6,
                  label="qiymətləndirmə")
        if failed.any():
            left.plot(iterations[failed], np.full(failed.sum(), ceiling),
                      "x", ms=5, color=PALETTE.danger, label="uğursuz")
        left.plot(iterations, np.minimum(best, ceiling), "-", lw=2,
                  color=PALETTE.accent, label="ən yaxşı")
        left.set_yscale("log")
        style_axes(left, "Konvergensiya", "Qiymətləndirmə", "Uyğunsuzluq")
        legend(left, fontsize=8)

        colours = RunComparisonRenderer.PALETTE_CYCLE
        for index, definition in enumerate(match_result.parameters.definitions):
            track = np.array([item.unit_values[index] for item in history])
            right.plot(iterations, track, lw=1.3,
                       color=colours[index % len(colours)],
                       label=definition.name)
        right.set_ylim(-0.05, 1.05)
        style_axes(right, "Parametrlər (normallaşdırılmış)",
                   "Qiymətləndirmə", "[0, 1]")
        legend(right, fontsize=8)


def MismatchReportPlaceholder():
    """Boş hesabat — müşahidə yüklənməyəndə qrafikləri təmizləmək üçün."""
    from ..history.mismatch import MismatchReport
    return MismatchReport()


class HistoryMatchRenderer:
    """Müşahidə və model əyrilərinin üst-üstə müqayisəsi.

    Müşahidə NÖQTƏ ilə, model XƏTT ilə çəkilir — bu, sənaye
    konvensiyasıdır və hansının ölçmə, hansının hesablama olduğunu
    bir baxışda göstərir.
    """

    def draw(self, figure, axes_list, report):
        for axes in axes_list:
            axes.clear()
        if not report.series:
            style_axes(axes_list[0], "Müşahidə tapılmadı", "Zaman, gün", "")
            for axes in axes_list[1:]:
                axes.set_visible(False)
            return

        ranked = sorted(report.series, key=lambda item: -item.nrmse)
        for axes, item in zip(axes_list, ranked):
            axes.set_visible(True)
            axes.plot(item.time, item.observed, "o", ms=5,
                      mfc="none", mec=PALETTE.oil, mew=1.4, label="ölçülmüş")
            axes.plot(item.time, item.simulated, "-", lw=1.8,
                      color=PALETTE.accent, label="model")
            style_axes(axes, f"{item.label}   NRMSE {item.nrmse:.3f}",
                       "Zaman, gün", item.quantity.unit)
            legend(axes, fontsize=8)

        for axes in axes_list[len(ranked):]:
            axes.set_visible(False)


class CrossPlotRenderer:
    """Ölçülmüş–hesablanmış çarpaz qrafiki.

    İdeal uyğunluqda bütün nöqtələr 45° xəttinin üzərindədir. Səpələnmə
    təsadüfi xətanı, xəttdən sistematik meyl isə model səhvini göstərir —
    zaman qrafikində bu ikisini ayırmaq çətindir.
    """

    def draw(self, axes, report):
        axes.clear()
        if not report.series:
            style_axes(axes, "Müşahidə tapılmadı", "", "")
            return

        colours = RunComparisonRenderer.PALETTE_CYCLE
        low = high = None
        for index, item in enumerate(report.series):
            # hər sıra öz miqyasında normallaşdırılır ki, bir qrafikə sığsın
            scale = item.scale
            observed = item.observed / scale
            simulated = item.simulated / scale
            axes.plot(observed, simulated, "o", ms=4, alpha=0.75,
                      color=colours[index % len(colours)], label=item.label)
            values = np.concatenate([observed, simulated])
            low = values.min() if low is None else min(low, values.min())
            high = values.max() if high is None else max(high, values.max())

        margin = 0.05 * max(high - low, 1e-9)
        limits = (low - margin, high + margin)
        axes.plot(limits, limits, "--", lw=1.2, color=PALETTE.text_dim)
        axes.set_xlim(limits)
        axes.set_ylim(limits)
        style_axes(axes, "Ölçülmüş — hesablanmış (normallaşdırılmış)",
                   "ölçülmüş", "hesablanmış")
        legend(axes, fontsize=8)


class PvtRenderer:
    """PVT əyriləri: Bo, μo, Rs, Bw/μw. Doyma təzyiqi şaquli xətlə."""

    def draw(self, axes, table):
        (a1, a2), (a3, a4) = axes
        for ax in (a1, a2, a3, a4):
            ax.clear()
        p = table.pressure
        pb = table.bubble_point

        a1.plot(p, table.oil_fvf, color=PALETTE.oil, lw=2)
        style_axes(a1, "Neftin həcm əmsalı", "Təzyiq, bar", "Bo, rm³/sm³")

        a2.plot(p, table.oil_viscosity, color=PALETTE.oil, lw=2)
        style_axes(a2, "Neftin lözlüyü", "Təzyiq, bar", "μo, cP")

        a3.plot(p, table.solution_gor, color=PALETTE.accent, lw=2)
        style_axes(a3, "Həll olmuş qaz nisbəti", "Təzyiq, bar", "Rs, sm³/sm³")

        a4.plot(p, table.water_fvf, color=PALETTE.water, lw=2, label="Bw")
        twin = a4.twinx()
        twin.plot(p, table.water_viscosity, color=PALETTE.text_dim, lw=1.2, ls=":")
        twin.set_ylabel("μw, cP", color=PALETTE.text_dim, fontsize=9)
        twin.tick_params(colors=PALETTE.text_dim, labelsize=8)
        for spine in twin.spines.values():
            spine.set_color(PALETTE.line)
        style_axes(a4, "Su xassələri", "Təzyiq, bar", "Bw, rm³/sm³")
        legend(a4)

        for ax in (a1, a2, a3):
            if p[0] <= pb <= p[-1]:
                ax.axvline(pb, color=PALETTE.danger, ls="--", lw=1)
                ax.annotate(f"Pb = {pb:.0f} bar", (pb, ax.get_ylim()[1]),
                            xytext=(6, -14), textcoords="offset points",
                            color=PALETTE.danger, fontsize=8)


class ValidationRenderer:
    """Analitik və ədədi həllin müqayisəsi."""

    def draw(self, ax, analytical, x_numeric, sw_numeric, time, ncell):
        ax.clear()
        ax.plot(analytical.distance, analytical.water_saturation,
                color=PALETTE.oil, lw=2.4, label="Analitik (Bukley-Leverett)")
        ax.step(x_numeric, sw_numeric, where="mid", color=PALETTE.water, lw=1.6,
                label=f"Ədədi IMPES ({ncell} hüceyrə)")
        ax.axvline(analytical.front_position, color=PALETTE.text_dim, ls="--", lw=1)
        ax.set_xlim(0, min(float(x_numeric[-1]), analytical.front_position * 2.2))
        style_axes(ax, f"1D sıxışdırma yoxlaması,  t = {time:.0f} gün",
                   "Məsafə, m", "Su doyumluluğu Sw")
        legend(ax, fontsize=9)
