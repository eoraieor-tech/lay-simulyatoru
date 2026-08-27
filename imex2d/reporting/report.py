"""PDF hesabat generatoru — B6.

Yeni asılılıq YOXDUR: matplotlib onsuz da tələb olunur, `PdfPages`
onun öz bir hissəsidir. Mövcud renderer-lər (`rendering/renderers.py`)
birbaşa işlədilir — hesabat qrafikləri interfeysdəki qrafiklərlə
EYNİ koddan çıxır, ikinci bir "hesabat üçün" versiyası yoxdur.

Struktur — hər bölmə öz səhifəsi (və ya səhifə qrupu):

    1. Başlıq        model adı, tarix, versiya
    2. Model xülasəsi  grid, süxur, flüid, SCAL, PVT, quyular, fay
    3. Diaqnostika     xəta/xəbərdarlıq siyahısı (varsa)
    4. Areal xəritə    məsaməlilik + keçiricilik
    5. SCAL əyriləri   (əgər model quruludursa)
    6. PVT əyriləri    (əgər cədvəl varsa)
    7. Nəticələr       debit/RF/WC qrafikləri (əgər simulyasiya varsa)
    8. Tarixçə uyğunluğu (əgər müşahidə varsa)

Hər bölmə İSTƏYƏ GÖRƏDİR (`None` verilsə keçilir) — hesabat model
tək başına da, tam nəticə dəstiylə də yaradıla bilər.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Optional

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure

from ..domain.reservoir_model import ReservoirModel
from ..history.mismatch import MismatchReport
from ..history.optimizer import MatchResult
from ..history.sensitivity import SensitivityReport
from ..logging_setup import get_logger
from ..rendering import renderers as R
from ..rendering.theme import PALETTE
from ..simulation.results import SimulationResult
from ..version import VERSION

LOG = get_logger(__name__)

PAGE_SIZE = (11.7, 8.3)      # A4 üfüqi, düym


@dataclass
class ReportSections:
    """Hansı bölmələr daxil edilsin — hamısı defolt aktivdir."""
    summary: bool = True
    diagnostics: bool = True
    maps: bool = True
    scal: bool = True
    pvt: bool = True
    results: bool = True
    history_match: bool = True


@dataclass
class ReportContext:
    """Hesabata daxil olan bütün material — nəyi versən, o görünür."""
    model: ReservoirModel
    result: Optional[SimulationResult] = None
    mismatch: Optional[MismatchReport] = None
    match_result: Optional[MatchResult] = None
    sensitivity: Optional[SensitivityReport] = None
    sections: ReportSections = field(default_factory=ReportSections)
    author: str = ""


class ReportGenerator:
    """Model və nəticələrdən çoxsəhifəli PDF qurur.

    Vəziyyətsizdir (stateless) — bölmə seçimi `ReportContext.sections`
    üzərindədir, ona görə eyni generator müxtəlif konfiqurasiyalı
    hesabatlar üçün təkrar işlədilə bilər.
    """

    def write(self, context: ReportContext, path: str) -> str:
        directory = os.path.dirname(os.path.abspath(path))
        if directory:
            os.makedirs(directory, exist_ok=True)

        pages = 0
        with PdfPages(path) as pdf:
            pages += self._title_page(pdf, context)
            if context.sections.summary:
                pages += self._summary_page(pdf, context)
            if context.sections.diagnostics:
                pages += self._diagnostics_page(pdf, context)
            if context.sections.maps:
                pages += self._map_pages(pdf, context)
            if context.sections.scal:
                pages += self._scal_page(pdf, context)
            if context.sections.pvt and context.model.pvt_table is not None:
                pages += self._pvt_page(pdf, context)
            if context.sections.results and context.result is not None:
                pages += self._results_page(pdf, context)
            if (context.sections.history_match
                    and context.mismatch is not None and context.mismatch.series):
                pages += self._history_match_page(pdf, context)

            info = pdf.infodict()
            info["Title"] = f"{context.model.name} — IMEX-2D hesabatı"
            info["Author"] = context.author or "IMEX-2D"
            info["Subject"] = "Rezervuar simulyasiyası hesabatı"
            info["CreationDate"] = datetime.now()

        LOG.info("PDF hesabat yazıldı: %s (%d səhifə)", path, pages)
        return path

    # ═══════════════════════════════════════════════════════ səhifələr
    def _new_figure(self) -> Figure:
        return Figure(figsize=PAGE_SIZE, facecolor="white")

    def _save(self, pdf: PdfPages, figure: Figure, layout: bool = True) -> int:
        self._print_friendly(figure)
        if layout:
            figure.tight_layout()
        pdf.savefig(figure, facecolor="white")
        return 1

    @staticmethod
    def _print_friendly(figure: Figure) -> None:
        """Renderer-lərin tünd UI teması ağ səhifədə oxunmur — məcburi düzəliş.

        `rendering/renderers.py`-dəki bütün çəkmə funksiyaları
        `PALETTE.text` kimi AÇIQ rəngləri TÜND fon üçün işlədir
        (interfeysin öz teması). Həmin qrafiki birbaşa ağ PDF
        səhifəsinə köçürsək, başlıqlar və oxlar demək olar görünməz
        olur. Hər renderer-i palitra qəbul edəcək şəkildə dəyişmək
        əvəzinə, çəkildikdən SONRA bütün oxları məcburi işıqlandırırıq
        — bu, hər renderer-in daxili detallarından asılı olmur və yeni
        renderer əlavə olunanda unudulma riski yaratmır.
        """
        DARK = "#111827"
        MID = "#4B5563"
        LIGHT_GRID = "#E5E7EB"

        for axes in figure.get_axes():
            axes.set_facecolor("white")
            for spine in axes.spines.values():
                spine.set_color(MID)
            axes.tick_params(colors=MID, labelsize=8)
            axes.grid(True, color=LIGHT_GRID, lw=0.6)
            if axes.get_title():
                axes.title.set_color(DARK)
            axes.xaxis.label.set_color(MID)
            axes.yaxis.label.set_color(MID)

            legend = axes.get_legend()
            if legend is not None:
                legend.get_frame().set_facecolor("white")
                legend.get_frame().set_edgecolor(LIGHT_GRID)
                for text in legend.get_texts():
                    text.set_color(DARK)

        if figure._suptitle is not None:
            figure._suptitle.set_color(DARK)

    def _title_page(self, pdf: PdfPages, context: ReportContext) -> int:
        figure = self._new_figure()
        axes = figure.add_axes([0, 0, 1, 1])
        axes.axis("off")
        stamp = datetime.now().strftime("%Y-%m-%d %H:%M")

        axes.text(0.5, 0.62, context.model.name, ha="center", fontsize=26,
                  weight="bold", color="#111827")
        axes.text(0.5, 0.54, "Rezervuar simulyasiyası hesabatı",
                  ha="center", fontsize=15, color="#374151")
        axes.text(0.5, 0.44, f"IMEX-2D v{VERSION}   ·   {stamp}", ha="center",
                  fontsize=11, color="#6B7280")
        if context.author:
            axes.text(0.5, 0.39, f"Hazırladı: {context.author}", ha="center",
                      fontsize=10, color="#6B7280")

        grid = context.model.grid
        highlights = [
            f"Grid: {grid.nx} × {grid.ny} × {grid.nz}  =  {context.model.ncell} hüceyrə",
            f"Quyu sayı: {len(context.model.active_wells())}",
        ]
        if context.result is not None:
            highlights.append(
                f"Recovery Factor: {context.result.final_recovery_factor:.2f} %")
        for index, line in enumerate(highlights):
            axes.text(0.5, 0.28 - index * 0.045, line, ha="center",
                      fontsize=11, color="#111827")
        return self._save(pdf, figure, layout=False)

    def _summary_page(self, pdf: PdfPages, context: ReportContext) -> int:
        model = context.model
        figure = self._new_figure()
        axes = figure.add_axes([0.06, 0.05, 0.9, 0.9])
        axes.axis("off")
        axes.text(0, 1.0, "Model xülasəsi", fontsize=17, weight="bold",
                  color="#111827")

        lines = self._summary_lines(model)
        text = "\n".join(lines)
        axes.text(0, 0.93, text, va="top", fontsize=10.5, color="#1F2937",
                  family="monospace", linespacing=1.6)
        return self._save(pdf, figure, layout=False)

    @staticmethod
    def _summary_lines(model: ReservoirModel) -> List[str]:
        grid, geometry = model.grid, model.geometry
        rock = model.rock
        lines = [
            "GRID",
            f"  Ölçü            {grid.nx} x {grid.ny} x {grid.nz}  "
            f"({model.ncell} hüceyrə)",
            f"  Hüceyrə ölçüsü  {geometry.dx:g} x {geometry.dy:g} x "
            f"{geometry.dz:g} m",
            f"  Dərinlik        {float(np.min(geometry.cell_depths())):.0f} – "
            f"{float(np.max(geometry.cell_depths())):.0f} m",
            "",
            "SÜXUR",
            f"  Orta məsaməlilik  {float(np.mean(rock.porosity.values)):.3f}",
            f"  Orta Kx           {float(np.mean(rock.permx.values)):.1f} mD",
            f"  Süxur sıxılması   {rock.compressibility:.2e}  1/bar",
            "",
            "FLÜİD",
            f"  Neft sıxlığı      {model.fluids.oil_density:.0f} kg/m³",
            f"  Su sıxlığı        {model.fluids.water_density:.0f} kg/m³",
            f"  PVT cədvəli       {'var (' + str(model.pvt_table.pressure.size) + ' sətir)' if model.pvt_table else 'yoxdur (sabit)'}",
            "",
            "SCAL",
            f"  Mənbə             {'laboratoriya cədvəli (' + str(len(model.scal_tables)) + ' region)' if model.scal_tables else 'Corey düsturu'}",
            f"  Swc / Sor         {model.scal_parameters.swc:.3f} / "
            f"{model.scal_parameters.sor:.3f}",
            f"  Kapilyar təzyiq   {'aktiv' if model.capillary_parameters.enabled else 'söndürülüb'}",
            "",
            "QUYULAR",
        ]
        for well in model.active_wells():
            perforations = well.open_perforations()
            k_range = (f"K{perforations[0].k + 1}-{perforations[-1].k + 1}"
                      if perforations else "—")
            lines.append(
                f"  {well.name:<10} {well.well_type.name:<9} "
                f"{well.control.mode.name:<5} {well.control.target:>8.1f}  {k_range}")
        lines += ["", "FAULTS"]
        if model.fault_references:
            for fault in model.fault_references:
                lines.append(f"  {fault.summary()}")
        else:
            lines.append("  Fault yoxdur.")
        return lines

    def _diagnostics_page(self, pdf: PdfPages, context: ReportContext) -> int:
        report = context.model.diagnose()
        if not report.items:
            return 0
        figure = self._new_figure()
        axes = figure.add_axes([0.06, 0.05, 0.9, 0.9])
        axes.axis("off")
        axes.text(0, 1.0, "Diaqnostika", fontsize=17, weight="bold",
                  color="#111827")

        lines = []
        for item in report.errors:
            lines.append(f"[XƏTA] {item.source}: {item.message}")
        for item in report.warnings:
            lines.append(f"[XƏBƏRDARLIQ] {item.source}: {item.message}")
        axes.text(0, 0.93, "\n".join(lines), va="top", fontsize=10.5,
                  color="#1F2937", family="monospace", linespacing=1.7)
        return self._save(pdf, figure, layout=False)

    def _map_pages(self, pdf: PdfPages, context: ReportContext) -> int:
        model = context.model
        keys = [key for key in (R.POROSITY, R.PERMEABILITY)
               if key is not None]
        if not keys:
            return 0
        figure = self._new_figure()
        axes_list = figure.subplots(1, len(keys))
        if len(keys) == 1:
            axes_list = [axes_list]
        renderer = R.MapRenderer()
        for axes, key in zip(axes_list, keys):
            renderer.draw(axes, figure, model, key)
        figure.suptitle("Areal xəritələr", fontsize=14, weight="bold")
        return self._save(pdf, figure)

    def _scal_page(self, pdf: PdfPages, context: ReportContext) -> int:
        model = context.model
        if model.scal_tables:
            return 0            # cədvəl əsaslı SCAL-ın öz qrafiki fərqlidir, B4-ün əhatəsi xaricində
        figure = self._new_figure()
        axes = figure.subplots(1, 2)
        R.ScalRenderer().draw(
            axes, model.scal_parameters, model.fluids.water_viscosity,
            model.fluids.oil_viscosity,
            model.capillary_parameters if model.capillary_parameters.enabled
            else None)
        figure.suptitle("Nisbi keçiricilik", fontsize=14, weight="bold")
        return self._save(pdf, figure)

    def _pvt_page(self, pdf: PdfPages, context: ReportContext) -> int:
        figure = self._new_figure()
        axes = figure.subplots(2, 2)
        R.PvtRenderer().draw(axes, context.model.pvt_table)
        figure.suptitle("PVT xassələri", fontsize=14, weight="bold")
        return self._save(pdf, figure)

    def _results_page(self, pdf: PdfPages, context: ReportContext) -> int:
        figure = self._new_figure()
        axes = figure.subplots(2, 2)
        R.ProductionCurveRenderer().draw(axes, context.result)
        figure.suptitle(
            f"Simulyasiya nəticələri — RF "
            f"{context.result.final_recovery_factor:.2f} %",
            fontsize=14, weight="bold")
        return self._save(pdf, figure)

    def _history_match_page(self, pdf: PdfPages, context: ReportContext) -> int:
        figure = self._new_figure()
        axes = figure.subplots(2, 2)
        R.HistoryMatchRenderer().draw(figure, list(np.ravel(axes)),
                                      context.mismatch)
        figure.suptitle(
            f"Tarixçə uyğunluğu — yekun {context.mismatch.total:.4f}",
            fontsize=14, weight="bold")
        return self._save(pdf, figure)
