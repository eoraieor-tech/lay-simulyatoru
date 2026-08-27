"""Vizual tema — rənglər və oxların formatlanması.

Rendering qatına aiddir, UI-yə yox: eyni palitra həm Qt interfeysində,
həm də başsız (headless) hesabat generasiyasında işlədilə bilər.
"""

from __future__ import annotations
from dataclasses import dataclass

from matplotlib.colors import LinearSegmentedColormap


@dataclass(frozen=True)
class Palette:
    background: str = "#141A1F"
    panel: str = "#1C242B"
    panel_alt: str = "#232D36"
    line: str = "#33414D"
    text: str = "#DCE4EA"
    text_dim: str = "#8A9AA6"
    oil: str = "#D98E2B"
    water: str = "#2AA7A0"
    accent: str = "#4B9FD6"
    danger: str = "#C0574B"


PALETTE = Palette()

SATURATION_CMAP = LinearSegmentedColormap.from_list(
    "sw_oil_water",
    ["#7A4A10", PALETTE.oil, "#B9A98C", "#5FBFB8", PALETTE.water, "#0E5F5C"])

PRESSURE_CMAP = "inferno"
PERMEABILITY_CMAP = "viridis"
POROSITY_CMAP = "cividis"


def style_axes(ax, title="", xlabel="", ylabel="", palette: Palette = PALETTE):
    ax.set_facecolor(palette.panel)
    for spine in ax.spines.values():
        spine.set_color(palette.line)
    ax.tick_params(colors=palette.text_dim, labelsize=8)
    ax.grid(True, color=palette.line, lw=0.5, alpha=0.55)
    ax.set_axisbelow(True)
    if title:
        ax.set_title(title, color=palette.text, fontsize=10, pad=8)
    if xlabel:
        ax.set_xlabel(xlabel, color=palette.text_dim, fontsize=9)
    if ylabel:
        ax.set_ylabel(ylabel, color=palette.text_dim, fontsize=9)


def legend(ax, palette: Palette = PALETTE, fontsize=8):
    ax.legend(fontsize=fontsize, facecolor=palette.panel_alt,
              edgecolor=palette.line, labelcolor=palette.text)
