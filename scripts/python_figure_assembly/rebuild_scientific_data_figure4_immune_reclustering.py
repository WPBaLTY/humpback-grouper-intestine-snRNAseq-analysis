from pathlib import Path
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from PIL import Image


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_ROOT = REPO_ROOT.parent / "zenodo_processed_data_record" if (REPO_ROOT.parent / "zenodo_processed_data_record").is_dir() else REPO_ROOT
SOURCE_DIR = DATA_ROOT / "source_data" / "figures"
OUT_DIR = REPO_ROOT / "outputs" / "figures" / "scientific_data"

UMAP_SOURCE = SOURCE_DIR / "Fig4_immune_umap_coordinates.csv"
DOT_SOURCE = SOURCE_DIR / "Fig4_immune_marker_dotplot_source.csv"

PNG_OUT = OUT_DIR / "Figure4_immune_reclustering_final.png"
PDF_OUT = OUT_DIR / "Figure4_immune_reclustering_final.pdf"
TIFF_OUT = OUT_DIR / "Figure4_immune_reclustering_final.tiff"
WRITE_TIFF = os.environ.get("GROUPER_WRITE_TIFF", "0") == "1"


SUBTYPE_ORDER = [
    "T cell (CCR7+)",
    "Activated T (RORA+)",
    "Activated lymphoid (CCL20hi)",
    "NK-like cytotoxic",
    "B cells",
    "cDC1-like (XCR1+)",
    "Monocytes/macrophages",
    "MoDC-like (CD209d+)",
    "Granulocyte-like",
    "Cycling (G2/M)",
]

PALETTE = {
    "T cell (CCR7+)": "#55B8AF",
    "Activated T (RORA+)": "#D8B900",
    "Activated lymphoid (CCL20hi)": "#8E7CC3",
    "NK-like cytotoxic": "#E64B35",
    "B cells": "#3B9BC2",
    "cDC1-like (XCR1+)": "#F39B34",
    "Monocytes/macrophages": "#76A82A",
    "MoDC-like (CD209d+)": "#CC5C99",
    "Granulocyte-like": "#999999",
    "Cycling (G2/M)": "#9B4AA3",
}

GENE_ORDER = [
    "CCR7",
    "BCL11B",
    "LCK",
    "roraa",
    "FKBP5",
    "ddit4.1",
    "Ccl20",
    "CCR6",
    "Prf1",
    "Gzmb.1",
    "TYROBP",
    "FCER1G",
    "XCR1",
    "ZNF366",
    "csf1r1",
    "Axl",
    "CMKLR1",
    "Cd209d.1",
    "FN1",
    "AOC3",
    "EPX",
    "Ncf4",
    "CYBB.1",
    "CD79A",
    "EBF1",
    "BLNK",
]


def draw_umap(ax, df):
    for subtype in SUBTYPE_ORDER:
        sub = df[df["immune_subtype"] == subtype]
        ax.scatter(
            sub["UMAP_1"],
            sub["UMAP_2"],
            s=1.1,
            c=PALETTE[subtype],
            label=subtype,
            rasterized=True,
            linewidths=0,
            alpha=0.95,
        )

    ax.set_title("Immune-lineage UMAP", loc="left", fontsize=17, fontweight="bold", pad=8)
    ax.set_xlabel("UMAP 1", fontsize=13)
    ax.set_ylabel("UMAP 2", fontsize=13)
    ax.set_xlim(-10.4, 12.6)
    ax.set_ylim(-13.2, 13.6)
    ax.set_xticks([-5, 0, 5, 10])
    ax.set_yticks([-10, -5, 0, 5, 10])
    ax.tick_params(labelsize=11, width=1.1, length=4)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)

    label_pos = {
        "T cell (CCR7+)": (-3.1, -0.8),
        "Activated T (RORA+)": (-2.6, -5.2),
        "Activated lymphoid (CCL20hi)": (-6.2, -3.5),
        "NK-like cytotoxic": (-4.8, 10.8),
        "B cells": (-0.9, -9.8),
        "cDC1-like (XCR1+)": (3.2, -11.4),
        "Monocytes/macrophages": (8.6, 2.6),
        "MoDC-like (CD209d+)": (9.4, 7.1),
        "Granulocyte-like": (-1.0, 1.4),
        "Cycling (G2/M)": (3.6, -4.8),
    }
    for subtype, (x, y) in label_pos.items():
        ax.text(
            x,
            y,
            subtype,
            ha="center",
            va="center",
            fontsize=11,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.65, pad=0.15),
        )

    handles = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=PALETTE[subtype], markeredgecolor="none", markersize=7, label=subtype)
        for subtype in SUBTYPE_ORDER
    ]
    ax.legend(
        handles=handles,
        title="Immune subtype",
        frameon=False,
        bbox_to_anchor=(0.99, 0.88),
        loc="upper left",
        fontsize=9,
        title_fontsize=10,
        labelspacing=0.42,
        handletextpad=0.4,
        borderaxespad=0,
    )


def draw_dotplot(ax, cax, legend_ax, dot):
    dot = dot.copy()
    dot["immune_subtype"] = pd.Categorical(dot["immune_subtype"], categories=SUBTYPE_ORDER, ordered=True)
    dot["marker_gene"] = pd.Categorical(dot["marker_gene"], categories=GENE_ORDER, ordered=True)
    dot = dot.dropna(subset=["immune_subtype", "marker_gene"])

    x = dot["marker_gene"].cat.codes
    y = dot["immune_subtype"].cat.codes
    sizes = 2 + 135 * (dot["percent_expressed"].clip(0, 75) / 75.0) ** 1.35
    vals = dot["average_expression_scaled"].clip(-1.5, 2.5)

    cmap = LinearSegmentedColormap.from_list("dotplot_blue", ["#D2D2D2", "#B79CE8", "#7A55E8", "#0000FF"])
    sc = ax.scatter(x, y, s=sizes, c=vals, cmap=cmap, vmin=-1.5, vmax=2.5, linewidths=0)
    ax.set_title("Immune-subtype marker-gene dot plot", loc="left", fontsize=16, fontweight="bold", pad=8)
    ax.set_xlabel("Marker genes", fontsize=13, labelpad=1)
    ax.set_ylabel("Immune subtype", fontsize=13)
    ax.set_xticks(range(len(GENE_ORDER)))
    ax.set_xticklabels(GENE_ORDER, rotation=50, ha="right", fontsize=9.3, fontstyle="italic")
    ax.set_yticks(range(len(SUBTYPE_ORDER)))
    ax.set_yticklabels(SUBTYPE_ORDER, fontsize=11)
    ax.set_xlim(-0.7, len(GENE_ORDER) - 0.3)
    ax.set_ylim(-0.6, len(SUBTYPE_ORDER) - 0.4)
    ax.grid(False)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    cbar = plt.colorbar(sc, cax=cax)
    cbar.ax.set_title("Scaled average\nexpression", fontsize=11, pad=8)
    cbar.ax.tick_params(labelsize=9, width=0.8, length=3)

    legend_ax.axis("off")
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.text(0.0, 0.98, "Percent Expressed", fontsize=10, ha="left", va="top")
    for yv, label, size in zip([0.73, 0.52, 0.31, 0.10], ["0", "25", "50", "75"], [4, 24, 58, 95]):
        legend_ax.scatter([0.18], [yv], s=size, color="black", linewidths=0)
        legend_ax.text(0.35, yv, label, fontsize=10, va="center", ha="left")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["pdf.fonttype"] = 42

    umap = pd.read_csv(UMAP_SOURCE)
    dot = pd.read_csv(DOT_SOURCE)

    fig = plt.figure(figsize=(11.5, 8.15), dpi=600, facecolor="white")

    ax_a = fig.add_axes([0.065, 0.535, 0.720, 0.405])
    draw_umap(ax_a, umap)
    ax_a.text(-0.055, 1.07, "a", transform=ax_a.transAxes, fontsize=16, fontweight="bold", va="top")

    ax_b = fig.add_axes([0.205, 0.120, 0.665, 0.310])
    size_ax = fig.add_axes([0.875, 0.305, 0.120, 0.115])
    cax = fig.add_axes([0.912, 0.135, 0.022, 0.095])
    draw_dotplot(ax_b, cax, size_ax, dot)
    ax_b.text(-0.180, 1.10, "b", transform=ax_b.transAxes, fontsize=16, fontweight="bold", va="top")

    fig.savefig(PNG_OUT, dpi=600, facecolor="white")
    fig.savefig(
        PDF_OUT,
        dpi=600,
        facecolor="white",
        metadata={"CreationDate": None},
    )
    plt.close(fig)

    if WRITE_TIFF:
        Image.open(PNG_OUT).convert("RGB").save(TIFF_OUT, dpi=(600, 600))
    print(PNG_OUT)


if __name__ == "__main__":
    main()
