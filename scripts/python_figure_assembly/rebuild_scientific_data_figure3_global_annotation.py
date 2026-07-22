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

CLUSTER_UMAP = SOURCE_DIR / "Fig3A_transcriptional_cluster_umap_source.png"
BROAD_UMAP = SOURCE_DIR / "Fig3B_broad_category_umap_source.png"
DOT_SOURCE = SOURCE_DIR / "Fig3_global_marker_dotplot_source.csv"

PNG_OUT = OUT_DIR / "Figure3_global_annotation_final.png"
PDF_OUT = OUT_DIR / "Figure3_global_annotation_final.pdf"
TIFF_OUT = OUT_DIR / "Figure3_global_annotation_final.tiff"
WRITE_TIFF = os.environ.get("GROUPER_WRITE_TIFF", "0") == "1"

UMAP_EXTENT = (-16, 10, -15.5, 12.5)

DOT_ROWS = [
    ("Enterocytes", "Enterocytes"),
    ("LREs", "LREs"),
    ("Goblet cells", "Goblet cells"),
    ("Tuft-like cells", "Tuft cells"),
    ("Best4+ cells", "Best4+ cell"),
    ("Enteroendocrine cells", "Enteroendocrine cells"),
    ("Neuronal cells", "Neuronal cells"),
    ("Leukocytes", "Leukocytes"),
    ("Fibroblasts", "Fibroblasts"),
    ("Endothelial cells", "Endothelial cells"),
    ("Smooth muscle cells", "Smooth muscle"),
    ("Acinar-like cells", "Acinar-like"),
]

DOT_GENES = [
    "fabp2", "cd36", "SI",
    "slc10a2", "CUBN", "LRP2",
    "muc2", "spdef", "FER1L6",
    "pou2f3", "avil", "Pik3ap1",
    "best4", "cftr", "slc20a1a",
    "neurod1", "scgn", "ISL1",
    "syt1", "elavl3", "phox2a",
    "ptprc", "lcp1", "BCL11B", "SATB1", "FCER1G", "XCR1", "Axl",
    "col1a1", "col1a2", "dcn",
    "pecam1", "cdh5", "kdrl",
    "tagln", "CNTNAP5", "Pld5.1",
    "cel", "Cela1.1", "cpa2",
]

def strip_dark_text(im: Image.Image) -> Image.Image:
    arr = np.asarray(im.convert("RGB")).copy()
    mean = arr.mean(axis=2)
    chroma = arr.max(axis=2) - arr.min(axis=2)
    mask = (mean < 118) & (chroma < 58)
    expanded = mask.copy()
    for dy in (-2, -1, 0, 1, 2):
        for dx in (-2, -1, 0, 1, 2):
            expanded |= np.roll(np.roll(mask, dy, axis=0), dx, axis=1)
    arr[expanded] = 255
    return Image.fromarray(arr)


def clean_umap_layer(path: Path, crop_box, strip_text=False) -> Image.Image:
    layer = Image.open(path).convert("RGB").crop(crop_box)
    return strip_dark_text(layer) if strip_text else layer


def source_to_umap(x, y, crop_box):
    x0, y0, x1, y1 = crop_box
    xmin, xmax, ymin, ymax = UMAP_EXTENT
    ux = xmin + (x - x0) / (x1 - x0) * (xmax - xmin)
    uy = ymax - (y - y0) / (y1 - y0) * (ymax - ymin)
    return ux, uy


def draw_clean_umap(ax, path, crop_box, title, labels=None, legend=None, strip_text=False, legend_anchor=(1.02, 0.5)):
    layer = clean_umap_layer(path, crop_box, strip_text=strip_text)
    xmin, xmax, ymin, ymax = UMAP_EXTENT
    ax.imshow(layer, extent=[xmin, xmax, ymin, ymax], origin="upper", aspect="auto", interpolation="nearest")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xticks([-10, 0, 10])
    ax.set_yticks([-15, -10, -5, 0, 5, 10])
    ax.set_xlabel("UMAP 1", fontsize=11)
    ax.set_ylabel("UMAP 2", fontsize=11)
    ax.set_title(title, loc="left", fontsize=13, fontweight="bold", pad=8)
    ax.tick_params(axis="both", labelsize=10, width=1.0, length=4)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)

    for item in labels or []:
        text, sx, sy, bw, bh = item
        ux, uy = source_to_umap(sx, sy, crop_box)
        ax.text(
            ux,
            uy,
            text,
            ha="center",
            va="center",
            fontsize=10,
            color="black",
            zorder=4,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.92, pad=0.25),
        )

    if legend:
        handles = [
            Line2D([0], [0], marker="o", linestyle="", markerfacecolor=color, markeredgecolor="none", markersize=6.0, label=label)
            for color, label in legend
        ]
        ax.legend(
            handles=handles,
            title="Cell type",
            frameon=False,
            bbox_to_anchor=legend_anchor,
            loc="center left",
            fontsize=9,
            title_fontsize=10,
            borderaxespad=0,
            handletextpad=0.55,
            labelspacing=0.45,
        )


def draw_dotplot(ax, cax, legend_ax):
    dot = pd.read_csv(DOT_SOURCE)
    dot["marker_gene_l"] = dot["marker_gene"].astype(str).str.lower()

    expected_rows = len(DOT_ROWS) * len(DOT_GENES)
    if len(dot) != expected_rows or dot.duplicated(["cell_type", "marker_gene_l"]).any():
        raise ValueError("Figure 3 dot-plot source must contain one row per cell type and marker")

    xs, ys, sizes, vals = [], [], [], []
    for yi, (_, source_label) in enumerate(DOT_ROWS):
        sub = dot[dot["cell_type"] == source_label]
        for xi, gene in enumerate(DOT_GENES):
            gsub = sub[sub["marker_gene_l"] == gene.lower()]
            if len(gsub) != 1:
                raise ValueError(f"Missing or duplicated Figure 3 dot value: {source_label} / {gene}")
            value = gsub.iloc[0]
            pct = float(value["percent_expressed"])
            val = float(value["average_expression_scaled"])
            xs.append(xi)
            ys.append(yi)
            sizes.append(2 + 135 * (min(max(pct, 0), 75) / 75.0) ** 1.35)
            vals.append(max(min(val, 2.5), -1.5))

    cmap = LinearSegmentedColormap.from_list("dotplot_blue", ["#CFCFCF", "#B79CE8", "#7A55E8", "#0000FF"])
    sc = ax.scatter(xs, ys, s=sizes, c=vals, cmap=cmap, vmin=-1.5, vmax=2.5, edgecolors="none")
    ax.set_xticks(range(len(DOT_GENES)))
    ax.set_xticklabels(DOT_GENES, rotation=58, ha="right", fontsize=9, fontstyle="italic")
    ax.set_yticks(range(len(DOT_ROWS)))
    ax.set_yticklabels([display_label for display_label, _ in DOT_ROWS], fontsize=11)
    ax.invert_yaxis()
    ax.set_xlabel("Marker genes", fontsize=12)
    ax.set_ylabel("Cell type", fontsize=12)
    ax.set_title("Representative marker-gene dot plot", loc="left", fontsize=14, fontweight="bold", pad=7)
    ax.set_xlim(-0.8, len(DOT_GENES) - 0.2)
    ax.set_ylim(len(DOT_ROWS) - 0.35, -0.65)
    ax.grid(False)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    cbar = plt.colorbar(sc, cax=cax)
    cbar.ax.set_title("Scaled average\nexpression", fontsize=10, pad=8)
    cbar.ax.tick_params(labelsize=9, width=0.8, length=3)

    legend_ax.axis("off")
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.text(0.0, 0.98, "Percent Expressed", fontsize=11, ha="left", va="top")
    legend_sizes = [2 + 135 * (value / 75.0) ** 1.35 for value in (0, 25, 50, 75)]
    for y, label, size in zip([0.73, 0.52, 0.31, 0.10], ["0", "25", "50", "75"], legend_sizes):
        legend_ax.scatter([0.18], [y], s=size, color="black", linewidths=0)
        legend_ax.text(0.35, y, label, fontsize=10, va="center", ha="left")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["pdf.fonttype"] = 42

    fig = plt.figure(figsize=(11.8, 10.2), dpi=600, facecolor="white")

    ax_a = fig.add_axes([0.070, 0.595, 0.360, 0.325])
    ax_b = fig.add_axes([0.510, 0.595, 0.360, 0.325])
    cluster_crop = (222, 141, 1966, 1605)
    broad_crop = (222, 141, 2140, 1605)
    legend = [
        ("#F8766D", "Epithelial"),
        ("#C49A00", "Leukocytes"),
        ("#53B400", "Endothelial"),
        ("#00C1A2", "Stromal"),
        ("#00A9D6", "Smooth muscle"),
        ("#9A7DFF", "Neural"),
        ("#F564C9", "Other"),
    ]

    draw_clean_umap(ax_a, CLUSTER_UMAP, cluster_crop, "UMAP by transcriptional cluster")
    draw_clean_umap(
        ax_b,
        BROAD_UMAP,
        broad_crop,
        "UMAP by broad annotated category",
        labels=[("Leukocytes", 915, 965, 3.6, 1.15), ("Smooth muscle", 1000, 1528, 6.4, 2.55)],
        legend=legend,
    )
    ax_a.text(-0.18, 1.08, "a", transform=ax_a.transAxes, fontsize=15, fontweight="bold", va="top")
    ax_b.text(-0.18, 1.08, "b", transform=ax_b.transAxes, fontsize=15, fontweight="bold", va="top")

    ax_c = fig.add_axes([0.160, 0.085, 0.650, 0.405])
    size_ax = fig.add_axes([0.852, 0.300, 0.120, 0.150])
    cax = fig.add_axes([0.895, 0.130, 0.018, 0.135])
    draw_dotplot(ax_c, cax, size_ax)
    ax_c.text(-0.060, 1.080, "c", transform=ax_c.transAxes, fontsize=15, fontweight="bold", va="top")

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
