from __future__ import annotations

import argparse
import gzip
import hashlib
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np
import pandas as pd
import seaborn as sns


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = REPO_ROOT.parent / "zenodo_processed_data_record" if (REPO_ROOT.parent / "zenodo_processed_data_record").is_dir() else REPO_ROOT
SOURCE_DIR = DATA_ROOT / "source_data" / "figures"
FIGURE_DIR = REPO_ROOT / "outputs" / "figures" / "supplementary"
SAMPLE_ORDER = ["CTL1", "CTL2", "Llac1", "Llac2", "Slim1", "Slim2", "COM1", "COM2"]
SAMPLE_PLOT_ORDER = ["CTL_1", "CTL_2", "Llac_1", "Llac_2", "Slim_1", "Slim_2", "COM_1", "COM_2"]
GROUP_ORDER = ["CTL", "Llac", "Slim", "COM"]
GROUP_COLORS = {
    "CTL": "#4477AA",
    "Llac": "#228833",
    "Slim": "#CCAA33",
    "COM": "#CC6677",
}
SAMPLE_COLORS = {
    "CTL_1": "#315F8C",
    "CTL_2": "#7EA6CF",
    "Llac_1": "#187A67",
    "Llac_2": "#63B8A5",
    "Slim_1": "#A67C16",
    "Slim_2": "#E1BD68",
    "COM_1": "#A84A5B",
    "COM_2": "#E39AA7",
}
CELL_TYPE_ORDER = [
    "Enterocytes",
    "LREs",
    "Goblet cells",
    "Tuft-like cells",
    "Best4+ cells",
    "Enteroendocrine cells",
    "Neuronal cells",
    "Leukocytes",
    "Fibroblasts",
    "Endothelial cells",
    "Smooth muscle cells",
    "Acinar-like cells",
]
CELL_TYPE_DISPLAY_MAP = {
    "Acinar-like": "Acinar-like cells",
    "Best4+ cell": "Best4+ cells",
    "Tuft cells": "Tuft-like cells",
    "Smooth muscle": "Smooth muscle cells",
}
CELL_TYPE_COLORS = {
    "Enterocytes": "#4E79A7",
    "LREs": "#A0CBE8",
    "Goblet cells": "#F28E2B",
    "Tuft-like cells": "#FFBE7D",
    "Best4+ cells": "#59A14F",
    "Enteroendocrine cells": "#8CD17D",
    "Neuronal cells": "#8C6D5A",
    "Leukocytes": "#B6992D",
    "Fibroblasts": "#B07AA1",
    "Endothelial cells": "#499894",
    "Smooth muscle cells": "#79706E",
    "Acinar-like cells": "#E15759",
}


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 5.5,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "svg.hashsalt": "grouper-submission-20260714",
    }
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def panel_label(ax: plt.Axes, text: str) -> None:
    ax.text(
        -0.12,
        1.05,
        text,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def save_figure(fig: plt.Figure, out_dir: Path, stem: str) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    svg = out_dir / f"{stem}.svg"
    pdf = out_dir / f"{stem}.pdf"
    png = out_dir / f"{stem}.png"
    fig.savefig(svg, dpi=600, bbox_inches="tight", metadata={"Date": None})
    fig.savefig(pdf, dpi=600, bbox_inches="tight", metadata={"CreationDate": None})
    fig.savefig(png, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return {path.suffix.lstrip("."): sha256(path) for path in (svg, pdf, png)}


def draw_s1(qc: pd.DataFrame, out_dir: Path) -> dict[str, str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2), constrained_layout=True)
    palette = [GROUP_COLORS[group] for group in ["CTL", "CTL", "Llac", "Llac", "Slim", "Slim", "COM", "COM"]]
    panels = [
        ("nFeature_RNA", "Detected genes per nucleus", "Detected genes"),
        ("nCount_RNA", "UMI counts per nucleus", "UMIs"),
        ("percent_mt", "Mitochondrial fraction", "Mitochondrial UMIs (%)"),
    ]
    for ax, (column, title, ylabel), label in zip(axes.flat[:3], panels, "abc"):
        sns.violinplot(
            data=qc,
            x="sample_plot",
            y=column,
            hue="sample_plot",
            order=SAMPLE_PLOT_ORDER,
            hue_order=SAMPLE_PLOT_ORDER,
            palette=palette,
            cut=0,
            inner=None,
            linewidth=0.45,
            saturation=0.82,
            legend=False,
            ax=ax,
        )
        sns.boxplot(
            data=qc,
            x="sample_plot",
            y=column,
            order=SAMPLE_PLOT_ORDER,
            width=0.14,
            showfliers=False,
            color="white",
            boxprops={"edgecolor": "#222222", "linewidth": 0.55, "facecolor": "white"},
            whiskerprops={"color": "#222222", "linewidth": 0.55},
            capprops={"color": "#222222", "linewidth": 0.55},
            medianprops={"color": "#111111", "linewidth": 0.8},
            ax=ax,
        )
        ax.set_title(title, fontweight="bold", pad=5)
        ax.set_xlabel("")
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x", rotation=45)
        panel_label(ax, label)
        ax.margins(x=0.02)

    ax = axes.flat[3]
    counts = qc.groupby(["sample_plot", "group"], sort=False).size().reset_index(name="retained_nuclei")
    counts["sample_plot"] = pd.Categorical(counts["sample_plot"], SAMPLE_PLOT_ORDER, ordered=True)
    counts = counts.sort_values("sample_plot")
    ax.bar(
        counts["sample_plot"].astype(str),
        counts["retained_nuclei"],
        color=[GROUP_COLORS[group] for group in counts["group"]],
        edgecolor="white",
        linewidth=0.5,
    )
    ax.set_title("Retained nuclei after QC", fontweight="bold", pad=5)
    ax.set_xlabel("")
    ax.set_ylabel("Retained nuclei")
    ax.tick_params(axis="x", rotation=45)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{int(value):,}"))
    ax.set_ylim(0, max(counts["retained_nuclei"]) * 1.12)
    panel_label(ax, "d")

    return save_figure(fig, out_dir, "Supplementary_Figure_S1_QC")


def scatter_umap(ax: plt.Axes, data: pd.DataFrame, column: str, order: list[str], colors: dict[str, str]) -> None:
    shuffled = data.sample(frac=1, random_state=20260714)
    for value in order:
        subset = shuffled.loc[shuffled[column] == value]
        ax.scatter(
            subset["UMAP_1"],
            subset["UMAP_2"],
            s=0.11,
            c=colors[value],
            linewidths=0,
            alpha=0.72,
            rasterized=True,
        )
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.xaxis.set_major_locator(MaxNLocator(3))
    ax.yaxis.set_major_locator(MaxNLocator(3))
    ax.set_aspect("equal", adjustable="datalim")


def legend_handles(order: list[str], colors: dict[str, str]) -> list[Line2D]:
    return [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=colors[value], markeredgewidth=0, markersize=3.5, label=value)
        for value in order
    ]


def draw_s2(umap: pd.DataFrame, out_dir: Path) -> dict[str, str]:
    umap = umap.copy()
    umap["cell_type_display"] = umap["global_cell_type"].replace(CELL_TYPE_DISPLAY_MAP)
    observed_types = set(umap["cell_type_display"].unique())
    if observed_types != set(CELL_TYPE_ORDER):
        raise RuntimeError(f"Unexpected global cell types: {sorted(observed_types)}")

    fig = plt.figure(figsize=(7.2, 3.35), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=[5.2, 1.0])
    axes = [fig.add_subplot(grid[0, index]) for index in range(3)]
    legend_axes = [fig.add_subplot(grid[1, index]) for index in range(3)]

    scatter_umap(axes[0], umap, "sample_plot", SAMPLE_PLOT_ORDER, SAMPLE_COLORS)
    scatter_umap(axes[1], umap, "group", GROUP_ORDER, GROUP_COLORS)
    scatter_umap(axes[2], umap, "cell_type_display", CELL_TYPE_ORDER, CELL_TYPE_COLORS)

    titles = ["UMAP by library", "UMAP by treatment", "UMAP by global cell type"]
    for index, (ax, title) in enumerate(zip(axes, titles)):
        ax.set_title(title, fontweight="bold", pad=4)
        panel_label(ax, "abc"[index])

    legend_specs = [
        (SAMPLE_PLOT_ORDER, SAMPLE_COLORS, 4),
        (GROUP_ORDER, GROUP_COLORS, 4),
        (CELL_TYPE_ORDER, CELL_TYPE_COLORS, 3),
    ]
    for legend_ax, (order, colors, columns) in zip(legend_axes, legend_specs):
        legend_ax.axis("off")
        legend_ax.legend(
            handles=legend_handles(order, colors),
            loc="center",
            ncol=columns,
            columnspacing=0.8,
            handletextpad=0.3,
            borderaxespad=0,
        )

    return save_figure(fig, out_dir, "Supplementary_Figure_S2_UMAP")


def draw_s3(composition: pd.DataFrame, correlation: pd.DataFrame, out_dir: Path) -> dict[str, str]:
    comp = composition.copy()
    comp["cell_type_display"] = comp["cell_type"].replace(CELL_TYPE_DISPLAY_MAP)
    observed_types = set(comp["cell_type_display"].unique())
    if observed_types != set(CELL_TYPE_ORDER):
        raise RuntimeError(f"Unexpected composition cell types: {sorted(observed_types)}")
    heat = (
        comp.pivot(index="cell_type_display", columns="sample_plot", values="prop")
        .reindex(index=CELL_TYPE_ORDER, columns=SAMPLE_PLOT_ORDER)
        * 100
    )
    if heat.isna().any().any():
        raise RuntimeError("Composition heatmap contains missing values")

    corr = correlation.set_index("sample_plot").reindex(index=SAMPLE_PLOT_ORDER, columns=SAMPLE_PLOT_ORDER)
    if corr.isna().any().any() or not np.allclose(corr.values, corr.values.T, atol=1e-12):
        raise RuntimeError("Pseudobulk correlation matrix is incomplete or asymmetric")

    fig = plt.figure(figsize=(7.2, 4.3))
    grid = fig.add_gridspec(
        1,
        5,
        width_ratios=[1.55, 0.07, 0.35, 1.25, 0.07],
        left=0.19,
        right=0.98,
        bottom=0.18,
        top=0.88,
        wspace=0.12,
    )
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 3])]
    colorbar_axes = [fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[0, 4])]

    sns.heatmap(
        heat,
        ax=axes[0],
        cmap=sns.light_palette("#245D74", as_cmap=True),
        vmin=0,
        vmax=float(np.ceil(heat.to_numpy().max() / 5) * 5),
        cbar_ax=colorbar_axes[0],
        linewidths=0.25,
        linecolor="white",
    )
    colorbar_axes[0].set_title("Nuclei\n(%)", fontsize=6, pad=4)
    axes[0].set_title("Cell-type representation by library", fontweight="bold", pad=5)
    axes[0].set_xlabel("")
    axes[0].set_ylabel("")
    axes[0].tick_params(axis="x", rotation=45)
    axes[0].tick_params(axis="y", rotation=0)
    panel_label(axes[0], "a")

    off_diagonal = corr.to_numpy()[~np.eye(len(corr), dtype=bool)]
    lower = max(0.0, np.floor(off_diagonal.min() * 1000) / 1000)
    sns.heatmap(
        corr,
        ax=axes[1],
        cmap="Blues",
        vmin=lower,
        vmax=1.0,
        annot=True,
        fmt=".3f",
        annot_kws={"fontsize": 5},
        cbar_ax=colorbar_axes[1],
        cbar_kws={"format": "%.3f"},
        linewidths=0.25,
        linecolor="white",
    )
    colorbar_axes[1].set_title("Pearson r", fontsize=6, pad=4)
    axes[1].set_title("Pseudobulk expression correlation", fontweight="bold", pad=5)
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")
    axes[1].tick_params(axis="x", rotation=45)
    axes[1].tick_params(axis="y", rotation=0)
    panel_label(axes[1], "b")

    return save_figure(fig, out_dir, "Supplementary_Figure_S3_Representation")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Scientific Data Supplementary Figures S1-S3.")
    parser.add_argument("--qc", type=Path, default=SOURCE_DIR / "SuppFigureS1_qc_per_cell.csv.gz")
    parser.add_argument(
        "--final-umap",
        type=Path,
        default=SOURCE_DIR / "Fig3_global_umap_coordinates_final_annotation.csv.gz",
    )
    parser.add_argument(
        "--composition",
        type=Path,
        default=SOURCE_DIR / "Fig5a_global_celltype_composition_leukocytes_merged.csv",
    )
    parser.add_argument(
        "--correlation",
        type=Path,
        default=SOURCE_DIR / "SuppFigureS3_pseudobulk_correlation.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=FIGURE_DIR)
    parser.add_argument(
        "--write-hash-table",
        action="store_true",
        help="Write supplementary_figure_hashes.csv for a standalone QA run.",
    )
    args = parser.parse_args()

    qc = pd.read_csv(args.qc, compression="gzip")
    umap = pd.read_csv(args.final_umap, compression="gzip")
    composition = pd.read_csv(args.composition)
    correlation = pd.read_csv(args.correlation)

    if len(qc) != 102_036 or len(umap) != 102_036:
        raise RuntimeError("QC and UMAP sources must each contain 102,036 nuclei")
    if set(qc["cell_id"]) != set(umap["cell_id"]):
        raise RuntimeError("QC and UMAP cell sets differ")

    results = {
        "S1": draw_s1(qc, args.output_dir),
        "S2": draw_s2(umap, args.output_dir),
        "S3": draw_s3(composition, correlation, args.output_dir),
    }
    rows = []
    for figure, formats in results.items():
        for file_format, digest in formats.items():
            rows.append({"figure": figure, "format": file_format, "sha256": digest})
    if args.write_hash_table:
        pd.DataFrame(rows).to_csv(args.output_dir / "supplementary_figure_hashes.csv", index=False)
    for figure, formats in results.items():
        print(figure, formats)


if __name__ == "__main__":
    main()
