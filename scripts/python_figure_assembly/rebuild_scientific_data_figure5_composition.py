from pathlib import Path
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "sans-serif"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_ROOT = REPO_ROOT.parent / "zenodo_processed_data_record" if (REPO_ROOT.parent / "zenodo_processed_data_record").is_dir() else REPO_ROOT


def env_path(name, default):
    value = os.environ.get(name)
    return Path(value) if value else Path(default)


GLOBAL_COMP = env_path(
    "GROUPER_GLOBAL_COMPOSITION_CSV",
    REPO_ROOT / "external_inputs" / "composition_celltype_bySample.csv",
)
IMMUNE_COMP = env_path(
    "GROUPER_IMMUNE_COMPOSITION_CSV",
    REPO_ROOT / "external_inputs" / "Fig3_immune_composition_long.csv",
)
SOURCE_OUT = env_path(
    "GROUPER_SOURCE_DATA_OUT",
    DATA_ROOT / "source_data" / "figures",
)
FIGURE_OUT = env_path(
    "GROUPER_FIGURE5_OUT_DIR",
    REPO_ROOT / "outputs" / "figures" / "scientific_data",
)

PNG_OUT = env_path("GROUPER_FIGURE5_PNG", FIGURE_OUT / "Figure5_composition_final.png")
PDF_OUT = env_path("GROUPER_FIGURE5_PDF", FIGURE_OUT / "Figure5_composition_final.pdf")
TIFF_OUT = env_path("GROUPER_FIGURE5_TIFF", FIGURE_OUT / "Figure5_composition_final.tiff")
WRITE_TIFF = os.environ.get("GROUPER_WRITE_TIFF", "0") == "1"
WRITE_NORMALIZED_SOURCE = os.environ.get("GROUPER_WRITE_NORMALIZED_SOURCE", "0") == "1"


MAJOR_PALETTE = {
    "Enterocytes": "#2C7FB8",
    "LREs": "#67A9CF",
    "Goblet cells": "#FF8C1A",
    "Tuft-like cells": "#F39C12",
    "Best4+ cells": "#33A02C",
    "Enteroendocrine cells": "#7BC043",
    "Leukocytes": "#9C7A00",
    "Fibroblasts": "#D95F5F",
    "Endothelial cells": "#9467BD",
    "Smooth muscle cells": "#9A80B9",
    "Neuronal cells": "#8C564B",
    "Acinar-like cells": "#D65AA5",
}
MAJOR_DISPLAY_MAP = {
    "Tuft cells": "Tuft-like cells",
    "Best4+ cell": "Best4+ cells",
    "Smooth muscle": "Smooth muscle cells",
    "Acinar-like": "Acinar-like cells",
}

IMMUNE_PALETTE = {
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

SAMPLE_ORDER_GLOBAL = ["CTL1", "CTL2", "Llac1", "Llac2", "Slim1", "Slim2", "COM1", "COM2"]
SAMPLE_ORDER_IMMUNE = ["CTL1", "CTL2", "Llac1", "Llac2", "Slim1", "Slim2", "COM1", "COM2"]
SAMPLE_LABELS = ["CTL_1", "CTL_2", "Llac_1", "Llac_2", "Slim_1", "Slim_2", "COM_1", "COM_2"]


def sample_plot_label(series):
    return (
        series.astype(str)
        .str.replace("^CON", "COM", regex=True)
        .str.replace(r"([A-Za-z]+)(\d+)$", r"\1_\2", regex=True)
    )


def ensure_sample_plot(df):
    df = df.copy()
    if "sample" in df.columns:
        df["sample"] = df["sample"].astype(str).str.replace("^CON", "COM", regex=True)
    if "group" in df.columns:
        df["group"] = df["group"].astype(str).str.replace("^CON$", "COM", regex=True)
    if "group_plot" in df.columns:
        df["group_plot"] = df["group_plot"].astype(str).str.replace("^CON$", "COM", regex=True)
    if "sample_plot" not in df.columns:
        df["sample_plot"] = sample_plot_label(df["sample"])
    else:
        df["sample_plot"] = df["sample_plot"].astype(str).str.replace("^CON", "COM", regex=True)
    return df


def build_global_source():
    fallback = SOURCE_OUT / "Fig5a_global_celltype_composition_leukocytes_merged.csv"
    if not GLOBAL_COMP.exists() and fallback.exists():
        return ensure_sample_plot(pd.read_csv(fallback))
    comp = pd.read_csv(GLOBAL_COMP)
    label_map = {
        "Enterocyte": "Enterocytes",
        "LRE": "LREs",
        "Goblet_cell": "Goblet cells",
        "Tuft_cell": "Tuft cells",
        "Best4_cell": "Best4+ cell",
        "EEC": "Enteroendocrine cells",
        "Macrophages_monocytes": "Leukocytes",
        "Dc_cell": "Leukocytes",
        "T_cell": "Leukocytes",
        "NK_cell": "Leukocytes",
        "Endothelial_cell": "Endothelial cells",
        "Fibroblasts": "Fibroblasts",
        "Smooth_muscle": "Smooth muscle",
        "Nerve_cell": "Neuronal cells",
        "Acinar_like": "Acinar-like",
    }
    comp["cell_type"] = comp["celltype"].map(label_map)
    if comp["cell_type"].isna().any():
        missing = sorted(comp.loc[comp["cell_type"].isna(), "celltype"].unique())
        raise ValueError(f"Unmapped global cell types: {missing}")
    agg = comp.groupby(["sample", "group", "replicate", "cell_type"], as_index=False)["n"].sum()
    totals = agg.groupby("sample")["n"].transform("sum")
    agg["sample_plot"] = sample_plot_label(agg["sample"])
    agg["group_plot"] = agg["group"].replace({"CON": "COM"})
    agg["prop"] = agg["n"] / totals
    return ensure_sample_plot(agg[["sample", "group", "replicate", "cell_type", "n", "sample_plot", "group_plot", "prop"]])


def build_immune_source():
    fallback = SOURCE_OUT / "Fig5b_immune_subtype_composition.csv"
    if not IMMUNE_COMP.exists() and fallback.exists():
        immune = ensure_sample_plot(pd.read_csv(fallback))
        immune["subtype"] = immune["subtype"].replace(
            {"B cell": "B cells", "Monocyte/Macrophage": "Monocytes/macrophages"}
        )
        return immune[["group", "sample", "subtype", "n", "prop", "sample_plot"]]
    fig4_umap = SOURCE_OUT / "Fig4_immune_umap_coordinates.csv"
    if not IMMUNE_COMP.exists() and fig4_umap.exists():
        coords = pd.read_csv(fig4_umap)
        immune = (
            coords.groupby(["group", "sample", "immune_subtype"], as_index=False)
            .size()
            .rename(columns={"immune_subtype": "subtype", "size": "n"})
        )
        immune["prop"] = immune["n"] / immune.groupby("sample")["n"].transform("sum")
        immune["sample_plot"] = sample_plot_label(immune["sample"])
        return immune[["group", "sample", "subtype", "n", "prop", "sample_plot"]]
    immune = pd.read_csv(IMMUNE_COMP)
    display_map = {
        "B cell": "B cells",
        "Monocyte/Macrophage": "Monocytes/macrophages",
    }
    immune["subtype"] = immune["subtype"].replace(display_map)
    immune["sample_plot"] = sample_plot_label(immune["sample"])
    return ensure_sample_plot(immune[["group", "sample", "subtype", "n", "prop", "sample_plot"]])


def draw_stacked(ax, df, value_col, category_col, category_order, sample_order, labels, palette, title, legend_title, sample_col="sample_plot"):
    x = np.arange(len(sample_order))
    bottom = np.zeros(len(sample_order))
    for category in category_order:
        vals = []
        for sample in sample_order:
            hit = df[(df[sample_col] == sample) & (df[category_col] == category)]
            vals.append(float(hit[value_col].iloc[0]) if not hit.empty else 0.0)
        ax.bar(
            x,
            vals,
            bottom=bottom,
            color=palette[category],
            edgecolor="white",
            linewidth=0.2,
            width=0.72,
            label=category,
        )
        bottom += np.array(vals)
    ax.set_ylim(0, 1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=9)
    ax.set_ylabel("Relative abundance", fontsize=10)
    ax.set_xlabel("Library", fontsize=10)
    ax.set_title(title, loc="left", fontsize=14, fontweight="bold")
    ax.legend(
        title=legend_title,
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=8.5 if legend_title == "Immune subtype" else 9,
        title_fontsize=9.5 if legend_title == "Immune subtype" else 10,
    )
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


def main():
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    PNG_OUT.parent.mkdir(parents=True, exist_ok=True)
    PDF_OUT.parent.mkdir(parents=True, exist_ok=True)
    if WRITE_TIFF:
        TIFF_OUT.parent.mkdir(parents=True, exist_ok=True)
    plt.rcParams["pdf.fonttype"] = 42

    global_source = build_global_source()
    immune_source = build_immune_source()

    if WRITE_NORMALIZED_SOURCE:
        global_source.to_csv(SOURCE_OUT / "Fig5a_global_celltype_composition_leukocytes_merged.csv", index=False)
        immune_source.to_csv(SOURCE_OUT / "Fig5b_immune_subtype_composition.csv", index=False)

    global_plot = global_source.copy()
    global_plot["cell_type_display"] = global_plot["cell_type"].replace(MAJOR_DISPLAY_MAP)
    observed_global_types = set(global_plot["cell_type_display"])
    if observed_global_types != set(MAJOR_PALETTE):
        raise RuntimeError(f"Unexpected global display labels: {sorted(observed_global_types)}")

    fig = plt.figure(figsize=(11.8, 10.2), dpi=600)
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.42)
    ax1 = fig.add_subplot(gs[0])
    draw_stacked(
        ax1,
        global_plot,
        "prop",
        "cell_type_display",
        list(MAJOR_PALETTE.keys()),
        SAMPLE_LABELS,
        SAMPLE_LABELS,
        MAJOR_PALETTE,
        "Major intestinal cell-type composition per library",
        "Cell type",
    )
    ax1.text(-0.065, 1.05, "a", transform=ax1.transAxes, fontsize=15, fontweight="bold")

    ax2 = fig.add_subplot(gs[1])
    draw_stacked(
        ax2,
        immune_source,
        "prop",
        "subtype",
        list(IMMUNE_PALETTE.keys()),
        SAMPLE_LABELS,
        SAMPLE_LABELS,
        IMMUNE_PALETTE,
        "Immune-subtype composition per library",
        "Immune subtype",
    )
    ax2.text(-0.065, 1.05, "b", transform=ax2.transAxes, fontsize=15, fontweight="bold")

    fig.savefig(PNG_OUT, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(
        PDF_OUT,
        bbox_inches="tight",
        facecolor="white",
        metadata={"CreationDate": None},
    )
    if WRITE_TIFF:
        fig.savefig(TIFF_OUT, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
