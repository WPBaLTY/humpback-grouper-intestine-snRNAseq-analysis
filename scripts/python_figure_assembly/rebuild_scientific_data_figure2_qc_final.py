from pathlib import Path
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from PIL import Image


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_ROOT = REPO_ROOT.parent / "zenodo_processed_data_record" if (REPO_ROOT.parent / "zenodo_processed_data_record").is_dir() else REPO_ROOT
METADATA_DIR = DATA_ROOT / "metadata"
OUT_DIR = REPO_ROOT / "outputs" / "figures" / "scientific_data"

SEQ_FILE = METADATA_DIR / "samples.sequence.stat.xls"
ALIGN_FILE = METADATA_DIR / "samples.align.stat.xls"

PNG_OUT = OUT_DIR / "Figure2_QC_final.png"
PDF_OUT = OUT_DIR / "Figure2_QC_final.pdf"
TIFF_OUT = OUT_DIR / "Figure2_QC_final.tiff"
WRITE_TIFF = os.environ.get("GROUPER_WRITE_TIFF", "0") == "1"


GROUP_COLORS = {
    "CTL": "#4C78A8",
    "Llac": "#54A24B",
    "Slim": "#E88900",
    "COM": "#C44E52",
}


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def numeric(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.replace(",", "", regex=False)
        .str.replace("%", "", regex=False)
        .astype(float)
    )


def sample_groups(samples: pd.Series) -> pd.Series:
    return samples.astype(str).str.split("_").str[0]


def add_panel_label(ax, label):
    ax.text(-0.17, 1.08, label, transform=ax.transAxes, fontsize=16, fontweight="bold", va="top")


def polish_axes(ax):
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.6)
    ax.spines["bottom"].set_linewidth(1.6)
    ax.tick_params(axis="both", labelsize=11, width=1.4, length=5, colors="#333333")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="white", font="Arial")
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["pdf.fonttype"] = 42

    seq = read_table(SEQ_FILE)
    aln = read_table(ALIGN_FILE)
    df = seq.merge(aln, on="Sample", how="inner")
    df["reads_million"] = numeric(df["Number of Reads"]) / 1e6
    df["valid_barcodes"] = numeric(df["Valid Barcodes"])
    df["saturation"] = numeric(df["Sequencing Saturation"])
    df["q30_rna"] = numeric(df["Q30 Bases in RNA Read"])
    df["q30_umi"] = numeric(df["Q30 Bases in UMI"])
    df["genome_map"] = numeric(df["Reads Mapped Confidently to Genome"])
    df["transcriptome_map"] = numeric(df["Reads Mapped Confidently to Transcriptome"])
    df["estimated_nuclei"] = numeric(df["Estimated Number of Cells"])
    df["fraction_reads_in_cells"] = numeric(df["Fraction Reads in Cells"])
    df["median_genes"] = numeric(df["Median Genes per Cell"])
    df["median_umis"] = numeric(df["Median UMI Counts per Cell"])

    samples = df["Sample"]
    x = np.arange(len(samples))
    palette = [GROUP_COLORS[g] for g in sample_groups(samples)]

    fig = plt.figure(figsize=(3205 / 300, 2388 / 300), dpi=300, facecolor="white")
    gs = fig.add_gridspec(2, 2, hspace=0.72, wspace=0.60, left=0.105, right=0.915, top=0.920, bottom=0.135)

    ax1 = fig.add_subplot(gs[0, 0])
    ax1.bar(x, df["reads_million"], color=palette, edgecolor="white", linewidth=1.0)
    ax1.set_title("Sequencing depth by library", loc="left", fontsize=15, fontweight="bold", pad=8)
    ax1.set_ylabel("Reads (millions)", fontsize=13)
    ax1.set_xticks(x)
    ax1.set_xticklabels(samples, rotation=50, ha="right", fontsize=10)
    ax1.set_ylim(0, 535)
    polish_axes(ax1)
    add_panel_label(ax1, "a")

    ax2 = fig.add_subplot(gs[0, 1])
    heatmap_df = (
        df[
            [
                "valid_barcodes",
                "saturation",
                "q30_rna",
                "q30_umi",
                "genome_map",
                "transcriptome_map",
            ]
        ]
        .rename(
            columns={
                "valid_barcodes": "Valid barcodes",
                "saturation": "Saturation",
                "q30_rna": "Q30 RNA",
                "q30_umi": "Q30 UMI",
                "genome_map": "Genome map",
                "transcriptome_map": "Transcriptome map",
            }
        )
        .T
    )
    sns.heatmap(
        heatmap_df,
        ax=ax2,
        cmap=sns.color_palette("YlGnBu", as_cmap=True),
        vmin=55,
        vmax=100,
        linewidths=0.8,
        linecolor="white",
        cbar_kws={"label": "%", "ticks": [55, 60, 65, 70, 75, 80, 85, 90, 95, 100]},
    )
    ax2.set_title("Library-level technical metrics", loc="left", fontsize=15, fontweight="bold", pad=8)
    ax2.set_xlabel("")
    ax2.set_ylabel("")
    ax2.set_xticklabels(samples, rotation=50, ha="right", fontsize=10)
    ax2.tick_params(axis="y", labelsize=11, rotation=0)
    ax2.collections[0].colorbar.ax.tick_params(labelsize=10, width=1.2)
    ax2.collections[0].colorbar.set_label("%", fontsize=12)
    add_panel_label(ax2, "b")

    ax3 = fig.add_subplot(gs[1, 0])
    ax3.bar(x, df["estimated_nuclei"] / 1000, color=palette, edgecolor="white", linewidth=1.0)
    ax3.set_title("Nucleus yield and read recovery", loc="left", fontsize=15, fontweight="bold", pad=8)
    ax3.set_ylabel("Estimated nuclei (thousands)", fontsize=13)
    ax3.set_xticks(x)
    ax3.set_xticklabels(samples, rotation=50, ha="right", fontsize=10)
    ax3.set_ylim(0, 14.8)
    polish_axes(ax3)
    ax3b = ax3.twinx()
    ax3b.plot(x, df["fraction_reads_in_cells"], color="#222222", marker="o", lw=2.0, markersize=6)
    ax3b.set_ylabel("Reads in cells (%)", fontsize=13)
    ax3b.set_ylim(62, 73.2)
    ax3b.grid(False)
    ax3b.spines["top"].set_visible(False)
    ax3b.spines["right"].set_linewidth(1.6)
    ax3b.tick_params(axis="y", labelsize=11, width=1.4, length=5, colors="#333333")
    add_panel_label(ax3, "c")

    ax4 = fig.add_subplot(gs[1, 1])
    width = 0.34
    ax4.bar(x - width / 2, df["median_genes"], width=width, color="#8064A2", edgecolor="white", linewidth=0.8, label="Median genes")
    ax4b = ax4.twinx()
    ax4b.bar(x + width / 2, df["median_umis"], width=width, color="#EF5A83", edgecolor="white", linewidth=0.8, label="Median UMIs")
    ax4.set_title("Molecular complexity per nucleus", loc="left", fontsize=15, fontweight="bold", pad=8)
    ax4.set_xticks(x)
    ax4.set_xticklabels(samples, rotation=50, ha="right", fontsize=10)
    ax4.set_ylabel("Median genes per nucleus", fontsize=13)
    ax4b.set_ylabel("Median UMI counts per nucleus", fontsize=13, labelpad=7)
    ax4.set_ylim(0, 1450)
    ax4b.set_ylim(0, 2700)
    polish_axes(ax4)
    ax4b.grid(False)
    ax4b.spines["top"].set_visible(False)
    ax4b.spines["right"].set_linewidth(1.6)
    ax4b.tick_params(axis="y", labelsize=11, width=1.4, length=5, colors="#333333")
    handles_a, labels_a = ax4.get_legend_handles_labels()
    handles_b, labels_b = ax4b.get_legend_handles_labels()
    ax4.legend(handles_a + handles_b, labels_a + labels_b, loc="upper center", bbox_to_anchor=(0.50, 1.02), ncol=2, frameon=False, fontsize=10)
    add_panel_label(ax4, "d")

    fig.savefig(PNG_OUT, dpi=300, facecolor="white")
    fig.savefig(
        PDF_OUT,
        dpi=300,
        facecolor="white",
        metadata={"CreationDate": None},
    )
    plt.close(fig)

    if WRITE_TIFF:
        Image.open(PNG_OUT).convert("RGB").save(TIFF_OUT, dpi=(300, 300))
    print(PNG_OUT)


if __name__ == "__main__":
    main()
