# Humpback Grouper Intestine snRNA-seq Analysis Scripts

This repository contains R scripts for the Scientific Data manuscript:

**A single-nucleus transcriptomic atlas of the humpback grouper (*Cromileptes altivelis*) intestine across dietary supplementation conditions**

The scripts inspect the deposited Seurat objects, export public metadata tables, and generate the main R-based submission figures from the processed data archive.

## Data Archive

Large processed data files are archived on Zenodo:

```text
https://doi.org/10.5281/zenodo.20632791
```

Download and unpack `humpback_grouper_intestine_snRNAseq_processed_data.zip`. The archive includes:

- `integrated_object/grouper_intestine_snRNAseq_seurat.rds`
- `integrated_object/grouper_intestine_immune_reclustering_seurat.rds`
- `cell_metadata.csv`
- `umap_coordinates.csv`
- cluster annotation summaries, marker-gene tables, figure source data and final figure files

By default, place the unpacked data folder next to this code repository:

```text
project_directory/
  analysis_code/
  humpback_grouper_intestine_snRNAseq_processed_data/
```

If the data folder is elsewhere, set `GROUPER_DATA_DIR` before running the scripts:

```r
Sys.setenv(GROUPER_DATA_DIR = "D:/path/to/humpback_grouper_intestine_snRNAseq_processed_data")
```

## Quick Start

Install missing R packages if needed:

```bash
Rscript 00_setup/00_install_required_packages.R
```

Run the complete analysis workflow from the repository root:

```bash
Rscript 00_run_analysis_workflow.R
```

The workflow writes a timestamped `analysis_run_*/` directory containing session information, RDS inventory summaries, metadata exports, generated figures, logs and `RUN_SUMMARY.txt`.

## Repository Structure

- `00_run_analysis_workflow.R`: one-command analysis workflow.
- `00_setup/`: package setup helper.
- `01_object_inventory/`: deposited Seurat object inspection.
- `02_export_metadata/`: public metadata and annotation table export.
- `03_make_figures/`: Figure 2, Figure 3 and Figure 4 generation.
- `04_environment/`: R session information and package notes.

See `WORKFLOW.md` for detailed inputs, outputs and figure-generation scope.

## Software

The manuscript analysis used:

- Cell Ranger v5.0.0 for primary count generation.
- R v4.4.0 for downstream analysis and visualization.
- Seurat v5.3.0 and Harmony v0.1 for integration and downstream analysis.

This repository was validated on 2026-06-10 using R 4.6.0 with Seurat 5.5.0, Harmony 2.0.4, ggplot2 4.0.3, patchwork 1.3.2 and scales 1.4.0. The workflow recovered 102,036 nuclei, 27 global clusters and 1,347 global-cluster-14 acinar-like nuclei.
