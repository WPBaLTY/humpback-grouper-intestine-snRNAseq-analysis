# Figure Data Provenance

This document records the data source and regeneration route for Scientific Data Figures 2-5. Figure 1 is submitted workflow artwork and is not regenerated here.

## Raw Starting Point

The analysis starts from eight per-library 10x Genomics matrices deposited with GEO `GSE326285`:

- `COM1`, `COM2`
- `CTL1`, `CTL2`
- `Llac1`, `Llac2`
- `Slim1`, `Slim2`

`metadata/tenx_sample_sheet.csv` records the original merge order, public names, delivery folder names such as `COM_1`, and legacy `CON1/CON2` aliases.

Each sample directory can be the matrix directory itself or contain `filtered_feature_bc_matrix/`, `outs/filtered_feature_bc_matrix/`, `raw_feature_bc_matrix/`, or `outs/raw_feature_bc_matrix/`.

This is the validated executable boundary. The package records hashes and parsed settings for all eight Cell Ranger reports, source fingerprints for the recovered project FASTA/GFF3 and annotation bundle, and the exact ordered 26,134-row Cell Ranger feature list. All 26,121 nuclear GFF gene IDs plus 13 `NC_021614.1` mitochondrial protein-coding IDs close exactly to that list. The original combined `cal_MT.fa` and exact 2 kb 3' UTR-extended `cal.3utr_MT.gtf` are still unavailable, so FASTQ-to-matrix realignment is documented but not claimed as byte-exact. See `REFERENCE_PROVENANCE_AUDIT.md`.

## Global Reconstruction

Install the validated core R packages and run:

```bash
Rscript environment/install_exact_r_packages.R
R_LIBS_USER=.r-library \
GROUPER_10X_ROOT=/path/to/10x_matrices \
Rscript scripts/R_analysis/00_reconstruct_global_atlas_from_10x.R
```

Validated versions are Seurat 5.3.0, SeuratObject 5.2.0, and harmony 1.2.3. The script uses `orig.ident` for Harmony, matching the original analysis.

The full rerun reproduced:

- 102,036 filtered nuclei
- 26,134 features in every delivered Cell Ranger matrix and 23,740 retained in the merged Seurat object after per-library `min.cells = 3`
- all 27 resolution-0.4 cluster assignments with 0 cell-level mismatches
- the 87 Figure 3 selected marker rows exactly
- all Figure 5a sample-by-cell-type counts exactly

Important generated outputs include:

- `outputs/global_atlas_reconstruction/GLOBAL_ATLAS_harmony_res0.4.rds`
- `outputs/global_atlas_reconstruction/GLOBAL_LEUKOCYTES_subset_for_reclustering.rds`
- `outputs/global_atlas_reconstruction/markers_res0.4_all.csv`
- `outputs/global_atlas_reconstruction/composition_celltype_bySample.csv`
- `source_data/figures/Fig3_global_markers_res0.4_all.csv`
- `source_data/figures/Fig3_global_marker_dotplot_source.csv`
- `source_data/figures/Fig5a_global_celltype_composition_leukocytes_merged.csv`

Large generated RDS files are intentionally excluded from GitHub.

The generated `GLOBAL_LEUKOCYTES_subset_for_reclustering.rds` is the default input to the immune branch below; no manually prepared composition table is used between the global and immune workflows.

## Final Global Processed Source

The historical final annotated object is exported by:

```bash
GROUPER_FINAL_GLOBAL_RDS=/path/to/seurat_FINAL_celltype_fixedMeta.rds \
Rscript scripts/R_analysis/00b_export_final_global_source_data.R
```

This writes the final Figure 3 UMAP source images, Figure 5a composition, and the compressed per-cell table:

`source_data/figures/Fig3_global_umap_coordinates_final_annotation.csv.gz`

The table contains cell ID, sample, treatment, replicate, final cluster, original cell type, public global cell type, and final UMAP coordinates for all nuclei. It provides the exact Figure 3 coordinate-level processed source even when UMAP layout differs across uwot/R versions.

Raw counts can be checked against the final object with:

```bash
GROUPER_10X_ROOT=/path/to/10x_matrices \
GROUPER_FINAL_GLOBAL_RDS=/path/to/seurat_FINAL_celltype_fixedMeta.rds \
Rscript scripts/R_analysis/00c_validate_final_object_against_10x.R
```

The completed audit confirmed exact sparse RNA count matrices.

## Figure 2

Figure 2 reads sequencing and Cell Ranger QC summaries directly:

- `metadata/samples.sequence.stat.xls`
- `metadata/samples.align.stat.xls`

These are library-level processing outputs, not manually entered plotting values.

```bash
python -B scripts/python_figure_assembly/rebuild_scientific_data_figure2_qc_final.py
```

## Figure 3

Figure 3 uses the final global annotation layer:

- panels a/b: final global UMAP source images, with final per-cell coordinates included in the CSV.GZ above
- panel c: `Fig3_global_marker_dotplot_source.csv`, containing the 480 direct cell-type-by-marker values generated from the normalized raw-derived object. Dot size uses `percent_expressed`; color uses `average_expression_scaled`. `Fig3_global_markers_res0.4_all.csv` remains included as upstream differential-marker evidence and is not relabelled as average expression.

The global immune-lineage category is `Leukocytes`. Detailed immune subtypes are not introduced at this layer.

```bash
python -B scripts/python_figure_assembly/rebuild_scientific_data_figure3_global_annotation.py
```

## Figure 4

Figure 4 starts from the `GLOBAL_LEUKOCYTES_subset_for_reclustering.rds` object generated directly by the eight-matrix global reconstruction. The immune script recreates the original sample-Harmony/resolution-0.2 immune clustering before applying the documented cleaning, reclustering, annotation, and B-cell split:

```bash
Rscript scripts/R_analysis/01_reconstruct_immune_final_object.R
Rscript scripts/R_analysis/02c_export_scientific_data_figure4_source_data.R
```

With Seurat 5.3.0, SeuratObject 5.2.0, and harmony 1.2.3, the raw-derived immune branch matched the historical intermediate and final branch exactly:

- 10,440 initial immune-lineage nuclei
- identical initial cell IDs/order, cluster assignments, sparse RNA count slots, and UMAP coordinates
- 5,061 final immune nuclei
- identical final subtype assignments and UMAP coordinates
- final subtype totals of 1,110 T cell (CCR7+), 586 Activated T (RORA+), 723 Activated lymphoid (CCL20hi), 1,027 NK-like cytotoxic, 108 B cells, 298 cDC1-like (XCR1+), 640 Monocytes/macrophages, 219 MoDC-like (CD209d+), 242 Granulocyte-like, and 108 Cycling (G2/M)

Included lightweight inputs are:

- `Fig4_immune_umap_coordinates.csv`
- `Fig4_immune_marker_dotplot_source.csv`
- `Fig4_immune_marker_features_used.csv`

The final dendritic-cell-like label is `cDC1-like (XCR1+)`.

```bash
python -B scripts/python_figure_assembly/rebuild_scientific_data_figure4_immune_reclustering.py
```

## Figure 5

Figure 5 uses:

- panel a: `Fig5a_global_celltype_composition_leukocytes_merged.csv`
- panel b: `Fig5b_immune_subtype_composition.csv`

Panel a comes from the final global annotation and keeps `Leukocytes` merged. Panel b comes from leukocyte reclustering and displays immune subtypes. Figure 5a counts were independently reproduced from raw 10x matrices with the validated package versions.

All 77 nonzero Figure 5b sample-by-subtype counts are exact direct counts from `Fig4_immune_umap_coordinates.csv`; proportions sum to one within each sample. Figure 5b is therefore derived from the same per-cell immune annotations used by Figure 4, not from an untracked plotting spreadsheet.

```bash
python -B scripts/python_figure_assembly/rebuild_scientific_data_figure5_composition.py
```

## Supplementary Figures S1-S3

Supplementary Figure S1 uses `SuppFigureS1_qc_per_cell.csv.gz`, which contains detected genes, UMI counts and mitochondrial fraction for every retained nucleus. `scripts/prepare_supplement_source_data.py` converts the historical `CON1/CON2` identifiers to public `COM1/COM2` identifiers and performs a one-to-one join against `Fig3_global_umap_coordinates_final_annotation.csv.gz`. The join passed for all 102,036 public cell IDs and reproduced the exact eight library counts.

Supplementary Figure S2 reads the final Figure 3 per-cell UMAP table directly. The library, treatment and global-cell-type panels therefore use the same coordinates and final annotation set; display-name expansion such as `Tuft cells` to `Tuft-like cells` does not alter the source annotation table.

Supplementary Figure S3a reads the checked Figure 5a global composition table. Panel b is generated from three deposited files:

- `SuppFigureS3_pseudobulk_logCPM.csv`: 23,738 expressed features by eight capture libraries.
- `SuppFigureS3_pseudobulk_correlation.csv`: eight-by-eight Pearson correlation matrix.
- `SuppFigureS3_pseudobulk_parameters.csv`: aggregation, normalization and transformation definitions.

The R export sums raw RNA counts by capture library, normalizes by total pseudobulk library size, applies `log2(CPM + 1)`, and computes Pearson correlation across the retained expressed features. The release validator recomputes every coefficient from the deposited logCPM matrix and requires a maximum absolute difference no greater than `1e-12`. This is library-level validation, not nucleus-level inference.

```bash
python -B scripts/python_figure_assembly/rebuild_scientific_data_supplementary_figures.py
```

The script writes deterministic SVG, PDF and 600 dpi PNG files. Consecutive clean runs produced identical hashes for all nine files; PDF creation dates are removed, SVG text remains editable, and a fixed SVG hashsalt prevents random internal IDs.

## GitHub Boundary

The lightweight repository excludes FASTQ files, 10x matrix directories, BAM files, and large RDS objects. It includes the scripts, exact package-version records, aggregate source tables, final per-cell UMAP/annotation metadata, and audit results needed to trace and redraw Figures 2-5 and Supplementary Figures S1-S3.
