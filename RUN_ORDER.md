# Run Order

This repository covers Scientific Data Figures 2-5, Supplementary Figures S1-S3, and Table 2/Table 3. Figure 1 is retained as submitted artwork and is not regenerated.

## 0. Reference Provenance Snapshot

`metadata/cellranger_run_metadata.csv`, `metadata/reference_source_manifest.csv`, `metadata/reference_feature_audit.csv` and `metadata/reference_features_26134.csv.gz` are validated by the release checker below. To regenerate them from the external source FASTA/GFF3, eight Cell Ranger reports and eight delivered `genes.tsv` files, run `scripts/audit_reference_provenance.py` with the inputs shown in `REFERENCE_PROVENANCE_AUDIT.md`.

The reference snapshot distinguishes 26,134 Cell Ranger reference features from the 23,740-feature Seurat union produced after per-library `CreateSeuratObject(min.cells = 3)` filtering. It does not claim byte-exact FASTQ realignment without the original combined FASTA and exact 2 kb 3' UTR-extended GTF.

## 1. Direct Figure 2-5 and Supplementary Figure S1-S3 Rebuild

From the repository root:

```bash
conda env create -f environment/figure_environment_win-64.yml
conda activate grouper-figure-win64
python -B scripts/python_figure_assembly/rebuild_scientific_data_figure2_qc_final.py
python -B scripts/python_figure_assembly/rebuild_scientific_data_figure3_global_annotation.py
python -B scripts/python_figure_assembly/rebuild_scientific_data_figure4_immune_reclustering.py
python -B scripts/python_figure_assembly/rebuild_scientific_data_figure5_composition.py
python -B scripts/python_figure_assembly/rebuild_scientific_data_supplementary_figures.py
```

These commands use included metadata and source data. Main PNG/PDF outputs are written under `outputs/figures/scientific_data/`; supplementary SVG, PDF, and 600 dpi PNG outputs are written under `outputs/figures/supplementary/`. TIFF output is disabled by default to keep the public tree lightweight; set `GROUPER_WRITE_TIFF=1` before running a main-figure script only when a TIFF is required. Add `--write-hash-table` to the supplementary command only for a standalone QA hash table.

Expected PNG canvas sizes:

- Figure 2: `3205 x 2388`
- Figure 3: `7080 x 6120`
- Figure 4: `6900 x 4890`
- Figure 5: `7301 x 5419`
- Supplementary Figure S1: `4389 x 3189`
- Supplementary Figure S2: `4389 x 2079`
- Supplementary Figure S3: `4238 x 2280`

The images are freshly assembled from source data. The win-64 conda file pins package build strings, including FreeType 2.13.3, and strict validation checks the Arial font fingerprints; under that render lock the PNG and PDF outputs are byte-identical to the deposited reference hashes. Main-figure PDFs use embedded CID TrueType fonts, contain no Type 3 fonts, and omit dynamic creation dates.
The assembly scripts do not modify included source-data files by default.

The deposited S1 per-cell source closes exactly to the final 102,036-cell UMAP table. To recreate that join from the historical QC inventory, run:

```bash
python -B scripts/prepare_supplement_source_data.py \
  --inventory /path/to/Seurat_Inventory_cell_metadata.csv \
  --final-umap source_data/figures/Fig3_global_umap_coordinates_final_annotation.csv.gz \
  --output-dir source_data/figures
```

To recreate the S3 pseudobulk source from the final global Seurat object, run this under the exact R environment described below:

```bash
Rscript scripts/R_analysis/03_export_supplementary_pseudobulk_source.R \
  /path/to/seurat_FINAL_celltype_fixedMeta.rds \
  source_data/figures/Fig3_global_umap_coordinates_final_annotation.csv.gz \
  source_data/figures
```

The deposited matrix contains 23,738 expressed features across eight library-level pseudobulks. The validator recomputes all 64 Pearson coefficients from that matrix and requires a maximum absolute difference of at most `1e-12` from the plotting table.

For cross-platform source-data checks or a visually equivalent redraw, `python -m pip install -r environment/python_requirements.txt` is also supported. Top-level pip versions alone do not lock the FreeType rasterizer, so that route is not a byte-identical PNG claim.

## 2. Exact Global Reconstruction From 10x Matrices

The exact raw reconstruction requires Seurat 5.3.0, SeuratObject 5.2.0, and harmony 1.2.3. Install them into a repository-local library:

```bash
Rscript environment/install_exact_r_packages.R
```

On Linux/macOS:

```bash
R_LIBS_USER=.r-library \
GROUPER_10X_ROOT=/path/to/10x_matrices \
Rscript scripts/R_analysis/00_reconstruct_global_atlas_from_10x.R
```

On Windows PowerShell:

```powershell
$env:R_LIBS_USER = "$PWD\.r-library"
$env:GROUPER_10X_ROOT = (Resolve-Path ".\external_inputs\10x_matrices").Path
Rscript scripts\R_analysis\00_reconstruct_global_atlas_from_10x.R
```

Use an ASCII-only repository and output path on Windows. Some R/Windows locale combinations cannot open or write paths containing non-ASCII characters.

The sample sheet accepts public folders (`COM1`, `CTL1`, and so on), delivery folders (`COM_1`, `CTL_1`, and so on), and legacy `CON1/CON2` aliases. `analysis_order` preserves the original merge order.

The script stops if core package versions differ. `GROUPER_ALLOW_VERSION_DRIFT=1` bypasses the guard only for an explicitly non-exact exploratory run.

Validated exact results:

- 102,036 filtered nuclei
- 27 clusters at resolution 0.4
- 0 cell-level cluster-assignment mismatches versus the final object
- exact Figure 5a counts
- an exact 10,440-cell global `Leukocytes` subset for the immune branch

See `audit/annotation/raw_10x_full_rerun_audit.md`.

## 3. Export Final Global Processed Source Data

To regenerate the deposited final per-cell UMAP/annotation table and Figure 5a from the historical final annotated object:

```bash
GROUPER_FINAL_GLOBAL_RDS=/path/to/seurat_FINAL_celltype_fixedMeta.rds \
Rscript scripts/R_analysis/00b_export_final_global_source_data.R
```

Key output:

`source_data/figures/Fig3_global_umap_coordinates_final_annotation.csv.gz`

This route preserves the exact final UMAP coordinates used for Figure 3. UMAP coordinates can drift across uwot/R versions even when counts and clusters are identical.

## 4. Validate Final Object Against Raw 10x Counts

```bash
GROUPER_10X_ROOT=/path/to/10x_matrices \
GROUPER_FINAL_GLOBAL_RDS=/path/to/seurat_FINAL_celltype_fixedMeta.rds \
Rscript scripts/R_analysis/00c_validate_final_object_against_10x.R
```

This reconstructs the filtered count layer and compares feature order, public cell IDs, and sparse matrix `Dim/i/p/x` slots. The completed local audit passed all checks.

## 5. Leukocyte-Reclustering Branch

After section 2, the default input is the generated global `Leukocytes` subset:

```bash
Rscript scripts/R_analysis/01_reconstruct_immune_final_object.R
Rscript scripts/R_analysis/02c_export_scientific_data_figure4_source_data.R
```

This branch first recreates the historical initial immune Harmony clustering from the global subset, then applies the documented contaminant-removal sequence, immune annotation, and B-cell split. It regenerates Figure 4 source data and Figure 5b composition.

Validated exact results:

- 10,440 initial immune-lineage nuclei with cell-level cluster assignments identical to the historical intermediate
- exact initial sparse RNA count matrix and UMAP coordinates
- 5,061 final immune nuclei
- exact final immune subtype assignments and UMAP coordinates
- all 77 nonzero Figure 5b sample-by-subtype counts equal direct counts from the Figure 4 per-cell table

For compatibility with the historical intermediate object:

```bash
GROUPER_IMMUNE_INPUT_RDS=/path/to/IMMUNEbaicell.rds \
GROUPER_IMMUNE_INPUT_STAGE=preclustered \
Rscript scripts/R_analysis/01_reconstruct_immune_final_object.R
```

The exact Seurat/SeuratObject/harmony version guard also applies to this branch.

## 6. Doublet Assessment

Run this after sections 2 and 5, because it uses the reconstructed global object and the saved 10,440-cell initial immune object. Install the validated doublet environment into a separate local library:

```bash
Rscript environment/install_doublet_r_packages.R
```

On Linux/macOS:

```bash
R_LIBS_USER=.r-library-doublet \
Rscript scripts/R_analysis/00d_assess_doublets_scDblFinder.R
```

On Windows PowerShell:

```powershell
$env:R_LIBS_USER = "$PWD\.r-library-doublet"
Rscript scripts\R_analysis\00d_assess_doublets_scDblFinder.R
```

The script analyzes all 102,036 nuclei separately by capture library in cluster-aware and random sensitivity modes. It writes one merged per-cell table plus summary tables under `audit/doublet_assessment/`. The seeds and `dbr.per1k` value are recorded in `doublet_analysis_parameters.csv`.

## 7. Audit the 10,440-to-5,061 Immune Filtering

Switch back to the exact Seurat reconstruction library and run:

```bash
R_LIBS_USER=.r-library \
Rscript scripts/R_analysis/01b_audit_immune_filtering_10440_to_5061.R
```

This reruns the four cluster exclusions, writes a cell-level fate table, marker and QC evidence, merges the scDblFinder calls, and requires the reconstructed 5,061-cell object to match the deposited Figure 4 cell set, annotations, and UMAP coordinates. See `DOUBLETS_AND_IMMUNE_FILTERING_AUDIT.md`.

## 8. Annotation Tables

```bash
Rscript scripts/R_analysis/02b_scientific_data_global_leukocyte_update.R
```

This writes the checked Table 2/Table 3 and Figure 3 marker-panel files. Figure 3c reads the raw-derived `percent_expressed` and `average_expression_scaled` summary directly. Global outputs use `Leukocytes`; immune-focused outputs use subtypes including `cDC1-like (XCR1+)`.

See `FIGURE_DATA_PROVENANCE.md` and `metadata/script_manifest.csv` for the full input/output map.

## 9. Automated Release Validation

Quick source-data and package validation uses only the Python standard library:

```bash
python -B scripts/validate_release_package.py
```

After creating and activating `grouper-figure-win64` and making `Rscript` available on `PATH`, run the full validation:

```bash
python -B scripts/validate_release_package.py --check-r --run-figures
```

Full mode parses retained R and Python scripts, validates Cell Ranger/reference provenance, the deposited doublet and filtering ledgers, closes S1 to the final cell set, recomputes S3 correlations, rebuilds Figures 2-5 and S1-S3, checks expected canvases and nonblank image content, verifies deterministic PNG/PDF hashes and main-figure PDF font/metadata properties, and confirms that figure assembly did not alter any included source-data file.
