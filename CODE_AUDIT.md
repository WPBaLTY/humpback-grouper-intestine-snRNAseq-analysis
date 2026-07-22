# Code Audit Notes

Date: 2026-07-14

## Summary

This lightweight package was audited for Scientific Data Figures 2-5, Supplementary Figures S1-S3, Table 2/Table 3, and the provenance chain from eight per-library 10x matrices to the final global and leukocyte-focused source data.

Figure 1 workflow artwork is out of scope and is retained as the submitted manuscript/Zenodo figure.

## Reproducibility Layers

Directly runnable from included metadata/source data:

- `scripts/python_figure_assembly/rebuild_scientific_data_figure2_qc_final.py`
- `scripts/python_figure_assembly/rebuild_scientific_data_figure3_global_annotation.py`
- `scripts/python_figure_assembly/rebuild_scientific_data_figure4_immune_reclustering.py`
- `scripts/python_figure_assembly/rebuild_scientific_data_figure5_composition.py`
- `scripts/python_figure_assembly/rebuild_scientific_data_supplementary_figures.py`

Raw 10x reconstruction requiring externally deposited matrices:

- `scripts/R_analysis/00_reconstruct_global_atlas_from_10x.R`

Final global-object export and raw-count validation requiring the historical final annotated RDS:

- `scripts/R_analysis/00b_export_final_global_source_data.R`
- `scripts/R_analysis/00c_validate_final_object_against_10x.R`

Leukocyte-reclustering branch runnable from the global reconstruction output:

- `scripts/R_analysis/01_reconstruct_immune_final_object.R`
- `scripts/R_analysis/02c_export_scientific_data_figure4_source_data.R`

Doublet and immune-filtering validation:

- `scripts/R_analysis/00d_assess_doublets_scDblFinder.R`
- `scripts/R_analysis/01b_audit_immune_filtering_10440_to_5061.R`

Source-data/table synchronization from included files:

- `scripts/R_analysis/02b_scientific_data_global_leukocyte_update.R`

Supplementary source exports with external historical inputs:

- `scripts/prepare_supplement_source_data.py`
- `scripts/R_analysis/03_export_supplementary_pseudobulk_source.R`

## Cell Ranger Reference Audit

The eight delivered Cell Ranger HTML reports were parsed structurally and locked by SHA256. All report Cell Ranger 5.0.0, `Single Cell 3' v3`, intronic-read counting, custom reference `Cal`, and transcriptome label `Cal-`; all eight estimated-cell values close to the Figure 2 alignment table.

The common delivered `genes.tsv` contains 26,134 reference features. The recovered project GFF3 contains 26,121 unique nuclear gene IDs, all present in that feature table with no GFF-only IDs. The remaining 13 IDs and symbols are exactly the 13 protein-coding genes from `NC_021614.1`. GFF sequence IDs and gene coordinates also close to the recovered nuclear FASTA. Evidence and source hashes are in `REFERENCE_PROVENANCE_AUDIT.md` and `metadata/reference_*`.

## Exact Raw 10x Audit

The eight delivered matrices were processed end to end on 2026-07-11. Exact reconstruction required:

- Seurat 5.3.0
- SeuratObject 5.2.0
- harmony 1.2.3
- Harmony batch variable `orig.ident`
- sample merge order recorded by `metadata/tenx_sample_sheet.csv`

The run was validated under R 4.6.1. The historical final object was serialized by R 4.4.0, but changing R itself was not required once the three core package versions were restored.

Results:

- Each delivered matrix contained 26,134 Cell Ranger reference features. `CreateSeuratObject(min.cells = 3)` retained a 23,740-feature union in the merged Seurat object.
- The 23,740 analysis features and 102,036 filtered nuclei matched the final annotated object.
- Sparse RNA count-matrix `Dim`, `i`, `p`, and `x` slots matched exactly.
- All public cell IDs matched.
- Cell-level `cluster_res04` comparison found 0 mismatches.
- All 27 resolution-0.4 cluster counts matched exactly.
- All 96 Figure 5a sample-by-cell-type counts matched; maximum proportion delta was `4.9960036108132044e-16`.
- All 87 marker rows selected for the Figure 3 marker panel matched exactly.

The exact checks are recorded in:

- `audit/annotation/raw_10x_full_rerun_audit.md`
- `audit/annotation/raw_10x_full_rerun_audit.csv`
- `audit/annotation/raw_10x_to_final_object_counts_validation.csv`
- `audit/annotation/raw_10x_to_final_object_sample_counts.csv`

## Exact Immune-Branch Audit

The immune workflow was rerun on 2026-07-11 from `GLOBAL_LEUKOCYTES_subset_for_reclustering.rds`, which is produced directly by the eight-matrix global reconstruction. It no longer requires an undocumented plotting table or preclustered object for the default path.

Using the same locked Seurat 5.3.0, SeuratObject 5.2.0, and harmony 1.2.3 environment:

- The raw-derived 10,440-cell initial immune object had the same public cell order, all 14 cluster assignments, sparse RNA count-matrix `Dim/i/p/x` slots, and UMAP coordinates as the historical `IMMUNEbaicell.rds` after `CON`/`COM` alias normalization.
- The final 5,061-cell object had identical cell order, all 10 immune subtype assignments, and UMAP coordinates as the exact-version historical branch.
- The final subtype counts were 1,110 T cell (CCR7+), 586 Activated T (RORA+), 723 Activated lymphoid (CCL20hi), 1,027 NK-like cytotoxic, 108 B cells, 298 cDC1-like (XCR1+), 640 Monocytes/macrophages, 219 MoDC-like (CD209d+), 242 Granulocyte-like, and 108 Cycling (G2/M).
- Every one of the 77 nonzero Figure 5b sample-by-subtype counts equals direct aggregation of the 5,061-row Figure 4 coordinate/annotation table.

An unpinned reconstruction had shifted one Slim1 cell from B cells to Granulocyte-like. The exact-version rerun restored B cells=18 and Granulocyte-like=58 for Slim1, matching the Figure 5b source and proving that this was package-version drift rather than an unexplained manual edit.

Audit files are under `audit/annotation/`, and the successful environment is recorded in `environment/immune_raw_branch_sessionInfo.txt`.

## Doublet Assessment Audit

Routine filtering on `nFeature_RNA`, `nCount_RNA`, and mitochondrial fraction was not misrepresented as doublet detection. A separate retrospective scDblFinder 1.26.7 analysis was run on raw RNA counts for all 102,036 nuclei, separately for each of the eight 10x libraries. The primary cluster-aware mode used the 27 global resolution-0.4 clusters; random artificial-doublet mode was used as a sensitivity analysis.

Results for nuclei in the 10,440-cell initial immune branch:

- Cluster-aware: 995/5,379 excluded nuclei (18.50%) versus 30/5,061 final nuclei (0.59%); odds ratio 38.06; Fisher exact P `4.64764115998859e-259`.
- Random sensitivity: 673/5,379 excluded nuclei (12.51%) versus 10/5,061 final nuclei (0.20%); odds ratio 72.19; Fisher exact P `4.37003961879297e-185`.
- The two modes jointly called 3,539 nuclei across the complete atlas.
- The merged table has 102,036 rows, 102,036 unique cell IDs, and no missing calls.

Calls were retained as a validation diagnostic and were not used for automatic deletion. The 5,379 excluded nuclei are not described as 5,379 doublets. Full parameters, per-cell calls, summaries, and session information are under `audit/doublet_assessment/`.

## Exact 10,440-to-5,061 Filtering Audit

The immune branch was independently reconstructed and audited as four cluster-level exclusions:

- Initial clusters 0 and 13: 10,440 to 6,782 nuclei; 3,658 excluded.
- Reclustered cluster 1: 6,782 to 6,006; 776 excluded.
- Reclustered clusters 8 and 11: 6,006 to 5,701; 305 excluded.
- Reclustered cluster 3: 5,701 to 5,061; 640 excluded.

These were marker- and cluster-identity exclusions, not a second low-gene-count filter. The excluded groups had equal or higher molecular complexity than the retained set. Evidence included a strong `CUBN`/`GBGT1`/`MYO5B` epithelial program, mixed immune/non-immune profiles with `LRP4`/`CPNE9`/`LAMA` expression, and low-`PTPRC`, `CD74`-high lineage-ambiguous clusters. Stage-level predicted-doublet rates were also substantially higher than in the final set.

The final reconstruction was hard-checked against the deposited Figure 4 table: 5,061 nuclei, exact cell set, exact annotation labels, and maximum absolute UMAP-coordinate delta `4.973799e-14`. The fate ledger, marker panel, QC table, stage-level doublet summary, flowchart, and UMAP audit are under `audit/immune_filtering/`. See `DOUBLETS_AND_IMMUNE_FILTERING_AUDIT.md` for the reviewer-facing interpretation.

## Version Guard

`00_reconstruct_global_atlas_from_10x.R`, `01_reconstruct_immune_final_object.R`, and `02c_export_scientific_data_figure4_source_data.R` stop by default when the validated reconstruction versions are absent. Run `environment/install_exact_r_packages.R` to create a repository-local R library. The scDblFinder script has a separate exact environment installed by `environment/install_doublet_r_packages.R` and rejects unvalidated versions unless `GROUPER_ALLOW_DOUBLET_VERSION_DRIFT=1` is explicitly set. `GROUPER_ALLOW_VERSION_DRIFT=1` remains limited to non-exact reconstruction runs.

Source URLs and SHA256 hashes for the three archived package releases are in `environment/legacy_r_package_sources.csv`. The successful run session is in `environment/raw_10x_rerun_sessionInfo.txt`.

## UMAP Boundary

The raw rerun reproduced filtered counts, cell identities, clusters, selected Figure 3 marker values, and Figure 5a exactly. UMAP coordinates were not numerically identical across the historical and current R/uwot stacks even though the graph partition was identical.

The actual final per-cell UMAP coordinates and annotations are therefore included as:

`source_data/figures/Fig3_global_umap_coordinates_final_annotation.csv.gz`

This is a processed cell-level source file, not a manually entered plotting table.

## Figure Assembly Verification

- R syntax parsing passed for all retained R scripts.
- Python syntax checks passed for all four Figure 2-5 assembly scripts, the release validator, and the final submission-artifact validator.
- Figure 2-5 assembly scripts ran successfully from a clean ASCII-only staging path.
- A fresh minimal conda environment created from `environment/figure_environment_win-64.yml` reproduced all four submitted PNG hashes. The strict preflight locks Python 3.13.9, direct package versions, FreeType 2.13.3, win-64, and the Arial regular/bold/italic font fingerprints.
- A separate pip-only clean-room test showed why build-level locking is necessary: the same top-level versions used a different embedded FreeType build and changed only rasterized text pixels, with Figure 5 tight-cropping differing by one pixel. Portable pip redraws remain supported but are not presented as byte-identical.
- PNG canvas sizes were Figure 2 `3205 x 2388`, Figure 3 `7080 x 6120`, Figure 4 `6900 x 4890`, and Figure 5 `7450 x 5434`.
- Figure 2-5 PDFs were regenerated twice consecutively with identical hashes, embedded `/CIDFontType2` TrueType fonts, no `/Type3` fonts, and no dynamic `/CreationDate` metadata.
- Relative to the prior Type 3 PDFs, Poppler 300 dpi differences affected only glyph contours (`0.73%-1.34%` of pixels); the submitted PNGs remained byte-identical and the quantitative graphics, axes, and layout were unchanged.
- The final Figure 2-5 PNGs were synchronized byte-for-byte into the final manuscript DOCX during release assembly; the manuscript is not tracked in this code repository.
- Figure 3 panel c now reads the 480-row percent-expressed/scaled-average-expression source directly. Differential-marker `avg_log2FC` values remain upstream evidence and are not displayed under an average-expression label.
- Figure 4 uses the exact 5,061-cell UMAP table and 260-row expression-summary table; `cDC1-like (XCR1+)` is the conservative cross-species label.
- Figure assembly leaves all included source-data files byte-for-byte unchanged by default.
- Figure 5a/Figure 5b proportions sum to 1.000000 by sample.
- Table 2 keeps global `Leukocytes` and contains no immune-subtype labels.
- Table 3, Figure 4, and Figure 5b use `cDC1-like (XCR1+)`.
- Public source data use `COM1/COM2`; `CON1/CON2` remain compatibility aliases only.
- `scripts/validate_release_package.py` automates release-wide source, count, label, manifest, syntax, canvas, PNG/PDF reference-hash, PDF font/metadata, and source-immutability checks.

Supplementary Figure verification added on 2026-07-14:

- The 102,036-row S1 QC table has unique public cell IDs and closes exactly to the final Figure 3 UMAP cell set and eight library counts.
- The S3 source contains 23,738 expressed features by eight library-level pseudobulks. Recalculation of all 64 Pearson coefficients from the deposited `log2(CPM + 1)` matrix differs from the plotting table by at most `1e-12`.
- S1-S3 SVG, PDF and 600 dpi PNG files were generated twice consecutively with identical hashes for all nine outputs.
- SVG text remains editable, PDF fonttype is 42, dynamic PDF creation dates are absent, and a fixed SVG hashsalt prevents random element IDs.
- Final PNG canvases are S1 `4389 x 3189`, S2 `4389 x 2079`, and S3 `4238 x 2280`.

## Excluded From The Lightweight Package

- FASTQ, 10x matrix directories, BAM files, and large Seurat/RDS objects.
- Generated validation figures under `outputs/figures/`.
- Internal planning, reviewer tokens, temporary scripts, and private links.
- Historical main/supplementary scripts outside the Figure 2-5 and S1-S3 scope.

## Upstream Reference Boundary

The verified end-to-end analysis begins with the eight GEO feature-barcode matrices. Cell Ranger 5.0.0 settings, source-report hashes, the custom-reference name, nuclear FASTA/GFF3 fingerprints, mitochondrial accession, feature-ID closure, 3' UTR extension statement and genome provenance are recorded. The original combined `cal_MT.fa` and exact 2 kb 3' UTR-extended `cal.3utr_MT.gtf` were not found, so FASTQ-to-matrix byte equivalence is not claimed. A future record version should add those derived files and SHA256 checksums, or a deterministic reference-build script validated against the historical reference and re-alignment output.

Large RDS hashes are retained in the audit document so an external object can be identified without placing it in GitHub.

## External-Input Environment Variables

- `GROUPER_10X_ROOT`
- `GROUPER_10X_SAMPLE_SHEET`
- `GROUPER_GLOBAL_CLUSTER_ANNOTATION_MAP`
- `GROUPER_HARMONY_BATCH_VAR`
- `GROUPER_ALLOW_VERSION_DRIFT`
- `GROUPER_GLOBAL_RECON_OUT`
- `GROUPER_SOURCE_DATA_OUT`
- `GROUPER_FINAL_GLOBAL_RDS`
- `GROUPER_FINAL_GLOBAL_AUDIT_OUT`
- `GROUPER_ANNOTATION_AUDIT_OUT`
- `GROUPER_IMMUNE_INPUT_RDS`
- `GROUPER_IMMUNE_INPUT_STAGE`
- `GROUPER_DOUBLET_OUT`
- `GROUPER_DOUBLET_SEED`
- `GROUPER_DOUBLET_DBR_PER_1K`
- `GROUPER_ALLOW_DOUBLET_VERSION_DRIFT`
- `GROUPER_DOUBLET_CALLS`
- `GROUPER_IMMUNE_FINAL_RDS`
- `GROUPER_GLOBAL_COMPOSITION_CSV`
- `GROUPER_IMMUNE_COMPOSITION_CSV`
- `GROUPER_FIGURE5_OUT_DIR`
- `GROUPER_WRITE_NORMALIZED_SOURCE`
