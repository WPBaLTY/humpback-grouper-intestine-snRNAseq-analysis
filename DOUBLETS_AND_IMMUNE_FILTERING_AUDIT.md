# Doublet assessment and immune-filtering audit

## Scope and interpretation

Filtering on detected genes (`nFeature_RNA`), UMI counts (`nCount_RNA`), and mitochondrial fraction is routine nucleus-level quality control. It is not a dedicated doublet assessment. The deposited workflow therefore reports routine QC, cluster-level immune cleanup, and retrospective doublet assessment as three separate operations.

The 5,379 nuclei excluded while refining the 10,440-nucleus immune branch must not be described as 5,379 doublets. They were excluded in four cluster-level steps because their expression profiles were epithelial-like, lineage-ambiguous, or mixed. Independent doublet assessment showed strong enrichment of predicted doublets in the excluded set and provides supporting evidence, not a replacement label for every excluded nucleus.

## Doublet method

- Input: raw RNA counts for all 102,036 nuclei in the final global object.
- Method: `scDblFinder` 1.26.7, applied separately to each of the eight 10x capture libraries.
- Primary mode: cluster-aware, using the 27 global resolution-0.4 clusters.
- Sensitivity mode: cluster-independent random artificial doublets.
- Parameters: `samples="sample"`, `multiSampleMode="split"`, `dbr.per1k=0.008`, and `BiocParallel::SerialParam`.
- Seeds: 20260711 for the primary mode and 20260712 for the sensitivity mode.
- Role in the analysis: retrospective validation. Calls were not used for automatic cell-level deletion from the submitted atlas.

## Doublet results

| Mode | Excluded from initial immune set | Final immune set | Odds ratio | Fisher exact P |
|---|---:|---:|---:|---:|
| Cluster-aware | 995/5,379 (18.50%) | 30/5,061 (0.59%) | 38.06 | 4.65 x 10^-259 |
| Random sensitivity | 673/5,379 (12.51%) | 10/5,061 (0.20%) | 72.19 | 4.37 x 10^-185 |

The two modes jointly called 3,539 nuclei as doublets across the 102,036-nucleus atlas. Within the final immune set, the largest primary-mode rates occurred in Cycling (G2/M), 7/108 (6.48%), and cDC1-like (XCR1+), 6/298 (2.01%). The sensitivity mode called 5/108 Cycling nuclei and 2/298 cDC1-like nuclei. These low residual calls were retained because the analysis is cluster- and marker-based, cycling programs can elevate RNA complexity, and isolated algorithmic calls do not by themselves justify deleting coherent biological states.

## Exact 10,440 to 5,061 transition

| Step | Cluster(s) excluded | N excluded | N retained | Cluster-aware doublet rate | Random doublet rate | Primary expression-based rationale |
|---|---|---:|---:|---:|---:|---|
| 1 | Initial 0 and 13 | 3,658 | 6,782 | 16.13% | 12.27% | Cluster 13 showed a strong lysosome-rich enterocyte/epithelial program (`CUBN`, `CUBN.1`, `CUBN.2`, `GBGT1.2`, `MYO5B`). Cluster 0 showed a broad mixed profile rather than a clean immune subtype. |
| 2 | Reclustered 1 | 776 | 6,006 | 19.33% | 11.60% | Mixed leukocyte/APC signal with discordant `LRP4`, `CPNE9`, and `LAMA3`-associated expression; no stable final lineage identity. |
| 3 | Reclustered 8 and 11 | 305 | 5,701 | 17.70% | 8.52% | `CD74`-high but very low-`PTPRC` lineage-ambiguous populations. Cluster 8 was APC-like (`CLEC10A.2`, `SPIC`); cluster 11 included `SLC6A4.1` and other discordant markers. |
| 4 | Reclustered 3 | 640 | 5,061 | 31.41% | 16.88% | Strongly mixed APC/myeloid and non-leukocyte-associated profile (`CD74`, `Axl`, `csf1r1`, `LRP4.1`, `CPNE9`, `LAMA5.1`) and the highest doublet enrichment among the four removed stages. |

These exclusions were not a second low-gene-count QC pass. Excluded clusters generally had equal or higher RNA complexity than the retained set. After the fourth step, reclustering of the former granulocyte/B-cell-containing cluster resolved 108 B cells and 242 granulocyte-like nuclei without changing the final total of 5,061.

## Final-object identity checks

The reconstructed final immune object was compared directly with the deposited Figure 4 per-cell source table:

- Final nuclei: 5,061.
- Cell-set match: true.
- Cell order after matching: true.
- Annotation mismatches: 0.
- Maximum absolute UMAP-coordinate difference: 4.97 x 10^-14.

## Reproduction files

- `scripts/R_analysis/00d_assess_doublets_scDblFinder.R`
- `scripts/R_analysis/01b_audit_immune_filtering_10440_to_5061.R`
- `environment/install_doublet_r_packages.R`
- `audit/doublet_assessment/`
- `audit/immune_filtering/`

The merged per-cell doublet table contains exactly 102,036 unique identifiers and no missing calls. The immune fate table contains exactly 10,440 identifiers, comprising 5,379 excluded and 5,061 retained nuclei.

## Reference

Germain, P.-L. et al. Doublet identification in single-cell sequencing data using scDblFinder. *F1000Research* 10, 979 (2022). https://doi.org/10.12688/f1000research.73600.2.
