# Cell Ranger Reference Provenance Audit

## Scope and Verdict

This audit separates three reproducibility levels that must not be conflated:

1. The delivered Cell Ranger feature-barcode matrices and report settings.
2. The source genome, nuclear annotation and mitochondrial accession used to construct the custom reference.
3. The exact derived `cal_MT.fa` and `cal.3utr_MT.gtf` files supplied to `cellranger mkref`.

Levels 1 and 2 are now supported by machine-readable hashes and exact feature-ID closure. Level 3 remains incomplete because the original combined FASTA and the exact 2 kb 3' UTR-extended GTF were not present in the available project or sequencing-delivery directories. The release therefore supports exact matrix-to-figure reconstruction but does not claim byte-exact FASTQ-to-matrix realignment.

## Cell Ranger Report Evidence

Eight delivered Cell Ranger HTML web summaries were parsed structurally. All eight report:

- Cell Ranger `5.0.0`;
- `Single Cell 3' v3` chemistry;
- intronic-read counting enabled;
- custom reference basename `Cal`;
- transcriptome label `Cal-`.

The report-derived estimated-cell values close exactly to `metadata/samples.align.stat.xls` for every library. `metadata/cellranger_run_metadata.csv` records the report filename, byte count and SHA256 for each source report. The full HTML reports are omitted from this lightweight archive because their extracted settings and hashes are sufficient for provenance and the reports are not Figure 2-5 plotting inputs.

The web summaries do not expose the original `cellranger count` command line. They verify chemistry, intronic counting, reference identity, Cell Ranger version and estimated-cell outputs, but they do not independently prove the historical `--force-cells` command text.

## Reference Feature Closure

All eight delivered `genes.tsv` files are byte-identical:

- SHA256: `677f20a78964349c7697faab4f6f6006e23c9f927b1f3d639ab92116d8214a12`
- rows: `26,134`

The project nuclear GFF3 contains `26,121` unique gene records. Every one of those `26,121` gene IDs occurs in the delivered Cell Ranger feature table, with no GFF-only IDs. The remaining `13` Cell Ranger IDs are exactly the protein-coding mitochondrial genes from RefSeq accession `NC_021614.1`, and all 13 gene symbols agree with the current NCBI GFF3 record.

The exact ordered 26,134-row feature list is included as `metadata/reference_features_26134.csv.gz`. The summary checks are in `metadata/reference_feature_audit.csv`.

## Source Fingerprints

The identified nuclear source FASTA contains 266 records and 1,066,743,859 sequence bases. The nuclear GFF3 uses 97 of those sequence records; every gene coordinate is valid and within its matching FASTA sequence.

The project annotation archive `Cal.zip` contains a GFF3 entry byte-identical to the audited nuclear GFF3 and an annotation workbook byte-identical to the local `Cal-annotation.xlsx`. File-level SHA256/MD5 values, sequence-level fingerprints, record counts, source identifiers and retrieval URLs are recorded in `metadata/reference_source_manifest.csv`.

The Figshare DOI in that manifest is the source identified by Sun et al. for the nuclear genome assembly. The project nuclear GFF3 was recovered from `Cal.zip`; a separate public download for that exact GFF3 was not independently located, so the manifest does not claim one.

## 26,134 Versus 23,740 Features

These numbers refer to different stages:

- `26,134`: all features in each delivered Cell Ranger reference matrix.
- `23,740`: the union retained in the reconstructed Seurat object after each library was created with `CreateSeuratObject(min.cells = 3)` and then merged.

No reference genes were silently relabelled. The difference is a documented Seurat feature-presence filter, not a different reference genome.

## Re-running the Audit

The audit uses only the Python standard library. Supply the eight delivered gene tables, eight report HTML files, project FASTA/GFF3, current NCBI `NC_021614.1` FASTA/GFF3 and optional project annotation bundle:

```bash
python -B scripts/audit_reference_provenance.py \
  --report-dir /path/to/CellRanger_Report \
  --alignment-table metadata/samples.align.stat.xls \
  --genes-tsv /path/to/COM_1/genes.tsv \
  --genes-tsv /path/to/COM_2/genes.tsv \
  --genes-tsv /path/to/CTL_1/genes.tsv \
  --genes-tsv /path/to/CTL_2/genes.tsv \
  --genes-tsv /path/to/Llac_1/genes.tsv \
  --genes-tsv /path/to/Llac_2/genes.tsv \
  --genes-tsv /path/to/Slim_1/genes.tsv \
  --genes-tsv /path/to/Slim_2/genes.tsv \
  --nuclear-fasta /path/to/Cromileptes_altivelis.genome.fasta \
  --nuclear-gff3 /path/to/Cromileptes_altivelis.gff3 \
  --mitochondrial-fasta /path/to/NC_021614.1.fasta \
  --mitochondrial-gff3 /path/to/NC_021614.1.gff3 \
  --annotation-bundle /path/to/Cal.zip \
  --annotation-workbook /path/to/Cal-annotation.xlsx \
  --output-dir /path/to/audit_output
```

An exact FASTQ-to-matrix claim should be made only after the original `cal_MT.fa` and `cal.3utr_MT.gtf` are deposited with hashes, or after a deterministic 3' UTR-extension build script is validated against the historical reference and re-alignment outputs.
