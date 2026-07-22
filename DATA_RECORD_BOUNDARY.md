# Data record boundary

This file defines what is and is not contained in this lightweight public release.

## Included here

- Code and environment records for downstream reconstruction from deposited 10x matrices.
- Metadata and source-data tables needed to rebuild Scientific Data Figures 2-5 and Supplementary Figures S1-S3.
- Compact cell-level ledgers for final annotation, doublet assessment and immune filtering.
- Reference provenance fingerprints and the ordered Cell Ranger feature table.
- Validators for package safety, source-data closure, figure regeneration and submission-artifact cross-checks.

## Hosted outside this lightweight tree

- Raw FASTQ files and processed per-library 10x matrices: GEO `GSE326285`.
- Large Seurat/RDS processed objects cited by the manuscript and Table 1: the existing Zenodo data record, or a new Zenodo version if the owner interface shows they are not inherited.
- Main manuscript DOCX, cover letter, reviewer links and administrative forms: journal submission package only.

## Zenodo upload rule

Do not treat the small Zenodo ZIP as a complete replacement for the full data record unless the Zenodo draft also contains the large processed objects promised in the manuscript. If the existing Zenodo record's large objects are inherited into the new draft, upload this ZIP as the reproducibility/code/source-data update. If they are not inherited, upload the required large processed objects separately or revise the manuscript, Table 1 and Data Availability wording together.

## Reproducibility boundary

The validated executable boundary starts from the deposited eight per-library 10x matrices. Byte-exact FASTQ-to-matrix realignment is not claimed because the historical combined FASTA and exact 2 kb 3-prime UTR-extended GTF are not available in this release.
