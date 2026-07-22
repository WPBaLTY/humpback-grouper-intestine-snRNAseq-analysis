# Humpback grouper intestine snRNA-seq analysis code

This repository contains the code and validation workflow that support the Scientific Data Data Descriptor on a humpback grouper (*Cromileptes altivelis*) intestine single-nucleus RNA-seq atlas.

## Repository role

This is the **code and executable reproducibility repository**. It is not the formal processed-data record.

- Raw FASTQ files and per-library 10x matrices are in GEO/SRA under `GSE326285`.
- Processed data and large reusable objects are in Zenodo DOI `10.5281/zenodo.21481391`.
- This repository supplies scripts, environment records, run order and validators.
- The validators automatically read data from a sibling `zenodo_processed_data_record/` directory when the Zenodo package is unpacked next to this repository.

## Main contents

- `scripts/`: R analysis, Python figure assembly, source-data preparation and validators.
- `environment/`: R/Python environment records and install helpers.
- `metadata/`: code-archive metadata template and script manifest.
- `manifest/files_manifest.csv`: file sizes and SHA256 values.

## Quick validation

```bash
python -B scripts/validate_release_package.py
```

Run this after unpacking `humpback-grouper-intestine-snRNAseq-processed-data-record.zip` as a sibling directory named `zenodo_processed_data_record`.

For the full deterministic figure check on the tested Windows stack:

```bash
conda env create -f environment/figure_environment_win-64.yml
conda activate grouper-figure-win64
python -B scripts/validate_release_package.py --check-r --run-figures
```

Generated figure outputs are written to `outputs/figures/` and are ignored by Git.

## Boundary

The executable reconstruction begins at the deposited eight per-library 10x feature-barcode matrices. Byte-exact FASTQ-to-matrix realignment is not claimed because the historical combined FASTA and exact 2 kb 3-prime UTR-extended GTF are not available in this release.

## Citation and metadata

The file `metadata/zenodo_code_archive_metadata_template.json` is a template for archiving this code release, not active metadata for the Zenodo processed-data record. Replace its creator list with the final manuscript author order, affiliations and verified ORCIDs before creating a citable code archive.

Code is MIT licensed. Source-data tables and processed-data records in Zenodo are CC BY 4.0 unless otherwise stated by the final Zenodo record.
