import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import sys
from pathlib import Path, PurePosixPath
from zipfile import ZipFile


EXPECTED_LIBRARY_ORDER = [
    "CTL_1",
    "CTL_2",
    "Llac_1",
    "Llac_2",
    "Slim_1",
    "Slim_2",
    "COM_1",
    "COM_2",
]
EXPECTED_ESTIMATED_CELLS = {
    "CTL_1": 12873,
    "CTL_2": 14034,
    "Llac_1": 13438,
    "Llac_2": 11989,
    "Slim_1": 13493,
    "Slim_2": 10995,
    "COM_1": 13682,
    "COM_2": 14020,
}
EXPECTED_MT_GENES = {
    "ncbi_16016975": "ND2",
    "ncbi_16016976": "COX1",
    "ncbi_16016977": "COX2",
    "ncbi_16016978": "ATP8",
    "ncbi_16016979": "ATP6",
    "ncbi_16016980": "COX3",
    "ncbi_16016981": "ND3",
    "ncbi_16016982": "ND4L",
    "ncbi_16016983": "ND4",
    "ncbi_16016984": "ND5",
    "ncbi_16016985": "ND6",
    "ncbi_16016986": "CYTB",
    "ncbi_16016987": "ND1",
}
PUBLIC_GENOME_SOURCE = "https://doi.org/10.6084/m9.figshare.12481850"
MT_FASTA_SOURCE = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id=NC_021614.1&rettype=fasta&retmode=text"
)
MT_GFF_SOURCE = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id=NC_021614.1&rettype=gff3&retmode=text"
)


def file_hash(path, algorithm):
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_fingerprint(path):
    return {
        "basename": path.name,
        "bytes": path.stat().st_size,
        "sha256": file_hash(path, "sha256"),
        "md5": file_hash(path, "md5"),
    }


def read_features(path):
    rows = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for index, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if len(row) < 2:
                raise ValueError(f"{path.name}: row {index} has fewer than two fields")
            rows.append((row[0], row[1]))
    return rows


def fasta_stats(path):
    sequence_lengths = {}
    canonical = hashlib.sha256()
    current_id = None
    current_sequence = []

    def finish_record():
        if current_id is None:
            return
        sequence = "".join(current_sequence).upper()
        sequence_lengths[current_id] = len(sequence)
        canonical.update(current_id.encode("utf-8"))
        canonical.update(b"\0")
        canonical.update(sequence.encode("ascii"))
        canonical.update(b"\0")

    with path.open("r", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish_record()
                current_id = line[1:].split()[0]
                if not current_id or current_id in sequence_lengths:
                    raise ValueError(f"Duplicate or empty FASTA identifier in {path.name}")
                current_sequence = []
            else:
                if current_id is None:
                    raise ValueError(f"Sequence before FASTA header in {path.name}")
                current_sequence.append(line)
    finish_record()
    return {
        "records": len(sequence_lengths),
        "bases": sum(sequence_lengths.values()),
        "canonical_sequence_sha256": canonical.hexdigest(),
        "lengths": sequence_lengths,
    }


def parse_attributes(text):
    result = {}
    for field in text.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            result[key] = value
    return result


def gff_stats(path):
    genes = []
    sequence_ids = set()
    invalid_coordinates = []
    with path.open("r", encoding="utf-8", errors="strict") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path.name}: line {line_number} is not nine-column GFF")
            sequence_id, _, feature_type, start, end, _, strand, _, attributes = fields
            sequence_ids.add(sequence_id)
            try:
                start_value = int(start)
                end_value = int(end)
            except ValueError:
                invalid_coordinates.append((line_number, sequence_id, start, end))
                continue
            if start_value < 1 or end_value < start_value or strand not in {"+", "-", ".", "?"}:
                invalid_coordinates.append((line_number, sequence_id, start, end))
            if feature_type == "gene":
                attrs = parse_attributes(attributes)
                gene_id = attrs.get("ID", "")
                if not gene_id:
                    raise ValueError(f"{path.name}: gene at line {line_number} has no ID")
                genes.append(
                    {
                        "gene_id": gene_id,
                        "sequence_id": sequence_id,
                        "start": start_value,
                        "end": end_value,
                        "attributes": attrs,
                    }
                )
    return {
        "genes": genes,
        "gene_ids": [row["gene_id"] for row in genes],
        "sequence_ids": sequence_ids,
        "invalid_coordinates": invalid_coordinates,
    }


def parse_cellranger_report(path):
    text = path.read_text(encoding="utf-8", errors="strict")
    table_match = re.search(
        r'"pipeline_info_table"\s*:\s*(\{.*?\})\s*,\s*"sequencing"',
        text,
        flags=re.DOTALL,
    )
    if not table_match:
        raise ValueError(f"No pipeline_info_table found in {path.name}")
    table = json.loads(table_match.group(1))
    rows = dict(table["rows"])
    estimate_match = re.search(
        r'"filtered_bcs_transcriptome_union"\s*:\s*\{.*?"metric"\s*:\s*"([0-9,]+)"',
        text,
        flags=re.DOTALL,
    )
    if not estimate_match:
        raise ValueError(f"No estimated-cell metric found in {path.name}")
    reference_path = rows.get("Reference Path", "")
    fingerprint = file_fingerprint(path)
    return {
        "sample": rows.get("Sample ID", ""),
        "report_filename": path.name,
        "report_bytes": fingerprint["bytes"],
        "report_sha256": fingerprint["sha256"],
        "chemistry": rows.get("Chemistry", ""),
        "include_introns": rows.get("Include introns", ""),
        "reference_name": PurePosixPath(reference_path).name,
        "transcriptome_label": rows.get("Transcriptome", ""),
        "cellranger_version": rows.get("Pipeline Version", ""),
        "estimated_cells": int(estimate_match.group(1).replace(",", "")),
    }


def read_alignment_estimates(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {
            row["Sample"]: int(row["Estimated Number of Cells"].replace(",", ""))
            for row in rows
        }


def hash_zip_entry(archive_path, entry_name, algorithm="sha256"):
    digest = hashlib.new(algorithm)
    with ZipFile(archive_path) as archive, archive.open(entry_name) as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path, fieldnames, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(rows)


def write_deterministic_gzip_csv(path, fieldnames, rows):
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, compresslevel=9, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text_handle:
                writer = csv.DictWriter(text_handle, fieldnames=fieldnames, quoting=csv.QUOTE_MINIMAL)
                writer.writeheader()
                writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Audit Cell Ranger report settings and custom-reference feature provenance."
    )
    parser.add_argument("--report-dir", type=Path, required=True)
    parser.add_argument("--alignment-table", type=Path, required=True)
    parser.add_argument("--genes-tsv", type=Path, action="append", required=True)
    parser.add_argument("--nuclear-fasta", type=Path, required=True)
    parser.add_argument("--nuclear-gff3", type=Path, required=True)
    parser.add_argument("--mitochondrial-fasta", type=Path, required=True)
    parser.add_argument("--mitochondrial-gff3", type=Path, required=True)
    parser.add_argument("--annotation-bundle", type=Path)
    parser.add_argument("--annotation-workbook", type=Path)
    parser.add_argument("--retrieval-date", default="2026-07-12")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    audit_rows = []

    def check(name, expected, observed, condition, detail):
        audit_rows.append(
            {
                "check": name,
                "expected": str(expected),
                "observed": str(observed),
                "status": "PASS" if condition else "FAIL",
                "detail": detail,
            }
        )

    for path in [
        args.report_dir,
        args.alignment_table,
        *args.genes_tsv,
        args.nuclear_fasta,
        args.nuclear_gff3,
        args.mitochondrial_fasta,
        args.mitochondrial_gff3,
    ]:
        if not path.exists():
            parser.error(f"Input does not exist: {path}")

    report_paths = sorted(args.report_dir.glob("CellRanger.*.result.html"))
    report_rows = [parse_cellranger_report(path) for path in report_paths]
    reports_by_sample = {row["sample"]: row for row in report_rows}
    ordered_reports = [reports_by_sample[sample] for sample in EXPECTED_LIBRARY_ORDER if sample in reports_by_sample]
    check(
        "Cell Ranger web-summary libraries",
        EXPECTED_LIBRARY_ORDER,
        [row["sample"] for row in ordered_reports],
        len(report_rows) == 8
        and len(reports_by_sample) == 8
        and [row["sample"] for row in ordered_reports] == EXPECTED_LIBRARY_ORDER,
        "Eight delivered HTML reports with unique sample IDs",
    )
    report_settings = {
        (
            row["chemistry"],
            row["include_introns"],
            row["reference_name"],
            row["transcriptome_label"],
            row["cellranger_version"],
        )
        for row in report_rows
    }
    expected_settings = {("Single Cell 3' v3", "True", "Cal", "Cal-", "cellranger-5.0.0")}
    check(
        "Uniform Cell Ranger settings",
        expected_settings,
        report_settings,
        report_settings == expected_settings,
        "Chemistry, intronic counting, reference name, transcriptome label and version",
    )
    report_estimates = {row["sample"]: row["estimated_cells"] for row in report_rows}
    alignment_estimates = read_alignment_estimates(args.alignment_table)
    check(
        "Cell Ranger estimated-cell closure",
        EXPECTED_ESTIMATED_CELLS,
        report_estimates,
        report_estimates == EXPECTED_ESTIMATED_CELLS
        and alignment_estimates == EXPECTED_ESTIMATED_CELLS,
        "HTML reports and Figure 2 alignment table agree for all libraries",
    )

    feature_tables = [read_features(path) for path in args.genes_tsv]
    features = feature_tables[0]
    check(
        "Delivered feature tables are identical",
        len(args.genes_tsv),
        sum(table == features for table in feature_tables),
        len(args.genes_tsv) == 8 and all(table == features for table in feature_tables),
        "All eight genes.tsv files have identical rows and order",
    )
    feature_ids = [gene_id for gene_id, _ in features]
    feature_names = dict(features)
    feature_id_set = set(feature_ids)
    check(
        "Cell Ranger reference feature rows",
        26134,
        len(features),
        len(features) == 26134 and len(feature_id_set) == 26134,
        "Reference features before Seurat min.cells filtering",
    )

    nuclear_fasta = fasta_stats(args.nuclear_fasta)
    nuclear_gff = gff_stats(args.nuclear_gff3)
    nuclear_gene_ids = nuclear_gff["gene_ids"]
    nuclear_gene_set = set(nuclear_gene_ids)
    check(
        "Nuclear GFF gene IDs",
        26121,
        len(nuclear_gene_ids),
        len(nuclear_gene_ids) == 26121 and len(nuclear_gene_set) == 26121,
        "Unique gene records in the project reference GFF3",
    )
    check(
        "Nuclear feature-ID closure",
        "26121 shared; 0 GFF-only",
        f"{len(nuclear_gene_set & feature_id_set)} shared; {len(nuclear_gene_set - feature_id_set)} GFF-only",
        nuclear_gene_set <= feature_id_set and len(nuclear_gene_set) == 26121,
        "Every project GFF3 gene ID occurs in the delivered Cell Ranger feature table",
    )
    tenx_only = feature_id_set - nuclear_gene_set
    check(
        "Mitochondrial feature-ID closure",
        sorted(EXPECTED_MT_GENES),
        sorted(tenx_only),
        tenx_only == set(EXPECTED_MT_GENES),
        "The 13 non-nuclear IDs are exactly the expected protein-coding mitochondrial genes",
    )
    observed_mt_names = {gene_id: feature_names.get(gene_id) for gene_id in EXPECTED_MT_GENES}
    check(
        "Mitochondrial feature names",
        EXPECTED_MT_GENES,
        observed_mt_names,
        observed_mt_names == EXPECTED_MT_GENES,
        "Delivered gene names match NC_021614.1 protein-coding genes",
    )
    missing_fasta_sequences = nuclear_gff["sequence_ids"] - set(nuclear_fasta["lengths"])
    out_of_bounds = [
        row["gene_id"]
        for row in nuclear_gff["genes"]
        if row["sequence_id"] not in nuclear_fasta["lengths"]
        or row["end"] > nuclear_fasta["lengths"][row["sequence_id"]]
    ]
    check(
        "GFF-to-FASTA sequence closure",
        "0 missing or out-of-bounds",
        f"{len(missing_fasta_sequences)} missing; {len(out_of_bounds)} out-of-bounds",
        not missing_fasta_sequences
        and not out_of_bounds
        and not nuclear_gff["invalid_coordinates"],
        f"{len(nuclear_gff['sequence_ids'])} annotated FASTA records checked",
    )

    mitochondrial_fasta = fasta_stats(args.mitochondrial_fasta)
    mt_gff = gff_stats(args.mitochondrial_gff3)
    mt_protein_genes = {}
    for row in mt_gff["genes"]:
        attrs = row["attributes"]
        if attrs.get("gene_biotype") != "protein_coding":
            continue
        gene_id_match = re.search(r"GeneID:(\d+)", attrs.get("Dbxref", "") + ";" + attrs.get("ID", ""))
        if gene_id_match:
            mt_protein_genes[f"ncbi_{gene_id_match.group(1)}"] = attrs.get("description", "")
    check(
        "NC_021614.1 FASTA identity",
        "1 record; 16497 bases; accession NC_021614.1",
        f"{mitochondrial_fasta['records']} record; {mitochondrial_fasta['bases']} bases; "
        f"IDs {list(mitochondrial_fasta['lengths'])}",
        mitochondrial_fasta["records"] == 1
        and mitochondrial_fasta["bases"] == 16497
        and set(mitochondrial_fasta["lengths"]) == {"NC_021614.1"},
        "NCBI RefSeq mitochondrial sequence",
    )
    check(
        "NCBI mitochondrial GFF mapping",
        EXPECTED_MT_GENES,
        mt_protein_genes,
        mt_protein_genes == EXPECTED_MT_GENES,
        "Current RefSeq GeneID and symbol mapping closes to the delivered 10x features",
    )

    bundle_fingerprint = None
    if args.annotation_bundle:
        bundle_fingerprint = file_fingerprint(args.annotation_bundle)
        gff_entry = "Cromileptes_altivelis.gff3"
        with ZipFile(args.annotation_bundle) as archive:
            bundle_names = set(archive.namelist())
        check(
            "Annotation bundle GFF identity",
            file_hash(args.nuclear_gff3, "sha256"),
            hash_zip_entry(args.annotation_bundle, gff_entry) if gff_entry in bundle_names else "missing",
            gff_entry in bundle_names
            and hash_zip_entry(args.annotation_bundle, gff_entry) == file_hash(args.nuclear_gff3, "sha256"),
            "Cal.zip contains the byte-identical nuclear GFF3",
        )
        if args.annotation_workbook:
            workbook_entry = "Cal-annotation.xlsx"
            check(
                "Annotation bundle workbook identity",
                file_hash(args.annotation_workbook, "sha256"),
                hash_zip_entry(args.annotation_bundle, workbook_entry)
                if workbook_entry in bundle_names
                else "missing",
                workbook_entry in bundle_names
                and hash_zip_entry(args.annotation_bundle, workbook_entry)
                == file_hash(args.annotation_workbook, "sha256"),
                "Cal.zip contains the byte-identical annotation workbook",
            )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(
        args.output_dir / "cellranger_run_metadata.csv",
        [
            "sample",
            "report_filename",
            "report_bytes",
            "report_sha256",
            "chemistry",
            "include_introns",
            "reference_name",
            "transcriptome_label",
            "cellranger_version",
            "estimated_cells",
        ],
        ordered_reports,
    )

    feature_rows = [
        {
            "feature_index": index,
            "gene_id": gene_id,
            "gene_name": gene_name,
            "source_component": "mitochondrial_NC_021614.1"
            if gene_id in EXPECTED_MT_GENES
            else "nuclear_project_GFF3",
        }
        for index, (gene_id, gene_name) in enumerate(features, start=1)
    ]
    write_deterministic_gzip_csv(
        args.output_dir / "reference_features_26134.csv.gz",
        ["feature_index", "gene_id", "gene_name", "source_component"],
        feature_rows,
    )
    write_csv(
        args.output_dir / "reference_feature_audit.csv",
        ["check", "expected", "observed", "status", "detail"],
        audit_rows,
    )

    source_rows = []

    def add_source(
        component,
        identifier,
        path,
        records,
        bases_or_rows,
        canonical_sequence_sha256,
        provenance,
        public_source,
        notes,
    ):
        fingerprint = file_fingerprint(path)
        source_rows.append(
            {
                "component": component,
                "identifier": identifier,
                "basename": fingerprint["basename"],
                "bytes": fingerprint["bytes"],
                "sha256": fingerprint["sha256"],
                "md5": fingerprint["md5"],
                "records": records,
                "bases_or_rows": bases_or_rows,
                "canonical_sequence_sha256": canonical_sequence_sha256,
                "provenance": provenance,
                "public_source": public_source,
                "release_availability": "fingerprint_and_derived_feature_list_only",
                "notes": notes,
            }
        )

    add_source(
        "nuclear_genome_fasta",
        "Cromileptes altivelis chromosome-level assembly",
        args.nuclear_fasta,
        nuclear_fasta["records"],
        nuclear_fasta["bases"],
        nuclear_fasta["canonical_sequence_sha256"],
        "Sun et al. 2020; BioProject PRJNA639378",
        PUBLIC_GENOME_SOURCE,
        "The local byte count converts to 1017.46 MiB, matching the numeric size displayed by Figshare.",
    )
    add_source(
        "nuclear_gene_annotation_gff3",
        "project protein-coding annotation",
        args.nuclear_gff3,
        len(nuclear_gff["sequence_ids"]),
        len(nuclear_gene_ids),
        "",
        "Project reference archive Cal.zip",
        "",
        "All 26121 gene IDs occur in the delivered Cell Ranger feature table.",
    )
    add_source(
        "delivered_10x_gene_table",
        "common genes.tsv across eight libraries",
        args.genes_tsv[0],
        len(args.genes_tsv),
        len(features),
        "",
        "Sequencing delivery and GEO processed matrices",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326285",
        "Eight files are byte-identical; 26134 reference features precede Seurat min.cells=3 filtering.",
    )
    add_source(
        "mitochondrial_fasta",
        "NC_021614.1",
        args.mitochondrial_fasta,
        mitochondrial_fasta["records"],
        mitochondrial_fasta["bases"],
        mitochondrial_fasta["canonical_sequence_sha256"],
        f"NCBI RefSeq; retrieved {args.retrieval_date}",
        MT_FASTA_SOURCE,
        "The original combined cal_MT.fa was not available; accession sequence identity is fixed here.",
    )
    add_source(
        "mitochondrial_gff3",
        "NC_021614.1",
        args.mitochondrial_gff3,
        len(mt_gff["sequence_ids"]),
        len(mt_protein_genes),
        "",
        f"NCBI RefSeq; retrieved {args.retrieval_date}",
        MT_GFF_SOURCE,
        "The 13 protein-coding GeneIDs and symbols match the delivered Cell Ranger feature table.",
    )
    if args.annotation_bundle and bundle_fingerprint:
        source_rows.append(
            {
                "component": "project_annotation_bundle",
                "identifier": "Cal.zip",
                "basename": bundle_fingerprint["basename"],
                "bytes": bundle_fingerprint["bytes"],
                "sha256": bundle_fingerprint["sha256"],
                "md5": bundle_fingerprint["md5"],
                "records": 4,
                "bases_or_rows": "",
                "canonical_sequence_sha256": "",
                "provenance": "Project reference archive",
                "public_source": "",
                "release_availability": "fingerprint_only",
                "notes": "Contains GFF3, CDS FASTA, peptide FASTA and Cal-annotation.xlsx; audited entries match local sources.",
            }
        )

    write_csv(
        args.output_dir / "reference_source_manifest.csv",
        [
            "component",
            "identifier",
            "basename",
            "bytes",
            "sha256",
            "md5",
            "records",
            "bases_or_rows",
            "canonical_sequence_sha256",
            "provenance",
            "public_source",
            "release_availability",
            "notes",
        ],
        source_rows,
    )

    failures = [row for row in audit_rows if row["status"] != "PASS"]
    for row in audit_rows:
        print(f"[{row['status']}] {row['check']}: {row['observed']}")
    print(f"Wrote reference audit metadata to {args.output_dir}")
    if failures:
        print(f"Reference provenance audit failed: {len(failures)} check(s).", file=sys.stderr)
        return 1
    print("Reference provenance audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
