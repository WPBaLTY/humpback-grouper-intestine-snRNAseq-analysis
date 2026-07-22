import argparse
import ast
import csv
import gzip
import hashlib
import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from pathlib import PurePosixPath
import shutil
import statistics
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT.parent / "zenodo_processed_data_record" if (ROOT.parent / "zenodo_processed_data_record").is_dir() else ROOT
DATA_METADATA = DATA_ROOT / "metadata"
ANNOTATION_AUDIT = DATA_ROOT / "audit" / "annotation"
SOURCE = DATA_ROOT / "source_data" / "figures"
DOUBLET = DATA_ROOT / "audit" / "doublet_assessment"
IMMUNE_AUDIT = DATA_ROOT / "audit" / "immune_filtering"

SAMPLES = ["COM1", "COM2", "CTL1", "CTL2", "Llac1", "Llac2", "Slim1", "Slim2"]
FIGURE_LIBRARY_IDS = ["CTL_1", "CTL_2", "Llac_1", "Llac_2", "Slim_1", "Slim_2", "COM_1", "COM_2"]
GLOBAL_SAMPLE_COUNTS = {
    "COM1": 13024,
    "COM2": 13824,
    "CTL1": 12647,
    "CTL2": 13767,
    "Llac1": 13286,
    "Llac2": 11803,
    "Slim1": 13140,
    "Slim2": 10545,
}
INTERNAL_SAMPLE_COUNTS = {
    ("CON" + sample[3:]) if sample.startswith("COM") else sample: count
    for sample, count in GLOBAL_SAMPLE_COUNTS.items()
}
GLOBAL_TYPES = {
    "Acinar-like",
    "Best4+ cell",
    "Endothelial cells",
    "Enterocytes",
    "Enteroendocrine cells",
    "Fibroblasts",
    "Goblet cells",
    "LREs",
    "Leukocytes",
    "Neuronal cells",
    "Smooth muscle",
    "Tuft cells",
}
IMMUNE_TYPES = {
    "T cell (CCR7+)",
    "Activated T (RORA+)",
    "Activated lymphoid (CCL20hi)",
    "NK-like cytotoxic",
    "B cells",
    "cDC1-like (XCR1+)",
    "Monocytes/macrophages",
    "MoDC-like (CD209d+)",
    "Granulocyte-like",
    "Cycling (G2/M)",
}
IMMUNE_TOTALS = {
    "T cell (CCR7+)": 1110,
    "Activated T (RORA+)": 586,
    "Activated lymphoid (CCL20hi)": 723,
    "NK-like cytotoxic": 1027,
    "B cells": 108,
    "cDC1-like (XCR1+)": 298,
    "Monocytes/macrophages": 640,
    "MoDC-like (CD209d+)": 219,
    "Granulocyte-like": 242,
    "Cycling (G2/M)": 108,
}
FINAL_UMAP_SHA256 = "07fda9640621eb290678815e13418b3fbbd5553cd5e4cbf95237ebd023190654"
EXPECTED_FIGURE_SHA256 = {
    "Figure2_QC_final.png": "e9664d83972b3034060ead4b1672e120020d5d08369db651592a9943ea3f5ed0",
    "Figure3_global_annotation_final.png": "ee9a9d9ed22e7eeee91a4c754bbb9a40e6e7f78117510447673daab2c5927292",
    "Figure4_immune_reclustering_final.png": "36c01c79a814a53108ed3c23f643f7d957dd6826c814d06120673a553a980ae9",
    "Figure5_composition_final.png": "dff85e9df8132666250221fa243515d2da4708550b2f6faba700b7ef93c3d329",
}
EXPECTED_FIGURE_PDF_SHA256 = {
    "Figure2_QC_final.pdf": "6114aaa6159e02ffb8b81c1018da13b19bffedeab136320803073620aa59764f",
    "Figure3_global_annotation_final.pdf": "b8d0dae355cfd959da04d95ced88f1b7d7d97ffa1d3838710e783a6618546d21",
    "Figure4_immune_reclustering_final.pdf": "5f327e1c456cc1f63676f33c220e209bbc8790f07e8129b8c0a3e66c94d9abea",
    "Figure5_composition_final.pdf": "cc1e2c97c4cfa777fa74cccd5508e06b78d6593198e5d5e4e09a2e648cf79437",
}
EXPECTED_SUPPLEMENTARY_SHA256 = {
    "Supplementary_Figure_S1_QC.svg": "002ffee1bbac7f105b2daa6911fb19de4736614122171fda8a1e86c0a4cd78a4",
    "Supplementary_Figure_S1_QC.pdf": "922a7f8cdc7284ed8def4af4c12e5411a83f35c882c40b9021436a8cb0d1572e",
    "Supplementary_Figure_S1_QC.png": "6d3b9cf869be2ebe194b10b7c9471d0a35d2f0713f509d7c0f20f8b714e9846b",
    "Supplementary_Figure_S2_UMAP.svg": "8e911ce772428ee1739421621917a9f3d1778f3bcd1b407be2dec7751e49e3aa",
    "Supplementary_Figure_S2_UMAP.pdf": "89c239e45de405a3b1e572510587a7095e23cceb91c08c5cccec7cc8d499ca88",
    "Supplementary_Figure_S2_UMAP.png": "e3113969381099b6c8b14208fe3aab9de0fb7cc9694d84ae3f2a3426ee78ffea",
    "Supplementary_Figure_S3_Representation.svg": "e4ba45c21f608813bc7ed4acb832f19a13b40357a953cde5476db40f91deb27d",
    "Supplementary_Figure_S3_Representation.pdf": "5dc8aff98a7561e6e8e40a60212cd08fe16eb7a0f544c97db6399273fcf3eef6",
    "Supplementary_Figure_S3_Representation.png": "a52ced6201ef3c9205d2bbfed68b004eca717d8ff4278323bd24655c62efcd17",
}
EXPECTED_FIGURE_PACKAGE_VERSIONS = {
    "numpy": "2.3.5",
    "pandas": "2.3.3",
    "matplotlib": "3.10.6",
    "seaborn": "0.13.2",
    "Pillow": "12.0.0",
}
EXPECTED_FIGURE_PYTHON_VERSION = (3, 13, 9)
EXPECTED_FREETYPE_VERSION = "2.13.3"
EXPECTED_ARIAL_SHA256 = {
    "regular": "b3658eadae55e682b5f69eb64c439c1ecc8f196c0bb8d4756d145d13bc86476a",
    "bold": "e8f4e3baf6cc35fed6fcce3a540e8b39e8f6cda1d22a28f2ec8f526fef7a43f5",
    "italic": "86b32db9a06f9694e2a3760c42e5117bcdc5cc1255bb5186ca8ce0305e22f288",
}

FORBIDDEN_PUBLIC_SUFFIXES = {
    ".rds", ".rda", ".rdata", ".fastq", ".fq", ".bam", ".bai", ".sam",
    ".cram", ".h5", ".h5ad", ".loom", ".mtx", ".tif", ".tiff", ".docx",
}
FORBIDDEN_PUBLIC_DIRECTORIES = {
    ".git", "__pycache__", "node_modules", ".r-library", ".r-library-doublet",
    ".venv", "venv", "tmp", "temp",
}
PUBLIC_TEXT_SUFFIXES = {
    ".md", ".txt", ".csv", ".json", ".py", ".r", ".gitignore", ".xls",
    ".yml", ".yaml",
}
LOCAL_PATH_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])(?:[A-Za-z]:[\\/]|\\\\[^\\/\s]+[\\/][^\\/\s]+)"
)
SECRET_PATTERNS = {
    "JWT": re.compile(
        r"\beyJ[A-Za-z0-9_-]{8,}\.[A-Za-z0-9_-]{8,}\.[A-Za-z0-9_-]{8,}\b"
    ),
    "credential query": re.compile(
        r"[?&](?:access[_-]?token|token|authkey|api[_-]?key)=[A-Za-z0-9._-]{10,}",
        re.IGNORECASE,
    ),
    "GitHub token": re.compile(r"\bgh[pousr]_[A-Za-z0-9_]{20,}\b"),
    "API secret": re.compile(r"\bsk-[A-Za-z0-9_-]{20,}\b"),
    "private key": re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----"),
}
EXPECTED_GEO_SRX = {
    "COM_1": ("GSM9627310", "SRX32711943"),
    "COM_2": ("GSM9627311", "SRX32711944"),
    "CTL_1": ("GSM9627312", "SRX32711945"),
    "CTL_2": ("GSM9627313", "SRX32711946"),
    "Llac_1": ("GSM9627314", "SRX32711947"),
    "Llac_2": ("GSM9627315", "SRX32711948"),
    "Slim_1": ("GSM9627316", "SRX32711949"),
    "Slim_2": ("GSM9627317", "SRX32711950"),
}
EXPECTED_CELLRANGER_REPORT_SHA256 = {
    "CTL_1": "362ef8ff59e77380a1c971b6047468e1bc7de165e31aa3cfcf955e996064b5c2",
    "CTL_2": "46ec3f05d58288a03e964cca660fd781aa2a279a0cdc9fe4e9f2ef99642b2f0c",
    "Llac_1": "8f1afdc92b498a51be333fab0cd7ee250b548178b70d5ece6f5134d7fa14d8fb",
    "Llac_2": "9d56166e280c774130b04b4c5689beeb99a688922751c3abbf88224d2557073f",
    "Slim_1": "6e73509d26c05ea2fa2f407c2e59dba43574c0e11df3ae909a11ac77be46d1be",
    "Slim_2": "30db24c77c76ae4538828bdb61dd3af664e4045c22d157ece3f62e5e1860f371",
    "COM_1": "5228e7ef6ae7b7dc8ac337971aafc9cac37d814d9eb243be2bcbc80cb774e5b3",
    "COM_2": "7d3e4a7c0ed1b0331e8c1d4fcb6cf3627c7660d8b912781bf7824d2006475e3c",
}
EXPECTED_REFERENCE_SOURCE_SHA256 = {
    "nuclear_genome_fasta": "d0bca3d294d2e1185ca84e64736bff0cee6cc919ddd517183d57688d4331956a",
    "nuclear_gene_annotation_gff3": "399371c9617667317eb70309bf1f86f6e1418c103ebc9d3298100ade52349161",
    "delivered_10x_gene_table": "677f20a78964349c7697faab4f6f6006e23c9f927b1f3d639ab92116d8214a12",
    "mitochondrial_fasta": "ccb1b2c38198370a6f925fb999f48fbfcf83cbbf07b133d070d0e0a426fd449f",
    "mitochondrial_gff3": "b91ca1f73e3cdfbccfc5b1f20713d0e75a38154d4c80870dbba9ffa945125316",
    "project_annotation_bundle": "730b2cab9aa9827b0e74649d737804e2d3f344931a375ee308e093b1fb72ec8a",
}
EXPECTED_REFERENCE_CANONICAL_SEQUENCE_SHA256 = {
    "nuclear_genome_fasta": "25a6bf3008f16594b5537c02595633ccf95091a018ba31b64abc332f5ed2e2eb",
    "mitochondrial_fasta": "57e35587112a5c127df0e82a8e8c39ea303686dbcaed22695f98c8dafb06c922",
}
EXPECTED_MT_FEATURES = {
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
IMMUNE_INTERNAL_TOTALS = {
    "T_cell_CCR7_like": 1110,
    "Activated_T_RORA_stress": 586,
    "Activated_lymphoid_CCL20_high": 723,
    "NK_like_cytotoxic": 1027,
    "B_cell": 108,
    "cDC1_XCR1_ZNF366": 298,
    "Macrophage_Monocyte_Axl_Csf1r1": 640,
    "MoDC_like_Cd209d_FN1_AOC3": 219,
    "Granulocyte_like_EPX_Ncf4_CYBB": 242,
    "Cycling_immune_G2M": 108,
}


class Audit:
    def __init__(self):
        self.failures = []
        self.records = []

    def check(self, condition, label, detail=""):
        self.records.append({
            "status": "PASS" if condition else "FAIL",
            "check": label,
            "detail": detail,
        })
        if condition:
            suffix = f": {detail}" if detail else ""
            print(f"[PASS] {label}{suffix}")
        else:
            self.failures.append(label)
            suffix = f": {detail}" if detail else ""
            print(f"[FAIL] {label}{suffix}")


def read_csv(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_tsv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def numeric_text(value):
    return float(str(value).replace(",", "").replace("%", ""))


def audit_file_manifest(audit, manifest_path, label):
    rows = read_csv(manifest_path)
    errors = []
    seen = set()
    for row in rows:
        relative_text = row.get("relative_path", "")
        relative = PurePosixPath(relative_text)
        if (
            not relative_text
            or relative.is_absolute()
            or ".." in relative.parts
            or "\\" in relative_text
            or relative_text in seen
        ):
            errors.append(f"unsafe/duplicate path: {relative_text}")
            continue
        seen.add(relative_text)
        target = ROOT.joinpath(*relative.parts)
        if not target.is_file():
            errors.append(f"missing: {relative_text}")
            continue
        try:
            expected_size = int(row["bytes"])
        except (KeyError, TypeError, ValueError):
            errors.append(f"invalid size: {relative_text}")
            continue
        if target.stat().st_size != expected_size:
            errors.append(f"size mismatch: {relative_text}")
        if sha256(target) != row.get("sha256", "").lower():
            errors.append(f"hash mismatch: {relative_text}")
    audit.check(
        bool(rows) and not errors,
        label,
        f"{len(rows)} listed files exact" if not errors else " | ".join(errors[:8]),
    )


def audit_public_release_safety(audit, manifest_path):
    rows = read_csv(manifest_path)
    targets = [(manifest_path.relative_to(ROOT).as_posix(), manifest_path)]
    for row in rows:
        relative_text = row.get("relative_path", "")
        relative = PurePosixPath(relative_text)
        if relative_text and not relative.is_absolute() and ".." not in relative.parts:
            targets.append((relative_text, ROOT.joinpath(*relative.parts)))

    problems = set()
    for relative_text, path in targets:
        if not path.is_file():
            continue
        relative = PurePosixPath(relative_text)
        suffix = path.suffix.lower()
        parts_lower = {part.lower() for part in relative.parts}
        if suffix in FORBIDDEN_PUBLIC_SUFFIXES:
            problems.add(f"forbidden extension: {relative_text}")
        if suffix == ".xlsx" and relative_text != "submission_support/Supplementary_Tables_S1-S4.xlsx":
            problems.add(f"unexpected workbook: {relative_text}")
        if parts_lower & FORBIDDEN_PUBLIC_DIRECTORIES:
            problems.add(f"forbidden directory: {relative_text}")
        if path.stat().st_size >= 100 * 1024 * 1024:
            problems.add(f"file at or above 100 MiB: {relative_text}")

        is_gzip_text = suffix == ".gz" and path.name.lower().endswith((".csv.gz", ".tsv.gz", ".txt.gz"))
        is_plain_text = suffix in PUBLIC_TEXT_SUFFIXES or path.name == ".gitignore"
        if not (is_gzip_text or is_plain_text):
            continue
        opener = gzip.open if is_gzip_text else open
        with opener(path, "rt", encoding="utf-8-sig", errors="ignore") as handle:
            for line in handle:
                if LOCAL_PATH_PATTERN.search(line):
                    problems.add(f"absolute local path: {relative_text}")
                for label, pattern in SECRET_PATTERNS.items():
                    if pattern.search(line):
                        problems.add(f"{label}: {relative_text}")

    audit.check(
        not problems,
        "Public release safety",
        (
            f"{len(targets)} manifest-controlled files; no forbidden data, local paths, "
            "credentials, caches, or >=100 MiB files"
            if not problems
            else " | ".join(sorted(problems)[:12])
        ),
    )


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def describe_path(path):
    for base in (ROOT, DATA_ROOT):
        try:
            return str(path.relative_to(base))
        except ValueError:
            pass
    return str(path)


def audit_release(audit):
    release_manifests = [ROOT / "manifest" / "files_manifest.csv"] if (ROOT / "manifest" / "files_manifest.csv").is_file() else []
    zenodo_manifests = [ROOT / "manifest" / "zenodo_files_manifest.csv"] if (ROOT / "manifest" / "zenodo_files_manifest.csv").is_file() else []
    audit.check(
        len(release_manifests) == 1,
        "Unique GitHub release manifest",
        release_manifests[0].name if len(release_manifests) == 1 else str([path.name for path in release_manifests]),
    )
    audit.check(
        len(zenodo_manifests) <= 1,
        "Unique optional Zenodo release manifest",
        zenodo_manifests[0].name if len(zenodo_manifests) == 1 else "not present in GitHub package",
    )
    if len(release_manifests) != 1 or len(zenodo_manifests) > 1:
        return
    release_manifest = release_manifests[0]
    zenodo_manifest = zenodo_manifests[0] if zenodo_manifests else None

    required = [
        ROOT / "metadata/zenodo_code_archive_metadata_template.json",
        release_manifest,
        ROOT / "DOUBLETS_AND_IMMUNE_FILTERING_AUDIT.md",
        ROOT / "REFERENCE_PROVENANCE_AUDIT.md",
        DATA_METADATA / "tenx_sample_sheet.csv",
        ROOT / "metadata" / "script_manifest.csv",
        DATA_METADATA / "figure_source_manifest.csv",
        DATA_METADATA / "geo_gsm_mapping.csv",
        DATA_METADATA / "sample_metadata.csv",
        DATA_METADATA / "samples.sequence.stat.xls",
        DATA_METADATA / "samples.align.stat.xls",
        DATA_METADATA / "sequencing_QC_metrics.csv",
        DATA_METADATA / "cellranger_run_metadata.csv",
        DATA_METADATA / "reference_source_manifest.csv",
        DATA_METADATA / "reference_feature_audit.csv",
        DATA_METADATA / "reference_features_26134.csv.gz",
        ROOT / "environment" / "install_doublet_r_packages.R",
        ROOT / "environment" / "figure_environment_win-64.yml",
        ROOT / "environment" / "python_requirements.txt",
        ROOT / "scripts" / "audit_reference_provenance.py",
        ROOT / "scripts" / "validate_submission_artifacts.py",
        ROOT / "scripts" / "R_analysis" / "00d_assess_doublets_scDblFinder.R",
        ROOT / "scripts" / "R_analysis" / "01b_audit_immune_filtering_10440_to_5061.R",
        SOURCE / "Fig3_global_umap_coordinates_final_annotation.csv.gz",
        SOURCE / "Fig3_global_markers_res0.4_all.csv",
        SOURCE / "Fig3_global_marker_dotplot_source.csv",
        SOURCE / "Fig3_global_marker_dotplot_gene_panel.csv",
        SOURCE / "Fig3_global_marker_dotplot_features_used.csv",
        SOURCE / "Fig4_immune_umap_coordinates.csv",
        SOURCE / "Fig4_immune_marker_dotplot_source.csv",
        SOURCE / "Fig5a_global_celltype_composition_leukocytes_merged.csv",
        SOURCE / "Fig5b_immune_subtype_composition.csv",
        SOURCE / "Table2_global_cell_type_annotations.csv",
        SOURCE / "Table3_leukocyte_reclustering_annotations.csv",
        DATA_ROOT / "source_data" / "celltype_counts_by_library_percentages.csv",
        DATA_ROOT / "supplementary_tables" / "Supplementary_Table_S1_sample_accession_mapping.csv",
        DATA_ROOT / "supplementary_tables" / "Supplementary_Table_S2_cell_metadata_schema.csv",
        DATA_ROOT / "supplementary_tables" / "Supplementary_Table_S3_annotation_evidence.csv",
        DATA_ROOT / "supplementary_tables" / "Supplementary_Table_S4_software_parameters.csv",
        DOUBLET / "doublet_assessment_per_cell.csv.gz",
        DOUBLET / "doublet_analysis_parameters.csv",
        DOUBLET / "doublet_enrichment_removed_5379_vs_final_5061.csv",
        IMMUNE_AUDIT / "immune_filtering_fate_10440_cells.csv.gz",
        IMMUNE_AUDIT / "immune_filtering_transition_ledger.csv",
        IMMUNE_AUDIT / "immune_filtering_final_object_validation.csv",
        IMMUNE_AUDIT / "immune_filtering_doublet_summary_by_removal_stage.csv",
        DATA_ROOT / "integrated_object" / "grouper_intestine_snRNAseq_seurat.rds",
        DATA_ROOT / "integrated_object" / "grouper_intestine_immune_reclustering_seurat.rds",
        ANNOTATION_AUDIT / "raw_10x_full_rerun_audit.csv",
        ANNOTATION_AUDIT / "final_global_export_summary.csv",
    ]
    missing = [describe_path(path) for path in required if not path.is_file()]
    audit.check(not missing, "Required release files", "none missing" if not missing else ", ".join(missing))
    if missing:
        return

    rds_expected = {
        "integrated_object/grouper_intestine_snRNAseq_seurat.rds": "cb5655900bcd2ac795e38e1371bf1d71053ee0f26a1ec9330651e82f605e7bd4",
        "integrated_object/grouper_intestine_immune_reclustering_seurat.rds": "b0a1d6bbf0aee048b8bfe1d960bc3041ce7e0954db871eeab86862ec2742867d",
    }
    rds_observed = {rel: sha256(DATA_ROOT.joinpath(*PurePosixPath(rel).parts)) for rel in rds_expected}
    audit.check(
        rds_observed == rds_expected,
        "Integrated Seurat RDS hashes",
        "two deposited RDS objects exact" if rds_observed == rds_expected else str(rds_observed),
    )

    figure_lock = (ROOT / "environment" / "figure_environment_win-64.yml").read_text(encoding="utf-8")
    required_render_specs = [
        "python=3.13.9=h260b955_100_cp313",
        "matplotlib=3.10.6=py313haa95532_1",
        "freetype=2.13.3=h0620614_0",
        "libpng=1.6.50=h46444df_0",
        "python-docx==1.2.0",
    ]
    audit.check(
        all(spec in figure_lock for spec in required_render_specs),
        "Figure render environment lock",
        "win-64 Python, matplotlib, FreeType, libpng and python-docx builds pinned",
    )

    audit_file_manifest(
        audit,
        release_manifest,
        "GitHub release file manifest",
    )
    if zenodo_manifest is not None:
        audit_file_manifest(audit, zenodo_manifest, "Zenodo release file manifest")
    public_manifest = zenodo_manifest if zenodo_manifest is not None else release_manifest
    audit_public_release_safety(audit, public_manifest)

    with (ROOT / "metadata/zenodo_code_archive_metadata_template.json").open("r", encoding="utf-8") as handle:
        zenodo_metadata = json.load(handle)
    audit.check(True, "Zenodo metadata JSON", "valid JSON")
    audit.check(
        zenodo_metadata.get("license") == "cc-by-4.0"
        and bool(zenodo_metadata.get("creators")),
        "Zenodo data license and creator field",
        "CC BY 4.0; creator list present",
    )

    sample_sheet = read_csv(DATA_METADATA / "tenx_sample_sheet.csv")
    audit.check(len(sample_sheet) == 8, "10x sample sheet row count", str(len(sample_sheet)))
    audit.check([row["sample"] for row in sample_sheet] == SAMPLES, "Original merge order", ", ".join(SAMPLES))
    public_sheet_values = [row["sample"] for row in sample_sheet]
    audit.check(not any(value.startswith("CON") for value in public_sheet_values), "Public sample names", "COM1/COM2")
    aliases = {row["folder_alias"] for row in sample_sheet if row["sample"].startswith("COM")}
    audit.check(aliases == {"CON1", "CON2"}, "Legacy aliases confined to sample sheet", ", ".join(sorted(aliases)))

    geo_rows = read_csv(DATA_METADATA / "geo_gsm_mapping.csv")
    observed_geo_srx = {
        row["library"]: (row["gsm_accession"], row["sra_experiment"])
        for row in geo_rows
    }
    audit.check(len(geo_rows) == 8, "GEO/SRX mapping row count", str(len(geo_rows)))
    audit.check(observed_geo_srx == EXPECTED_GEO_SRX, "GEO/GSM/SRX mapping", "all eight exact")

    sample_metadata = read_csv(DATA_METADATA / "sample_metadata.csv")
    sample_metadata_fields = set(sample_metadata[0]) if sample_metadata else set()
    audit.check(
        "Pooled_library_replicate" in sample_metadata_fields
        and "Biological_library_replicate" not in sample_metadata_fields,
        "Pooled-library replicate terminology",
        "no unsupported biological-replicate claim",
    )

    sequence_rows = read_tsv(DATA_METADATA / "samples.sequence.stat.xls")
    alignment_rows = read_tsv(DATA_METADATA / "samples.align.stat.xls")
    qc_rows = read_csv(DATA_METADATA / "sequencing_QC_metrics.csv")
    sequence_by_sample = {row["Sample"]: row for row in sequence_rows}
    alignment_by_sample = {row["Sample"]: row for row in alignment_rows}
    qc_by_sample = {row["Sample"]: row for row in qc_rows}
    audit.check(
        [row["Sample"] for row in sequence_rows]
        == [row["Sample"] for row in alignment_rows]
        == [row["Sample"] for row in qc_rows]
        == FIGURE_LIBRARY_IDS,
        "Figure 2 library order",
        ", ".join(FIGURE_LIBRARY_IDS),
    )
    figure2_field_map = {
        "Number_of_reads": ("sequence", "Number of Reads"),
        "Valid_barcodes_percent": ("sequence", "Valid Barcodes"),
        "Sequencing_saturation_percent": ("sequence", "Sequencing Saturation"),
        "Q30_RNA_read_percent": ("sequence", "Q30 Bases in RNA Read"),
        "Q30_UMI_percent": ("sequence", "Q30 Bases in UMI"),
        "Estimated_nuclei_CellRanger": ("alignment", "Estimated Number of Cells"),
        "Fraction_reads_in_cells_percent": ("alignment", "Fraction Reads in Cells"),
        "Median_genes_per_nucleus": ("alignment", "Median Genes per Cell"),
        "Median_UMI_counts_per_nucleus": ("alignment", "Median UMI Counts per Cell"),
        "Genome_mapping_percent": ("alignment", "Reads Mapped Confidently to Genome"),
        "Transcriptome_mapping_percent": ("alignment", "Reads Mapped Confidently to Transcriptome"),
    }
    figure2_values_match = True
    try:
        for sample in FIGURE_LIBRARY_IDS:
            qc_row = qc_by_sample[sample]
            for qc_field, (source_name, source_field) in figure2_field_map.items():
                source_row = sequence_by_sample[sample] if source_name == "sequence" else alignment_by_sample[sample]
                if not math.isclose(
                    numeric_text(qc_row[qc_field]),
                    numeric_text(source_row[source_field]),
                    rel_tol=0,
                    abs_tol=1e-9,
                ):
                    figure2_values_match = False
            if int(qc_row["Pooled_fish_per_library"]) != 3:
                figure2_values_match = False
    except (KeyError, TypeError, ValueError):
        figure2_values_match = False
    audit.check(
        len(sequence_rows) == len(alignment_rows) == len(qc_rows) == 8 and figure2_values_match,
        "Figure 2 source-table field closure",
        "all plotted metrics exact",
    )

    cellranger_rows = read_csv(DATA_METADATA / "cellranger_run_metadata.csv")
    report_hashes = {row["sample"]: row["report_sha256"] for row in cellranger_rows}
    report_estimates = {
        row["sample"]: int(row["estimated_cells"]) for row in cellranger_rows
    }
    alignment_estimates = {
        sample: int(numeric_text(row["Estimated Number of Cells"]))
        for sample, row in alignment_by_sample.items()
    }
    audit.check(
        len(cellranger_rows) == 8
        and [row["sample"] for row in cellranger_rows] == FIGURE_LIBRARY_IDS,
        "Cell Ranger report metadata libraries",
        ", ".join(row["sample"] for row in cellranger_rows),
    )
    audit.check(
        report_hashes == EXPECTED_CELLRANGER_REPORT_SHA256,
        "Cell Ranger source-report SHA256",
        "all eight locked",
    )
    audit.check(
        all(
            row["chemistry"] == "Single Cell 3' v3"
            and row["include_introns"] == "True"
            and row["reference_name"] == "Cal"
            and row["transcriptome_label"] == "Cal-"
            and row["cellranger_version"] == "cellranger-5.0.0"
            for row in cellranger_rows
        ),
        "Cell Ranger report settings",
        "chemistry, introns, Cal reference and version exact",
    )
    audit.check(
        report_estimates == alignment_estimates,
        "Cell Ranger report-to-Figure 2 closure",
        "all eight estimated-cell values exact",
    )

    reference_source_rows = read_csv(DATA_METADATA / "reference_source_manifest.csv")
    reference_source_by_component = {
        row["component"]: row for row in reference_source_rows
    }
    observed_reference_hashes = {
        component: row["sha256"]
        for component, row in reference_source_by_component.items()
    }
    audit.check(
        len(reference_source_rows) == 6
        and observed_reference_hashes == EXPECTED_REFERENCE_SOURCE_SHA256,
        "Reference source fingerprints",
        "six source components locked by SHA256",
    )
    observed_canonical_hashes = {
        component: reference_source_by_component.get(component, {}).get(
            "canonical_sequence_sha256", ""
        )
        for component in EXPECTED_REFERENCE_CANONICAL_SEQUENCE_SHA256
    }
    audit.check(
        observed_canonical_hashes == EXPECTED_REFERENCE_CANONICAL_SEQUENCE_SHA256,
        "Reference canonical sequence fingerprints",
        "nuclear and mitochondrial FASTA exact",
    )
    expected_reference_shapes = {
        "nuclear_genome_fasta": ("266", "1066743859"),
        "nuclear_gene_annotation_gff3": ("97", "26121"),
        "delivered_10x_gene_table": ("8", "26134"),
        "mitochondrial_fasta": ("1", "16497"),
        "mitochondrial_gff3": ("1", "13"),
        "project_annotation_bundle": ("4", ""),
    }
    observed_reference_shapes = {
        component: (row["records"], row["bases_or_rows"])
        for component, row in reference_source_by_component.items()
    }
    audit.check(
        observed_reference_shapes == expected_reference_shapes,
        "Reference record and feature counts",
        "266 FASTA records; 26,121 nuclear plus 13 mitochondrial genes",
    )

    reference_feature_audit = read_csv(
        DATA_METADATA / "reference_feature_audit.csv"
    )
    audit.check(
        len(reference_feature_audit) == 14
        and all(row["status"] == "PASS" for row in reference_feature_audit),
        "Reference provenance audit ledger",
        "14/14 checks passed",
    )
    reference_features = read_csv(
        DATA_METADATA / "reference_features_26134.csv.gz"
    )
    reference_feature_ids = {row["gene_id"] for row in reference_features}
    reference_origins = Counter(row["source_component"] for row in reference_features)
    reference_mt_features = {
        row["gene_id"]: row["gene_name"]
        for row in reference_features
        if row["source_component"] == "mitochondrial_NC_021614.1"
    }
    audit.check(
        len(reference_features) == len(reference_feature_ids) == 26134
        and [int(row["feature_index"]) for row in reference_features]
        == list(range(1, 26135)),
        "Ordered Cell Ranger reference features",
        "26,134 unique rows in exact source order",
    )
    audit.check(
        reference_origins
        == Counter({"nuclear_project_GFF3": 26121, "mitochondrial_NC_021614.1": 13}),
        "Reference feature source partition",
        str(dict(reference_origins)),
    )
    audit.check(
        reference_mt_features == EXPECTED_MT_FEATURES,
        "Reference mitochondrial IDs and symbols",
        "all 13 NC_021614.1 protein-coding genes exact",
    )

    raw_rerun_rows = read_csv(
        ANNOTATION_AUDIT / "raw_10x_full_rerun_audit.csv"
    )
    final_export_summary = {
        row["metric"]: row["value"]
        for row in read_csv(
            ANNOTATION_AUDIT / "final_global_export_summary.csv"
        )
    }
    global_reconstruction_script = (
        ROOT / "scripts" / "R_analysis" / "00_reconstruct_global_atlas_from_10x.R"
    ).read_text(encoding="utf-8")
    audit.check(
        len(raw_rerun_rows) == 8
        and {int(row["matrix_features"]) for row in raw_rerun_rows} == {26134}
        and final_export_summary.get("n_features") == "23740"
        and "min.cells = 3" in global_reconstruction_script,
        "26,134-to-23,740 feature-stage distinction",
        "Cell Ranger reference features versus Seurat min.cells=3 union",
    )
    reference_audit_text = (ROOT / "REFERENCE_PROVENANCE_AUDIT.md").read_text(
        encoding="utf-8"
    )
    audit.check(
        "does not claim byte-exact FASTQ-to-matrix realignment" in reference_audit_text
        and "`26,134`" in reference_audit_text
        and "`23,740`" in reference_audit_text,
        "Reference reproducibility boundary wording",
        "upstream limitation and feature stages explicit",
    )

    script_manifest = read_csv(ROOT / "metadata" / "script_manifest.csv")
    figure_manifest = read_csv(DATA_METADATA / "figure_source_manifest.csv")
    audit.check(len(script_manifest) == 20, "Script manifest row count", str(len(script_manifest)))
    audit.check(len(figure_manifest) == 17, "Figure source manifest row count", str(len(figure_manifest)))
    manifest_missing = [row["script"] for row in script_manifest if not (ROOT / row["script"]).is_file()]
    audit.check(not manifest_missing, "Manifest script paths", "all included" if not manifest_missing else ", ".join(manifest_missing))

    final_umap = SOURCE / "Fig3_global_umap_coordinates_final_annotation.csv.gz"
    audit.check(sha256(final_umap) == FINAL_UMAP_SHA256, "Final per-cell UMAP SHA256", FINAL_UMAP_SHA256)

    global_rows = read_csv(final_umap)
    global_ids = {row["cell_id"] for row in global_rows}
    global_samples = Counter(row["sample"] for row in global_rows)
    global_clusters = {row["cluster_res04"] for row in global_rows}
    global_types = {row["global_cell_type"] for row in global_rows}
    global_counts = Counter((row["sample"], row["global_cell_type"]) for row in global_rows)
    audit.check(len(global_rows) == 102036, "Final global cell rows", str(len(global_rows)))
    audit.check(len(global_ids) == 102036, "Unique final global cell IDs", str(len(global_ids)))
    audit.check(global_samples == Counter(GLOBAL_SAMPLE_COUNTS), "Global per-sample counts", str(dict(global_samples)))
    audit.check(global_clusters == {str(value) for value in range(27)}, "Resolution-0.4 clusters", str(len(global_clusters)))
    audit.check(global_types == GLOBAL_TYPES, "Global annotation categories", str(len(global_types)))
    audit.check(not any(row["sample"].startswith("CON") or row["cell_id"].startswith("CON") for row in global_rows), "Global public identifiers", "no CON-prefixed public IDs")

    fig5a = read_csv(SOURCE / "Fig5a_global_celltype_composition_leukocytes_merged.csv")
    fig5a_counts = {(row["sample"], row["cell_type"]): int(row["n"]) for row in fig5a}
    fig5a_sums = defaultdict(float)
    for row in fig5a:
        fig5a_sums[row["sample"]] += float(row["prop"])
    audit.check(len(fig5a) == 96, "Figure 5a source rows", str(len(fig5a)))
    audit.check(fig5a_counts == dict(global_counts), "Figure 5a counts versus final per-cell table", "all 96 exact")
    audit.check(all(math.isclose(value, 1.0, abs_tol=1e-12) for value in fig5a_sums.values()), "Figure 5a sample proportions", "all sum to 1")

    immune_rows = read_csv(SOURCE / "Fig4_immune_umap_coordinates.csv")
    immune_ids = {row["cell_id"] for row in immune_rows}
    immune_counts = Counter((row["sample"], row["immune_subtype"]) for row in immune_rows)
    immune_totals = Counter(row["immune_subtype"] for row in immune_rows)
    audit.check(len(immune_rows) == 5061, "Final immune cell rows", str(len(immune_rows)))
    audit.check(len(immune_ids) == 5061, "Unique final immune cell IDs", str(len(immune_ids)))
    audit.check(set(immune_totals) == IMMUNE_TYPES, "Immune subtype categories", str(len(immune_totals)))
    audit.check(immune_totals == Counter(IMMUNE_TOTALS), "Immune subtype totals", str(dict(immune_totals)))
    audit.check(not any(row["sample"].startswith("CON") or row["cell_id"].startswith("CON") for row in immune_rows), "Immune public identifiers", "no CON-prefixed public IDs")

    fig5b = read_csv(SOURCE / "Fig5b_immune_subtype_composition.csv")
    fig5b_counts = {(row["sample"], row["subtype"]): int(row["n"]) for row in fig5b}
    fig5b_sums = defaultdict(float)
    for row in fig5b:
        fig5b_sums[row["sample"]] += float(row["prop"])
    audit.check(len(fig5b) == 77, "Figure 5b source rows", str(len(fig5b)))
    audit.check(fig5b_counts == dict(immune_counts), "Figure 5b counts versus Figure 4 per-cell table", "all 77 exact")
    audit.check(all(math.isclose(value, 1.0, abs_tol=1e-12) for value in fig5b_sums.values()), "Figure 5b sample proportions", "all sum to 1")

    fig3_dot = read_csv(SOURCE / "Fig3_global_marker_dotplot_source.csv")
    fig3_panel = read_csv(SOURCE / "Fig3_global_marker_dotplot_gene_panel.csv")
    fig3_features = read_csv(SOURCE / "Fig3_global_marker_dotplot_features_used.csv")
    fig3_markers = read_csv(SOURCE / "Fig3_global_markers_res0.4_all.csv")
    fig4_dot = read_csv(SOURCE / "Fig4_immune_marker_dotplot_source.csv")
    fig4_features = read_csv(SOURCE / "Fig4_immune_marker_features_used.csv")
    fig3_dot_pairs = {(row["cell_type"], row["marker_gene"].lower()) for row in fig3_dot}
    fig3_dot_genes = {row["marker_gene"].lower() for row in fig3_dot}
    fig3_feature_genes = {row["feature"].lower() for row in fig3_features}
    try:
        fig3_numeric_valid = all(
            0.0 <= float(row["percent_expressed"]) <= 100.0
            and math.isfinite(float(row["average_expression"]))
            and math.isfinite(float(row["average_expression_scaled"]))
            for row in fig3_dot
        )
    except (KeyError, TypeError, ValueError):
        fig3_numeric_valid = False
    audit.check(len(fig3_dot) == 480, "Figure 3 expression-summary rows", "12 categories x 40 genes")
    audit.check(
        len(fig3_dot_pairs) == 480
        and {row["cell_type"] for row in fig3_dot} == GLOBAL_TYPES
        and len(fig3_dot_genes) == 40,
        "Figure 3 expression-summary dimensions",
        "unique cell-type/marker combinations",
    )
    audit.check(fig3_dot_genes == fig3_feature_genes, "Figure 3 marker feature list", "all 40 exact")
    audit.check(fig3_numeric_valid, "Figure 3 expression-summary values", "finite and percent in [0, 100]")
    audit.check(len(fig3_panel) == 480, "Figure 3 marker-panel documentation", "12 categories x 40 genes")
    audit.check(len(fig3_markers) == 13387, "Figure 3 full marker rows", str(len(fig3_markers)))
    fig4_dot_pairs = {(row["immune_subtype"], row["marker_gene"].lower()) for row in fig4_dot}
    fig4_dot_genes = {row["marker_gene"].lower() for row in fig4_dot}
    fig4_feature_genes = {row["feature"].lower() for row in fig4_features}
    try:
        fig4_numeric_valid = all(
            0.0 <= float(row["percent_expressed"]) <= 100.0
            and math.isfinite(float(row["average_expression"]))
            and math.isfinite(float(row["average_expression_scaled"]))
            for row in fig4_dot
        )
    except (KeyError, TypeError, ValueError):
        fig4_numeric_valid = False
    audit.check(
        len(fig4_dot) == len(fig4_dot_pairs) == 260,
        "Figure 4 dot-plot source",
        "10 subtypes x 26 unique genes",
    )
    audit.check({row["immune_subtype"] for row in fig4_dot} == IMMUNE_TYPES, "Figure 4 dot-plot subtype labels", "exact")
    audit.check(fig4_dot_genes == fig4_feature_genes, "Figure 4 marker feature list", "all 26 exact")
    audit.check(fig4_numeric_valid, "Figure 4 expression-summary values", "finite and percent in [0, 100]")

    figure3_script = (
        ROOT / "scripts" / "python_figure_assembly" / "rebuild_scientific_data_figure3_global_annotation.py"
    ).read_text(encoding="utf-8")
    audit.check(
        "Fig3_global_marker_dotplot_source.csv" in figure3_script
        and "avg_log2FC" not in figure3_script
        and '["pct.1"]' not in figure3_script,
        "Figure 3 plotting variable semantics",
        "direct percent-expressed and scaled-average-expression source",
    )

    table2 = read_csv(SOURCE / "Table2_global_cell_type_annotations.csv")
    table3 = read_csv(SOURCE / "Table3_leukocyte_reclustering_annotations.csv")
    table2_types = {row["Cell type"] for row in table2}
    table3_types = {row["Cell type"] for row in table3}
    audit.check("Leukocytes" in table2_types and not (table2_types & IMMUNE_TYPES), "Table 2 annotation layer", "global Leukocytes only")
    audit.check(table3_types == IMMUNE_TYPES, "Table 3 annotation layer", "10 immune subtypes")

    doublet_rows = read_csv(DOUBLET / "doublet_assessment_per_cell.csv.gz")
    doublet_ids = {row["cell_id"] for row in doublet_rows}
    doublet_samples = Counter(row["sample"] for row in doublet_rows)
    analysis_sets = Counter(row["analysis_set"] for row in doublet_rows)
    audit.check(len(doublet_rows) == 102036, "Doublet per-cell rows", str(len(doublet_rows)))
    audit.check(len(doublet_ids) == 102036, "Unique doublet cell IDs", str(len(doublet_ids)))
    audit.check(doublet_samples == Counter(INTERNAL_SAMPLE_COUNTS), "Doublet per-library rows", str(dict(doublet_samples)))
    audit.check(
        analysis_sets
        == Counter({
            "outside_initial_immune": 91596,
            "removed_from_initial_immune_5379": 5379,
            "final_immune_5061": 5061,
        }),
        "Doublet analysis-set partition",
        str(dict(analysis_sets)),
    )
    audit.check(
        all(
            row["cluster_aware_class"] in {"singlet", "doublet"}
            and row["random_class"] in {"singlet", "doublet"}
            and row["cluster_aware_score"] != ""
            and row["random_score"] != ""
            for row in doublet_rows
        ),
        "Complete doublet calls and scores",
        "no missing values",
    )
    cluster_doublets = Counter(
        row["analysis_set"] for row in doublet_rows if row["cluster_aware_class"] == "doublet"
    )
    random_doublets = Counter(
        row["analysis_set"] for row in doublet_rows if row["random_class"] == "doublet"
    )
    audit.check(
        cluster_doublets
        == Counter({
            "outside_initial_immune": 4962,
            "removed_from_initial_immune_5379": 995,
            "final_immune_5061": 30,
        }),
        "Cluster-aware doublet totals",
        str(dict(cluster_doublets)),
    )
    audit.check(
        random_doublets
        == Counter({
            "outside_initial_immune": 4689,
            "removed_from_initial_immune_5379": 673,
            "final_immune_5061": 10,
        }),
        "Random-mode doublet totals",
        str(dict(random_doublets)),
    )
    concordant_doublets = sum(
        row["cluster_aware_class"] == "doublet" and row["random_class"] == "doublet"
        for row in doublet_rows
    )
    audit.check(concordant_doublets == 3539, "Concordant doublet calls", str(concordant_doublets))

    fate_rows = read_csv(IMMUNE_AUDIT / "immune_filtering_fate_10440_cells.csv.gz")
    fate_ids = {row["cell_id"] for row in fate_rows}
    fate_public_ids = {row["public_cell_id"] for row in fate_rows}
    disposition = Counter(row["disposition"] for row in fate_rows)
    removal_stages = Counter(row["removal_stage"] for row in fate_rows)
    final_internal = Counter(
        row["final_immune_annotation"]
        for row in fate_rows
        if row["disposition"] == "retained_final_5061"
    )
    audit.check(len(fate_rows) == 10440, "Immune fate rows", str(len(fate_rows)))
    audit.check(len(fate_ids) == 10440 and len(fate_public_ids) == 10440, "Unique immune fate IDs", "internal and public")
    audit.check(
        disposition == Counter({"excluded": 5379, "retained_final_5061": 5061}),
        "Immune fate disposition",
        str(dict(disposition)),
    )
    audit.check(
        removal_stages
        == Counter({
            "NA": 5061,
            "step1_drop_initial_clusters_0_13": 3658,
            "step2_drop_reclustered_cluster_1": 776,
            "step3_drop_reclustered_clusters_8_11": 305,
            "step4_drop_reclustered_cluster_3": 640,
        }),
        "Immune removal-stage ledger",
        str(dict(removal_stages)),
    )
    audit.check(final_internal == Counter(IMMUNE_INTERNAL_TOTALS), "Final internal subtype totals", str(dict(final_internal)))
    retained_fate = {
        row["public_cell_id"]: row["final_immune_annotation"]
        for row in fate_rows
        if row["disposition"] == "retained_final_5061"
    }
    deposited_immune = {
        row["cell_id"]: row["immune_subtype_internal"] for row in immune_rows
    }
    audit.check(retained_fate == deposited_immune, "Fate ledger versus Figure 4 cells and labels", "all 5,061 exact")

    transition_rows = read_csv(IMMUNE_AUDIT / "immune_filtering_transition_ledger.csv")
    observed_transitions = [
        (
            int(row["input_nuclei"]),
            int(row["excluded_nuclei"]),
            int(row["retained_nuclei"]),
        )
        for row in transition_rows
    ]
    audit.check(
        observed_transitions
        == [(10440, 3658, 6782), (6782, 776, 6006), (6006, 305, 5701), (5701, 640, 5061)],
        "10,440-to-5,061 transition ledger",
        str(observed_transitions),
    )

    final_validation = {
        row["check"]: row["value"]
        for row in read_csv(IMMUNE_AUDIT / "immune_filtering_final_object_validation.csv")
    }
    audit.check(final_validation.get("final_nuclei") == "5061", "Final immune validation count", "5,061")
    audit.check(
        final_validation.get("cell_set_match") == "TRUE"
        and final_validation.get("cell_order_after_matching") == "TRUE"
        and final_validation.get("annotation_mismatches") == "0",
        "Final immune identity assertions",
        "cell set, order, and labels exact",
    )
    audit.check(
        float(final_validation.get("maximum_absolute_umap_delta", "inf")) < 1e-12,
        "Final immune UMAP delta",
        final_validation.get("maximum_absolute_umap_delta", "missing"),
    )

    python_files = sorted((ROOT / "scripts").rglob("*.py"))
    parse_errors = []
    for path in python_files:
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except SyntaxError as error:
            parse_errors.append(f"{path.relative_to(ROOT)}: {error}")
    audit.check(not parse_errors, "Python syntax", f"{len(python_files)} files parsed" if not parse_errors else "; ".join(parse_errors))


def check_r_syntax(audit):
    rscript = shutil.which("Rscript")
    if not rscript:
        audit.check(False, "R syntax", "Rscript was not found on PATH")
        return
    errors = []
    r_files = sorted((ROOT / "scripts" / "R_analysis").glob("*.R")) + sorted(
        (ROOT / "environment").glob("install*_r_packages.R")
    )
    for path in r_files:
        expression = f"parse(file={json.dumps(path.as_posix())})"
        result = subprocess.run([rscript, "-e", expression], capture_output=True, text=True)
        if result.returncode:
            errors.append(f"{path.relative_to(ROOT)}: {result.stderr.strip()}")
    audit.check(not errors, "R syntax", f"{len(r_files)} files parsed" if not errors else " | ".join(errors))


def pearson(left, right):
    left_mean = statistics.fmean(left)
    right_mean = statistics.fmean(right)
    numerator = sum((a - left_mean) * (b - right_mean) for a, b in zip(left, right))
    left_ss = sum((value - left_mean) ** 2 for value in left)
    right_ss = sum((value - right_mean) ** 2 for value in right)
    return numerator / math.sqrt(left_ss * right_ss)


def audit_supplementary_sources(audit):
    qc_path = SOURCE / "SuppFigureS1_qc_per_cell.csv.gz"
    qc_rows = read_csv(qc_path)
    qc_ids = [row["cell_id"] for row in qc_rows]
    audit.check(len(qc_rows) == 102036, "Supplement S1 QC rows", f"{len(qc_rows):,}")
    audit.check(len(set(qc_ids)) == 102036, "Supplement S1 unique cell IDs", f"{len(set(qc_ids)):,}")
    qc_counts = Counter(row["sample"] for row in qc_rows)
    audit.check(qc_counts == Counter(GLOBAL_SAMPLE_COUNTS), "Supplement S1 per-library counts", str(dict(qc_counts)))
    qc_numeric_ok = all(
        float(row["nFeature_RNA"]) > 0
        and float(row["nCount_RNA"]) > 0
        and 0 <= float(row["percent_mt"]) <= 2
        for row in qc_rows
    )
    audit.check(qc_numeric_ok, "Supplement S1 retained-QC numeric ranges", "positive genes/UMIs; mitochondrial fraction 0-2%")

    global_rows = read_csv(SOURCE / "Fig3_global_umap_coordinates_final_annotation.csv.gz")
    audit.check(
        set(qc_ids) == {row["cell_id"] for row in global_rows},
        "Supplement S1 versus final UMAP cell closure",
        "102,036 exact public cell IDs",
    )
    summary_rows = read_csv(SOURCE / "SuppFigureS1_qc_summary.csv")
    summary_counts = {row["sample"]: int(row["retained_nuclei"]) for row in summary_rows}
    audit.check(
        len(summary_rows) == 8 and summary_counts == GLOBAL_SAMPLE_COUNTS,
        "Supplement S1 summary closure",
        "eight libraries exact",
    )

    log_cpm_rows = read_csv(SOURCE / "SuppFigureS3_pseudobulk_logCPM.csv")
    audit.check(len(log_cpm_rows) == 23738, "Supplement S3 expressed features", f"{len(log_cpm_rows):,}")
    log_columns_ok = bool(log_cpm_rows) and set(log_cpm_rows[0]) == {"feature_id", *FIGURE_LIBRARY_IDS}
    audit.check(log_columns_ok, "Supplement S3 pseudobulk columns", str(FIGURE_LIBRARY_IDS))
    vectors = {
        sample: [float(row[sample]) for row in log_cpm_rows]
        for sample in FIGURE_LIBRARY_IDS
    }
    finite_ok = all(math.isfinite(value) and value >= 0 for values in vectors.values() for value in values)
    audit.check(finite_ok, "Supplement S3 logCPM numeric values", "finite and nonnegative")

    correlation_rows = read_csv(SOURCE / "SuppFigureS3_pseudobulk_correlation.csv")
    observed_correlation = {
        row["sample_plot"]: {sample: float(row[sample]) for sample in FIGURE_LIBRARY_IDS}
        for row in correlation_rows
    }
    matrix_shape_ok = set(observed_correlation) == set(FIGURE_LIBRARY_IDS)
    audit.check(matrix_shape_ok, "Supplement S3 correlation matrix shape", "8 x 8")
    max_difference = 0.0
    if matrix_shape_ok:
        for left in FIGURE_LIBRARY_IDS:
            for right in FIGURE_LIBRARY_IDS:
                recalculated = pearson(vectors[left], vectors[right])
                max_difference = max(max_difference, abs(recalculated - observed_correlation[left][right]))
    audit.check(
        matrix_shape_ok and max_difference <= 1e-12,
        "Supplement S3 correlation recomputation",
        f"maximum absolute difference {max_difference:.3e}",
    )
    off_diagonal = [
        observed_correlation[left][right]
        for left in FIGURE_LIBRARY_IDS
        for right in FIGURE_LIBRARY_IDS
        if left != right
    ]
    audit.check(
        min(off_diagonal) > 0.97 and max(off_diagonal) < 1.0,
        "Supplement S3 off-diagonal correlation range",
        f"{min(off_diagonal):.6f}-{max(off_diagonal):.6f}",
    )

    parameters = {
        row["parameter"]: row["value"]
        for row in read_csv(SOURCE / "SuppFigureS3_pseudobulk_parameters.csv")
    }
    expected_parameters = {
        "input_features": "23740",
        "input_nuclei": "102036",
        "retained_expressed_features": "23738",
        "aggregation": "sum raw RNA counts by capture library",
        "normalization": "counts per million using total pseudobulk library size",
        "transformation": "log2(CPM + 1)",
        "correlation": "Pearson correlation across retained expressed features",
    }
    audit.check(parameters == expected_parameters, "Supplement S3 parameter record", "exact")


def check_figure_environment(audit):
    try:
        import matplotlib
        import matplotlib.ft2font as ft2font
        from matplotlib import font_manager
        import numpy
        import pandas
        import PIL
        import seaborn
    except Exception as error:
        audit.check(False, "Figure render imports", repr(error))
        return False

    python_ok = sys.version_info[:3] == EXPECTED_FIGURE_PYTHON_VERSION
    audit.check(
        python_ok,
        "Figure Python version",
        ".".join(map(str, sys.version_info[:3])),
    )

    observed_versions = {
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
        "matplotlib": matplotlib.__version__,
        "seaborn": seaborn.__version__,
        "Pillow": PIL.__version__,
    }
    packages_ok = observed_versions == EXPECTED_FIGURE_PACKAGE_VERSIONS
    audit.check(packages_ok, "Figure package versions", str(observed_versions))

    freetype_version = ft2font.__freetype_version__
    freetype_build = getattr(ft2font, "__freetype_build_type__", "unknown")
    freetype_ok = freetype_version == EXPECTED_FREETYPE_VERSION and freetype_build == "system"
    audit.check(
        freetype_ok,
        "Figure FreeType build",
        f"{freetype_version}; {freetype_build}",
    )

    platform_ok = sys.platform == "win32"
    audit.check(platform_ok, "Byte-identical figure platform", sys.platform)

    font_specs = {
        "regular": font_manager.FontProperties(family="Arial", style="normal", weight="normal"),
        "bold": font_manager.FontProperties(family="Arial", style="normal", weight="bold"),
        "italic": font_manager.FontProperties(family="Arial", style="italic", weight="normal"),
    }
    observed_fonts = {}
    try:
        for label, properties in font_specs.items():
            font_path = Path(font_manager.findfont(properties, fallback_to_default=False))
            observed_fonts[label] = {
                "file": font_path.name,
                "sha256": sha256(font_path),
            }
    except Exception as error:
        observed_fonts["error"] = repr(error)
    fonts_ok = all(
        observed_fonts.get(label, {}).get("sha256") == expected
        for label, expected in EXPECTED_ARIAL_SHA256.items()
    )
    audit.check(fonts_ok, "Figure Arial font fingerprints", str(observed_fonts))

    return python_ok and packages_ok and freetype_ok and platform_ok and fonts_ok


def run_figures(audit):
    if not check_figure_environment(audit):
        return
    source_hashes_before = {path: sha256(path) for path in SOURCE.iterdir() if path.is_file()}
    scripts = [
        ROOT / "scripts" / "python_figure_assembly" / "rebuild_scientific_data_figure2_qc_final.py",
        ROOT / "scripts" / "python_figure_assembly" / "rebuild_scientific_data_figure3_global_annotation.py",
        ROOT / "scripts" / "python_figure_assembly" / "rebuild_scientific_data_figure4_immune_reclustering.py",
        ROOT / "scripts" / "python_figure_assembly" / "rebuild_scientific_data_figure5_composition.py",
        ROOT / "scripts" / "python_figure_assembly" / "rebuild_scientific_data_supplementary_figures.py",
    ]
    all_scripts_passed = True
    for script in scripts:
        result = subprocess.run([sys.executable, "-B", str(script)], cwd=ROOT)
        audit.check(result.returncode == 0, f"Run {script.name}", f"exit {result.returncode}")
        all_scripts_passed = all_scripts_passed and result.returncode == 0

    if not all_scripts_passed:
        return

    from PIL import Image, ImageStat

    expected_sizes = {
        "Figure2_QC_final.png": (3205, 2388),
        "Figure3_global_annotation_final.png": (7080, 6120),
        "Figure4_immune_reclustering_final.png": (6900, 4890),
        "Figure5_composition_final.png": (7301, 5419),
    }
    figure_dir = ROOT / "outputs" / "figures" / "scientific_data"
    for name, expected in expected_sizes.items():
        with Image.open(figure_dir / name) as image:
            audit.check(image.size == expected, f"{name} canvas", f"{image.size[0]} x {image.size[1]}")
            variance = ImageStat.Stat(image.convert("L")).var[0]
            audit.check(variance > 1.0, f"{name} nonblank pixels", f"variance {variance:.2f}")
            audit.check(
                sha256(figure_dir / name) == EXPECTED_FIGURE_SHA256[name],
                f"{name} reference hash",
                "byte-identical to final manuscript PNG",
            )

    for name, expected_hash in EXPECTED_FIGURE_PDF_SHA256.items():
        path = figure_dir / name
        content = path.read_bytes() if path.is_file() else b""
        audit.check(
            path.is_file() and sha256(path) == expected_hash,
            f"{name} reference hash",
            "byte-identical deterministic PDF",
        )
        audit.check(
            b"/Type3" not in content and b"/CIDFontType2" in content and b"/FontFile2" in content,
            f"{name} embedded TrueType fonts",
            "CIDFontType2 embedded; no Type3 font",
        )
        audit.check(
            b"/CreationDate" not in content,
            f"{name} deterministic metadata",
            "no CreationDate",
        )

    supplementary_sizes = {
        "Supplementary_Figure_S1_QC.png": (4389, 3189),
        "Supplementary_Figure_S2_UMAP.png": (4389, 2079),
        "Supplementary_Figure_S3_Representation.png": (4238, 2280),
    }
    supplementary_dir = ROOT / "outputs" / "figures" / "supplementary"
    for name, expected in supplementary_sizes.items():
        with Image.open(supplementary_dir / name) as image:
            audit.check(image.size == expected, f"{name} canvas", f"{image.size[0]} x {image.size[1]}")
            variance = ImageStat.Stat(image.convert("L")).var[0]
            audit.check(variance > 1.0, f"{name} nonblank pixels", f"variance {variance:.2f}")
    for name, expected_hash in EXPECTED_SUPPLEMENTARY_SHA256.items():
        path = supplementary_dir / name
        audit.check(
            path.is_file() and sha256(path) == expected_hash,
            f"{name} reference hash",
            "byte-identical deterministic supplementary output",
        )

    source_hashes_after = {path: sha256(path) for path in SOURCE.iterdir() if path.is_file()}
    audit.check(source_hashes_before == source_hashes_after, "Figure assembly source immutability", "all source hashes unchanged")

    zenodo_manifests = [ROOT / "manifest" / "zenodo_files_manifest.csv"] if (ROOT / "manifest" / "zenodo_files_manifest.csv").is_file() else []
    if len(zenodo_manifests) == 1:
        audit_file_manifest(
            audit,
            zenodo_manifests[0],
            "Zenodo post-generation manifest",
        )


def main():
    parser = argparse.ArgumentParser(description="Validate the lightweight reproducibility release.")
    parser.add_argument("--check-r", action="store_true", help="Parse retained R scripts with Rscript.")
    parser.add_argument(
        "--run-figures",
        action="store_true",
        help="Run Figure 2-5 and Supplementary Figure S1-S3 assembly and verify outputs.",
    )
    parser.add_argument("--report", type=Path, help="Optional CSV path for the validation results.")
    args = parser.parse_args()

    audit = Audit()
    audit_release(audit)
    audit_supplementary_sources(audit)
    if args.check_r:
        check_r_syntax(audit)
    if args.run_figures:
        run_figures(audit)

    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with args.report.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["status", "check", "detail"])
            writer.writeheader()
            writer.writerows(audit.records)
        print(f"Report: {args.report}")

    if audit.failures:
        print(f"\nRelease validation failed: {len(audit.failures)} check(s).")
        return 1
    print("\nRelease validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
