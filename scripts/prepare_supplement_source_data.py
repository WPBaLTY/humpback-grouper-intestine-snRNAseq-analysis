from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd


SAMPLE_ORDER = ["CTL1", "CTL2", "Llac1", "Llac2", "Slim1", "Slim2", "COM1", "COM2"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def public_cell_id(value: str) -> str:
    if value.startswith("CON1_"):
        return "COM1_" + value[len("CON1_") :]
    if value.startswith("CON2_"):
        return "COM2_" + value[len("CON2_") :]
    return value


def public_sample(value: str) -> str:
    return {"CON1": "COM1", "CON2": "COM2"}.get(value, value)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare Scientific Data supplementary figure sources.")
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--final-umap", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    inventory = pd.read_csv(args.inventory)
    first_column = inventory.columns[0]
    inventory = inventory.rename(columns={first_column: "historical_cell_id"})
    required_inventory = {
        "historical_cell_id",
        "nCount_RNA",
        "nFeature_RNA",
        "sample",
        "percent.mt.true",
    }
    missing = required_inventory - set(inventory.columns)
    if missing:
        raise RuntimeError(f"QC inventory is missing columns: {sorted(missing)}")

    inventory["cell_id"] = inventory["historical_cell_id"].astype(str).map(public_cell_id)
    inventory["sample"] = inventory["sample"].astype(str).map(public_sample)
    inventory = inventory.rename(columns={"percent.mt.true": "percent_mt"})

    final_umap = pd.read_csv(args.final_umap, compression="gzip")
    required_final = {"cell_id", "sample", "sample_plot", "group"}
    missing = required_final - set(final_umap.columns)
    if missing:
        raise RuntimeError(f"Final UMAP table is missing columns: {sorted(missing)}")

    if len(final_umap) != 102_036 or final_umap["cell_id"].nunique() != 102_036:
        raise RuntimeError("Final UMAP table must contain 102,036 unique nuclei")
    if inventory["cell_id"].nunique() != len(inventory):
        raise RuntimeError("QC inventory contains duplicate public cell IDs")

    qc = final_umap[["cell_id", "sample", "sample_plot", "group"]].merge(
        inventory[["cell_id", "nCount_RNA", "nFeature_RNA", "percent_mt"]],
        on="cell_id",
        how="left",
        validate="one_to_one",
    )
    if qc[["nCount_RNA", "nFeature_RNA", "percent_mt"]].isna().any().any():
        raise RuntimeError("Final cell set does not close to the QC inventory")
    if list(qc.groupby("sample", sort=False).size().reindex(SAMPLE_ORDER).index) != SAMPLE_ORDER:
        raise RuntimeError("Unexpected sample order after QC join")

    expected_counts = {
        "CTL1": 12_647,
        "CTL2": 13_767,
        "Llac1": 13_286,
        "Llac2": 11_803,
        "Slim1": 13_140,
        "Slim2": 10_545,
        "COM1": 13_024,
        "COM2": 13_824,
    }
    observed_counts = qc.groupby("sample").size().astype(int).to_dict()
    if observed_counts != expected_counts:
        raise RuntimeError(f"Per-library cell counts differ: {observed_counts}")

    qc_path = args.output_dir / "SuppFigureS1_qc_per_cell.csv.gz"
    qc.to_csv(
        qc_path,
        index=False,
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
    )

    summary = (
        qc.groupby(["sample", "sample_plot", "group"], sort=False)
        .agg(
            retained_nuclei=("cell_id", "size"),
            median_detected_genes=("nFeature_RNA", "median"),
            median_umis=("nCount_RNA", "median"),
            median_percent_mt=("percent_mt", "median"),
        )
        .reset_index()
    )
    summary_path = args.output_dir / "SuppFigureS1_qc_summary.csv"
    summary.to_csv(summary_path, index=False)

    manifest = pd.DataFrame(
        [
            {
                "role": "input_qc_inventory",
                "file": args.inventory.name,
                "rows": len(inventory),
                "sha256": sha256(args.inventory),
            },
            {
                "role": "input_final_umap",
                "file": args.final_umap.name,
                "rows": len(final_umap),
                "sha256": sha256(args.final_umap),
            },
            {
                "role": "output_qc_per_cell",
                "file": qc_path.name,
                "rows": len(qc),
                "sha256": sha256(qc_path),
            },
            {
                "role": "output_qc_summary",
                "file": summary_path.name,
                "rows": len(summary),
                "sha256": sha256(summary_path),
            },
        ]
    )
    manifest.to_csv(args.output_dir / "supplement_source_manifest_stage1.csv", index=False)

    print(f"QC rows: {len(qc):,}")
    print(f"QC per-library counts: {observed_counts}")
    print(f"Output: {qc_path}")


if __name__ == "__main__":
    main()
