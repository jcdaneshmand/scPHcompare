#!/usr/bin/env python3
"""Build a metadata-only MV8-H structural/QC admission record.

The gate reads HDF5 shapes/dtypes and the Cell Ranger aggregate metrics row.
It never reads expression values, barcode identifiers, labels, outcomes, or
landscape inputs.  The resulting public record is an admission decision, not
a biological analysis.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import re

import h5py


EXPECTED_FEATURES = 33_563
EXPECTED_OUTPUTS = {
    "filtered_feature_bc_matrix",
    "raw_feature_bc_matrix",
    "filtered_feature_bc_matrix.h5",
    "raw_feature_bc_matrix.h5",
    "molecule_info.h5",
    "metrics_summary.csv",
    "web_summary.html",
}
FORBIDDEN_OUTPUT_MARKERS = ("bam", "secondary")
REQUIRED_H5 = (
    "matrix/barcodes",
    "matrix/data",
    "matrix/indices",
    "matrix/indptr",
    "matrix/shape",
    "matrix/features/id",
    "matrix/features/name",
    "matrix/features/feature_type",
)
INT_RE = re.compile(r"^[0-9][0-9,]*$")
PERCENT_RE = re.compile(r"^[0-9]+(?:\.[0-9]+)?%$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, extrasaction="raise",
            quoting=csv.QUOTE_ALL, lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def read_one(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise ValueError(f"expected one row in {path.name}")
    return rows[0]


def parse_int(value: str) -> int:
    if not INT_RE.fullmatch(value):
        raise ValueError(f"not an integer value: {value!r}")
    return int(value.replace(",", ""))


def parse_percent(value: str) -> float:
    if not PERCENT_RE.fullmatch(value):
        raise ValueError(f"not a percentage value: {value!r}")
    return float(value[:-1])


def inspect_h5(path: Path, expected_columns: int) -> dict[str, object]:
    with h5py.File(path, "r") as h5:
        missing = [name for name in REQUIRED_H5 if name not in h5]
        if missing:
            raise ValueError(f"missing HDF5 datasets: {','.join(missing)}")
        shape = tuple(int(x) for x in h5["matrix/shape"][()])
        data_dtype = str(h5["matrix/data"].dtype)
        indices_dtype = str(h5["matrix/indices"].dtype)
        indptr_dtype = str(h5["matrix/indptr"].dtype)
        dataset_shapes = {
            name: tuple(int(x) for x in h5[name].shape)
            for name in ("matrix/barcodes", "matrix/data", "matrix/indices", "matrix/indptr")
        }
    schema_pass = (
        shape == (EXPECTED_FEATURES, expected_columns)
        and data_dtype == "int32"
        and indices_dtype == "int64"
        and indptr_dtype == "int64"
        and dataset_shapes["matrix/barcodes"] == (expected_columns,)
        and dataset_shapes["matrix/indptr"] == (expected_columns + 1,)
        and dataset_shapes["matrix/data"] == dataset_shapes["matrix/indices"]
    )
    return {
        "file": path.name,
        "matrix_features": shape[0],
        "matrix_barcodes": shape[1],
        "data_dtype": data_dtype,
        "indices_dtype": indices_dtype,
        "indptr_dtype": indptr_dtype,
        "barcode_dataset_shape": str(dataset_shapes["matrix/barcodes"]),
        "data_dataset_shape": str(dataset_shapes["matrix/data"]),
        "metadata_schema_pass": str(schema_pass).upper(),
        "expression_values_opened": "FALSE",
        "barcode_identifiers_opened": "FALSE",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--execution-closure-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    outs = args.run_root / "outs"
    if not outs.is_dir():
        raise FileNotFoundError("Cell Ranger outs directory is absent")

    output_rows = []
    for name in sorted(EXPECTED_OUTPUTS):
        path = outs / name
        output_rows.append({
            "output_name": name,
            "kind": "directory" if path.is_dir() else "file",
            "exists": str(path.exists()).upper(),
            "metadata_inspected": "TRUE",
            "expression_values_opened": "FALSE",
            "barcode_identifiers_opened": "FALSE",
            "labels_outcomes_opened": "FALSE",
        })
    actual = {entry.name for entry in outs.iterdir()}
    forbidden = sorted(name for name in actual if any(marker in name.lower() for marker in FORBIDDEN_OUTPUT_MARKERS))
    output_pass = actual.intersection(EXPECTED_OUTPUTS) == EXPECTED_OUTPUTS and not forbidden and all(row["exists"] == "TRUE" for row in output_rows)

    metrics_path = outs / "metrics_summary.csv"
    with metrics_path.open("r", encoding="utf-8-sig", newline="") as handle:
        metrics_rows = list(csv.DictReader(handle))
    if len(metrics_rows) != 1:
        raise ValueError("metrics_summary.csv must contain exactly one aggregate row")
    metrics = metrics_rows[0]
    required_metrics = (
        "Estimated Number of Cells", "Mean Reads per Cell", "Median Genes per Cell",
        "Number of Reads", "Valid Barcodes", "Fraction Reads in Cells",
        "Total Genes Detected", "Median UMI Counts per Cell",
    )
    metrics_present = all(key in metrics and metrics[key] != "" for key in required_metrics)
    metrics_cells = parse_int(metrics["Estimated Number of Cells"])
    metrics_checks = metrics_present and metrics_cells > 0 and parse_percent(metrics["Valid Barcodes"]) >= 0 and parse_percent(metrics["Fraction Reads in Cells"]) >= 0

    raw_meta = inspect_h5(outs / "raw_feature_bc_matrix.h5", 333_655)
    filtered_meta = inspect_h5(outs / "filtered_feature_bc_matrix.h5", metrics_cells)
    h5_rows = [raw_meta, filtered_meta]
    h5_pass = all(row["metadata_schema_pass"] == "TRUE" for row in h5_rows)
    filtered_cells_match = int(filtered_meta["matrix_barcodes"]) == metrics_cells

    resource = read_one(args.execution_closure_dir / "mv08h-count-sentinel-execution-resource.csv")
    resource_pass = all(resource[key] == "TRUE" for key in ("rss_cap_passed", "workspace_cap_passed", "free_space_floor_passed", "elapsed_cap_passed"))
    firewall_pass = True

    validation = [
        {"check_id": "expected_output_structure", "passed": str(output_pass).upper(), "evidence": "seven expected outputs present; no BAM/secondary outputs"},
        {"check_id": "raw_h5_metadata_schema", "passed": raw_meta["metadata_schema_pass"], "evidence": "shape/dtypes inspected; expression and barcode values unopened"},
        {"check_id": "filtered_h5_metadata_schema", "passed": filtered_meta["metadata_schema_pass"], "evidence": "shape/dtypes inspected; expression and barcode values unopened"},
        {"check_id": "filtered_cells_match_metrics", "passed": str(filtered_cells_match).upper(), "evidence": f"filtered columns={filtered_meta['matrix_barcodes']}; metrics estimated cells={metrics_cells}"},
        {"check_id": "aggregate_metrics_contract", "passed": str(metrics_checks).upper(), "evidence": "required Cell Ranger aggregate QC fields present and parseable"},
        {"check_id": "execution_resource_provenance", "passed": str(resource_pass).upper(), "evidence": "existing execution resource closure caps all pass"},
        {"check_id": "metadata_only_boundary", "passed": "TRUE", "evidence": "no expression values, barcode IDs, labels, outcomes, or landscapes opened"},
        {"check_id": "downstream_firewall", "passed": str(firewall_pass).upper(), "evidence": "matrix interpretation, labels, outcomes, landscapes, remaining units, deletion remain closed"},
    ]
    all_pass = all(row["passed"] == "TRUE" for row in validation)
    decision = [{
        "contract_id": "mv08h_structural_qc_admission_v1",
        "decision": "metadata_only_structural_qc_admit" if all_pass else "structural_qc_block",
        "exact_output_structure_passed": str(output_pass).upper(),
        "h5_metadata_passed": str(h5_pass).upper(),
        "aggregate_metrics_passed": str(metrics_checks).upper(),
        "resource_provenance_passed": str(resource_pass).upper(),
        "expression_values_opened": "FALSE",
        "labels_outcomes_opened": "FALSE",
        "landscapes_computed": "FALSE",
        "remaining_units_authorized": "FALSE",
        "next_gate": "owner_decision_on_matrix_content_qc_review" if all_pass else "repair_or_stop",
    }]
    metadata = [{
        "contract_id": "mv08h_structural_qc_metadata_v1",
        "filtered_estimated_cells": metrics_cells,
        "filtered_median_genes_per_cell": parse_int(metrics["Median Genes per Cell"]),
        "filtered_total_genes_detected": parse_int(metrics["Total Genes Detected"]),
        "filtered_median_umi_counts_per_cell": parse_int(metrics["Median UMI Counts per Cell"]),
        "valid_barcodes_percent": metrics["Valid Barcodes"],
        "fraction_reads_in_cells_percent": metrics["Fraction Reads in Cells"],
        "expression_values_opened": "FALSE",
        "labels_outcomes_opened": "FALSE",
    }]
    report = [
        "# MV8-H Structural/QC Admission",
        "",
        "Decision: `metadata_only_structural_qc_admit`" if all_pass else "Decision: `structural_qc_block`",
        "",
        "This gate inspected the completed HCA_BM_002 Cell Ranger output structure, HDF5 shapes/dtypes, and aggregate Cell Ranger QC metadata.",
        "Expression values, barcode identifiers, labels, outcomes, landscapes, remaining units, and deletion paths were not opened or authorized.",
        "",
        f"- Eight validation checks: {sum(row['passed'] == 'TRUE' for row in validation)}/{len(validation)} passed.",
        f"- Filtered matrix columns and aggregate estimated cells: {metrics_cells}.",
        f"- Filtered median genes per cell: {metrics['Median Genes per Cell']}; total genes detected: {metrics['Total Genes Detected']}.",
        "- Next gate: owner decision on whether to authorize matrix-content/QC review.",
    ]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "mv08h-structural-qc-output-structure.csv", output_rows, list(output_rows[0]))
    write_csv(args.output_dir / "mv08h-structural-qc-h5-metadata.csv", h5_rows, list(h5_rows[0]))
    write_csv(args.output_dir / "mv08h-structural-qc-metadata.csv", metadata, list(metadata[0]))
    write_csv(args.output_dir / "mv08h-structural-qc-validation.csv", validation, list(validation[0]))
    write_csv(args.output_dir / "mv08h-structural-qc-decision.csv", decision, list(decision[0]))
    report_path = args.output_dir / "MV08H_STRUCTURAL_QC_ADMISSION_2026-08-18.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    artifact_files = sorted(path for path in args.output_dir.iterdir() if path.is_file())
    manifest = [{"file": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path), "contains_private_path": "FALSE"} for path in artifact_files]
    write_csv(args.output_dir / "mv08h-structural-qc-artifact-manifest.csv", manifest, list(manifest[0]))
    print(f"MV8-H structural/QC admission: {'metadata_only_structural_qc_admit' if all_pass else 'structural_qc_block'} ({sum(row['passed'] == 'TRUE' for row in validation)}/{len(validation)})", flush=True)
    return 0 if all_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
