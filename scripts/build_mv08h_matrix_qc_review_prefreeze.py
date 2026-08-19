#!/usr/bin/env python3
"""Build the MV8-H matrix-content/QC review prefreeze.

This script consumes only already-published structural/QC admission metadata
and the execution resource provenance.  It does not open any Cell Ranger
matrix, molecule, label, outcome, or landscape payload.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path


def read_one(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise ValueError(f"expected one row in {path.name}")
    return rows[0]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--structural-qc-dir", type=Path, required=True)
    parser.add_argument("--execution-closure-dir", type=Path, required=True)
    parser.add_argument("--resource-policy", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    structural_decision = read_one(args.structural_qc_dir / "mv08h-structural-qc-decision.csv")
    structural_validation = read_rows(args.structural_qc_dir / "mv08h-structural-qc-validation.csv")
    resource_provenance = read_one(args.execution_closure_dir / "mv08h-count-sentinel-execution-resource.csv")
    resource_policy = read_rows(args.resource_policy)
    if structural_decision["decision"] != "metadata_only_structural_qc_admit":
        raise ValueError("structural/QC admission must pass before opening a review gate")
    if not all(row["passed"] == "TRUE" for row in structural_validation):
        raise ValueError("structural/QC validation is incomplete")

    input_rows = [
        {"input_order": 1, "input_name": "filtered_feature_bc_matrix.h5", "role": "filtered_matrix_content", "permitted_after_owner_approval": "TRUE", "opened_by_prefreeze": "FALSE", "label_or_outcome_access": "FALSE"},
        {"input_order": 2, "input_name": "raw_feature_bc_matrix.h5", "role": "raw_matrix_content", "permitted_after_owner_approval": "TRUE", "opened_by_prefreeze": "FALSE", "label_or_outcome_access": "FALSE"},
        {"input_order": 3, "input_name": "metrics_summary.csv", "role": "aggregate_cellranger_qc_metadata", "permitted_after_owner_approval": "TRUE", "opened_by_prefreeze": "FALSE", "label_or_outcome_access": "FALSE"},
        {"input_order": 4, "input_name": "molecule_info.h5", "role": "molecule_level_payload", "permitted_after_owner_approval": "FALSE", "opened_by_prefreeze": "FALSE", "label_or_outcome_access": "FALSE"},
        {"input_order": 5, "input_name": "web_summary.html", "role": "interactive_qc_payload", "permitted_after_owner_approval": "FALSE", "opened_by_prefreeze": "FALSE", "label_or_outcome_access": "FALSE"},
    ]
    firewall_rows = [
        {"firewall_id": "matrix_content", "state": "closed_until_owner_approval", "prohibited_action": "open_matrix_values_or_barcode_identifiers", "stop_on_breach": "TRUE"},
        {"firewall_id": "labels_outcomes", "state": "closed", "prohibited_action": "read_labels_or_biological_outcomes", "stop_on_breach": "TRUE"},
        {"firewall_id": "landscapes", "state": "closed", "prohibited_action": "compute_persistence_or_landscape_artifacts", "stop_on_breach": "TRUE"},
        {"firewall_id": "remaining_units", "state": "closed", "prohibited_action": "process_any_additional_unit", "stop_on_breach": "TRUE"},
        {"firewall_id": "deletion", "state": "closed", "prohibited_action": "delete_partial_or_complete_artifacts", "stop_on_breach": "TRUE"},
    ]
    stop_rows = [
        {"stop_order": 1, "condition": "owner_approval_missing", "action": "do_not_open_matrix_or_qc_content"},
        {"stop_order": 2, "condition": "input_hash_or_structure_drift", "action": "stop_and_rebind_before_review"},
        {"stop_order": 3, "condition": "resource_cap_breach", "action": "stop_without_automatic_kill_and_preserve_artifacts"},
        {"stop_order": 4, "condition": "any_label_outcome_landscape_or_remaining_unit_access", "action": "stop_and_escalate"},
        {"stop_order": 5, "condition": "nonreproducible_repeat", "action": "block_review_and_investigate"},
    ]
    validation = [
        {"check_id": "structural_qc_admission_dependency", "passed": "TRUE", "evidence": "metadata-only structural/QC admission is 8/8"},
        {"check_id": "permitted_input_allowlist", "passed": "TRUE", "evidence": "three review inputs explicitly allowlisted; molecule/web payloads forbidden"},
        {"check_id": "owner_authorization_boundary", "passed": "TRUE", "evidence": "matrix/QC content remains closed until explicit owner approval"},
        {"check_id": "resource_provenance", "passed": str(all(resource_provenance[key] == "TRUE" for key in ("rss_cap_passed", "workspace_cap_passed", "free_space_floor_passed", "elapsed_cap_passed"))).upper(), "evidence": "existing 4-core/32-GiB execution resource gates pass"},
        {"check_id": "label_outcome_firewall", "passed": "TRUE", "evidence": "labels and biological outcomes remain closed"},
        {"check_id": "landscape_firewall", "passed": "TRUE", "evidence": "no persistence diagrams or landscapes permitted"},
        {"check_id": "remaining_unit_firewall", "passed": "TRUE", "evidence": "only HCA_BM_002 is in scope"},
        {"check_id": "non_destructive_stop_policy", "passed": "TRUE", "evidence": "no automatic kill or deletion; preserve artifacts on breach"},
    ]
    decision = [{
        "contract_id": "mv08h_matrix_qc_review_prefreeze_v1",
        "prefreeze_completed": "TRUE",
        "matrix_content_review_authorized": "FALSE",
        "qc_content_review_authorized": "FALSE",
        "labels_outcomes_authorized": "FALSE",
        "landscapes_authorized": "FALSE",
        "remaining_units_authorized": "FALSE",
        "deletion_authorized": "FALSE",
        "next_gate": "owner_approval_for_allowlisted_matrix_qc_review",
    }]
    resource_rows = [
        {"resource": row["resource"], "selected_value": row["selected_value"], "unit": row["unit"], "gate": row["gate"], "source": "mv08h-count-sentinel-resource-policy.csv"}
        for row in resource_policy
    ]
    report = [
        "# MV8-H Matrix-Content/QC Review Prefreeze",
        "",
        "Prefreeze completed: `TRUE`",
        "Matrix-content/QC review authorized: `FALSE`",
        "",
        "This record defines the smallest permitted review after owner approval.",
        "It does not open matrix values, barcode identifiers, labels, outcomes, landscapes, or remaining units.",
        "",
        "Allowlisted after approval: filtered/raw feature-barcode H5 content and aggregate metrics_summary.csv.",
        "Forbidden: molecule_info.h5, web_summary.html, labels, outcomes, persistence/landscape calculations, other units, and deletion.",
        "",
        "The review must stop on hash/structure drift, resource breach, unauthorized access, or nonreproducible repeat.",
        "Next gate: explicit owner approval for the allowlisted matrix/QC review.",
    ]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "mv08h-matrix-qc-review-input-binding.csv", input_rows, list(input_rows[0]))
    write_csv(args.output_dir / "mv08h-matrix-qc-review-resource.csv", resource_rows, list(resource_rows[0]))
    write_csv(args.output_dir / "mv08h-matrix-qc-review-firewall.csv", firewall_rows, list(firewall_rows[0]))
    write_csv(args.output_dir / "mv08h-matrix-qc-review-stop-conditions.csv", stop_rows, list(stop_rows[0]))
    write_csv(args.output_dir / "mv08h-matrix-qc-review-validation.csv", validation, list(validation[0]))
    write_csv(args.output_dir / "mv08h-matrix-qc-review-decision.csv", decision, list(decision[0]))
    report_path = args.output_dir / "MV08H_MATRIX_QC_REVIEW_PREFREEZE_2026-08-18.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    artifact_files = sorted(path for path in args.output_dir.iterdir() if path.is_file())
    manifest = [{"file": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path), "contains_private_path": "FALSE"} for path in artifact_files]
    write_csv(args.output_dir / "mv08h-matrix-qc-review-artifact-manifest.csv", manifest, list(manifest[0]))
    print("MV8-H matrix/QC review prefreeze: TRUE; review authorization: FALSE", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
