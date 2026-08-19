#!/usr/bin/env python3
"""Independently validate exact-500 feasibility from an existing Cell Ranger sentinel.

This reads only feature/barcode axes, the already-approved aggregate QC summary,
and execution metadata. It does not publish expression values, barcodes, labels,
outcomes, or matrix contents.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import re
import sys
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_gzip_first_column(path: Path) -> list[str]:
    values: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            values.append(line.rstrip("\n").split("\t", 1)[0])
    return values


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def add(checks: list[dict[str, object]], check_id: str, passed: bool, evidence: str) -> None:
    checks.append({"check_id": check_id, "passed": bool(passed), "evidence": evidence})


def main(argv: list[str]) -> int:
    if len(argv) != 6:
        raise SystemExit(
            "usage: validate_mv08h_exact500_sentinel.py <count-root> <exact500-panel.csv> "
            "<common475-panel.csv> <matrix-qc-summary.csv> <execution-audit-dir> <output-dir>"
        )
    count_root = Path(argv[0]).resolve()
    exact_panel_path = Path(argv[1]).resolve()
    common_panel_path = Path(argv[2]).resolve()
    qc_summary_path = Path(argv[3]).resolve()
    execution_dir = Path(argv[4]).resolve()
    output_dir = Path(argv[5]).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    outs = count_root / "outs"
    filtered = outs / "filtered_feature_bc_matrix"
    required = [filtered / "features.tsv.gz", filtered / "barcodes.tsv.gz",
                outs / "filtered_feature_bc_matrix.h5", outs / "raw_feature_bc_matrix.h5",
                outs / "molecule_info.h5", outs / "metrics_summary.csv"]
    if not all(path.is_file() for path in required):
        missing = [str(path.name) for path in required if not path.is_file()]
        raise SystemExit(f"missing sentinel outputs: {missing}")

    feature_ids = [re.sub(r"\\.[0-9]+$", "", value) for value in read_gzip_first_column(filtered / "features.tsv.gz")]
    barcode_count = sum(1 for _ in gzip.open(filtered / "barcodes.tsv.gz", "rt", encoding="utf-8"))
    exact_rows = read_rows(exact_panel_path)
    common_rows = read_rows(common_panel_path)
    exact_ids = [re.search(r"(ENSG[0-9]+)", row["feature_id"]).group(1) for row in exact_rows]
    common_ids = [row["ensembl_stable_id"] for row in common_rows]
    feature_set = set(feature_ids)

    qc_rows = read_rows(qc_summary_path)
    if len(qc_rows) != 1:
        raise SystemExit("matrix QC summary must contain exactly one sentinel row")
    qc = qc_rows[0]
    resource_rows = read_rows(execution_dir / "mv08h-count-sentinel-execution-resource.csv")
    decision_rows = read_rows(execution_dir / "mv08h-count-sentinel-execution-decision.csv")
    if len(resource_rows) != 1 or len(decision_rows) != 1:
        raise SystemExit("execution audit must contain one resource and one decision row")
    resource = resource_rows[0]
    decision = decision_rows[0]

    checks: list[dict[str, object]] = []
    add(checks, "feature_axis_count", len(feature_ids) == 33563, f"unique feature rows={len(feature_ids)}")
    add(checks, "feature_axis_unique", len(feature_ids) == len(set(feature_ids)), "feature IDs unique after version stripping")
    add(checks, "exact500_panel_present", len(exact_ids) == 500 and len(set(exact_ids)) == 500 and set(exact_ids) <= feature_set,
        f"exact500 present={sum(value in feature_set for value in exact_ids)}/500")
    add(checks, "common475_panel_present", len(common_ids) == 475 and len(set(common_ids)) == 475 and set(common_ids) <= feature_set,
        f"common475 present={sum(value in feature_set for value in common_ids)}/475")
    add(checks, "common475_subset_exact500", set(common_ids) <= set(exact_ids), "common475 is an exact ordered subset of exact500")
    add(checks, "filtered_cell_axis", barcode_count == 5037, f"filtered barcode rows={barcode_count}")
    post_qc = int(float(qc["post_qc_cells"]))
    add(checks, "qc_384_eligibility", post_qc >= 384 and qc["depth_384_pass"].upper() == "TRUE", f"post-QC eligible cells={post_qc}")
    add(checks, "execution_success", decision["execution_performed"].upper() == "TRUE" and decision["run_success"].upper() == "TRUE",
        "existing Cell Ranger sentinel completed successfully")
    add(checks, "resource_closure", resource["rss_cap_passed"].upper() == "TRUE" and resource["workspace_cap_passed"].upper() == "TRUE" and resource["free_space_floor_passed"].upper() == "TRUE",
        f"elapsed={resource['elapsed_seconds']} s; peak RSS={resource['peak_rss_bytes']} B; cores={resource['local_cores']}; memory={resource['local_memory_gib']} GiB")
    add(checks, "downstream_firewall", all(decision[key].upper() == "FALSE" for key in ("matrix_access_authorized", "qc_authorized", "pca_ph_landscape_authorized", "label_access_authorized", "biological_outcomes_authorized", "remaining_units_authorized", "deletion_authorized")),
        "labels, outcomes, landscapes, remaining units, and deletion remain closed")

    artifact_rows = []
    for path in (outs / "filtered_feature_bc_matrix.h5", outs / "raw_feature_bc_matrix.h5", outs / "molecule_info.h5", outs / "metrics_summary.csv"):
        artifact_rows.append({"artifact": path.name, "bytes": path.stat().st_size, "sha256": sha256(path), "private_content": "TRUE"})
    with (output_dir / "mv08h-exact500-sentinel-artifact-identity.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["artifact", "bytes", "sha256", "private_content"])
        writer.writeheader(); writer.writerows(artifact_rows)
    with (output_dir / "mv08h-exact500-sentinel-validation.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check_id", "passed", "evidence"])
        writer.writeheader(); writer.writerows(checks)
    report = [
        "# MV8-H exact-500 sentinel feasibility",
        "",
        "This audit reopens only the existing HCA_BM_002 Cell Ranger 8.0.1 feature/barcode axes, approved aggregate QC summary, and execution metadata.",
        "",
        f"- Feature axis: {len(feature_ids):,} unique features; exact500 present {sum(value in feature_set for value in exact_ids)}/500; common475 present {sum(value in feature_set for value in common_ids)}/475.",
        f"- Filtered cell axis: {barcode_count:,} cells; existing frozen QC review reports {post_qc:,} cells passing the 384-cell minimum.",
        f"- Existing run: {resource['elapsed_seconds']} seconds, peak RSS {resource['peak_rss_bytes']} bytes, {resource['local_cores']} cores, {resource['local_memory_gib']} GiB.",
        "- Exact500 is technically feasible for this sentinel; this does not authorize remaining units, PH, landscapes, labels, outcomes, fusion, manuscript work, or deletion.",
        "",
        f"Independent checks: {sum(bool(row['passed']) for row in checks)}/{len(checks)} passed.",
    ]
    (output_dir / "MV08H_EXACT500_SENTINEL_FEASIBILITY_2026-08-19.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    if not all(row["passed"] for row in checks):
        return 1
    print(f"MV8-H exact500 sentinel feasibility passed {len(checks)}/{len(checks)} checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
