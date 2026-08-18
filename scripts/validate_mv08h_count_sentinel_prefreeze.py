#!/usr/bin/env python3
"""Independently validate the MV8-H count-sentinel prefreeze."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import re


EXPECTED_REFERENCE_SHA256 = (
    "5e2aff9e7154e6b02f98552a4419bd48edce66e617e579ae562e714f79199f1c"
)
EXPECTED_RUNTIME_SHA256 = (
    "aafd39e293e0ba9d14dba3896a6aeda077304531a2702d26bda0c62c4688fdf3"
)
EXPECTED_LAUNCHER_SHA256 = (
    "4ee3a1670b4f14c826004fe8e17b4759e1edc701b15ff2e9623753bf1b34d4d6"
)
UNIT = "HCA_BM_002"
UNIT_BYTES = 11_249_623_632
RUN_ID = "mv08h_count_sentinel_hca_bm_002"
PRIVATE_PATH = re.compile(r"(^|[^A-Za-z])(?:[A-Za-z]:[\\/]|/mnt/|/home/)")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError("refusing empty validation")
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), extrasaction="raise",
            quoting=csv.QUOTE_ALL, lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows({
            key: "TRUE" if value is True else "FALSE" if value is False else value
            for key, value in row.items()
        } for row in rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def tree_identity(root: Path) -> tuple[int, int, int, str]:
    entries: list[tuple[str, str, Path]] = []
    for directory, dirnames, filenames in os.walk(root, followlinks=False):
        base = Path(directory)
        retained: list[str] = []
        for name in dirnames:
            path = base / name
            if path.is_symlink():
                entries.append(("L", path.relative_to(root).as_posix(), path))
            else:
                retained.append(name)
        dirnames[:] = retained
        for name in filenames:
            path = base / name
            if path.is_symlink():
                entries.append(("L", path.relative_to(root).as_posix(), path))
            elif path.is_file():
                entries.append(("F", path.relative_to(root).as_posix(), path))
    entries.sort(key=lambda item: (item[1], item[0]))
    digest = hashlib.sha256()
    regular = symlinks = regular_bytes = 0
    for kind, relative, path in entries:
        if kind == "F":
            size = path.stat().st_size
            file_sha = sha256_file(path)
            regular += 1
            regular_bytes += size
            record = f"F\0{relative}\0{size}\0{file_sha}\n"
        else:
            symlinks += 1
            record = f"L\0{relative}\0{os.readlink(path)}\n"
        digest.update(record.encode("utf-8"))
    return regular, symlinks, regular_bytes, digest.hexdigest()


def artifact_manifest(output: Path, names: list[str]) -> list[dict[str, object]]:
    return [{
        "contract_id": "mv08h_count_sentinel_artifact_manifest_v1",
        "artifact_order": order,
        "file": name,
        "bytes": (output / name).stat().st_size,
        "sha256": sha256_file(output / name),
        "contains_expression": False,
        "contains_cell_barcode": False,
        "contains_absolute_private_path": False,
        "contains_donor_attribute": False,
        "contains_qc_value": False,
        "contains_outcome_label": False,
    } for order, name in enumerate(names, 1)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--unit-manifest", type=Path, required=True)
    parser.add_argument("--fastq-manifest", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repeat-dir", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    args = parser.parse_args()

    output = args.output_dir.resolve()
    repeat = args.repeat_dir.resolve()
    selection = read_csv(output / "mv08h-count-sentinel-selection.csv")
    fastq_evidence = read_csv(output / "mv08h-count-sentinel-fastqs.csv")
    reference_evidence = read_csv(output / "mv08h-count-sentinel-reference-binding.csv")
    command = read_csv(output / "mv08h-count-sentinel-command.csv")
    resources = read_csv(output / "mv08h-count-sentinel-resource-policy.csv")
    validation_contract = read_csv(output / "mv08h-count-sentinel-validation-contract.csv")
    firewall = read_csv(output / "mv08h-count-sentinel-firewall.csv")
    gates = read_csv(output / "mv08h-count-sentinel-gate-status.csv")
    decision = read_csv(output / "mv08h-count-sentinel-decision.csv")
    units = read_csv(args.unit_manifest)
    fastqs = read_csv(args.fastq_manifest)

    ranked = sorted(
        units, key=lambda row: (-int(row["fastq_bytes"]), int(row["unit_order"]))
    )
    selection_exact = all([
        len(selection) == 8,
        [row["unit_id"] for row in selection] == [row["unit_id"] for row in ranked],
        [int(row["burden_rank"]) for row in selection] == list(range(1, 9)),
        sum(row["selected"] == "TRUE" for row in selection) == 1,
        selection[0]["unit_id"] == UNIT,
        int(selection[0]["fastq_bytes"]) == UNIT_BYTES,
        all(row["label_fields_consulted"] == "FALSE" for row in selection),
        all(row["expression_or_qc_consulted"] == "FALSE" for row in selection),
    ])

    selected_manifest = sorted(
        (row for row in fastqs if row["unit_id"] == UNIT),
        key=lambda row: int(row["file_order"]),
    )
    fastq_exact = len(fastq_evidence) == len(selected_manifest) == 6
    for evidence, manifest in zip(fastq_evidence, selected_manifest):
        path = args.cache_dir / UNIT / manifest["file_name"]
        fastq_exact = fastq_exact and all([
            evidence["file_order"] == manifest["file_order"],
            evidence["file_name"] == manifest["file_name"],
            evidence["file_uuid"] == manifest["file_uuid"],
            evidence["expected_bytes"] == manifest["file_size_bytes"],
            evidence["expected_sha256"] == manifest["file_sha256"],
            evidence["cache_identity_passed"] == "TRUE",
            evidence["private_path_published"] == "FALSE",
            path.is_file(),
            path.stat().st_size == int(manifest["file_size_bytes"]),
            sha256_file(path) == manifest["file_sha256"],
        ])
    fastq_exact = fastq_exact and sum(
        int(row["expected_bytes"]) for row in fastq_evidence
    ) == UNIT_BYTES

    regular, symlinks, reference_bytes, reference_sha = tree_identity(args.reference_dir)
    reference_exact = len(reference_evidence) == 1 and all([
        regular == 19,
        symlinks == 0,
        reference_bytes == 20_765_871_518,
        reference_sha == EXPECTED_REFERENCE_SHA256,
        reference_evidence[0]["tree_sha256"] == reference_sha,
        reference_evidence[0]["runtime_tree_sha256"] == EXPECTED_RUNTIME_SHA256,
        reference_evidence[0]["launcher_sha256"] == EXPECTED_LAUNCHER_SHA256,
        reference_evidence[0]["binding_passed"] == "TRUE",
    ])

    by_parameter = {row["parameter"]: row["frozen_value"] for row in command}
    expected_parameters = {
        "id": RUN_ID, "sample": "MantonBM2_HiSeq_9", "chemistry": "SC3Pv2",
        "expect-cells": "7000", "include-introns": "false",
        "create-bam": "false", "nosecondary": "true", "localcores": "4",
        "localmem": "32", "disable-ui": "true",
    }
    command_exact = len(command) == 13 and all(
        by_parameter.get(key) == value for key, value in expected_parameters.items()
    )
    public_command = by_parameter.get("public_command", "")
    for term in (
        "--id=mv08h_count_sentinel_hca_bm_002",
        "--transcriptome=<verified_custom_reference>",
        "--fastqs=<verified_hca_bm_002_fastq_directory>",
        "--sample=MantonBM2_HiSeq_9", "--chemistry=SC3Pv2",
        "--expect-cells=7000", "--include-introns=false",
        "--create-bam=false", "--nosecondary", "--localcores=4",
        "--localmem=32", "--disable-ui",
    ):
        command_exact = command_exact and term in public_command
    command_exact = command_exact and not PRIVATE_PATH.search(public_command)

    by_resource = {row["resource"]: row for row in resources}
    resource_exact = len(resources) == 9 and all([
        by_resource["local_cores"]["selected_value"] == "4",
        by_resource["local_memory"]["selected_value"] == "32",
        by_resource["process_tree_rss_absolute_cap"]["selected_value"] == str(80 * 1024**3),
        by_resource["run_workspace_cap"]["selected_value"] == str(200 * 1024**3),
        by_resource["minimum_free_space"]["selected_value"] == str(1024**4),
        by_resource["elapsed_observation_cap"]["selected_value"] == str(96 * 60 * 60),
        by_resource["concurrency"]["selected_value"] == "1",
        by_resource["automatic_termination"]["selected_value"] == "false",
        by_resource["automatic_deletion"]["selected_value"] == "false",
    ])

    validation_exact = len(validation_contract) == 12 and all(
        row["count_execution_performed_by_prefreeze"] == "FALSE"
        for row in validation_contract
    )
    validation_text = " ".join(row["frozen_requirement"] for row in validation_contract)
    for term in ("HDF5", "33,563", "exact500", "common475", "barcodes unique and private",
                 "persistence landscape", "remaining-unit decision"):
        validation_exact = validation_exact and term in validation_text

    firewall_exact = len(firewall) == 11
    public_true = [row["field_class"] for row in firewall
                   if row["public_release_permitted"] == "TRUE"]
    firewall_exact = firewall_exact and public_true == [
        "FASTQ byte sizes and SHA-256",
        "reference/runtime relative identities",
        "aggregate resource samples and run status",
        "matrix dimensions and feature-axis identities",
    ]
    forbidden = [
        "expression or UMI values", "cell barcodes", "donor attributes or identifiers",
        "QC values or eligibility decisions", "study/tissue/approach labels",
        "biological outcomes",
    ]
    firewall_exact = firewall_exact and all(
        any(row["field_class"] == field and row["public_release_permitted"] == "FALSE"
            for row in firewall) for field in forbidden
    )

    decision_exact = len(decision) == 1 and all([
        decision[0]["decision"] == "count_sentinel_prefreeze_exact_await_execution_authorization",
        decision[0]["count_sentinel_prefreeze_completed"] == "TRUE",
        decision[0]["count_sentinel_execution_authorized"] == "FALSE",
        decision[0]["matrix_access_authorized"] == "FALSE",
        decision[0]["qc_authorized"] == "FALSE",
        decision[0]["remaining_units_authorized"] == "FALSE",
        decision[0]["pca_ph_landscape_authorized"] == "FALSE",
        decision[0]["label_access_authorized"] == "FALSE",
        decision[0]["biological_outcomes_authorized"] == "FALSE",
        decision[0]["deletion_authorized"] == "FALSE",
    ])
    gates_exact = len(gates) == 6 and all(row["status"] == "pass" for row in gates)

    builder_names = sorted([
        "MV08H_CELLRANGER8_COUNT_SENTINEL_PREFREEZE_2026-08-18.md",
        "mv08h-count-sentinel-command.csv",
        "mv08h-count-sentinel-decision.csv",
        "mv08h-count-sentinel-fastqs.csv",
        "mv08h-count-sentinel-firewall.csv",
        "mv08h-count-sentinel-gate-status.csv",
        "mv08h-count-sentinel-reference-binding.csv",
        "mv08h-count-sentinel-resource-policy.csv",
        "mv08h-count-sentinel-selection.csv",
        "mv08h-count-sentinel-validation-contract.csv",
    ])
    repeat_exact = all(
        (output / name).is_file()
        and (repeat / name).is_file()
        and sha256_file(output / name) == sha256_file(repeat / name)
        for name in builder_names
    )

    report = (output / builder_names[0]).read_text(encoding="utf-8")
    report_exact = all(term in report for term in (
        "does **not** execute or authorize count", "4 cores", "32 GiB",
        "96 hours", "does not kill or delete", "separate typed observation views",
        "H0 and H1 remain separate", "every consecutive active level",
        "no fixed grid", "no universal level cap",
    ))
    runner_text = args.runner.read_text(encoding="utf-8")
    runner_exact = all(term in runner_text for term in (
        "--dry-run", "EXECUTION_TOKEN", "EXPECTED_REFERENCE_TREE_SHA256",
        "FASTQS", "competing_cellranger_processes", "process_tree_rss_kib",
        "RSS_CAP_BYTES", "WORKSPACE_CAP_BYTES",
        "ELAPSED_CAP_SECONDS", "FREE_SPACE_FLOOR_BYTES", "resource_breach_detected",
        "automatic_kill_used", "deletion_used", "--create-bam=false",
        "--nosecondary", "--localcores=4", "--localmem=32",
    ))

    public_files = [output / name for name in builder_names]
    public_firewall = all(
        not PRIVATE_PATH.search(path.read_text(encoding="utf-8", errors="replace"))
        for path in public_files
    )
    count_not_executed = not any(
        child.name == RUN_ID
        for child in (output.parent.parent.parent / "tmp").rglob(RUN_ID)
    )

    checks = [
        ("selection_reconstructed", selection_exact, "unique largest complete unit selected without labels/QC"),
        ("six_fastq_live_identity", fastq_exact, "6/6 files independently rehashed; exact bytes"),
        ("reference_tree_live_identity", reference_exact, "19-file reference independently rehashed"),
        ("command_contract", command_exact, "all Cell Ranger count controls and placeholders exact"),
        ("resource_contract", resource_exact, "4 cores/32 GiB; 80/200-GiB caps; 1-TiB floor; 96 h"),
        ("future_validation_contract", validation_exact, "12 structural/resource/firewall/stop checks frozen"),
        ("firewall_contract", firewall_exact, "no public expression/barcode/donor/QC/label/outcome fields"),
        ("gate_status", gates_exact, "6/6 prefreeze gates pass"),
        ("authorization_boundary", decision_exact, "count and every downstream stage remain unauthorized"),
        ("deterministic_builder_repeat", repeat_exact, "all ten builder artifacts byte-identical"),
        ("human_readable_report", report_exact, "resources, stops, and landscape contract explicit"),
        ("prospective_runner_bound", runner_exact, "dry-run default and exact non-destructive monitor controls present"),
        ("public_firewall", public_firewall, "no absolute private paths in public artifacts"),
        ("count_not_executed", count_not_executed, "no sentinel output directory detected by prefreeze"),
    ]
    rows = [{
        "contract_id": "mv08h_count_sentinel_independent_validation_v1",
        "check_order": order, "check_id": check_id, "passed": passed,
        "evidence": evidence,
    } for order, (check_id, passed, evidence) in enumerate(checks, 1)]
    if not all(row["passed"] for row in rows):
        failed = [row["check_id"] for row in rows if not row["passed"]]
        raise RuntimeError("independent prefreeze validation failed: " + ", ".join(failed))
    validation_name = "mv08h-count-sentinel-independent-validation.csv"
    write_csv(output / validation_name, rows)
    artifact_names = sorted(builder_names + [validation_name])
    write_csv(
        output / "mv08h-count-sentinel-artifact-manifest.csv",
        artifact_manifest(output, artifact_names),
    )
    print("MV8-H count-sentinel prefreeze independent validation passed: 14/14 checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
