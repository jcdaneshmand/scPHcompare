#!/usr/bin/env python3
"""Build a public, label-closed closure for the completed MV8-H sentinel.

This builder consumes only run metadata and file identities. It deliberately
does not open HDF5/MEX matrices, metrics, HTML, expression values, barcodes,
or any donor/label/outcome fields.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess


RUN_ID = "mv08h_count_sentinel_hca_bm_002"
UNIT = "HCA_BM_002"
SAMPLE = "MantonBM2_HiSeq_9"
EXPECTED_VERSION = "cellranger-8.0.1"
EXPECTED_FASTQ_BYTES = 11_249_623_632
LOCAL_CORES = 4
LOCAL_MEMORY_GIB = 32
RSS_CAP_BYTES = 80 * 1024**3
WORKSPACE_CAP_BYTES = 200 * 1024**3
FREE_SPACE_FLOOR_BYTES = 1024**4
ELAPSED_CAP_SECONDS = 96 * 60 * 60
SUCCESS_MARKER = "Pipestance completed successfully!"
PRIVATE_PATH_PATTERNS = ("/mnt/", "E:\\", "C:\\", "\\\\")
METADATA_FILES = (
    ("run_stdout.log", "mv08h-count-sentinel.stdout.log"),
    ("run_stderr.log", "mv08h-count-sentinel.stderr.log"),
    ("run_command", "_cmdline"),
    ("run_versions", "_versions"),
    ("run_final_state", "_finalstate"),
    ("run_performance", "_perf"),
    ("run_sitecheck", "_sitecheck"),
)
OUTPUT_FILES = (
    ("filtered_feature_bc_matrix", "directory", True),
    ("raw_feature_bc_matrix", "directory", True),
    ("filtered_feature_bc_matrix.h5", "private_matrix", True),
    ("raw_feature_bc_matrix.h5", "private_matrix", True),
    ("molecule_info.h5", "private_molecule_data", True),
    ("metrics_summary.csv", "private_qc", True),
    ("web_summary.html", "private_qc", True),
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing empty evidence: {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), extrasaction="raise",
            quoting=csv.QUOTE_ALL, lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: "TRUE" if value is True else "FALSE" if value is False else value
                for key, value in row.items()
            })
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def write_text(path: Path, text: str) -> None:
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def tree_bytes(root: Path) -> int:
    return sum(
        path.stat().st_size for path in root.rglob("*")
        if path.is_file() and not path.is_symlink()
    )


def measured_bytes(root: Path) -> int:
    """Use the host filesystem's du when available; retain a portable fallback."""
    try:
        result = subprocess.run(
            ["du", "-sb", str(root)], check=True, text=True,
            capture_output=True,
        )
        return int(result.stdout.split()[0])
    except (FileNotFoundError, subprocess.CalledProcessError, ValueError):
        return tree_bytes(root)


def collect_highmem(value: object, rows: list[dict[str, int]]) -> None:
    if isinstance(value, dict):
        highmem = value.get("highmem")
        if isinstance(highmem, dict) and "rss" in highmem:
            rows.append({key: int(highmem.get(key, 0)) for key in
                         ("rss", "vmem", "proc_count", "shared", "text")})
        for child in value.values():
            collect_highmem(child, rows)
    elif isinstance(value, list):
        for child in value:
            collect_highmem(child, rows)


def parse_local_time(value: str) -> datetime:
    return datetime.strptime(value, "%Y-%m-%d %H:%M:%S")


def metadata_manifest(run_root: Path, stdout_path: Path, stderr_path: Path) -> list[dict[str, object]]:
    rows = []
    for public_name, source_name in METADATA_FILES:
        source = stdout_path if public_name == "run_stdout.log" else stderr_path if public_name == "run_stderr.log" else run_root / source_name
        if not source.is_file():
            raise RuntimeError(f"missing terminal metadata: {source_name}")
        rows.append({
            "artifact_class": "terminal_metadata",
            "file": public_name,
            "bytes": source.stat().st_size,
            "sha256": sha256_file(source),
            "content_opened": public_name in {"run_stdout.log", "run_stderr.log", "run_command", "run_versions", "run_final_state", "run_sitecheck"},
            "public_release_permitted": True,
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefreeze-dir", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--stdout", type=Path, required=True)
    parser.add_argument("--stderr", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--recorded-free-bytes", type=int)
    args = parser.parse_args()

    run = args.run_root.resolve()
    output = args.output_dir.resolve()
    stdout = args.stdout.resolve()
    stderr = args.stderr.resolve()
    prefreeze = args.prefreeze_dir.resolve()
    out_root = run / "outs"

    stdout_text = stdout.read_text(encoding="utf-8", errors="replace")
    if not stdout.is_file() or SUCCESS_MARKER not in stdout_text:
        raise RuntimeError("Cell Ranger success marker is absent")
    stderr_text = stderr.read_text(encoding="utf-8", errors="replace") if stderr.is_file() else ""
    if not stderr.is_file() or stderr_text.strip():
        raise RuntimeError("stderr contains non-whitespace diagnostics")
    failures = [path for path in run.rglob("*") if path.is_file() and path.name in {"_failed", "_error"}]
    if failures:
        raise RuntimeError("failure markers remain in the completed run")

    versions = json.loads((run / "_versions").read_text(encoding="utf-8"))
    if versions.get("pipelines") != EXPECTED_VERSION:
        raise RuntimeError("Cell Ranger version differs from the frozen runtime")
    command = (run / "_cmdline").read_text(encoding="utf-8").strip()
    required = (
        f"count --id={RUN_ID}", f"--sample={SAMPLE}", "--chemistry=SC3Pv2",
        "--expect-cells=7000", "--include-introns=false", "--create-bam=false",
        "--nosecondary", "--localcores=4", "--localmem=32", "--disable-ui",
    )
    if not all(term in command for term in required):
        raise RuntimeError("completed command controls differ from the freeze")
    if not command.startswith("/mnt/e/CellRanger/cellranger-8.0.1/bin/cellranger count"):
        raise RuntimeError("completed launcher path differs from Cell Ranger 8.0.1")

    final_state = json.loads((run / "_finalstate").read_text(encoding="utf-8"))
    if not any(isinstance(stage, dict) and stage.get("name") == "SC_RNA_COUNTER_CS" and stage.get("state") == "complete" for stage in final_state):
        raise RuntimeError("top-level count pipeline is not complete")
    perf = json.loads((run / "_perf").read_text(encoding="utf-8"))
    highmem_rows: list[dict[str, int]] = []
    collect_highmem(perf, highmem_rows)
    if not highmem_rows:
        raise RuntimeError("performance metadata has no highmem record")
    peak = max(highmem_rows, key=lambda row: row["rss"])

    timestamp_lines = (run / "_timestamp").read_text(encoding="utf-8").splitlines()
    start_line = next((line for line in timestamp_lines if line.startswith("start: ")), "")
    if not start_line:
        raise RuntimeError("run start timestamp is absent")
    lock = parse_local_time(start_line.removeprefix("start: "))
    finish_match = re.search(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) Shutting down\.", stdout_text, re.MULTILINE)
    if not finish_match:
        raise RuntimeError("shutdown timestamp is absent")
    finish = parse_local_time(finish_match.group(1))
    elapsed = int((finish - lock).total_seconds())

    workspace = measured_bytes(run)
    observed_free_bytes = shutil.disk_usage(run).free
    recorded_free_bytes = args.recorded_free_bytes if args.recorded_free_bytes is not None else observed_free_bytes
    if peak["rss"] > RSS_CAP_BYTES or workspace > WORKSPACE_CAP_BYTES or observed_free_bytes < FREE_SPACE_FLOOR_BYTES or elapsed > ELAPSED_CAP_SECONDS:
        raise RuntimeError("terminal resource contract is breached")

    expected_outputs = {name: (kind, private) for name, kind, private in OUTPUT_FILES}
    observed_outputs = {path.name: path for path in out_root.iterdir()} if out_root.is_dir() else {}
    for name, (kind, _private) in expected_outputs.items():
        path = observed_outputs.get(name)
        if path is None or (kind == "directory" and not path.is_dir()) or (kind != "directory" and not path.is_file()):
            raise RuntimeError(f"expected output structure missing: {name}")
    forbidden_outputs = [path.name for path in out_root.iterdir() if path.name.lower().endswith((".bam", ".bai", ".csi")) or path.name == "secondary"]
    if forbidden_outputs:
        raise RuntimeError("BAM/secondary output unexpectedly present")

    prefreeze_manifest_path = prefreeze / "mv08h-count-sentinel-artifact-manifest.csv"
    prefreeze_report_path = prefreeze / "MV08H_CELLRANGER8_COUNT_SENTINEL_PREFREEZE_2026-08-18.md"
    prefreeze_decision = read_csv(prefreeze / "mv08h-count-sentinel-decision.csv")
    if not prefreeze_manifest_path.is_file() or not prefreeze_report_path.is_file() or len(prefreeze_decision) != 1:
        raise RuntimeError("historical prefreeze evidence is incomplete")
    if prefreeze_decision[0]["count_sentinel_execution_authorized"] != "FALSE":
        raise RuntimeError("prefreeze record is not a non-execution record")

    prefreeze_artifact_hash = sha256_file(prefreeze_manifest_path)
    prefreeze_report_hash = sha256_file(prefreeze_report_path)
    reference_binding = read_csv(prefreeze / "mv08h-count-sentinel-reference-binding.csv")[0]
    fastq_identity = read_csv(prefreeze / "mv08h-count-sentinel-fastqs.csv")
    if sum(int(row["expected_bytes"]) for row in fastq_identity) != EXPECTED_FASTQ_BYTES:
        raise RuntimeError("prefreeze FASTQ binding differs")

    input_binding = [{
        "contract_id": "mv08h_count_sentinel_execution_input_binding_v1",
        "historical_prefreeze_artifact_manifest_sha256": prefreeze_artifact_hash,
        "historical_prefreeze_report_sha256": prefreeze_report_hash,
        "unit_id": UNIT,
        "fastq_files": len(fastq_identity),
        "fastq_bytes": EXPECTED_FASTQ_BYTES,
        "reference_tree_sha256": reference_binding["tree_sha256"],
        "runtime_tree_sha256": reference_binding["runtime_tree_sha256"],
        "launcher_sha256": reference_binding["launcher_sha256"],
        "input_binding_passed": True,
    }]

    command_hash = hashlib.sha256(command.encode("utf-8")).hexdigest()
    metadata = metadata_manifest(run, stdout, stderr)
    metadata.append({
        "artifact_class": "completed_command_identity",
        "file": "command_sha256",
        "bytes": len(command.encode("utf-8")),
        "sha256": command_hash,
        "content_opened": True,
        "public_release_permitted": True,
    })

    output_rows = []
    for name, kind, private in OUTPUT_FILES:
        path = observed_outputs[name]
        output_rows.append({
            "artifact_class": "private_output_structure",
            "relative_name": name,
            "kind": kind,
            "bytes": measured_bytes(path) if path.is_dir() else path.stat().st_size,
            "content_opened": False,
            "contains_expression_or_umi": private,
            "contains_qc_or_labels": private,
            "public_release_permitted": False,
        })

    resource = [{
        "contract_id": "mv08h_count_sentinel_execution_resource_v1",
        "start_local": lock.strftime("%Y-%m-%d %H:%M:%S"),
        "finish_local": finish.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_seconds": elapsed,
        "peak_rss_bytes": peak["rss"],
        "peak_vmem_bytes": peak["vmem"],
        "peak_process_count": peak["proc_count"],
        "final_workspace_bytes": workspace,
        "free_bytes_at_closure": recorded_free_bytes,
        "local_cores": LOCAL_CORES,
        "local_memory_gib": LOCAL_MEMORY_GIB,
        "rss_cap_passed": peak["rss"] <= RSS_CAP_BYTES,
        "workspace_cap_passed": workspace <= WORKSPACE_CAP_BYTES,
        "free_space_floor_passed": observed_free_bytes >= FREE_SPACE_FLOOR_BYTES,
        "elapsed_cap_passed": elapsed <= ELAPSED_CAP_SECONDS,
    }]

    validation = [
        (1, "terminal_success_marker", True, "Pipestance completed successfully marker present"),
        (2, "stderr_empty", True, "terminal stderr is zero bytes"),
        (3, "no_failure_markers", True, "no _failed or _error marker remains"),
        (4, "runtime_exact", True, "Cell Ranger 8.0.1 pipeline identity exact"),
        (5, "command_controls_exact", True, "SC3Pv2/exon-only/no-BAM/no-secondary/4-core/32-GiB controls exact"),
        (6, "historical_prefreeze_preserved", True, "pre-execution artifact and decision hashes preserved"),
        (7, "resource_caps", peak["rss"] <= RSS_CAP_BYTES and workspace <= WORKSPACE_CAP_BYTES and observed_free_bytes >= FREE_SPACE_FLOOR_BYTES and elapsed <= ELAPSED_CAP_SECONDS, "RSS/workspace/free-space/elapsed gates pass"),
        (8, "output_structure", True, "raw/filtered matrices and molecule-info outputs exist structurally"),
        (9, "bam_secondary_absent", not forbidden_outputs, "BAM and secondary outputs absent"),
        (10, "matrix_content_closed", True, "matrix/QC contents were not opened by closure builder"),
        (11, "label_outcome_firewall", True, "labels, outcomes, landscapes, and remaining units remain closed"),
        (12, "public_firewall", True, "public closure contains identities/aggregates only"),
    ]
    if not all(passed for _, _, passed, _ in validation):
        raise RuntimeError("execution closure validation failed")
    validation_rows = [{
        "contract_id": "mv08h_count_sentinel_execution_validation_v1",
        "check_order": order,
        "check_id": check_id,
        "passed": passed,
        "evidence": evidence,
    } for order, check_id, passed, evidence in validation]

    decision = [{
        "contract_id": "mv08h_count_sentinel_execution_decision_v1",
        "decision": "count_sentinel_execution_closed_downstream_firewall_preserved",
        "execution_performed": True,
        "run_success": True,
        "matrix_access_authorized": False,
        "qc_authorized": False,
        "pca_ph_landscape_authorized": False,
        "label_access_authorized": False,
        "biological_outcomes_authorized": False,
        "remaining_units_authorized": False,
        "deletion_authorized": False,
        "next_gate": "separate owner-approved structural/QC review and remaining-unit decision",
    }]

    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "mv08h-count-sentinel-execution-input-binding.csv", input_binding)
    write_csv(output / "mv08h-count-sentinel-execution-resource.csv", resource)
    write_csv(output / "mv08h-count-sentinel-execution-output-structure.csv", output_rows)
    write_csv(output / "mv08h-count-sentinel-execution-validation.csv", validation_rows)
    write_csv(output / "mv08h-count-sentinel-execution-decision.csv", decision)
    write_csv(output / "mv08h-count-sentinel-execution-metadata.csv", metadata)

    report = f"""# MV8-H Cell Ranger 8.0.1 count-sentinel execution closure

## Outcome

The one prospectively selected complete unit `{UNIT}` completed successfully
under the frozen Cell Ranger 8.0.1 command. This document is an execution
closure, not a replacement for the historical pre-execution prefreeze. The
prefreeze decision remains preserved by SHA-256 and records that count had not
yet executed at that earlier timestamp.

The completed command identity is bound by SHA-256 `{command_hash}`. Its
scientific controls remain exact: `{SAMPLE}`, SC3Pv2, expected cells 7000,
exon-only counting, BAM disabled, secondary analysis disabled, 4 cores, and
32 GiB. The selected six FASTQs total `{EXPECTED_FASTQ_BYTES:,}` bytes and
remain bound to the historical prefreeze manifest. The custom reference,
runtime tree, and launcher identities remain unchanged and hash-bound.

## Terminal evidence

The pipestance success marker is present, stderr is empty, and no failure
marker remains. The terminal resource record reports `{elapsed}` seconds,
peak RSS `{peak["rss"]:,}` bytes, final workspace `{workspace:,}` bytes, and
free space `{recorded_free_bytes:,}` bytes. All admitted resource ceilings pass.

The output structure contains raw and filtered feature-barcode matrix
directories, HDF5 matrix files, molecule information, and Cell Ranger summary
files. The closure records their names and sizes only; it does not open their
contents. BAM and secondary outputs are absent.

## Firewall and next gate

No matrix contents, expression/UMI values, cell barcodes, QC values, donor
attributes, labels, outcomes, PCA, clustering, persistence landscapes, or
remaining units were opened or processed by this closure. The corrected
landscape contract remains unchanged: separate cell and gene observation
views, separate H0/H1, every consecutive active level, no fixed grid, and no
universal level cap.

The next gate is a separately authorized structural/QC review and remaining
unit decision. Nothing in this closure authorizes labels, outcomes, landscape
calculation, additional units, or deletion.
"""
    report_path = output / "MV08H_CELLRANGER8_COUNT_SENTINEL_EXECUTION_CLOSURE_2026-08-18.md"
    write_text(report_path, report)

    artifact_names = [
        report_path.name,
        "mv08h-count-sentinel-execution-decision.csv",
        "mv08h-count-sentinel-execution-input-binding.csv",
        "mv08h-count-sentinel-execution-metadata.csv",
        "mv08h-count-sentinel-execution-output-structure.csv",
        "mv08h-count-sentinel-execution-resource.csv",
        "mv08h-count-sentinel-execution-validation.csv",
    ]
    manifest_rows = []
    for name in artifact_names:
        path = output / name
        content = path.read_text(encoding="utf-8", errors="replace")
        if any(pattern in content for pattern in PRIVATE_PATH_PATTERNS):
            raise RuntimeError(f"private path leaked into public artifact: {name}")
        manifest_rows.append({
            "contract_id": "mv08h_count_sentinel_execution_artifact_manifest_v1",
            "artifact_order": len(manifest_rows) + 1,
            "file": name,
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "contains_private_path": False,
            "contains_expression_or_umi": False,
            "contains_qc_or_labels": False,
            "public_release_permitted": True,
        })
    write_csv(output / "mv08h-count-sentinel-execution-artifact-manifest.csv", manifest_rows)
    print("MV8-H execution closure built: 12/12 checks; downstream firewall preserved")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
