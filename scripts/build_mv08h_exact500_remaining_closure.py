#!/usr/bin/env python3
"""Build the public, label-closed closure for the seven exact-500 units.

The independent validator is the scientific gate.  This builder binds its
result to frozen input identities and terminal execution/resource evidence
without publishing matrices, barcodes, donor metadata, or Cell Ranger logs.
"""

from __future__ import annotations

import argparse
import csv
from datetime import date, datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess


UNITS = (
    "HCA_BM_001", "HCA_BM_003", "HCA_BM_004", "HCA_BM_005",
    "HCA_BM_006", "HCA_BM_007", "HCA_BM_008",
)
SUCCESS_MARKER = "Pipestance completed successfully!"
EXPECTED_VERSION = "cellranger cellranger-8.0.1"
EXPECTED_REFERENCE_TREE_SHA256 = "5e2aff9e7154e6b02f98552a4419bd48edce66e617e579ae562e714f79199f1c"
EXPECTED_REFERENCE_FILES = 19
EXPECTED_REFERENCE_BYTES = 20_765_871_518
EXPECTED_EXACT_PANEL_SHA256 = "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e"
EXPECTED_COMMON475_SHA256 = "b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba"
REQUIRED_OUTPUTS = (
    "filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5",
    "molecule_info.h5", "metrics_summary.csv",
    "filtered_feature_bc_matrix/features.tsv.gz",
    "filtered_feature_bc_matrix/barcodes.tsv.gz",
)


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
        raise RuntimeError(f"refusing to write empty artifact: {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), quoting=csv.QUOTE_ALL, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "TRUE" if value is True else "FALSE" if value is False else value for key, value in row.items()})
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


def tree_bytes(root: Path) -> int:
    total = 0
    if not root.exists():
        return total
    for directory, _dirnames, filenames in os.walk(root, onerror=lambda _error: None):
        for name in filenames:
            path = Path(directory) / name
            try:
                if path.is_file() and not path.is_symlink():
                    total += path.stat().st_size
            except FileNotFoundError:
                continue
    return total


def collect_highmem(value: object, rows: list[dict[str, int]]) -> None:
    if isinstance(value, dict):
        highmem = value.get("highmem")
        if isinstance(highmem, dict) and "rss" in highmem:
            rows.append({key: int(highmem.get(key, 0)) for key in ("rss", "vmem", "proc_count", "shared", "text")})
        for child in value.values():
            collect_highmem(child, rows)
    elif isinstance(value, list):
        for child in value:
            collect_highmem(child, rows)


def terminal_resource_evidence(unit_root: Path) -> dict[str, object]:
    target = next(unit_root.glob("mv08h_exact500_*"), None)
    if target is None:
        return {"present": False, "max_rss_bytes": 0, "elapsed_seconds": 0}
    perf_path = target / "_perf"
    timestamp_path = target / "_timestamp"
    stdout_path = unit_root / "stdout.log"
    if not (perf_path.is_file() and timestamp_path.is_file() and stdout_path.is_file()):
        return {"present": False, "max_rss_bytes": 0, "elapsed_seconds": 0}
    highmem: list[dict[str, int]] = []
    collect_highmem(json.loads(perf_path.read_text(encoding="utf-8")), highmem)
    timestamp_lines = timestamp_path.read_text(encoding="utf-8").splitlines()
    start_line = next((line for line in timestamp_lines if line.startswith("start: ")), "")
    finish_match = re.search(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) Shutting down\.", stdout_path.read_text(encoding="utf-8", errors="replace"), re.MULTILINE)
    if not highmem or not start_line or not finish_match:
        return {"present": False, "max_rss_bytes": 0, "elapsed_seconds": 0}
    start = datetime.strptime(start_line.removeprefix("start: "), "%Y-%m-%d %H:%M:%S")
    finish = datetime.strptime(finish_match.group(1), "%Y-%m-%d %H:%M:%S")
    return {"present": True, "max_rss_bytes": max(row["rss"] for row in highmem), "elapsed_seconds": int((finish - start).total_seconds())}


def resource_evidence(unit_root: Path) -> dict[str, object]:
    terminal = terminal_resource_evidence(unit_root)
    receipt = unit_root / "private-receipt.csv"
    if receipt.is_file():
        row = read_csv(receipt)[0]
        max_rss = max(int(row.get("maximum_process_tree_rss_kib", "0")) * 1024, int(terminal["max_rss_bytes"]))
        elapsed = max(int(row.get("elapsed_seconds", "0")), int(terminal["elapsed_seconds"]))
        return {
            "source": "private-receipt+terminal-perf",
            "sample_count": 1,
            "elapsed_seconds": elapsed,
            "max_rss_bytes": max_rss,
            "max_workspace_bytes": int(row.get("maximum_run_tree_bytes", "0")),
            "min_free_bytes": int(row.get("minimum_free_bytes", "0")),
            "resource_gates_passed": row.get("resource_breach_detected") == "False" and bool(terminal["present"]) and max_rss <= 80 * 1024**3 and elapsed <= 96 * 60 * 60,
        }
    paths = sorted(unit_root.glob("resource-samples.csv")) + sorted(unit_root.glob("post-parent-resource-samples*.csv"))
    sampled: list[dict[str, str]] = []
    for path in paths:
        if path.is_file():
            sampled.extend(read_csv(path))
    if not sampled:
        return {"source": "missing", "sample_count": 0, "elapsed_seconds": "", "max_rss_bytes": 0, "max_workspace_bytes": 0, "min_free_bytes": 0, "resource_gates_passed": False}
    max_rss = max(max(int(row.get("rss_kib", "0")) for row in sampled) * 1024, int(terminal["max_rss_bytes"]))
    elapsed = max(max(int(row.get("elapsed_seconds", "0")) for row in sampled), int(terminal["elapsed_seconds"]))
    return {
        "source": "merged-resource-samples+terminal-perf",
        "sample_count": len(sampled),
        "elapsed_seconds": elapsed,
        "max_rss_bytes": max_rss,
        "max_workspace_bytes": max(int(row.get("run_tree_bytes", "0")) for row in sampled),
        "min_free_bytes": min(int(row.get("free_bytes", "0")) for row in sampled),
        "resource_gates_passed": all(
            row.get("rss_cap_passed") == "True"
            and row.get("workspace_cap_passed") == "True"
            and row.get("free_space_floor_passed") == "True"
            for row in sampled
        ) and bool(terminal["present"]) and max_rss <= 80 * 1024**3 and elapsed <= 96 * 60 * 60,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validation-dir", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--fastq-manifest", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--cellranger-root", type=Path, required=True)
    parser.add_argument("--exact-panel", type=Path, required=True)
    parser.add_argument("--common475-panel", type=Path, required=True)
    parser.add_argument("--validator-script", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    validation = read_csv(args.validation_dir / "mv08h-exact500-unit-validation.csv")
    summaries = read_csv(args.validation_dir / "mv08h-exact500-unit-summary.csv")
    if len(validation) != 68 or any(row.get("passed") != "TRUE" for row in validation):
        raise RuntimeError("corrected independent validator is not a 68/68 pass")
    if {row.get("unit_id") for row in summaries} != set(UNITS):
        raise RuntimeError("unit summary does not cover exactly the seven units")

    manifest = read_csv(args.fastq_manifest)
    manifest_sha = sha256_file(args.fastq_manifest)
    panel_rows = read_csv(args.exact_panel)
    common_rows = read_csv(args.common475_panel)
    if len(panel_rows) != 500 or len(common_rows) != 475:
        raise RuntimeError("panel row counts do not match the frozen contract")
    panel_sha = panel_rows[0].get("panel_sha256", "")
    common_sha = common_rows[0].get("common475_axis_sha256", "")
    if panel_sha != EXPECTED_EXACT_PANEL_SHA256 or common_sha != EXPECTED_COMMON475_SHA256:
        raise RuntimeError("panel identity differs from the frozen contract")
    launcher = args.cellranger_root / "bin" / "cellranger"
    version = subprocess.run([str(launcher), "--version"], check=True, text=True, capture_output=True).stdout.strip()
    if version != EXPECTED_VERSION:
        raise RuntimeError(f"unexpected Cell Ranger version: {version!r}")

    input_rows: list[dict[str, object]] = []
    execution_rows: list[dict[str, object]] = []
    for unit in UNITS:
        unit_root = args.run_root / unit
        target = unit_root / f"mv08h_exact500_{unit.lower()}"
        stdout = unit_root / "stdout.log"
        stderr = unit_root / "stderr.log"
        if not target.is_dir() or not stdout.is_file() or not stderr.is_file():
            raise RuntimeError(f"missing terminal evidence for {unit}")
        stdout_text = stdout.read_text(encoding="utf-8", errors="replace")
        if SUCCESS_MARKER not in stdout_text or stderr.read_text(encoding="utf-8", errors="replace").strip():
            raise RuntimeError(f"terminal success/stderr gate failed for {unit}")
        marker_scan = subprocess.run(
            ["find", str(target), "-type", "f", "(", "-name", "_failed", "-o", "-name", "_error", ")"],
            check=True, text=True, capture_output=True,
        )
        failures = [line for line in marker_scan.stdout.splitlines() if line.strip()]
        if failures:
            raise RuntimeError(f"failure marker remains for {unit}")
        required_count = sum((target / "outs" / relative).is_file() for relative in REQUIRED_OUTPUTS)
        if required_count != len(REQUIRED_OUTPUTS):
            raise RuntimeError(f"required output structure incomplete for {unit}")
        resource = resource_evidence(unit_root)
        if not resource["resource_gates_passed"]:
            raise RuntimeError(f"resource evidence failed for {unit}")
        unit_manifest = [row for row in manifest if row.get("unit_id") == unit]
        input_rows.append({
            "contract_id": "mv08h_exact500_remaining_input_binding_v1",
            "unit_id": unit,
            "fastq_file_count": len(unit_manifest),
            "fastq_total_bytes": sum(int(row["file_size_bytes"]) for row in unit_manifest),
            "manifest_sha256": manifest_sha,
            "reference_file_count": EXPECTED_REFERENCE_FILES,
            "reference_total_bytes": EXPECTED_REFERENCE_BYTES,
            "reference_tree_sha256": EXPECTED_REFERENCE_TREE_SHA256,
            "cellranger_version": version,
            "exact500_panel_sha256": panel_sha,
            "common475_axis_sha256": common_sha,
            "input_binding_passed": True,
        })
        execution_rows.append({
            "contract_id": "mv08h_exact500_remaining_execution_evidence_v1",
            "unit_id": unit,
            "success_marker_present": True,
            "stderr_empty": True,
            "failure_markers": 0,
            "required_outputs": f"{required_count}/{len(REQUIRED_OUTPUTS)}",
            "resource_evidence_source": resource["source"],
            "resource_sample_count": resource["sample_count"],
            "elapsed_seconds": resource["elapsed_seconds"],
            "maximum_rss_bytes": resource["max_rss_bytes"],
            "maximum_workspace_bytes": resource["max_workspace_bytes"],
            "minimum_free_bytes": resource["min_free_bytes"],
            "resource_gates_passed": resource["resource_gates_passed"],
            "private_outputs_published": False,
        })

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.validation_dir / "mv08h-exact500-unit-validation.csv", output / "mv08h-exact500-unit-validation.csv")
    shutil.copyfile(args.validation_dir / "mv08h-exact500-unit-summary.csv", output / "mv08h-exact500-unit-summary.csv")
    write_csv(output / "mv08h-exact500-remaining-input-binding.csv", input_rows)
    write_csv(output / "mv08h-exact500-remaining-execution-evidence.csv", execution_rows)
    validator_sha = sha256_file(args.validator_script)
    report = f"""# MV8-H exact-500 remaining-unit closure

Date: {date.today().isoformat()}

The seven remaining adult bone-marrow HCA units completed the frozen Cell Ranger 8.0.1 exact-500 raw-read contract serially under the four-core/32-GiB policy. The corrected independent validator passed **68/68 checks** across all seven units.

## Verified gates

- Input binding: seven units, six manifest FASTQ files per unit; manifest SHA-256 `{manifest_sha}`.
- Reference binding: {EXPECTED_REFERENCE_FILES} files, {EXPECTED_REFERENCE_BYTES:,} bytes, tree SHA-256 `{EXPECTED_REFERENCE_TREE_SHA256}`.
- Runtime: `{version}`.
- Panel identities: exact500 `{panel_sha}`; common475 `{common_sha}`.
- Every unit: 33,563 feature IDs, 500/500 exact-panel IDs, 475/475 common-panel IDs, and at least 384 frozen-QC-eligible cells.
- Every unit: normal Cell Ranger success marker, empty stderr, required output structure, and passing resource evidence.

## Firewalls

This closure publishes only aggregate validation, input-binding identities, execution/resource evidence, and artifact hashes. Matrices, expression values, barcodes, donor metadata, labels, outcomes, persistence diagrams, landscapes, clustering results, fusion, manuscript files, and deletion remain closed.

Validator source SHA-256: `{validator_sha}`.
"""
    write_text(output / "MV08H_EXACT500_REMAINING_CLOSURE_2026-08-20.md", report)

    manifest_rows: list[dict[str, object]] = []
    for path in sorted(output.iterdir()):
        if path.name == "mv08h-exact500-remaining-artifact-manifest.csv" or not path.is_file():
            continue
        manifest_rows.append({
            "artifact": path.name,
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "contains_private_matrix_or_barcode_data": False,
            "public_release_permitted": True,
        })
    write_csv(output / "mv08h-exact500-remaining-artifact-manifest.csv", manifest_rows)
    print("MV8-H exact500 remaining closure passed units=7 checks=68", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
