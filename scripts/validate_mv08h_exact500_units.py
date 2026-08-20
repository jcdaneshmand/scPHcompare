#!/usr/bin/env python3
"""Independently validate the seven MV08-H exact-500 Cell Ranger outputs.

The validator opens only frozen FASTQ identities, runtime/reference identities,
feature axes, aggregate matrix/QC values, and execution/resource evidence.  It
never reads barcode identifiers, donor metadata, labels, outcomes, PCA,
persistence, landscapes, or clustering artifacts, and it writes aggregate
public evidence only.
"""

from __future__ import annotations

import argparse
import csv
from datetime import date, datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

import h5py
import numpy as np
from scipy import sparse


UNITS = (
    "HCA_BM_001", "HCA_BM_003", "HCA_BM_004", "HCA_BM_005",
    "HCA_BM_006", "HCA_BM_007", "HCA_BM_008",
)
FEATURES = 33_563
MIN_ENTRY_FEATURES = 200
MIN_ENTRY_CELLS = 3
MIN_NFEATURE = 500
MAX_NFEATURE = 9_000
MAX_MITO_PERCENT = 20.0
MIN_RIBO_PERCENT = 5.0
MIN_POST_QC_CELLS = 384
EXPECTED_VERSION = "cellranger cellranger-8.0.1"
EXPECTED_LAUNCHER_BYTES = 19_924_984
EXPECTED_LAUNCHER_SHA256 = "4ee3a1670b4f14c826004fe8e17b4759e1edc701b15ff2e9623753bf1b34d4d6"
EXPECTED_REFERENCE_FILES = 19
EXPECTED_REFERENCE_BYTES = 20_765_871_518
EXPECTED_REFERENCE_TREE_SHA256 = "5e2aff9e7154e6b02f98552a4419bd48edce66e617e579ae562e714f79199f1c"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def exact_file(path: Path, size: int, digest: str) -> bool:
    return path.is_file() and path.stat().st_size == size and sha256_file(path) == digest


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
            regular += 1
            regular_bytes += size
            digest.update(f"F\0{relative}\0{size}\0{sha256_file(path)}\n".encode())
        else:
            symlinks += 1
            digest.update(f"L\0{relative}\0{os.readlink(path)}\n".encode())
    return regular, symlinks, regular_bytes, digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_panels(exact_path: Path, common_path: Path) -> tuple[list[str], list[str]]:
    exact_rows = rows(exact_path)
    common_rows = rows(common_path)
    exact_ids = []
    for row in exact_rows:
        match = re.search(r"(ENSG[0-9]+)", row.get("feature_id", ""))
        if not match:
            raise ValueError("exact panel contains a row without an Ensembl stable ID")
        exact_ids.append(match.group(1))
    common_ids = [row["ensembl_stable_id"] for row in common_rows]
    return exact_ids, common_ids


def load_matrix(path: Path) -> tuple[sparse.csc_matrix, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5:
        shape = tuple(int(value) for value in h5["matrix/shape"][()])
        if shape[0] != FEATURES or shape[1] <= 0:
            raise ValueError(f"unexpected shape in {path.name}: {shape}")
        data = np.asarray(h5["matrix/data"][()])
        indices = np.asarray(h5["matrix/indices"][()]).astype(np.int64, copy=False)
        indptr = np.asarray(h5["matrix/indptr"][()]).astype(np.int64, copy=False)
        ids_raw = np.asarray(h5["matrix/features/id"][()])
        names_raw = np.asarray(h5["matrix/features/name"][()])
        feature_ids = np.asarray([
            value.decode("utf-8") if isinstance(value, (bytes, np.bytes_)) else str(value)
            for value in ids_raw
        ], dtype=object)
        names = np.asarray([
            value.decode("utf-8") if isinstance(value, (bytes, np.bytes_)) else str(value)
            for value in names_raw
        ], dtype=object)
        feature_types = h5["matrix/features/feature_type"][()]
        if data.dtype != np.int32 or indices.dtype != np.int64 or indptr.dtype != np.int64:
            raise ValueError(f"unexpected dtypes in {path.name}")
        if len(feature_ids) != shape[0] or len(names) != shape[0] or len(feature_types) != shape[0] or len(indptr) != shape[1] + 1:
            raise ValueError(f"axis mismatch in {path.name}")
        if indptr[-1] != len(data) or np.any(data < 0):
            raise ValueError(f"invalid sparse matrix values in {path.name}")
    return sparse.csc_matrix((data, indices, indptr), shape=shape), feature_ids, names


def qc_summary(matrix: sparse.csc_matrix, names: np.ndarray) -> dict[str, object]:
    text = names.astype(str)
    mito = np.char.startswith(text, "MT-")
    ribo = np.char.startswith(text, "RPS") | np.char.startswith(text, "RPL")
    nfeature = np.asarray(matrix.getnnz(axis=0)).ravel().astype(float)
    counts = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    if np.any(counts <= 0):
        raise ValueError("filtered matrix contains a zero-count cell")
    percent_mito = 100.0 * np.asarray(matrix[mito, :].sum(axis=0)).ravel() / counts
    percent_ribo = 100.0 * np.asarray(matrix[ribo, :].sum(axis=0)).ravel() / counts
    entry_cells = nfeature >= MIN_ENTRY_FEATURES
    entry_features = np.asarray(matrix[:, entry_cells].getnnz(axis=1)).ravel() >= MIN_ENTRY_CELLS
    qc_cells = (
        entry_cells & (nfeature >= MIN_NFEATURE) & (nfeature <= MAX_NFEATURE)
        & (percent_mito <= MAX_MITO_PERCENT) & (percent_ribo > MIN_RIBO_PERCENT)
    )
    final_features = int(np.sum(np.asarray(matrix[:, qc_cells].getnnz(axis=1)).ravel() > 3))
    return {
        "filtered_cells": int(matrix.shape[1]), "filtered_features": int(matrix.shape[0]),
        "filtered_nnz": int(matrix.nnz), "entry_cells_min_200_features": int(entry_cells.sum()),
        "entry_features_min_3_cells": int(entry_features.sum()), "post_qc_cells": int(qc_cells.sum()),
        "final_features_gt_3_cells": final_features, "nfeature_min": float(nfeature.min()),
        "nfeature_median": float(np.median(nfeature)), "nfeature_max": float(nfeature.max()),
        "count_min": float(counts.min()), "count_median": float(np.median(counts)),
        "count_max": float(counts.max()), "percent_mito_median": float(np.median(percent_mito)),
        "percent_ribo_median": float(np.median(percent_ribo)),
    }


def add(checks: list[dict[str, object]], unit: str, check_id: str, passed: bool, evidence: str) -> None:
    checks.append({"unit_id": unit, "check_id": check_id, "passed": str(bool(passed)).upper(), "evidence": evidence})


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
    return {
        "present": True,
        "max_rss_bytes": max(row["rss"] for row in highmem),
        "elapsed_seconds": int((finish - start).total_seconds()),
    }


def resource_evidence(unit_root: Path) -> tuple[bool, dict[str, object]]:
    terminal = terminal_resource_evidence(unit_root)
    receipt = unit_root / "private-receipt.csv"
    if receipt.is_file():
        row = rows(receipt)[0]
        sample_passed = (
            row.get("resource_breach_detected", "True") == "False"
            and int(row.get("maximum_process_tree_rss_kib", "0")) * 1024 <= 80 * 1024**3
            and int(row.get("maximum_run_tree_bytes", "0")) <= 200 * 1024**3
            and int(row.get("minimum_free_bytes", "0")) >= 1024**4
        )
        max_rss = max(int(row.get("maximum_process_tree_rss_kib", "0")) * 1024, int(terminal["max_rss_bytes"]))
        elapsed = max(int(row.get("elapsed_seconds", "0")), int(terminal["elapsed_seconds"]))
        passed = sample_passed and bool(terminal["present"]) and max_rss <= 80 * 1024**3 and elapsed <= 96 * 60 * 60
        return passed, {**row, "source": "private-receipt+terminal-perf", "maximum_process_tree_rss_kib": max_rss // 1024, "elapsed_seconds": elapsed}
    monitor_paths = sorted(unit_root.glob("resource-samples.csv")) + sorted(
        unit_root.glob("post-parent-resource-samples*.csv")
    )
    monitor_paths = [path for path in monitor_paths if path.is_file()]
    if not monitor_paths:
        return False, {"source": "missing-resource-evidence"}
    sampled: list[dict[str, str]] = []
    for path in monitor_paths:
        sampled.extend(rows(path))
    if not sampled:
        return False, {"source": "merged-resource-samples", "rows": 0}
    sample_passed = all(
        row["rss_cap_passed"] == "True" and row["workspace_cap_passed"] == "True"
        and row["free_space_floor_passed"] == "True" for row in sampled
    )
    max_rss = max(max(int(row["rss_kib"]) for row in sampled) * 1024, int(terminal["max_rss_bytes"]))
    elapsed = max(max(int(row.get("elapsed_seconds", "0")) for row in sampled), int(terminal["elapsed_seconds"]))
    passed = sample_passed and bool(terminal["present"]) and max_rss <= 80 * 1024**3 and elapsed <= 96 * 60 * 60
    return passed, {
        "source": "merged-resource-samples+terminal-perf", "rows": len(sampled),
        "maximum_process_tree_rss_kib": max_rss // 1024,
        "elapsed_seconds": elapsed,
        "maximum_run_tree_bytes": max(int(row["run_tree_bytes"]) for row in sampled),
        "minimum_free_bytes": min(int(row["free_bytes"]) for row in sampled),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--fastq-cache", type=Path, required=True)
    parser.add_argument("--fastq-manifest", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--cellranger-root", type=Path, required=True)
    parser.add_argument("--exact-panel", type=Path, required=True)
    parser.add_argument("--common475-panel", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    manifest_rows = rows(args.fastq_manifest)
    by_unit = {unit: [row for row in manifest_rows if row.get("unit_id") == unit] for unit in UNITS}
    launcher = args.cellranger_root / "bin" / "cellranger"
    launcher_ok = exact_file(launcher, EXPECTED_LAUNCHER_BYTES, EXPECTED_LAUNCHER_SHA256)
    version = ""
    if launcher_ok:
        version = subprocess.run([str(launcher), "--version"], check=True, text=True, capture_output=True).stdout.strip()
    add(checks, "ALL", "runtime_launcher_identity", launcher_ok and version == EXPECTED_VERSION, f"version={version}; sha256={sha256_file(launcher) if launcher_ok else 'missing'}")
    ref_files, ref_links, ref_bytes, ref_sha = tree_identity(args.reference_dir)
    add(checks, "ALL", "reference_tree_identity", (ref_files, ref_links, ref_bytes, ref_sha) == (EXPECTED_REFERENCE_FILES, 0, EXPECTED_REFERENCE_BYTES, EXPECTED_REFERENCE_TREE_SHA256), f"files={ref_files}; bytes={ref_bytes}; tree_sha256={ref_sha}")
    exact_ids, common_ids = read_panels(args.exact_panel, args.common475_panel)
    add(checks, "ALL", "exact500_panel_contract", len(exact_ids) == 500 and len(set(exact_ids)) == 500, f"rows={len(exact_ids)}")
    add(checks, "ALL", "common475_panel_contract", len(common_ids) == 475 and len(set(common_ids)) == 475, f"rows={len(common_ids)}")
    add(checks, "ALL", "common475_subset_exact500", set(common_ids) <= set(exact_ids), "ordered common475 IDs are contained in exact500 IDs")
    for unit in UNITS:
        unit_root = args.run_root / unit
        target = unit_root / f"mv08h_exact500_{unit.lower()}"
        outs = target / "outs"
        unit_rows = by_unit[unit]
        input_ok = len(unit_rows) == 6
        observed = sorted(path.name for path in (args.fastq_cache / unit).glob("*.fastq.gz"))
        expected = sorted(row["file_name"] for row in unit_rows)
        input_ok = input_ok and observed == expected
        for row in unit_rows:
            input_ok = input_ok and exact_file(args.fastq_cache / unit / row["file_name"], int(row["file_size_bytes"]), row["file_sha256"])
        add(checks, unit, "input_binding", input_ok, f"manifest_files={len(unit_rows)}; observed_files={len(observed)}")
        required = [outs / "filtered_feature_bc_matrix.h5", outs / "raw_feature_bc_matrix.h5", outs / "molecule_info.h5", outs / "metrics_summary.csv", outs / "filtered_feature_bc_matrix" / "features.tsv.gz", outs / "filtered_feature_bc_matrix" / "barcodes.tsv.gz"]
        output_ok = all(path.is_file() for path in required)
        add(checks, unit, "output_structure", output_ok, f"required_files={sum(path.is_file() for path in required)}/{len(required)}")
        if not output_ok:
            continue
        filtered, feature_ids_raw, names = load_matrix(outs / "filtered_feature_bc_matrix.h5")
        raw, _, _ = load_matrix(outs / "raw_feature_bc_matrix.h5")
        feature_ids = [re.sub(r"\.[0-9]+$", "", value) for value in feature_ids_raw.astype(str)]
        feature_set = set(feature_ids)
        add(checks, unit, "feature_axis", len(feature_ids) == FEATURES and len(feature_set) == FEATURES, f"unique_features={len(feature_set)}")
        add(checks, unit, "exact500_present", set(exact_ids) <= feature_set, f"present={sum(value in feature_set for value in exact_ids)}/500")
        add(checks, unit, "common475_present", set(common_ids) <= feature_set, f"present={sum(value in feature_set for value in common_ids)}/475")
        add(checks, unit, "raw_filtered_shapes", raw.shape[0] == FEATURES and raw.shape[1] >= filtered.shape[1] > 0, f"raw={raw.shape}; filtered={filtered.shape}")
        qc = qc_summary(filtered, names)
        add(checks, unit, "qc_384_eligibility", qc["post_qc_cells"] >= MIN_POST_QC_CELLS, f"post_qc_cells={qc['post_qc_cells']}")
        stdout = (unit_root / "stdout.log").read_text(encoding="utf-8", errors="replace") if (unit_root / "stdout.log").is_file() else ""
        stderr = (unit_root / "stderr.log").read_text(encoding="utf-8", errors="replace") if (unit_root / "stderr.log").is_file() else "missing"
        add(checks, unit, "execution_success_marker", "Pipestance completed successfully!" in stdout and not stderr.strip(), "Cell Ranger success marker present and stderr empty")
        resource_ok, resource = resource_evidence(unit_root)
        add(checks, unit, "resource_closure", resource_ok, f"source={resource.get('source')}; rows={resource.get('rows', 'receipt')}; max_rss_kib={resource.get('maximum_process_tree_rss_kib', '')}")
        summaries.append({"unit_id": unit, "filtered_cells": qc["filtered_cells"], "post_qc_cells": qc["post_qc_cells"], "filtered_nnz": qc["filtered_nnz"], "raw_cells": raw.shape[1], "exact500_present": sum(value in feature_set for value in exact_ids), "common475_present": sum(value in feature_set for value in common_ids), "resource_evidence": resource.get("source", "missing"), "resource_passed": str(resource_ok).upper()})
    check_fields = ["unit_id", "check_id", "passed", "evidence"]
    with (args.output_dir / "mv08h-exact500-unit-validation.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=check_fields, extrasaction="raise", quoting=csv.QUOTE_ALL); writer.writeheader(); writer.writerows(checks)
    if summaries:
        with (args.output_dir / "mv08h-exact500-unit-summary.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summaries[0]), extrasaction="raise", quoting=csv.QUOTE_ALL); writer.writeheader(); writer.writerows(summaries)
    passed = sum(row["passed"] == "TRUE" for row in checks)
    report = [
        "# MV8-H exact-500 remaining-unit validation",
        "",
        f"Validation date: {date.today().isoformat()}",
        "",
        "This independent closure opens only frozen input identities, the Cell Ranger 8.0.1/reference identities, feature axes, aggregate matrix/QC values, and execution/resource evidence.",
        "Barcode IDs, donor metadata, labels, outcomes, PCA, persistence, landscapes, clustering, and manuscript artifacts remain closed.",
        "",
        f"Checks passed: {passed}/{len(checks)}.",
        f"Units with summaries: {len(summaries)}/{len(UNITS)}.",
        "",
        "A complete pass requires every check to be TRUE and every unit summary to be present.",
    ]
    (args.output_dir / "MV08H_EXACT500_REMAINING_UNIT_VALIDATION.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"MV8-H exact500 remaining validation checks={passed}/{len(checks)} units={len(summaries)}/{len(UNITS)}", flush=True)
    return 0 if passed == len(checks) and len(summaries) == len(UNITS) else 2


if __name__ == "__main__":
    raise SystemExit(main())
