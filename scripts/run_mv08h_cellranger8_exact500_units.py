#!/usr/bin/env python3
"""Run the remaining MV8-H HCA units under the frozen exact-500 contract.

This runner is intentionally serial and fail-closed.  It never kills a
process, deletes an artifact, overwrites a target, resumes a partial run, or
starts a second Cell Ranger process.  Private logs, resource samples, and
receipts remain under ``tmp``; only identity/status summaries should be
published by a later independent validation step.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import time


UNITS = (
    "HCA_BM_001", "HCA_BM_003", "HCA_BM_004", "HCA_BM_005",
    "HCA_BM_006", "HCA_BM_007", "HCA_BM_008",
)
EXPECTED_VERSION = "cellranger cellranger-8.0.1"
EXPECTED_LAUNCHER_BYTES = 19_924_984
EXPECTED_LAUNCHER_SHA256 = (
    "4ee3a1670b4f14c826004fe8e17b4759e1edc701b15ff2e9623753bf1b34d4d6"
)
EXPECTED_REFERENCE_FILES = 19
EXPECTED_REFERENCE_BYTES = 20_765_871_518
EXPECTED_REFERENCE_TREE_SHA256 = (
    "5e2aff9e7154e6b02f98552a4419bd48edce66e617e579ae562e714f79199f1c"
)
LOCAL_CORES = 4
LOCAL_MEMORY_GIB = 32
RSS_CAP_BYTES = 80 * 1024**3
WORKSPACE_CAP_BYTES = 200 * 1024**3
ELAPSED_CAP_SECONDS = 96 * 60 * 60
FREE_SPACE_FLOOR_BYTES = 1024**4
MONITOR_INTERVAL_SECONDS = 30
EXECUTION_TOKEN = "MV08H_EXACT500_REMAINING_UNITS_AUTHORIZED"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_file(path: Path, expected_bytes: int, expected_sha: str) -> bool:
    return (
        path.is_file()
        and path.stat().st_size == expected_bytes
        and sha256_file(path) == expected_sha
    )


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


def tree_bytes(root: Path) -> int:
    if not root.exists():
        return 0
    return sum(
        path.stat().st_size for path in root.rglob("*")
        if path.is_file() and not path.is_symlink()
    )


def process_tree_rss_kib(root_pid: int) -> tuple[int, int]:
    result = subprocess.run(
        ["ps", "-e", "-o", "pid=,ppid=,rss="],
        check=True, text=True, capture_output=True,
    )
    rows: dict[int, tuple[int, int]] = {}
    children: dict[int, list[int]] = {}
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) != 3:
            continue
        pid, parent, rss = map(int, fields)
        rows[pid] = (parent, rss)
        children.setdefault(parent, []).append(pid)
    descendants: set[int] = set()
    pending = [root_pid]
    while pending:
        pid = pending.pop()
        if pid in descendants:
            continue
        descendants.add(pid)
        pending.extend(children.get(pid, []))
    return len(descendants), sum(rows.get(pid, (0, 0))[1] for pid in descendants)


def competing_cellranger_processes() -> list[tuple[int, str]]:
    result = subprocess.run(
        ["ps", "-e", "-o", "pid=,comm="],
        check=True, text=True, capture_output=True,
    )
    competing: list[tuple[int, str]] = []
    for line in result.stdout.splitlines():
        fields = line.split(maxsplit=1)
        if len(fields) != 2:
            continue
        pid, command = int(fields[0]), fields[1]
        if command.lower() in {"cellranger", "martian", "star"}:
            competing.append((pid, command))
    return competing


def write_receipt(path: Path, row: dict[str, object]) -> None:
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerow(row)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def load_manifest(path: Path) -> dict[str, list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    result: dict[str, list[dict[str, str]]] = {unit: [] for unit in UNITS}
    for row in rows:
        unit = row.get("unit_id", "")
        if unit in result:
            result[unit].append(row)
    for unit in UNITS:
        if len(result[unit]) != 6:
            raise RuntimeError(f"{unit} does not have exactly six frozen FASTQ rows")
    return result


def unit_number(unit: str) -> int:
    return int(unit.rsplit("_", 1)[1])


def validate_inputs(
    launcher: Path, fastq_cache: Path, reference: Path, manifest: Path,
    run_root: Path,
) -> dict[str, list[dict[str, str]]]:
    if not exact_file(launcher, EXPECTED_LAUNCHER_BYTES, EXPECTED_LAUNCHER_SHA256):
        raise RuntimeError("Cell Ranger launcher identity differs from the freeze")
    version = subprocess.run(
        [str(launcher), "--version"], check=True, text=True, capture_output=True,
    ).stdout.strip()
    if version != EXPECTED_VERSION:
        raise RuntimeError(f"Cell Ranger version differs: {version!r}")
    regular, symlinks, reference_bytes, reference_sha = tree_identity(reference)
    if (regular, symlinks, reference_bytes, reference_sha) != (
        EXPECTED_REFERENCE_FILES, 0, EXPECTED_REFERENCE_BYTES,
        EXPECTED_REFERENCE_TREE_SHA256,
    ):
        raise RuntimeError("reference tree identity differs from the freeze")
    rows_by_unit = load_manifest(manifest)
    for unit, rows in rows_by_unit.items():
        fastq_dir = fastq_cache / unit
        observed = sorted(path.name for path in fastq_dir.glob("*.fastq.gz"))
        expected = sorted(row["file_name"] for row in rows)
        if observed != expected:
            raise RuntimeError(f"{unit} FASTQ directory has missing or extra files")
        if list(fastq_dir.glob("*.partial")):
            raise RuntimeError(f"{unit} FASTQ directory contains a partial file")
        for row in rows:
            path = fastq_dir / row["file_name"]
            if not exact_file(path, int(row["file_size_bytes"]), row["file_sha256"]):
                raise RuntimeError(f"{unit} FASTQ identity differs: {row['file_name']}")
    run_root.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(run_root).free < FREE_SPACE_FLOOR_BYTES:
        raise RuntimeError("1-TiB free-space floor is not met")
    competing = competing_cellranger_processes()
    if competing:
        summary = ", ".join(f"pid={pid}:{command}" for pid, command in competing)
        raise RuntimeError("another Cell Ranger process is active; refusing competition: " + summary)
    return rows_by_unit


def public_command(unit: str) -> str:
    number = unit_number(unit)
    run_id = f"mv08h_exact500_{unit.lower()}"
    sample = f"MantonBM{number}_HiSeq_9"
    return (
        f"cellranger count --id={run_id} --transcriptome=<verified_custom_reference> "
        f"--fastqs=<verified_{unit.lower()}_fastq_directory> --sample={sample} "
        "--chemistry=SC3Pv2 --expect-cells=7000 --include-introns=false "
        "--create-bam=false --nosecondary --localcores=4 --localmem=32 --disable-ui"
    )


def run_unit(
    unit: str, rows: list[dict[str, str]], launcher: Path, fastq_cache: Path,
    reference: Path, run_root: Path,
) -> dict[str, object]:
    number = unit_number(unit)
    run_id = f"mv08h_exact500_{unit.lower()}"
    unit_root = run_root / unit
    target = unit_root / run_id
    if unit_root.exists() or target.exists():
        raise RuntimeError(f"{unit} run target already exists; refusing overwrite or resume")
    unit_root.mkdir(parents=False)
    fastq_dir = fastq_cache / unit
    command = [
        str(launcher), "count", f"--id={run_id}",
        f"--transcriptome={reference}", f"--fastqs={fastq_dir}",
        f"--sample=MantonBM{number}_HiSeq_9", "--chemistry=SC3Pv2",
        "--expect-cells=7000", "--include-introns=false", "--create-bam=false",
        "--nosecondary", "--localcores=4", "--localmem=32", "--disable-ui",
    ]
    stdout_path = unit_root / "stdout.log"
    stderr_path = unit_root / "stderr.log"
    samples_path = unit_root / "resource-samples.csv"
    receipt_path = unit_root / "private-receipt.csv"
    started = utc_now()
    start_monotonic = time.monotonic()
    initial_free = shutil.disk_usage(unit_root).free
    max_rss_kib = max_tree = 0
    min_free = initial_free
    resource_breach = False
    with (
        stdout_path.open("wb") as stdout_handle,
        stderr_path.open("wb") as stderr_handle,
        samples_path.open("w", encoding="utf-8", newline="") as samples_handle,
    ):
        fields = [
            "sample_time_utc", "elapsed_seconds", "processes", "rss_kib",
            "run_tree_bytes", "free_bytes", "rss_cap_passed",
            "workspace_cap_passed", "elapsed_cap_passed", "free_space_floor_passed",
        ]
        writer = csv.DictWriter(samples_handle, fieldnames=fields, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        samples_handle.flush()
        process = subprocess.Popen(
            command, cwd=unit_root, stdout=stdout_handle, stderr=stderr_handle,
            start_new_session=True,
        )
        while True:
            elapsed = int(time.monotonic() - start_monotonic)
            try:
                processes, rss_kib = process_tree_rss_kib(process.pid)
            except subprocess.CalledProcessError:
                processes, rss_kib = 0, 0
            current_tree = tree_bytes(target)
            current_free = shutil.disk_usage(unit_root).free
            max_rss_kib = max(max_rss_kib, rss_kib)
            max_tree = max(max_tree, current_tree)
            min_free = min(min_free, current_free)
            gates = {
                "rss_cap_passed": rss_kib * 1024 <= RSS_CAP_BYTES,
                "workspace_cap_passed": current_tree <= WORKSPACE_CAP_BYTES,
                "elapsed_cap_passed": elapsed <= ELAPSED_CAP_SECONDS,
                "free_space_floor_passed": current_free >= FREE_SPACE_FLOOR_BYTES,
            }
            resource_breach = resource_breach or not all(gates.values())
            writer.writerow({
                "sample_time_utc": utc_now(), "elapsed_seconds": elapsed,
                "processes": processes, "rss_kib": rss_kib,
                "run_tree_bytes": current_tree, "free_bytes": current_free, **gates,
            })
            samples_handle.flush()
            print(
                f"{unit} monitor elapsed={elapsed}s processes={processes} "
                f"rss_gib={rss_kib / 1024**2:.2f} tree_gib={current_tree / 1024**3:.2f} "
                f"free_tib={current_free / 1024**4:.3f} "
                f"resource_breach={str(resource_breach).lower()}", flush=True,
            )
            return_code = process.poll()
            if return_code is not None:
                break
            time.sleep(MONITOR_INTERVAL_SECONDS)
    elapsed = int(time.monotonic() - start_monotonic)
    receipt = {
        "contract_id": "mv08h_exact500_private_receipt_v1",
        "unit_id": unit, "started_utc": started, "finished_utc": utc_now(),
        "elapsed_seconds": elapsed, "public_command": public_command(unit),
        "return_code": return_code, "run_directory_present": target.is_dir(),
        "initial_free_bytes": initial_free, "minimum_free_bytes": min_free,
        "final_free_bytes": shutil.disk_usage(unit_root).free,
        "maximum_process_tree_rss_kib": max_rss_kib,
        "maximum_run_tree_bytes": max_tree, "final_run_tree_bytes": tree_bytes(target),
        "resource_breach_detected": resource_breach,
        "automatic_kill_used": False, "deletion_used": False,
        "fastq_count": len(rows), "fastq_bytes": sum(int(row["file_size_bytes"]) for row in rows),
    }
    write_receipt(receipt_path, receipt)
    if return_code != 0 or not target.is_dir() or resource_breach:
        raise RuntimeError(f"{unit} did not close cleanly; private logs and artifacts were preserved")
    print(f"{unit} completed elapsed_seconds={elapsed}", flush=True)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cellranger-root", type=Path, required=True)
    parser.add_argument("--fastq-cache", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument(
        "--start-after", choices=UNITS, default=None,
        help="skip units through this already-closed unit; never resumes its target",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--execution-token")
    args = parser.parse_args()
    if not args.dry_run and args.execution_token != EXECUTION_TOKEN:
        raise RuntimeError("remaining-unit count execution is not authorized")
    launcher = args.cellranger_root.resolve() / "bin" / "cellranger"
    fastq_cache = args.fastq_cache.resolve()
    reference = args.reference_dir.resolve()
    manifest = args.manifest.resolve()
    run_root = args.run_root.resolve()
    rows_by_unit = validate_inputs(launcher, fastq_cache, reference, manifest, run_root)
    start_index = 0
    if args.start_after is not None:
        start_index = UNITS.index(args.start_after) + 1
    if args.dry_run:
        for unit in UNITS[start_index:]:
            rows = rows_by_unit[unit]
            print(public_command(unit))
            print(f"{unit} fastq_files={len(rows)} fastq_bytes={sum(int(row['file_size_bytes']) for row in rows)}")
        print(f"reference_tree_sha256={EXPECTED_REFERENCE_TREE_SHA256}")
        print(f"initial_free_bytes={shutil.disk_usage(run_root).free}")
        print("count_executed=false")
        return 0
    receipts: list[dict[str, object]] = []
    for unit in UNITS[start_index:]:
        receipts.append(run_unit(unit, rows_by_unit[unit], launcher, fastq_cache, reference, run_root))
    print(f"MV8-H exact500 remaining units completed count={len(receipts)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
