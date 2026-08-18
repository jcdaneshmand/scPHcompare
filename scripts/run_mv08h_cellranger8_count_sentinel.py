#!/usr/bin/env python3
"""Fail-closed launcher for the frozen MV8-H count sentinel.

The default and currently authorized mode is ``--dry-run``. Actual execution
also requires an explicit authorization token in a later owner-authorized
sprint. The launcher never kills a process and never deletes an artifact.

Non-destructive sentinel enforcement defaults are enforced with
``--localcores=4`` and ``--localmem=32``.
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


UNIT = "HCA_BM_002"
SAMPLE = "MantonBM2_HiSeq_9"
RUN_ID = "mv08h_count_sentinel_hca_bm_002"
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
FASTQS = (
    ("MantonBM2_HiSeq_9_S2_L002_I1_001.fastq.gz", 418_604_700,
     "644ebfcb0e60e81cd2f6912582cf39a82dc80fac615be643c8d7cb2296c088eb"),
    ("MantonBM2_HiSeq_9_S2_L002_R1_001.fastq.gz", 1_299_741_789,
     "c25a5807ceb616073c862253f90ee35e2b913bbb9ef81d7e3729b61bd8e881ab"),
    ("MantonBM2_HiSeq_9_S2_L002_R2_001.fastq.gz", 3_886_709_372,
     "4c223a1fb9d6ef44c166d7cf22223982287ebbdd3e55fb2cf33994796646051c"),
    ("MantonBM2_HiSeq_9_S2_L003_I1_001.fastq.gz", 419_885_066,
     "0dcea8a7b5c40c42fe77ed2ab9a086f1ce40456169d3b1cce6f3247de2b3514f"),
    ("MantonBM2_HiSeq_9_S2_L003_R1_001.fastq.gz", 1_306_937_357,
     "3560974c562e79b3ba2fbddf4e553f32d964383793911d60fbe38f1583fad020"),
    ("MantonBM2_HiSeq_9_S2_L003_R2_001.fastq.gz", 3_917_745_348,
     "bc867c9574a2d6fc98f79caba2d97ed5f7b19e2ca7136fb20bd27103b305f2cc"),
)
LOCAL_CORES = 4
LOCAL_MEMORY_GIB = 32
LOCAL_CORES_COMMAND = "--localcores=4"
LOCAL_MEMORY_COMMAND = "--localmem=32"
RSS_CAP_BYTES = 80 * 1024**3
WORKSPACE_CAP_BYTES = 200 * 1024**3
ELAPSED_CAP_SECONDS = 96 * 60 * 60
FREE_SPACE_FLOOR_BYTES = 1024**4
MONITOR_INTERVAL_SECONDS = 30
EXECUTION_TOKEN = "MV08H_COUNT_SENTINEL_EXECUTION_AUTHORIZED"


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


def public_command() -> str:
    return (
        f"cellranger count --id={RUN_ID} "
        "--transcriptome=<verified_custom_reference> "
        "--fastqs=<verified_hca_bm_002_fastq_directory> "
        f"--sample={SAMPLE} --chemistry=SC3Pv2 --expect-cells=7000 "
        "--include-introns=false --create-bam=false --nosecondary "
        f"--localcores={LOCAL_CORES} --localmem={LOCAL_MEMORY_GIB} --disable-ui"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cellranger-root", type=Path, required=True)
    parser.add_argument("--fastq-cache", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--execution-token")
    args = parser.parse_args()

    if not args.dry_run and args.execution_token != EXECUTION_TOKEN:
        raise RuntimeError(
            "count execution is not authorized; use --dry-run or obtain the exact later authorization"
        )
    launcher = args.cellranger_root.resolve() / "bin" / "cellranger"
    fastq_dir = args.fastq_cache.resolve() / UNIT
    reference = args.reference_dir.resolve()
    run_root = args.run_root.resolve()
    target = run_root / RUN_ID
    if target.exists():
        raise RuntimeError("sentinel run target already exists; refusing overwrite or resume")
    if not exact_file(launcher, EXPECTED_LAUNCHER_BYTES, EXPECTED_LAUNCHER_SHA256):
        raise RuntimeError("Cell Ranger launcher identity differs from the freeze")
    version = subprocess.run(
        [str(launcher), "--version"], check=True, text=True, capture_output=True
    ).stdout.strip()
    if version != EXPECTED_VERSION:
        raise RuntimeError(f"Cell Ranger version differs: {version!r}")
    for name, expected_bytes, expected_sha in FASTQS:
        if not exact_file(fastq_dir / name, expected_bytes, expected_sha):
            raise RuntimeError(f"sentinel FASTQ identity differs: {name}")
    observed_fastqs = sorted(path.name for path in fastq_dir.glob("*.fastq.gz"))
    if observed_fastqs != sorted(row[0] for row in FASTQS):
        raise RuntimeError("sentinel FASTQ directory has missing or extra FASTQ files")
    if list(fastq_dir.glob("*.partial")):
        raise RuntimeError("sentinel FASTQ directory contains a partial file")
    regular, symlinks, reference_bytes, reference_sha = tree_identity(reference)
    if (regular, symlinks, reference_bytes, reference_sha) != (
        EXPECTED_REFERENCE_FILES, 0, EXPECTED_REFERENCE_BYTES,
        EXPECTED_REFERENCE_TREE_SHA256,
    ):
        raise RuntimeError("reference tree identity differs from the freeze")
    if not run_root.is_dir():
        raise RuntimeError("run root must already exist")
    initial_free = shutil.disk_usage(run_root).free
    if initial_free < FREE_SPACE_FLOOR_BYTES:
        raise RuntimeError("1-TiB free-space floor is not met")
    competing = competing_cellranger_processes()
    if competing:
        summary = ", ".join(f"pid={pid}:{command}" for pid, command in competing)
        raise RuntimeError("another Cell Ranger process is active; refusing competition: " + summary)

    command = [
        str(launcher), "count", f"--id={RUN_ID}",
        f"--transcriptome={reference}", f"--fastqs={fastq_dir}",
        f"--sample={SAMPLE}", "--chemistry=SC3Pv2", "--expect-cells=7000",
        "--include-introns=false", "--create-bam=false", "--nosecondary",
        LOCAL_CORES_COMMAND, LOCAL_MEMORY_COMMAND,
        "--disable-ui",
    ]
    if args.dry_run:
        print(public_command())
        print(f"sentinel_fastq_files={len(FASTQS)}")
        print(f"sentinel_fastq_bytes={sum(row[1] for row in FASTQS)}")
        print(f"reference_tree_sha256={reference_sha}")
        print(f"initial_free_bytes={initial_free}")
        print("count_executed=false")
        return 0

    stdout_path = run_root / "mv08h-count-sentinel.stdout.log"
    stderr_path = run_root / "mv08h-count-sentinel.stderr.log"
    samples_path = run_root / "mv08h-count-sentinel-resource-samples.csv"
    receipt_path = run_root / "mv08h-count-sentinel-private-receipt.csv"
    for path in (stdout_path, stderr_path, samples_path, receipt_path):
        if path.exists():
            raise RuntimeError(f"private run artifact already exists: {path.name}")
    started = utc_now()
    start_monotonic = time.monotonic()
    max_rss_kib = max_tree_bytes = 0
    min_free_bytes = initial_free
    resource_breach = False

    with (
        stdout_path.open("wb") as stdout_handle,
        stderr_path.open("wb") as stderr_handle,
        samples_path.open("w", encoding="utf-8", newline="") as samples_handle,
    ):
        fields = [
            "sample_time_utc", "elapsed_seconds", "processes", "rss_kib",
            "run_tree_bytes", "free_bytes", "rss_cap_passed",
            "workspace_cap_passed", "elapsed_cap_passed",
            "free_space_floor_passed",
        ]
        writer = csv.DictWriter(samples_handle, fieldnames=fields, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        samples_handle.flush()
        process = subprocess.Popen(
            command, cwd=run_root, stdout=stdout_handle, stderr=stderr_handle,
            start_new_session=True,
        )
        while True:
            elapsed = int(time.monotonic() - start_monotonic)
            try:
                processes, rss_kib = process_tree_rss_kib(process.pid)
            except subprocess.CalledProcessError:
                processes, rss_kib = 0, 0
            current_tree = tree_bytes(target)
            current_free = shutil.disk_usage(run_root).free
            max_rss_kib = max(max_rss_kib, rss_kib)
            max_tree_bytes = max(max_tree_bytes, current_tree)
            min_free_bytes = min(min_free_bytes, current_free)
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
                "run_tree_bytes": current_tree, "free_bytes": current_free,
                **gates,
            })
            samples_handle.flush()
            print(
                f"count monitor elapsed={elapsed}s processes={processes} "
                f"rss_gib={rss_kib / 1024**2:.2f} "
                f"tree_gib={current_tree / 1024**3:.2f} "
                f"free_tib={current_free / 1024**4:.3f} "
                f"resource_breach={str(resource_breach).lower()}", flush=True,
            )
            return_code = process.poll()
            if return_code is not None:
                break
            time.sleep(MONITOR_INTERVAL_SECONDS)

    elapsed = int(time.monotonic() - start_monotonic)
    receipt = {
        "contract_id": "mv08h_count_sentinel_private_receipt_v1",
        "started_utc": started, "finished_utc": utc_now(),
        "elapsed_seconds": elapsed, "public_command": public_command(),
        "return_code": return_code, "run_directory_present": target.is_dir(),
        "initial_free_bytes": initial_free, "minimum_free_bytes": min_free_bytes,
        "final_free_bytes": shutil.disk_usage(run_root).free,
        "maximum_process_tree_rss_kib": max_rss_kib,
        "maximum_run_tree_bytes": max_tree_bytes,
        "final_run_tree_bytes": tree_bytes(target),
        "resource_breach_detected": resource_breach,
        "automatic_kill_used": False, "deletion_used": False,
    }
    write_receipt(receipt_path, receipt)
    if return_code != 0 or not target.is_dir() or resource_breach:
        raise RuntimeError(
            "count sentinel did not close cleanly; private logs and artifacts were preserved"
        )
    print(f"MV8-H count sentinel completed elapsed_seconds={elapsed}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
