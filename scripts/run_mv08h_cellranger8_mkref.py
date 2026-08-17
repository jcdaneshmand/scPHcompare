#!/usr/bin/env python3
"""Run the prospectively frozen MV8-H Cell Ranger 8.0.1 reference build.

The launcher fails closed before execution when an input, runtime, target, or
storage identity differs. During execution it records only private, Git-ignored
logs and resource samples. It never kills a process or deletes an artifact.
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


GENOME = "GRCh38_Ensembl93_target_complete_33563"
EXPECTED_VERSION = "cellranger cellranger-8.0.1"
EXPECTED_LAUNCHER_BYTES = 19_924_984
EXPECTED_LAUNCHER_SHA256 = (
    "4ee3a1670b4f14c826004fe8e17b4759e1edc701b15ff2e9623753bf1b34d4d6"
)
EXPECTED_FASTA_BYTES = 3_151_425_857
EXPECTED_FASTA_SHA256 = (
    "78777b0886e8dfa5e14e4957fbbaa53736fcbaa5668d59e09b6b7945fca93d8c"
)
EXPECTED_GTF_BYTES = 1_099_737_654
EXPECTED_GTF_SHA256 = (
    "e28e4c4faf0dd76884d5e94c481fce2db43ad303968067c1276092a234727182"
)
NTHREADS = 4
MEM_GB = 32
REFERENCE_CAP_BYTES = 50 * 1024**3
FREE_SPACE_FLOOR_BYTES = 1024**4
MONITOR_INTERVAL_SECONDS = 30


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_file(path: Path, expected_bytes: int, expected_sha256: str) -> bool:
    return (
        path.is_file()
        and path.stat().st_size == expected_bytes
        and sha256_file(path) == expected_sha256
    )


def tree_bytes(root: Path) -> int:
    if not root.exists():
        return 0
    total = 0
    for directory, _, filenames in os.walk(root, followlinks=False):
        base = Path(directory)
        for name in filenames:
            path = base / name
            if path.is_file() and not path.is_symlink():
                total += path.stat().st_size
    return total


def process_tree_rss_kib(root_pid: int) -> tuple[int, int]:
    result = subprocess.run(
        ["ps", "-e", "-o", "pid=,ppid=,rss="],
        check=True,
        text=True,
        capture_output=True,
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
    rss = sum(rows.get(pid, (0, 0))[1] for pid in descendants)
    return len(descendants), rss


def write_receipt(path: Path, row: dict[str, object]) -> None:
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerow(row)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cellranger-root", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    launcher = args.cellranger_root.resolve() / "bin" / "cellranger"
    fasta = args.fasta.resolve()
    genes = args.genes.resolve()
    run_root = args.run_root.resolve()
    reference = run_root / GENOME

    if run_root.exists() or reference.exists():
        raise RuntimeError("run root or reference target already exists; refusing overwrite")
    if not exact_file(launcher, EXPECTED_LAUNCHER_BYTES, EXPECTED_LAUNCHER_SHA256):
        raise RuntimeError("Cell Ranger launcher identity differs from the freeze")
    version = subprocess.run(
        [str(launcher), "--version"], check=True, text=True, capture_output=True
    ).stdout.strip()
    if version != EXPECTED_VERSION:
        raise RuntimeError(f"Cell Ranger version differs: {version!r}")
    if not exact_file(fasta, EXPECTED_FASTA_BYTES, EXPECTED_FASTA_SHA256):
        raise RuntimeError("uncompressed FASTA identity differs from the freeze")
    if not exact_file(genes, EXPECTED_GTF_BYTES, EXPECTED_GTF_SHA256):
        raise RuntimeError("target-complete GTF identity differs from the freeze")
    initial_free = shutil.disk_usage(run_root.parent).free
    if initial_free < FREE_SPACE_FLOOR_BYTES:
        raise RuntimeError("free-space floor is not met")

    command = [
        str(launcher), "mkref", f"--genome={GENOME}",
        f"--fasta={fasta}", f"--genes={genes}",
        f"--nthreads={NTHREADS}", f"--memgb={MEM_GB}",
        f"--localcores={NTHREADS}", f"--localmem={MEM_GB}", "--disable-ui",
    ]
    public_command = (
        f"cellranger mkref --genome={GENOME} "
        "--fasta=<verified_uncompressed_ensembl93_primary_assembly> "
        "--genes=<verified_target_complete_33563_gtf> "
        f"--nthreads={NTHREADS} --memgb={MEM_GB} "
        f"--localcores={NTHREADS} --localmem={MEM_GB} --disable-ui"
    )
    if args.dry_run:
        print(public_command)
        print(f"initial_free_bytes={initial_free}")
        return 0

    run_root.mkdir(parents=True, exist_ok=False)
    monitor = run_root / "mv08h-mkref-resource-samples.csv"
    stdout_path = run_root / "mv08h-mkref.stdout.log"
    stderr_path = run_root / "mv08h-mkref.stderr.log"
    receipt_path = run_root / "mv08h-mkref-private-receipt.csv"
    started = utc_now()
    start_monotonic = time.monotonic()
    max_rss_kib = 0
    max_tree_bytes = 0
    min_free_bytes = initial_free
    resource_breach = False

    with (
        stdout_path.open("wb") as stdout_handle,
        stderr_path.open("wb") as stderr_handle,
        monitor.open("w", encoding="utf-8", newline="") as monitor_handle,
    ):
        fields = [
            "sample_time_utc", "elapsed_seconds", "processes", "rss_kib",
            "run_tree_bytes", "free_bytes", "reference_cap_passed",
            "free_space_floor_passed",
        ]
        writer = csv.DictWriter(monitor_handle, fieldnames=fields, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        monitor_handle.flush()
        process = subprocess.Popen(
            command,
            cwd=run_root,
            stdout=stdout_handle,
            stderr=stderr_handle,
            start_new_session=True,
        )
        sample_index = 0
        while True:
            sample_index += 1
            elapsed = int(time.monotonic() - start_monotonic)
            try:
                processes, rss_kib = process_tree_rss_kib(process.pid)
            except subprocess.CalledProcessError:
                processes, rss_kib = 0, 0
            current_tree = tree_bytes(run_root)
            current_free = shutil.disk_usage(run_root).free
            max_rss_kib = max(max_rss_kib, rss_kib)
            max_tree_bytes = max(max_tree_bytes, current_tree)
            min_free_bytes = min(min_free_bytes, current_free)
            cap_passed = current_tree <= REFERENCE_CAP_BYTES
            floor_passed = current_free >= FREE_SPACE_FLOOR_BYTES
            resource_breach = resource_breach or not cap_passed or not floor_passed
            writer.writerow({
                "sample_time_utc": utc_now(),
                "elapsed_seconds": elapsed,
                "processes": processes,
                "rss_kib": rss_kib,
                "run_tree_bytes": current_tree,
                "free_bytes": current_free,
                "reference_cap_passed": cap_passed,
                "free_space_floor_passed": floor_passed,
            })
            monitor_handle.flush()
            if sample_index == 1 or sample_index % 2 == 0 or resource_breach:
                print(
                    f"mkref monitor elapsed={elapsed}s processes={processes} "
                    f"rss_gib={rss_kib / 1024**2:.2f} "
                    f"tree_gib={current_tree / 1024**3:.2f} "
                    f"free_tib={current_free / 1024**4:.3f} "
                    f"resource_breach={str(resource_breach).lower()}",
                    flush=True,
                )
            return_code = process.poll()
            if return_code is not None:
                break
            time.sleep(MONITOR_INTERVAL_SECONDS)

    elapsed_seconds = int(time.monotonic() - start_monotonic)
    final_tree_bytes = tree_bytes(run_root)
    final_free_bytes = shutil.disk_usage(run_root).free
    reference_exists = reference.is_dir()
    receipt = {
        "contract_id": "mv08h_cellranger8_mkref_private_receipt_v1",
        "started_utc": started,
        "finished_utc": utc_now(),
        "elapsed_seconds": elapsed_seconds,
        "public_command": public_command,
        "return_code": return_code,
        "reference_directory_present": reference_exists,
        "initial_free_bytes": initial_free,
        "minimum_free_bytes": min_free_bytes,
        "final_free_bytes": final_free_bytes,
        "maximum_process_tree_rss_kib": max_rss_kib,
        "maximum_run_tree_bytes": max_tree_bytes,
        "final_run_tree_bytes": final_tree_bytes,
        "resource_breach_detected": resource_breach,
        "automatic_kill_used": False,
        "deletion_used": False,
    }
    write_receipt(receipt_path, receipt)
    if return_code != 0 or not reference_exists or resource_breach:
        raise RuntimeError(
            "mkref did not close cleanly; private logs and partial artifacts were preserved"
        )
    print(
        f"MV8-H mkref completed elapsed_seconds={elapsed_seconds} "
        f"final_run_tree_bytes={final_tree_bytes}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
