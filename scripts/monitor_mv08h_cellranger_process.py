#!/usr/bin/env python3
"""Read-only monitor for an already-running MV08-H Cell Ranger process.

This utility never signals, pauses, resumes, or deletes a process.  It is
used only when the parent runner is lost while the Cell Ranger child remains
active, so the scientific run can finish with auditable resource evidence.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
import subprocess
import time


RSS_CAP_BYTES = 80 * 1024**3
WORKSPACE_CAP_BYTES = 200 * 1024**3
FREE_SPACE_FLOOR_BYTES = 1024**4


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def process_tree(root_pid: int) -> tuple[int, int, bool]:
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
    present = root_pid in rows
    descendants: set[int] = set()
    pending = [root_pid]
    while pending:
        pid = pending.pop()
        if pid in descendants:
            continue
        descendants.add(pid)
        pending.extend(children.get(pid, []))
    return len(descendants), sum(rows.get(pid, (0, 0))[1] for pid in descendants), present


def tree_bytes(root: Path) -> int:
    if not root.exists():
        return 0
    return sum(
        path.stat().st_size for path in root.rglob("*")
        if path.is_file() and not path.is_symlink()
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pid", type=int, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--interval", type=int, default=30)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample_time_utc", "root_pid", "root_present", "processes", "rss_kib",
        "run_tree_bytes", "free_bytes", "rss_cap_passed",
        "workspace_cap_passed", "free_space_floor_passed",
    ]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        handle.flush()
        while True:
            processes, rss_kib, present = process_tree(args.pid)
            run_tree = tree_bytes(args.workspace)
            free_bytes = shutil.disk_usage(args.workspace).free
            row = {
                "sample_time_utc": utc_now(), "root_pid": args.pid,
                "root_present": present, "processes": processes,
                "rss_kib": rss_kib, "run_tree_bytes": run_tree,
                "free_bytes": free_bytes,
                "rss_cap_passed": rss_kib * 1024 <= RSS_CAP_BYTES,
                "workspace_cap_passed": run_tree <= WORKSPACE_CAP_BYTES,
                "free_space_floor_passed": free_bytes >= FREE_SPACE_FLOOR_BYTES,
            }
            writer.writerow(row)
            handle.flush()
            print(
                f"monitor pid={args.pid} present={str(present).lower()} "
                f"processes={processes} rss_gib={rss_kib / 1024**2:.2f} "
                f"tree_gib={run_tree / 1024**3:.2f} "
                f"free_tib={free_bytes / 1024**4:.3f}", flush=True,
            )
            if not present:
                return 0
            time.sleep(args.interval)


if __name__ == "__main__":
    raise SystemExit(main())
