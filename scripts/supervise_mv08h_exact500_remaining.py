#!/usr/bin/env python3
"""Wait for the orphaned BM_001 count, then continue units 003--008.

The supervisor is deliberately conservative: it never signals or resumes the
active process.  It continues only when the original Cell Ranger process has
exited, the process log has the normal successful-completion marker, required
outputs exist, and all post-parent monitor gates passed.  Otherwise it leaves
all artifacts in place and exits for human review.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import time


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def process_present(pid: int) -> bool:
    result = subprocess.run(["ps", "-p", str(pid), "-o", "pid="], text=True, capture_output=True)
    return bool(result.stdout.strip())


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_state(path: Path, state: str, evidence: str) -> None:
    partial = path.with_name(path.name + ".partial")
    partial.write_text(
        f"state={state}\nupdated_utc={utc_now()}\nevidence={evidence}\n",
        encoding="utf-8",
    )
    partial.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pid", type=int, required=True)
    parser.add_argument("--unit-root", type=Path, required=True)
    parser.add_argument("--monitor-csv", type=Path, required=True)
    parser.add_argument("--state-file", type=Path, required=True)
    parser.add_argument("--runner-cwd", type=Path, required=True)
    parser.add_argument("--runner-args", nargs=argparse.REMAINDER, required=True)
    parser.add_argument("--interval", type=int, default=60)
    args = parser.parse_args()
    args.state_file.parent.mkdir(parents=True, exist_ok=True)
    write_state(args.state_file, "waiting_for_bm001", f"pid={args.pid}")
    while process_present(args.pid):
        time.sleep(args.interval)
    stdout_path = args.unit_root / "stdout.log"
    stderr_path = args.unit_root / "stderr.log"
    target = args.unit_root / "mv08h_exact500_hca_bm_001"
    outs = target / "outs"
    required = [
        outs / "filtered_feature_bc_matrix.h5",
        outs / "raw_feature_bc_matrix.h5",
        outs / "molecule_info.h5",
        outs / "metrics_summary.csv",
    ]
    stdout = stdout_path.read_text(encoding="utf-8", errors="replace") if stdout_path.is_file() else ""
    stderr = stderr_path.read_text(encoding="utf-8", errors="replace") if stderr_path.is_file() else "missing"
    monitors = csv_rows(args.monitor_csv) if args.monitor_csv.is_file() else []
    gates_pass = bool(monitors) and all(
        row.get("rss_cap_passed") == "True"
        and row.get("workspace_cap_passed") == "True"
        and row.get("free_space_floor_passed") == "True"
        for row in monitors
    )
    success = (
        "Pipestance completed successfully!" in stdout
        and not stderr.strip()
        and all(path.is_file() for path in required)
        and gates_pass
    )
    if not success:
        write_state(
            args.state_file, "bm001_failed_or_incomplete",
            f"success_marker={('Pipestance completed successfully!' in stdout)};"
            f"stderr_empty={not stderr.strip()};required_outputs={all(path.is_file() for path in required)};"
            f"monitor_rows={len(monitors)};resource_gates_passed={gates_pass}",
        )
        return 2
    write_state(args.state_file, "bm001_closed_continue_authorized", f"monitor_rows={len(monitors)}")
    command = ["python3", "scripts/run_mv08h_cellranger8_exact500_units.py", *args.runner_args]
    result = subprocess.run(command, cwd=args.runner_cwd)
    write_state(args.state_file, "remaining_runner_finished", f"return_code={result.returncode}")
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
