#!/usr/bin/env python3
"""Resumable exact all-level landscape chunks for MV5-C2 query-training pairs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import resource
import time
from collections import defaultdict

import numpy as np
from persim import PersLandscapeExact


ENGINE_ID = "persim_0.3.8_mv05c2_exact_chunked_v1"
LANDSCAPE_CONTRACT_ID = "mv05c_full_l2_exact_all_active_levels_v1"
REQUEST_CONTRACT_ID = "mv05c2_query_training_pair_scope_v1"


def read_csv(path: str) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def atomic_write(path: str, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("Refusing to write an empty CSV artifact.")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    temporary = path + f".partial.{os.getpid()}"
    try:
        with open(temporary, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        if os.path.exists(path):
            raise FileExistsError(f"Refusing to overwrite {path}")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def file_sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def truth(value: str) -> bool:
    return value.strip().lower() in {"true", "t", "1"}


def build_landscape(values: list[tuple[float, float]], dimension: int):
    if not values:
        return None
    diagrams = [np.empty((0, 2)) for _ in range(dimension + 1)]
    diagrams[dimension] = np.asarray(values, dtype=float).reshape((-1, 2))
    return PersLandscapeExact(dgms=diagrams, hom_deg=dimension)


def integrate_critical_pairs(levels) -> tuple[float, int, int]:
    squared = 0.0
    segments = 0
    points = 0
    for level in levels:
        points += len(level)
        for (x0, y0), (x1, y1) in zip(level, level[1:]):
            width = x1 - x0
            if width < 0:
                raise ValueError("Critical-pair filtration values are unsorted.")
            squared += width * (y0 * y0 + y0 * y1 + y1 * y1) / 3.0
            segments += 1
    if squared < -1e-12:
        raise ValueError(f"Negative squared landscape distance: {squared}")
    return max(0.0, squared), segments, points


def corrected_distance(first, second):
    if first is None and second is None:
        levels = []
    elif first is None:
        levels = second.critical_pairs
    elif second is None:
        levels = first.critical_pairs
    else:
        levels = (first - second).critical_pairs
    squared, segments, points = integrate_critical_pairs(levels)
    return math.sqrt(squared), squared, len(levels), segments, points


def subset_sha256(rows: list[dict[str, str]]) -> str:
    payload = [
        {key: row[key] for key in sorted(row)}
        for row in sorted(rows, key=lambda item: item["pair_request_id"])
    ]
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def safe_chunk_name(chunk_id: str) -> str:
    return chunk_id.replace(":", "_")


def validate_requests(rows: list[dict[str, str]], max_pairs: int) -> None:
    required = {
        "contract_id", "pair_request_id", "source_pair_id", "chunk_id",
        "chunk_offset", "stratum_id", "fold_id", "seed", "representation",
        "view_id", "homology_dimension", "first_sample_id",
        "second_sample_id", "first_diagram_id", "second_diagram_id",
        "pair_scope", "exact", "all_active_levels", "outcome_label_state",
        "biological_outcomes_computed",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Pair request chunk has an invalid schema.")
    if len(rows) > max_pairs:
        raise ValueError("Pair request chunk exceeds the resource guard.")
    if len({row["pair_request_id"] for row in rows}) != len(rows):
        raise ValueError("Pair request IDs are duplicated within a chunk.")
    for row in rows:
        if (
            row["contract_id"] != REQUEST_CONTRACT_ID
            or row["pair_scope"] != "held_out_query_to_training_reference"
            or row["outcome_label_state"] != "closed"
            or truth(row["biological_outcomes_computed"])
            or not truth(row["exact"])
            or not truth(row["all_active_levels"])
            or row["homology_dimension"] not in {"H0", "H1"}
        ):
            raise ValueError("Pair request violates the MV5-C2 boundary.")


def validate_resume(
    output_path: str, status_path: str, rows: list[dict[str, str]]
) -> bool:
    output_exists = os.path.exists(output_path)
    status_exists = os.path.exists(status_path)
    if not output_exists and not status_exists:
        return False
    if output_exists != status_exists:
        raise RuntimeError("Chunk has an incomplete output/status pair.")
    status = read_csv(status_path)
    output = read_csv(output_path)
    expected_ids = sorted(row["pair_request_id"] for row in rows)
    observed_ids = sorted(row["pair_request_id"] for row in output)
    if (
        len(status) != 1
        or status[0]["status"] != "completed"
        or status[0]["request_subset_sha256"] != subset_sha256(rows)
        or status[0]["output_sha256"] != file_sha256(output_path)
        or expected_ids != observed_ids
    ):
        raise RuntimeError("Existing chunk output failed resume validation.")
    return True


def execute_chunk(
    chunk_id: str,
    rows: list[dict[str, str]],
    intervals: dict[tuple[str, int], list[tuple[float, float]]],
    landscapes: dict,
    output_dir: str,
    status_dir: str,
    max_pairs: int,
    max_seconds: float,
) -> str:
    validate_requests(rows, max_pairs)
    safe_name = safe_chunk_name(chunk_id)
    output_path = os.path.join(output_dir, safe_name + ".csv")
    status_path = os.path.join(status_dir, safe_name + "-status.csv")
    if validate_resume(output_path, status_path, rows):
        return "reused"
    started = time.perf_counter()
    output_rows = []
    for row in sorted(rows, key=lambda item: int(item["chunk_offset"])):
        if time.perf_counter() - started > max_seconds:
            raise TimeoutError("Chunk exceeded max_seconds before completion.")
        dimension = int(row["homology_dimension"][1:])
        first_key = (row["first_diagram_id"], dimension)
        second_key = (row["second_diagram_id"], dimension)
        for key in (first_key, second_key):
            if key not in landscapes:
                landscapes[key] = build_landscape(intervals[key], dimension)
        pair_started = time.perf_counter()
        distance, squared, levels, segments, points = corrected_distance(
            landscapes[first_key], landscapes[second_key]
        )
        output_rows.append({
            "engine_id": ENGINE_ID,
            "contract_id": LANDSCAPE_CONTRACT_ID,
            "pair_request_id": row["pair_request_id"],
            "source_pair_id": row["source_pair_id"],
            "chunk_id": chunk_id,
            "stratum_id": row["stratum_id"],
            "fold_id": row["fold_id"],
            "seed": row["seed"],
            "representation": row["representation"],
            "view_id": row["view_id"],
            "homology_dimension": row["homology_dimension"],
            "first_sample_id": row["first_sample_id"],
            "second_sample_id": row["second_sample_id"],
            "first_diagram_id": row["first_diagram_id"],
            "second_diagram_id": row["second_diagram_id"],
            "difference_active_levels": levels,
            "difference_segments": segments,
            "difference_critical_points": points,
            "distance": distance,
            "squared_distance": squared,
            "pair_seconds": time.perf_counter() - pair_started,
            "exact": True,
            "all_active_levels": True,
            "absolute_error_estimate": 0.0,
            "status": "completed",
            "outcome_label_state": "closed",
            "biological_outcomes_computed": False,
        })
    atomic_write(output_path, output_rows)
    elapsed = time.perf_counter() - started
    status_rows = [{
        "contract_id": "mv05c2_pair_chunk_status_v1",
        "engine_id": ENGINE_ID,
        "chunk_id": chunk_id,
        "request_count": len(rows),
        "completed_count": len(output_rows),
        "request_subset_sha256": subset_sha256(rows),
        "output_file": os.path.basename(output_path),
        "output_sha256": file_sha256(output_path),
        "elapsed_seconds": elapsed,
        "peak_process_rss_bytes": resource.getrusage(
            resource.RUSAGE_SELF
        ).ru_maxrss * 1024,
        "max_pairs_guard": max_pairs,
        "max_seconds_guard": max_seconds,
        "status": "completed",
        "outcome_label_state": "closed",
        "biological_outcomes_computed": False,
    }]
    atomic_write(status_path, status_rows)
    return "completed"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--requests", required=True)
    parser.add_argument("--intervals", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--status-dir", required=True)
    parser.add_argument("--max-pairs", type=int, default=250)
    parser.add_argument("--max-seconds", type=float, default=1800.0)
    args = parser.parse_args()
    if args.max_pairs < 1 or args.max_seconds <= 0:
        raise ValueError("Resource guards must be positive.")

    requests = read_csv(args.requests)
    groups = defaultdict(list)
    for row in requests:
        groups[row["chunk_id"]].append(row)
    if not groups:
        raise ValueError("Pair request manifest is empty.")
    interval_rows = read_csv(args.intervals)
    intervals = defaultdict(list)
    for row in interval_rows:
        intervals[(row["diagram_id"], int(row["homology_dimension"]))].append(
            (float(row["birth"]), float(row["death"]))
        )

    completed = 0
    reused = 0
    ordered_groups = sorted(
        groups.items(),
        key=lambda item: (
            item[1][0]["stratum_id"],
            item[1][0]["homology_dimension"],
            int(item[1][0]["chunk_index"]),
        ),
    )
    landscape_cache = {}
    active_locality = None
    for chunk_id, chunk_rows in ordered_groups:
        locality = (
            chunk_rows[0]["stratum_id"],
            chunk_rows[0]["homology_dimension"],
        )
        if locality != active_locality:
            landscape_cache.clear()
            active_locality = locality
        disposition = execute_chunk(
            chunk_id, chunk_rows, intervals, landscape_cache, args.output_dir,
            args.status_dir, args.max_pairs, args.max_seconds
        )
        completed += disposition == "completed"
        reused += disposition == "reused"
    print(
        f"MV5-C2 chunks complete: built={completed} reused={reused} "
        f"total={len(groups)}"
    )


if __name__ == "__main__":
    main()
