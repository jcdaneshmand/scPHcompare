#!/usr/bin/env python3
"""Immutable exact all-level landscape chunks for one MV5-I fold-seed group."""

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


ENGINE_ID = "persim_0.3.8_mv05i_exact_critical_pairs_v1"
CONTRACT_ID = "mv05i_all_active_exact_landscape_distance_v1"
REQUEST_ID = "mv05i_cell_landscape_pair_v1"


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def subset_sha256(rows):
    payload = [
        {key: row[key] for key in sorted(row)}
        for row in sorted(rows, key=lambda item: item["pair_request_id"])
    ]
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def atomic_csv(path, rows):
    if not rows:
        raise ValueError("Refusing to write an empty MV5-I artifact.")
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


def safe_name(value):
    return "".join(char if char.isalnum() or char in "_.-" else "_"
                   for char in value)


def build_landscape(values, dimension):
    if not values:
        return None
    diagrams = [np.empty((0, 2)) for _ in range(dimension + 1)]
    diagrams[dimension] = np.asarray(values, dtype=float).reshape((-1, 2))
    return PersLandscapeExact(dgms=diagrams, hom_deg=dimension)


def exact_segment_integral(levels):
    squared = 0.0
    segments = 0
    points = 0
    for level in levels:
        points += len(level)
        for (x0, y0), (x1, y1) in zip(level, level[1:]):
            width = x1 - x0
            if width < 0:
                raise ValueError("Critical-pair coordinates are unsorted.")
            squared += width * (y0*y0 + y0*y1 + y1*y1) / 3.0
            segments += 1
    if squared < -1e-12:
        raise ValueError("Exact landscape squared distance is negative.")
    return max(0.0, squared), segments, points


def landscape_distance(first, second):
    if first is None and second is None:
        levels = []
    elif first is None:
        levels = second.critical_pairs
    elif second is None:
        levels = first.critical_pairs
    else:
        levels = (first - second).critical_pairs
    squared, segments, points = exact_segment_integral(levels)
    return math.sqrt(squared), squared, len(levels), segments, points


def validate_requests(rows, manifest_sha, implementation_sha, max_pairs):
    required = {
        "contract_id", "pair_request_id", "group_id", "fold_id", "seed",
        "homology_dimension", "query_sample_id", "training_sample_id",
        "query_record_cache_key", "training_record_cache_key",
        "query_diagram_sha256", "training_diagram_sha256", "chunk_id",
        "chunk_offset", "pair_scope", "landscape_definition_id",
        "essential_h0_policy", "level_policy", "integration_policy",
        "representation", "view_id",
        "pair_manifest_sha256", "outcome_label_state",
        "biological_outcomes_computed",
    }
    if not rows or not required.issubset(rows[0]) or len(rows) > max_pairs:
        raise ValueError("MV5-I request chunk has an invalid schema or size.")
    if len({row["pair_request_id"] for row in rows}) != len(rows):
        raise ValueError("MV5-I request identifiers are duplicated.")
    for row in rows:
        if (
            row["contract_id"] != REQUEST_ID
            or row["pair_manifest_sha256"] != manifest_sha
            or row["homology_dimension"] not in {"H0", "H1"}
            or row["query_sample_id"] == row["training_sample_id"]
            or row["pair_scope"] != "held_out_query_to_training_reference"
            or row["landscape_definition_id"] !=
                "all_active_exact_critical_pairs_v1"
            or row["essential_h0_policy"] != "exclude"
            or row["level_policy"] != "all_consecutive_active_levels"
            or row["integration_policy"] !=
                "exact_linear_critical_pair_segments"
            or row["representation"] != "inductive_integrated"
            or row["view_id"] != "cell_topology_v1"
            or row["outcome_label_state"] != "closed"
            or row["biological_outcomes_computed"].lower() == "true"
        ):
            raise ValueError("MV5-I request violates the frozen boundary.")
    if len(implementation_sha) != 64:
        raise ValueError("MV5-I implementation identity is invalid.")


def validate_resume(output_path, status_path, rows, manifest_sha,
                    implementation_sha):
    exists = (os.path.exists(output_path), os.path.exists(status_path))
    if exists == (False, False):
        return False
    if exists[0] != exists[1]:
        raise RuntimeError("MV5-I chunk has a partial output/status pair.")
    output = read_csv(output_path)
    status = read_csv(status_path)
    if (
        len(status) != 1
        or status[0]["status"] != "completed"
        or status[0]["engine_id"] != ENGINE_ID
        or status[0]["pair_manifest_sha256"] != manifest_sha
        or status[0]["implementation_sha256"] != implementation_sha
        or status[0]["request_subset_sha256"] != subset_sha256(rows)
        or status[0]["output_sha256"] != file_sha256(output_path)
        or sorted(row["pair_request_id"] for row in output) !=
            sorted(row["pair_request_id"] for row in rows)
    ):
        raise RuntimeError("Existing MV5-I chunk failed identity validation.")
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--requests", required=True)
    parser.add_argument("--intervals", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--status-dir", required=True)
    parser.add_argument("--pair-manifest-sha256", required=True)
    parser.add_argument("--implementation-sha256", required=True)
    parser.add_argument("--max-pairs", type=int, default=250)
    parser.add_argument("--max-seconds", type=float, default=900.0)
    args = parser.parse_args()

    requests = read_csv(args.requests)
    grouped = defaultdict(list)
    for row in requests:
        grouped[row["chunk_id"]].append(row)
    intervals = defaultdict(list)
    interval_hashes = {}
    for row in read_csv(args.intervals):
        key = (row["record_cache_key"], row["homology_dimension"])
        birth, death = float(row["birth"]), float(row["death"])
        if not math.isfinite(birth) or not math.isfinite(death) or death <= birth:
            raise ValueError("MV5-I interval input is not finite positive-persistence.")
        intervals[key].append((birth, death))
        interval_hashes[key] = row["diagram_sha256"]

    landscapes = {}
    built = reused = 0
    for chunk_id, rows in sorted(
        grouped.items(), key=lambda item: int(item[1][0]["chunk_index"])
    ):
        validate_requests(rows, args.pair_manifest_sha256,
                          args.implementation_sha256, args.max_pairs)
        stem = safe_name(chunk_id)
        output_path = os.path.join(args.output_dir, stem + ".csv")
        status_path = os.path.join(args.status_dir, stem + "__status.csv")
        if validate_resume(output_path, status_path, rows,
                           args.pair_manifest_sha256,
                           args.implementation_sha256):
            reused += 1
            continue
        started = time.perf_counter()
        output_rows = []
        pair_operation_seconds = 0.0
        for row in sorted(rows, key=lambda item: int(item["chunk_offset"])):
            if time.perf_counter() - started > args.max_seconds:
                raise TimeoutError("MV5-I chunk exceeded its elapsed guard.")
            dimension = row["homology_dimension"]
            degree = int(dimension[1:])
            first_key = (row["query_record_cache_key"], dimension)
            second_key = (row["training_record_cache_key"], dimension)
            if first_key not in intervals or second_key not in intervals:
                raise ValueError("MV5-I request lacks finite interval input.")
            if interval_hashes[first_key] != row["query_diagram_sha256"] or \
               interval_hashes[second_key] != row["training_diagram_sha256"]:
                raise ValueError("MV5-I diagram identity differs from intervals.")
            for key in (first_key, second_key):
                if key not in landscapes:
                    landscapes[key] = build_landscape(intervals[key], degree)
            pair_started = time.perf_counter()
            distance, squared, levels, segments, points = landscape_distance(
                landscapes[first_key], landscapes[second_key]
            )
            query_levels = len(landscapes[first_key].critical_pairs)
            training_levels = len(landscapes[second_key].critical_pairs)
            pair_operation_seconds += time.perf_counter() - pair_started
            output_rows.append({
                "contract_id": CONTRACT_ID, "engine_id": ENGINE_ID,
                "pair_request_id": row["pair_request_id"],
                "group_id": row["group_id"], "fold_id": row["fold_id"],
                "seed": row["seed"], "homology_dimension": dimension,
                "query_sample_id": row["query_sample_id"],
                "training_sample_id": row["training_sample_id"],
                "query_record_cache_key": row["query_record_cache_key"],
                "training_record_cache_key": row["training_record_cache_key"],
                "query_diagram_sha256": row["query_diagram_sha256"],
                "training_diagram_sha256": row["training_diagram_sha256"],
                "chunk_id": chunk_id, "difference_active_levels": levels,
                "query_active_levels": query_levels,
                "training_active_levels": training_levels,
                "source_active_levels_processed":
                    query_levels + training_levels,
                "difference_segments": segments,
                "difference_critical_points": points,
                "distance": distance, "squared_distance": squared,
                "exact": True, "all_active_levels": True,
                "level_cap_applied": False,
                "absolute_error_estimate": 0.0,
                "pair_manifest_sha256": args.pair_manifest_sha256,
                "implementation_sha256": args.implementation_sha256,
                "status": "completed", "outcome_label_state": "closed",
                "biological_outcomes_computed": False,
                "retrieval_jobs_executed": 0,
                "clustering_jobs_executed": 0,
                "gene_view_jobs_executed": 0,
                "fusion_jobs_executed": 0,
                "new_data_jobs_executed": 0,
            })
        atomic_csv(output_path, output_rows)
        status_rows = [{
            "contract_id": "mv05i_chunk_status_v1", "engine_id": ENGINE_ID,
            "chunk_id": chunk_id, "group_id": rows[0]["group_id"],
            "homology_dimension": rows[0]["homology_dimension"],
            "request_count": len(rows), "completed_count": len(output_rows),
            "request_subset_sha256": subset_sha256(rows),
            "pair_manifest_sha256": args.pair_manifest_sha256,
            "implementation_sha256": args.implementation_sha256,
            "output_file": os.path.basename(output_path),
            "output_sha256": file_sha256(output_path),
            "elapsed_seconds": time.perf_counter() - started,
            "pair_operation_seconds": pair_operation_seconds,
            "peak_process_rss_bytes":
                resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
            "max_pairs_guard": args.max_pairs,
            "max_seconds_guard": args.max_seconds, "status": "completed",
            "outcome_label_state": "closed",
            "biological_outcomes_computed": False,
        }]
        atomic_csv(status_path, status_rows)
        built += 1
    print(f"MV5-I group complete: built={built} reused={reused} total={len(grouped)}")


if __name__ == "__main__":
    main()
