#!/usr/bin/env python3
"""Exact all-active-level MV5-O landscape production for one frozen group."""

from __future__ import annotations

import argparse
import math
import os
import time
from collections import defaultdict

from mv05n_landscape_admission import (
    CONTRACT_ID as ADMISSION_CONTRACT_ID,
    atomic_csv,
    build_landscape,
    exact_distance,
    file_sha256,
    peak_process_rss_bytes,
    read_csv,
    rowset_sha256,
    safe_name,
)

CONTRACT_ID = "mv05o_all_active_exact_landscape_group_v1"
ENGINE_ID = "persim_exact_critical_pairs_all_active_v1"


def validate_requests(rows, implementation_sha256, max_rows_per_chunk):
    required = {
        "production_group_id", "production_chunk_id", "pair_request_id",
        "fold_id", "seed", "representation", "homology_dimension",
        "first_sample_id", "second_sample_id", "first_record_cache_key",
        "second_record_cache_key", "first_diagram_sha256",
        "second_diagram_sha256", "pair_scope", "landscape_definition_id",
        "essential_h0_policy", "level_policy", "integration_policy",
        "source_freeze_sha256", "implementation_sha256",
        "outcome_label_state", "biological_outcomes_computed",
        "production_executed", "clustering_jobs_executed",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError("MV5-O group request schema is incomplete.")
    forbidden = {"tissue", "approach", "class", "label", "outcome"}
    if forbidden.intersection(key.lower() for key in rows[0]):
        raise ValueError("MV5-O group requests contain prohibited labels.")
    if len({row["pair_request_id"] for row in rows}) != len(rows):
        raise ValueError("MV5-O group request identifiers are duplicated.")
    if len({row["production_group_id"] for row in rows}) != 1:
        raise ValueError("MV5-O runner accepts exactly one production group.")
    chunks = defaultdict(list)
    for row in rows:
        chunks[row["production_chunk_id"]].append(row)
        if (
            row["implementation_sha256"] != implementation_sha256
            or row["representation"] not in {"sct_whole", "inductive_integrated"}
            or row["homology_dimension"] not in {"H0", "H1"}
            or row["first_sample_id"] >= row["second_sample_id"]
            or row["pair_scope"] != "training_training_unordered"
            or row["landscape_definition_id"] != "all_active_exact_critical_pairs_v1"
            or row["essential_h0_policy"] != "exclude"
            or row["level_policy"] != "all_consecutive_active_levels"
            or row["integration_policy"] != "exact_linear_critical_pair_segments"
            or row["outcome_label_state"] != "closed"
            or row["biological_outcomes_computed"].lower() == "true"
            or row["production_executed"].lower() == "true"
            or int(row["clustering_jobs_executed"]) != 0
        ):
            raise ValueError("MV5-O request violates the prospective freeze.")
    if any(len(chunk) > max_rows_per_chunk for chunk in chunks.values()):
        raise ValueError("MV5-O landscape chunk exceeds its 250-row cap.")
    return chunks


def validate_resume(output_path, status_path, rows, request_sha, implementation_sha):
    exists = os.path.exists(output_path), os.path.exists(status_path)
    if exists == (False, False):
        return False
    if exists[0] != exists[1]:
        raise RuntimeError("MV5-O found a partial output/status pair.")
    output = read_csv(output_path)
    status = read_csv(status_path)
    if (
        len(status) != 1
        or status[0]["status"] != "completed"
        or status[0]["engine_id"] != ENGINE_ID
        or status[0]["request_file_sha256"] != request_sha
        or status[0]["implementation_sha256"] != implementation_sha
        or status[0]["request_subset_sha256"] != rowset_sha256(rows)
        or status[0]["output_sha256"] != file_sha256(output_path)
        or sorted(row["pair_request_id"] for row in output)
        != sorted(row["pair_request_id"] for row in rows)
    ):
        raise RuntimeError("MV5-O existing chunk failed immutable resume validation.")
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--requests", required=True)
    parser.add_argument("--intervals", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--status-dir", required=True)
    parser.add_argument("--implementation-sha256", required=True)
    parser.add_argument("--max-rows-per-chunk", type=int, default=250)
    parser.add_argument("--max-seconds-per-group", type=float, default=900.0)
    args = parser.parse_args()
    if file_sha256(__file__) != args.implementation_sha256:
        raise ValueError("MV5-O landscape implementation hash is stale.")
    requests = read_csv(args.requests)
    chunks = validate_requests(requests, args.implementation_sha256,
                               args.max_rows_per_chunk)
    intervals = defaultdict(list)
    interval_hashes = {}
    for row in read_csv(args.intervals):
        key = row["record_cache_key"], row["homology_dimension"]
        birth, death = float(row["birth"]), float(row["death"])
        if not math.isfinite(birth) or not math.isfinite(death) or not death > birth:
            raise ValueError("MV5-O intervals must have positive finite persistence.")
        intervals[key].append((birth, death))
        interval_hashes[key] = row["diagram_sha256"]
    request_sha = file_sha256(args.requests)
    landscapes = {}
    group_started = time.perf_counter()
    built = reused = 0
    for chunk_id, rows in sorted(chunks.items()):
        stem = safe_name(chunk_id)
        output_path = os.path.join(args.output_dir, stem + ".csv")
        status_path = os.path.join(args.status_dir, stem + "__status.csv")
        if validate_resume(output_path, status_path, rows, request_sha,
                           args.implementation_sha256):
            reused += 1
            continue
        if time.perf_counter() - group_started > args.max_seconds_per_group:
            raise TimeoutError("MV5-O landscape group exceeded 900 seconds.")
        chunk_started = time.perf_counter()
        operation_seconds = 0.0
        output = []
        for row in sorted(rows, key=lambda item: item["pair_request_id"]):
            dimension = row["homology_dimension"]
            degree = int(dimension[1:])
            first = row["first_record_cache_key"], dimension
            second = row["second_record_cache_key"], dimension
            if first not in intervals or second not in intervals:
                raise ValueError("MV5-O request lacks finite staged intervals.")
            if (interval_hashes[first] != row["first_diagram_sha256"] or
                    interval_hashes[second] != row["second_diagram_sha256"]):
                raise ValueError("MV5-O staged diagram identity is stale.")
            for key in (first, second):
                if key not in landscapes:
                    landscapes[key] = build_landscape(intervals[key], degree)
            started = time.perf_counter()
            distance, squared, levels, segments, points = exact_distance(
                landscapes[first], landscapes[second]
            )
            operation_seconds += time.perf_counter() - started
            output.append({
                "contract_id": CONTRACT_ID, "engine_id": ENGINE_ID,
                "production_group_id": row["production_group_id"],
                "production_chunk_id": chunk_id,
                "pair_request_id": row["pair_request_id"],
                "fold_id": row["fold_id"], "seed": row["seed"],
                "representation": row["representation"],
                "homology_dimension": dimension,
                "first_sample_id": row["first_sample_id"],
                "second_sample_id": row["second_sample_id"],
                "distance": distance, "squared_distance": squared,
                "difference_active_levels": levels,
                "difference_segments": segments,
                "difference_critical_points": points,
                "exact": True, "all_active_levels": True,
                "level_cap_applied": False, "absolute_error_estimate": 0.0,
                "source_freeze_sha256": row["source_freeze_sha256"],
                "implementation_sha256": args.implementation_sha256,
                "request_file_sha256": request_sha, "status": "completed",
                "outcome_label_state": "closed",
                "biological_outcomes_computed": False,
                "clustering_jobs_executed": 0,
            })
        atomic_csv(output_path, output)
        elapsed = time.perf_counter() - chunk_started
        atomic_csv(status_path, [{
            "contract_id": "mv05o_landscape_chunk_status_v1",
            "engine_id": ENGINE_ID, "production_chunk_id": chunk_id,
            "status": "completed", "request_count": len(rows),
            "elapsed_seconds": elapsed,
            "pair_operation_seconds": operation_seconds,
            "peak_process_rss_bytes": peak_process_rss_bytes(),
            "request_file_sha256": request_sha,
            "request_subset_sha256": rowset_sha256(rows),
            "implementation_sha256": args.implementation_sha256,
            "output_sha256": file_sha256(output_path),
            "output_size_bytes": os.path.getsize(output_path),
            "outcome_label_state": "closed",
            "biological_outcomes_computed": False,
            "clustering_jobs_executed": 0,
        }])
        built += 1
    print({"group": requests[0]["production_group_id"], "built": built,
           "reused": reused, "chunks": len(chunks)})


if __name__ == "__main__":
    main()
