#!/usr/bin/env python3
"""Bounded label-closed exact all-level landscape admission for MV5-N."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import time
from collections import defaultdict

import numpy as np
from persim import PersLandscapeExact

try:
    import resource

    def peak_process_rss_bytes():
        value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # Linux reports KiB; macOS reports bytes.
        return value if value > 10**9 else value * 1024
except ImportError:  # pragma: no cover - exercised by the Windows admission.
    import ctypes
    from ctypes import wintypes

    class PROCESS_MEMORY_COUNTERS(ctypes.Structure):
        _fields_ = [
            ("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
            ("PeakWorkingSetSize", ctypes.c_size_t),
            ("WorkingSetSize", ctypes.c_size_t),
            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
            ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
            ("PagefileUsage", ctypes.c_size_t),
            ("PeakPagefileUsage", ctypes.c_size_t),
        ]

    def peak_process_rss_bytes():
        counters = PROCESS_MEMORY_COUNTERS()
        counters.cb = ctypes.sizeof(counters)
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        psapi = ctypes.WinDLL("psapi", use_last_error=True)
        kernel32.GetCurrentProcess.restype = wintypes.HANDLE
        query = psapi.GetProcessMemoryInfo
        query.argtypes = [wintypes.HANDLE,
                          ctypes.POINTER(PROCESS_MEMORY_COUNTERS), wintypes.DWORD]
        query.restype = wintypes.BOOL
        ok = query(kernel32.GetCurrentProcess(), ctypes.byref(counters), counters.cb)
        if not ok:
            raise OSError("Could not read Windows process memory counters.")
        return int(counters.PeakWorkingSetSize)


ENGINE_ID = "persim_0.3.8_mv05n_exact_critical_pairs_admission_v1"
CONTRACT_ID = "mv05n_all_active_exact_landscape_admission_v1"


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rowset_sha256(rows):
    payload = [{key: row[key] for key in sorted(row)} for row in sorted(
        rows, key=lambda item: item["pair_request_id"]
    )]
    return hashlib.sha256(json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()


def atomic_csv(path, rows):
    if not rows:
        raise ValueError("Refusing to write an empty MV5-N artifact.")
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


def build_landscape(values, degree):
    if not values:
        return None
    diagrams = [np.empty((0, 2)) for _ in range(degree + 1)]
    diagrams[degree] = np.asarray(values, dtype=float).reshape((-1, 2))
    return PersLandscapeExact(dgms=diagrams, hom_deg=degree)


def exact_distance(first, second):
    if first is None and second is None:
        levels = []
    elif first is None:
        levels = second.critical_pairs
    elif second is None:
        levels = first.critical_pairs
    else:
        levels = (first - second).critical_pairs
    squared = 0.0
    segments = points = 0
    for level in levels:
        points += len(level)
        for (x0, y0), (x1, y1) in zip(level, level[1:]):
            width = x1 - x0
            if width < 0:
                raise ValueError("MV5-N critical-pair coordinates are unsorted.")
            squared += width * (y0*y0 + y0*y1 + y1*y1) / 3.0
            segments += 1
    if squared < -1e-12:
        raise ValueError("MV5-N exact squared landscape distance is negative.")
    squared = max(0.0, squared)
    return math.sqrt(squared), squared, len(levels), segments, points


def validate_requests(rows, max_rows):
    required = {
        "admission_contract_id", "admission_group_id", "pair_request_id",
        "profile", "fold_id", "seed", "representation",
        "homology_dimension", "first_sample_id", "second_sample_id",
        "first_record_cache_key", "second_record_cache_key",
        "first_diagram_sha256", "second_diagram_sha256", "pair_scope",
        "landscape_definition_id", "essential_h0_policy", "level_policy",
        "integration_policy", "outcome_label_state",
        "biological_outcomes_computed", "full_production_authorized",
    }
    if not rows or not required.issubset(rows[0]) or len(rows) > max_rows:
        raise ValueError("MV5-N admission request schema or size is invalid.")
    if len({row["pair_request_id"] for row in rows}) != len(rows):
        raise ValueError("MV5-N admission request identifiers are duplicated.")
    forbidden = {"tissue", "approach", "class", "label", "outcome"}
    if forbidden.intersection(key.lower() for key in rows[0]):
        raise ValueError("MV5-N admission requests contain outcome columns.")
    for row in rows:
        if (
            row["admission_contract_id"] != "mv05n_label_closed_admission_request_v1"
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
            or row["full_production_authorized"].lower() == "true"
        ):
            raise ValueError("MV5-N admission request violates its frozen boundary.")


def validate_resume(output_path, status_path, rows, request_sha, implementation_sha):
    exists = os.path.exists(output_path), os.path.exists(status_path)
    if exists == (False, False):
        return False
    if exists[0] != exists[1]:
        raise RuntimeError("MV5-N admission has a partial output/status pair.")
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
        or sorted(row["pair_request_id"] for row in output) !=
           sorted(row["pair_request_id"] for row in rows)
    ):
        raise RuntimeError("Existing MV5-N admission output failed resume validation.")
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--requests", required=True)
    parser.add_argument("--intervals", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--status-dir", required=True)
    parser.add_argument("--implementation-sha256", required=True)
    parser.add_argument("--max-rows-per-group", type=int, default=32)
    parser.add_argument("--max-seconds-per-group", type=float, default=900.0)
    args = parser.parse_args()
    if len(args.implementation_sha256) != 64:
        raise ValueError("MV5-N implementation SHA-256 is invalid.")

    requests = read_csv(args.requests)
    groups = defaultdict(list)
    for row in requests:
        groups[row["admission_group_id"]].append(row)
    if len(groups) != 12:
        raise ValueError("MV5-N admission must contain exactly 12 profile groups.")
    intervals = defaultdict(list)
    interval_hashes = {}
    for row in read_csv(args.intervals):
        key = row["record_cache_key"], row["homology_dimension"]
        birth, death = float(row["birth"]), float(row["death"])
        if not math.isfinite(birth) or not math.isfinite(death) or death <= birth:
            raise ValueError("MV5-N interval input is not finite positive-persistence.")
        intervals[key].append((birth, death))
        interval_hashes[key] = row["diagram_sha256"]

    request_sha = file_sha256(args.requests)
    landscapes = {}
    built = reused = 0
    for group_id, rows in sorted(groups.items()):
        validate_requests(rows, args.max_rows_per_group)
        stem = safe_name(group_id)
        output_path = os.path.join(args.output_dir, stem + ".csv")
        status_path = os.path.join(args.status_dir, stem + "__status.csv")
        if validate_resume(output_path, status_path, rows, request_sha,
                           args.implementation_sha256):
            reused += 1
            continue
        started = time.perf_counter()
        operation_seconds = 0.0
        output_rows = []
        for row in sorted(rows, key=lambda item: item["pair_request_id"]):
            if time.perf_counter() - started > args.max_seconds_per_group:
                raise TimeoutError("MV5-N admission group exceeded elapsed guard.")
            dimension = row["homology_dimension"]
            degree = int(dimension[1:])
            first_key = row["first_record_cache_key"], dimension
            second_key = row["second_record_cache_key"], dimension
            if first_key not in intervals or second_key not in intervals:
                raise ValueError("MV5-N request lacks staged finite intervals.")
            if (interval_hashes[first_key] != row["first_diagram_sha256"] or
                    interval_hashes[second_key] != row["second_diagram_sha256"]):
                raise ValueError("MV5-N staged diagram identity is stale.")
            for key in (first_key, second_key):
                if key not in landscapes:
                    landscapes[key] = build_landscape(intervals[key], degree)
            operation_started = time.perf_counter()
            distance, squared, levels, segments, points = exact_distance(
                landscapes[first_key], landscapes[second_key]
            )
            operation_seconds += time.perf_counter() - operation_started
            output_rows.append({
                "contract_id": CONTRACT_ID, "engine_id": ENGINE_ID,
                "pair_request_id": row["pair_request_id"],
                "admission_group_id": group_id, "profile": row["profile"],
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
                "implementation_sha256": args.implementation_sha256,
                "request_file_sha256": request_sha, "status": "completed",
                "outcome_label_state": "closed",
                "biological_outcomes_computed": False,
                "clustering_jobs_executed": 0,
            })
        atomic_csv(output_path, output_rows)
        status_rows = [{
            "contract_id": "mv05n_admission_group_status_v1",
            "engine_id": ENGINE_ID, "admission_group_id": group_id,
            "profile": rows[0]["profile"],
            "representation": rows[0]["representation"],
            "homology_dimension": rows[0]["homology_dimension"],
            "request_count": len(rows), "completed_count": len(output_rows),
            "request_subset_sha256": rowset_sha256(rows),
            "request_file_sha256": request_sha,
            "implementation_sha256": args.implementation_sha256,
            "output_file": os.path.basename(output_path),
            "output_sha256": file_sha256(output_path),
            "elapsed_seconds": time.perf_counter() - started,
            "pair_operation_seconds": operation_seconds,
            "peak_process_rss_bytes": peak_process_rss_bytes(),
            "max_rows_guard": args.max_rows_per_group,
            "max_seconds_guard": args.max_seconds_per_group,
            "status": "completed", "outcome_label_state": "closed",
            "biological_outcomes_computed": False,
            "clustering_jobs_executed": 0,
        }]
        atomic_csv(status_path, status_rows)
        built += 1
    print(f"MV5-N admission complete: built={built} reused={reused} total={len(groups)}")


if __name__ == "__main__":
    main()
