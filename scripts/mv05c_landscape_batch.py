#!/usr/bin/env python3
"""Exact all-level persistence-landscape distances for MV5-C strata."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import resource
import sys
import time
from collections import defaultdict
from importlib import metadata

import numpy as np
from persim import PersLandscapeExact


ENGINE_ID = "persim_0.3.8_mv05c_exact_critical_pairs_v1"
CONTRACT_ID = "mv05c_full_l2_exact_all_active_levels_v1"


def write_rows(path: str, rows: list[dict]):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


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


def stable_pair_id(stratum_id, first, second, dimension):
    encoded = json.dumps(
        {
            "contract_id": CONTRACT_ID,
            "stratum_id": stratum_id,
            "first_diagram_id": first["diagram_id"],
            "second_diagram_id": second["diagram_id"],
            "first_diagram_sha256": first["diagram_sha256"],
            "second_diagram_sha256": second["diagram_sha256"],
            "homology_dimension": f"H{dimension}",
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return "mv05c_landscape_pair_v1:" + hashlib.sha256(encoded).hexdigest()


def self_test():
    first = build_landscape([(0.0, 2.0)], 0)
    observed, _, _, _, _ = corrected_distance(first, None)
    expected = math.sqrt(2.0 / 3.0)
    error = abs(observed - expected)
    if error > 1e-12:
        raise RuntimeError("Exact landscape self-test failed.")
    return [{"case": "single_tent", "expected": expected,
             "observed": observed, "absolute_error": error, "passed": True}]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--intervals", required=True)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--build-metrics", required=True)
    parser.add_argument("--self-test", required=True)
    parser.add_argument("--environment", required=True)
    args = parser.parse_args()

    started = time.perf_counter()
    write_rows(args.self_test, self_test())
    with open(args.manifest, newline="", encoding="utf-8") as handle:
        manifest = list(csv.DictReader(handle))
    if not manifest or len({row["diagram_id"] for row in manifest}) != len(manifest):
        raise ValueError("MV5-C diagram manifest is empty or has duplicate IDs.")
    with open(args.intervals, newline="", encoding="utf-8") as handle:
        interval_rows = list(csv.DictReader(handle))
    intervals = defaultdict(list)
    for row in interval_rows:
        intervals[(row["diagram_id"], int(row["homology_dimension"]))].append(
            (float(row["birth"]), float(row["death"]))
        )
    strata = defaultdict(list)
    for row in manifest:
        strata[row["stratum_id"]].append(row)
    if any(len(rows) != 6 for rows in strata.values()):
        raise ValueError("Every completed MV5-C stratum must contain six diagrams.")

    pair_rows = []
    build_rows = []
    for stratum_id in sorted(strata):
        diagrams = sorted(strata[stratum_id], key=lambda row: row["sample_id"])
        for dimension in (0, 1):
            landscapes = {}
            for row in diagrams:
                values = intervals[(row["diagram_id"], dimension)]
                build_started = time.perf_counter()
                landscape = build_landscape(values, dimension)
                elapsed = time.perf_counter() - build_started
                critical = [] if landscape is None else landscape.critical_pairs
                build_rows.append({
                    "engine_id": ENGINE_ID, "contract_id": CONTRACT_ID,
                    "stratum_id": stratum_id, "diagram_id": row["diagram_id"],
                    "sample_id": row["sample_id"],
                    "homology_dimension": f"H{dimension}",
                    "finite_intervals": len(values),
                    "active_levels": len(critical),
                    "critical_points": sum(len(level) for level in critical),
                    "build_seconds": elapsed,
                    "peak_process_rss_bytes": resource.getrusage(
                        resource.RUSAGE_SELF
                    ).ru_maxrss * 1024,
                    "outcome_label_state": "closed",
                    "biological_outcomes_computed": False,
                })
                landscapes[row["diagram_id"]] = landscape
            for left in range(len(diagrams)):
                for right in range(left + 1, len(diagrams)):
                    first = diagrams[left]
                    second = diagrams[right]
                    pair_started = time.perf_counter()
                    distance, squared, levels, segments, points = corrected_distance(
                        landscapes[first["diagram_id"]],
                        landscapes[second["diagram_id"]],
                    )
                    pair_rows.append({
                        "engine_id": ENGINE_ID, "contract_id": CONTRACT_ID,
                        "pair_id": stable_pair_id(
                            stratum_id, first, second, dimension
                        ),
                        "stratum_id": stratum_id, "fold_id": first["fold_id"],
                        "seed": first["seed"],
                        "representation": first["representation"],
                        "view_id": first["view_id"],
                        "homology_dimension": f"H{dimension}",
                        "first_sample_id": first["sample_id"],
                        "second_sample_id": second["sample_id"],
                        "first_diagram_id": first["diagram_id"],
                        "second_diagram_id": second["diagram_id"],
                        "difference_active_levels": levels,
                        "difference_segments": segments,
                        "difference_critical_points": points,
                        "distance": distance, "squared_distance": squared,
                        "pair_seconds": time.perf_counter() - pair_started,
                        "peak_process_rss_bytes": resource.getrusage(
                            resource.RUSAGE_SELF
                        ).ru_maxrss * 1024,
                        "exact": True, "all_active_levels": True,
                        "absolute_error_estimate": 0.0, "status": "completed",
                        "outcome_label_state": "closed",
                        "biological_outcomes_computed": False,
                    })
                    write_rows(args.pairs, pair_rows)
            write_rows(args.build_metrics, build_rows)

    environment = [{
        "engine_id": ENGINE_ID, "contract_id": CONTRACT_ID,
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(), "numpy": metadata.version("numpy"),
        "scipy": metadata.version("scipy"), "persim": metadata.version("persim"),
        "elapsed_seconds": time.perf_counter() - started,
        "peak_process_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        "diagram_count": len(manifest), "stratum_count": len(strata),
        "dimension_pair_rows": len(pair_rows),
        "outcome_label_state": "closed",
        "biological_outcomes_computed": False,
    }]
    write_rows(args.environment, environment)
    if len(pair_rows) != len(strata) * 30:
        raise RuntimeError("MV5-C pair-row count is incomplete.")
    print(f"completed {len(pair_rows)} exact distances across {len(strata)} strata")


if __name__ == "__main__":
    main()
