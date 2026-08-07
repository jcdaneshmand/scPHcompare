#!/usr/bin/env python3
"""Persistent corrected-Persim batch engine for eligible MV-04 diagrams."""

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


ENGINE_ID = "persim_0.3.8_corrected_critical_pairs_batch_v1"
CONTRACT_ID = "full_l2_exact_critical_pairs_v1"


def as_array(values: list[tuple[float, float]]) -> np.ndarray:
    return np.asarray(values, dtype=float).reshape((-1, 2))


def build_landscape(values: list[tuple[float, float]], dimension: int):
    if not values:
        return None
    diagrams = [np.empty((0, 2)) for _ in range(dimension + 1)]
    diagrams[dimension] = as_array(values)
    return PersLandscapeExact(dgms=diagrams, hom_deg=dimension)


def integrate_critical_pairs(critical_pairs) -> tuple[float, int, int]:
    squared = 0.0
    segments = 0
    points = 0
    for level in critical_pairs:
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


def corrected_distance(first, second) -> tuple[float, float, int, int, int]:
    if first is None and second is None:
        return 0.0, 0.0, 0, 0, 0
    if first is None:
        critical_pairs = second.critical_pairs
    elif second is None:
        critical_pairs = first.critical_pairs
    else:
        critical_pairs = (first - second).critical_pairs
    squared, segments, points = integrate_critical_pairs(critical_pairs)
    return math.sqrt(squared), squared, len(critical_pairs), segments, points


def self_test() -> list[dict]:
    cases = {
        "single_tent": ([(0.0, 2.0)], [], 0, math.sqrt(2.0 / 3.0)),
        "translated_tents": (
            [(0.0, 2.0)], [(3.0, 5.0)], 0, math.sqrt(4.0 / 3.0)
        ),
        "sign_crossing": (
            [(0.0, 2.0)], [(0.25, 2.25)], 0, math.sqrt(7.0 / 64.0)
        ),
        "narrow_feature": (
            [(0.499, 0.501)], [], 0, math.sqrt((0.002**3) / 12.0)
        ),
    }
    rows = []
    for case, (first, second, dimension, expected) in cases.items():
        observed, _, _, _, _ = corrected_distance(
            build_landscape(first, dimension),
            build_landscape(second, dimension),
        )
        error = abs(observed - expected)
        passed = error <= max(1e-12, 1e-12 * abs(expected))
        rows.append(
            {
                "case": case,
                "expected": expected,
                "observed": observed,
                "absolute_error": error,
                "passed": passed,
            }
        )
        if not passed:
            raise RuntimeError(f"Corrected critical-pair self-test failed: {case}")
    return rows


def read_manifest(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 56 or len({row["diagram_id"] for row in rows}) != 56:
        raise ValueError("MV-04 manifest must contain 56 unique diagrams.")
    return rows


def read_intervals(path: str):
    values = defaultdict(list)
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            values[(row["diagram_id"], int(row["homology_dimension"]))].append(
                (float(row["birth"]), float(row["death"]))
            )
    return values


def stable_pair_id(stratum_id, first, second, dimension):
    payload = {
        "contract_id": CONTRACT_ID,
        "stratum_id": stratum_id,
        "first_diagram_id": first["diagram_id"],
        "second_diagram_id": second["diagram_id"],
        "first_diagram_sha256": first["diagram_sha256"],
        "second_diagram_sha256": second["diagram_sha256"],
        "homology_dimension": f"H{dimension}",
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return "mv04_landscape_pair_v1:" + hashlib.sha256(encoded).hexdigest()


def write_rows(path: str, rows: list[dict]):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--intervals", required=True)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--build-metrics", required=True)
    parser.add_argument("--self-test", required=True)
    parser.add_argument("--environment", required=True)
    parser.add_argument("--repetition", type=int, required=True)
    parser.add_argument(
        "--stratum", action="append", default=[],
        help="Optional frozen stratum ID to select; repeat for multiple strata.",
    )
    args = parser.parse_args()
    if args.repetition < 1:
        parser.error("--repetition must be positive")

    started = time.perf_counter()
    self_test_rows = self_test()
    write_rows(args.self_test, self_test_rows)
    manifest = read_manifest(args.manifest)
    if args.stratum:
        requested = set(args.stratum)
        known = {row["stratum_id"] for row in manifest}
        if not requested <= known:
            raise ValueError(f"Unknown stratum selection: {sorted(requested - known)}")
        manifest = [row for row in manifest if row["stratum_id"] in requested]
    intervals = read_intervals(args.intervals)
    strata = defaultdict(list)
    for row in manifest:
        strata[row["stratum_id"]].append(row)

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
                build_seconds = time.perf_counter() - build_started
                critical_pairs = [] if landscape is None else landscape.critical_pairs
                build_rows.append(
                    {
                        "repetition": args.repetition,
                        "engine_id": ENGINE_ID,
                        "contract_id": CONTRACT_ID,
                        "stratum_id": stratum_id,
                        "diagram_id": row["diagram_id"],
                        "sample_id": row["sample_id"],
                        "homology_dimension": f"H{dimension}",
                        "finite_intervals": len(values),
                        "active_levels": len(critical_pairs),
                        "critical_points": sum(len(level) for level in critical_pairs),
                        "build_seconds": build_seconds,
                        "peak_process_rss_bytes": resource.getrusage(
                            resource.RUSAGE_SELF
                        ).ru_maxrss
                        * 1024,
                    }
                )
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
                    elapsed = time.perf_counter() - pair_started
                    pair_rows.append(
                        {
                            "repetition": args.repetition,
                            "engine_id": ENGINE_ID,
                            "contract_id": CONTRACT_ID,
                            "pair_id": stable_pair_id(
                                stratum_id, first, second, dimension
                            ),
                            "stratum_id": stratum_id,
                            "cohort": first["cohort"],
                            "representation": first["representation"],
                            "view_id": first["view_id"],
                            "homology_dimension": f"H{dimension}",
                            "first_sample_id": first["sample_id"],
                            "second_sample_id": second["sample_id"],
                            "first_diagram_id": first["diagram_id"],
                            "second_diagram_id": second["diagram_id"],
                            "first_finite_intervals": len(
                                intervals[(first["diagram_id"], dimension)]
                            ),
                            "second_finite_intervals": len(
                                intervals[(second["diagram_id"], dimension)]
                            ),
                            "difference_active_levels": levels,
                            "difference_segments": segments,
                            "difference_critical_points": points,
                            "distance": distance,
                            "squared_distance": squared,
                            "pair_seconds": elapsed,
                            "peak_process_rss_bytes": resource.getrusage(
                                resource.RUSAGE_SELF
                            ).ru_maxrss
                            * 1024,
                            "exact": True,
                            "all_active_levels": True,
                            "absolute_error_estimate": 0.0,
                            "status": "completed",
                        }
                    )
                    write_rows(args.pairs, pair_rows)
            write_rows(args.build_metrics, build_rows)

    environment = {
        "engine_id": ENGINE_ID,
        "contract_id": CONTRACT_ID,
        "repetition": args.repetition,
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "numpy": metadata.version("numpy"),
        "scipy": metadata.version("scipy"),
        "persim": metadata.version("persim"),
        "psutil": metadata.version("psutil"),
        "batch_elapsed_seconds": time.perf_counter() - started,
        "peak_process_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        * 1024,
        "diagram_count": len(manifest),
        "stratum_count": len(strata),
        "dimension_pair_rows": len(pair_rows),
        "process_id": os.getpid(),
    }
    write_rows(args.environment, [environment])
    print(
        f"completed {len(pair_rows)} dimension-pair distances across "
        f"{len(strata)} strata in {environment['batch_elapsed_seconds']:.3f}s"
    )


if __name__ == "__main__":
    main()
