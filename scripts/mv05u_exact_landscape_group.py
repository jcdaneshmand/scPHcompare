#!/usr/bin/env python3
"""Exact all-active landscape staging for one MV5-U admission unit."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import sys
from collections import defaultdict

import numpy as np
import scipy
import persim

from mv05n_landscape_admission import (
    atomic_csv,
    build_landscape,
    exact_distance,
    file_sha256,
    read_csv,
)


CONTRACT_ID = "mv05u_exact_all_active_landscape_admission_v1"
ENGINE_ID = "persim_exact_critical_pairs_mv05u_v1"


def canonical_sha(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")).hexdigest()


def validate_square_oracle():
    persistence = math.sqrt(2.0) - 1.0
    landscape = build_landscape([(1.0, math.sqrt(2.0))], 1)
    distance, squared, levels, _, _ = exact_distance(landscape, None)
    expected = persistence ** 3 / 12.0
    if (levels != 1 or abs(squared - expected) > 1e-12 or
            abs(distance - math.sqrt(expected)) > 1e-12):
        raise RuntimeError("MV5-U analytic H1 landscape oracle failed.")


def validate_inputs(views, intervals, pairs, admission_unit_id):
    view_required = {
        "admission_unit_id", "sample_id", "diagram_sha256", "h0_intervals",
        "h1_intervals", "outcome_label_state", "outcomes_computed",
    }
    interval_required = {
        "admission_unit_id", "sample_id", "homology_dimension", "birth",
        "death", "diagram_sha256", "outcome_label_state", "outcomes_computed",
    }
    pair_required = {
        "pair_request_id", "pair_ordinal", "pair_scope", "first_sample_id",
        "second_sample_id", "outcome_label_state",
        "biological_outcomes_computed",
    }
    if (len(views) != 90 or not view_required.issubset(views[0])):
        raise ValueError("MV5-U view inventory is incomplete.")
    if not intervals or not interval_required.issubset(intervals[0]):
        raise ValueError("MV5-U finite interval staging is incomplete.")
    if len(pairs) != 32 or not pair_required.issubset(pairs[0]):
        raise ValueError("MV5-U pair coverage is incomplete.")
    sample_ids = [row["sample_id"] for row in views]
    if len(set(sample_ids)) != 90:
        raise ValueError("MV5-U sample identities are duplicated.")
    if len({row["pair_request_id"] for row in pairs}) != 32:
        raise ValueError("MV5-U pair identities are duplicated.")
    forbidden = {"tissue", "approach", "class", "label", "outcome"}
    for rows in (views, intervals, pairs):
        if forbidden.intersection(key.lower() for key in rows[0]):
            raise ValueError("MV5-U private staging contains outcome columns.")
        if any(row.get("outcome_label_state") != "closed" for row in rows):
            raise ValueError("MV5-U label state is not closed.")
        for row in rows:
            value = row.get("outcomes_computed",
                            row.get("biological_outcomes_computed", "false"))
            if value.lower() == "true":
                raise ValueError("MV5-U outcome computation is prohibited.")
    if any(row["admission_unit_id"] != admission_unit_id for row in views):
        raise ValueError("MV5-U view unit identity drifted.")
    if any(row["admission_unit_id"] != admission_unit_id for row in intervals):
        raise ValueError("MV5-U interval unit identity drifted.")
    for row in pairs:
        if row["first_sample_id"] not in sample_ids or row["second_sample_id"] not in sample_ids:
            raise ValueError("MV5-U pair references an absent sample.")
    return sample_ids


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--views")
    parser.add_argument("--intervals")
    parser.add_argument("--pairs")
    parser.add_argument("--summary-output")
    parser.add_argument("--pair-output")
    parser.add_argument("--admission-unit-id")
    parser.add_argument("--script-sha256")
    parser.add_argument("--oracle-only", action="store_true")
    parser.add_argument("--environment-only", action="store_true")
    args = parser.parse_args()

    if args.script_sha256 is None or file_sha256(__file__) != args.script_sha256:
        raise ValueError("MV5-U landscape implementation hash is stale.")
    validate_square_oracle()
    if args.environment_only:
        print("python_version\t" + platform.python_version())
        print("persim_version\t" + persim.__version__)
        print("numpy_version\t" + np.__version__)
        print("scipy_version\t" + scipy.__version__)
        print("python_implementation\t" + platform.python_implementation())
        print("python_major_minor\t" + f"{sys.version_info.major}.{sys.version_info.minor}")
        return
    if args.oracle_only:
        print({"analytic_square_h1_landscape_oracle": "passed"})
        return
    required = (args.views, args.intervals, args.pairs, args.summary_output,
                args.pair_output, args.admission_unit_id)
    if any(value is None for value in required):
        raise ValueError("MV5-U landscape execution arguments are incomplete.")
    views = read_csv(args.views)
    intervals_rows = read_csv(args.intervals)
    pairs = read_csv(args.pairs)
    sample_ids = validate_inputs(
        views, intervals_rows, pairs, args.admission_unit_id
    )

    expected_hashes = {row["sample_id"]: row["diagram_sha256"] for row in views}
    intervals = defaultdict(list)
    observed_hashes = {}
    for row in intervals_rows:
        dimension = row["homology_dimension"]
        if dimension not in {"H0", "H1"}:
            raise ValueError("MV5-U interval dimension is invalid.")
        birth, death = float(row["birth"]), float(row["death"])
        if not math.isfinite(birth) or not math.isfinite(death) or death <= birth:
            raise ValueError("MV5-U intervals must have finite positive persistence.")
        key = row["sample_id"], dimension
        intervals[key].append((birth, death))
        observed_hashes[row["sample_id"]] = row["diagram_sha256"]
    if any(observed_hashes.get(sample_id) != expected_hashes[sample_id]
           for sample_id in sample_ids):
        raise ValueError("MV5-U staged diagram identity is incomplete or stale.")

    landscapes = {}
    summary = []
    for sample_id in sorted(sample_ids):
        for dimension in ("H0", "H1"):
            degree = int(dimension[1:])
            values = sorted(intervals.get((sample_id, dimension), []))
            landscape = build_landscape(values, degree)
            landscapes[sample_id, dimension] = landscape
            critical_pairs = [] if landscape is None else landscape.critical_pairs
            canonical_pairs = [
                [[float(x), float(y)] for x, y in level]
                for level in critical_pairs
            ]
            summary.append({
                "contract_id": CONTRACT_ID,
                "engine_id": ENGINE_ID,
                "admission_unit_id": args.admission_unit_id,
                "sample_id": sample_id,
                "homology_dimension": dimension,
                "finite_intervals": len(values),
                "active_levels": len(critical_pairs),
                "critical_points": sum(len(level) for level in critical_pairs),
                "landscape_sha256": canonical_sha(canonical_pairs),
                "essential_h0_policy": "exclude",
                "level_policy": "all_consecutive_active_levels",
                "integration_policy": "exact_linear_critical_pair_segments",
                "outcome_label_state": "closed",
                "outcomes_computed": False,
            })

    output = []
    for pair in sorted(pairs, key=lambda row: int(row["pair_ordinal"])):
        for dimension in ("H0", "H1"):
            first = landscapes[pair["first_sample_id"], dimension]
            second = landscapes[pair["second_sample_id"], dimension]
            distance, squared, levels, segments, points = exact_distance(first, second)
            output.append({
                "contract_id": CONTRACT_ID,
                "engine_id": ENGINE_ID,
                "admission_unit_id": args.admission_unit_id,
                "pair_request_id": pair["pair_request_id"],
                "pair_ordinal": pair["pair_ordinal"],
                "pair_scope": pair["pair_scope"],
                "first_sample_id": pair["first_sample_id"],
                "second_sample_id": pair["second_sample_id"],
                "homology_dimension": dimension,
                "distance": distance,
                "squared_distance": squared,
                "difference_active_levels": levels,
                "difference_segments": segments,
                "difference_critical_points": points,
                "exact": True,
                "all_active_levels": True,
                "level_cap_applied": False,
                "absolute_error_estimate": 0.0,
                "outcome_label_state": "closed",
                "outcomes_computed": False,
            })
    if len(summary) != 180 or len(output) != 64:
        raise RuntimeError("MV5-U landscape output cardinality is invalid.")
    atomic_csv(args.summary_output, summary)
    atomic_csv(args.pair_output, output)
    print({"unit": args.admission_unit_id, "landscapes": len(summary),
           "pair_distances": len(output), "square_oracle": "passed"})


if __name__ == "__main__":
    main()
