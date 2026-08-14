#!/usr/bin/env python3
"""Independent stratified exact-landscape oracles for MV5-AF."""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict

import numpy as np
from persim.landscapes import PersLandscapeExact


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def safe(value):
    return "".join(char if char.isalnum() or char in "_.-" else "_"
                   for char in value)


def exact_distance(values_a, values_b, degree):
    def build(values):
        if not values:
            return None
        diagrams = [np.empty((0, 2)) for _ in range(degree + 1)]
        diagrams[degree] = np.asarray(values, dtype=float).reshape((-1, 2))
        return PersLandscapeExact(dgms=diagrams, hom_deg=degree)

    first, second = build(values_a), build(values_b)
    if first is None and second is None:
        levels = []
    elif first is None:
        levels = second.critical_pairs
    elif second is None:
        levels = first.critical_pairs
    else:
        levels = (first - second).critical_pairs
    squared = 0.0
    for level in levels:
        for (x0, y0), (x1, y1) in zip(level, level[1:]):
            width = x1 - x0
            if width < 0:
                raise RuntimeError("Independent critical pairs are unsorted")
            squared += width * (y0 * y0 + y0 * y1 + y1 * y1) / 3.0
    return math.sqrt(max(0.0, squared))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", required=True)
    parser.add_argument("--result-root", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--tolerance", required=True, type=float)
    args = parser.parse_args()
    if os.path.exists(args.output):
        raise RuntimeError("Refusing to overwrite independent oracle output")

    queue = read_csv(args.queue)
    fold_sizes = {}
    for row in queue:
        fold_sizes[row["fold_id"]] = int(row["training_samples"])
    folds = sorted(fold_sizes, key=lambda value: (fold_sizes[value], value))
    selected = {folds[0], folds[(len(folds) - 1) // 2], folds[-1]}
    output = []
    for unit in queue:
        if unit["fold_id"] not in selected:
            continue
        root = os.path.join(args.result_root,
                            safe(unit["robustness_group_id"]))
        pairs = sorted(read_csv(os.path.join(root, "pair_scope.csv")),
                       key=lambda row: int(row["pair_ordinal"]))
        pair = pairs[0]
        intervals = read_csv(os.path.join(root, "finite_intervals.csv"))
        observed_rows = read_csv(os.path.join(root, "landscape_pairs.csv"))
        by_key = defaultdict(list)
        for row in intervals:
            by_key[row["sample_id"], row["homology_dimension"]].append(
                (float(row["birth"]), float(row["death"])))
        observed = {
            row["homology_dimension"]: float(row["distance"])
            for row in observed_rows
            if row["pair_request_id"] == pair["pair_request_id"]
        }
        for dimension, degree in (("H0", 0), ("H1", 1)):
            expected = exact_distance(
                sorted(by_key[pair["first_sample_id"], dimension]),
                sorted(by_key[pair["second_sample_id"], dimension]), degree)
            actual = observed.get(dimension)
            error = abs(expected - actual) if actual is not None else math.inf
            if not math.isfinite(error) or error > args.tolerance * max(
                    1.0, abs(expected), abs(actual)):
                raise RuntimeError("Independent landscape oracle mismatch")
            output.append({
                "contract_id": "mv05af_independent_landscape_oracle_v1",
                "robustness_group_id": unit["robustness_group_id"],
                "fold_id": unit["fold_id"],
                "seed": unit["seed"],
                "representation": unit["representation"],
                "pair_request_id": pair["pair_request_id"],
                "homology_dimension": dimension,
                "expected_distance": format(expected, ".17g"),
                "observed_distance": format(actual, ".17g"),
                "absolute_error": format(error, ".17g"),
                "exact_match": "true",
                "production_scientific_helpers_called": "false",
                "labels_opened": "false",
                "rankings_computed": "false",
                "outcomes_computed": "false",
            })
    if len(output) != 60:
        raise RuntimeError("Independent landscape oracle count is not 60")
    with open(args.output, "x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    print("MV5-AF independently reproduced 60 stratified landscape distances.")


if __name__ == "__main__":
    main()
