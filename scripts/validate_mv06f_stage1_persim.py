#!/usr/bin/env python3
"""Validate frozen MV6-F stage-one oracle rows with grouped exact Persim."""

import argparse
import csv
import importlib.util
import math
import os


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_engine(path):
    spec = importlib.util.spec_from_file_location("mv05d4_engine", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def atomic_csv(path, rows):
    if not rows:
        raise ValueError("Refusing to publish empty MV6-F Persim evidence.")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    partial = f"{path}.partial.{os.getpid()}"
    try:
        with open(partial, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        if os.path.exists(path):
            raise FileExistsError(f"Refusing to overwrite {path}")
        os.replace(partial, path)
    finally:
        if os.path.exists(partial):
            os.unlink(partial)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--intervals", required=True)
    parser.add_argument("--r-oracles", required=True)
    parser.add_argument("--engine-source", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    engine = load_engine(args.engine_source)
    grouped = {}
    for row in read_csv(args.intervals):
        grouped.setdefault(row["diagram_id"], []).append(
            (float(row["birth"]), float(row["death"]))
        )
    landscapes = {}
    output = []
    for row in read_csv(args.r_oracles):
        degree = int(row["homology_dimension"][1:])
        first_id, second_id = row["first_diagram_id"], row["second_diagram_id"]
        keys = ((first_id, degree), (second_id, degree))
        for key in keys:
            if key not in landscapes:
                landscapes[key] = engine.build_landscape(grouped[key[0]], key[1])
        distance, squared, levels, segments, points = engine.landscape_distance(
            landscapes[keys[0]], landscapes[keys[1]]
        )
        rust_squared = float(row["rust_squared_distance"])
        r_squared = float(row["r_squared_distance"])
        tolerance = 1e-10 + 1e-10 * abs(squared)
        rust_error = abs(squared - rust_squared)
        r_tolerance = float(row["acceptance_tolerance"])
        r_error = abs(squared - r_squared)
        passed = (math.isfinite(squared) and squared >= 0 and
                  rust_error <= tolerance and
                  r_error <= max(tolerance, r_tolerance))
        output.append({
            "contract_id": "mv06f_stage1_persim_oracle_v1",
            "pair_id": row["pair_id"], "view_id": row["view_id"],
            "homology_dimension": row["homology_dimension"],
            "selection_stratum": row["selection_stratum"],
            "persim_squared_distance": format(squared, ".17g"),
            "rust_squared_distance": row["rust_squared_distance"],
            "r_squared_distance": row["r_squared_distance"],
            "absolute_error_vs_rust": format(rust_error, ".17g"),
            "absolute_error_vs_r": format(r_error, ".17g"),
            "acceptance_tolerance": format(max(tolerance, r_tolerance), ".17g"),
            "active_levels": levels, "event_segments": segments,
            "critical_points": points, "exact": "TRUE",
            "all_active_levels": "TRUE", "level_cap_applied": "FALSE",
            "passed": "TRUE" if passed else "FALSE",
            "outcome_label_state": "closed",
            "biological_outcomes_computed": "FALSE",
        })
    failures = [row for row in output if row["passed"] != "TRUE"]
    if len(output) != 12 or failures:
        for row in failures:
            print("FAILED", row)
        raise RuntimeError("MV6-F stage-one Persim oracle validation failed.")
    atomic_csv(args.output, output)
    print("Validated MV6-F stage 1: 12/12 grouped Persim oracle rows.")


if __name__ == "__main__":
    main()
