#!/usr/bin/env python3
"""MV6-E grouped exact Persim benchmark over frozen private intervals."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import resource
import sys
import time
from collections import defaultdict


CONTRACT_ID = "mv06e_grouped_persim_result_v1"
ENGINE_ID = "persim_0.3.8_mv05d4_exact_critical_pairs_v1"


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def file_sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_csv(path, rows):
    if not rows:
        raise ValueError("Refusing to publish empty MV6-E output.")
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


def load_engine(path):
    spec = importlib.util.spec_from_file_location("mv05d4_engine", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def result_row(row_type, pair_id, view_id, dimension, first_id, second_id,
               first_count, second_count, result):
    distance, squared, levels, segments, points = result
    return {
        "contract_id": CONTRACT_ID,
        "engine_id": ENGINE_ID,
        "row_type": row_type,
        "pair_id": pair_id,
        "view_id": view_id,
        "homology_dimension": dimension,
        "first_diagram_id": first_id,
        "second_diagram_id": second_id,
        "first_finite_intervals": first_count,
        "second_finite_intervals": second_count,
        "squared_distance": format(squared, ".17g"),
        "distance": format(distance, ".17g"),
        "active_levels": levels,
        "event_segments": segments,
        "critical_points": points,
        "exact": "TRUE",
        "all_active_levels": "TRUE",
        "level_cap_applied": "FALSE",
        "outcome_label_state": "closed",
        "biological_outcomes_computed": "FALSE",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--intervals", required=True)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--references", required=True)
    parser.add_argument("--engine-source", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--metrics", required=True)
    args = parser.parse_args()
    if os.path.exists(args.output) or os.path.exists(args.metrics):
        if not (os.path.exists(args.output) and os.path.exists(args.metrics)):
            raise RuntimeError("MV6-E grouped Persim has an incomplete artifact pair.")
        rows = read_csv(args.output)
        metrics = read_csv(args.metrics)
        if (len(rows) == 240 and
                all(row.get("contract_id") == CONTRACT_ID for row in rows) and
                len(metrics) == 1 and
                metrics[0].get("contract_id") ==
                "mv06e_grouped_persim_metrics_v1"):
            print("Reused validated MV6-E grouped Persim output.")
            return
        raise RuntimeError("MV6-E grouped Persim has partial or stale output.")

    started = time.perf_counter()
    engine = load_engine(args.engine_source)
    intervals_rows = read_csv(args.intervals)
    pairs = read_csv(args.pairs)
    references = read_csv(args.references)
    if len(pairs) != 180 or len(references) != 20:
        raise ValueError("MV6-E pair/reference cardinality differs from prefreeze.")
    intervals = defaultdict(list)
    for row in intervals_rows:
        birth, death = float(row["birth"]), float(row["death"])
        if not death > birth:
            raise ValueError("MV6-E interval is not positive persistence.")
        intervals[(row["diagram_id"], row["homology_dimension"])].append(
            (birth, death)
        )
    for key in intervals:
        intervals[key].sort()
    loaded_seconds = time.perf_counter() - started

    built_started = time.perf_counter()
    landscapes = {}
    for key in sorted(intervals):
        dimension = int(key[1][1])
        landscapes[key] = engine.build_landscape(intervals[key], dimension)
    build_seconds = time.perf_counter() - built_started

    primary_started = time.perf_counter()
    output = []
    primary_by_id = {}
    for row in sorted(pairs, key=lambda value: value["pair_id"]):
        first_key = (row["first_diagram_id"], row["homology_dimension"])
        second_key = (row["second_diagram_id"], row["homology_dimension"])
        result = engine.landscape_distance(
            landscapes.get(first_key), landscapes.get(second_key)
        )
        item = result_row(
            "primary", row["pair_id"], row["view_id"],
            row["homology_dimension"], row["first_diagram_id"],
            row["second_diagram_id"], len(intervals[first_key]),
            len(intervals[second_key]), result
        )
        output.append(item)
        primary_by_id[row["pair_id"]] = item
    primary_seconds = time.perf_counter() - primary_started

    checks_started = time.perf_counter()
    for row in sorted(references, key=lambda value: value["pair_id"]):
        primary = primary_by_id[row["pair_id"]]
        first_key = (primary["second_diagram_id"], primary["homology_dimension"])
        second_key = (primary["first_diagram_id"], primary["homology_dimension"])
        output.append(result_row(
            "reverse", row["pair_id"], primary["view_id"],
            primary["homology_dimension"], primary["second_diagram_id"],
            primary["first_diagram_id"], len(intervals[first_key]),
            len(intervals[second_key]), engine.landscape_distance(
                landscapes.get(first_key), landscapes.get(second_key)
            )
        ))
    for key in sorted(intervals):
        diagram_id, dimension = key
        output.append(result_row(
            "self", "self:" + diagram_id + ":" + dimension,
            next(row["view_id"] for row in pairs
                 if row["first_diagram_id"] == diagram_id or
                 row["second_diagram_id"] == diagram_id),
            dimension, diagram_id, diagram_id, len(intervals[key]),
            len(intervals[key]), engine.landscape_distance(
                landscapes[key], landscapes[key]
            )
        ))
    check_seconds = time.perf_counter() - checks_started
    output.sort(key=lambda row: (
        row["row_type"], row["view_id"], row["homology_dimension"],
        row["pair_id"], row["first_diagram_id"], row["second_diagram_id"]
    ))
    identity = {
        "intervals_sha256": file_sha(args.intervals),
        "pairs_sha256": file_sha(args.pairs),
        "references_sha256": file_sha(args.references),
        "engine_source_sha256": file_sha(args.engine_source),
        "rows": len(output),
    }
    identity_sha = hashlib.sha256(json.dumps(
        identity, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    for row in output:
        row["execution_identity_sha256"] = identity_sha
    atomic_csv(args.output, output)
    total_seconds = time.perf_counter() - started
    metrics = [{
        "contract_id": "mv06e_grouped_persim_metrics_v1",
        "engine_id": ENGINE_ID,
        "diagram_dimensions_built": len(landscapes),
        "primary_pair_rows": len(pairs),
        "reverse_rows": len(references),
        "self_rows": len(intervals),
        "load_seconds": format(loaded_seconds, ".17g"),
        "landscape_build_seconds": format(build_seconds, ".17g"),
        "primary_pair_seconds": format(primary_seconds, ".17g"),
        "reverse_self_seconds": format(check_seconds, ".17g"),
        "total_seconds": format(total_seconds, ".17g"),
        "peak_process_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        "python_version": sys.version.split()[0],
        "persim_version": getattr(sys.modules.get("persim"), "__version__", "0.3.8"),
        "execution_identity_sha256": identity_sha,
        "outcome_label_state": "closed",
        "biological_outcomes_computed": "FALSE",
    }]
    atomic_csv(args.metrics, metrics)
    print(f"Completed MV6-E grouped Persim: {len(output)} scientific rows.")


if __name__ == "__main__":
    main()
