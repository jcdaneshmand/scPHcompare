#!/usr/bin/env python3
"""Corrected-Persim equivalence and timing against accepted MV5-AY values."""

from __future__ import annotations

import argparse
import csv
import ctypes
import gc
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


ENGINE_ID = "persim_0.3.8_corrected_segment_integral_mv05ba_v1"


def release_memory():
    gc.collect()
    try:
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except (AttributeError, OSError):
        pass


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path, rows):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    temporary = path + ".tmp"
    with open(temporary, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def build_landscape(values, degree):
    values = np.asarray(values, dtype=float).reshape((-1, 2))
    if len(values) == 0:
        return None
    diagrams = [np.empty((0, 2)) for _ in range(degree + 1)]
    diagrams[degree] = values
    return PersLandscapeExact(dgms=diagrams, hom_deg=degree)


def corrected_squared(first, second):
    if first is None and second is None:
        levels = []
    elif first is None:
        levels = second.critical_pairs
    elif second is None:
        levels = first.critical_pairs
    else:
        levels = (first - second).critical_pairs
    squared = 0.0
    segments = 0
    points = 0
    for level in levels:
        points += len(level)
        for (x0, y0), (x1, y1) in zip(level, level[1:]):
            width = x1 - x0
            if width < 0:
                raise ValueError("Corrected critical pairs are unsorted.")
            squared += width * (y0 * y0 + y0 * y1 + y1 * y1) / 3.0
            segments += 1
    if squared < -1e-10:
        raise ValueError("Negative corrected squared distance.")
    return max(0.0, squared), len(levels), segments, points


def fixture_rows():
    cases = {
        "single_tent": ([(0.0, 2.0)], [], 2.0 / 3.0),
        "sign_changing_tents": ([(0.0, 2.0)], [(0.25, 2.25)], 7.0 / 64.0),
        "narrow_feature": ([(0.499, 0.501)], [], 0.002**3 / 12.0),
    }
    rows = []
    for case, (first, second, expected) in cases.items():
        observed, levels, segments, points = corrected_squared(
            build_landscape(first, 0), build_landscape(second, 0)
        )
        rows.append({
            "contract_id": "mv05ba_fixture_v1", "case": case,
            "expected_squared_distance": format(expected, ".17g"),
            "observed_squared_distance": format(observed, ".17g"),
            "absolute_error": format(abs(observed - expected), ".17g"),
            "passed": abs(observed - expected) <= 1e-12,
            "active_levels": levels, "segments": segments,
            "critical_points": points,
        })
    return rows


def stable_identity(row, squared):
    value = {
        "engine_id": ENGINE_ID, "stratum_id": row["stratum_id"],
        "pair_order": int(row["pair_order"]),
        "pair_cache_key": row["pair_cache_key"],
        "dimension": row["dimension"],
        "squared_distance": format(squared, ".17g"),
    }
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return "mv05ba_persim_result_v1:" + hashlib.sha256(encoded).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--intervals", required=True)
    parser.add_argument("--references", required=True)
    parser.add_argument("--speed-panel", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    intervals = read_rows(args.intervals)
    references = read_rows(args.references)
    panel = read_rows(args.speed_panel)
    values = defaultdict(list)
    for row in intervals:
        values[(row["diagram_id"], int(row["dimension"]))].append(
            (float(row["birth"]), float(row["death"]))
        )

    build_rows = []
    result_rows = []
    # Retain critical-pair objects for only one sample pair. All-strata and
    # one-stratum prototypes exceeded the frozen 2 GiB gate on deep gene H1.
    pair_groups = defaultdict(list)
    for row in references:
        pair_groups[(row["stratum_id"], int(row["pair_order"]))].append(row)
    for pair_key in sorted(pair_groups):
        pair_rows = pair_groups[pair_key]
        required = sorted({
            (diagram_id, int(row["dimension"][1:]))
            for row in pair_rows
            for diagram_id in (row["first_diagram_id"], row["second_diagram_id"])
        })
        landscapes = {}
        for diagram_id, degree in required:
            build_started = time.perf_counter()
            landscape = build_landscape(values[(diagram_id, degree)], degree)
            elapsed = time.perf_counter() - build_started
            landscapes[(diagram_id, degree)] = landscape
            critical = [] if landscape is None else landscape.critical_pairs
            build_rows.append({
                "contract_id": "mv05ba_private_build_v1", "engine_id": ENGINE_ID,
                "stratum_id": pair_key[0], "pair_order": pair_key[1],
                "diagram_id": diagram_id,
                "dimension": f"H{degree}",
                "finite_intervals": len(values[(diagram_id, degree)]),
                "active_levels": len(critical),
                "critical_points": sum(len(level) for level in critical),
                "build_seconds": format(elapsed, ".17g"),
            })
        for row in pair_rows:
            degree = int(row["dimension"][1:])
            pair_started = time.perf_counter()
            squared, levels, segments, points = corrected_squared(
                landscapes[(row["first_diagram_id"], degree)],
                landscapes[(row["second_diagram_id"], degree)],
            )
            pair_seconds = time.perf_counter() - pair_started
            reference = float(row["reference_squared_distance"])
            error = abs(squared - reference)
            if row["reference_exact"].lower() == "true":
                threshold = 1e-10 + 1e-10 * abs(reference)
                comparison_class = "exact_reference"
            else:
                threshold = float(row["achieved_absolute_error_estimate"]) + \
                    100.0 * sys.float_info.epsilon * max(1.0, abs(reference))
                comparison_class = "adaptive_certificate"
            result_rows.append({
                "contract_id": "mv05ba_private_equivalence_v1",
                "engine_id": ENGINE_ID, "stratum_id": row["stratum_id"],
                "pair_order": row["pair_order"],
                "pair_cache_key": row["pair_cache_key"],
                "dimension": row["dimension"],
                "comparison_class": comparison_class,
                "reference_squared_distance": row["reference_squared_distance"],
                "candidate_squared_distance": format(squared, ".17g"),
                "absolute_error": format(error, ".17g"),
                "acceptance_threshold": format(threshold, ".17g"),
                "equivalent": error <= threshold,
                "active_levels": levels, "segments": segments,
                "critical_points": points,
                "pair_seconds": format(pair_seconds, ".17g"),
                "result_identity": stable_identity(row, squared),
            })
        del landscapes
        release_memory()

    panel_rows = []
    reference_by_pair = defaultdict(dict)
    for row in references:
        reference_by_pair[(row["stratum_id"], row["pair_order"])][row["dimension"]] = row
    for item in panel:
        key = (item["stratum_id"], item["pair_order"])
        fresh_started = time.perf_counter()
        outputs = {}
        for degree in (0, 1):
            ref = reference_by_pair[key][f"H{degree}"]
            first = build_landscape(values[(ref["first_diagram_id"], degree)], degree)
            second = build_landscape(values[(ref["second_diagram_id"], degree)], degree)
            outputs[f"H{degree}"] = corrected_squared(first, second)[0]
        elapsed = time.perf_counter() - fresh_started
        panel_rows.append({
            "contract_id": "mv05ba_private_python_speed_v1",
            "panel_order": item["panel_order"], "stratum_id": item["stratum_id"],
            "pair_order": item["pair_order"],
            "candidate_h0_squared": format(outputs["H0"], ".17g"),
            "candidate_h1_squared": format(outputs["H1"], ".17g"),
            "fresh_pair_seconds": format(elapsed, ".17g"),
        })

    fixtures = fixture_rows()
    os.makedirs(args.output_dir, exist_ok=True)
    write_rows(os.path.join(args.output_dir, "private-build-metrics.csv"), build_rows)
    write_rows(os.path.join(args.output_dir, "private-equivalence.csv"), result_rows)
    write_rows(os.path.join(args.output_dir, "private-python-speed.csv"), panel_rows)
    write_rows(os.path.join(args.output_dir, "fixture-validation.csv"), fixtures)
    environment = [{
        "contract_id": "mv05ba_environment_v1", "engine_id": ENGINE_ID,
        "python": sys.version.replace("\n", " "), "platform": platform.platform(),
        "persim": metadata.version("persim"), "numpy": metadata.version("numpy"),
        "scipy": metadata.version("scipy"),
        "elapsed_seconds": format(time.perf_counter() - started, ".17g"),
        "max_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        "dimension_results": len(result_rows), "speed_panel_pairs": len(panel_rows),
        "retention_policy": "one_pair_at_a_time_malloc_trim_v1",
        "labels_opened": False, "outcomes_computed": False,
    }]
    write_rows(os.path.join(args.output_dir, "environment.csv"), environment)
    if not all(row["passed"] for row in fixtures):
        raise RuntimeError("Corrected-Persim fixture gate failed.")
    if not all(row["equivalent"] for row in result_rows):
        raise RuntimeError("Corrected-Persim MV5-AY equivalence gate failed.")
    print(f"MV5-BA corrected Persim passed {len(result_rows)} dimension results")


if __name__ == "__main__":
    main()
