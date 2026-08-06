#!/usr/bin/env python3
"""Independent Persim exact oracle and candidate scaling profiler."""

import argparse
import csv
import json
import math
import os
import platform
import resource
import subprocess
import sys
import time
from importlib import metadata

import numpy as np
import psutil
from gudhi.representations import Landscape
from persim import PersLandscapeExact
from scipy.integrate import quad


def array(values):
    return np.asarray(values, dtype=float).reshape((-1, 2))


def corpus():
    deep = [[i / 100, 3 - i / 100] for i in range(1, 13)]
    return {
        "single_tent": ({0: [[0, 2]], 1: []}, {0: [], 1: []}),
        "translated_tents": ({0: [[0, 2]], 1: []}, {0: [[3, 5]], 1: []}),
        "overlapping_dimensions": (
            {0: [[0, 4], [1, 3]], 1: [[0, 2]]},
            {0: [[0.5, 3.5], [1.25, 2.75]], 1: [[0.25, 1.75]]},
        ),
        "narrow_feature": ({0: [[0.499, 0.501]], 1: []}, {0: [], 1: []}),
        "deep_landscape": ({0: deep, 1: []}, {0: [], 1: []}),
        "both_dimensions": (
            {0: [[0, 2]], 1: [[0, 1]]},
            {0: [[0.25, 2.25]], 1: [[0.2, 0.9]]},
        ),
    }


def persim_exact_distance(first, second, homology_dimension=0):
    first = array(first)
    second = array(second)
    if len(first) == 0 and len(second) == 0:
        return 0.0

    def landscape(diagram):
        diagrams = [np.empty((0, 2)) for _ in range(homology_dimension + 1)]
        diagrams[homology_dimension] = diagram
        return PersLandscapeExact(dgms=diagrams, hom_deg=homology_dimension)

    if len(first) == 0:
        return float(landscape(second).p_norm(p=2))
    if len(second) == 0:
        return float(landscape(first).p_norm(p=2))
    return float((landscape(first) - landscape(second)).p_norm(p=2))


def persim_corrected_distance(first, second, homology_dimension=0):
    first = array(first)
    second = array(second)
    if len(first) == 0 and len(second) == 0:
        return 0.0

    def landscape(diagram):
        diagrams = [np.empty((0, 2)) for _ in range(homology_dimension + 1)]
        diagrams[homology_dimension] = diagram
        return PersLandscapeExact(dgms=diagrams, hom_deg=homology_dimension)

    if len(first) == 0:
        critical_pairs = landscape(second).critical_pairs
    elif len(second) == 0:
        critical_pairs = landscape(first).critical_pairs
    else:
        critical_pairs = (landscape(first) - landscape(second)).critical_pairs
    squared = 0.0
    for level in critical_pairs:
        for (x0, y0), (x1, y1) in zip(level, level[1:]):
            squared += (x1 - x0) * (y0 * y0 + y0 * y1 + y1 * y1) / 3
    return math.sqrt(max(0.0, squared))


def gudhi_grid_distance(first, second, resolution):
    first = array(first)
    second = array(second)
    if len(first) == 0 and len(second) == 0:
        return 0.0
    nonempty = [x for x in (first, second) if len(x)]
    lower = min(float(np.min(x[:, 0])) for x in nonempty)
    upper = max(float(np.max(x[:, 1])) for x in nonempty)
    levels = max(len(first), len(second))
    transformer = Landscape(
        num_landscapes=levels,
        resolution=resolution,
        sample_range=[lower, upper],
        keep_endpoints=True,
    )
    vectors = transformer.fit_transform([first, second])
    step = (upper - lower) / (resolution - 1)
    # GUDHI multiplies sampled landscapes by sqrt(2). Multiplication by
    # sqrt(step / 2) restores the uniform-grid L2 quadrature scale.
    return float(np.linalg.norm(vectors[0] - vectors[1]) * math.sqrt(step / 2))


def scipy_quad_distance(first, second):
    first = array(first)
    second = array(second)
    if len(first) == 0 and len(second) == 0:
        return 0.0
    nonempty = [x for x in (first, second) if len(x)]
    breaks = sorted(
        {
            float(value)
            for diagram in nonempty
            for value in np.concatenate(
                [diagram[:, 0], np.mean(diagram, axis=1), diagram[:, 1]]
            )
        }
    )

    def levels(diagram, location):
        if len(diagram) == 0:
            return np.empty(0)
        values = np.minimum(location - diagram[:, 0], diagram[:, 1] - location)
        return np.sort(values[values > 0])[::-1]

    def integrand(location):
        a = levels(first, location)
        b = levels(second, location)
        depth = max(len(a), len(b))
        difference = np.zeros(depth)
        difference[: len(a)] += a
        difference[: len(b)] -= b
        return float(np.dot(difference, difference))

    squared = sum(
        quad(integrand, left, right, epsabs=1e-12, epsrel=1e-12, limit=200)[0]
        for left, right in zip(breaks, breaks[1:])
    )
    return math.sqrt(max(0.0, squared))


def write_corpus(path):
    rows = []
    for case, (first, second) in corpus().items():
        for dimension in (0, 1):
            a = first[dimension]
            b = second[dimension]
            start = time.perf_counter()
            exact = persim_exact_distance(a, b, dimension)
            rows.append(
                {
                    "case": case,
                    "dimension": f"H{dimension}",
                    "engine": "persim_exact_0.3.8",
                    "resolution": "",
                    "distance": exact,
                    "elapsed_seconds": time.perf_counter() - start,
                }
            )
            start = time.perf_counter()
            corrected = persim_corrected_distance(a, b, dimension)
            rows.append(
                {
                    "case": case,
                    "dimension": f"H{dimension}",
                    "engine": "persim_0.3.8_corrected_pnorm",
                    "resolution": "",
                    "distance": corrected,
                    "elapsed_seconds": time.perf_counter() - start,
                }
            )
            start = time.perf_counter()
            quadrature = scipy_quad_distance(a, b)
            rows.append(
                {
                    "case": case,
                    "dimension": f"H{dimension}",
                    "engine": "scipy_quad_1.15.3",
                    "resolution": "",
                    "distance": quadrature,
                    "elapsed_seconds": time.perf_counter() - start,
                }
            )
            for resolution in (250, 500, 1000):
                start = time.perf_counter()
                distance = gudhi_grid_distance(a, b, resolution)
                rows.append(
                    {
                        "case": case,
                        "dimension": f"H{dimension}",
                        "engine": "gudhi_grid_3.12.0",
                        "resolution": resolution,
                        "distance": distance,
                        "elapsed_seconds": time.perf_counter() - start,
                    }
                )
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def scaling_diagrams(size):
    first = [[i / (20 * size), 3 - i / (20 * size)] for i in range(1, size + 1)]
    second = [[b + 0.01, d + 0.015] for b, d in first]
    return first, second


def python_worker(candidate, size):
    first, second = scaling_diagrams(size)
    started = time.perf_counter()
    if candidate == "persim_exact":
        distance = persim_exact_distance(first, second, 0)
    elif candidate == "persim_corrected_pnorm":
        distance = persim_corrected_distance(first, second, 0)
    elif candidate == "scipy_quad":
        distance = scipy_quad_distance(first, second)
    elif candidate.startswith("gudhi_grid_"):
        resolution = int(candidate.rsplit("_", 1)[1])
        distance = gudhi_grid_distance(first, second, resolution)
    else:
        raise ValueError(f"Unknown Python candidate: {candidate}")
    elapsed = time.perf_counter() - started
    maximum_rss_bytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    print(f"{distance}\t{elapsed}\t{maximum_rss_bytes}")


def profile_process(command, timeout_seconds):
    started = time.perf_counter()
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    peak_rss = 0
    timed_out = False
    while process.poll() is None:
        try:
            root = psutil.Process(process.pid)
            members = [root] + root.children(recursive=True)
            peak_rss = max(
                peak_rss,
                sum(member.memory_info().rss for member in members if member.is_running()),
            )
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        if time.perf_counter() - started > timeout_seconds:
            timed_out = True
            process.kill()
            break
        time.sleep(0.01)
    stdout, stderr = process.communicate()
    elapsed = time.perf_counter() - started
    if timed_out:
        return None, elapsed, None, peak_rss, "timeout", stderr.strip()
    if process.returncode != 0:
        return None, elapsed, None, peak_rss, "failed", stderr.strip()
    fields = stdout.strip().splitlines()[-1].split("\t")
    return (
        float(fields[0]),
        elapsed,
        float(fields[1]),
        max(peak_rss, int(float(fields[2]))),
        "completed",
        stderr.strip(),
    )


def write_benchmark(path, repo, repetitions, timeout_seconds, sizes):
    rows = []
    python_candidates = (
        "persim_exact", "persim_corrected_pnorm", "scipy_quad",
        "gudhi_grid_250", "gudhi_grid_1000"
    )
    r_candidates = ("exact", "adaptive")
    for size in sizes:
        for candidate in python_candidates:
            for repetition in range(1, repetitions + 1):
                command = [sys.executable, __file__, "--worker", candidate, str(size)]
                distance, elapsed, kernel_elapsed, peak, status, error = profile_process(
                    command, timeout_seconds
                )
                rows.append(
                    {
                        "candidate": candidate,
                        "implementation": "python",
                        "intervals": size,
                        "maximum_active_levels": size,
                        "repetition": repetition,
                        "timeout_seconds": timeout_seconds,
                        "status": status,
                        "distance": distance,
                        "elapsed_seconds": elapsed,
                        "kernel_elapsed_seconds": kernel_elapsed,
                        "peak_process_tree_rss_bytes": peak,
                        "error": error[:500],
                    }
                )
        for candidate in r_candidates:
            for repetition in range(1, repetitions + 1):
                command = [
                    "Rscript",
                    "--vanilla",
                    os.path.join(repo, "scripts", "landscape_reference_worker.R"),
                    candidate,
                    str(size),
                ]
                distance, elapsed, kernel_elapsed, peak, status, error = profile_process(
                    command, timeout_seconds
                )
                rows.append(
                    {
                        "candidate": f"r_{candidate}",
                        "implementation": "R",
                        "intervals": size,
                        "maximum_active_levels": size,
                        "repetition": repetition,
                        "timeout_seconds": timeout_seconds,
                        "status": status,
                        "distance": distance,
                        "elapsed_seconds": elapsed,
                        "kernel_elapsed_seconds": kernel_elapsed,
                        "peak_process_tree_rss_bytes": peak,
                        "error": error[:500],
                    }
                )
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(path, repetitions, timeout_seconds):
    packages = (
        "persim", "gudhi", "numpy", "scipy", "scikit-learn", "psutil"
    )
    row = {
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "logical_cpus": os.cpu_count(),
        "benchmark_repetitions": repetitions,
        "timeout_seconds": timeout_seconds,
    }
    for package in packages:
        row[f"package_{package.replace('-', '_')}"] = metadata.version(package)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=row.keys())
        writer.writeheader()
        writer.writerow(row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus")
    parser.add_argument("--benchmark")
    parser.add_argument("--repo")
    parser.add_argument("--manifest")
    parser.add_argument("--repetitions", type=int, default=2)
    parser.add_argument("--timeout-seconds", type=float, default=120)
    parser.add_argument("--sizes", default="10,25,50,100,200")
    parser.add_argument("--worker", nargs=2, metavar=("CANDIDATE", "SIZE"))
    args = parser.parse_args()
    if args.worker:
        python_worker(args.worker[0], int(args.worker[1]))
    elif args.corpus:
        write_corpus(args.corpus)
    elif args.benchmark:
        if not args.repo:
            parser.error("--repo is required with --benchmark")
        sizes = tuple(int(value) for value in args.sizes.split(","))
        if not sizes or any(value < 1 for value in sizes):
            parser.error("--sizes must contain positive comma-separated integers")
        write_benchmark(
            args.benchmark, args.repo, args.repetitions, args.timeout_seconds, sizes
        )
    elif args.manifest:
        write_manifest(args.manifest, args.repetitions, args.timeout_seconds)
    else:
        parser.error("choose --corpus, --benchmark, --manifest, or --worker")


if __name__ == "__main__":
    main()
