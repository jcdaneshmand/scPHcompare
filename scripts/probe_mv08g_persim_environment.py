#!/usr/bin/env python3
"""Emit a path-safe dependency fingerprint for the MV8-G Persim oracle."""

import argparse
import csv
import hashlib
import importlib.metadata
import importlib.util
import os
import sys


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine-source", required=True)
    args = parser.parse_args()

    import numpy
    import persim

    spec = importlib.util.spec_from_file_location("mv08g_probe_engine", args.engine_source)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load the corrected Persim engine.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "landscape_distance"):
        raise RuntimeError("Corrected Persim engine lacks landscape_distance.")

    fields = [
        "python_version",
        "environment_name",
        "persim_version",
        "persim_init_sha256",
        "numpy_version",
        "engine_import_passed",
    ]
    writer = csv.DictWriter(sys.stdout, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerow({
        "python_version": sys.version.split()[0],
        "environment_name": os.path.basename(sys.prefix),
        "persim_version": importlib.metadata.version("persim"),
        "persim_init_sha256": sha256(persim.__file__),
        "numpy_version": numpy.__version__,
        "engine_import_passed": "TRUE",
    })


if __name__ == "__main__":
    main()
