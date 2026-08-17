#!/usr/bin/env python3
"""Validate and freeze the authorized MV8-H Ensembl-93 FASTA input."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
from pathlib import Path
import subprocess


FASTA_NAME = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
FASTA_URL = (
    "https://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/"
    + FASTA_NAME
)
EXPECTED_FASTA_BYTES = 881_214_682
EXPECTED_FASTA_SHA256 = (
    "2a27436d44f0d6350f86894fbe5edec56faa5467028879784508041562406aa0"
)
EXPECTED_BSD_SUM = "07119"
EXPECTED_BSD_BLOCKS = 860_562
EXPECTED_CHECKSUMS_SHA256 = (
    "d0cd8313040eebcb42db28360477278579ed1e7b14aef0d6b119626f53c12c11"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), extrasaction="raise",
            quoting=csv.QUOTE_ALL, lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def official_checksum_line(path: Path) -> str:
    matches = [
        line.strip() for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip().endswith(" " + FASTA_NAME)
    ]
    if matches != [f"{EXPECTED_BSD_SUM} {EXPECTED_BSD_BLOCKS} {FASTA_NAME}"]:
        raise RuntimeError("official Ensembl CHECKSUMS record differs from the freeze")
    return matches[0]


def local_bsd_sum(path: Path) -> tuple[str, int]:
    result = subprocess.run(
        ["sum", str(path)], check=True, text=True, capture_output=True,
    )
    fields = result.stdout.strip().split()
    if len(fields) < 2:
        raise RuntimeError("BSD sum output is incomplete")
    return fields[0], int(fields[1])


def validate_gzip(path: Path) -> int:
    uncompressed_bytes = 0
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            uncompressed_bytes += len(block)
    return uncompressed_bytes


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--checksums", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.fasta.name != FASTA_NAME or not args.fasta.is_file():
        raise RuntimeError("exact Ensembl-93 primary-assembly FASTA is absent")
    if not args.checksums.is_file():
        raise RuntimeError("official Ensembl CHECKSUMS file is absent")

    observed_bytes = args.fasta.stat().st_size
    observed_sha256 = sha256_file(args.fasta)
    checksums_sha256 = sha256_file(args.checksums)
    checksum_line = official_checksum_line(args.checksums)
    observed_bsd_sum, observed_bsd_blocks = local_bsd_sum(args.fasta)
    uncompressed_bytes = validate_gzip(args.fasta)

    gates = [
        observed_bytes == EXPECTED_FASTA_BYTES,
        observed_sha256 == EXPECTED_FASTA_SHA256,
        checksums_sha256 == EXPECTED_CHECKSUMS_SHA256,
        checksum_line == f"{EXPECTED_BSD_SUM} {EXPECTED_BSD_BLOCKS} {FASTA_NAME}",
        observed_bsd_sum == EXPECTED_BSD_SUM,
        observed_bsd_blocks == EXPECTED_BSD_BLOCKS,
        uncompressed_bytes > observed_bytes,
    ]
    if not all(gates):
        raise RuntimeError("one or more Ensembl-93 reference-input gates failed")

    identity = [{
        "contract_id": "mv08h_ensembl93_reference_input_identity_v1",
        "source_url": FASTA_URL,
        "file_name": FASTA_NAME,
        "expected_bytes": EXPECTED_FASTA_BYTES,
        "observed_bytes": observed_bytes,
        "expected_sha256": EXPECTED_FASTA_SHA256,
        "observed_sha256": observed_sha256,
        "official_bsd_sum": EXPECTED_BSD_SUM,
        "observed_bsd_sum": observed_bsd_sum,
        "official_bsd_blocks": EXPECTED_BSD_BLOCKS,
        "observed_bsd_blocks": observed_bsd_blocks,
        "official_checksums_sha256": checksums_sha256,
        "gzip_integrity": "TRUE",
        "uncompressed_bytes": uncompressed_bytes,
        "all_identity_gates_passed": "TRUE",
    }]
    write_csv(args.output_dir / "mv08h-reference-input-identity.csv", identity)

    decision = [{
        "contract_id": "mv08h_reference_input_decision_v1",
        "decision": "ensembl93_reference_input_exact_stop_before_runtime_or_mkref",
        "reference_input_identity_closed": "TRUE",
        "cellranger_acquisition_authorized": "FALSE",
        "cellranger_runtime_validated": "FALSE",
        "mkref_authorized": "FALSE",
        "raw_reprocessing_authorized": "FALSE",
        "label_access_authorized": "FALSE",
        "biological_outcomes_authorized": "FALSE",
        "next_gate": "owner_provided_exact_cellranger_3_0_0_runtime_then_custom_reference_validation",
    }]
    write_csv(args.output_dir / "mv08h-reference-input-decision.csv", decision)

    artifact_names = sorted(
        path.name for path in args.output_dir.iterdir()
        if path.is_file() and not path.name.endswith("artifact-manifest.csv")
    )
    artifacts = [{
        "contract_id": "mv08h_reference_input_artifact_manifest_v1",
        "file": name,
        "bytes": (args.output_dir / name).stat().st_size,
        "sha256": sha256_file(args.output_dir / name),
        "contains_expression": "FALSE",
        "contains_cell_barcode": "FALSE",
        "contains_absolute_private_path": "FALSE",
        "contains_donor_attribute": "FALSE",
        "contains_outcome_label": "FALSE",
    } for name in artifact_names]
    write_csv(args.output_dir / "mv08h-reference-input-artifact-manifest.csv", artifacts)
    print("MV8-H Ensembl-93 reference-input validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
