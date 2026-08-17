#!/usr/bin/env python3
"""Independently validate an MV8-H sentinel or complete FASTQ cache."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import shutil


EXPECTED_MANIFEST_SHA256 = "b315e097ef9b2f8cf1cad31cda08c4a440919641f9fff1674cf09e86928b2136"
EXPECTED_TOTAL_BYTES = 85_034_239_918
RESERVE_BYTES = 1_649_267_441_664
CACHE_CAP_BYTES = 92_274_688_000


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--expected-files", type=int, required=True)
    parser.add_argument("--stage", choices=("sentinel", "complete"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if sha256_file(args.manifest) != EXPECTED_MANIFEST_SHA256:
        raise RuntimeError("MV8-H manifest SHA-256 changed")
    manifest = read_csv(args.manifest)
    if len(manifest) != 48 or sum(int(row["file_size_bytes"]) for row in manifest) != EXPECTED_TOTAL_BYTES:
        raise RuntimeError("MV8-H manifest cardinality or total changed")
    if args.expected_files not in (1, 48):
        raise ValueError("only the one-file sentinel or complete 48-file cache may validate")
    if args.stage == "sentinel" and args.expected_files != 1:
        raise ValueError("sentinel stage requires exactly one expected file")
    if args.stage == "complete" and args.expected_files != 48:
        raise ValueError("complete stage requires exactly 48 expected files")

    receipt = read_csv(args.receipt)
    if len(receipt) != args.expected_files:
        raise RuntimeError("private receipt cardinality differs from the requested gate")
    receipt_by_order = {int(row["file_order"]): row for row in receipt}
    if sorted(receipt_by_order) != list(range(1, args.expected_files + 1)):
        raise RuntimeError("private receipt order is incomplete or duplicated")

    final_files = sorted(args.cache_dir.rglob("*.fastq.gz"))
    partial_files = sorted(args.cache_dir.rglob("*.partial"))
    if len(final_files) != args.expected_files or partial_files:
        raise RuntimeError("cache final/partial file counts differ from the requested gate")

    validation_rows: list[dict[str, object]] = []
    for row in manifest[:args.expected_files]:
        order = int(row["file_order"])
        target = args.cache_dir / row["unit_id"] / row["file_name"]
        if not target.is_file() or target not in final_files:
            raise RuntimeError(f"expected final file is absent: order {order}")
        observed_bytes = target.stat().st_size
        observed_sha = sha256_file(target)
        with target.open("rb") as handle:
            gzip_magic = handle.read(2) == b"\x1f\x8b"
        receipt_row = receipt_by_order[order]
        receipt_exact = all([
            receipt_row["contract_id"] == "mv08h_private_download_receipt_v1",
            receipt_row["file_uuid"] == row["file_uuid"],
            receipt_row["file_version"] == row["file_version"],
            receipt_row["expected_bytes"] == row["file_size_bytes"],
            receipt_row["observed_bytes"] == row["file_size_bytes"],
            receipt_row["expected_sha256"] == row["file_sha256"],
            receipt_row["observed_sha256"] == row["file_sha256"],
            receipt_row["verified"] == "TRUE",
        ])
        passed = observed_bytes == int(row["file_size_bytes"]) and observed_sha == row["file_sha256"] and gzip_magic and receipt_exact
        validation_rows.append({
            "contract_id": f"mv08h_{args.stage}_file_validation_v1",
            "file_order": order,
            "unit_id": row["unit_id"],
            "file_name": row["file_name"],
            "file_uuid": row["file_uuid"],
            "expected_bytes": int(row["file_size_bytes"]),
            "observed_bytes": observed_bytes,
            "expected_sha256": row["file_sha256"],
            "observed_sha256": observed_sha,
            "gzip_magic": str(gzip_magic).upper(),
            "receipt_exact": str(receipt_exact).upper(),
            "passed": str(passed).upper(),
        })
    if not all(row["passed"] == "TRUE" for row in validation_rows):
        raise RuntimeError("one or more MV8-H cache files failed independent identity validation")

    observed_cache_bytes = sum(path.stat().st_size for path in final_files)
    expected_cache_bytes = sum(int(row["file_size_bytes"]) for row in manifest[:args.expected_files])
    remaining_bytes = EXPECTED_TOTAL_BYTES - expected_cache_bytes
    free_bytes = shutil.disk_usage(args.cache_dir).free
    storage_pass = observed_cache_bytes == expected_cache_bytes and observed_cache_bytes <= CACHE_CAP_BYTES
    storage_pass = storage_pass and free_bytes - remaining_bytes >= RESERVE_BYTES
    if not storage_pass:
        raise RuntimeError("MV8-H cache validation failed the byte cap or free-space reserve")

    prefix = f"mv08h-{args.stage}"
    write_csv(args.output_dir / f"{prefix}-file-validation.csv", validation_rows)
    summary_rows = [{
        "contract_id": f"mv08h_{args.stage}_summary_v1",
        "expected_files": args.expected_files,
        "verified_files": len(validation_rows),
        "expected_cache_bytes": expected_cache_bytes,
        "observed_cache_bytes": observed_cache_bytes,
        "remaining_manifest_bytes": remaining_bytes,
        "partial_files": len(partial_files),
        "free_bytes_at_validation": free_bytes,
        "post_remaining_free_bytes": free_bytes - remaining_bytes,
        "required_reserve_bytes": RESERVE_BYTES,
        "storage_gate_passed": str(storage_pass).upper(),
        "all_file_gates_passed": "TRUE",
    }]
    write_csv(args.output_dir / f"{prefix}-summary.csv", summary_rows)
    if args.stage == "sentinel":
        decision = "sentinel_exact_authorize_remaining_47_file_download"
        remaining_authorized = "TRUE"
        reference_next = "FALSE"
        next_gate = "complete_48_file_checksum_closure"
    else:
        decision = "complete_download_exact_authorize_reference_input_validation_only"
        remaining_authorized = "FALSE"
        reference_next = "TRUE"
        next_gate = "exact_cellranger_3_0_0_and_custom_reference_validation"
    decision_rows = [{
        "contract_id": f"mv08h_{args.stage}_decision_v1",
        "decision": decision,
        "remaining_fastq_download_authorized": remaining_authorized,
        "reference_input_validation_authorized": reference_next,
        "mkref_authorized": "FALSE",
        "raw_reprocessing_authorized": "FALSE",
        "label_access_authorized": "FALSE",
        "biological_outcomes_authorized": "FALSE",
        "next_gate": next_gate,
    }]
    write_csv(args.output_dir / f"{prefix}-decision.csv", decision_rows)

    artifact_names = sorted(
        path.name for path in args.output_dir.iterdir()
        if path.is_file() and not path.name.endswith("artifact-manifest.csv")
    )
    artifact_rows = [{
        "contract_id": f"mv08h_{args.stage}_artifact_manifest_v1",
        "file": name,
        "bytes": (args.output_dir / name).stat().st_size,
        "sha256": sha256_file(args.output_dir / name),
        "contains_expression": "FALSE",
        "contains_cell_barcode": "FALSE",
        "contains_absolute_private_path": "FALSE",
        "contains_donor_attribute": "FALSE",
        "contains_outcome_label": "FALSE",
    } for name in artifact_names]
    write_csv(args.output_dir / f"{prefix}-artifact-manifest.csv", artifact_rows)
    print(f"MV8-H {args.stage} cache validation passed: {len(validation_rows)}/{args.expected_files} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
