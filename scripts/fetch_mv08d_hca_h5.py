#!/usr/bin/env python3
"""Fetch the exact MV8-C HCA H5 manifest into a private cache.

The script is fail-closed: a present but incompatible file is never replaced,
and a download is atomically published only after byte and SHA-256 validation.
Transient redirected URLs are intentionally never logged.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import sys
import urllib.request


EXPECTED_CONTRACT = "mv08c_hca_primary_file_manifest_v1"
EXPECTED_TOTAL_BYTES = 202_770_089
EXPECTED_ROWS = 8
AUTHORIZATION_SENTENCE = (
    "Exact eight-file structural/QC download authorized"
)


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {
        "contract_id", "file_order", "unit_id", "file_name", "file_uuid",
        "file_version", "file_size_bytes", "file_sha256",
        "azul_download_url",
    }
    if len(rows) != EXPECTED_ROWS or not rows or not required.issubset(rows[0]):
        raise ValueError("MV8-C exact-file manifest has an invalid schema or row count")
    rows.sort(key=lambda row: int(row["file_order"]))
    if [int(row["file_order"]) for row in rows] != list(range(1, 9)):
        raise ValueError("MV8-C file order is not exactly 1 through 8")
    if any(row["contract_id"] != EXPECTED_CONTRACT for row in rows):
        raise ValueError("MV8-C manifest contract identity changed")
    if sum(int(row["file_size_bytes"]) for row in rows) != EXPECTED_TOTAL_BYTES:
        raise ValueError("MV8-C manifest total byte count changed")
    if len({row["unit_id"] for row in rows}) != 8:
        raise ValueError("MV8-C biological unit identities are not unique")
    if len({row["file_uuid"] for row in rows}) != 8:
        raise ValueError("MV8-C file identities are not unique")
    for row in rows:
        if Path(row["file_name"]).name != row["file_name"]:
            raise ValueError("Manifest file_name is not a safe basename")
        if not row["file_name"].endswith("_raw_feature_bc_matrix.h5"):
            raise ValueError("Manifest filename is outside the frozen H5 contract")
        if len(row["file_sha256"]) != 64:
            raise ValueError("Manifest SHA-256 is malformed")
    return rows


def write_receipt(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "contract_id", "file_order", "unit_id", "file_name", "file_uuid",
        "file_version", "expected_bytes", "observed_bytes",
        "expected_sha256", "observed_sha256", "disposition", "verified",
    ]
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def fetch_one(row: dict[str, str], cache_dir: Path) -> dict[str, object]:
    target = cache_dir / row["unit_id"] / row["file_name"]
    target.parent.mkdir(parents=True, exist_ok=True)
    expected_bytes = int(row["file_size_bytes"])
    expected_sha = row["file_sha256"].lower()

    if target.exists():
        observed_bytes = target.stat().st_size
        observed_sha = sha256_file(target)
        if observed_bytes != expected_bytes or observed_sha != expected_sha:
            raise RuntimeError(
                f"existing cache file is incompatible and was not replaced: {row['unit_id']}"
            )
        disposition = "verified_cache_hit"
    else:
        partial = target.with_name(target.name + ".partial")
        if partial.exists():
            raise RuntimeError(
                f"pre-existing partial file requires manual review: {row['unit_id']}"
            )
        request = urllib.request.Request(
            row["azul_download_url"],
            headers={"User-Agent": "scPHcompare-MV8D/1.0"},
        )
        digest = hashlib.sha256()
        observed_bytes = 0
        try:
            with urllib.request.urlopen(request, timeout=120) as response, partial.open("xb") as handle:
                while chunk := response.read(1024 * 1024):
                    handle.write(chunk)
                    digest.update(chunk)
                    observed_bytes += len(chunk)
                handle.flush()
                os.fsync(handle.fileno())
        except Exception:
            # Preserve a partial payload for forensic review; never open or publish it.
            raise
        observed_sha = digest.hexdigest()
        if observed_bytes != expected_bytes or observed_sha != expected_sha:
            raise RuntimeError(
                f"download identity failed and partial was retained: {row['unit_id']}"
            )
        os.replace(partial, target)
        disposition = "downloaded_verified_atomic"

    return {
        "contract_id": "mv08d_hca_private_download_receipt_v1",
        "file_order": int(row["file_order"]),
        "unit_id": row["unit_id"],
        "file_name": row["file_name"],
        "file_uuid": row["file_uuid"],
        "file_version": row["file_version"],
        "expected_bytes": expected_bytes,
        "observed_bytes": observed_bytes,
        "expected_sha256": expected_sha,
        "observed_sha256": observed_sha,
        "disposition": disposition,
        "verified": "TRUE",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--authorization-contract", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()

    authorization = args.authorization_contract.read_text(encoding="utf-8")
    if AUTHORIZATION_SENTENCE not in authorization:
        raise RuntimeError("owner download authorization is not frozen in the supplied contract")
    rows = read_manifest(args.manifest)
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    receipts: list[dict[str, object]] = []
    for row in rows:
        print(f"MV8-D verified acquisition {row['file_order']}/8: {row['unit_id']}", flush=True)
        receipts.append(fetch_one(row, args.cache_dir))
        write_receipt(args.receipt, receipts)
    print(f"MV8-D exact-file acquisition complete: {len(receipts)} files", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"MV8-D acquisition stopped: {error}", file=sys.stderr)
        raise
