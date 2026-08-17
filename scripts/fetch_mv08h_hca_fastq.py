#!/usr/bin/env python3
"""Resumably acquire the exact MV8-H HCA FASTQ manifest.

Files are downloaded serially, resumed only with a validated HTTP Range
response, SHA-256 checked, and atomically published. Signed redirect URLs are
never printed or persisted.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import shutil
import sys
import time
import urllib.error
import urllib.request


EXPECTED_CONTRACT = "mv08h_exact_hca_fastq_manifest_v1"
EXPECTED_MANIFEST_SHA256 = "b315e097ef9b2f8cf1cad31cda08c4a440919641f9fff1674cf09e86928b2136"
EXPECTED_FILES = 48
EXPECTED_BYTES = 85_034_239_918
MINIMUM_POST_DOWNLOAD_FREE_BYTES = 1_649_267_441_664
CACHE_CAP_BYTES = 92_274_688_000
MAX_ATTEMPTS = 5
AUTHORIZATION_SENTENCE = (
    "MV8-H exact 48-file HCA FASTQ download authorized by the project owner on 2026-08-17."
)


class IdentityError(RuntimeError):
    """A non-retryable content or manifest identity failure."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path: Path) -> list[dict[str, str]]:
    if sha256_file(path) != EXPECTED_MANIFEST_SHA256:
        raise IdentityError("MV8-H public FASTQ manifest SHA-256 changed")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {
        "contract_id", "file_order", "unit_id", "file_name", "file_uuid",
        "file_version", "file_size_bytes", "file_sha256", "azul_download_url",
        "download_authorized",
    }
    if len(rows) != EXPECTED_FILES or not rows or not required.issubset(rows[0]):
        raise IdentityError("MV8-H manifest schema or cardinality changed")
    rows.sort(key=lambda row: int(row["file_order"]))
    if [int(row["file_order"]) for row in rows] != list(range(1, 49)):
        raise IdentityError("MV8-H file order is not exactly 1 through 48")
    if any(row["contract_id"] != EXPECTED_CONTRACT for row in rows):
        raise IdentityError("MV8-H manifest contract identity changed")
    if any(row["download_authorized"] != "TRUE" for row in rows):
        raise IdentityError("MV8-H manifest contains an unauthorized row")
    if sum(int(row["file_size_bytes"]) for row in rows) != EXPECTED_BYTES:
        raise IdentityError("MV8-H manifest total byte count changed")
    if len({row["file_uuid"] for row in rows}) != EXPECTED_FILES:
        raise IdentityError("MV8-H file UUIDs are not unique")
    for row in rows:
        if Path(row["file_name"]).name != row["file_name"]:
            raise IdentityError("manifest file_name is not a safe basename")
        if not row["file_name"].endswith(".fastq.gz"):
            raise IdentityError("manifest filename is outside the FASTQ contract")
        if len(row["file_sha256"]) != 64:
            raise IdentityError("manifest SHA-256 is malformed")
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


def cache_bytes(cache_dir: Path) -> int:
    return sum(path.stat().st_size for path in cache_dir.rglob("*") if path.is_file())


def local_state(row: dict[str, str], cache_dir: Path) -> tuple[Path, Path, int]:
    target = cache_dir / row["unit_id"] / row["file_name"]
    partial = target.with_name(target.name + ".partial")
    if target.exists() and partial.exists():
        raise IdentityError(f"both final and partial files exist: order {row['file_order']}")
    expected = int(row["file_size_bytes"])
    present = target.stat().st_size if target.exists() else partial.stat().st_size if partial.exists() else 0
    if present > expected:
        raise IdentityError(f"local payload exceeds expected bytes: order {row['file_order']}")
    return target, partial, present


def preflight(rows: list[dict[str, str]], cache_dir: Path) -> tuple[int, int]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    remaining = 0
    for row in rows:
        target, partial, present = local_state(row, cache_dir)
        if target.exists() and present == int(row["file_size_bytes"]):
            continue
        remaining += int(row["file_size_bytes"]) - present
    used = cache_bytes(cache_dir)
    if used > CACHE_CAP_BYTES or used + remaining > CACHE_CAP_BYTES:
        raise RuntimeError("MV8-H cache cap would be exceeded")
    free = shutil.disk_usage(cache_dir).free
    if free - remaining < MINIMUM_POST_DOWNLOAD_FREE_BYTES:
        raise RuntimeError("MV8-H 1.5-TiB post-download free-space reserve would be violated")
    return remaining, free


def validate_or_download(row: dict[str, str], cache_dir: Path) -> dict[str, object]:
    target, partial, present = local_state(row, cache_dir)
    target.parent.mkdir(parents=True, exist_ok=True)
    expected_bytes = int(row["file_size_bytes"])
    expected_sha = row["file_sha256"].lower()
    disposition = ""

    if target.exists():
        observed_sha = sha256_file(target)
        if target.stat().st_size != expected_bytes or observed_sha != expected_sha:
            raise IdentityError(f"incompatible final cache file retained: order {row['file_order']}")
        disposition = "verified_cache_hit"
    else:
        if partial.exists() and partial.stat().st_size == expected_bytes:
            observed_sha = sha256_file(partial)
            if observed_sha != expected_sha:
                raise IdentityError(f"complete partial has wrong SHA-256 and was retained: order {row['file_order']}")
            os.replace(partial, target)
            disposition = "verified_complete_partial_atomic"
        else:
            for attempt in range(1, MAX_ATTEMPTS + 1):
                offset = partial.stat().st_size if partial.exists() else 0
                headers = {"User-Agent": "scPHcompare-MV8H/1.0"}
                if offset:
                    headers["Range"] = f"bytes={offset}-"
                request = urllib.request.Request(row["azul_download_url"], headers=headers)
                try:
                    with urllib.request.urlopen(request, timeout=180) as response:
                        status = getattr(response, "status", response.getcode())
                        if offset and status != 206:
                            raise IdentityError(
                                f"server did not honor HTTP Range; partial retained: order {row['file_order']}"
                            )
                        if offset:
                            content_range = response.headers.get("Content-Range", "")
                            if not content_range.startswith(f"bytes {offset}-"):
                                raise IdentityError(
                                    f"server returned incompatible Content-Range: order {row['file_order']}"
                                )
                        mode = "ab" if offset else "xb"
                        with partial.open(mode) as handle:
                            while block := response.read(8 * 1024 * 1024):
                                handle.write(block)
                                if handle.tell() > expected_bytes:
                                    raise IdentityError(
                                        f"download exceeded expected size; partial retained: order {row['file_order']}"
                                    )
                            handle.flush()
                            os.fsync(handle.fileno())
                    observed_bytes = partial.stat().st_size
                    if observed_bytes != expected_bytes:
                        raise urllib.error.ContentTooShortError(
                            f"incomplete file order {row['file_order']}", None
                        )
                    break
                except IdentityError:
                    raise
                except Exception as error:
                    if attempt == MAX_ATTEMPTS:
                        raise RuntimeError(
                            f"bounded network attempts exhausted; partial retained: order {row['file_order']}"
                        ) from error
                    print(
                        f"MV8-H transient acquisition interruption: order {row['file_order']}; "
                        f"attempt {attempt}/{MAX_ATTEMPTS}; partial retained",
                        flush=True,
                    )
                    time.sleep(2 ** attempt)
            observed_sha = sha256_file(partial)
            if observed_sha != expected_sha:
                raise IdentityError(f"download SHA-256 failed; partial retained: order {row['file_order']}")
            os.replace(partial, target)
            disposition = "downloaded_verified_atomic"

    return {
        "contract_id": "mv08h_private_download_receipt_v1",
        "file_order": int(row["file_order"]),
        "unit_id": row["unit_id"],
        "file_name": row["file_name"],
        "file_uuid": row["file_uuid"],
        "file_version": row["file_version"],
        "expected_bytes": expected_bytes,
        "observed_bytes": target.stat().st_size,
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
    parser.add_argument("--max-files", type=int, default=EXPECTED_FILES)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    authorization = args.authorization_contract.read_text(encoding="utf-8")
    if AUTHORIZATION_SENTENCE not in authorization:
        raise RuntimeError("owner download authorization is not frozen in the supplied contract")
    rows = read_manifest(args.manifest)
    if args.max_files < 1 or args.max_files > EXPECTED_FILES:
        raise ValueError("--max-files must be between 1 and 48")
    remaining, free = preflight(rows, args.cache_dir)
    print(
        f"MV8-H preflight: {len(rows)} files; remaining_bytes={remaining}; "
        f"free_bytes={free}; reserve_bytes={MINIMUM_POST_DOWNLOAD_FREE_BYTES}",
        flush=True,
    )
    if args.dry_run:
        print("MV8-H dry-run complete: no network payload requested", flush=True)
        return 0

    receipts: list[dict[str, object]] = []
    if args.receipt.exists():
        with args.receipt.open("r", encoding="utf-8-sig", newline="") as handle:
            prior = list(csv.DictReader(handle))
        if any(row.get("contract_id") != "mv08h_private_download_receipt_v1" for row in prior):
            raise IdentityError("existing private receipt has incompatible identity")
    for row in rows[:args.max_files]:
        print(
            f"MV8-H acquisition {row['file_order']}/{EXPECTED_FILES}: "
            f"{row['unit_id']} {row['file_name']}",
            flush=True,
        )
        receipts.append(validate_or_download(row, args.cache_dir))
        write_receipt(args.receipt, receipts)
    print(f"MV8-H acquisition gate complete: {len(receipts)} verified files", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"MV8-H acquisition stopped: {error}", file=sys.stderr)
        raise
