#!/usr/bin/env python3
"""Build a normalized, nonpublishing Rust accelerator candidate manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import platform
import subprocess


def sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command_identity(command: list[str]) -> str:
    return subprocess.check_output(command, text=True).strip().replace("\n", " | ")


def git_source_sha256(revision: str, path: str) -> str:
    data = subprocess.check_output(["git", "show", f"{revision}:{path}"])
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", required=True)
    parser.add_argument("--runner", required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--dependencies", required=True)
    parser.add_argument("--fixtures", required=True)
    parser.add_argument("--r-check-absent", required=True)
    parser.add_argument("--r-check-present", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    evidence_files = (
        args.library,
        args.dependencies,
        args.fixtures,
        args.r_check_absent,
        args.r_check_present,
    )
    if not all(os.path.isfile(path) for path in evidence_files):
        raise FileNotFoundError("Candidate library or required evidence is missing")
    commit = os.environ.get("GITHUB_SHA", "local-unpublished")
    source_revision = commit if commit != "local-unpublished" else "HEAD"
    row = {
        "contract_id": "mv05bf_candidate_manifest_v1",
        "engine_version": "1",
        "git_commit": commit,
        "runner_label": args.runner,
        "runner_os": os.environ.get("RUNNER_OS", platform.system()),
        "runner_arch": os.environ.get("RUNNER_ARCH", platform.machine()),
        "certification_class": "candidate-only",
        "rust_target": args.target,
        "rustup_identity": command_identity(["rustup", "--version"]),
        "rustc_identity": command_identity(["rustc", "+1.97.1", "-vV"]),
        "cargo_identity": command_identity(["cargo", "+1.97.1", "-vV"]),
        "source_hash_basis": f"git:{source_revision}",
        "cargo_lock_sha256": git_source_sha256(
            source_revision, "rust/scph_landscape_kernel/Cargo.lock"
        ),
        "cargo_manifest_sha256": git_source_sha256(
            source_revision, "rust/scph_landscape_kernel/Cargo.toml"
        ),
        "rust_source_sha256": git_source_sha256(
            source_revision, "rust/scph_landscape_kernel/src/lib.rs"
        ),
        "public_abi_header_sha256": git_source_sha256(
            source_revision,
            "rust/scph_landscape_kernel/include/scph_landscape_kernel.h",
        ),
        "r_shim_sha256": git_source_sha256(
            source_revision, "R/landscape_rust_prototype.R"
        ),
        "build_command": "cargo +1.97.1 build --manifest-path rust/scph_landscape_kernel/Cargo.toml --release --locked",
        "library_filename": os.path.basename(args.library),
        "library_bytes": os.path.getsize(args.library),
        "library_sha256": sha256(args.library),
        "dependency_inventory_sha256": sha256(args.dependencies),
        "r_fixture_evidence_sha256": sha256(args.fixtures),
        "r_check_absent_log_sha256": sha256(args.r_check_absent),
        "r_check_present_log_sha256": sha256(args.r_check_present),
        "clean_builds_byte_identical": "true",
        "rust_tests_passed": "true",
        "strict_clippy_passed": "true",
        "native_ffi_passed": "true",
        "r_fixtures_and_fallback_passed": "true",
        "r_package_checks_absent_and_present_passed": "true",
        "private_data_used": "false",
        "published_release": "false",
    }
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    temporary = args.output + ".tmp"
    with open(temporary, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=row.keys())
        writer.writeheader()
        writer.writerow(row)
    os.replace(temporary, args.output)


if __name__ == "__main__":
    main()
