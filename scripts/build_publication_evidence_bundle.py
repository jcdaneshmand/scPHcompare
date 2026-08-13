#!/usr/bin/env python3
"""Build or verify a deterministic DOI-archive evidence bundle from Git blobs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath


CONTRACT_ID = "scphcompare_publication_evidence_v1"
FIXED_ZIP_TIME = (1980, 1, 1, 0, 0, 0)


@dataclass(frozen=True)
class EvidenceFile:
    path: str
    category: str
    media_type: str
    size: int
    git_blob: str


def run_git(repo: Path, *args: str, text: bool = True) -> str | bytes:
    result = subprocess.run(
        ["git", "-C", str(repo), *args],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=text,
    )
    return result.stdout


def resolve_revision(repo: Path, revision: str) -> str:
    return str(run_git(repo, "rev-parse", f"{revision}^{{commit}}")).strip()


def classify(path: str) -> tuple[str, str] | None:
    if path.startswith("results/"):
        return "derived_r_result", "application/x-r-rds"
    if path.startswith("docs/audits/") and path.endswith(".csv.gz"):
        return "generated_audit_table", "application/gzip"
    if path.startswith("docs/audits/") and path.endswith(".csv"):
        return "generated_audit_table", "text/csv"
    return None


def changed_paths(repo: Path, base: str, revision: str) -> set[str]:
    raw = run_git(
        repo,
        "diff",
        "--name-only",
        "--diff-filter=ACMR",
        "-z",
        f"{base}..{revision}",
        text=False,
    )
    return {item.decode("utf-8") for item in raw.split(b"\0") if item}


def inventory(repo: Path, base: str, revision: str) -> list[EvidenceFile]:
    changed = changed_paths(repo, base, revision)
    raw = run_git(repo, "ls-tree", "-r", "-l", "-z", revision, text=False)
    rows: list[EvidenceFile] = []
    for record in raw.split(b"\0"):
        if not record:
            continue
        metadata, encoded_path = record.split(b"\t", 1)
        mode, kind, blob, encoded_size = metadata.decode("ascii").split()
        path = encoded_path.decode("utf-8")
        classification = classify(path)
        if path not in changed or classification is None:
            continue
        if kind != "blob" or mode == "160000":
            raise RuntimeError(f"Unsupported evidence object: {path}")
        category, media_type = classification
        rows.append(EvidenceFile(path, category, media_type, int(encoded_size), blob))
    return sorted(rows, key=lambda item: item.path.encode("utf-8"))


class GitBlobReader:
    def __init__(self, repo: Path) -> None:
        self.process = subprocess.Popen(
            ["git", "-C", str(repo), "cat-file", "--batch"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    def read(self, oid: str, expected_size: int) -> bytes:
        assert self.process.stdin is not None
        assert self.process.stdout is not None
        self.process.stdin.write(oid.encode("ascii") + b"\n")
        self.process.stdin.flush()
        header = self.process.stdout.readline().decode("ascii").strip().split()
        if len(header) != 3 or header[1] != "blob":
            raise RuntimeError(f"Unable to read Git blob {oid}: {' '.join(header)}")
        size = int(header[2])
        data = self.process.stdout.read(size)
        terminator = self.process.stdout.read(1)
        if size != expected_size or len(data) != size or terminator != b"\n":
            raise RuntimeError(f"Git blob framing or size mismatch for {oid}")
        return data

    def close(self) -> None:
        if self.process.stdin is not None:
            self.process.stdin.close()
        return_code = self.process.wait()
        if return_code != 0:
            assert self.process.stderr is not None
            raise RuntimeError(self.process.stderr.read().decode("utf-8", errors="replace"))

    def __enter__(self) -> "GitBlobReader":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()


MANIFEST_FIELDS = [
    "contract_id",
    "archive_member",
    "source_path",
    "category",
    "media_type",
    "bytes",
    "sha256",
    "git_blob",
    "source_revision",
    "base_revision",
]


def manifest_bytes(rows: list[dict[str, str]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=MANIFEST_FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def zip_info(name: str) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo(name, FIXED_ZIP_TIME)
    info.compress_type = zipfile.ZIP_DEFLATED
    info.create_system = 3
    info.external_attr = 0o100644 << 16
    return info


def archive_readme(source_revision: str, base_revision: str, count: int) -> bytes:
    return (
        "scPHcompare derived evidence bundle\n"
        "====================================\n\n"
        f"Contract: {CONTRACT_ID}\n"
        f"Source revision: {source_revision}\n"
        f"Base revision: {base_revision}\n"
        f"Evidence files: {count}\n\n"
        "This archive contains generated audit tables and derived R result objects.\n"
        "Source code, tests, specifications, and human-readable reports remain in\n"
        "the scPHcompare Git repository. Verify every member against MANIFEST.csv.\n\n"
        "DOI: pending Zenodo draft and author-team metadata review.\n"
    ).encode("utf-8")


def metadata_bytes(source_revision: str, count: int, total_bytes: int) -> bytes:
    metadata = {
        "schema": "scphcompare_zenodo_metadata_draft_v1",
        "status": "draft_owner_and_author_team_review_required",
        "resource_type": "dataset",
        "title": "scPHcompare derived computational evidence",
        "description": (
            "Generated audit tables and derived persistence-homology result "
            "objects supporting the scPHcompare modernization and validation."
        ),
        "source_revision": source_revision,
        "evidence_file_count": count,
        "uncompressed_evidence_bytes": total_bytes,
        "license": "license_review_required_before_deposit",
        "creators": [
            {"name": "Jonah Daneshmand", "order": 1},
            {"name": "Julia H. Chariker", "order": 2},
            {"name": "Akshitkumar Mistry", "order": 3},
            {"name": "Eric C. Rouchka", "order": 4},
        ],
        "creator_order_status": "provisional_author_team_review_required",
        "credit_note": (
            "Provisional creator order follows the public preprint record and "
            "preserves project-owner direction that Dr. Eric C. Rouchka and "
            "Dr. Akshitkumar Mistry receive credit; final creator/contributor "
            "classification, roles, ORCIDs, affiliations, and CRediT statements "
            "require author-team approval."
        ),
        "related_identifiers": [
            {
                "identifier": "https://github.com/jcdaneshmand/scPHcompare",
                "relation": "isSupplementTo",
                "scheme": "url",
            },
            {
                "identifier": "software_version_doi_pending",
                "relation": "isSupplementTo",
                "scheme": "doi",
            },
        ],
        "doi": "pending_reservation",
    }
    return (json.dumps(metadata, indent=2, sort_keys=True) + "\n").encode("utf-8")


def build(args: argparse.Namespace) -> None:
    repo = Path(args.repo).resolve()
    base = resolve_revision(repo, args.base)
    revision = resolve_revision(repo, args.revision)
    if subprocess.run(
        ["git", "-C", str(repo), "merge-base", "--is-ancestor", base, revision]
    ).returncode != 0:
        raise RuntimeError("The publication base is not an ancestor of the source revision.")

    files = inventory(repo, base, revision)
    if not files:
        raise RuntimeError("No publication evidence matched the frozen classification.")

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    archive_path = output_dir / f"scphcompare-derived-evidence-{revision[:8]}.zip"
    output_manifest = output_dir / "MANIFEST.csv"
    metadata_path = output_dir / "ZENODO_METADATA_DRAFT.json"
    for path in (archive_path, output_manifest, metadata_path):
        if path.exists() and not args.force:
            raise FileExistsError(f"Refusing to overwrite {path}; use --force explicitly.")

    payloads: dict[str, bytes] = {}
    rows: list[dict[str, str]] = []
    with GitBlobReader(repo) as reader:
        for item in files:
            data = reader.read(item.git_blob, item.size)
            member = f"evidence/{item.path}"
            payloads[member] = data
            rows.append(
                {
                    "contract_id": CONTRACT_ID,
                    "archive_member": member,
                    "source_path": item.path,
                    "category": item.category,
                    "media_type": item.media_type,
                    "bytes": str(item.size),
                    "sha256": hashlib.sha256(data).hexdigest(),
                    "git_blob": item.git_blob,
                    "source_revision": revision,
                    "base_revision": base,
                }
            )

    manifest = manifest_bytes(rows)
    metadata = metadata_bytes(revision, len(files), sum(item.size for item in files))
    readme = archive_readme(revision, base, len(files))
    with zipfile.ZipFile(
        archive_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as archive:
        archive.writestr(zip_info("README.txt"), readme)
        archive.writestr(zip_info("MANIFEST.csv"), manifest)
        archive.writestr(zip_info("ZENODO_METADATA_DRAFT.json"), metadata)
        for member in sorted(payloads, key=lambda value: value.encode("utf-8")):
            archive.writestr(zip_info(member), payloads[member])

    output_manifest.write_bytes(manifest)
    metadata_path.write_bytes(metadata)
    archive_sha = hashlib.sha256(archive_path.read_bytes()).hexdigest()
    (output_dir / f"{archive_path.name}.sha256").write_text(
        f"{archive_sha}  {archive_path.name}\n", encoding="ascii", newline="\n"
    )
    if args.git_manifest:
        git_manifest = (repo / args.git_manifest).resolve()
        if repo not in git_manifest.parents:
            raise RuntimeError("The Git manifest must remain inside the repository.")
        git_manifest.parent.mkdir(parents=True, exist_ok=True)
        git_manifest.write_bytes(manifest)

    if args.lean_tree_dir:
        lean_dir = Path(args.lean_tree_dir).resolve()
        if lean_dir.exists() and any(lean_dir.iterdir()):
            raise RuntimeError(f"Lean tree directory is not empty: {lean_dir}")
        lean_dir.mkdir(parents=True, exist_ok=True)
        external = {item.path for item in files}
        raw_tree = run_git(repo, "ls-tree", "-r", "-l", "-z", revision, text=False)
        tree_files: list[EvidenceFile] = []
        for record in raw_tree.split(b"\0"):
            if not record:
                continue
            metadata_row, encoded_path = record.split(b"\t", 1)
            mode, kind, blob, encoded_size = metadata_row.decode("ascii").split()
            path = encoded_path.decode("utf-8")
            if path in external or kind != "blob" or mode == "160000":
                continue
            tree_files.append(EvidenceFile(path, "source", "", int(encoded_size), blob))
        with GitBlobReader(repo) as reader:
            for item in tree_files:
                destination = lean_dir.joinpath(*PurePosixPath(item.path).parts)
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(reader.read(item.git_blob, item.size))

    print(f"contract_id={CONTRACT_ID}")
    print(f"source_revision={revision}")
    print(f"base_revision={base}")
    print(f"evidence_files={len(files)}")
    print(f"uncompressed_evidence_bytes={sum(item.size for item in files)}")
    print(f"archive={archive_path}")
    print(f"archive_bytes={archive_path.stat().st_size}")
    print(f"archive_sha256={archive_sha}")


def verify(args: argparse.Namespace) -> None:
    archive_path = Path(args.archive).resolve()
    external_manifest = Path(args.manifest).read_bytes()
    with zipfile.ZipFile(archive_path, "r") as archive:
        if archive.read("MANIFEST.csv") != external_manifest:
            raise RuntimeError("External and archived manifests differ.")
        rows = list(csv.DictReader(io.TextIOWrapper(
            io.BytesIO(external_manifest), encoding="utf-8", newline=""
        )))
        expected = {"README.txt", "MANIFEST.csv", "ZENODO_METADATA_DRAFT.json"}
        for row in rows:
            member = row["archive_member"]
            data = archive.read(member)
            if len(data) != int(row["bytes"]):
                raise RuntimeError(f"Size mismatch: {member}")
            if hashlib.sha256(data).hexdigest() != row["sha256"]:
                raise RuntimeError(f"SHA-256 mismatch: {member}")
            expected.add(member)
        if set(archive.namelist()) != expected:
            raise RuntimeError("Archive member set differs from the manifest contract.")
    print(f"verified_files={len(rows)}")
    print(f"archive_sha256={hashlib.sha256(archive_path.read_bytes()).hexdigest()}")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)
    build_parser = subparsers.add_parser("build")
    build_parser.add_argument("--repo", default=Path(__file__).resolve().parents[1])
    build_parser.add_argument("--base", required=True)
    build_parser.add_argument("--revision", required=True)
    build_parser.add_argument("--output-dir", required=True)
    build_parser.add_argument("--git-manifest")
    build_parser.add_argument("--lean-tree-dir")
    build_parser.add_argument("--force", action="store_true")
    build_parser.set_defaults(function=build)
    verify_parser = subparsers.add_parser("verify")
    verify_parser.add_argument("--archive", required=True)
    verify_parser.add_argument("--manifest", required=True)
    verify_parser.set_defaults(function=verify)
    return result


def main() -> None:
    args = parser().parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
