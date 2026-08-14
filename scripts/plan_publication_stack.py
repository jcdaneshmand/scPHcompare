#!/usr/bin/env python3
"""Plan dependency-ordered Git publication slices after evidence externalization."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


EXCLUDES = (
    ":(exclude,glob)docs/audits/*.csv",
    ":(exclude,glob)docs/audits/*.csv.gz",
    ":(exclude,glob)docs/audits/**/*.csv",
    ":(exclude,glob)docs/audits/**/*.csv.gz",
    ":(exclude,glob)results/**",
)


def git(repo: Path, *args: str, text: bool = True) -> str | bytes:
    return subprocess.run(
        ["git", "-C", str(repo), *args],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=text,
    ).stdout


def resolve(repo: Path, revision: str) -> str:
    return str(git(repo, "rev-parse", f"{revision}^{{commit}}")).strip()


def metrics(repo: Path, base: str, head: str) -> tuple[int, int, int, int]:
    pathspec = ("--", ".", *EXCLUDES)
    names = str(git(repo, "diff", "--name-only", base, head, *pathspec)).splitlines()
    insertions = deletions = 0
    for line in str(git(repo, "diff", "--numstat", base, head, *pathspec)).splitlines():
        added, removed, _ = line.split("\t", 2)
        if added.isdigit():
            insertions += int(added)
        if removed.isdigit():
            deletions += int(removed)
    raw_bytes = len(git(repo, "diff", "--no-ext-diff", "--binary", base, head, *pathspec, text=False))
    return len(names), insertions, deletions, raw_bytes


def subject(repo: Path, revision: str) -> str:
    return str(git(repo, "show", "-s", "--format=%s", revision)).strip()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", default=Path(__file__).resolve().parents[1])
    parser.add_argument("--base", required=True)
    parser.add_argument("--revision", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-files", type=int, default=250)
    parser.add_argument("--max-lines", type=int, default=15000)
    parser.add_argument("--max-raw-bytes", type=int, default=750000)
    args = parser.parse_args()
    repo = Path(args.repo).resolve()
    base = resolve(repo, args.base)
    revision = resolve(repo, args.revision)
    commits = str(git(repo, "rev-list", "--reverse", f"{base}..{revision}")).splitlines()
    if not commits:
        raise RuntimeError("No commits are available for publication planning.")

    slices: list[dict[str, str | int | bool]] = []
    cursor = base
    index = 0
    while index < len(commits):
        start_index = index
        accepted = index
        accepted_metrics = metrics(repo, cursor, commits[index])
        index += 1
        while index < len(commits):
            candidate = metrics(repo, cursor, commits[index])
            files, added, removed, raw_bytes = candidate
            if (
                files > args.max_files
                or added + removed > args.max_lines
                or raw_bytes > args.max_raw_bytes
            ):
                break
            accepted = index
            accepted_metrics = candidate
            index += 1
        head = commits[accepted]
        files, added, removed, raw_bytes = accepted_metrics
        slices.append(
            {
                "contract_id": "scphcompare_publication_stack_v1",
                "slice_id": f"P{len(slices) + 1:02d}",
                "base_revision": cursor,
                "head_revision": head,
                "first_source_commit": commits[start_index],
                "last_source_commit": head,
                "source_commit_count": accepted - start_index + 1,
                "files_changed": files,
                "insertions": added,
                "deletions": removed,
                "raw_diff_bytes": raw_bytes,
                "first_subject": subject(repo, commits[start_index]),
                "last_subject": subject(repo, head),
                "within_file_target": files <= args.max_files,
                "within_line_target": added + removed <= args.max_lines,
                "within_raw_diff_target": raw_bytes <= args.max_raw_bytes,
                "publication_order": len(slices) + 1,
            }
        )
        cursor = head
        index = accepted + 1

    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(slices[0]))
        writer.writeheader()
        writer.writerows(slices)
    print(f"slices={len(slices)}")
    print(f"max_files={max(int(row['files_changed']) for row in slices)}")
    print(f"max_lines={max(int(row['insertions']) + int(row['deletions']) for row in slices)}")
    print(f"max_raw_diff_bytes={max(int(row['raw_diff_bytes']) for row in slices)}")


if __name__ == "__main__":
    main()
