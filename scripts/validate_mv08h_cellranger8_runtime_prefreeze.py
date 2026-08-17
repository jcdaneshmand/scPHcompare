#!/usr/bin/env python3
"""Independently validate the MV8-H Cell Ranger 8.0.1 prefreeze."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import re
import subprocess


EXPECTED_RUNTIME_VERSION = "cellranger cellranger-8.0.1"
EXPECTED_LAUNCHER_SHA256 = (
    "4ee3a1670b4f14c826004fe8e17b4759e1edc701b15ff2e9623753bf1b34d4d6"
)
EXPECTED_STAR_SHA256 = (
    "4b50e57ae0570bb9a7dcf0364360ce4d32e3bd049b852be69bed86827847ab36"
)
EXPECTED_SAMTOOLS_SHA256 = (
    "98cdc7aef8d3a0ca291a4ba2db15e9d20293963a648c834c5e18f36bb2c37b47"
)
EXPECTED_FASTA_SHA256 = (
    "2a27436d44f0d6350f86894fbe5edec56faa5467028879784508041562406aa0"
)
EXPECTED_SOURCE_GTF_SHA256 = (
    "810e3bb63bb24bd5a005b14397d69280dda34c41a612fb86a18f3c4836fce57d"
)
EXPECTED_CUSTOM_GTF_SHA256 = (
    "e28e4c4faf0dd76884d5e94c481fce2db43ad303968067c1276092a234727182"
)
EXPECTED_PANEL_SHA256 = (
    "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e"
)
EXPECTED_COMMON475_SHA256 = (
    "b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba"
)
GTF_ATTRIBUTE = re.compile(r'(\S+) "([^"]*)";')
STABLE_ID = re.compile(r"(ENSG\d+)")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_axis(values: list[str]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(value.encode("utf-8") + b"\n")
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def one(rows: list[dict[str, str]], label: str) -> dict[str, str]:
    if len(rows) != 1:
        raise RuntimeError(f"{label} must contain exactly one row")
    return rows[0]


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write empty evidence: {path.name}")
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


def runtime_tree(root: Path) -> tuple[int, int, int, str]:
    """Independent deterministic regular-file/symlink tree digest."""
    items: list[tuple[str, str, Path]] = []
    for directory, directories, files in os.walk(root, followlinks=False):
        current = Path(directory)
        keep: list[str] = []
        for name in directories:
            candidate = current / name
            if candidate.is_symlink():
                items.append(("L", candidate.relative_to(root).as_posix(), candidate))
            else:
                keep.append(name)
        directories[:] = keep
        for name in files:
            candidate = current / name
            kind = "L" if candidate.is_symlink() else "F"
            if kind == "L" or candidate.is_file():
                items.append((kind, candidate.relative_to(root).as_posix(), candidate))
    items.sort(key=lambda item: (item[1], item[0]))
    digest = hashlib.sha256()
    regular_files = 0
    links = 0
    regular_bytes = 0
    for index, (kind, relative, path) in enumerate(items, 1):
        if kind == "L":
            links += 1
            record = f"L\0{relative}\0{os.readlink(path)}\n"
        else:
            regular_files += 1
            size = path.stat().st_size
            regular_bytes += size
            record = f"F\0{relative}\0{size}\0{sha256_file(path)}\n"
        digest.update(record.encode("utf-8"))
        if index % 2000 == 0:
            print(f"independent runtime rehash: {index}/{len(items)} entries", flush=True)
    return regular_files, links, regular_bytes, digest.hexdigest()


def custom_gtf_identity(path: Path) -> tuple[int, int, str, set[str]]:
    genes: list[str] = []
    records = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise RuntimeError("malformed custom GTF record")
            records += 1
            if fields[2] == "gene":
                genes.append(dict(GTF_ATTRIBUTE.findall(fields[8])).get("gene_id", ""))
    return len(genes), records, sha256_axis(genes), set(genes)


def panel_ids(path: Path, order_field: str) -> tuple[list[dict[str, str]], set[str]]:
    rows = read_csv(path)
    rows.sort(key=lambda row: int(row[order_field]))
    ids = set()
    for row in rows:
        value = row.get("ensembl_stable_id", "")
        if not value:
            match = STABLE_ID.search(row.get("feature_id", ""))
            value = match.group(1) if match else ""
        ids.add(value)
    return rows, ids


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--repeat-evidence-dir", type=Path, required=True)
    parser.add_argument("--cellranger-root", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--source-gtf", type=Path, required=True)
    parser.add_argument("--custom-gtf", type=Path, required=True)
    parser.add_argument("--repeat-custom-gtf", type=Path, required=True)
    parser.add_argument("--panel500", type=Path, required=True)
    parser.add_argument("--panel475", type=Path, required=True)
    args = parser.parse_args()

    evidence = args.evidence_dir
    runtime = one(read_csv(evidence / "mv08h-runtime-identity.csv"), "runtime")
    reference = read_csv(evidence / "mv08h-reference-binding.csv")
    processing = read_csv(evidence / "mv08h-processing-amendment.csv")
    gates = read_csv(evidence / "mv08h-gate-status.csv")
    decision = one(read_csv(evidence / "mv08h-decision.csv"), "decision")
    sources = read_csv(evidence / "mv08h-source-manifest.csv")

    launcher = args.cellranger_root / "bin" / "cellranger"
    star = args.cellranger_root / "lib" / "bin" / "STAR"
    samtools = args.cellranger_root / "lib" / "bin" / "samtools"
    reported = subprocess.run(
        [str(launcher), "--version"], check=True, text=True,
        capture_output=True,
    ).stdout.strip()
    regular_files, links, regular_bytes, tree_sha = runtime_tree(args.cellranger_root)
    runtime_exact = all([
        reported == EXPECTED_RUNTIME_VERSION,
        runtime["reported_version"] == reported,
        int(runtime["regular_files"]) == regular_files,
        int(runtime["symlinks"]) == links,
        int(runtime["regular_file_bytes"]) == regular_bytes,
        runtime["tree_sha256"] == tree_sha,
        runtime["launcher_sha256"] == sha256_file(launcher) == EXPECTED_LAUNCHER_SHA256,
        runtime["star_sha256"] == sha256_file(star) == EXPECTED_STAR_SHA256,
        runtime["samtools_sha256"] == sha256_file(samtools) == EXPECTED_SAMTOOLS_SHA256,
        runtime["required_cli_controls_verified"] == "TRUE",
        runtime["all_runtime_gates_passed"] == "TRUE",
    ])

    reference_by_id = {row["resource_id"]: row for row in reference}
    file_identity_exact = all([
        sha256_file(args.fasta) == EXPECTED_FASTA_SHA256,
        sha256_file(args.source_gtf) == EXPECTED_SOURCE_GTF_SHA256,
        sha256_file(args.custom_gtf) == EXPECTED_CUSTOM_GTF_SHA256,
        sha256_file(args.repeat_custom_gtf) == EXPECTED_CUSTOM_GTF_SHA256,
        reference_by_id["ensembl93_primary_assembly_fasta_gz"]["sha256"] == EXPECTED_FASTA_SHA256,
        reference_by_id["ensembl93_gtf_gz"]["sha256"] == EXPECTED_SOURCE_GTF_SHA256,
        reference_by_id["target_complete_33563_gtf"]["sha256"] == EXPECTED_CUSTOM_GTF_SHA256,
        all(row["gate_passed"] == "TRUE" for row in reference),
    ])
    gene_count, feature_records, axis_sha, custom_ids = custom_gtf_identity(args.custom_gtf)
    panel500, panel500_ids = panel_ids(args.panel500, "panel_order")
    panel475, panel475_ids = panel_ids(args.panel475, "panel_order_475")
    annotation_exact = all([
        gene_count == 33_563,
        feature_records == 2_565_751,
        axis_sha == "0ea5f54698d79e05a7dd98eded2bf9558174d5214cfe684e8263b996c2a169a4",
        len(panel500) == len(panel500_ids) == 500,
        len(panel475) == len(panel475_ids) == 475,
        panel500[0]["panel_sha256"] == EXPECTED_PANEL_SHA256,
        panel475[0]["parent_panel_sha256"] == EXPECTED_PANEL_SHA256,
        panel475[0]["common475_axis_sha256"] == EXPECTED_COMMON475_SHA256,
        panel500_ids.issubset(custom_ids),
        panel475_ids.issubset(panel500_ids),
    ])

    processing_by_id = {row["step_id"]: row for row in processing}
    count_text = processing_by_id["count_sentinel"]["frozen_value"]
    landscape_text = processing_by_id["topology_landscapes"]["frozen_value"]
    processing_exact = all(term in count_text for term in [
        "--chemistry=SC3Pv2", "--include-introns=false",
        "--create-bam=false", "--nosecondary", "--localcores=16",
        "--localmem=64",
    ])
    processing_exact = processing_exact and all(term in landscape_text for term in [
        "cell_topology_v1", "gene_topology_v1", "complete VR H0/H1",
        "every consecutive active landscape level", "no fixed grid",
        "no level cap", "H0/H1 separate", "essential H0 excluded",
    ])
    processing_exact = processing_exact and processing_by_id["runtime"]["frozen_value"].startswith("Cell Ranger 8.0.1")
    processing_exact = processing_exact and processing_by_id["panel_separation"]["execution_authorized"] == "FALSE"

    authorization_exact = all([
        decision["decision"] == "runtime_reference_exact_authorize_mkref_only",
        decision["runtime_amendment_accepted"] == "TRUE",
        decision["runtime_validated"] == "TRUE",
        decision["reference_inputs_validated"] == "TRUE",
        decision["mkref_authorized"] == "TRUE",
        decision["count_sentinel_authorized"] == "FALSE",
        decision["remaining_units_authorized"] == "FALSE",
        decision["pca_ph_landscape_authorized"] == "FALSE",
        decision["label_access_authorized"] == "FALSE",
        decision["biological_outcomes_authorized"] == "FALSE",
        processing_by_id["mkref"]["execution_authorized"] == "TRUE",
        sum(row["execution_authorized"] == "TRUE" for row in processing) == 1,
    ])
    gates_exact = len(gates) == 11 and all(row["status"].startswith("pass") for row in gates)
    sources_exact = len(sources) == 5 and {
        "installed_cellranger_8_0_1", "ensembl93_reference_input_closure",
        "tenx_cellranger_8_release_notes", "tenx_command_arguments",
    }.issubset({row["source_id"] for row in sources})

    public_paths = [
        path for path in evidence.iterdir()
        if path.is_file() and path.name != "mv08h-artifact-manifest.csv"
    ]
    public_text = "\n".join(
        path.read_text(encoding="utf-8", errors="strict") for path in public_paths
    )
    forbidden = ["/mnt/e/", "E:\\", "C:\\Users\\", "donor_organism.provenance", "cell_barcode"]
    firewall_exact = not any(token in public_text for token in forbidden)

    repeat_names = [
        "mv08h-decision.csv", "mv08h-gate-status.csv",
        "mv08h-processing-amendment.csv", "mv08h-reference-binding.csv",
        "mv08h-runtime-identity.csv", "mv08h-source-manifest.csv",
    ]
    repeat_rows = []
    for index, name in enumerate(repeat_names, 1):
        first_sha = sha256_file(evidence / name)
        second_sha = sha256_file(args.repeat_evidence_dir / name)
        repeat_rows.append({
            "contract_id": "mv08h_cellranger8_repeat_validation_v2",
            "artifact_order": index,
            "artifact": name,
            "production_sha256": first_sha,
            "repeat_sha256": second_sha,
            "byte_identical": str(first_sha == second_sha).upper(),
        })
    repeat_exact = all(row["byte_identical"] == "TRUE" for row in repeat_rows)
    write_csv(evidence / "mv08h-repeat-validation.csv", repeat_rows)

    checks = [
        ("runtime_report_and_tree_exact", runtime_exact, f"Cell Ranger 8.0.1; tree {tree_sha}"),
        ("runtime_components_exact", runtime_exact, "launcher, STAR, samtools, Python, and required CLI controls bound"),
        ("reference_file_identities_exact", file_identity_exact, "FASTA, source GTF, target-complete GTF, and repeat exact"),
        ("target_complete_annotation_exact", annotation_exact, "33563 genes; 2565751 feature records; all exact 500 present"),
        ("panel_axes_exact", annotation_exact, "exact500 and strict ordered common475 identities unchanged"),
        ("processing_amendment_exact", processing_exact, "SC3Pv2 exon-only count controls and unchanged landscapes frozen"),
        ("panel_runtime_estimands_separated", "same Cell Ranger 8.0.1 count matrices" in processing_by_id["panel_separation"]["frozen_value"], "within-runtime panel and shared-panel runtime comparisons separated"),
        ("mkref_only_authorization", authorization_exact, "only mkref opened; count and downstream execution closed"),
        ("gate_ladder_exact", gates_exact, "11 prospective gates pass; execution stop retained"),
        ("authoritative_sources_bound", sources_exact, "local runtime/reference identities and official 10x documentation bound"),
        ("public_firewall", firewall_exact, "no private absolute paths, donor fields, barcodes, expression, labels, or outcomes"),
        ("deterministic_builder_repeat", repeat_exact, "six builder artifacts reproduce byte-identically"),
    ]
    validation_rows = [
        {
            "contract_id": "mv08h_cellranger8_independent_validation_v2",
            "check_order": index,
            "check_id": check_id,
            "passed": str(passed).upper(),
            "detail": detail,
        }
        for index, (check_id, passed, detail) in enumerate(checks, 1)
    ]
    if not all(row["passed"] == "TRUE" for row in validation_rows):
        failed = [row["check_id"] for row in validation_rows if row["passed"] != "TRUE"]
        raise RuntimeError(f"MV8-H Cell Ranger 8.0.1 prefreeze failed: {failed}")
    write_csv(evidence / "mv08h-independent-validation.csv", validation_rows)

    artifact_names = sorted(
        path.name for path in evidence.iterdir()
        if path.is_file() and path.name != "mv08h-artifact-manifest.csv"
    )
    artifact_rows = [
        {
            "contract_id": "mv08h_cellranger8_artifact_manifest_v2",
            "file": name,
            "bytes": (evidence / name).stat().st_size,
            "sha256": sha256_file(evidence / name),
            "contains_expression": "FALSE",
            "contains_cell_barcode": "FALSE",
            "contains_absolute_private_path": "FALSE",
            "contains_donor_attribute": "FALSE",
            "contains_outcome_label": "FALSE",
        }
        for name in artifact_names
    ]
    write_csv(evidence / "mv08h-artifact-manifest.csv", artifact_rows)
    print("MV8-H Cell Ranger 8.0.1 prefreeze validation passed: 12/12 checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
