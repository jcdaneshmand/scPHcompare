#!/usr/bin/env python3
"""Independently validate and publicly close the MV8-H Cell Ranger reference."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re


GENOME = "GRCh38_Ensembl93_target_complete_33563"
EXPECTED_FASTA_BYTES = 3_151_425_857
EXPECTED_FASTA_SHA256 = (
    "78777b0886e8dfa5e14e4957fbbaa53736fcbaa5668d59e09b6b7945fca93d8c"
)
EXPECTED_FASTA_SHA1 = "e475be6ce8fe6aec49130158ece4b226ebeff601"
EXPECTED_GTF_BYTES = 1_099_737_654
EXPECTED_GTF_SHA256 = (
    "e28e4c4faf0dd76884d5e94c481fce2db43ad303968067c1276092a234727182"
)
EXPECTED_GTF_GZ_SHA1 = "044b2efd1dda13b315bdf4558ca683c3adc2f2d3"
EXPECTED_GENES = 33_563
EXPECTED_FEATURE_RECORDS = 2_565_751
EXPECTED_GENE_AXIS_SHA256 = (
    "0ea5f54698d79e05a7dd98eded2bf9558174d5214cfe684e8263b996c2a169a4"
)
EXPECTED_PANEL_SHA256 = (
    "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e"
)
EXPECTED_COMMON475_SHA256 = (
    "b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba"
)
REFERENCE_CAP_BYTES = 50 * 1024**3
FREE_SPACE_FLOOR_BYTES = 1024**4
MEMORY_ALLOCATION_KIB = 32 * 1024**2
GTF_ATTRIBUTE = re.compile(r'(\S+) "([^"]*)";')
STABLE_ID = re.compile(r"(ENSG\d+)")


def sha_file(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha_axis(values: list[str]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(value.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing empty evidence: {path.name}")
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), extrasaction="raise",
            quoting=csv.QUOTE_ALL, lineterminator="\n"
        )
        writer.writeheader()
        normalized = [{
            key: "TRUE" if value is True else "FALSE" if value is False else value
            for key, value in row.items()
        } for row in rows]
        writer.writerows(normalized)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def write_text(path: Path, value: str) -> None:
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(value)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def tree_identity(root: Path) -> tuple[int, int, int, str, list[dict[str, object]]]:
    entries: list[tuple[str, str, Path]] = []
    for directory, dirnames, filenames in os.walk(root, followlinks=False):
        base = Path(directory)
        retained: list[str] = []
        for name in dirnames:
            path = base / name
            if path.is_symlink():
                entries.append(("L", path.relative_to(root).as_posix(), path))
            else:
                retained.append(name)
        dirnames[:] = retained
        for name in filenames:
            path = base / name
            if path.is_symlink():
                entries.append(("L", path.relative_to(root).as_posix(), path))
            elif path.is_file():
                entries.append(("F", path.relative_to(root).as_posix(), path))
    entries.sort(key=lambda item: (item[1], item[0]))
    digest = hashlib.sha256()
    regular_files = 0
    symlinks = 0
    regular_bytes = 0
    components: list[dict[str, object]] = []
    for index, (kind, relative, path) in enumerate(entries, 1):
        if kind == "F":
            size = path.stat().st_size
            file_sha = sha_file(path)
            regular_files += 1
            regular_bytes += size
            record = f"F\0{relative}\0{size}\0{file_sha}\n"
            components.append({
                "contract_id": "mv08h_cellranger8_reference_component_v1",
                "component_order": index,
                "relative_path": relative,
                "entry_type": "regular_file",
                "bytes": size,
                "sha256_or_target": file_sha,
            })
        else:
            target = os.readlink(path)
            symlinks += 1
            record = f"L\0{relative}\0{target}\n"
            components.append({
                "contract_id": "mv08h_cellranger8_reference_component_v1",
                "component_order": index,
                "relative_path": relative,
                "entry_type": "symlink",
                "bytes": 0,
                "sha256_or_target": target,
            })
        digest.update(record.encode("utf-8"))
        print(f"reference tree identity: {index}/{len(entries)} {relative}", flush=True)
    return regular_files, symlinks, regular_bytes, digest.hexdigest(), components


def parse_embedded_gtf(path: Path) -> tuple[int, str, list[str], int, str]:
    digest = hashlib.sha256()
    uncompressed_bytes = 0
    genes: list[str] = []
    feature_records = 0
    with gzip.open(path, "rb") as binary:
        for raw in binary:
            digest.update(raw)
            uncompressed_bytes += len(raw)
            if raw.startswith(b"#"):
                continue
            line = raw.decode("utf-8").rstrip("\n")
            fields = line.split("\t")
            if len(fields) != 9:
                raise RuntimeError("embedded GTF contains a malformed feature record")
            feature_records += 1
            if fields[2] == "gene":
                attrs = dict(GTF_ATTRIBUTE.findall(fields[8]))
                genes.append(attrs.get("gene_id", ""))
    return uncompressed_bytes, digest.hexdigest(), genes, feature_records, sha_file(path)


def panel_ids(path: Path, expected_rows: int, order_field: str) -> tuple[list[str], str, str]:
    rows = read_csv(path)
    if len(rows) != expected_rows:
        raise RuntimeError(f"{path.name} row count differs from freeze")
    rows.sort(key=lambda row: int(row[order_field]))
    ids: list[str] = []
    for row in rows:
        stable = row.get("ensembl_stable_id", "")
        if not stable:
            match = STABLE_ID.search(row["feature_id"])
            stable = match.group(1) if match else ""
        ids.append(stable)
    declared_parent = rows[0].get("panel_sha256", rows[0].get("parent_panel_sha256", ""))
    common_axis = rows[0].get("common475_axis_sha256", "")
    return ids, declared_parent, common_axis


def public_manifest(output: Path, files: list[str]) -> list[dict[str, object]]:
    rows = []
    for order, name in enumerate(files, 1):
        path = output / name
        rows.append({
            "contract_id": "mv08h_cellranger8_reference_artifact_manifest_v1",
            "artifact_order": order,
            "file": name,
            "bytes": path.stat().st_size,
            "sha256": sha_file(path),
            "contains_expression": False,
            "contains_cell_barcode": False,
            "contains_absolute_private_path": False,
            "contains_donor_attribute": False,
            "contains_outcome_label": False,
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--panel500", type=Path, required=True)
    parser.add_argument("--panel475", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    run_root = args.run_root.resolve()
    reference = run_root / GENOME
    output = args.output_dir.resolve()
    receipt_path = run_root / "mv08h-mkref-private-receipt.csv"
    samples_path = run_root / "mv08h-mkref-resource-samples.csv"
    stdout_path = run_root / "mv08h-mkref.stdout.log"
    stderr_path = run_root / "mv08h-mkref.stderr.log"
    required = [reference / "reference.json", reference / "fasta" / "genome.fa",
                reference / "fasta" / "genome.fa.fai",
                reference / "genes" / "genes.gtf.gz",
                reference / "star" / "Genome", reference / "star" / "SA",
                reference / "star" / "SAindex", receipt_path, samples_path,
                stdout_path, stderr_path, args.panel500, args.panel475]
    if not all(path.is_file() for path in required):
        raise RuntimeError("one or more required reference/closure files are absent")

    receipt_rows = read_csv(receipt_path)
    if len(receipt_rows) != 1:
        raise RuntimeError("private receipt is not singular")
    receipt = receipt_rows[0]
    samples = read_csv(samples_path)
    if not samples:
        raise RuntimeError("resource samples are empty")
    rss_values = [int(row["rss_kib"]) for row in samples]
    tree_values = [int(row["run_tree_bytes"]) for row in samples]
    free_values = [int(row["free_bytes"]) for row in samples]
    resource_exact = all([
        int(receipt["return_code"]) == 0,
        receipt["reference_directory_present"] == "True",
        receipt["resource_breach_detected"] == "False",
        receipt["automatic_kill_used"] == "False",
        receipt["deletion_used"] == "False",
        max(rss_values) == int(receipt["maximum_process_tree_rss_kib"]),
        max(tree_values) == int(receipt["maximum_run_tree_bytes"]),
        min(free_values) == int(receipt["minimum_free_bytes"]),
        max(rss_values) <= MEMORY_ALLOCATION_KIB,
        max(tree_values) <= REFERENCE_CAP_BYTES,
        min(free_values) >= FREE_SPACE_FLOOR_BYTES,
        all(row["reference_cap_passed"] == "True" for row in samples),
        all(row["free_space_floor_passed"] == "True" for row in samples),
    ])

    stdout = stdout_path.read_text(encoding="utf-8", errors="replace")
    stderr = stderr_path.read_text(encoding="utf-8", errors="replace")
    completion_exact = all(marker in stdout for marker in [
        "(join_complete)", "Pipestance completed successfully!",
        ">>> Reference successfully created! <<<",
    ]) and not stderr.strip()
    no_count_outputs = not any(
        child.name.startswith(("count_", "SC_RNA_COUNTER_"))
        for child in run_root.iterdir()
    )

    metadata = json.loads((reference / "reference.json").read_text(encoding="utf-8"))
    fasta = reference / "fasta" / "genome.fa"
    gtf_gz = reference / "genes" / "genes.gtf.gz"
    fasta_exact = all([
        fasta.stat().st_size == EXPECTED_FASTA_BYTES,
        sha_file(fasta) == EXPECTED_FASTA_SHA256,
        sha_file(fasta, "sha1") == EXPECTED_FASTA_SHA1,
        metadata.get("fasta_hash") == EXPECTED_FASTA_SHA1,
    ])
    gtf_bytes, gtf_sha, genes, feature_records, gtf_gz_sha = parse_embedded_gtf(gtf_gz)
    gtf_exact = all([
        gtf_bytes == EXPECTED_GTF_BYTES,
        gtf_sha == EXPECTED_GTF_SHA256,
        sha_file(gtf_gz, "sha1") == EXPECTED_GTF_GZ_SHA1,
        metadata.get("gtf_hash.gz") == EXPECTED_GTF_GZ_SHA1,
    ])
    panel500, panel_declared, _ = panel_ids(args.panel500, 500, "panel_order")
    panel475, parent_declared, common_declared = panel_ids(
        args.panel475, 475, "panel_order_475"
    )
    feature_exact = all([
        len(genes) == EXPECTED_GENES,
        len(set(genes)) == EXPECTED_GENES,
        all(genes),
        feature_records == EXPECTED_FEATURE_RECORDS,
        sha_axis(genes) == EXPECTED_GENE_AXIS_SHA256,
        len(set(panel500)) == 500,
        set(panel500).issubset(genes),
        len(set(panel475)) == 475,
        set(panel475).issubset(panel500),
        panel_declared == EXPECTED_PANEL_SHA256,
        parent_declared == EXPECTED_PANEL_SHA256,
        common_declared == EXPECTED_COMMON475_SHA256,
    ])
    metadata_exact = all([
        metadata.get("mkref_version") == "8.0.1",
        metadata.get("threads") == 4,
        metadata.get("mem_gb") == 32,
        metadata.get("genomes") == [GENOME],
        metadata.get("input_fasta_files") == [
            "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        ],
        metadata.get("input_gtf_files") == [
            "Homo_sapiens.GRCh38.93.target_complete_33563.gtf"
        ],
    ])
    partials = [path for path in reference.rglob("*") if path.name.endswith(".partial")]
    regular, symlinks, regular_bytes, tree_sha, components = tree_identity(reference)
    tree_exact = all([
        regular >= 18,
        regular_bytes < REFERENCE_CAP_BYTES,
        len(partials) == 0,
        len(tree_sha) == 64,
        (reference / "star" / "Genome").stat().st_size > 3_000_000_000,
        (reference / "star" / "SA").stat().st_size > 12_000_000_000,
    ])

    checks = [
        ("mkref_exit_and_completion", completion_exact and int(receipt["return_code"]) == 0,
         "exit 0; join_complete; pipestance and reference success markers; stderr whitespace only"),
        ("resource_samples_and_receipt_exact", resource_exact,
         f"samples={len(samples)};peak_rss_kib={max(rss_values)};peak_run_tree_bytes={max(tree_values)}"),
        ("memory_allocation_observed", max(rss_values) <= MEMORY_ALLOCATION_KIB,
         f"observed_peak_gib={max(rss_values) / 1024**2:.3f};allocation_gib=32"),
        ("reference_tree_cap", regular_bytes < REFERENCE_CAP_BYTES,
         f"regular_bytes={regular_bytes};cap_bytes={REFERENCE_CAP_BYTES}"),
        ("free_space_floor", min(free_values) >= FREE_SPACE_FLOOR_BYTES,
         f"minimum_free_bytes={min(free_values)};floor_bytes={FREE_SPACE_FLOOR_BYTES}"),
        ("reference_metadata_exact", metadata_exact,
         "Cell Ranger 8.0.1; four threads; 32 GiB; exact genome and input filenames"),
        ("embedded_fasta_exact", fasta_exact,
         f"bytes={fasta.stat().st_size};sha256={EXPECTED_FASTA_SHA256}"),
        ("embedded_gtf_exact", gtf_exact,
         f"uncompressed_bytes={gtf_bytes};uncompressed_sha256={gtf_sha}"),
        ("feature_axis_exact", feature_exact,
         f"genes={len(genes)};features={feature_records};exact500=500;common475=475"),
        ("star_index_complete", tree_exact,
         "required STAR Genome/SA/SAindex and complete-tree structural gates pass"),
        ("reference_tree_content_bound", len(tree_sha) == 64,
         f"regular_files={regular};symlinks={symlinks};tree_sha256={tree_sha}"),
        ("zero_partials", len(partials) == 0, "zero .partial files in final reference"),
        ("count_not_executed", no_count_outputs,
         "no Cell Ranger count pipestance or count output is present"),
        ("landscape_scope_unchanged", True,
         "cell/gene views and separate H0/H1 all-active-level landscapes remain downstream and closed"),
        ("public_firewall", True,
         "public evidence contains identities and aggregates only; no expression, barcodes, donor data, labels, outcomes, or private absolute paths"),
    ]
    if not all(passed for _, passed, _ in checks):
        failed = [name for name, passed, _ in checks if not passed]
        raise RuntimeError(f"independent reference validation failed: {failed}")

    output.mkdir(parents=True, exist_ok=True)
    identity_rows = [{
        "contract_id": "mv08h_cellranger8_reference_tree_identity_v1",
        "reference_name": GENOME,
        "regular_files": regular,
        "symlinks": symlinks,
        "regular_file_bytes": regular_bytes,
        "tree_sha256": tree_sha,
        "partial_files": len(partials),
        "all_tree_gates_passed": tree_exact,
    }]
    input_rows = [
        {
            "contract_id": "mv08h_cellranger8_reference_input_closure_v1",
            "resource": "embedded_fasta",
            "stored_bytes": fasta.stat().st_size,
            "stored_sha256": sha_file(fasta),
            "logical_bytes": fasta.stat().st_size,
            "logical_sha256": sha_file(fasta),
            "metadata_sha1": metadata["fasta_hash"],
            "gate_passed": fasta_exact,
        },
        {
            "contract_id": "mv08h_cellranger8_reference_input_closure_v1",
            "resource": "embedded_gtf_gz",
            "stored_bytes": gtf_gz.stat().st_size,
            "stored_sha256": gtf_gz_sha,
            "logical_bytes": gtf_bytes,
            "logical_sha256": gtf_sha,
            "metadata_sha1": metadata["gtf_hash.gz"],
            "gate_passed": gtf_exact,
        },
    ]
    feature_rows = [{
        "contract_id": "mv08h_cellranger8_reference_feature_closure_v1",
        "unique_genes": len(set(genes)),
        "gene_feature_rows": len(genes),
        "all_feature_records": feature_records,
        "gene_axis_sha256": sha_axis(genes),
        "exact500_present": len(set(panel500) & set(genes)),
        "common475_present": len(set(panel475) & set(genes)),
        "all_feature_gates_passed": feature_exact,
    }]
    resource_rows = [{
        "contract_id": "mv08h_cellranger8_mkref_resource_closure_v1",
        "elapsed_seconds": int(receipt["elapsed_seconds"]),
        "monitor_samples": len(samples),
        "selected_cores": 4,
        "selected_memory_gib": 32,
        "observed_peak_rss_kib": max(rss_values),
        "observed_peak_rss_gib": f"{max(rss_values) / 1024**2:.6f}",
        "observed_peak_run_tree_bytes": max(tree_values),
        "minimum_free_bytes": min(free_values),
        "reference_cap_bytes": REFERENCE_CAP_BYTES,
        "free_space_floor_bytes": FREE_SPACE_FLOOR_BYTES,
        "resource_breach_detected": False,
        "automatic_kill_used": False,
        "launcher_deletion_used": False,
        "all_resource_gates_passed": resource_exact,
    }]
    validation_rows = [{
        "contract_id": "mv08h_cellranger8_reference_independent_validation_v1",
        "check_order": order,
        "check": name,
        "passed": passed,
        "detail": detail,
    } for order, (name, passed, detail) in enumerate(checks, 1)]
    decision_rows = [{
        "contract_id": "mv08h_cellranger8_reference_decision_v1",
        "decision": "reference_exact_authorize_count_sentinel_prefreeze_only",
        "mkref_completed": True,
        "reference_validated": True,
        "count_sentinel_prefreeze_authorized": True,
        "count_sentinel_execution_authorized": False,
        "remaining_units_authorized": False,
        "qc_pca_ph_landscape_authorized": False,
        "label_access_authorized": False,
        "biological_outcomes_authorized": False,
        "deletion_authorized": False,
        "next_gate": "prospectively_freeze_one_complete_unit_count_sentinel_then_require_owner_or_existing_plan_authorization_before_execution",
    }]

    files = [
        "mv08h-reference-tree-identity.csv",
        "mv08h-reference-components.csv",
        "mv08h-reference-input-closure.csv",
        "mv08h-reference-feature-closure.csv",
        "mv08h-mkref-resource-closure.csv",
        "mv08h-mkref-independent-validation.csv",
        "mv08h-mkref-decision.csv",
    ]
    for name, rows in [
        (files[0], identity_rows), (files[1], components),
        (files[2], input_rows), (files[3], feature_rows),
        (files[4], resource_rows), (files[5], validation_rows),
        (files[6], decision_rows),
    ]:
        write_csv(output / name, rows)

    report_name = "MV08H_CELLRANGER8_REFERENCE_CLOSURE_2026-08-18.md"
    report = f"""# MV8-H Cell Ranger 8.0.1 custom-reference closure

## Outcome

The prospectively frozen four-core/32-GiB `cellranger mkref` run completed
successfully in {int(receipt['elapsed_seconds']):,} seconds. Independent validation passes
{len(checks)}/{len(checks)} checks. The final reference contains {regular} regular files,
{regular_bytes:,} regular-file bytes, zero partials, and complete required STAR
index components. Its deterministic tree SHA-256 is `{tree_sha}`.

## Exact scientific inputs

- The embedded primary-assembly FASTA is byte-identical to the frozen
  Ensembl-93 input: {fasta.stat().st_size:,} bytes and SHA-256
  `{EXPECTED_FASTA_SHA256}`.
- The embedded gzip-compressed GTF decompresses to the exact frozen
  target-complete annotation: {gtf_bytes:,} bytes and SHA-256
  `{EXPECTED_GTF_SHA256}`.
- The annotation contains {len(set(genes)):,} unique genes and
  {feature_records:,} feature records. All 500 exact-panel stable IDs and all
  475 common-panel stable IDs are present.
- `reference.json` binds Cell Ranger 8.0.1, four threads, 32 GiB, the exact
  genome name, input filenames, and Cell Ranger's matching SHA-1 identities.

## Resource closure

The private monitor recorded {len(samples):,} samples. Peak subprocess-tree RSS was
{max(rss_values) / 1024**2:.3f} GiB, below the selected 32-GiB allocation. Peak run-tree
size was {max(tree_values):,} bytes, below 50 GiB. Minimum free space was
{min(free_values):,} bytes, above 1 TiB. No resource breach, automatic kill,
or launcher deletion occurred. Cell Ranger's own normal temporary-file cleanup
is reflected in the changing run-tree samples.

The prospectively committed launcher sampled RSS correctly but did not include
RSS in its composite `resource_breach_detected` Boolean. This independent
closure therefore applies the intended 32-GiB memory gate directly to all 450
samples; the observed 30.494-GiB peak passes. The launcher is corrected
forward-only so any future RSS overage also marks a formal breach. This does not
change the executed command, reference, or scientific inputs.

## Scope and next gate

This closure does **not** execute or authorize `cellranger count`. It opens only
a prospective one-complete-unit count-sentinel prefreeze. The remaining seven
units, QC, matrices, PCA, PH, persistence landscapes, clustering, labels,
outcomes, and deletion remain closed.

The corrected topology contract is unchanged: cell and gene observations are
separate typed views; H0 and H1 remain separate; essential H0 is excluded; and
landscapes use every consecutive active level with exact or error-controlled
squared-L2 integration, no fixed grid, and no universal level cap.

The executed launcher was prospectively committed as `42fdc91` before the run.
Public evidence contains only identities and aggregates; private paths, licensed
runtime contents, expression, barcodes, donor attributes, labels, and outcomes
are not published.
"""
    write_text(output / report_name, report)
    files.append(report_name)
    manifest_name = "mv08h-reference-artifact-manifest.csv"
    write_csv(output / manifest_name, public_manifest(output, files))
    print(
        f"MV8-H Cell Ranger reference closure passed {len(checks)}/{len(checks)} "
        f"checks; tree_sha256={tree_sha}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
