#!/usr/bin/env python3
"""Build the public MV8-H Cell Ranger 8.0.1 runtime/reference prefreeze.

The builder reads private local runtime/reference inputs but publishes only
content identities, structural summaries, frozen commands, and gate states.
It never runs ``cellranger mkref`` or ``cellranger count``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import re
import subprocess


EXPECTED_RUNTIME_VERSION = "cellranger cellranger-8.0.1"
EXPECTED_FASTA_BYTES = 881_214_682
EXPECTED_FASTA_SHA256 = (
    "2a27436d44f0d6350f86894fbe5edec56faa5467028879784508041562406aa0"
)
EXPECTED_SOURCE_GTF_BYTES = 42_921_931
EXPECTED_SOURCE_GTF_SHA256 = (
    "810e3bb63bb24bd5a005b14397d69280dda34c41a612fb86a18f3c4836fce57d"
)
EXPECTED_CUSTOM_GTF_BYTES = 1_099_737_654
EXPECTED_CUSTOM_GTF_SHA256 = (
    "e28e4c4faf0dd76884d5e94c481fce2db43ad303968067c1276092a234727182"
)
EXPECTED_CUSTOM_GENE_COUNT = 33_563
EXPECTED_CUSTOM_FEATURE_RECORDS = 2_565_751
EXPECTED_CUSTOM_AXIS_SHA256 = (
    "0ea5f54698d79e05a7dd98eded2bf9558174d5214cfe684e8263b996c2a169a4"
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
        digest.update(value.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty evidence: {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            extrasaction="raise",
            quoting=csv.QUOTE_ALL,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def run_text(command: list[str]) -> str:
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    return (result.stdout + result.stderr).strip()


def distribution_tree_identity(root: Path) -> tuple[int, int, int, str]:
    """Hash regular files and symlink targets in relative-path order."""
    regular: list[Path] = []
    links: list[Path] = []
    for directory, dirnames, filenames in os.walk(root, followlinks=False):
        directory_path = Path(directory)
        retained_dirs: list[str] = []
        for name in dirnames:
            candidate = directory_path / name
            if candidate.is_symlink():
                links.append(candidate)
            else:
                retained_dirs.append(name)
        dirnames[:] = retained_dirs
        for name in filenames:
            candidate = directory_path / name
            if candidate.is_symlink():
                links.append(candidate)
            elif candidate.is_file():
                regular.append(candidate)

    digest = hashlib.sha256()
    regular_bytes = 0
    entries = sorted(
        [("F", path.relative_to(root).as_posix(), path) for path in regular]
        + [("L", path.relative_to(root).as_posix(), path) for path in links],
        key=lambda item: (item[1], item[0]),
    )
    for index, (kind, relative, path) in enumerate(entries, 1):
        if kind == "F":
            size = path.stat().st_size
            file_sha = sha256_file(path)
            regular_bytes += size
            record = f"F\0{relative}\0{size}\0{file_sha}\n"
        else:
            target = os.readlink(path)
            record = f"L\0{relative}\0{target}\n"
        digest.update(record.encode("utf-8"))
        if index % 2000 == 0:
            print(f"runtime tree identity: {index}/{len(entries)} entries", flush=True)
    return len(regular), len(links), regular_bytes, digest.hexdigest()


def parse_custom_gtf(path: Path) -> tuple[list[str], int]:
    genes: list[str] = []
    feature_records = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise RuntimeError("custom GTF contains a malformed feature record")
            feature_records += 1
            if fields[2] == "gene":
                attrs = dict(GTF_ATTRIBUTE.findall(fields[8]))
                genes.append(attrs.get("gene_id", ""))
    return genes, feature_records


def panel_ids(path: Path, expected_rows: int, order_field: str) -> tuple[list[str], str]:
    rows = read_csv(path)
    if len(rows) != expected_rows:
        raise RuntimeError(f"{path.name} has {len(rows)} rows, expected {expected_rows}")
    rows.sort(key=lambda row: int(row[order_field]))
    if [int(row[order_field]) for row in rows] != list(range(1, expected_rows + 1)):
        raise RuntimeError(f"{path.name} order is not consecutive")
    ids = []
    for row in rows:
        stable = row.get("ensembl_stable_id", "")
        if not stable:
            match = STABLE_ID.search(row["feature_id"])
            stable = match.group(1) if match else ""
        ids.append(stable)
    if len(set(ids)) != expected_rows or any(not value for value in ids):
        raise RuntimeError(f"{path.name} stable-ID axis is not unique and complete")
    declared = rows[0].get("panel_sha256", rows[0].get("parent_panel_sha256", ""))
    return ids, declared


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cellranger-root", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--source-gtf", type=Path, required=True)
    parser.add_argument("--custom-gtf", type=Path, required=True)
    parser.add_argument("--repeat-custom-gtf", type=Path, required=True)
    parser.add_argument("--panel500", type=Path, required=True)
    parser.add_argument("--panel475", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    root = args.cellranger_root.resolve()
    launcher = root / "bin" / "cellranger"
    star = root / "lib" / "bin" / "STAR"
    samtools = root / "lib" / "bin" / "samtools"
    python = root / "external" / "anaconda" / "bin" / "python3.10"
    version_file = root / ".version"
    required = [launcher, star, samtools, python, version_file, args.fasta,
                args.source_gtf, args.custom_gtf, args.repeat_custom_gtf,
                args.panel500, args.panel475]
    if not all(path.is_file() for path in required):
        missing = [str(path) for path in required if not path.is_file()]
        raise RuntimeError(f"required runtime/reference inputs are absent: {missing}")

    reported_version = run_text([str(launcher), "--version"])
    count_help = run_text([str(launcher), "count", "--help"])
    mkref_help = run_text([str(launcher), "mkref", "--help"])
    required_count_flags = [
        "--chemistry", "--include-introns", "--create-bam", "--nosecondary",
        "--expect-cells", "--localcores", "--localmem", "--transcriptome",
        "--fastqs", "--sample",
    ]
    required_mkref_flags = ["--genome", "--fasta", "--genes", "--nthreads", "--memgb"]
    cli_controls = all(flag in count_help for flag in required_count_flags)
    cli_controls = cli_controls and all(flag in mkref_help for flag in required_mkref_flags)
    if reported_version != EXPECTED_RUNTIME_VERSION or not cli_controls:
        raise RuntimeError("Cell Ranger version or required CLI controls differ from the freeze")

    regular_files, symlinks, regular_bytes, tree_sha = distribution_tree_identity(root)
    star_version = run_text([str(star), "--version"]).splitlines()[0]
    samtools_version = run_text([str(samtools), "--version"]).splitlines()[0]
    python_version = run_text([str(python), "--version"]).splitlines()[0]
    version_text = version_file.read_text(encoding="utf-8").strip()
    if version_text != "cellranger-8.0.1":
        raise RuntimeError("Cell Ranger .version file differs from the runtime report")

    fasta_exact = (
        args.fasta.stat().st_size == EXPECTED_FASTA_BYTES
        and sha256_file(args.fasta) == EXPECTED_FASTA_SHA256
    )
    source_gtf_exact = (
        args.source_gtf.stat().st_size == EXPECTED_SOURCE_GTF_BYTES
        and sha256_file(args.source_gtf) == EXPECTED_SOURCE_GTF_SHA256
    )
    custom_sha = sha256_file(args.custom_gtf)
    repeat_custom_sha = sha256_file(args.repeat_custom_gtf)
    genes, feature_records = parse_custom_gtf(args.custom_gtf)
    panel500_ids, declared_panel_sha = panel_ids(args.panel500, 500, "panel_order")
    panel475_ids, declared_parent_sha = panel_ids(args.panel475, 475, "panel_order_475")
    common475_declared = read_csv(args.panel475)[0]["common475_axis_sha256"]
    custom_exact = all([
        args.custom_gtf.stat().st_size == EXPECTED_CUSTOM_GTF_BYTES,
        custom_sha == EXPECTED_CUSTOM_GTF_SHA256,
        repeat_custom_sha == custom_sha,
        len(genes) == EXPECTED_CUSTOM_GENE_COUNT,
        len(set(genes)) == EXPECTED_CUSTOM_GENE_COUNT,
        feature_records == EXPECTED_CUSTOM_FEATURE_RECORDS,
        sha256_axis(genes) == EXPECTED_CUSTOM_AXIS_SHA256,
        set(panel500_ids).issubset(genes),
    ])
    panels_exact = all([
        declared_panel_sha == EXPECTED_PANEL_SHA256,
        declared_parent_sha == EXPECTED_PANEL_SHA256,
        common475_declared == EXPECTED_COMMON475_SHA256,
        set(panel475_ids).issubset(panel500_ids),
    ])
    if not all([fasta_exact, source_gtf_exact, custom_exact, panels_exact]):
        raise RuntimeError("one or more reference/panel identity gates failed")

    runtime_rows = [{
        "contract_id": "mv08h_cellranger8_runtime_identity_v2",
        "runtime": "Cell Ranger 8.0.1",
        "reported_version": reported_version,
        "distribution_label": root.name,
        "regular_files": regular_files,
        "symlinks": symlinks,
        "regular_file_bytes": regular_bytes,
        "tree_sha256": tree_sha,
        "launcher_relative_path": "bin/cellranger",
        "launcher_bytes": launcher.stat().st_size,
        "launcher_sha256": sha256_file(launcher),
        "version_file_sha256": sha256_file(version_file),
        "star_version": star_version,
        "star_sha256": sha256_file(star),
        "samtools_version": samtools_version,
        "samtools_sha256": sha256_file(samtools),
        "python_version": python_version,
        "python_sha256": sha256_file(python),
        "required_cli_controls_verified": "TRUE",
        "all_runtime_gates_passed": "TRUE",
    }]
    reference_rows = [
        {
            "contract_id": "mv08h_cellranger8_reference_binding_v2",
            "resource_order": 1,
            "resource_id": "ensembl93_primary_assembly_fasta_gz",
            "public_locator": "Ensembl release-93 primary-assembly FASTA archive",
            "bytes": args.fasta.stat().st_size,
            "sha256": sha256_file(args.fasta),
            "structural_detail": "gzip identity previously closed; uncompressed_bytes=3151425857",
            "gate_passed": "TRUE",
        },
        {
            "contract_id": "mv08h_cellranger8_reference_binding_v2",
            "resource_order": 2,
            "resource_id": "ensembl93_gtf_gz",
            "public_locator": "Ensembl release-93 Homo_sapiens GRCh38.93 GTF archive",
            "bytes": args.source_gtf.stat().st_size,
            "sha256": sha256_file(args.source_gtf),
            "structural_detail": "source for deterministic target-complete annotation",
            "gate_passed": "TRUE",
        },
        {
            "contract_id": "mv08h_cellranger8_reference_binding_v2",
            "resource_order": 3,
            "resource_id": "target_complete_33563_gtf",
            "public_locator": "private deterministic derivative; identity only",
            "bytes": args.custom_gtf.stat().st_size,
            "sha256": custom_sha,
            "structural_detail": (
                f"genes={len(genes)};feature_records={feature_records};"
                f"gene_axis_sha256={sha256_axis(genes)};all_exact500_present=true;"
                "independent_repeat_byte_identical=true"
            ),
            "gate_passed": "TRUE",
        },
        {
            "contract_id": "mv08h_cellranger8_reference_binding_v2",
            "resource_order": 4,
            "resource_id": "exact500_panel",
            "public_locator": "docs/audits/mv07h-prefreeze-evidence-v4/mv07h-panel.csv",
            "bytes": args.panel500.stat().st_size,
            "sha256": sha256_file(args.panel500),
            "structural_detail": f"features=500;panel_sha256={EXPECTED_PANEL_SHA256}",
            "gate_passed": "TRUE",
        },
        {
            "contract_id": "mv08h_cellranger8_reference_binding_v2",
            "resource_order": 5,
            "resource_id": "common475_panel",
            "public_locator": "docs/audits/mv08e-reference-reconciliation-evidence/mv08e-common475-panel.csv",
            "bytes": args.panel475.stat().st_size,
            "sha256": sha256_file(args.panel475),
            "structural_detail": f"features=475;axis_sha256={EXPECTED_COMMON475_SHA256}",
            "gate_passed": "TRUE",
        },
    ]
    processing_rows = [
        {
            "contract_id": "mv08h_cellranger8_processing_amendment_v2",
            "step_order": 1, "step_id": "runtime",
            "frozen_value": f"Cell Ranger 8.0.1 exact installed tree_sha256={tree_sha}",
            "rationale": "prospective modern raw-read validation; historical 3.0.0 matrices retained as comparator",
            "execution_authorized": "FALSE",
        },
        {
            "contract_id": "mv08h_cellranger8_processing_amendment_v2",
            "step_order": 2, "step_id": "mkref",
            "frozen_value": "cellranger mkref --genome=GRCh38_Ensembl93_target_complete_33563 --fasta=<verified uncompressed Ensembl93 primary assembly> --genes=<verified target-complete GTF> --nthreads=16 --memgb=64",
            "rationale": "same Ensembl93 genome and minimal 33563-gene annotation; restore all exact 500 targets",
            "execution_authorized": "TRUE",
        },
        {
            "contract_id": "mv08h_cellranger8_processing_amendment_v2",
            "step_order": 3, "step_id": "count_sentinel",
            "frozen_value": "one complete HCA unit: cellranger count --id=<unit_id> --transcriptome=<verified custom reference> --fastqs=<verified unit FASTQs> --sample=MantonBM<unit>_HiSeq_9 --chemistry=SC3Pv2 --expect-cells=7000 --include-introns=false --create-bam=false --nosecondary --localcores=16 --localmem=64",
            "rationale": "explicit exon-only historical alignment target; omit BAM/secondary outputs without changing feature-barcode counts",
            "execution_authorized": "FALSE",
        },
        {
            "contract_id": "mv08h_cellranger8_processing_amendment_v2",
            "step_order": 4, "step_id": "panel_separation",
            "frozen_value": "derive exact ordered 500 and ordered common 475 from the same Cell Ranger 8.0.1 count matrices; compare 8.0.1 common475 to historical 3.0.0 common475 separately",
            "rationale": "isolate panel effect within runtime and quantify runtime effect on a shared panel",
            "execution_authorized": "FALSE",
        },
        {
            "contract_id": "mv08h_cellranger8_processing_amendment_v2",
            "step_order": 5, "step_id": "qc_cell_selection",
            "frozen_value": "reuse mv08d-hca-qc-contract-v1.csv exactly; sorted eligible barcodes; five deterministic 384-cell selections using seeds 20260805:20260809",
            "rationale": "no data-dependent QC or selection amendment",
            "execution_authorized": "FALSE",
        },
        {
            "contract_id": "mv08h_cellranger8_processing_amendment_v2",
            "step_order": 6, "step_id": "topology_landscapes",
            "frozen_value": "cell_topology_v1 and gene_topology_v1 separately; complete VR H0/H1; every consecutive active landscape level; exact or error-controlled squared L2; no fixed grid; no level cap; H0/H1 separate; essential H0 excluded",
            "rationale": "preserve the corrected dissertation-modernization landscape contract unchanged",
            "execution_authorized": "FALSE",
        },
    ]
    gate_rows = [
        {"contract_id": "mv08h_cellranger8_gate_status_v2", "gate_order": 1, "gate": "owner_runtime_amendment", "status": "pass", "evidence": "owner identified the installed Cell Ranger runtime and authorized this goal"},
        {"contract_id": "mv08h_cellranger8_gate_status_v2", "gate_order": 2, "gate": "fastq_acquisition", "status": "pass", "evidence": "48/48 files; 85034239918 bytes; independent complete checksum closure"},
        {"contract_id": "mv08h_cellranger8_gate_status_v2", "gate_order": 3, "gate": "runtime_identity", "status": "pass", "evidence": f"Cell Ranger 8.0.1 tree/executable/components bound; tree_sha256={tree_sha}"},
        {"contract_id": "mv08h_cellranger8_gate_status_v2", "gate_order": 4, "gate": "runtime_controls", "status": "pass", "evidence": "local CLI exposes SC3Pv2, include-introns, create-bam, nosecondary, compute, and mkref controls"},
        {"contract_id": "mv08h_cellranger8_gate_status_v2", "gate_order": 5, "gate": "reference_inputs", "status": "pass", "evidence": "Ensembl93 FASTA/source GTF/custom GTF identities exact; independent GTF repeat exact"},
        {"contract_id": "mv08h_cellranger8_gate_status_v2", "gate_order": 6, "gate": "exact500_annotation", "status": "pass", "evidence": "33563 unique genes; all exact 500 present; common475 is a strict ordered subset"},
        {"contract_id": "mv08h_cellranger8_gate_status_v2", "gate_order": 7, "gate": "panel_runtime_separation", "status": "pass", "evidence": "within-8.0.1 500-vs-475 and shared-475 historical-runtime comparisons frozen separately"},
        {"contract_id": "mv08h_cellranger8_gate_status_v2", "gate_order": 8, "gate": "resource_policy", "status": "pass_for_mkref", "evidence": "16 cores; 64 GiB request; 50-GiB reference cap; 1-TiB free-space floor; no deletion"},
        {"contract_id": "mv08h_cellranger8_gate_status_v2", "gate_order": 9, "gate": "scientific_scope", "status": "pass_for_prefreeze", "evidence": "QC, 384-cell seeds, typed views, H0/H1, and all-active-level landscapes unchanged"},
        {"contract_id": "mv08h_cellranger8_gate_status_v2", "gate_order": 10, "gate": "label_firewall", "status": "pass", "evidence": "no donor attributes, cell barcodes, expression, labels, or outcomes read or published"},
        {"contract_id": "mv08h_cellranger8_gate_status_v2", "gate_order": 11, "gate": "execution_stop", "status": "pass", "evidence": "mkref not executed; count sentinel and all downstream work remain closed"},
    ]
    decision_rows = [{
        "contract_id": "mv08h_cellranger8_prefreeze_decision_v2",
        "decision": "runtime_reference_exact_authorize_mkref_only",
        "runtime_amendment_accepted": "TRUE",
        "runtime_validated": "TRUE",
        "reference_inputs_validated": "TRUE",
        "mkref_authorized": "TRUE",
        "count_sentinel_authorized": "FALSE",
        "remaining_units_authorized": "FALSE",
        "pca_ph_landscape_authorized": "FALSE",
        "label_access_authorized": "FALSE",
        "biological_outcomes_authorized": "FALSE",
        "next_gate": "build_and_independently_validate_custom_reference_then_prefreeze_one_complete_unit_count_sentinel",
    }]
    source_rows = [
        {"contract_id": "mv08h_cellranger8_source_manifest_v2", "source_order": 1, "source_id": "installed_cellranger_8_0_1", "locator": "private local installation; relative identities only", "sha256": tree_sha, "access_class": "licensed_runtime_identity", "published_content": "identity_only"},
        {"contract_id": "mv08h_cellranger8_source_manifest_v2", "source_order": 2, "source_id": "ensembl93_reference_input_closure", "locator": "docs/audits/mv08h-reference-input-v1/", "sha256": EXPECTED_FASTA_SHA256, "access_class": "public_reference_identity", "published_content": "identity_only"},
        {"contract_id": "mv08h_cellranger8_source_manifest_v2", "source_order": 3, "source_id": "mv08h_acquisition_prefreeze_v1", "locator": "docs/audits/mv08h-raw-read-prefreeze-v1/", "sha256": "historical_v1_artifact_bundle_retained", "access_class": "public_prior_contract", "published_content": "contract_only"},
        {"contract_id": "mv08h_cellranger8_source_manifest_v2", "source_order": 4, "source_id": "tenx_cellranger_8_release_notes", "locator": "https://www.10xgenomics.com/support/software/cell-ranger/8.0/release-notes/cr-release-notes", "sha256": "not_local_payload", "access_class": "public_authoritative_documentation", "published_content": "citation_only"},
        {"contract_id": "mv08h_cellranger8_source_manifest_v2", "source_order": 5, "source_id": "tenx_command_arguments", "locator": "https://www.10xgenomics.com/support/software/cell-ranger/latest/resources/cr-command-line-arguments", "sha256": "not_local_payload", "access_class": "public_authoritative_documentation", "published_content": "citation_only"},
    ]

    outputs = {
        "mv08h-runtime-identity.csv": runtime_rows,
        "mv08h-reference-binding.csv": reference_rows,
        "mv08h-processing-amendment.csv": processing_rows,
        "mv08h-gate-status.csv": gate_rows,
        "mv08h-decision.csv": decision_rows,
        "mv08h-source-manifest.csv": source_rows,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, rows in outputs.items():
        write_csv(args.output_dir / name, rows)
    print(
        f"MV8-H Cell Ranger 8.0.1 prefreeze built: tree {tree_sha}; "
        f"{len(genes)} genes; mkref only opened",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
