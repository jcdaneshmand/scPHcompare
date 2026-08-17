#!/usr/bin/env python3
"""Build deterministic public MV8-H raw-read acquisition evidence.

Private compact HCA manifests are reduced to public file identities and unit
ownership. Donor identifiers, biological labels, signed URLs, expression
values, and cell barcodes are never emitted.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
from pathlib import Path
import re


EXPECTED_FASTQ_FILES = 48
EXPECTED_FASTQ_BYTES = 85_034_239_918
EXPECTED_GTF_SHA256 = "810e3bb63bb24bd5a005b14397d69280dda34c41a612fb86a18f3c4836fce57d"
EXPECTED_PANEL_SHA256 = "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e"
EXPECTED_COMMON475_SHA256 = "b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba"
EXPECTED_UNITS = [f"HCA_BM_{index:03d}" for index in range(1, 9)]
FASTA_URL = (
    "https://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/"
    "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
)
FASTA_BYTES = 881_214_682
FASTA_BSD_SUM = "07119"
FASTA_BSD_BLOCKS = 860_562
DOWNLOAD_RESERVE_BYTES = 1_649_267_441_664  # 1.5 TiB
RAW_CACHE_CAP_BYTES = 92_274_688_000  # safely above the 85.0-GB manifest
GTF_ATTRIBUTE = re.compile(r'(\S+) "([^"]*)";')
FASTQ_NAME = re.compile(
    r"^(MantonBM([1-8]))_HiSeq_9_S([1-8])_L00([23])_(I1|R1|R2)_001\.fastq\.gz$"
)
ALLOWED_BIOTYPES = {
    "protein_coding", "lincRNA", "antisense", "IG_LV_gene", "IG_V_gene",
    "IG_V_pseudogene", "IG_D_gene", "IG_J_gene", "IG_J_pseudogene",
    "IG_C_gene", "IG_C_pseudogene", "TR_V_gene", "TR_V_pseudogene",
    "TR_D_gene", "TR_J_gene", "TR_J_pseudogene", "TR_C_gene",
}


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


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty evidence: {path.name}")
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


def attributes(line: str) -> dict[str, str]:
    fields = line.rstrip("\n").split("\t")
    return dict(GTF_ATTRIBUTE.findall(fields[8])) if len(fields) == 9 else {}


def build_custom_gtf(
    source: Path, output: Path, target_stable_ids: set[str]
) -> tuple[list[str], int]:
    included: set[str] = set()
    ordered_ids: list[str] = []
    with gzip.open(source, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = attributes(line)
            gene_id = attrs.get("gene_id", "")
            biotype = attrs.get("gene_biotype", attrs.get("gene_type", ""))
            if gene_id and (biotype in ALLOWED_BIOTYPES or gene_id in target_stable_ids):
                included.add(gene_id)
                ordered_ids.append(gene_id)
    if len(ordered_ids) != 33_563 or len(included) != 33_563:
        raise ValueError("custom reference must contain exactly 33,563 unique genes")
    if not target_stable_ids.issubset(included):
        raise ValueError("custom reference does not retain all 500 target IDs")

    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_name(output.name + ".partial")
    included_lines = 0
    with gzip.open(source, "rt", encoding="utf-8") as source_handle, partial.open(
        "w", encoding="utf-8", newline=""
    ) as output_handle:
        for line in source_handle:
            if line.startswith("#"):
                output_handle.write(line)
                continue
            attrs = attributes(line)
            if attrs.get("gene_id") in included:
                output_handle.write(line)
                included_lines += 1
        output_handle.flush()
        os.fsync(output_handle.fileno())
    os.replace(partial, output)
    return ordered_ids, included_lines


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq-manifest", type=Path, required=True)
    parser.add_argument("--h5-compact-manifest", type=Path, required=True)
    parser.add_argument("--h5-public-manifest", type=Path, required=True)
    parser.add_argument("--unit-resource", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--missing-panel", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--ensembl-checksums", type=Path, required=True)
    parser.add_argument("--custom-gtf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if sha256_file(args.gtf) != EXPECTED_GTF_SHA256:
        raise ValueError("Ensembl-93 GTF identity changed")
    panel = read_csv(args.panel)
    if len(panel) != 500 or {row["panel_sha256"] for row in panel} != {EXPECTED_PANEL_SHA256}:
        raise ValueError("exact ordered 500-gene panel identity changed")
    panel.sort(key=lambda row: int(row["panel_order"]))
    if [int(row["panel_order"]) for row in panel] != list(range(1, 501)):
        raise ValueError("panel order is not exactly 1 through 500")
    target_ids = {re.search(r"(ENSG\d+)\.\d+$", row["feature_id"]).group(1) for row in panel}
    if len(target_ids) != 500:
        raise ValueError("panel does not contain 500 unique stable IDs")
    missing = read_csv(args.missing_panel)
    if len(missing) != 25 or {row["ensembl_stable_id"] for row in missing} - target_ids:
        raise ValueError("missing-25 identity changed")

    h5_private = read_csv(args.h5_compact_manifest, delimiter="\t")
    h5_public = read_csv(args.h5_public_manifest)
    if len(h5_private) != 8 or len(h5_public) != 8:
        raise ValueError("HCA H5 cohort is not exactly eight units")
    unit_by_h5_uuid = {row["file_uuid"]: row["unit_id"] for row in h5_public}
    donor_to_unit: dict[str, str] = {}
    for row in h5_private:
        unit = unit_by_h5_uuid.get(row["file_uuid"])
        donor = row["donor_organism.provenance.document_id"]
        if not unit or not donor or donor in donor_to_unit:
            raise ValueError("private-to-public H5 unit join is not one-to-one")
        donor_to_unit[donor] = unit
    if sorted(donor_to_unit.values()) != EXPECTED_UNITS:
        raise ValueError("H5 donor join does not reproduce the frozen unit axis")

    fastq = read_csv(args.fastq_manifest, delimiter="\t")
    if len(fastq) != EXPECTED_FASTQ_FILES:
        raise ValueError("official FASTQ manifest does not contain exactly 48 rows")
    public_fastq: list[dict[str, object]] = []
    seen_names: set[str] = set()
    seen_uuids: set[str] = set()
    role_order = {"I1": 1, "R1": 2, "R2": 3}
    read_index = {"I1": "index1", "R1": "read1", "R2": "read2"}
    prepared: list[tuple[int, int, int, dict[str, str], str]] = []
    for row in fastq:
        unit = donor_to_unit.get(row["donor_organism.provenance.document_id"])
        match = FASTQ_NAME.fullmatch(row["file_name"])
        if unit is None or match is None:
            raise ValueError("FASTQ row cannot be assigned to a frozen safe unit/name")
        sample_number = int(match.group(2))
        lane = int(match.group(4))
        role = match.group(5)
        if unit != f"HCA_BM_{sample_number:03d}" or int(match.group(3)) != sample_number:
            raise ValueError("FASTQ sample token disagrees with H5 donor ownership")
        if row["read_index"] != read_index[role] or row["file_format"] != "fastq.gz":
            raise ValueError("FASTQ role or format differs from the frozen read structure")
        if row["file_name"] in seen_names or row["file_uuid"] in seen_uuids:
            raise ValueError("FASTQ file names and UUIDs must be unique")
        seen_names.add(row["file_name"])
        seen_uuids.add(row["file_uuid"])
        prepared.append((sample_number, lane, role_order[role], row, role))
    prepared.sort(key=lambda item: item[:3])
    for file_order, (unit_order, lane, _, row, role) in enumerate(prepared, 1):
        version = row["file_version"]
        public_fastq.append({
            "contract_id": "mv08h_exact_hca_fastq_manifest_v1",
            "file_order": file_order,
            "unit_order": unit_order,
            "unit_id": f"HCA_BM_{unit_order:03d}",
            "lane": lane,
            "read_role": role,
            "read_index": row["read_index"],
            "file_name": row["file_name"],
            "file_uuid": row["file_uuid"],
            "file_version": version,
            "file_size_bytes": int(row["file_size"]),
            "file_sha256": row["file_sha256"].lower(),
            "file_crc32c": row["file_crc32c"].lower(),
            "file_content_type": row["file_content_type"],
            "drs_uri": row["file_drs_uri"],
            "azul_download_url": row["file_azul_url"],
            "mirror_uri": row["file_mirror_uri"],
            "download_authorized": "TRUE",
        })
    if sum(int(row["file_size_bytes"]) for row in public_fastq) != EXPECTED_FASTQ_BYTES:
        raise ValueError("official FASTQ total differs from frozen aggregate")

    unit_resource = {row["unit_id"]: row for row in read_csv(args.unit_resource)}
    unit_rows: list[dict[str, object]] = []
    for unit_order, unit_id in enumerate(EXPECTED_UNITS, 1):
        rows = [row for row in public_fastq if row["unit_id"] == unit_id]
        if len(rows) != 6 or {(row["lane"], row["read_role"]) for row in rows} != {
            (lane, role) for lane in (2, 3) for role in ("I1", "R1", "R2")
        }:
            raise ValueError(f"{unit_id} does not have the exact six-file read structure")
        total = sum(int(row["file_size_bytes"]) for row in rows)
        if total != int(unit_resource[unit_id]["fastq_bytes"]):
            raise ValueError(f"{unit_id} total disagrees with MV8-E")
        h5_row = next(row for row in h5_public if row["unit_id"] == unit_id)
        unit_rows.append({
            "contract_id": "mv08h_exact_hca_unit_manifest_v1",
            "unit_order": unit_order,
            "unit_id": unit_id,
            "fastq_files": 6,
            "fastq_bytes": total,
            "lanes": "L002;L003",
            "read_roles": "I1;R1;R2",
            "library_construction_approach": "10x 3' v2",
            "paired_end": "FALSE",
            "expected_cells": 7000,
            "source_h5_uuid": h5_row["file_uuid"],
            "biological_unit_preserved": "TRUE",
            "donor_identifier_published": "FALSE",
        })

    custom_axis, custom_lines = build_custom_gtf(args.gtf, args.custom_gtf, target_ids)
    checksums_text = args.ensembl_checksums.read_text(encoding="utf-8")
    checksum_line = f"{FASTA_BSD_SUM} {FASTA_BSD_BLOCKS} Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    if checksum_line not in checksums_text:
        raise ValueError("Ensembl FASTA checksum record changed")
    custom_sha = sha256_file(args.custom_gtf)
    reference_rows = [
        {
            "contract_id": "mv08h_custom_reference_contract_v1", "resource_order": 1,
            "resource_id": "ensembl93_gtf", "source_locator": "https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz",
            "source_bytes": args.gtf.stat().st_size, "source_sha256": EXPECTED_GTF_SHA256,
            "identity_state": "locally_sha256_verified", "derived_bytes": "", "derived_sha256": "", "gate": "pass",
        },
        {
            "contract_id": "mv08h_custom_reference_contract_v1", "resource_order": 2,
            "resource_id": "ensembl93_primary_assembly_fasta", "source_locator": FASTA_URL,
            "source_bytes": FASTA_BYTES, "source_sha256": "pending_local_acquisition",
            "identity_state": f"official_Ensembl_BSD_sum_{FASTA_BSD_SUM}_blocks_{FASTA_BSD_BLOCKS}",
            "derived_bytes": "", "derived_sha256": "", "gate": "stop_before_mkref_until_local_sha256_frozen",
        },
        {
            "contract_id": "mv08h_custom_reference_contract_v1", "resource_order": 3,
            "resource_id": "target_complete_33563_gtf", "source_locator": "derived_private_cache_from_ensembl93_gtf",
            "source_bytes": args.gtf.stat().st_size, "source_sha256": EXPECTED_GTF_SHA256,
            "identity_state": f"cellranger_3_0_0_allowlist_union_exact_missing25;genes={len(custom_axis)};lines={custom_lines};axis_sha256={sha256_axis(custom_axis)}",
            "derived_bytes": args.custom_gtf.stat().st_size, "derived_sha256": custom_sha, "gate": "pass_exact500_present",
        },
        {
            "contract_id": "mv08h_custom_reference_contract_v1", "resource_order": 4,
            "resource_id": "cellranger_mkref_output", "source_locator": "private_runtime_output",
            "source_bytes": "", "source_sha256": "pending_exact_runtime",
            "identity_state": "Cell_Ranger_3.0.0_mkref_name=GRCh38_Ensembl93_target_complete_33563",
            "derived_bytes": "", "derived_sha256": "pending_post_build_tree_manifest", "gate": "stop_before_count_until_independently_validated",
        },
    ]
    processing_rows = [
        {"contract_id": "mv08h_processing_contract_v1", "step_order": 1, "step_id": "counter_runtime", "frozen_value": "Cell Ranger 3.0.0 exact binary distribution", "status": "selected_not_installed", "execution_authorized": "FALSE"},
        {"contract_id": "mv08h_processing_contract_v1", "step_order": 2, "step_id": "chemistry", "frozen_value": "SC3Pv2; paired_end=false; lanes L002,L003; I1/R1/R2 retained", "status": "frozen_from_official_metadata", "execution_authorized": "FALSE"},
        {"contract_id": "mv08h_processing_contract_v1", "step_order": 3, "step_id": "mkref", "frozen_value": "cellranger mkref --genome=GRCh38_Ensembl93_target_complete_33563 --fasta=<verified Ensembl93 primary assembly> --genes=<verified target-complete GTF> --nthreads=16 --memgb=64", "status": "blocked_on_fasta_and_runtime_validation", "execution_authorized": "FALSE"},
        {"contract_id": "mv08h_processing_contract_v1", "step_order": 4, "step_id": "count", "frozen_value": "one unit at a time: cellranger count --id=<unit_id> --transcriptome=<verified custom reference> --fastqs=<verified unit FASTQs> --sample=MantonBM<unit>_HiSeq_9 --chemistry=SC3Pv2 --expect-cells=7000 --localcores=16 --localmem=64 --nosecondary", "status": "blocked_until_download_closure_and_runtime_reference_validation", "execution_authorized": "FALSE"},
        {"contract_id": "mv08h_processing_contract_v1", "step_order": 5, "step_id": "qc", "frozen_value": "reuse docs/specifications/mv08d-hca-qc-contract-v1.csv exactly", "status": "frozen", "execution_authorized": "FALSE"},
        {"contract_id": "mv08h_processing_contract_v1", "step_order": 6, "step_id": "cell_selection", "frozen_value": "sort post-QC barcodes; select 384 without replacement for seeds 20260805:20260809 using select_matched_cells", "status": "frozen", "execution_authorized": "FALSE"},
        {"contract_id": "mv08h_processing_contract_v1", "step_order": 7, "step_id": "topology", "frozen_value": "exact ordered 500; cell_topology_v1 and gene_topology_v1; complete VR H0/H1; every active landscape level; exact/error-controlled squared L2; no fixed grid or level cap; H0/H1 separate; essential H0 excluded", "status": "frozen_unchanged", "execution_authorized": "FALSE"},
    ]
    resource_rows = [
        {"contract_id": "mv08h_resource_caps_v1", "resource_order": 1, "resource": "FASTQ_manifest_bytes", "cap_bytes": EXPECTED_FASTQ_BYTES, "policy": "exact equality", "execution_scope": "download"},
        {"contract_id": "mv08h_resource_caps_v1", "resource_order": 2, "resource": "FASTQ_cache_including_partials", "cap_bytes": RAW_CACHE_CAP_BYTES, "policy": "hard stop above cap", "execution_scope": "download"},
        {"contract_id": "mv08h_resource_caps_v1", "resource_order": 3, "resource": "minimum_free_space_after_remaining_FASTQs", "cap_bytes": DOWNLOAD_RESERVE_BYTES, "policy": "hard stop below 1.5 TiB reserve", "execution_scope": "download"},
        {"contract_id": "mv08h_resource_caps_v1", "resource_order": 4, "resource": "download_concurrency", "cap_bytes": 1, "policy": "one file at a time; five bounded network attempts; resume only by HTTP Range", "execution_scope": "download"},
        {"contract_id": "mv08h_resource_caps_v1", "resource_order": 5, "resource": "custom_reference_tree", "cap_bytes": 53_687_091_200, "policy": "hard stop above 50 GiB", "execution_scope": "future_reference"},
        {"contract_id": "mv08h_resource_caps_v1", "resource_order": 6, "resource": "cellranger_per_unit_peak_rss", "cap_bytes": 85_899_345_920, "policy": "hard stop above 80 GiB; units serialized", "execution_scope": "future_count"},
        {"contract_id": "mv08h_resource_caps_v1", "resource_order": 7, "resource": "cellranger_per_unit_workspace", "cap_bytes": 214_748_364_800, "policy": "hard stop above 200 GiB; no deletion without owner approval", "execution_scope": "future_count"},
        {"contract_id": "mv08h_resource_caps_v1", "resource_order": 8, "resource": "cellranger_per_unit_elapsed", "cap_bytes": 86_400, "policy": "hard stop after 24 h; no scientific retry", "execution_scope": "future_count"},
        {"contract_id": "mv08h_resource_caps_v1", "resource_order": 9, "resource": "minimum_free_space_during_processing", "cap_bytes": 1_099_511_627_776, "policy": "hard stop below 1 TiB reserve", "execution_scope": "future_count"},
    ]
    gate_rows = [
        {"contract_id": "mv08h_gate_status_v1", "gate_order": 1, "gate": "owner_resource_authorization", "status": "pass", "evidence": "project owner: I authorize including the download (2026-08-17)"},
        {"contract_id": "mv08h_gate_status_v1", "gate_order": 2, "gate": "immutable_fastq_manifest", "status": "pass", "evidence": "48 UUID/version/size/SHA256 rows and eight-unit ownership frozen"},
        {"contract_id": "mv08h_gate_status_v1", "gate_order": 3, "gate": "custom_reference_identity", "status": "prefrozen", "evidence": "Ensembl93 GTF and target-complete GTF exact; FASTA local SHA required before mkref"},
        {"contract_id": "mv08h_gate_status_v1", "gate_order": 4, "gate": "exact500_feature_oracle", "status": "pass_for_annotation", "evidence": "all 500 ordered stable IDs occur in the 33,563-gene custom GTF; post-mkref and post-count oracles still required"},
        {"contract_id": "mv08h_gate_status_v1", "gate_order": 5, "gate": "aligner_counter_freeze", "status": "selected_not_installed", "evidence": "Cell Ranger 3.0.0/SC3Pv2 commands frozen; count remains unauthorized"},
        {"contract_id": "mv08h_gate_status_v1", "gate_order": 6, "gate": "biological_unit_preservation", "status": "pass", "evidence": "exact eight units; six FASTQs per unit; no pooling"},
        {"contract_id": "mv08h_gate_status_v1", "gate_order": 7, "gate": "qc_and_384_cell_rule", "status": "pass_for_prefreeze", "evidence": "MV8-D QC and five deterministic 384-cell seeds reused unchanged"},
        {"contract_id": "mv08h_gate_status_v1", "gate_order": 8, "gate": "resource_and_storage_caps", "status": "pass_for_download", "evidence": "exact manifest cap, serialized acquisition, and 1.5-TiB reserve"},
        {"contract_id": "mv08h_gate_status_v1", "gate_order": 9, "gate": "label_firewall", "status": "pass", "evidence": "no donor attributes or biological outcomes published/accessed; topology labels closed"},
        {"contract_id": "mv08h_gate_status_v1", "gate_order": 10, "gate": "independent_validation_and_stop_rule", "status": "prefrozen", "evidence": "validate download before reference; validate reference/runtime before count; stop before outcomes"},
    ]
    decision_rows = [{
        "contract_id": "mv08h_prefreeze_decision_v1",
        "decision": "authorize_exact_48_file_download_only",
        "biological_units": 8,
        "fastq_files": EXPECTED_FASTQ_FILES,
        "fastq_bytes": EXPECTED_FASTQ_BYTES,
        "fastq_gib": f"{EXPECTED_FASTQ_BYTES / 1024 ** 3:.6f}",
        "sentinel_files_first": 1,
        "download_authorized": "TRUE",
        "reference_download_authorized": "TRUE",
        "mkref_authorized": "FALSE",
        "raw_reprocessing_authorized": "FALSE",
        "pca_ph_landscape_authorized": "FALSE",
        "label_access_authorized": "FALSE",
        "biological_outcomes_authorized": "FALSE",
        "next_gate": "verified_48_file_download_closure_then_exact_runtime_reference_validation",
    }]
    source_rows = [
        {"contract_id": "mv08h_source_manifest_v1", "source_order": 1, "source_id": "official_hca_fastq_compact_manifest", "locator": "https://service.azul.data.humancellatlas.org/fetch/manifest/files; dcp60; frozen MV8-H filters", "bytes": args.fastq_manifest.stat().st_size, "sha256": sha256_file(args.fastq_manifest), "access_class": "private_local_metadata_hash_only", "contains_expression": "FALSE", "contains_donor_identifiers": "TRUE", "published": "FALSE"},
        {"contract_id": "mv08h_source_manifest_v1", "source_order": 2, "source_id": "published_fastq_manifest", "locator": "mv08h-fastq-manifest.csv", "bytes": "pending_after_write", "sha256": "pending_after_write", "access_class": "public_file_identity_metadata", "contains_expression": "FALSE", "contains_donor_identifiers": "FALSE", "published": "TRUE"},
        {"contract_id": "mv08h_source_manifest_v1", "source_order": 3, "source_id": "ensembl93_gtf", "locator": "official archive URL in reference contract", "bytes": args.gtf.stat().st_size, "sha256": EXPECTED_GTF_SHA256, "access_class": "public_reference", "contains_expression": "FALSE", "contains_donor_identifiers": "FALSE", "published": "FALSE"},
        {"contract_id": "mv08h_source_manifest_v1", "source_order": 4, "source_id": "ensembl93_fasta_checksums", "locator": "https://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/CHECKSUMS", "bytes": args.ensembl_checksums.stat().st_size, "sha256": sha256_file(args.ensembl_checksums), "access_class": "public_reference_metadata", "contains_expression": "FALSE", "contains_donor_identifiers": "FALSE", "published": "FALSE"},
    ]

    outputs = {
        "mv08h-fastq-manifest.csv": public_fastq,
        "mv08h-unit-manifest.csv": unit_rows,
        "mv08h-reference-contract.csv": reference_rows,
        "mv08h-processing-contract.csv": processing_rows,
        "mv08h-resource-caps.csv": resource_rows,
        "mv08h-gate-status.csv": gate_rows,
        "mv08h-decision.csv": decision_rows,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, rows in outputs.items():
        write_csv(args.output_dir / name, rows)
    source_rows[1]["bytes"] = (args.output_dir / "mv08h-fastq-manifest.csv").stat().st_size
    source_rows[1]["sha256"] = sha256_file(args.output_dir / "mv08h-fastq-manifest.csv")
    write_csv(args.output_dir / "mv08h-source-manifest.csv", source_rows)
    print(
        f"MV8-H prefreeze built: {len(public_fastq)} FASTQs; "
        f"{EXPECTED_FASTQ_BYTES} bytes; custom GTF {custom_sha}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
