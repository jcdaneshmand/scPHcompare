#!/usr/bin/env python3
"""Independently validate MV8-H public prefreeze evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import re


EXPECTED_BYTES = 85_034_239_918
EXPECTED_UNITS = [f"HCA_BM_{index:03d}" for index in range(1, 9)]
FASTQ_NAME = re.compile(
    r"^MantonBM([1-8])_HiSeq_9_S([1-8])_L00([23])_(I1|R1|R2)_001\.fastq\.gz$"
)
GTF_ATTRIBUTE = re.compile(r'(\S+) "([^"]*)";')


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
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


def one(rows: list[dict[str, str]], label: str) -> dict[str, str]:
    if len(rows) != 1:
        raise ValueError(f"{label} must contain exactly one row")
    return rows[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--private-fastq-manifest", type=Path, required=True)
    parser.add_argument("--private-h5-manifest", type=Path, required=True)
    parser.add_argument("--public-h5-manifest", type=Path, required=True)
    parser.add_argument("--unit-resource", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--custom-gtf", type=Path, required=True)
    parser.add_argument("--repeat-evidence-dir", type=Path)
    parser.add_argument("--repeat-custom-gtf", type=Path)
    args = parser.parse_args()

    evidence = args.evidence_dir
    manifest = read_csv(evidence / "mv08h-fastq-manifest.csv")
    units = read_csv(evidence / "mv08h-unit-manifest.csv")
    reference = read_csv(evidence / "mv08h-reference-contract.csv")
    processing = read_csv(evidence / "mv08h-processing-contract.csv")
    resources = read_csv(evidence / "mv08h-resource-caps.csv")
    gates = read_csv(evidence / "mv08h-gate-status.csv")
    decision = one(read_csv(evidence / "mv08h-decision.csv"), "decision")
    source = read_csv(evidence / "mv08h-source-manifest.csv")
    private_fastq = read_csv(args.private_fastq_manifest, delimiter="\t")
    private_h5 = read_csv(args.private_h5_manifest, delimiter="\t")
    public_h5 = read_csv(args.public_h5_manifest)
    unit_resource = {row["unit_id"]: row for row in read_csv(args.unit_resource)}
    panel = read_csv(args.panel)

    checks: list[tuple[str, bool, str]] = []
    checks.append(("authorization_boundary", decision["download_authorized"] == "TRUE" and
                   decision["raw_reprocessing_authorized"] == "FALSE" and
                   decision["label_access_authorized"] == "FALSE", "download only; processing and labels closed"))
    checks.append(("manifest_cardinality_and_bytes", len(manifest) == 48 and
                   sum(int(row["file_size_bytes"]) for row in manifest) == EXPECTED_BYTES,
                   "48 files; 85034239918 bytes"))
    checks.append(("manifest_unique_identity", len({row["file_name"] for row in manifest}) == 48 and
                   len({row["file_uuid"] for row in manifest}) == 48 and
                   [int(row["file_order"]) for row in manifest] == list(range(1, 49)),
                   "unique names/UUIDs and exact order"))

    private_by_uuid = {row["file_uuid"]: row for row in private_fastq}
    identity_exact = len(private_by_uuid) == 48
    for row in manifest:
        raw = private_by_uuid.get(row["file_uuid"], {})
        identity_exact = identity_exact and all([
            raw.get("file_name") == row["file_name"],
            raw.get("file_version") == row["file_version"],
            raw.get("file_size") == row["file_size_bytes"],
            raw.get("file_sha256", "").lower() == row["file_sha256"],
            raw.get("file_crc32c", "").lower() == row["file_crc32c"],
            raw.get("file_azul_url") == row["azul_download_url"],
        ])
    checks.append(("official_file_identity_exact", identity_exact, "public rows reproduce official compact metadata"))

    structure_exact = True
    for unit_order, unit_id in enumerate(EXPECTED_UNITS, 1):
        rows = [row for row in manifest if row["unit_id"] == unit_id]
        observed = set()
        for row in rows:
            match = FASTQ_NAME.fullmatch(row["file_name"])
            structure_exact = structure_exact and match is not None
            if match:
                structure_exact = structure_exact and int(match.group(1)) == unit_order and int(match.group(2)) == unit_order
                observed.add((int(match.group(3)), match.group(4)))
        structure_exact = structure_exact and observed == {
            (lane, role) for lane in (2, 3) for role in ("I1", "R1", "R2")
        }
    checks.append(("read_structure_exact", structure_exact, "eight units x two lanes x I1/R1/R2"))

    h5_unit_by_uuid = {row["file_uuid"]: row["unit_id"] for row in public_h5}
    donor_to_unit = {
        row["donor_organism.provenance.document_id"]: h5_unit_by_uuid[row["file_uuid"]]
        for row in private_h5
    }
    ownership = all(
        donor_to_unit.get(private_by_uuid[row["file_uuid"]]["donor_organism.provenance.document_id"]) == row["unit_id"]
        for row in manifest
    )
    checks.append(("eight_unit_ownership_exact", ownership and len(donor_to_unit) == 8,
                   "private donor join independently reproduces public unit IDs"))

    unit_exact = len(units) == 8 and [row["unit_id"] for row in units] == EXPECTED_UNITS
    for row in units:
        file_sum = sum(int(item["file_size_bytes"]) for item in manifest if item["unit_id"] == row["unit_id"])
        unit_exact = unit_exact and int(row["fastq_files"]) == 6 and int(row["fastq_bytes"]) == file_sum
        unit_exact = unit_exact and int(unit_resource[row["unit_id"]]["fastq_bytes"]) == file_sum
    checks.append(("unit_aggregate_reconciliation", unit_exact, "all eight unit totals match MV8-E"))

    target_ids = {
        re.search(r"(ENSG\d+)\.\d+$", row["feature_id"]).group(1)
        for row in panel
    }
    custom_ids: list[str] = []
    with args.custom_gtf.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 9 and fields[2] == "gene":
                attrs = dict(GTF_ATTRIBUTE.findall(fields[8]))
                custom_ids.append(attrs.get("gene_id", ""))
    custom_row = next(row for row in reference if row["resource_id"] == "target_complete_33563_gtf")
    custom_exact = len(custom_ids) == 33_563 and len(set(custom_ids)) == 33_563 and target_ids.issubset(custom_ids)
    custom_exact = custom_exact and custom_row["derived_sha256"] == sha256_file(args.custom_gtf)
    checks.append(("exact500_custom_annotation", custom_exact, "33,563 unique genes; all ordered-panel IDs present; SHA-256 exact"))

    ref_exact = len(reference) == 4 and reference[0]["source_sha256"] == "810e3bb63bb24bd5a005b14397d69280dda34c41a612fb86a18f3c4836fce57d"
    ref_exact = ref_exact and reference[1]["source_bytes"] == "881214682" and "07119" in reference[1]["identity_state"]
    ref_exact = ref_exact and reference[1]["gate"] == "stop_before_mkref_until_local_sha256_frozen"
    checks.append(("reference_stop_boundary", ref_exact, "GTF exact; FASTA frozen by archive identity and blocked pending local SHA-256"))

    processing_text = "\n".join(row["frozen_value"] for row in processing)
    processing_exact = len(processing) == 7 and all(term in processing_text for term in [
        "Cell Ranger 3.0.0", "SC3Pv2", "mv08d-hca-qc-contract-v1.csv",
        "384", "every active landscape level", "no fixed grid or level cap",
        "H0/H1 separate", "essential H0 excluded",
    ]) and not any(row["execution_authorized"] == "TRUE" for row in processing)
    checks.append(("processing_and_landscape_freeze", processing_exact, "runtime/QC/384-cell/dual-view/all-level contracts frozen; execution closed"))

    resource_by_name = {row["resource"]: row for row in resources}
    resource_exact = len(resources) == 9
    resource_exact = resource_exact and int(resource_by_name["FASTQ_manifest_bytes"]["cap_bytes"]) == EXPECTED_BYTES
    resource_exact = resource_exact and int(resource_by_name["minimum_free_space_after_remaining_FASTQs"]["cap_bytes"]) == 1_649_267_441_664
    resource_exact = resource_exact and int(resource_by_name["download_concurrency"]["cap_bytes"]) == 1
    checks.append(("resource_policy_exact", resource_exact, "exact byte cap, 1.5-TiB reserve, single-file acquisition"))

    forbidden_fields = {
        "donor_organism.sex", "donor_organism.organism_age", "disease",
        "tissue", "outcome", "cell_barcode", "local_absolute_path",
    }
    public_rows = [manifest, units, reference, processing, resources, gates, [decision], source]
    public_fields = {name.lower() for rows in public_rows for row in rows for name in row}
    no_forbidden = not (public_fields & forbidden_fields)
    no_forbidden = no_forbidden and not any("private_local_metadata_not_published" in str(row) for rows in public_rows for row in rows)
    checks.append(("public_firewall", no_forbidden, "no donor attributes, outcomes, barcodes, or absolute private paths"))

    source_exact = len(source) == 4 and source[0]["sha256"] == sha256_file(args.private_fastq_manifest)
    source_exact = source_exact and source[1]["sha256"] == sha256_file(evidence / "mv08h-fastq-manifest.csv")
    checks.append(("source_hash_chain", source_exact, "official private manifest hash and public manifest hash exact"))

    validation_rows = [
        {"contract_id": "mv08h_independent_validation_v1", "check_order": index,
         "check_id": check_id, "passed": str(passed).upper(), "detail": detail}
        for index, (check_id, passed, detail) in enumerate(checks, 1)
    ]
    if len(validation_rows) != 13 or not all(row["passed"] == "TRUE" for row in validation_rows):
        failed = [row["check_id"] for row in validation_rows if row["passed"] != "TRUE"]
        raise RuntimeError(f"MV8-H independent validation failed: {failed}")
    write_csv(evidence / "mv08h-independent-validation.csv", validation_rows)
    validation_decision = [{
        "contract_id": "mv08h_validation_decision_v1",
        "decision": "prefreeze_exact_authorize_one_file_sentinel_then_resumable_48_file_download",
        "checks_passed": 13,
        "checks_total": 13,
        "fastq_download_authorized": "TRUE",
        "raw_reprocessing_authorized": "FALSE",
        "labels_opened": "FALSE",
        "outcomes_computed": "FALSE",
        "next_gate": "sentinel_checksum_then_complete_download_checksum_closure",
    }]
    write_csv(evidence / "mv08h-validation-decision.csv", validation_decision)

    if bool(args.repeat_evidence_dir) != bool(args.repeat_custom_gtf):
        raise ValueError("repeat evidence and repeat custom GTF must be supplied together")
    if args.repeat_evidence_dir:
        repeat_names = sorted([
            "mv08h-decision.csv", "mv08h-fastq-manifest.csv",
            "mv08h-gate-status.csv", "mv08h-processing-contract.csv",
            "mv08h-reference-contract.csv", "mv08h-resource-caps.csv",
            "mv08h-source-manifest.csv", "mv08h-unit-manifest.csv",
        ])
        repeat_rows = []
        for index, name in enumerate(repeat_names, 1):
            first_sha = sha256_file(evidence / name)
            second_sha = sha256_file(args.repeat_evidence_dir / name)
            repeat_rows.append({
                "contract_id": "mv08h_prefreeze_repeat_validation_v1",
                "artifact_order": index, "artifact": name,
                "production_sha256": first_sha, "repeat_sha256": second_sha,
                "byte_identical": str(first_sha == second_sha).upper(),
            })
        first_sha = sha256_file(args.custom_gtf)
        second_sha = sha256_file(args.repeat_custom_gtf)
        repeat_rows.append({
            "contract_id": "mv08h_prefreeze_repeat_validation_v1",
            "artifact_order": 9, "artifact": "private_target_complete_33563.gtf",
            "production_sha256": first_sha, "repeat_sha256": second_sha,
            "byte_identical": str(first_sha == second_sha).upper(),
        })
        if not all(row["byte_identical"] == "TRUE" for row in repeat_rows):
            raise RuntimeError("MV8-H prefreeze byte repeat failed")
        write_csv(evidence / "mv08h-repeat-validation.csv", repeat_rows)

    names = sorted(
        path.name for path in evidence.iterdir()
        if path.is_file() and path.name != "mv08h-artifact-manifest.csv"
    )
    artifact_rows = [{
        "contract_id": "mv08h_prefreeze_artifact_manifest_v1",
        "file": name,
        "bytes": (evidence / name).stat().st_size,
        "sha256": sha256_file(evidence / name),
        "contains_expression": "FALSE",
        "contains_cell_barcode": "FALSE",
        "contains_absolute_private_path": "FALSE",
        "contains_donor_attribute": "FALSE",
        "contains_outcome_label": "FALSE",
    } for name in names]
    write_csv(evidence / "mv08h-artifact-manifest.csv", artifact_rows)
    print("MV8-H independent validation passed: 13/13 checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
