#!/usr/bin/env python3
"""Build public MV8-E annotation/reference reconciliation evidence.

The Ensembl GTF, HCA H5 count matrix, HCA search response, and current-Ensembl
lookup response are read-only private/ignored inputs.  Outputs contain only
annotation facts, aggregate file sizes, stable public identifiers, hashes, and
prospective decisions; no expression values or cell barcodes are exported.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
from typing import Iterable

import h5py


EXPECTED_GTF_SHA256 = (
    "810e3bb63bb24bd5a005b14397d69280dda34c41a612fb86a18f3c4836fce57d"
)
EXPECTED_PANEL_SHA256 = (
    "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e"
)
EXPECTED_HCA_FEATURE_AXIS_SHA256 = (
    "a6a914db1d218e2c3ae4c6680ac8fc546ba52ee4cfaedf423eb635c826c085be"
)
EXPECTED_FASTQ_BYTES = 85_034_239_918
EXPECTED_FASTQ_FILES = 48
EXPECTED_UNITS = [f"HCA_BM_{index:03d}" for index in range(1, 9)]
ENSEMBL_AT_END = re.compile(r"(ENSG\d+)\.(\d+)$")
GTF_ATTRIBUTE = re.compile(r'(\S+) "([^"]*)";')

# Exact Cell Ranger 3.0.0 GRCh38 mkgtf biotype allow-list documented by 10x.
ALLOWED_BIOTYPES = {
    "protein_coding", "lincRNA", "antisense", "IG_LV_gene", "IG_V_gene",
    "IG_V_pseudogene", "IG_D_gene", "IG_J_gene", "IG_J_pseudogene",
    "IG_C_gene", "IG_C_pseudogene", "TR_V_gene", "TR_V_pseudogene",
    "TR_D_gene", "TR_J_gene", "TR_J_pseudogene", "TR_C_gene",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def axis_sha256(values: Iterable[str]) -> str:
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
            handle, fieldnames=list(rows[0]), extrasaction="raise",
            quoting=csv.QUOTE_ALL, lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def parse_gtf(path: Path) -> list[dict[str, object]]:
    genes: list[dict[str, object]] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = dict(GTF_ATTRIBUTE.findall(fields[8]))
            gene_id = attributes.get("gene_id", "")
            gene_version = attributes.get("gene_version", "")
            gene_name = attributes.get("gene_name", "")
            gene_biotype = attributes.get(
                "gene_biotype", attributes.get("gene_type", "")
            )
            if not all((gene_id, gene_version, gene_name, gene_biotype)):
                raise ValueError("Ensembl-93 gene row lacks required annotation")
            genes.append({
                "gene_id": gene_id,
                "gene_version": gene_version,
                "gene_name": gene_name,
                "gene_biotype": gene_biotype,
            })
    if len(genes) != 58_395:
        raise ValueError(f"unexpected Ensembl-93 gene count: {len(genes)}")
    return genes


def decode_h5_axis(dataset: h5py.Dataset) -> list[str]:
    values = dataset[()]
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def donor_count(hit: dict[str, object]) -> int:
    organisms = hit.get("donorOrganisms", [])
    if not isinstance(organisms, list) or len(organisms) != 1:
        return -1
    return int(organisms[0].get("donorCount", -1))


def file_summary(hit: dict[str, object], file_format: str) -> dict[str, object]:
    matches = [
        value for value in hit.get("fileTypeSummaries", [])
        if value.get("format") == file_format
    ]
    if len(matches) != 1:
        raise ValueError(
            f"entry {hit.get('entryId')} has {len(matches)} {file_format} summaries"
        )
    return matches[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--h5", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--missing-summary", type=Path, required=True)
    parser.add_argument("--current-crosswalk", type=Path, required=True)
    parser.add_argument("--hca-sample-metadata", type=Path, required=True)
    parser.add_argument("--hca-file-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if sha256_file(args.gtf) != EXPECTED_GTF_SHA256:
        raise ValueError("Ensembl-93 GTF SHA-256 differs from the frozen source")
    genes = parse_gtf(args.gtf)
    by_id = {str(row["gene_id"]): row for row in genes}
    if len(by_id) != len(genes):
        raise ValueError("Ensembl-93 gene IDs are duplicated")
    filtered = [
        row for row in genes if str(row["gene_biotype"]) in ALLOWED_BIOTYPES
    ]
    gtf_ids = [str(row["gene_id"]) for row in filtered]
    gtf_names = [str(row["gene_name"]) for row in filtered]

    with h5py.File(args.h5, "r") as handle:
        h5_ids = decode_h5_axis(handle["matrix/features/id"])
        h5_names = decode_h5_axis(handle["matrix/features/name"])
    if len(h5_ids) != 33_538 or len(h5_names) != 33_538:
        raise ValueError("HCA feature axis is not the expected 33,538 genes")
    if axis_sha256(h5_ids) != EXPECTED_HCA_FEATURE_AXIS_SHA256:
        raise ValueError("HCA feature axis differs from the MV8-D exact axis")
    if gtf_ids != h5_ids or gtf_names != h5_names:
        raise ValueError("HCA axis is not exact Cell Ranger 3.0.0/Ensembl-93")

    panel = read_csv(args.panel)
    panel.sort(key=lambda row: int(row["panel_order"]))
    if (
        len(panel) != 500
        or [int(row["panel_order"]) for row in panel] != list(range(1, 501))
        or {row["panel_sha256"] for row in panel} != {EXPECTED_PANEL_SHA256}
    ):
        raise ValueError("MV7-H panel identity differs from the accepted 500")

    parsed_panel: list[tuple[dict[str, str], str, str]] = []
    for row in panel:
        match = ENSEMBL_AT_END.search(row["feature_id"])
        if not match:
            raise ValueError(f"panel feature lacks versioned Ensembl ID: {row}")
        parsed_panel.append((row, match.group(1), match.group(2)))
    panel_ids = [stable_id for _, stable_id, _ in parsed_panel]
    filtered_set = set(gtf_ids)
    if sum(stable_id in filtered_set for stable_id in panel_ids) != 475:
        raise ValueError("accepted panel does not reproduce the exact 475 intersection")
    common475_axis = [
        row["feature_id"] for row, stable_id, _ in parsed_panel
        if stable_id in filtered_set
    ]
    common475_sha256 = axis_sha256(common475_axis)
    common475_rows: list[dict[str, object]] = []
    common_order = 0
    for row, stable_id, _ in parsed_panel:
        if stable_id not in filtered_set:
            continue
        common_order += 1
        common475_rows.append({
            "contract_id": "mv08e_common475_panel_v1",
            "panel_order_475": common_order,
            "panel_order_500": int(row["panel_order"]),
            "feature_id": row["feature_id"],
            "gene": row["gene"],
            "ensembl_stable_id": stable_id,
            "parent_panel_sha256": EXPECTED_PANEL_SHA256,
            "common475_axis_sha256": common475_sha256,
            "selection_rule": "ordered_exact_stable_id_intersection_with_hca_cellranger_3_0_0_reference",
            "hca_expression_used_for_selection": "FALSE",
            "biological_outcomes_used_for_selection": "FALSE",
        })
    if common_order != 475:
        raise ValueError("common panel construction did not emit 475 genes")

    missing_source = read_csv(args.missing_summary)
    missing_source.sort(key=lambda row: int(row["panel_order"]))
    missing_ids = [row["ensembl_stable_id"] for row in missing_source]
    expected_missing = [value for value in panel_ids if value not in filtered_set]
    if len(missing_source) != 25 or missing_ids != expected_missing:
        raise ValueError("MV8-D missing-gene axis differs from Ensembl-93 filtering")

    current = read_csv(args.current_crosswalk)
    current.sort(key=lambda row: int(row["panel_order"]))
    if (
        [row["reference_ensembl_stable_id"] for row in current] != missing_ids
        or not all(row["current_lookup_present"].lower() == "true" for row in current)
    ):
        raise ValueError("current Ensembl lookup does not resolve all 25 stable IDs")
    current_by_id = {
        row["reference_ensembl_stable_id"]: row for row in current
    }

    missing_rows: list[dict[str, object]] = []
    for source in missing_source:
        annotation = by_id[source["ensembl_stable_id"]]
        now = current_by_id[source["ensembl_stable_id"]]
        biotype = str(annotation["gene_biotype"])
        if biotype in ALLOWED_BIOTYPES:
            raise ValueError("a missing gene unexpectedly passes the 10x filter")
        missing_rows.append({
            "contract_id": "mv08e_missing_biotype_reconciliation_v1",
            "panel_order": int(source["panel_order"]),
            "panel_feature_id": source["reference_feature_id"],
            "panel_gene": source["reference_gene"],
            "ensembl_stable_id": source["ensembl_stable_id"],
            "ensembl93_gene_version": annotation["gene_version"],
            "ensembl93_gene_name": annotation["gene_name"],
            "ensembl93_gene_biotype": biotype,
            "cellranger_3_0_0_filter_status": "excluded_by_documented_biotype_filter",
            "current_ensembl_lookup_present": "TRUE",
            "current_ensembl_display_name": now["current_lookup_display_name"],
            "crosswalk_can_restore_matrix_count": "FALSE",
        })

    category_counts: dict[str, int] = {}
    for row in missing_rows:
        key = str(row["ensembl93_gene_biotype"])
        category_counts[key] = category_counts.get(key, 0) + 1
    biotype_rows = [
        {
            "contract_id": "mv08e_missing_biotype_summary_v1",
            "ensembl93_gene_biotype": key,
            "genes": category_counts[key],
            "cellranger_3_0_0_allowed": "FALSE",
        }
        for key in sorted(category_counts)
    ]

    version_matches = sum(
        by_id[stable_id]["gene_version"] == version
        for _, stable_id, version in parsed_panel
    )
    name_matches = sum(
        by_id[stable_id]["gene_name"] == row["gene"]
        for row, stable_id, _ in parsed_panel
    )
    fingerprint_rows = [{
        "contract_id": "mv08e_reference_fingerprint_v1",
        "reference_name": "Cell Ranger 3.0.0 GRCh38",
        "source_annotation": "Ensembl release 93 GRCh38 GTF",
        "source_gtf_sha256": EXPECTED_GTF_SHA256,
        "unfiltered_ensembl93_genes": len(genes),
        "documented_allowed_biotypes": len(ALLOWED_BIOTYPES),
        "filtered_reference_genes": len(filtered),
        "hca_h5_feature_genes": len(h5_ids),
        "filtered_gtf_id_axis_sha256": axis_sha256(gtf_ids),
        "hca_h5_id_axis_sha256": axis_sha256(h5_ids),
        "id_axis_byte_exact": "TRUE",
        "name_axis_byte_exact": "TRUE",
        "panel_genes": 500,
        "panel_in_unfiltered_ensembl93": 500,
        "panel_in_filtered_cellranger_reference": 475,
        "common475_axis_sha256": common475_sha256,
        "panel_versions_matching_ensembl93": version_matches,
        "panel_names_matching_ensembl93": name_matches,
        "panel_excluded_by_reference_filter": 25,
    }]

    with args.hca_sample_metadata.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    one_donor_hits = [hit for hit in metadata["hits"] if donor_count(hit) == 1]
    manifest = read_csv(args.hca_file_manifest)
    manifest.sort(key=lambda row: int(row["file_order"]))
    if len(one_donor_hits) != 8 or [row["unit_id"] for row in manifest] != EXPECTED_UNITS:
        raise ValueError("HCA one-donor axis differs from MV8-C")
    hits_by_h5_bytes: dict[int, dict[str, object]] = {}
    for hit in one_donor_hits:
        h5_summary = file_summary(hit, "h5")
        size = int(h5_summary["totalSize"])
        if size in hits_by_h5_bytes:
            raise ValueError("one-donor H5 size is not unique for metadata join")
        hits_by_h5_bytes[size] = hit

    fastq_rows: list[dict[str, object]] = []
    for row in manifest:
        h5_bytes = int(row["file_size_bytes"])
        hit = hits_by_h5_bytes.get(h5_bytes)
        if hit is None:
            raise ValueError(f"no one-donor metadata hit for {row['unit_id']}")
        fastq = file_summary(hit, "fastq.gz")
        fastq_rows.append({
            "contract_id": "mv08e_hca_raw_read_resource_v1",
            "unit_order": int(row["file_order"]),
            "unit_id": row["unit_id"],
            "hca_sample_entry_id": hit["entryId"],
            "h5_file_uuid": row["file_uuid"],
            "h5_bytes": h5_bytes,
            "fastq_files": int(fastq["count"]),
            "fastq_bytes": int(fastq["totalSize"]),
            "fastq_gib": f"{int(fastq['totalSize']) / (1024 ** 3):.6f}",
            "raw_reads_downloaded": "FALSE",
        })
    if (
        sum(int(row["fastq_files"]) for row in fastq_rows) != EXPECTED_FASTQ_FILES
        or sum(int(row["fastq_bytes"]) for row in fastq_rows) != EXPECTED_FASTQ_BYTES
    ):
        raise ValueError("HCA raw-read totals differ from frozen metadata")
    fastq_summary_rows = [{
        "contract_id": "mv08e_hca_raw_read_resource_summary_v1",
        "biological_units": 8,
        "fastq_files": EXPECTED_FASTQ_FILES,
        "fastq_bytes": EXPECTED_FASTQ_BYTES,
        "fastq_gib": f"{EXPECTED_FASTQ_BYTES / (1024 ** 3):.6f}",
        "download_status": "not_downloaded",
        "reprocessing_status": "deferred_until_500_vs_475_sensitivity_decision",
    }]

    source_rows = [
        {
            "contract_id": "mv08e_authoritative_reference_source_v1",
            "source_order": 1,
            "source_id": "ensembl_release_93_gtf",
            "url": "https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz",
            "evidence_role": "unfiltered historical gene annotation",
            "accessed_date": "2026-08-16",
        },
        {
            "contract_id": "mv08e_authoritative_reference_source_v1",
            "source_order": 2,
            "source_id": "cellranger_3_0_0_reference_build_steps",
            "url": "https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps",
            "evidence_role": "Ensembl-93 provenance and exact allowed-biotype policy",
            "accessed_date": "2026-08-16",
        },
        {
            "contract_id": "mv08e_authoritative_reference_source_v1",
            "source_order": 3,
            "source_id": "hca_bone_marrow_publication",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6296228/",
            "evidence_role": "Cell Ranger and standard 10x GRCh38 processing provenance",
            "accessed_date": "2026-08-16",
        },
        {
            "contract_id": "mv08e_authoritative_reference_source_v1",
            "source_order": 4,
            "source_id": "ensembl_rest_lookup",
            "url": "https://rest.ensembl.org/documentation/info/lookup",
            "evidence_role": "current stable-ID existence check; not count recovery",
            "accessed_date": "2026-08-16",
        },
    ]
    decision_rows = [{
        "contract_id": "mv08e_reference_reconciliation_decision_v1",
        "decision": "exact_reference_identified_crosswalk_cannot_restore_counts",
        "hca_reference": "Cell Ranger 3.0.0 GRCh38 Ensembl93 filtered",
        "common_panel_genes": 475,
        "filtered_panel_genes": 25,
        "missing_due_identifier_retirement": "FALSE",
        "missing_due_documented_reference_filter": "TRUE",
        "symbol_or_stable_id_crosswalk_authorized": "FALSE",
        "zero_imputation_authorized": "FALSE",
        "next_authorized_stage": "reference_only_500_vs_475_sensitivity_prefreeze",
        "hca_fastq_download_authorized_now": "FALSE",
        "raw_read_reprocessing_preference_recorded": "TRUE",
        "raw_read_decision_boundary": "after_reference_sensitivity_and_custom_reference_prefreeze",
        "outcome_label_state": "closed",
        "biological_outcomes_computed": "FALSE",
    }]

    outputs = {
        "mv08e-reference-fingerprint.csv": fingerprint_rows,
        "mv08e-common475-panel.csv": common475_rows,
        "mv08e-missing-biotypes.csv": missing_rows,
        "mv08e-missing-biotype-summary.csv": biotype_rows,
        "mv08e-hca-fastq-resource.csv": fastq_rows,
        "mv08e-hca-fastq-resource-summary.csv": fastq_summary_rows,
        "mv08e-authoritative-sources.csv": source_rows,
        "mv08e-decision.csv": decision_rows,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, rows in outputs.items():
        write_csv(args.output_dir / name, rows)
    manifest_rows = [
        {
            "contract_id": "mv08e_reference_reconciliation_artifact_manifest_v1",
            "file": name,
            "bytes": (args.output_dir / name).stat().st_size,
            "sha256": sha256_file(args.output_dir / name),
            "contains_expression": "FALSE",
            "contains_cell_barcode": "FALSE",
            "contains_local_absolute_path": "FALSE",
        }
        for name in sorted(outputs)
    ]
    write_csv(args.output_dir / "mv08e-artifact-manifest.csv", manifest_rows)
    print("MV8-E exact reference reconciliation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
