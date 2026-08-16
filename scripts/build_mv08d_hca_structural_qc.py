#!/usr/bin/env python3
"""Build deterministic public structural/QC evidence for MV8-D.

Only aggregate structural, QC, mapping, and digest evidence is emitted. Raw H5
payloads, expression values, and cell barcodes remain in the private cache.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import re
import sys

import h5py
import numpy as np
import scipy
from scipy import sparse


PANEL_SHA256 = "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e"
EXPECTED_TOTAL_BYTES = 202_770_089
EXPECTED_UNITS = [f"HCA_BM_{index:03d}" for index in range(1, 9)]
ENSEMBL_AT_END = re.compile(r"(ENSG\d+)(?:\.\d+)?$")
FINAL_VERSION = re.compile(r"\.\d+$")


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def axis_sha256(values: list[str] | np.ndarray) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(str(value).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def decode_array(dataset: h5py.Dataset) -> np.ndarray:
    raw = dataset[()]
    decoded = np.array([
        value.decode("utf-8") if isinstance(value, (bytes, np.bytes_)) else str(value)
        for value in raw
    ], dtype=object)
    if np.any(decoded == ""):
        raise ValueError(f"empty identifier in {dataset.name}")
    return decoded


def stable_ensembl(value: str) -> str:
    match = ENSEMBL_AT_END.search(value)
    return match.group(1) if match else ""


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, extrasaction="raise",
            quoting=csv.QUOTE_ALL, lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    rows.sort(key=lambda row: int(row["file_order"]))
    if len(rows) != 8 or [row["unit_id"] for row in rows] != EXPECTED_UNITS:
        raise ValueError("manifest biological-unit axis differs from MV8-C")
    if [int(row["file_order"]) for row in rows] != list(range(1, 9)):
        raise ValueError("manifest order differs from MV8-C")
    if sum(int(row["file_size_bytes"]) for row in rows) != EXPECTED_TOTAL_BYTES:
        raise ValueError("manifest total bytes differ from MV8-C")
    return rows


def read_panel(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    rows.sort(key=lambda row: int(row["panel_order"]))
    if len(rows) != 500 or [int(row["panel_order"]) for row in rows] != list(range(1, 501)):
        raise ValueError("ordered panel axis is not exactly 1 through 500")
    if {row["panel_sha256"] for row in rows} != {PANEL_SHA256}:
        raise ValueError("panel identity differs from MV7-FP/MV7-H")
    stable_ids = [stable_ensembl(row["feature_id"]) for row in rows]
    if any(not value for value in stable_ids) or len(set(stable_ids)) != 500:
        raise ValueError("panel Ensembl stable identifiers are incomplete or duplicated")
    for row, stable_id in zip(rows, stable_ids):
        row["ensembl_stable_id"] = stable_id
    return rows


def finite_stat(values: np.ndarray, function: str) -> str:
    if values.size == 0:
        return ""
    if function == "min":
        result = float(np.min(values))
    elif function == "median":
        result = float(np.median(values))
    elif function == "max":
        result = float(np.max(values))
    else:
        raise ValueError(function)
    return f"{result:.6f}"


def inspect_unit(
    manifest_row: dict[str, str], panel: list[dict[str, str]], cache_dir: Path
) -> tuple[dict[str, object], dict[str, object], dict[str, object], list[dict[str, object]]]:
    path = cache_dir / manifest_row["unit_id"] / manifest_row["file_name"]
    expected_bytes = int(manifest_row["file_size_bytes"])
    expected_sha = manifest_row["file_sha256"].lower()
    if not path.is_file():
        raise FileNotFoundError(f"verified H5 is absent for {manifest_row['unit_id']}")
    observed_bytes = path.stat().st_size
    observed_sha = sha256_file(path)
    if observed_bytes != expected_bytes or observed_sha != expected_sha:
        raise ValueError(f"byte/checksum identity failed for {manifest_row['unit_id']}")

    verification = {
        "contract_id": "mv08d_hca_exact_file_verification_v1",
        "file_order": int(manifest_row["file_order"]),
        "unit_id": manifest_row["unit_id"],
        "file_name": manifest_row["file_name"],
        "file_uuid": manifest_row["file_uuid"],
        "file_version": manifest_row["file_version"],
        "expected_bytes": expected_bytes,
        "observed_bytes": observed_bytes,
        "expected_sha256": expected_sha,
        "observed_sha256": observed_sha,
        "verified": "TRUE",
    }

    with h5py.File(path, "r") as h5:
        required = [
            "matrix/barcodes", "matrix/data", "matrix/indices",
            "matrix/indptr", "matrix/shape", "matrix/features/id",
            "matrix/features/name", "matrix/features/feature_type",
        ]
        absent = [name for name in required if name not in h5]
        if absent:
            raise ValueError(f"missing required H5 datasets: {', '.join(absent)}")
        shape_raw = np.asarray(h5["matrix/shape"][()]).astype(np.int64)
        if shape_raw.shape != (2,) or np.any(shape_raw <= 0):
            raise ValueError("H5 matrix shape is not two positive dimensions")
        n_features, n_barcodes = (int(shape_raw[0]), int(shape_raw[1]))
        barcodes = decode_array(h5["matrix/barcodes"])
        feature_ids = decode_array(h5["matrix/features/id"])
        feature_names = decode_array(h5["matrix/features/name"])
        feature_types = decode_array(h5["matrix/features/feature_type"])
        data = np.asarray(h5["matrix/data"][()])
        indices = np.asarray(h5["matrix/indices"][()]).astype(np.int64, copy=False)
        indptr = np.asarray(h5["matrix/indptr"][()]).astype(np.int64, copy=False)
        data_dtype = str(data.dtype)
        index_dtype = str(h5["matrix/indices"].dtype)
        pointer_dtype = str(h5["matrix/indptr"].dtype)
        chemistry_description = str(h5.attrs.get("chemistry_description", ""))
        h5_filetype = str(h5.attrs.get("filetype", ""))
        h5_format_version = str(h5.attrs.get("version", ""))
        genomes = decode_array(h5["matrix/features/genome"]) if (
            "matrix/features/genome" in h5
        ) else np.array([], dtype=object)

    if len(barcodes) != n_barcodes or len(feature_ids) != n_features:
        raise ValueError("H5 axis lengths differ from shape")
    if len(feature_names) != n_features or len(feature_types) != n_features:
        raise ValueError("H5 feature annotation lengths differ from shape")
    if len(indptr) != n_barcodes + 1 or indptr[0] != 0:
        raise ValueError("H5 CSC pointer axis is invalid")
    if np.any(np.diff(indptr) < 0) or indptr[-1] != len(data) or len(indices) != len(data):
        raise ValueError("H5 CSC pointer/data lengths are inconsistent")
    if len(data) == 0 or np.min(indices) < 0 or np.max(indices) >= n_features:
        raise ValueError("H5 sparse indices are empty or out of bounds")
    if not np.all(np.isfinite(data)) or np.any(data < 0) or not np.all(data == np.floor(data)):
        raise ValueError("H5 counts are not finite nonnegative integers")
    if len(set(barcodes.tolist())) != n_barcodes:
        raise ValueError("H5 barcode axis is duplicated")

    ge_mask = feature_types == "Gene Expression"
    ge_original = np.flatnonzero(ge_mask)
    if ge_original.size == 0:
        raise ValueError("H5 has no exact Gene Expression feature rows")
    ge_stable = np.array([stable_ensembl(value) for value in feature_ids[ge_mask]], dtype=object)
    nonempty_stable = ge_stable[ge_stable != ""]
    if len(set(nonempty_stable.tolist())) != len(nonempty_stable):
        raise ValueError("Gene Expression Ensembl stable identifiers are duplicated")

    matrix = sparse.csc_matrix((data, indices, indptr), shape=(n_features, n_barcodes))
    matrix.eliminate_zeros()
    gene_expression = matrix[ge_mask, :].tocsc()
    entry_cells = np.asarray(gene_expression.getnnz(axis=0)).ravel() >= 200
    if not np.any(entry_cells):
        raise ValueError("no cell passes the frozen 200-feature entry rule")
    entry_feature_detection = np.asarray(
        gene_expression[:, entry_cells].getnnz(axis=1)
    ).ravel()
    entry_features = entry_feature_detection >= 3
    entry_matrix = gene_expression[entry_features, :][:, entry_cells].tocsc()
    entry_names = feature_names[ge_mask][entry_features]
    nfeature = np.asarray(entry_matrix.getnnz(axis=0)).ravel()
    total_counts = np.asarray(entry_matrix.sum(axis=0)).ravel().astype(np.float64)
    if np.any(total_counts <= 0):
        raise ValueError("entry-filtered cells include a nonpositive total count")
    mito_mask = np.char.startswith(entry_names.astype(str), "MT-")
    ribo_names = entry_names.astype(str)
    ribo_mask = np.char.startswith(ribo_names, "RPS") | np.char.startswith(ribo_names, "RPL")
    mito_counts = np.asarray(entry_matrix[mito_mask, :].sum(axis=0)).ravel()
    ribo_counts = np.asarray(entry_matrix[ribo_mask, :].sum(axis=0)).ravel()
    percent_mito = 100.0 * mito_counts / total_counts
    percent_ribo = 100.0 * ribo_counts / total_counts
    explicit_cells = (
        (nfeature >= 500) & (nfeature <= 9000) &
        (percent_mito <= 20.0) & (percent_ribo > 5.0)
    )
    post_matrix = entry_matrix[:, explicit_cells].tocsc()
    final_entry_features = np.asarray(post_matrix.getnnz(axis=1)).ravel() > 3

    entry_original = ge_original[entry_features]
    final_original = entry_original[final_entry_features]
    entry_original_set = set(entry_original.tolist())
    final_original_set = set(final_original.tolist())
    # Frozen mapping is stable-Ensembl-only: no symbol-only rescue is allowed.
    stable_to_original = {
        stable_id: int(original)
        for stable_id, original in zip(ge_stable.tolist(), ge_original.tolist())
        if stable_id
    }
    symbol_to_originals: dict[str, list[int]] = {}
    for original in ge_original.tolist():
        symbol_to_originals.setdefault(str(feature_names[original]), []).append(int(original))

    mapping_rows: list[dict[str, object]] = []
    mapped = 0
    final_mapped = 0
    symbol_matches = 0
    for panel_row in panel:
        stable_id = panel_row["ensembl_stable_id"]
        original = stable_to_original.get(stable_id)
        mapping_status = "missing"
        h5_id = ""
        h5_name = ""
        symbol_match = False
        entry_retained = False
        final_retained = False
        symbol_candidates: list[int] = []
        if original is not None:
            mapped += 1
            h5_id = str(feature_ids[original])
            h5_name = str(feature_names[original])
            symbol_match = h5_name == panel_row["gene"]
            symbol_matches += int(symbol_match)
            entry_retained = original in entry_original_set
            final_retained = original in final_original_set
            final_mapped += int(final_retained)
            mapping_status = "mapped_final_retained" if final_retained else (
                "mapped_entry_retained_final_failed" if entry_retained
                else "mapped_entry_failed"
            )
        else:
            # Diagnostic only; these candidates never rescue the stable-ID gate.
            symbol_candidates = symbol_to_originals.get(panel_row["gene"], [])
        mapping_rows.append({
            "contract_id": "mv08d_hca_panel_mapping_v1",
            "unit_id": manifest_row["unit_id"],
            "panel_order": int(panel_row["panel_order"]),
            "panel_sha256": PANEL_SHA256,
            "reference_feature_id": panel_row["feature_id"],
            "reference_gene": panel_row["gene"],
            "ensembl_stable_id": stable_id,
            "h5_feature_id": h5_id,
            "h5_feature_name": h5_name,
            "symbol_match": str(symbol_match).upper(),
            "entry_feature_retained": str(entry_retained).upper(),
            "final_feature_retained": str(final_retained).upper(),
            "symbol_only_candidate_count": len(symbol_candidates),
            "symbol_only_candidate_ids": "|".join(
                str(feature_ids[index]) for index in symbol_candidates
            ),
            "symbol_only_rescue_applied": "FALSE",
            "mapping_status": mapping_status,
        })

    post_barcodes = barcodes[entry_cells][explicit_cells]
    post_count = int(explicit_cells.sum())
    structure = {
        "contract_id": "mv08d_hca_h5_structure_v1",
        "file_order": int(manifest_row["file_order"]),
        "unit_id": manifest_row["unit_id"],
        "matrix_features": n_features,
        "matrix_barcodes": n_barcodes,
        "gene_expression_features": int(ge_mask.sum()),
        "matrix_nnz": int(matrix.nnz),
        "data_dtype": data_dtype,
        "indices_dtype": index_dtype,
        "indptr_dtype": pointer_dtype,
        "feature_type_count": len(set(feature_types.tolist())),
        "genome": "|".join(sorted(set(genomes.tolist()))),
        "chemistry_description": chemistry_description,
        "h5_filetype": h5_filetype,
        "h5_format_version": h5_format_version,
        "barcode_axis_sha256": axis_sha256(barcodes),
        "feature_id_axis_sha256": axis_sha256(feature_ids),
        "schema_pass": "TRUE",
        "count_type_pass": "TRUE",
        "unique_barcode_pass": "TRUE",
        "unique_ensembl_pass": "TRUE",
    }
    qc = {
        "contract_id": "mv08d_hca_legacy_qc_depth_v1",
        "file_order": int(manifest_row["file_order"]),
        "unit_id": manifest_row["unit_id"],
        "raw_barcodes": n_barcodes,
        "entry_cells_min_200_features": int(entry_cells.sum()),
        "entry_features_min_3_cells": int(entry_features.sum()),
        "post_qc_cells": post_count,
        "final_features_gt_3_cells": int(final_entry_features.sum()),
        "panel_features_mapped": mapped,
        "panel_features_symbol_match": symbol_matches,
        "panel_features_final_retained": final_mapped,
        "post_qc_barcode_axis_sha256": axis_sha256(post_barcodes),
        "post_qc_nfeature_min": finite_stat(nfeature[explicit_cells], "min"),
        "post_qc_nfeature_median": finite_stat(nfeature[explicit_cells], "median"),
        "post_qc_nfeature_max": finite_stat(nfeature[explicit_cells], "max"),
        "post_qc_count_min": finite_stat(total_counts[explicit_cells], "min"),
        "post_qc_count_median": finite_stat(total_counts[explicit_cells], "median"),
        "post_qc_count_max": finite_stat(total_counts[explicit_cells], "max"),
        "post_qc_percent_mito_median": finite_stat(percent_mito[explicit_cells], "median"),
        "post_qc_percent_ribo_median": finite_stat(percent_ribo[explicit_cells], "median"),
        "depth_384_pass": str(post_count >= 384).upper(),
        "ordered_panel_500_mapped_pass": str(mapped == 500).upper(),
        "ordered_panel_500_final_retained_pass": str(final_mapped == 500).upper(),
    }
    return verification, structure, qc, mapping_rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    manifest = read_manifest(args.manifest)
    panel = read_panel(args.panel)
    verification_rows: list[dict[str, object]] = []
    structure_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    mapping_rows: list[dict[str, object]] = []
    for row in manifest:
        print(f"MV8-D structural/QC audit {row['file_order']}/8: {row['unit_id']}", flush=True)
        verification, structure, qc, mappings = inspect_unit(row, panel, args.cache_dir)
        verification_rows.append(verification)
        structure_rows.append(structure)
        qc_rows.append(qc)
        mapping_rows.extend(mappings)

    write_csv(args.output_dir / "mv08d-file-verification.csv", verification_rows,
        list(verification_rows[0]))
    write_csv(args.output_dir / "mv08d-h5-structure.csv", structure_rows,
        list(structure_rows[0]))
    write_csv(args.output_dir / "mv08d-qc-depth.csv", qc_rows, list(qc_rows[0]))
    write_csv(args.output_dir / "mv08d-panel-mapping.csv", mapping_rows,
        list(mapping_rows[0]))
    software_rows = [{
        "contract_id": "mv08d_hca_software_v1",
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "h5py_version": h5py.__version__,
        "matrix_orientation": "features_by_barcodes_csc",
        "qc_contract": "mv08d_hca_legacy_qc_v1",
        "outcome_label_state": "closed",
        "biological_outcomes_computed": "FALSE",
    }]
    write_csv(args.output_dir / "mv08d-software.csv", software_rows,
        list(software_rows[0]))

    all_pass = (
        len(verification_rows) == 8 and
        all(row["schema_pass"] == "TRUE" and row["count_type_pass"] == "TRUE"
            for row in structure_rows) and
        all(row["depth_384_pass"] == "TRUE" and
            row["ordered_panel_500_mapped_pass"] == "TRUE" and
            row["ordered_panel_500_final_retained_pass"] == "TRUE"
            for row in qc_rows)
    )
    next_gate = (
        "immutable_transform_runtime_validation" if all_pass
        else "panel_annotation_reconciliation_or_cohort_block"
    )
    gate_rows = [{
        "contract_id": "mv08d_hca_structural_qc_gate_v1",
        "decision": "structural_qc_pass" if all_pass else "structural_qc_block",
        "exact_files_verified": len(verification_rows),
        "structural_files_passed": sum(row["schema_pass"] == "TRUE" for row in structure_rows),
        "units_depth_384_passed": sum(row["depth_384_pass"] == "TRUE" for row in qc_rows),
        "units_panel_500_mapped_passed": sum(row["ordered_panel_500_mapped_pass"] == "TRUE" for row in qc_rows),
        "units_panel_500_final_retained_passed": sum(row["ordered_panel_500_final_retained_pass"] == "TRUE" for row in qc_rows),
        "pca_coordinates_computed": "FALSE",
        "ph_computed": "FALSE",
        "landscapes_computed": "FALSE",
        "outcomes_computed": "FALSE",
        "next_gate": next_gate,
    }]
    write_csv(args.output_dir / "mv08d-structural-qc-gate.csv", gate_rows,
        list(gate_rows[0]))
    print(f"MV8-D structural/QC gate: {gate_rows[0]['decision']}", flush=True)
    return 0 if all_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
