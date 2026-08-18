#!/usr/bin/env python3
"""Build the public, label-closed MV8-H count-sentinel prefreeze.

This builder inspects only immutable file identities and resource metadata. It
does not invoke Cell Ranger, open a count matrix, or inspect donor attributes,
cell barcodes, QC values, labels, or biological outcomes.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path


EXPECTED_UNIT_MANIFEST_SHA256 = (
    "765788d28087a16ce0e6e73fb4a532285809b96f187a08ed09621c7ebe996ef2"
)
EXPECTED_FASTQ_MANIFEST_SHA256 = (
    "b315e097ef9b2f8cf1cad31cda08c4a440919641f9fff1674cf09e86928b2136"
)
EXPECTED_RUNTIME_TREE_SHA256 = (
    "aafd39e293e0ba9d14dba3896a6aeda077304531a2702d26bda0c62c4688fdf3"
)
EXPECTED_LAUNCHER_SHA256 = (
    "4ee3a1670b4f14c826004fe8e17b4759e1edc701b15ff2e9623753bf1b34d4d6"
)
EXPECTED_REFERENCE_TREE_SHA256 = (
    "5e2aff9e7154e6b02f98552a4419bd48edce66e617e579ae562e714f79199f1c"
)
EXPECTED_REFERENCE_FILES = 19
EXPECTED_REFERENCE_BYTES = 20_765_871_518
EXPECTED_TOTAL_FASTQ_BYTES = 85_034_239_918
SENTINEL_UNIT = "HCA_BM_002"
SENTINEL_UNIT_ORDER = 2
SENTINEL_FASTQ_BYTES = 11_249_623_632
SENTINEL_SAMPLE = "MantonBM2_HiSeq_9"
RUN_ID = "mv08h_count_sentinel_hca_bm_002"
LOCAL_CORES = 4
LOCAL_MEMORY_GIB = 32
ABSOLUTE_RSS_CAP_BYTES = 80 * 1024**3
WORKSPACE_CAP_BYTES = 200 * 1024**3
ELAPSED_CAP_SECONDS = 96 * 60 * 60
FREE_SPACE_FLOOR_BYTES = 1024**4


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing empty evidence: {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), extrasaction="raise",
            quoting=csv.QUOTE_ALL, lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows({
            key: "TRUE" if value is True else "FALSE" if value is False else value
            for key, value in row.items()
        } for row in rows)
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


def tree_identity(root: Path) -> tuple[int, int, int, str]:
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
    regular = symlinks = regular_bytes = 0
    for kind, relative, path in entries:
        if kind == "F":
            size = path.stat().st_size
            file_sha = sha256_file(path)
            regular += 1
            regular_bytes += size
            record = f"F\0{relative}\0{size}\0{file_sha}\n"
        else:
            target = os.readlink(path)
            symlinks += 1
            record = f"L\0{relative}\0{target}\n"
        digest.update(record.encode("utf-8"))
    return regular, symlinks, regular_bytes, digest.hexdigest()


def public_command() -> str:
    return (
        f"cellranger count --id={RUN_ID} "
        "--transcriptome=<verified_custom_reference> "
        "--fastqs=<verified_hca_bm_002_fastq_directory> "
        f"--sample={SENTINEL_SAMPLE} --chemistry=SC3Pv2 "
        "--expect-cells=7000 --include-introns=false --create-bam=false "
        f"--nosecondary --localcores={LOCAL_CORES} "
        f"--localmem={LOCAL_MEMORY_GIB} --disable-ui"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--unit-manifest", type=Path, required=True)
    parser.add_argument("--fastq-manifest", type=Path, required=True)
    parser.add_argument("--runtime-identity", type=Path, required=True)
    parser.add_argument("--reference-identity", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if sha256_file(args.unit_manifest) != EXPECTED_UNIT_MANIFEST_SHA256:
        raise RuntimeError("unit manifest identity differs from the admitted freeze")
    if sha256_file(args.fastq_manifest) != EXPECTED_FASTQ_MANIFEST_SHA256:
        raise RuntimeError("FASTQ manifest identity differs from the admitted freeze")
    units = read_csv(args.unit_manifest)
    fastqs = read_csv(args.fastq_manifest)
    runtime = read_csv(args.runtime_identity)
    reference = read_csv(args.reference_identity)
    if len(units) != 8 or len(fastqs) != 48 or len(runtime) != 1 or len(reference) != 1:
        raise RuntimeError("an admitted input has unexpected cardinality")
    if sum(int(row["fastq_bytes"]) for row in units) != EXPECTED_TOTAL_FASTQ_BYTES:
        raise RuntimeError("unit manifest total differs from the admitted freeze")
    if runtime[0]["tree_sha256"] != EXPECTED_RUNTIME_TREE_SHA256:
        raise RuntimeError("runtime tree identity differs")
    if runtime[0]["launcher_sha256"] != EXPECTED_LAUNCHER_SHA256:
        raise RuntimeError("Cell Ranger launcher identity differs")
    if reference[0]["tree_sha256"] != EXPECTED_REFERENCE_TREE_SHA256:
        raise RuntimeError("reference closure identity differs")

    ranked = sorted(
        units, key=lambda row: (-int(row["fastq_bytes"]), int(row["unit_order"]))
    )
    chosen = ranked[0]
    if not all([
        chosen["unit_id"] == SENTINEL_UNIT,
        int(chosen["unit_order"]) == SENTINEL_UNIT_ORDER,
        int(chosen["fastq_bytes"]) == SENTINEL_FASTQ_BYTES,
        int(chosen["fastq_files"]) == 6,
        chosen["expected_cells"] == "7000",
    ]):
        raise RuntimeError("deterministic maximum-burden sentinel selection changed")

    selection_rows: list[dict[str, object]] = []
    for rank, row in enumerate(ranked, 1):
        selection_rows.append({
            "contract_id": "mv08h_count_sentinel_selection_v1",
            "burden_rank": rank,
            "unit_order": int(row["unit_order"]),
            "unit_id": row["unit_id"],
            "fastq_files": int(row["fastq_files"]),
            "fastq_bytes": int(row["fastq_bytes"]),
            "selection_rule": "maximum_complete_unit_fastq_bytes_then_unit_order",
            "selected": row["unit_id"] == SENTINEL_UNIT,
            "label_fields_consulted": False,
            "expression_or_qc_consulted": False,
        })

    selected_fastqs = sorted(
        (row for row in fastqs if row["unit_id"] == SENTINEL_UNIT),
        key=lambda row: int(row["file_order"]),
    )
    if len(selected_fastqs) != 6:
        raise RuntimeError("sentinel FASTQ cardinality differs")
    fastq_rows: list[dict[str, object]] = []
    for row in selected_fastqs:
        path = args.cache_dir / SENTINEL_UNIT / row["file_name"]
        observed_bytes = path.stat().st_size if path.is_file() else -1
        observed_sha = sha256_file(path) if path.is_file() else ""
        passed = (
            observed_bytes == int(row["file_size_bytes"])
            and observed_sha == row["file_sha256"]
        )
        if not passed:
            raise RuntimeError(f"sentinel FASTQ identity failed: order {row['file_order']}")
        fastq_rows.append({
            "contract_id": "mv08h_count_sentinel_fastq_v1",
            "file_order": int(row["file_order"]),
            "unit_id": SENTINEL_UNIT,
            "lane": int(row["lane"]),
            "read_role": row["read_role"],
            "file_name": row["file_name"],
            "file_uuid": row["file_uuid"],
            "expected_bytes": int(row["file_size_bytes"]),
            "expected_sha256": row["file_sha256"],
            "cache_identity_passed": passed,
            "private_path_published": False,
        })
    if sum(int(row["expected_bytes"]) for row in fastq_rows) != SENTINEL_FASTQ_BYTES:
        raise RuntimeError("sentinel FASTQ byte total differs")

    regular, symlinks, regular_bytes, reference_sha = tree_identity(args.reference_dir)
    reference_exact = all([
        regular == EXPECTED_REFERENCE_FILES,
        symlinks == 0,
        regular_bytes == EXPECTED_REFERENCE_BYTES,
        reference_sha == EXPECTED_REFERENCE_TREE_SHA256,
    ])
    if not reference_exact:
        raise RuntimeError("live reference tree no longer matches its closure")
    reference_rows = [{
        "contract_id": "mv08h_count_sentinel_reference_binding_v1",
        "reference_name": "GRCh38_Ensembl93_target_complete_33563",
        "regular_files": regular,
        "symlinks": symlinks,
        "regular_file_bytes": regular_bytes,
        "tree_sha256": reference_sha,
        "runtime_tree_sha256": EXPECTED_RUNTIME_TREE_SHA256,
        "launcher_sha256": EXPECTED_LAUNCHER_SHA256,
        "binding_passed": reference_exact,
    }]

    command_rows = [
        {"contract_id": "mv08h_count_sentinel_command_v1", "parameter_order": 1,
         "parameter": "id", "frozen_value": RUN_ID, "scientific_role": "unique private run directory"},
        {"contract_id": "mv08h_count_sentinel_command_v1", "parameter_order": 2,
         "parameter": "transcriptome", "frozen_value": "verified_custom_reference_tree_sha256=" + EXPECTED_REFERENCE_TREE_SHA256, "scientific_role": "exact Ensembl93 target-complete reference"},
        {"contract_id": "mv08h_count_sentinel_command_v1", "parameter_order": 3,
         "parameter": "fastqs", "frozen_value": "verified_complete_HCA_BM_002_six_file_directory", "scientific_role": "one preserved biological unit"},
        {"contract_id": "mv08h_count_sentinel_command_v1", "parameter_order": 4,
         "parameter": "sample", "frozen_value": SENTINEL_SAMPLE, "scientific_role": "exact FASTQ prefix"},
        {"contract_id": "mv08h_count_sentinel_command_v1", "parameter_order": 5,
         "parameter": "chemistry", "frozen_value": "SC3Pv2", "scientific_role": "official 10x 3-prime-v2 metadata"},
        {"contract_id": "mv08h_count_sentinel_command_v1", "parameter_order": 6,
         "parameter": "expect-cells", "frozen_value": 7000, "scientific_role": "official expected-cell target"},
        {"contract_id": "mv08h_count_sentinel_command_v1", "parameter_order": 7,
         "parameter": "include-introns", "frozen_value": "false", "scientific_role": "historical exon-only comparison target"},
        {"contract_id": "mv08h_count_sentinel_command_v1", "parameter_order": 8,
         "parameter": "create-bam", "frozen_value": "false", "scientific_role": "omit unneeded BAM without changing feature-barcode counts"},
        {"contract_id": "mv08h_count_sentinel_command_v1", "parameter_order": 9,
         "parameter": "nosecondary", "frozen_value": "true", "scientific_role": "defer Cell Ranger secondary analysis"},
        {"contract_id": "mv08h_count_sentinel_command_v1", "parameter_order": 10,
         "parameter": "localcores", "frozen_value": LOCAL_CORES, "scientific_role": "downward-only compute allocation"},
        {"contract_id": "mv08h_count_sentinel_command_v1", "parameter_order": 11,
         "parameter": "localmem", "frozen_value": LOCAL_MEMORY_GIB, "scientific_role": "downward-only GiB allocation"},
        {"contract_id": "mv08h_count_sentinel_command_v1", "parameter_order": 12,
         "parameter": "disable-ui", "frozen_value": "true", "scientific_role": "noninteractive local execution"},
        {"contract_id": "mv08h_count_sentinel_command_v1", "parameter_order": 13,
         "parameter": "public_command", "frozen_value": public_command(), "scientific_role": "exact prospective command without private paths"},
    ]

    resource_rows = [
        {"contract_id": "mv08h_count_sentinel_resource_policy_v1", "resource_order": 1,
         "resource": "local_cores", "selected_value": LOCAL_CORES, "unit": "cores", "gate": "exact", "rationale": "quarter original 16-core maximum; longer runtime accepted"},
        {"contract_id": "mv08h_count_sentinel_resource_policy_v1", "resource_order": 2,
         "resource": "local_memory", "selected_value": LOCAL_MEMORY_GIB, "unit": "GiB", "gate": "exact", "rationale": "half original 64-GiB allocation"},
        {"contract_id": "mv08h_count_sentinel_resource_policy_v1", "resource_order": 3,
         "resource": "process_tree_rss_absolute_cap", "selected_value": ABSOLUTE_RSS_CAP_BYTES, "unit": "bytes", "gate": "invalidate_above", "rationale": "retain admitted 80-GiB hard scientific-acceptance ceiling"},
        {"contract_id": "mv08h_count_sentinel_resource_policy_v1", "resource_order": 4,
         "resource": "run_workspace_cap", "selected_value": WORKSPACE_CAP_BYTES, "unit": "bytes", "gate": "invalidate_above", "rationale": "retain admitted 200-GiB per-unit ceiling"},
        {"contract_id": "mv08h_count_sentinel_resource_policy_v1", "resource_order": 5,
         "resource": "minimum_free_space", "selected_value": FREE_SPACE_FLOOR_BYTES, "unit": "bytes", "gate": "must_remain_at_or_above", "rationale": "retain 1-TiB reserve"},
        {"contract_id": "mv08h_count_sentinel_resource_policy_v1", "resource_order": 6,
         "resource": "elapsed_observation_cap", "selected_value": ELAPSED_CAP_SECONDS, "unit": "seconds", "gate": "invalidate_above", "rationale": "96 hours accommodates quarter-core execution without scientific retry"},
        {"contract_id": "mv08h_count_sentinel_resource_policy_v1", "resource_order": 7,
         "resource": "concurrency", "selected_value": 1, "unit": "unit", "gate": "exact", "rationale": "serialized sentinel only"},
        {"contract_id": "mv08h_count_sentinel_resource_policy_v1", "resource_order": 8,
         "resource": "automatic_termination", "selected_value": "false", "unit": "boolean", "gate": "non_destructive", "rationale": "breach invalidates admission and requires owner review; launcher never kills"},
        {"contract_id": "mv08h_count_sentinel_resource_policy_v1", "resource_order": 9,
         "resource": "automatic_deletion", "selected_value": "false", "unit": "boolean", "gate": "non_destructive", "rationale": "preserve all partial and complete artifacts for diagnosis"},
    ]

    validation_rows = [
        (1, "prelaunch_runtime", "Cell Ranger version, full-tree closure, and launcher SHA-256 exact"),
        (2, "prelaunch_reference", "19-file reference tree bytes and SHA-256 exact"),
        (3, "prelaunch_fastqs", "six selected FASTQs independently size/SHA-256 exact; zero partials"),
        (4, "prelaunch_scope", "run target absent; only one process; command and four-core/32-GiB allocation exact"),
        (5, "execution_completion", "exit zero, Cell Ranger pipestance success marker, and nonempty outs directory"),
        (6, "resource_closure", "RSS/workspace/elapsed caps and 1-TiB floor pass from complete monitor receipt"),
        (7, "output_structure", "raw and filtered feature-barcode HDF5 plus molecule-info outputs present; BAM absent"),
        (8, "feature_axis", "33,563 unique reference genes; ordered exact500 and common475 targets present exactly once"),
        (9, "matrix_integrity", "HDF5 sparse dimensions/indptr/indices/data consistent; barcodes unique and private"),
        (10, "qc_firewall", "no QC threshold, cell selection, normalization, PCA, PH, persistence landscape, cluster, label, or outcome calculation"),
        (11, "public_firewall", "publish identities and aggregates only; no paths, expression values, barcodes, donor attributes, labels, or outcomes"),
        (12, "authorization_stop", "successful sentinel opens only separate structural/QC review and remaining-unit decision; nothing automatic"),
    ]
    validation = [{
        "contract_id": "mv08h_count_sentinel_validation_contract_v1",
        "check_order": order, "check_id": check_id, "frozen_requirement": requirement,
        "checked_before_execution": order <= 4,
        "count_execution_performed_by_prefreeze": False,
    } for order, check_id, requirement in validation_rows]

    firewall_rows = [
        (1, "FASTQ byte sizes and SHA-256", True, True),
        (2, "reference/runtime relative identities", True, True),
        (3, "aggregate resource samples and run status", True, True),
        (4, "matrix dimensions and feature-axis identities", True, True),
        (5, "absolute private paths", False, False),
        (6, "expression or UMI values", False, False),
        (7, "cell barcodes", False, False),
        (8, "donor attributes or identifiers", False, False),
        (9, "QC values or eligibility decisions", False, False),
        (10, "study/tissue/approach labels", False, False),
        (11, "biological outcomes", False, False),
    ]
    firewall = [{
        "contract_id": "mv08h_count_sentinel_firewall_v1",
        "field_order": order, "field_class": field_class,
        "private_validation_permitted": private_permitted,
        "public_release_permitted": public_permitted,
    } for order, field_class, private_permitted, public_permitted in firewall_rows]

    gate_rows = [
        (1, "deterministic_label_closed_selection", "pass", "unique maximum FASTQ-byte unit; stable tie-break by unit order"),
        (2, "six_fastq_live_identity", "pass", "6/6 exact size and SHA-256; 11249623632 bytes"),
        (3, "reference_live_identity", "pass", "19 files; 20765871518 bytes; exact frozen tree SHA-256"),
        (4, "runtime_reference_contract", "pass", "Cell Ranger 8.0.1 and exact custom reference closures rebound"),
        (5, "command_and_resources", "pass", "SC3Pv2; exon-only; no BAM/secondary; 4 cores; 32 GiB; serialized"),
        (6, "firewall_and_stop", "pass", "count/QC/matrix/label/outcome execution remain closed"),
    ]
    gates = [{
        "contract_id": "mv08h_count_sentinel_gate_status_v1",
        "gate_order": order, "gate": gate, "status": status, "evidence": evidence,
    } for order, gate, status, evidence in gate_rows]
    decision = [{
        "contract_id": "mv08h_count_sentinel_decision_v1",
        "decision": "count_sentinel_prefreeze_exact_await_execution_authorization",
        "count_sentinel_prefreeze_completed": True,
        "count_sentinel_execution_authorized": False,
        "matrix_access_authorized": False,
        "qc_authorized": False,
        "remaining_units_authorized": False,
        "pca_ph_landscape_authorized": False,
        "label_access_authorized": False,
        "biological_outcomes_authorized": False,
        "deletion_authorized": False,
        "next_gate": "owner_or_existing_plan_authorization_for_exact_frozen_one_unit_count_execution",
    }]

    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    files: dict[str, list[dict[str, object]]] = {
        "mv08h-count-sentinel-selection.csv": selection_rows,
        "mv08h-count-sentinel-fastqs.csv": fastq_rows,
        "mv08h-count-sentinel-reference-binding.csv": reference_rows,
        "mv08h-count-sentinel-command.csv": command_rows,
        "mv08h-count-sentinel-resource-policy.csv": resource_rows,
        "mv08h-count-sentinel-validation-contract.csv": validation,
        "mv08h-count-sentinel-firewall.csv": firewall,
        "mv08h-count-sentinel-gate-status.csv": gates,
        "mv08h-count-sentinel-decision.csv": decision,
    }
    for name, rows in files.items():
        write_csv(output / name, rows)

    report = f"""# MV8-H Cell Ranger 8.0.1 count-sentinel prefreeze

## Outcome

One complete unit is prospectively frozen without running `cellranger count`.
The deterministic, label-closed selection rule chooses the complete unit with
the greatest total compressed FASTQ bytes, breaking any tie by original unit
order. The unique result is `{SENTINEL_UNIT}`: six exact files totaling
{SENTINEL_FASTQ_BYTES:,} bytes. All six live cache files pass independent size
and SHA-256 checks. Selection did not inspect expression, QC, donor attributes,
labels, or outcomes.

## Exact execution contract

```text
{public_command()}
```

The live 19-file custom reference independently retains its exact
{EXPECTED_REFERENCE_BYTES:,}-byte tree and SHA-256
`{EXPECTED_REFERENCE_TREE_SHA256}`. Cell Ranger 8.0.1 remains bound by runtime
tree SHA-256 `{EXPECTED_RUNTIME_TREE_SHA256}` and launcher SHA-256
`{EXPECTED_LAUNCHER_SHA256}`.

The fixed contract is **4 cores** and **32 GiB** (four-core/32-GiB), a
downward-only resource amendment from the earlier 16-core/64-GiB maximum.
The existing 80-GiB process-tree RSS ceiling, 200-GiB workspace ceiling, and
1-TiB free-space floor remain. The elapsed observation ceiling is 96 hours so the
deliberately smaller core allocation may take longer. Monitoring is
non-destructive: a breach invalidates scientific admission and preserves artifacts
for owner review; it does not kill or delete.

## Scientific and privacy stop

This prefreeze does **not** execute or authorize count. It does not open a
matrix, evaluate QC, select cells, normalize data, run PCA, compute PH or
persistence landscapes, cluster, access labels/outcomes, process the remaining
seven units, or delete anything. Public closure may contain only file/runtime/
reference identities, matrix shapes and feature-axis identities, aggregate
resource measurements, and run status—not private paths, expression values,
barcodes, donor attributes, labels, or outcomes.

The corrected downstream landscape definition is unchanged: cells and genes
are separate typed observation views; H0 and H1 remain separate; essential H0
is excluded; and landscapes use every consecutive active level with exact or
error-controlled squared-L2 integration, no fixed grid, and no universal level cap.
A successful future sentinel can open only a separate structural/QC review
and remaining-unit decision.
"""
    write_text(output / "MV08H_CELLRANGER8_COUNT_SENTINEL_PREFREEZE_2026-08-18.md", report)
    print(
        f"MV8-H count sentinel prefreeze built: {SENTINEL_UNIT}; "
        f"{len(fastq_rows)}/6 FASTQs; count_executed=false"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
