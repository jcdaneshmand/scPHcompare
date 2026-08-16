# MV8-D HCA structural/QC prefreeze v1

## Document control

| Field | Frozen value |
|---|---|
| Status | Prospectively frozen before any HCA H5 expression payload is opened |
| Date | 2026-08-16 |
| Parent decision | MV8-C `metadata_closed_download_candidate` |
| Owner authorization | Exact eight-file structural/QC download authorized; later sprints may continue only through already-prespecified gates |
| Project | HCA `cc95ff89-2e68-4a08-a234-480eca21ce79` |
| Primary units | Exactly `HCA_BM_001` through `HCA_BM_008` from the published MV8-C manifest |
| Outcomes | Closed; no biological label, retrieval, clustering, or manuscript endpoint calculation |
| Current scope | Download, byte/checksum verification, H5 structure, legacy-comparable QC, ordered-panel mapping/depth, and immutable-transform interface admission |
| Explicitly excluded | PCA coordinates, PH, landscapes, topology distances, clustering, endpoint evaluation, claims, donor substitution, or reference refit |

## 1. Exact inputs and fail-closed acquisition

The only admissible payloads are the eight rows of
`docs/audits/mv08c-hca-admission-evidence/mv08c-primary-file-manifest.csv`, in
`file_order`. File UUID, version, basename, byte count, and SHA-256 are jointly
binding. Downloads are sequential and published to a private ignored cache only
after byte count and SHA-256 pass. A verified existing file is a cache hit. A
partial, changed, extra, substituted, or checksum-incompatible file is never
opened and stops the sprint.

The downloader may follow the HCA stable locator's redirect, but it must not log
or publish any transient signed redirect URL. Public evidence contains the
stable file UUID/version and cryptographic receipt, never local absolute paths.

## 2. H5 structural contract

Each verified file must expose the 10x sparse CSC count structure under
`/matrix`: `barcodes`, `data`, `indices`, `indptr`, `shape`, and
`features/{id,name,feature_type}`. The following are hard gates:

- `shape` is exactly two positive dimensions interpreted as features by barcodes;
- `indptr` has `barcodes + 1` entries, begins at zero, is monotone, and ends at
  the common `data`/`indices` length;
- feature and barcode axis lengths equal `shape`;
- sparse row indices are in bounds;
- counts are finite, nonnegative, and integer-valued;
- feature IDs, feature names, feature types, and barcodes decode without empty
  identifiers;
- the analysis uses only rows whose type is exactly `Gene Expression`; and
- duplicated Gene Expression Ensembl stable IDs or duplicated barcodes stop
  admission rather than being repaired after observation.

These checks establish matrix semantics; the filename alone does not.

## 3. Frozen legacy-comparable QC

MV8-D reproduces the actual historical `R/PH_Calculation.R` behavior rather
than choosing thresholds after viewing the HCA results. In the locked accepted
environment this corresponds to Seurat 5.3.0 / SeuratObject 5.0.2 and the
following ordered operations on Gene Expression counts:

1. retain barcodes with at least 200 detected features (`min.features = 200`);
2. on those barcodes, retain features detected in at least three cells
   (`min.cells = 3`);
3. recompute per-cell `nFeature_RNA`, total counts, mitochondrial count fraction
   from symbols matching `^MT-`, and ribosomal count fraction from symbols
   matching `^RP[SL]`;
4. retain cells with `500 <= nFeature_RNA <= 9000`, mitochondrial percentage
   `<= 20`, and ribosomal percentage `> 5`; and
5. retain features detected in more than three post-QC cells.

The 200-feature entry rule is not an EmptyDrops-style statistical droplet caller;
it is a deterministic legacy heuristic. MV8-D labels it accordingly and
does not claim that it is an independently optimized modern QC procedure. A
modern cell-calling/QC sensitivity may be designed later, prospectively and
separately; it may not replace this primary comparability gate after outcomes
are opened.

Every biological unit must retain at least 384 cells after step 4. Failure is a
technical cohort-level block: no donor is replaced, and the eight-unit set is
not narrowed after seeing results.

## 4. Ordered 500-gene mapping

The immutable panel is the ordered 500 rows of
`docs/audits/mv07h-prefreeze-evidence-v4/mv07h-panel.csv`, with panel digest
`48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e`.

Mapping is identifier-led and fixed before H5 inspection:

1. parse the panel Ensembl identifier from `feature_id` and remove only a final
   numeric version suffix;
2. remove only a final numeric version suffix from H5 feature IDs;
3. require exactly one Gene Expression H5 row for every panel Ensembl stable ID;
4. record symbol agreement, but use the stable-ID match when an annotation-era
   symbol differs; and
5. prohibit symbol-only rescue, zero filling, duplicate aggregation, panel
   replacement, or order changes.

All 500 mapped features must also survive the final post-QC `>3`-cell gene gate
in every donor. Any missing, ambiguous, or post-QC-ineligible panel member
blocks both corrected views for this cohort.

## 5. Frozen cells, normalization, and projection boundary

MV8-D may prove that deterministic selection is feasible but does not run
SCTransform or calculate coordinates. The next admitted calculation must:

- sort post-QC barcode IDs and select exactly 384 without replacement under
  seeds `20260805` through `20260809`, using the existing
  `select_matched_cells()` contract;
- run per-donor, per-seed `SCTransform(variable.features.n = 3000,
  return.only.var.genes = FALSE)` on the selected raw counts, with the accepted
  locked runtime and one computational thread;
- extract the exact ordered panel without substitution;
- apply the corresponding accepted seed-specific 124-sample reference center
  and scale; and
- project cell coordinates by multiplying the transposed standardized
  `500 x 384` query matrix by the immutable `500 x 30` reference rotation.

The same standardized query matrix supplies the gene view. The reference panel,
fit samples, center, scale, PCA rotation, and source topology remain immutable.
No HCA value contributes to a reference fit. The five accepted source-bundle
SHA-256 values are frozen in the MV8-D reference-transform ledger.

## 6. Admission decision and next authorization

MV8-D passes only if all eight exact files pass every checksum, structural,
count, identifier, panel, post-QC depth, reproducibility, and immutable-interface
gate. Production and clean-repeat public ledgers must be byte-identical, and an
independent validator must pass before the decision is published.

A pass authorizes only a separately auditable label-closed normalization and
immutable-projection sprint. It does not authorize PCA refitting, PH,
landscapes, distances, clustering, label opening, endpoint claims, or manuscript
promotion. Those remain later gates. If eventually reached, landscapes remain
all finite positive intervals, essential H0 excluded, every consecutive active
level, H0 and H1 separate, exact or error-controlled squared-L2, grid-free as
the primary definition, and streamed/chunked.

## 7. Public/private boundary

The eight H5 payloads, partial downloads, cell barcodes, selected-cell lists,
and expression values remain private ignored artifacts. Public evidence may
contain official file identity/checksum receipts, aggregate dimensions, QC
counts and distribution summaries, panel mapping/retention status, hashes of
ordered private axes, reference-transform identities, software versions,
validation results, and the admission decision. No local absolute path,
transient signed URL, expression matrix, or barcode is published.
