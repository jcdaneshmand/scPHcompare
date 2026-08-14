# MV6-B matched gene-view scale-up and resource gate v1

| Field | Frozen value |
|---|---|
| Contract | `mv06b_matched_gene_view_scaleup_resource_gate_v1` |
| Date | 2026-08-14 |
| Accepted base | `53c2a75` |
| Existing-data axis | 90 samples, 15 leave-one-study-out folds, five seeds |
| Representations | SCT and inductive integration audited separately |
| Views | Existing accepted cell views plus candidate `gene_topology_v1` |
| Gene geometry | Pearson correlation-chord distance over the exact 384 selected cells |
| Outcomes and biological labels | Prohibited |
| Authorized computation | Identity, structural eligibility, reconstructability, and resource inventory only |
| Required stop | Before gene-view construction, PH, landscapes, fusion, endpoints, or clustering |

## Purpose

Determine whether the accepted 90-sample MV-05 blocked benchmark can support a
complete matched gene-view axis without changing its scientific unit,
held-out-study design, feature identities, selected cells, or accepted cell
views. This gate must distinguish three questions:

1. whether every required sample-fold-seed instance has the full frozen gene
   panel with finite nonzero within-sample variance;
2. whether the accepted private caches retain enough expression information to
   reconstruct the candidate gene view without inventing an integrated gene
   space; and
3. whether a later bounded profiling sprint can estimate full H0/H1 diagram and
   all-active-level landscape cost safely.

MV6-B cannot establish biological complementarity, select a fusion weight,
exclude inconvenient views, or authorize a full calculation.

## Frozen evidence sources

The inventory may read only the accepted, label-closed MV5-D1 and MV5-G cache
families and their accepted resource ledgers:

| Evidence family | Required records | Accepted ledger SHA-256 |
|---|---:|---|
| MV5-D1 SCT fold caches | 75 | `205c78dd33d509627fea32517ff2c03326c04f88a0f4ec34e5034f12cd1aca71` |
| MV5-G integrated coordinate caches | 75 | `8ac47aaaaa0a9ee16dc749f0bb507c999c88c7d2346654951fc0447f3025a1dd` |

The ledgers bind 75 fold-seed groups, 6,750 cell-view instances per
representation, the exact D1 500-gene panel per group, 384 selected cells per
sample, training/query roles, cache identities, private result hashes, and
zero downstream counters. Source hashes must verify before structural fields
are summarized. Missing, stale, duplicated, or extra records stop the gate.

## Strict matched-view eligibility

The primary scale-up candidate inherits the corrected dual-view contract. For
each representation, fold, seed, and sample:

- the candidate gene view must use the same 384 selected cells and the exact
  ordered 500-gene D1 panel associated with the accepted cell view;
- every panel gene must be observed, finite, and nonconstant over those cells;
- the gene distance must be the explicit Pearson correlation-chord metric;
- all 90 sample instances in a fold-seed group must pass; and
- no absent or constant gene may be zero-filled, imputed, removed per sample,
  or replaced after seeing the held-out sample.

Training-mean zero filling remains valid only for the accepted cell-coordinate
projection. It is prospectively invalid for gene topology because it creates a
constant gene point. Query-specific active feature subsets are also
prospectively invalid for the strict candidate because they change point
identity and count across samples.

An instance with incomplete feature presence is recorded as structurally
ineligible. An instance whose variance has not been checked from an accepted
expression cache is recorded as `variance_unresolved`, never as eligible.

## Representation-specific reconstructability

### SCT

The audit may use MV5-D1 panels, centers, scales, missing-feature counts, and
accepted MV5-D0 identities to determine whether a strict standardized
gene-by-cell source is reconstructable. Training instances already admitted by
the D1 typed-source constructor may be credited only for the exact panel and
fit scope it validated. Held-out instances require complete feature presence;
variance remains unresolved unless the accepted matrix cache is directly
checked.

### Inductive integration

The accepted MV5-G artifact is a 384-cell by 30-PC cell-coordinate cache. Its
query `active_features` provenance may diagnose incomplete D1 panels, but cell
coordinates are not a gene-by-cell expression matrix and must never be
transposed or relabeled as gene topology.

The audit must determine whether an accepted, immutable, fold-safe corrected
gene-expression payload exists. If none exists, integrated gene topology is
`representation_undefined_from_accepted_artifacts`. Rerunning SCTransform or
creating a new corrected expression assay is a new scientific and resource
contract and is not authorized by this gate.

## Required public-safe outputs

The deterministic inventory must report:

1. ledger hashes, expected and observed cache counts, and cache-hash matches;
2. group/view totals by representation and training/query role;
3. full-panel, incomplete-panel, variance-unresolved, and undefined-expression
   counts without sample, cell, gene, tissue, or outcome identifiers;
4. affected group counts and the distribution of absent-feature counts;
5. whether the original 6,750-view cell axis can be matched without exclusions;
6. what expression reconstruction would be required and whether it is already
   accepted;
7. a lower-bound workload inventory for views, PH dimensions, directed
   query-to-training pairs, scales, and fusion components; and
8. a deterministic decision record and artifact manifest.

Group-level output may use already-public held-out study IDs and seeds, but it
must not expose private paths, feature names, cell IDs, matrices, coordinates,
tissue, approach, or biological endpoints.

## Decision rule

The decision is `profile_candidate_gene_ph` only if at least one
representation has all 6,750 matched instances structurally eligible, its
fold-safe expression source is already accepted and reconstructable, and no
sample/view exclusion or point-identity change is required. That decision
would authorize only a separately prefrozen, bounded, resource-capped gene-PH
profile.

Otherwise the decision is `stop_contract_revision_required`. The audit must
name each failed invariant and enumerate possible future contracts without
selecting one. It must not automatically:

- remove affected samples, studies, folds, seeds, or genes;
- form a panel from held-out availability;
- accept sample-specific point sets;
- reuse cell PCs as a gene space;
- treat SCT expression as an integrated gene representation;
- launch PH, landscapes, fusion, outcomes, or clustering; or
- proceed to MV-07 as if G-MV6 had passed.

## Landscape and fusion boundary

If a future candidate is authorized, its diagrams must retain the accepted
dissertation-aligned definition: finite positive-persistence intervals, H0 and
H1 separately, every active consecutive landscape level, and exact or
error-controlled L2 integration with no universal level or grid cap. Fusion
must refit component scales on training partitions only and preserve all four
cell/gene H0/H1 components.

This gate does not calculate a persistence diagram or landscape. It cannot
modify the accepted landscape definition.

## Prohibited scope

Biological labels and endpoints, result-informed exclusions, new data,
Similarity Network Fusion, learned weights, clustering, Rust adoption, package
defaults, manuscript claims, PDFs, reviewer correspondence, release, DOI,
binaries, and `example_run.r` remain outside this sprint.
