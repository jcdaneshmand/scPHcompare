# MV6-C global-core matched SCT feasibility prefreeze v1

| Field | Frozen value |
|---|---|
| Contract | `mv06c_global_core_matched_sct_feasibility_v1` |
| Date | 2026-08-14 |
| Accepted base | `ff6cb80` |
| Owner decision | Pursue MV6-B option A before abandoning full fusion |
| Existing-data axis | 90 samples × five seeds = 450 accepted SCT caches |
| Cells per cache | Exact accepted 384 selected cells |
| Candidate panel | One fixed 500-gene global core |
| Harmonization scope | All existing samples/seeds; label-closed but explicitly transductive |
| Later fold transforms | Training-only standardization and PCA remain mandatory |
| Biological labels/outcomes | Prohibited |
| Authorized computation | Eligibility, ranking, panel stability, and bounded resource inventory only |
| Required stop | Before source reconstruction, PCA, PH, landscapes, fusion, clustering, or outcomes |

## Purpose

Determine whether the existing 90-sample/five-seed SCT corpus contains a
single scientifically valid 500-gene panel that can support matched cell and
gene topology for every accepted sample-seed instance. This is the first
bounded feasibility step for the project-owner-selected global-core direction.
It does not authorize a rerun or establish that fusion improves an endpoint.

The global panel is intentionally a different estimand from the MV5-D1
training-derived panels. It uses feature availability and within-sample
variance from the complete existing-data corpus, but no tissue, study outcome,
approach, clustering, retrieval result, or manuscript endpoint. Any later
claim must disclose this as transductive technical harmonization rather than
future-study-inductive feature discovery.

## Frozen input identity

Only the accepted MV5-D0 v2 SCT matrix caches may be read:

| Input | Frozen value |
|---|---|
| Resource ledger | `mv05d0-sct-cache-resources-v2-2026-08-07.csv` |
| Ledger SHA-256 | `73f757a91c202a8e38dfa746fc8816f8e272caf912c33b4b959ef55102c68308` |
| Expected records | 450 |
| Expected biological samples | 90 |
| Expected seeds | `20260805:20260809` |
| Payload | `mv05d0_sct_data_matrix_v1` |
| Expected shape | named genes/features × 384 selected cells |
| Outcome state | closed |

Every private file SHA-256, cache key, payload identity, selected-cell identity,
shape, seed, and sample axis must agree with the accepted ledger. Missing,
extra, duplicated, unreadable, or stale records abort the sprint.

## Frozen global-core rule

The rule reuses the accepted MV-03/MV5 feature semantics:

1. Form the exact feature-ID intersection across all 450 matrices.
2. Canonicalize feature IDs with `canonical_mv03_gene_ids()`.
3. Classify canonical genes with `technical_feature_rules_v1` and exclude
   mitochondrial `^MT-`, ribosomal-protein `^RP[SL]`, and hemoglobin
   `^HB(?!P)` categories from the primary panel.
4. For every common feature in every cache, calculate the sample variance over
   the exact accepted 384 cells with `.mv03_row_variance()`.
5. A feature is eligible only when every one of its 450 variances is finite and
   strictly greater than `.Machine$double.eps`.
6. Within each cache, rank common features by descending variance using
   `rank(-variance, ties.method = "min", na.last = "keep")`.
7. Calculate the median of the 450 within-cache ranks for each feature.
8. Order eligible retained candidates by median rank, canonical gene, then
   feature ID using radix order; retain only the first feature for a duplicated
   canonical gene.
9. Select the first 500 rows as the frozen candidate global-core panel.

There is no result-dependent alternative threshold, panel size, category rule,
imputation, query-specific active set, or sample/view exclusion. If fewer than
500 unique canonical retained candidates pass, the decision is
`stop_global_core_insufficient`.

## Frozen diagnostics

The feasibility output must report:

1. verified record/sample/seed counts and accepted ledger/file hashes;
2. union and exact-intersection feature counts;
3. counts excluded for technical category, nonfinite variance, nonpositive
   variance, or duplicate canonical identity;
4. eligible unique canonical gene count and margin above/below 500;
5. the exact ordered 500-gene candidate panel, its SHA-256, feature category,
   median rank, and label-free rank summaries;
6. for each seed, the panel eligible under that seed alone, its top-500 overlap
   with the global panel, Jaccard index, and rank correlation;
7. confirmation that every selected panel feature is present, finite, and
   nonconstant in all 450 caches;
8. elapsed time, peak process-tree RSS when measurable, bytes read, and output
   size for the inventory itself;
9. lower-bound future axes for matched cell/gene views, H0/H1 diagram
   components, directed landscape pairs, training-fitted scales, and the frozen
   five-weight fusion grid; and
10. a deterministic artifact manifest and decision record.

The compact ordered panel may be public because it is required method
configuration. Per-sample expression values, variances, ranks, cell IDs,
matrices, coordinates, and private paths remain prohibited from public output.

## Independent validation and repeat

An independent validator must not call the production panel selector. It must:

- rehash the accepted ledger and all output artifacts;
- reconstruct the record/sample/seed axes and all aggregate eligibility counts;
- independently verify the exact ordered panel from sufficient aggregate
  candidate evidence;
- confirm all 500 selected features pass the all-cache presence/variance flags;
- reconstruct seed-overlap and workload summaries;
- reject biological-label or private-value columns; and
- confirm all prohibited execution counters remain zero.

A clean run must reproduce every production output byte-for-byte.

## Decision rule and next boundary

The decision is `go_bounded_matched_sct_profile` only when:

- all 450 accepted cache identities verify;
- at least 500 unique canonical retained candidates pass the frozen rule;
- the selected 500 all pass every-cache presence, finite, and variance checks;
- the exact panel and all diagnostics independently reconstruct;
- the clean repeat is byte-identical; and
- the inventory remains within 8 GiB RSS and 1,800 seconds elapsed.

A `go` authorizes only a separately prefrozen bounded MV6-D profile that
reconstructs matched global-core SCT cell/gene sources for prespecified
fold-seed sentinels and measures cell/gene PH and landscape cost. It does not
authorize all 6,750 views, blocked outcomes, clustering, or fusion evaluation.

## Landscape, integration, and scope boundary

Any later profile must retain the accepted dissertation-aligned landscape
definition: separate H0/H1, finite positive-persistence intervals, all active
consecutive levels, and exact or error-controlled L2 integration without a
universal level or grid cap.

Integrated gene topology remains undefined and outside MV6-C. The candidate
direction is matched SCT cell/gene topology first; accepted integrated cell
results remain a separate cell-view analysis unless a later contract defines a
fold-safe integrated gene-expression representation.

MV6-C prohibits PCA, PH, landscapes, distances, fusion rows, labels, endpoints,
clusters, learned weights, new data, Rust adoption, package defaults,
manuscript claims, PDFs, reviewer material, release, DOI, binaries, and
`example_run.r`.
