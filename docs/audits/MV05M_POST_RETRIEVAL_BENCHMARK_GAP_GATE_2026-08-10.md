# MV5-M post-retrieval benchmark-gap and gate audit

## Outcome

MV5-M is complete. A deterministic, no-outcome gate selected an MV5-N
label-closed clustering contract and complete-matrix resource sprint as the
highest-value next action. It did not authorize clustering outcomes or full
distance production.

The selection used six frozen criteria and did not consult MV5-L tissue-level
heterogeneity. All seven production artifacts were independently reconstructed
and repeated byte for byte.

## What is now complete

Biological-conservation retrieval is complete for the existing-data LOSO
benchmark:

- MV5-E: SCT H0/H1 did not outperform matched SCT energy;
- MV5-K: integrated H0/H1 did not outperform matched integrated energy; and
- MV5-L: integration did not favorably change either topology-minus-energy
  contrast.

These null/negative results are retained. No future sensitivity or clustering
result may retroactively replace them.

## Remaining-axis audit

| Axis | Current disposition | Core reason |
|---|---|---|
| Technical mixing | Blocked pending identifiability design | Study is nested within tissue; held-out query has no same-study training comparator |
| Label-free clustering | Selected for contract/resource gate | High value; diagrams/query distances/stable-k helper exist; training matrices missing |
| Robustness | Eligible second | Scientifically strong, but grid is unfrozen and post-outcome sensitivity cannot rescue primary result |
| Integration expansion | Deferred | Avoid broad method search before one complete benchmark path |
| Gene topology | Blocked | SCT incomplete and integrated gene view structurally unavailable in pilot folds |
| Cell/gene fusion | Blocked | Both components have not independently passed |
| External validation/new data | Deferred | Important but existing-data-first rule still holds |
| Optimization/Rust | Deferred | Current measured gates pass; no fresh bottleneck trigger |

## Scored decision

| Candidate | Score | Validity | Disposition |
|---|---:|---|---|
| Clustering contract/resource gate | 45 | Pass | Selected MV5-N |
| Retrieval robustness | 43 | Pass | Next eligible |
| External validation/new data | 40 | Pass | Deferred existing-data-first |
| Technical mixing | 38 | **Fail** | Identifiability-blocked |
| Integration-method expansion | 37 | Pass | Deferred |
| Optimization/Rust | 35 | Pass | Deferred |
| Gene topology | 30 | Pass | Readiness-blocked |
| Cell/gene fusion | 18 | **Fail** | Component-blocked |

The weights were scientific value 3, reviewer relevance 2,
identifiability/validity 3, artifact readiness 2, resource feasibility 1, and
outcome-selection safety 2. Scores used a frozen 0-to-4 scale.

## Why technical mixing is not next

Technical mixing is central to modern integration benchmarking and remains
necessary for the eventual paper. It is not being skipped. It is blocked
because every currently eligible study belongs to one tissue, and the accepted
LOSO query-to-training design deliberately removes the query's study. A
same-study versus cross-study held-out contrast is therefore absent. An
in-sample training-only or cell-neighborhood metric could be built, but its
estimand and overcorrection safeguards must be specified separately before it
can carry scientific weight.

## Why clustering is next

Sample-level clustering is central to the dissertation and preprint, is a
major generalized reviewer workstream, and now has a corrected same-unit path.
The valid inductive task is not to cluster all samples transductively. It is to:

1. build clusters from each outer fold's training samples only;
2. choose k without labels from five-seed partition stability;
3. freeze training medoids; and
4. assign held-out samples using immutable query-to-training distances before
   opening biological or technical labels.

PAM/k-medoids is primary because it accepts arbitrary validated distances and
has a natural held-out assignment. Average linkage is a sensitivity only after
an out-of-sample rule is frozen. Existing spectral code remains ineligible.

## Exact missing workload

The existing artifacts contain 35,350 query-to-training pairs per
representation. Across the same 15 folds and five seeds, training-only PAM
requires 262,675 additional training-training pairs per representation, or
525,350 H0/H1 component rows per representation. Held-out query-query pairs are
not required.

Landscape-only linear lower bounds are 8.655 worker-hours for SCT and 4.280 for
integrated. They exclude baselines, validation, I/O, repeats, and clustering,
so full production is not authorized. MV5-N must profile bounded fold strata
and set honest caps first.

## Independent validation

The independent validator reconstructed:

- all six criteria and weights;
- all eight candidate totals and the unique score-45 selection;
- all nine axis dispositions;
- 7,070 existing query-training and 52,535 missing training pairs per seed;
- the 262,675/525,350 all-seed pair/component totals;
- the confidential-source boundary; and
- zero biological outcomes and zero downstream jobs.

All seven validation categories passed, and a clean rebuild reproduced all
seven production artifacts byte for byte. The focused test file passed 24
expectations. The complete package test suite passed with only the two existing
CRAN-gated skips. A fresh prefiltered source tarball built successfully, and
`R CMD check --no-manual` completed with `Status: OK` on R 4.4.1 under Ubuntu
22.04.4.

## MV5-N boundary

MV5-N may freeze and technically validate the training-only PAM/held-out-medoid
contract, k=2:10 five-seed one-SE stability, exact pair requests, cache/resume,
and bounded real-data resource admission. It must stop before full matrix
production, label opening, clustering ARI/NMI, spectral promotion, method or
tissue selection, gene/fusion work, new data, or optimization.

## Public-safety boundary

The audit used generalized reviewer workstreams only. No confidential wording
or reviewer source file is included. The dissertation and preprint PDFs,
reviewer material, private caches, and `example_run.r` remain untracked.
Nothing was pushed.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | Approve MV5-N contract/resource gate only |
| Correctness demonstrated? | MV5-M prioritization passes independent and byte-repeat validation |
| Computation feasible? | Unknown for full clustering; bounded admission required first |
| Biological interpretation permitted? | Prohibited; MV5-M computed no outcome |
| Next action | Execute MV5-N label-closed clustering contract and complete-matrix resource gate |
