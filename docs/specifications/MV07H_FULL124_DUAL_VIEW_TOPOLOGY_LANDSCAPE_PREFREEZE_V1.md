# MV7-H full-124 dual-view topology and landscape prefreeze v1

Date: 2026-08-15

MV7-H constructs the complete label-closed descriptive topology over the exact
124 retained samples, five deterministic cell-subsampling seeds, and the
MV7-FP 500-gene panel. It does not change the panel, selected cells, global
standardization, global PCA, topology orientation, filtration, or persistence-
landscape definition accepted by MV7-E through MV7-G.

## Scientific axes

For each seed, MV7-H reuses the immutable MV7-G center, scale, 30-component PCA
rotation, standardization identity, and input-cache identities. It does not
refit those transformations. It reconstructs all 124 typed sample views and
requires the six sentinel view pairs to reproduce their MV7-G cache and payload
identities exactly.

- `cell_topology_v1`: 384 cells are points in the shared 30-PC Euclidean space.
- `gene_topology_v1`: 500 panel genes are points under explicit Pearson-
  correlation chord distance over the same 384 matched cells.

Every view receives complete Vietoris-Rips H0/H1 persistence with field 2 and
no filtration truncation. Essential H0 is excluded from landscapes. All finite
H1 intervals must have positive persistence, and every finite H0 death multiset
must reproduce the view's minimum-spanning-tree edge multiset.

The complete PH workload is five source bundles, 1,240 typed views, and 1,240
PH records: 124 samples by five seeds by two views. Source and PH artifacts are
written atomically and resumed only by exact hash, size, and identity checks.

## Landscape definition and grouping

Landscapes use finite positive-persistence intervals, every consecutive active
level, zero-padding where a paired diagram has less depth, and exact squared-L2
integration. H0 and H1 are calculated and stored separately. There is no fixed
grid, no fixed level cap, and no dense landscape-function matrix. The accepted
Rust kernel is execution-only and must match its frozen SHA-256; the R reference
engine remains the canonical numerical oracle.

Within each seed and view, all 7,626 unordered sample pairs are evaluated. The
landscape workload is split into 20 independent groups (five seeds by two views
by H0/H1), each containing exactly 7,626 component rows. No cross-seed pair is
computed. The group with the greatest MV7-G sentinel interval burden runs first
as the prospective stress gate. Only that group is authorized until resource,
repeat, R-oracle, independent, and immutable-resume gates pass; the remaining
19 groups require a separate stage-two authorization.

## Resource and validation policy

Execution is serial with one worker and zero retries. Source groups have
1,800-second/8-GiB caps; cell-PH records 600-second/4-GiB caps; gene-PH records
1,800-second/8-GiB caps; and landscape groups 3,600-second/12-GiB caps. The
complete sprint has a 48-worker-hour, 12-GiB peak, and 4-GiB retained-private-
state boundary. A cap breach stops without retry or scientific exclusion.

One full seed source bundle and the twelve already frozen sentinel PH records
for that seed must repeat byte-for-byte. Full PH validation must independently
reopen all 1,240 records, reconstruct axis counts, check MST/H1 validity, and
perform a prespecified new-sample Ripserr/GUDHI corroboration. The stress
landscape group must repeat byte-for-byte and pass minimum-, quantile-, and
maximum-depth R exact/adaptive oracle checks before stage two.

## Firewall and stop boundary

Production inputs are limited to frozen sample IDs, seeds, cache identities,
selected-cell identities, feature/panel identities, and accepted transform
identities. Tissue, approach, study, class, labels, outcomes, and previous
benchmark results are prohibited from execution and stopping decisions.
Clustering, distance combination, label access, outcome analysis, method
selection, and manuscript claims remain zero or closed throughout MV7-H.
