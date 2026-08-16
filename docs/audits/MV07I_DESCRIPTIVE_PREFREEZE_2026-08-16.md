# MV7-I descriptive topology and clustering prefreeze

Date: 2026-08-16

Status: complete; label-closed matrix and clustering production authorized

## Frozen representations and summaries

MV7-I freezes six representations on the exact 124-sample, five-seed axis:

- primary cell-H0 and cell-H1;
- secondary cell `sqrt(H0^2 + H1^2)`;
- primary gene-H0 and gene-H1; and
- secondary gene `sqrt(H0^2 + H1^2)`.

The composite cannot replace either component. H1 squared-distance contribution
is descriptive only. Combined distances and H1 fractions are calculated within
each seed before aggregation. Every pair reports all five seed values plus
median, minimum, maximum, IQR, and raw MAD (`constant = 1`). No seed may be
selected for favorable results.

## Frozen label-closed clustering

For each representation, fit PAM independently to all five seed matrices for
every `k = 2:10`, yielding 270 candidate fits. Select `k` from the ten pairwise
seed ARIs using the five-seed delete-one-seed jackknife and the smallest-
`k`-within-one-SE rule. If the complete grid or matched axes fail, report
`no_stable_k`. Average linkage is the sole algorithm sensitivity and uses the
PAM-selected `k`, for 30 selected fits. Cluster IDs use sorted member
signatures. Spectral clustering remains deferred pending an affinity/eigengap
prefreeze; Ward and direct distance-matrix k-means remain prohibited.

## Population and metadata boundaries

The full 124 samples are descriptive across eight tissues. The prior 90-sample,
five-tissue population remains the cross-study primary context. The 34 added
samples from three single-study tissues cannot support new cross-study tissue
claims. The three below-threshold samples remain outside this analysis.

Tissue and study come only from the frozen MV7-D reconciliation. Approach comes
only from the MV7-E canonical field; historical heuristic approach fields are
prohibited. These metadata may not be loaded by matrix reconstruction,
clustering, or `k` selection. A separate outcome prefreeze is required after
partitions are immutable.

## Validation and authorization

The builder passes 13/13 categories. An independent rebuild is byte-identical
for all 11 artifacts, and independent validation passes 13/13. The first
validator attempt failed closed because a recomputed SHA vector retained R path
names while the CSV values were unnamed; removing names from both sides fixed
the comparison without changing a hash or scientific contract. No failed
attempt output was published.

Execution is limited to one worker, 3,600 seconds, 4 GiB process-tree RSS,
atomic private artifacts, deterministic repeat, and immutable resume. The next
sprint may reconstruct 30 seed-specific matrices, derive six descriptive
aggregates, fit the frozen label-free PAM grid, select `k`, and fit average
linkage. Labels, outcomes, claims, cross-view fusion, new distances, and
external data remain closed.
