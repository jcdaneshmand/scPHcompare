# MV5-D5 label-closed SCT cell retrieval-input completion audit

## Outcome

MV5-D5 is complete. The accepted 75 fold-seed bundles contain five immutable,
label-closed query-to-training retrieval methods for every one of the 35,350
frozen biological pairs. All 176,750 ranked rows independently validate and no
method failed. Tissue, approach, biological outcomes, clustering, integration,
gene topology, and fusion remained closed.

This audit authorizes only the next prediction-locked label-opening stage. It
does not report biological performance or select a method.

## Inputs and environment

| Input | Accepted identity |
|---|---|
| Cell fold coordinates | 75 MV5-D1 caches; 6,750 typed 384-cell/30-PC views |
| SCT source values | 450 MV5-D0 runtime-complete matrix caches |
| Topological components | MV5-D4 35,350 query-to-training rows with separate H0/H1 |
| Seeds | `20260805:20260809` |
| Runtime | Ubuntu WSL, R 4.4.1, locked project library |
| Workers | At most two |
| Group guards | 600 seconds and 4 GiB process-tree RSS |
| Stage guard | 7,200 seconds |

No dependency was installed or changed.

## Scientific contract

The production panel is:

1. raw exact all-active-level H0 landscape L2;
2. raw exact all-active-level H1 landscape L2;
3. raw `sqrt(H0^2 + H1^2)` as descriptive secondary only;
4. energy distance on the same cells and training-fit shared coordinates as
   the matched cell baseline; and
5. training-standardized same-panel pseudobulk Euclidean distance as context.

H0 and H1 are separate confirmatory methods. Their within-method rankings are
invariant to positive scalar multiplication. The raw combination is not
promoted to a primary method.

## Component-scale audit

The accepted MV5-D4 table contains zero within-training topology pairs. Using
its held-out-query rows to fit a scale would cross the frozen fit boundary.
Computing the missing within-training scope would add 262,675 biological pairs
and 525,350 H0/H1 rows, 7.431 times the completed MV5-D4 workload.

MV5-D5 therefore records, rather than hides, the following dispositions:

- H0: no scale required for its own ranking;
- H1: no scale required for its own ranking; and
- H0/H1 composite: no training scale fitted; raw Euclidean composite remains
  descriptive secondary.

No held-out query contributes to a fitted distance scale, and no additional
topology job was launched.

## Staging and admission

All 450 accepted SCT files were rehashed and read once. Their named row means
were reduced into five seed-specific private bundles. Each bundle binds all 90
normalization cache keys, its source-manifest hash, and the implementation hash.
The staging pass took 410 seconds on the Windows-mounted volume.

The representative `SRA550660` / `20260805` group contains five queries, 85
training references, 425 biological pairs, and 2,125 ranked rows. Two clean
computations produced byte-identical RDS files. Independent admission checks
found:

- exact H0, H1, and raw-composite agreement with MV5-D4;
- seven energy-oracle checks within tolerance;
- seven independently reconstructed pseudobulk checks within tolerance;
- canonical rankings and tie fields; and
- no outcome fields.

An initial successful admission emitted a shutdown warning from an explicitly
opened gzip connection. The reader was corrected to let `read.csv()` own the
compressed-file lifecycle, the implementation identity changed, and admission
was rerun into fresh private paths. The earlier private result was not
published or overwritten.

## Production completion

| Measure | Accepted value |
|---|---:|
| Fold-seed groups | 75 / 75 |
| Fixed seeds | 5 |
| Biological query-to-training pairs | 35,350 |
| Methods | 5 |
| Ranked rows | 176,750 |
| Completed method groups | 375 / 375 |
| Failed method groups | 0 |
| Aggregate worker elapsed | 1,461.448 seconds |
| Aggregate scientific operation time | 561.532 seconds |
| Median group elapsed | 17.985 seconds |
| Maximum group elapsed | 35.368 seconds |
| Peak process-tree RSS | 280,870,912 bytes (267.9 MiB) |
| Private bundle/audit bytes | 161,447,925 |
| Exact distance-tie rows | 0 |

The monitored queue completed in 873.8 seconds wall time with two workers. No
time, memory, or stage guard fired.

## Independent production validation

Every bundle was reopened independently. The validator checked bundle and file
hashes, fold-cache identities, mean-profile identities, exact query/training
partitions, common method axes, canonical ranks, unique pair IDs, completion
records, and the label firewall.

All H0, H1, and raw-composite rows were compared with MV5-D4. Three spread-out
energy pairs and three spread-out pseudobulk pairs were independently
recomputed in every group, for 225 checks per baseline. The maximum difference
over all topological and baseline checks was `1.136868e-13`, below the frozen
`1e-12` validation tolerance.

The result contains 176,750 unique pair IDs, 75 unique bundle keys, 75 closed
partitions, 75 outcome-absent groups, and zero downstream job counters.

## Determinism and resume

- The two representative admission bundles are byte-identical.
- A completed-queue resume reopened and validated all 75 bundles without
  launching replacement work.
- All 75 private bundle files and 75 private audit files retained identical
  paths, SHA-256 hashes, and byte sizes after resume.
- A second independent public assembly produced byte-identical copies of all
  six public retrieval artifacts, including timestamp-free gzip output.

The first resume verifier result was a false negative caused by comparing CSV-
restored integer sizes with in-memory numeric sizes using strict R type
identity. Direct diagnosis showed all 150 paths, hashes, and values were equal.
The verifier now compares numeric size values and passes. No production file
changed.

## Public evidence

| Artifact | Purpose |
|---|---|
| `mv05d5-cell-retrieval-rankings-2026-08-08.csv.gz` | 176,750 immutable distance/rank rows |
| `mv05d5-method-completion-2026-08-08.csv` | Complete method/failure ledger |
| `mv05d5-group-bundle-index-2026-08-08.csv` | Public identities and hashes for 75 private bundles |
| `mv05d5-method-registry-2026-08-08.csv` | Frozen five-method roles and policies |
| `mv05d5-component-scale-disposition-2026-08-08.csv` | Explicit scale decision for each topology group/method |
| `mv05d5-public-assembly-summary-2026-08-08.csv` | Counts, public hashes, and stop counters |
| `mv05d5-mean-profile-staging-2026-08-08.csv` | Five seed-bundle source identities |
| `mv05d5-production-group-resources-2026-08-08.csv` | Per-group completion and resource evidence |
| `mv05d5-independent-group-validation-2026-08-08.csv` | Per-group independent checks |
| `mv05d5-independent-validation-summary-2026-08-08.csv` | Full validation summary |
| `mv05d5-admission-validation-2026-08-08.csv` | Admission oracle and replay evidence |
| `mv05d5-resume-validation-2026-08-08.csv` | 150-file immutable resume proof |
| `mv05d5-public-repeat-validation-2026-08-08.csv` | Six-file public assembly repeat proof |

Private fold caches, SCT caches, mean-profile bundles, production RDS bundles,
logs, PDFs, reviewer material, and `example_run.r` are not tracked.

## Package verification

The correctly package-loaded test suite passed 578 expectations with zero
failures, warnings, or skips. The first source-package check completed with
zero errors and warnings but one namespace NOTE for undeclared `ave()`. The
precise `stats::ave` import was added to `NAMESPACE`; this did not change the
production helper or its bundle-bound implementation hash
`887b38f222946bfb3e2db31f738fbcd6c6ba4c4c9f304faa8fdcd6ddd7dd91f9`.

The initial build also exposed an operational inefficiency: `R CMD build`
copied the 37-GiB ignored private cache tree before applying `.Rbuildignore`.
It completed safely, but subsequent checks should use a pre-filtered source
staging area. The clean-staged rerun completed in 301.6 seconds with
`Status: OK` and zero errors, warnings, or notes.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | Approve SCT cell retrieval inputs; retain separate H0/H1 and descriptive raw composite |
| Correctness demonstrated? | Pass |
| Computation feasible? | Yes; all groups completed far below guards |
| Biological interpretation permitted? | Prohibited in MV5-D5; no labels were opened |
| Next action | Advance to prediction-locked MV5-E retrieval endpoint evaluation only |

## Next-stage boundary

MV5-E must verify the public ranking hash and 375/375 completion records before
opening labels. It may calculate only prespecified held-out retrieval endpoints
and blocked uncertainty from the immutable rankings. It may not refit,
rescale, rerank, select methods, cluster, integrate, construct gene topology,
fuse views, or substitute new data.
