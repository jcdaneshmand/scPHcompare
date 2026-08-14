# MV5-H label-closed integrated cell-PH production audit

## Outcome

MV5-H is complete. All 6,750 corrected cells-as-observations integrated H0/H1
records are valid, deterministic, resumable, and within the frozen resource
caps. A separately specified integrated landscape-distance sprint is authorized.
No landscape, distance, retrieval, clustering, or biological outcome ran.

## Input and orientation

The complete public manifest binds 75 MV5-G groups across 15 LOSO studies and
five seeds. Each view uses the accepted ordered 384-by-30 integrated coordinate
matrix: rows are cells, columns are reference-fitted inductive integrated
coordinates, and Euclidean distance defines the point cloud. Integrated views
have their own identity and are never mislabeled as SCT or gene topology.

| Axis | Result |
|---|---:|
| Groups | 75 / 75 |
| Views | 6,750 / 6,750 |
| Held-out / training views | 450 / 6,300 |
| H0 intervals | 2,592,000 |
| H1 intervals | 1,545,943 |
| Failures | 0 |

Every view retains 384 H0 intervals: 383 finite MST merges and one explicit
essential class. All retained H1 intervals are finite and have positive
persistence. The essential H0 class remains stored and is reserved for explicit
exclusion when the later landscape input is formed.

## Admission and production resources

Two separate admission groups completed 180 views with zero failures at about
53.2 seconds and 266 MiB peak RSS per group. Full production then completed with
two heavy workers:

| Measure | Result | Guard |
|---|---:|---:|
| Worker time | 3,952.673 s (1.098 h) | 43,200 s |
| Median group | 52.887 s | 900 s |
| Maximum group | 54.641 s | 900 s |
| Peak process-tree RSS | 288,251,904 B (274.9 MiB) | 4 GiB |
| Private PH storage | 179,901,005 B (179.9 MB) | 1 GB |

## Independent correctness and reproducibility

Independent validation reconstructed 6,750 coordinate identities, record
identities, file hashes, diagram invariants, stored MST checks, and scope checks.
All passed. It also recomputed one complete Euclidean 384-cell MST in every
group; all 75 passed.

The separate admission and production copies of
`mv05h_group__SRA550660__20260805` match in all 90 identities, diagram hashes,
R objects, and serialized files. A completed-queue rerun rebuilt zero groups;
all 75 view-audit hashes, result-set hashes, record counts, and byte totals were
unchanged.

## Post-PH projection

| Component | Worker hours | State |
|---|---:|---|
| Integrated coordinates | 9.107 | Measured |
| Integrated cell PH | 1.098 | Measured |
| Integrated landscapes | 1.456 | Projected with 25% reserve |
| Integrated baseline/retrieval | 0.507 | Projected with 25% reserve |
| **Total** | **12.169** | **Passes 21.6-hour cap** |

Projected complete storage is 1,269,456,100 bytes (1.269 GB). The conservative
combined peak remains 3,303,829,504 bytes and the maximum upstream group remains
790.95 seconds. Every time, RSS, storage, and group-duration cap passes.

## Software and boundary gates

The focused MV5-H contract suite passes 12/12 expectations. The complete
source-loaded suite passes 590/590 expectations with zero failures and zero
warnings; two existing CRAN-only tests are skipped. A public-only temporary
source build and installed-package `R CMD check --no-manual` pass with
`Status: OK`, zero errors, zero warnings, and zero notes. No dependency was
installed and no lockfile or dependency declaration changed.

Landscape, distance, retrieval, clustering, gene-view, fusion, new-data, and
biological-outcome counters are all zero. Labels remained closed. Private
diagrams remain under `tmp/`; PDFs, reviewer material, private caches, and
`example_run.r` remain untracked. Nothing was pushed.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | Approve integrated cells-as-points H0/H1 |
| Correctness demonstrated? | Pass independently, deterministically, and under resume |
| Computation feasible? | Yes under all frozen caps |
| Biological interpretation permitted? | Prohibited |
| Next action | Specify exact all-active-level integrated H0/H1 landscape distances, then stop and reproject |

MV5-H establishes corrected integrated diagrams only. It does not show that
integration improves topology, retrieval, clustering, or biology.
