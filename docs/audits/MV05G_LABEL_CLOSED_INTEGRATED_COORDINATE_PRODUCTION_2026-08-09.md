# MV5-G label-closed integrated-coordinate production audit

## Outcome

MV5-G is complete. All 75 fixed-panel label-closed integrated coordinate groups
are valid, deterministic, resumable, and within the frozen resource caps. The
measured post-coordinate projection authorizes a separate integrated cell-PH
sprint. No PH or biological evaluation occurred in MV5-G.

## Production axis and result

The manifest fixes 15 leave-one-study-out folds crossed with seeds
`20260805:20260809`. Each group reconstructs a training-only reference on its
exact D1 500-gene panel and exact D0-selected cells, maps every held-out query
without label transfer, and emits 90 ordered 384-by-30 coordinate matrices.

| Measure | Result |
|---|---:|
| Groups | 75 / 75 |
| Studies / seeds | 15 / 5 |
| Coordinate views | 6,750 |
| Held-out query mappings | 450 |
| Failures | 0 |
| Active query features | 496–500 |
| Anchors | 145–907 |
| Reference mutations | 0 |

Every unavailable query gene was omitted only from that query's active anchor
features; none was replaced and the full reference PCA remained fixed.

## Resource evidence

The sequential one-worker stage used 32,786.38 worker-seconds (9.107 hours).
Maximum process-tree RSS was 3,303,829,504 bytes (3.077 GiB), the slowest group
used 790.950 seconds, and private results used 806,423,810 bytes (769.1 MiB).
These pass the 12-hour coordinate-stage, 8-GiB per-group RSS, 1,800-second
per-group, and 10-GB storage guards.

## Independent, deterministic, and resume validation

The independent validator reconstructed all identities, axes, coordinates,
mapping records, active features, reference immutability, hashes, and scope.
All 675 checks in nine validation categories passed across all 75 groups.

A clean admission build and the production result for
`mv05g_group__SRA550660__20260805` are byte identical at SHA-256
`eacde0ed761d81a6a347876b1b48fbdb5b0adfb14cef845ac5ae6cdb2926b7d7`.
Re-running the completed production queue exited successfully, rebuilt zero
groups, and preserved all 75 cache keys, payload hashes, coordinate hashes, and
file hashes.

## Post-coordinate projection

| Component | Worker hours | Evidence basis |
|---|---:|---|
| Integrated coordinates | 9.107 | Complete MV5-G measurement |
| Integrated cell PH | 1.308 | Accepted D3 rate + 25% reserve |
| Integrated landscapes | 1.456 | Accepted D4 rate + 25% reserve |
| Integrated baseline/retrieval | 0.507 | Accepted D5 rate + 25% reserve |
| **Projected total** | **12.379** | **Measured + reserved future work** |

Projected total storage is 1,334,921,626 bytes (1.335 GB). Measured peak RSS
and maximum group duration remain the applicable conservative guards. The
projection passes the 21.6-worker-hour, 10-GB, 8-GiB, and 1,800-second caps.

## Software gates

The focused MV5-G test context passes 10/10 expectations. The complete
source-loaded suite passes 578/578 expectations with zero failures and zero
warnings; two existing CRAN-only tests are skipped. Shell syntax and Git
whitespace checks pass.

A public-only temporary staging tree was built into a source tarball and checked
as an installed package with `R CMD check --no-manual`. The result is
`Status: OK` with zero errors, zero warnings, and zero notes. The first
in-worktree build wrapper exceeded ten minutes while R copied Git-ignored
private artifacts before applying `.Rbuildignore`; the bounded public-only
staging route excludes those artifacts before the source build. No dependency
was installed and no lockfile or dependency declaration changed.

## Boundary audit

Every label-transfer, PH, landscape, distance, retrieval, clustering, gene-view,
fusion, new-data, and biological-outcome counter is zero. Tissue labels remained
closed. Private coordinate bundles remain under `tmp/`. PDFs, confidential
reviewer material, private caches, and `example_run.r` remain untracked. Nothing
was pushed.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | Approve fixed-D1-panel label-closed coordinates |
| Correctness demonstrated? | Pass independently, deterministically, and under resume |
| Computation feasible? | Yes under all frozen caps |
| Biological interpretation permitted? | Prohibited |
| Next action | Specify and execute separate integrated cell H0/H1 PH, then stop and reproject |

MV5-G establishes a trustworthy integrated coordinate substrate. It does not
show that integration improves topology, retrieval, clustering, or biology.
