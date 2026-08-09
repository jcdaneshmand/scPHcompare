# MV5-F label-closed integration-induction resource gate audit

## Outcome

MV5-F is complete. The corrected D1-panel-preserving real-data mapping route is
valid, deterministic, and feasible under the frozen caps. The decision permits
a separately authorized full label-closed integrated execution sprint. No full
run or biological evaluation occurred in MV5-F.

## Corrected reconstruction

The accepted D0 caches contain matrix-only SCT data and selected-cell identities,
not reusable Seurat models. The earlier MV5-C pilot also selected a new
query-compatible panel. MV5-F instead validates all accepted source identities,
recovers the exact 384 cells, and reconstructs reference and query objects on the
fixed D1 panel. Missing held-out genes are dropped only from that query's anchor
feature set; no gene is replaced and the reference PCA is never refit.

## Pilot manifest and measurements

| Role | Study | Training/query | Aggregate missing | Elapsed (s) | Peak RSS (GiB) | Active genes |
|---|---|---:|---:|---:|---:|---:|
| Minimum query | SRA713577 | 89 / 1 | 0 | 325.764 | 3.176 | 500 |
| Maximum query | SRA779509 | 65 / 25 | 5 | 840.886 | 2.676 | 499–500 |
| Maximum missing | SRA826293 | 82 / 8 | 35 | 532.277 | 2.671 | 498–500 |
| Median nonzero missing | SRA749327 | 86 / 4 | 2 | 393.305 | 2.867 | 499–500 |

Selection used only D1 resource evidence. Tissue labels and MV5-E method results
were not read. All four groups completed under the caps; the stage took 2,118.7
seconds. All 360 coordinate views and 38 query mappings completed. Every view is
an ordered finite 384-by-30 matrix and every reference before/after hash matched.

Input validation took approximately 210 seconds per group; it conservatively
rehashes and reads all 90 raw and SCT inputs. This is projected separately from
reference fitting and query mapping.

## Independent and deterministic validation

The independent validator reconstructed payload digest, group identity,
coordinate axes, mapping identities, reference immutability, file hashes, and
scope. All 28 group-category checks passed without calling the production group
validator.

Two clean builds under the finalized implementation produced the same cache key,
payload hash, coordinate-set hash, query-mapping identities, and exact RDS
SHA-256 `aa7d7c46ef0c4fa19855fb51c57529a4d33658990ddfa38b00d91ced2d98cd39`.

## Complete workload projection

The D1 axis contains 6,300 training and 450 held-out sample instances across 75
groups. With a 25% reserve:

| Component | Worker hours |
|---|---:|
| Input validation | 5.524 |
| Reference SCT/PCA | 1.497 |
| Query SCT | 0.225 |
| Held-out mapping | 4.443 |
| Coordinate assembly | 0.002 |
| Cell PH | 1.308 |
| Landscape distances | 1.456 |
| Baseline/retrieval | 0.507 |
| **Total** | **14.963** |

The maximum MV5-C integrated-to-SCT finite-interval ratio was 0.963; the
projection conservatively used 1.0 before reserve. Projected storage is 1.537
GB, peak RSS 3.97 GiB, and worst-group time 1,051 seconds. All pass the 21.6-hour,
10-GB, 8-GiB, and 1,800-second caps.

## Software gates

The focused MV5-F test file passes 13/13 expectations. The complete source test
suite passes 614/614 expectations with no failures, warnings, or skips. A clean
source build and installed-package `R CMD check --no-manual` pass with
`Status: OK`, zero errors, zero warnings, and zero notes.

The first combined build/check wrapper expired after 20 minutes while its log
remained clean through namespace-isolation checks. The already built tarball was
then checked directly to completion in 217.8 seconds with `Status: OK`. No
dependency was installed and no project dependency or lockfile changed.

## Boundary audit

Label transfer, PH, landscapes, distances, clustering, gene views, fusion, new
data, biological outcomes, and full integrated production each have count zero.
Private coordinate bundles remain under `tmp/`. PDFs, reviewer material, caches,
and `example_run.r` remain untracked. Nothing was pushed.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | Approve fixed-D1-panel label-closed induction |
| Correctness demonstrated? | Pass independently and deterministically |
| Computation feasible? | Yes under frozen caps |
| Biological interpretation permitted? | Prohibited |
| Next action | Separately authorize all 75 integrated coordinate groups, then stop and reproject before PH |

MV5-F supports only technical feasibility and leakage-resistant execution. It
does not show that integration improves topology, retrieval, clustering, or
biology.
