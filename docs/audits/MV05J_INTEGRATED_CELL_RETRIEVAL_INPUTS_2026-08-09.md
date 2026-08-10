# MV5-J label-closed integrated cell retrieval-input completion audit

## Outcome

MV5-J is complete. All 75 frozen fold-seed groups contain five immutable,
label-closed integrated retrieval-input methods for all 35,350 directed
held-out-query-to-training-reference biological pairs. The 176,750 public
ranking rows and 375 method-completion records independently validate with no
method failure, exact-distance tie, held-out-derived component-scale fit, or
prohibited downstream job.

This sprint did not open tissue or other biological outcomes and did not compare
the integrated view with the accepted SCT result. It authorizes only a later,
separately specified prediction-locked integrated retrieval evaluation.

## Frozen inputs and method panel

Each group binds one accepted MV5-D1 fold, one read-only MV5-D5 mean-profile
bundle, one accepted MV5-G integrated-coordinate record, the matching MV5-I
component subset, and the implementation hash.

The five methods are:

1. raw exact all-active-level integrated H0 landscape L2;
2. raw exact all-active-level integrated H1 landscape L2;
3. raw `sqrt(H0^2 + H1^2)` as descriptive secondary only;
4. empirical V-statistic energy divergence on the exact same 384-cell by
   30-coordinate MV5-G matrices; and
5. the same MV5-D1 training-selected, training-standardized 500-gene
   pseudobulk context baseline used by the SCT contract.

No SCT cell distance or ranking was reused. H0 and H1 remain separate primary
topology methods. Because MV5-I contains no within-training topology scope,
MV5-J fitted no component scale from held-out queries and launched no topology
expansion.

## Prospective admission

The specification and ADR were frozen before running real groups. Admission
covered the representative `SRA550660` / `20260805` group (425 biological
pairs) and the prospectively identified maximum `SRA779509` / `20260805` group
(1,625 pairs).

Both groups were calculated twice into fresh paths and produced byte-identical
RDS bundles. Seven spread-out direct energy and seven pseudobulk oracles passed
for each group. H0 and H1 agreed exactly with MV5-I; the maximum raw-composite
roundoff was `9.379165e-13`, below the frozen `1e-12` tolerance. All partitions,
ranks, identities, tie fields, and label-firewall checks passed.

## Production resources

| Measure | Accepted value |
|---|---:|
| Fold-seed groups | 75 / 75 |
| Fixed seeds | 5 |
| Biological pairs | 35,350 |
| Methods | 5 |
| Ranked rows | 176,750 |
| Completed method groups | 375 / 375 |
| Failed method groups | 0 |
| Aggregate worker elapsed | 1,896.369 seconds |
| Aggregate scientific operation time | 737.767 seconds |
| Observed queue wall time | 1,062.4 seconds |
| Median group elapsed | 23.645 seconds |
| Maximum group elapsed | 45.440 seconds |
| Peak process-tree RSS | 301,789,184 bytes (287.8 MiB) |
| Private bundle/audit bytes | 209,877,213 |

Production used two workers. Every group remained below the 600-second and
4-GiB guards; aggregate worker time remained below 7,200 seconds; and new
private output remained below 2 GiB. No guard fired and no process was killed.

## Independent validation

Every private bundle, audit, MV5-D1 fold, MV5-D5 mean-profile bundle, and MV5-G
coordinate record was reopened and rehashed. The validator checked all query
and training partitions, method axes, completion records, pair and ranking IDs,
canonical ranks, tie fields, label closure, and public/private identity joins.

All 35,350 H0 rows and all 35,350 H1 rows agree exactly with MV5-I. The raw
composite agrees within `9.805490e-13`. Across three spread-out pairs per group,
225 direct integrated-energy oracles agree within `2.220446e-15`, and all 225
independently reconstructed pseudobulk oracles agree exactly. All 176,750 pair
IDs and 176,750 ranking IDs are unique. There are zero exact-distance tie rows.

## Determinism and resume

- Representative and maximum complete-group repeats are byte-identical.
- A completed-queue resume reopened and validated every group without launching
  a replacement calculation.
- All 75 private bundles and 75 audits retained identical paths, hashes, and
  sizes after resume.
- An independent public rebuild produced byte-identical copies of all six
  public artifacts, including timestamp-free gzip rankings.

## Public evidence

| Artifact | Purpose |
|---|---|
| `mv05j-integrated-cell-retrieval-rankings-2026-08-09.csv.gz` | 176,750 immutable label-closed distance/rank rows |
| `mv05j-method-completion-2026-08-09.csv` | Five-method completion/failure ledger per group |
| `mv05j-group-bundle-index-2026-08-09.csv` | MV5-D1/D5/G/I and private-bundle identities |
| `mv05j-method-registry-2026-08-09.csv` | Frozen roles, formulas, scales, and representations |
| `mv05j-component-scale-disposition-2026-08-09.csv` | Explicit no-leakage H0/H1 scale decision |
| `mv05j-public-assembly-summary-2026-08-09.csv` | Counts, hashes, and zero-job counters |
| `mv05j-production-group-resources-2026-08-09.csv` | Per-group guard and resource measurements |
| `mv05j-independent-group-validation-2026-08-09.csv` | Per-group independent numerical checks |
| `mv05j-independent-validation-summary-2026-08-09.csv` | Complete all-group validation summary |
| `mv05j-resume-validation-2026-08-09.csv` | 150-file immutable resume proof |
| `mv05j-public-repeat-validation-2026-08-09.csv` | Six-file public repeat proof |
| `mv05j-integrated-retrieval-evaluation-authorization-2026-08-09.csv` | Mandatory next-stage decision and boundary |

Private folds, coordinates, mean profiles, production bundles, logs, PDFs,
reviewer material, and `example_run.r` remain untracked.

## Package verification

The focused MV5-J file passed all 10 expectations. The complete package-loaded
suite passed 654 expectations with zero failures or warnings; two longstanding
CRAN-gated integration checks were skipped. A fresh prefiltered source build and
`R CMD check --no-manual` completed with `Status: OK` and zero errors, warnings,
or notes. No dependency was installed or changed.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | Approve the five-method integrated retrieval-input panel; H0/H1 separate and raw composite descriptive |
| Correctness demonstrated? | Pass |
| Computation feasible? | Yes; complete production is far below all frozen guards |
| Biological interpretation permitted? | Prohibited in MV5-J; labels remained closed |
| Next action | Separately specify prediction-locked integrated retrieval evaluation |

## Next-stage boundary

The next sprint may open frozen labels only after verifying the public ranking
hash and all 375 completion records. It may calculate prespecified held-out
retrieval endpoints and blocked uncertainty. It may not refit, rescale, rerank,
tune, select, or replace methods; compare integrated results with SCT until the
new evaluation is complete and locked; cluster; construct gene topology; fuse
views; introduce new data; or promote manuscript claims.
