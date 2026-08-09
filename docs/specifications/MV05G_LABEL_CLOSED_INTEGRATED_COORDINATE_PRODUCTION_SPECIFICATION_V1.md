# MV5-G label-closed integrated-coordinate production specification v1

## Document control

| Field | Value |
|---|---|
| Contract ID | `mv05g_label_closed_integrated_coordinate_production_v1` |
| Dates | 2026-08-08 to 2026-08-09 |
| Status | Executed and independently validated |
| Source revision | `1ea094bc7ad78cc140a402646814ab79bd355258` |
| Existing-data axis | 90 samples, 15 studies, five seeds |
| Mapping axis | Reference-only SCT transfer projection, 500 genes, 384 cells, 30 PCs |
| Outcome state | Closed |
| Integrated PH execution | Not executed; separately authorized after measured reprojection |

## 1. Purpose and stop boundary

MV5-G executes the complete 75-group integrated-coordinate cache authorized by
MV5-F. It preserves the accepted D1-panel label-closed induction contract and
measures the real cost of all integrated coordinate groups before any integrated
persistence calculation.

This sprint permits source validation, reference/query reconstruction, held-out
mapping, ordered coordinate assembly, immutable private caching, resource
measurement, independent validation, deterministic repeat, resume proof, and
post-coordinate projection. It prohibits label transfer, persistence homology,
landscapes, distances, retrieval, clustering, gene topology, fusion, new data,
biological outcomes, and claim promotion.

## 2. Frozen axes

Production must preserve:

- all 15 leave-one-study-out folds and seeds `20260805:20260809`;
- the exact D1 training-selected 500-gene panel for each fold-seed;
- the exact D0-selected 384 cells for every sample-seed;
- joint training-only SCT and exact 30-PC reference PCA;
- independent held-out query SCT and SCT transfer projection;
- query-specific ordered intersections with the frozen panel, without feature
  replacement;
- 90 ordered 384-by-30 coordinate matrices per group;
- no tissue-label reads or label transfer; and
- the accepted MV5-F cache and provenance identity.

MV5-E results cannot select folds, features, exclusions, seeds, mapping settings,
or downstream authorization thresholds.

## 3. Execution and cache contract

The public manifest contains exactly 75 unique fold-seed groups, ordered by
study and seed. One heavy worker executes groups sequentially. A completed group
is accepted only when its private RDS validates against its manifest row and
all input identities. Valid groups are reused without modification on resume.

Each private payload records coordinate and mapping provenance but excludes
timing so a clean same-implementation rebuild can be byte identical. Every
group identity binds fold, seed, sample axes, panel and dimension settings, D1
payload identity, raw and D0 source identities, mapping implementation identity,
runtime versions, reference identity, query mapping identities, and ordered
coordinate hashes.

## 4. Resource guards

The frozen guards are one heavy worker, 1,800 seconds and 8 GiB process-tree RSS
per group, 43,200 seconds measured worker time for the coordinate stage, and
10 GB private result storage. A group failure, reference mutation, malformed
coordinate view, resource breach, or prohibited counter stops acceptance.

After all coordinates complete, a new projection must combine measured MV5-G
coordinate cost with the accepted D3, D4, and D5 rates, each carrying the frozen
25% reserve. The complete integrated cell-primary route must remain below 21.6
worker-hours, 8 GiB RSS, 1,800 seconds per group, and 10 GB storage.

## 5. Required validation

An independent validator must reconstruct, without calling the production group
validator:

- all group and private-file identities;
- ordered coordinate axes and finite 384-by-30 matrices;
- all query mapping identities, active feature sets, and anchor counts;
- reference immutability;
- cache/payload/coordinate/file hashes;
- exact study, seed, sample, and mapping counts; and
- absence of all prohibited jobs and outcomes.

One complete production group must match a clean same-implementation admission
build byte for byte. A completed-queue rerun must rebuild zero groups and leave
all 75 cache, payload, coordinate, and file hashes unchanged.

## 6. Decision rule

Integrated cell PH is authorized only as a separate sprint if all 75 groups
pass, repeat and resume evidence pass, no prohibited work occurs, and the
measured reserved total passes every cap. Authorization does not execute PH and
does not authorize landscapes, distances, retrieval, outcomes, clustering,
gene topology, fusion, new data, or manuscript claims.
