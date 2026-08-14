# MV5-C2 resource-safe execution audit

| Field | Value |
|---|---|
| Date | 2026-08-06 |
| Contract | `mv05c2_resource_safe_execution_v1` |
| Real-data validation | Frozen six-sample/two-study MV5-C artifacts |
| Full planning set | 90 samples, 15 studies, 5 seeds |
| Outcome-label state | Closed |
| Endpoints computed | None |
| Gate | SCT cell-primary conditional go; all-view MV5-D prohibited |

## Outcome

MV5-C2 implemented and tested the two resource changes required by ADR-006:
sample–seed normalization reuse and exact query-to-training landscape chunks.
The implementation preserves the dissertation-aligned landscape definition:
all active levels, exact integration, and separate H0/H1 results.

The sprint does not make the full method panel feasible under the 24-hour cap.
It narrows the next executable stage to label-closed SCT cell-primary
precomputation and leaves gene, integration, clustering, technical within-study
contrasts, and all outcomes gated.

## Implemented contracts

- immutable `mv05c2_sample_seed_sct_cache_v1` records;
- deterministic 90-sample/15-fold/5-seed execution plan;
- `mv05c2_query_training_pair_scope_v1` with explicit unsupported tasks;
- fold/seed/representation/view/dimension-local chunks of at most 250 pairs;
- exact Persim critical-pair execution with stratum-bounded landscape reuse;
- atomic output/status publication and hash-validated resume;
- stale, corrupt, partial, label-open, and oversized request refusal; and
- scenario projection with a 10% planning reserve.

## Real-data equivalence

The normalization pilot built one seed for all six MV5-C samples. The six
private cache records total 368,898,488 bytes and remain ignored. Every cached
500-gene-by-384-cell SCT slice exactly matches the immutable MV5-C fold bundle.

Both LOSO directions were then reconstructed from cache. All eligible cell and
gene diagrams and all matched baseline matrices have maximum absolute
difference zero from MV5-C. Fold computation took 17.718 seconds for the
five-query direction and 10.382 seconds for the one-query direction, compared
with about 301–306 seconds for the original SCT phase that included repeated
normalization.

The exact landscape pilot selected only held-out-query/training-reference pairs
from the immutable MV5-C distances. It produced 250/250 exact distances across
50 local chunks. Every distance equals its MV5-C reference exactly. Persim
emitted its known internal zero-width slope warning on some landscape
construction, but all outputs remained finite and passed exact reference
comparison. A second run built zero chunks, reused all 50, and left all 100
output/status files byte-identical.

## Resources

| Measurement | Result |
|---|---:|
| Six-cache build | 514.4 s; 4,722,900,992 B peak RSS |
| Six-cache validated resume | 67.41 s; 490,651,648 B peak RSS |
| Cache build/resume speedup | 7.63× |
| Projected 450-cache storage | 27,667,386,600 B |
| Two cached SCT folds | 60.14 s process; 1,392,304,128 B peak RSS |
| Cached fold calculation | 10.382–17.718 s/fold |
| Exact landscape pilot | 250 pairs; 87.31 chunk-seconds |
| Largest landscape chunk | 18.33 s; 323,506,176 B peak RSS |

The installed Seurat environment again lacked `glmGamPoi`; the existing native
fallback was retained, no dependency was installed, and the lockfile was not
changed.

## Full-plan reduction and gate

The frozen 90 candidate samples all have at least 421 filtered cells, so the
384-cell contract does not exclude a candidate before execution. The new plan
contains 450 normalization entries, 75 fold/seed jobs, 212,100 exact
query-to-training landscape requests, and 1,080 bounded pair chunks.

| Scenario | Lower bound | Gate |
|---|---:|---|
| Naïve full panel | 242.52 h | Prohibited |
| Resource-safe full panel | 25.82 h before integration mapping | Prohibited |
| Resource-safe SCT cell + gene | 22.25 h | Prohibited; no reserve |
| Resource-safe SCT cell primary | 18.68 h | Conditional go |

The conditional go authorizes only the label-closed SCT cell-primary artifact
stage. After the 450 normalization caches are built, actual runtime, storage,
failures, and sample coverage must be reprojected before the 75 fold/seed jobs.

## Boundaries and next action

- No tissue/approach endpoint join, method winner, fusion, or biological claim
  was created.
- Private counts, cache objects, chunk outputs, PDFs, reviewer material, and
  `example_run.r` remain outside Git.
- Query-to-training pairs cannot support complete-matrix clustering or
  within-study contrasts; those remain deferred rather than approximated.
- Integrated mapping has no measured full-scale lower bound and remains
  prohibited.
- G-MV5 remains open.

The next sprint is MV5-D0 stage 1: build the 450-entry SCT normalization cache
under the frozen guards, publish only its audit/hashes, and stop for resource
reprojection before any full fold execution or outcome calculation.
