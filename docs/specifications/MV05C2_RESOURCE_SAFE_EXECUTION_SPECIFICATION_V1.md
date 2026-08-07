# MV5-C2 resource-safe execution specification v1

## Document control

| Field | Value |
|---|---|
| Contract | `mv05c2_resource_safe_execution_v1` |
| Date frozen | 2026-08-06 |
| Status | Accepted for label-closed SCT cell-primary precomputation only |
| Samples | Frozen 90-sample, 15-study existing-data candidate set |
| Seeds | `20260805`–`20260809` |
| Landscape definition | Exact all-active-level L2; H0 and H1 separate |
| Outcome-label state | Closed |
| Biological/technical endpoints | Not computed |

## 1. Purpose and boundary

This specification replaces the prohibited naïve MV5-D execution architecture.
It does not change the scientific cell-topology definition, persistence
landscape, seeds, cells per sample, feature count, PCA dimension, or LOSO split.
It changes execution units and requested pair scope so that reusable work is
not repeated once per outer fold.

The only current conditional authorization is label-closed SCT cell-primary
precomputation. Gene topology, inductive integration, full sample-distance
matrices, clustering, within-study pair contrasts, endpoint calculation,
method ranking, fusion, and manuscript claims remain outside this contract.

## 2. Frozen execution graph

1. Select the same 384 cells per sample and seed from the immutable manifest.
2. Run `SCTransform` once per sample–seed, yielding 450 private cache entries
   instead of 6,750 fold-repeated normalizations.
3. Bind each cache entry to the raw-cache SHA-256, selected-cell SHA-256,
   normalization parameters, seed, sample ID, and Seurat version.
4. For each LOSO fold and seed, load validated cache records, fit the 500-gene
   panel, training standardization, and 30-PC model using training samples only,
   then project held-out cells without refitting.
5. Generate complete H0/H1 diagrams with structured failures and immutable
   provenance.
6. Request exact landscape distances only from each held-out query to every
   training reference within the same fold/seed/representation/view/dimension.
7. Execute pair requests in deterministic fold-local chunks of at most 250
   pairs. Atomically publish one output and one status artifact per chunk.
8. Resume only when request-subset and output hashes validate; stale, partial,
   mismatched, or unreadable artifacts stop rather than overwrite.

## 3. Normalization cache contract

`mv05c2_sample_seed_sct_cache_v1` is private data. Its public audit row records
only stable IDs, cache keys, filenames, sizes, hashes, timing, versions, and
label state. A cache identity changes if any of these change:

- sample ID or seed;
- ordered selected-cell digest;
- frozen private raw-cache digest;
- normalization method or parameters; or
- Seurat version.

The cache is reusable across LOSO folds because none of those inputs depends on
the held-out study. Feature selection, scaling, PCA, and integration fitting are
not reusable across folds and remain training-only.

Existing cache files are never silently replaced. Missing entries may be built
atomically. Valid entries may be reused. Stale or unreadable entries require an
explicitly new cache location or contract rather than mutation in place.

## 4. Pair-scope contract

`mv05c2_query_training_pair_scope_v1` contains only
`held_out_query_to_training_reference` requests. This scope is sufficient for:

- cross-study tissue mean reciprocal rank after labels are eventually opened;
- fixed 1-nearest-neighbor retrieval after labels are eventually opened; and
- label-closed held-out neighbor rankings.

It is not sufficient for:

- PAM, hierarchical, or spectral clustering over a complete sample matrix;
- same-study versus cross-study distance contrasts;
- within-held-out-study sample-pair contrasts; or
- distance normalization that requires a complete fold matrix.

Those tasks are deferred, not approximated. They require a separately approved
pair-scope extension after the primary retrieval bundle is immutable.

Pair scope never changes the landscape calculation itself. Every requested H0
or H1 distance uses every active landscape level and exact critical-pair
integration with zero numerical error estimate.

## 5. Chunk and resource guards

| Guard | Value |
|---|---:|
| Heavy workers | At most 2 |
| Sample/cache or fold job elapsed cap | 1,800 seconds |
| Process RSS cap | 8 GiB |
| Landscape requests per chunk | At most 250 |
| Nominal aggregate cap | 24 worker-hours |
| Planning cap with reserve | 21.6 worker-hours |
| Projected private SCT cache budget | 40 GB maximum |

Landscape objects may be retained only across adjacent chunks in the same
stratum and homology dimension. They are cleared at a stratum boundary so a
long queue does not accumulate every landscape in memory.

## 6. Measured evidence

- Six real sample–seed cache entries built in 514.4 seconds at 4.72 GB peak
  RSS; validated resume took 67.41 seconds at 0.49 GB.
- All six cached SCT matrices were exactly identical to immutable MV5-C.
- Both cached LOSO directions reproduced every eligible cell/gene persistence
  diagram and baseline matrix exactly; fold computation took 10.382–17.718
  seconds after cache loading.
- The exact chunk engine reproduced 250/250 MV5-C landscape distances with
  maximum absolute difference zero.
- A second queue run built zero chunks, reused all 50, and left all 100
  output/status files byte-identical.
- The full plan reduces normalization operations 15-fold and landscape rows
  from 1,802,250 to 212,100 (8.497-fold).

## 7. Scenario gate

| Scenario | Lower-bound worker-hours | Disposition |
|---|---:|---|
| Naïve all-view MV5-D | 242.52 | Prohibited |
| Resource-safe all planned views | 25.82 before integrated mapping | Prohibited |
| Resource-safe SCT cell + gene | 22.25 | Prohibited; insufficient reserve |
| Resource-safe SCT cell primary | 18.68 | Conditional label-closed go |

The 18.68-hour projection does not authorize outcomes. It authorizes a staged
SCT cell-primary precomputation only, with a stop and re-projection after the
450-entry normalization cache and before all fold/seed diagrams are launched.

## 8. Acceptance and stop conditions

The next execution stage must stop if:

- projected cache storage exceeds 40 GB;
- any sample–seed normalization exceeds 30 minutes or 8 GiB;
- the completed cache projection exceeds 21.6 aggregate worker-hours;
- cached-fold equivalence no longer holds;
- labels enter any fit, cache, request, chunk, or method-selection artifact;
- a stale/partial artifact would need overwrite; or
- the full SCT cell-primary projection no longer retains the 10% reserve.

G-MV5 remains open. Outcome labels remain closed until all primary predictions
and requested distance artifacts are immutable and a separate MV5-E gate is
approved.
