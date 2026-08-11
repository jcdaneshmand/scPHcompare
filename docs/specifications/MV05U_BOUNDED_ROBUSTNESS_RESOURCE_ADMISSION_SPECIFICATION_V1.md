# MV5-U bounded robustness resource-admission specification v1

Date frozen: 2026-08-10

Accepted gate base: `47235d8`

Biological or technical outcomes: prohibited

## Purpose

MV5-U executes only the 24 label-closed admission units frozen by MV5-T. Its
purpose is to determine whether the admitted one-factor-at-a-time coordinate
perturbations can be generated, validated, streamed, and resumed within the
frozen resource envelope. It does not estimate robustness effects, calculate
retrieval or clustering outcomes, compare configurations, rank methods, or
authorize full robustness production.

## Frozen unit axis

The exact axis is three folds by two representations by four configurations at
seed `20260805`:

- minimum training: `large_loso_v1:SRA779509`, n=65;
- median training: `large_loso_v1:SRA703206`, n=86;
- maximum training: `large_loso_v1:SRA713577`, n=89;
- representations: `sct_whole` and `inductive_integrated`;
- configurations: nested 192 cells/30 PCs/Euclidean, nested 256 cells/30
  PCs/Euclidean, 384 cells/20 PCs/Euclidean, and 384 cells/30 PCs/cosine chord.

Every unit constructs all 90 sample views. Nested subsets use the frozen
SHA-256 cell order shared across representations. The 20-PC condition truncates
accepted coordinates. Cosine chord row-normalizes the accepted 30-PC cell
vectors before applying Euclidean PH and energy calculations.

## Computation within one unit

Each unit must:

1. verify the committed MV5-T queue, all 150 private source hashes, and the
   selected SCT/integrated source pair before reading coordinates;
2. construct and validate all 90 transformed typed views;
3. run corrected `ripserr` Vietoris-Rips PH through H1 on every view with
   threshold `-1` and field 2;
4. compare every finite H0 death multiset with an independently calculated
   Euclidean MST;
5. exclude the single essential H0 interval and stage exact, all-active,
   consecutive-level H0 and H1 persistence landscapes for all 90 views;
6. calculate exact landscape L2 distances for a frozen 32-pair coverage set in
   each homology dimension; and
7. calculate matched empirical energy distance for the same 32 pairs.

The 32-pair coverage set contains 16 training-training unordered pairs and 16
held-out-training directed pairs. Pair selection is the canonical first 16
under SHA-256 order of fold, seed, scope, and sample identities. It exercises
both later retrieval scopes while remaining a resource admission rather than a
complete robustness analysis.

## Landscape definition

The dissertation-aligned corrected definition is immutable in MV5-U:

- H0 and H1 remain separate;
- the essential H0 interval is excluded;
- every consecutive active landscape level is retained;
- no fixed level cap or uniform grid is used; and
- L2 integration is exact over linear critical-pair segments.

Landscape summaries and pair distances are deterministic private artifacts.
No dense landscape matrix is saved.

## Artifact and resume contract

Every unit publishes one immutable directory only after all stages validate.
It contains deterministic interval, view, energy, landscape-summary, and
landscape-pair artifacts; a resource-only timing record; and a status manifest
binding every file hash and size to the exact queue, sources, implementation,
and prospective commit. A missing directory is buildable. A complete valid
directory is reusable. Any partial, stale, or hash-invalid directory is a hard
failure and is never overwritten automatically.

A clean repeat uses an independent private root and must reproduce every
deterministic scientific artifact byte-for-byte. Resource timings may vary but
must satisfy the same schema and caps. Resume validation must reuse all 24 units
without changing any path, hash, size, or modification timestamp.

## Independent validation

Validation is separate from the execution runner and must cover:

- all committed and private source identities;
- exact 24-unit queue axes and configuration isolation;
- 90 transformed views and expected point/coordinate counts per unit;
- nested cell membership across the 192/256 conditions;
- exact first-20 coordinate identity;
- cosine row norms and absence of zero/nonfinite values;
- all PH interval invariants and all-view H0 MST agreement;
- an analytic square H1 interval `[1, sqrt(2)]` and its exact landscape
  squared norm `(sqrt(2)-1)^3/12`;
- independent direct-loop energy recalculation on frozen sampled pairs;
- exact artifact hashes, clean repeat, immutable resume, resource caps, and
  public label safety.

## Resource and abort gates

Execution uses one worker. Abort the queue on any unit above 600 seconds or 4
GiB process-tree RSS, total worker time above two hours, or new private storage
above 2 GiB. Abort also on source drift, label access, outcome calculation,
shape/cell/configuration leakage, zero-norm cosine rows, PH/MST/landscape/energy
disagreement, partial publication, repeat mismatch, or resume mutation.

Completed immutable units are preserved after an abort. No automatic retry or
overwrite is allowed.

## Decision boundary

Even a fully passing admission sets `full_robustness_authorized=FALSE`. Its only
positive next action is to prospectively freeze a separate streamed full-
robustness execution gate using measured unit resources. A failed admission
stops for review.

MV5-U cannot calculate labels or outcomes, retrieval endpoints, clusters,
p-values, ranks, or winners; expand the perturbation grid; alter accepted
MV5-D through MV5-T artifacts; track private caches, labels, PDFs, reviewer
material, or `example_run.r`; push; add spectral, gene, fusion, new-data,
optimization, or Rust work; change package defaults; or make manuscript claims.
