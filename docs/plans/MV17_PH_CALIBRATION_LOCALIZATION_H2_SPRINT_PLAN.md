# MV17 PH calibration, localization, and H2 extension sprint plan

Date: 2026-08-26
Status: proposed execution branch; no scientific execution authorized by this document alone

## Goal

Extend the corrected all-QC persistent-homology program in three directions:
calibrate H0/H1 against prospectively frozen nulls, localize the observations
and features responsible for persistent structure, and determine whether H2 is
a reproducible, interpretable, computationally feasible, and incrementally
useful representation. A negative result at any gate is a complete result.

## Starting evidence

- MV13-D/E closes 636 all-QC cell views and 1,272 H0/H1 PH records after exact
  independent PCA/view/diagram reconstruction.
- MV14 closes 76,372 all-active-level cell-landscape distances with 14 exact-R
  oracles and no grid or universal level cap.
- MV15/MV16 closes and reviews 36 cell-distance comparisons: cell topology is
  stable across seeds and panels but complementary, not interchangeable, with
  the Pearson-residual gene topology.
- The historical fusion result is weight-sensitive and does not authorize a
  selected fusion weight.

## Non-negotiable boundaries

- H0, H1, and H2 remain separate primary objects. No dimension is silently
  pooled, weighted, or selected.
- The corrected complete-VR/all-finite-positive-interval contract remains the
  H0/H1 reference. Essential H0 is excluded from finite summaries.
- Existing closed PH, landscape, and comparison artifacts are immutable.
- Null families, seeds, sentinels, metrics, tolerances, stopping rules, and
  reporting tables are frozen before their values are inspected.
- Sample IDs and unit-level assignments remain private; public evidence is
  aggregate-only unless an already accepted policy explicitly permits more.
- Labels, outcomes, tissues, biology, clustering, view ranking, fusion-weight
  tuning, and manuscript claims remain closed through MV17-A to MV17-F.
- H2 is not a rescue analysis. It advances only by passing meaning, resource,
  stability, redundancy, and incremental-value gates.
- Execution is isolated from the main line by distinct scripts, private/public
  roots, manifests, receipts, and tests. Heavy jobs use one worker unless a
  later prospective resource amendment explicitly authorizes otherwise.

## Branch and merge policy

Implement MV17 in a dedicated worktree/branch from the immutable D-250 head.
Do not edit active MV13-MV16 production scripts or reuse their output roots.
Inputs are referenced by exact SHA-256, byte size, schema, and closure audit.
Each execution sprint requires a committed prefreeze, fresh atomic roots, an
independent reconstruction, resource reconciliation, and an explicit next-step
decision. Merge only closed audit artifacts and reusable implementation; never
merge private paths, identifiers, caches, partials, or exploratory outputs.

## MV17-A — source inventory and estimand freeze

Purpose: bind the exact eligible H0/H1 cell and gene corpora and define the
scientific questions without reading new result values.

Tasks:

1. Rehash the accepted MV8 gene and MV13/MV14 cell closure artifacts.
2. Inventory exact diagram, distance, selected-cell, PCA, and feature axes.
3. Define separate cell and gene calibration estimands for H0 and H1.
4. Define localization outputs for H0 mergers and H1 representative cycles.
5. Define H2-positive, H2-negative, and shuffled-null fixtures.
6. Publish only identities, counts, hashes, schemas, and prospective gates.

Acceptance: independent source reconstruction; zero null, localization, H2,
label, outcome, clustering, fusion, or claim computation.

## MV17-B — null-model qualification

Purpose: admit null generators before applying them to the real corpus.

Prospectively compare at least four families where mathematically compatible:

- coordinate-wise permutation;
- covariance-preserving Gaussian simulation;
- density/sample-size-matched geometric clouds; and
- within-unit feature or observation shuffles that preserve the frozen axis.

Each family must state exactly which structure it destroys and preserves. Use
fixed seeds and analytic or synthetic fixtures. Validate cardinality, marginal
properties, determinism, separation of positive/negative controls, and privacy.
Reject nulls whose preserved structure does not match the intended question.

Acceptance: primary and repeat artifacts agree scientifically; an independent
validator reconstructs every admitted null; no real-corpus calibration runs.

## MV17-C — bounded H0/H1 calibration sentinel

Purpose: measure whether null calibration is feasible and well behaved before
whole-corpus production.

Freeze label-blind minimum-, median-, and maximum-burden cell and gene units,
with H0/H1 separate. For every admitted null, prospectively fix replicate
counts, persistence summaries, empirical-tail handling, Monte Carlo precision,
time/RSS/storage caps, and abort rules. Candidate summaries include finite
feature count, total persistence, maximum persistence, Betti-curve integrals,
landscape norm, and observed-to-null distance; the prefreeze selects the exact
set before execution.

Acceptance: deterministic source identities, independent null regeneration,
valid empirical calibration behavior, complete resource receipts, and a
measured projection. A separate decision is required before full calibration.

## MV17-D — H0/H1 localization qualification

Purpose: connect persistent structure to influential cells, genes, samples, or
neighborhoods without selecting examples by labels or outcomes.

Tasks:

1. For H0, recover the MST/component merger edges corresponding to deaths.
2. For H1, compare at least two representative-cycle/localization algorithms
   on analytic circles, loops with noise, and negative controls.
3. Freeze tie handling, canonicalization, coefficient field, scale selection,
   and unit/feature influence summaries.
4. Test invariance to row order and stability across fixed seeds, cell-count
   perturbations, PC counts, and metric sensitivities.
5. Select real sentinels only by prospective burden/identity criteria.

Acceptance: known fixtures localize correctly; repeats and independent checks
agree; unstable or non-unique representatives are reported as such. Full-corpus
localization requires its own production prefreeze.

## MV17-E — H2 meaning and implementation gate

Purpose: prove that the implementation detects known voids and rejects known
non-voids before any real-data H2 interpretation.

Frozen fixture classes:

- sphere: H2 positive;
- torus: H1 and H2 positive with known Betti structure;
- circle: H1 positive and H2 negative;
- Gaussian/density-matched clouds: negative or null-calibrated controls; and
- shuffled versions of each admitted positive fixture.

Compare the primary backend with an independent engine on bounded fixtures.
Require correct dimension, interval validity, determinism, expected qualitative
ordering, and prospectively fixed numerical tolerances. Record simplex counts,
elapsed time, peak RSS, and failure behavior.

Acceptance: every positive/negative fixture gate passes twice and independently.
Failure closes real-data H2 while leaving H0/H1 calibration/localization open.

## MV17-F — H2 resource sentinel

Purpose: establish the real-data domain in which complete-VR H2 is credible.

1. Profile a deterministic synthetic grid over point count and dimension.
2. Freeze label-blind minimum-, median-, and maximum-burden real cell and gene
   sentinels from the already closed all-QC axes.
3. Record simplex growth, diagram size, time, peak RSS, storage, and stderr.
4. Predefine timeout/RSS/storage stops and one exact fallback attempt, if any.
5. If complete VR is infeasible, stop and prospectively evaluate a sparse or
   subsampled filtration as a different estimand; do not silently substitute it.

Acceptance: primary/repeat agreement, independent bounded-engine agreement,
and a conservative full-corpus projection. Full H2 remains closed.

## MV17-G — full H0/H1 calibration and localization

Eligibility: MV17-C and MV17-D independently close and a separate production
prefreeze is committed.

Run only the admitted nulls and localization methods over the exact frozen
corpora. Publish complete aggregate summaries, including failures and unstable
representatives. Independently reconstruct a prospectively fixed subset plus
all aggregate calculations. Do not rank dimensions, views, units, or tissues.

## MV17-H — real-data H2 stability and redundancy

Eligibility: MV17-E/F close, resource projection passes, and an exact queue is
prospectively frozen.

Compute H2 separately for cell and gene views. Evaluate only fixed perturbations
of seed, point count, PC count, metric, and preprocessing. Compare H2 with H0,
H1, energy, and pseudobulk using prespecified global and local agreement metrics.
No labels or outcomes are opened. Retain H2 only if it is repeatable and not
merely a numerical proxy for an existing representation.

## MV17-I — prediction lock and incremental evaluation

Eligibility: MV17-G/H independently close and the owner explicitly authorizes
an outcome-bearing evaluation.

Freeze predictions and all comparisons before labels are opened. Compare H0/H1,
H2 alone, H0/H1/H2, and matched conventional baselines under the established
sample-blocked framework. No representation or weight is selected from outcome
performance. H2 advances only with reproducible incremental value under the
prospective rule; null or negative results remain primary evidence.

## Parallel execution schedule

After MV17-A, the following may proceed concurrently in separate roots:

- lane 1: MV17-B then MV17-C;
- lane 2: MV17-D;
- lane 3: MV17-E then synthetic profiling in MV17-F.

Real-data H2 sentinels in MV17-F must use immutable D-250-era source closures.
MV17-G waits for lanes 1 and 2. MV17-H waits for lane 3. MV17-I waits for both
MV17-G and MV17-H plus an explicit owner decision. CPU-light implementation and
fixture work may overlap; memory-heavy real-data children must be serialized.

## Stop rules

Stop and preserve evidence on identity drift, partial or ambiguous state,
unowned stderr, resource breach, nondeterminism, independent-engine disagreement,
invalid null behavior, unstable localization, or failed H2 controls. Do not
retry automatically, delete evidence, loosen thresholds after inspection, or
replace an estimand without a new committed amendment.

## Completion criteria

MV17 is complete only when the project can state, with auditable evidence:

1. whether observed H0/H1 exceeds appropriate matched null structure;
2. which observations or features stably drive that structure;
3. whether H2 detects known voids under the implemented contract;
4. where real-data H2 is computationally and statistically stable; and
5. whether H2 adds information beyond H0/H1 and conventional baselines.

Publication value does not require a positive H2 result. A rigorously calibrated
negative result, coupled with interpretable H0/H1 localization and explicit
computational boundaries, is an acceptable and potentially useful contribution.
