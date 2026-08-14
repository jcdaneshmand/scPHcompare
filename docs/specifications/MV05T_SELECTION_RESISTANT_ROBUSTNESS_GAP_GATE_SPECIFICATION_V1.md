# MV5-T selection-resistant robustness/gap gate specification v1

Date frozen: 2026-08-10

Accepted outcome base: `64d84d9`

New robustness outcomes: not computed

## Purpose and boundary

MV5-T converts the already planned robustness work into a bounded,
selection-resistant sensitivity program after completion of the retrieval and
clustering benchmarks. MV5-S outcomes establish that robustness remains useful;
their values do not choose a representation, distance, algorithm, fold, cell
count, coordinate count, or metric.

This gate may inspect artifact identities, coordinate schemas, measured
resources, and previously completed public outcomes. It does not compute a new
diagram, distance, retrieval endpoint, cluster, or biological/technical
outcome. Any later robustness result remains secondary sensitivity evidence and
cannot rescue, replace, or promote a primary/null result.

## Already covered axes

The following do not justify new execution:

- five fixed cell-subsample seeds are already complete;
- H0, H1, and the unweighted raw H0/H1 combination are already separate;
- exact all-active-level landscapes, matched energy, and pseudobulk context are
  already reported;
- SCT and inductive integrated representations are already paired;
- PAM primary and average-linkage sensitivity are already complete;
- label-free `k` selection is complete, while oracle `k` remains prohibited.

Repeating these axes or adding favorable seeds/algorithms after MV5-S would be
post-outcome expansion rather than robustness.

## Admitted perturbation families

Only three outcome-independent families enter a resource admission:

1. Nested cell count: 192 and 256 cells versus the accepted 384-cell baseline,
   with 30 PCs and Euclidean geometry.
2. Coordinate dimension: the first 20 coordinates versus the accepted first 30,
   with all 384 cells and Euclidean geometry.
3. Point geometry: row-unit-normalized 30-PC coordinates with Euclidean chord
   distance versus accepted Euclidean geometry, with all 384 cells.

These are one-factor-at-a-time settings. No factorial interaction grid is
authorized. The 192/256 subsets are nested within each accepted 384-cell view.
Order cells by SHA-256 of the contract ID, sample ID, technical seed, and cell
ID; use the same order for both representations. Cell and coordinate identities
must agree across paired SCT/integrated cache views before admission.

The 20-PC analysis truncates the existing training-fit decomposition and does
not refit PCA. Cosine sensitivity uses Euclidean chord distance after
unit-normalizing each nonzero cell vector. A zero-norm cell is a hard failure.

## Deferred or rejected families

- 250/1,000-feature panels require new training-only panel/PCA/integration
  fits and remain deferred.
- 50 PCs are unavailable in the accepted 30-PC caches and remain blocked until
  a separately justified refit.
- additional technical seeds are redundant with the completed five-seed plan.
- spectral, oracle-`k`, or broad clustering searches are prohibited or deferred;
  they cannot be introduced to improve MV5-S outcomes.
- small-study exclusion requires a complete fold rerun, not post-hoc removal
  from the held-out aggregate, and remains deferred.
- additional integration methods require native inductive contracts and a
  separate method-expansion gate.

## Scientific estimand

The eventual estimand is change in the already frozen cross-study retrieval
effects under each perturbation, using held-out biological samples blocked by
study. Reuse the MV5-E/MV5-K endpoint, aggregation, uncertainty, and complete-
reporting rules. Recompute H0/H1 landscapes and matched energy in the perturbed
coordinate space. The raw H0/H1 combination remains descriptive.

The next admission does not calculate that estimand. It validates only
coordinate construction, PH/landscape correctness, exact source identity,
runtime, memory, storage, deterministic repeat, and resume.

## Bounded admission queue

Use one seed (`20260805`) and three outcome-independent folds:

- minimum training size: `large_loso_v1:SRA779509`, n=65;
- median training size with lexical tie break: `large_loso_v1:SRA703206`, n=86;
- maximum training size with lexical tie break: `large_loso_v1:SRA713577`, n=89.

Cross three folds, two representations, and four configurations for exactly 24
admission groups. Labels remain closed. Every group must construct all 90 views,
run full Euclidean VR H0/H1, stage exact all-active-level landscapes, and
exercise the matched-energy calculation. No retrieval outcome is calculated.

## Resource rationale

Accepted baseline measurements imply approximately 7,720.432 worker-seconds
for SCT plus integrated PH and 6,267.045 worker-seconds for their landscape
distance stages per complete setting: 13,987.477 seconds (3.885 hours). Four
full settings therefore project to 15.54 worker-hours before assembly and
validation.

Observed diagram/output/staging evidence conservatively projects above 10 GB
when repeated four times. Full execution is therefore not authorized by this
gate. The 24-group admission is capped at one worker, two worker-hours, 600
seconds and 4 GiB RSS per group, and 2 GiB of new private storage. It must
produce a measured streaming/staging plan before a later full queue can pass.

## Validation and abort rules

The admission must validate:

- accepted D1 and G manifest hashes and all 75 private files per representation;
- exact paired 384-by-30 shapes and cell IDs in the three selected fold/seed
  groups before any perturbation;
- nested 192/256 membership and representation-neutral cell order;
- first-20-coordinate identity;
- finite row normalization and Euclidean chord invariants;
- independent H0 MST and small exact H1/landscape oracles;
- atomic artifact/status pairs, deterministic repeat, immutable resume, and
  complete failure reporting.

Abort on source/hash/shape drift, cell-ID mismatch, missing coordinate, zero
norm, leakage, non-nested subsets, PH/landscape oracle failure, partial/stale
artifact, cap breach, or any attempt to open labels or calculate outcomes.

## Stop boundary

MV5-T stops at a committed no-outcome gate and 24-unit admission queue. The
next sprint may execute only that resource admission. Full robustness outcomes,
method ranking, spectral promotion, gene topology, fusion, new data,
optimization/Rust, package-default changes, manuscript claims, PDFs/reviewer
material/private-cache tracking, `example_run.r`, and pushing remain closed.
