# MV5-AA selection-resistant robustness continuation gate

Date: 2026-08-11

Accepted PC20 evidence base: `d4c332a`

Prospective decision-contract commit: `45c7685`

New robustness calculation or outcome access: zero

## Decision

MV5-AA authorizes exactly one later label-closed calculation:
`cells384_pc30_cosine_chord_v1`. It is the unique second configuration in the
canonical MV5-V execution order after completed PC20. The two nested-cell
configurations remain closed and cannot be substituted if cosine fails.

This gate does not execute cosine transformation, PH, landscapes, energy,
rankings, outcomes, or clustering. All 150 authorized queue rows retain false
labels-opened, rankings-computed, outcomes-computed, and execution-completed
states.

## Why continuation is scientifically justified

PC20 and cosine-chord address distinct named alternatives. PC20 asks whether
retaining 20 rather than 30 accepted training-fit coordinates changes the
result while Euclidean geometry remains fixed. Cosine-chord asks whether radial
coordinate magnitude drives the geometry while 384 cells and all 30 accepted
coordinates remain fixed. The complete PC20 result cannot resolve that radial-
scale question.

The gate binds all 24 PC20 estimands and intervals and all four primary tests.
It does not use any representation, H0/H1 component, tissue, endpoint, seed,
estimate sign, interval, raw p-value, or adjusted p-value to choose the next
axis. Heterogeneous PC20 evidence is reported completely, but cosine was next
before those values were evaluated.

## Exact later calculation scope

| Object | Count |
|---|---:|
| Atomic groups | 150 |
| LOSO folds | 15 |
| Seeds | 5 |
| Representations | 2 |
| Transformed views | 13,500 |
| Heldout-training biological pairs | 70,700 |
| H0/H1 landscape requests | 141,400 |
| Deterministic landscape subchunks | 720 |
| Matched-energy rows | 70,700 |
| Four-method rows | 282,800 |

Each row is inherited unchanged from the frozen MV5-V queue and differs from
the accepted reference only through row-unit normalization followed by
Euclidean chord distance. Zero or nonfinite cell-coordinate norms remain hard
failures.

## Resources and execution requirements

The accepted PC20 precedent completed 150 groups in 11,163.624 worker-seconds
(3.101 hours), with a 185.348-second maximum group and 638,365,696-byte peak
process-tree RSS. This is feasibility evidence, not a promise that cosine has
the same runtime.

The later cosine sprint retains one heavy worker, 600 seconds and 4 GiB RSS per
group, 8 worker-hours for the configuration, and 4 GiB new private storage. It
must bind all coordinate/source/runtime/implementation hashes, publish groups
atomically, preserve deterministic subchunks, pass typed-PH and H0-MST checks,
retain all-active exact H0/H1 landscapes and matched energy, clean-repeat
prospectively selected artifacts, prove an immutable full resume, and pass an
implementation-independent reconstruction.

## Validation and reproducibility

The decision contract passed 12 focused expectations, including adversarial
order drift, incomplete PC20 evidence, failed selection-firewall criteria, and
outcome-contaminated queue cases. A standalone validator that does not call
MV5-AA production scientific helpers passed 12 categories, independently
rehashing all 19 frozen sources, recovering the exact order, reconstructing the
scope and resource precedent, and verifying every closed-state field.

Two clean assemblies, each independently validated, reproduced all 11 ledgers
byte-for-byte. Twelve later-execution validation classes and ten hard-abort
classes are frozen.

## Boundary and next action

The next sprint may bind a cosine-only execution engine and perform the same
real-runner launch-readiness discipline used before PC20, or—if the already
accepted MV5-W runner is proven semantically identical for cosine—bind and
execute only the 150 frozen cosine groups under the MV5-AA caps. It must stop
before labels, rankings, outcomes, clustering, nested-cell configurations, or
scientific comparison. A separate prediction-locked outcome prefreeze remains
mandatory after a successful label-closed calculation.

Gene topology, fusion, new data, optimization/Rust, package defaults,
manuscript claims, PDFs/reviewer/private material, `example_run.r`, and pushing
remain outside this sprint.
