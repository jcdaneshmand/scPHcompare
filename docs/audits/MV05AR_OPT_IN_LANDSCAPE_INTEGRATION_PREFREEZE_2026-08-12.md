# MV5-AR opt-in corrected-landscape integration prefreeze

Date: 2026-08-12

Authorization: MV5-AP-R1 completion `f6df036`

Prospective contract: `aac9ef9`

Validator correction: `1da7fe5`

## Outcome

MV5-AR passes and authorizes only MV5-AS: implementation of an additive,
default-off corrected-landscape artifact stream. It does not authorize
corrected-matrix consumption by clustering, visualization, Betti analysis, or
cross-iteration analysis; a workflow default change; a legacy artifact rewrite;
project-data calculation; or optimization.

## Why the boundary is additive

The actual workflow currently has two incompatible concepts behind the word
“landscape.” `process_iteration_calculate_matrices()` writes the historical
level-1, 100-point unit-grid landscape list and one combined matrix. Clustering
and visualization consume that single matrix, while modular and cross-iteration
analysis consume historical grid curves.

The corrected dissertation-aligned API instead produces separate all-active
H0 and H1 squared-L2 matrices plus a descriptive combined matrix. Silently
placing that object or its combined matrix into the legacy field would erase
the H0/H1 distinction and could make legacy curve consumers misinterpret a
distance object. MV5-AR therefore freezes corrected calculation as a new
sidecar artifact tree only.

## Caller boundary

A later implementation may add `corrected_landscape_control = NULL` to
`run_postprocessing_pipeline()` and pass it through
`run_unified_pipeline()`. `NULL` performs no corrected work and preserves the
historical workflow. The other four exported/internal exposure points remain
unchanged.

An enabled v1 control is artifacts-only and binds explicit wall/pair budgets,
at least 1.5-GiB RSS capacity, one worker, resume-or-fail behavior, scientific
`auto` routing, exact guard 500, and strict absolute/relative `1e-8`
tolerances. Scientific fields cannot be loosened through `...`.

## Versioned artifact and resume design

Each iteration/input-set hash receives a non-colliding directory under
`results/corrected_landscape_v1/`. Seven public artifact classes are frozen:

1. resource plan;
2. input manifest;
3. one versioned RDS pair shard per canonical pair;
4. pair index;
5. versioned H0/H1/combined matrix object;
6. provenance table; and
7. hash-bound completion marker written last.

Every write is create-only and atomic. Existing pair shards are reused only
after class, input identity, cache key, contract, certificate, size, and hash
verification. A mismatch hard-fails. The final matrix is reconstructed from
and cross-checked against all shards before atomic publication. Interruption
therefore leaves verified resumable work but no apparently complete result.

## Resource admission

MV5-AP-R1's accepted maximum unit was 567.94 wall seconds and peak RSS was
990,363,648 bytes. The frozen v1 planner uses one worker and a conservative
planning model:

- 30 seconds per exact/exact pair;
- 240 seconds per pair with any adaptive dimension;
- 30 seconds per iteration for startup/finalization;
- minimum caller RSS budget 1.5 GiB; default 2 GiB;
- explicit caller `max_pairs` and `max_wall_seconds` ceilings.

H0/H1 interval counts outside the observed 383--499 / 79--2802 envelope
require profiling and are refused. Pair or time estimates beyond the caller
budget are also refused. These are execution admission rules: they never cap
intervals, levels, or filtration range and never loosen tolerance.

## Legacy coexistence and rollback

All legacy functions, filenames, iteration fields, custom override keys,
clustering inputs, Betti summaries, cross-iteration curves, and defaults remain
unchanged. Corrected artifacts never reuse a legacy filename and cannot be
silently converted. Rollback is disabling the explicit control; legacy outputs
remain available and incomplete corrected work is retained for audited resume
or explicit owner disposition, never automatically deleted.

## Validation

- Six exported workflow exposure points were rescanned with exact formals.
- Four workflow sources and `NAMESPACE` equal MV5-AP-R1 completion `f6df036`.
- Twelve control fields, seven artifact schemas, ten atomic/resume steps,
  fifteen required implementation validations, fourteen abort rules, and five
  migration stages are frozen.
- Fourteen public prefreeze ledgers reproduce byte-for-byte across clean builds.
- Sixteen independent validation categories pass twice; validator outputs are
  byte-identical.
- Fourteen prohibited-change counters are zero.
- No artifact or project calculation was performed.

## Next sprint

MV5-AS may implement only the additive default-off artifact producer and its
tests. It must demonstrate default-off equivalence, analytic oracles,
pair-shard/matrix equivalence, strict certificates, atomic interruption/resume,
immutable completed resume, cache identity, serialization, legacy coexistence,
resource admission/rejection, and a bounded synthetic/realistic smoke.

Downstream use of H0, H1, or descriptive combined matrices requires its own
later prefreeze. Optimization, including Rust, also remains separate so that
performance work cannot change scientific semantics.
