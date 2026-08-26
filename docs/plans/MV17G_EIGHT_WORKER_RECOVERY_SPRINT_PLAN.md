# MV17-G eight-worker recovery sprint plan

Date: 2026-08-26
Status: 24/24 recovery-prefreeze gates passed at implementation head `7b95335`; parallel execution remains closed until that audit commits and pushes

## Purpose

Resume the prospectively frozen MV17-G primary and repeat calibration with eight
single-threaded children while preserving every completed serial child exactly.
This changes execution concurrency only. It does not change a source matrix,
unit ordering, null model, seed, replicate count, PH dimension, summary,
probability formula, privacy rule, or downstream scientific firewall.

## Capacity basis

- 40 logical CPUs were visible in WSL.
- The observed load was approximately 7.7 at the decision checkpoint.
- Approximately 1.4 TiB RAM was visible, with about 23 GiB used.
- Completed MV17-G children used approximately 160 MiB peak RSS.
- Eight observed children therefore project to roughly 1--2 GiB aggregate RSS.
- The conservative hard boundary remains eight times the existing 8-GiB
  per-child cap, or 64 GiB, with one thread per child.

## Immutable scientific contract

- Primary: 1,188 grouped children and 91,740 scientific runs.
- Repeat: 27 grouped children and 2,085 scientific runs.
- Ninety-nine fixed replicates per compatible unit/null-family group.
- H0 and H1 remain separate; four summaries remain separate.
- Existing seed registry, private locator, queues, matrices, and grouped worker
  remain unchanged.
- Completed serial children are adopted only after exact byte/hash, empty-stream,
  payload, MST-oracle, and resource validation. They are never recomputed.
- Eight workers, one thread per child, zero retries, and atomic per-child outputs.
- Public outputs remain aggregate-only.
- Real localization, labels, outcomes, clustering, fusion, biology, and manuscript
  claims remain closed.

## Prospective transition sequence

1. Commit the isolated recovery helper, prefreeze builder, parallel controller,
   parallel-aware closure, tests, this plan, and the PROJECT_PLAN decision.
2. Verify focused and complete tests without touching the live serial process.
3. Stop the serial controller at an atomic boundary:
   pause only the R parent while its current child finishes, require no active
   child and a consecutive complete prefix, then send signal 9 to that still-
   paused R parent only. Do not continue the R parent. Require the outer
   GNU-time receipt plus empty outer streams.
4. Build a fresh private/public recovery prefreeze at the exact implementation
   head. Privately bind every prefix artifact and all 264 matrix hashes; publish
   only counts, aggregate resources, capacity, implementation hashes, and
   firewalls.
5. Commit and push the recovery prefreeze before launching parallel work.
6. Resume the first pending job through the end of the primary in waves of at
   most eight distinct queue orders. Require one-thread environment variables,
   zero retries, atomic outputs, and ordered parent-side validation.
7. Validate the first complete eight-child wave before allowing continued waves.
8. After 1,188/1,188 primary children close, run the frozen 27-child repeat with
   the same controller and eight-worker contract.
9. Build the parallel-aware immutable MV17-G closure, run focused and complete
   tests, review privacy and git state, then commit and push.

## Stop and recovery rules

- Never run the serial and parallel controllers concurrently.
- Never stop or signal an active grouped child.
- Signal 15 followed by `CONT` is forbidden. A disposable R-process test showed
  that this sequence can execute an additional post-child R expression before
  termination, so it cannot guarantee that another queue item will not launch.
- Parent-only signal 9 after the active child exits while the parent remains
  stopped is required; it permits no further R code to execute and preserves
  the child's already atomically promoted four-artifact result.
- Any partial four-artifact child set is a terminal stop; preserve it and do not
  retry without a new committed recovery prefreeze.
- Any duplicate queue order, nonempty stream, nonzero child exit, hash drift,
  non-finite metric, H0 MST-oracle failure, timeout, RSS breach, aggregate cap,
  or public privacy breach is terminal.
- Completion order may differ within a wave, but all ledgers and scientific
  metrics are reconstructed in frozen queue order.
- A controller interruption may adopt only complete validated parallel waves;
  recovery after any incomplete wave requires a new prospective record.

## Acceptance evidence

- Exact committed implementation head and exact committed recovery-prefreeze head.
- Consecutive serial prefix with four artifacts per child and no partials.
- Signal-9 serial controller receipt and empty outer streams.
- Eight-worker/one-thread/zero-retry status for primary and repeat.
- Exact primary/repeat queue and scientific-run cardinalities.
- Exact maximum-burden repeat agreement.
- Independent reconstruction of 7,392 private empirical rows and 56 public
  aggregate calibration strata.
- Resource inequalities appropriate for parallel execution rather than the
  obsolete serial wall-time sum rule.
- Complete privacy, source-binding, focused-test, and full-suite closure.
