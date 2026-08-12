# MV5-AQ Numerical Persistence-Landscape Engine Remediation

Date: 2026-08-12

Branch: `codex/phase-0-audit-foundation`

Starting revision: `6d28da2`

Implementation revision: `ecc4957`

Validation-runner revision: `a8d8e89`

## Decision

MV5-AQ passes. The strict `1e-8` adaptive calculation that failed in MV5-AP
now returns a certified result and agrees with the exact critical-breakpoint
calculation. A representative high-H1-depth pressure pair also certifies within
the frozen 180-second bound.

This decision authorizes **only MV5-AP-R1**, a clean rerun of the previously
stopped realistic compatibility/resource gate. It does not authorize workflow
integration, a workflow default change, a legacy artifact rewrite, a project
data PH recomputation, or a manuscript claim.

## Frozen scientific contract

The remediation preserves the dissertation-aligned corrected definition:

- every finite persistence interval is retained;
- every active consecutive landscape level is used;
- H0 and H1 are calculated and reported separately;
- distances are squared-L2 integrals, with the descriptive combined value
  derived from the two dimension-specific values;
- exact critical-breakpoint integration is used where computationally
  appropriate;
- adaptive results must carry and pass an explicit global error certificate;
- no fixed grid, universal level cap, interval deletion, or silent tolerance
  relaxation is permitted;
- a failed certificate produces an error, not a numerical result.

The `auto`/500 policy selects a numerical engine; it does not change which
topological information enters the distance.

## Reproduced failure and root cause

The frozen MV5-AP sentinel contains 383 H0 intervals in each sample and
130/137 H1 intervals. The original adaptive method partitioned the integration
domain at feature endpoints but not at every landscape-order crossing. It then
distributed the global absolute budget in proportion to very narrow partition
widths. On ordinary failing partitions this demanded local absolute tolerances
as small as approximately `1.739e-12`, while local relative demands further
overconstrained a result whose global error could have been acceptable.

The error was therefore numerical budget allocation, not evidence that H1,
all active levels, or the dissertation-aligned definition should be removed.

## Repair

`landscape_reference_v2` adds
`adaptive_quadpack_partitioned_v2`. A deterministic midpoint pilot estimates
where work is needed and establishes one global allocation scale, after which
the requested absolute budget is divided equally across the deterministic
partitions. Pilot values never supply the reported integral.

The final certificate is deliberately conservative:

`achieved error = summed fine QUADPACK error + absolute coarse/fine delta`.

The result is accepted only when that value is below the caller's global
absolute/relative threshold. The public pair and matrix APIs now default to
`method = "auto"` with `exact_max_intervals = 500`. Counts at or below the
guard use exact integration; larger counts use the same-contract adaptive
engine. The accepted MV-04 corpus supporting this routing decision spans H0
counts 383--499 and H1 counts 79--2802.

## Frozen-sentinel results

| Quantity | Exact | Adaptive `1e-8` | Absolute difference |
|---|---:|---:|---:|
| H0 | 54.1735947138941 | 54.1735947138941 | 7.1054e-15 |
| H1 | 6.23387432164545 | 6.23387432146666 | 1.7879e-10 |
| Combined | 54.5310879525003 | 54.5310879524799 | 2.0435e-11 |

The repaired adaptive run completed in 21.994 seconds; exact completed in
5.307 seconds. For H1, the summed fine quadrature error was
`1.04225628852233e-8`, the independent refinement delta was
`2.72352451702318e-9`, and the conservative achieved estimate was
`1.31460874022465e-8`, below the relative-aware final threshold
`3.88611890558414e-7`. H0 and H1 both certified.

The `auto` sentinel result routed both dimensions to exact integration and had
the exact result's cache identity. Serialization/reload preserved the object,
distances, and cache key exactly.

## High-depth pressure result

The bounded pressure pair contains 499 H0 intervals in both samples and
1,206/1,471 H1 intervals. `auto` routed H0 to exact and H1 to adaptive. It
completed in 92.436 seconds, below the frozen 180-second limit, and returned:

- H0: 1.10205038028309;
- H1: 0.0661068137224299;
- combined: 1.1040313181711;
- H1 achieved error estimate: `1.86863835169085e-10`.

Both dimension-specific results certified.

## Verification

- Focused reference tests: 31 expectations passed.
- Focused public-API tests: 62 expectations passed.
- Full repository suite: passed, with only two expected CRAN-guard skips.
- Independent validation: 15/15 categories passed in each of two runs.
- Clean repeat: all scientific fields and certificates were identical; only
  measured wall time and serialization byte-size were excluded where stated.
- Legacy/workflow source blobs match MV5-AP completion `6d28da2`.
- All eight prohibited-change counters are zero.
- Final exact-staged package check is recorded in the completion commit.

## Public evidence

The tracked `docs/audits/mv05aq-*` tables contain the sentinel execution,
strict exact/adaptive agreement, error certificates, routing policy,
serialization check, pressure execution, root-cause findings, clean-repeat
ledger, source freeze, prohibited-change counters, independent validation, and
continuation decision. Sample paths and private serialized result objects
remain only under ignored `tmp/` and are not published.

## Next sprint

Run MV5-AP-R1 from the same frozen 24-diagram manifest. It must use the repaired
engine and routing policy, rerun the complete realistic compatibility/resource
evaluation, reproduce independently, and stop again if any strict certificate,
resource bound, serialization identity, or legacy/workflow immutability check
fails. Only MV5-AP-R1 may decide whether later opt-in integration planning can
resume.
