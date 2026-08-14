# MV5-BD Rust full-equivalence and adoption gate v1

Date: 2026-08-13
Authorization: automatic safe continuation after accepted MV5-BC `c1de547`

## Purpose and fixed corpus

MV5-BD evaluates the accepted bounded Rust prototype against every accepted
MV5-AY dimension reference before any production-adoption decision. The fixed
private corpus is all 408 dimension rows from the 56 raw-hash-verified existing
diagrams: Tier D contains 318 exact R references and Tier E contains 90
adaptive-certified H1 references. No new diagram, seed, partition, label,
outcome, or biological claim is authorized.

## Scientific gates

- Tier D must pass 318/318 with squared-distance absolute error no larger than
  `1e-10 + 1e-10 * abs(reference)`.
- Tier E must pass 90/90 within the accepted R achieved absolute-error estimate
  plus `100 * machine epsilon * max(1, abs(reference))`.
- Every result must have status zero, engine version one, finite nonnegative
  squared distance, and finite-interval counts equal to the accepted manifest.
- Reversing all 408 pairs must produce bit-identical squared distance, swap the
  two input counts, and preserve levels/segments.
- All 112 diagram/dimension self-comparisons must be exactly zero.
- Two clean executions must yield byte-identical normalized scientific output.
- Whole-run peak RSS must remain at or below 1 GiB.

## Adoption boundary

Passing the numerical corpus can establish scientific kernel equivalence only.
Production adoption additionally requires an auditable build/distribution plan
for every supported R platform, artifact provenance and verification, installed
package fallback tests, and an explicit later default/API decision. This sprint
must not embed a binary, enable Rust by default, change cache identities, or
remove the canonical R path. If numerical equivalence passes but portable
distribution is not yet certified, the decision is scientific acceptance with
engineering adoption deferred.
