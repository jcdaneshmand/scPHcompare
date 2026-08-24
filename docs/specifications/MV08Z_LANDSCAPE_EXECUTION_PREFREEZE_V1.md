# MV8-Z production-landscape execution prefreeze v1

MV8-Z is the first execution contract after independent full-PH closure and
Rust-kernel admission. It does not change the landscape estimand.

## Scientific object

- Use every finite positive-persistence interval and exclude essential H0.
- Compute every consecutive active landscape level.
- Integrate exact, dimension-specific squared L2 over the complete support.
- Do not impose a universal grid or a landscape-level cap.
- Keep H0 and H1 as separate primary outputs. Any later combined summary is
  secondary and requires a separate comparison contract.

## Immutable execution queue

The 28-group queue contains 152,744 dimension-specific unordered pairs. A
private, label-blind axis orders each group's PH records by decreasing interval
burden and then diagram SHA-256. The exact unit and pair axes are hash-bound.
The queue is partitioned into 628 chunks of at most 250 pairs. Each chunk binds
its inclusive pair-ordinal range and the SHA-256 of the ordered pair identities.

Completed chunks publish private distances and a private status receipt
atomically. Resume is allowed only when both files, the execution head, admitted
candidate, pair subset, output hash, and ordered identities validate exactly.
Missing, partial, ambiguous, stale, or resource-stopped evidence is preserved
and fails closed; it is never deleted or retried automatically.

## Engine and fallback policy

The only sentinel engine is the explicit private Linux candidate admitted by
MV8-Y and verified by exact SHA-256 before execution. A Rust error or invalid
diagnostic stops the chunk. Canonical R exact integration is used for a
maximum-burden sentinel pair with at most 500 intervals per diagram; larger
pairs use the admitted adaptive route with an explicit error certificate. The
grouped Persim implementation is retained
as a separately gated recovery option; neither R nor Persim may silently mix
with Rust results inside a production chunk.

## Bounded sentinel

After the prefreeze commit, authorize only the first 250-pair chunk of the
label-blind group with the greatest top-two interval burden. Run it once into a
fresh primary root and once into a separate fresh repeat root, then run the
canonical-R exact or error-certified oracle on the maximum-burden pair. Require byte-identical
scientific results, oracle agreement, exact active-depth diagnostics, empty
unexpected stderr, resource-cap compliance, and independent rehash closure.

The full 152,744-pair execution, comparisons, clustering, fusion, labels,
outcomes, default adoption, and manuscript claims remain closed.
