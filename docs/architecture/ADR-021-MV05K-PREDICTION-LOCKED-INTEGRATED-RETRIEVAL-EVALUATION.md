# ADR-021: Evaluate integrated cell retrieval only after a new immutable prediction lock

- Status: Accepted and executed
- Date: 2026-08-10
- Decision scope: MV5-K only

## Context

MV5-J produced 176,750 label-closed rankings for five frozen integrated methods
across 75 fold-seed bundles. The accepted MV5-E SCT evaluation provides a
prospectively reusable statistical contract, but its rankings and observed
results cannot be consulted to alter the integrated evaluation. Opening labels
before a new MV5-J lock would permit outcome-informed reranking, scaling,
selection, or replacement.

The five subsample seeds are technical repetitions. Studies, not cells, pairs,
or seeds, are the independent generalization blocks.

## Decision

1. Require base commit `6836d1c`, the exact MV5-J ranking hash, 75 bundle
   identities, 375/375 completion rows, the five-method registry, the explicit
   scale disposition, and the frozen metadata hash to pass before label access.
2. Reuse the MV5-E endpoint and inference design without modification:
   cross-study tissue MRR as primary, fixed 1-NN balanced accuracy as
   supportive, seed-within-sample averaging, equal-tissue macro aggregation,
   tissue-stratified study-block bootstrap, two topology-minus-energy MRR sign
   flips, and Holm correction.
3. Keep integrated H0 and H1 separate. Retain the raw H0/H1 composite as
   descriptive and pseudobulk as contextual.
4. Use immutable neighbor ranks; prohibit refitting, rescaling, reranking,
   tuning, method selection, tissue selection, or replacement after labels
   open.
5. Do not read or compare the accepted SCT outcome during MV5-K execution or
   validation. Complete and lock the integrated result first.
6. Retain every null, negative, heterogeneous, tied, failed, and non-estimable
   disposition.
7. Require independent endpoint and inference reconstruction plus a
   byte-identical public repeat before accepting the result.

## Consequences

MV5-K can answer the same cross-study retrieval question for the integrated
view without changing the scientific target or borrowing information from the
SCT result. Any later SCT-versus-integrated comparison must use both already
locked results under its own contract.

MV5-K does not authorize post hoc method tuning, selective tissue reporting,
SCT comparison, clustering, gene topology, fusion, new data, or manuscript
claim promotion.

## Execution evidence

The prospective lock passed against base commit `6836d1c` and the frozen
MV5-J identities before labels opened. All 2,250 query-method-seed endpoints
were estimable. Integrated H0 minus matched integrated energy had macro-MRR
difference -0.005159 (95% blocked-bootstrap CI -0.140977 to 0.164341), and H1
minus energy had difference -0.040203 (-0.165125 to 0.062618); both Holm
adjusted p-values were 1. The accepted SCT outcome was not read.

Independent reconstruction passed all 11 validation categories, and a clean
repeat reproduced all 15 result artifacts byte for byte. The initial
randomization output was rejected when validation exposed a floating-point
boundary ambiguity; the final contract uniformly counts 64-epsilon boundary
ties as exceedances and was regenerated and revalidated. The complete record
is `docs/audits/MV05K_PREDICTION_LOCKED_INTEGRATED_RETRIEVAL_EVALUATION_2026-08-10.md`.
