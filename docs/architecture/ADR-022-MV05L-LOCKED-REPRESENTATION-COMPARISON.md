# ADR-022: Compare SCT and integrated retrieval through locked paired difference-in-differences

- Status: Accepted and executed after pre-join freeze
- Date: 2026-08-10
- Decision scope: MV5-L only

## Context

MV5-E and MV5-K are independently accepted retrieval evaluations over the same
90 samples, 15 LOSO studies, five tissues, five seeds, and endpoint definitions.
Their marginal aggregate outcomes are already known. Directly subtracting H0
or H1 retrieval across representations would mix a change in topology with the
broader change in coordinate geometry caused by integration.

Each representation already has a matched same-coordinate energy baseline.
That permits a paired difference-in-differences estimand that asks whether
integration changes topology's incremental retrieval contrast relative to the
non-topological geometry available in the same representation.

## Decision

1. Bind the accepted MV5-E and MV5-K query endpoint, prediction-lock,
   validation, repeat, and artifact-manifest hashes before reading endpoints.
2. Disclose that marginal aggregate results were known; freeze the joint
   estimand before any sample-level cross-representation join.
3. Pair only exact common fold-study-seed-sample-tissue endpoints under a fixed
   five-family role map.
4. Use H0 and H1 topology-minus-energy representation
   difference-in-differences on macro MRR as the two-test primary family.
5. Keep fixed-1-NN analogues supportive and direct representation contrasts
   descriptive; never let either replace the primary estimands.
6. Require exact pseudobulk endpoint equality as an identity control.
7. Preserve seed-within-sample and equal-tissue aggregation, paired
   tissue-stratified study bootstrap, paired study sign flips, conservative
   numerical-boundary handling, and two-test Holm correction.
8. Prohibit reranking, recomputation, tuning, selection, tissue filtering, or
   modification of either accepted source evaluation.
9. Require independent reconstruction and a byte-identical clean repeat.

## Consequences

MV5-L can estimate whether integration changes topology's incremental
retrieval contrast beyond a matched within-representation geometric baseline.
It cannot recover full outcome blindness because the two marginal aggregate
results were necessarily known when the comparison was authorized.

A positive difference-in-differences means only that the topology-minus-energy
contrast became more favorable after integration; it does not establish that
topology beats energy in either representation. A negative value means the
relative topology contrast became less favorable. Direct representation
differences remain diagnostic because integration changes the underlying
coordinate space.

MV5-L does not authorize clustering, gene topology, fusion, new data,
optimization, package-default changes, or manuscript claim promotion.

## Execution evidence

The contract was committed at `b3f7e28` before endpoint pairing. All 2,250
paired rows and common denominators reconciled, and the 450 shared-pseudobulk
endpoints were exactly identical. H0 DID macro MRR was -0.001195 (95% paired
blocked-bootstrap CI -0.136724 to 0.167949; Holm p=1), and H1 DID was -0.032016
(-0.114822 to 0.052640; Holm p=1). All 11 independent reconstruction
categories and all 13 byte-repeat comparisons passed. The complete evidence is
`docs/audits/MV05L_LOCKED_REPRESENTATION_COMPARISON_2026-08-10.md`.
