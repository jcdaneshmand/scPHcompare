# MV7-I descriptive outcome prefreeze

Date: 2026-08-16

## Decision

The prospective outcome design passes 15/15 builder checks and 15/15
independent checks. All 13 canonical evidence artifacts are byte-identical to a
clean repeat. This authorizes only the fixed MV7-I descriptive outcome
execution; no association was computed during prefreeze.

## Frozen evaluation

The immutable full-124 cluster assignments are evaluated without refitting for
six representations, five seeds, and two algorithms. PAM remains primary and
average linkage remains sensitivity. Five population/label endpoints are
scheduled:

- full 124: tissue, study, and canonical approach;
- primary-90 contextual restriction: tissue and study.

The primary-90 canonical-approach endpoint is explicitly retained as
`structurally_not_estimable_single_class` and is not queued because all 90
samples are scRNA-seq. The primary-90 restriction uses assignments fit on all
124 samples and therefore cannot replace the accepted blocked MV5-S analysis.

Each scheduled endpoint uses adjusted Rand index and max-normalized mutual
information. The queue contains 120 fixed summary units and expects 600 seed
rows. Every five-seed unit reports all values, mean, median, minimum, maximum,
and delete-one-seed jackknife technical SE. P-values, multiplicity adjustment,
winner selection, ranking, favorable-seed selection, algorithm pooling, and
representation pooling are prohibited.

## Confounding boundary

The six canonical snRNA-seq samples are all substantia nigra samples from study
SRA850958. Full-124 approach alignment is therefore a perfectly nested,
confounded descriptive proxy and cannot identify a technology effect. The
three added tissues are also each single-study and cannot support cross-study
generalization. Tissue and study alignment remain descriptive associations;
H1 alignment does not establish a biological cycle or mechanism.

## Validation and privacy

The independent validator reconstructed all metadata category counts, the 12
partition identities, all 120 queue IDs and their order, the approach nesting,
and the non-estimable endpoint. It verified exact input hashes, the public
schema firewall, all 15 prospective checks, and a 13/13 byte-identical repeat.
No cluster/metadata join, metric estimate, p-value, or method selection was
present in the prefreeze evidence.

Public execution outputs may contain only fixed identifiers, metric estimates,
uncertainty, counts, and status. Sample-level metadata joins and contingency
tables remain private; only their hashes and counts may be public.

## Durable evidence

- Specification:
  `docs/specifications/MV07I_DESCRIPTIVE_OUTCOME_PREFREEZE_SPECIFICATION_V1.md`
- Canonical evidence: `docs/audits/mv07i-outcome-prefreeze-evidence/`
- Independent validation: `docs/audits/mv07i-outcome-prefreeze-validation/`
- Prospective implementation: commit `f9659e1`
- Manifest-name validation correction: commit `3f5c225`

## Next authorized action

Implement the one-worker resource-capped descriptive outcome runner, execute
all 120 units without selection, produce a clean deterministic repeat, verify
immutable resume, and independently reconstruct all 600 seed metrics and 120
summaries before any result-dependent manuscript claim or figure map.
