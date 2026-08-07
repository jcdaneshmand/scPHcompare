# MV-05 statistical-plan freeze and design-feasibility audit

| Field | Value |
|---|---|
| Date | 2026-08-06 |
| Contract | `mv05_statistical_benchmark_contract_v1` |
| Scope | Design, feasibility, label boundaries, and executable validators only |
| Biological outcomes | Not computed |
| Pilot role | Technical smoke testing only |
| New data | Not triggered |
| Gate | Plan frozen; G-MV5 remains open |

## Outcome

MV-05 now has a prospective, executable statistical plan rather than a list of
candidate metrics. The plan freezes the sample as the inferential unit, study
as the validation and uncertainty block, five repeated cell subsamples as
technical realizations, and H0/H1 as separate results. It also freezes endpoint,
method, label-use, and multiplicity registries before biological outcomes are
calculated.

The design-only builder reproduced the existing-data inventory:

- large cohort: 124 samples, 18 studies, eight tissues, two approaches;
- cross-study candidate subset: 90 samples across 15 studies and five tissues;
- eligible tissues: bone marrow, PBMC, liver, testis, and colon;
- bone-marrow cohort: 25 samples in one study, suitable only for technical
  approach diagnostics; and
- MV-04 pilot: zero tissues eligible for cross-study confirmation.

It generated all 18 leave-one-study-out fold summaries with no study overlap.
These counts establish design feasibility only; no distance was associated with
a tissue, study, or approach outcome.

## Legacy-code disposition

The audit records why legacy evaluation cannot be reused as primary inference:

- class labels directly determine tissue/study/approach cluster counts;
- sample clusters are copied to cell rows, creating cell-count weighting and
  pseudoreplication;
- the random-null fraction can produce exact zero p-values and does not respect
  study blocks; and
- Ward and k-means are not generally valid for arbitrary sample-distance
  matrices, while the current spectral helper lacks frozen randomness and
  eigengap provenance.

PAM with label-free repeated-subsample stability is the primary clustering
contract. Average linkage is a sensitivity. Spectral clustering requires a new
deterministic implementation gate, and oracle `k` is historical-only.

## Integration boundary

`R/Integration_flexible.R` contains `FindIntegrationAnchors()` and
`IntegrateData()` but no `FindTransferAnchors()`, `MapQuery()`, or
`TransferData()` path. Current all-sample integrated objects are therefore
transductive and descriptive. Confirmatory integrated LOSO analysis must return
`held_out_mapping_unavailable` until a tested training-reference/held-out-query
contract exists; it may not silently substitute an all-sample fit.

## Statistical boundary

The primary biological endpoint is cross-study tissue retrieval macro mean
reciprocal rank. Fixed 1-NN balanced accuracy and study-pair-reduced distance
contrasts are supportive. Residual study and approach effects remain separate
technical endpoints. Integration effects are shown on biological and technical
axes, never collapsed into a weighted score.

The matched cell baseline is energy distance over the same 384 cells in the
same training-fit 30-PC coordinates. Pseudobulk is contextual, and RMS
Frobenius distance between matched gene-correlation matrices is the gene-view
baseline. Composition is conditional on an externally frozen harmonized cell
label vocabulary.

Primary intervals use 2,000 paired tissue-stratified study-block bootstrap
replicates. Optional randomization uses 9,999 paired study-level sign flips and
`(b+1)/(B+1)`. Holm correction is restricted to the three frozen families. No
inferential claim is allowed from fewer than four independent study blocks.

## Verification

- The design builder completed and emitted cohort, tissue, pilot, LOSO,
  integration, legacy-evaluation, and registry evidence tables.
- Focused MV-05 package tests passed 26/26 expectations, and the complete
  package suite passed 425/425 expectations.
- The final `R CMD build` completed in 453.3 seconds, and
  `R CMD check --no-manual`
  finished with `Status: OK` (zero errors, warnings, or notes).
- The contract cache key and input SHA-256 digests are recorded.
- `mv05-artifact-manifest-2026-08-06.csv` hashes the specification, code, tests,
  builder, and generated evidence.
- The private reviewer material and source PDFs remain ignored, and
  `example_run.r` remains preserved and untracked.

## Gate disposition

The statistical-plan-freeze subgate passes. G-MV5 does not pass because no
single-view biological benchmark has been executed and neither view has a
works/fails/uncertain disposition. The next authorized sprint is MV5-A/MV5-B:
implement fold-manifest and matched-baseline fixtures, then test inductive
integration on a bounded two-study technical case before any label opening.
