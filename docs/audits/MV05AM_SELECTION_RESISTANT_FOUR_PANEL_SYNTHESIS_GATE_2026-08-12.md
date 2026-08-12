# MV5-AM selection-resistant four-panel synthesis gate

Date: 2026-08-12

Accepted base: `d5c4b7e`

Prospective contract: `9319b35`

Validator numeric-storage correction: `5b082b2`

## Result

MV5-AM completed the prespecified synthesis gate over the full PC20,
cosine-chord, nested-192, and nested-256 retrieval-robustness panels. It binds
all 96 macro estimands, all 96 blocked intervals, and all 16 primary contrasts.
No result, representation, dimension, endpoint, tissue, or seed was filtered.

The decision helper received zero numerical result rows. It used only the
immutable four-position sequence, panel completeness, shared-axis contracts,
clustering identifiability, the frozen landscape definition, and the
prerequisite-ordered evidence-gap registry. The four-position robustness
sequence is exhausted; no fifth configuration is authorized.

## Complete primary panel

The topology-increment MRR contrasts are:

| Contrast | Representation | Dimension | Estimate | 95% blocked interval | Raw p | Within-panel Holm p |
|---|---|---:|---:|---:|---:|---:|
| PC20 minus PC30 | Integrated | H0 | -0.02036 | [-0.08386, 0.02537] | 0.5893 | 1.0000 |
| PC20 minus PC30 | Integrated | H1 | -0.06444 | [-0.08203, -0.03374] | 0.0880 | 0.3520 |
| PC20 minus PC30 | SCT | H0 | -0.02567 | [-0.05116, -0.00042] | 0.1046 | 0.3520 |
| PC20 minus PC30 | SCT | H1 | 0.00537 | [-0.03537, 0.02673] | 0.8102 | 1.0000 |
| Cosine minus Euclidean | Integrated | H0 | -0.07216 | [-0.21895, 0.00552] | 0.2734 | 0.2734 |
| Cosine minus Euclidean | Integrated | H1 | -0.10106 | [-0.21066, -0.02869] | 0.0526 | 0.1052 |
| Cosine minus Euclidean | SCT | H0 | -0.11833 | [-0.20089, -0.06917] | 0.0070 | 0.0224 |
| Cosine minus Euclidean | SCT | H1 | -0.10088 | [-0.15610, -0.07495] | 0.0056 | 0.0224 |
| 192 minus 384 cells | Integrated | H0 | -0.01080 | [-0.03676, 0.05266] | 0.7511 | 1.0000 |
| 192 minus 384 cells | Integrated | H1 | -0.01115 | [-0.04821, 0.05495] | 0.7139 | 1.0000 |
| 192 minus 384 cells | SCT | H0 | 0.00131 | [-0.01601, 0.02958] | 0.9225 | 1.0000 |
| 192 minus 384 cells | SCT | H1 | -0.01366 | [-0.04051, 0.01803] | 0.5091 | 1.0000 |
| 256 minus 384 cells | Integrated | H0 | -0.00951 | [-0.03399, 0.06169] | 0.7198 | 0.7758 |
| 256 minus 384 cells | Integrated | H1 | -0.04579 | [-0.06398, -0.00569] | 0.0149 | 0.0596 |
| 256 minus 384 cells | SCT | H0 | 0.00808 | [-0.00340, 0.03157] | 0.3140 | 0.7758 |
| 256 minus 384 cells | SCT | H1 | -0.01386 | [-0.04671, 0.00354] | 0.2586 | 0.7758 |

These are four distinct intervention panels, each with its own prospectively
defined four-test family. They are not pooled into one effect and MV5-AM does
not introduce a post-hoc 16-test family. The complete table shows that geometry
has the clearest detected sensitivity in the current benchmark: cosine chord
reduces topology's MRR increment for SCT H0 and H1 after the prespecified
within-panel Holm adjustment. Cell-depth results are generally smaller; the
nested-256 integrated-H1 result remains suggestive but not confirmatory after
its frozen adjustment. Coordinate truncation shows negative interval shifts in
two cells of the table without adjusted detection. These observations are
descriptive synthesis, not configuration selection or default-setting evidence.

## Landscape interpretation

All four panels preserve the same corrected scientific definition: all finite
intervals, essential H0 excluded, all consecutive active levels, exact
critical-pair squared-L2 integration, H0/H1 separate, and the Euclidean H0/H1
composite descriptive only. No fixed grid, level cap, post-hoc normalization,
or level weighting enters the synthesis.

The four panels therefore support an important methodological conclusion: the
landscape calculation itself is held stable while coordinate number, geometry,
or cell realization changes. The observed cosine sensitivity is attributable
to the named geometry intervention, not to a changing landscape approximation.

## Validation and reproducibility

The source freeze contains 38 exact sources, including 12 complete numerical
tables. Independent validation rehashes every source and reconstructs all
96/96/16 rows, axes, estimates, intervals, p-values, and identities without
calling the production synthesis or decision helper. All 12 validation
categories pass. Two clean builds reproduce all 15 production/validation
ledgers byte-for-byte.

The first validation pass exposed only an integer-versus-double storage class
difference for MV5-AH Holm values equal to one. The independent validator was
corrected to demand exact numeric equality while ignoring integer/double
storage class, prospectively committed at `5b082b2`, and both accepted builds
were generated fresh. No scientific value or decision changed.

## Decision and next boundary

Complete within-training distance matrices remain absent for every panel, so
clustering is not identifiable from the accepted directed retrieval rows. Gene
topology, fusion, new data, Rust, optimization, and manuscript claims remain
closed.

The sole authorized next action is MV5-AN: a prospective public-landscape-
contract reconciliation prefreeze. It should map every current public and
internal landscape API, legacy grid/cap behavior, exact production engine,
tests, documentation, compatibility risk, and migration path. MV5-AM does not
authorize changing a default; that requires MV5-AN evidence and a separately
committed implementation sprint.
