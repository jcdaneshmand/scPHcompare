# MV5-AE selection-resistant nested-cell continuation gate audit

Date: 2026-08-11

Accepted base: `0b32d76`

Status: complete; nested-192 calculation eligible but not executed

## Decision

MV5-AE authorizes exactly one later label-closed calculation:
`nested_cells_192_pc30_euclidean_v1`, 150 groups. Nested 256 remains closed.

No nested-192 view, persistence diagram, landscape, energy distance, ranking,
label join, endpoint, or outcome was calculated in this sprint.

## Why continuation is scientifically distinct

The completed PC20 analysis changed coordinate count at 384 cells under
Euclidean geometry. The completed cosine analysis changed point geometry at
384 cells and 30 coordinates. Neither answers how sensitive the topological
sample comparison is to the number of cells used to represent each sample.

Nested 192 holds 30 coordinates and Euclidean geometry fixed and takes exactly
the first 192 cells from each frozen deterministic 384-cell realization. It is
therefore a one-factor-at-a-time cell-subsampling sensitivity rather than an
unrelated resampling experiment.

## Selection firewall

Both prior analyses were bound completely and unsliced: 24 estimands, 24
intervals, and four primary tests each. The decision helper received only the
analysis identities and complete-panel states. It received no representation,
homology dimension, tissue, endpoint, seed, estimate sign/magnitude, interval
exclusion, p-value, or method ranking.

The canonical MV5-V order fixes nested 192 at position three after completed
PC20 and cosine, with nested 256 at position four. The prior outcomes therefore
cannot reorder the remaining configurations.

## Authorized scope

| Item | Frozen value |
|---|---:|
| Groups | 150 |
| Folds | 15 |
| Seeds | 5 |
| Representations | 2 |
| Cells per view | 192 |
| Coordinates | 30 |
| Views | 13,500 |
| Held-out/training biological pairs | 70,700 |
| H0/H1 landscape rows | 141,400 |
| Landscape subchunks | 720 |
| Energy rows | 70,700 |
| Four-method rows | 282,800 |

H0 and H1 remain separate exact all-active-level landscape distances. The raw
H0/H1 composite is descriptive; energy remains the matched same-coordinate
baseline. The later calculation must stop before ranking or labels.

## Resources

Six real nested-192 admissions completed 540 views in 153.880 seconds total,
with 27.052 seconds maximum per group and 419,532,800 bytes peak RSS. The
complete PC20 and cosine precedents required 3.101 and 2.608 worker-hours.

The conservative later envelope is one worker, 600 seconds and 4 GiB RSS per
group, six aggregate worker-hours, and 4 GiB private storage. These are
admission limits, not a runtime promise.

## Validation and reproducibility

Ten helper-independent categories pass:

- 30/30 source hashes;
- exact canonical configuration order;
- both complete unsliced evidence panels;
- nine mandatory criteria and nine prohibited selection inputs;
- the nested-192-only decision;
- all 150 queue axes;
- exact scientific cardinalities;
- real-admission and full-precedent resource evidence;
- 12 later validation requirements and 10 abort rules; and
- zero new calculation, ranking, label, or outcome activity.

Two clean assemblies reproduced all 11 generated gate ledgers byte-for-byte.
The public repeat and independent-validation ledgers record this evidence.

## Next action and boundary

The next sprint should bind the existing streamed runner/runtime to only these
150 nested-192 groups, execute the complete label-closed calculation, and
independently reconstruct nested inclusion, PH/MST, landscapes, energy, and all
method rows. It must stop before any ranking, tissue access, outcome, clustering,
or nested-256 calculation.

This gate does not imply that 192 cells will be robust, better, or worse. It
does not authorize equivalence, a new default, gene/fusion/new-data/Rust work,
manuscript claims, private/PDF/reviewer tracking, pushing, or `example_run.r`.
