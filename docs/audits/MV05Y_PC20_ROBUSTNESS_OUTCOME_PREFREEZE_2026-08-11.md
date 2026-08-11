# MV5-Y PC20 robustness-outcome prefreeze audit

Date: 2026-08-11

Accepted calculation base: `f69c6e8`

Decision: approve a later, separately bound prediction-locked PC20 retrieval
robustness execution; do not execute outcomes in MV5-Y.

## What was inspected

The audit inspected the actual accepted MV5-X queue, all 150 private group
manifests and structural method-row axes, the accepted MV5-D5/J outcome-closed
ranking inputs, and the accepted MV5-E/K/L retrieval and MV5-R/S clustering
contracts. The dissertation, preprint, and confidential reviewer material were
consulted privately only for generalized design requirements. No PDF, reviewer
text, label file, private cache, or `example_run.r` was added to Git.

Accepted outcome-bearing MV5-E/K query files were hash-bound but not opened.
The external metadata SHA-256 was verified; only sample and study key columns
were parsed, while tissue and approach values were skipped. No label was joined
to a PC20 distance or prediction.

## Compatibility finding

All eight representation/family axes pair exactly:

| Representation | Family | Baseline rows | PC20 rows | Missing | Excess |
|---|---|---:|---:|---:|---:|
| SCT | H0 | 35,350 | 35,350 | 0 | 0 |
| SCT | H1 | 35,350 | 35,350 | 0 | 0 |
| SCT | raw H0/H1 composite | 35,350 | 35,350 | 0 | 0 |
| SCT | energy | 35,350 | 35,350 | 0 | 0 |
| Inductive integrated | H0 | 35,350 | 35,350 | 0 | 0 |
| Inductive integrated | H1 | 35,350 | 35,350 | 0 | 0 |
| Inductive integrated | raw H0/H1 composite | 35,350 | 35,350 | 0 | 0 |
| Inductive integrated | energy | 35,350 | 35,350 | 0 | 0 |

Equality was tested on representation, fold, seed, held-out query sample,
training reference sample, and mapped method. Distances and outcomes were not
compared. The resulting execution queue contains 150 groups, 70,700 biological
pairs, 282,800 method rows, 3,600 expected query/method results, and 7,200
expected query/method/endpoint rows. Every execution flag is false.

## Frozen analysis

The contract retains the accepted MRR primary endpoint and fixed-1-NN balanced
accuracy supportive endpoint. Five technical seeds are averaged within sample,
then samples within tissue, then the five tissues equally. Biological sample is
the inferential unit and held-out study is the resampling block.

The complete registry contains:

- 16 direct PC20-minus-PC30 sensitivity estimands: two representations by four
  method families by two endpoints; and
- 8 topology-increment difference-in-differences: two representations by H0/H1
  by two endpoints.

All 24 receive paired 2,000-replicate tissue-stratified held-out-study bootstrap
intervals. The only p-value family contains the four MRR H0/H1
difference-in-differences, using 9,999 paired study-block sign flips and Holm
adjustment. No p-values are authorized for the other 20 estimands.

Signed positive values mean the named PC20 quantity increased, negatives mean
it decreased, and zero means unchanged. No equivalence or noninferiority claim
is authorized because no defensible practical margin has been frozen.

## Clustering disposition

PC20 clustering is not identifiable from MV5-X. The accepted artifacts contain
only directed held-out-query to training-reference distances and zero
within-training PC20 pairs. MV5-R/S clustering requires complete training
distance matrices, frozen k selection, training partitions, and held-out
assignments. The missing label-closed scope is 262,675 within-training
biological pairs per representation, 525,350 total before method/component
expansion.

No imputation, symmetrization of incomplete rows, or reuse of incompatible
30-coordinate cluster artifacts is permitted. A future PC20 clustering
sensitivity would need its own label-closed calculation and later prefreeze.

## Audit artifacts and guards

Fourteen public ledgers bind 178 source identities, the method map, two
endpoints, 24 estimands, exact axis checks, structural label join, 150-group
queue, clustering disposition, 14 independent validations, 12 abort rules,
resource/repeat/resume controls, complete-reporting rules, and the final
decision. A clean second prefreeze assembly reproduced all 13 generated
contract ledgers byte for byte.

A later runner is limited to one worker, 300 seconds and 4 GiB per group, two
aggregate worker-hours, and 1 GiB public output. Four prospectively selected
private groups, a clean full public repeat, a full immutable resume, and an
implementation-independent reconstruction are mandatory.

## Stop state

- PC20 ranks computed: 0
- labels joined to PC20 predictions: 0
- PC20 retrieval outcomes computed: 0
- PC20 clustering outcomes computed: 0
- method selection performed: 0
- other robustness configurations authorized: 0
- execution authorized by MV5-Y: 0

The next eligible sprint is MV5-Z: implement, bind, and execute only this exact
prediction-locked PC20 retrieval robustness contract, with all 24 estimands and
complete validation. It must not calculate PC20 clustering or begin another
robustness configuration.
