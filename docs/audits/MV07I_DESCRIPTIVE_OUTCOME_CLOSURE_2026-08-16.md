# MV7-I descriptive outcome closure

Date: 2026-08-16

## Decision

MV7-I descriptive outcome execution is complete and independently passes
15/15 checks. All six private artifacts are byte-identical across production
and clean-repeat roots; immutable resume preserves every production hash, byte
size, and modification time. Complete public seed and summary tables are now
available. This authorizes MV7-J claim-map and figure planning only. It does not
yet authorize manuscript claims.

## Execution result

Production completed in 4.206 seconds with 79,429,632 bytes peak process-tree
RSS. The repeat completed in 5.779 seconds with 79,376,384 bytes peak RSS. Both
were far below the frozen 900-second and 2-GiB caps. The run produced all 600
seed metrics, all 120 five-seed summaries, 5,620 private contingency cells,
and the one prespecified non-estimable endpoint. No p-value, ranking, tuning,
or method selection was computed.

The independent validator reconstructed the 124-row canonical metadata join,
all ARI and max-NMI values from contingency definitions, every five-seed
summary, and every private contingency cell. It also reconfirmed that all six
snRNA-seq samples are nested in substantia nigra and SRA850958 and that the
primary-90 approach endpoint contains only one class.

## Complete-pattern synopsis

All six representations used the label-free selected `k=2`. The cluster
alignment metrics are therefore coarse two-cluster descriptions against eight
tissues, 18 studies, or two nested approach classes; they are not classifiers
and were assigned no inferential p-values.

For primary PAM partitions in the full 124 samples, across the four separate
cell/gene H0/H1 representations:

- tissue mean ARI ranges from 0.0476 to 0.0927 and mean max-NMI from 0.1096 to
  0.1514;
- study mean ARI ranges from 0.0552 to 0.0868 and mean max-NMI from 0.1403 to
  0.1484; and
- canonical-approach mean ARI ranges from -0.0557 to 0.0475 and mean max-NMI
  from 0.0242 to 0.0654.

For the contextual primary-90 restriction of those same full-124 PAM
partitions:

- tissue mean ARI ranges from 0.0838 to 0.1101 and mean max-NMI from 0.1180 to
  0.1439; and
- study mean ARI ranges from 0.0506 to 0.1044 and mean max-NMI from 0.1362 to
  0.1610.

These are modest descriptive alignments, not evidence of tissue prediction or
technology causation. Tissue and study values are similar enough that their
confounding remains material. The full-124 approach values are especially
non-identifiable as technology effects because the six snRNA-seq samples are
exactly the substantia-nigra/SRA850958 group. The 90-sample approach result is
correctly absent rather than replaced with a numeric metric.

Average linkage supplies the prespecified algorithm sensitivity and sometimes
differs substantially from PAM. Consequently, any later figure or claim must
show PAM as primary and average linkage as sensitivity; it may not present one
algorithm as a favorable replacement. The public tables retain every
representation, endpoint, algorithm, metric, seed, and uncertainty summary.

## Interpretation boundary

The results support only the statement that the corrected dual-view topology
contains coarse, seed-stable structure with modest descriptive alignment to
tissue and study labels, whose magnitude and membership can depend on the
clustering algorithm and topology component. They do not establish
cross-study prediction, a causal batch or assay effect, a biological H1 cycle,
external generalization, or superiority of a view, homology dimension,
composite, or algorithm.

## Durable evidence

- Outcome prefreeze:
  `docs/audits/MV07I_DESCRIPTIVE_OUTCOME_PREFREEZE_2026-08-16.md`
- Complete public results and validation:
  `docs/audits/mv07i-outcome-validation/`
- Outcome runner: commit `9cc3357`
- Independent validator: commit `1b199ae`
- Contingency serialization correction: commit `c79360d`

## Next authorized action

Build an auditable MV7-J claim and figure map that uses the complete validated
tables, quantifies H0/H1 contribution and algorithm sensitivity without
selection, links every proposed statement to its population and limitation,
and separates corrected new results from the prior blocked MV5-S evidence.
Literature and reference updates must be searched prospectively and cannot be
used to retrofit the computational result.
