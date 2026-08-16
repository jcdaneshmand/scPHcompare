# MV7-J claim, figure, and literature map prefreeze specification v1

Date: 2026-08-16

## Purpose and stop boundary

MV7-J converts the completed, independently validated MV5--MV7 evidence into
an auditable manuscript claim map, figure plan, and current-literature
positioning statement. It is an evidence-synthesis stage, not a new biological
analysis. It may derive complete descriptive summaries from already validated
artifacts, but it may not refit a model, select a favorable result, recompute
persistent homology, add a dataset, or promote a manuscript claim beyond its
source estimand.

Passing MV7-J authorizes corrected figure implementation and a read-only
external-dataset admission audit. It does not authorize external-data download,
new PH, confirmatory inference, causal interpretation, journal submission, or
replacement of the accepted landscape definition.

## Scientific hierarchy

Every proposed statement is assigned one of five statuses:

1. `supported_method`: validated implementation or mathematical contract;
2. `supported_descriptive`: complete existing-corpus descriptive result;
3. `conditional_context`: supported only under a named population,
   representation, or sensitivity condition;
4. `hypothesis_only`: scientifically motivated but not tested by the current
   estimand; or
5. `retire`: incompatible with the corrected method or evidence.

No result may be described as predictive, causal, biologically mechanistic,
externally generalizable, or superior unless a separately frozen estimand
directly supports that word.

## Landscape and observation contract

The manuscript and all figures must use the corrected primary definition:

- one biological sample is the unit compared and clustered;
- within each sample, the cell view treats cells as observations in the fixed
  cell coordinate space;
- within each sample, the gene view treats genes from the fixed global-core
  panel as observations in the fixed gene coordinate space;
- all finite positive-persistence intervals are retained;
- the essential H0 interval is excluded;
- all consecutive active landscape levels are used;
- H0 and H1 landscapes and squared-L2 distances are computed separately;
- integration is exact or error-controlled over the support;
- no universal grid count or landscape-level cap is used; and
- the within-view unweighted H0/H1 Euclidean composite is secondary and
  descriptive only.

The prior first-level, 100-point-grid, pointwise H0/H1 aggregation is a legacy
method and must not be presented as the current scPHcompare estimand.

## Frozen result synthesis

The following result summaries are scheduled and must be reported completely:

- per-view distribution of the median-across-seed H1 squared-distance
  contribution for all 7,626 sample pairs;
- H0 versus secondary H0/H1-composite partition agreement for every view,
  algorithm, and seed;
- PAM versus average-linkage partition agreement for every representation and
  seed;
- the complete selected-k and stability curve for all six representations;
- the complete 120-unit ARI/max-NMI outcome table by population, label axis,
  representation, algorithm, and metric; and
- a compact range synopsis that never replaces the complete tables.

The five seeds are technical realizations. Their variation is not biological
replication. H1 contribution is an energy decomposition of the secondary
Euclidean composite; it is not evidence that an H1 feature is biologically a
cell cycle, regulatory loop, or mechanism.

## Claim boundaries fixed before synthesis

- The corrected pipeline and its exact/error-controlled landscape distance are
  supportable method claims.
- The 124-sample results may be called coarse, seed-stable, and descriptive.
- Tissue and study alignment must be described together because their
  confounding is material.
- Canonical approach is a confounded proxy: all six full-corpus snRNA-seq
  samples are nested in substantia nigra and SRA850958; the primary 90 contains
  only scRNA-seq.
- Cell and gene views are complementary. Neither is globally superior.
- Equal-weight fusion remains a negative result, not a failed implementation.
- Existing integration evidence is corpus-specific and negative or null; PH
  measures geometry but does not remove batch effects.
- H0 and H1 stay separate. A combined result cannot replace either component.
- Algorithm sensitivity must be visible wherever clustering membership is
  interpreted.
- No current result supports external generalization, a causal technology
  effect, biological H1-cycle interpretation, or tissue prediction.

## Literature audit rules

The literature registry is frozen from primary sources verified through
2026-08-16. Each entry records publication state, DOI or stable URL, relevance,
and its effect on novelty or interpretation. At minimum it covers:

- persistence-landscape definition and stability;
- robust topological inference and outlier sensitivity;
- atlas-level single-cell integration benchmarking;
- contemporary topology-aware integration evaluation;
- single-cell methods using topology on cell and gene structures;
- topology-regularized single-cell dimensionality reduction; and
- current geometry-aware single-cell representation analysis.

The project must not claim that PH in single-cell analysis, topology-aware
integration evaluation, or dual cell/gene topology is unprecedented. Its
defensible contribution is the corrected, auditable, sample-comparison
workflow with explicit dual views, all-level landscapes, uncertainty,
confounding boundaries, negative results, and resource-validated execution.

## Figure rules

Each planned figure must identify its exact artifact, population, estimand,
uncertainty, claim status, and limitation. The main figures are planned as:

1. corrected pipeline and observation-level schematic;
2. cohort/accountability and confounding structure;
3. persistence-diagram to all-level H0/H1-landscape explanation;
4. H1 contribution and H0/composite relationship;
5. label-free stability and selected k;
6. complete descriptive tissue/study alignment;
7. PAM/average-linkage algorithm sensitivity; and
8. accepted prior integration and dual-view results, visibly separated from
   the full-124 descriptive analysis.

No old image may be reused if it encodes the first-level/100-grid method.
Figures must show complete prespecified families rather than a favorable
subset. The saved conceptual diagram may be regenerated only after its labels
match this specification and its data panels come from validated artifacts.

## Legacy manuscript and confidential-review policy

The dissertation, preprint PDF, and confidential review record are external
audit inputs and are not added to Git. Public evidence may paraphrase general
scientific risks, but may not quote, reproduce, identify, or hash-locate the
confidential reviewer text. The revised manuscript must explicitly correct:

- first-level/100-grid landscape language;
- label-chosen cluster counts;
- combined-only H0/H1 reporting;
- claims that PH removes batch effects;
- causal technology language;
- biological interpretation of loops without validation;
- claims of universal robustness or deformation invariance; and
- any implication that independently fitted PCA spaces are directly
  comparable.

## External-data decision gate

Existing data are sufficient for a methods-focused paper with a descriptive
case study, transparent negative findings, and bounded claims. An independent
external dataset is necessary before any of these named claims:

- performance on unseen studies or tissues;
- external reproducibility of cell/gene complementarity;
- general superiority to conventional or integration-specific baselines;
- robustness of the fixed 500-gene transductive panel outside this corpus; or
- general conclusions about how integration changes topology.

MV7-J may therefore authorize a read-only dataset and metadata audit, but not
new computation. Dataset admission must be prospective and based on label,
study, tissue, raw-data, gene-overlap, and resource suitability--never on a
favorable preliminary PH result.

## Authorship and credit

The author-team registry retains Jonah Daneshmand, Julia H. Chariker,
Akshitkumar Mistry, and Eric C. Rouchka. Dr. Mistry and Dr. Rouchka must remain
credited. Author order, corresponding-author responsibility, final CRediT
roles, and approval of the revised scientific narrative remain human author-
team decisions before manuscript or release submission.

## Acceptance criteria

MV7-J passes only if:

- all computational sources match their accepted hashes and decisions;
- the landscape contract is reproduced verbatim in machine-readable form;
- all scheduled result families are complete and independently recomputed;
- every claim has evidence, population, uncertainty, and limitation fields;
- every figure maps to complete validated artifacts;
- current literature entries have stable primary-source locators;
- confidential review content is absent from public outputs;
- no favorable selection, p-value, new PH, or new data is introduced;
- the external-data gate distinguishes methods sufficiency from
  generalization sufficiency; and
- all public outputs rebuild byte-for-byte and pass an independent validator.
