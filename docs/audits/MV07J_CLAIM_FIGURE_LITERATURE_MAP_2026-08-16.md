# MV7-J claim, figure, and literature map audit

Date: 2026-08-16

## Decision

MV7-J is complete. The evidence builder passes all 14 acceptance criteria, a
clean repeat reproduces all 18 output files byte-for-byte, and independent
validation passes 17/17 checks. This authorizes corrected figure implementation
and a read-only external-dataset admission audit. It does not authorize new
data download, new persistent homology, confirmatory claim promotion,
manuscript submission, or public release of confidential review material.

## Corrected manuscript estimand

The revised work compares biological samples. Each sample yields two separate
point-cloud views: cells are the observations in the cell view, and fixed
global-core genes are the observations in the gene view. Persistent homology is
computed within each sample and view; landscape distances then compare samples,
and the resulting sample-by-sample dissimilarities are clustered.

The primary persistence-landscape definition is now fixed throughout the
claim and figure map:

- retain every finite positive-persistence interval;
- exclude the essential H0 interval;
- use every consecutive active landscape level;
- calculate and report H0 and H1 separately;
- integrate squared-L2 distance exactly or with controlled numerical error;
- use no universal grid count or landscape-level cap; and
- treat the within-view unweighted H0/H1 Euclidean composite as secondary and
  descriptive only.

The legacy first-level, 100-point-grid, pointwise H0/H1 aggregation is retired.
It cannot appear as the current scPHcompare method.

## H0/H1 synthesis

All 7,626 sample pairs per view were summarized without selection. For the
median-across-seed fraction of the secondary composite's squared distance
contributed by H1:

- cell view: mean 0.004535, median 0.002053, 90th percentile 0.010744,
  maximum 0.119180; 10.73% of pairs exceed 0.01, 0.47% exceed 0.05, 0.026%
  exceed 0.10, and none exceed 0.50;
- gene view: mean 0.001271, median 0.000683, 90th percentile 0.003090,
  maximum 0.019597; 0.41% of pairs exceed 0.01 and none exceed 0.05.

The secondary composite is therefore strongly H0-dominated in this corpus.
That does not make H1 redundant: H1-only clustering can be distinct, and one
cell-view average-linkage seed changed substantially between H0 and the
composite (ARI 0.13846). PAM produced identical H0 and composite partitions in
all ten view-seed comparisons; gene average linkage was also identical in all
five seeds, while cell average linkage was identical in four of five.

The correct interpretation is consequently narrow: H1 should remain separate
and visible, but current sample-level H1 does not establish biological cycles,
pathways, or mechanisms.

## Clustering and descriptive outcomes

All six representations selected k=2 using the frozen label-free five-seed
stability rule. This is a coarse two-cluster description, not recovery or
prediction of eight tissues or 18 studies.

Algorithm sensitivity is material and representation dependent. PAM versus
average-linkage ARI ranges from essentially zero to 0.391 for cell H0 and is
near zero for cell H1. It is high for gene H0 (0.830--0.964) but low for gene
H1 (0.016--0.076). Neither algorithm may be selected as a favorable
replacement; PAM stays primary and average linkage stays a complete
sensitivity.

The prior MV7-I conclusion remains unchanged: tissue and study alignment is
modest and confounded; the six snRNA-seq samples are nested in substantia nigra
and SRA850958; and the primary-90 approach endpoint contains only scRNA-seq.
No technology effect, external generalization, or biological mechanism is
identified.

## Literature consequence

The current literature changes the novelty and robustness language materially:

- Bubenik's persistence-landscape work supports the functional, stable
  landscape representation, but not a blanket claim that a raw point-cloud
  Rips analysis is outlier robust.
- Chazal and colleagues show that an empirical distance-function persistence
  analysis can be destroyed by even a single outlier and motivate robust
  alternatives such as distance-to-a-measure. The revised paper must therefore
  retire claims of universal or inherent outlier robustness.
- The scIB benchmark separates batch removal from biological conservation and
  documents their trade-off. PH describes geometry; it does not itself correct
  batch effects.
- scGeom already analyzes topology on both cell and gene networks, so the
  project cannot claim that dual cell/gene topology is unprecedented.
- CONCORD now uses persistent homology and Betti curves alongside geometric and
  biological metrics to evaluate integration fidelity. Topology-aware
  integration evaluation is therefore not itself a sufficient novelty claim.
- TopoMetry and related current work reinforce the need to position PH as one
  component of a broader geometry-aware evaluation rather than as a complete
  replacement for conventional metrics.

The defensible contribution is the reproducible sample-comparison workflow:
explicit dual observation views, corrected all-level landscapes, separate H0
and H1, complete uncertainty and sensitivity reporting, negative fusion and
integration findings, confounding boundaries, and audited cross-engine and
resource-safe execution.

## Figure plan

The map defines eight main and three supplementary figures. The main sequence
is:

1. corrected sample-comparison and dual-view pipeline;
2. cohort accountability and confounding structure;
3. persistence diagrams to all-level H0/H1 landscapes;
4. H1 contribution and H0/composite concordance;
5. label-free stability and selected k;
6. complete descriptive tissue/study alignment;
7. PAM versus average-linkage sensitivity; and
8. prior blocked integration/dual-view evidence, visibly separated from the
   full-124 descriptive analysis.

The saved explanatory diagram should be regenerated now that the corrected
pipeline is complete. It must use validated artifacts and must not reproduce
the legacy first-level/100-grid method.

## External-data gate

The existing corpus is sufficient for the corrected method definition and a
methods-focused descriptive paper. It is not sufficient for four named claim
families:

- unseen-study or unseen-tissue generalization;
- external reproducibility of cell/gene complementarity;
- integration-method superiority; and
- inductive robustness of the fixed gene-panel strategy.

An external dataset is now the scientifically preferred expansion if any of
those claims is important to the resubmission. The next authorized step is
only a read-only dataset and metadata admission audit. Selection must precede
PH and use study structure, tissue replication, raw-data availability,
gene-overlap feasibility, label quality, and resource requirements--not a
favorable topology result.

## Credit and confidentiality

The author-team registry retains Jonah Daneshmand, Julia H. Chariker,
Akshitkumar Mistry, and Eric C. Rouchka. Dr. Mistry and Dr. Rouchka remain
explicitly protected in the credit plan. Final order, CRediT roles,
corresponding-author responsibility, and manuscript approval remain human
author-team decisions.

The public evidence contains only generalized scientific requirements. It
contains no reviewer quotation, confidential-review locator, or confidential
source text. The dissertation and manuscript PDFs remain external and are not
tracked in Git.

## Durable evidence

- Specification:
  `docs/specifications/MV07J_CLAIM_FIGURE_LITERATURE_MAP_PREFREEZE_V1.md`
- Complete public evidence:
  `docs/audits/mv07j-synthesis-evidence/`
- Claim map:
  `docs/audits/mv07j-synthesis-evidence/mv07j-claim-map.csv`
- Figure map:
  `docs/audits/mv07j-synthesis-evidence/mv07j-figure-map.csv`
- Literature registry:
  `docs/audits/mv07j-synthesis-evidence/mv07j-literature-registry.csv`
- External-data gate:
  `docs/audits/mv07j-synthesis-evidence/mv07j-external-data-gate.csv`
- Independent validation:
  `docs/audits/mv07j-synthesis-evidence/mv07j-independent-validation.csv`
- Prospective implementation commit: `f2deaf99375091b2030d00ef4f2cbd4be6861785`

## Next authorized action

Implement the corrected figures from the frozen map and, in parallel only at
the planning level, conduct a read-only external-dataset admission audit. Do
not download or compute on new data until a candidate and a named external
estimand pass a separate prospective gate.
