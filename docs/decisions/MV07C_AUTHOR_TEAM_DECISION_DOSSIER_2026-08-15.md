# MV7-C author-team decision dossier

Date: 2026-08-15

Status: owner and author-team decision required

## Recommended direction

Prepare a methods-focused existing-data manuscript and figures before
authorizing more PH. The project now has a defensible methodological core: an
orientation-safe cell/gene framework, the dissertation-aligned corrected
landscape definition, exact or certified integration with all active levels,
study-blocked evaluation, multiple cell-view robustness panels, a complete
gene-view calculation, and a prospectively negative fusion result.

That is publishable scientific work even though it does not produce a simple
winner. The strongest honest story is that cell and gene topologies expose
different aspects of sample geometry; equal-weight fusion does not reliably
improve retrieval; and apparent relative performance depends materially on
study and tissue. This is more informative than presenting one clustering as
universally correct.

## What the current evidence supports

- The corrected landscape definition is stable and should remain unchanged:
  finite positive-persistence intervals, essential H0 excluded, all
  consecutive active levels, separate H0/H1, exact or error-controlled
  squared-L2, no universal grid, and no universal level cap.
- The software and computation are reproducible, resumable, independently
  validated, and sufficiently fast for the accepted 90-sample analysis.
- Cell-view robustness has been tested across 192/256/384 cells, PC20/PC30,
  Euclidean/cosine-chord geometry, five seeds, H0/H1, SCT/integrated
  representations, and two secondary clustering algorithms.
- The gene view is technically feasible and distinct, but conditional on a
  fixed transductive 500-gene panel, one gene metric, SCT, and 384 cells.
- Equal-weight fusion failed its required improvement over both component
  views. The descriptively higher 0.25 gene weight must not be promoted.
- Study/tissue influence is material. The cell-minus-gene composite MRR order
  reverses when one study (`SRA716608`) or colon is removed.
- Current retained-cell and mixed-study approach diagnostics do not trigger
  their fixed alerts, but their uncertainty is too broad to prove absence of
  confounding.

## What must not be claimed

Do not claim that cell or gene topology is globally superior, that integration
is generally beneficial or harmful, that sequencing technology causes the
observed differences, that the gene panel is inductive, that topology alone
identifies biological mechanism, or that current performance generalizes to
unseen datasets.

## The decision

The author team should choose one of two scientific ambitions:

1. **Methods-focused paper now (recommended).** Draft the revised paper around
   the corrected method, dual-view framework, transparent negative fusion
   result, robustness, and limits. During drafting, identify only analyses
   necessary to support a concrete sentence or figure.
2. **Generalization-focused paper.** First conduct a dataset/literature audit
   and prospectively select an external validation corpus. This is more
   valuable than running many additional gene-panel/metric variants on the
   same confounded data. A small gene robustness family may still be included,
   but only after its manuscript role is explicit.

Broad fusion searches, tissue-specific weights, and large clustering method
enumerations should remain closed. They would be difficult to distinguish from
post-outcome search without independent validation.

## Specific approvals needed

The project owner, Dr. Rouchka, and Dr. Mistry should review:

- the methods-focused versus generalization-focused target;
- whether gene topology is central enough to justify a small panel/metric
  sensitivity family;
- the intended journal and claim level;
- whether external data are required before drafting or can follow a first
  manuscript outline;
- authorship order and CRediT roles. Both Dr. Rouchka and Dr. Mistry remain in
  the credit plan, but exact roles and order must be agreed by the people
  involved rather than inferred from repository history.

## Proposed next sprint after approval

If option 1 is approved, the next sprint is a manuscript/figure claim map:
every proposed sentence and figure linked to a code artifact, data artifact,
result, limitation, and literature source. It should also regenerate the saved
persistence-diagram/landscape explanatory figure from corrected pipeline
artifacts.

If option 2 is approved, the next sprint is a read-only external dataset and
literature suitability audit. No dataset is downloaded and no analysis is run
until sample/study/tissue/technology coverage, gene compatibility, licensing,
compute, and a frozen validation estimand are approved.

## Gate disposition

G-MV7 is `owner_author_team_decision_required`. Existing data are sufficient
for a methods-focused paper, but insufficient for external-generalization or
causal technology claims. No further PH, new data, or manuscript claim
promotion is authorized by this dossier alone.

Repository verification at this gate passes 1,605 tests with zero failures or
warnings and four established skips (two CRAN-only, the optional Rust library,
and public audit documents intentionally excluded from R builds).
