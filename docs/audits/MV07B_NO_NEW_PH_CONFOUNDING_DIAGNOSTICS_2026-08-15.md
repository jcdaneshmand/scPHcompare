# MV7-B no-new-PH confounding diagnostics audit

## Disposition

`narrow_claims_before_gene_rerun_or_new_data`

The exact-commit MV7-B diagnostic completed without recomputing PH,
landscapes, distances, rankings, clustering, or fusion. Independent
reconstruction passes 7/7 categories and all 12 production artifacts repeat
byte-for-byte.

## Design

The accepted axis contains 90 samples, 15 studies, and five tissues. Retained
cell counts range from 421 to 9,071 (median 2,676). Eighty samples are scRNA-seq
and ten are snRNA-seq; only three studies contain both approaches. Library size
and harmonized cell-type composition are unavailable.

The six fixed methods are cell H0/H1, gene H0/H1, and the two descriptive
H0/H1 composites. Fusion is excluded. MRR is the primary diagnostic and 1-NN
balanced accuracy supportive. No p-values or outcome-driven selection were
computed.

## Influence result

Study and tissue influence exceed the prospectively frozen 0.05 threshold.
For composite MRR, deleting study `SRA716608` changes the gene composite from
0.36943 to 0.43896 (change +0.06953). Deleting `SRA779509` changes it by
+0.05174. The cell composite's largest study deletion is +0.03581.

Tissue deletion is larger, as expected from an equal-five-tissue estimand:
removing PBMC changes cell/gene composite MRR by -0.11227/-0.12165, while
removing liver changes the cell composite by +0.09444. These rows are evidence
that corpus-wide performance is materially tissue-composition dependent, not a
reason to remove a tissue.

The full cell-minus-gene composite MRR contrast is +0.03706. It changes sign
when `SRA716608` is deleted (-0.03003) and when colon is deleted (-0.03266).
Cell-H0 minus gene-H0 also changes sign after deleting `SRA716608`. Therefore
the relative ordering of the two views is not stable enough to support a
general claim that one view is superior.

## Retained-cell and approach diagnostics

The retained-cell flag does not trigger. The within-study rank correlation is
0.254 for the cell composite and -0.151 for the gene composite, but their
study-block intervals are broad and cross zero: [-0.252, 0.454] and
[-0.531, 0.464]. The result does not establish absence of cell-count
confounding; it says the current diagnostic does not detect a stable monotone
within-study association.

The approach flag also does not trigger. Across the three mixed studies, the
snRNA-minus-scRNA composite MRR differences are +0.016 for cell and +0.119 for
gene, with intervals [-0.126, 0.174] and [-0.192, 0.522]. With only three mixed
studies and no snRNA samples in colon or liver, these are descriptive and do
not identify a causal technology effect.

## Scientific consequence

The final manuscript should report cell and gene topology as complementary,
separate views whose relative retrieval performance varies by study and tissue.
It should not claim that cell topology globally outperforms gene topology, that
gene topology is uniformly stronger, or that sequencing approach causes the
observed differences. The negative equal-weight fusion result remains intact.

The diagnostics do not automatically justify a costly gene-panel or gene-
metric rerun, and they do not automatically trigger new data. MV7-C should now
synthesize the existing-data evidence, define the narrowest defensible claims,
and present the author team with an explicit choice between a methods-focused
existing-data paper and additional validation/robustness computation.

## Reproducibility

- exact execution commit: `758780adec641f23d62e5b322224bfdf9e679ba9`;
- prefreeze validation: 8/8;
- outcome reconstruction: 7/7;
- production repeat: 12/12 byte-identical;
- first execution attempt was interrupted by the orchestration timeout after
  deterministic pre-bootstrap tables; its partial output is preserved under
  ignored `tmp/` and the accepted run started in a clean directory.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | approve as post-outcome diagnostic |
| Correctness demonstrated? | pass; independent 7/7 and repeat 12/12 |
| Computation feasible? | yes; no new PH and bounded block bootstraps |
| Biological interpretation permitted? | descriptive/secondary with study, tissue, and approach limits |
| Next action | MV7-C existing-data synthesis and author decision dossier |
