# MV7-B no-new-PH confounding diagnostics prefreeze v1

Date: 2026-08-15

MV7-B is a post-outcome, prespecified diagnostic analysis. It reuses the locked
MV6-H sample-level summaries and authoritative metadata; it computes no PH,
landscapes, distances, rankings, clustering, or fusion. The six fixed methods
are cell H0/H1, gene H0/H1, and the two already-defined descriptive composites.
Fusion weights are excluded. MRR is primary diagnostic and 1-NN balanced
accuracy supportive.

For each method/endpoint, recompute the equal-five-tissue macro estimate and
delete each of 15 studies in turn while retaining all five tissues. Also delete
each tissue and average the remaining four tissues. Report every deletion and
the maximum absolute change; do not suppress influential groups. Fixed
cell-minus-gene contrasts are H0, H1, and composite, and their sign stability is
reported under every deletion.

Retained-cell association uses within-study centered ranks of log2 retained
cell count and the outcome, pooled across informative studies. Its interval is
a 2,000-replicate study-block bootstrap with seed 20260817. This is descriptive
association, not causal adjustment. Library size is unavailable and retained
cells are not called a proxy for it.

Approach association is restricted to the three studies containing both
scRNA-seq and snRNA-seq. Within each mixed study, calculate snRNA minus scRNA;
average the three study effects equally and use a 2,000-replicate study
bootstrap with seed 20260818. Also publish the design inventory. The result
cannot identify a causal technology effect and cannot be generalized to colon
or liver, which lack snRNA-seq in the accepted 90-sample axis.

Prospective flags are: absolute delete-one change at least 0.05; any sign change
in a fixed cell-minus-gene contrast; absolute within-study count correlation at
least 0.30 with a percentile interval excluding zero; or absolute mixed-study
approach difference at least 0.10 with an interval excluding zero. Flags narrow
claims; they do not choose a method, view, gene panel, metric, or fusion weight.
No p-values or multiplicity-driven discoveries are produced.

All input hashes, method/endpoint/contrast registries, seeds, thresholds,
metadata fields, output schemas, and prohibitions must be committed before the
diagnostic runner executes. A label-access receipt is written before parsing
metadata. A standalone validator must reconstruct every result and a second
exact-commit run must reproduce every production artifact byte-for-byte.
