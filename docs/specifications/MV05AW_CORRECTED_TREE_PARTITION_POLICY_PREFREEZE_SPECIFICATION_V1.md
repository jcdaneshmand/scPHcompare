# MV5-AW corrected-tree partition-policy prefreeze v1

Date: 2026-08-12

Authorization: MV5-AV `e26e04c`.

## Decision

No corrected-tree partition is currently identifiable from the MV5-AV smoke
artifacts. Each of the eight strata has exactly three samples. The admissible
grid `2:min(10,n-1)` therefore contains only `k = 2`, making selection vacuous.
Each stratum also has one matrix rather than the five matched seed matrices
required by the existing label-closed stability rule.

The trees remain valid descriptive, label-free hierarchical objects. Cutting
them at two clusters would be mechanically possible but methodologically
unjustified and is not authorized.

## Requirements before reconsideration

A later gate requires, within each independently analyzed H0/H1 and cell/gene
stratum: a larger complete sample matrix; at least two candidate values of `k`;
matched resampled or seed-specific matrices on identical sample axes; a frozen
label-closed stability statistic and uncertainty rule; and prospective
resource evidence. Labels, outcomes, the descriptive combined matrix, and
cross-view fusion remain prohibited during selection.

No default, legacy path, tree, source artifact, or scientific result changes in
this sprint.
