# MV11-L threshold-free common-K cross-view review

## Result

The matched historical cell and gene views are not redundant. Under the
prospectively primary PAM method, H1 agreement is numerically higher than H0 at
both common K values: mean ARI is 0.3206 at K=2 and 0.1974 at K=3 for H1,
versus 0.0570 and 0.0611 for H0. The H1 PAM K=2 seed range is 0.1559 to 0.4116;
the H1 K=3 range is narrower at 0.1700 to 0.2194.

Most hierarchical sensitivity means are near zero, and several seed ranges
cross zero. Single linkage at K=3 is especially seed-variable: H0 spans
-0.0354 to 0.2643 and H1 spans -0.0195 to 0.4897, while both medians remain
negative. No method/dimension/K combination has exact partition agreement in
any seed.

The fixed K=3-minus-K=2 contrast is -0.1232 for primary H1 PAM but +0.0040 for
primary H0 PAM. The other contrasts are retained in full without a selection
threshold. These are symmetric partition-agreement diagnostics: they do not
show that either cell or gene topology is superior and do not by themselves
authorize fusion or biological claims.

## Disposition

The evidence supports treating cell and gene topology as complementary views,
with a reproducible shared H1 partition signal under primary PAM and much less
shared structure under most hierarchical sensitivities. The next step is a
major scientific choice between a prospectively constrained fusion feasibility
test and constructing a newer all-QC cell topology to remove the historical
transform limitation before fusion.

Validation: 22/22 checks pass.
