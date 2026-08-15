# MV7-E metadata provenance and full-124 descriptive-estimand prefreeze v1

Date: 2026-08-15

## Purpose and boundary

MV7-E resolves the 16 sequencing-approach disagreements and freezes the
label-closed mathematical contract for a separate 124-sample descriptive
analysis. It does not alter or enlarge the accepted 90-sample cross-study LOSO
estimand. It computes no PCA, PH, landscapes, distances, clustering, retrieval,
or biological outcomes.

## Approach provenance

The historical retained table is a join between `filtered_cells.csv` and the
public Panglao-derived candidate metadata. `Approach.x` came from an expression
heuristic in the historical pipeline, while `Approach.y` came from public
metadata. `Approach.y` matches the public 127-row candidate metadata for all
124 retained samples. `Approach.x` differs for exactly 16.

The fixed accession table links those 16 SRS records to official NCBI GEO
sample records and method descriptions. All 16 official records agree with
the public value: 14 are single-cell and two substantia-nigra samples are
single-nucleus. Future scientific approach labels therefore use the public
metadata, confirmed by GEO where disputed. Expression characteristics never
classify assay technology.

This correction makes the accepted 90-sample axis 90/90 scRNA-seq with zero
mixed-approach studies. MV7-B's mixed-study approach association is therefore
not estimable and is superseded. Its topology, landscape, study, tissue,
retained-cell, and influence results are unchanged because they did not depend
on the approach field. No causal technology claim is authorized. In the full
124, snRNA-seq remains nested in the substantia-nigra study/tissue, so the
descriptive join cannot identify a technology effect.

## Sample, cell, panel, and transform contract

The descriptive population is the exact 124 retained samples. The five seeds
are `20260805` through `20260809`, with 384 deterministically selected cells per
sample and seed. There are 620 sample-seed states.

Availability of the accepted MV6-C 500-gene panel is tested using only feature
names from the 34 added sparse sources after reproducing Seurat's historical
underscore-to-dash normalization. No expression value, topology, or biological
label enters the decision. One sample, `SRA701877_SRS3279688`, lacks
`KLF2-ENSG00000127528.5`; the other 33 contain all 500. The prospectively
specified fallback is therefore mandatory.

After MV7-F creates the 170 missing SCT caches, apply the unchanged
`mv06c_global_core_panel_v1` algorithm to all 620 caches: exact feature
intersection, finite and positive within-cache variance, the fixed technical
gene exclusions, deterministic canonical de-duplication, median variance rank,
then deterministic feature tie-breaking. Select exactly 500 genes and report
overlap with the accepted 90-derived panel. No label or topology result may be
read. The exact resulting panel must be committed and independently validated
before any fitted coordinates or PH.

For each seed independently, pool the equal-depth cells from all 124 samples,
fit one global gene center/scale and one conventional 30-component cell PCA,
then apply them to those same 124 samples. This is explicitly transductive.
It must not be described as the same estimand as training-only LOSO PCA.

## Typed views and topology

`cell_topology_v1` contains 384 selected cells as points in the 30 global
descriptive PCs with Euclidean geometry. `gene_topology_v1` contains the same
500 panel genes as points described by the same 384 cells; its distance is
correlation chord, `sqrt(2 * (1 - r))`, with numerical clamping to `[0, 2]`.
Nonfinite or zero-variance inputs stop; there is no imputation.

Both views use complete Vietoris-Rips PH in H0 and H1. Essential H0 is excluded.
The landscape definition remains: all finite positive-persistence intervals,
all consecutive active levels, H0/H1 separate, exact or error-controlled
squared-L2 integration on dimension-specific support, no universal grid, no
universal level cap, and streamed/chunked pair calculation.

Pairs are unordered lexicographic sample pairs within seed only: 7,626 per
seed, 38,130 per view over five seeds, 76,260 view-specific pairs, and 152,520
H0/H1 component-distance rows. Cross-seed distances are prohibited.

## Resource, atomicity, and label firewall

MV7-F is limited to 34 raw shards, 170 SCT caches, and their deterministic cell
identities. Each child has an 1,800-second and 8-GiB cap; aggregate production
has a four-hour worker-time and 4-GiB storage gate. Writes are atomic,
validated existing entries are immutable, resume is identity-checked, and no
partial state is publishable.

The post-MV7-F measured projection must set the MV7-G sentinel cap. MV7-G must
set the MV7-H cap. No full PH is authorized by this prefreeze.

Panel fitting, standardization, PCA, typed-view construction, PH, and landscape
production remain label-closed. Sample IDs may route artifacts; tissue, study,
approach, and biological outcomes may not select features, fits, methods,
sentinels beyond the already frozen design strata, stopping rules, or retained
results.

## Decision

Passing MV7-E authorizes MV7-F upstream production only. Because the fallback
was triggered by availability, the exact 124-derived panel must be frozen in a
separate post-MV7-F evidence commit before MV7-G. PCA, PH, landscapes,
clustering, outcomes, external data, default changes, and primary-90
recalculation remain closed.
