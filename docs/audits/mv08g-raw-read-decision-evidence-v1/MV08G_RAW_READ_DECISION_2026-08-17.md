# MV8-G raw-read decision audit

Date: 2026-08-17

## Decision

The frozen 500-versus-475 reference sensitivity independently reproduces as `material_panel_sensitivity`. Exact-500 HCA raw-read reprocessing is therefore recommended before making an external topology claim. This audit authorizes prospective planning only; it does not authorize FASTQ download, raw processing, label access, or biological claims.

## Why

The HCA processed matrices exactly match the historical Cell Ranger 3.0.0 Ensembl93 filtered reference. That reference excludes 25 of the ordered 500 target genes through its documented biotype filter, so identifier crosswalks or zero imputation cannot recover their counts.

The corrected landscape calculation used every active finite positive-persistence level, excluded essential H0, kept H0 and H1 separate, used exact or error-controlled squared L2 integration, and applied no grid or level cap. All 40 landscape groups, eight repeats, and 24 R/Persim oracles passed.

The paired comparison passed all 13 independent checks. Cell H1 is the decisive component: median Spearman is 0.916 (threshold 0.95) and median top-10 overlap is 0.702 (threshold 0.80), while fixed-k PAM ARI is 0.814. Cell H0 and both gene components pass all three prospective classification thresholds.

## Scope of the recommendation

The common-475 results remain useful as a transparent secondary harmonization sensitivity. They are not an adequate sole basis for claiming that external cell-H1 topology replicates the exact-500 reference analysis.

## Required next gate

Before download, MV8-H must prospectively freeze the 48-FASTQ/eight-unit manifest; the Ensembl93 custom reference retaining all 500 ordered stable IDs; software and commands; unit-level QC and 384-cell rule; resource/storage caps; the label firewall; and independent count, topology, PH, all-level H0/H1 landscape, distance, and clustering validation. The current FASTQ inventory is 85,034,239,918 bytes (79.194 GiB), so owner resource approval is required.
