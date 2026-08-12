# MV5-AT broader realistic corrected-landscape workflow smoke

Date: 2026-08-12

Authorization: MV5-AS `4d427ce`; prospective specification `9444c1d`; bound
scope `818d505`; independent validator `7684bb3`.

## Outcome

MV5-AT passes. The actual corrected-only postprocessing workflow processed all
eight frozen MV5-AP-R1 strata as independent, sequential, one-worker units:
24 existing persistence diagrams and all 24 within-stratum pairs. No new data,
labels, outcomes, cross-stratum pairs, or downstream analysis were used.

The revised landscape definition remained intact: all finite intervals and all
active levels, H0 and H1 separate, exact or strict error-controlled integration,
and descriptive combined distance only. H0 was exact in every unit. H1 was
exact in six units and used `adaptive_quadpack_partitioned_v2` in the two
high-depth gene units. All 48 dimension results certified at `1e-8`.

The prospective scope ledger described the adaptive route with the planning
label `adaptive_global_error_v1`; the engine's actual provenance identifier is
`adaptive_quadpack_partitioned_v2`. The validator was corrected to the actual
identifier after its first safe stop. This is a nomenclature correction, not a
scientific or numerical change.

## Artifact, resume, and resource validation

All eight completion markers verify. The independent validator loaded every
pair shard and exactly reconstructed every H0, H1, and combined matrix entry.
Completed resume preserved every artifact path, size, modification time, and
SHA-256. All 24 source diagram files remained unchanged, and no unit populated
either legacy landscape field.

Twelve independent validation categories pass; 15 prohibited-change counters
are zero. Process wall time ranged from 28.11 to 578.21 seconds and peak RSS
from 858,382,336 to 936,157,184 bytes. Every unit stayed below 750 seconds and
2 GiB. Scientific times for the adaptive units were 432.755 and 564.856 seconds.

Private diagrams, shards, RDS matrices, and process logs remain ignored. Only
public-safe aggregate CSVs are tracked.

## Decision

MV5-AT authorizes MV5-AU only: a prospective corrected-matrix consumer
prefreeze. Downstream consumption remains off. Defaults, legacy rewrites, new
data, biological claims, clustering execution, and Rust/optimization remain
unauthorized.
