# MV6-D matched SCT source/PH/landscape profile

| Field | Result |
|---|---|
| Date | 2026-08-14 |
| Prefreeze commit | `88937d9` |
| Implementation commit | `9a8e42a` |
| Fold sentinels | 5, one per seed and held-out-size rank 1/4/8/12/15 |
| Profiled samples | 10, one held-out and one training sample per sentinel |
| Source / PH / landscape units | 5 / 20 / 10 |
| Stage-1 exact repeat | 7/7 files byte-identical |
| Independent validation | 12/12 categories pass |
| Labels/outcomes | Closed / zero |
| Decision | `revise_for_targeted_acceleration` |

## Outcome

Option A remains viable. The fixed MV6-C 500-gene global core supports the real
fold-local corrected pipeline in both orientations. Every sentinel rebuilt its
training-only scaling and 30-PC cell transform from the accepted caches, and
every selected sample produced both a 384-cell view and a 500-gene view without
zero filling, panel shrinking, or sample exclusion.

All 20 complete Vietoris–Rips H0/H1 jobs passed typed provenance, interval,
essential-H0, and independent metric-MST checks. All ten held-out/training
landscape comparisons reproduced the dissertation-aligned H0/H1 definition.
The bounded profile therefore rejects neither gene topology nor future fusion.

The prospective 72-worker-hour full-run gate does not pass. The reason is not
PH. It is evaluating every landscape pair as a separate public-R-API process,
especially for gene H1 diagrams with deep interval sets. The correct next
action is a production-style batched/streamed landscape acceleration gate, not
abandoning fusion and not rewriting Ripser.

## Frozen sentinel coverage

The deterministic selector sorted studies by held-out sample count and study
ID, selected ranks 1, 4, 8, 12, and 15, and paired them with the five seeds.
This covered held-out study sizes 1, 2, 4, 8, and 25. Within each fold, the
largest accepted cache in each role was selected, with sample ID as the tie
breaker. No tissue, assay approach, endpoint, outcome, PH, landscape, or fusion
result entered selection.

Stage 1 used the rank-8/seed-20260807 fold. It completed one source bundle, four
PH jobs, and two landscape pairs, then rebuilt and reran all seven artifacts.
All seven repeated files matched byte count and SHA-256 exactly. Both gene-PH
jobs remained far below the stricter 900-second and 6-GiB continuation gate, so
the other four sentinels executed as prefrozen.

## Correctness and topology profile

| View | PH jobs | H1 intervals per diagram | Median PH | Maximum PH | Maximum process-tree RSS |
|---|---:|---:|---:|---:|---:|
| Cell topology | 10 | 196–486 | 1.042 s | 1.058 s | 86,319,104 B |
| Gene topology | 10 | 105–2,540 | 1.046 s | 2.068 s | 109,096,960 B |

For cell views, finite H0 deaths were compared with a Prim MST over Euclidean
30-PC coordinates. For gene views, the independent oracle operated on the
actual 500-gene Pearson-correlation chord-distance matrix. All 20 multisets
match within the prospective tolerance. Every diagram contains one essential
H0 class; all H1 intervals are finite and have positive persistence.

The independent validator separately reconstructed the ten-row sentinel
manifest, rehashed all private outputs, validated all five typed source
bundles, recomputed all 20 metric MSTs, reran all ten landscape calculations,
reconstructed both resource projections and the prospective decision, and
checked the public schema. All 12 categories pass.

The complete source-loaded R suite passes after the final accounting
correction, including all 25 MV6-D expectations. Its two skips are the
established optional Rust-library and public-audit-in-build exclusions, not
MV6-D failures. The separately documented installed-package PH child-process
boundary remains outside this sprint.

## Source and landscape resources

Fold source/PCA construction includes validating and deserializing all 90
accepted caches for the seed, extracting the fixed panel, fitting scaling on
training cells only, fitting PCA on training samples only, and constructing
the two matched sentinel views.

| Unit | Measured units | Median | Maximum | Maximum process-tree RSS |
|---|---:|---:|---:|---:|
| Fold source/PCA | 5 | 123.983 s | 125.518 s | 2,145,976,320 B |
| Cell landscape pair process | 5 | 5.664 s | 8.249 s | 113,840,128 B |
| Gene landscape pair process | 5 | 28.113 s | 76.211 s | 122,781,696 B |

All five cell pairs used the exact R route for H0 and H1. All gene H0
components were exact; two of five gene H1 components were exact and three
used the error-controlled adaptive route because at least one diagram exceeded
the 500-interval exact-engine routing guard. No interval or landscape level was
discarded, and no grid or level cap was introduced.

The bounded run, including the full stage-1 repeat, used 1,002.898 monitored
worker-seconds and 16,136,889 private bytes. No unit approached its hard
elapsed or RSS cap.

## Full-workload projection

| Component | Median scenario | P90 scenario | Maximum scenario |
|---|---:|---:|---:|
| 75 fold source/PCA builds | 2.583 h | 2.615 h | 2.615 h |
| 6,750 cell PH jobs | 1.954 h | 1.984 h | 1.984 h |
| 6,750 gene PH jobs | 1.961 h | 3.360 h | 3.877 h |
| 35,350 cell landscape-pair processes | 55.618 h | 81.000 h | 81.000 h |
| 35,350 gene landscape-pair processes | 276.053 h | 748.352 h | 748.352 h |
| **Total** | **338.169 h** | **837.311 h** | **837.828 h** |

This projection is deliberately faithful to the bounded runner: every pair
starts an R process, loads two records, and invokes the public reference API.
It is not a claim that a production-batched implementation must take 838 hours.
The accepted MV5-D4 cell calculation already demonstrated that grouped exact
critical-pair execution can compute all 35,350 cell pairs in 1.165 worker-hours.
Therefore the large pairwise projection identifies a missing production route
for the new matched cell/gene diagrams; it does not establish an inherent
mathematical cost of the landscape definition.

The reconciled private-storage projection is 8,738,555,180 bytes, under the
10-GiB cap. Gene-view distance objects account for 6,961,389,750 projected
bytes. The initial sidecar mixed serialized bundle bytes with in-memory object
sizes; it was rejected and regenerated on a single serialized-byte basis before
acceptance. A full implementation should construct each gene
view, emit its diagram, and discard the distance object rather than retaining
6,750 gene views. That streaming change does not alter any scientific value.

## Landscape definition and acceleration boundary

The revised dissertation-aligned definition remains unchanged:

- H0 and H1 are separate;
- only finite positive-persistence intervals enter landscapes;
- every active consecutive level is retained;
- missing levels are zero-padded;
- L2 integration is exact or error-controlled on dimension-specific support;
- there is no universal level cap or fixed uniform grid; and
- any combined distance remains secondary to recorded H0 and H1 components.

Acceleration must reproduce that contract. The next gate should compare a
production-style batched exact critical-pair route and the already certified
Rust landscape candidate against the canonical R results on the frozen MV6-D
cell and gene diagrams. It should measure grouped throughput, memory, resume,
and streaming storage. Rust is relevant to the landscape kernel only if it
adds a material advantage over the existing grouped implementation; it is not
currently justified for PH or source/PCA reconstruction.

## Gate disposition

MV6-D completes with `revise_for_targeted_acceleration`. It authorizes only a
separately prefrozen, correctness-first MV6-E batched landscape admission and
resource sprint. It does not authorize the 75-fold matched rerun, blocked
fusion, clustering, outcomes, G-MV6, new data, integrated gene topology,
default changes, manuscript claims, release actions, binaries, PDFs, reviewer
material, or `example_run.r`.

The result is constructive: Option A solved the matched-panel problem; both PH
orientations are inexpensive and correct; and the remaining obstacle is a
well-localized landscape execution strategy with existing candidate engines.
