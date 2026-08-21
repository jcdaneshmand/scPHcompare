# MV8-L exact-500 view-specific observation-scope closure

## Decision

**Do not advance a view-specific exact-500 gene-topology contract.**

MV8-L tested the strongest non-gene-aware version of the proposed rescue: gene-side observations were expanded from the fixed 384-cell subset to all 4,614 cells passing the frozen QC rule, in deterministic barcode order. The existing cell-side 384-cell/immutable-30-PC comparator was not changed or rerun.

Both tested SCT configurations retained every exact-panel gene and produced finite standardized values. Nevertheless, each yielded one exactly zero-variance standardized gene. Correlation—and therefore the frozen correlation-chord gene distance—is undefined for that vector. No feature was dropped, imputed, padded, or identified publicly to evade the gate.

## Aggregate result

| Configuration | Exact-500 retained | Finite standardized values | Zero-variance genes | Correlation-chord gene distance |
| --- | ---: | --- | ---: | --- |
| Standard SCT (`min_cells=5`) | 500 / 500 | Yes | 1 | Not admissible |
| Low-count SCT (`min_cells=1`) | 500 / 500 | Yes | 1 | Not admissible |

The least-detected exact-panel gene was observed in 33 of the 4,614 QC-eligible cells in the raw matrix. Thus the failure is not missing panel coverage and was not repaired by a lower SCT inclusion threshold; it is a limitation of the tested transformed representation for the frozen correlation-based gene view.

Two fresh sentinel runs matched on all scientific aggregate fields. Resource use was 17:26 / 6,097,884 KiB and 19:08 / 6,166,200 KiB peak RSS, respectively; both stayed below the 30-minute and 12-GiB prospective caps. The intentional exit code `2` means an admission check declined, not that the calculation crashed.

## Method rationale evaluated—not adopted

The prefreeze was motivated by the fact that SCTransform performs a regularized negative-binomial transformation and, by default, drops genes expressed in fewer than five cells. [Seurat reference](https://satijalab.org/seurat/reference/sctransform) documents this behavior and the SCT output layers. Larger numbers of cells can improve single-cell gene-network reconstruction, but association metric choice also materially affects such networks, as demonstrated by [Skinnider, Squair, and Foster (2019)](https://doi.org/10.1038/s41592-019-0612-7). Those considerations justified the feasibility test; they do not justify changing the established correlation-chord metric after the result failed.

## Consequences and firewalls

Common-475 remains the admitted external dual-view HCA analysis. Exact-500 raw-read recovery remains useful panel-coverage evidence, but it is not an external exact-500 gene-topology route under the current representation.

No other HCA unit was opened at matrix level, and no persistence diagrams, landscapes, labels, outcomes, fusion, clustering, manuscript claims, or deletions were performed. Any future exact-500 strategy would require a new owner-approved scientific contract that changes either the representation or the gene association geometry; it is not authorized by MV8-L.
