# MV8-K exact-500 transform-contract closure

## Decision

**Do not admit exact-500 HCA topology under the current frozen contract.**

The result is a closed, label-free negative gate rather than a data or runtime
failure.  The recovered HCA_BM_002 raw-read matrix contains the whole
exact-500 panel, but the frozen 384-cell transform cannot simultaneously
preserve all 500 features *and* supply valid gene-side correlation-chord
geometry.

## Fixed comparison

Both configurations used the same HCA_BM_002 filtered matrix, ordered
exact-500 panel, frozen 384-cell selection (seed `20260805`), global
center/scale vectors, and immutable 30-PC reference rotation.  No reference
model was re-fit.

| Configuration | Panel retention | Finite standardized values | Frozen 30-PC projection | Effective zero-variance features | Separate views admitted |
| --- | ---: | --- | --- | ---: | --- |
| Standard SCTransform (`min_cells=5`) | 497 / 500 | Not evaluated after failed retention | Not evaluated | Not evaluated | No |
| Low-count-preserving SCTransform (`min_cells=1`) | 500 / 500 | Yes | Yes, 384 x 30 | 3 | No |

The standard configuration drops three panel genes.  Lowering only the SCT
inclusion threshold retains the full panel and yields finite standardized
values and a finite immutable PCA projection, but three standardized genes
are exactly constant across the frozen selected cells.  A constant feature has
undefined correlation; therefore the required 500-gene correlation-chord view
cannot be constructed.  The existing `new_dual_view_source()` nonzero-variance
guard was retained and made explicit in the audit runner; no genes were
dropped, padded, or silently altered.

## Deterministic and resource evidence

Two fresh low-count-preserving runs produced the same selected-cell identity,
panel retention, zero-variance count, and standardized-source SHA-256
`6a74c707fa61ebbe5dbc26b14d7d7bed565cdfb60776b1005501fbd915385ee0`.
The repeat resource envelopes were 143.32 seconds / 1,148,080 KiB and 142.71
seconds / 1,158,876 KiB peak resident memory.  Their nonzero exit status is
intentional: the runner exits `2` whenever an admission check fails, after
publishing aggregate diagnostics.

## Firewalls and next boundary

No persistence diagrams, landscapes, fusion, clustering, labels, outcomes,
or manuscript claims were computed.  The seven exact-500 Cell Ranger recovery
units remain useful raw-read evidence, but they are not authorized for
exact-500 topology by this result.

Resolving this gate would require a **new scientific contract**, not an
implementation workaround—for example, a prospectively justified
coverage-aware cell-selection rule or a separately evaluated transform/metric
definition.  Neither alternative is adopted here.
