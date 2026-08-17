# MV8-G common-475 paired reference sensitivity prefreeze v2

Status: prospective replacement for the stopped v1 source execution.

This contract incorporates every scientific, landscape, comparison, resource,
label-firewall, and HCA stop rule in
`MV08G_COMMON475_PAIRED_REFERENCE_SENSITIVITY_PREFREEZE_V1.md`. It changes only
the typed source-shape implementation needed to represent the prospectively
frozen 475-gene panel.

## Reason for replacement

The first v1 source job stopped under its zero-retry policy before producing a
source bundle. The input was correctly shaped as 475 genes by 384 cells, but
the reusable dual-view constructor recognized only the original 500-gene
scientific profile and rejected the input as `475 x 384` rather than
`500 x 384`. The job used 1,950,502,912 peak process-tree RSS bytes, ran for
151.452103376389 seconds, and did not approach either resource cap. This was an
implementation-contract defect, not a data, landscape, or resource failure.

## Corrected typed profile

The replacement implementation adds an explicit
`scientific_common475` profile with exactly:

- 475 ordered genes;
- 384 matched cells per sample;
- 30 shared cell-view principal components; and
- scientific eligibility for both typed topology views.

The original `scientific` profile remains exactly 500 genes by 384 cells with
30 principal components. The `analytical_fixture` profile remains
non-scientific. No shape is inferred from dimensions, and neither scientific
profile permits caller-selected dimensions.

The source record validator must confirm that the PCA model and every cell and
gene view carry `scientific_common475` and remain scientifically eligible.

## Execution boundary

The stopped v1 ledger is immutable and is not resumed, adopted, deleted, or
overwritten. A new commit, new prefreeze evidence, and new private execution
root are required. The v2 source gate again authorizes only five source bundles
and one exact repeat with zero retries. PH, landscapes, comparisons, labels,
HCA FASTQ download, and raw-read reprocessing remain unauthorized until their
respective independent gates close.

The persistence-landscape definition is unchanged: finite positive-persistence
intervals, essential H0 excluded, every consecutive active level, H0 and H1 separate,
exact or error-controlled squared-L2 integration, no fixed grid, no level cap,
and streamed or chunked computation.
