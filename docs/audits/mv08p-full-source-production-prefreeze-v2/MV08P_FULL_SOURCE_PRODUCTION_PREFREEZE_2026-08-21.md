# MV8-P full Pearson-residual source-production prefreeze

## Decision

After this prefreeze is committed, exactly 129 remaining source fits may run serially in ascending frozen-QC cell count. The three completed MV8-O primary source fits are not rerun. Each new fit uses the MV8-N all-QC standard-SCTransform/Pearson-residual contract and the unchanged selected-384 axes.

## Resource policy

The policy remains one worker, zero automatic retries, 1,800 seconds, and 12 GiB per source process. The maximum completed sentinel used 12.708 GB, leaving 176.7 MB (about 1.37%) below the cap; the largest remaining source has 9,071 rather than 11,475 QC cells. Runs remain serial and stop safely on a resource/precondition failure.

## Scope firewall

This authorizes only source fitting and private cache construction. Persistence diagrams, landscapes, pairwise comparisons, clustering, fusion, labels, biological outcomes, manuscript claims, and default adoption remain closed.

## Next gate

A source-production closure must rehash every private cache and publish aggregate resource/geometry evidence before any PH or landscape execution is considered.
