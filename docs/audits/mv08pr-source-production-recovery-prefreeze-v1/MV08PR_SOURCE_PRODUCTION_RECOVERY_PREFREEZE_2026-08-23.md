# MV8-PR source-production recovery prefreeze

## Preserved stop

MV8-P stopped safely after 123 completed fits when job 124 attempted to export 654.06 MiB through R `future` while the runtime default allowed 500 MiB. The failed child exited after 23.13 seconds at 2.12 GB peak RSS, so neither the 1,800-second nor 12 GiB source-process cap was approached. The original private and public v1 roots, including the failed ledger row and child logs, remain immutable.

## Approved recovery

The project owner approved a fresh overlay for original queue jobs 124 through 129. The overlay permits exactly one explicit second attempt for job 124 and first attempts for jobs 125 through 129. It raises only `future.globals.maxSize` to a bounded 2 GiB; one worker, the 1,800-second cap, the 12 GiB process cap, source hashes, selected-cell axes, SCTransform parameters, Pearson-residual definition, and geometry gates remain unchanged.

## Closure

MV8-Q must independently rehash the 123 accepted v1 caches plus six completed overlay caches, preserve and hash-bind the stopped v1 evidence, and then bind those 129 sources to the three MV8-O primary sources for 132/132 coverage.

## Scope firewall

This amendment authorizes source fitting only. Persistence diagrams, landscapes, comparisons, clustering, fusion, labels, outcomes, manuscript claims, default adoption, cleanup, and deletion remain closed.
