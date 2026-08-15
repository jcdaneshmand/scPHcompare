# MV7-F omitted-34 upstream production prefreeze v1

Date: 2026-08-15

MV7-F extends only the accepted individual-source/raw/SCT path to the 34
retained single-study-tissue samples. It creates 34 atomic raw shards, 170 SCT
caches across seeds `20260805` through `20260809`, and 170 deterministic
384-cell identities. It does not fit the 124-sample panel, standardization,
PCA, typed views, PH, landscapes, distances, clustering, labels, or outcomes.

The queue is exact and serial: all 34 lexicographically ordered raw jobs first,
then 170 sample-seed SCT jobs ordered by seed and sample. There is one worker,
no retry, and stop-on-first-failure. A previously created cache may be reused
only after its complete identity and payload validate; stale state is an abort.

Each raw and SCT child has an 1,800-second and 8-GiB process-tree RSS cap. The
aggregate worker-time cap is four hours and the final unique private-state cap
is 4 GiB. These exceed the MV7-D linear projection by the frozen safety margin.
Crossing a cap stops publication and triggers a new prospective resource gate.

All cache writes are temporary-file plus atomic rename. Partial files, missing
audits, duplicate identities, unexpected source counts, nonfinite SCT values,
or cell-selection drift stop the stage. The final public evidence contains
only sample IDs, hashes, dimensions, versions, dispositions, and resource
summaries; private paths and expression values are excluded.

An independent validator must reopen and validate all 34 raw and 170 SCT
caches, reconstruct all selected-cell hashes, verify the 204-job queue, enforce
the resource and label gates, and reject partial state. A clean representative
repeat uses the lexicographically first added sample at seed `20260809` and
must match raw counts, selected cells, SCT payload, and cache bytes. A complete
resume pass must preserve every accepted cache hash, size, and modification
time.

Passing MV7-F authorizes only the availability-only global-core inventory over
the complete 620-cache corpus and the exact 124-derived 500-gene panel lock.
PCA and PH remain closed until that panel is durable and independently
validated.
