# MV7-F omitted-34 upstream production audit

Date: 2026-08-15

Status: complete; exact 124-derived panel lock authorized

## Scope and boundary

MV7-F extended only the accepted private upstream layer to the 34 retained
single-study-tissue samples. It created or identity-validated 34 raw shards,
170 sample-seed SCT caches, and 170 deterministic 384-cell selections across
seeds `20260805` through `20260809`. It did not fit a gene panel, standardize
features, fit PCA, construct typed views, calculate PH or landscapes, open
labels, or calculate outcomes.

The implementation was frozen at `7097dcc`, and the exact 204-job queue was
durably frozen at `c71cdea`. The prefreeze independently passed 7/7 categories
and bound 19 public-safe source identities, 34 uniquely resolved private
sources, one worker, zero retries, 1,800-second and 8-GiB child caps, a
four-hour aggregate worker cap, and a 4-GiB private cache cap.

## Production result

The accepted run completed every frozen job:

- raw shards: 34/34;
- SCT caches: 170/170;
- deterministic selection identities: 170/170;
- raw dispositions: 24 built atomically and 10 validated reuses;
- SCT dispositions: 164 built atomically and 6 validated reuses;
- partial or unexpected cache files: zero; and
- panel, PCA, PH, landscape, label, and outcome jobs: zero.

The accepted execution used 917.069 raw worker-seconds and 6,034.980 SCT
worker-seconds, for 6,952.050 total worker-seconds. Peak monitored process-tree
RSS was 1,697,538,048 bytes and the complete private cache state was
1,618,456,042 bytes. All three measurements passed their prospective caps.

The first launcher wrapper had a 10-second client timeout. Its WSL process
completed five atomic raw entries before exiting without public evidence. The
accepted relaunch reopened those artifacts under the same frozen queue; this
accounts for four additional validated raw reuses beyond the six MV7-D
predecessors. No partial cache was accepted and no scientific contract changed.

## Independent and deterministic validation

The versioned independent validator reopened all caches and passed 9/9
categories:

- exact 204-job queue;
- 34/34 raw content identities;
- 170/170 SCT content identities;
- 170/170 reconstructed cell-selection hashes;
- 204/204 cache hashes and byte sizes;
- zero unexpected cache files;
- resource gates;
- upstream-only production summary; and
- execution receipt before source parsing.

The fixed representative sample-seed was rebuilt from a clean private root.
Raw counts, selected cells, and SCT payload identities matched, and both raw
and SCT serialized caches were byte-identical to the accepted artifacts.

The complete repository R test suite exited zero with four established skips.

The immutable check found 0 hash mismatches, 0 byte-size mismatches, and 0
manifest mismatches across 204 artifacts. CSV serialization rounded stored
mtime values by at most 5.0068 microseconds; the prospectively amended checker
uses a 100-microsecond tolerance and passed 204/204 artifacts. The private and
public execution receipts are byte-identical.

## Validator amendments

Two validation-only false negatives were retained rather than erased:

1. `vapply()` attached path names to recomputed hashes, so strict `identical()`
   rejected a vector whose 204 content hashes and sizes all matched. The
   correction was committed at `066e3dd`; amendment evidence at `4c79690`
   records zero cache mismatches and no production-runner change.
2. Exact in-memory comparison rejected a byte-identical CSV receipt after type
   round-tripping, and exact floating-point mtime comparison rejected
   microsecond-level CSV precision loss. The tested corrections were committed
   at `79a2a13`; amendment evidence at `9691120` records byte-identical receipts,
   the explicit 100-microsecond tolerance, zero cache mismatches, and no cache
   mutation.

Neither correction changes a cache, normalization, cell identity, resource
result, scientific estimand, or downstream interpretation. The accepted main
run itself exercised validated reuse for 10 raw and 6 SCT artifacts.

## Public and private separation

Public evidence contains only sample identifiers, dimensions, package/runtime
versions, hashes, byte sizes, dispositions, Boolean firewalls, and aggregate
resource measurements. It contains no expression matrices, private source
paths, home or WSL paths, PDFs, generated binaries, labels, or biological
outcomes. Raw/SCT caches, logs, receipts, and clean-repeat artifacts remain
under ignored private `tmp/` state.

## Decision

MV7-F passes. The next authorized stage is availability-only application of the
unchanged MV6-C global-core algorithm across the complete 620 SCT-cache corpus,
followed by an exact, independently validated 124-derived 500-gene panel lock.
PCA, typed-view construction, PH, landscapes, labels, and outcomes remain
closed until that panel is durable. The revised landscape definition remains
unchanged: finite positive intervals, essential H0 excluded, all consecutive
active levels, H0/H1 separate, exact or error-controlled squared-L2 distance,
no universal grid or level cap, and streamed/chunked computation.
