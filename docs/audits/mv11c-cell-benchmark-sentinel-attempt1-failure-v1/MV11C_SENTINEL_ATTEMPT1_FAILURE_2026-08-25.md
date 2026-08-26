# MV11-C sentinel attempt 1 failure audit

The first prospectively authorized sentinel stopped before clustering output.
The worker compared the in-memory catalog with its CSV-frozen queue using
R's type-strict `identical()`: `unordered_pairs` was numeric in memory and
integer after CSV parsing, although both values were exactly 7,626. The worker
therefore failed closed with `MV11 worker source catalog drift`.

No partition, quality, status, label, outcome, cross-view, fusion, biological,
or manuscript artifact was produced. Peak process-tree RSS was 74,981,376
bytes, elapsed time was 1.467 seconds, and all partial evidence remains in the
private ignored execution roots. The correction compares frozen catalog fields
by their canonical character values; it does not change any scientific value,
method, K, source, cardinality, or resource contract.

A second execution is not an automatic retry. It requires a newly committed
implementation and a new prospective prefreeze binding this failure record.
