# MV8-H Cell Ranger 8.0.1 count-sentinel execution closure

## Outcome

The one prospectively selected complete unit `HCA_BM_002` completed successfully
under the frozen Cell Ranger 8.0.1 command. This document is an execution
closure, not a replacement for the historical pre-execution prefreeze. The
prefreeze decision remains preserved by SHA-256 and records that count had not
yet executed at that earlier timestamp.

The completed command identity is bound by SHA-256 `0b2e322bf2b3a86eb9eeb65ec6a783d31e9ae9cb7b5a840b0cd126085cf0dfdf`. Its
scientific controls remain exact: `MantonBM2_HiSeq_9`, SC3Pv2, expected cells 7000,
exon-only counting, BAM disabled, secondary analysis disabled, 4 cores, and
32 GiB. The selected six FASTQs total `11,249,623,632` bytes and
remain bound to the historical prefreeze manifest. The custom reference,
runtime tree, and launcher identities remain unchanged and hash-bound.

## Terminal evidence

The pipestance success marker is present, stderr is empty, and no failure
marker remains. The terminal resource record reports `14348` seconds,
peak RSS `20,141,236,224` bytes, final workspace `161,962,028` bytes, and
free space `1,903,052,161,024` bytes. All admitted resource ceilings pass.

The output structure contains raw and filtered feature-barcode matrix
directories, HDF5 matrix files, molecule information, and Cell Ranger summary
files. The closure records their names and sizes only; it does not open their
contents. BAM and secondary outputs are absent.

## Firewall and next gate

No matrix contents, expression/UMI values, cell barcodes, QC values, donor
attributes, labels, outcomes, PCA, clustering, persistence landscapes, or
remaining units were opened or processed by this closure. The corrected
landscape contract remains unchanged: separate cell and gene observation
views, separate H0/H1, every consecutive active level, no fixed grid, and no
universal level cap.

The next gate is a separately authorized structural/QC review and remaining
unit decision. Nothing in this closure authorizes labels, outcomes, landscape
calculation, additional units, or deletion.
