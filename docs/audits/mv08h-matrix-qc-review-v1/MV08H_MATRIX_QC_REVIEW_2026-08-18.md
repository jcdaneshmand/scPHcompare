# MV8-H Matrix/QC Review

Decision: `matrix_qc_review_pass`
Unit: `HCA_BM_002`

The approved review opened only the filtered/raw feature-barcode H5 matrices and aggregate Cell Ranger metrics.
Barcode identifiers, labels, outcomes, PCA, persistence homology, landscapes, other units, and deletion paths remained closed.

- Validation: 8/8 passed.
- Filtered matrix: 5037 cells, 33563 features, 4517625 nonzero entries.
- Frozen legacy QC depth rule: 4614 cells pass the 384-cell minimum.
- Median detected features/cell: 801.000000; median counts/cell: 2545.000000.
- No topology or biological interpretation was performed.
- Next gate: owner decision on topology review or stop.
