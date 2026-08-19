# MV8-H Structural/QC Admission

Decision: `metadata_only_structural_qc_admit`

This gate inspected the completed HCA_BM_002 Cell Ranger output structure, HDF5 shapes/dtypes, and aggregate Cell Ranger QC metadata.
Expression values, barcode identifiers, labels, outcomes, landscapes, remaining units, and deletion paths were not opened or authorized.

- Eight validation checks: 8/8 passed.
- Filtered matrix columns and aggregate estimated cells: 5037.
- Filtered median genes per cell: 801; total genes detected: 18,845.
- Next gate: owner decision on whether to authorize matrix-content/QC review.
