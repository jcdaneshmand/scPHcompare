# MV8-U full-PH production prefreeze

**Date:** 2026-08-23

**Result:** 23/23 prospective checks pass; no PH or landscape was executed.

## Exact scope

MV8-T already closed 23 of MV8-R's 1,280 PH records. The exact remaining queue is 1,257 gene-view records: 1,236 internal and 21 external. All eight cell-view records are already closed. The remaining queue contains 625 SCT-data and 632 Pearson-residual records across 131 accepted source caches.

## Evidence-based execution order and resources

Execution is serial and label-blind: 14 common-475 records first, then 625 exact-500 residual records, then 618 exact-500 SCT-data records. Every remaining stratum was exercised by the closed sentinel. The measured-max projection is 15.25 hours; the two-times planning projection is 30.49 hours inside a 72-hour aggregate stop. The measured-max output projection is 164 MiB inside a 1-GiB private-output cap.

Ripserr remains primary with an 1,800-second and 8-GiB child cap. One exact GUDHI attempt is permitted only after a Ripserr RSS-cap stop, with an 1,800-second and 12-GiB cap. A timeout or ordinary child failure stops the run without retry.

## Resume and closure

Outputs are atomic and resume accepts only a fully rehashed completed prefix. Ambiguous mid-job evidence is preserved for a separate recovery prefreeze. MV8-W must reopen and reconstruct all 1,257 new records, rerun every H0 MST oracle, independently rehash them, and bind the 23 closed MV8-T records for exact 1,280/1,280 coverage.

## Scientific firewall

H0 and H1 remain separate. Landscapes retain the dissertation-aligned all-active-level, no-fixed-grid, no-level-cap, exact/error-controlled streamed definition, but no landscape is authorized here. Comparisons, clustering, fusion, labels, outcomes, adoption, and manuscript claims remain closed.
