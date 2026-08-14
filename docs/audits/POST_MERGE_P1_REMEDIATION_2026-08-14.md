# Post-merge P1 remediation audit

Date: 2026-08-14

Repository: `jcdaneshmand/scPHcompare`

Base: `d0192d35a4ab52006aa83b0ad3b0ad6a19f066cb` (`main`)

Branch: `agent/post-merge-p1-remediation`

Source review: <https://github.com/jcdaneshmand/scPHcompare/pull/112>

## Scope

Codex submitted four P1 review threads at `2026-08-14T02:44:37Z`, after
PR 112 had merged. This remediation treats all four threads as actionable and
keeps them ahead of MV6-F fusion stage 2.

| Review finding | Remediation | Regression evidence |
|---|---|---|
| Raw cache identity was accepted from size plus expected-hash syntax without hashing the file | `verify_frozen_file_identity()` now validates size and actual SHA-256 before `readRDS()`; the verified digest is propagated into the job bundle | A different same-size file is rejected by SHA-256 |
| SCT queue checked elapsed and RSS caps only after worker exit | `mv05d0_enforce_live_process_caps_v1()` samples every live poll, closes admission, and calls `kill_tree()` immediately on the first cap breach | Fake live workers are terminated independently for elapsed and RSS breaches; a below-cap worker remains alive |
| Analytical-fixture eligibility could be changed to scientific without invalidating the view | `validate_topology_view()` derives the only permitted eligibility state from the validated contract profile | A fixture view forged to `scientific_eligible = TRUE` is rejected before PH |
| Shared PCA accepted sources with different fitted standardization identities | `standardization_id` is now required across PCA sources, carried into PCA provenance and identity, and checked again at projection | Mixed-standardization fitting and projection are both rejected |

## Additional baseline repair

The complete source-tree test suite exposed one unrelated merged-main omission:
`test-mv05p-production-validation.R` required the public, label-closed
`mv05o-validation-plan-2026-08-10.csv`, but the CSV was absent from `main`.
The exact 15-record plan from historical prefreeze commit `46765a7` is restored.
It contains no private data or biological outcomes.

## Local verification

- All changed R and test files parsed successfully.
- Focused remediation, dual-view topology, resource-safe execution,
  provenance, and MV5-P production-validation tests passed.
- The complete source-tree R suite passed in 164.9 seconds.
- Two established optional skips remained: the absent optional Rust prototype
  library and public audit documents intentionally excluded from R builds.
- No unexpected test failures or warnings remained.

The complete suite used the established Ubuntu R 4.4 dependency library with
renv autoloading disabled for the isolated worktree. Hosted GitHub Actions are
the remaining independent environment check.

## Scientific boundary

This remediation changes integrity validation and process safety only. It does
not change the persistence-landscape definition: H0 and H1 remain separate,
all active consecutive landscape levels remain included, and exact or
error-controlled integration remains required. No outcome labels were opened,
no clustering jobs were executed, and no biological conclusions were produced.
