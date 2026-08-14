# MV6-F stage-two remediated-root admission

| Field | Result |
|---|---|
| Date | 2026-08-14 |
| Prospective rebind commit | `5762d0e` |
| Queue root | `f5471633e21d229eeabecadf12989dece2a3a7ab5b5d09f4584b0c3b6410bb5d` |
| Remediated implementation root | `5a1258e87eb30c367648daa4bb02d1aec1f4b40dd3799a91fbb7067e9558d292` |
| Rust SHA-256 | `51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d` |
| Rebind validation | 8/8 categories pass |
| Admission | `stage2_authorized = TRUE` |
| Labels/outcomes | Closed / zero |

## Result

The corrected-root maximum group completed in 337.429 seconds with
3,100,360,704 bytes peak process-tree RSS. Its clean repeat completed in
349.147 seconds with 3,029,950,464 bytes peak RSS. Both remain well inside the
prospective 1,200-second and 6-GiB stage-one caps.

The old and remediated roots are scientifically equivalent:

- all 180 cell/gene PH diagram matrices are identical;
- manifest sample/view roles, diagram hashes, and interval counts are
  identical;
- all 6,500 separate cell-H0/cell-H1/gene-H0/gene-H1 exact landscape rows are
  identical;
- the three remediated scientific files repeat byte-for-byte;
- 12/12 canonical R oracles pass, with maximum error
  `5.00222085975111e-12`;
- 12/12 grouped-Persim oracles pass, with maximum Rust error
  `8.185452315956354e-12` and R error `3.183231456205249e-12`;
- oracle depth spans 25 through 1,856 active levels without a cap; and
- all five remediated group files survive zero-rebuild resume with unchanged
  hashes, sizes, and modification times.

The added PCA standardization identity therefore changes provenance/cache
identity as intended but does not change the maximum-group scientific values.
The landscape contract remains H0/H1 separate, finite-positive intervals,
all active consecutive levels, zero padding, exact or error-controlled L2, and
no universal grid or level cap.

## Caught pre-admission defect

The first rebind attempt exposed an unquoted WSL path in the new resume checker
before its child runner started. No production artifact changed. That complete
root was invalidated and quarantined; the checker was corrected and committed,
and every primary/repeat/oracle/resume gate was rerun under the replacement
root. This admission refers only to `5a1258e8…8d292`.

## Decision boundary

The exact remaining 74 queue rows may now execute serially under the committed
live-cap/checkpoint monitor. Each group remains label closed and contains no
fusion, clustering, or outcome job. Any identity, partial-state, numerical,
resource, or storage failure stops new launches and preserves validated groups
without automatic retry.

This admission does not authorize blocked fusion, clustering, outcome access,
public Rust adoption, binary release, new data, or manuscript claims. Complete
75-group structural/numerical validation and immutable resume must pass before
the next MV-06 sprint.

Machine-readable evidence is in
`docs/audits/mv06f-stage2-rebind-evidence/`; private diagrams, logs, interval
staging, caches, binaries, and quarantined outputs remain excluded from Git.
