# MV7-H landscape-stress closure

Date: 2026-08-16

Status: complete; remaining 19 landscape groups authorized serially

## Production and repeat

The prospectively selected maximum-burden stress group was seed 20260807,
gene topology, H1. It contains 124 samples and all 7,626 unordered pairs. The
accepted all-level Rust kernel completed the primary calculation in 663.772
seconds at 199,913,472 bytes peak process-tree RSS. Its independent repeat
completed in 659.259 seconds at 205,836,288 bytes peak RSS.

Both distance files contain 7,626 finite nonnegative squared-L2 values and are
byte-identical:

- bytes: 4,556,982;
- SHA-256:
  `77eee28eb99cadd58559ab2dcb13c5df1ecef8ea246f837ea5355ea4056b372b`.

Both runs remained far below the frozen 3,600-second and 12-GiB per-group
limits. Every row reports finite intervals only, all active landscape levels,
no level cap, and H1 kept separate from H0.

## Independent validation

After the two documented failed-closed validator rebinds, the clean v3 attempt
passed all 8/8 prospective categories:

- exact frozen group and pair axes;
- complete distance contract;
- accepted Rust-library identity;
- byte-identical scientific repeat;
- minimum-, median-, and maximum-depth independent R oracles;
- closed-form analytical H1 oracle;
- immutable resume with hash, byte, and mtime equality; and
- closed clustering, label, and outcome firewall.

The minimum-depth exact oracle agreed with Rust to
`1.62630325872826e-19`. The median adaptive oracle agreed to
`8.10170832331314e-13`, below its `8.49574431001189e-11` acceptance tolerance.
The maximum adaptive oracle agreed to `8.96927931232083e-13`, below its
`6.18823484598098e-12` acceptance tolerance. The maximum pair contained 3,172
and 2,635 finite H1 intervals and exercised the error-budget-preserving
roundoff bisection path.

## Authorization

The public validation decision authorizes the exact remaining 19 frozen MV7-H
landscape groups serially. It does not authorize dimension combination,
clustering, biological labels, outcomes, result-dependent claims, or any change
to the frozen landscape estimand. Those stages remain closed until the complete
152,520-row H0/H1 corpus passes its own independent validation, representative
repeats, immutable-resume, and resource reconciliation gates.
