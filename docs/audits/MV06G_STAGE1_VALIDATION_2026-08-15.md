# MV6-G maximum-group stage-one validation and admission audit

## Disposition

The corrected-root maximum-group sentinel passes every prospective scientific,
deterministic, numerical, resource, projection, and immutable-resume gate.
MV6-G stage one completes as `pass_stage2_label_closed_only`. A separately
committed full-production policy may authorize the remaining 74 groups, but
metadata reading, label opening, outcome evaluation, advanced fusion,
clustering, defaults, release, and claims remain closed.

## Accepted identities and scope

| Item | Accepted value |
|---|---|
| Corrected launch commit | `4f90d2e` |
| Parent fusion contract | `f72326c0411c6c954ffb570fbf3e019adb7b09f82b74cd96984409e762077f8b` |
| Queue root | `f5471633e21d229eeabecadf12989dece2a3a7ab5b5d09f4584b0c3b6410bb5d` |
| Stage-one implementation root | `6a76a11d1b211fcf89fddcd67b7591161619950023420c0bfccbbdc65b76ce82` |
| Rust library | `51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d` |
| Group | `mv06f_group_v1:fbed6ad0…b6ebf1` |
| Fold / seed | `large_loso_v1:SRA779509` / `20260807` |

The runner consumed the accepted 180 PH records and 6,500 query component
rows. It did not rerun PH. It generated exactly 2,080 unordered training pairs,
8,320 cell/gene H0/H1 component rows, four training medians, and 14,625
label-closed ranking rows for 1,625 query pairs and nine fixed methods.

## Corrected execution and resources

| Run | Elapsed | Peak process-tree RSS | Private root | Result |
|---|---:|---:|---:|---|
| Primary | 221.397 s | 166,313,984 B | 8,064,765 B | Complete |
| Clean repeat | 227.835 s | 166,301,696 B | 8,064,765 B | Complete |

Both runs are far below the 1,800-second, 12-GiB RSS, and 5-GiB private caps.
Using the slower repeat for every one of 75 groups yields a conservative
4.746553-worker-hour projection. Linear private-storage projection is about
0.605 GB. Both pass their frozen 12-hour and 5-GiB gates.

The first execution attempt under root `883bbe32…16e2` had already passed the
same scientific/resource/repeat/R/Persim gates, but its resume checker split a
spaced WSL path before the child runner could validate reuse. All five files
were unchanged. That complete attempt is preserved in ignored quarantine and
is not part of this acceptance. The corrected attempt above was regenerated
from clean primary and repeat roots after commit `4f90d2e`.

## Training-only scales and rankings

| Component | Training median scale | Reconstruction error |
|---|---:|---:|
| Cell H0 | 30.4780478412656 | 0 |
| Cell H1 | 1.40300829721247 | 0 |
| Gene H0 | 1.39737575116798 | 5.11e-15 |
| Gene H1 | 0.0301507224656368 | 0 |

Every scale used all 2,080 training-to-training distances in its component and
zero query rows or labels. The independent validator reconstructed all 14,625
normalized distances and deterministic ranks across all nine methods. Maximum
formula error was `1.47e-14`, below the maximum `7.09e-14` tolerance; every
65-sample rank block matched ascending distance followed by canonical training
sample ID.

The three scientific artifacts repeat byte-identically:

| Artifact | SHA-256 |
|---|---|
| Training distances | `a97f38c82dacfb86380aafc9cec0406ffa110079b5046b2e9a9db473d24ba095` |
| Four scales | `e5458c6875990ec9349d2baac61aeccc6e5aa21c13d0b1970d3d015e031afac4` |
| Nine-method rankings | `ee84d39d3f9c6db359406a741cc205d7d14bfa95cd5e5a77ab856335481f5dd7` |

## Landscape numerical oracles

The prospectively selected 12 rows span minimum, median, and maximum combined
interval depth separately in cell H0, cell H1, gene H0, and gene H1. They
exercise 239–5,696 combined finite intervals and 19–1,856 active landscape
levels. No level cap or uniform grid was used.

- Canonical R passes 12/12; maximum Rust-versus-R squared-distance error is
  `5.46e-12`, below the largest certified tolerance `4.38e-7`.
- Grouped Persim passes 12/12; maximum squared error is `7.28e-12` versus Rust
  and `9.09e-12` versus R.

The Persim oracle used an isolated Ubuntu Python 3.10.12 environment with the
same accepted scientific package versions: Persim 0.3.8, NumPy 2.2.6, and
SciPy 1.15.2. The prior MV6-F oracle used Python 3.12.11; this interpreter-minor
difference is disclosed, and direct agreement with both production Rust and
canonical R is the acceptance criterion.

## Independent validation and resume

All 10 independent categories pass: source/implementation identity, primary
and repeat artifacts, scientific byte repeat, scale reconstruction, ranking
reconstruction, R oracles, resource caps, full projection, and the
label/downstream firewall. The corrected runner then validated reuse and
returned zero; SHA-256, byte size, and modification time remained unchanged
for all five files.

| Public evidence | SHA-256 |
|---|---|
| Independent validation | `18c74376b8fd96c9e54d1b10845a1ccb1adc1349f7e659ee5e625cfbe1a5c48f` |
| Scientific repeat | `d3076b2728f708a8180c3d78df007bb2aa6c0ffaea5623bfe33d5b9e09caf393` |
| Scale reconstruction | `1f743cc54c58c7d3652ca976fa34fd71cd933e5850f62671901785e16897e726` |
| Ranking reconstruction | `ba47098837d4d5aac69d5e9ffce54115563f1256169e3b6ac9049e0c61b43aec` |
| R oracles | `e23711e09f6da8b3c13a0efc850afff472d7b18f409c1754a44d236918f753cf` |
| Persim oracles | `0f22720a7a0049f1682171bd9b9dc528cd0cf036d885c257275e65f93ab5d2db` |
| Immutable resume | `5140795201777b45039d0bebdef942ad544f2f1be143930cbdea139810ae308f` |

No metadata or outcome label was read, biological outcomes remain zero, and
no fusion evaluation occurred. The next sprint may only prefreeze the serial
74-group label-closed production policy and its complete validation/resume
gates before execution.
