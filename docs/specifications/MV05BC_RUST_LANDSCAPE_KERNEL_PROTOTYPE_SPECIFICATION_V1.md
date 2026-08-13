# MV5-BC bounded Rust landscape-kernel prototype v1

Date: 2026-08-13
Authorization: project owner, following MV5-BB `49fb442`

## Authorized scope

MV5-BC may install the pinned, isolated Rust toolchain below and implement only
the MV5-BB pair-bounded exact squared-L2 persistence-landscape kernel. The
accepted R engine remains canonical and every Rust failure falls back to R.
Prototype success does not authorize production adoption, additional seeds,
partition selection, labels/outcomes, or a change to any public default.

## Toolchain freeze

- Host: Ubuntu 22.04, `x86_64-unknown-linux-gnu`.
- Rust toolchain: `1.97.1-x86_64-unknown-linux-gnu`.
- Installer: official `rustup-init` 1.29.0 from `static.rust-lang.org`.
- Installer SHA-256: `4acc9acc76d5079515b46346a485974457b5a79893cfb01112423c89aeb5aa10`.
- Installation: minimal profile with dedicated ignored MV5-BC `CARGO_HOME` and
  `RUSTUP_HOME`; no distro package, home-directory, or shell-profile mutation.
- Build command: `cargo +1.97.1 build --release --locked`.
- External Rust crates: none. `Cargo.lock` and a source-hash manifest must be
  frozen before the first compilation.

## Prototype boundary and gates

The scientific, ABI, ownership, fallback, determinism, memory, equivalence, and
speed boundaries are exactly those frozen in MV5-BB. The implementation is one
dimension and one pair per call, uses all finite intervals and consecutive
active levels with zero padding, and integrates the exact piecewise-linear
squared difference without a grid or level cap.

MV5-BC evaluates only equivalence Tiers A-C: three analytical fixtures, 20
tractable exact R-oracle results, and the 12 frozen MV5-BA worst-depth
certificates. It must also demonstrate stable error/fallback behavior, two
clean deterministic builds/runs, peak RSS at or below 1 GiB, at least 3x median
speedup over R on the six-pair panel with no slower pair, and exact-staged
package checks with the prototype absent and available. Tiers D-E and a
separate adoption decision remain mandatory before production use.
