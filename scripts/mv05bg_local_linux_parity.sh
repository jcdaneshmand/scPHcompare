#!/usr/bin/env bash
set -euo pipefail

repo=$(git rev-parse --show-toplevel)
cd "${repo}"
export CARGO_HOME="${repo}/tmp/mv05bc/toolchain/cargo-home"
export RUSTUP_HOME="${repo}/tmp/mv05bc/toolchain/rustup-home"
export CARGO_TARGET_DIR="${repo}/tmp/mv05bg-local/target"
export PATH="${CARGO_HOME}/bin:${PATH}"
evidence="${repo}/tmp/mv05bg-local/evidence"
mkdir -p "${evidence}"

bash scripts/mv05bg_local_guard.sh
Rscript scripts/mv05bg_verify_locked_library.R \
  "${evidence}/dependency-bootstrap-manifest.csv"
rustc +1.97.1 -vV > "${evidence}/rustc-identity.txt"
cargo +1.97.1 fmt --manifest-path \
  rust/scph_landscape_kernel/Cargo.toml -- --check
cargo +1.97.1 test --manifest-path \
  rust/scph_landscape_kernel/Cargo.toml --locked
cargo +1.97.1 clippy --manifest-path \
  rust/scph_landscape_kernel/Cargo.toml --locked --all-targets -- -D warnings
cargo +1.97.1 clean --manifest-path rust/scph_landscape_kernel/Cargo.toml
cargo +1.97.1 build --manifest-path \
  rust/scph_landscape_kernel/Cargo.toml --release --locked
library="${CARGO_TARGET_DIR}/release/libscph_landscape_kernel.so"
first=$(sha256sum "${library}" | cut -d' ' -f1)
cargo +1.97.1 clean --manifest-path rust/scph_landscape_kernel/Cargo.toml
cargo +1.97.1 build --manifest-path \
  rust/scph_landscape_kernel/Cargo.toml --release --locked
second=$(sha256sum "${library}" | cut -d' ' -f1)
test "${first}" = "${second}"
printf '%s' "${first}" > "${evidence}/clean-build-a.sha256"
printf '%s' "${second}" > "${evidence}/clean-build-b.sha256"
candidate="${evidence}/scph-landscape-kernel-v1_x86_64-unknown-linux-gnu.so"
cp "${library}" "${candidate}"

cc -std=c11 -Wall -Wextra -Werror \
  -I rust/scph_landscape_kernel/include \
  tests/native/mv05bc_ffi_sanitizer.c \
  -L "${CARGO_TARGET_DIR}/release" -lscph_landscape_kernel \
  -Wl,-rpath,"${CARGO_TARGET_DIR}/release" -lm \
  -o "${evidence}/ffi-harness"
"${evidence}/ffi-harness"
ldd "${library}" > "${evidence}/native-dependencies-raw.txt"
python3 scripts/mv05bg_normalize_dependencies.py \
  --platform Linux \
  --input "${evidence}/native-dependencies-raw.txt" \
  --output "${evidence}/native-dependencies.txt"
Rscript scripts/mv05bf_rust_ci_fixture.R \
  "${library}" "${evidence}/r-fixtures.csv"

for mode in absent present; do
  source_log="${repo}/tmp/mv05bf-local/evidence/r-check-${mode}.log"
  test -f "${source_log}"
  grep -F 'Status: OK' "${source_log}" >/dev/null
  cp "${source_log}" "${evidence}/r-check-${mode}.log"
done

python3 scripts/mv05bf_build_candidate_manifest.py \
  --target x86_64-unknown-linux-gnu \
  --runner local-ubuntu-22.04 \
  --library "${candidate}" \
  --dependencies "${evidence}/native-dependencies.txt" \
  --fixtures "${evidence}/r-fixtures.csv" \
  --r-check-absent "${evidence}/r-check-absent.log" \
  --r-check-present "${evidence}/r-check-present.log" \
  --output "${evidence}/candidate-manifest.csv"

python3 - "${evidence}/candidate-manifest.csv" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="", encoding="utf-8") as handle:
    row = next(csv.DictReader(handle))
required = (
    "cargo_manifest_sha256",
    "public_abi_header_sha256",
    "r_shim_sha256",
)
assert all(len(row[field]) == 64 for field in required)
assert row["certification_class"] == "candidate-only"
assert row["private_data_used"] == "false"
assert row["published_release"] == "false"
print("MV5-BG local Linux command parity: PASS")
PY
