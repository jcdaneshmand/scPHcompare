#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
  echo "usage: run_mv08x_rust_rebuild.sh TOOLCHAIN_ROOT EVIDENCE_ROOT EXECUTION_HEAD CARGO_TOML_SHA CARGO_LOCK_SHA RUST_SOURCE_SHA ABI_HEADER_SHA" >&2
  exit 64
fi

toolchain_root=$1
evidence=$2
execution_head=$3
expected_manifest=$4
expected_lock=$5
expected_source=$6
expected_header=$7
if [[ ! "${execution_head}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "invalid execution head" >&2
  exit 65
fi
for value in "${expected_manifest}" "${expected_lock}" "${expected_source}" "${expected_header}"; do
  if [[ ! "${value}" =~ ^[0-9a-f]{64}$ ]]; then
    echo "invalid expected source SHA-256" >&2
    exit 65
  fi
done
if [[ -e "${evidence}" ]]; then
  echo "refusing to overwrite Rust rebuild evidence" >&2
  exit 66
fi

repo=$(pwd -P)
manifest="rust/scph_landscape_kernel/Cargo.toml"
lock="rust/scph_landscape_kernel/Cargo.lock"
source_file="rust/scph_landscape_kernel/src/lib.rs"
header="rust/scph_landscape_kernel/include/scph_landscape_kernel.h"
observed_manifest=$(sha256sum "${manifest}" | awk '{print $1}')
observed_lock=$(sha256sum "${lock}" | awk '{print $1}')
observed_source=$(sha256sum "${source_file}" | awk '{print $1}')
observed_header=$(sha256sum "${header}" | awk '{print $1}')
if [[ "${observed_manifest}" != "${expected_manifest}" ||
      "${observed_lock}" != "${expected_lock}" ||
      "${observed_source}" != "${expected_source}" ||
      "${observed_header}" != "${expected_header}" ]]; then
  echo "Rust implementation differs from MV8-X prefreeze" >&2
  exit 67
fi

export CARGO_HOME="${toolchain_root}/toolchain/cargo-home"
export RUSTUP_HOME="${toolchain_root}/toolchain/rustup-home"
export PATH="${CARGO_HOME}/bin:${PATH}"
export CARGO_BUILD_JOBS=1
export SOURCE_DATE_EPOCH=0
if [[ ! -x "${CARGO_HOME}/bin/cargo" ]]; then
  echo "isolated Cargo is absent" >&2
  exit 68
fi
rustc_identity=$(rustc +1.97.1 -vV)
host=$(printf '%s\n' "${rustc_identity}" | awk '/^host:/ {print $2}')
release=$(printf '%s\n' "${rustc_identity}" | awk '/^release:/ {print $2}')
if [[ "${host}" != "x86_64-unknown-linux-gnu" || "${release}" != "1.97.1" ]]; then
  echo "isolated Rust identity differs from MV8-X" >&2
  exit 69
fi

mkdir -p "${evidence}"
printf '%s\n' "${rustc_identity}" > "${evidence}/rustc-identity.txt"
cargo +1.97.1 --version > "${evidence}/cargo-identity.txt"
rustup --version > "${evidence}/rustup-identity.txt"

export CARGO_TARGET_DIR="${evidence}/check-target"
timeout 3600 cargo +1.97.1 fmt --manifest-path "${manifest}" -- --check \
  > "${evidence}/cargo-fmt.log" 2>&1
timeout 3600 cargo +1.97.1 test --manifest-path "${manifest}" --locked -j 1 \
  > "${evidence}/cargo-test.log" 2>&1
timeout 3600 cargo +1.97.1 clippy --manifest-path "${manifest}" --locked \
  --all-targets -j 1 -- -D warnings > "${evidence}/cargo-clippy.log" 2>&1

build_one() {
  local label=$1
  local target="${evidence}/target-${label}"
  CARGO_TARGET_DIR="${target}" /usr/bin/time -v -o "${evidence}/build-${label}-time.txt" \
    timeout 3600 cargo +1.97.1 build --manifest-path "${manifest}" \
      --release --locked -j 1 > "${evidence}/build-${label}.log" 2>&1
  local library="${target}/release/libscph_landscape_kernel.so"
  test -f "${library}"
  sha256sum "${library}" | awk '{print $1}' > "${evidence}/build-${label}.sha256"
  stat --printf='%s\n' "${library}" > "${evidence}/build-${label}.bytes"
}
build_one a
build_one b
hash_a=$(cat "${evidence}/build-a.sha256")
hash_b=$(cat "${evidence}/build-b.sha256")
bytes_a=$(cat "${evidence}/build-a.bytes")
bytes_b=$(cat "${evidence}/build-b.bytes")
if [[ "${hash_a}" != "${hash_b}" || "${bytes_a}" != "${bytes_b}" ]]; then
  echo "clean Rust builds are not byte-identical" >&2
  exit 70
fi
candidate="${evidence}/scph-landscape-kernel-v1_x86_64-unknown-linux-gnu.so"
cp "${evidence}/target-b/release/libscph_landscape_kernel.so" "${candidate}"

cc -std=c11 -Wall -Wextra -Werror \
  -I rust/scph_landscape_kernel/include \
  tests/native/mv05bc_ffi_sanitizer.c \
  -L "${evidence}/target-b/release" -lscph_landscape_kernel \
  -Wl,-rpath,"${evidence}/target-b/release" -lm \
  -o "${evidence}/ffi-harness"
"${evidence}/ffi-harness" > "${evidence}/ffi-harness.log" 2>&1
ldd "${candidate}" > "${evidence}/native-dependencies.txt"

cat > "${evidence}/source-bindings.csv" <<EOF
"role","bytes","sha256"
"cargo_manifest","$(stat --printf='%s' "${manifest}")","${observed_manifest}"
"cargo_lock","$(stat --printf='%s' "${lock}")","${observed_lock}"
"rust_kernel","$(stat --printf='%s' "${source_file}")","${observed_source}"
"c_abi_header","$(stat --printf='%s' "${header}")","${observed_header}"
EOF
prior_hash="51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d"
cat > "${evidence}/build-validation.csv" <<EOF
"contract_id","execution_head","rustc_release","rustc_host","build_jobs","external_crates","format_passed","unit_tests_passed","clippy_passed","clean_builds_byte_identical","build_a_sha256","build_b_sha256","candidate_sha256","candidate_bytes","prior_accepted_sha256","matches_prior_accepted_binary","native_c_abi_passed","published_release","public_default_changed"
"mv08x_rust_rebuild_validation_v1","${execution_head}","${release}","${host}",1,0,TRUE,TRUE,TRUE,TRUE,"${hash_a}","${hash_b}","${hash_b}",${bytes_b},"${prior_hash}",$([[ "${hash_b}" == "${prior_hash}" ]] && echo TRUE || echo FALSE),TRUE,FALSE,FALSE
EOF

find "${evidence}" -type f ! -name 'artifact-manifest.csv' ! -name '*.partial' -printf '%P\n' | sort | \
  while IFS= read -r relative; do
    file="${evidence}/${relative}"
    printf '"%s","%s","%s"\n' "${relative}" \
      "$(stat --printf='%s' "${file}")" "$(sha256sum "${file}" | awk '{print $1}')"
  done | { printf '"artifact","bytes","sha256"\n'; cat; } \
  > "${evidence}/artifact-manifest.csv.partial"
mv "${evidence}/artifact-manifest.csv.partial" "${evidence}/artifact-manifest.csv"
printf 'MV8-X Rust rebuild complete: sha256=%s bytes=%s\n' "${hash_b}" "${bytes_b}"
