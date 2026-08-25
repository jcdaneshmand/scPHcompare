#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "usage: run_mv08zq_rust_acceptance.sh TOOLCHAIN_ROOT CANDIDATE EXPECTED_CANDIDATE_SHA EXPECTED_SOURCE_SHA EXECUTION_HEAD OUTPUT_ROOT" >&2
  exit 64
fi

toolchain_root=$1
candidate=$2
expected_candidate_sha=$3
expected_source_sha=$4
execution_head=$5
output_root=$6
if [[ ! "${expected_candidate_sha}" =~ ^[0-9a-f]{64}$ ||
      ! "${expected_source_sha}" =~ ^[0-9a-f]{64}$ ||
      ! "${execution_head}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "invalid MV8-ZQ hash binding" >&2
  exit 65
fi
if [[ -e "${output_root}" ]]; then
  echo "refusing to overwrite MV8-ZQ Rust acceptance evidence" >&2
  exit 66
fi

manifest=rust/scph_landscape_kernel/Cargo.toml
source_file=rust/scph_landscape_kernel/src/lib.rs
candidate_sha=$(sha256sum "${candidate}" | awk '{print $1}')
source_sha=$(sha256sum "${source_file}" | awk '{print $1}')
if [[ "${candidate_sha}" != "${expected_candidate_sha}" ||
      "${source_sha}" != "${expected_source_sha}" ]]; then
  echo "MV8-ZQ candidate or corrected source drift" >&2
  exit 67
fi

export CARGO_HOME="${toolchain_root}/toolchain/cargo-home"
export RUSTUP_HOME="${toolchain_root}/toolchain/rustup-home"
export PATH="${CARGO_HOME}/bin:${PATH}"
export CARGO_BUILD_JOBS=1
if [[ ! -x "${CARGO_HOME}/bin/cargo" ]]; then
  echo "isolated Cargo is absent" >&2
  exit 68
fi
rustc_release=$(rustc +1.97.1 -vV | awk '/^release:/ {print $2}')
rustc_host=$(rustc +1.97.1 -vV | awk '/^host:/ {print $2}')
if [[ "${rustc_release}" != "1.97.1" ||
      "${rustc_host}" != "x86_64-unknown-linux-gnu" ]]; then
  echo "isolated Rust identity drift" >&2
  exit 69
fi

mkdir -p "${output_root}"
export CARGO_TARGET_DIR="${output_root}/target"
timeout 3600 cargo +1.97.1 fmt --manifest-path "${manifest}" -- --check \
  > "${output_root}/cargo-fmt.log" 2>&1
timeout 3600 cargo +1.97.1 test --manifest-path "${manifest}" --locked -j 1 \
  > "${output_root}/cargo-test.log" 2>&1
timeout 3600 cargo +1.97.1 clippy --manifest-path "${manifest}" --locked \
  --all-targets --all-features -j 1 -- -D warnings \
  > "${output_root}/cargo-clippy.log" 2>&1

cat > "${output_root}/rust-validation.csv" <<EOF
"contract_id","execution_head","rustc_release","rustc_host","build_jobs","format_passed","unit_tests_passed","unit_tests_total","clippy_passed","candidate_sha256","corrected_source_sha256","scientific_engine_version","production_landscape_jobs","comparison_jobs","clustering_jobs","fusion_jobs","label_jobs","outcome_jobs","outcome_label_state","biological_outcomes_computed"
"mv08zq_rust_validation_v1","${execution_head}","${rustc_release}","${rustc_host}",1,TRUE,TRUE,6,TRUE,"${candidate_sha}","${source_sha}",2,0,0,0,0,0,0,"closed",FALSE
EOF

find "${output_root}" -type f ! -path '*/target/*' ! -name 'artifact-manifest.csv' ! -name '*.partial' -printf '%P\n' | sort |
  while IFS= read -r relative; do
    file="${output_root}/${relative}"
    printf '"%s","%s","%s"\n' "${relative}" \
      "$(stat --printf='%s' "${file}")" "$(sha256sum "${file}" | awk '{print $1}')"
  done | { printf '"artifact","bytes","sha256"\n'; cat; } \
  > "${output_root}/artifact-manifest.csv.partial"
mv "${output_root}/artifact-manifest.csv.partial" "${output_root}/artifact-manifest.csv"
printf 'MV8-ZQ Rust acceptance: format=pass tests=6/6 clippy=pass engine=2\n'
