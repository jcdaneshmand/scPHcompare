#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: acquire_mv08x_rust_toolchain.sh PRIVATE_ROOT INSTALLER_URL INSTALLER_SHA256" >&2
  exit 64
fi

private_root=$1
installer_url=$2
expected_installer_sha=$3
case "${expected_installer_sha}" in
  [0-9a-f][0-9a-f][0-9a-f][0-9a-f]*) ;;
  *) echo "invalid installer SHA-256" >&2; exit 65 ;;
esac
if [[ ${#expected_installer_sha} -ne 64 ]]; then
  echo "invalid installer SHA-256 length" >&2
  exit 65
fi

mkdir -p "${private_root}"
installer="${private_root}/rustup-init"
partial="${installer}.partial"
cargo_home="${private_root}/toolchain/cargo-home"
rustup_home="${private_root}/toolchain/rustup-home"
receipt="${private_root}/toolchain-acquisition-receipt.txt"

if [[ -e "${partial}" ]]; then
  echo "refusing stale partial installer" >&2
  exit 66
fi
if [[ ! -f "${installer}" ]]; then
  curl --fail --location --silent --show-error --retry 3 \
    --output "${partial}" "${installer_url}"
  observed=$(sha256sum "${partial}" | awk '{print $1}')
  if [[ "${observed}" != "${expected_installer_sha}" ]]; then
    echo "installer SHA-256 mismatch" >&2
    exit 67
  fi
  mv "${partial}" "${installer}"
  chmod 0700 "${installer}"
fi
observed=$(sha256sum "${installer}" | awk '{print $1}')
if [[ "${observed}" != "${expected_installer_sha}" ]]; then
  echo "cached installer SHA-256 mismatch" >&2
  exit 67
fi

if [[ ! -x "${cargo_home}/bin/rustup" ]]; then
  CARGO_HOME="${cargo_home}" RUSTUP_HOME="${rustup_home}" \
    "${installer}" -y --no-modify-path --profile minimal \
    --default-toolchain 1.97.1
fi
export CARGO_HOME="${cargo_home}"
export RUSTUP_HOME="${rustup_home}"
export PATH="${CARGO_HOME}/bin:${PATH}"
rustup component add --toolchain 1.97.1 rustfmt clippy

rustc_identity=$(rustc +1.97.1 -vV)
cargo_identity=$(cargo +1.97.1 --version)
rustup_identity=$(rustup --version | head -n 1)
host=$(printf '%s\n' "${rustc_identity}" | awk '/^host:/ {print $2}')
release=$(printf '%s\n' "${rustc_identity}" | awk '/^release:/ {print $2}')
if [[ "${host}" != "x86_64-unknown-linux-gnu" || "${release}" != "1.97.1" ]]; then
  echo "installed Rust identity differs from MV8-X" >&2
  exit 68
fi

tmp_receipt="${receipt}.partial"
{
  printf 'installer_sha256=%s\n' "${observed}"
  printf 'rustc_release=%s\n' "${release}"
  printf 'rustc_host=%s\n' "${host}"
  printf 'cargo_identity=%s\n' "${cargo_identity}"
  printf 'rustup_identity=%s\n' "${rustup_identity}"
  printf 'profile=minimal_plus_rustfmt_clippy\n'
  printf 'home_mutation=false\n'
  printf 'system_package_mutation=false\n'
} > "${tmp_receipt}"
mv "${tmp_receipt}" "${receipt}"
printf 'MV8-X isolated Rust toolchain ready: %s %s\n' "${release}" "${host}"
