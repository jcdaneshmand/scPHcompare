#!/usr/bin/env bash
set -euo pipefail

workflow=.github/workflows/rust-accelerator-certification.yml
candidate_section=$(mktemp /tmp/mv05bg-guard.XXXXXX)
case "${candidate_section}" in
  /tmp/mv05bg-guard.*) ;;
  *) exit 91 ;;
esac
trap 'rm -f -- "${candidate_section}"' EXIT

sed '/^  no-release-guard:/,$d' "${workflow}" > "${candidate_section}"
! grep -Eq '(^|/)release-action|gh release|create-release|contents: write' \
  "${candidate_section}"
! grep -Eq 'tmp/mv05|private-|docs/private|\.rds|\.pdf' \
  "${candidate_section}"
forbidden=$(git ls-files | grep -E \
  '(^docs/private/|\.pdf$|^example_run\.r$|^tmp/|^rust-candidate-evidence/)' \
  || true)
test -z "${forbidden}"
printf 'Hardened no-release/private-data guard: PASS\n'
