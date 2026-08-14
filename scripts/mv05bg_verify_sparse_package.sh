#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  printf 'Usage: mv05bg_verify_sparse_package.sh <new-output-directory>\n' >&2
  exit 64
fi
repo=$(git rev-parse --show-toplevel)
output=$1
case "${output}" in
  "${repo}"/tmp/mv05bg-sparse-*) ;;
  *)
    printf 'Output must be a new absolute path under repo/tmp/mv05bg-sparse-*\n' >&2
    exit 65
    ;;
esac
test ! -e "${output}"
mkdir -p "${output}/source" "${output}/build" \
  "${output}/accepted" "${output}/current"

paths=(
  .Rbuildignore .Rprofile DESCRIPTION LICENSE NAMESPACE README.md
  R inst man tests
)
git archive --format=tar HEAD -- "${paths[@]}" |
  tar -xf - -C "${output}/source"
cd "${output}/build"
R CMD build "${output}/source" >/dev/null
tarball=$(find . -maxdepth 1 -type f -name 'scPHcompare_*.tar.gz' \
  -print -quit)
test -n "${tarball}"
accepted_tar="${repo}/tmp/mv05bc/package-check-v3/scPHcompare_0.1.0.tar.gz"
test -f "${accepted_tar}"
tar -xzf "${accepted_tar}" -C "${output}/accepted"
tar -xzf "${tarball}" -C "${output}/current"
sed -i '/^Packaged:/d' \
  "${output}/accepted/scPHcompare/DESCRIPTION" \
  "${output}/current/scPHcompare/DESCRIPTION"
find "${output}/accepted/scPHcompare" -type f -exec chmod 0644 {} +
find "${output}/current/scPHcompare" -type f -exec chmod 0644 {} +
file_count=$(find "${output}/current/scPHcompare" -type f | wc -l)
git -C "${output}/accepted/scPHcompare" init -q
git -C "${output}/accepted/scPHcompare" add -A
accepted_tree=$(git -C "${output}/accepted/scPHcompare" write-tree)
git -C "${output}/current/scPHcompare" init -q
git -C "${output}/current/scPHcompare" add -A
current_tree=$(git -C "${output}/current/scPHcompare" write-tree)
test "${accepted_tree}" = "${current_tree}"
test "${file_count}" = 219
printf 'contract_id,accepted_tree,current_tree,file_count,identical\n' \
  > "${output}/sparse-package-invariance.csv"
printf 'mv05bg_sparse_package_v1,%s,%s,%s,true\n' \
  "${accepted_tree}" "${current_tree}" "${file_count}" \
  >> "${output}/sparse-package-invariance.csv"
printf 'MV5-BG sparse R package invariance: PASS (%s; %s files)\n' \
  "${current_tree}" "${file_count}"
