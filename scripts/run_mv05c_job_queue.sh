#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: run_mv05c_job_queue.sh REPOSITORY_ROOT" >&2
  exit 2
fi

repo_root="$1"
cd "$repo_root"

cache="tmp/mv05c/mv05c-existing-data-raw-cache.rds"
sample_manifest="docs/audits/mv05c-label-closed-sample-manifest-2026-08-06.csv"
cell_manifest="docs/audits/mv05c-matched-cell-manifest-2026-08-06.csv"
private_dir="tmp/mv05c/jobs"
audit_dir="docs/audits/mv05c-jobs"
resource_dir="tmp/mv05c/resources"
mkdir -p "$private_dir" "$audit_dir" "$resource_dir"

folds=("SRA716608" "SRA760933")
seeds=(20260805 20260806 20260807 20260808 20260809)
for study in "${folds[@]}"; do
  fold_id="mv05c_loso_v1:${study}"
  safe_job_id="mv05c_loso_v1_${study}"
  for seed in "${seeds[@]}"; do
    status_path="${audit_dir}/mv05c-job-status-${safe_job_id}__${seed}.csv"
    if [[ -f "$status_path" ]]; then
      echo "skip completed artifact: ${fold_id} seed ${seed}"
      continue
    fi
    echo "start: ${fold_id} seed ${seed}"
    /usr/bin/time -v -o "${resource_dir}/${safe_job_id}__${seed}.txt" \
      Rscript scripts/run_mv05c_existing_data_job.R \
      "$cache" "$sample_manifest" "$cell_manifest" "$fold_id" "$seed" \
      "$private_dir" "$audit_dir"
    echo "finish: ${fold_id} seed ${seed}"
  done
done
