#!/usr/bin/env bash
set -euo pipefail

Rscript --vanilla scripts/monitor_mv05f_mapping_pilot.R \
  docs/audits/mv05g-integrated-coordinate-manifest-2026-08-08.csv \
  tmp/mv05d1/production-cache-v2 \
  docs/audits/mv05d0-individual-source-raw-shards-2026-08-07.csv \
  tmp/mv05d0/raw-sample-cache-v2 \
  docs/audits/mv05d0-sct-cache-resources-v2-2026-08-07.csv \
  tmp/mv05d0/sct-cache-v2-production \
  tmp/mv05g/production/results \
  tmp/mv05g/production/group-audits \
  tmp/mv05g/production/logs \
  docs/audits/mv05g-integrated-coordinate-resources-2026-08-08.csv
