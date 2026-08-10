#!/usr/bin/env bash
set -euo pipefail

Rscript scripts/monitor_mv05h_ph_groups.R \
  docs/audits/mv05h-integrated-ph-manifest-2026-08-09.csv \
  tmp/mv05g/production/results \
  tmp/mv05h/production/results \
  tmp/mv05h/production/view_audits \
  tmp/mv05h/production/logs \
  docs/audits/mv05h-integrated-ph-resources-2026-08-09.csv \
  75 2 900 4294967296 43200
