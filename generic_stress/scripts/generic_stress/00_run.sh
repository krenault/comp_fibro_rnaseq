#!/usr/bin/env bash
# 00 — run the generic-stress path.
# Ingest (01) is skipped unless --ingest is passed (slow; already on disk).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
if [[ -f "$ROOT/.env" ]]; then
  set -a
  # shellcheck disable=SC1091
  source "$ROOT/.env"
  set +a
fi
export FLEX_ROOT="${FLEX_ROOT:-$ROOT}"

if [[ "${1:-}" == "--ingest" ]]; then
  Rscript scripts/generic_stress/01_ingest.R
fi
Rscript scripts/generic_stress/02_public_hmp.R
Rscript scripts/generic_stress/03_drive_holdout.R
Rscript scripts/generic_stress/04_plots.R
Rscript scripts/generic_stress/05_species_hmp.R
