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
  shift || true
fi

# Reviewer leave-one-series: 4-day tendon hypoxia off. Does not overwrite primary.
if [[ "${1:-}" == "--sens-328685" ]]; then
  export HMP_OUT="${FLEX_ROOT}/results/mammalian_generic_stress/sens_no_GSE328685"
  export DROP_GSE=GSE328685
  export RESUME_MIXED=1
  export LONG_RDS="${FLEX_ROOT}/results/mammalian_generic_stress/cell_meta_ec_long.rds"
  export SKIP_GSEA=1
  unset RESUME_SIGS SKIP_FIBRO FORCE_GSEA GENERIC_RES
  mkdir -p "$HMP_OUT"
  Rscript scripts/generic_stress/02_public_hmp.R
  export GENERIC_RES="$HMP_OUT"
  Rscript scripts/generic_stress/03_drive_holdout.R
  unset GENERIC_RES
  Rscript scripts/generic_stress/07_sensitivity_compare.R
  exit 0
fi

Rscript scripts/generic_stress/02_public_hmp.R
Rscript scripts/generic_stress/03_drive_holdout.R
Rscript scripts/generic_stress/04_plots.R
# 02 refits mixed models unless RESUME_MIXED=1 (and RESUME_SIGS=1 to reuse class rma).
# 05 recomputes unless SKIP_SPECIES_HMP=1.
Rscript scripts/generic_stress/05_species_hmp.R
# Univariate REML + cell-type Jaccard. Genome-wide rma.mv only if RUN_RMAMV=1.
if [[ "${SKIP_REML:-}" != "1" ]]; then
  Rscript scripts/generic_stress/06_reml_celltype.R
fi
