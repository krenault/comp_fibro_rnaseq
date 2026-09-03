#!/usr/bin/env bash
# Copy StringTie GTFs from Unity (read-only) and run official prepDE.py.
# Never modifies files on the cluster.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
: "${UNITY_HOST:=unity}"
: "${UNITY_BASE:=/project/pi_rachel_daniels_uri_edu/flexible_homeostasis}"
: "${DATA_DIR:?Set DATA_DIR to the local prepDE output directory}"
: "${PREPDE:=$SCRIPT_DIR/prepDE.py}"
: "${PREPDE_LENGTH:=75}"   # StringTie prepDE default; production Drive run omitted -l
: "${SPECIES:=bactrian_camel bat dolphin dromedary_camel gelada honey_badger human rat rhino rousette seal squirrel whale}"
: "${EXPERIMENTS:=glucose temperature hypoxia}"

LOG="${DATA_DIR}/processing_log_$(date +%Y%m%d_%H%M%S).txt"
TMP="${DATA_DIR}/temp_gtfs"
mkdir -p "$DATA_DIR" "$TMP"

log() { local m="[$(date '+%Y-%m-%d %H:%M:%S')] $*"; echo "$m"; echo "$m" >> "$LOG"; }

pattern_for() {
  case "$1" in
    glucose)     echo "2-5mM|2-5-mM|2p5mM|8mM|8Mm|30mM" ;;
    temperature) echo "32C|41C|-32_|-41_|-32[^0-9]|-41[^0-9]" ;;
    hypoxia)     echo "6H|24H|-6_|-6[^0-9H]|-24_|-24[^0-9H]" ;;
    *) echo "unknown experiment: $1" >&2; exit 1 ;;
  esac
}

normalize_sample() {
  local s
  s="$(basename "$1" | sed 's/_stringtie\.gtf$//')"
  s="$(echo "$s" | sed 's/2-5mM/2.5mM/g; s/2p5mM/2.5mM/g; s/2-5-mM/2.5mM/g')"
  echo "$s" | sed 's/-32_/-32C/g; s/-41_/-41C/g; s/-32$/-32C/g; s/-41$/-41C/g'
}

process_one() {
  local species="$1" experiment="$2"
  local remote="${UNITY_BASE}/${species}/stringtie_november2025_one2one"
  local out="${DATA_DIR}/${species}/${experiment}"
  local list="${out}/sample_list.txt"
  mkdir -p "$out"
  rm -rf "$TMP"; mkdir -p "$TMP"

  log "Processing ${species} / ${experiment}"
  local pat files file_count downloaded filename local_file sample
  pat="$(pattern_for "$experiment")"
  files="$(ssh "$UNITY_HOST" "ls ${remote}/*.gtf 2>/dev/null | grep -iE '${pat}'" 2>/dev/null || true)"
  if [[ -z "$files" ]]; then
    log "  no GTF files; skip"
    return 0
  fi
  file_count="$(echo "$files" | wc -l | tr -d ' ')"
  log "  found ${file_count} GTFs"
  : > "$list"
  downloaded=0
  for remote_file in $files; do
    filename="$(basename "$remote_file")"
    local_file="${TMP}/${filename}"
    sample="$(normalize_sample "$filename")"
    if scp -q "${UNITY_HOST}:${remote_file}" "$local_file" 2>/dev/null; then
      printf "%s\t%s\n" "$sample" "$local_file" >> "$list"
      downloaded=$((downloaded + 1))
    else
      log "  WARNING: failed ${filename}"
    fi
  done
  log "  downloaded ${downloaded}/${file_count}"
  [[ "$downloaded" -eq 0 ]] && return 0

  if python3 "$PREPDE" -i "$list" -g "${out}/gene_count_matrix.csv" \
      -t "${out}/transcript_count_matrix.csv" -l "$PREPDE_LENGTH" >> "$LOG" 2>&1; then
    log "  SUCCESS: $(wc -l < "${out}/gene_count_matrix.csv" | tr -d ' ') genes"
  else
    log "  ERROR: prepDE.py failed"
  fi
  rm -rf "$TMP"; mkdir -p "$TMP"
}

log "Unity ${UNITY_HOST}:${UNITY_BASE}"
log "DATA_DIR=${DATA_DIR}  PREPDE_LENGTH=${PREPDE_LENGTH}"
ssh -o ConnectTimeout=10 "$UNITY_HOST" "echo Connection_OK" >/dev/null
log "Unity connection OK"

for species in $SPECIES; do
  for experiment in $EXPERIMENTS; do
    process_one "$species" "$experiment"
  done
done
rm -rf "$TMP"
log "DONE. Matrices under ${DATA_DIR}/<species>/<experiment>/gene_count_matrix.csv"
