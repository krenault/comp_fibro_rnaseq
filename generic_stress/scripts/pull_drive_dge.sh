#!/usr/bin/env bash
# Sync DGE CSVs from the shared Drive folder into drive_dge/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REMOTE="${RCLONE_REMOTE:-gdrive_fh}"
FOLDER_ID="1-cCvd4b5SyubpE2uWA3Kl3Eg4sMb3huo"
OUT="${ROOT}/drive_dge"
mkdir -p "$OUT"

if ! rclone listremotes | grep -qx "${REMOTE}:"; then
  echo "Remote ${REMOTE}: not found. Run: bash scripts/setup_gdrive_rclone.sh" >&2
  exit 1
fi

echo "Syncing Drive folder ${FOLDER_ID} -> ${OUT}"
rclone sync "${REMOTE}:" "$OUT" \
  --drive-root-folder-id "$FOLDER_ID" \
  --progress \
  --exclude ".DS_Store"

echo "Downloaded:"
find "$OUT" -type f | head -50
echo "n_files=$(find "$OUT" -type f | wc -l | tr -d ' ')"
