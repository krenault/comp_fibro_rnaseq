#!/usr/bin/env bash
# One-time Google Drive OAuth for Flexible_homeostasis GEO comparison.
# Opens a browser; sign in with the Google account that owns/shared the DGE folder.
set -euo pipefail
REMOTE="${RCLONE_REMOTE:-gdrive_fh}"
if rclone listremotes | grep -qx "${REMOTE}:"; then
  echo "Remote ${REMOTE}: already configured."
else
  echo "Creating rclone remote ${REMOTE}: (browser OAuth)..."
  rclone config create "$REMOTE" drive \
    scope=drive.readonly \
    config_is_local=false
fi
echo "Testing access..."
rclone about "${REMOTE}:" | head -5
echo "OK. Next: bash scripts/pull_drive_dge.sh"
