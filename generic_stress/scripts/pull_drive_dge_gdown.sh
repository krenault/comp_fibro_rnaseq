#!/usr/bin/env bash
# Alternative if folder is shared as "Anyone with the link" (Viewer).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${ROOT}/drive_dge"
VENV="${ROOT}/.venv_gdown"
mkdir -p "$OUT"
[[ -x "$VENV/bin/gdown" ]] || { python3 -m venv "$VENV"; "$VENV/bin/pip" install -q -U gdown; }
"$VENV/bin/gdown" --folder "https://drive.google.com/drive/folders/1-cCvd4b5SyubpE2uWA3Kl3Eg4sMb3huo" -O "$OUT"
find "$OUT" -type f | head -50
