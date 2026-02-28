#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
SNAP_DIR="$ROOT_DIR/.snapshots/src_checkpoint_20260227_213739"

if [[ ! -d "$SNAP_DIR" ]]; then
  echo "snapshot not found: $SNAP_DIR" >&2
  exit 1
fi

# Restore source/config files from the checkpoint.
# Data outputs are intentionally kept.
rsync -a --delete \
  --exclude 'Data/' \
  --exclude 'BasicData/' \
  --exclude '.snapshots/' \
  "$SNAP_DIR"/ "$ROOT_DIR"/

echo "restored source checkpoint: $SNAP_DIR"
