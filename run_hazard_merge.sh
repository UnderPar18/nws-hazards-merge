#!/usr/bin/env bash
# Cron wrapper: runs the merge with a lock so overlapping runs never collide.
# Usage: run_hazard_merge.sh [extra hazard_merge.py args]
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${HAZARD_MERGE_PYTHON:-$DIR/.venv/bin/python}"
LOGDIR="${HAZARD_MERGE_LOGDIR:-$DIR/logs}"
mkdir -p "$LOGDIR"
exec 9>"$DIR/.hazard_merge.lock"
if ! flock -n 9; then
  echo "$(date -u +%FT%TZ) previous run still active, skipping" >> "$LOGDIR/hazard_merge.log"
  exit 0
fi
"$PY" "$DIR/hazard_merge.py" "$@" >> "$LOGDIR/hazard_merge.log" 2>&1
# keep the log from growing forever (~5 MB)
if [ "$(stat -c %s "$LOGDIR/hazard_merge.log")" -gt 5000000 ]; then
  mv "$LOGDIR/hazard_merge.log" "$LOGDIR/hazard_merge.log.1"
fi
