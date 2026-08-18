#!/bin/bash
# Run AlphaGenome feature extraction for Grimm2015 benchmark in parallel (20 runs).
# Each run uses a different API key (RUN_ID 0-19).

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
cd "$PROJECT_ROOT"

INPUT_BASE="external_validation/processing/chunks/external_grimm2015"
OUTPUT_BASE="external_validation/processing/features/grimm2015"
LOG_DIR="external_validation/logs/feature_extraction/grimm2015"

mkdir -p "$LOG_DIR"

PYTHON=$(which python3)
SCRIPT="external_validation/scripts/code_external_validation.py"

pids=()
for RUN_ID in $(seq 0 19); do
    LOG_FILE="$LOG_DIR/run_$((RUN_ID+1)).log"
    echo "Launching RUN_ID=$RUN_ID -> $LOG_FILE"
    (
        export RUN_ID
        export INPUT_BASE
        export OUTPUT_BASE
        "$PYTHON" "$SCRIPT" > "$LOG_FILE" 2>&1
    ) &
    pids+=($!)
done

echo "Launched ${#pids[@]} parallel feature extraction runs"
echo "PIDs: ${pids[@]}"

# Wait for all background jobs
wait

echo "All Grimm2015 feature extraction runs complete"
