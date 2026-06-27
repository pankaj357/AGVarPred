#!/bin/bash
# Launch parallel feature extraction jobs for DVD benchmark

set -e

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$PROJECT_ROOT"

N_JOBS=12

echo "🚀 Launching $N_JOBS parallel DVD feature extraction jobs..."
echo ""

for i in $(seq 0 $((N_JOBS-1))); do
    export RUN_ID=$i
    LOG="external_validation/logs/feature_extraction/dvd/run_${i}.log"
    mkdir -p "external_validation/logs/feature_extraction/dvd"
    
    echo "   Starting RUN_ID=$i -> $LOG"
    python external_validation/scripts/run_dvd_fe.py > "$LOG" 2>&1 &
done

echo ""
echo "✅ All $N_JOBS jobs launched in background."
echo "   Monitor with: tail -f external_validation/logs/feature_extraction/dvd/run_*.log"
echo "   Wait for completion with: wait"
wait
echo ""
echo "🎉 All DVD feature extraction jobs completed!"
