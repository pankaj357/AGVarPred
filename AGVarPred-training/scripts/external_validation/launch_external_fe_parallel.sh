#!/bin/bash
# Launch 20 parallel feature extraction jobs for external validation
# Each RUN_ID uses a different API key

set -e

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$PROJECT_ROOT"

echo "🚀 Launching 20 parallel feature extraction jobs..."
echo ""

for i in $(seq 0 19); do
    export RUN_ID=$i
    LOG="external_validation/logs/feature_extraction/parallel/run_${i}.log"
    mkdir -p "external_validation/logs/feature_extraction/parallel"
    
    echo "   Starting RUN_ID=$i -> $LOG"
    python external_validation/scripts/run_all_external_fe.py > "$LOG" 2>&1 &
done

echo ""
echo "✅ All 20 jobs launched in background."
echo "   Monitor with: tail -f external_validation/logs/feature_extraction/parallel/run_*.log"
echo "   Wait for completion with: wait"
wait
echo ""
echo "🎉 All feature extraction jobs completed!"
