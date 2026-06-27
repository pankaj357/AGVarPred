#!/bin/bash
# Launch all 20 AlphaGenome feature extraction runs for Humsavar benchmark
# Usage: bash run_humsavar_feature_extraction.sh

set -e

echo "🚀 Launching Humsavar feature extraction (20 parallel runs)"
echo "   Output: external_validation/processing/features/humsavar/"
echo ""

mkdir -p external_validation/summaries
mkdir -p external_validation/processing/features/humsavar

for i in $(seq 0 19); do
    echo "  Starting RUN_ID=$i ..."
    RUN_ID=$i nohup /home/info_lab/miniconda3/envs/alphagenome/bin/python -u external_validation/scripts/code_external_humsavar.py >> external_validation/summaries/log_run_$((i+1)).log 2>&1 &
done

echo ""
echo "✅ All 20 runs launched in background."
echo "   Monitor progress: tail -f external_validation/summaries/log_run_*.log"
echo "   Check summaries:  external_validation/summaries/summary_humsavar_run_*.csv"
