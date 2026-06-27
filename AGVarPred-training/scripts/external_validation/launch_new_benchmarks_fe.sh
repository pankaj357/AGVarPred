#!/bin/bash
for RUN_ID in $(seq 0 19); do
  export RUN_ID
  python3 external_validation/scripts/run_new_benchmarks_fe.py > "external_validation/logs/feature_extraction/new_benchmarks/run_${RUN_ID}.log" 2>&1 &
done
wait
echo "All new benchmark feature extraction jobs completed!"
