#!/bin/bash
# Launch AlphaGenome feature extraction for PanelApp benign benchmark
mkdir -p external_validation/logs/feature_extraction/panelapp_benign

for RUN_ID in $(seq 0 19); do
  export RUN_ID
  python3 external_validation/scripts/run_panelapp_fe.py > "external_validation/logs/feature_extraction/panelapp_benign/run_${RUN_ID}.log" 2>&1 &
done
wait
echo "All PanelApp benign feature extraction jobs completed!"
