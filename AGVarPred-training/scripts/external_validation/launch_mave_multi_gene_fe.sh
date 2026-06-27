#!/bin/bash
# Launch feature extraction for remaining MAVE multi-gene chunks

cd .
PYTHON=/home/info_lab/miniconda3/envs/alphagenome/bin/python
CODE=external_validation/scripts/code_external_generic.py

export INPUT_BASE=external_validation/processing/chunks/mave_multi_gene
export OUTPUT_BASE=external_validation/processing/features/mave_multi_gene

for RUN_ID in $(seq 1 19); do
    LOG="logs_mave_multi_gene/run_$((RUN_ID+1)).log"
    echo "Launching RUN_ID=$RUN_ID"
    nohup env RUN_ID=$RUN_ID PYTHONUNBUFFERED=1 $PYTHON $CODE > "$LOG" 2>&1 &
    sleep 2
done

echo "All 19 workers launched."
