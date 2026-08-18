#!/usr/bin/env python3
"""
Wrapper to run feature extraction for external validation datasets
for a single RUN_ID. Processes datasets sequentially.
"""
import os
import sys
import subprocess

RUN_ID = os.environ.get("RUN_ID", "0")
PYTHON = sys.executable
CODE = os.path.join(os.path.dirname(__file__), "code_external_validation.py")

# (name, input_chunks_dir, output_features_dir)
DATASETS = [
    ("humsavar", "external_validation/processing/chunks/external_humsavar", "external_validation/processing/features/humsavar"),
    ("mave_independent", "external_validation/processing/chunks/mave_independent", "external_validation/processing/features/mave_independent"),
    ("gnomad_benign", "external_validation/processing/chunks/gnomad_benign", "external_validation/processing/features/gnomad_benign"),
]

for name, inp, out in DATASETS:
    # Skip if input chunks don't exist yet
    if not os.path.exists(inp):
        print(f"\n⚠️  Skipping {name}: input chunks not found at {inp}")
        continue

    log_dir = f"external_validation/logs/feature_extraction/{name}"
    os.makedirs(log_dir, exist_ok=True)
    log_file = f"{log_dir}/run_{int(RUN_ID)+1}.log"
    
    print(f"\n{'='*60}")
    print(f"RUN_ID={RUN_ID} | Dataset={name}")
    print(f"Input: {inp}")
    print(f"Output: {out}")
    print(f"Log: {log_file}")
    print(f"{'='*60}")
    
    env = os.environ.copy()
    env["INPUT_BASE"] = inp
    env["OUTPUT_BASE"] = out
    env["RUN_ID"] = RUN_ID
    
    with open(log_file, "w") as fh:
        result = subprocess.run(
            [PYTHON, CODE],
            env=env,
            stdout=fh,
            stderr=subprocess.STDOUT,
        )
    
    print(f"Completed {name} with exit code {result.returncode}")

print(f"\nRUN_ID={RUN_ID} ALL DONE")
