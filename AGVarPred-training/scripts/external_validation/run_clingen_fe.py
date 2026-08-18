#!/usr/bin/env python3
"""
Run AlphaGenome feature extraction for ClinGen benchmark.
Single RUN_ID job — processes one run directory.
"""
import os
import sys

RUN_ID = os.environ.get("RUN_ID", "0")
PYTHON = sys.executable
CODE = os.path.join(os.path.dirname(__file__), "code_external_validation.py")

LOG_DIR = "external_validation/logs/feature_extraction/clingen"
os.makedirs(LOG_DIR, exist_ok=True)
log_file = f"{LOG_DIR}/run_{int(RUN_ID)+1}.log"

print(f"\n{'='*60}")
print(f"RUN_ID={RUN_ID} | Dataset=clingen")
print(f"Log: {log_file}")
print(f"{'='*60}")

env = os.environ.copy()
env["INPUT_BASE"] = "external_validation/processing/chunks/external_clingen"
env["OUTPUT_BASE"] = "external_validation/processing/features/clingen"
env["RUN_ID"] = RUN_ID

with open(log_file, "w") as fh:
    import subprocess
    result = subprocess.run(
        [PYTHON, CODE],
        env=env,
        stdout=fh,
        stderr=subprocess.STDOUT,
    )

print(f"Completed clingen with exit code {result.returncode}")
