#!/usr/bin/env python3
"""
Wrapper to run AlphaGenome feature extraction for the PanelApp benign benchmark.
Processes a single RUN_ID.
"""
import os
import sys
import subprocess

RUN_ID = os.environ.get("RUN_ID", "0")
PYTHON = sys.executable
CODE = os.path.join(os.path.dirname(__file__), "code_external_generic.py")

INPUT_BASE = "external_validation/processing/chunks/panelapp_benign"
OUTPUT_BASE = "external_validation/processing/features/panelapp_benign"

if not os.path.exists(INPUT_BASE):
    print(f"\n⚠️  Input chunks not found at {INPUT_BASE}")
    sys.exit(1)

log_dir = "external_validation/logs/feature_extraction/panelapp_benign"
os.makedirs(log_dir, exist_ok=True)
log_file = f"{log_dir}/run_{int(RUN_ID)+1}.log"

print(f"\n{'='*60}")
print(f"RUN_ID={RUN_ID} | Dataset=panelapp_benign")
print(f"Input: {INPUT_BASE}")
print(f"Output: {OUTPUT_BASE}")
print(f"Log: {log_file}")
print(f"{'='*60}")

env = os.environ.copy()
env["INPUT_BASE"] = INPUT_BASE
env["OUTPUT_BASE"] = OUTPUT_BASE
env["RUN_ID"] = RUN_ID

with open(log_file, "w") as fh:
    result = subprocess.run(
        [PYTHON, CODE],
        env=env,
        stdout=fh,
        stderr=subprocess.STDOUT,
    )

print(f"Completed panelapp_benign RUN_ID={RUN_ID} with exit code {result.returncode}")
