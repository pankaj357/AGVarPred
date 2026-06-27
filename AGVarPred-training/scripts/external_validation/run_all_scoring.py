#!/usr/bin/env python3
"""Run external validation scoring for all benchmarks with BOTH models."""

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
SCRIPTS_DIR = ROOT / "external_validation/scripts"

BENCHMARKS = [
    "humsavar",
    "mave_independent",
    "gnomad_benign",
    "grimm2015",
    "clingen",
]

MODELS = [
    ("regularized", ""),
    ("regularized_no_af", "_no_af"),
]

def run_script(script_path, bench, model_name):
    print("\n" + "="*80)
    print(f"RUNNING: {script_path.name}  |  Benchmark: {bench}  |  Model: {model_name}")
    print("="*80)
    result = subprocess.run(
        [sys.executable, str(script_path)],
        cwd=str(ROOT)
    )
    return result.returncode == 0

failed = []
successful = []

for bench in BENCHMARKS:
    for model_name, suffix in MODELS:
        script = SCRIPTS_DIR / f"score_{bench}{suffix}.py"
        if not script.exists():
            print(f"\n⚠ Script not found: {script}")
            failed.append(f"{bench} ({model_name})")
            continue

        ok = run_script(script, bench, model_name)
        if ok:
            print(f"\n✅ SUCCESS: {bench} — {model_name}")
            successful.append(f"{bench} ({model_name})")
        else:
            print(f"\n❌ FAILED: {bench} — {model_name}")
            failed.append(f"{bench} ({model_name})")

print("\n" + "="*80)
print("SCORING RUN COMPLETE")
print("="*80)
print(f"\nSuccessful ({len(successful)}):")
for s in successful:
    print(f"  ✅ {s}")

if failed:
    print(f"\nFailed ({len(failed)}):")
    for f in failed:
        print(f"  ❌ {f}")
else:
    print("\nAll benchmarks scored successfully!")
