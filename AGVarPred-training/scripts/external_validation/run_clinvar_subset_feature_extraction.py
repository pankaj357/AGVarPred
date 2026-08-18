#!/usr/bin/env python3
"""
Run AlphaGenome feature extraction for ClinVar 3+ star SUBSET benchmark.
Spawns 20 parallel processes, one per RUN_ID.
"""

import os
import sys
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "external_validation/scripts/code_external_generic.py"
INPUT_BASE = ROOT / "external_validation/processing/chunks/clinvar_3star_subset"
OUTPUT_BASE = ROOT / "external_validation/processing/features/clinvar_3star_subset"
LOG_DIR = ROOT / "external_validation/logs/clinvar_subset"

N_RUNS = 20


def run_single(run_id: int):
    env = os.environ.copy()
    env["RUN_ID"] = str(run_id)
    env["INPUT_BASE"] = str(INPUT_BASE)
    env["OUTPUT_BASE"] = str(OUTPUT_BASE)

    log_file = LOG_DIR / f"run_{run_id}.log"
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    cmd = [sys.executable, str(SCRIPT)]

    with open(log_file, "w") as lf:
        lf.write(f"RUN_ID={run_id}\n")
        lf.write(f"CMD={' '.join(cmd)}\n")
        lf.write("=" * 60 + "\n")
        lf.flush()

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
            cwd=str(ROOT),
        )

        for line in proc.stdout:
            lf.write(line)
            lf.flush()

        proc.wait()

    return {
        "run_id": run_id,
        "returncode": proc.returncode,
        "log": str(log_file),
    }


def main():
    print(f"🚀 Starting ClinVar SUBSET feature extraction")
    print(f"   Input:  {INPUT_BASE}")
    print(f"   Output: {OUTPUT_BASE}")
    print(f"   Runs:   {N_RUNS}")
    print(f"   Logs:   {LOG_DIR}")
    print()

    start = time.time()

    results = []
    with ProcessPoolExecutor(max_workers=20) as executor:
        futures = {executor.submit(run_single, i): i for i in range(N_RUNS)}
        for future in as_completed(futures):
            run_id = futures[future]
            try:
                res = future.result()
                status = "✅" if res["returncode"] == 0 else "❌"
                print(f"{status} Run {run_id:2d} finished (code={res['returncode']})")
                results.append(res)
            except Exception as e:
                print(f"❌ Run {run_id:2d} crashed: {e}")
                results.append({"run_id": run_id, "returncode": -1, "error": str(e)})

    elapsed = time.time() - start
    print(f"\n⏱ Total time: {elapsed/60:.1f} minutes")

    failed = [r for r in results if r.get("returncode", -1) != 0]
    if failed:
        print(f"\n⚠️ {len(failed)} runs failed:")
        for r in failed:
            print(f"   Run {r['run_id']}: code={r.get('returncode')}")
    else:
        print("\n🎉 All runs completed successfully!")


if __name__ == "__main__":
    main()
