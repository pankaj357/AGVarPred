#!/usr/bin/env python3
"""Run all external validation scoring for both regularized and no_af models."""

import subprocess
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

SCRIPT_DIR = Path(__file__).parent
ROOT = SCRIPT_DIR.parent.parent

BENCHMARKS = [
    "humsavar",
    "mave_independent",
    "gnomad_benign",
    "vip",
    "grimm2015",
    "dvd",
]

def run_scoring(benchmark, model_type):
    script_name = f"score_{benchmark}_{model_type}.py" if model_type == "no_af" else f"score_{benchmark}.py"
    script_path = SCRIPT_DIR / script_name
    
    if not script_path.exists():
        return (benchmark, model_type, "MISSING_SCRIPT", "")
    
    cmd = [sys.executable, str(script_path)]
    
    log_dir = ROOT / "external_validation/logs/scoring"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{benchmark}_{model_type}.log"
    
    try:
        with open(log_file, "w") as f:
            result = subprocess.run(
                cmd, stdout=f, stderr=subprocess.STDOUT, timeout=1800, cwd=str(ROOT)
            )
        return (benchmark, model_type, result.returncode, str(log_file))
    except subprocess.TimeoutExpired:
        return (benchmark, model_type, "TIMEOUT", str(log_file))
    except Exception as e:
        return (benchmark, model_type, f"ERROR: {e}", str(log_file))

def main():
    tasks = []
    for bench in BENCHMARKS:
        for model in ["regularized", "no_af"]:
            tasks.append((bench, model))
    
    print(f"Running {len(tasks)} scoring tasks with max 3 concurrent...")
    
    completed = 0
    with ProcessPoolExecutor(max_workers=3) as executor:
        futures = {executor.submit(run_scoring, bench, model): (bench, model) 
                   for bench, model in tasks}
        
        for future in as_completed(futures):
            bench, model = futures[future]
            try:
                result = future.result()
                completed += 1
                status = result[2]
                if status == 0:
                    print(f"  ✅ [{completed}/{len(tasks)}] {bench}/{model}")
                else:
                    print(f"  ❌ [{completed}/{len(tasks)}] {bench}/{model}: {status}")
            except Exception as e:
                completed += 1
                print(f"  ❌ [{completed}/{len(tasks)}] {bench}/{model}: {e}")
    
    print("\nAll scoring complete!")

if __name__ == "__main__":
    main()
