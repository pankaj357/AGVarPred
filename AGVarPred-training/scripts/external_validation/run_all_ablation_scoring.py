#!/usr/bin/env python3
"""Run ablation model scoring across all benchmarks."""

import subprocess
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

SCRIPT_DIR = Path(__file__).parent
ROOT = SCRIPT_DIR.parent.parent

MODELS = {
    "af_only": "ablation_feature_groups_output/AF_only_pipeline.pkl",
    "af_plus_vep": "ablation_feature_groups_output/AF_plus_VEP_pipeline.pkl",
    "alphagenome_only": "ablation_feature_groups_output/AlphaGenome_only_pipeline.pkl",
    "vep_only": "ablation_feature_groups_output/VEP_only_pipeline.pkl",
}

BENCHMARKS = [
    "humsavar",
    "mave_independent",
    "gnomad_benign",
    "vip",
    "grimm2015",
    "dvd",
]

# Benchmarks where gnomAD AF should be queried for subgroup analysis
HAS_AF_BENCHES = {"humsavar", "mave_independent", "vip", "grimm2015", "clingen", "dvd"}


def run_scoring(benchmark, model_name, model_path):
    cmd = [
        sys.executable,
        str(SCRIPT_DIR / "score_ablation.py"),
        "--model", str(ROOT / model_path),
        "--benchmark", benchmark,
        "--model-name", model_name,
    ]
    if benchmark in HAS_AF_BENCHES:
        cmd.append("--has-af")

    log_dir = ROOT / "external_validation/logs/scoring"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{benchmark}_{model_name}.log"

    try:
        with open(log_file, "w") as f:
            result = subprocess.run(
                cmd, stdout=f, stderr=subprocess.STDOUT, timeout=1800, cwd=str(ROOT)
            )
        return (benchmark, model_name, result.returncode, str(log_file))
    except subprocess.TimeoutExpired:
        return (benchmark, model_name, "TIMEOUT", str(log_file))
    except Exception as e:
        return (benchmark, model_name, f"ERROR: {e}", str(log_file))


def main():
    tasks = []
    for bench in BENCHMARKS:
        for model_name, model_path in MODELS.items():
            tasks.append((bench, model_name, model_path))

    print(f"Running {len(tasks)} ablation scoring tasks with max 3 concurrent...")

    completed = 0
    with ProcessPoolExecutor(max_workers=3) as executor:
        futures = {executor.submit(run_scoring, bench, model_name, model_path): (bench, model_name)
                   for bench, model_name, model_path in tasks}

        for future in as_completed(futures):
            bench, model_name = futures[future]
            try:
                result = future.result()
                completed += 1
                status = result[2]
                if status == 0:
                    print(f"  ✅ [{completed}/{len(tasks)}] {bench}/{model_name}")
                else:
                    print(f"  ❌ [{completed}/{len(tasks)}] {bench}/{model_name}: {status}")
            except Exception as e:
                completed += 1
                print(f"  ❌ [{completed}/{len(tasks)}] {bench}/{model_name}: {e}")

    print("\nAll ablation scoring complete!")


if __name__ == "__main__":
    main()
