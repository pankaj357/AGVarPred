#!/usr/bin/env python3
"""Score MAVE Independent benchmark with the no-AF regularized model."""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from score_mave_independent import score_mave

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
EXT_DIR = ROOT / "external_validation"
MODEL_PATH = ROOT / "model_6_minus_af_output/Model_1_no_AF_pipeline.pkl"
OUTDIR = EXT_DIR / "results/mave_independent/regularized_no_af"

if __name__ == "__main__":
    score_mave(MODEL_PATH, OUTDIR)
