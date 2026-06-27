#!/usr/bin/env python3
"""
Finish fetching missing VEP for benchmarks that were interrupted or skipped.
"""

import sys
sys.path.insert(0, str(__file__).rsplit('/', 1)[0])
from fetch_vep_rest_api import process_benchmark

benchmarks = [
    "humsavar_vep.parquet",
    "mave_independent_vep.parquet",
    "gnomad_benign_vep.parquet",
]

for fname in benchmarks:
    process_benchmark(fname)

print("\n✅ Finished fetching VEP for remaining benchmarks!")
