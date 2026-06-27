import os
from pathlib import Path
import shutil

INPUT_DIR = "alphagenome_input/test"
OUTPUT_BASE = "chunks_input_test"

N_CHUNKS = 20

files = list(Path(INPUT_DIR).glob("*.txt"))

os.makedirs(OUTPUT_BASE, exist_ok=True)

for i in range(N_CHUNKS):
    os.makedirs(f"{OUTPUT_BASE}/run_{i+1}", exist_ok=True)

for i, f in enumerate(files):
    chunk_id = i % N_CHUNKS
    dest = f"{OUTPUT_BASE}/run_{chunk_id+1}/{f.name}"
    shutil.copy(f, dest)

print("✅ TEST chunking complete")
