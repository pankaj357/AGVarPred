#!/usr/bin/env python3
"""
Sequential feature extraction for MAVE multi-gene benchmark.
Runs all remaining chunks one at a time to avoid the hanging issue
with code_external_generic.py.
"""

import os
import sys
import time
import pandas as pd
from pathlib import Path
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

API_KEYS = [
    "AIzaSyDomOPmToBNr-y06A6bY2Pmqh6UJ5_7VLA",
    "AIzaSyDHHygGFg5LyMU00iPaKblvRz-bD2vMvmo",
    "AIzaSyBQbzMnScSPrTb1MI7920bj-xgKPlBChEk",
    "AIzaSyCbl9UbnFcqS7nQlL-UVZdTmbAkdUSbzgU",
    "AIzaSyCbZlR0wdaxILkPrWPjJltPsN1drrKJTCw",
    "AIzaSyD6n7hWocSn-L43DNF0jb2UMjbVOkNQgkg",
    "AIzaSyBm5oPZjvvhSW_4GzbTJ3AEry0TGBClN7U",
    "AIzaSyAr0iOqyy4Qvp_Lzq2KBVzjSJWQEla0QnE",
    "AIzaSyBy-bqASkFAa6N0-sFgRSE9wqBfltaDQi0",
    "AIzaSyCpr_wTZHd_v0oYLr70C_fhWCgLuBfQ8Ec",
    "AIzaSyCSqFoSs5W6k61eN8LGOG7jfBJrExRmppw",
    "AIzaSyAPaS9SE7WO4uw1PIG8GqX-s7BsAPflzkc",
    "AIzaSyCW5SltHTqkjSFxjnH3uuGO6jxUN9I9b8M",
    "AIzaSyBD_g6eVWLR9RcZa97l49n2stN3ZeHL664",
    "AIzaSyA-zr5rVAOPDxQduxFsDRFENIz8Nq_95F4",
    "AIzaSyC2JQLmQznLovT6Xc7snyJH0IuMuzoctXY",
    "AIzaSyBOVUkFsuq0t7ACiYxbb2dXSjg-5OrztgE",
    "AIzaSyC_sf5FZRxkPK20DfcauZM7qPyA39sL-Lw",
    "AIzaSyDU6s91wEB0C5dY4GJNt4IJk2lfRuHApUs",
    "AIzaSyCWSEt1652KiBkHJiLmSsCCQuC3f1_eBdc"
]

ROOT = Path(".")
INPUT_BASE = ROOT / "external_validation/processing/chunks/mave_independent"
OUTPUT_BASE = ROOT / "external_validation" / "processing/features/mave_independent"

def process_chunk(run_id):
    api_key = API_KEYS[run_id]
    input_path = INPUT_BASE / f"run_{run_id + 1}"
    files = [f for f in input_path.glob("*.txt") if not f.name.startswith("._")]
    
    if not files:
        print(f"RUN_ID={run_id}: No files found")
        return True
    
    print(f"\n{'='*60}")
    print(f"RUN_ID={run_id} | API_KEY={api_key[:20]}...")
    print(f"{'='*60}")
    
    print("Creating client...")
    model = dna_client.create(api_key)
    print("Client created")
    
    output_dir = OUTPUT_BASE / f"output_test_run_{run_id + 1}"
    output_dir.mkdir(exist_ok=True, parents=True)
    
    all_scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
    seq_len = dna_client.SEQUENCE_LENGTH_1MB
    
    success = True
    for f in files:
        gene_name = f.stem
        outdir = output_dir / f"{gene_name}_VEP"
        outdir.mkdir(exist_ok=True, parents=True)
        outfile = outdir / f"{gene_name}_ALL_VEP_RAW_SCORE_MATRIX.parquet"
        
        if outfile.exists():
            print(f"  SKIP: {outfile.name} already exists")
            continue
        
        print(f"\n  Processing: {gene_name}")
        variants = pd.read_csv(f, sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "gene"])
        print(f"  Variants: {len(variants)}")
        
        all_rows = []
        all_columns = set()
        errors = 0
        
        for idx, row in variants.iterrows():
            if idx % 50 == 0:
                print(f"    {idx}/{len(variants)}")
            
            try:
                chrom = str(row["chrom"])
                if chrom == "chrMT":
                    chrom = "chrM"
                pos = int(row["pos"])
                ref = str(row["ref"])
                alt = str(row["alt"])
                variant_id = f"{chrom}_{pos}_{ref}_{alt}"
                
                variant_obj = genome.Variant(
                    chromosome=chrom, position=pos, reference_bases=ref,
                    alternate_bases=alt, name=variant_id,
                )
                interval = variant_obj.reference_interval.resize(seq_len)
                
                scores = model.score_variant(
                    interval=interval, variant=variant_obj, variant_scorers=all_scorers,
                    organism=dna_client.Organism.HOMO_SAPIENS,
                )
                df_scores = variant_scorers.tidy_scores(scores)
                
                vector = {"variant_id": variant_id}
                for _, s in df_scores.iterrows():
                    col = f"{s.get('output_type','unk')}__{s.get('biosample_name','unk')}__{s.get('track_name','unk')}"
                    vector[col] = s.get("raw_score", None)
                    all_columns.add(col)
                
                all_rows.append(vector)
                
            except Exception as e:
                errors += 1
                if errors <= 3:
                    print(f"    ERROR: {e}")
                continue
        
        if all_rows:
            df = pd.DataFrame(all_rows).set_index("variant_id")
            df = df.reindex(columns=sorted(all_columns))
            df.to_parquet(outfile)
            df.to_csv(str(outfile).replace(".parquet", ".csv"))
            print(f"  DONE: {len(df)} variants ({errors} errors)")
        else:
            print(f"  FAILED: no variants scored ({errors} errors)")
            success = False
    
    return success

def main():
    print("=" * 70)
    print("MAVE MULTI-GENE FEATURE EXTRACTION (SEQUENTIAL)")
    print("=" * 70)
    
    # Find which chunks are missing
    todo = []
    for run_id in range(20):
        outfile = OUTPUT_BASE / f"output_test_run_{run_id + 1}" / f"chunk_{run_id + 1}_VEP" / f"chunk_{run_id + 1}_ALL_VEP_RAW_SCORE_MATRIX.parquet"
        if outfile.exists():
            print(f"Chunk {run_id + 1}: ALREADY DONE")
        else:
            todo.append(run_id)
    
    print(f"\nChunks to process: {len(todo)} / 20")
    print(f"Run IDs: {todo}")
    
    for run_id in todo:
        start = time.time()
        try:
            ok = process_chunk(run_id)
            elapsed = time.time() - start
            print(f"\n✅ Chunk {run_id + 1} complete in {elapsed:.0f}s")
        except Exception as e:
            elapsed = time.time() - start
            print(f"\n❌ Chunk {run_id + 1} failed after {elapsed:.0f}s: {e}")
    
    print("\n" + "=" * 70)
    print("ALL CHUNKS PROCESSED")
    print("=" * 70)

if __name__ == "__main__":
    main()
