#!/usr/bin/env python3
"""Download per-gene variant CSVs from DVD in parallel with timeout."""
import os
import sys
import time
import argparse
import urllib.request
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


def download_gene(gene, out_dir, base_url="https://deafnessvariationdatabase.org", timeout=180):
    url = f"{base_url}/variants/export-csv?gene={gene}&build=GRCH38"
    out_path = out_dir / f"{gene}.csv"
    if out_path.exists() and out_path.stat().st_size > 100:
        return True, f"{gene}: already exists"
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=timeout) as response:
            data = response.read()
        with open(out_path, "wb") as f:
            f.write(data)
        size_kb = len(data) / 1024
        return True, f"{gene}: downloaded ({size_kb:.1f} KB)"
    except Exception as e:
        if out_path.exists():
            try:
                out_path.unlink()
            except Exception:
                pass
        return False, f"{gene}: ERROR {type(e).__name__}: {e}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-list", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=180)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.gene_list) as f:
        genes = [line.strip() for line in f if line.strip()]

    print(f"Downloading {len(genes)} genes with {args.workers} workers, timeout={args.timeout}s...", flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures = {ex.submit(download_gene, g, out_dir, timeout=args.timeout): g for g in genes}
        for future in as_completed(futures):
            ok, msg = future.result()
            print(msg, flush=True)

    n_done = sum(1 for p in out_dir.glob("*.csv") if p.stat().st_size > 100)
    print(f"\nDone: {n_done}/{len(genes)} downloaded.", flush=True)


if __name__ == "__main__":
    main()
