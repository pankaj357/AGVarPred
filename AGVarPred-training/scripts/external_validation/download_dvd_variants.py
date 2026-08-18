#!/usr/bin/env python3
"""Download per-gene variant CSVs from the Deafness Variation Database (DVD)."""
import os
import time
import argparse
import urllib.request
from pathlib import Path


def download_gene(gene, out_dir, base_url="https://deafnessvariationdatabase.org", delay=0.5):
    url = f"{base_url}/variants/export-csv?gene={gene}&build=GRCH38"
    out_path = out_dir / f"{gene}.csv"
    if out_path.exists():
        return True, f"{gene}: already exists"
    try:
        urllib.request.urlretrieve(url, out_path)
        time.sleep(delay)
        return True, f"{gene}: downloaded"
    except Exception as e:
        return False, f"{gene}: ERROR {e}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-list", required=True, help="File with one gene symbol per line")
    parser.add_argument("--out-dir", required=True, help="Output directory for CSVs")
    parser.add_argument("--delay", type=float, default=0.5, help="Seconds between requests")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.gene_list) as f:
        genes = [line.strip() for line in f if line.strip()]

    for gene in genes:
        ok, msg = download_gene(gene, out_dir, delay=args.delay)
        print(msg)


if __name__ == "__main__":
    main()
