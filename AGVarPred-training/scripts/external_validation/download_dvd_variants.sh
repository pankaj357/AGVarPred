#!/bin/bash
# Download DVD per-gene CSVs using curl with timeout/retry.
GENE_LIST="$1"
OUT_DIR="$2"
DELAY="${3:-0.5}"

mkdir -p "$OUT_DIR"

total=$(wc -l < "$GENE_LIST")
count=0
while IFS= read -r gene; do
    count=$((count + 1))
    out="$OUT_DIR/${gene}.csv"
    if [[ -s "$out" ]]; then
        echo "[$count/$total] $gene: already exists"
        continue
    fi
    echo "[$count/$total] Downloading $gene ..."
    curl -s --max-time 60 --retry 2 \
        "https://deafnessvariationdatabase.org/variants/export-csv?gene=${gene}&build=GRCH38" \
        -o "$out"
    sleep "$DELAY"
done < "$GENE_LIST"

echo "Done. Downloaded $(ls -1 "$OUT_DIR"/*.csv 2>/dev/null | wc -l) files."
