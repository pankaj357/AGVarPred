#!/usr/bin/env bash
set -euo pipefail

RELEASE_URL="https://github.com/pankaj357/AGVarPred/releases/download/runtime-gnomad-af-vep-r2.1.1-grch38-v1"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="${GNOMAD_DATA_DIR:-$REPO_ROOT/data/gnomad}"

TSV="gnomad.exomes.r2.1.1.sites.liftover_grch38.af_vep.tsv.gz"
INDEX="${TSV}.tbi"

mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

curl --fail --location --retry 3 --retry-delay 2 \
  --silent --show-error \
  "${RELEASE_URL}/SHA256SUMS.txt" \
  --output SHA256SUMS.txt.part

mv SHA256SUMS.txt.part SHA256SUMS.txt

# Reuse existing files only if both match the published checksums.
if [[ -s "$TSV" && -s "$INDEX" ]] &&
   sha256sum --check --status SHA256SUMS.txt; then
    echo "Existing gnomAD files pass checksum verification."
    sha256sum --check SHA256SUMS.txt
    exit 0
fi

echo "Downloading reduced gnomAD TSV..."
curl --fail --location --retry 3 --retry-delay 2 \
  --retry-all-errors --continue-at - \
  --silent --show-error \
  "${RELEASE_URL}/${TSV}" --output "${TSV}.part"

echo "Downloading tabix index..."
curl --fail --location --retry 3 --retry-delay 2 \
  --retry-all-errors --continue-at - \
  --silent --show-error \
  "${RELEASE_URL}/${INDEX}" --output "${INDEX}.part"

mv "${TSV}.part" "$TSV"
mv "${INDEX}.part" "$INDEX"

echo "Verifying downloaded files..."
sha256sum --check SHA256SUMS.txt

echo "gnomAD data download and verification completed."
