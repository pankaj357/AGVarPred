#!/usr/bin/env python3
"""
Local HGVS cDNA to genomic coordinate mapper using refGene.

Maps HGVS notations like NM_000546.6:c.818G>A to genomic coordinates.
Only handles simple substitutions (SNVs) for now.
"""

import pandas as pd
import re
import gzip
from pathlib import Path

REFGENE = Path("external_validation/raw_data/refGene.txt.gz")

# Cache for transcript info
_TRANSCRIPT_CACHE = {}


def load_transcript_info(nm_id):
    """Load exon info for a transcript from refGene."""
    if nm_id in _TRANSCRIPT_CACHE:
        return _TRANSCRIPT_CACHE[nm_id]
    
    # Strip version if present
    nm_base = nm_id.split('.')[0]
    
    df = pd.read_csv(REFGENE, sep="\t", header=None,
                     names=["bin", "name", "chrom", "strand", "txStart", "txEnd",
                            "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
                            "score", "gene", "cdsStartStat", "cdsEndStat", "exonFrames"])
    
    # Match by base name (without version)
    row = df[df["name"] == nm_base]
    if len(row) == 0:
        return None
    
    row = row.iloc[0]
    
    # Parse exons
    starts = [int(x) for x in row["exonStarts"].rstrip(",").split(",") if x]
    ends = [int(x) for x in row["exonEnds"].rstrip(",").split(",") if x]
    
    info = {
        "chrom": row["chrom"].replace("chr", ""),
        "strand": row["strand"],
        "txStart": int(row["txStart"]),
        "txEnd": int(row["txEnd"]),
        "cdsStart": int(row["cdsStart"]),
        "cdsEnd": int(row["cdsEnd"]),
        "exons": list(zip(starts, ends)),
    }
    
    _TRANSCRIPT_CACHE[nm_id] = info
    return info


def map_cdna_to_genomic(nm_id, cdna_pos):
    """Map a cDNA position to genomic coordinate."""
    info = load_transcript_info(nm_id)
    if info is None:
        return None
    
    exons = info["exons"]
    strand = info["strand"]
    
    if strand == "+":
        # Walk left to right through exons
        cum_len = 0
        for start, end in exons:
            exon_len = end - start
            if cum_len + exon_len >= cdna_pos:
                # Position is in this exon
                offset = cdna_pos - cum_len - 1  # 0-based offset within exon
                return start + offset + 1  # Convert to 1-based VCF position
            cum_len += exon_len
    else:
        # Negative strand: walk right to left through exons
        # cDNA position 1 = txEnd (rightmost exon end)
        cum_len = 0
        for start, end in reversed(exons):
            exon_len = end - start
            if cum_len + exon_len >= cdna_pos:
                offset = cdna_pos - cum_len - 1
                return end - offset  # Position decreases going left
            cum_len += exon_len
    
    return None


def parse_hgvs(hgvs_str):
    """Parse HGVS cDNA notation. Returns (nm_id, cdna_pos, ref, alt) or None."""
    # Pattern: NM_12345.6:c.123A>G
    pattern = r"^(NM_\d+(?:\.\d+)?):c\.(\d+)([ACGT])>([ACGT])$"
    m = re.match(pattern, hgvs_str)
    if m:
        return m.group(1), int(m.group(2)), m.group(3), m.group(4)
    return None


def map_hgvs(hgvs_str):
    """Map HGVS notation to genomic coordinates."""
    parsed = parse_hgvs(hgvs_str)
    if parsed is None:
        return None
    
    nm_id, cdna_pos, ref, alt = parsed
    genomic_pos = map_cdna_to_genomic(nm_id, cdna_pos)
    if genomic_pos is None:
        return None
    
    info = load_transcript_info(nm_id)
    return {
        "chrom": info["chrom"],
        "pos": genomic_pos,
        "ref": ref,
        "alt": alt,
        "strand": info["strand"],
    }


if __name__ == "__main__":
    # Test with known variants
    test_cases = [
        "NM_000546.6:c.818G>A",  # TP53 p.R273H, chr17:7577120
        "NM_000314.8:c.800A>G",   # PTEN
        "NM_007294.3:c.5565A>T",  # BRCA1
    ]
    
    for hgvs in test_cases:
        result = map_hgvs(hgvs)
        print(f"{hgvs} -> {result}")
