#!/usr/bin/env python3
"""
Build ClinVar 3+ Star Holdout Benchmark

Filters ClinVar VCF for:
- 3+ star review status (multiple submitters, expert panel, practice guideline)
- Pathogenic/Likely_pathogenic vs Benign/Likely_benign only
- Single nucleotide variants only
- Genes NOT in training data
"""

import gzip
import re
from pathlib import Path
from collections import Counter
import argparse

# Review statuses that qualify as 3+ stars
THREE_PLUS_STAR = {
    "criteria_provided,_multiple_submitters,_no_conflicts",
    "reviewed_by_expert_panel",
    "practice_guideline",
}

# Pathogenic classifications
PATHO_SIGS = {
    "Pathogenic",
    "Likely_pathogenic",
    "Pathogenic/Likely_pathogenic",
    "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance",
}

# Benign classifications
BENIGN_SIGS = {
    "Benign",
    "Likely_benign",
    "Benign/Likely_benign",
}


def load_genes(gene_file):
    with open(gene_file) as f:
        return set(line.strip() for line in f if line.strip())


def parse_info(info_str):
    """Parse VCF INFO field into dict."""
    result = {}
    for pair in info_str.split(";"):
        if "=" in pair:
            k, v = pair.split("=", 1)
            result[k] = v
    return result


def is_snv(ref, alt):
    """Check if ref and alt are single nucleotide."""
    return len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT"


def classify(clnsig):
    """Return 'pathogenic', 'benign', or None."""
    # Handle composite classifications (pipe or slash separated)
    parts = set(re.split(r"[|/]", clnsig))
    # Check for pathogenic first
    for p in parts:
        if p in PATHO_SIGS:
            return "pathogenic"
    for p in parts:
        if p in BENIGN_SIGS:
            return "benign"
    return None


def build_benchmark(vcf_path, gene_file, output_path):
    train_genes = load_genes(gene_file)
    print(f"Loaded {len(train_genes)} training genes")

    pathogenic = []
    benign = []
    skipped = Counter()

    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            info = parse_info(parts[7])

            # Filter SNVs
            if not is_snv(ref, alt):
                skipped["not_snv"] += 1
                continue

            # Filter 3+ star
            revstat = info.get("CLNREVSTAT", "")
            if revstat not in THREE_PLUS_STAR:
                skipped[f"revstat:{revstat}"] += 1
                continue

            # Filter pathogenic/benign only
            clnsig = info.get("CLNSIG", "")
            label = classify(clnsig)
            if label is None:
                skipped[f"clnsig:{clnsig}"] += 1
                continue

            # Get gene
            gene_info = info.get("GENEINFO", "")
            if not gene_info:
                skipped["no_gene"] += 1
                continue

            # GENEINFO format: SYMBOL:ID|SYMBOL2:ID2
            gene = gene_info.split(":")[0].split("|")[0]

            # Exclude training genes
            if gene in train_genes:
                skipped["in_training"] += 1
                continue

            variant = {
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "gene": gene,
                "clnsig": clnsig,
                "clnrevstat": revstat,
            }

            if label == "pathogenic":
                pathogenic.append(variant)
            else:
                benign.append(variant)

    print(f"\nPathogenic: {len(pathogenic)}")
    print(f"Benign: {len(benign)}")
    print(f"\nSkipped reasons (top 20):")
    for reason, count in skipped.most_common(20):
        print(f"  {reason}: {count}")

    # Write output
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        f.write("chrom\tpos\tref\talt\tgene\tlabel\tclnsig\tclnrevstat\n")
        for v in pathogenic + benign:
            f.write(f"{v['chrom']}\t{v['pos']}\t{v['ref']}\t{v['alt']}\t{v['gene']}\t")
            label_str = "1" if v["clnsig"] in PATHO_SIGS or any(p in PATHO_SIGS for p in re.split(r"[|/]", v["clnsig"])) else "0"
            f.write(f"{label_str}\t{v['clnsig']}\t{v['clnrevstat']}\n")

    print(f"\nWrote {len(pathogenic) + len(benign)} variants to {output_path}")
    return len(pathogenic), len(benign)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", default="external_data/clinvar_GRCh38_latest.vcf.gz")
    parser.add_argument("--genes", default="external_validation/processing/all_genes.txt")
    parser.add_argument("--output", default="external_validation/source_datasets/clinvar_3star_holdout.tsv")
    args = parser.parse_args()

    build_benchmark(args.vcf, args.genes, args.output)
