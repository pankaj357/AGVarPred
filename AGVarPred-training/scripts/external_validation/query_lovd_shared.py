#!/usr/bin/env python3
"""
Query the LOVD3 shared instance for expert-curated pathogenic variants
in genes NOT present in the training set.

LOVD3 shared API: https://databases.lovd.nl/shared/api
"""

import json
import requests
import pandas as pd
from pathlib import Path
from time import sleep

ROOT = Path(__file__).resolve().parents[3]
TRAIN_GENES_PATH = ROOT / "external_validation/train_genes_upper.json"
OUTDIR = ROOT / "external_validation/lovd_data"
OUTDIR.mkdir(exist_ok=True, parents=True)

LOVD_API = "https://databases.lovd.nl/shared/api/rest.php/genes"
LOVD_VARIANTS_API = "https://databases.lovd.nl/shared/api/rest.php/variants"
DELAY = 0.5


def load_train_genes():
    with open(TRAIN_GENES_PATH) as f:
        return set(json.load(f))


def fetch_lovd_genes():
    """Fetch all genes available in the LOVD shared instance."""
    print("🔄 Fetching LOVD gene list...")
    try:
        resp = requests.get(LOVD_API, timeout=60)
        resp.raise_for_status()
        data = resp.json()
        genes = []
        for entry in data:
            gene_symbol = entry.get("id", "").upper()
            if gene_symbol:
                genes.append({
                    "gene": gene_symbol,
                    "name": entry.get("name", ""),
                    "chromosome": entry.get("chromosome", ""),
                    "url": entry.get("url", ""),
                })
        print(f"✅ LOVD genes: {len(genes):,}")
        return genes
    except Exception as e:
        print(f"❌ Error fetching LOVD genes: {e}")
        return []


def fetch_variants_for_gene(gene_symbol):
    """Fetch variants for a specific gene from LOVD."""
    url = f"{LOVD_VARIANTS_API}/{gene_symbol}"
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        data = resp.json()
        variants = []
        for v in data:
            # LOVD returns various fields; extract what we can
            variant_id = v.get("id", "")
            dna_change = v.get("VariantOnGenome/DNA", "")
            classification = v.get("VariantOnGenome/Classification", "")
            # Only keep clearly pathogenic
            if classification and "pathogenic" in classification.lower():
                variants.append({
                    "gene": gene_symbol,
                    "variant_id": variant_id,
                    "dna_change": dna_change,
                    "classification": classification,
                    "raw": json.dumps(v),
                })
        return variants
    except Exception as e:
        # 404 or empty means no variants for this gene
        return []


def main():
    train_genes = load_train_genes()
    print(f"Training genes: {len(train_genes):,}")

    all_genes = fetch_lovd_genes()
    if not all_genes:
        print("❌ No genes fetched. Exiting.")
        return

    df_genes = pd.DataFrame(all_genes)
    df_genes["in_training"] = df_genes["gene"].isin(train_genes)
    novel_genes_df = df_genes[~df_genes["in_training"]].copy()
    novel_genes = novel_genes_df["gene"].tolist()

    print(f"\n📊 LOVD GENE ANALYSIS")
    print(f"   Total LOVD genes: {len(df_genes):,}")
    print(f"   In training: {df_genes['in_training'].sum():,}")
    print(f"   Novel genes: {len(novel_genes_df):,}")

    df_genes.to_csv(OUTDIR / "lovd_all_genes.csv", index=False)
    novel_genes_df.to_csv(OUTDIR / "lovd_novel_genes.csv", index=False)

    # Fetch variants for novel genes (limit to first 200 to avoid timeout)
    print(f"\n🔄 Fetching pathogenic variants for novel genes (limit: 200 genes)...")
    all_variants = []
    limit = min(200, len(novel_genes))
    for i, gene in enumerate(novel_genes[:limit]):
        variants = fetch_variants_for_gene(gene)
        all_variants.extend(variants)
        if (i + 1) % 20 == 0:
            print(f"   Processed {i+1}/{limit} genes | Variants found: {len(all_variants)}")
        sleep(DELAY)

    if all_variants:
        df_var = pd.DataFrame(all_variants)
        df_var.to_csv(OUTDIR / "lovd_pathogenic_variants_novel_genes.csv", index=False)
        print(f"✅ Saved: {OUTDIR / 'lovd_pathogenic_variants_novel_genes.csv'}")
        print(f"   Total pathogenic variants: {len(df_var):,}")
        print(f"   Genes with variants: {df_var['gene'].nunique():,}")
    else:
        print("⚠ No pathogenic variants found in queried genes.")

    print(f"\n📝 NEXT STEPS:")
    print(f"   1. Review {OUTDIR / 'lovd_novel_genes.csv'} for relevant genes")
    print(f"   2. Increase the 'limit' variable in this script to query more genes")
    print(f"   3. Map DNA changes (e.g., c.100A>G) to GRCh38 coordinates using VEP")
    print(f"   4. Filter for SNVs only and apply variant-level deduplication")


if __name__ == "__main__":
    main()
