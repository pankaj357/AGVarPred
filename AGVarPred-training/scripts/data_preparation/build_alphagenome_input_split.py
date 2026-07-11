import pandas as pd
import os
import re

def process_file(input_csv, output_dir):

    print(f"\n🔄 Processing {input_csv}...")

    df = pd.read_csv(input_csv, low_memory=False)

    df = df[[
        'Chromosome',
        'PositionVCF',
        'ReferenceAlleleVCF',
        'AlternateAlleleVCF',
        'GeneSymbol'
    ]]

    df.columns = ['chrom','pos','ref','alt','gene']

    # remove invalid rows
    df = df.dropna()
    df = df[df['chrom'] != 'na']

    # keep only valid DNA
    df = df[
        df["ref"].astype(str).str.match("^[ACGTN]+$") &
        df["alt"].astype(str).str.match("^[ACGTN]+$")
    ]

    # add chr prefix
    df['chrom'] = "chr" + df['chrom'].astype(str)

    # remove duplicates
    df['vid'] = df['chrom']+"_"+df['pos'].astype(str)+"_"+df['ref']+"_"+df['alt']
    df = df.drop_duplicates(subset=['vid'])

    print(f"  Input variants: {len(df):,}")
    print(f"  Unique genes: {df['gene'].nunique()}")

    os.makedirs(output_dir, exist_ok=True)

    groups = df.groupby("gene")

    processed = 0
    total_variants_output = 0

    for gene, gdf in groups:

        if len(gdf) < 10:
            continue

        gene_clean = re.sub(r'[^a-zA-Z0-9_\-]', '_', str(gene))[:50]

        # Cap at 2000 variants per gene
        gdf_sampled = gdf.sample(n=min(len(gdf), 2000), random_state=42)

        gdf_sampled[['chrom','pos','ref','alt','gene']].to_csv(
            f"{output_dir}/{gene_clean}.txt",
            sep="\t", index=False, header=False
        )

        processed += 1
        total_variants_output += len(gdf_sampled)

        if processed % 500 == 0:
            print(f"  {input_csv}: {processed} genes processed...")

    print(f"  Output genes: {processed}")
    print(f"  Output variants: {total_variants_output:,}")
    print(f"✅ Done: {input_csv}")


# =========================
# RUN FOR ALL SPLITS (USING NEW GENE-BASED SPLITS)
# =========================
print("\n" + "="*60)
print("CREATING ALPHAGENOME INPUT FILES FROM PROPER GENE SPLITS")
print("="*60)

process_file("train.csv", "alphagenome_input/train")
process_file("cal.csv",   "alphagenome_input/cal")
process_file("test.csv",  "alphagenome_input/test")

print("\n" + "="*60)
print("🎯 ALL INPUT FILES READY - PROPER GENE SPLITS!")
print("="*60)