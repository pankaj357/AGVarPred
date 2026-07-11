import pandas as pd

print("🔄 Loading ClinVar file...")

df = pd.read_csv(
    "variant_summary.txt",
    sep="\t",
    dtype={"Chromosome": str},
    low_memory=False
)
print("Total variants:", len(df))

# =========================
# STEP 1: Keep only GRCh38
# =========================
df = df[df["Assembly"] == "GRCh38"]
print("After GRCh38 filter:", len(df))

# =========================
# STEP 2: Keep germline only ← NEW
# =========================
df = df[df["OriginSimple"] == "germline"]
print("After germline filter:", len(df))

# =========================
# STEP 3: Keep high review status ← NEW
# =========================
high_quality_review = [
    "criteria provided, multiple submitters, no conflicts",
    "reviewed by expert panel",
    "practice guideline"
]
df = df[df["ReviewStatus"].isin(high_quality_review)]
print("After review status filter:", len(df))

# =========================
# STEP 4: Keep strict labels only
# =========================
df = df[df["ClinicalSignificance"].isin(["Pathogenic", "Benign"])]
print("After strict label filter:", len(df))

# =========================
# STEP 5: Remove conflicting entries ← NEW
# =========================
df = df[~df["ClinicalSignificance"].str.contains("conflicting", case=False, na=False)]
print("After conflict filter:", len(df))

# =========================
# STEP 6: Remove missing VCF info
# =========================
df = df.dropna(subset=[
    "Chromosome",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF"
])
print("After missing VCF filter:", len(df))

# =========================
# STEP 7: Create label column
# =========================
df["label"] = df["ClinicalSignificance"].map({
    "Pathogenic": 1,
    "Benign": 0
})

# =========================
# STEP 8: Create variant_id
# =========================
df["variant_id"] = (
    df["Chromosome"] + ":" +
    df["PositionVCF"].astype(str) + ":" +
    df["ReferenceAlleleVCF"] + ":" +
    df["AlternateAlleleVCF"]
)

# =========================
# STEP 9: Remove duplicates
# =========================
before = len(df)
df = df.drop_duplicates(subset=["variant_id"])
after = len(df)
print(f"Removed duplicates: {before - after}")

# =========================
# STEP 10: Check class balance
# =========================
print("\n📊 Class distribution:")
print(df["label"].value_counts())
print(f"Pathogenic: {df['label'].sum():,}")
print(f"Benign: {(df['label'] == 0).sum():,}")
ratio = df['label'].sum() / (df['label'] == 0).sum()
print(f"Imbalance ratio: {ratio:.2f}")

# =========================
# STEP 11: Save
# =========================
df.to_csv("clinvar_gold_grch38_clean.csv", index=False)
print(f"\n✅ Clean germline gold dataset saved! Final size: {len(df):,}")