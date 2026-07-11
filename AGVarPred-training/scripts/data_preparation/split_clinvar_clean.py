import pandas as pd
from sklearn.model_selection import train_test_split

print("🔄 Loading clean dataset...")

df = pd.read_csv("clinvar_gold_grch38_clean.csv", low_memory=False)

print("Total variants:", len(df))
print("Unique genes:", df['GeneSymbol'].nunique())

# =========================
# STEP 1: GET UNIQUE GENES AS LIST
# =========================
unique_genes = df['GeneSymbol'].dropna().unique().tolist()
print(f"\nUnique genes available: {len(unique_genes)}")

# =========================
# STEP 2: SPLIT GENES (80/10/10)
# =========================
# First split: 80% train, 20% temp (cal + test)
train_genes, temp_genes = train_test_split(
    unique_genes,
    test_size=0.2,
    random_state=42
)

# Second split: 50% of temp = 10% cal, 50% = 10% test
cal_genes, test_genes = train_test_split(
    temp_genes,
    test_size=0.5,
    random_state=42
)

print(f"\n🎯 GENE-LEVEL SPLIT:")
print(f"  Train genes: {len(train_genes)} ({len(train_genes)/len(unique_genes)*100:.1f}%)")
print(f"  Cal genes:   {len(cal_genes)} ({len(cal_genes)/len(unique_genes)*100:.1f}%)")
print(f"  Test genes:  {len(test_genes)} ({len(test_genes)/len(unique_genes)*100:.1f}%)")

# =========================
# STEP 3: CREATE DATAFRAMES BASED ON GENES
# =========================
train_df = df[df['GeneSymbol'].isin(train_genes)]
cal_df = df[df['GeneSymbol'].isin(cal_genes)]
test_df = df[df['GeneSymbol'].isin(test_genes)]

print(f"\n📊 VARIANT COUNTS:")
print(f"  Train: {len(train_df):,} variants")
print(f"  Cal:   {len(cal_df):,} variants")
print(f"  Test:  {len(test_df):,} variants")

# =========================
# STEP 4: VERIFY NO GENE OVERLAP
# =========================
train_genes_set = set(train_genes)
cal_genes_set = set(cal_genes)
test_genes_set = set(test_genes)

print(f"\n✅ VERIFICATION (No gene overlap expected):")
print(f"  Train ∩ Cal: {len(train_genes_set & cal_genes_set)} genes")
print(f"  Train ∩ Test: {len(train_genes_set & test_genes_set)} genes")
print(f"  Cal ∩ Test: {len(cal_genes_set & test_genes_set)} genes")

# =========================
# STEP 5: CHECK LABEL DISTRIBUTION
# =========================
print(f"\n📊 LABEL DISTRIBUTION:")
print(f"  Train - Benign: {train_df['label'].value_counts(normalize=True).get(0, 0)*100:.1f}%, Pathogenic: {train_df['label'].value_counts(normalize=True).get(1, 0)*100:.1f}%")
print(f"  Cal   - Benign: {cal_df['label'].value_counts(normalize=True).get(0, 0)*100:.1f}%, Pathogenic: {cal_df['label'].value_counts(normalize=True).get(1, 0)*100:.1f}%")
print(f"  Test  - Benign: {test_df['label'].value_counts(normalize=True).get(0, 0)*100:.1f}%, Pathogenic: {test_df['label'].value_counts(normalize=True).get(1, 0)*100:.1f}%")

# =========================
# STEP 6: SAVE FILES
# =========================
train_df.to_csv("train.csv", index=False)
cal_df.to_csv("cal.csv", index=False)
test_df.to_csv("test.csv", index=False)

print("\n✅ Split saved: train.csv, cal.csv, test.csv")
print("\n📁 These splits are GENE-BASED with NO overlap!")