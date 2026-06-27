# ============================================================
# FINAL UNIFIED PIPELINE
# ============================================================
#
# This script:
#
# 1. Loads train/cal/test ClinVar splits
# 2. Builds:
#       - variant → label
#       - variant → gene
# 3. Extracts from gnomAD VCF:
#       - gnomAD_AF
#       - VEP annotations (SIFT, PolyPhen, LoF, Consequence, IMPACT)
# 4. Processes AlphaGenome parquet files
# 5. Adds:
#       - label
#       - GeneSymbol
#       - gnomAD_AF, log10_gnomAD_AF
#       - AF_missing, is_ultra_rare
#       - VEP features
# 6. Saves final parquet datasets
#
# ============================================================

import pandas as pd
import numpy as np
from glob import glob
import os
import gc
import pysam
import re
from tqdm import tqdm

# ============================================================
# PATHS
# ============================================================

GNOMAD_VCF = os.environ.get(
    "GNOMAD_VCF",
    "external_data/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz",
)

if not os.path.exists(GNOMAD_VCF):
    raise FileNotFoundError(
        f"gnomAD VCF not found: {GNOMAD_VCF}\n"
        "Set GNOMAD_VCF to the correct path or download the file "
        "(see data_manifest.json)."
    )

# ============================================================
# VEP PARSING
# ============================================================

VEP_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
    "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
    "Codons", "Existing_variation", "ALLELE_NUM", "DISTANCE", "STRAND",
    "FLAGS", "VARIANT_CLASS", "MINIMISED", "SYMBOL_SOURCE", "HGNC_ID",
    "CANONICAL", "TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT", "TREMBL",
    "UNIPARC", "GENE_PHENO", "SIFT", "PolyPhen", "DOMAINS", "HGVS_OFFSET",
    "GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF", "EUR_MAF", "SAS_MAF", "AA_MAF",
    "EA_MAF", "ExAC_MAF", "ExAC_Adj_MAF", "ExAC_AFR_MAF", "ExAC_AMR_MAF",
    "ExAC_EAS_MAF", "ExAC_FIN_MAF", "ExAC_NFE_MAF", "ExAC_OTH_MAF",
    "ExAC_SAS_MAF", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME",
    "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "LoF", "LoF_filter",
    "LoF_flags", "LoF_info"
]

VEP_IDX = {f: i for i, f in enumerate(VEP_FIELDS)}


def parse_vep(vep_string):
    """
    Parse a single VEP annotation string.
    Returns dict with extracted features.
    gnomAD VEP may contain multiple transcripts separated by commas.
    We take the first (canonical) annotation.
    """
    if not vep_string:
        return {}

    # Take first transcript annotation
    first_ann = vep_string.split(",")[0]
    parts = first_ann.split("|")

    if len(parts) < 65:
        return {}

    def get(idx):
        val = parts[idx] if idx < len(parts) else ""
        return val if val != "" else None

    # SIFT: "tolerated(0.05)" or "deleterious(0.01)"
    sift_raw = get(VEP_IDX["SIFT"])
    sift_score = None
    sift_pred = None
    if sift_raw:
        m = re.search(r"([^(]+)\(([\d.]+)\)", sift_raw)
        if m:
            try:
                sift_score = float(m.group(2))
            except ValueError:
                pass
            sift_pred = m.group(1)

    # PolyPhen: "probably_damaging(0.999)" or "benign(0.123)"
    polyphen_raw = get(VEP_IDX["PolyPhen"])
    polyphen_score = None
    polyphen_pred = None
    if polyphen_raw:
        m = re.search(r"([^(]+)\(([\d.]+)\)", polyphen_raw)
        if m:
            try:
                polyphen_score = float(m.group(2))
            except ValueError:
                pass
            polyphen_pred = m.group(1)

    consequence = get(VEP_IDX["Consequence"])
    impact = get(VEP_IDX["IMPACT"])
    lof = get(VEP_IDX["LoF"])
    protein_pos = get(VEP_IDX["Protein_position"])
    amino_acids = get(VEP_IDX["Amino_acids"])

    # Encode IMPACT as numeric
    impact_map = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
    impact_score = impact_map.get(impact, 0) if impact else 0

    # Parse protein position (may be range like "123-125")
    protein_pos_num = None
    if protein_pos:
        m = re.search(r"^(\d+)", protein_pos)
        if m:
            try:
                protein_pos_num = int(m.group(1))
            except ValueError:
                pass

    cons_list = consequence.split("&") if consequence else []

    return {
        "vep_SIFT_score": sift_score,
        "vep_SIFT_pred": sift_pred,
        "vep_PolyPhen_score": polyphen_score,
        "vep_PolyPhen_pred": polyphen_pred,
        "vep_IMPACT": impact,
        "vep_IMPACT_score": impact_score,
        "vep_Consequence": consequence,
        "vep_is_missense": 1 if "missense_variant" in cons_list else 0,
        "vep_is_synonymous": 1 if "synonymous_variant" in cons_list else 0,
        "vep_is_stop_gained": 1 if "stop_gained" in cons_list else 0,
        "vep_is_frameshift": 1 if "frameshift_variant" in cons_list else 0,
        "vep_is_splice": 1 if any(c in cons_list for c in ["splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"]) else 0,
        "vep_LoF": lof,
        "vep_is_LoF_HC": 1 if lof == "HC" else 0,
        "vep_Protein_position": protein_pos_num,
        "vep_Amino_acids": amino_acids,
        "vep_has_SIFT": 1 if sift_score is not None else 0,
        "vep_has_PolyPhen": 1 if polyphen_score is not None else 0,
    }


# ============================================================
# LOAD SPLIT CSVs
# ============================================================

print("\n🔄 Loading split CSVs...")


def load_split(csv_file, split_name):

    print(f"\n📂 Loading {split_name}: {csv_file}")

    df = pd.read_csv(csv_file, low_memory=False)

    df['variant_id'] = (
        "chr" + df['Chromosome'].astype(str) + "_" +
        df['PositionVCF'].astype(str) + "_" +
        df['ReferenceAlleleVCF'].astype(str) + "_" +
        df['AlternateAlleleVCF'].astype(str)
    )

    df['label'] = df['ClinicalSignificance'].map({
        'Pathogenic': 1,
        'Benign': 0
    })

    df = df.dropna(subset=['label'])

    label_map = dict(zip(df['variant_id'], df['label']))
    gene_map = dict(zip(df['variant_id'], df['GeneSymbol']))

    variant_set = set(df['variant_id'])

    print(f"✅ {split_name}: {len(df):,} variants")

    return label_map, gene_map, variant_set


train_labels, train_genes, train_vars = load_split(
    "train.csv",
    "TRAIN"
)

cal_labels, cal_genes, cal_vars = load_split(
    "cal.csv",
    "CAL"
)

test_labels, test_genes, test_vars = load_split(
    "test.csv",
    "TEST"
)

# ============================================================
# MERGE ALL LOOKUPS
# ============================================================

print("\n🔄 Merging mappings...")

label_dict = {
    **train_labels,
    **cal_labels,
    **test_labels
}

gene_dict = {
    **train_genes,
    **cal_genes,
    **test_genes
}

all_variants = (
    train_vars |
    cal_vars |
    test_vars
)

print(f"✅ Total unique variants: {len(all_variants):,}")

# ============================================================
# EXTRACT GNOMAD AF + VEP ANNOTATIONS
# ============================================================

print("\n🔄 Extracting gnomAD AFs and VEP annotations...")

vcf = pysam.VariantFile(GNOMAD_VCF)

af_map = {}
vep_map = {}

for rec in tqdm(vcf.fetch()):

    try:

        chrom = rec.chrom

        if not chrom.startswith("chr"):
            chrom = "chr" + chrom

        pos = rec.pos
        ref = rec.ref

        if rec.alts is None:
            continue

        alt = rec.alts[0]

        variant_id = f"{chrom}_{pos}_{ref}_{alt}"

        if variant_id not in all_variants:
            continue

        # AF
        af = rec.info.get("AF")
        if af is None:
            af_value = 0.0
        else:
            af_value = float(af[0])
        af_map[variant_id] = af_value

        # VEP
        vep = rec.info.get("vep")
        if vep:
            vep_map[variant_id] = parse_vep(vep[0])

    except Exception:
        continue

print(f"✅ AFs found for {len(af_map):,} variants")
print(f"✅ VEP annotations found for {len(vep_map):,} variants")

# ============================================================
# BUILD GLOBAL ORDINAL MAPS FOR CATEGORICAL VEP COLUMNS
# ============================================================

print("\n🔄 Building global ordinal encodings for VEP categoricals...")

ORDINAL_COLS = [
    'vep_SIFT_pred', 'vep_PolyPhen_pred', 'vep_IMPACT', 'vep_LoF'
]
ONEHOT_COLS = [
    'vep_Consequence'
]
CAT_VEP_COLS = ORDINAL_COLS + ONEHOT_COLS

# Custom ordinal orders (least severe -> most severe / most benign -> most damaging)
CUSTOM_ORDERS = {
    'vep_IMPACT': ["MODIFIER", "LOW", "MODERATE", "HIGH"],
    'vep_SIFT_pred': ["tolerated", "tolerated_low_confidence", "deleterious_low_confidence", "deleterious"],
    'vep_PolyPhen_pred': ["benign", "possibly_damaging", "probably_damaging", "unknown"],
    'vep_LoF': ["LC", "HC"],
}

# Collect all unique non-null categories across all variants
global_cat_values = {col: set() for col in CAT_VEP_COLS}
for vep_dict in vep_map.values():
    if isinstance(vep_dict, dict):
        for col in CAT_VEP_COLS:
            val = vep_dict.get(col)
            if val is not None:
                global_cat_values[col].add(val)

# Build ordinal mappings
ordinal_maps = {}
for col in CAT_VEP_COLS:
    if col in CUSTOM_ORDERS:
        # Use custom order, keeping only categories that actually exist
        ordered = [c for c in CUSTOM_ORDERS[col] if c in global_cat_values[col]]
        # Append any unexpected categories alphabetically
        extra = sorted(global_cat_values[col] - set(ordered))
        all_cats = ordered + extra
    else:
        # Alphabetical for Consequence and Amino_acids
        all_cats = sorted(global_cat_values[col])

    ordinal_maps[col] = {cat: i for i, cat in enumerate(all_cats)}
    print(f"  {col}: {len(all_cats)} categories")

print("✅ Ordinal maps ready")

# ============================================================
# PROCESS FUNCTION
# ============================================================


def process_split(
    input_pattern,
    output_dir,
    split_name
):

    print(f"\n🚀 PROCESSING {split_name.upper()}")

    os.makedirs(output_dir, exist_ok=True)

    old_files = glob(f"{output_dir}/*.parquet")

    for f in old_files:
        os.remove(f)

    folders = glob(input_pattern)

    folders = [
        f for f in folders
        if not f.endswith("/-_VEP")
    ]

    print(f"📂 Found {len(folders)} folders")

    file_id = 0
    total_variants = 0
    total_matched = 0

    for idx, folder in enumerate(folders):

        gene_name = os.path.basename(folder).replace(
            "_VEP",
            ""
        )

        parquet_files = glob(
            os.path.join(folder, "*.parquet")
        )

        if not parquet_files:
            continue

        try:

            df = pd.read_parquet(parquet_files[0])

            # ====================================================
            # HANDLE VARIANT ID
            # ====================================================

            df = df.reset_index()

            if 'variant_id' not in df.columns:

                if df.columns[0] == 'index':

                    df = df.rename(
                        columns={'index': 'variant_id'}
                    )

                else:

                    df = df.rename(
                        columns={
                            df.columns[0]: 'variant_id'
                        }
                    )

            df['variant_id'] = (
                df['variant_id']
                .astype(str)
            )

            total_variants += len(df)

            # ====================================================
            # FILTER MATCHED VARIANTS
            # ====================================================

            matched_df = df[
                df['variant_id'].isin(label_dict)
            ].copy()

            if len(matched_df) == 0:
                continue

            # ====================================================
            # ADD LABEL
            # ====================================================

            matched_df['label'] = (
                matched_df['variant_id']
                .map(label_dict)
                .astype('int8')
            )

            # ====================================================
            # ADD GENE SYMBOL
            # ====================================================

            matched_df['GeneSymbol'] = (
                matched_df['variant_id']
                .map(gene_dict)
                .astype('category')
            )

            # ====================================================
            # ADD AF + log10 AF
            # ====================================================

            matched_df['gnomAD_AF'] = (
                matched_df['variant_id']
                .map(af_map)
                .fillna(0)
                .astype('float32')
            )

            matched_df['log10_gnomAD_AF'] = np.log10(
                matched_df['gnomAD_AF'] + 1e-8
            ).astype('float32')

            # ====================================================
            # AF FEATURES
            # ====================================================

            matched_df['AF_missing'] = (
                matched_df['gnomAD_AF'] == 0
            ).astype('int8')

            matched_df['is_ultra_rare'] = (
                matched_df['gnomAD_AF'] < 0.0001
            ).astype('int8')

            # ====================================================
            # ADD VEP ANNOTATIONS
            # ====================================================

            def _get_vep_value(vid, col, default):
                d = vep_map.get(vid)
                if isinstance(d, dict):
                    return d.get(col, default)
                return default

            vid_series = matched_df['variant_id']

            # Numeric scores
            matched_df['vep_SIFT_score'] = vid_series.apply(
                lambda v: _get_vep_value(v, 'vep_SIFT_score', np.nan)
            ).astype('float32')

            matched_df['vep_PolyPhen_score'] = vid_series.apply(
                lambda v: _get_vep_value(v, 'vep_PolyPhen_score', np.nan)
            ).astype('float32')

            matched_df['vep_Protein_position'] = vid_series.apply(
                lambda v: _get_vep_value(v, 'vep_Protein_position', np.nan)
            ).astype('float32')

            matched_df['vep_IMPACT_score'] = vid_series.apply(
                lambda v: _get_vep_value(v, 'vep_IMPACT_score', 0)
            ).astype('int8')

            # Binary flags
            for flag_col in [
                'vep_is_missense', 'vep_is_synonymous', 'vep_is_stop_gained',
                'vep_is_frameshift', 'vep_is_splice', 'vep_is_LoF_HC',
                'vep_has_SIFT', 'vep_has_PolyPhen'
            ]:
                matched_df[flag_col] = vid_series.apply(
                    lambda v: _get_vep_value(v, flag_col, 0)
                ).astype('int8')

            # Ordinal-encoded categorical columns (int16, -1 = missing)
            for col in ORDINAL_COLS:
                raw_vals = vid_series.apply(
                    lambda v: _get_vep_value(v, col, None)
                )
                matched_df[col] = raw_vals.map(
                    ordinal_maps[col]
                ).fillna(-1).astype('int16')

            # One-hot encoded categorical columns (int8, 0/1)
            for col in ONEHOT_COLS:
                raw_vals = vid_series.apply(
                    lambda v: _get_vep_value(v, col, None)
                )
                dummies = {
                    f"{col}_{cat_val}": (raw_vals == cat_val).astype('int8')
                    for cat_val in global_cat_values[col]
                }
                if dummies:
                    matched_df = pd.concat(
                        [matched_df, pd.DataFrame(dummies)],
                        axis=1
                    )

            # ====================================================
            # OPTIONAL GENE COLUMN
            # ====================================================

            matched_df['gene'] = gene_name
            matched_df['gene'] = matched_df['gene'].astype('category')

            # ====================================================
            # SAVE
            # ====================================================

            out_path = f"{output_dir}/part_{file_id}.parquet"

            matched_df.to_parquet(
                out_path,
                compression="snappy"
            )

            total_matched += len(matched_df)

            file_id += 1

            del df
            del matched_df

            gc.collect()

            if idx % 100 == 0:

                print(
                    f"{split_name}: "
                    f"{idx}/{len(folders)} | "
                    f"Matched: {total_matched:,}"
                )

        except Exception as e:

            print(f"❌ Error in {folder}")
            print(e)

    print(f"\n📊 {split_name.upper()} SUMMARY")
    print(f"Saved parts: {file_id}")
    print(f"Matched variants: {total_matched:,}")


# ============================================================
# RUN ALL SPLITS
# ============================================================

process_split(
    "feature_extraction_train_dataset/output_train_run_*/*_VEP",
    "final_dataset_parts_train",
    "train"
)

process_split(
    "feature_extraction_cal_dataset/output_cal_run_*/*_VEP",
    "final_dataset_parts_cal",
    "cal"
)

process_split(
    "feature_extraction_test_dataset/output_test_run_*/*_VEP",
    "final_dataset_parts_test",
    "test"
)

print("\n🎯 ALL DATASETS PROCESSED SUCCESSFULLY")
