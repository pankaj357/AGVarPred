import os
import pandas as pd
from pathlib import Path

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

# -------------------------------------------------------
# API KEYS (same as main pipeline)
# -------------------------------------------------------
API_KEYS = [
    "AIzaSyDomOPmToBNr-y06A6bY2Pmqh6UJ5_7VLA",
    "AIzaSyDHHygGFg5LyMU00iPaKblvRz-bD2vMvmo",
    "AIzaSyBQbzMnScSPrTb1MI7920bj-xgKPlBChEk",
    "AIzaSyCbl9UbnFcqS7nQlL-UVZdTmbAkdUSbzgU",
    "AIzaSyCbZlR0wdaxILkPrWPjJltPsN1drrKJTCw",
    "AIzaSyD6n7hWocSn-L43DNF0jb2UMjbVOkNQgkg",
    "AIzaSyBm5oPZjvvhSW_4GzbTJ3AEry0TGBClN7U",
    "AIzaSyAr0iOqyy4Qvp_Lzq2KBVzjSJWQEla0QnE",
    "AIzaSyBy-bqASkFAa6N0-sFgRSE9wqBfltaDQi0",
    "AIzaSyCpr_wTZHd_v0oYLr70C_fhWCgLuBfQ8Ec",
    "AIzaSyCSqFoSs5W6k61eN8LGOG7jfBJrExRmppw",
    "AIzaSyAPaS9SE7WO4uw1PIG8GqX-s7BsAPflzkc",
    "AIzaSyCW5SltHTqkjSFxjnH3uuGO6jxUN9I9b8M",
    "AIzaSyBD_g6eVWLR9RcZa97l49n2stN3ZeHL664",
    "AIzaSyA-zr5rVAOPDxQduxFsDRFENIz8Nq_95F4",
    "AIzaSyC2JQLmQznLovT6Xc7snyJH0IuMuzoctXY",
    "AIzaSyBOVUkFsuq0t7ACiYxbb2dXSjg-5OrztgE",
    "AIzaSyC_sf5FZRxkPK20DfcauZM7qPyA39sL-Lw",
    "AIzaSyDU6s91wEB0C5dY4GJNt4IJk2lfRuHApUs",
    "AIzaSyCWSEt1652KiBkHJiLmSsCCQuC3f1_eBdc"
]

# -------------------------------------------------------
# PROCESS FUNCTION
# -------------------------------------------------------
def process_variant_file(input_file_path, dna_model, output_dir):
    gene_name = Path(input_file_path).stem
    print(f"\n==== Processing: {gene_name} ====")

    try:
        variants = pd.read_csv(
            input_file_path,
            sep="\t",
            header=None,
            names=["chrom", "pos", "ref", "alt", "gene"]
        )

        if variants.empty:
            print("⚠ Empty file")
            return {"gene": gene_name, "status": "empty"}

        print(f"Total variants in input: {len(variants)}")

        gene_outdir = os.path.join(output_dir, f"{gene_name}_VEP")
        os.makedirs(gene_outdir, exist_ok=True)

        outfile = os.path.join(
            gene_outdir, f"{gene_name}_ALL_VEP_RAW_SCORE_MATRIX.parquet"
        )

        if os.path.exists(outfile):
            print(f"⏩ SKIPPING {gene_name}: already exists")
            return {"gene": gene_name, "status": "skipped"}

        all_scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
        SEQ_LEN = dna_client.SEQUENCE_LENGTH_1MB

        all_rows = []
        all_columns = set()

        for idx, row in variants.iterrows():
            if idx % 100 == 0:
                print(f"  {idx}/{len(variants)}")

            try:
                chrom = str(row["chrom"])
                if chrom == "chrMT":
                    chrom = "chrM"
                pos = int(row["pos"])
                ref = str(row["ref"])
                alt = str(row["alt"])

                variant_id = f"{chrom}_{pos}_{ref}_{alt}"

                variant_obj = genome.Variant(
                    chromosome=chrom,
                    position=pos,
                    reference_bases=ref,
                    alternate_bases=alt,
                    name=variant_id,
                )

                interval = variant_obj.reference_interval.resize(SEQ_LEN)

                scores = dna_model.score_variant(
                    interval=interval,
                    variant=variant_obj,
                    variant_scorers=all_scorers,
                    organism=dna_client.Organism.HOMO_SAPIENS,
                )

                df_scores = variant_scorers.tidy_scores(scores)

                vector = {"variant_id": variant_id}

                for _, s in df_scores.iterrows():
                    col = f"{s.get('output_type','unk')}__{s.get('biosample_name','unk')}__{s.get('track_name','unk')}"
                    vector[col] = s.get("raw_score", None)
                    all_columns.add(col)

                all_rows.append(vector)

            except Exception as e:
                print(f"⚠ Variant error: {str(e)[:80]}")
                continue

        if all_rows:
            df = pd.DataFrame(all_rows).set_index("variant_id")
            df = df.reindex(columns=sorted(all_columns))

            df.to_parquet(outfile)
            df.to_csv(outfile.replace(".parquet", ".csv"))

            print(f"✅ DONE: {gene_name} (Total variants: {len(df)})")
            return {"gene": gene_name, "status": "success", "variants": len(df)}

        return {"gene": gene_name, "status": "failed"}

    except Exception as e:
        print(f"❌ ERROR in {gene_name}: {e}")
        return {"gene": gene_name, "status": "error"}


# -------------------------------------------------------
# MAIN
# -------------------------------------------------------
def main():
    RUN_ID = int(os.environ.get("RUN_ID", 0))

    if RUN_ID >= len(API_KEYS):
        raise ValueError("RUN_ID exceeds API_KEYS")

    API_KEY = API_KEYS[RUN_ID]
    print(f"🚀 RUN_ID: {RUN_ID}")

    dna_model = dna_client.create(API_KEY)

    input_path = f"external_validation/processing/chunks/external_humsavar/run_{RUN_ID+1}"
    files = list(Path(input_path).glob("*.txt"))

    print(f"📂 Files: {len(files)}")

    output_dir = f"external_validation/processing/features/humsavar/output_run_{RUN_ID+1}"
    os.makedirs(output_dir, exist_ok=True)

    summary_dir = "external_validation/summaries"
    os.makedirs(summary_dir, exist_ok=True)

    results = []
    for f in files:
        res = process_variant_file(str(f), dna_model, output_dir)
        results.append(res)

    pd.DataFrame(results).to_csv(
        f"{summary_dir}/summary_humsavar_run_{RUN_ID+1}.csv", index=False
    )

    print("🎉 RUN COMPLETE")


if __name__ == "__main__":
    main()
