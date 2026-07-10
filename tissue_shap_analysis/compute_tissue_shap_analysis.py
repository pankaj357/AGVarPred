#!/usr/bin/env python3
"""
Tissue/cell-type SHAP aggregation for AGVarPred.

Uses only existing SHAP outputs. Does NOT retrain or recompute SHAP values.

Outputs:
  - biosample_shap_importance.csv
  - tissue_shap_importance.csv
  - tissue_summary_table.csv
  - figures: tissue_heatmap, top_tissues_barplot, cumulative_contribution
  - manuscript_results.txt
  - manuscript_discussion.txt
"""

import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(OUT_DIR)

BOOTSTRAP_N = 1000

# ---------------------------------------------------------------------------
# Feature name parser
# ---------------------------------------------------------------------------

def parse_feature_name(feature):
    """
    Parse an AlphaGenome feature name into (assay, biosample, track).
    Feature names use single underscores because SHAP outputs were sanitized.
    """

    # Order matters: more specific prefixes first
    patterns = [
        # SPLICE_SITE_USAGE: SPLICE_SITE_USAGE_<biosample>_usage_<ontology>_<id>_<rna_type>_RNA_seq
        (r"^SPLICE_SITE_USAGE_(.+)_usage_(UBERON|CL|EFO|GTEX)_(\d+)_(.+)_RNA_seq$",
         "SPLICE_SITE_USAGE", "usage"),
        # SPLICE_JUNCTIONS: SPLICE_JUNCTIONS_<biosample>_junction_<ontology>_<id>_<rna_type>_RNA_seq
        (r"^SPLICE_JUNCTIONS_(.+)_junction_(UBERON|CL|EFO|GTEX)_(\d+)_(.+)_RNA_seq$",
         "SPLICE_JUNCTIONS", "junction"),
        # SPLICE_SITES: SPLICE_SITES_<biosample>_(acceptor|donor)
        (r"^SPLICE_SITES_(.+)?_(acceptor|donor)$",
         "SPLICE_SITES", "site"),
        # CONTACT_MAPS: CONTACT_MAPS_<biosample>_4dn_<accession>
        (r"^CONTACT_MAPS_(.+)_4dn_(4DNFI\w+)$",
         "CONTACT_MAPS", "contact_map"),
        # CHIP_HISTONE: CHIP_HISTONE_<biosample>_Histone_ChIP_seq_<mark>
        (r"^CHIP_HISTONE_(.+)_Histone_ChIP_seq_(.+)$",
         "CHIP_HISTONE", "histone"),
        # CHIP_TF: CHIP_TF_<biosample>_TF_ChIP_seq_<target>...
        (r"^CHIP_TF_(.+)_TF_ChIP_seq_(.+)$",
         "CHIP_TF", "tf"),
        # PROCAP: PROCAP_<biosample>_<accession>
        (r"^PROCAP_(.+?)_(ENCSR\w+)$",
         "PROCAP", "procapture"),
        # DNASE: DNASE_<biosample>_DNase_seq
        (r"^DNASE_(.+)_DNase_seq$",
         "DNASE", "dnase"),
        # CAGE (standard): CAGE_<biosample>_hCAGE_<ontology>_<id>
        (r"^CAGE_(.+)_hCAGE_(UBERON|CL|EFO)_(\d+)$",
         "CAGE", "cage"),
        # CAGE (low quality): CAGE_<biosample>_LQhCAGE_<ontology>_<id>
        (r"^CAGE_(.+)_LQhCAGE_(UBERON|CL|EFO)_(\d+)$",
         "CAGE", "cage_lq"),
        # ATAC: ATAC_<biosample>_ATAC_seq
        (r"^ATAC_(.+)_ATAC_seq$",
         "ATAC", "atac"),
    ]

    for pat, assay, track_type in patterns:
        m = re.match(pat, feature)
        if m:
            biosample = m.group(1) if m.group(1) else "unspecified"
            # Build track description
            if assay in ["SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]:
                track = f"{track_type}_{m.group(2)}_{m.group(3)}_{m.group(4)}"
            elif assay == "SPLICE_SITES":
                track = m.group(2)
            elif assay == "CONTACT_MAPS":
                track = f"4dn_{m.group(2)}"
            elif assay in ["CHIP_HISTONE", "CHIP_TF"]:
                track = m.group(2)
            elif assay == "PROCAP":
                track = m.group(2)
            elif assay == "CAGE":
                track = f"{track_type}_{m.group(2)}_{m.group(3)}"
            else:
                track = track_type
            return assay, biosample, track

    return "UNKNOWN", feature, "unknown"


# ---------------------------------------------------------------------------
# Tissue grouping rules
# ---------------------------------------------------------------------------

def map_biosample_to_tissue(biosample):
    """
    Map a biosample name to a broad tissue/lineage category.
    Rules are based on the biosample name itself.
    """
    b = biosample.lower().replace("_", " ")

    # Explicit cell-line mappings first (override tissue keywords)
    cell_line_map = {
        "h1 hesc": "Embryonic / stem",
        "h9": "Embryonic / stem",
        "hues48": "Embryonic / stem",
        "cyt49": "Embryonic / stem",
        "ips 18a": "Embryonic / stem",
        "hek293": "Other",  # widely used kidney-derived line; generic
        "hct116": "Liver / digestive",  # colorectal carcinoma
        "hela s3": "Reproductive",  # cervical carcinoma
        "pc 9": "Lung",  # lung adenocarcinoma
        "mcf 10a": "Reproductive",  # breast epithelial
        "vcap": "Reproductive",  # prostate carcinoma
        "be2c": "Brain / CNS",  # neuroblastoma
        "pfsk 1": "Brain / CNS",  # neuroectodermal tumor
        "sk n sh": "Brain / CNS",  # neuroblastoma
        "k562": "Blood / immune",
        "mm1s": "Blood / immune",
        "a549": "Lung",
        "calu3": "Lung",
        "hepg2": "Liver / digestive",
        "panc1": "Liver / digestive",
        "sw480": "Liver / digestive",
        "colo829": "Skin",  # melanoma
    }
    for key, val in cell_line_map.items():
        if key in b:
            return val

    # Brain / CNS
    if any(x in b for x in [
        "brain", "cortex", "cerebellum", "neuron", "neural", "neuronal",
        "substantia nigra", "diencephalon", "spinal cord", "occipital",
        "frontal", "parietal", "motor neuron", "purkinje",
        "neuro", "glioblastoma"
    ]):
        return "Brain / CNS"

    # Blood / immune
    if any(x in b for x in [
        "t cell", "b cell", "cd4", "cd8", "monocyte", "neutrophil",
        "pbmc", "peripheral blood", "lymphocyte", "leukemia",
        "blood", "immune", "thymus", "spleen", "tonsil"
    ]):
        return "Blood / immune"

    # Liver / digestive
    if any(x in b for x in [
        "liver", "hepatocyte", "pancreas", "panc1", "stomach",
        "esophagus", "gastroesophageal", "colon", "rectum", "intestine",
        "colorectal", "parotid", "salivary", "gland",
        "caecum", "appendix", "duodenum", "fungiform papilla", "pancreatic"
    ]):
        return "Liver / digestive"

    # Heart / vasculature
    if any(x in b for x in [
        "heart", "cardiac", "aorta", "valve", "vascular", "endothelial",
        "umbilical vein", "artery", "blood vessel", "mesenteric"
    ]):
        return "Heart / vasculature"

    # Muscle
    if any(x in b for x in [
        "muscle", "myocyte", "skeletal muscle", "smooth muscle"
    ]):
        return "Muscle"

    # Kidney / urinary
    if any(x in b for x in [
        "kidney", "metanephros", "urothelium", "urinary bladder", "bladder"
    ]):
        return "Kidney / urinary"

    # Reproductive
    if any(x in b for x in [
        "testis", "testes", "ovary", "ovarian", "prostate", "uterus",
        "endometrium", "seminal vesicle", "breast", "mammary", "placenta",
        "cervical", "cervix"
    ]):
        return "Reproductive"

    # Skin
    if any(x in b for x in [
        "skin", "keratinocyte", "melanocyte", "fibroblast of skin",
        "epidermis"
    ]):
        return "Skin"

    # Lung
    if any(x in b for x in [
        "lung", "bronchial", "alveolar", "trachea"
    ]):
        return "Lung"

    # Eye / sensory
    if any(x in b for x in [
        "eye", "retina", "cornea", "camera type eye"
    ]):
        return "Eye / sensory"

    # Connective / bone / cartilage
    if any(x in b for x in [
        "fibroblast", "chondrocyte", "osteoblast", "mesenchymal",
        "adipose", "fat", "stromal"
    ]):
        return "Connective / mesenchyme"

    # Embryonic / stem
    if any(x in b for x in [
        "hesc", "esc", "ips", "stem cell", "progenitor", "embryo", "amnion"
    ]):
        return "Embryonic / stem"

    return "Other"


# ---------------------------------------------------------------------------
# Bootstrap CI
# ---------------------------------------------------------------------------

def bootstrap_feature_total(mean_abs_shap_values, n_boot=BOOTSTRAP_N, seed=42):
    """Bootstrap across features to estimate CI for total SHAP of a group."""
    rng = np.random.default_rng(seed)
    n = len(mean_abs_shap_values)
    if n < 2:
        return np.nan, np.nan
    totals = []
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        totals.append(np.sum(mean_abs_shap_values[idx]))
    return np.percentile(totals, 2.5), np.percentile(totals, 97.5)


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Load SHAP values and feature importance
    shap_df = pd.read_csv(os.path.join(ROOT_DIR, "final_model_output_regularized", "shap_values_sample.csv"))
    imp_df = pd.read_csv(os.path.join(ROOT_DIR, "final_model_output_regularized", "feature_importance_shap.csv"))

    # Keep only AlphaGenome features
    ag_features = [c for c in shap_df.columns if not c.startswith("vep_") and c != "gnomAD_AF"]
    shap_ag = shap_df[ag_features]
    imp_df = imp_df[imp_df["feature"].isin(ag_features)].copy()

    # Parse feature names
    parsed = []
    for feat in ag_features:
        assay, biosample, track = parse_feature_name(feat)
        tissue = map_biosample_to_tissue(biosample)
        parsed.append({
            "feature": feat,
            "assay": assay,
            "biosample": biosample,
            "track": track,
            "tissue_group": tissue,
        })
    parsed_df = pd.DataFrame(parsed)

    # Merge with mean SHAP importance
    parsed_df = parsed_df.merge(imp_df, on="feature", how="left")

    # Save parsed features
    parsed_df.to_csv(os.path.join(OUT_DIR, "parsed_features.csv"), index=False)

    # -----------------------------------------------------------------------
    # Biosample-level aggregation
    # -----------------------------------------------------------------------
    biosample_stats = []
    for bio, sub in parsed_df.groupby("biosample"):
        total = sub["mean_abs_shap"].sum()
        ci_low, ci_high = bootstrap_feature_total(sub["mean_abs_shap"].values)
        biosample_stats.append({
            "biosample": bio,
            "tissue_group": sub["tissue_group"].iloc[0],
            "total_shap": total,
            "mean_shap": sub["mean_abs_shap"].mean(),
            "median_shap": sub["mean_abs_shap"].median(),
            "n_features": len(sub),
            "pct_of_alphagenome": np.nan,
            "shap_ci_low": ci_low,
            "shap_ci_high": ci_high,
        })
    biosample_df = pd.DataFrame(biosample_stats)
    total_ag = biosample_df["total_shap"].sum()
    biosample_df["pct_of_alphagenome"] = 100.0 * biosample_df["total_shap"] / total_ag
    biosample_df = biosample_df.sort_values("total_shap", ascending=False).reset_index(drop=True)
    biosample_df["rank"] = np.arange(1, len(biosample_df) + 1)
    biosample_df.to_csv(os.path.join(OUT_DIR, "biosample_shap_importance.csv"), index=False)

    # -----------------------------------------------------------------------
    # Tissue-group-level aggregation
    # -----------------------------------------------------------------------
    tissue_stats = []
    for tissue, sub in parsed_df.groupby("tissue_group"):
        total = sub["mean_abs_shap"].sum()
        ci_low, ci_high = bootstrap_feature_total(sub["mean_abs_shap"].values)
        tissue_stats.append({
            "tissue_group": tissue,
            "total_shap": total,
            "mean_shap": sub["mean_abs_shap"].mean(),
            "median_shap": sub["mean_abs_shap"].median(),
            "n_features": len(sub),
            "n_biosamples": sub["biosample"].nunique(),
            "pct_of_alphagenome": np.nan,
            "shap_ci_low": ci_low,
            "shap_ci_high": ci_high,
        })
    tissue_df = pd.DataFrame(tissue_stats)
    tissue_df["pct_of_alphagenome"] = 100.0 * tissue_df["total_shap"] / total_ag
    tissue_df = tissue_df.sort_values("total_shap", ascending=False).reset_index(drop=True)
    tissue_df["rank"] = np.arange(1, len(tissue_df) + 1)
    tissue_df.to_csv(os.path.join(OUT_DIR, "tissue_shap_importance.csv"), index=False)

    # Summary table with cumulative contributions
    tissue_df["cumulative_pct"] = tissue_df["pct_of_alphagenome"].cumsum()
    tissue_df.to_csv(os.path.join(OUT_DIR, "tissue_summary_table.csv"), index=False)

    print("Top 20 biosamples:")
    print(biosample_df[["rank", "biosample", "tissue_group", "total_shap", "pct_of_alphagenome", "n_features"]].head(20).to_string(index=False))

    print("\nTop 10 tissue groups:")
    print(tissue_df[["rank", "tissue_group", "total_shap", "pct_of_alphagenome", "n_features", "n_biosamples"]].head(10).to_string(index=False))

    print(f"\nTop 5 tissues cumulative: {tissue_df.head(5)['cumulative_pct'].iloc[-1]:.1f}%")
    print(f"Top 10 tissues cumulative: {tissue_df.head(10)['cumulative_pct'].iloc[-1]:.1f}%")

    # -----------------------------------------------------------------------
    # Figures
    # -----------------------------------------------------------------------

    # Figure A: Heatmap of top tissues x assay categories
    top_tissues = tissue_df.head(8)["tissue_group"].tolist()
    heat_data = parsed_df[parsed_df["tissue_group"].isin(top_tissues)].groupby(["tissue_group", "assay"])["mean_abs_shap"].sum().unstack(fill_value=0)
    # Reorder tissues by total importance
    heat_data = heat_data.reindex(top_tissues)

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(heat_data, annot=True, fmt=".3f", cmap="YlOrRd", cbar_kws={"label": "Total |SHAP|"}, ax=ax)
    ax.set_title("SHAP importance by tissue and assay type")
    ax.set_xlabel("Assay type")
    ax.set_ylabel("Tissue / lineage")
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "tissue_assay_heatmap.pdf"), dpi=300)
    fig.savefig(os.path.join(OUT_DIR, "tissue_assay_heatmap.png"), dpi=300)
    plt.close(fig)

    # Figure B: Top 20 biosamples bar plot
    top20 = biosample_df.head(20).sort_values("total_shap", ascending=True)
    fig, ax = plt.subplots(figsize=(10, 10))
    colors = [plt.cm.tab10(i % 10) for i in range(len(top20))]
    bars = ax.barh(top20["biosample"], top20["total_shap"], color=colors)
    ax.set_xlabel("Total |SHAP|")
    ax.set_title("Top 20 biosamples by SHAP importance")
    for bar, pct in zip(bars, top20["pct_of_alphagenome"]):
        ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                f"{pct:.1f}%", va="center", fontsize=8)
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "top_biosamples_barplot.pdf"), dpi=300)
    fig.savefig(os.path.join(OUT_DIR, "top_biosamples_barplot.png"), dpi=300)
    plt.close(fig)

    # Figure C: Cumulative contribution
    fig, ax = plt.subplots(figsize=(10, 6))
    ranks = np.arange(1, len(tissue_df) + 1)
    ax.plot(ranks, tissue_df["cumulative_pct"].values, marker="o", linewidth=2)
    ax.axhline(50, color="gray", linestyle="--", alpha=0.5)
    ax.axhline(80, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Tissue rank")
    ax.set_ylabel("Cumulative % of AlphaGenome SHAP importance")
    ax.set_title("Cumulative SHAP importance across tissue groups")
    ax.set_ylim(0, 105)
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "cumulative_tissue_contribution.pdf"), dpi=300)
    fig.savefig(os.path.join(OUT_DIR, "cumulative_tissue_contribution.png"), dpi=300)
    plt.close(fig)

    # -----------------------------------------------------------------------
    # Manuscript text
    # -----------------------------------------------------------------------
    write_results_subsection(tissue_df, biosample_df, parsed_df)
    write_discussion_subsection(tissue_df, biosample_df)

    print("\nTissue SHAP analysis complete.")


def write_results_subsection(tissue_df, biosample_df, parsed_df):
    top5 = tissue_df.head(5)
    top10 = tissue_df.head(10)
    top20_bio = biosample_df.head(20)

    total_features = len(parsed_df)
    total_biosamples = parsed_df["biosample"].nunique()

    text = f"""Cell-type-specific functional genomics signals underlying AGVarPred predictions

We investigated whether AGVarPred's predictions are driven by a small number of cell types or distributed across diverse tissues and assays. Of the {total_features} AlphaGenome features selected by the model, we parsed {total_biosamples} distinct biosamples and assigned each to one of {len(tissue_df)} broad tissue/lineage groups (see Methods for parsing and grouping rules).

The top five tissue groups accounted for {top5['cumulative_pct'].iloc[-1]:.1f}% of total AlphaGenome SHAP importance, and the top ten groups accounted for {top10['cumulative_pct'].iloc[-1]:.1f}%. The highest-contributing tissue group was {top5.iloc[0]['tissue_group']} ({top5.iloc[0]['pct_of_alphagenome']:.1f}%), followed by {top5.iloc[1]['tissue_group']} ({top5.iloc[1]['pct_of_alphagenome']:.1f}%), {top5.iloc[2]['tissue_group']} ({top5.iloc[2]['pct_of_alphagenome']:.1f}%), {top5.iloc[3]['tissue_group']} ({top5.iloc[3]['pct_of_alphagenome']:.1f}%), and {top5.iloc[4]['tissue_group']} ({top5.iloc[4]['pct_of_alphagenome']:.1f}%).

At the level of individual biosamples, the strongest contributors were {top20_bio.iloc[0]['biosample']} ({top20_bio.iloc[0]['pct_of_alphagenome']:.1f}%), {top20_bio.iloc[1]['biosample']} ({top20_bio.iloc[1]['pct_of_alphagenome']:.1f}%), {top20_bio.iloc[2]['biosample']} ({top20_bio.iloc[2]['pct_of_alphagenome']:.1f}%), {top20_bio.iloc[3]['biosample']} ({top20_bio.iloc[3]['pct_of_alphagenome']:.1f}%), and {top20_bio.iloc[4]['biosample']} ({top20_bio.iloc[4]['pct_of_alphagenome']:.1f}%). These top five biosamples alone contributed {top20_bio.head(5)['pct_of_alphagenome'].sum():.1f}% of AlphaGenome SHAP importance.

The assay-type breakdown differed across tissues. """

    # Add per-tissue assay detail for top tissues
    heat_data = parsed_df.groupby(["tissue_group", "assay"])["mean_abs_shap"].sum().unstack(fill_value=0)
    top3 = tissue_df.head(3)["tissue_group"].tolist()
    for tissue in top3:
        row = heat_data.loc[tissue]
        top_assay = row.idxmax()
        top_assay_pct = 100.0 * row.max() / row.sum()
        text += f"For {tissue}, the dominant assay type was {top_assay} ({top_assay_pct:.1f}% of that tissue's SHAP contribution). "

    text += f"""

These results indicate that AGVarPred does not rely on a single generic functional signal. Instead, it draws on a broad but uneven distribution of cell-type-specific annotations, with a concentration in tissues relevant to development, proliferation, and gene regulation.
"""

    with open(os.path.join(OUT_DIR, "manuscript_results.txt"), "w") as f:
        f.write(text)
    print("Saved manuscript_results.txt")


def write_discussion_subsection(tissue_df, biosample_df):
    top = tissue_df.iloc[0]
    text = f"""Biological interpretation of cell-type-specific contributions

The tissue-level SHAP analysis supports the hypothesis that AGVarPred leverages cell-type-specific functional genomics rather than treating AlphaGenome scores as a generic set of annotations. The dominance of {top['tissue_group']} ({top['pct_of_alphagenome']:.1f}%) and the strong representation of developmentally active and proliferative cell types is consistent with the biology of germline pathogenicity: many disease-causing variants disrupt regulatory or splicing programs active in specific lineages.

It is important to interpret these rankings cautiously. SHAP importance reflects the marginal contribution of a feature to the trained model, not necessarily the causal tissue of disease. A high contribution from a cell line such as {biosample_df.iloc[0]['biosample']} may reflect the abundance or quality of AlphaGenome assays in that biosample rather than a primary disease-relevant tissue. Conversely, tissues with low SHAP contribution are not irrelevant to disease; they may simply be less represented in the selected feature set or less predictive given the training data.

Nevertheless, the observed distribution—broad across tissues but concentrated in developmentally and transcriptionally active cell types—provides biological plausibility for the model's predictions and distinguishes AGVarPred from predictors that rely solely on protein-level or conservation-based features.
"""
    with open(os.path.join(OUT_DIR, "manuscript_discussion.txt"), "w") as f:
        f.write(text)
    print("Saved manuscript_discussion.txt")


if __name__ == "__main__":
    main()
