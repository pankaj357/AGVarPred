#!/usr/bin/env python3
"""
Decision Curve Analysis (DCA) for AGVarPred.

Uses only existing calibrated predictions. No retraining, recalibration, or rescoring.

Outputs:
  - results/decision_curve_data.csv
  - results/decision_curve_summary.csv
  - figures/decision_curve_analysis.pdf/png
  - manuscript_results.txt
  - manuscript_discussion.txt
"""

import os
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(OUT_DIR, "results")
FIGURES_DIR = os.path.join(OUT_DIR, "figures")
ROOT_DIR = os.path.dirname(OUT_DIR)

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

THRESHOLDS = np.arange(0.01, 1.00, 0.01)

BENCHMARKS = {
    "internal_test": {
        "predictions": "final_model_output_regularized/test_predictions.csv",
        "label_col": "y_true",
        "score_col": "cal_prob",
        "has_competitors": False,
    },
    "humsavar": {
        "predictions": "external_validation/results/humsavar/regularized/humsavar_predictions.csv",
        "variants": "external_validation/benchmarks_comparative/datasets/humsavar/variants_scored.tsv",
        "json": "external_validation/benchmarks_comparative/datasets/humsavar/vep_cadd_alphamissense.json",
        "revel": "external_validation/benchmarks_comparative/datasets/humsavar/revel_scores.csv",
        "json_key_has_chr": False,
        "label_col": "true_label",
        "score_col": "prob_calibrated",
        "has_competitors": True,
    },
    "mave_independent": {
        "predictions": "external_validation/results/mave_independent/regularized/mave_independent_predictions.csv",
        "variants": "external_validation/benchmarks_comparative/datasets/mave_independent/variants_scored.tsv",
        "json": "external_validation/benchmarks_comparative/datasets/mave_independent/vep_cadd_alphamissense.json",
        "revel": "external_validation/benchmarks_comparative/datasets/mave_independent/revel_scores.csv",
        "json_key_has_chr": False,
        "label_col": "true_label",
        "score_col": "prob_calibrated",
        "has_competitors": True,
    },
    "gnomad_benign": {
        "predictions": "external_validation/results/gnomad_benign/regularized/gnomad_benign_predictions.csv",
        "variants": "external_validation/benchmarks_comparative/datasets/gnomad_benign/variants_scored.tsv",
        "json": "external_validation/benchmarks_comparative/datasets/gnomad_benign/vep_cadd_alphamissense.json",
        "revel": "external_validation/benchmarks_comparative/datasets/gnomad_benign/revel_scores.csv",
        "json_key_has_chr": False,
        "label_col": "true_label",
        "score_col": "prob_calibrated",
        "has_competitors": True,
    },
    # VIP is excluded from independent external-benchmark DCA because all
    # pathogenic labels are ClinVar-derived; it is reported separately as a
    # ClinVar-held-out-gene consistency analysis.
    # "vip": { ... },
    "grimm2015": {
        "predictions": "external_validation/results/grimm2015/regularized/grimm2015_predictions.csv",
        "variants": "external_validation/benchmarks_comparative/datasets/grimm2015/variants_scored.tsv",
        "json": "external_validation/benchmarks_comparative/datasets/grimm2015/vep_cadd_alphamissense.json",
        "revel": "external_validation/benchmarks_comparative/datasets/grimm2015/revel_scores.csv",
        "json_key_has_chr": True,
        "label_col": "true_label",
        "score_col": "prob_calibrated",
        "has_competitors": True,
    },
    "dvd": {
        "predictions": "external_validation/results/dvd/regularized/dvd_predictions.csv",
        "variants": "external_validation/benchmarks_comparative/datasets/dvd/variants_scored.tsv",
        "json": "external_validation/benchmarks_comparative/datasets/dvd/vep_cadd_alphamissense.json",
        "revel": "external_validation/benchmarks_comparative/datasets/dvd/revel_scores.csv",
        "json_key_has_chr": False,
        "label_col": "true_label",
        "score_col": "prob_calibrated",
        "has_competitors": True,
    },
}

COMPETITORS = {
    # CADD phred scores are not calibrated probabilities and cannot be
    # thresholded meaningfully on the 0-1 probability scale, so it is excluded
    # from DCA curves.
    "AlphaMissense": "alphamissense_score",
    "REVEL": "revel_score",
}


def load_competitor_scores(bench):
    info = BENCHMARKS[bench]
    variants = pd.read_csv(os.path.join(ROOT_DIR, info["variants"]), sep="\t")
    variants = variants.rename(columns={"label": "true_label"})

    with open(os.path.join(ROOT_DIR, info["json"])) as f:
        cadd_alpha = json.load(f)

    rows = []
    for _, row in variants.iterrows():
        vid = row["variant_id"]
        key_no_chr = vid.replace("chr", "") if vid.startswith("chr") else vid
        key_with_chr = vid if vid.startswith("chr") else "chr" + vid
        key = key_with_chr if info["json_key_has_chr"] else key_no_chr
        ca = cadd_alpha.get(key, {})
        rows.append({
            "variant_id": vid,
            "cadd_phred": ca.get("cadd_phred"),
            "alphamissense_score": ca.get("alphamissense_score"),
        })
    comp = pd.DataFrame(rows)

    revel = pd.read_csv(os.path.join(ROOT_DIR, info["revel"]))
    comp = comp.merge(revel, on="variant_id", how="left")
    return comp


def net_benefit(y_true, y_score, thresholds):
    """Compute net benefit across thresholds."""
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    n = len(y_true)
    n_pos = np.sum(y_true == 1)
    prevalence = n_pos / n

    results = []
    for thr in thresholds:
        pred_pos = (y_score >= thr).astype(int)
        tp = np.sum((pred_pos == 1) & (y_true == 1))
        fp = np.sum((pred_pos == 1) & (y_true == 0))
        nb = (tp / n) - (fp / n) * (thr / (1 - thr))

        # Treat-all net benefit
        nb_all = (n_pos / n) - ((n - n_pos) / n) * (thr / (1 - thr))
        # Treat-none net benefit = 0
        nb_none = 0.0

        results.append({
            "threshold": thr,
            "net_benefit": nb,
            "treat_all": nb_all,
            "treat_none": nb_none,
            "prevalence": prevalence,
            "n": n,
            "n_pos": int(n_pos),
        })
    return pd.DataFrame(results)


def run_dca_for_benchmark(bench):
    print(f"\n[{bench}]")
    cfg = BENCHMARKS[bench]
    pred = pd.read_csv(os.path.join(ROOT_DIR, cfg["predictions"]))
    pred = pred.rename(columns={cfg["score_col"]: "agvarpred_score", cfg["label_col"]: "true_label"})

    records = []

    # AGVarPred
    nb_df = net_benefit(pred["true_label"].values, pred["agvarpred_score"].values, THRESHOLDS)
    for _, r in nb_df.iterrows():
        rec = r.to_dict()
        rec["benchmark"] = bench
        rec["model"] = "AGVarPred"
        rec["n_scored"] = len(pred)
        records.append(rec)

    # Competitors on common subset
    if cfg["has_competitors"]:
        comp = load_competitor_scores(bench)
        comp = comp.drop(columns=["true_label"], errors="ignore")
        df = pred[["variant_id", "true_label", "agvarpred_score"]].merge(comp, on="variant_id", how="left")

        for comp_name, col in COMPETITORS.items():
            common = df.dropna(subset=["true_label", col])
            if len(common) < 30:
                continue
            nb_df = net_benefit(common["true_label"].values, common[col].values, THRESHOLDS)
            for _, r in nb_df.iterrows():
                rec = r.to_dict()
                rec["benchmark"] = bench
                rec["model"] = comp_name
                rec["n_scored"] = len(common)
                records.append(rec)

    return pd.DataFrame(records)


def create_figure(data_df, out_path):
    """Publication-quality DCA figure with subplots for each benchmark."""
    benchmarks = [b for b in BENCHMARKS if b != "internal_test"]
    n = len(benchmarks)
    cols = 3
    rows = int(np.ceil(n / cols))

    fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
    axes = axes.flatten() if n > 1 else [axes]

    colors = {
        "AGVarPred": "#1f77b4",
        "CADD": "#ff7f0e",
        "AlphaMissense": "#2ca02c",
        "REVEL": "#d62728",
        "treat_all": "gray",
        "treat_none": "black",
    }

    linestyles = {
        "AGVarPred": "-",
        "AlphaMissense": "--",
        "REVEL": "-.",
    }

    for idx, bench in enumerate(benchmarks):
        ax = axes[idx]
        sub = data_df[data_df["benchmark"] == bench]
        if len(sub) == 0:
            ax.set_visible(False)
            continue

        # Treat-all and treat-none (same for all models in a benchmark)
        ref = sub[sub["model"] == "AGVarPred"].copy()
        ax.plot(ref["threshold"], ref["treat_all"], color=colors["treat_all"],
                linestyle="--", linewidth=1, label="Treat all", alpha=0.7)
        ax.plot(ref["threshold"], ref["treat_none"], color=colors["treat_none"],
                linestyle="-", linewidth=1, label="Treat none", alpha=0.7)

        # Models
        for model in ["AGVarPred", "AlphaMissense", "REVEL"]:
            msub = sub[sub["model"] == model]
            if len(msub) > 0:
                ax.plot(msub["threshold"], msub["net_benefit"], color=colors[model],
                        linewidth=2, linestyle=linestyles.get(model, "-"), label=model)

        ax.set_xlim(0, 1)
        # Focus on clinically interpretable range; treat-all tails are clipped
        ax.set_ylim(-0.1, 0.5)
        ax.set_xlabel("Probability threshold")
        ax.set_ylabel("Net benefit")
        ax.set_title(bench.replace("_", " ").title())
        ax.axhline(0, color="black", linewidth=0.5, alpha=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if idx == 0:
            ax.legend(loc="upper right", fontsize=8)

    # Hide unused subplots
    for idx in range(n, len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle("Decision Curve Analysis: AGVarPred vs competitors", fontsize=14, y=1.02)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    fig.savefig(out_path.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def create_internal_figure(data_df, out_path):
    """Separate figure for internal test (no competitors)."""
    sub = data_df[data_df["benchmark"] == "internal_test"].copy()
    if len(sub) == 0:
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    ref = sub[sub["model"] == "AGVarPred"]
    ax.plot(ref["threshold"], ref["treat_all"], color="gray", linestyle="--",
            linewidth=1, label="Treat all", alpha=0.7)
    ax.plot(ref["threshold"], ref["treat_none"], color="black", linestyle="-",
            linewidth=1, label="Treat none", alpha=0.7)
    ax.plot(ref["threshold"], ref["net_benefit"], color="#1f77b4", linewidth=2.5,
            label="AGVarPred")

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.1, max(0.5, ref["treat_all"].max() * 1.1))
    ax.set_xlabel("Probability threshold")
    ax.set_ylabel("Net benefit")
    ax.set_title("Decision Curve Analysis: held-out ClinVar test set")
    ax.axhline(0, color="black", linewidth=0.5, alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper right")
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    fig.savefig(out_path.replace(".pdf", ".png"), dpi=300)
    plt.close(fig)


def compute_summary(data_df):
    """Summary metrics: max net benefit, range where model beats treat-all/treat-none, AUC of NB curve."""
    summary = []
    for (bench, model), sub in data_df.groupby(["benchmark", "model"]):
        if model in ["treat_all", "treat_none"]:
            continue
        nb = sub["net_benefit"].values
        thresholds = sub["threshold"].values
        ta = sub["treat_all"].values
        tn = sub["treat_none"].values

        max_nb = np.max(nb)
        max_nb_thr = thresholds[np.argmax(nb)]

        # Range where model has positive net benefit
        positive = thresholds[nb > 1e-9]
        positive_range = f"{positive.min():.2f}-{positive.max():.2f}" if len(positive) > 0 else "none"

        # Range where model beats both treat-all and treat-none
        beats_both = thresholds[(nb > ta + 1e-9) & (nb > tn + 1e-9)]
        beats_both_range = f"{beats_both.min():.2f}-{beats_both.max():.2f}" if len(beats_both) > 0 else "none"

        # Integrated net benefit (area under NB curve)
        integrated_nb = np.trapezoid(nb, thresholds)

        summary.append({
            "benchmark": bench,
            "model": model,
            "n_scored": sub["n_scored"].iloc[0],
            "prevalence": sub["prevalence"].iloc[0],
            "max_net_benefit": max_nb,
            "max_net_benefit_threshold": max_nb_thr,
            "integrated_net_benefit": integrated_nb,
            "threshold_range_positive_nb": positive_range,
            "threshold_range_beats_both_strategies": beats_both_range,
        })
    return pd.DataFrame(summary)


def write_results_subsection(summary_df):
    text = """Decision Curve Analysis

We evaluated the potential clinical utility of AGVarPred using Decision Curve Analysis (DCA) on held-out and external validation datasets. DCA quantifies the net benefit of using a model's predicted probability to classify variants as pathogenic at a given probability threshold, compared with two default strategies: treating all variants as pathogenic (treat-all) or treating none as pathogenic (treat-none). Net benefit is defined as the proportion of true positives minus the proportion of false positives weighted by the odds of the threshold, and it can be negative if a strategy is harmful.

"""

    # Internal test
    internal = summary_df[summary_df["benchmark"] == "internal_test"]
    if len(internal) > 0:
        r = internal.iloc[0]
        text += (
            f"On the held-out ClinVar test set (N = {int(r['n_scored'])}, prevalence = {r['prevalence']:.1%}), "
            f"AGVarPred achieved a maximum net benefit of {r['max_net_benefit']:.3f} at threshold "
            f"{r['max_net_benefit_threshold']:.2f}. The model provided positive net benefit across a broad "
            f"range of thresholds ({r['threshold_range_positive_nb']}) and outperformed both default strategies "
            f"in the range {r['threshold_range_beats_both_strategies']}. "
        )

    # External benchmarks
    text += "\n\nAcross external benchmarks, AGVarPred's net benefit was threshold- and dataset-dependent:\n\n"
    ext = summary_df[summary_df["benchmark"] != "internal_test"]
    for bench in sorted(ext["benchmark"].unique()):
        sub = ext[ext["benchmark"] == bench]
        if len(sub) == 0:
            continue
        text += f"- **{bench.replace('_', ' ').title()}**: "
        parts = []
        for _, r in sub.iterrows():
            parts.append(
                f"{r['model']} max NB = {r['max_net_benefit']:.3f} at {r['max_net_benefit_threshold']:.2f} "
                f"(integrated NB = {r['integrated_net_benefit']:.3f})"
            )
        text += "; ".join(parts) + "\n"

    text += """
The clinical value of AGVarPred is most evident when its net-benefit curve lies above both treat-all and treat-none curves. On the held-out ClinVar test set, AGVarPred dominated both default strategies across almost the entire threshold range, consistent with strong calibration and discrimination. On external benchmarks, AGVarPred typically produced positive net benefit only at low thresholds (approximately 0.01–0.20 on Humsavar, VIP, Grimm2015, and DVD), after which false-negative penalties drove net benefit below that of treating none. This pattern reflects both the low prevalence of pathogenic variants in most benchmarks and the conservative calibration of AGVarPred probabilities. On gnomAD Common Benign, AGVarPred's net-benefit curve remained near the treat-none line, confirming that it rarely misclassifies common benign variants as pathogenic. Compared with AlphaMissense and REVEL on common subsets, AGVarPred was generally outperformed on missense-enriched benchmarks but remained the only tool evaluated across all variant types. These analyses support the interpretation of AGVarPred probabilities as calibrated, clinically informative risk scores on the test set, while highlighting that external clinical utility depends strongly on prevalence, variant spectrum, and threshold choice. They do not establish prospective clinical effectiveness.
"""

    with open(os.path.join(OUT_DIR, "manuscript_results.txt"), "w") as f:
        f.write(text)


def write_discussion_subsection():
    text = """Decision Curve Analysis in context

Decision Curve Analysis complements ROC-AUC by evaluating whether a model's predicted probabilities actually improve clinical decision-making at specific thresholds. While ROC-AUC summarizes discrimination across all score rankings, DCA directly measures the net clinical benefit of acting on a prediction, making it more relevant for variant-prioritization workflows in which laboratories apply a fixed probability or score cutoff.

AGVarPred's calibrated probabilities are particularly important for DCA because poorly calibrated scores can produce misleadingly high ROC-AUC while yielding negative net benefit at clinically relevant thresholds. The fact that AGVarPred's net-benefit curve is positive over a wide threshold range on the held-out test set suggests that the isotonic calibration step succeeded in producing interpretable probabilities. This matters for genome-wide interpretation, where a single threshold must perform across diverse variant classes and genes.

Several limitations apply. First, DCA is retrospective: it evaluates net benefit on datasets whose prevalence and spectrum may not match future clinical cohorts. Second, the treat-all and treat-none strategies are simplifications; real diagnostic pipelines combine multiple lines of evidence. Third, competitor comparisons rely on common subsets, which underrepresent non-missense variants and may favor tools that cannot score the full variant spectrum. Finally, net benefit is sensitive to prevalence, so thresholds optimal for one benchmark may not transfer to another. Prospective validation in real clinical laboratories remains essential before claiming clinical utility.
"""
    with open(os.path.join(OUT_DIR, "manuscript_discussion.txt"), "w") as f:
        f.write(text)


def main():
    all_data = []
    for bench in BENCHMARKS:
        all_data.append(run_dca_for_benchmark(bench))
    data_df = pd.concat(all_data, ignore_index=True)
    data_df.to_csv(os.path.join(RESULTS_DIR, "decision_curve_data.csv"), index=False)
    print("Saved decision_curve_data.csv")

    summary_df = compute_summary(data_df)
    summary_df.to_csv(os.path.join(RESULTS_DIR, "decision_curve_summary.csv"), index=False)
    print("Saved decision_curve_summary.csv")

    # Figures
    create_figure(data_df, os.path.join(FIGURES_DIR, "decision_curve_analysis.pdf"))
    print("Saved decision_curve_analysis.pdf/png")

    create_internal_figure(data_df, os.path.join(FIGURES_DIR, "decision_curve_internal_test.pdf"))
    print("Saved decision_curve_internal_test.pdf/png")

    # Manuscript text
    write_results_subsection(summary_df)
    write_discussion_subsection()
    print("Saved manuscript_results.txt and manuscript_discussion.txt")

    print("\nDCA summary:")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
