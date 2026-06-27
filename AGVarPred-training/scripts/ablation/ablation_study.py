import pandas as pd
import numpy as np
import joblib
from glob import glob
import gc
import os
import re

from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    confusion_matrix,
    roc_curve,
    precision_recall_curve,
    brier_score_loss,
)
from sklearn.calibration import calibration_curve
from sklearn.isotonic import IsotonicRegression
from sklearn.preprocessing import RobustScaler
from sklearn.dummy import DummyClassifier

from lightgbm import LGBMClassifier

import matplotlib.pyplot as plt

print("🔄 Loading data...")

# =========================
# OUTPUT DIR
# =========================

OUT_DIR = "ablation_output"
os.makedirs(OUT_DIR, exist_ok=True)

# =========================
# CLEAN NAME
# =========================

def clean_name(s):
    s = str(s)
    s = s.replace(":", "_").replace("-", "_").replace(" ", "_")
    s = re.sub(r'[^A-Za-z0-9_]', '', s)
    s = re.sub(r'_+', '_', s)
    return s

# =========================
# LOAD DATA
# =========================

def load_dataset(path):
    files = glob(f"{path}/*.parquet")
    dfs = []
    for i, f in enumerate(files):
        df = pd.read_parquet(f)
        df.columns = [clean_name(c) for c in df.columns]
        dfs.append(df)
        if i % 500 == 0:
            print(f"{path}: {i}/{len(files)}")
    df = pd.concat(dfs, ignore_index=True)
    del dfs
    gc.collect()
    return df

train_df = load_dataset("final_dataset_parts_train")
cal_df   = load_dataset("final_dataset_parts_cal")
test_df  = load_dataset("final_dataset_parts_test")

# =========================
# LOAD FINALIZED REGULARIZED LGBM PIPELINE
# =========================

pipeline = joblib.load("final_model_output_regularized/final_pipeline.pkl")

base_model = pipeline["model"]
selected_features = pipeline["features"]

print(f"✅ Loaded finalized regularized LGBM model")
print(f"Features in pipeline: {len(selected_features)}")

# =========================
# LOAD SHAP RANKING FROM REGULARIZED MODEL
# =========================

shap_df = pd.read_csv("final_model_output_regularized/feature_importance_shap.csv")

ranked_features = [
    clean_name(f)
    for f in shap_df["feature"].tolist()
]

ranked_features = [
    f for f in ranked_features
    if f in train_df.columns
]

print(f"✅ Using {len(ranked_features)} SHAP-ranked features")

# =========================
# ECE
# =========================

def compute_ece(y_true, y_prob, n_bins=10):
    bins = np.linspace(0, 1, n_bins + 1)
    binids = np.digitize(y_prob, bins) - 1
    ece = 0
    for i in range(n_bins):
        mask = binids == i
        if np.sum(mask) > 0:
            acc = np.mean(y_true[mask])
            conf = np.mean(y_prob[mask])
            ece += np.abs(acc - conf) * np.sum(mask) / len(y_true)
    return ece

# =========================
# ABLATION SETTINGS
# =========================

steps = [0, 5, 10, 15, 20, 25, 30, 35]
results = []

# =========================
# ABLATION LOOP
# =========================

for remove_n in steps:

    print(f"\n{'='*60}")
    print(f"🚀 ABLATION STEP: Removing top {remove_n} SHAP features")
    print(f"{'='*60}")

    remaining_features = ranked_features[remove_n:]
    if len(remaining_features) == 0:
        print("⚠️ No remaining features. Skipping.")
        continue

    print(f"Remaining features: {len(remaining_features)}")

    # =========================
    # PREPARE DATA FOR THIS STEP
    # =========================

    X_train = train_df[remaining_features]
    y_train = train_df["label"]

    X_cal = cal_df[remaining_features]
    y_cal = cal_df["label"]

    X_test = test_df[remaining_features]
    y_test = test_df["label"]

    # =========================
    # IMPUTATION + SCALING (full train, consistent with final model training)
    # =========================

    imputer_full = X_train.median()

    X_train_imp = X_train.fillna(imputer_full)
    X_cal_imp   = X_cal.fillna(imputer_full)
    X_test_imp  = X_test.fillna(imputer_full)

    scaler_full = RobustScaler()
    X_train_f = pd.DataFrame(
        scaler_full.fit_transform(X_train_imp),
        columns=X_train_imp.columns,
        index=X_train_imp.index
    )
    X_cal_f = pd.DataFrame(
        scaler_full.transform(X_cal_imp),
        columns=X_cal_imp.columns,
        index=X_cal_imp.index
    )
    X_test_f = pd.DataFrame(
        scaler_full.transform(X_test_imp),
        columns=X_test_imp.columns,
        index=X_test_imp.index
    )

    # =========================
    # BASELINE (dummy)
    # =========================

    dummy = DummyClassifier(strategy="stratified")
    dummy.fit(X_train_f, y_train)
    baseline_probs = dummy.predict_proba(X_test_f)[:, 1]
    baseline_auc = roc_auc_score(y_test, baseline_probs)
    print(f"Baseline AUC: {baseline_auc:.4f}")

    # =========================
    # FINAL MODEL — SAME HYPERPARAMS, FEWER FEATURES
    # =========================

    model = LGBMClassifier(**base_model.get_params())

    print(f"\nTraining LGBM with fixed hyperparams (remove {remove_n})...")
    model.fit(X_train_f, y_train)

    # =========================
    # CALIBRATION (retuned per step)
    # =========================

    probs_cal = model.predict_proba(X_cal_f)[:, 1]

    iso = IsotonicRegression(out_of_bounds='clip')
    iso.fit(probs_cal, y_cal)

    # =========================
    # THRESHOLD (retuned on cal set)
    # =========================

    probs_cal_calibrated = iso.transform(probs_cal)

    thresholds_grid = np.arange(0.05, 0.95, 0.01)
    final_thresh = max(
        thresholds_grid,
        key=lambda t: f1_score(y_cal, (probs_cal_calibrated > t).astype(int))
    )

    np.savetxt(f"{OUT_DIR}/threshold_remove_{remove_n}.txt", [final_thresh])
    print(f"Threshold tuned on CAL set: {final_thresh:.2f}")

    # =========================
    # TEST
    # =========================

    probs_test_raw = model.predict_proba(X_test_f)[:, 1]
    probs_test = iso.transform(probs_test_raw)

    preds = (probs_test > final_thresh).astype(int)

    roc_auc = roc_auc_score(y_test, probs_test)
    auprc   = average_precision_score(y_test, probs_test)
    f1      = f1_score(y_test, preds)
    precision = precision_score(y_test, preds)
    recall  = recall_score(y_test, preds)
    mcc     = matthews_corrcoef(y_test, preds)
    brier   = brier_score_loss(y_test, probs_test)
    ece     = compute_ece(y_test.values, probs_test)

    print(f"\n📊 TEST PERFORMANCE (remove {remove_n})")
    print(f"AUC:    {roc_auc:.4f}")
    print(f"AUPRC:  {auprc:.4f}")
    print(f"F1:     {f1:.4f}")

    # =========================
    # BOOTSTRAP CI
    # =========================

    rng = np.random.RandomState(42)
    auc_scores = []

    for _ in range(200):
        idx = rng.choice(len(y_test), len(y_test), replace=True)
        auc_scores.append(roc_auc_score(y_test.values[idx], probs_test[idx]))

    ci_low = np.percentile(auc_scores, 2.5)
    ci_high = np.percentile(auc_scores, 97.5)

    # =========================
    # SAVE PREDICTIONS
    # =========================

    pd.DataFrame({
        "y_true": y_test.values,
        "raw_prob": probs_test_raw,
        "cal_prob": probs_test,
    }).to_csv(
        f"{OUT_DIR}/predictions_remove_{remove_n}.csv",
        index=False
    )

    # =========================
    # CONFUSION MATRIX
    # =========================

    np.savetxt(
        f"{OUT_DIR}/confusion_matrix_remove_{remove_n}.txt",
        confusion_matrix(y_test, preds),
        fmt='%d'
    )

    # =========================
    # ROC CURVE
    # =========================

    fpr, tpr, _ = roc_curve(y_test, probs_test)

    plt.figure()
    plt.plot(fpr, tpr)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC Curve (Remove {remove_n})")
    plt.savefig(f"{OUT_DIR}/roc_remove_{remove_n}.png")
    plt.close()

    # =========================
    # PR CURVE
    # =========================

    prec_curve, rec_curve, _ = precision_recall_curve(y_test, probs_test)

    plt.figure()
    plt.plot(rec_curve, prec_curve)
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(f"PR Curve (Remove {remove_n})")
    plt.savefig(f"{OUT_DIR}/pr_remove_{remove_n}.png")
    plt.close()

    # =========================
    # CALIBRATION CURVE
    # =========================

    prob_true, prob_pred = calibration_curve(y_test, probs_test, n_bins=10)

    plt.figure()
    plt.plot(prob_pred, prob_true, marker='o')
    plt.plot([0, 1], [0, 1], '--')
    plt.xlabel("Predicted Probability")
    plt.ylabel("Observed Frequency")
    plt.title(f"Calibration Curve (Remove {remove_n})")
    plt.savefig(f"{OUT_DIR}/calibration_remove_{remove_n}.png")
    plt.close()

    # =========================
    # FEATURE IMPORTANCE
    # =========================

    importance_df = pd.DataFrame({
        "feature": remaining_features,
        "importance": model.feature_importances_,
    }).sort_values("importance", ascending=False)

    importance_df.to_csv(
        f"{OUT_DIR}/feature_importance_remove_{remove_n}.csv",
        index=False
    )

    # =========================
    # SAVE MODEL PIPELINE
    # =========================

    joblib.dump({
        "model": model,
        "calibrator": iso,
        "imputer": imputer_full,
        "scaler": scaler_full,
        "features": remaining_features,
        "threshold": final_thresh,
    }, f"{OUT_DIR}/pipeline_remove_{remove_n}.pkl")

    # =========================
    # SAVE METRICS
    # =========================

    pd.DataFrame({
        "metric": [
            "AUC", "AUPRC", "F1", "Precision",
            "Recall", "MCC", "Brier", "ECE",
            "Baseline_AUC",
            "AUC_CI_low", "AUC_CI_high",
        ],
        "value": [
            roc_auc, auprc, f1, precision,
            recall, mcc, brier, ece,
            baseline_auc,
            ci_low, ci_high,
        ],
    }).to_csv(
        f"{OUT_DIR}/metrics_remove_{remove_n}.csv",
        index=False
    )

    # =========================
    # ACCUMULATE RESULTS
    # =========================

    results.append({
        "features_removed": remove_n,
        "features_remaining": len(remaining_features),
        "AUC": roc_auc,
        "AUPRC": auprc,
        "F1": f1,
        "Precision": precision,
        "Recall": recall,
        "MCC": mcc,
        "Brier": brier,
        "ECE": ece,
        "Baseline_AUC": baseline_auc,
        "AUC_CI_low": ci_low,
        "AUC_CI_high": ci_high,
    })

    # cleanup
    del X_train, X_cal, X_test, model, iso, scaler_full
    gc.collect()

# =========================
# FINAL RESULTS TABLE
# =========================

df_res = pd.DataFrame(results)

df_res.to_csv(
    f"{OUT_DIR}/ablation_results.csv",
    index=False
)

# =========================
# MAIN ABLATION PLOT
# =========================

plt.figure(figsize=(7, 5))

plt.plot(
    df_res["features_remaining"],
    df_res["AUC"],
    marker='o',
    label="AUC"
)

plt.plot(
    df_res["features_remaining"],
    df_res["AUPRC"],
    marker='s',
    label="AUPRC"
)

plt.xlabel("Number of Features Remaining")
plt.ylabel("Performance")
plt.title("SHAP-based Ablation Study (Regularized LGBM)")
plt.legend()
plt.grid()

plt.savefig(
    f"{OUT_DIR}/ablation_curve.png",
    dpi=200
)

plt.close()

print("\n📊 Ablation results:")
print(df_res)

print("\n🎯 Ablation study complete")
