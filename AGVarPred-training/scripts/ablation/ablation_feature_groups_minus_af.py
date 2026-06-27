import pandas as pd
import numpy as np
import gc, re, joblib, os
from glob import glob

from sklearn.metrics import (
    roc_auc_score, f1_score, precision_score, recall_score,
    roc_curve, average_precision_score,
    brier_score_loss, confusion_matrix, matthews_corrcoef,
    precision_recall_curve
)
from sklearn.isotonic import IsotonicRegression
from sklearn.calibration import calibration_curve
from sklearn.preprocessing import RobustScaler
from sklearn.dummy import DummyClassifier

from lightgbm import LGBMClassifier
import matplotlib.pyplot as plt

# =========================
# OUTPUT DIR
# =========================

OUTDIR = "model_6_minus_af_output"
os.makedirs(OUTDIR, exist_ok=True)

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
        dfs.append(pd.read_parquet(f))

        if i % 500 == 0:
            print(f"{path}: {i}/{len(files)}")

    df = pd.concat(dfs, ignore_index=True)
    del dfs
    gc.collect()

    df.columns = [clean_name(c) for c in df.columns]

    return df

print("🔄 Loading datasets...")

train_df = load_dataset("final_dataset_parts_train")
cal_df   = load_dataset("final_dataset_parts_cal")
test_df  = load_dataset("final_dataset_parts_test")

# =========================
# FEATURES
# =========================

with open("feature_selection_output_nested/selected_features.txt") as f:
    selected_features = [clean_name(x.strip()) for x in f]

selected_features = [f for f in selected_features if f in train_df.columns]

print("Selected features:", len(selected_features))

# =========================
# FIXED BEST PARAMETERS
# (from regularized_final_model.log)
# =========================

BEST_PARAMS = {
    'n_estimators': 485,
    'learning_rate': 0.06525555612936938,
    'num_leaves': 73,
    'max_depth': 8,
    'min_child_samples': 91,
    'subsample': 0.6012817574242904,
    'colsample_bytree': 0.6756617006329491,
    'reg_alpha': 0.03731217662652025,
    'reg_lambda': 0.013979283607874253
}

SCALE_POS_WEIGHT = (train_df['label'] == 0).sum() / (train_df['label'] == 1).sum()

# =========================
# FEATURE GROUPS
# =========================

af_features = [f for f in selected_features if f.lower() in ['gnomad_af']]
model_6_features = [f for f in selected_features if f not in af_features]

print(f"AF features: {len(af_features)}")
print(f"Model 6 features (All minus AF): {len(model_6_features)}")

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
# MODEL 6: ALL FEATURES MINUS AF
# =========================

results = []

model_name = "Model_1_no_AF"
features = model_6_features

print(f"\n{'='*60}")
print(f"🚀 {model_name} — {len(features)} features")
print(f"{'='*60}")

if len(features) == 0:
    print("⚠️ No features. Skipping.")
else:
    X_train = train_df[features]
    y_train = train_df['label']
    X_cal = cal_df[features]
    y_cal = cal_df['label']
    X_test = test_df[features]
    y_test = test_df['label']

    # Train from scratch
    imputer = X_train.median()

    X_train_imp = X_train.fillna(imputer)
    X_cal_imp   = X_cal.fillna(imputer)
    X_test_imp  = X_test.fillna(imputer)

    scaler = RobustScaler()
    X_train_f = pd.DataFrame(
        scaler.fit_transform(X_train_imp),
        columns=X_train_imp.columns,
        index=X_train_imp.index
    )
    X_cal_f = pd.DataFrame(
        scaler.transform(X_cal_imp),
        columns=X_cal_imp.columns,
        index=X_cal_imp.index
    )
    X_test_f = pd.DataFrame(
        scaler.transform(X_test_imp),
        columns=X_test_imp.columns,
        index=X_test_imp.index
    )

    model = LGBMClassifier(
        **BEST_PARAMS,
        scale_pos_weight=SCALE_POS_WEIGHT,
        random_state=42,
        verbose=-1,
        n_jobs=-1
    )

    print("Training...")
    model.fit(X_train_f, y_train)

    # Calibration
    probs_cal = model.predict_proba(X_cal_f)[:, 1]
    iso = IsotonicRegression(out_of_bounds='clip')
    iso.fit(probs_cal, y_cal)

    # Threshold
    probs_cal_calibrated = iso.transform(probs_cal)
    thresholds_grid = np.arange(0.05, 0.95, 0.01)
    final_thresh = max(
        thresholds_grid,
        key=lambda t: f1_score(y_cal, (probs_cal_calibrated > t).astype(int))
    )

    print(f"Threshold tuned on CAL set: {final_thresh:.2f}")

    # =========================
    # BASELINE
    # =========================

    dummy = DummyClassifier(strategy="stratified")
    dummy.fit(X_train_f, y_train)
    baseline_probs = dummy.predict_proba(X_test_f)[:, 1]
    baseline_auc = roc_auc_score(y_test, baseline_probs)
    print(f"Baseline AUC: {baseline_auc:.4f}")

    # =========================
    # TEST
    # =========================

    probs_test_raw = model.predict_proba(X_test_f)[:, 1]
    probs_test = iso.transform(probs_test_raw)

    preds = (probs_test > final_thresh).astype(int)

    roc_auc = roc_auc_score(y_test, probs_test)
    auprc = average_precision_score(y_test, probs_test)
    f1 = f1_score(y_test, preds)
    precision = precision_score(y_test, preds)
    recall = recall_score(y_test, preds)
    mcc = matthews_corrcoef(y_test, preds)
    brier = brier_score_loss(y_test, probs_test)
    ece = compute_ece(y_test.values, probs_test)

    print("\n📊 TEST PERFORMANCE")
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
    # SAVE OUTPUTS
    # =========================

    safe_name = model_name

    pd.DataFrame({
        "y_true": y_test.values,
        "raw_prob": probs_test_raw,
        "cal_prob": probs_test
    }).to_csv(f"{OUTDIR}/{safe_name}_predictions.csv", index=False)

    np.savetxt(f"{OUTDIR}/{safe_name}_confusion_matrix.txt",
               confusion_matrix(y_test, preds), fmt='%d')

    np.savetxt(f"{OUTDIR}/{safe_name}_threshold.txt", [final_thresh])

    # ROC
    nfpr, ntpr, _ = roc_curve(y_test, probs_test)
    plt.figure()
    plt.plot(nfpr, ntpr)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC Curve — {safe_name}")
    plt.savefig(f"{OUTDIR}/{safe_name}_roc_curve.png")
    plt.close()

    # PR
    prec, rec, _ = precision_recall_curve(y_test, probs_test)
    plt.figure()
    plt.plot(rec, prec)
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(f"Precision-Recall Curve — {safe_name}")
    plt.savefig(f"{OUTDIR}/{safe_name}_pr_curve.png")
    plt.close()

    # Calibration
    prob_true, prob_pred = calibration_curve(y_test, probs_test, n_bins=10)
    plt.figure()
    plt.plot(prob_pred, prob_true, marker='o')
    plt.plot([0,1],[0,1],'--')
    plt.xlabel("Predicted Probability")
    plt.ylabel("Observed Frequency")
    plt.title(f"Calibration Curve — {safe_name}")
    plt.savefig(f"{OUTDIR}/{safe_name}_calibration_curve.png")
    plt.close()

    # Feature importance
    pd.DataFrame({
        "feature": features,
        "importance": model.feature_importances_
    }).sort_values("importance", ascending=False).to_csv(
        f"{OUTDIR}/{safe_name}_feature_importance.csv", index=False
    )

    joblib.dump({
        "model": model,
        "calibrator": iso,
        "imputer": imputer,
        "scaler": scaler,
        "features": features,
        "threshold": final_thresh
    }, f"{OUTDIR}/{safe_name}_pipeline.pkl")

    # Metrics
    pd.DataFrame({
        "metric": [
            "AUC","AUPRC","F1","Precision","Recall","MCC",
            "Brier","ECE","Baseline_AUC",
            "AUC_CI_low","AUC_CI_high"
        ],
        "value": [
            roc_auc, auprc, f1, precision, recall, mcc,
            brier, ece, baseline_auc,
            ci_low, ci_high
        ]
    }).to_csv(f"{OUTDIR}/{safe_name}_metrics.csv", index=False)

    # Accumulate
    results.append({
        "model": model_name,
        "n_features": len(features),
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
    del X_train, X_cal, X_test, model, iso, scaler
    gc.collect()

# =========================
# SUMMARY
# =========================

df_res = pd.DataFrame(results)
df_res.to_csv(f"{OUTDIR}/ablation_summary.csv", index=False)

print("\n📊 MODEL 6 SUMMARY")
print(df_res.to_string(index=False))

print(f"\n🎯 All results saved to: {OUTDIR}/")
